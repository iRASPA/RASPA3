module;

module interactions_polarization_derivatives;

import std;

import int3;
import double3;
import double3x3;
import atom;
import molecule;
import component;
import forcefield;
import framework;
import simulationbox;
import units;
import system;
import potential_pair_derivatives;
import potential_pair_coulomb;

namespace
{
// Row-major 3x3 helpers.
using Mat3 = std::array<double, 9>;

Mat3 identityScaled(double s)
{
  return Mat3{s, 0.0, 0.0, 0.0, s, 0.0, 0.0, 0.0, s};
}

Mat3 outer(const double3& u, const double3& v)
{
  const std::array<double, 3> a = {u.x, u.y, u.z};
  const std::array<double, 3> b = {v.x, v.y, v.z};
  Mat3 m{};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) m[i * 3 + j] = a[i] * b[j];
  return m;
}

Mat3 matmul(const Mat3& a, const Mat3& b)
{
  Mat3 c{};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j)
    {
      double sum = 0.0;
      for (std::size_t k = 0; k < 3; ++k) sum += a[i * 3 + k] * b[k * 3 + j];
      c[i * 3 + j] = sum;
    }
  return c;
}

// Field-point separation from a source, d = r_field - r_source (minimum image).
struct RealSource
{
  std::size_t globalIndex{};
  double3 d{};
  int3 latticeVector{};  // image cell of the source relative to the (home-cell) field point
  double3 offset{};      // sigma = r_source - COM_source (nonzero only for rigid-molecule sources)
  double charge{};       // scalingCoulomb * charge of the source
  double firstDerivativeFactor{};
  double secondDerivativeFactor{};
  double thirdDerivativeFactor{};
};

// One reciprocal wave vector contribution needed for the on-site (A,A) blocks.
struct ReciprocalTerm
{
  double3 k{};
  double weight{};  // 2 * temp
  double f{};       // Im[e^{i k.r_A}] Re[S] - Re[e^{i k.r_A}] Im[S]
  double hValue{};  // Re[e^{i k.r_A}] Re[S] + Im[e^{i k.r_A}] Im[S] = d f / d(k.r_A)
};

// Unit Cartesian axis vector.
double3 unitAxis(std::size_t axis)
{
  double3 e(0.0, 0.0, 0.0);
  e[axis] = 1.0;
  return e;
}

// Reciprocal-space data that does not depend on the field point: a surviving wave vector, its Gaussian
// weight, and the fixed-framework structure factor S(k) = sum_f q_f exp(i k.r_f).
struct PrecomputedReciprocal
{
  double3 k{};
  double weight{};  // 2 * temp
  double structureReal{};
  double structureImaginary{};
};

// Enumerates the Ewald wave vectors and precomputes S(k) once; the per-field-point loops then only
// evaluate their own phase against the cached structure factor.
std::vector<PrecomputedReciprocal> precomputeFixedFrameworkReciprocal(const ForceField& forceField,
                                                                      const SimulationBox& box,
                                                                      std::span<const Atom> frameworkAtoms)
{
  const double alphaSquared = forceField.EwaldAlpha * forceField.EwaldAlpha;
  const std::size_t reciprocalIntegerCutOffSquared = forceField.reciprocalIntegerCutOffSquared;
  const double reciprocalCutOffSquared = forceField.reciprocalCutOffSquared;
  const double3x3 inverseCell = box.inverseCell;
  const double3 aStar(inverseCell.ax, inverseCell.bx, inverseCell.cx);
  const double3 bStar(inverseCell.ay, inverseCell.by, inverseCell.cy);
  const double3 cStar(inverseCell.az, inverseCell.bz, inverseCell.cz);
  constexpr double twoPi = 2.0 * std::numbers::pi;
  const double prefactor = Units::CoulombicConversionFactor * (twoPi / box.volume);
  const int kxMax = static_cast<int>(forceField.numberOfWaveVectors.x);
  const int kyMax = static_cast<int>(forceField.numberOfWaveVectors.y);
  const int kzMax = static_cast<int>(forceField.numberOfWaveVectors.z);

  std::vector<PrecomputedReciprocal> precomputedReciprocal;
  for (int kx = 0; kx <= kxMax; ++kx)
  {
    const double factor = (kx == 0) ? prefactor : 2.0 * prefactor;
    for (int ky = -kyMax; ky <= kyMax; ++ky)
    {
      for (int kz = -kzMax; kz <= kzMax; ++kz)
      {
        const std::size_t integerSquared = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if (integerSquared == 0uz || integerSquared > reciprocalIntegerCutOffSquared) continue;
        const double3 kVector = twoPi * (static_cast<double>(kx) * aStar + static_cast<double>(ky) * bStar +
                                         static_cast<double>(kz) * cStar);
        const double kSquared = double3::dot(kVector, kVector);
        if (kSquared >= reciprocalCutOffSquared) continue;
        const double temp = factor * std::exp(-0.25 * kSquared / alphaSquared) / kSquared;

        double structureReal = 0.0;
        double structureImaginary = 0.0;
        for (std::size_t f = 0; f < frameworkAtoms.size(); ++f)
        {
          const double charge = frameworkAtoms[f].scalingCoulomb * frameworkAtoms[f].charge;
          const double argument = double3::dot(kVector, frameworkAtoms[f].position);
          structureReal += charge * std::cos(argument);
          structureImaginary += charge * std::sin(argument);
        }

        precomputedReciprocal.push_back(PrecomputedReciprocal{
            .k = kVector, .weight = 2.0 * temp, .structureReal = structureReal, .structureImaginary = structureImaginary});
      }
    }
  }
  return precomputedReciprocal;
}
}  // namespace

Interactions::PolarizationDerivatives Interactions::computePolarizationDerivatives(
    const System& system, std::span<const std::uint8_t> movable, std::span<const double3x3> strainBases,
    bool computeHessian, bool molecularCenterOfMassStrain)
{
  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();

  const std::size_t numberOfFrameworkAtoms = frameworkAtoms.size();
  const std::size_t numberOfMoleculeAtoms = moleculeAtoms.size();
  const std::size_t numberOfAtoms = numberOfFrameworkAtoms + numberOfMoleculeAtoms;

  const std::size_t numberOfStrainBases = strainBases.size();
  const bool computeStrain = numberOfStrainBases > 0;
  const bool wantHessian = computeHessian;
  const bool wantGradient = computeHessian;
  const bool wantStrainSecond = computeStrain && computeHessian;

  PolarizationDerivatives result;
  result.numberOfFrameworkAtoms = numberOfFrameworkAtoms;
  result.numberOfMoleculeAtoms = numberOfMoleculeAtoms;
  result.gradient.assign(numberOfAtoms, double3(0.0, 0.0, 0.0));
  result.numberOfStrainBases = numberOfStrainBases;
  if (computeStrain)
  {
    result.cellGradient.assign(numberOfStrainBases, 0.0);
    if (wantStrainSecond)
    {
      result.positionStrain.assign(numberOfAtoms * numberOfStrainBases, double3(0.0, 0.0, 0.0));
      result.strainStrain.assign(numberOfStrainBases * numberOfStrainBases, 0.0);
    }
  }

  if (movable.size() != numberOfAtoms)
  {
    throw std::runtime_error("computePolarizationDerivatives: movable mask size does not match the number of atoms");
  }
  if (!forceField.computePolarization || !forceField.useCharge || numberOfMoleculeAtoms == 0) return result;

  // Under cell strain, molecule atoms deform via their center of mass (y_i = F COM + sigma_i, with the
  // body offset sigma_i strain-independent). Precompute sigma_i so the affine separation-arm can be
  // corrected from the atom-atom vector to the COM-COM vector, and so the reciprocal structure-factor
  // phase (which is only invariant for the COM part) carries the extra body-offset terms.
  // Minimization uses stored rigid-molecule COMs; molecular-pressure sampling recomputes a mass-weighted
  // COM from current atom positions for every molecule (matching the COM-scaling volume move).
  std::vector<double3> moleculeAtomOffset;
  if (computeStrain)
  {
    moleculeAtomOffset.assign(numberOfMoleculeAtoms, double3(0.0, 0.0, 0.0));
    for (const Molecule& molecule : system.moleculeData)
    {
      if (molecularCenterOfMassStrain)
      {
        double totalMass = 0.0;
        double3 com(0.0, 0.0, 0.0);
        for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
        {
          const std::size_t index = molecule.atomIndex + localAtom;
          const double mass = forceField.pseudoAtoms[static_cast<std::size_t>(moleculeAtoms[index].type)].mass;
          com += mass * moleculeAtoms[index].position;
          totalMass += mass;
        }
        com = com / totalMass;
        for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
        {
          const std::size_t index = molecule.atomIndex + localAtom;
          moleculeAtomOffset[index] = moleculeAtoms[index].position - com;
        }
      }
      else
      {
        if (!system.components[molecule.componentId].rigid) continue;
        for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
        {
          const std::size_t index = molecule.atomIndex + localAtom;
          moleculeAtomOffset[index] = moleculeAtoms[index].position - molecule.centerOfMassPosition;
        }
      }
    }
  }

  // Precompute the (symmetric) strain generators B_a, their traces, and the second-order log-strain
  // generators C_ab = 1/2 (B_a B_b + B_b B_a) = d^2 F / ds_a ds_b (matching cellStrainSecondDerivative).
  std::vector<double3x3> strainGenerator(strainBases.begin(), strainBases.end());
  std::vector<double> strainTrace(numberOfStrainBases, 0.0);
  std::vector<double3x3> strainSecond(numberOfStrainBases * numberOfStrainBases);
  for (std::size_t a = 0; a < numberOfStrainBases; ++a)
  {
    strainTrace[a] = strainGenerator[a].trace();
    for (std::size_t b = 0; b < numberOfStrainBases; ++b)
    {
      strainSecond[a * numberOfStrainBases + b] =
          0.5 * (strainGenerator[a] * strainGenerator[b] + strainGenerator[b] * strainGenerator[a]);
    }
  }

  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const bool includeInter = !forceField.omitInterInteractions && !forceField.omitInterPolarization;

  // Fixed-framework reciprocal structure factor sources (rigid framework only, matching the model's
  // computeEwaldFourierElectricField which uses the fixed framework structure factor).
  const bool useReciprocal = forceField.usesEwaldFourier() && system.framework.has_value() &&
                             system.framework->rigid && numberOfFrameworkAtoms > 0;

  const double alpha = forceField.EwaldAlpha;
  const double alphaSquared = alpha * alpha;

  // Accumulate a Cartesian block into both the minimum-image-folded Hessian (Gamma / normal modes) and
  // the image-resolved Hessian (phonon force constants). The field point sits in the home cell, so the
  // lattice vector is the source image relative to it: block(field, source) carries +R, block(source,
  // field) carries -R, and on-site / same-image blocks carry R = 0.
  auto addBlock = [&](std::size_t i, std::size_t j, const int3& latticeVector, const Mat3& contribution)
  {
    Mat3& destination = result.hessianBlocks[{i, j}];
    for (std::size_t index = 0; index < 9; ++index) destination[index] += contribution[index];

    Mat3& imageDestination = result.imageHessianBlocks[latticeVector][{i, j}];
    for (std::size_t index = 0; index < 9; ++index) imageDestination[index] += contribution[index];
  };

  const std::vector<PrecomputedReciprocal> precomputedReciprocal =
      useReciprocal ? precomputeFixedFrameworkReciprocal(forceField, box, frameworkAtoms)
                    : std::vector<PrecomputedReciprocal>{};

  std::vector<RealSource> sources;
  std::vector<ReciprocalTerm> reciprocalTerms;

  for (std::size_t a = 0; a < numberOfMoleculeAtoms; ++a)
  {
    const std::size_t globalA = numberOfFrameworkAtoms + a;
    const std::size_t typeA = static_cast<std::size_t>(moleculeAtoms[a].type);
    const double alphaA = forceField.pseudoAtoms[typeA].polarizability / Units::CoulombicConversionFactor;
    if (alphaA == 0.0) continue;
    if (!movable[globalA])
    {
      throw std::runtime_error("computePolarizationDerivatives: a polarizable atom must be movable");
    }

    const double3 positionA = moleculeAtoms[a].position;
    const std::size_t moleculeIdA = static_cast<std::size_t>(moleculeAtoms[a].moleculeId);
    // Rigid body offset of the field point (zero for flexible molecules); used by the strain derivatives.
    const double3 sigmaA = computeStrain ? moleculeAtomOffset[a] : double3(0.0, 0.0, 0.0);

    sources.clear();
    reciprocalTerms.clear();
    double3 field(0.0, 0.0, 0.0);
    Mat3 fieldGradientAA{};  // dE_A/dr_A (row = field component, column = displacement component)

    const auto gatherRealSource = [&](std::size_t globalIndex, const double3& sourcePosition, double sourceCharge)
    {
      const double3 rawSeparation = positionA - sourcePosition;
      double3 d = box.applyPeriodicBoundaryConditions(rawSeparation);
      const double rr = double3::dot(d, d);
      if (rr >= cutOffChargeSquared) return;
      const double r = std::sqrt(rr);
      // The third-derivative factor is only consumed by the Hessian / second-strain blocks; the
      // molecular-pressure path (computeHessian == false) needs derivative orders <= 2 only.
      double firstDerivativeFactor;
      double secondDerivativeFactor;
      double thirdDerivativeFactor;
      if (wantHessian)
      {
        const Potentials::PairDerivatives<3> factors =
            Potentials::potentialCoulomb<3>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        firstDerivativeFactor = factors.firstDerivativeFactor;
        secondDerivativeFactor = factors.secondDerivativeFactor;
        thirdDerivativeFactor = factors.thirdDerivativeFactor;
      }
      else
      {
        const Potentials::PairDerivatives<2> factors =
            Potentials::potentialCoulomb<2>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        firstDerivativeFactor = factors.firstDerivativeFactor;
        secondDerivativeFactor = factors.secondDerivativeFactor;
        thirdDerivativeFactor = 0.0;
      }

      // Lattice vector of the source image nearest the (home-cell) field point: R = round((raw - d) / cell).
      const double3 latticeCartesian = rawSeparation - d;
      const double3 fractional = box.inverseCell * latticeCartesian;
      const int3 latticeVector(static_cast<std::int32_t>(std::llround(fractional.x)),
                               static_cast<std::int32_t>(std::llround(fractional.y)),
                               static_cast<std::int32_t>(std::llround(fractional.z)));

      // E_A += -q_source * firstDerivativeFactor * d.
      field -= sourceCharge * firstDerivativeFactor * d;

      // T = q_source * (firstDerivativeFactor * I + secondDerivativeFactor * d (x) d) = q_source * H_phi.
      if (wantHessian || wantGradient)
      {
        Mat3 tensor = identityScaled(sourceCharge * firstDerivativeFactor);
        const Mat3 dd = outer(d, d);
        for (std::size_t index = 0; index < 9; ++index)
          tensor[index] += sourceCharge * secondDerivativeFactor * dd[index];

        // dE_A/dr_A gets -T from every source.
        for (std::size_t index = 0; index < 9; ++index) fieldGradientAA[index] -= tensor[index];
      }

      const double3 sourceOffset = (computeStrain && globalIndex >= numberOfFrameworkAtoms)
                                       ? moleculeAtomOffset[globalIndex - numberOfFrameworkAtoms]
                                       : double3(0.0, 0.0, 0.0);

      sources.push_back(RealSource{.globalIndex = globalIndex,
                                   .d = d,
                                   .latticeVector = latticeVector,
                                   .offset = sourceOffset,
                                   .charge = sourceCharge,
                                   .firstDerivativeFactor = firstDerivativeFactor,
                                   .secondDerivativeFactor = secondDerivativeFactor,
                                   .thirdDerivativeFactor = thirdDerivativeFactor});
    };

    if (includeInter)
    {
      for (std::size_t b = 0; b < numberOfMoleculeAtoms; ++b)
      {
        if (b == a) continue;
        if (static_cast<std::size_t>(moleculeAtoms[b].moleculeId) == moleculeIdA) continue;
        gatherRealSource(numberOfFrameworkAtoms + b, moleculeAtoms[b].position,
                         moleculeAtoms[b].scalingCoulomb * moleculeAtoms[b].charge);
      }
    }
    for (std::size_t f = 0; f < numberOfFrameworkAtoms; ++f)
    {
      gatherRealSource(f, frameworkAtoms[f].position, frameworkAtoms[f].scalingCoulomb * frameworkAtoms[f].charge);
    }

    if (useReciprocal)
    {
      for (const PrecomputedReciprocal& entry : precomputedReciprocal)
      {
        const double argumentA = double3::dot(entry.k, positionA);
        const double cosA = std::cos(argumentA);
        const double sinA = std::sin(argumentA);
        const double fValue = sinA * entry.structureReal - cosA * entry.structureImaginary;
        const double hValue = cosA * entry.structureReal + sinA * entry.structureImaginary;

        field += entry.weight * fValue * entry.k;
        if (wantHessian || wantGradient)
        {
          const Mat3 kk = outer(entry.k, entry.k);
          for (std::size_t index = 0; index < 9; ++index) fieldGradientAA[index] += entry.weight * hValue * kk[index];
        }

        reciprocalTerms.push_back(ReciprocalTerm{.k = entry.k, .weight = entry.weight, .f = fValue, .hValue = hValue});
      }
    }

    result.energy -= 0.5 * alphaA * double3::dot(field, field);

    if (wantHessian || wantGradient)
    {
      // Term A on-site: -alpha_A (dE_A/dr_A)^T (dE_A/dr_A). fieldGradientAA is symmetric.
      Mat3 blockAA{};
      {
        const Mat3 product = matmul(fieldGradientAA, fieldGradientAA);
        for (std::size_t index = 0; index < 9; ++index) blockAA[index] += -alphaA * product[index];
      }

      // Reciprocal Term B (on-site): +alpha_A sum_k weight f (E.k) (k (x) k).
      for (const ReciprocalTerm& term : reciprocalTerms)
      {
        const double eDotK = double3::dot(field, term.k);
        const Mat3 kk = outer(term.k, term.k);
        const double scalar = alphaA * term.weight * term.f * eDotK;
        for (std::size_t index = 0; index < 9; ++index) blockAA[index] += scalar * kk[index];
      }

      // gradient contribution for the field point itself: -alpha_A E_A^T (dE_A/dr_A).
      {
        double3 g(0.0, 0.0, 0.0);
        for (std::size_t column = 0; column < 3; ++column)
        {
          const double value = field.x * fieldGradientAA[0 * 3 + column] + field.y * fieldGradientAA[1 * 3 + column] +
                               field.z * fieldGradientAA[2 * 3 + column];
          (&g.x)[column] = -alphaA * value;
        }
        result.gradient[globalA] += g;
      }

      // Per-source contributions (Term B on-site over all sources; cross/off-diagonal only for movable sources).
      const double3 fieldVector = field;
      for (const RealSource& source : sources)
      {
        // T = q (f1 I + f2 d(x)d); D_A(source) = +T.
        Mat3 tensor = identityScaled(source.charge * source.firstDerivativeFactor);
        const Mat3 dd = outer(source.d, source.d);
        for (std::size_t index = 0; index < 9; ++index)
          tensor[index] += source.charge * source.secondDerivativeFactor * dd[index];

        // W(a,b) = -q [ f2 (E_a d_b + E_b d_a + delta_ab (E.d)) + f3 (E.d) d_a d_b ]  (field-Hessian contracted E).
        const double eDotD = double3::dot(fieldVector, source.d);
        const std::array<double, 3> e = {fieldVector.x, fieldVector.y, fieldVector.z};
        const std::array<double, 3> d = {source.d.x, source.d.y, source.d.z};
        Mat3 w{};
        for (std::size_t i = 0; i < 3; ++i)
          for (std::size_t j = 0; j < 3; ++j)
          {
            double value = source.secondDerivativeFactor * (e[i] * d[j] + e[j] * d[i] + (i == j ? eDotD : 0.0));
            value += source.thirdDerivativeFactor * eDotD * d[i] * d[j];
            w[i * 3 + j] = -source.charge * value;
          }

        // Term B on-site (all sources contribute to d^2E_A/dr_A dr_A).
        for (std::size_t index = 0; index < 9; ++index) blockAA[index] += -alphaA * w[index];

        if (!movable[source.globalIndex]) continue;

        // Term A cross: block(A,src) = -alpha_A D_AA^T T ; block(src,A) transpose.
        const Mat3 crossAsrcTermA = matmul(fieldGradientAA, tensor);
        Mat3 blockAsrc{};
        Mat3 blockSrcA{};
        for (std::size_t index = 0; index < 9; ++index)
        {
          blockAsrc[index] += -alphaA * crossAsrcTermA[index];
        }
        const Mat3 crossSrcATermA = matmul(tensor, fieldGradientAA);
        for (std::size_t index = 0; index < 9; ++index)
        {
          blockSrcA[index] += -alphaA * crossSrcATermA[index];
        }

        // Term B cross (+alpha_A W on both off-diagonal blocks) and on-site source block (-alpha_A W).
        for (std::size_t index = 0; index < 9; ++index)
        {
          blockAsrc[index] += alphaA * w[index];
          blockSrcA[index] += alphaA * w[index];
        }
        // block(field, source) sits at +R (source image), block(source, field) at -R; the source self-block
        // couples the same image to itself, so it stays in the home image R = 0.
        addBlock(globalA, source.globalIndex, source.latticeVector, blockAsrc);
        addBlock(source.globalIndex, globalA, -source.latticeVector, blockSrcA);

        Mat3 blockSrcSrcTermB{};
        for (std::size_t index = 0; index < 9; ++index) blockSrcSrcTermB[index] += -alphaA * w[index];
        addBlock(source.globalIndex, source.globalIndex, int3(0, 0, 0), blockSrcSrcTermB);

        // gradient for the source atom: -alpha_A E_A^T T.
        double3 g(0.0, 0.0, 0.0);
        for (std::size_t column = 0; column < 3; ++column)
        {
          const double value =
              fieldVector.x * tensor[0 * 3 + column] + fieldVector.y * tensor[1 * 3 + column] +
              fieldVector.z * tensor[2 * 3 + column];
          (&g.x)[column] = -alphaA * value;
        }
        result.gradient[source.globalIndex] += g;
      }

      addBlock(globalA, globalA, int3(0, 0, 0), blockAA);

      // Term A three-center: block(B,C) += -alpha_A T_AB^T T_AC over movable source pairs (including B == C).
      for (std::size_t bIndex = 0; bIndex < sources.size(); ++bIndex)
      {
        if (!movable[sources[bIndex].globalIndex]) continue;
        Mat3 tensorB = identityScaled(sources[bIndex].charge * sources[bIndex].firstDerivativeFactor);
        const Mat3 ddB = outer(sources[bIndex].d, sources[bIndex].d);
        for (std::size_t index = 0; index < 9; ++index)
          tensorB[index] += sources[bIndex].charge * sources[bIndex].secondDerivativeFactor * ddB[index];

        for (std::size_t cIndex = 0; cIndex < sources.size(); ++cIndex)
        {
          if (!movable[sources[cIndex].globalIndex]) continue;
          Mat3 tensorC = identityScaled(sources[cIndex].charge * sources[cIndex].firstDerivativeFactor);
          const Mat3 ddC = outer(sources[cIndex].d, sources[cIndex].d);
          for (std::size_t index = 0; index < 9; ++index)
            tensorC[index] += sources[cIndex].charge * sources[cIndex].secondDerivativeFactor * ddC[index];

          const Mat3 product = matmul(tensorB, tensorC);
          Mat3 contribution{};
          for (std::size_t index = 0; index < 9; ++index) contribution[index] = -alphaA * product[index];
          // Both sources are referenced to the same home-cell field point, so their relative image is
          // R_C - R_B (the field-point cell cancels); B == C then lands in the home image R = 0.
          const int3 latticeVector = sources[cIndex].latticeVector - sources[bIndex].latticeVector;
          addBlock(sources[bIndex].globalIndex, sources[cIndex].globalIndex, latticeVector, contribution);
        }
      }
    }

    if (!computeStrain) continue;

    // Cell-strain derivatives via homogeneous deformation. All derivatives follow from
    //   U_A = -1/2 alpha_A |E_A|^2,   dU/dm = -alpha_A E_A . dE_A/dm,
    //   d^2U/dm dm' = -alpha_A [ (dE_A/dm) . (dE_A/dm') + E_A . (d^2 E_A/dm dm') ],
    // where a strain mode s_a deforms the cell as F = exp(sum_a s_a B_a). Rigid-molecule atoms move via
    // their center of mass, so a separation d = r_A - r_s deforms with the COM-COM arm
    //   delta = d - sigma_A + sigma_s      (sigma = r - COM, zero for flexible/framework atoms),
    //   dd/ds_a = B_a delta,   d^2 d/ds_a ds_b = C_ab delta   with C_ab = 1/2 (B_a B_b + B_b B_a).
    // Reciprocal k -> F^{-T} k and volume -> det F as usual, but for a rigid field point the phase
    // arg_A = k.r_A is no longer strain invariant: only the COM part deforms, giving
    //   d(arg_A)/ds_a = -(B_a k).sigma_A,   d^2(arg_A)/ds_a ds_b = (C_ab k).sigma_A.
    //
    // The position-strain block below is the *pure* force change d(grad_j)/ds_a (atoms move via their COM
    // under strain). The kinematic term Sum_j grad_j . d^2 y_j/dq ds_a (nonzero only for translational
    // degrees of freedom, i.e. B_a . grad projected onto the center-of-mass / Cartesian DOF) is added by
    // the caller so that rigid-orientation DOFs, which have no such term, stay correct.

    // (de/dd)(v) = -T_s v.
    const auto firstFieldApply = [](const RealSource& s, const double3& v) -> double3
    {
      return -s.charge * (s.firstDerivativeFactor * v + s.secondDerivativeFactor * double3::dot(s.d, v) * s.d);
    };
    // (d^2 e/dd^2)(u, v).
    const auto secondFieldApply = [](const RealSource& s, const double3& u, const double3& v) -> double3
    {
      const double du = double3::dot(s.d, u);
      const double dv = double3::dot(s.d, v);
      const double uv = double3::dot(u, v);
      return -s.charge *
             (s.secondDerivativeFactor * (dv * u + du * v + uv * s.d) + s.thirdDerivativeFactor * du * dv * s.d);
    };

    // First strain derivative of the field, dE_A/ds_a, and the stress dU/ds_a.
    std::vector<double3> dFieldStrain(numberOfStrainBases, double3(0.0, 0.0, 0.0));
    for (std::size_t a = 0; a < numberOfStrainBases; ++a)
    {
      double3 accumulated(0.0, 0.0, 0.0);
      for (const RealSource& source : sources)
      {
        // dd/ds_a = B_a delta; dE contribution = (de/dd)(B_a delta).
        const double3 delta = source.d - sigmaA + source.offset;
        accumulated += firstFieldApply(source, strainGenerator[a] * delta);
      }
      for (const ReciprocalTerm& term : reciprocalTerms)
      {
        const double kSquared = double3::dot(term.k, term.k);
        const double3 Bk = strainGenerator[a] * term.k;
        const double kBk = double3::dot(term.k, Bk);
        const double beta = -strainTrace[a] + (0.25 / alphaSquared + 1.0 / kSquared) * (2.0 * kBk);
        const double phaseVel = -double3::dot(Bk, sigmaA);  // d(arg_A)/ds_a
        accumulated += term.weight * ((beta * term.f + term.hValue * phaseVel) * term.k - term.f * Bk);
      }
      dFieldStrain[a] = accumulated;
      result.cellGradient[a] += -alphaA * double3::dot(field, accumulated);
    }

    if (!wantStrainSecond) continue;

    // Strain-strain block, d^2U/ds_a ds_b.
    for (std::size_t a = 0; a < numberOfStrainBases; ++a)
    {
      for (std::size_t b = 0; b < numberOfStrainBases; ++b)
      {
        const double3x3& secondGenerator = strainSecond[a * numberOfStrainBases + b];
        double3 secondField(0.0, 0.0, 0.0);
        for (const RealSource& source : sources)
        {
          const double3 delta = source.d - sigmaA + source.offset;
          const double3 Bad = strainGenerator[a] * delta;
          const double3 Bbd = strainGenerator[b] * delta;
          // (d^2 e/dd^2)(B_a delta, B_b delta) + (de/dd)(C_ab delta).
          secondField += secondFieldApply(source, Bad, Bbd);
          secondField += firstFieldApply(source, secondGenerator * delta);
        }
        for (const ReciprocalTerm& term : reciprocalTerms)
        {
          const double kSquared = double3::dot(term.k, term.k);
          const double phi = 0.25 / alphaSquared + 1.0 / kSquared;
          const double3 Bak = strainGenerator[a] * term.k;
          const double3 Bbk = strainGenerator[b] * term.k;
          const double3 Cabk = secondGenerator * term.k;
          const double kBak = double3::dot(term.k, Bak);
          const double kBbk = double3::dot(term.k, Bbk);
          const double kCabk = double3::dot(term.k, Cabk);
          const double betaA = -strainTrace[a] + 2.0 * phi * kBak;
          const double betaB = -strainTrace[b] + 2.0 * phi * kBbk;
          const double dBetaAB = 4.0 * kBak * kBbk / (kSquared * kSquared) - 4.0 * phi * kCabk;
          const double weightSecond = betaA * betaB + dBetaAB;  // d^2 weight / ds_a ds_b / weight
          // Phase (rigid field point only; sigma_A = 0 otherwise). f'' = -f, so f_ab picks up -f phaseVel^2.
          const double phaseVelA = -double3::dot(Bak, sigmaA);
          const double phaseVelB = -double3::dot(Bbk, sigmaA);
          const double phaseAccAB = double3::dot(Cabk, sigmaA);
          const double fa = term.hValue * phaseVelA;
          const double fb = term.hValue * phaseVelB;
          const double fab = -term.f * phaseVelA * phaseVelB + term.hValue * phaseAccAB;
          const double3 ga = -Bak;   // dg/ds_a
          const double3 gb = -Bbk;   // dg/ds_b
          const double3 gab = Cabk;  // d^2 g/ds_a ds_b
          secondField += term.weight * (weightSecond * term.f * term.k + fab * term.k + term.f * gab +
                                        betaA * fb * term.k + betaB * fa * term.k + betaA * term.f * gb +
                                        betaB * term.f * ga + fa * gb + fb * ga);
        }
        result.strainStrain[a * numberOfStrainBases + b] +=
            -alphaA * (double3::dot(dFieldStrain[a], dFieldStrain[b]) + double3::dot(field, secondField));
      }
    }

    // Position-strain block for the field point A (pure force change): d(grad_A)/ds_a per axis c.
    for (std::size_t a = 0; a < numberOfStrainBases; ++a)
    {
      double3 positionStrainA(0.0, 0.0, 0.0);
      for (std::size_t c = 0; c < 3; ++c)
      {
        const double3 axis = unitAxis(c);
        // dE_A/dr_{A,c} is column c of dE_A/dr_A (includes the reciprocal on-site term).
        const double3 dFieldPositionA(fieldGradientAA[0 * 3 + c], fieldGradientAA[1 * 3 + c], fieldGradientAA[2 * 3 + c]);
        double3 secondField(0.0, 0.0, 0.0);
        for (const RealSource& source : sources)
        {
          // d(dE_A/dr_{A,c})/ds_a = (d^2 e/dd^2)(e_c, B_a delta)  (D_AA = -sum_s T_s).
          const double3 delta = source.d - sigmaA + source.offset;
          secondField += secondFieldApply(source, axis, strainGenerator[a] * delta);
        }
        for (const ReciprocalTerm& term : reciprocalTerms)
        {
          const double kSquared = double3::dot(term.k, term.k);
          const double3 Bak = strainGenerator[a] * term.k;
          const double kBak = double3::dot(term.k, Bak);
          const double betaA = -strainTrace[a] + (0.25 / alphaSquared + 1.0 / kSquared) * 2.0 * kBak;
          const double phaseVelA = -double3::dot(Bak, sigmaA);
          // d(D_AA^recip column c)/ds_a = weight[ (beta hValue - f phaseVel) k_c k - hValue (B_a k)_c k
          //                                       - hValue k_c (B_a k) ].
          secondField += term.weight * ((betaA * term.hValue - term.f * phaseVelA) * term.k[c] * term.k -
                                        term.hValue * Bak[c] * term.k - term.hValue * term.k[c] * Bak);
        }
        positionStrainA[c] = -alphaA * (double3::dot(dFieldPositionA, dFieldStrain[a]) + double3::dot(field, secondField));
      }
      result.positionStrain[globalA * numberOfStrainBases + a] += positionStrainA;
    }

    // Position-strain block for each movable real-space source (pure force change): d(grad_source)/ds_a.
    for (const RealSource& source : sources)
    {
      if (!movable[source.globalIndex]) continue;
      const double3 delta = source.d - sigmaA + source.offset;
      for (std::size_t a = 0; a < numberOfStrainBases; ++a)
      {
        double3 positionStrainSource(0.0, 0.0, 0.0);
        for (std::size_t c = 0; c < 3; ++c)
        {
          const double3 axis = unitAxis(c);
          // dE_A/dr_{source,c} = column c of T_s.
          const double3 dFieldPositionSource =
              source.charge * (source.firstDerivativeFactor * axis + source.secondDerivativeFactor * source.d[c] * source.d);
          // d(T_s e_c)/ds_a = -(d^2 e/dd^2)(e_c, B_a delta).
          const double3 secondField = -secondFieldApply(source, axis, strainGenerator[a] * delta);
          positionStrainSource[c] =
              -alphaA * (double3::dot(dFieldPositionSource, dFieldStrain[a]) + double3::dot(field, secondField));
        }
        result.positionStrain[source.globalIndex * numberOfStrainBases + a] += positionStrainSource;
      }
    }
  }

  return result;
}

std::pair<double, double3x3> Interactions::computePolarizationMolecularPressureStrain(
    const System& system, std::span<double3> field, std::span<std::array<double3, 9>> fieldStrain,
    std::span<const double3> centerOfMassOffset, std::span<const double> polarizability)
{
  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  const std::size_t numberOfFrameworkAtoms = frameworkAtoms.size();
  const std::size_t numberOfMoleculeAtoms = moleculeAtoms.size();

  double energy = 0.0;
  double3x3 strainDerivative{};

  if (!forceField.computePolarization || !forceField.useCharge || numberOfMoleculeAtoms == 0)
  {
    return {energy, strainDerivative};
  }

  // Ewald reciprocal field of the rigid fixed framework (matching computeEwaldFourierElectricField and the
  // reciprocal part of computePolarizationDerivatives). The real-space contributions were already gathered
  // by the Coulomb strain loops.
  const bool useReciprocal = forceField.usesEwaldFourier() && system.framework.has_value() &&
                             system.framework->rigid && numberOfFrameworkAtoms > 0;
  if (useReciprocal)
  {
    const double alphaSquared = forceField.EwaldAlpha * forceField.EwaldAlpha;
    const std::vector<PrecomputedReciprocal> precomputedReciprocal =
        precomputeFixedFrameworkReciprocal(forceField, box, frameworkAtoms);

    for (std::size_t a = 0; a < numberOfMoleculeAtoms; ++a)
    {
      if (polarizability[a] == 0.0) continue;
      const double3 positionA = moleculeAtoms[a].position;
      const double3 sigmaA = centerOfMassOffset[a];
      const std::array<double, 3> sigmaComponents = {sigmaA.x, sigmaA.y, sigmaA.z};
      std::array<double3, 9>& tensor = fieldStrain[a];

      for (const PrecomputedReciprocal& entry : precomputedReciprocal)
      {
        const double argumentA = double3::dot(entry.k, positionA);
        const double cosA = std::cos(argumentA);
        const double sinA = std::sin(argumentA);
        const double fValue = sinA * entry.structureReal - cosA * entry.structureImaginary;
        const double hValue = cosA * entry.structureReal + sinA * entry.structureImaginary;

        field[a] += entry.weight * fValue * entry.k;

        // Strain response of the reciprocal field under F = exp(s B): the weight carries
        // beta = -tr B + (1/(4 alpha^2) + 1/k^2) 2 k.Bk, the wave vector deforms as -Bk, and the phase of
        // a COM-scaled field point moves by -(Bk).sigma_A. Linearized in B this gives
        //   dE_A[i]/dF[j][m] += w [ f (2 phi k_j k_m - delta_jm) k_i - h k_i sigma_j k_m - f delta_ij k_m ].
        const double kSquared = double3::dot(entry.k, entry.k);
        const double twoPhi = 2.0 * (0.25 / alphaSquared + 1.0 / kSquared);
        const std::array<double, 3> kComponents = {entry.k.x, entry.k.y, entry.k.z};
        for (std::size_t i = 0; i < 3; ++i)
        {
          for (std::size_t j = 0; j < 3; ++j)
          {
            const double scalar = fValue * twoPhi * kComponents[i] * kComponents[j] -
                                  hValue * kComponents[i] * sigmaComponents[j] - (i == j ? fValue : 0.0);
            double3 contribution = scalar * entry.k;
            contribution[j] -= fValue * kComponents[i];
            tensor[3 * i + j] += entry.weight * contribution;
          }
        }
      }
    }
  }

  // Contract the completed field and its strain response:
  //   U = -1/2 sum_A alpha_A |E_A|^2,   dU/dF[j][m] = -sum_A alpha_A sum_i E_A[i] dE_A[i]/dF[j][m].
  for (std::size_t a = 0; a < numberOfMoleculeAtoms; ++a)
  {
    const double alphaA = polarizability[a];
    if (alphaA == 0.0) continue;

    const double3 electricField = field[a];
    energy -= 0.5 * alphaA * double3::dot(electricField, electricField);

    const std::array<double3, 9>& tensor = fieldStrain[a];
    const double3 row0 =
        -alphaA * (electricField.x * tensor[0] + electricField.y * tensor[3] + electricField.z * tensor[6]);
    const double3 row1 =
        -alphaA * (electricField.x * tensor[1] + electricField.y * tensor[4] + electricField.z * tensor[7]);
    const double3 row2 =
        -alphaA * (electricField.x * tensor[2] + electricField.y * tensor[5] + electricField.z * tensor[8]);

    strainDerivative.ax += row0.x;
    strainDerivative.ay += row0.y;
    strainDerivative.az += row0.z;

    strainDerivative.bx += row1.x;
    strainDerivative.by += row1.y;
    strainDerivative.bz += row1.z;

    strainDerivative.cx += row2.x;
    strainDerivative.cy += row2.y;
    strainDerivative.cz += row2.z;
  }

  return {energy, strainDerivative};
}
