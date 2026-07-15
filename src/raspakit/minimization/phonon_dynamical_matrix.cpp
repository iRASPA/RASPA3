module;

module phonon_dynamical_matrix;

import std;

import units;
import int3;
import double3;
import double3x3;
import atom;
import molecule;
import component;
import framework;
import forcefield;
import simulationbox;
import system;
import phonon_force_constants;
import phonon_kpath;
import hermitian_eigensolver;
import generalized_hessian;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_rigid_kinematics;
import minimization_mass_metric;
import minimization_evaluate_derivatives;

namespace
{
double checkedInverseSqrtMass(double mass, std::string_view what)
{
  if (!(mass > 0.0) || !std::isfinite(mass))
  {
    throw std::runtime_error(std::format("Phonon analysis requires positive masses, but {} has mass {}", what, mass));
  }
  return 1.0 / std::sqrt(mass);
}

/** A charged site participating in the reciprocal sum, in force-constant site order. */
struct ChargedSite
{
  double3 position{};
  double charge{};  // scalingCoulomb * charge
};

/**
 * Collect the charged force-constant sites in the same order as computeRealSpaceForceConstants
 * (flexible framework atoms first, then molecule atoms). Atoms of a rigid framework carry no degrees
 * of freedom and are therefore not sites; their (fixed) charges are returned separately by
 * gatherStaticFrameworkCharges and enter only the molecule self term.
 */
std::vector<ChargedSite> gatherReciprocalSites(const System& system)
{
  std::vector<ChargedSite> sites;

  const bool flexibleFramework = system.framework.has_value() && !system.framework->rigid;
  if (flexibleFramework)
  {
    for (const Atom& atom : system.spanOfFrameworkAtoms())
    {
      sites.push_back(ChargedSite{.position = atom.position, .charge = atom.scalingCoulomb * atom.charge});
    }
  }

  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      throw std::runtime_error("computeEwaldReciprocalDynamicalMatrix: rigid molecules are not supported");
    }

    // Charged multi-atom molecules are handled: the intramolecular exclusion correction is folded into
    // the real-space force constants by computeRealSpaceForceConstants.
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      const Atom& atom = moleculeAtoms[molecule.atomIndex + localAtom];
      sites.push_back(ChargedSite{.position = atom.position, .charge = atom.scalingCoulomb * atom.charge});
    }
  }
  return sites;
}

/**
 * Fixed charges of a rigid framework. They do not move, so they contribute no matrix rows/columns, but
 * they add a static structure factor S_fixed(G) that shifts every molecule self term (each molecule
 * atom still interacts with all framework images). Empty for a flexible or uncharged framework.
 */
std::vector<ChargedSite> gatherStaticFrameworkCharges(const System& system)
{
  std::vector<ChargedSite> staticCharges;
  if (!system.framework.has_value() || !system.framework->rigid) return staticCharges;
  for (const Atom& atom : system.spanOfFrameworkAtoms())
  {
    const double charge = atom.scalingCoulomb * atom.charge;
    if (charge != 0.0) staticCharges.push_back(ChargedSite{.position = atom.position, .charge = charge});
  }
  return staticCharges;
}

bool systemHasRigidMolecules(const System& system)
{
  return std::ranges::any_of(system.moleculeData,
                             [&](const Molecule& molecule) { return system.components[molecule.componentId].rigid; });
}

/**
 * Unweighted Cartesian reciprocal-space Ewald dynamical matrix D_recip(k), dimension 3 * (number of
 * charged sites), in force-constant site order (flexible framework atoms first, then molecule atoms). No
 * mass weighting is applied: the flexible path multiplies by M^{-1/2} afterwards, while the rigid
 * (generalized) path adds this matrix to the Cartesian force constants and projects it onto the
 * center-of-mass/orientation coordinates before mass weighting with the generalized metric. Returns an
 * empty vector when electrostatics or the Fourier sum are disabled.
 */
std::vector<std::complex<double>> ewaldReciprocalCartesianMatrix(const System& system, double3 kFractional)
{
  const ForceField& forceField = system.forceField;
  if (!forceField.useCharge || !forceField.usesEwaldFourier()) return {};

  const std::vector<ChargedSite> sites = gatherReciprocalSites(system);
  const std::size_t numberOfAtoms = sites.size();
  const std::size_t dimension = 3 * numberOfAtoms;
  std::vector<std::complex<double>> matrix(dimension * dimension, std::complex<double>{});
  const std::vector<ChargedSite> staticCharges = gatherStaticFrameworkCharges(system);

  const SimulationBox& simulationBox = system.simulationBox;
  const double alpha = forceField.EwaldAlpha;
  const double alphaSquared = alpha * alpha;
  const double reciprocalCutOffSquared = forceField.reciprocalCutOffSquared;
  const std::size_t reciprocalIntegerCutOffSquared = forceField.reciprocalIntegerCutOffSquared;

  // Reciprocal-lattice basis (rows of the inverse cell), matching computeEwaldFourierHessian.
  const double3x3 inverseCell = simulationBox.inverseCell;
  const double3 aStar(inverseCell.ax, inverseCell.bx, inverseCell.cx);
  const double3 bStar(inverseCell.ay, inverseCell.by, inverseCell.cy);
  const double3 cStar(inverseCell.az, inverseCell.bz, inverseCell.cz);

  constexpr double twoPi = 2.0 * std::numbers::pi;
  // Full-space (both +/- G) Hessian coefficient per wave vector: 2 * (2 pi C / V) exp(-|q|^2/4a^2)/|q|^2.
  // The leading 2 is the factor arising from the two complex-conjugate terms of d^2|S|^2.
  const double coefficientPrefactor = 2.0 * Units::CoulombicConversionFactor * (twoPi / simulationBox.volume);

  const double3 kCartesian = twoPi * (kFractional.x * aStar + kFractional.y * bStar + kFractional.z * cStar);

  const int kxMax = static_cast<int>(forceField.numberOfWaveVectors.x);
  const int kyMax = static_cast<int>(forceField.numberOfWaveVectors.y);
  const int kzMax = static_cast<int>(forceField.numberOfWaveVectors.z);

  std::vector<std::complex<double>> phase(numberOfAtoms);  // e_a = q_a exp(i q . r_a)
  std::vector<double3> selfDiagonal(numberOfAtoms, double3(0.0, 0.0, 0.0));  // diag(G_a G_a) contributions
  std::vector<double> selfDiagonalXY(numberOfAtoms, 0.0);                    // off-diagonal xy/xz/yz
  std::vector<double> selfDiagonalXZ(numberOfAtoms, 0.0);
  std::vector<double> selfDiagonalYZ(numberOfAtoms, 0.0);

  // Loop the full reciprocal lattice; +G and -G partners keep D(0) real and give the acoustic sum rule.
  for (int kx = -kxMax; kx <= kxMax; ++kx)
  {
    for (int ky = -kyMax; ky <= kyMax; ++ky)
    {
      for (int kz = -kzMax; kz <= kzMax; ++kz)
      {
        const std::size_t integerSquared = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if (integerSquared == 0uz || integerSquared > reciprocalIntegerCutOffSquared)
        {
          continue;
        }

        const double3 gVector = twoPi * (static_cast<double>(kx) * aStar + static_cast<double>(ky) * bStar +
                                         static_cast<double>(kz) * cStar);

        // --- Self term (uses G at k = 0); its G-set matches the Fourier sum, keeping the sum rule exact.
        const double gSquared = double3::dot(gVector, gVector);
        if (gSquared < reciprocalCutOffSquared)
        {
          const double weight = coefficientPrefactor * std::exp(-0.25 * gSquared / alphaSquared) / gSquared;
          std::complex<double> structureFactor(0.0, 0.0);
          for (std::size_t site = 0; site < numberOfAtoms; ++site)
          {
            const double argument = double3::dot(gVector, sites[site].position);
            structureFactor += sites[site].charge * std::polar(1.0, argument);
          }
          // Rigid-framework charges are fixed but still shift every molecule self term.
          for (const ChargedSite& fixedCharge : staticCharges)
          {
            const double argument = double3::dot(gVector, fixedCharge.position);
            structureFactor += fixedCharge.charge * std::polar(1.0, argument);
          }
          for (std::size_t site = 0; site < numberOfAtoms; ++site)
          {
            const double argument = double3::dot(gVector, sites[site].position);
            const std::complex<double> conjugatePhase = std::polar(1.0, -argument);
            // sum_l q_l cos(G.(r_site - r_l)) = Re[exp(-iG.r_site) S(G)].
            const double neighbourSum = (conjugatePhase * structureFactor).real();
            const double scalar = weight * sites[site].charge * neighbourSum;
            selfDiagonal[site].x -= scalar * gVector.x * gVector.x;
            selfDiagonal[site].y -= scalar * gVector.y * gVector.y;
            selfDiagonal[site].z -= scalar * gVector.z * gVector.z;
            selfDiagonalXY[site] -= scalar * gVector.x * gVector.y;
            selfDiagonalXZ[site] -= scalar * gVector.x * gVector.z;
            selfDiagonalYZ[site] -= scalar * gVector.y * gVector.z;
          }
        }

        // --- Pair term (uses q = k + G).
        const double3 qVector = kCartesian + gVector;
        const double qSquared = double3::dot(qVector, qVector);
        if (qSquared < 1.0e-12 || qSquared >= reciprocalCutOffSquared)
        {
          continue;
        }
        const double weight = coefficientPrefactor * std::exp(-0.25 * qSquared / alphaSquared) / qSquared;
        for (std::size_t site = 0; site < numberOfAtoms; ++site)
        {
          const double argument = double3::dot(qVector, sites[site].position);
          phase[site] = sites[site].charge * std::polar(1.0, argument);
        }
        const std::array<double, 3> q = {qVector.x, qVector.y, qVector.z};
        for (std::size_t rowSite = 0; rowSite < numberOfAtoms; ++rowSite)
        {
          for (std::size_t columnSite = 0; columnSite < numberOfAtoms; ++columnSite)
          {
            const std::complex<double> pairPhase = phase[rowSite] * std::conj(phase[columnSite]);
            for (std::size_t alpha2 = 0; alpha2 < 3; ++alpha2)
            {
              for (std::size_t beta = 0; beta < 3; ++beta)
              {
                const std::size_t row = 3 * rowSite + alpha2;
                const std::size_t column = 3 * columnSite + beta;
                matrix[row * dimension + column] += weight * q[alpha2] * q[beta] * pairPhase;
              }
            }
          }
        }
      }
    }
  }

  // Add the (real, k-independent) self blocks on the diagonal sites.
  for (std::size_t site = 0; site < numberOfAtoms; ++site)
  {
    const std::size_t base = 3 * site;
    matrix[(base + 0) * dimension + (base + 0)] += selfDiagonal[site].x;
    matrix[(base + 1) * dimension + (base + 1)] += selfDiagonal[site].y;
    matrix[(base + 2) * dimension + (base + 2)] += selfDiagonal[site].z;
    matrix[(base + 0) * dimension + (base + 1)] += selfDiagonalXY[site];
    matrix[(base + 1) * dimension + (base + 0)] += selfDiagonalXY[site];
    matrix[(base + 0) * dimension + (base + 2)] += selfDiagonalXZ[site];
    matrix[(base + 2) * dimension + (base + 0)] += selfDiagonalXZ[site];
    matrix[(base + 1) * dimension + (base + 2)] += selfDiagonalYZ[site];
    matrix[(base + 2) * dimension + (base + 1)] += selfDiagonalYZ[site];
  }
  return matrix;
}

/**
 * Maps a Cartesian force-constant site displacement to the generalized DOFs that drive it: a single
 * contribution `direction * dq` per generalized DOF. Flexible atoms and rigid center-of-mass DOFs use
 * unit directions; rigid orientation DOFs use the rotation Jacobian dVec = e_axis x (r - r_com).
 */
struct SiteDofContribution
{
  std::size_t dof{};
  double3 direction{};
};

/**
 * Generalized dynamical matrix along a k-path for systems containing rigid molecules.
 *
 * The energy Hessian in generalized (center-of-mass + orientation) coordinates is H_gen = J^T H_cart J +
 * Sum_i grad_i . ddVec_i, where J is the rigid-body kinematic Jacobian (columns dVec) and the second term
 * is the geometric (gradient-curvature) contribution that RASPA's orientation-orientation scatter folds
 * into each molecule's own block. Applying this relation image by image gives
 *
 *   D_gen(k) = J^T [ Sum_R Phi_cart(R) e^{2 pi i k.R} + D_recip_cart(k) ] J + Xi,
 *
 * with Xi (the gradient-curvature term) k-independent because it is on-site. The reciprocal-space Ewald
 * term D_recip_cart(k) is a position-dependent potential like any other, so it enters the projection the
 * same way (and its force is already included in the Cartesian gradient that drives Xi). After mass
 * weighting with the generalized metric (mass for translations, inertia tensor for rotations) the
 * Hermitian matrix is diagonalized. At k = 0 this reproduces computeNormalModes exactly.
 */
std::vector<PhononModes> computeGeneralizedDispersion(const System& system, std::span<const double3> kPath)
{
  // Fully-flexible copy: rigid molecule atoms become independent Cartesian sites (no intramolecular
  // stiffness, since rigid internal geometry is a constraint), which is what the Cartesian force
  // constants and Cartesian gradients require.
  System cartesianSystem = system;
  cartesianSystem.cellMinimizationType = CellMinimizationType::Fixed;
  for (Component& component : cartesianSystem.components)
  {
    component.rigid = false;
    component.intraMolecularPotentials = {};
  }

  const ForceConstants forceConstants = computeRealSpaceForceConstants(cartesianSystem);
  const std::size_t numberOfSites = forceConstants.numberOfAtoms();
  const std::size_t cartesianDimension = forceConstants.dimension();

  const std::size_t numberOfFlexibleFrameworkAtoms =
      (system.framework && !system.framework->rigid) ? system.spanOfFrameworkAtoms().size() : 0;

  // Cartesian atomic gradients: with an all-flexible layout the generalized gradient equals the Cartesian
  // per-atom gradient dU/dr (RASPA convention grad_A = f1 dr), ordered as the force-constant sites.
  const MinimizationDofLayout cartesianLayout = buildMinimizationDofLayout(
      cartesianSystem.moleculeData, cartesianSystem.components, numberOfFlexibleFrameworkAtoms, 0);
  if (cartesianLayout.numDofs() != cartesianDimension)
  {
    throw std::runtime_error("computeGeneralizedDispersion: Cartesian DOF layout does not match the force constants");
  }
  std::vector<double> cartesianGradient(cartesianDimension, 0.0);
  {
    // The intermolecular / framework-molecule gradient is only accumulated alongside the Hessian, so the
    // Hessian capability must be enabled for the Cartesian atomic gradient to be complete. The Hessian
    // itself is discarded here; the image-resolved force constants provide the k-dependent matrix.
    GeneralizedHessian scratchHessian(cartesianDimension, 0);
    DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
    DerivativeResults results{.gradient = cartesianGradient, .hessian = scratchHessian};
    evaluateDerivatives(cartesianSystem, cartesianLayout, capabilities, results);
  }

  // Generalized (rigid-aware) layout and rigid-body kinematics for the projection.
  const MinimizationDofLayout generalizedLayout =
      buildMinimizationDofLayout(system.moleculeData, system.components, numberOfFlexibleFrameworkAtoms, 0);
  const std::size_t numberOfGeneralizedDofs = generalizedLayout.numDofs();

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, system.spanOfMoleculeAtoms());

  const std::array<double3, 3> unit = {double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0)};

  // Columns of the kinematic Jacobian J: per Cartesian site, the generalized DOFs that move it.
  std::vector<std::vector<SiteDofContribution>> siteContributions(numberOfSites);
  for (std::size_t atom = 0; atom < numberOfFlexibleFrameworkAtoms; ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      siteContributions[atom].push_back(
          {*generalizedLayout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis)), unit[axis]});
    }
  }
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const bool rigid = system.components[molecule.componentId].rigid;
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      const std::size_t site = numberOfFlexibleFrameworkAtoms + molecule.atomIndex + localAtom;
      if (rigid)
      {
        const std::size_t comBase = *generalizedLayout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
        for (std::size_t axis = 0; axis < 3; ++axis) siteContributions[site].push_back({comBase + axis, unit[axis]});
        const std::size_t orientationBase = *generalizedLayout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
        const Minimization::RigidAtomDerivatives& derivatives = rigidCache.atom(moleculeIndex, localAtom);
        siteContributions[site].push_back({orientationBase + 0, derivatives.dVecX});
        siteContributions[site].push_back({orientationBase + 1, derivatives.dVecY});
        siteContributions[site].push_back({orientationBase + 2, derivatives.dVecZ});
      }
      else
      {
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          siteContributions[site].push_back(
              {*generalizedLayout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis)),
               unit[axis]});
        }
      }
    }
  }

  // Gradient-curvature term Xi (real, k-independent): only rigid orientation blocks receive it.
  std::vector<double> gradientCurvature(numberOfGeneralizedDofs * numberOfGeneralizedDofs, 0.0);
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    if (!system.components[molecule.componentId].rigid) continue;
    const std::size_t orientationBase = *generalizedLayout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      const std::size_t site = numberOfFlexibleFrameworkAtoms + molecule.atomIndex + localAtom;
      const double3 gradient(cartesianGradient[3 * site + 0], cartesianGradient[3 * site + 1],
                             cartesianGradient[3 * site + 2]);
      const Minimization::RigidAtomDerivatives& derivatives = rigidCache.atom(moleculeIndex, localAtom);
      const std::array<std::array<double3, 3>, 3> ddVec = {
          {{derivatives.ddVecAX, derivatives.ddVecAY, derivatives.ddVecAZ},
           {derivatives.ddVecAY, derivatives.ddVecBY, derivatives.ddVecBZ},
           {derivatives.ddVecAZ, derivatives.ddVecBZ, derivatives.ddVecCZ}}};
      for (std::size_t alpha = 0; alpha < 3; ++alpha)
      {
        for (std::size_t beta = 0; beta < 3; ++beta)
        {
          gradientCurvature[(orientationBase + alpha) * numberOfGeneralizedDofs + (orientationBase + beta)] +=
              double3::dot(gradient, ddVec[alpha][beta]);
        }
      }
    }
  }

  const MassMetric metric = buildMassMetric(system, generalizedLayout);
  constexpr double twoPi = 2.0 * std::numbers::pi;

  // Charged rigid molecules: the long-ranged reciprocal-space Ewald term is a function of the atomic
  // positions like any other pair potential, so the rigid-body identity D_gen = J^T D_cart J + Xi holds
  // for it as well. Its unweighted Cartesian matrix is added to the short-ranged force constants before
  // projection; the intramolecular exclusion is already folded into the force constants and the total
  // Cartesian gradient (which drives Xi) already includes the reciprocal force.
  const bool useEwald = system.forceField.useCharge && system.forceField.usesEwaldFourier();

  std::vector<PhononModes> dispersion;
  dispersion.reserve(kPath.size());
  for (const double3 kFractional : kPath)
  {
    // Un-weighted Cartesian dynamical matrix D_cart(k) = Sum_R Phi(R) exp(2 pi i k.R).
    std::vector<std::complex<double>> cartesian(cartesianDimension * cartesianDimension, std::complex<double>{});
    for (const auto& [latticeVector, block] : forceConstants.blocks())
    {
      const double argument = twoPi * (kFractional.x * static_cast<double>(latticeVector.x) +
                                       kFractional.y * static_cast<double>(latticeVector.y) +
                                       kFractional.z * static_cast<double>(latticeVector.z));
      const std::complex<double> phase = std::polar(1.0, argument);
      for (std::size_t index = 0; index < cartesian.size(); ++index) cartesian[index] += block[index] * phase;
    }

    if (useEwald)
    {
      const std::vector<std::complex<double>> reciprocal = ewaldReciprocalCartesianMatrix(cartesianSystem, kFractional);
      if (reciprocal.size() != cartesian.size())
      {
        throw std::runtime_error(
            "computeGeneralizedDispersion: reciprocal-space matrix extent does not match the Cartesian force constants");
      }
      for (std::size_t index = 0; index < cartesian.size(); ++index) cartesian[index] += reciprocal[index];
    }

    // Project onto generalized coordinates: D_gen = J^T D_cart J.
    std::vector<std::complex<double>> generalized(numberOfGeneralizedDofs * numberOfGeneralizedDofs,
                                                  std::complex<double>{});
    for (std::size_t rowSite = 0; rowSite < numberOfSites; ++rowSite)
    {
      if (siteContributions[rowSite].empty()) continue;
      for (std::size_t columnSite = 0; columnSite < numberOfSites; ++columnSite)
      {
        if (siteContributions[columnSite].empty()) continue;
        std::array<std::complex<double>, 9> pairBlock{};
        for (std::size_t alpha = 0; alpha < 3; ++alpha)
        {
          for (std::size_t beta = 0; beta < 3; ++beta)
          {
            pairBlock[alpha * 3 + beta] = cartesian[(3 * rowSite + alpha) * cartesianDimension + (3 * columnSite + beta)];
          }
        }
        for (const SiteDofContribution& rowContribution : siteContributions[rowSite])
        {
          const double3 u = rowContribution.direction;
          std::array<std::complex<double>, 3> uTimesBlock{};
          for (std::size_t beta = 0; beta < 3; ++beta)
          {
            uTimesBlock[beta] =
                u.x * pairBlock[0 * 3 + beta] + u.y * pairBlock[1 * 3 + beta] + u.z * pairBlock[2 * 3 + beta];
          }
          for (const SiteDofContribution& columnContribution : siteContributions[columnSite])
          {
            const double3 v = columnContribution.direction;
            const std::complex<double> value = uTimesBlock[0] * v.x + uTimesBlock[1] * v.y + uTimesBlock[2] * v.z;
            generalized[rowContribution.dof * numberOfGeneralizedDofs + columnContribution.dof] += value;
          }
        }
      }
    }

    for (std::size_t index = 0; index < generalized.size(); ++index) generalized[index] += gradientCurvature[index];

    // Mass weighting D <- W D W with the real symmetric metric, applied to real and imaginary parts.
    std::vector<double> realPart(generalized.size());
    std::vector<double> imaginaryPart(generalized.size());
    for (std::size_t index = 0; index < generalized.size(); ++index)
    {
      realPart[index] = generalized[index].real();
      imaginaryPart[index] = generalized[index].imag();
    }
    const std::vector<double> weightedReal = applyMassMetric(realPart, numberOfGeneralizedDofs, metric);
    const std::vector<double> weightedImaginary = applyMassMetric(imaginaryPart, numberOfGeneralizedDofs, metric);
    std::vector<std::complex<double>> weighted(generalized.size());
    for (std::size_t index = 0; index < generalized.size(); ++index)
    {
      weighted[index] = std::complex<double>(weightedReal[index], weightedImaginary[index]);
    }

    const HermitianEigenSystem eigensystem = diagonalizeHermitian(weighted, numberOfGeneralizedDofs);
    dispersion.push_back(PhononModes{.kFractional = kFractional, .eigenvalues = eigensystem.eigenvalues});
  }
  return dispersion;
}
}  // namespace

std::vector<double> phononInverseSqrtMasses(const System& system)
{
  std::vector<double> weights;

  const bool flexibleFramework = system.framework.has_value() && !system.framework->rigid;
  if (flexibleFramework)
  {
    const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
    for (std::size_t atom = 0; atom < frameworkAtoms.size(); ++atom)
    {
      const std::size_t type = static_cast<std::size_t>(frameworkAtoms[atom].type);
      weights.push_back(checkedInverseSqrtMass(
          system.forceField.pseudoAtoms[type].mass,
          std::format("framework atom {} ({})", atom, system.forceField.pseudoAtoms[type].name)));
    }
  }

  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      throw std::runtime_error("phononInverseSqrtMasses: rigid molecules are not supported");
    }
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      const std::size_t type = static_cast<std::size_t>(moleculeAtoms[molecule.atomIndex + localAtom].type);
      weights.push_back(checkedInverseSqrtMass(
          system.forceField.pseudoAtoms[type].mass,
          std::format("atom {} of molecule {} ({})", localAtom, moleculeIndex, component.name)));
    }
  }
  return weights;
}

std::vector<std::complex<double>> computeDynamicalMatrix(const ForceConstants& forceConstants,
                                                         std::span<const double> inverseSqrtMass, double3 kFractional)
{
  const std::size_t dimension = forceConstants.dimension();
  if (inverseSqrtMass.size() != forceConstants.numberOfAtoms())
  {
    throw std::invalid_argument("computeDynamicalMatrix: mass vector length does not match the number of atoms");
  }

  std::vector<std::complex<double>> matrix(dimension * dimension, std::complex<double>{});
  constexpr double twoPi = 2.0 * std::numbers::pi;
  for (const auto& [latticeVector, block] : forceConstants.blocks())
  {
    const double phaseArgument = twoPi * (kFractional.x * static_cast<double>(latticeVector.x) +
                                          kFractional.y * static_cast<double>(latticeVector.y) +
                                          kFractional.z * static_cast<double>(latticeVector.z));
    const std::complex<double> phase = std::polar(1.0, phaseArgument);
    for (std::size_t index = 0; index < matrix.size(); ++index)
    {
      matrix[index] += block[index] * phase;
    }
  }

  // Symmetric mass weighting D <- M^{-1/2} D M^{-1/2}; each atom carries a scalar 1/sqrt(m).
  for (std::size_t row = 0; row < dimension; ++row)
  {
    const double rowWeight = inverseSqrtMass[row / 3];
    for (std::size_t column = 0; column < dimension; ++column)
    {
      matrix[row * dimension + column] *= rowWeight * inverseSqrtMass[column / 3];
    }
  }
  return matrix;
}

std::vector<std::complex<double>> computeEwaldReciprocalDynamicalMatrix(const System& system,
                                                                       std::span<const double> inverseSqrtMass,
                                                                       double3 kFractional)
{
  const std::size_t numberOfAtoms = inverseSqrtMass.size();
  const std::size_t dimension = 3 * numberOfAtoms;

  std::vector<std::complex<double>> matrix = ewaldReciprocalCartesianMatrix(system, kFractional);
  if (matrix.empty())
  {
    // Electrostatics or the Fourier sum are disabled: the reciprocal contribution vanishes.
    return std::vector<std::complex<double>>(dimension * dimension, std::complex<double>{});
  }
  if (matrix.size() != dimension * dimension)
  {
    throw std::invalid_argument(
        "computeEwaldReciprocalDynamicalMatrix: mass vector length does not match the number of charged sites");
  }

  // Symmetric mass weighting D <- M^{-1/2} D M^{-1/2}.
  for (std::size_t row = 0; row < dimension; ++row)
  {
    const double rowWeight = inverseSqrtMass[row / 3];
    for (std::size_t column = 0; column < dimension; ++column)
    {
      matrix[row * dimension + column] *= rowWeight * inverseSqrtMass[column / 3];
    }
  }
  return matrix;
}

PhononModes computePhononModes(const ForceConstants& forceConstants, std::span<const double> inverseSqrtMass,
                               double3 kFractional)
{
  const std::vector<std::complex<double>> dynamicalMatrix =
      computeDynamicalMatrix(forceConstants, inverseSqrtMass, kFractional);
  const HermitianEigenSystem eigensystem = diagonalizeHermitian(dynamicalMatrix, forceConstants.dimension());
  return PhononModes{.kFractional = kFractional, .eigenvalues = eigensystem.eigenvalues};
}

PhononModes computePhononModes(const System& system, const ForceConstants& forceConstants,
                               std::span<const double> inverseSqrtMass, double3 kFractional)
{
  std::vector<std::complex<double>> dynamicalMatrix =
      computeDynamicalMatrix(forceConstants, inverseSqrtMass, kFractional);
  const std::vector<std::complex<double>> reciprocal =
      computeEwaldReciprocalDynamicalMatrix(system, inverseSqrtMass, kFractional);
  for (std::size_t index = 0; index < dynamicalMatrix.size(); ++index)
  {
    dynamicalMatrix[index] += reciprocal[index];
  }
  const HermitianEigenSystem eigensystem = diagonalizeHermitian(dynamicalMatrix, forceConstants.dimension());
  return PhononModes{.kFractional = kFractional, .eigenvalues = eigensystem.eigenvalues};
}

std::vector<PhononModes> computePhononDispersion(const System& system, std::span<const double3> kPath)
{
  if (systemHasRigidMolecules(system))
  {
    return computeGeneralizedDispersion(system, kPath);
  }

  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);

  std::vector<PhononModes> dispersion;
  dispersion.reserve(kPath.size());
  for (const double3 k : kPath)
  {
    dispersion.push_back(computePhononModes(system, forceConstants, inverseSqrtMass, k));
  }
  return dispersion;
}

std::vector<double> phononFrequenciesWavenumber(std::span<const double> eigenvalues)
{
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
  // Internal angular frequency (per internal time unit) to cm^-1.
  const double conversion =
      reduced ? 1.0 : 1.0 / (2.0 * std::numbers::pi * 100.0 * Units::SpeedOfLight * Units::TimeConversionFactor);

  std::vector<double> frequencies;
  frequencies.reserve(eigenvalues.size());
  for (const double eigenvalue : eigenvalues)
  {
    const double magnitude = std::sqrt(std::abs(eigenvalue)) * conversion;
    frequencies.push_back(eigenvalue < 0.0 ? -magnitude : magnitude);
  }
  return frequencies;
}

PhononDispersionResult computePhononDispersionAlongPath(const System& system, std::span<const PhononPathNode> nodes,
                                                        std::size_t pointsPerSegment)
{
  PhononDispersionResult result;
  result.path = buildPhononKPath(nodes, pointsPerSegment, system.simulationBox);

  std::vector<double3> kFractionalPoints;
  kFractionalPoints.reserve(result.path.size());
  for (const PhononKPoint& point : result.path) kFractionalPoints.push_back(point.kFractional);

  result.modes = computePhononDispersion(system, kFractionalPoints);
  return result;
}

std::string writePhononDispersion(const PhononDispersionResult& result)
{
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
  const std::size_t numberOfBands = result.modes.empty() ? 0 : result.modes.front().eigenvalues.size();

  std::string output = std::format(
      "\nPhonon dispersion\n"
      "Number of k-points: {}\n"
      "Number of bands:    {}\n"
      "Frequencies are wavenumbers [{}]; negative values denote imaginary (unstable) modes.\n\n",
      result.path.size(), numberOfBands, reduced ? "reduced" : "cm^-1");

  for (std::size_t index = 0; index < result.path.size(); ++index)
  {
    const PhononKPoint& point = result.path[index];
    output += std::format("k-point {:5d}  {:<4}  ({: .5f}, {: .5f}, {: .5f})  path={: .6f}\n", index,
                          point.label.empty() ? "-" : point.label, point.kFractional.x, point.kFractional.y,
                          point.kFractional.z, point.pathCoordinate);
    const std::vector<double> frequencies = phononFrequenciesWavenumber(result.modes[index].eigenvalues);
    for (std::size_t band = 0; band < frequencies.size(); ++band)
    {
      output += std::format("    band {:5d}: {: 18.8f}\n", band, frequencies[band]);
    }
  }
  return output;
}
