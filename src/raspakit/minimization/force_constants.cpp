module;

module phonon_force_constants;

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import molecule;
import component;
import forcefield;
import framework;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import bend_potential_gradient_hessian_strain;
import torsion_potential_gradient_hessian_strain;
import generalized_hessian;
import minimization_hessian_scatter;
import minimization_dof_layout;
import system;
import interactions_pair_kernel;
import interactions_hessian_intramolecular;
import interactions_polarization_derivatives;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;

namespace
{
// Cartesian pair Hessian block B_{alpha,beta} = f2 dr_alpha dr_beta + f1 delta_{alpha,beta}.
std::array<double, 9> pairBlock(double f1, double f2, const double3& dr)
{
  const std::array<double, 3> d = {dr.x, dr.y, dr.z};
  std::array<double, 9> block{};
  for (std::size_t alpha = 0; alpha < 3; ++alpha)
  {
    for (std::size_t beta = 0; beta < 3; ++beta)
    {
      block[alpha * 3 + beta] = f2 * d[alpha] * d[beta] + (alpha == beta ? f1 : 0.0);
    }
  }
  return block;
}

// Integer lattice vector connecting the raw separation to the minimum-image separation.
// R_cart = (posA - posB) - dr is an exact lattice translation, so the rounding is unambiguous.
int3 latticeVectorFromImage(const SimulationBox& box, const double3& rawSeparation, const double3& dr)
{
  const double3 latticeCartesian = rawSeparation - dr;
  const double3 fractional = box.inverseCell * latticeCartesian;
  return int3(static_cast<std::int32_t>(std::llround(fractional.x)),
              static_cast<std::int32_t>(std::llround(fractional.y)),
              static_cast<std::int32_t>(std::llround(fractional.z)));
}

// Distribute a local dense Cartesian Hessian (built with atom bases {0, 3, 6, ...}) into the
// image-resolved force constants. By lattice translational invariance a term whose atom i sits in
// image cell shifts[i] couples atoms i and j through Phi_{ids[i], ids[j]}(shifts[j] - shifts[i]).
// Summing over all terms/images reproduces the folded (Gamma-point) framework Hessian.
void distributeLocalHessian(ForceConstants& forceConstants, std::span<const std::size_t> ids,
                            std::span<const int3> shifts, const GeneralizedHessian& local)
{
  const std::size_t n = ids.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      std::array<double, 9> block{};
      for (std::size_t alpha = 0; alpha < 3; ++alpha)
      {
        for (std::size_t beta = 0; beta < 3; ++beta)
        {
          block[alpha * 3 + beta] = local(3 * i + alpha, 3 * j + beta);
        }
      }
      forceConstants.addBlock(ids[i], ids[j], shifts[j] - shifts[i], block);
    }
  }
}

// Framework-framework intramolecular force constants (bonds, bends, torsions, improper torsions and
// the scaled/nonbonded intra van der Waals and Coulomb pairs) for a flexible framework, resolved by
// periodic image. Mirrors Interactions::computeFrameworkIntraMolecularHessian exactly, but records
// each term's local Cartesian block against the lattice vector between its atoms instead of folding.
void addFrameworkIntramolecular(const System& system, ForceConstants& forceConstants)
{
  const Framework& framework = *system.framework;
  if (framework.rigid) return;

  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> atoms = system.spanOfFrameworkAtoms();

  const auto shiftedPosition = [&](std::size_t atom, const int3& shift)
  {
    return atoms[atom].position + box.cell * double3(static_cast<double>(shift.x), static_cast<double>(shift.y),
                                                     static_cast<double>(shift.z));
  };

  const auto& potentials = framework.intraMolecularPotentials;
  const auto& images = framework.intraMolecularImageShifts;

  // 1-2/1-3 exclusions and 1-4 scaled pairs, reconstructed from the connectivity table exactly as in
  // the analytic framework-intra Hessian.
  std::set<std::array<std::size_t, 2>> pairs12And13;
  std::set<std::array<std::size_t, 2>> pairs14;
  if (!framework.connectivityTable.table.empty())
  {
    for (const std::array<std::size_t, 2>& bond : framework.connectivityTable.findAllBonds())
    {
      pairs12And13.insert({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
    }
    for (const std::array<std::size_t, 3>& bend : framework.connectivityTable.findAllBends())
    {
      pairs12And13.insert({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
    }
    for (const std::array<std::size_t, 4>& torsion : framework.connectivityTable.findAllTorsions())
    {
      const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
      if (!pairs12And13.contains(pair)) pairs14.insert(pair);
    }
  }

  const auto distributeRadial = [&](std::size_t A, std::size_t B, const std::array<int3, 2>& shifts, double f1,
                                    double f2, const double3& dr)
  {
    GeneralizedHessian local(6, 0);
    Minimization::scatterAtomicPositionPositionByDof(local, 0, 3, f1, f2, dr);
    const std::array<std::size_t, 2> ids{A, B};
    distributeLocalHessian(forceConstants, ids, shifts, local);
  };

  for (std::size_t index = 0; index < potentials.bonds.size(); ++index)
  {
    const BondPotential& bond = potentials.bonds[index];
    const std::size_t A = bond.identifiers[0];
    const std::size_t B = bond.identifiers[1];
    const double3 positionA = shiftedPosition(A, images.bonds[index][0]);
    const double3 positionB = shiftedPosition(B, images.bonds[index][1]);
    auto [energy, gradient, strain, f1, f2] = bond.potentialEnergyGradientHessianStrain(positionA, positionB);
    distributeRadial(A, B, images.bonds[index], f1, f2, positionA - positionB);
  }

  for (std::size_t index = 0; index < potentials.bends.size(); ++index)
  {
    const BendPotential& bend = potentials.bends[index];
    const auto& ids = bend.identifiers;
    const auto& shifts = images.bends[index];
    const std::array<double3, 3> shifted = {shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                            shiftedPosition(ids[2], shifts[2])};
    auto [energy, gradient, strain, geometry] = Potentials::Internal::bendPotentialEnergyGradientHessianStrain(
        bend.type, bend.parameters, shifted[0], shifted[1], shifted[2]);
    GeneralizedHessian local(9, 0);
    Minimization::scatterBendHessianByDof(local, {0, 3, 6}, geometry);
    distributeLocalHessian(forceConstants, ids, shifts, local);
  }

  const auto distributeTorsions = [&](const std::vector<TorsionPotential>& torsions,
                                      const std::vector<std::array<int3, 4>>& termImages)
  {
    for (std::size_t index = 0; index < torsions.size(); ++index)
    {
      const TorsionPotential& torsion = torsions[index];
      const auto& ids = torsion.identifiers;
      const auto& shifts = termImages[index];
      const std::array<double3, 4> shifted = {shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                              shiftedPosition(ids[2], shifts[2]), shiftedPosition(ids[3], shifts[3])};
      auto [energy, gradient, strain, geometry] = Potentials::Internal::torsionPotentialEnergyGradientHessianStrain(
          torsion.type, torsion.parameters, shifted[0], shifted[1], shifted[2], shifted[3]);
      GeneralizedHessian local(12, 0);
      Minimization::scatterTorsionHessianByDof(local, {0, 3, 6, 9}, geometry);
      distributeLocalHessian(forceConstants, ids, shifts, local);
    }
  };
  distributeTorsions(potentials.torsions, images.torsions);
  distributeTorsions(potentials.improperTorsions, images.improperTorsions);

  for (std::size_t index = 0; index < potentials.vanDerWaals.size(); ++index)
  {
    const VanDerWaalsPotential& potential = potentials.vanDerWaals[index];
    const std::size_t A = potential.identifiers[0];
    const std::size_t B = potential.identifiers[1];
    const double3 dr =
        shiftedPosition(A, images.vanDerWaals[index][0]) - shiftedPosition(B, images.vanDerWaals[index][1]);
    const double rr = double3::dot(dr, dr);
    const bool isScaled14 = pairs14.contains({std::min(A, B), std::max(A, B)}) && potential.scaling != 1.0;
    if (isScaled14 || rr < forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW)
    {
      const Potentials::PairDerivatives<2> factors = Potentials::potentialVDW<2>(
          forceField, 1.0, 1.0, rr, static_cast<std::size_t>(atoms[A].type), static_cast<std::size_t>(atoms[B].type));
      distributeRadial(A, B, images.vanDerWaals[index], potential.scaling * factors.firstDerivativeFactor,
                       potential.scaling * factors.secondDerivativeFactor, dr);
    }
  }

  for (std::size_t index = 0; index < potentials.coulombs.size(); ++index)
  {
    const CoulombPotential& potential = potentials.coulombs[index];
    const std::size_t A = potential.identifiers[0];
    const std::size_t B = potential.identifiers[1];
    const double3 dr = shiftedPosition(A, images.coulombs[index][0]) - shiftedPosition(B, images.coulombs[index][1]);
    const double rr = double3::dot(dr, dr);
    if (pairs14.contains({std::min(A, B), std::max(A, B)}) && potential.scaling != 1.0)
    {
      const double r = std::sqrt(rr);
      const double k = potential.scaling * Units::CoulombicConversionFactor * potential.chargeA * potential.chargeB;
      distributeRadial(A, B, images.coulombs[index], -k / (r * rr), 3.0 * k / (r * rr * rr), dr);
    }
    else if (rr < forceField.cutOffCoulomb * forceField.cutOffCoulomb)
    {
      const double r = std::sqrt(rr);
      const Potentials::PairDerivatives<2> factors =
          Potentials::potentialCoulomb<2>(forceField, potential.scaling, 1.0, r, potential.chargeA, potential.chargeB);
      distributeRadial(A, B, images.coulombs[index], factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
    }
  }
}

// Molecule intramolecular terms (bonds, bends, urey-bradley, torsions, cross terms and the intra
// van der Waals / Coulomb pairs). These are evaluated from the direct atom positions without periodic
// image shifts, so every contribution lands in the home image Phi(0). Assembled by reusing the analytic
// intramolecular Hessian routines with a full-system DOF layout, whose position-position block indexes
// atoms exactly as the force-constant sites (flexible framework atoms first, then molecule atoms).
void addMoleculeIntramolecular(const System& system, ForceConstants& forceConstants,
                               std::size_t numberOfFlexibleFrameworkAtoms)
{
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  if (moleculeAtoms.empty()) return;

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, numberOfFlexibleFrameworkAtoms);
  const std::size_t dimension = layout.numDofs();
  if (dimension != 3 * forceConstants.numberOfAtoms())
  {
    throw std::runtime_error("computeRealSpaceForceConstants: molecule intramolecular DOF layout mismatch");
  }

  GeneralizedHessian local(dimension, 0);
  std::vector<AtomDynamics> dynamicsStorage(moleculeAtoms.size());
  const std::span<const Molecule> moleculeData = system.moleculeData;
  const std::span<const Component> components = system.components;
  const std::span<AtomDynamics> dynamics = dynamicsStorage;

  Interactions::computeIntraMolecularBondHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularBendHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularUreyBradleyHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                        dynamics);
  Interactions::computeIntraMolecularTorsionHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularVanDerWaalsHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                        dynamics);
  Interactions::computeIntraMolecularCoulombHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularBondBondHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularBondBendHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularBendBendHessian(moleculeData, moleculeAtoms, components, layout, local, dynamics);
  Interactions::computeIntraMolecularBondTorsionHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                        dynamics);
  Interactions::computeIntraMolecularBendTorsionHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                        dynamics);
  Interactions::computeIntraMolecularInversionBendHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                          dynamics);
  Interactions::computeIntraMolecularOutOfPlaneBendHessian(moleculeData, moleculeAtoms, components, layout, local,
                                                           dynamics);

  forceConstants.addToImage(int3(0, 0, 0), local.positionPosition());
}

// Ewald intramolecular exclusion correction. The reciprocal Fourier sum counts every intramolecular
// charge pair; to obtain the intended electrostatics that spurious contribution is removed by the
// pairwise term U(r) = -C q_i q_j erf(alpha r) / r for all pairs within a molecule. This mirrors the
// molecule branch of Interactions::computeEwaldFourierHessian exactly, but records each pair's
// Cartesian block against the lattice vector between the two atoms instead of folding to Gamma.
void addMoleculeEwaldExclusion(const System& system, ForceConstants& forceConstants,
                               std::size_t numberOfFrameworkSites)
{
  const ForceField& forceField = system.forceField;
  if (!forceField.useCharge || !forceField.usesEwaldFourier()) return;

  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  const double alpha = forceField.EwaldAlpha;
  const double alphaSquared = alpha * alpha;

  for (const Molecule& molecule : system.moleculeData)
  {
    if (molecule.numberOfAtoms < 2) continue;
    const std::span<const Atom> atoms = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
    for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
    {
      for (std::size_t j = i + 1; j < atoms.size(); ++j)
      {
        const double3 rawSeparation = atoms[i].position - atoms[j].position;
        const double3 dr = box.applyPeriodicBoundaryConditions(rawSeparation);
        const double rr = double3::dot(dr, dr);
        const double r = std::sqrt(rr);

        const double chargeProduct = Units::CoulombicConversionFactor * atoms[i].scalingCoulomb *
                                     atoms[j].scalingCoulomb * atoms[i].charge * atoms[j].charge;
        if (chargeProduct == 0.0) continue;

        const double erfTerm = std::erf(alpha * r);
        const double gaussTerm = 2.0 * alpha * std::numbers::inv_sqrtpi * std::exp(-alphaSquared * rr);
        // U(r) = -chargeProduct erf(alpha r)/r; RASPA convention f1 = U'/r, f2 = (U'' - U'/r)/r^2.
        const double f1 = -chargeProduct * (gaussTerm / rr - erfTerm / (r * rr));
        const double f2 = chargeProduct * (2.0 * alphaSquared * gaussTerm / rr + 3.0 * gaussTerm / (rr * rr) -
                                           3.0 * erfTerm / (r * rr * rr));

        const std::array<double, 9> block = pairBlock(f1, f2, dr);
        std::array<double, 9> negativeBlock{};
        for (std::size_t index = 0; index < 9; ++index) negativeBlock[index] = -block[index];
        const int3 latticeVector = latticeVectorFromImage(box, rawSeparation, dr);

        const std::size_t siteI = numberOfFrameworkSites + molecule.atomIndex + i;
        const std::size_t siteJ = numberOfFrameworkSites + molecule.atomIndex + j;
        forceConstants.addBlock(siteI, siteI, int3(0, 0, 0), block);
        forceConstants.addBlock(siteJ, siteJ, int3(0, 0, 0), block);
        forceConstants.addBlock(siteI, siteJ, latticeVector, negativeBlock);
        forceConstants.addBlock(siteJ, siteI, -latticeVector, negativeBlock);
      }
    }
  }
}

// Ewald intramolecular exclusion correction for a flexible framework. Unlike molecules, a framework
// keeps its nonbonded intraframework Coulomb, so only the bonded/excluded pairs (bonds, 1-3 bends and
// 1-4 torsions that are absent from the Coulomb list or scaled away) receive the -C q_i q_j erf/r
// correction. Mirrors the framework branch of Interactions::computeEwaldFourierHessian; framework atoms
// are the leading force-constant sites, so their site index equals the atom index.
void addFrameworkEwaldExclusion(const System& system, ForceConstants& forceConstants)
{
  const ForceField& forceField = system.forceField;
  if (!forceField.useCharge || !forceField.usesEwaldFourier()) return;
  if (!system.framework.has_value() || system.framework->rigid) return;

  const Framework& framework = *system.framework;
  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  const double alpha = forceField.EwaldAlpha;
  const double alphaSquared = alpha * alpha;

  std::set<std::array<std::size_t, 2>> excludedPairs;
  std::map<std::array<std::size_t, 2>, double> coulombScaling;
  for (const CoulombPotential& potential : framework.intraMolecularPotentials.coulombs)
  {
    coulombScaling[{std::min(potential.identifiers[0], potential.identifiers[1]),
                    std::max(potential.identifiers[0], potential.identifiers[1])}] = potential.scaling;
  }
  const auto excludeIfAbsentOrScaled = [&](const std::array<std::size_t, 2>& pair)
  {
    const auto scaling = coulombScaling.find(pair);
    if (scaling == coulombScaling.end() || scaling->second != 1.0) excludedPairs.insert(pair);
  };
  if (!framework.connectivityTable.table.empty())
  {
    for (const std::array<std::size_t, 2>& bond : framework.connectivityTable.findAllBonds())
    {
      excludeIfAbsentOrScaled({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
    }
    for (const std::array<std::size_t, 3>& bend : framework.connectivityTable.findAllBends())
    {
      excludeIfAbsentOrScaled({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
    }
    for (const std::array<std::size_t, 4>& torsion : framework.connectivityTable.findAllTorsions())
    {
      excludeIfAbsentOrScaled({std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])});
    }
  }

  for (const std::array<std::size_t, 2>& pair : excludedPairs)
  {
    const std::size_t i = pair[0];
    const std::size_t j = pair[1];
    const double3 rawSeparation = frameworkAtoms[i].position - frameworkAtoms[j].position;
    const double3 dr = box.applyPeriodicBoundaryConditions(rawSeparation);
    const double rr = double3::dot(dr, dr);
    const double r = std::sqrt(rr);

    const double chargeProduct = Units::CoulombicConversionFactor * frameworkAtoms[i].scalingCoulomb *
                                 frameworkAtoms[j].scalingCoulomb * frameworkAtoms[i].charge * frameworkAtoms[j].charge;
    if (chargeProduct == 0.0) continue;

    const double erfTerm = std::erf(alpha * r);
    const double gaussTerm = 2.0 * alpha * std::numbers::inv_sqrtpi * std::exp(-alphaSquared * rr);
    const double f1 = -chargeProduct * (gaussTerm / rr - erfTerm / (r * rr));
    const double f2 = chargeProduct * (2.0 * alphaSquared * gaussTerm / rr + 3.0 * gaussTerm / (rr * rr) -
                                       3.0 * erfTerm / (r * rr * rr));

    const std::array<double, 9> block = pairBlock(f1, f2, dr);
    std::array<double, 9> negativeBlock{};
    for (std::size_t index = 0; index < 9; ++index) negativeBlock[index] = -block[index];
    const int3 latticeVector = latticeVectorFromImage(box, rawSeparation, dr);

    forceConstants.addBlock(i, i, int3(0, 0, 0), block);
    forceConstants.addBlock(j, j, int3(0, 0, 0), block);
    forceConstants.addBlock(i, j, latticeVector, negativeBlock);
    forceConstants.addBlock(j, i, -latticeVector, negativeBlock);
  }
}
}  // namespace

ForceConstants::ForceConstants(std::size_t numberOfAtoms)
    : _numberOfAtoms(numberOfAtoms), _dimension(3 * numberOfAtoms)
{
}

void ForceConstants::addBlock(std::size_t i, std::size_t j, int3 latticeVector, const std::array<double, 9>& block)
{
  std::vector<double>& matrix = _blocks.try_emplace(latticeVector, _dimension * _dimension, 0.0).first->second;
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      matrix[(3 * i + row) * _dimension + (3 * j + column)] += block[row * 3 + column];
    }
  }
}

void ForceConstants::addToImage(int3 latticeVector, std::span<const double> matrix)
{
  if (matrix.size() != _dimension * _dimension)
  {
    throw std::runtime_error("ForceConstants::addToImage: matrix extent does not match the force-constant dimension");
  }
  std::vector<double>& block = _blocks.try_emplace(latticeVector, _dimension * _dimension, 0.0).first->second;
  for (std::size_t index = 0; index < block.size(); ++index) block[index] += matrix[index];
}

std::vector<double> ForceConstants::foldToGamma() const
{
  std::vector<double> folded(_dimension * _dimension, 0.0);
  for (const auto& [latticeVector, matrix] : _blocks)
  {
    for (std::size_t index = 0; index < folded.size(); ++index)
    {
      folded[index] += matrix[index];
    }
  }
  return folded;
}

ForceConstants computeRealSpaceForceConstants(const System& system)
{
  for (const Component& component : system.components)
  {
    if (component.rigid)
    {
      throw std::runtime_error(
          "computeRealSpaceForceConstants: rigid molecules are not supported in this increment (flexible atoms only)");
    }
  }
  const bool flexibleFramework = system.framework.has_value() && !system.framework->rigid;
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();

  // Global atom ordering matches the minimization DOF layout: flexible framework atoms first, then
  // molecule atoms. A rigid framework contributes no degrees of freedom, so its atoms are not sites.
  const std::size_t numberOfFrameworkSites = flexibleFramework ? frameworkAtoms.size() : 0;
  const std::size_t numberOfSites = numberOfFrameworkSites + moleculeAtoms.size();
  ForceConstants forceConstants(numberOfSites);

  // Nonbonded (van der Waals + real-space Coulomb) pair: atom i sits in the home cell, atom j in the
  // periodic image R. When j carries no degrees of freedom (rigid framework) only the home self-block
  // of atom i is recorded, matching scatterFrameworkMoleculeHessian.
  const auto addNonbondedPair = [&](std::size_t i, std::size_t j, const Potentials::PairDerivatives<2>& factors,
                                    const double3& rawSeparation, const double3& dr, bool jHasDof)
  {
    const std::array<double, 9> block = pairBlock(factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
    forceConstants.addBlock(i, i, int3(0, 0, 0), block);
    if (!jHasDof) return;

    std::array<double, 9> negativeBlock{};
    for (std::size_t index = 0; index < 9; ++index) negativeBlock[index] = -block[index];
    const int3 latticeVector = latticeVectorFromImage(system.simulationBox, rawSeparation, dr);
    forceConstants.addBlock(j, j, int3(0, 0, 0), block);
    forceConstants.addBlock(i, j, latticeVector, negativeBlock);
    forceConstants.addBlock(j, i, -latticeVector, negativeBlock);
  };

  // Molecule-molecule intermolecular pairs.
  Interactions::forEachMoleculeMoleculePair<2>(
      system.forceField, system.simulationBox, moleculeAtoms,
      [&](std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
          const Potentials::PairDerivatives<2>& factors, const double3& dr)
      {
        addNonbondedPair(numberOfFrameworkSites + indexA, numberOfFrameworkSites + indexB, factors,
                         atomA.position - atomB.position, dr, true);
      },
      [&](std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
          const Potentials::PairDerivatives<2>& factors, const double3& dr)
      {
        addNonbondedPair(numberOfFrameworkSites + indexA, numberOfFrameworkSites + indexB, factors,
                         atomA.position - atomB.position, dr, true);
      });

  // Framework-molecule pairs. dr = moleculeAtom - frameworkAtom, so the molecule atom is the home site.
  if (!frameworkAtoms.empty() && !moleculeAtoms.empty())
  {
    for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeAtoms.size(); ++moleculeIndex)
    {
      const Atom& moleculeAtom = moleculeAtoms[moleculeIndex];
      const std::size_t moleculeSite = numberOfFrameworkSites + moleculeIndex;
      Interactions::forEachFrameworkMoleculePair<2>(
          system.forceField, system.simulationBox, moleculeAtom, frameworkAtoms,
          [&](std::size_t frameworkIndex, const Atom& frameworkAtom, const Potentials::PairDerivatives<2>& factors,
              const double3& dr)
          {
            addNonbondedPair(moleculeSite, frameworkIndex, factors, moleculeAtom.position - frameworkAtom.position, dr,
                             flexibleFramework);
          },
          [&](std::size_t frameworkIndex, const Atom& frameworkAtom, const Potentials::PairDerivatives<2>& factors,
              const double3& dr)
          {
            addNonbondedPair(moleculeSite, frameworkIndex, factors, moleculeAtom.position - frameworkAtom.position, dr,
                             flexibleFramework);
          });
    }
  }

  // Framework-framework intramolecular terms (flexible framework only).
  if (flexibleFramework)
  {
    addFrameworkIntramolecular(system, forceConstants);
    addFrameworkEwaldExclusion(system, forceConstants);
  }

  // Molecule intramolecular terms (image-free, home cell).
  addMoleculeIntramolecular(system, forceConstants, numberOfFrameworkSites);

  // Ewald intramolecular exclusion correction (removes the reciprocal-counted intramolecular charge
  // pairs). Its short-ranged Hessian folds directly into the image-resolved force constants.
  addMoleculeEwaldExclusion(system, forceConstants, numberOfFrameworkSites);

  // Polarization (non-self-consistent induced-dipole) Hessian, resolved by periodic image so the Bloch
  // phase is applied per lattice vector. The real-space field-gradient blocks connecting a home-cell
  // field point to a source in image R are recorded against R (many-body Term A three-center blocks land
  // at R_C - R_B); the reciprocal-space (fixed framework) contribution is an on-site R = 0 term that is
  // k-independent. Summing over R reproduces the folded Gamma-point Hessian.
  if (system.forceField.computePolarization)
  {
    std::vector<std::uint8_t> movable(frameworkAtoms.size() + moleculeAtoms.size(), 0);
    if (flexibleFramework)
    {
      for (std::size_t atom = 0; atom < frameworkAtoms.size(); ++atom) movable[atom] = 1;
    }
    for (std::size_t atom = 0; atom < moleculeAtoms.size(); ++atom) movable[frameworkAtoms.size() + atom] = 1;

    const Interactions::PolarizationDerivatives derivatives =
        Interactions::computePolarizationDerivatives(system, movable);

    const auto siteOfGlobal = [&](std::size_t global) -> std::size_t
    {
      if (global < frameworkAtoms.size()) return global;  // flexible framework atom -> leading site
      return numberOfFrameworkSites + (global - frameworkAtoms.size());
    };

    for (const auto& [latticeVector, blocks] : derivatives.imageHessianBlocks)
    {
      for (const auto& [pair, block] : blocks)
      {
        const std::size_t i = pair[0];
        const std::size_t j = pair[1];
        if (!movable[i] || !movable[j]) continue;
        forceConstants.addBlock(siteOfGlobal(i), siteOfGlobal(j), latticeVector, block);
      }
    }
  }

  return forceConstants;
}
