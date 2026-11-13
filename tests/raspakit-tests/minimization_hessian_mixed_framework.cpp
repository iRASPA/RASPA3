#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import component;
import framework;
import connectivity_table;
import intra_molecular_potentials;
import pseudo_atom;
import forcefield;
import system;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import normal_modes;
import phonon_dynamical_matrix;
import phonon_force_constants;
import minimization_generalized_coordinates;
import running_energy;
import interactions_hessian_ewald;
import interactions_hessian_intramolecular;

namespace
{
class TemporaryFile
{
 public:
  TemporaryFile(std::string name, std::string_view contents)
      : path(std::filesystem::temp_directory_path() / std::move(name))
  {
    std::ofstream stream(path);
    stream << contents;
  }
  ~TemporaryFile()
  {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
  }
  std::filesystem::path path;
};

Atom frameworkAtom(double3 fractionalPosition, double charge, std::size_t type)
{
  return Atom(fractionalPosition, charge, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

simd_quatd quatFromRotationVector(const double3& omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

// A six-atom carbon chain host partitioned into Fixed{0} | Rigid{1,2,3} (bent, nonlinear) |
// Flexible{4,5}, with harmonic bonds/bends and a TraPPE torsion whose junction terms cross the group
// boundaries (terms wholly inside the Fixed or Rigid groups are skipped by the reader), plus one
// flexible single-site adsorbate close to the flexible end. The charged variant places a net -0.4 on
// the framework balanced by a +0.4 adsorbate, so the mixed Ewald path (fixed-prefix stored structure
// factor + live rigid/flexible sites + net-charge correction) is exercised.
System makeMixedFrameworkSystem(bool charged)
{
  ForceField forceField(
      {{"C", true, 12.0, 0.0, 0.0, 6, true}, {"M", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}, {60.0, 2.5}},
      ForceField::MixingRule::Lorentz_Berthelot, 5.0, 5.0, 5.0, true, false, charged);
  const SimulationBox box(14.0, 14.0, 14.0);
  const std::size_t typeC = *forceField.findPseudoAtom("C");
  const std::size_t typeM = *forceField.findPseudoAtom("M");

  const std::array<double, 6> charges =
      charged ? std::array<double, 6>{0.3, -0.1, 0.2, -0.1, -0.35, -0.35} : std::array<double, 6>{};

  // ~1.43-1.49 A consecutive spacing (inside covalent-radius bond detection), zigzag in y so the
  // rigid group is nonlinear and all three rotations carry inertia.
  std::vector<Atom> atoms{frameworkAtom({4.0 / 14.0, 7.0 / 14.0, 7.0 / 14.0}, charges[0], typeC),
                          frameworkAtom({5.4 / 14.0, 7.3 / 14.0, 7.0 / 14.0}, charges[1], typeC),
                          frameworkAtom({6.8 / 14.0, 6.8 / 14.0, 7.0 / 14.0}, charges[2], typeC),
                          frameworkAtom({8.2 / 14.0, 7.3 / 14.0, 7.0 / 14.0}, charges[3], typeC),
                          frameworkAtom({9.6 / 14.0, 7.0 / 14.0, 7.0 / 14.0}, charges[4], typeC),
                          frameworkAtom({11.0 / 14.0, 7.2 / 14.0, 7.0 / 14.0}, charges[5], typeC)};
  Framework framework(forceField, "mixed-chain", box, 1, atoms, atoms, int3(1, 1, 1));

  TemporaryFile definition("raspa3-mixed-framework-hessian.json",
                           R"({"Type": "Flexible",
                               "Groups": [
                                 {"Type": "Fixed", "Atoms": [0]},
                                 {"Type": "Rigid", "Atoms": [1, 2, 3]},
                                 {"Type": "Flexible", "Atoms": [4, 5]}
                               ],
                               "Intra14VanDerWaalsScalingValue": 0.0,
                               "Intra14ChargeChargeScalingValue": 0.0,
                               "Bonds": [[["C", "C"], "HARMONIC", [1200.0, 1.45]]],
                               "Bends": [[["C", "C", "C"], "HARMONIC", [250.0, 160.0]]],
                               "Torsions": [[["C", "C", "C", "C"], "TRAPPE", [10.0, 20.0, 30.0, 40.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path.string());
  EXPECT_TRUE(framework.isMixed());
  EXPECT_EQ(framework.numberOfFixedAtoms(), 1u);
  EXPECT_EQ(framework.numberOfRigidGroups(), 1u);
  // Bonds 1-2 and 2-3 (inside the Rigid group) are skipped; junction bonds 0-1, 3-4, 4-5 remain.
  EXPECT_EQ(framework.intraMolecularPotentials.bonds.size(), 3u);

  Component molecule(forceField, "single-site", 16.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, charged ? 0.4 : 0.0, 1.0, 0, static_cast<std::uint16_t>(typeM), 0, false,
                           false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(11.6, 9.8, 7.6);
  return system;
}

struct MixedBaseState
{
  std::vector<double3> frameworkPositions;
  std::vector<double3> moleculePositions;
  std::vector<double3> rigidGroupCenters;  ///< Mass-weighted center per framework group index.
};

MixedBaseState captureBaseState(System& system)
{
  MixedBaseState base{};
  for (const Atom& atom : system.spanOfFrameworkAtoms()) base.frameworkPositions.push_back(atom.position);
  for (const Atom& atom : system.spanOfMoleculeAtoms()) base.moleculePositions.push_back(atom.position);
  base.rigidGroupCenters.resize(system.framework->groups.size());
  for (std::size_t groupIndex = 0; groupIndex < system.framework->groups.size(); ++groupIndex)
  {
    const FrameworkGroup& group = system.framework->groups[groupIndex];
    if (!group.isRigidBody()) continue;
    double3 centerOfMass{};
    for (std::size_t k = 0; k != group.atoms.size(); ++k)
    {
      centerOfMass += group.atomMasses[k] * base.frameworkPositions[group.atoms[k]];
    }
    base.rigidGroupCenters[groupIndex] = centerOfMass / group.mass;
  }
  return base;
}

// Rebuild all positions from the base state under a generalized displacement: rigid framework groups
// translate their center of mass and rotate their base offsets by the exponential map of the
// orientation-tangent components; flexible framework and molecule atoms move Cartesian; fixed
// framework atoms never move.
void applyGeneralizedDisplacement(System& system, const MinimizationDofLayout& layout, const MixedBaseState& base,
                                  std::span<const double> displacement)
{
  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  for (std::size_t groupIndex = 0; groupIndex < system.framework->groups.size(); ++groupIndex)
  {
    const FrameworkGroup& group = system.framework->groups[groupIndex];
    if (!group.isRigidBody()) continue;
    const std::size_t comBase = *layout.frameworkAtomRigidComDof(group.atoms.front());
    const double3 comDisplacement(displacement[comBase + 0], displacement[comBase + 1], displacement[comBase + 2]);
    const double3 omega(displacement[comBase + 3], displacement[comBase + 4], displacement[comBase + 5]);
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(quatFromRotationVector(omega));
    const double3 center = base.rigidGroupCenters[groupIndex];
    for (const std::size_t atom : group.atoms)
    {
      frameworkAtoms[atom].position = center + comDisplacement + rotation * (base.frameworkPositions[atom] - center);
    }
  }
  for (std::size_t atom = 0; atom < frameworkAtoms.size(); ++atom)
  {
    if (const auto dofBase = layout.frameworkAtomDof(atom, MinimizationDofAxis::X))
    {
      frameworkAtoms[atom].position =
          base.frameworkPositions[atom] +
          double3(displacement[*dofBase + 0], displacement[*dofBase + 1], displacement[*dofBase + 2]);
    }
  }
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t atom = 0; atom < moleculeAtoms.size(); ++atom)
  {
    if (const auto dofBase = layout.flexibleAtomDof(0, atom, MinimizationDofAxis::X))
    {
      moleculeAtoms[atom].position =
          base.moleculePositions[atom] +
          double3(displacement[*dofBase + 0], displacement[*dofBase + 1], displacement[*dofBase + 2]);
    }
  }
}

// Total energy along the same code path that assembles the analytic Hessian (the energy-only branch
// of evaluateDerivatives skips the pairwise and Ewald routines).
double energyAtDisplacement(System& system, const MinimizationDofLayout& layout, const MixedBaseState& base,
                            std::span<const double> displacement)
{
  applyGeneralizedDisplacement(system, layout, base, displacement);
  GeneralizedHessian scratch(layout.numDofs(), 0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = false, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = {}, .hessian = scratch};
  evaluateDerivatives(system, layout, capabilities, results);
  return results.energy;
}

MinimizationDofLayout makeMixedLayout(const System& system)
{
  return buildMinimizationDofLayout(system.moleculeData, system.components, system.framework->atoms.size(), 0,
                                    &*system.framework);
}

void expectHessianMatchesFiniteDifference(System& system, double absoluteTolerance, double relativeTolerance,
                                          double delta)
{
  const MinimizationDofLayout layout = makeMixedLayout(system);
  // Fixed atom: 0 DOFs, rigid group: 6, two flexible framework atoms: 6, single-site molecule: 3.
  ASSERT_EQ(layout.numDofs(), 15u);
  ASSERT_FALSE(layout.frameworkAtomDof(0, MinimizationDofAxis::X).has_value());
  ASSERT_TRUE(layout.frameworkAtomRigidComDof(1).has_value());
  ASSERT_EQ(layout.frameworkAtomRigidComDof(1), layout.frameworkAtomRigidComDof(3));

  const MixedBaseState base = captureBaseState(system);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, results);

  std::vector<double> displacement(layout.numDofs(), 0.0);

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    for (std::size_t column = row; column < layout.numDofs(); ++column)
    {
      auto finiteDifference = [&](double step)
      {
        std::ranges::fill(displacement, 0.0);
        displacement[row] += step;
        displacement[column] += step;
        const double ePP = energyAtDisplacement(system, layout, base, displacement);

        std::ranges::fill(displacement, 0.0);
        displacement[row] += step;
        displacement[column] -= step;
        const double ePM = energyAtDisplacement(system, layout, base, displacement);

        std::ranges::fill(displacement, 0.0);
        displacement[row] -= step;
        displacement[column] += step;
        const double eMP = energyAtDisplacement(system, layout, base, displacement);

        std::ranges::fill(displacement, 0.0);
        displacement[row] -= step;
        displacement[column] -= step;
        const double eMM = energyAtDisplacement(system, layout, base, displacement);
        return (ePP - ePM - eMP + eMM) / (4.0 * step * step);
      };
      const double numerical = (4.0 * finiteDifference(delta) - finiteDifference(2.0 * delta)) / 3.0;
      const double analytic = hessian(row, column);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale)
          << "row=" << row << " column=" << column;
      EXPECT_NEAR(hessian(column, row), analytic, 1e-12 * scale) << "asymmetric row=" << row << " column=" << column;
    }
  }

  // Restore base positions.
  std::ranges::fill(displacement, 0.0);
  applyGeneralizedDisplacement(system, layout, base, displacement);
}
}  // namespace

// The rigid groups of a mixed framework move with a generalized displacement whether or not the system keeps
// authoritative group states: without them the state is derived from the atoms, as it is for semi-flexible molecules.
TEST(minimization_hessian_mixed_framework, generalized_displacement_moves_the_rigid_groups_without_group_states)
{
  System withGroupStates = makeMixedFrameworkSystem(false);
  withGroupStates.initializeFrameworkGroupData();
  ASSERT_FALSE(withGroupStates.frameworkGroupData.empty());
  System withoutGroupStates = makeMixedFrameworkSystem(false);
  withoutGroupStates.frameworkGroupData.clear();
  const std::vector<double3> original = captureBaseState(withoutGroupStates).frameworkPositions;

  const MinimizationDofLayout layout = makeMixedLayout(withGroupStates);
  std::vector<double> displacement(layout.numDofs(), 0.0);
  const std::size_t comBase = *layout.frameworkAtomRigidComDof(1);
  displacement[comBase + 0] = 0.10;
  displacement[comBase + 1] = -0.05;
  displacement[comBase + 4] = 0.20;
  ::applyGeneralizedDisplacement(withGroupStates, layout, displacement);
  ::applyGeneralizedDisplacement(withoutGroupStates, layout, displacement);

  const std::span<const Atom> expected = withGroupStates.spanOfFrameworkAtoms();
  const std::span<const Atom> actual = withoutGroupStates.spanOfFrameworkAtoms();
  for (const std::size_t atom : {1uz, 2uz, 3uz})
  {
    EXPECT_GT((expected[atom].position - original[atom]).length(), 0.05) << "atom " << atom;
    EXPECT_NEAR(actual[atom].position.x, expected[atom].position.x, 1e-10) << "atom " << atom;
    EXPECT_NEAR(actual[atom].position.y, expected[atom].position.y, 1e-10) << "atom " << atom;
    EXPECT_NEAR(actual[atom].position.z, expected[atom].position.z, 1e-10) << "atom " << atom;
  }
}

// A mixed framework carries generalized degrees of freedom (none for its fixed atoms, six per rigid group), so the
// Cartesian force constants refuse it like they refuse rigid molecules.
TEST(minimization_hessian_mixed_framework, real_space_force_constants_refuse_a_mixed_framework)
{
  System system = makeMixedFrameworkSystem(false);
  EXPECT_THROW(computeRealSpaceForceConstants(system), std::runtime_error);
}

TEST(minimization_hessian_mixed_framework, generalized_gradient_matches_finite_difference)
{
  System system = makeMixedFrameworkSystem(false);
  const MinimizationDofLayout layout = makeMixedLayout(system);
  ASSERT_EQ(layout.numDofs(), 15u);
  const MixedBaseState base = captureBaseState(system);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, results);

  const double delta = 1.0e-5;
  std::vector<double> displacement(layout.numDofs(), 0.0);
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    std::ranges::fill(displacement, 0.0);
    displacement[dof] = delta;
    const double plus = energyAtDisplacement(system, layout, base, displacement);
    displacement[dof] = -delta;
    const double minus = energyAtDisplacement(system, layout, base, displacement);
    const double numerical = (plus - minus) / (2.0 * delta);
    EXPECT_NEAR(gradient[dof], numerical, 1e-7 + 1e-6 * std::max(1.0, std::abs(numerical))) << "dof=" << dof;
  }
}

TEST(minimization_hessian_mixed_framework, generalized_hessian_matches_finite_difference)
{
  System system = makeMixedFrameworkSystem(false);
  expectHessianMatchesFiniteDifference(system, 1.0e-5, 1.0e-5, 3.0e-4);
}

TEST(minimization_hessian_mixed_framework, ewald_generalized_hessian_matches_finite_difference)
{
  // The total Ewald energy is large, so the four-point finite-difference cancellation leaves more
  // roundoff noise; a larger step and a matching absolute floor keep the comparison meaningful.
  System system = makeMixedFrameworkSystem(true);
  expectHessianMatchesFiniteDifference(system, 1.0e-4, 1.0e-5, 1.2e-3);
}

TEST(minimization_hessian_mixed_framework, phonon_gamma_matches_normal_modes)
{
  System system = makeMixedFrameworkSystem(false);

  // Mixed framework -> generalized dispersion route: the Cartesian force constants of the flattened
  // (fully flexible) copy are projected onto the fixed/rigid/flexible DOFs and the gradient-curvature
  // term is added; at Gamma the spectrum must reproduce the analytic generalized normal modes.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 15u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 15u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }

  // Time-reversal symmetry away from Gamma: omega^2(k) = omega^2(-k).
  const std::array<double3, 2> kPair = {double3(0.3, 0.0, 0.0), double3(-0.3, 0.0, 0.0)};
  const std::vector<PhononModes> offGamma = computePhononDispersion(system, kPair);
  ASSERT_EQ(offGamma.size(), 2u);
  for (std::size_t mode = 0; mode < 15u; ++mode)
  {
    EXPECT_NEAR(offGamma[0].eigenvalues[mode], offGamma[1].eigenvalues[mode], 1e-9 * scale + 1e-12)
        << "mode=" << mode;
  }
}

TEST(minimization_hessian_mixed_framework, phonon_gamma_matches_normal_modes_charged)
{
  System system = makeMixedFrameworkSystem(true);

  // Charged variant: the mixed Ewald Hessian (fixed-prefix stored structure factor + live sites) feeds
  // computeNormalModes, while the phonon route adds the reciprocal-space Ewald matrix of the flattened
  // copy before the projection. Both must agree at Gamma.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 15u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 15u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(minimization_hessian_mixed_framework, writes_normal_mode_movies)
{
  System system = makeMixedFrameworkSystem(false);
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 15u);

  const std::filesystem::path directory = "mixed_framework_normal_modes_test_output";
  std::filesystem::remove_all(directory);
  ASSERT_NO_THROW(writeNormalModeMovies(system, modes, 0, 1, 4, 0.4, directory));

  std::size_t fileCount = 0;
  for (const auto& entry : std::filesystem::directory_iterator(directory))
  {
    if (entry.path().extension() == ".pdb") ++fileCount;
  }
  EXPECT_EQ(fileCount, modes.numberOfModes);
  std::filesystem::remove_all(directory);
}

// The system under the molecular strain: the Fixed and Flexible atoms scale affinely, a Rigid group through its center
// of mass only.
System scaledMixedSystem(const System& system, const MixedBaseState& base, double scale)
{
  System copy = system;
  copy.simulationBox = system.simulationBox.scaled(scale);
  copy.forceField.initializeEwaldParameters(copy.simulationBox);
  std::span<Atom> frameworkAtoms = copy.spanOfFrameworkAtoms();
  for (std::size_t i = 0; i < frameworkAtoms.size(); ++i) frameworkAtoms[i].position = scale * base.frameworkPositions[i];
  for (std::size_t g = 0; g < copy.framework->groups.size(); ++g)
  {
    const FrameworkGroup& group = copy.framework->groups[g];
    if (!group.isRigidBody()) continue;
    for (const std::size_t atom : group.atoms)
    {
      frameworkAtoms[atom].position =
          scale * base.rigidGroupCenters[g] + (base.frameworkPositions[atom] - base.rigidGroupCenters[g]);
    }
  }
  std::span<Atom> moleculeAtoms = copy.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < moleculeAtoms.size(); ++i) moleculeAtoms[i].position = scale * base.moleculePositions[i];
  copy.precomputeTotalRigidEnergy();
  return copy;
}

// The Ewald strain first derivative of a mixed framework is the derivative of the Ewald energy under the molecular
// strain: the Fixed and Flexible atoms scale affinely and a Rigid group scales through its center of mass only, so
// the exclusion pairs inside the group do not strain.
TEST(minimization_hessian_mixed_framework, ewald_strain_gradient_matches_finite_difference)
{
  System system = makeMixedFrameworkSystem(true);
  const MinimizationDofLayout layout = makeMixedLayout(system);
  const MixedBaseState base = captureBaseState(system);

  const auto evaluate = [&](double scale, double3x3& strain)
  {
    System copy = scaledMixedSystem(system, base, scale);

    GeneralizedHessian hessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> moleculeDynamics(copy.spanOfMoleculeAtoms().size());
    std::vector<AtomDynamics> frameworkDynamics(copy.spanOfFrameworkAtoms().size());
    const RunningEnergy energy = Interactions::computeEwaldFourierHessian(
        copy.forceField, copy.simulationBox, copy.framework, copy.fixedFrameworkStoredEik, copy.netChargeFramework,
        copy.moleculeData, copy.components, copy.spanOfFrameworkAtoms(), copy.spanOfMoleculeAtoms(), layout, hessian,
        moleculeDynamics, frameworkDynamics);
    strain = hessian.strainGradient();
    return energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion;
  };

  double3x3 strain{};
  double3x3 ignored{};
  (void)evaluate(1.0, strain);
  const double delta = 1.0e-4;
  const double fd = (evaluate(std::exp(delta), ignored) - evaluate(std::exp(-delta), ignored)) / (2.0 * delta);
  ASSERT_GT(std::abs(fd), 1.0e-3);
  EXPECT_NEAR(strain.trace(), fd, 1.0e-5 * std::abs(fd) + 1.0e-4);
}

// The same for the framework intramolecular terms, whose junction bonds, bends and torsions reach into the Rigid group.
TEST(minimization_hessian_mixed_framework, framework_intramolecular_strain_gradient_matches_finite_difference)
{
  System system = makeMixedFrameworkSystem(true);
  const MinimizationDofLayout layout = makeMixedLayout(system);
  const MixedBaseState base = captureBaseState(system);

  const auto evaluate = [&](double scale, double3x3& strain)
  {
    System copy = scaledMixedSystem(system, base, scale);
    GeneralizedHessian hessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> frameworkDynamics(copy.spanOfFrameworkAtoms().size());
    const RunningEnergy energy = Interactions::computeFrameworkIntraMolecularHessian(
        copy.forceField, *copy.framework, copy.simulationBox, copy.spanOfFrameworkAtoms(), layout, hessian,
        frameworkDynamics);
    strain = hessian.strainGradient();
    return energy.potentialEnergy();
  };

  double3x3 strain{};
  double3x3 ignored{};
  (void)evaluate(1.0, strain);
  const double delta = 1.0e-5;
  const double fd = (evaluate(std::exp(delta), ignored) - evaluate(std::exp(-delta), ignored)) / (2.0 * delta);
  ASSERT_GT(std::abs(fd), 1.0e-3);
  EXPECT_NEAR(strain.trace(), fd, 1.0e-5 * std::abs(fd) + 1.0e-4);
}

// The molecular pressure of a framework with a Rigid group is the derivative of the total energy under the molecular
// strain. The group scales through its center of mass, so every force on its atoms enters the center-of-mass virial,
// including the real-space forces of the adsorbate.
TEST(minimization_hessian_mixed_framework, molecular_pressure_matches_finite_difference_with_rigid_group)
{
  for (const bool charged : {false, true})
  {
    System system = makeMixedFrameworkSystem(charged);
    const MixedBaseState base = captureBaseState(system);

    const auto energyAtScale = [&](double scale)
    {
      System copy = scaledMixedSystem(system, base, scale);
      return copy.computeTotalEnergies().potentialEnergy();
    };
    const double delta = 1.0e-5;
    const double scalePlus = std::cbrt(1.0 + delta);
    const double scaleMinus = std::cbrt(1.0 - delta);
    const double volumePlus = system.simulationBox.scaled(scalePlus).volume;
    const double volumeMinus = system.simulationBox.scaled(scaleMinus).volume;
    const double fd = -(energyAtScale(scalePlus) - energyAtScale(scaleMinus)) / (volumePlus - volumeMinus);
    ASSERT_GT(std::abs(fd), 1.0e-6);

    const double3x3 virial = system.computeMolecularPressure().second;
    EXPECT_NEAR(virial.trace() / (3.0 * system.simulationBox.volume), fd, 1.0e-4 * std::abs(fd) + 1.0e-8)
        << "charged " << charged;
  }
}
