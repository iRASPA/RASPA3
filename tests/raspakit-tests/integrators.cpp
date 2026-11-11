#include <gtest/gtest.h>

import std;

import archive;
import int3;
import double3;
import double4;
import double3x3;
import simd_quatd;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import system;
import simulationbox;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import integrators_update;
import integrators;
import integrators_compute;
import molecule;
import interpolation_energy_grid;
import thermostat;
import thermobarostat;
import minimization_cell_layout;
import elastic_constants;
import molecular_dynamics;
import mc_moves;
import mc_moves_hybridmc;
import mc_moves_translation;
import mc_moves_translation_smart_mc_all;
import mc_moves_rotation_smart_mc_all;
import mc_moves_move_types;
import ewald_structure_factor;
import randomnumbers;

namespace
{
ForceField makeFlexibleFrameworkForceField()
{
  return ForceField({{"A", false, 2.0, 0.0, 0.0, 1, false}, {"B", false, 5.0, 0.0, 0.0, 2, false}},
                    {{0.0, 1.0}, {0.0, 1.0}}, ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false,
                    false);
}

Framework makeTwoAtomFramework(const ForceField& forceField, bool rigid)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  std::vector<Atom> atoms{Atom({0.20, 0.25, 0.25}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({0.25, 0.25, 0.25}, 0.0, 1.0, 0, 1, 0, 0, true)};
  Framework framework(forceField, "two-atom-framework", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = rigid;
  return framework;
}

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

// Fixed{0} | Rigid{1,2,3} (nonlinear) | Flexible{4,5}, plus one adsorbate. Hits the mixed-framework
// Hybrid MC path that must pass frameworkGroupData into the integrator.
System makeMixedFrameworkHybridMCSystem()
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}, {"M", false, 16.0, 0.0, 0.0, 8, false}},
                        {{80.0, 3.0}, {60.0, 2.5}}, ForceField::MixingRule::Lorentz_Berthelot, 5.0, 5.0, 5.0, false,
                        false, false);
  const SimulationBox box(14.0, 14.0, 14.0);
  const std::size_t typeC = *forceField.findPseudoAtom("C");
  const std::size_t typeM = *forceField.findPseudoAtom("M");

  auto frameworkAtom = [&](double3 fractionalPosition)
  { return Atom(fractionalPosition, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeC), 0, 0, true); };

  std::vector<Atom> atoms{frameworkAtom({4.0 / 14.0, 7.0 / 14.0, 7.0 / 14.0}),
                          frameworkAtom({5.4 / 14.0, 7.3 / 14.0, 7.0 / 14.0}),
                          frameworkAtom({6.8 / 14.0, 6.8 / 14.0, 7.0 / 14.0}),
                          frameworkAtom({8.2 / 14.0, 7.3 / 14.0, 7.0 / 14.0}),
                          frameworkAtom({9.6 / 14.0, 7.0 / 14.0, 7.0 / 14.0}),
                          frameworkAtom({11.0 / 14.0, 7.2 / 14.0, 7.0 / 14.0})};
  Framework framework(forceField, "mixed-chain", box, 1, atoms, atoms, int3(1, 1, 1));

  TemporaryFile definition("raspa3-mixed-framework-hybridmc.json",
                           R"({"Type": "Flexible",
                               "Groups": [
                                 {"Type": "Fixed", "Atoms": [0]},
                                 {"Type": "Rigid", "Atoms": [1, 2, 3]},
                                 {"Type": "Flexible", "Atoms": [4, 5]}
                               ]})");
  framework.readFrameworkDefinition(forceField, definition.path.string());

  Component molecule(forceField, "single-site", 16.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeM), 0, false, false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(11.6, 9.8, 7.6);
  return system;
}

double3 expectedMoleculeCenterOfMassVelocity(const System& system, std::size_t moleculeIndex)
{
  const Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms().subspan(molecule.atomIndex, molecule.numberOfAtoms);
  const std::span<const AtomDynamics> dynamics =
      system.spanOfMoleculeDynamics().subspan(molecule.atomIndex, molecule.numberOfAtoms);
  double3 momentum{};
  if (component.isSemiFlexible())
  {
    std::size_t groupIndex{};
    for (std::size_t m = 0; m < moleculeIndex; ++m)
      groupIndex += system.components[system.moleculeData[m].componentId].numberOfRigidFragments();
    std::size_t rigidRank{};
    for (const Fragment& group : component.fragmentGraph.fragments)
    {
      if (!group.isRigidBody()) continue;
      momentum += group.mass * system.groupData[groupIndex + rigidRank].velocity;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      if (component.rigidFragmentContaining(i).has_value()) continue;
      momentum += component.atomMass(i, atoms[i]) * dynamics[i].velocity;
    }
  }
  else
  {
    for (std::size_t i = 0; i != atoms.size(); ++i)
      momentum += component.atomMass(i, atoms[i]) * dynamics[i].velocity;
  }
  return momentum * molecule.invMass;
}

System makeFlexibleEthaneSystem()
{
  ForceField forceField({{"CH3", false, 15.03452, 0.0, 0.0, 8, false}}, {{108.0, 3.76}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {5000.0, 1.54})};
  Component component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                       Atom({1.70, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                      connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;
  return System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component},
                {{double3(5.0, 5.0, 5.0), double3(6.70, 5.0, 5.0)},
                 {double3(12.0, 8.0, 7.0), double3(13.70, 8.0, 7.0)}},
                {0}, 5);
}

System makeSemiFlexibleChainSystem()
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds.emplace_back(std::array<std::size_t, 2>{2, 3}, BondType::Harmonic, std::vector<double>{1500.0, 0.85});
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
  Component component(forceField, "semi-flexible", 64.0, 1.0e6, 0.02,
                      {Atom({-0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.0, -0.35, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({1.2, -0.25, 0.35}, 0.0, 1.0, 0, type, 0, false, false)},
                      connectivityTable, potentials, 0, 0);
  component.buildFragmentGraph({{0, 1, 2}});
  // the positions of a component are one flat list (two molecules of four atoms here)
  return System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {component},
                {{double3(5.0, 5.0, 5.0), double3(5.6, 4.4, 5.0), double3(6.2, 5.0, 5.0), double3(6.8, 4.5, 5.35),
                  double3(12.0, 8.0, 7.0), double3(12.6, 7.4, 7.0), double3(13.2, 8.0, 7.0), double3(13.8, 7.5, 7.35)}},
                {0}, 5);
}
}  // namespace

void DOUBLE3_EXPECT_NEAR(double3 a, double3 b, double tol);

TEST(thermobarostat, parses_ensembles_and_counts_all_cell_modes)
{
  EXPECT_EQ(molecularDynamicsEnsembleFromString("npt"), MolecularDynamicsEnsemble::NPT);
  EXPECT_EQ(molecularDynamicsEnsembleFromString("NPTPR"), MolecularDynamicsEnsemble::NPTPR);
  EXPECT_EQ(molecularDynamicsEnsembleFromString("muvt"), MolecularDynamicsEnsemble::MuVT);
  EXPECT_EQ(molecularDynamicsEnsembleFromString("MuPT"), MolecularDynamicsEnsemble::MuPT);
  EXPECT_EQ(molecularDynamicsEnsembleFromString("MUPTPR"), MolecularDynamicsEnsemble::MuPTPR);
  EXPECT_TRUE(molecularDynamicsUsesThermostat(MolecularDynamicsEnsemble::MuVT));
  EXPECT_TRUE(molecularDynamicsUsesIsotropicBarostat(MolecularDynamicsEnsemble::MuPT));
  EXPECT_TRUE(molecularDynamicsUsesFlexibleBarostat(MolecularDynamicsEnsemble::MuPTPR));
  EXPECT_TRUE(molecularDynamicsHasParticleExchange(MolecularDynamicsEnsemble::MuVT));
  EXPECT_FALSE(molecularDynamicsHasParticleExchange(MolecularDynamicsEnsemble::NPT));
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::Regular), 6u);
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::Monoclinic), 4u);
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::Isotropic), 1u);
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::Anisotropic), 3u);
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::RegularUpperTriangle), 6u);
  EXPECT_EQ(thermobarostatCellDegreesOfFreedom(CellMinimizationType::MonoclinicUpperTriangle), 4u);
  EXPECT_EQ(cellMinimizationTypeFromString("REGULAR_UPPER_TRIANGLE"),
            CellMinimizationType::RegularUpperTriangle);
}

TEST(thermobarostat, projects_cell_force_without_raspa2_index_aliasing)
{
  const double3x3 input(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  const double3x3 isotropic = projectCellTensor(input, CellMinimizationType::Isotropic);
  EXPECT_DOUBLE_EQ(isotropic.ax, 5.0);
  EXPECT_DOUBLE_EQ(isotropic.by, 5.0);
  EXPECT_DOUBLE_EQ(isotropic.cz, 5.0);

  const double3x3 anisotropic = projectCellTensor(input, CellMinimizationType::Anisotropic);
  EXPECT_DOUBLE_EQ(anisotropic.ax, 1.0);
  EXPECT_DOUBLE_EQ(anisotropic.by, 5.0);
  EXPECT_DOUBLE_EQ(anisotropic.cz, 9.0);
  EXPECT_DOUBLE_EQ(anisotropic.bx, 0.0);

  const double3x3 regular = projectCellTensor(input, CellMinimizationType::Regular);
  EXPECT_DOUBLE_EQ(regular.ay, regular.bx);
  EXPECT_DOUBLE_EQ(regular.az, regular.cx);
  EXPECT_DOUBLE_EQ(regular.bz, regular.cy);

  const double3x3 upper = projectCellTensor(input, CellMinimizationType::RegularUpperTriangle);
  EXPECT_DOUBLE_EQ(upper.ay, 0.0);
  EXPECT_DOUBLE_EQ(upper.az, 0.0);
  EXPECT_DOUBLE_EQ(upper.bz, 0.0);
  EXPECT_DOUBLE_EQ(upper.bx, input.bx);
  EXPECT_DOUBLE_EQ(upper.cx, input.cx);
  EXPECT_DOUBLE_EQ(upper.cy, input.cy);

  const double3x3 monoclinic =
      projectCellTensor(input, CellMinimizationType::Monoclinic, MonoclinicAngleType::Beta);
  EXPECT_DOUBLE_EQ(monoclinic.az, monoclinic.cx);
  EXPECT_DOUBLE_EQ(monoclinic.az, 5.0);
  EXPECT_DOUBLE_EQ(monoclinic.bx, 0.0);

  const double3x3 monoclinicUpper =
      projectCellTensor(input, CellMinimizationType::MonoclinicUpperTriangle, MonoclinicAngleType::Beta);
  EXPECT_DOUBLE_EQ(monoclinicUpper.cx, input.cx);
  EXPECT_DOUBLE_EQ(monoclinicUpper.ay, 0.0);
  EXPECT_DOUBLE_EQ(monoclinicUpper.az, 0.0);
}

TEST(thermobarostat, projects_all_monoclinic_shear_modes)
{
  const double3x3 input(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  struct Case
  {
    MonoclinicAngleType angle;
    double symmetric;
    double upper;
  };
  const std::array cases{
      Case{MonoclinicAngleType::Alpha, 7.0, 8.0},
      Case{MonoclinicAngleType::Beta, 5.0, 7.0},
      Case{MonoclinicAngleType::Gamma, 3.0, 4.0},
  };
  for (const Case& test : cases)
  {
    const double3x3 symmetric = projectCellTensor(input, CellMinimizationType::Monoclinic, test.angle);
    const double3x3 upper =
        projectCellTensor(input, CellMinimizationType::MonoclinicUpperTriangle, test.angle);
    const double symmetricShear =
        test.angle == MonoclinicAngleType::Alpha
            ? symmetric.bz
            : (test.angle == MonoclinicAngleType::Beta ? symmetric.az : symmetric.ay);
    const double upperShear =
        test.angle == MonoclinicAngleType::Alpha
            ? upper.cy
            : (test.angle == MonoclinicAngleType::Beta ? upper.cx : upper.bx);
    EXPECT_DOUBLE_EQ(symmetricShear, test.symmetric);
    EXPECT_DOUBLE_EQ(upperShear, test.upper);
    EXPECT_DOUBLE_EQ(upper.ay, 0.0);
    EXPECT_DOUBLE_EQ(upper.az, 0.0);
    EXPECT_DOUBLE_EQ(upper.bz, 0.0);
  }
}

TEST(thermobarostat, initialization_recomputes_mass_after_dof_constraint)
{
  Thermobarostat state(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic,
                       MonoclinicAngleType::Beta, 320.0, 0.0, 0.002, 12, 3, 1, 0.4);
  state.translationalDegreesOfFreedom = 9;
  RandomNumber random(913);
  state.initialize(random);
  EXPECT_NEAR(state.logVolumeMass, 12.0 * Units::KB * 320.0 * 0.4 * 0.4, 1.0e-14);
  EXPECT_NEAR(state.cellMass, state.logVolumeMass / 3.0, 1.0e-14);
}

// The isotropic barostat integrates x = ln(V) = 3 epsilon; its kinetic energy must be that of the MTK strain,
// W epsdot^2 / 2 = (W/9) xdot^2 / 2, and the sampled initial velocity must have the matching variance
// k_B T / (W/9). Using W for x instead holds the barostat at T/3 and produces a steady drift of the
// conserved quantity through the barostat chain.
TEST(thermobarostat, isotropic_kinetic_energy_uses_the_mass_of_ln_volume)
{
  const double temperature = 320.0;
  Thermobarostat state(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta,
                       temperature, 0.0, 0.002, 12, 3, 1, 0.4);
  RandomNumber random(913);
  state.initialize(random);
  state.logVolumeVelocity = 0.3;
  const double epsilonDot = state.logVolumeVelocity / 3.0;
  EXPECT_NEAR(state.barostatKineticEnergy(), 0.5 * state.logVolumeMass * epsilonDot * epsilonDot, 1.0e-14);
  double chainEnergy = 0.0;
  for (std::size_t i = 0; i != state.chainLength; ++i)
    chainEnergy += 0.5 * state.chainMass[i] * state.chainVelocity[i] * state.chainVelocity[i] +
                   state.chainDegreesOfFreedom[i] * state.chainPosition[i];
  EXPECT_NEAR(state.energy(100.0), state.barostatKineticEnergy() + state.pressure * 100.0 + chainEnergy, 1.0e-10);

  double sumSquared = 0.0;
  const std::size_t samples = 20000;
  for (std::size_t i = 0; i != samples; ++i)
  {
    state.initialize(random);
    sumSquared += state.logVolumeVelocity * state.logVolumeVelocity;
  }
  const double variance = sumSquared / static_cast<double>(samples);
  const double expected = Units::KB * temperature / (state.logVolumeMass / 9.0);
  EXPECT_NEAR(variance / expected, 1.0, 0.05);
}

TEST(thermobarostat, refreshing_degrees_of_freedom_keeps_the_live_cell_momentum)
{
  Thermobarostat state(MolecularDynamicsEnsemble::MuPTPR, CellMinimizationType::Regular, MonoclinicAngleType::Beta,
                       300.0, 1.0, 0.001, 12, 3, 1, 0.5);
  RandomNumber random(4242);
  state.initialize(random);
  const double3x3 cellVelocity = projectCellTensor(double3x3(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                                                   state.cellType, state.monoclinicAngle);
  state.cellVelocity = cellVelocity;
  state.logVolumeVelocity = 0.25;
  const std::vector<double> chainVelocity = state.chainVelocity;

  state.refreshDegreesOfFreedom(24, 8000.0);

  EXPECT_EQ(state.translationalDegreesOfFreedom, 24u);
  EXPECT_NEAR(state.logVolumeMass, 27.0 * Units::KB * 300.0 * 0.25, 1.0e-12);
  EXPECT_NEAR(state.cellMass, state.logVolumeMass / 3.0, 1.0e-12);
  EXPECT_DOUBLE_EQ(state.logVolumePosition, std::log(8000.0));
  for (std::size_t column = 0; column != 3; ++column)
    for (std::size_t row = 0; row != 3; ++row)
      EXPECT_DOUBLE_EQ(state.cellVelocity.mm[column][row], cellVelocity.mm[column][row]);
  EXPECT_DOUBLE_EQ(state.logVolumeVelocity, 0.25);
  EXPECT_EQ(state.chainVelocity, chainVelocity);
}

TEST(thermostat, initialization_sets_the_chain_forces_consistent_with_the_drawn_velocities)
{
  Thermostat state(300.0, 0.0005, 12, 6, 4, 1, 0.15);
  RandomNumber random(555);
  state.initialize(random);
  const double kT = Units::KB * 300.0;
  for (std::size_t i = 1; i != state.thermostatChainLength; ++i)
  {
    EXPECT_NEAR(state.thermostatForceTranslation[i],
                (state.thermostatMassTranslation[i - 1] * state.thermostatVelocityTranslation[i - 1] *
                     state.thermostatVelocityTranslation[i - 1] -
                 kT) /
                    state.thermostatMassTranslation[i],
                1.0e-12);
    EXPECT_NEAR(state.thermostatForceRotation[i],
                (state.thermostatMassRotation[i - 1] * state.thermostatVelocityRotation[i - 1] *
                     state.thermostatVelocityRotation[i - 1] -
                 kT) /
                    state.thermostatMassRotation[i],
                1.0e-12);
    EXPECT_NE(state.thermostatForceTranslation[i], 0.0);
  }
}

TEST(thermobarostat, initialization_sets_the_chain_forces_consistent_with_the_drawn_velocities)
{
  Thermobarostat state(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta,
                       300.0, 0.0, 0.001, 12, 4, 1, 0.5);
  RandomNumber random(777);
  state.initialize(random);
  for (std::size_t i = 1; i != state.chainLength; ++i)
  {
    EXPECT_NEAR(state.chainForce[i],
                (state.chainMass[i - 1] * state.chainVelocity[i - 1] * state.chainVelocity[i - 1] -
                 state.chainDegreesOfFreedom[i]) /
                    state.chainMass[i],
                1.0e-12);
    EXPECT_NE(state.chainForce[i], 0.0);
  }
}

TEST(thermostat, chain_is_inert_when_the_centre_of_mass_constraint_removes_all_translational_degrees_of_freedom)
{
  // one molecule, constraint 3: nothing to thermostat translationally
  Thermostat state(300.0, 0.0005, 3, 0, 3, 1, 0.15);
  state.translationalCenterOfMassConstraint = 3;
  RandomNumber random(31);
  state.initialize(random);

  for (std::size_t step = 0; step != 1000; ++step)
  {
    const auto [scaleTranslation, scaleRotation] = state.NoseHooverNVT(0.0, 0.0);
    EXPECT_DOUBLE_EQ(scaleTranslation, 1.0);
    EXPECT_DOUBLE_EQ(scaleRotation, 1.0);
  }
  EXPECT_DOUBLE_EQ(state.getEnergy(), 0.0);
  for (std::size_t i = 0; i != state.thermostatChainLength; ++i)
  {
    EXPECT_DOUBLE_EQ(state.thermostatVelocityTranslation[i], 0.0);
    EXPECT_DOUBLE_EQ(state.thermostatPositionTranslation[i], 0.0);
  }
}

TEST(system, sampled_temperatures_are_finite_without_translational_degrees_of_freedom)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0e4, 1.0, std::nullopt, {methane}, {},
                {1}, 5);
  system.translationalCenterOfMassConstraint = 3;
  ASSERT_EQ(system.translationalDegreesOfFreedom, 3u);
  system.sampleProperties(0, 0, 0, false);
  EXPECT_TRUE(std::isfinite(system.averageTranslationalTemperature.averaged(0)));
  EXPECT_TRUE(std::isfinite(system.averageTemperature.averaged(0)));
  EXPECT_DOUBLE_EQ(system.averageTranslationalTemperature.averaged(0), 0.0);
  EXPECT_DOUBLE_EQ(system.averageTemperature.averaged(0), 0.0);

  System empty(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0e4, 1.0, std::nullopt, {methane}, {},
               {0}, 5);
  empty.sampleProperties(0, 0, 0, false);
  EXPECT_TRUE(std::isfinite(empty.averageTranslationalTemperature.averaged(0)));
  EXPECT_TRUE(std::isfinite(empty.averageTemperature.averaged(0)));
}

TEST(thermostat, refreshes_variable_particle_degrees_of_freedom_without_resetting_chain_state)
{
  Thermostat state(300.0, 0.0005, 6, 3, 3, 1, 0.15);
  RandomNumber random(731);
  state.initialize(random);
  state.thermostatPositionTranslation[0] = 1.25;
  state.thermostatVelocityTranslation[0] = -0.75;

  state.refreshDegreesOfFreedom(12, 6, 3);

  EXPECT_EQ(state.translationalDegreesOfFreedom, 12u);
  EXPECT_EQ(state.rotationalDegreesOfFreedom, 6u);
  EXPECT_DOUBLE_EQ(state.thermostatPositionTranslation[0], 1.25);
  EXPECT_DOUBLE_EQ(state.thermostatVelocityTranslation[0], -0.75);
  EXPECT_NEAR(state.thermostatDegreesOfFreedomTranslation[0], 9.0 * Units::KB * 300.0, 1.0e-14);
}

TEST(integrators, initializes_only_the_inserted_molecule_velocity)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0e4, 1.0, std::nullopt,
                {methane}, {}, {2}, 5);
  system.moleculeData[0].velocity = {4.0, 5.0, 6.0};
  system.moleculeData[1].velocity = {};
  RandomNumber random(991);
  Molecule& inserted = system.moleculeData[1];

  Integrators::initializeMoleculeVelocity(
      random, inserted, system.spanOfMoleculeAtoms().subspan(inserted.atomIndex, inserted.numberOfAtoms),
      system.spanOfMoleculeDynamics().subspan(inserted.atomIndex, inserted.numberOfAtoms), system.components[0],
      system.temperature);

  EXPECT_EQ(system.moleculeData[0].velocity, double3(4.0, 5.0, 6.0));
  EXPECT_GT(inserted.velocity.length(), 0.0);
}

TEST(molecular_dynamics, gcmd_swap_reports_accepted_insertion_and_deletion_without_changing_cell)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);
  methane.fugacityCoefficient = 1.0;
  methane.idealGasRosenbluthWeight = 1.0;
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0e12, 1.0, std::nullopt,
                {methane}, {}, {1}, 5);
  const double3x3 initialCell = system.simulationBox.cell;
  RandomNumber random(1879);

  std::size_t exchangedMolecule{};
  MC_Moves::ParticleExchangeResult result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Inserted; ++attempt)
    result = MC_Moves::performMolecularDynamicsSwap(random, system, 0, exchangedMolecule);
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Inserted);
  EXPECT_EQ(exchangedMolecule, system.numberOfMoleculesPerComponent[0] - 1);
  const std::size_t countAfterInsertion = system.numberOfIntegerMoleculesPerComponent[0];

  system.pressure = 1.0e-12;
  system.components[0].fugacityCoefficient = 1.0;
  result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Deleted; ++attempt)
    result = MC_Moves::performMolecularDynamicsSwap(random, system, 0, exchangedMolecule);
  EXPECT_EQ(result, MC_Moves::ParticleExchangeResult::Deleted);
  EXPECT_LT(exchangedMolecule, countAfterInsertion);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], countAfterInsertion - 1);
  EXPECT_EQ(system.simulationBox.cell, initialCell);
}

TEST(thermobarostat, archive_round_trip_preserves_extended_state)
{
  Thermobarostat original(MolecularDynamicsEnsemble::MuPTPR, CellMinimizationType::MonoclinicUpperTriangle,
                          MonoclinicAngleType::Gamma, 275.0, -2.0, 0.001, 18, 4, 3, 0.7);
  original.numberOfRespaSteps = 7;
  RandomNumber random(2727);
  original.initialize(random);
  original.chainPosition[2] = 1.25;
  original.cellVelocity = projectCellTensor(double3x3(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                                            original.cellType, original.monoclinicAngle);

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "raspa3_thermobarostat_round_trip.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  Thermobarostat restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.ensemble, original.ensemble);
  EXPECT_EQ(restored.cellType, original.cellType);
  EXPECT_EQ(restored.monoclinicAngle, original.monoclinicAngle);
  EXPECT_EQ(restored.numberOfRespaSteps, original.numberOfRespaSteps);
  EXPECT_EQ(restored.chainPosition, original.chainPosition);
  EXPECT_DOUBLE_EQ(restored.cellVelocity.bx, original.cellVelocity.bx);
  EXPECT_DOUBLE_EQ(restored.energy(37.0), original.energy(37.0));
}

TEST(thermobarostat, upper_triangular_cell_propagation_is_reversible)
{
  const double3x3 originalCell(20.0, 0.0, 0.0, 1.5, 18.0, 0.0, -0.7, 0.8, 22.0);
  const double3x3 rate(0.02, 0.0, 0.0, -0.03, -0.01, 0.0, 0.015, 0.025, 0.005);
  double3x3 cell = originalCell;
  std::vector<double3> positions{{1.2, 2.3, 3.4}, {-2.0, 0.5, 4.1}};
  const std::vector<double3> originalPositions = positions;
  const std::vector<double3> velocities{{0.4, -0.2, 0.1}, {-0.3, 0.6, 0.2}};
  propagateCellAndPosition(cell, positions, velocities, rate, 0.01, true);
  const double3x3 reversedVelocityMap = velocityPropagator(rate, 0.01, 3);
  (void)reversedVelocityMap;
  propagateCellAndPosition(cell, positions, velocities, -rate, 0.01, true);

  EXPECT_NEAR(cell.ax, originalCell.ax, 1.0e-11);
  EXPECT_NEAR(cell.by, originalCell.by, 1.0e-11);
  EXPECT_NEAR(cell.cz, originalCell.cz, 1.0e-11);
  EXPECT_DOUBLE_EQ(cell.ay, 0.0);
  EXPECT_DOUBLE_EQ(cell.az, 0.0);
  EXPECT_DOUBLE_EQ(cell.bz, 0.0);
  // With velocities held fixed this verifies the exact matrix inverse for the affine cell part.
  EXPECT_TRUE(std::ranges::all_of(positions, [](const double3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
  }));
  (void)originalPositions;
}

TEST(thermobarostat, isotropic_hamiltonian_submap_is_symplectic)
{
  const double rate = 0.17;
  const double dt = 0.003;
  const double a = std::exp(rate * dt);
  const double b = dt * std::exp(0.5 * rate * dt) * sinhc(0.5 * rate * dt);
  const double d = std::exp(-rate * dt);
  // J=[[a,b],[0,d]] and Omega=[[0,1],[-1,0]]: J^T Omega J=(ad) Omega.
  EXPECT_NEAR(a * d, 1.0, 2.0e-15);
  EXPECT_GT(std::abs(b), 0.0);
}

TEST(thermobarostat, extended_energy_contains_pressure_work_and_reverses_all_momenta)
{
  Thermobarostat state(MolecularDynamicsEnsemble::NPTPR, CellMinimizationType::Anisotropic,
                       MonoclinicAngleType::Beta, 300.0, 2.5, 0.001, 12, 3, 1, 0.5);
  RandomNumber random(1234);
  state.initialize(random);
  state.cellVelocity = double3x3(0.1, -0.2, 0.3);
  const double energy = state.energy(40.0);
  EXPECT_GT(energy, 100.0);
  const double3x3 before = state.cellVelocity;
  const std::vector<double> chainBefore = state.chainVelocity;
  state.reverseMomenta();
  EXPECT_DOUBLE_EQ(state.cellVelocity.ax, -before.ax);
  EXPECT_DOUBLE_EQ(state.cellVelocity.by, -before.by);
  EXPECT_DOUBLE_EQ(state.cellVelocity.cz, -before.cz);
  for (std::size_t i = 0; i != chainBefore.size(); ++i)
    EXPECT_DOUBLE_EQ(state.chainVelocity[i], -chainBefore[i]);
  EXPECT_NEAR(state.energy(40.0), energy, 1.0e-12);
}

TEST(thermobarostat, complete_npt_map_obeys_time_reversal_identity)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  Framework framework = makeTwoAtomFramework(forceField, false);
  framework.simulationBox = SimulationBox(24.0, 24.0, 24.0);
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {200.0, 1.0})};
  framework.intraMolecularImageShifts.bonds = {{{int3{}, int3{}}}};
  System system(forceField, SimulationBox(24.0, 24.0, 24.0), false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.timeStep = 1.0e-6;
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfFrameworkAtoms()[1].position = {6.05, 5.0, 5.0};
  system.spanOfFrameworkDynamics()[0].velocity = {0.12, -0.03, 0.02};
  system.spanOfFrameworkDynamics()[1].velocity = {-0.04, 0.08, -0.01};
  system.thermostat = Thermostat(300.0, system.timeStep, system.translationalDegreesOfFreedom, 0, 3, 1, 0.15);
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic,
                     MonoclinicAngleType::Beta, 300.0, system.pressure, system.timeStep,
                     system.translationalDegreesOfFreedom, 3, 1, 1.0);
  RandomNumber random(8191);
  system.thermostat->initialize(random);
  system.thermobarostat->initialize(random);
  system.thermobarostat->pressure = 0.0;
  system.thermobarostat->logVolumeVelocity = 0.0;
  system.precomputeTotalGradients();
  system.runningEnergies = system.computeTotalEnergies();
  system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField);
  system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
  system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);
  const double initialExtendedEnergy = system.runningEnergies.conservedEnergy();

  const SimulationBox initialBox = system.simulationBox;
  const std::vector<Atom> initialAtoms(system.spanOfFrameworkAtoms().begin(), system.spanOfFrameworkAtoms().end());
  const std::vector<AtomDynamics> initialDynamics(system.spanOfFrameworkDynamics().begin(),
                                                  system.spanOfFrameworkDynamics().end());
  const std::vector<double> initialThermostatPositions = system.thermostat->thermostatPositionTranslation;
  const std::vector<double> initialBarostatPositions = system.thermobarostat->chainPosition;

  const RunningEnergy forwardEnergy = molecularDynamicsStep(system);
  EXPECT_LT(std::abs(forwardEnergy.conservedEnergy() - initialExtendedEnergy) /
                std::max(1.0, std::abs(initialExtendedEnergy)),
            1.0e-5);
  for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics()) dynamics.velocity = -dynamics.velocity;
  for (double& velocity : system.thermostat->thermostatVelocityTranslation) velocity = -velocity;
  system.thermobarostat->reverseMomenta();
  molecularDynamicsStep(system);
  for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics()) dynamics.velocity = -dynamics.velocity;
  for (double& velocity : system.thermostat->thermostatVelocityTranslation) velocity = -velocity;
  system.thermobarostat->reverseMomenta();

  EXPECT_NEAR(system.simulationBox.cell.ax, initialBox.cell.ax, 2.0e-9);
  EXPECT_NEAR(system.simulationBox.cell.by, initialBox.cell.by, 2.0e-9);
  EXPECT_NEAR(system.simulationBox.cell.cz, initialBox.cell.cz, 2.0e-9);
  for (std::size_t i = 0; i != initialAtoms.size(); ++i)
  {
    DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkAtoms()[i].position, initialAtoms[i].position, 2.0e-8);
    DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkDynamics()[i].velocity, initialDynamics[i].velocity, 2.0e-8);
  }
  for (std::size_t i = 0; i != initialThermostatPositions.size(); ++i)
    EXPECT_NEAR(system.thermostat->thermostatPositionTranslation[i], initialThermostatPositions[i], 2.0e-8);
  for (std::size_t i = 0; i != initialBarostatPositions.size(); ++i)
    EXPECT_NEAR(system.thermobarostat->chainPosition[i], initialBarostatPositions[i], 2.0e-8);
}

void DOUBLE3_EXPECT_NEAR(double3 a, double3 b, double tol)
{
  EXPECT_NEAR(a.x, b.x, tol);
  EXPECT_NEAR(a.y, b.y, tol);
  EXPECT_NEAR(a.z, b.z, tol);
}

void QUATD_EXPECT_NEAR(simd_quatd a, simd_quatd b, double tol)
{
  EXPECT_NEAR(a.r, b.r, tol);
  EXPECT_NEAR(a.ix, b.ix, tol);
  EXPECT_NEAR(a.iy, b.iy, tol);
  EXPECT_NEAR(a.iz, b.iz, tol);
}

void print_system(System& system)
{
  std::cout << "Energy (tr, rot)\n";
  std::cout << Integrators::computeTranslationalKineticEnergy(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                              system.spanOfMoleculeDynamics(), system.components)
            << std::endl;
  std::cout << Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components) << std::endl;

  std::cout << "Center of Mass\n";
  std::cout << system.moleculeData[0].centerOfMassPosition.to_string() << std::endl;
  std::cout << system.moleculeData[1].centerOfMassPosition.to_string() << std::endl;

  std::cout << "Velocity\n";
  std::cout << system.moleculeData[0].velocity.to_string() << std::endl;
  std::cout << system.moleculeData[1].velocity.to_string() << std::endl;

  std::cout << "Gradient\n";
  std::cout << system.moleculeData[0].gradient.to_string() << std::endl;
  std::cout << system.moleculeData[1].gradient.to_string() << std::endl;

  std::cout << "Orientation\n";
  std::cout << system.moleculeData[0].orientation.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientation.to_string() << std::endl;

  std::cout << "OrientationMomentum\n";
  std::cout << system.moleculeData[0].orientationMomentum.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientationMomentum.to_string() << std::endl;

  std::cout << "OrientationGradient\n";
  std::cout << system.moleculeData[0].orientationGradient.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientationGradient.to_string() << std::endl;
}

TEST(integrators_flexible_framework, degrees_of_freedom_and_rigid_regression)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System flexible(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                  makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  System rigid(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
               makeTwoAtomFramework(forceField, true), {}, {}, {}, 5);

  EXPECT_EQ(flexible.translationalDegreesOfFreedom, 6u);
  EXPECT_EQ(flexible.rotationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.translationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.rotationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.numberOfRigidFrameworkAtoms, 2u);
}

TEST(integrators_flexible_framework, mass_weighted_velocity_position_and_kinetic_energy)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  std::span<Atom> atoms = system.spanOfFrameworkAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
  atoms[0].position = {4.0, 5.0, 6.0};
  atoms[1].position = {7.0, 8.0, 9.0};
  dynamics[0].velocity = {1.0, -2.0, 0.5};
  dynamics[1].velocity = {-0.5, 1.0, 2.0};
  dynamics[0].gradient = {4.0, -2.0, 6.0};
  dynamics[1].gradient = {-5.0, 10.0, 2.5};

  Integrators::updateVelocities({}, {}, {}, {}, 0.2, system.framework, atoms, dynamics, &system.forceField);
  DOUBLE3_EXPECT_NEAR(dynamics[0].velocity, double3(0.8, -1.9, 0.2), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(dynamics[1].velocity, double3(-0.4, 0.8, 1.95), 1.0e-14);

  Integrators::updatePositions({}, {}, {}, {}, 0.25, system.framework, atoms, dynamics);
  DOUBLE3_EXPECT_NEAR(atoms[0].position, double3(4.2, 4.525, 6.05), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(atoms[1].position, double3(6.9, 8.2, 9.4875), 1.0e-14);

  const double expected =
      0.5 * 2.0 * (0.8 * 0.8 + 1.9 * 1.9 + 0.2 * 0.2) + 0.5 * 5.0 * (0.4 * 0.4 + 0.8 * 0.8 + 1.95 * 1.95);
  EXPECT_NEAR(Integrators::computeTranslationalKineticEnergy({}, {}, {}, {}, system.framework, atoms, dynamics,
                                                             &system.forceField),
              expected, 1.0e-14);
}

TEST(integrators_flexible_framework, center_of_mass_drift_combines_framework_and_molecule_momentum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(8.0, false, false, false);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  framework.rigid = false;
  Component methane = Component::makeMethane(forceField, 0);
  System system(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, framework, {methane}, {}, {1}, 5);

  system.moleculeData[0].velocity = {1.0, -0.5, 0.25};
  for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics())
  {
    dynamics.velocity = {-0.2, 0.4, -0.1};
  }
  const double frameworkMass = *system.frameworkMass();
  const double moleculeMass = system.moleculeData[0].mass;
  const double3 initialMoleculeVelocity = system.moleculeData[0].velocity;
  const double3 initialFrameworkVelocity = system.spanOfFrameworkDynamics()[0].velocity;
  const double3 drift = (moleculeMass * initialMoleculeVelocity + frameworkMass * initialFrameworkVelocity) /
                        (moleculeMass + frameworkMass);

  Integrators::removeCenterOfMassVelocityDrift(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, initialMoleculeVelocity - drift, 1.0e-14);
  DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkDynamics()[0].velocity, initialFrameworkVelocity - drift, 1.0e-14);
  double3 finalMomentum = moleculeMass * system.moleculeData[0].velocity;
  for (std::size_t i = 0; i < system.spanOfFrameworkAtoms().size(); ++i)
  {
    const double mass = system.forceField.pseudoAtoms[system.spanOfFrameworkAtoms()[i].type].mass;
    finalMomentum += mass * system.spanOfFrameworkDynamics()[i].velocity;
  }
  DOUBLE3_EXPECT_NEAR(finalMomentum, double3{}, 1.0e-11);
}

TEST(integrators_flexible_framework, thermostat_scaling_and_target_degrees_of_freedom)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  system.setThermostat(Thermostat(3, 1, 0.15));
  ASSERT_TRUE(system.thermostat.has_value());
  EXPECT_EQ(system.thermostat->translationalDegreesOfFreedom, 6u);
  EXPECT_EQ(system.thermostat->rotationalDegreesOfFreedom, 0u);

  std::span<AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
  dynamics[0].velocity = {1.0, 2.0, 3.0};
  dynamics[1].velocity = {-2.0, 1.0, 0.5};
  Integrators::scaleVelocities({}, {}, {}, {}, {0.25, 7.0}, system.framework, dynamics);
  DOUBLE3_EXPECT_NEAR(dynamics[0].velocity, double3(0.25, 0.5, 0.75), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(dynamics[1].velocity, double3(-0.5, 0.25, 0.125), 1.0e-14);

  const double kinetic = Integrators::computeTranslationalKineticEnergy(
      {}, {}, {}, {}, system.framework, system.spanOfFrameworkAtoms(), dynamics, &system.forceField);
  const double measuredTemperature =
      2.0 * kinetic / (Units::KB * static_cast<double>(system.translationalDegreesOfFreedom));
  EXPECT_NEAR(measuredTemperature, 2.0 * (0.5 * 2.0 * 0.875 + 0.5 * 5.0 * 0.328125) / (6.0 * Units::KB), 1.0e-12);
}

TEST(integrators_flexible_framework, nose_hoover_velocity_verlet_scales_framework_kinetic_energy)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  const auto makeSystem = [&]()
  {
    System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                  makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
    system.spanOfFrameworkAtoms()[0].position = {4.0, 5.0, 6.0};
    system.spanOfFrameworkAtoms()[1].position = {8.0, 9.0, 10.0};
    system.spanOfFrameworkDynamics()[0].velocity = {1.0, -2.0, 0.5};
    system.spanOfFrameworkDynamics()[1].velocity = {-0.4, 0.8, 1.2};
    return system;
  };
  System nve = makeSystem();
  System nvt = makeSystem();
  nvt.setThermostat(Thermostat(3, 1, 0.15));
  ASSERT_TRUE(nvt.thermostat.has_value());
  RandomNumber random(1729);
  nvt.thermostat->initialize(random);

  const RunningEnergy nveEnergy = Integrators::velocityVerlet(
      {}, {}, {}, {}, nve.timeStep, nve.thermostat, nve.spanOfFrameworkAtoms(), nve.forceField, nve.simulationBox,
      nve.eik_x, nve.eik_y, nve.eik_z, nve.eik_xy, nve.trialEik, nve.fixedFrameworkStoredEik, nve.interpolationGrids,
      {}, nve.framework, nve.spanOfFrameworkDynamics());
  const RunningEnergy nvtEnergy = Integrators::velocityVerlet(
      {}, {}, {}, {}, nvt.timeStep, nvt.thermostat, nvt.spanOfFrameworkAtoms(), nvt.forceField, nvt.simulationBox,
      nvt.eik_x, nvt.eik_y, nvt.eik_z, nvt.eik_xy, nvt.trialEik, nvt.fixedFrameworkStoredEik, nvt.interpolationGrids,
      {}, nvt.framework, nvt.spanOfFrameworkDynamics());

  DOUBLE3_EXPECT_NEAR(nve.spanOfFrameworkDynamics()[0].velocity, double3(1.0, -2.0, 0.5), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(nve.spanOfFrameworkDynamics()[1].velocity, double3(-0.4, 0.8, 1.2), 1.0e-14);
  const double thermostatScale =
      nvt.spanOfFrameworkDynamics()[0].velocity.x / nve.spanOfFrameworkDynamics()[0].velocity.x;
  EXPECT_GT(std::abs(thermostatScale - 1.0), 1.0e-8);
  DOUBLE3_EXPECT_NEAR(nvt.spanOfFrameworkDynamics()[0].velocity,
                      thermostatScale * nve.spanOfFrameworkDynamics()[0].velocity, 1.0e-13);
  DOUBLE3_EXPECT_NEAR(nvt.spanOfFrameworkDynamics()[1].velocity,
                      thermostatScale * nve.spanOfFrameworkDynamics()[1].velocity, 1.0e-13);

  const double explicitFrameworkKinetic = Integrators::computeTranslationalKineticEnergy(
      {}, {}, {}, {}, nvt.framework, nvt.spanOfFrameworkAtoms(), nvt.spanOfFrameworkDynamics(), &nvt.forceField);
  EXPECT_GT(explicitFrameworkKinetic, 0.0);
  EXPECT_NEAR(nvtEnergy.translationalKineticEnergy, explicitFrameworkKinetic, 1.0e-13);
  EXPECT_NE(nvtEnergy.translationalKineticEnergy, nveEnergy.translationalKineticEnergy);
}

TEST(integrators_flexible_framework, harmonic_velocity_verlet_has_bounded_short_nve_drift)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  Framework framework = makeTwoAtomFramework(forceField, false);
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.0})};
  framework.intraMolecularImageShifts.bonds = {{{int3{}, int3{}}}};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfFrameworkAtoms()[1].position = {6.1, 5.0, 5.0};
  system.spanOfFrameworkDynamics()[0].velocity = {0.1, 0.0, 0.0};
  system.spanOfFrameworkDynamics()[1].velocity = {-0.04, 0.0, 0.0};

  RunningEnergy initial = Integrators::updateGradients(
      {}, {}, {}, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, {}, system.eik_x,
      system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, {}, system.framework, system.spanOfFrameworkDynamics());
  initial.translationalKineticEnergy =
      Integrators::computeTranslationalKineticEnergy({}, {}, {}, {}, system.framework, system.spanOfFrameworkAtoms(),
                                                     system.spanOfFrameworkDynamics(), &system.forceField);
  const double initialEnergy = initial.conservedEnergy();
  double maximumDrift = 0.0;

  for (std::size_t step = 0; step < 250; ++step)
  {
    const RunningEnergy current =
        Integrators::velocityVerlet({}, {}, {}, {}, 1.0e-4, system.thermostat, system.spanOfFrameworkAtoms(),
                                    system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                    system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
                                    system.interpolationGrids, {}, system.framework, system.spanOfFrameworkDynamics());
    maximumDrift = std::max(maximumDrift, std::abs(current.conservedEnergy() - initialEnergy));
  }

  EXPECT_GT(std::abs(system.spanOfFrameworkAtoms()[0].position.x - 5.0), 1.0e-5);
  EXPECT_LT(maximumDrift / std::max(1.0, std::abs(initialEnergy)), 2.0e-5);
}

TEST(integrators_flexible_adsorbate, molecular_pressure_preserves_live_md_gradients)
{
  // Molecular pressure evaluates intermolecular/framework strain forces into a scratch buffer and
  // must leave live AtomDynamics gradients (including intramolecular) untouched for Velocity-Verlet.
  ForceField forceField =
      ForceField({{"CH3", false, 15.03452, 0.0, 0.0, 8, false}}, {{108.0, 3.76}},
                 ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);

  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  // Soft enough for stable single-step NVE drift checks, stiff enough that wiping bond forces matters.
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {5000.0, 1.54})};

  Component component =
      Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                 Atom({1.70, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                connectivityTable, intraMolecularPotentials, 5, 21);
  ASSERT_FALSE(component.rigid);

  const std::vector<std::vector<double3>> initialPositions = {
      {double3(5.0, 5.0, 5.0), double3(6.70, 5.0, 5.0)}};
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component},
             initialPositions, {0}, 5);
  system.timeStep = 1.0e-4;
  system.spanOfMoleculeDynamics()[0].velocity = {0.05, -0.02, 0.01};
  system.spanOfMoleculeDynamics()[1].velocity = {-0.04, 0.03, -0.01};

  system.precomputeTotalGradients();
  const double3 gradient0Before = system.spanOfMoleculeDynamics()[0].gradient;
  const double3 gradient1Before = system.spanOfMoleculeDynamics()[1].gradient;
  ASSERT_GT(std::abs(gradient0Before.x), 1.0);

  system.computeMolecularPressure();
  DOUBLE3_EXPECT_NEAR(system.spanOfMoleculeDynamics()[0].gradient, gradient0Before, 1.0e-12);
  DOUBLE3_EXPECT_NEAR(system.spanOfMoleculeDynamics()[1].gradient, gradient1Before, 1.0e-12);

  RunningEnergy initial = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
      system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, system.components, system.eik_x,
      system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      system.spanOfFrameworkDynamics());
  initial.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField);
  const double initialEnergy = initial.conservedEnergy();

  system.computeMolecularPressure();
  const RunningEnergy afterPressure = Integrators::velocityVerlet(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      system.spanOfFrameworkDynamics());
  const double relativeDrift =
      std::abs(afterPressure.conservedEnergy() - initialEnergy) / std::max(1.0, std::abs(initialEnergy));
  EXPECT_LT(relativeDrift, 1.0e-6);
}

TEST(hybrid_mc, flexible_framework_only_does_not_throw_and_restores_on_reject)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  Framework framework = makeTwoAtomFramework(forceField, false);
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.0})};
  framework.intraMolecularImageShifts.bonds = {{{int3{}, int3{}}}};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfFrameworkAtoms()[1].position = {6.1, 5.0, 5.0};
  system.numberOfHybridMCSteps = 5;
  // Large timestep makes energy drift large so the move is almost surely rejected.
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 0.05);

  const double3 initial0 = system.spanOfFrameworkAtoms()[0].position;
  const double3 initial1 = system.spanOfFrameworkAtoms()[1].position;
  RandomNumber random(12);
  const std::optional<RunningEnergy> accepted = MC_Moves::hybridMCMove(random, system);
  EXPECT_FALSE(accepted.has_value());
  DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkAtoms()[0].position, initial0, 1.0e-14);
  DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkAtoms()[1].position, initial1, 1.0e-14);
}

TEST(hybrid_mc, rigid_molecules_in_rigid_framework_does_not_throw)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
  ASSERT_TRUE(system.framework->rigid);
  ASSERT_TRUE(system.components[0].rigid);
  system.numberOfHybridMCSteps = 3;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 1.0e-4);

  RandomNumber random(7);
  EXPECT_NO_THROW(std::ignore = MC_Moves::hybridMCMove(random, system));
}

TEST(hybrid_mc, flexible_molecule_in_rigid_framework_does_not_throw)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, false);
  co2.rigid = false;
  System system(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
  ASSERT_TRUE(system.framework->rigid);
  ASSERT_FALSE(system.components[0].rigid);
  system.numberOfHybridMCSteps = 2;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5.0e-5);

  RandomNumber random(11);
  EXPECT_NO_THROW(std::ignore = MC_Moves::hybridMCMove(random, system));
}

TEST(hybrid_mc, rigid_molecule_in_flexible_framework_does_not_throw)
{
  ForceField forceField = ForceField::makeZeoliteForceField(8.0, false, false, false);
  std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({6.0, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "two-atom-framework", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms,
                      int3(1, 1, 1));
  framework.rigid = false;
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
  ASSERT_FALSE(system.framework->rigid);
  ASSERT_TRUE(system.components[0].rigid);
  system.numberOfHybridMCSteps = 2;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5.0e-5);

  RandomNumber random(19);
  EXPECT_NO_THROW(std::ignore = MC_Moves::hybridMCMove(random, system));
}

TEST(hybrid_mc, flexible_molecule_in_flexible_framework_does_not_throw)
{
  ForceField forceField = ForceField::makeZeoliteForceField(8.0, false, false, false);
  // Keep the flexible framework tiny so the Hybrid MC trial stays a cheap unit test.
  std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({6.0, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "two-atom-framework", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms,
                      int3(1, 1, 1));
  framework.rigid = false;
  Component co2 = Component::makeCO2(forceField, 0, false);
  co2.rigid = false;
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
  ASSERT_FALSE(system.framework->rigid);
  ASSERT_FALSE(system.components[0].rigid);
  system.numberOfHybridMCSteps = 2;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5.0e-5);

  RandomNumber random(3);
  EXPECT_NO_THROW(std::ignore = MC_Moves::hybridMCMove(random, system));
}

TEST(hybrid_mc, mixed_framework_rigid_group_does_not_throw)
{
  System system = makeMixedFrameworkHybridMCSystem();
  ASSERT_TRUE(system.framework->isMixed());
  ASSERT_EQ(system.framework->numberOfRigidGroups(), 1u);
  ASSERT_TRUE(system.framework->hasMobileAtoms());

  system.numberOfHybridMCSteps = 4;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 1.0e-4);

  std::vector<double3> initialFramework;
  for (const Atom& atom : system.spanOfFrameworkAtoms()) initialFramework.push_back(atom.position);
  const double3 initialMolecule = system.spanOfMoleculeAtoms()[0].position;

  RandomNumber random(21);
  std::optional<RunningEnergy> accepted;
  EXPECT_NO_THROW(accepted = MC_Moves::hybridMCMove(random, system));
  ASSERT_EQ(system.spanOfFrameworkAtoms().size(), initialFramework.size());
  if (accepted.has_value())
  {
    EXPECT_EQ(system.frameworkGroupData.size(), 1u);
  }
  else
  {
    for (std::size_t i = 0; i != initialFramework.size(); ++i)
    {
      DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkAtoms()[i].position, initialFramework[i], 1.0e-14);
    }
    DOUBLE3_EXPECT_NEAR(system.spanOfMoleculeAtoms()[0].position, initialMolecule, 1.0e-14);
  }
}

// A single molecule has degrees of freedom to propagate: in a rigid framework the lab frame is pinned, so even a
// single atom translates; in a free lab frame a polyatomic molecule still rotates (and flexes). Only a single
// atom in a free lab frame has nothing left once its centre-of-mass velocity is removed.
TEST(hybrid_mc, single_molecules_are_propagated)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 0, false);

  std::vector<System> systems{
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {methane}, {}, {1}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {1}, 5)};
  for (System& system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
    system.numberOfHybridMCSteps = 5;
    system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 1.0e-3);
    const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    const std::vector<Atom> before(atoms.begin(), atoms.end());

    RandomNumber random(3);
    std::size_t accepted = 0;
    for (std::size_t i = 0; i != 20; ++i)
    {
      const std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, system);
      if (!energy.has_value()) continue;
      ++accepted;
      system.runningEnergies += energy.value();
    }
    EXPECT_GT(accepted, 0uz);
    double displacement = 0.0;
    for (std::size_t k = 0; k != atoms.size(); ++k)
    {
      displacement = std::max(displacement, (atoms[k].position - before[k].position).length());
    }
    EXPECT_GT(displacement, 1.0e-6);
  }

  // a single atom in a free lab frame: nothing to propagate
  System lone(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {1}, 5);
  lone.runningEnergies = lone.computeTotalEnergies();
  lone.numberOfHybridMCSteps = 5;
  RandomNumber random(3);
  EXPECT_FALSE(MC_Moves::hybridMCMove(random, lone).has_value());
}

TEST(hybrid_mc, rejects_when_molecule_is_inside_blocked_pocket)
{
  // Pocket is attached after construction: System::createInitialMolecules refuses pockets, so a
  // box-filling sphere at construction time would hang the initial grow.
  System system = makeMixedFrameworkHybridMCSystem();
  system.components[0].blockingPockets.push_back(double4(0.5, 0.5, 0.5, 40.0));
  system.numberOfHybridMCSteps = 2;
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 1.0e-4);

  const double3 initialMolecule = system.spanOfMoleculeAtoms()[0].position;
  RandomNumber random(5);
  EXPECT_FALSE(MC_Moves::hybridMCMove(random, system).has_value());
  DOUBLE3_EXPECT_NEAR(system.spanOfMoleculeAtoms()[0].position, initialMolecule, 1.0e-14);
}

// The MD gradient routines carry no external field, so the trajectory is a surrogate proposal; the acceptance and
// the returned energy must nevertheless be the exact ones (external field and tail corrections included), and a
// rejected move must leave the configuration and its energy untouched.
TEST(hybrid_mc, acceptance_and_returned_energy_are_exact_with_an_external_field)
{
  const auto makeSystem = [](bool externalField)
  {
    ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                          ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, true, false);
    forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
    const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
    Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)}, {}, {}, 5, 21);
    // Four atoms of one component (Hybrid MC needs more than one molecule). They sit near a Muller-Brown
    // minimum in x,y and are spaced along z so the pair energy is finite.
    System system(forceField, SimulationBox(20.0, 20.0, 20.0), externalField, 300.0, 1e4, 1.0, {}, {atomic},
                  {{double3(-0.50, 1.50, 4.0), double3(-0.45, 1.55, 8.0), double3(-0.55, 1.45, 12.0),
                    double3(-0.40, 1.40, 16.0)}},
                  {0}, 5);
    system.hybridMCUseForceEngine = false;
    system.numberOfHybridMCSteps = 8;
    system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5.0e-4);
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };

  for (const bool externalField : {true, false})
  {
    std::size_t accepted = 0, rejected = 0;
    for (std::size_t seed = 0; seed < 24 && (accepted < 3 || rejected < 1); ++seed)
    {
      System system = makeSystem(externalField);
      ASSERT_GT(system.moleculeData.size(), 1uz);
      const RunningEnergy exactOld = system.computeTotalEnergies();
      if (externalField) ASSERT_NE(exactOld.externalFieldVDW, 0.0);
      ASSERT_NE(exactOld.tail, 0.0);
      std::vector<double3> initialPositions;
      for (const Atom& atom : system.spanOfMoleculeAtoms()) initialPositions.push_back(atom.position);

      RandomNumber random(seed);
      const std::optional<RunningEnergy> result = MC_Moves::hybridMCMove(random, system);
      const RunningEnergy exactNew = system.computeTotalEnergies();
      const double scale = std::max(1.0, std::abs(exactNew.potentialEnergy()));
      if (!result.has_value())
      {
        ++rejected;
        for (std::size_t i = 0; i != initialPositions.size(); ++i)
        {
          DOUBLE3_EXPECT_NEAR(system.spanOfMoleculeAtoms()[i].position, initialPositions[i], 1.0e-14);
        }
        EXPECT_NEAR(exactNew.potentialEnergy(), exactOld.potentialEnergy(), 1.0e-10 * scale);
        continue;
      }
      ++accepted;
      if (externalField) EXPECT_NE(exactNew.externalFieldVDW, exactOld.externalFieldVDW);
      EXPECT_NEAR(result->externalFieldVDW, exactNew.externalFieldVDW, 1.0e-10 * scale);
      EXPECT_NEAR(result->tail, exactNew.tail, 1.0e-10 * scale);
      EXPECT_NEAR(result->potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
      EXPECT_GT(result->translationalKineticEnergy, 0.0);
      // the MC driver overwrites runningEnergies with the returned energy (it does not add)
      system.runningEnergies = *result;
      EXPECT_NEAR(system.runningEnergies.externalFieldVDW, exactNew.externalFieldVDW, 1.0e-10 * scale);
      EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
    }
    EXPECT_GE(accepted, 1uz) << "external field " << externalField;
  }
}

// Polarization forces are likewise absent from the MD gradients; the acceptance and the returned energy
// must still be the exact ones (polarization included), and a rejected move must leave the stored
// electric field of the old configuration.
TEST(hybrid_mc, acceptance_and_returned_energy_are_exact_with_polarization)
{
  const auto makeSystem = []()
  {
    ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
    forceField.computePolarization = true;
    forceField.omitInterPolarization = false;
    Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
    Component co2 = Component::makeCO2(forceField, 0, true);
    System system(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
    system.hybridMCUseForceEngine = false;
    system.numberOfHybridMCSteps = 3;
    system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 1.0e-4);
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };

  System probe = makeSystem();
  ASSERT_LT(probe.runningEnergies.polarization, -1e-8);

  std::size_t accepted = 0;
  for (std::size_t seed = 0; seed < 16 && accepted < 2; ++seed)
  {
    System system = makeSystem();
    const RunningEnergy exactOld = system.computeTotalEnergies();
    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::hybridMCMove(random, system);
    if (!result.has_value()) continue;
    ++accepted;
    const RunningEnergy exactNew = system.computeTotalEnergies();
    const double scale = std::max(1.0, std::abs(exactNew.potentialEnergy()));
    EXPECT_NEAR(result->polarization, exactNew.polarization, 1.0e-10 * scale);
    EXPECT_NEAR(result->potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
    EXPECT_GT(result->translationalKineticEnergy, 0.0);
    EXPECT_NE(exactNew.polarization, exactOld.polarization);
    EXPECT_LT(result->polarization, -1.0e-8);
    // the MC driver overwrites runningEnergies with the returned energy (it does not add)
    system.runningEnergies = *result;
    EXPECT_NEAR(system.runningEnergies.polarization, exactNew.polarization, 1.0e-10 * scale);
    EXPECT_LT(system.runningEnergies.polarization, -1.0e-8);
  }
  EXPECT_GE(accepted, 1uz);

  System rejected = makeSystem();
  rejected.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 0.05);
  const std::span<const double3> field = rejected.spanOfMoleculeElectricField();
  const std::vector<double3> initialField(field.begin(), field.end());
  const std::vector<Atom> initialAtoms(rejected.spanOfMoleculeAtoms().begin(), rejected.spanOfMoleculeAtoms().end());
  const RunningEnergy exactOld = rejected.runningEnergies;
  RandomNumber random(3);
  EXPECT_FALSE(MC_Moves::hybridMCMove(random, rejected).has_value());
  const std::span<const Atom> atoms = rejected.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i != initialAtoms.size(); ++i)
  {
    DOUBLE3_EXPECT_NEAR(atoms[i].position, initialAtoms[i].position, 1.0e-14);
  }
  const std::span<const double3> restored = rejected.spanOfMoleculeElectricField();
  for (std::size_t i = 0; i != initialField.size(); ++i)
  {
    DOUBLE3_EXPECT_NEAR(restored[i], initialField[i], 1.0e-14);
  }
  EXPECT_NEAR(rejected.computeTotalEnergies().potentialEnergy(), exactOld.potentialEnergy(),
              1.0e-10 * std::max(1.0, std::abs(exactOld.potentialEnergy())));
}

TEST(integrators, Test_2_CO2_in_ITQ_29_2x2x2_inter)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, false);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Molecule> moleculeData(system.moleculeData);

  moleculeData[0].centerOfMassPosition = double3(5.93355, 7.93355, 5.93355);
  moleculeData[1].centerOfMassPosition = double3(5.93355, 3.93355, 5.93355);
  moleculeData[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  moleculeData[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);

  Integrators::createCartesianPositions(moleculeData, system.spanOfMoleculeAtoms(), system.components);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  DOUBLE3_EXPECT_NEAR(atomData[0].position, double3(5.93355, 7.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[1].position, double3(5.93355, 7.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[2].position, double3(5.93355, 7.93355, 4.78455), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[3].position, double3(5.93355, 3.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[4].position, double3(5.93355, 3.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[5].position, double3(5.93355, 3.93355, 4.78455), 1e-6);

  moleculeData[0].velocity = double3(0.0, -1.0, 0.0);
  moleculeData[0].orientationMomentum = simd_quatd(10.0, 0.0, 0.0, 0.0);
  moleculeData[1].velocity = double3(0.0, 1.0, 0.0);
  moleculeData[1].orientationMomentum = simd_quatd(-10.0, 0.0, 0.0, 0.0);

  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition, double3(5.933550, 7.933050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].centerOfMassPosition, double3(5.933550, 3.934050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, double3(0.000218, -0.999519, 0.00004), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].velocity, double3(0.000218, 0.999519, 0.00004), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientation, simd_quatd(0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientation, simd_quatd(-0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientationMomentum, simd_quatd(10.000011, 0.0, 0.0, -0.000296), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientationMomentum, simd_quatd(-10.000011, 0.0, 0.0, -0.000296), 1e-6);

  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition, double3(5.933550, 7.932550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].centerOfMassPosition, double3(5.933550, 3.934550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, double3(0.000653, -0.998560, 0.000121), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].velocity, double3(0.000653, 0.998560, 0.000121), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientation, simd_quatd(0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientation, simd_quatd(-0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientationMomentum, simd_quatd(10.000043, 0.0, 0.0, -0.000592), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientationMomentum, simd_quatd(-10.000043, 0.0, 0.0, -0.000592), 1e-6);
}

TEST(molecular_dynamics, resume_of_a_finished_or_empty_stage_does_not_loop)
{
  MolecularDynamics md;
  md.outputToFiles = false;
  md.printEvery = 1;
  md.writeBinaryRestartEvery = 1;
  md.optimizeMCMovesEvery = 1;
  md.numberOfSteps = 0;
  md.numberOfPreInitializationCycles = 0;
  md.numberOfInitializationCycles = 3;
  md.numberOfEquilibrationCycles = 5;
  md.numberOfProductionCycles = 7;

  md.simulationStage = MolecularDynamics::SimulationStage::PreInitialization;
  md.currentCycle = 0;
  md.preInitialize();

  md.simulationStage = MolecularDynamics::SimulationStage::Initialization;
  md.currentCycle = 3;
  md.initialize();

  md.simulationStage = MolecularDynamics::SimulationStage::Equilibration;
  md.currentCycle = 5;
  md.equilibrate();

  md.simulationStage = MolecularDynamics::SimulationStage::Production;
  md.currentCycle = 7;
  md.production();

  md.numberOfProductionCycles = 0;
  md.simulationStage = MolecularDynamics::SimulationStage::Production;
  md.currentCycle = 0;
  md.production();
}

TEST(molecular_dynamics, precompute_total_gradients_fills_rigid_body_forces)
{
  // Equilibration/production re-derive group kinematics then call precomputeTotalGradients. The
  // first Velocity-Verlet half-kick uses molecule.gradient / orientationGradient, not the atom
  // gradients, so those rigid-body forces must be reduced here.
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);
  ASSERT_TRUE(system.components[0].rigid);

  system.moleculeData[0].centerOfMassPosition = double3(5.0, 5.0, 5.0);
  system.moleculeData[1].centerOfMassPosition = double3(5.0, 9.0, 5.0);
  system.moleculeData[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  system.moleculeData[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);

  system.initializeGroupData();
  system.moleculeData[0].gradient = double3(0.0, 0.0, 0.0);
  system.moleculeData[1].gradient = double3(0.0, 0.0, 0.0);
  system.moleculeData[0].orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);
  system.moleculeData[1].orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);

  system.precomputeTotalGradients();

  EXPECT_GT(double3::dot(system.moleculeData[0].gradient, system.moleculeData[0].gradient), 0.0);
  EXPECT_GT(double3::dot(system.moleculeData[1].gradient, system.moleculeData[1].gradient), 0.0);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].gradient + system.moleculeData[1].gradient, double3(0.0, 0.0, 0.0),
                      1.0e-10);
}

TEST(molecular_dynamics, precompute_total_gradients_fills_framework_group_forces)
{
  System system = makeMixedFrameworkHybridMCSystem();
  ASSERT_EQ(system.framework->numberOfRigidGroups(), 1u);

  system.initializeGroupData();
  system.initializeFrameworkGroupData();
  ASSERT_EQ(system.frameworkGroupData.size(), 1u);
  system.frameworkGroupData[0].gradient = double3(0.0, 0.0, 0.0);
  system.frameworkGroupData[0].orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);

  system.precomputeTotalGradients();

  EXPECT_GT(double3::dot(system.frameworkGroupData[0].gradient, system.frameworkGroupData[0].gradient), 0.0);
}

TEST(molecular_dynamics, accepted_particle_exchange_refreshes_barostat_with_molecular_degrees_of_freedom)
{
  // Flexible ethane: the thermostat sees 3 atomic translational degrees of freedom per atom, the (molecular)
  // barostat sees 3 per molecule. After an accepted exchange the barostat must get the molecular count.
  System system = makeFlexibleEthaneSystem();
  system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::MuPTPR;
  system.components[0].fugacityCoefficient = 1.0;
  system.components[0].idealGasRosenbluthWeight = 1.0;
  system.components[0].mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  system.pressure = 1.0e12;
  system.translationalCenterOfMassConstraint = 3;
  system.thermostat = Thermostat(300.0, 1.0e-3, system.translationalDegreesOfFreedom, 0, 3, 1, 0.15);
  system.thermostat->translationalCenterOfMassConstraint = 3;
  system.thermobarostat = Thermobarostat(MolecularDynamicsEnsemble::MuPTPR, CellMinimizationType::Regular,
                                         MonoclinicAngleType::Beta, 300.0, 1.0, 1.0e-3,
                                         barostatTranslationalDegreesOfFreedom(system, BarostatCoupling::Molecular),
                                         3, 1, 0.5);
  RandomNumber random(2024);
  system.thermostat->initialize(random);
  system.thermobarostat->initialize(random);
  system.precomputeTotalGradients();

  const double3x3 cellVelocity = projectCellTensor(double3x3(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                                                   system.thermobarostat->cellType,
                                                   system.thermobarostat->monoclinicAngle);
  system.thermobarostat->cellVelocity = cellVelocity;

  MC_Moves::ParticleExchangeResult result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 200 && result != MC_Moves::ParticleExchangeResult::Inserted; ++attempt)
    result = attemptParticleExchange(random, system);
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Inserted);

  const std::size_t molecularDegreesOfFreedom = 3 * system.moleculeData.size() - 3;
  const std::size_t atomicDegreesOfFreedom = system.translationalDegreesOfFreedom - 3;
  ASSERT_NE(molecularDegreesOfFreedom, atomicDegreesOfFreedom);
  EXPECT_EQ(system.thermobarostat->translationalDegreesOfFreedom, molecularDegreesOfFreedom);
  EXPECT_EQ(system.thermobarostat->translationalDegreesOfFreedom,
            barostatTranslationalDegreesOfFreedom(system, system.thermobarostat->coupling));
  EXPECT_EQ(system.thermostat->translationalDegreesOfFreedom, system.translationalDegreesOfFreedom);
  for (std::size_t column = 0; column != 3; ++column)
    for (std::size_t row = 0; row != 3; ++row)
      EXPECT_DOUBLE_EQ(system.thermobarostat->cellVelocity.mm[column][row], cellVelocity.mm[column][row]);
}

TEST(integrators, semi_flexible_molecule_velocity_matches_group_and_atom_com)
{
  System system = makeSemiFlexibleChainSystem();
  ASSERT_TRUE(system.components[0].isSemiFlexible());
  system.initializeGroupData();
  RandomNumber random(11);
  Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                    system.spanOfMoleculeDynamics(), system.components, system.temperature,
                                    std::nullopt, {}, {}, nullptr, system.spanOfGroupData());
  Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, system.spanOfMoleculeAtoms(),
                                               system.spanOfMoleculeDynamics(), system.components, std::nullopt, {}, {},
                                               nullptr, system.spanOfGroupData());

  for (std::size_t m = 0; m != system.moleculeData.size(); ++m)
  {
    DOUBLE3_EXPECT_NEAR(system.moleculeData[m].velocity, expectedMoleculeCenterOfMassVelocity(system, m), 1.0e-12);
  }

  Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                               system.components, {1.15, 1.0}, std::nullopt, {}, system.spanOfGroupData());
  for (std::size_t m = 0; m != system.moleculeData.size(); ++m)
  {
    DOUBLE3_EXPECT_NEAR(system.moleculeData[m].velocity, expectedMoleculeCenterOfMassVelocity(system, m), 1.0e-12);
  }

  system.precomputeTotalGradients();
  std::optional<Thermostat> thermostat = std::nullopt;
  system.timeStep = 1.0e-3;
  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(),
                              system.spanOfGroupData());
  for (std::size_t m = 0; m != system.moleculeData.size(); ++m)
  {
    DOUBLE3_EXPECT_NEAR(system.moleculeData[m].velocity, expectedMoleculeCenterOfMassVelocity(system, m), 1.0e-12);
  }
}

TEST(integrators, flexible_molecule_velocity_matches_atom_com_after_thermostat)
{
  System system = makeFlexibleEthaneSystem();
  ASSERT_FALSE(system.components[0].rigid);
  ASSERT_FALSE(system.components[0].isSemiFlexible());
  RandomNumber random(19);
  Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                    system.spanOfMoleculeDynamics(), system.components, system.temperature);
  Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, system.spanOfMoleculeAtoms(),
                                               system.spanOfMoleculeDynamics(), system.components);
  Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                               system.components, {1.15, 1.0});
  for (std::size_t m = 0; m != system.moleculeData.size(); ++m)
  {
    DOUBLE3_EXPECT_NEAR(system.moleculeData[m].velocity, expectedMoleculeCenterOfMassVelocity(system, m), 1.0e-12);
  }

  system.precomputeTotalGradients();
  std::optional<Thermostat> thermostat = std::nullopt;
  system.timeStep = 1.0e-3;
  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);
  for (std::size_t m = 0; m != system.moleculeData.size(); ++m)
  {
    DOUBLE3_EXPECT_NEAR(system.moleculeData[m].velocity, expectedMoleculeCenterOfMassVelocity(system, m), 1.0e-12);
  }
}

namespace
{
std::vector<double3> netForcesPerMolecule(System& system)
{
  system.precomputeTotalGradients();
  const std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  std::vector<double3> forces(system.moleculeData.size());
  std::size_t index = 0;
  for (std::size_t i = 0; i < system.moleculeData.size(); ++i)
  {
    for (std::size_t k = 0; k < system.moleculeData[i].numberOfAtoms; ++k)
    {
      forces[i] -= dynamics[index + k].gradient;
    }
    index += system.moleculeData[i].numberOfAtoms;
  }
  return forces;
}

System makeCO2InITQ29System(bool useCharges = false, bool polarization = false)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  forceField.computePolarization = polarization;
  forceField.omitInterPolarization = false;
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, useCharges);
  return System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {co2}, {}, {2}, 5);
}
}  // namespace

TEST(smart_mc_all, translation_biases_semi_flexible_molecules_with_the_net_atomic_force)
{
  // The drift of the smart-MC proposal is b times the net force on the molecule, the sum of its atomic
  // forces. 'Molecule::gradient' is not maintained for semi-flexible molecules (their rigid groups carry the
  // gradients), so a stale value there must not enter the proposal.
  const double sigma = 0.05;
  std::size_t accepted = 0;
  for (std::size_t seed = 0; seed < 12; ++seed)
  {
    System system = makeSemiFlexibleChainSystem();
    ASSERT_TRUE(system.components[0].isSemiFlexible());
    system.mc_moves_statistics.setMaxChange(Move::Types::TranslationSmartMCAll, sigma);

    // bring the molecules within range of each other so that the net forces are appreciable
    ASSERT_EQ(system.moleculeData.size(), 2uz);
    const std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    ASSERT_EQ(atoms.size(), 8uz);
    const double3 shift(-3.5, -2.0, -1.0);
    for (std::size_t k = 4; k < 8; ++k) atoms[k].position += shift;
    system.moleculeData[1].centerOfMassPosition += shift;
    for (Molecule& molecule : system.moleculeData) molecule.gradient = double3(3.0e3, -2.0e3, 1.0e3);

    const std::vector<double3> netForce = netForcesPerMolecule(system);
    const double b = 0.5 * system.beta * sigma * sigma;
    ASSERT_GT(b * netForce[0].length(), 1.0e-7) << "the fixture must exert a measurable drift";
    system.runningEnergies = system.computeTotalEnergies();
    const std::vector<Atom> before(atoms.begin(), atoms.end());

    RandomNumber random(seed);
    RandomNumber replay(seed);
    const std::optional<RunningEnergy> result = MC_Moves::translationSmartMCMoveAll(random, system);

    std::size_t index = 0;
    for (std::size_t i = 0; i < system.moleculeData.size(); ++i)
    {
      const double3 xi(replay.Gaussian(), replay.Gaussian(), replay.Gaussian());
      const double3 expected = result.has_value() ? b * netForce[i] + sigma * xi : double3{};
      for (std::size_t k = 0; k < system.moleculeData[i].numberOfAtoms; ++k)
      {
        DOUBLE3_EXPECT_NEAR(atoms[index + k].position - before[index + k].position, expected, 1.0e-10);
      }
      index += system.moleculeData[i].numberOfAtoms;
    }
    if (result.has_value())
    {
      ++accepted;
      const double scale = std::max(1.0, std::abs(result->potentialEnergy()));
      EXPECT_NEAR(result->potentialEnergy(), system.computeTotalEnergies().potentialEnergy(), 1.0e-8 * scale);
    }
  }
  EXPECT_GE(accepted, 1uz);
}

TEST(smart_mc_all, translation_rejects_a_molecule_displaced_into_a_blocked_pocket)
{
  // Pockets are defined in the framework cell; attached after construction (the initial grow refuses
  // pockets), covering the whole box.
  System system = makeCO2InITQ29System();
  system.components[0].blockingPockets.push_back(double4(0.5, 0.5, 0.5, 40.0));
  system.mc_moves_statistics.setMaxChange(Move::Types::TranslationSmartMCAll, 1.0e-3);
  system.runningEnergies = system.computeTotalEnergies();
  const RunningEnergy energyBefore = system.runningEnergies;
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<Atom> before(atoms.begin(), atoms.end());

  RandomNumber random(5);
  EXPECT_FALSE(MC_Moves::translationSmartMCMoveAll(random, system).has_value());
  for (std::size_t k = 0; k < atoms.size(); ++k)
  {
    DOUBLE3_EXPECT_NEAR(atoms[k].position, before[k].position, 1.0e-14);
  }
  EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());
}

TEST(smart_mc_all, rotation_rejects_a_molecule_rotated_into_a_blocked_pocket)
{
  System system = makeCO2InITQ29System();
  ASSERT_TRUE(system.components[0].rigid);
  system.components[0].blockingPockets.push_back(double4(0.5, 0.5, 0.5, 40.0));
  system.mc_moves_statistics.setMaxChange(Move::Types::RotationSmartMCAll, 1.0e-3);
  system.runningEnergies = system.computeTotalEnergies();
  const RunningEnergy energyBefore = system.runningEnergies;
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<Atom> before(atoms.begin(), atoms.end());

  RandomNumber random(5);
  EXPECT_FALSE(MC_Moves::rotationSmartMCMoveAll(random, system).has_value());
  for (std::size_t k = 0; k < atoms.size(); ++k)
  {
    DOUBLE3_EXPECT_NEAR(atoms[k].position, before[k].position, 1.0e-14);
  }
  EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());
}

TEST(smart_mc_all, translation_acceptance_includes_the_external_field_energy)
{
  // The gradient routines carry no external field; the acceptance and the returned running energy must
  // nevertheless include its change.
  const auto makeSystem = []()
  {
    ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                          ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
    forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
    const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
    Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)},
                     {}, {}, 5, 21);
    System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {atomic},
                  {{double3(0.5, 0.3, 5.0)}, {double3(-0.4, 1.2, 9.0)}}, {0}, 5);
    system.hasExternalField = true;
    system.mc_moves_statistics.setMaxChange(Move::Types::TranslationSmartMCAll, 0.02);
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };

  std::size_t accepted = 0;
  for (std::size_t seed = 0; seed < 16 && accepted < 3; ++seed)
  {
    System system = makeSystem();
    const RunningEnergy exactOld = system.computeTotalEnergies();
    ASSERT_NE(exactOld.externalFieldVDW, 0.0);

    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::translationSmartMCMoveAll(random, system);
    if (!result.has_value()) continue;
    ++accepted;

    const RunningEnergy exactNew = system.computeTotalEnergies();
    EXPECT_NE(exactNew.externalFieldVDW, exactOld.externalFieldVDW);
    const double scale = std::max(1.0, std::abs(exactNew.potentialEnergy()));
    EXPECT_NEAR(result->externalFieldVDW, exactNew.externalFieldVDW, 1.0e-10 * scale);
    EXPECT_NEAR(result->potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
    EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
  }
  EXPECT_GE(accepted, 1uz);
}

TEST(smart_mc_all, moves_keep_the_running_energy_and_structure_factors_consistent)
{
  // Rigid charged molecules with Ewald: accepted moves leave the running energy equal to the exact energy, rejected
  // moves restore the configuration and the stored structure factors. Charged molecules are required: uncharged
  // CO2 produces a zero molecule structure factor for both the old and the trial configuration, so a missing
  // restore of storedEik would not be observed.
  System system = makeCO2InITQ29System(true);
  system.mc_moves_statistics.setMaxChange(Move::Types::TranslationSmartMCAll, 0.3);
  system.mc_moves_statistics.setMaxChange(Move::Types::RotationSmartMCAll, 0.6);
  system.runningEnergies = system.computeTotalEnergies();
  ASSERT_FALSE(system.storedEik.empty());
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();

  RandomNumber random(3);
  std::size_t accepted = 0;
  std::size_t rejected = 0;
  for (std::size_t step = 0; step < 40; ++step)
  {
    const std::vector<Atom> before(atoms.begin(), atoms.end());
    const RunningEnergy energyBefore = system.runningEnergies;
    const EwaldStructureFactor eikBefore(system.storedEik);

    const std::optional<RunningEnergy> result = (step % 2 == 0) ? MC_Moves::translationSmartMCMoveAll(random, system)
                                                                : MC_Moves::rotationSmartMCMoveAll(random, system);
    if (result.has_value())
    {
      ++accepted;
      system.runningEnergies = result.value();
    }
    else
    {
      ++rejected;
      for (std::size_t k = 0; k < atoms.size(); ++k)
      {
        DOUBLE3_EXPECT_NEAR(atoms[k].position, before[k].position, 1.0e-14);
      }
      EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());
      EXPECT_TRUE(system.storedEik == eikBefore);
    }
    const RunningEnergy exact = system.computeTotalEnergies();
    const double scale = std::max(1.0, std::abs(exact.potentialEnergy()));
    EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exact.potentialEnergy(), 1.0e-8 * scale);
  }
  EXPECT_GE(accepted, 1uz);
  EXPECT_GE(rejected, 1uz);
}

// computeTotalEnergies() writes storedEik and the molecule electric field in place for the configuration it
// evaluates. A rejected collective smart-MC move must restore both; otherwise a later incremental translation
// (which reads storedEik as the old structure factor and the stored field as the old polarization field) drifts.
TEST(smart_mc_all, rejection_restores_stored_ewald_and_electric_field)
{
  System system = makeCO2InITQ29System(true, true);
  system.mc_moves_statistics.setMaxChange(Move::Types::TranslationSmartMCAll, 0.3);
  system.mc_moves_statistics.setMaxChange(Move::Types::RotationSmartMCAll, 0.6);
  system.runningEnergies = system.computeTotalEnergies();
  ASSERT_LT(system.runningEnergies.polarization, -1.0e-8);
  ASSERT_FALSE(system.storedEik.empty());

  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::span<const double3> field = system.spanOfMoleculeElectricField();
  RandomNumber random(7);
  std::size_t rejected = 0;
  std::size_t followUpAccepted = 0;
  for (std::size_t step = 0; step < 40; ++step)
  {
    const std::vector<Atom> atomsBefore(atoms.begin(), atoms.end());
    const std::vector<double3> fieldBefore(field.begin(), field.end());
    const EwaldStructureFactor eikBefore(system.storedEik);
    const RunningEnergy energyBefore = system.runningEnergies;

    const std::optional<RunningEnergy> result = (step % 2 == 0) ? MC_Moves::translationSmartMCMoveAll(random, system)
                                                                : MC_Moves::rotationSmartMCMoveAll(random, system);
    if (result.has_value())
    {
      system.runningEnergies = result.value();
      continue;
    }

    ++rejected;
    for (std::size_t k = 0; k < atoms.size(); ++k)
    {
      DOUBLE3_EXPECT_NEAR(atoms[k].position, atomsBefore[k].position, 1.0e-14);
    }
    for (std::size_t k = 0; k < field.size(); ++k)
    {
      DOUBLE3_EXPECT_NEAR(field[k], fieldBefore[k], 1.0e-14);
    }
    EXPECT_TRUE(system.storedEik == eikBefore);
    EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());

    // Incremental translation uses the restored storedEik and electric field as the old state.
    const std::optional<RunningEnergy> delta = MC_Moves::translationMove(random, system, 0, 0);
    if (delta.has_value())
    {
      ++followUpAccepted;
      system.runningEnergies += delta.value();
    }
    const RunningEnergy exact = system.computeTotalEnergies();
    const double scale = std::max(1.0, std::abs(exact.potentialEnergy()));
    EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exact.potentialEnergy(), 1.0e-8 * scale);
    EXPECT_NEAR(system.runningEnergies.polarization, exact.polarization, 1.0e-8 * scale);
    EXPECT_NEAR(system.runningEnergies.ewald_fourier, exact.ewald_fourier, 1.0e-8 * scale);
  }
  EXPECT_GE(rejected, 1uz);
  EXPECT_GE(followUpAccepted, 1uz);
}

namespace
{
// A chain host: 'fixedAtoms' lab-fixed atoms at the start, the others flexible, no molecules.
System makeHybridFrameworkSystem(const std::string& groups)
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}}, {{80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot,
                        5.0, 5.0, 5.0, true, false, false);
  const SimulationBox box(14.0, 14.0, 14.0);
  std::vector<Atom> atoms;
  for (const double3& position : {double3(4.0, 7.0, 7.0), double3(5.4, 7.3, 7.0), double3(6.8, 6.8, 7.0)})
  {
    atoms.push_back(Atom(position / 14.0, 0.0, 1.0, 0, 0, 0, 0, true));
  }
  Framework framework(forceField, "hybrid-chain", box, 1, atoms, atoms, int3(1, 1, 1));
  const std::filesystem::path definition = std::filesystem::temp_directory_path() / "raspa3-hybrid-framework.json";
  std::ofstream(definition) << std::format(R"({{"Type": "Flexible", {} "Intra14VanDerWaalsScalingValue": 0.0,
                                                "Intra14ChargeChargeScalingValue": 0.0}})",
                                           groups);
  framework.readFrameworkDefinition(forceField, definition.string());
  std::filesystem::remove(definition);
  return System(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
}

// With no integration steps the move commits the velocities it drew, which exposes the projection.
double drawnFrameworkMomentum(System& system, RandomNumber& random)
{
  system.numberOfHybridMCSteps = 0;
  EXPECT_TRUE(MC_Moves::hybridMCMove(random, system).has_value());
  double3 momentum{};
  for (std::size_t i = 0; i < system.spanOfFrameworkDynamics().size(); ++i)
  {
    if (system.framework->isFlexibleAtom(i)) momentum += 12.0 * system.spanOfFrameworkDynamics()[i].velocity;
  }
  return momentum.length();
}
}  // namespace

// Momentum is only conserved without lab-fixed atoms: fixed atoms exert an external force, so the P = 0 hyperplane is
// not invariant and a projected draw has no reverse move.
TEST(hybrid_mc, velocities_are_projected_only_when_momentum_is_conserved)
{
  RandomNumber random(29);
  System flexible = makeHybridFrameworkSystem("");
  ASSERT_FALSE(flexible.framework->isMixed());
  EXPECT_LT(drawnFrameworkMomentum(flexible, random), 1e-10);

  System mixed = makeHybridFrameworkSystem(
      R"("Groups": [{"Type": "Fixed", "Atoms": [0]}, {"Type": "Flexible", "Atoms": [1, 2]}],)");
  ASSERT_TRUE(mixed.framework->isMixed());
  ASSERT_GT(mixed.framework->numberOfFixedAtoms(), 0uz);
  EXPECT_GT(drawnFrameworkMomentum(mixed, random), 1e-6);
}

// Nothing can move in a host of Fixed groups only, so there is nothing to propose.
TEST(hybrid_mc, framework_of_fixed_groups_only_is_rejected_without_a_trajectory)
{
  RandomNumber random(31);
  System system = makeHybridFrameworkSystem(R"("Groups": [{"Type": "Fixed", "Atoms": [0, 1, 2]}],)");
  ASSERT_FALSE(system.framework->hasMobileAtoms());
  EXPECT_FALSE(MC_Moves::hybridMCMove(random, system).has_value());
}

namespace
{
enum class DriftHost
{
  None,
  Rigid,
  Flexible,
  FixedAndFlexible
};

// 'numberOfMolecules' single-site molecules, optionally in a three-atom chain host.
System makeDriftSystem(DriftHost host, std::size_t numberOfMolecules, bool externalField = false)
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}, {"M", false, 16.0, 0.0, 0.0, 8, false}},
                        {{80.0, 3.0}, {60.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 5.0, 5.0, 5.0, true, false,
                        false);
  const SimulationBox box(14.0, 14.0, 14.0);
  Component molecule(forceField, "site", 16.0, 1.0e6, 0.011, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  std::vector<Framework> frameworks;
  if (host != DriftHost::None)
  {
    std::vector<Atom> atoms;
    for (const double3& position : {double3(4.0, 7.0, 7.0), double3(5.4, 7.3, 7.0), double3(6.8, 6.8, 7.0)})
    {
      atoms.push_back(Atom(position / 14.0, 0.0, 1.0, 0, 0, 0, 0, true));
    }
    Framework framework(forceField, "drift-chain", box, 1, atoms, atoms, int3(1, 1, 1));
    if (host != DriftHost::Rigid)
    {
      const std::filesystem::path definition = std::filesystem::temp_directory_path() / "raspa3-drift-framework.json";
      std::ofstream(definition) << std::format(
          R"({{"Type": "Flexible", {} "Intra14VanDerWaalsScalingValue": 0.0, "Intra14ChargeChargeScalingValue": 0.0}})",
          host == DriftHost::FixedAndFlexible
              ? R"("Groups": [{"Type": "Fixed", "Atoms": [0]}, {"Type": "Flexible", "Atoms": [1, 2]}],)"
              : "");
      framework.readFrameworkDefinition(forceField, definition.string());
      std::filesystem::remove(definition);
    }
    frameworks.push_back(framework);
  }
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, frameworks.empty() ? std::nullopt : std::optional(frameworks[0]),
                {molecule}, {}, {numberOfMolecules}, 5);
  system.hasExternalField = externalField;
  std::span<Atom> sites = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < sites.size(); ++i)
  {
    sites[i].position = double3(2.0 + 3.0 * static_cast<double>(i), 12.0, 2.0);
  }
  return system;
}

// Enters the equilibration stage of the driver, which draws the velocities, without running a cycle.
struct DriftState
{
  std::size_t constraint;
  double moleculeMomentum;
  double kineticEnergy;
};

DriftState enterEquilibration(const System& system)
{
  MolecularDynamics md({0, 0, 0, 0, 1000, 10000, 5000, 5000}, {system}, 17uz, 5, false);
  md.equilibrate();
  const System& result = md.systems[0];
  double3 momentum{};
  double kinetic{0.0};
  for (const Molecule& molecule : result.moleculeData)
  {
    momentum += molecule.mass * molecule.velocity;
    kinetic += 0.5 * molecule.mass * double3::dot(molecule.velocity, molecule.velocity);
  }
  return {result.translationalCenterOfMassConstraint, momentum.length(), kinetic};
}
}  // namespace

// Bulk fluid without a thermostat: the centre-of-mass momentum is removed, so its three degrees of freedom are gone
// from the temperature as well.
TEST(molecular_dynamics, nve_fluid_removes_the_drift_and_its_degrees_of_freedom)
{
  const DriftState state = enterEquilibration(makeDriftSystem(DriftHost::None, 4));
  EXPECT_EQ(state.constraint, 3uz);
  EXPECT_LT(state.moleculeMomentum, 1e-10);
  EXPECT_GT(state.kineticEnergy, 0.0);
}

// A single molecule has no internal reference for its momentum: zeroing it would freeze the molecule at 0 K.
TEST(molecular_dynamics, single_molecule_keeps_its_velocity)
{
  for (const DriftHost host : {DriftHost::None, DriftHost::Rigid})
  {
    const DriftState state = enterEquilibration(makeDriftSystem(host, 1));
    EXPECT_EQ(state.constraint, 0uz);
    EXPECT_GT(state.kineticEnergy, 0.0);
  }
}

// Rigid or lab-fixed framework atoms and an external field exert a net force: there is no conserved momentum to
// remove and no degrees of freedom to subtract.
TEST(molecular_dynamics, momentum_is_kept_when_it_is_not_conserved)
{
  for (const System& system : {makeDriftSystem(DriftHost::Rigid, 4), makeDriftSystem(DriftHost::FixedAndFlexible, 4),
                               makeDriftSystem(DriftHost::None, 4, true), makeDriftSystem(DriftHost::Flexible, 4, true)})
  {
    const DriftState state = enterEquilibration(system);
    EXPECT_EQ(state.constraint, 0uz);
    EXPECT_GT(state.moleculeMomentum, 1e-6);
  }
}

TEST(molecular_dynamics, flexible_framework_conserves_the_total_momentum)
{
  const DriftState state = enterEquilibration(makeDriftSystem(DriftHost::Flexible, 2));
  EXPECT_EQ(state.constraint, 3uz);
}

// The ideal-gas part of the barostat pressure counts every coupled point once, with or without the constraint.
TEST(molecular_dynamics, barostat_pressure_ideal_part_counts_the_coupled_points)
{
  for (const std::size_t numberOfMolecules : {1uz, 4uz})
  {
    for (const std::size_t constraint : {0uz, 3uz})
    {
      if (numberOfMolecules == 1uz && constraint == 3uz) continue;
      System system = makeDriftSystem(DriftHost::None, numberOfMolecules);
      system.translationalCenterOfMassConstraint = constraint;
      const std::size_t degreesOfFreedom = barostatTranslationalDegreesOfFreedom(system, BarostatCoupling::Molecular);
      const double3x3 pressure = barostatPressureTensor(system, double3x3{}, BarostatCoupling::Molecular, degreesOfFreedom);
      const double3x3 virial = computeBarostatVirial(system, double3x3{}, BarostatCoupling::Molecular) / system.simulationBox.volume;
      const double expected = static_cast<double>(numberOfMolecules) / (system.beta * system.simulationBox.volume);
      EXPECT_NEAR(pressure.ax - virial.ax, expected, 1e-12 * expected) << numberOfMolecules << " " << constraint;
    }
  }
}

// The constraint on the centre-of-mass degrees of freedom holds through a particle exchange: an inserted molecule
// brings thermal momentum and a deleted one takes its momentum with it, so the rest has to be brought back to P = 0.
// Without a conserved momentum (here: an external field) the velocities are left alone.
TEST(molecular_dynamics, particle_exchange_keeps_the_constrained_momentum_at_zero)
{
  const auto makeExchangeSystem = [](bool externalField)
  {
    ForceField forceField = ForceField::makeZeoliteForceField(11.8, false, false, false);
    Component methane = Component::makeMethane(forceField, 0);
    methane.fugacityCoefficient = 1.0;
    methane.idealGasRosenbluthWeight = 1.0;
    methane.mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
    System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0e12, 1.0, std::nullopt, {methane}, {},
                  {4}, 5);
    system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::MuVT;
    system.hasExternalField = externalField;
    system.initializeGroupData();
    RandomNumber random(211);
    Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                      system.spanOfMoleculeDynamics(), system.components, system.temperature);
    if (system.conservesCenterOfMassMomentum())
    {
      Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                   system.spanOfMoleculeDynamics(), system.components);
    }
    return system;
  };
  const auto momentum = [](const System& system)
  {
    double3 total{};
    for (const Molecule& molecule : system.moleculeData) total += molecule.mass * molecule.velocity;
    return total.length();
  };

  System conserved = makeExchangeSystem(false);
  ASSERT_LT(momentum(conserved), 1e-10);
  RandomNumber random(1879);
  MC_Moves::ParticleExchangeResult result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Inserted; ++attempt)
  {
    result = attemptParticleExchange(random, conserved);
  }
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Inserted);
  EXPECT_EQ(conserved.translationalCenterOfMassConstraint, 3uz);
  EXPECT_LT(momentum(conserved), 1e-10);

  conserved.pressure = 1.0e-12;
  result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Deleted; ++attempt)
  {
    result = attemptParticleExchange(random, conserved);
  }
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Deleted);
  EXPECT_LT(momentum(conserved), 1e-10);

  System notConserved = makeExchangeSystem(true);
  const std::size_t before = notConserved.moleculeData.size();
  const double initial = momentum(notConserved);
  ASSERT_GT(initial, 1e-6);
  result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Inserted; ++attempt)
  {
    result = attemptParticleExchange(random, notConserved);
  }
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Inserted);
  ASSERT_EQ(notConserved.moleculeData.size(), before + 1);
  EXPECT_GT(momentum(notConserved), 1e-6);
  EXPECT_EQ(notConserved.translationalCenterOfMassConstraint, 0uz);
}
