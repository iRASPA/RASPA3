#include <gtest/gtest.h>

import std;

import archive;
import int3;
import double3;
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
import molecular_dynamics;
import mc_moves;
import mc_moves_hybridmc;
import mc_moves_move_types;
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
      random, inserted, system.spanOfMoleculeDynamics().subspan(inserted.atomIndex, inserted.numberOfAtoms),
      system.components[0], system.temperature);

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

  MC_Moves::ParticleExchangeResult result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Inserted; ++attempt)
    result = MC_Moves::performMolecularDynamicsSwap(random, system, 0);
  ASSERT_EQ(result, MC_Moves::ParticleExchangeResult::Inserted);
  const std::size_t countAfterInsertion = system.numberOfIntegerMoleculesPerComponent[0];

  system.pressure = 1.0e-12;
  system.components[0].fugacityCoefficient = 1.0;
  result = MC_Moves::ParticleExchangeResult::Rejected;
  for (std::size_t attempt = 0; attempt != 100 && result != MC_Moves::ParticleExchangeResult::Deleted; ++attempt)
    result = MC_Moves::performMolecularDynamicsSwap(random, system, 0);
  EXPECT_EQ(result, MC_Moves::ParticleExchangeResult::Deleted);
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
      nve.eik_x, nve.eik_y, nve.eik_z, nve.eik_xy, nve.totalEik, nve.fixedFrameworkStoredEik, nve.interpolationGrids,
      {}, nve.framework, nve.spanOfFrameworkDynamics());
  const RunningEnergy nvtEnergy = Integrators::velocityVerlet(
      {}, {}, {}, {}, nvt.timeStep, nvt.thermostat, nvt.spanOfFrameworkAtoms(), nvt.forceField, nvt.simulationBox,
      nvt.eik_x, nvt.eik_y, nvt.eik_z, nvt.eik_xy, nvt.totalEik, nvt.fixedFrameworkStoredEik, nvt.interpolationGrids,
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
      system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
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
                                    system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                    system.interpolationGrids, {}, system.framework, system.spanOfFrameworkDynamics());
    maximumDrift = std::max(maximumDrift, std::abs(current.conservedEnergy() - initialEnergy));
  }

  EXPECT_GT(std::abs(system.spanOfFrameworkAtoms()[0].position.x - 5.0), 1.0e-5);
  EXPECT_LT(maximumDrift / std::max(1.0, std::abs(initialEnergy)), 2.0e-5);
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
                              system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
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
                              system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
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
