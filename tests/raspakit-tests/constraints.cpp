#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

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
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import system;
import simulationbox;
import running_energy;
import integrators_update;
import integrators;
import integrators_compute;
import molecule;
import thermostat;
import thermobarostat;
import minimization_cell_layout;
import molecular_dynamics;
import randomnumbers;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;
import spatial_decomposition_integrator;
import force_engine;

// Molecular dynamics holds the FIXED bonds of a flexible molecule with SHAKE/RATTLE: the constraints follow from
// the component, take their degrees of freedom off, and are kept by every integration path.

namespace
{
constexpr double bondLength = 1.54;

ForceField makeForceField()
{
  return ForceField({{"CH3", false, 15.03452, 0.0, 0.0, 6, false}}, {{108.0, 3.76}},
                    ForceField::MixingRule::Lorentz_Berthelot, 6.5, 6.5, 6.5, true, false, false);
}

ConnectivityTable chain(std::size_t numberOfAtoms)
{
  ConnectivityTable table(numberOfAtoms);
  for (std::size_t i = 0; i + 1 < numberOfAtoms; ++i)
  {
    table[i, i + 1] = true;
    table[i + 1, i] = true;
  }
  return table;
}

Atom site(double3 position) { return Atom(position, 0.0, 1.0, 0, 0, 0, false, false); }

// A dimer: with a FIXED bond (flexible atoms, one constraint) or as a rigid body (no bond at all).
Component makeDimer(const ForceField& forceField, bool constrained)
{
  Potentials::IntraMolecularPotentials potentials{};
  if (constrained) potentials.bonds = {BondPotential({0, 1}, BondType::Fixed, {bondLength})};
  return Component(forceField, "dimer", 305.33, 4871800.0, 0.0993,
                   {site({-0.5 * bondLength, 0.0, 0.0}), site({0.5 * bondLength, 0.0, 0.0})}, chain(2), potentials, 5,
                   21);
}

// A bent trimer held by two FIXED bonds and the FIXED angle between them.
Component makeTrimer(const ForceField& forceField, BondType armType = BondType::Fixed)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bonds = {BondPotential({0, 1}, BondType::Fixed, {1.0}),
                      armType == BondType::Fixed ? BondPotential({1, 2}, BondType::Fixed, {1.0})
                                                 : BondPotential({1, 2}, BondType::Harmonic, {50000.0, 1.0})};
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Fixed, {109.47})};
  const double half = 0.5 * 109.47 * std::numbers::pi / 180.0;
  return Component(forceField, "trimer", 305.33, 4871800.0, 0.0993,
                   {site({std::sin(half), std::cos(half), 0.0}), site({0.0, 0.0, 0.0}),
                    site({-std::sin(half), std::cos(half), 0.0})},
                   chain(3), potentials, 5, 21);
}

// The molecules on a lattice of spacing 'spacing', each with the geometry of the component turned to a random
// orientation: the constraints are satisfied, the molecules are within range of one another.
std::vector<double3> latticePositions(const Component& component, std::size_t perAxis, double spacing,
                                      RandomNumber& random)
{
  std::vector<double3> positions;
  for (std::size_t ix = 0; ix != perAxis; ++ix)
  {
    for (std::size_t iy = 0; iy != perAxis; ++iy)
    {
      for (std::size_t iz = 0; iz != perAxis; ++iz)
      {
        const double3 center = spacing * double3(static_cast<double>(ix) + 0.5, static_cast<double>(iy) + 0.5,
                                                 static_cast<double>(iz) + 0.5);
        const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
        for (const Atom& atom : component.atoms) positions.push_back(center + rotation * atom.position);
      }
    }
  }
  return positions;
}

System makeSystem(const ForceField& forceField, const Component& component, std::size_t perAxis, double spacing,
                  std::size_t seed)
{
  RandomNumber random(seed);
  const double length = spacing * static_cast<double>(perAxis);
  return System(forceField, SimulationBox(length, length, length), false, 300.0, 1e5, 1.0, {}, {component},
                {latticePositions(component, perAxis, spacing, random)}, {0}, 5);
}

RunningEnergy step(System& system, std::optional<Thermostat>& thermostat)
{
  return Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(),
                                     system.spanOfMoleculeDynamics(), system.components, system.timeStep, thermostat,
                                     system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
                                     system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
                                     system.fixedFrameworkStoredEik, system.interpolationGrids,
                                     system.numberOfMoleculesPerComponent);
}

void drawVelocities(System& system, std::size_t seed)
{
  RandomNumber random(seed);
  Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                    system.spanOfMoleculeDynamics(), system.components, system.temperature);
}

// The largest relative error of a constrained length, and the largest velocity along a constraint relative to the
// largest atomic velocity.
struct Residuals
{
  double length{};
  double velocity{};
};

Residuals residuals(const System& system)
{
  Residuals result{};
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  double scale{};
  for (const AtomDynamics& d : dynamics) scale = std::max(scale, d.velocity.length());
  for (const Molecule& molecule : system.moleculeData)
  {
    for (const Component::DistanceConstraint& c : system.components[molecule.componentId].distanceConstraints)
    {
      const double3 separation = atoms[molecule.atomIndex + c.a].position - atoms[molecule.atomIndex + c.b].position;
      const double3 relative = dynamics[molecule.atomIndex + c.a].velocity - dynamics[molecule.atomIndex + c.b].velocity;
      result.length = std::max(result.length, std::abs(separation.length() - c.length) / c.length);
      if (scale > 0.0)
      {
        result.velocity =
            std::max(result.velocity, std::abs(double3::dot(separation, relative)) / (separation.length() * scale));
      }
    }
  }
  return result;
}
}  // namespace

TEST(constraints, component_lists_fixed_bonds_and_takes_their_degrees_of_freedom_off)
{
  const ForceField forceField = makeForceField();

  const Component dimer = makeDimer(forceField, true);
  ASSERT_FALSE(dimer.rigid);
  ASSERT_EQ(dimer.distanceConstraints.size(), 1uz);
  EXPECT_DOUBLE_EQ(dimer.distanceConstraints[0].length, bondLength);
  EXPECT_TRUE(dimer.unsupportedDynamicsConstraints.empty());
  EXPECT_EQ(dimer.translationalDegreesOfFreedom, 5uz);
  EXPECT_EQ(dimer.rotationalDegreesOfFreedom, 0uz);

  // the rigid dimer has the same five: three translations and the two rotations of a linear body
  const Component rigid = makeDimer(forceField, false);
  ASSERT_TRUE(rigid.rigid);
  EXPECT_TRUE(rigid.distanceConstraints.empty());
  EXPECT_EQ(rigid.translationalDegreesOfFreedom + rigid.rotationalDegreesOfFreedom, 5uz);

  // two fixed bonds and the fixed angle between them: the 1-3 distance closes the triangle
  const Component trimer = makeTrimer(forceField);
  ASSERT_EQ(trimer.distanceConstraints.size(), 3uz);
  EXPECT_NEAR(trimer.distanceConstraints[2].length, 2.0 * std::sin(0.5 * 109.47 * std::numbers::pi / 180.0), 1e-12);
  EXPECT_TRUE(trimer.unsupportedDynamicsConstraints.empty());
  EXPECT_EQ(trimer.translationalDegreesOfFreedom, 6uz);
  EXPECT_NO_THROW(trimer.requireSupportedDynamicsConstraints());

  // the system sums the components
  System system = makeSystem(forceField, dimer, 3, 5.0, 3);
  EXPECT_EQ(system.translationalDegreesOfFreedom, 27uz * 5uz);
}

TEST(constraints, a_fixed_angle_with_a_flexible_arm_is_refused_by_the_dynamics)
{
  const ForceField forceField = makeForceField();
  const Component trimer = makeTrimer(forceField, BondType::Harmonic);
  EXPECT_EQ(trimer.distanceConstraints.size(), 1uz);
  ASSERT_EQ(trimer.unsupportedDynamicsConstraints.size(), 1uz);
  EXPECT_THROW(trimer.requireSupportedDynamicsConstraints(), std::runtime_error);

  System system = makeSystem(forceField, trimer, 3, 6.0, 5);
  system.precomputeTotalGradients();
  std::optional<Thermostat> thermostat{};
  EXPECT_THROW(step(system, thermostat), std::runtime_error);
}

TEST(constraints, drawn_velocities_are_stationary_and_carry_the_temperature_of_the_remaining_freedom)
{
  const ForceField forceField = makeForceField();
  for (const Component& component : {makeDimer(forceField, true), makeTrimer(forceField)})
  {
    System system = makeSystem(forceField, component, 3, 6.0, 11);
    RandomNumber random(29);
    double kinetic{};
    const std::size_t draws = 400;
    for (std::size_t draw = 0; draw != draws; ++draw)
    {
      Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                        system.spanOfMoleculeDynamics(), system.components, system.temperature);
      EXPECT_LT(residuals(system).velocity, 1e-7);
      kinetic += Integrators::computeTranslationalKineticEnergy(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                                system.spanOfMoleculeDynamics(), system.components);
    }
    // equipartition over the degrees of freedom the constraints leave (27 molecules, 400 draws: 1% is 3 sigma)
    const double expected =
        0.5 * static_cast<double>(system.translationalDegreesOfFreedom) * Units::KB * system.temperature;
    EXPECT_NEAR(kinetic / static_cast<double>(draws), expected, 0.012 * expected) << component.name;
  }
}

TEST(constraints, velocity_verlet_keeps_the_lengths_and_conserves_the_energy)
{
  const ForceField forceField = makeForceField();
  for (const Component& component : {makeDimer(forceField, true), makeTrimer(forceField)})
  {
    System system = makeSystem(forceField, component, 3, 5.5, 17);
    system.timeStep = 0.001;
    drawVelocities(system, 41);
    system.precomputeTotalGradients();
    std::optional<Thermostat> thermostat{};

    const double first = step(system, thermostat).conservedEnergy();
    double kinetic{}, last{first}, largest{};
    for (std::size_t i = 0; i != 500; ++i)
    {
      const RunningEnergy energy = step(system, thermostat);
      last = energy.conservedEnergy();
      largest = std::max(largest, std::abs(last - first));
      kinetic = energy.translationalKineticEnergy;
      const Residuals r = residuals(system);
      ASSERT_LT(r.length, 1e-7) << component.name << " step " << i;
      ASSERT_LT(r.velocity, 1e-7) << component.name << " step " << i;
    }
    // no drift: the energy stays within a small fraction of the kinetic energy over 0.5 ps
    EXPECT_LT(largest, 2e-4 * kinetic) << component.name;
  }
}

// The constrained dimer and the rigid dimer are the same mechanical system, integrated by two unrelated schemes
// (RATTLE on atoms, the quaternion free rotor on the body): from the same start they follow the same trajectory,
// and the difference between the two falls with the square of the time step.
TEST(constraints, constrained_dimer_follows_the_rigid_dimer)
{
  const ForceField forceField = makeForceField();
  const auto difference = [&](double timeStep, std::size_t numberOfSteps)
  {
    System constrained = makeSystem(forceField, makeDimer(forceField, true), 3, 4.6, 23);
    System rigid = makeSystem(forceField, makeDimer(forceField, false), 3, 4.6, 23);
    std::optional<Thermostat> thermostat{};
    double result{};
    for (System* system : {&constrained, &rigid})
    {
      system->timeStep = timeStep;
      system->precomputeTotalGradients();
    }
    EXPECT_NEAR(constrained.computeTotalEnergies().potentialEnergy(), rigid.computeTotalEnergies().potentialEnergy(),
                1e-9 * std::abs(rigid.computeTotalEnergies().potentialEnergy()));
    for (std::size_t i = 0; i != numberOfSteps; ++i)
    {
      step(constrained, thermostat);
      step(rigid, thermostat);
    }
    const std::span<const Atom> a = constrained.spanOfMoleculeAtoms();
    const std::span<const Atom> b = rigid.spanOfMoleculeAtoms();
    double moved{};
    for (std::size_t i = 0; i != a.size(); ++i) result = std::max(result, (a[i].position - b[i].position).length());
    // the start is at rest: the molecules have to have moved for the comparison to mean anything
    System start = makeSystem(forceField, makeDimer(forceField, true), 3, 4.6, 23);
    for (std::size_t i = 0; i != a.size(); ++i)
      moved = std::max(moved, (a[i].position - start.spanOfMoleculeAtoms()[i].position).length());
    EXPECT_GT(moved, 0.05);
    return result;
  };

  const double coarse = difference(0.002, 200);
  const double fine = difference(0.001, 400);
  EXPECT_LT(fine, 1e-4);
  EXPECT_LT(fine, 0.3 * coarse);
}

TEST(constraints, nose_hoover_thermostat_holds_the_temperature_of_the_constrained_fluid)
{
  const ForceField forceField = makeForceField();
  System system = makeSystem(forceField, makeDimer(forceField, true), 4, 5.5, 31);
  system.timeStep = 0.002;
  drawVelocities(system, 43);
  system.precomputeTotalGradients();
  std::optional<Thermostat> thermostat = Thermostat(system.temperature, system.timeStep,
                                                    system.translationalDegreesOfFreedom, 0, 3, 1, 0.15);
  RandomNumber random(47);
  thermostat->initialize(random);

  double kinetic{};
  std::size_t samples{};
  for (std::size_t i = 0; i != 6000; ++i)
  {
    const RunningEnergy energy = step(system, thermostat);
    if (i < 1000) continue;
    kinetic += energy.translationalKineticEnergy;
    ++samples;
  }
  EXPECT_LT(residuals(system).length, 1e-7);
  const double temperature = 2.0 * kinetic / (static_cast<double>(samples) * Units::KB *
                                              static_cast<double>(system.translationalDegreesOfFreedom));
  EXPECT_NEAR(temperature, system.temperature, 0.03 * system.temperature);
}

TEST(constraints, barostat_step_keeps_the_lengths_and_needs_molecular_coupling)
{
  const ForceField forceField = makeForceField();
  System system = makeSystem(forceField, makeTrimer(forceField), 3, 6.0, 37);
  system.timeStep = 0.001;
  drawVelocities(system, 53);
  system.precomputeTotalGradients();
  system.thermostat =
      Thermostat(system.temperature, system.timeStep, system.translationalDegreesOfFreedom, 0, 3, 1, 0.15);
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta,
                     system.temperature, system.pressure, system.timeStep, 3 * system.moleculeData.size(), 3, 1, 1.0);
  RandomNumber random(59);
  system.thermostat->initialize(random);
  system.thermobarostat->initialize(random);

  const double volume = system.simulationBox.volume;
  for (std::size_t i = 0; i != 200; ++i)
  {
    molecularDynamicsStep(system);
    const Residuals r = residuals(system);
    ASSERT_LT(r.length, 1e-7) << "step " << i;
    ASSERT_LT(r.velocity, 1e-7) << "step " << i;
  }
  EXPECT_NE(system.simulationBox.volume, volume);

  system.thermobarostat->coupling = BarostatCoupling::Atomic;
  EXPECT_THROW(molecularDynamicsStep(system), std::runtime_error);
}

TEST(constraints, domain_decomposed_step_follows_the_serial_step)
{
  const ForceField forceField = makeForceField();
  System serial = makeSystem(forceField, makeTrimer(forceField), 5, 6.0, 61);
  System decomposed = serial;
  for (System* system : {&serial, &decomposed})
  {
    system->timeStep = 0.001;
    drawVelocities(*system, 67);
  }
  serial.precomputeTotalGradients();

  SpatialDecompositionSettings settings;
  settings.numberOfThreads = 4;
  settings.verletSkin = 1.0;
  ForceEngine engine(settings);
  engine.initialize(decomposed);
  engine.computeGradients(decomposed, true);

  std::optional<Thermostat> thermostat{};
  for (std::size_t i = 0; i != 50; ++i)
  {
    step(serial, thermostat);
    SpatialDecompositionIntegrator::engineVelocityVerlet(decomposed, engine);
  }
  EXPECT_LT(residuals(decomposed).length, 1e-7);
  EXPECT_LT(residuals(decomposed).velocity, 1e-7);
  const std::span<const Atom> a = serial.spanOfMoleculeAtoms();
  const std::span<const Atom> b = decomposed.spanOfMoleculeAtoms();
  double difference{};
  for (std::size_t i = 0; i != a.size(); ++i)
    difference = std::max(difference, (a[i].position - b[i].position).length());
  EXPECT_LT(difference, 1e-9);
}
