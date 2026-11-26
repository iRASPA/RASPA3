#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import double3;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import fragment;
import simulationbox;
import system;
import thermobarostat;
import minimization_cell_layout;
import simd_quatd;
import molecule;
import mc_moves_probabilities;
import mc_moves_move_types;
import minimization;
import minimization_options;
import monte_carlo;
import simulation_schedule;
import molecular_dynamics;
import thermostat;
import property_elastic_constants_fluctuation;
import write_lammps_data;

// A light rigid body whose moments lie between the zero thresholds of the two counts: the inverse inertia zeroes an
// axis below max(1e-2, S) * 1e-5, and the degrees of freedom were counted against S * 1e-5, or against 1e-5 when the
// sum S of the moments is below 1e-5. A diatomic with I = 1e-6 about its two transverse axes rotated about both
// (their inverse inertia is not zero) while it counted no rotational degree of freedom at all.
TEST(SIMULATION_SETUP_FIXES, rotational_degrees_of_freedom_count_the_axes_the_integrator_rotates_about)
{
  const double a = std::sqrt(5.0e-7);
  const std::vector<double3> positions{double3(-a, 0.0, 0.0), double3(a, 0.0, 0.0)};
  const std::vector<double> masses{1.0, 1.0};

  Fragment fragment(std::vector<std::size_t>{0, 1});
  fragment.computeRigidProperties(positions, masses);
  const std::size_t rotatedAxes = (fragment.inverseInertiaVector.x != 0.0) +
                                  (fragment.inverseInertiaVector.y != 0.0) + (fragment.inverseInertiaVector.z != 0.0);
  EXPECT_EQ(rotatedAxes, 2uz);
  EXPECT_EQ(fragment.rotationalDegreesOfFreedom, rotatedAxes);

  const ForceField forceField({{"H", false, 1.0, 0.0, 0.0, 1, false}}, {VDWParameters(10.0, 2.0)},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  const Component diatomic(forceField, "diatomic", 33.0, 1.3e6, 0.0,
                           {Atom(positions[0], 0.0, 1.0, 0, 0, 0, false, false),
                            Atom(positions[1], 0.0, 1.0, 0, 0, 0, false, false)},
                           {}, {}, 5, 21);
  const std::size_t componentAxes = (diatomic.inverseInertiaVector.x != 0.0) +
                                    (diatomic.inverseInertiaVector.y != 0.0) +
                                    (diatomic.inverseInertiaVector.z != 0.0);
  EXPECT_EQ(componentAxes, 2uz);
  EXPECT_EQ(diatomic.rotationalDegreesOfFreedom, componentAxes);
}

// An automatic cut-off follows the cell: under a barostat the energy then depends on the volume through the cut-off,
// which the virial leaves out. Refused, except the automatic Coulomb cut-off of Ewald with automatic parameters (whose
// energy does not depend on the cut-off within the precision).
TEST(SIMULATION_SETUP_FIXES, a_barostat_refuses_cut_offs_that_follow_the_cell)
{
  ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}}, {VDWParameters(119.8, 3.4)},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
  const Component argon(forceField, "argon", 150.86, 4.898e6, -0.00219,
                        {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);

  const auto reason = [&](const ForceField& field)
  {
    System system(field, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {argon}, {}, {0}, 5);
    system.thermobarostat =
        Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta,
                       300.0, system.pressure, 0.0005, system.translationalDegreesOfFreedom);
    return system.thermobarostatUnsupportedReason();
  };

  EXPECT_FALSE(reason(forceField).has_value());

  ForceField automaticVDW = forceField;
  automaticVDW.cutOffMoleculeVDWAutomatic = true;
  EXPECT_TRUE(reason(automaticVDW).has_value());

  ForceField automaticEwald = forceField;
  automaticEwald.cutOffCoulombAutomatic = true;
  automaticEwald.automaticEwald = true;
  EXPECT_FALSE(reason(automaticEwald).has_value());

  ForceField automaticWolf = automaticEwald;
  automaticWolf.chargeMethod = ForceField::ChargeMethod::Wolf;
  EXPECT_TRUE(reason(automaticWolf).has_value());
}

// United-atom butane with a rigid, linear CH2-CH2 middle: the orientation of a linear group is fixed by its atoms up
// to a rotation about its axis.
constexpr std::string_view kButaneWithLinearMiddleJson = R"({
  "CriticalTemperature" : 425.12,
  "CriticalPressure" : 3796000.0,
  "AcentricFactor" : 0.2,
  "pseudoAtoms" : [["CH3", [0.0, 0.0, 0.0]], ["CH2", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH3", [0.0, 0.0, 0.0]]],
  "Connectivity" : [[0, 1], [1, 2], [2, 3]],
  "RigidBodies" : [[1, 2]],
  "Bonds" : [[["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]], [["CH2", "CH2"], "FIXED", [1.54]]],
  "Bends" : [[["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114]]],
  "Torsions" : [[["CH3", "CH2", "CH2", "CH3"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]],
  "VanDerWaals" : "auto"
}
)";

// The group states are fitted to the atoms again when production starts; the stored orientation of a linear group
// can differ from the fit by a rotation about its axis, and its momentum has to follow the fitted orientation, or
// the angular velocity in the laboratory frame turns about the axis of the group by that angle.
TEST(SIMULATION_SETUP_FIXES, refitted_linear_group_keeps_its_angular_momentum)
{
  const ForceField forceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                              {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                              12.0, true, false, false);
  TemporaryFile file("butane-linear-middle.json", kButaneWithLinearMiddleJson);
  const Component butane(Component::Type::Adsorbate, 0, forceField, "butane", file.stemPath().string(), 5, 21,
                         MCMoveProbabilities(), std::nullopt, false);
  ASSERT_TRUE(butane.isSemiFlexible());
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {butane}, {}, {1}, 5);
  system.initializeGroupData();
  ASSERT_EQ(system.groupData.size(), 1uz);

  const auto rotate = [](const simd_quatd& q, const double3& v)
  {
    const simd_quatd r = q * simd_quatd(0.0, v) * simd_quatd(-q.ix, -q.iy, -q.iz, q.r);
    return double3(r.ix, r.iy, r.iz);
  };
  // the angular momentum in the laboratory frame of a group state, from L_body = Im(p~ q) / 2 with p~ = (-p_r, p)
  const auto laboratoryAngularMomentum = [&](const GroupState& state)
  {
    const simd_quatd& p = state.orientationMomentum;
    const simd_quatd pq = simd_quatd(p.ix, p.iy, p.iz, -p.r) * state.orientation;
    return rotate(state.orientation, 0.5 * double3(pq.ix, pq.iy, pq.iz));
  };

  // the axis of the group is its body axis without moment of inertia
  const Fragment& group = *std::ranges::find_if(system.components[0].fragmentGraph.fragments,
                                                [](const Fragment& f) { return f.isRigidBody(); });
  const double3 inertia = group.inertiaVector;
  const std::size_t axis = inertia.x < inertia.y ? (inertia.x < inertia.z ? 0 : 2) : (inertia.y < inertia.z ? 1 : 2);
  double3 axialDirection{};
  axialDirection[axis] = 1.0;
  double3 transverse{};
  transverse[(axis + 1) % 3] = 1.0;

  // a stored orientation that differs from the fit by a rotation of one radian about the axis, spinning about a
  // transverse body axis
  GroupState& state = system.groupData[0];
  const simd_quatd fitted = state.orientation;
  const double angle = 1.0;
  const simd_quatd gauge(std::cos(0.5 * angle), -std::sin(0.5 * angle) * axialDirection);
  state.orientation = fitted * gauge;
  const double3 bodyAngularMomentum = 25.0 * transverse;
  state.orientationMomentum = 2.0 * (state.orientation * simd_quatd(0.0, bodyAngularMomentum));
  const double3 before = laboratoryAngularMomentum(state);

  system.initializeGroupData();
  const double3 after = laboratoryAngularMomentum(system.groupData[0]);
  EXPECT_NEAR((after - before).length(), 0.0, 1.0e-9 * before.length());
  EXPECT_NEAR(system.groupData[0].orientation.r * fitted.r + system.groupData[0].orientation.ix * fitted.ix +
                  system.groupData[0].orientation.iy * fitted.iy + system.groupData[0].orientation.iz * fitted.iz,
              1.0, 1.0e-12);
}

// A fractional molecule pinned at a lambda is created switched off, where it feels nothing and may sit on top of
// another molecule; at the pinned lambda it is a full molecule there (above lambda = 0.5 without soft core). It was
// grown anew only from inside a blocked pocket, and an overlap stayed.
TEST(SIMULATION_SETUP_FIXES, pinned_fractional_molecule_is_regrown_out_of_an_overlap)
{
  const ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}}, {VDWParameters(119.8, 3.4)},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  Component argon(forceField, "argon", 150.86, 4.898e6, -0.00219,
                  {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities);
  argon.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::GC;
  argon.fixedLambdaBin = 18;
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {argon}, {}, {20}, 5);

  const std::size_t fractional = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(0);
  const std::size_t other = fractional == 0 ? 1 : 0;
  std::span<Atom> molecule = system.spanOfMolecule(0, fractional);
  ASSERT_GT(molecule[0].scalingVDW, 0.5);
  const double3 occupied = system.spanOfMolecule(0, other)[0].position;
  molecule[0].position = occupied + double3(0.1, 0.0, 0.0);

  system.initializeFixedLambdaFractionalMolecules();
  const double3 separation = system.simulationBox.applyPeriodicBoundaryConditions(
      system.spanOfMolecule(0, fractional)[0].position - occupied);
  EXPECT_GT(separation.length(), 2.0);
  EXPECT_GT(system.spanOfMolecule(0, fractional)[0].scalingVDW, 0.5);
}

// A FIXED bond between flexible atoms has no energy, gradient or Hessian: the dynamics holds it with SHAKE/RATTLE.
// The minimizer moved its atoms freely, pulled the bond apart along the other forces, and the normal modes counted it
// as a mode of zero frequency; refused.
TEST(SIMULATION_SETUP_FIXES, minimization_refuses_fixed_bonds_between_flexible_atoms)
{
  const ForceField forceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                              {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                              12.0, true, false, false);
  TemporaryFile file("propane-fixed-bonds.json", molecule_fixtures::kPropaneJson);
  const Component propane(Component::Type::Adsorbate, 0, forceField, "propane", file.stemPath().string(), 5, 21,
                          MCMoveProbabilities(), std::nullopt, false);
  ASSERT_FALSE(propane.distanceConstraints.empty());
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {propane}, {}, {2}, 5);

  Minimization minimization(MinimizationOptions{}, {system});
  EXPECT_THROW(minimization.setup(), std::runtime_error);
}

// The stress-fluctuation property keeps blocks as the averages of the system do: one made with another number of
// blocks was taken as it was, and the sampling of a block past its last wrote past its end; the system now reports
// it as not sized and sizes it with the rest.
TEST(SIMULATION_SETUP_FIXES, elastic_constants_property_follows_the_blocks_of_the_system)
{
  const ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}}, {VDWParameters(119.8, 3.4)},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component argon(forceField, "argon", 150.86, 4.898e6, -0.00219,
                        {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 8, 21);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {argon}, {}, {2}, 8);
  ASSERT_TRUE(system.hasNumberOfBlocks(8));
  system.propertyElasticConstantsFluctuation = PropertyElasticConstantsFluctuation(5);
  EXPECT_FALSE(system.hasNumberOfBlocks(8));
  system.setNumberOfBlocks(8);
  EXPECT_EQ(system.propertyElasticConstantsFluctuation->numberOfBlocks, 8uz);
  EXPECT_TRUE(system.hasNumberOfBlocks(8));
}

// The thermostat, the barostat and the correlation functions keep a copy of the time step: set on the system alone,
// the chains were integrated with another time step than the particles.
TEST(SIMULATION_SETUP_FIXES, the_time_step_of_a_system_is_that_of_its_thermostat_and_barostat)
{
  const ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}}, {VDWParameters(119.8, 3.4)},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component argon(forceField, "argon", 150.86, 4.898e6, -0.00219,
                        {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {argon}, {}, {2}, 5);
  system.setThermostat(Thermostat(300.0, system.timeStep, 3, 0, 5, 5, 0.15));
  system.setThermobarostat(Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic,
                                          MonoclinicAngleType::Beta, 300.0, 1e5, system.timeStep, 3));
  system.setTimeStep(0.002);
  EXPECT_EQ(system.timeStep, 0.002);
  EXPECT_EQ(system.thermostat->timeStep, 0.002);
  EXPECT_EQ(system.thermobarostat->timeStep, 0.002);
  EXPECT_THROW(system.setTimeStep(0.0), std::invalid_argument);
}

// A driver without systems reads the components of the first one at every stage, and a LAMMPS writer every zero
// cycles divides by zero: both are refused when they are made.
TEST(SIMULATION_SETUP_FIXES, drivers_without_systems_and_writers_without_an_interval_are_refused)
{
  const SimulationSchedule schedule{20, 0, 3, 3, 1000, 10000, 5000, 5000};
  EXPECT_THROW((MonteCarlo(schedule, std::vector<System>{}, 47uz, 5, false)), std::invalid_argument);
  EXPECT_THROW((MolecularDynamics(schedule, std::vector<System>{}, 47uz, 5, false)), std::invalid_argument);
  EXPECT_THROW(WriteLammpsData(0, 0), std::invalid_argument);
}
