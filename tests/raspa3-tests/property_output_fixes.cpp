#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import double3x3;
import simd_quatd;
import units;
import atom;
import molecule;
import forcefield;
import component;
import connectivity_table;
import intra_molecular_potentials;
import simulationbox;
import system;
import running_energy;
import averages;
import integrators_update;
import lammps_topology;
import property_molecule_properties;
import property_number_of_molecules_evolution;
import property_conserved_energy_evolution;
import property_energy;

namespace
{
class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::current_path(path);
  }
  ~ScopedCurrentPath() { std::filesystem::current_path(original_); }

 private:
  std::filesystem::path original_;
};

std::vector<std::string> linesOf(const std::filesystem::path& path)
{
  std::ifstream stream(path);
  std::vector<std::string> lines;
  for (std::string line; std::getline(stream, line);) lines.push_back(line);
  return lines;
}
}  // namespace

TEST(PROPERTY_OUTPUT_FIXES, lammps_velocities_of_a_rigid_molecule_include_its_rotation)
{
  // The exported atom velocities of a rigid molecule are those of the rigid body, v_com + omega x (r - r_com): the
  // displacement of its atoms over a short free rotation (and translation) of the body.
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}, {"H", false, 1.0, 0.0, 0.0, 1, false}},
                        {{80.0, 3.0}, {10.0, 2.0}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false,
                        false, false);
  Component water(
      forceField, "water", 647.0, 22.064e6, 0.344,
      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false), Atom({0.8, 0.6, 0.0}, 0.0, 1.0, 0, 1, 0, false, false),
       Atom({-0.8, 0.6, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
      ConnectivityTable(3), Potentials::IntraMolecularPotentials{}, 5, 21);
  water.rigid = true;
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {water}, {}, {1}, 5);

  Molecule& molecule = system.moleculeData[0];
  molecule.centerOfMassPosition = double3(12.0, 11.0, 13.0);
  molecule.orientation = simd_quatd::fromAxisAngle(0.7, double3(1.0, 2.0, 3.0).normalized());
  molecule.velocity = double3(0.3, -0.2, 0.5);
  // the quaternion momentum of a body-frame angular velocity (as Integrators::initializeVelocities sets it)
  const double3 I = system.components[0].inertiaVector;
  const double3 w(1.3, -0.7, 2.1);
  const simd_quatd q = molecule.orientation;
  molecule.orientationMomentum.ix = 2.0 * (q.r * (I.x * w.x) - q.iz * (I.y * w.y) + q.iy * (I.z * w.z));
  molecule.orientationMomentum.iy = 2.0 * (q.iz * (I.x * w.x) + q.r * (I.y * w.y) - q.ix * (I.z * w.z));
  molecule.orientationMomentum.iz = 2.0 * (-q.iy * (I.x * w.x) + q.ix * (I.y * w.y) + q.r * (I.z * w.z));
  molecule.orientationMomentum.r = 2.0 * (-q.ix * (I.x * w.x) - q.iy * (I.y * w.y) - q.iz * (I.z * w.z));
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);

  const LAMMPS::Topology topology = LAMMPS::buildTopology(
      system.components, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.moleculeData,
      system.simulationBox, system.forceField, system.numberOfIntegerMoleculesPerComponent, system.framework);
  ASSERT_EQ(topology.atoms.size(), 3uz);

  constexpr double dt = 1.0e-7;
  std::vector<Molecule> moved(system.moleculeData.begin(), system.moleculeData.end());
  Integrators::noSquishFreeRotorOrderTwo(moved, system.components, dt);
  moved[0].centerOfMassPosition += dt * moved[0].velocity;
  std::vector<Atom> movedAtoms(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
  Integrators::createCartesianPositions(moved, movedAtoms, system.components);

  const double toAngstromPerFemtosecond = Units::VelocityConversionFactor * 1e-5;
  for (std::size_t i = 0; i < 3; ++i)
  {
    const double3 expected =
        (movedAtoms[i].position - system.spanOfMoleculeAtoms()[i].position) / dt * toAngstromPerFemtosecond;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      EXPECT_NEAR(topology.atoms[i].velocity[axis], expected[axis], 1e-6 * std::max(1.0, std::abs(expected[axis])))
          << "atom " << i << ", axis " << axis;
    }
  }
}

TEST(PROPERTY_OUTPUT_FIXES, molecule_property_errors_use_the_blocks_with_samples)
{
  // Five blocks of which three hold samples: the confidence interval is t(2) s / sqrt(3) of those three blocks.
  PropertyMoleculeProperties property;
  property.numberOfBlocks = 5;
  property.numberOfComponents = 1;
  property.numberOfCounts = {{1.0}, {1.0}, {1.0}, {0.0}, {0.0}};
  const std::vector<std::vector<std::vector<std::vector<double>>>> histogram = {
      {{{1.0}}}, {{{2.0}}}, {{{4.0}}}, {{{0.0}}}, {{{0.0}}}};

  const auto [bins, average, error] = property.result(histogram, 0, 0, 1.0, 0.0);
  ASSERT_EQ(average.size(), 1uz);
  EXPECT_NEAR(average[0], 7.0 / 3.0, 1e-12);
  const double sumOfSquares =
      std::pow(1.0 - 7.0 / 3.0, 2) + std::pow(2.0 - 7.0 / 3.0, 2) + std::pow(4.0 - 7.0 / 3.0, 2);
  EXPECT_NEAR(error[0], studentTCriticalValue(2) * std::sqrt(sumOfSquares / 2.0) / std::sqrt(3.0), 1e-12);
}

TEST(PROPERTY_OUTPUT_FIXES, evolution_files_hold_every_sample_in_kelvin_and_survive_zero_components)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  // ten cycles, sampled every third (cycles 0, 3, 6, 9) and written only at multiples of ten: during the run the
  // file is never complete, the last write holds all four samples
  PropertyNumberOfMoleculesEvolution numberOfMolecules(10, 1, 3, 10);
  PropertyConservedEnergyEvolution conservedEnergy(10, 3, 10);
  RunningEnergy energy;
  energy.moleculeMoleculeVDW = 2.0;
  for (std::size_t cycle = 0; cycle < 10; ++cycle)
  {
    numberOfMolecules.addSample(cycle, {cycle});
    conservedEnergy.addSample(cycle, energy);
  }
  numberOfMolecules.writeOutput(0, 10, true);
  conservedEnergy.writeOutput(0, 10, true);

  const std::vector<std::string> molecules =
      linesOf("number_of_molecules_evolution/number_of_molecules_evolution.s0.txt");
  ASSERT_EQ(molecules.size(), 4uz);
  EXPECT_EQ(molecules[3], "9 ");

  const std::vector<std::string> energies = linesOf("conserved_energy/conserved_energy.s0.txt");
  ASSERT_EQ(energies.size(), 4uz);
  std::istringstream row(energies[0]);
  double index{}, conserved{}, potential{};
  row >> index >> conserved >> potential;
  EXPECT_NEAR(potential, 2.0 * Units::EnergyToKelvin, 1e-9);

  // no components: nothing to write, and no crash
  PropertyNumberOfMoleculesEvolution noComponents(10, 0, 3, 10);
  EXPECT_NO_THROW(noComponents.writeOutput(1, 10, true));
}

TEST(PROPERTY_OUTPUT_FIXES, energy_repr_returns_its_text)
{
  const PropertyEnergy energy(5, 0, 0, 1);
  EXPECT_FALSE(energy.repr().empty());
}
