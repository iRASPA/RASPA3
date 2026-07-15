#include <gtest/gtest.h>

import std;

import units;
import double3;
import double3x3;
import simd_quatd;
import int3;
import atom;
import molecule;
import component;
import framework;
import bond_potential;
import forcefield;
import system;
import simulationbox;
import minimization_options;
import minimization;
import normal_modes;

namespace
{
void setRigidMoleculePositions(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[localAtom].position;
  }
}
}  // namespace

TEST(normal_modes, harmonic_framework_dimer_matches_analytic_frequency)
{
  const double mass = 12.0;
  const double forceConstantKelvin = 1000.0;
  ForceField forceField({{"C", false, mass, 0.0, 0.0, 6, false}}, {{1.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, false, false, false);
  const SimulationBox box(10.0, 10.0, 10.0);
  const Atom atomA({0.20, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  const Atom atomB({0.40, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "flexible-dimer", box, 1, {atomA, atomB}, {atomA, atomB}, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds.emplace_back(std::array<std::size_t, 2>{0, 1}, BondType::Harmonic,
                                                        std::vector<double>{forceConstantKelvin, 1.5});
  framework.intraMolecularImageShifts.bonds.push_back({int3{}, int3{}});
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {}, {}, {}, 5);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 100;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1.0e-8;
  options.maxGradientTolerance = 1.0e-7;
  options.printEvery = 100;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());
  ASSERT_TRUE(minimization.results[0].converged);

  const NormalModesResult modes = computeNormalModes(minimization.systems[0]);
  ASSERT_EQ(modes.numberOfModes, 6u);
  ASSERT_EQ(modes.eigenvalues.size(), 6u);
  EXPECT_EQ(modes.negativeModes, 0u);
  EXPECT_EQ(modes.discardedRotationalDofs, 0u);

  // Three translations plus two rotations of the linear dimer are zero modes.
  EXPECT_EQ(modes.zeroModes, 5u);

  // U = (1/2) k (r - r0)^2 gives a single vibration with omega^2 = k / mu, mu = m/2.
  const double forceConstant = forceConstantKelvin * Units::KelvinToEnergy;
  const double expectedEigenvalue = forceConstant * (2.0 / mass);
  EXPECT_NEAR(modes.eigenvalues.back(), expectedEigenvalue, 1.0e-8 * expectedEigenvalue);
}

TEST(normal_modes, rigid_methane_pair_lennard_jones_stretch_frequency)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(10.0, 12.5, 12.5);
  system.moleculeData[1].centerOfMassPosition = double3(13.0, 12.5, 12.5);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-8;
  options.maxGradientTolerance = 1e-7;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());
  ASSERT_TRUE(minimization.results[0].converged);

  const NormalModesResult modes = computeNormalModes(minimization.systems[0]);
  ASSERT_EQ(modes.numberOfModes, 12u);
  EXPECT_EQ(modes.negativeModes, 0u);

  // Single-bead rigid molecules have zero inertia in all directions.
  EXPECT_EQ(modes.discardedRotationalDofs, 6u);

  // 6 discarded orientation DOFs + 3 translations + 2 transverse modes are zero;
  // one radial Lennard-Jones stretch mode remains.
  EXPECT_EQ(modes.zeroModes, 11u);

  // At the minimum r* = 2^(1/6) sigma: U''(r*) = 72 epsilon / r*^2, omega^2 = U'' (2/m).
  const double epsilon = 158.5 * Units::KelvinToEnergy;
  const double sigma = 3.72;
  const double rStar = std::pow(2.0, 1.0 / 6.0) * sigma;
  const double moleculeMass = minimization.systems[0].moleculeData[0].mass;
  const double expectedEigenvalue = (72.0 * epsilon / (rStar * rStar)) * (2.0 / moleculeMass);
  EXPECT_NEAR(modes.eigenvalues.back(), expectedEigenvalue, 1.0e-6 * expectedEigenvalue);

  const std::vector<double> frequencies = normalModeFrequencies(modes);
  ASSERT_EQ(frequencies.size(), 12u);
  EXPECT_GT(frequencies.back(), 0.0);
}

TEST(normal_modes, writes_one_pdb_movie_per_mode)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(10.0, 12.5, 12.5);
  system.moleculeData[1].centerOfMassPosition = double3(13.0, 12.5, 12.5);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-8;
  options.maxGradientTolerance = 1e-7;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());
  ASSERT_TRUE(minimization.results[0].converged);

  const NormalModesResult modes = computeNormalModes(minimization.systems[0]);

  const std::filesystem::path directory = "normal_modes_test_output";
  std::filesystem::remove_all(directory);
  const std::size_t periods = 2;
  const std::size_t pointsPerPeriod = 8;
  ASSERT_NO_THROW(writeNormalModeMovies(minimization.systems[0], modes, 0, periods, pointsPerPeriod, 0.4, directory));

  std::size_t fileCount = 0;
  for (const auto& entry : std::filesystem::directory_iterator(directory))
  {
    if (entry.path().extension() == ".pdb") ++fileCount;
  }
  EXPECT_EQ(fileCount, modes.numberOfModes);

  std::ifstream movie(directory / "mode_0011.s0.pdb");
  ASSERT_TRUE(movie.is_open());
  const std::string contents((std::istreambuf_iterator<char>(movie)), std::istreambuf_iterator<char>());
  EXPECT_TRUE(contents.contains("MODEL"));
  EXPECT_TRUE(contents.contains("ATOM"));
  EXPECT_TRUE(contents.contains("ENDMDL"));
  EXPECT_EQ(std::ranges::count(contents, 'M') > 0, true);
  // One MODEL/ENDMDL pair per frame.
  const std::string endmdl = "ENDMDL";
  std::size_t frames = 0;
  for (std::size_t pos = contents.find(endmdl); pos != std::string::npos; pos = contents.find(endmdl, pos + 1))
  {
    ++frames;
  }
  EXPECT_EQ(frames, periods * pointsPerPeriod);

  std::filesystem::remove_all(directory);
}
