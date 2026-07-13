#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;
import forcefield;
import system;
import simulationbox;
import symmetric_eigensolver;
import minimization_options;
import minimization_baker;
import minimization_dof_layout;
import minimization_generalized_coordinates;
import minimization;
import sample_movies;

namespace
{
void setRigidMoleculePositions(System &system, std::size_t moleculeIndex)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  const Component &component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[localAtom].position;
  }
}
}  // namespace

TEST(minimization_driver, lapack_symmetric_eigensystem_has_small_residual)
{
  constexpr std::size_t n = 3;
  const std::array<double, n * n> matrix = {4.0, 1.0, -2.0, 1.0, 2.0, 0.5, -2.0, 0.5, 3.0};
  const SymmetricEigenSystem result = diagonalizeSymmetric(matrix, n);

  ASSERT_EQ(result.eigenvalues.size(), n);
  EXPECT_LE(result.eigenvalues[0], result.eigenvalues[1]);
  EXPECT_LE(result.eigenvalues[1], result.eigenvalues[2]);
  for (std::size_t mode = 0; mode < n; ++mode)
  {
    double norm = 0.0;
    for (std::size_t row = 0; row < n; ++row)
    {
      norm += result.eigenvector(row, mode) * result.eigenvector(row, mode);
      double residual = -result.eigenvalues[mode] * result.eigenvector(row, mode);
      for (std::size_t column = 0; column < n; ++column)
      {
        residual += matrix[row * n + column] * result.eigenvector(column, mode);
      }
      EXPECT_NEAR(residual, 0.0, 1e-12);
    }
    EXPECT_NEAR(norm, 1.0, 1e-12);
  }
}

TEST(minimization_driver, baker_step_filters_zero_modes_and_clips)
{
  MinimizationOptions options{};
  options.maximumStepLength = 0.25;
  options.convergenceFactor = 0.0;
  options.minimumEigenvalue = 1e-3;

  const std::array<double, 4> hessian = {0.0, 0.0, 0.0, 2.0};
  const std::array<double, 2> gradient = {1.0, 2.0};
  const BakerStep step = computeBakerStep(hessian, gradient, options);

  EXPECT_EQ(step.zeroModes, 1u);
  EXPECT_EQ(step.negativeModes, 0u);
  EXPECT_NEAR(step.displacement[0], 0.0, 1e-14);
  EXPECT_LT(step.displacement[1], 0.0);
  EXPECT_NEAR(step.stepNorm, 0.25, 1e-12);
}

TEST(minimization_driver, baker_positive_definite_step_and_convergence_classification)
{
  MinimizationOptions options{};
  options.maximumStepLength = 10.0;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-6;
  options.maxGradientTolerance = 1e-6;
  const std::array<double, 4> hessian = {2.0, 0.0, 0.0, 4.0};
  const std::array<double, 2> gradient = {0.2, -0.4};
  const BakerStep step = computeBakerStep(hessian, gradient, options);
  EXPECT_FALSE(step.converged);
  EXPECT_EQ(step.negativeModes, 0u);
  EXPECT_LT(step.displacement[0] * gradient[0] + step.displacement[1] * gradient[1], 0.0);

  const std::array<double, 2> convergedGradient = {1e-8, -1e-8};
  const BakerStep converged = computeBakerStep(hessian, convergedGradient, options);
  EXPECT_TRUE(converged.converged);
  EXPECT_EQ(converged.stepNorm, 0.0);
}

TEST(minimization_driver, baker_step_follows_negative_mode_downhill)
{
  MinimizationOptions options{};
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  const std::array<double, 1> hessian = {-1.0};
  const std::array<double, 1> gradient = {0.3};
  const BakerStep step = computeBakerStep(hessian, gradient, options);

  EXPECT_EQ(step.negativeModes, 1u);
  EXPECT_LT(step.shift, -1.0);
  EXPECT_LT(step.displacement[0], 0.0);
  EXPECT_LE(step.stepNorm, options.maximumStepLength);
}

TEST(minimization_driver, generalized_displacement_updates_rigid_quaternion_chart)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {1}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(8.0, 9.0, 10.0);
  system.moleculeData[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  setRigidMoleculePositions(system, 0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  std::vector<double> displacement(layout.numDofs(), 0.0);
  const std::size_t com = *layout.rigidMoleculeDof(0, RigidDof::ComX);
  const std::size_t orientation = *layout.rigidMoleculeDof(0, RigidDof::OriX);
  displacement[com] = 0.2;
  displacement[orientation + 1] = 0.1;
  applyGeneralizedDisplacement(system, layout, displacement);

  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.x, 8.2, 1e-14);
  const double3 oxygenOffset = system.spanOfMoleculeAtoms()[0].position -
                               system.moleculeData[0].centerOfMassPosition;
  EXPECT_NEAR(std::sqrt(double3::dot(oxygenOffset, oxygenOffset)), 1.149, 1e-12);
  EXPECT_GT(std::abs(oxygenOffset.x), 0.05);
}

TEST(minimization_driver, rigid_methane_pair_reaches_local_minimum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                         {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(10.0, 12.5, 12.5);
  system.moleculeData[1].centerOfMassPosition = double3(13.0, 12.5, 12.5);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-7;
  options.maxGradientTolerance = 1e-6;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());

  ASSERT_EQ(minimization.results.size(), 1u);
  EXPECT_TRUE(minimization.results[0].converged);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);
  EXPECT_EQ(minimization.results[0].negativeModes, 0u);
}

TEST(minimization_driver, writes_pdb_movie_during_minimization)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                         {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(10.0, 12.5, 12.5);
  system.moleculeData[1].centerOfMassPosition = double3(13.0, 12.5, 12.5);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);
  system.forceField.pseudoAtoms[static_cast<std::size_t>(system.components[0].atoms[0].type)].printToPDB = true;
  system.setSamplePDBMovie(SampleMovie(0, 1, false, std::nullopt));

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-7;
  options.maxGradientTolerance = 1e-6;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());

  std::ifstream movie("movies/movie.s0.pdb");
  ASSERT_TRUE(movie.is_open());
  const std::string contents((std::istreambuf_iterator<char>(movie)), std::istreambuf_iterator<char>());
  EXPECT_TRUE(contents.contains("ATOM"));
  EXPECT_TRUE(contents.contains("ENDMDL"));
  EXPECT_GT(std::ranges::count(contents, '\n'), 5u);
}

TEST(minimization_driver, flexible_methane_pair_reaches_local_minimum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                         {2}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(10.0, 12.5, 12.5);
  system.spanOfMoleculeAtoms()[1].position = double3(13.0, 12.5, 12.5);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-7;
  options.maxGradientTolerance = 1e-6;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());

  ASSERT_EQ(minimization.results.size(), 1u);
  EXPECT_TRUE(minimization.results[0].converged);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);
  EXPECT_EQ(minimization.results[0].negativeModes, 0u);
}
