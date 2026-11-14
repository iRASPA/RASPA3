#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;
import framework;
import int3;
import bond_potential;
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
import cross_links;

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

// Below a negative lowest eigenvalue the secular residual sum_i g_i^2/(lambda - e_i) - lambda falls monotonically,
// with its root arbitrarily close below the eigenvalue for a small gradient along that mode and far below it for a
// large one. The shift has to be the root, not an endpoint of a fixed bracket.
TEST(minimization_driver, baker_shift_solves_the_secular_equation_for_small_and_large_gradients)
{
  MinimizationOptions options{};
  options.maximumStepLength = 1.0e6;
  options.convergenceFactor = 0.0;
  const std::array<double, 4> hessian = {-1.0, 0.0, 0.0, 2.0};
  for (const double gradientAlongTheNegativeMode : {0.01, 1000.0})
  {
    SCOPED_TRACE(gradientAlongTheNegativeMode);
    const std::array<double, 2> gradient = {gradientAlongTheNegativeMode, 1.0};
    const BakerStep step = computeBakerStep(hessian, gradient, options);
    EXPECT_LT(step.shift, -1.0);
    const double residual = gradient[0] * gradient[0] / (step.shift + 1.0) + 1.0 / (step.shift - 2.0) - step.shift;
    EXPECT_NEAR(residual, 0.0, 1e-8 * std::max(1.0, std::abs(step.shift)));
  }
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
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {1}, 5);
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
  const double3 oxygenOffset = system.spanOfMoleculeAtoms()[0].position - system.moleculeData[0].centerOfMassPosition;
  EXPECT_NEAR(std::sqrt(double3::dot(oxygenOffset, oxygenOffset)), 1.149, 1e-12);
  EXPECT_GT(std::abs(oxygenOffset.x), 0.05);
}

TEST(minimization_driver, rigid_methane_pair_reaches_local_minimum)
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
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
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

TEST(minimization_driver, flexible_framework_bond_reaches_local_minimum)
{
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{1.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, false, false, false);
  const SimulationBox box(10.0, 10.0, 10.0);
  const Atom atomA({0.20, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  const Atom atomB({0.40, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "flexible-dimer", box, 1, {atomA, atomB}, {atomA, atomB}, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds.emplace_back(std::array<std::size_t, 2>{0, 1}, BondType::Harmonic,
                                                        std::vector<double>{1000.0, 1.5});
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

  ASSERT_EQ(minimization.results.size(), 1u);
  EXPECT_TRUE(minimization.results[0].converged);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);
  const std::span<const Atom> minimizedAtoms = minimization.systems[0].spanOfFrameworkAtoms();
  const double3 dr = minimizedAtoms[0].position - minimizedAtoms[1].position;
  EXPECT_NEAR(std::sqrt(double3::dot(dr, dr)), 1.5, 1.0e-8);
}

TEST(minimization_driver, flexible_methane_pair_reaches_local_minimum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
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

// The phonon primitive-cell reduction keyed the symmetry species on the pseudo-atom type alone: a CsCl-like cell of one
// type with charges +0.5 and -0.5 reduced to a one-atom cell (the body centre looked like a lattice translation) with
// the charge of the first atom. With (type, charge) species the two atoms stay inequivalent with their own charges.
TEST(minimization_driver, primitive_reduction_keeps_atoms_of_one_type_with_distinct_charges)
{
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, false, false, false);
  const SimulationBox box(6.0, 6.0, 6.0);
  const Atom atomA({0.0, 0.0, 0.0}, 0.5, 1.0, 0, 0, 0, 0, true);
  const Atom atomB({0.5, 0.5, 0.5}, -0.5, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "cscl", box, 1, {atomA, atomB}, {atomA, atomB}, int3(1, 1, 1));
  framework.rigid = false;
  framework.frameworkDefinitionName = "cscl";
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {}, {}, {}, 5);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 20;
  options.printEvery = 100;
  options.computePhononDispersion = true;
  options.phononDispersionPointsPerSegment = 2;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());
  ASSERT_TRUE(minimization.results[0].converged);

  const std::span<const Atom> atoms = minimization.systems[0].spanOfFrameworkAtoms();
  ASSERT_EQ(atoms.size(), 2u);
  EXPECT_NEAR(atoms[0].charge + atoms[1].charge, 0.0, 1e-12);
  EXPECT_NEAR(std::abs(atoms[0].charge - atoms[1].charge), 1.0, 1e-12);
}

// Minimization has no cross-link or polymer-link terms and refuses them instead of pulling linked atoms apart.
TEST(minimization_driver, minimization_refuses_cross_links)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
  system.crossLinks.bondTypes.push_back(CrossLinkBondType{});
  system.crossLinks.bondTypes.back().bond = BondPotential({0, 0}, BondType::Harmonic, {1000.0, 3.0});
  system.crossLinks.links.push_back(CrossLink{{0, 0, 0}, {0, 1, 0}, 0});
  system.crossLinks.rebuildIndex();
  MinimizationOptions options{};
  options.printEvery = 100;
  Minimization minimization(options, {system}, false);
  EXPECT_THROW(minimization.run(), std::runtime_error);
}

// Three methanes on a line are a symmetric saddle of a fixed-cell minimization: the gradient lies along the line, the
// bending mode is negative with no gradient along it, and the rational-function step shrinks to zero there. The
// negative-mode escape (before: only for a variable cell) bends the chain off the hilltop into a minimum.
TEST(minimization_driver, fixed_cell_minimization_escapes_a_symmetric_saddle)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {3}, 5);
  for (std::size_t i = 0; i < 3; ++i)
  {
    system.moleculeData[i].centerOfMassPosition = double3(8.3 + 4.2 * static_cast<double>(i), 12.5, 12.5);
    setRigidMoleculePositions(system, i);
  }

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
  EXPECT_EQ(minimization.results[0].negativeModes, 0u);

  // off the line: the middle methane left the axis through the outer two
  const std::span<const Molecule> molecules = minimization.systems[0].moleculeData;
  const double3 axis = molecules[2].centerOfMassPosition - molecules[0].centerOfMassPosition;
  const double3 middle = molecules[1].centerOfMassPosition - molecules[0].centerOfMassPosition;
  const double3 offAxis = middle - (double3::dot(middle, axis) / double3::dot(axis, axis)) * axis;
  EXPECT_GT(std::sqrt(double3::dot(offAxis, offAxis)), 0.5);
}

// The derivative evaluation has no external-field term, so a minimization with an external field is refused for a
// fixed cell as it was for a variable cell.
TEST(minimization_driver, minimization_refuses_an_external_field_for_a_fixed_cell)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(10.0, 12.5, 12.5);
  system.moleculeData[1].centerOfMassPosition = double3(13.0, 12.5, 12.5);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);
  system.hasExternalField = true;

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 5;
  options.printEvery = 5;
  Minimization minimization(options, {system}, false);
  try
  {
    minimization.run();
    ADD_FAILURE() << "the minimization ran";
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string(error.what()).find("external field"), std::string::npos) << error.what();
  }
}

// A soft positive mode under a large gradient (|g_1| / e_1 = 10^4): the shift is the root of the secular equation
// below the lowest eigenvalue, as for a negative one. The fixed-point iteration used for positive eigenvalues settled
// into a two-cycle here and returned whichever iterate it stopped on.
TEST(minimization_driver, baker_shift_solves_the_secular_equation_on_a_soft_positive_mode)
{
  MinimizationOptions options{};
  options.maximumStepLength = 1.0e6;
  const std::array<double, 4> hessian = {0.01, 0.0, 0.0, 2.0};
  const std::array<double, 2> gradient = {100.0, 1.0};
  const BakerStep step = computeBakerStep(hessian, gradient, options);
  EXPECT_LT(step.shift, 0.01);
  const double residual =
      gradient[0] * gradient[0] / (step.shift - 0.01) + gradient[1] * gradient[1] / (step.shift - 2.0) - step.shift;
  EXPECT_NEAR(residual, 0.0, 1e-8 * std::max(1.0, std::abs(step.shift)));
}
