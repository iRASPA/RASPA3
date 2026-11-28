#include <gtest/gtest.h>

import std;

import finite_value;
import int3;
import double3;
import double3x3;
import units;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import mc_moves_move_types;
import mc_moves_probabilities;
import generalized_hessian;
import minimization_hessian_scatter;
import torsion_potential;
import skdefinitions;
import skrotationmatrix;
import skseitzintegermatrix;
import skseitzmatrix;
import skspacegroupsetting;
import skrotationalchangeofbasis;
import skspacegroupdatabase;
import skspacegroup;

namespace
{
Component makeMethane(const ForceField& forceField)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  return Component(forceField, "methane", 190.564, 4599200, 0.01142,
                   {Atom({0, 0, 0}, 0.0, 1.0, 0, 6, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                   false);
}

// the operations of a tabulated setting (translations in fractions of the cell)
std::vector<SKSeitzMatrix> operationsOfSetting(std::size_t hallNumber)
{
  std::vector<SKSeitzMatrix> operations;
  for (const SKSeitzIntegerMatrix& m : SKSpaceGroup(hallNumber).spaceGroupSetting().fullSeitzMatrices().operations)
  {
    operations.emplace_back(m.rotation, double3(m.translation.x, m.translation.y, m.translation.z) / 24.0);
  }
  return operations;
}
}  // namespace

// An atom and its own periodic image (the image lists of a cell narrower than twice the cutoff) move together: the
// separation is the lattice vector whatever the displacement, and the pair adds nothing to the position-position
// Hessian. One +B on the diagonal gave the uniform translation of the crystal a curvature.
TEST(BUG_REPORT_FIXES, periodic_self_image_pair_adds_no_position_position_hessian)
{
  GeneralizedHessian hessian(6, 0);
  Minimization::scatterAtomicPositionPositionByDof(hessian, 3, 3, -0.3, 0.7, double3(1.0, 2.0, 3.0));
  for (std::size_t i = 0; i < 6; ++i)
  {
    for (std::size_t j = 0; j < 6; ++j) EXPECT_EQ(hessian(i, j), 0.0);
  }

  // two atoms: the blocks +B, +B, -B, -B, whose sum over the columns of a row (a uniform translation) is zero
  Minimization::scatterAtomicPositionPositionByDof(hessian, 0, 3, -0.3, 0.7, double3(1.0, 2.0, 3.0));
  for (std::size_t i = 0; i < 6; ++i)
  {
    for (std::size_t axis = 0; axis < 3; ++axis) EXPECT_NEAR(hessian(i, axis) + hessian(i, 3 + axis), 0.0, 1.0e-14);
  }
  EXPECT_NE(hessian(0, 0), 0.0);
}

// Three atoms of a torsion on a line: the angle is undefined. The gradients divided by the length of the
// perpendicular arm, zero, and came out as NaN; they are zero, and the energy is that of the energy routine.
TEST(BUG_REPORT_FIXES, collinear_torsion_has_finite_energy_and_zero_gradient)
{
  const TorsionPotential torsion({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32});
  const double3 a(0.0, 0.0, 0.0), b(1.5, 0.0, 0.0), c(3.0, 0.0, 0.0), d(4.0, 1.0, 0.0);

  const auto [energy, gradient, strain] = torsion.potentialEnergyGradientStrain(a, b, c, d);
  EXPECT_TRUE(isFiniteValue(energy));
  EXPECT_NEAR(energy, torsion.calculateEnergy(a, b, c, d), 1.0e-10 * std::max(1.0, std::abs(energy)));
  for (const double3& g : gradient)
  {
    EXPECT_EQ(g.x, 0.0);
    EXPECT_EQ(g.y, 0.0);
    EXPECT_EQ(g.z, 0.0);
  }
  EXPECT_TRUE(isFiniteValue(strain.ax) && isFiniteValue(strain.by) && isFiniteValue(strain.cz));

  // a bent geometry is untouched: a gradient that sums to zero and is not zero
  const auto [bentEnergy, bentGradient, bentStrain] =
      torsion.potentialEnergyGradientStrain(double3(0.0, 1.0, 0.0), b, c, double3(4.0, 0.5, 0.8));
  const double3 sum = bentGradient[0] + bentGradient[1] + bentGradient[2] + bentGradient[3];
  EXPECT_NEAR(sum.x, 0.0, 1.0e-9);
  EXPECT_NEAR(sum.y, 0.0, 1.0e-9);
  EXPECT_NEAR(sum.z, 0.0, 1.0e-9);
  EXPECT_GT(double3::dot(bentGradient[0], bentGradient[0]), 0.0);
  EXPECT_TRUE(isFiniteValue(bentEnergy));
}

// The ideal-gas Rosenbluth weight of a long chain is not a double: exp(-800) is zero, and the acceptance rules took
// its logarithm, minus infinity. A computed weight keeps the logarithm it was computed as.
TEST(BUG_REPORT_FIXES, ideal_gas_rosenbluth_weight_is_used_through_its_logarithm)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  Component methane = makeMethane(forceField);

  EXPECT_EQ(methane.logIdealGasRosenbluthWeight(), 0.0);

  methane.idealGasRosenbluthWeight = 0.5;
  EXPECT_NEAR(methane.logIdealGasRosenbluthWeight(), std::log(0.5), 1.0e-14);

  methane.idealGasRosenbluthWeight = std::exp(-800.0);
  methane.idealGasRosenbluthWeightEstimate = Component::IdealGasRosenbluthWeightEstimate{
      .temperature = 300.0, .logRosenbluthWeight = -800.0, .relativeError = 0.01, .numberOfGrows = 100, .converged = true};
  EXPECT_EQ(methane.idealGasRosenbluthWeight.value(), 0.0);
  EXPECT_EQ(methane.logIdealGasRosenbluthWeight(), -800.0);

  // a weight set by hand after one was computed is the weight
  methane.idealGasRosenbluthWeight = 0.25;
  EXPECT_NEAR(methane.logIdealGasRosenbluthWeight(), std::log(0.25), 1.0e-14);
}

// A system built in code: no initial molecule counts is no molecules (the counts were read past the end of the
// empty list), fewer of them than components and a temperature that is not positive are refused.
TEST(BUG_REPORT_FIXES, system_checks_the_arguments_the_input_reader_checks_for_a_file)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  const Component methane = makeMethane(forceField);
  const SimulationBox box(25.0, 25.0, 25.0);

  const System system(forceField, box, false, 300.0, 1.0e4, 1.0, {}, {methane}, {}, {}, 5);
  EXPECT_EQ(system.numberOfMolecules(), 0uz);

  EXPECT_THROW(System(forceField, box, false, 300.0, 1.0e4, 1.0, {}, {methane, methane}, {}, {1}, 5),
               std::runtime_error);
  EXPECT_THROW(System(forceField, box, false, 0.0, 1.0e4, 1.0, {}, {methane}, {}, {0}, 5), std::runtime_error);
  EXPECT_THROW(System(forceField, box, false, -10.0, 1.0e4, 1.0, {}, {methane}, {}, {0}, 5), std::runtime_error);
}

// A precision in Angstrom is a distance: with c = 30 Angstrom the screw translation 1/6 of P 6_1 is 5 Angstrom, and
// compared as a fraction with a precision of 0.2 it passed for no translation, P 6.
TEST(BUG_REPORT_FIXES, screw_translation_is_measured_in_the_lattice)
{
  const std::size_t p6 = SKSpaceGroupDataBase::spaceGroupHallData[168].front();
  const std::size_t p61 = SKSpaceGroupDataBase::spaceGroupHallData[169].front();
  const std::vector<SKSeitzMatrix> operations = operationsOfSetting(p61);
  ASSERT_EQ(operations.size(), 6uz);

  const double3x3 hexagonal(double3(10.0, 0.0, 0.0), double3(-5.0, 5.0 * std::sqrt(3.0), 0.0),
                            double3(0.0, 0.0, 30.0));
  EXPECT_FALSE(SKSpaceGroup::getOriginShift(p6, Centring::primitive, SKRotationalChangeOfBasis::identity, operations,
                                            hexagonal, 0.2)
                   .has_value());
  EXPECT_TRUE(SKSpaceGroup::getOriginShift(p61, Centring::primitive, SKRotationalChangeOfBasis::identity, operations,
                                           hexagonal, 0.2)
                  .has_value());
  EXPECT_FALSE(SKSpaceGroup::matchSpaceGroup(p6, hexagonal, Centring::primitive, operations, 0.2).has_value());
  EXPECT_TRUE(SKSpaceGroup::matchSpaceGroup(p61, hexagonal, Centring::primitive, operations, 0.2).has_value());
}
