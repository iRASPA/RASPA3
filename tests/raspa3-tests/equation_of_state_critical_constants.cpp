#include <gtest/gtest.h>

import std;

import atom;
import component;
import equation_of_states;
import forcefield;
import isotherm_bet;
import mc_moves_move_types;
import mc_moves_probabilities;
import parallel_tmmc;
import simulationbox;
import system;

// A cubic equation of state is built from a component's critical temperature and pressure. They are optional
// in a molecule file and nought when left out, and a component without them has no equation of state of its
// own: it is an ideal gas, with a fugacity coefficient of one unless the input gives it another. Dividing by
// the noughts instead made the coefficient, the compressibility and everything computed with them not a
// number, in an ordinary Monte Carlo run and in the drivers that recompute the coefficient at every state
// point alike.

namespace
{

constexpr double methaneCriticalTemperature = 190.564;
constexpr double methaneCriticalPressure = 4599200.0;
constexpr double methaneAcentricFactor = 0.01142;

// One swappable particle in an empty box, with the critical constants given, and with a fugacity coefficient
// of its own or none.
System makeSystem(double criticalTemperature, double criticalPressure, std::optional<double> fugacityCoefficient)
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 1, false}}, {{150.0, 3.7}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  const Component particle(forceField, "particle", criticalTemperature, criticalPressure, methaneAcentricFactor,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities,
                           fugacityCoefficient);
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {particle}, {}, {0}, 5);
}

// A number, and one a fugacity coefficient or a compressibility of a gas can be. A comparison is what tells
// a value that is not a number here: the build is one that may assume there are none.
bool sensible(double value) { return value > 0.1 && value < 10.0; }

}  // namespace

TEST(EQUATION_OF_STATE, a_component_without_critical_constants_is_an_ideal_gas)
{
  const std::vector<EquationOfState::FluidInput> alone = {{0.0, 0.0, 0.0, 1.0, true}};
  for (double pressure : {1.0, 1.0e5, 1.0e7})
  {
    const std::vector<EquationOfState::FluidResult> result = EquationOfState::computeFluidProperties(
        300.0, pressure, alone, EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals);
    ASSERT_EQ(result.size(), 1uz);
    ASSERT_TRUE(result.front().fugacityCoefficient.has_value());
    EXPECT_DOUBLE_EQ(result.front().fugacityCoefficient.value(), 1.0) << "at " << pressure << " Pa";
    EXPECT_DOUBLE_EQ(result.front().compressibility, 1.0) << "at " << pressure << " Pa";
  }

  // One constant without the other is no equation of state either.
  for (const auto& [criticalTemperature, criticalPressure] :
       std::vector<std::pair<double, double>>{{methaneCriticalTemperature, 0.0}, {0.0, methaneCriticalPressure}})
  {
    const std::vector<EquationOfState::FluidInput> half = {{criticalTemperature, criticalPressure, 0.0, 1.0, true}};
    const std::vector<EquationOfState::FluidResult> result = EquationOfState::computeFluidProperties(
        300.0, 1.0e5, half, EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals);
    ASSERT_TRUE(result.front().fugacityCoefficient.has_value());
    EXPECT_DOUBLE_EQ(result.front().fugacityCoefficient.value(), 1.0);
  }
}

// In a mixture it takes no part in the attraction and has no volume of its own, and leaves the component
// that has an equation of state with a coefficient that is a number.
TEST(EQUATION_OF_STATE, an_ideal_component_does_not_spoil_the_mixture_it_is_in)
{
  const std::vector<EquationOfState::FluidInput> mixture = {
      {methaneCriticalTemperature, methaneCriticalPressure, methaneAcentricFactor, 0.5, true},
      {0.0, 0.0, 0.0, 0.5, true}};
  const std::vector<EquationOfState::FluidResult> result = EquationOfState::computeFluidProperties(
      300.0, 1.0e6, mixture, EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals);
  ASSERT_EQ(result.size(), 2uz);
  for (const EquationOfState::FluidResult& component : result)
  {
    ASSERT_TRUE(component.fugacityCoefficient.has_value());
    EXPECT_TRUE(sensible(component.fugacityCoefficient.value())) << component.fugacityCoefficient.value();
    EXPECT_TRUE(sensible(component.compressibility)) << component.compressibility;
  }

  // And methane by itself is what it was: a little under one at room temperature and a bar.
  const std::vector<EquationOfState::FluidInput> methane = {
      {methaneCriticalTemperature, methaneCriticalPressure, methaneAcentricFactor, 1.0, true}};
  const std::vector<EquationOfState::FluidResult> pure = EquationOfState::computeFluidProperties(
      300.0, 1.0e5, methane, EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals);
  EXPECT_NEAR(pure.front().fugacityCoefficient.value(), 0.998, 0.001);
}

// A system is given its fugacity coefficients as it is built. Without critical constants the coefficient is
// one, or what the input gave; with them it is Peng-Robinson's, or what the input gave.
TEST(EQUATION_OF_STATE, a_system_without_critical_constants_has_a_coefficient_that_is_a_number)
{
  const System ideal = makeSystem(0.0, 0.0, std::nullopt);
  ASSERT_TRUE(ideal.components.front().swappable);
  ASSERT_TRUE(ideal.components.front().fugacityCoefficient.has_value());
  EXPECT_DOUBLE_EQ(ideal.components.front().fugacityCoefficient.value(), 1.0);
  EXPECT_DOUBLE_EQ(ideal.components.front().compressibility, 1.0);
  EXPECT_GT(ideal.components.front().bulkFluidDensity, 0.0);

  const System given = makeSystem(0.0, 0.0, 0.9);
  EXPECT_DOUBLE_EQ(given.components.front().fugacityCoefficient.value(), 0.9);
  EXPECT_DOUBLE_EQ(given.components.front().compressibility, 1.0);

  const System real = makeSystem(methaneCriticalTemperature, methaneCriticalPressure, std::nullopt);
  EXPECT_NEAR(real.components.front().fugacityCoefficient.value(), 0.998, 0.001);
}

// The drivers that run one system at many state points take the coefficient from Peng-Robinson at each of
// them, and put aside one the input gives: it would be right at one state point only. That holds where
// there is a Peng-Robinson to take it from. A component without critical constants has none, and the
// coefficient it was given, or the ideal gas's one, is all there is: it is kept, the same at every state
// point, where it used to be thrown away and replaced by something that was not a number.
TEST(EQUATION_OF_STATE, a_driver_keeps_the_coefficient_of_a_component_it_cannot_recompute)
{
  for (std::optional<double> givenCoefficient : {std::optional<double>{}, std::optional<double>{0.9}})
  {
    const double expected = givenCoefficient.value_or(1.0);

    System pinned = makeSystem(0.0, 0.0, givenCoefficient);
    for (double pressure : {1.0e3, 1.0e6})
    {
      pinSystemPengRobinsonPressure(pinned, pressure);
      ASSERT_TRUE(pinned.components.front().fugacityCoefficient.has_value());
      EXPECT_DOUBLE_EQ(pinned.components.front().fugacityCoefficient.value(), expected);
    }

    System walkerTemplate = makeSystem(0.0, 0.0, givenCoefficient);
    walkerTemplate.tmmc.minMacrostate = 0;
    walkerTemplate.tmmc.maxMacrostate = 20;
    ParallelTMMCParameters parameters;
    parameters.numberOfWindows = 2;
    parameters.randomSeed = 7;
    const ParallelTMMC simulation(std::move(walkerTemplate), {250.0, 350.0}, parameters);
    ASSERT_EQ(simulation.systems.size(), 4uz);
    for (const System& walker : simulation.systems)
    {
      ASSERT_TRUE(walker.components.front().fugacityCoefficient.has_value());
      EXPECT_DOUBLE_EQ(walker.components.front().fugacityCoefficient.value(), expected);
    }
  }

  // With critical constants the coefficient is recomputed, whatever was given: 0.9 is not what methane has
  // at a bar, and a megapascal gives another than a kilopascal does.
  System real = makeSystem(methaneCriticalTemperature, methaneCriticalPressure, 0.9);
  EXPECT_DOUBLE_EQ(real.components.front().fugacityCoefficient.value(), 0.9);
  pinSystemPengRobinsonPressure(real, 1.0e3);
  const double dilute = real.components.front().fugacityCoefficient.value();
  pinSystemPengRobinsonPressure(real, 1.0e6);
  const double dense = real.components.front().fugacityCoefficient.value();
  EXPECT_NEAR(dilute, 1.0, 1.0e-4);
  EXPECT_LT(dense, 0.99);
  EXPECT_GT(dense, 0.95);
}

// In the three-root regime the mixture is one bulk phase: every component takes the same root, the one with the
// lower molar Gibbs energy of the mixture. Methane with n-butane at 260 K and 5 bar: methane alone would prefer
// the vapor root and butane the liquid root, and the mixture's choice depends on the composition.
TEST(EQUATION_OF_STATE, a_mixture_takes_one_compressibility_root_for_all_its_components)
{
  constexpr double butaneCriticalTemperature = 425.125;
  constexpr double butaneCriticalPressure = 3796000.0;
  constexpr double butaneAcentricFactor = 0.2002;
  constexpr double temperature = 260.0;
  constexpr double pressure = 5.0e5;

  auto mixture = [&](double methaneFraction)
  {
    const std::vector<EquationOfState::FluidInput> inputs = {
        {methaneCriticalTemperature, methaneCriticalPressure, methaneAcentricFactor, methaneFraction, true},
        {butaneCriticalTemperature, butaneCriticalPressure, butaneAcentricFactor, 1.0 - methaneFraction, true}};
    return EquationOfState::computeFluidProperties(temperature, pressure, inputs, EquationOfState::Type::PengRobinson,
                                                   EquationOfState::MixingRules::VanDerWaals);
  };

  // butane-rich: the liquid root for both
  {
    const std::vector<EquationOfState::FluidResult> results = mixture(0.1);
    ASSERT_EQ(results.size(), 2uz);
    EXPECT_DOUBLE_EQ(results[0].compressibility, results[1].compressibility);
    EXPECT_EQ(results[0].fluidState, results[1].fluidState);
    EXPECT_EQ(results[0].fluidState, EquationOfState::FluidState::Liquid);
    EXPECT_LT(results[0].compressibility, 0.1);
  }
  // half methane: the vapor root for both
  {
    const std::vector<EquationOfState::FluidResult> results = mixture(0.5);
    ASSERT_EQ(results.size(), 2uz);
    EXPECT_DOUBLE_EQ(results[0].compressibility, results[1].compressibility);
    EXPECT_EQ(results[0].fluidState, results[1].fluidState);
    EXPECT_EQ(results[0].fluidState, EquationOfState::FluidState::Vapor);
    EXPECT_GT(results[0].compressibility, 0.5);
  }
}
