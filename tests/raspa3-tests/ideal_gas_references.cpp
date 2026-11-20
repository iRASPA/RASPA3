#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import units;
import forcefield;
import component;
import simulationbox;
import system;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;
import mc_moves;
import randomnumbers;
import running_energy;
import reweighted_histogram;

// The ideal-gas references of a flexible molecule: the Rosenbluth weight of the temperature of each replica, and the
// ratio of the intramolecular configuration integrals of two temperatures that an exchange of configurations with
// different molecule counts (and a reweighting) between them needs.

namespace
{

constexpr std::string_view kHarmonicDimerJson =
    R"({
  "CriticalTemperature" : 305.32,
  "CriticalPressure" : 4872000.0,
  "AcentricFactor" : 0.099,
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1]
  ],
  "Bonds" : [
    [["CH3", "CH3"], "HARMONIC", [96500.0, 1.54]]
  ]
}
)";

ForceField makeForceField()
{
  return ForceField({{"CH3", false, 15.0, 0.0, 0.0, 6, false}, {"CH2", false, 14.0, 0.0, 0.0, 6, false}},
                    {{108.0, 3.76}, {56.0, 3.96}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

System makeSystem(const ForceField& forceField, std::string_view name, std::string_view json, double temperature,
                  std::size_t numberOfMolecules = 0)
{
  const TemporaryFile file{std::format("ideal-gas-{}.json", name), json};
  Component component{Component::Type::Adsorbate, 0, forceField, std::string(name), file.stemPath().string(), 5, 21,
                      MCMoveProbabilities(), std::nullopt, false};
  component.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  component.mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, temperature, 1e5, 1.0, {}, {component}, {},
                {numberOfMolecules}, 5);
  system.components[0].swappable = true;
  return system;
}

// ln int w(x) exp(-u(x) / T) dx by the midpoint rule
double logIntegral(double lower, double upper, auto&& integrand)
{
  constexpr std::size_t steps = 200000;
  const double width = (upper - lower) / static_cast<double>(steps);
  double sum = 0.0;
  for (std::size_t i = 0; i < steps; ++i)
  {
    sum += integrand(lower + (static_cast<double>(i) + 0.5) * width);
  }
  return std::log(sum * width);
}

}  // namespace

// One harmonic bond: the growth draws its length from the normalized r^2 exp(-beta u(r)), so that the Rosenbluth
// weight is one at every temperature, while the configuration integral is int r^2 exp(-beta u(r)) dr.
TEST(IDEAL_GAS_REFERENCES, harmonic_dimer_has_unit_weight_and_the_integral_of_its_bond)
{
  const ForceField forceField = makeForceField();
  const System system = makeSystem(forceField, "dimer", kHarmonicDimerJson, 300.0);
  const std::vector<double> temperatures{600.0, 300.0, 400.0};
  const std::vector<MC_Moves::IdealGasReference> references = MC_Moves::idealGasReferences(system, 0, temperatures);

  auto exact = [](double T)
  {
    return logIntegral(0.5, 3.0, [&](double r)
                       { return r * r * std::exp(-0.5 * 96500.0 * (r - 1.54) * (r - 1.54) / T); });
  };
  ASSERT_EQ(references.size(), 3uz);
  for (std::size_t i = 0; i < 3; ++i)
  {
    EXPECT_EQ(references[i].temperature, temperatures[i]);
    EXPECT_NEAR(references[i].logRosenbluthWeight, 0.0, 1e-6);
    // relative to the lowest temperature
    EXPECT_NEAR(references[i].logPartitionFunction, exact(temperatures[i]) - exact(300.0), 0.02);
    EXPECT_LT(references[i].logPartitionFunctionError, 0.05);
  }
  EXPECT_EQ(references[1].logPartitionFunction, 0.0);
}

// Propane with fixed bonds: one harmonic bend, int sin(theta) exp(-beta u(theta)) dtheta.
TEST(IDEAL_GAS_REFERENCES, propane_has_the_integral_of_its_bend)
{
  const ForceField forceField = makeForceField();
  const System system = makeSystem(forceField, "propane", molecule_fixtures::kPropaneJson, 300.0);
  const std::vector<double> temperatures{300.0, 500.0};
  const std::vector<MC_Moves::IdealGasReference> references = MC_Moves::idealGasReferences(system, 0, temperatures);

  auto exact = [](double T)
  {
    const double theta0 = 114.0 * std::numbers::pi / 180.0;
    return logIntegral(0.0, std::numbers::pi, [&](double theta)
                       { return std::sin(theta) * std::exp(-0.5 * 62500.0 * (theta - theta0) * (theta - theta0) / T); });
  };
  EXPECT_NEAR(references[1].logPartitionFunction - references[0].logPartitionFunction, exact(500.0) - exact(300.0),
              0.02);
}

// The replicas of a ladder: each gets the weight of its temperature and the normalization of its insertions, and
// the exchange of configurations with different molecule counts has the ratio of the two in its acceptance.
TEST(IDEAL_GAS_REFERENCES, replicas_of_a_ladder_exchange_a_flexible_molecule)
{
  const ForceField forceField = makeForceField();
  const std::vector<double> temperatures{300.0, 400.0};
  std::vector<System> systems;
  systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0, 2));
  systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 400.0, 5));
  for (System& system : systems) system.buildConformationReservoirs();

  // not prepared: the ratio of the configuration integrals is not known, and the pair is refused
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(systems[0], false).has_value());
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(systems[0], true).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(systems[0], systems[1]).has_value());

  const std::vector<double> other{350.0};
  const std::vector<std::vector<double>> others = MC_Moves::prepareIdealGasReferences(systems, temperatures, other);
  const std::vector<double> all{300.0, 400.0, 350.0};
  const std::vector<MC_Moves::IdealGasReference> references = MC_Moves::idealGasReferences(systems[0], 0, all);

  for (std::size_t i = 0; i < 2; ++i)
  {
    const Component& component = systems[i].components[0];
    ASSERT_TRUE(component.idealGasRosenbluthWeight.has_value());
    ASSERT_TRUE(component.logInsertionNormalization.has_value());
    // no weight was given: the computed one, and the normalization is the ln of the configuration integral
    EXPECT_NEAR(std::log(component.idealGasRosenbluthWeight.value()), references[i].logRosenbluthWeight, 1e-12);
    EXPECT_NEAR(component.logInsertionNormalization.value(), references[i].logPartitionFunction, 1e-12);
  }
  EXPECT_NEAR(others[0][0], references[2].logPartitionFunction, 1e-12);
  // between the two, as the integral grows with the temperature
  EXPECT_GT(others[0][0], references[0].logPartitionFunction);
  EXPECT_LT(others[0][0], references[1].logPartitionFunction);

  const std::optional<double> logAcceptance = MC_Moves::ParallelTemperingLogAcceptance(systems[0], systems[1]);
  ASSERT_TRUE(logAcceptance.has_value());
  auto logActivity = [](const System& system)
  {
    const Component& component = system.components[0];
    return std::log(system.beta * component.molFraction * component.fugacityCoefficient.value_or(1.0) *
                    system.pressure) -
           component.logInsertionNormalization.value();
  };
  const double expected = (systems[1].beta - systems[0].beta) * (systems[1].runningEnergies.potentialEnergy() -
                                                                 systems[0].runningEnergies.potentialEnergy()) +
                          3.0 * (logActivity(systems[0]) - logActivity(systems[1]));
  EXPECT_NEAR(logAcceptance.value(), expected, 1e-9);
}

// A list of weights is one per temperature of the ladder; one number for a flexible molecule on several
// temperatures is the weight of one of them.
TEST(IDEAL_GAS_REFERENCES, weights_of_a_ladder_are_given_per_temperature)
{
  const ForceField forceField = makeForceField();
  const std::vector<double> temperatures{300.0, 400.0};
  auto makeLadder = [&]()
  {
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 400.0));
    for (System& system : systems) system.buildConformationReservoirs();
    return systems;
  };
  {
    std::vector<System> systems = makeLadder();
    for (System& system : systems) system.components[0].idealGasRosenbluthWeightLadder = {0.11, 0.22};
    MC_Moves::prepareIdealGasReferences(systems, temperatures);
    EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeight.value(), 0.11);
    EXPECT_EQ(systems[1].components[0].idealGasRosenbluthWeight.value(), 0.22);
    // ln(W_used Z_b): the normalization follows the weight the replica is run with
    std::vector<System> computed = makeLadder();
    MC_Moves::prepareIdealGasReferences(computed, temperatures);
    for (std::size_t i = 0; i < 2; ++i)
    {
      EXPECT_NEAR(systems[i].components[0].logInsertionNormalization.value() -
                      std::log(systems[i].components[0].idealGasRosenbluthWeight.value()),
                  computed[i].components[0].logInsertionNormalization.value() -
                      std::log(computed[i].components[0].idealGasRosenbluthWeight.value()),
                  1e-12);
    }
  }
  {
    std::vector<System> systems = makeLadder();
    for (System& system : systems) system.components[0].idealGasRosenbluthWeightLadder = {0.11, 0.22, 0.33};
    EXPECT_THROW(MC_Moves::prepareIdealGasReferences(systems, temperatures), std::runtime_error);
  }
  {
    std::vector<System> systems = makeLadder();
    for (System& system : systems) system.components[0].idealGasRosenbluthWeight = 0.11;
    EXPECT_THROW(MC_Moves::prepareIdealGasReferences(systems, temperatures), std::runtime_error);
    // systems that were declared one by one have a weight each
    EXPECT_NO_THROW(MC_Moves::prepareIdealGasReferences(systems));
  }
  {
    // at one temperature nothing is computed
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    const std::vector<double> one{300.0};
    MC_Moves::prepareIdealGasReferences(systems, one);
    EXPECT_FALSE(systems[0].components[0].idealGasRosenbluthWeight.has_value());
    EXPECT_FALSE(systems[0].components[0].logInsertionNormalization.has_value());
  }
}

// The exchange itself: two replicas of harmonic dimers at 300 and 600 K, each sampling its own grand-canonical
// distribution with insertions and deletions. Exchanges of their configurations leave both distributions as they
// are when the acceptance has the ratio of the configuration integrals, here sqrt(T_B / T_A) per molecule of
// difference to a good approximation; without it the mean loading of the two replicas is pulled apart.
TEST(IDEAL_GAS_REFERENCES, exchanges_between_temperatures_leave_the_loadings_as_they_are)
{
  const ForceField forceField = makeForceField();
  const std::vector<double> temperatures{300.0, 600.0};

  enum class Exchange { None, Prepared, WithoutTheIntegrals };
  auto meanLoadings = [&](Exchange exchange) -> std::array<double, 2>
  {
    std::vector<System> systems;
    for (const double temperature : temperatures)
    {
      systems.push_back(makeSystem(forceField, "dimer", kHarmonicDimerJson, temperature));
      systems.back().buildConformationReservoirs();
      systems.back().precomputeTotalRigidEnergy();
      systems.back().runningEnergies = systems.back().computeTotalEnergies();
    }
    MC_Moves::prepareIdealGasReferences(systems, temperatures);
    if (exchange == Exchange::WithoutTheIntegrals)
    {
      // the weights alone (both are one): what the acceptance had before
      for (System& system : systems) system.components[0].logInsertionNormalization = 0.0;
    }

    RandomNumber random(11);
    std::size_t fractionalMoleculeSystem = 0;
    std::array<double, 2> sum{0.0, 0.0};
    constexpr std::size_t numberOfSteps = 400000;
    for (std::size_t step = 0; step < numberOfSteps; ++step)
    {
      for (std::size_t i = 0; i < 2; ++i)
      {
        MC_Moves::performRandomMoveProduction(random, systems[i], systems[i], 0, fractionalMoleculeSystem, 0);
      }
      if (exchange != Exchange::None && step % 2 == 0)
      {
        (void)MC_Moves::ParallelTemperingSwap(random, systems[0], systems[1]);
      }
      for (std::size_t i = 0; i < 2; ++i)
      {
        sum[i] += static_cast<double>(systems[i].numberOfIntegerMoleculesPerComponent[0]);
      }
    }
    return {sum[0] / static_cast<double>(numberOfSteps), sum[1] / static_cast<double>(numberOfSteps)};
  };

  const std::array<double, 2> reference = meanLoadings(Exchange::None);
  const std::array<double, 2> prepared = meanLoadings(Exchange::Prepared);
  const std::array<double, 2> without = meanLoadings(Exchange::WithoutTheIntegrals);
  EXPECT_NEAR(prepared[0], reference[0], 0.02);
  EXPECT_NEAR(prepared[1], reference[1], 0.02);
  EXPECT_GT(std::abs(without[0] - reference[0]) + std::abs(without[1] - reference[1]), 0.05);
}

// The reweighting driver takes a flexible adsorbate across temperatures (it was refused): the replicas carry the
// normalization of their temperature, which the reduced potentials and the configuration swaps use.
TEST(IDEAL_GAS_REFERENCES, reweighted_histogram_takes_a_flexible_adsorbate_across_temperatures)
{
  const ForceField forceField = makeForceField();
  System system = makeSystem(forceField, "dimer", kHarmonicDimerJson, 300.0);
  system.tmmc.maxMacrostate = 20;
  ReweightedHistogramParameters parameters;
  parameters.numberOfInitializationCycles = 0;
  parameters.numberOfProductionCycles = 10;
  parameters.writeBinaryRestartEvery = 0;
  const ReweightedHistogram driver(system, {300.0, 600.0}, {1e4}, parameters);
  ASSERT_EQ(driver.systems.size(), 2uz);
  const std::optional<double> cold = driver.systems[0].components[0].logInsertionNormalization;
  const std::optional<double> hot = driver.systems[1].components[0].logInsertionNormalization;
  ASSERT_TRUE(cold.has_value() && hot.has_value());
  // one harmonic bond: the configuration integral goes as sqrt(T) to a good approximation
  EXPECT_NEAR(hot.value() - cold.value(), 0.5 * std::log(2.0), 0.03);

  // one number for the two temperatures is the weight of one of them
  system.components[0].idealGasRosenbluthWeight = 0.5;
  EXPECT_THROW(ReweightedHistogram(system, {300.0, 600.0}, {1e4}, parameters), std::runtime_error);
  system.components[0].idealGasRosenbluthWeightLadder = {0.5, 0.25};
  EXPECT_NO_THROW(ReweightedHistogram(system, {300.0, 600.0}, {1e4}, parameters));
}
