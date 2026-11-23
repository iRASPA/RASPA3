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
import system_preparation;
import mc_moves;
import randomnumbers;
import threadpool;
import archive;
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
  const std::vector<Preparation::IdealGasReference> references = Preparation::idealGasReferences(system, 0, temperatures);

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
  const std::vector<Preparation::IdealGasReference> references = Preparation::idealGasReferences(system, 0, temperatures);

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
  const std::vector<std::vector<double>> others = Preparation::prepareIdealGasReferences(systems, temperatures, other);
  const std::vector<double> all{300.0, 400.0, 350.0};
  const std::vector<Preparation::IdealGasReference> references = Preparation::idealGasReferences(systems[0], 0, all);

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
    Preparation::prepareIdealGasReferences(systems, temperatures);
    EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeight.value(), 0.11);
    EXPECT_EQ(systems[1].components[0].idealGasRosenbluthWeight.value(), 0.22);
    // ln(W_used Z_b): the normalization follows the weight the replica is run with
    std::vector<System> computed = makeLadder();
    Preparation::prepareIdealGasReferences(computed, temperatures);
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
    EXPECT_THROW(Preparation::prepareIdealGasReferences(systems, temperatures), std::runtime_error);
  }
  {
    std::vector<System> systems = makeLadder();
    for (System& system : systems) system.components[0].idealGasRosenbluthWeight = 0.11;
    EXPECT_THROW(Preparation::prepareIdealGasReferences(systems, temperatures), std::runtime_error);
    // systems that were declared one by one have a weight each
    EXPECT_NO_THROW(Preparation::prepareIdealGasReferences(systems));
  }
  {
    // at one temperature there is no configuration integral to compare; the replicas share the computed weight
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    const std::vector<double> one{300.0};
    Preparation::prepareIdealGasReferences(systems, one);
    ASSERT_TRUE(systems[0].components[0].idealGasRosenbluthWeight.has_value());
    EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeight.value(),
              systems[1].components[0].idealGasRosenbluthWeight.value());
    EXPECT_FALSE(systems[0].components[0].logInsertionNormalization.has_value());
  }
}

// A weight that was given is in the activity of a replica only where the insertions divide by it: a molecule that
// is inserted without growth samples beta f over its configuration integral (beta f when it is rigid), whatever
// the weight, and the exchange compares the replicas with that.
TEST(IDEAL_GAS_REFERENCES, given_weight_is_in_the_exchange_only_where_the_insertions_divide_by_it)
{
  const ForceField forceField = makeForceField();
  const std::vector<double> temperatures{300.0, 400.0};
  {
    auto makeLadder = [&](bool grown)
    {
      std::vector<System> systems;
      for (const double temperature : temperatures)
      {
        systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, temperature));
        systems.back().buildConformationReservoirs();
        if (!grown)
        {
          systems.back().components[0].mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 0.0);
          systems.back().components[0].mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
        }
      }
      return systems;
    };
    std::vector<System> computed = makeLadder(false);
    Preparation::prepareIdealGasReferences(computed, temperatures);
    std::vector<System> given = makeLadder(false);
    for (System& system : given) system.components[0].idealGasRosenbluthWeightLadder = {0.11, 0.22};
    Preparation::prepareIdealGasReferences(given, temperatures);
    std::vector<System> grown = makeLadder(true);
    for (System& system : grown) system.components[0].idealGasRosenbluthWeightLadder = {0.11, 0.22};
    Preparation::prepareIdealGasReferences(grown, temperatures);
    for (std::size_t i = 0; i < 2; ++i)
    {
      EXPECT_FALSE(Preparation::insertionsDivideByIdealGasRosenbluthWeight(given[i], 0));
      EXPECT_TRUE(Preparation::insertionsDivideByIdealGasRosenbluthWeight(grown[i], 0));
      // without growth: the ln of the configuration integral, as when no weight is given
      EXPECT_NEAR(given[i].components[0].logInsertionNormalization.value(),
                  computed[i].components[0].logInsertionNormalization.value(), 1e-12);
      // with growth: the weight that was given in the place of the computed one
      EXPECT_NEAR(grown[i].components[0].logInsertionNormalization.value() -
                      std::log(grown[i].components[0].idealGasRosenbluthWeight.value()),
                  computed[i].components[0].logInsertionNormalization.value() -
                      std::log(computed[i].components[0].idealGasRosenbluthWeight.value()),
                  1e-12);
    }
  }
  {
    // a rigid molecule: no normalization is prepared, and the weight divides beta f in the grown insertions only
    constexpr std::string_view rigidDimerJson =
        R"({
  "CriticalTemperature" : 305.32,
  "CriticalPressure" : 4872000.0,
  "AcentricFactor" : 0.099,
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH3", [1.54, 0.0, 0.0]]
    ]
}
)";
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "rigid-dimer", rigidDimerJson, 300.0, 2));
    systems.push_back(makeSystem(forceField, "rigid-dimer", rigidDimerJson, 400.0, 5));
    ASSERT_TRUE(systems[0].components[0].rigid);
    systems[0].components[0].idealGasRosenbluthWeight = 0.5;
    systems[1].components[0].idealGasRosenbluthWeight = 0.25;
    Preparation::prepareIdealGasReferences(systems);
    EXPECT_FALSE(systems[0].components[0].logInsertionNormalization.has_value());

    auto logBetaFugacity = [](const System& system)
    {
      const Component& component = system.components[0];
      return std::log(system.beta * component.molFraction * component.fugacityCoefficient.value_or(1.0) *
                      system.pressure);
    };
    const double withoutWeights =
        (systems[1].beta - systems[0].beta) *
            (systems[1].runningEnergies.potentialEnergy() - systems[0].runningEnergies.potentialEnergy()) +
        3.0 * (logBetaFugacity(systems[0]) - logBetaFugacity(systems[1]));

    const std::optional<double> grown = MC_Moves::ParallelTemperingLogAcceptance(systems[0], systems[1]);
    ASSERT_TRUE(grown.has_value());
    EXPECT_NEAR(grown.value(), withoutWeights - 3.0 * (std::log(0.5) - std::log(0.25)), 1e-9);

    for (System& system : systems)
    {
      system.components[0].mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 0.0);
      system.components[0].mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
    }
    const std::optional<double> placed = MC_Moves::ParallelTemperingLogAcceptance(systems[0], systems[1]);
    ASSERT_TRUE(placed.has_value());
    EXPECT_NEAR(placed.value(), withoutWeights, 1e-9);
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
    Preparation::prepareIdealGasReferences(systems, temperatures);
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

// A system on its own: the weight of a flexible molecule that is inserted with configurational bias and was given
// none is computed at the temperature of the system, to the tolerance, and a second estimate from other grows
// agrees with it.
TEST(IDEAL_GAS_REFERENCES, weight_of_a_single_system_is_computed_when_not_given)
{
  const ForceField forceField = makeForceField();
  std::vector<System> systems;
  systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
  Preparation::prepareSystems(systems);

  const Component& component = systems[0].components[0];
  ASSERT_TRUE(component.idealGasRosenbluthWeight.has_value());
  ASSERT_TRUE(component.idealGasRosenbluthWeightEstimate.has_value());
  EXPECT_EQ(component.idealGasRosenbluthWeightEstimate->temperature, 300.0);
  EXPECT_TRUE(component.idealGasRosenbluthWeightEstimate->converged);
  EXPECT_LE(component.idealGasRosenbluthWeightEstimate->relativeError, 5.0e-3);
  EXPECT_GE(component.idealGasRosenbluthWeightEstimate->numberOfGrows, 20000uz);
  // nothing is compared across temperatures
  EXPECT_FALSE(component.logInsertionNormalization.has_value());
  // said in the status of the component
  EXPECT_NE(systems[0].writeComponentStatus().find("grows of the isolated molecule"), std::string::npos);

  // other chunks are other grows: an independent estimate
  const Preparation::IdealGasReference other = Preparation::idealGasRosenbluthWeight(
      systems[0], 0, {.growsPerChunk = 1999, .minimumNumberOfGrows = 80000, .maximumNumberOfGrows = 80000});
  const double error = std::hypot(component.idealGasRosenbluthWeightEstimate->relativeError, other.relativeError);
  EXPECT_NEAR(std::log(component.idealGasRosenbluthWeight.value()), other.logRosenbluthWeight, 5.0 * error);
  EXPECT_NE(std::log(component.idealGasRosenbluthWeight.value()), other.logRosenbluthWeight);

  // prepared again, it is left as it is
  const double weight = component.idealGasRosenbluthWeight.value();
  Preparation::prepareSystems(systems);
  EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeight.value(), weight);
}

// What is not computed: a weight that was given, a molecule without conformations, and a flexible molecule that no
// move grows against its ideal gas.
TEST(IDEAL_GAS_REFERENCES, weight_is_not_computed_where_it_is_given_or_not_needed)
{
  const ForceField forceField = makeForceField();
  {
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems[0].components[0].idealGasRosenbluthWeight = 0.25;
    Preparation::prepareSystems(systems);
    EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeight.value(), 0.25);
    EXPECT_FALSE(systems[0].components[0].idealGasRosenbluthWeightEstimate.has_value());
  }
  {
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems[0].components[0].mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 0.0);
    EXPECT_FALSE(Preparation::needsIdealGasRosenbluthWeight(systems[0], 0));
    Preparation::prepareSystems(systems);
    EXPECT_FALSE(systems[0].components[0].idealGasRosenbluthWeight.has_value());

    // a Widom probe divides by the weight as an insertion does
    systems[0].components[0].mc_moves_probabilities.setProbability(Move::Types::Widom, 1.0);
    EXPECT_TRUE(Preparation::needsIdealGasRosenbluthWeight(systems[0], 0));

    // the pair and the group swap grow their molecules (the plain ones too)
    for (const Move::Types move : {Move::Types::GroupSwap, Move::Types::PairSwap})
    {
      systems[0].components[0].mc_moves_probabilities.setProbability(Move::Types::Widom, 0.0);
      EXPECT_FALSE(Preparation::needsIdealGasRosenbluthWeight(systems[0], 0));
      systems[0].components[0].mc_moves_probabilities.setProbability(move, 1.0);
      EXPECT_TRUE(Preparation::needsIdealGasRosenbluthWeight(systems[0], 0));
      systems[0].components[0].mc_moves_probabilities.setProbability(move, 0.0);
    }
  }
  {
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems[0].components[0].rigid = true;
    EXPECT_FALSE(Preparation::needsIdealGasRosenbluthWeight(systems[0], 0));
  }
}

// A computed weight belongs to its temperature: at another state point it is computed again, where a given one is
// kept.
TEST(IDEAL_GAS_REFERENCES, computed_weight_follows_the_temperature)
{
  const ForceField forceField = makeForceField();
  std::vector<System> systems;
  systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
  Preparation::prepareSystems(systems);
  const double cold = systems[0].components[0].idealGasRosenbluthWeight.value();

  systems[0].temperature = 500.0;
  systems[0].beta = 1.0 / (Units::KB * 500.0);
  Preparation::prepareStatePoint(systems[0]);
  ASSERT_TRUE(systems[0].components[0].idealGasRosenbluthWeightEstimate.has_value());
  EXPECT_EQ(systems[0].components[0].idealGasRosenbluthWeightEstimate->temperature, 500.0);
  const double hot = systems[0].components[0].idealGasRosenbluthWeight.value();
  // the torsion and the 1-4 pair are less of a hindrance at the higher temperature
  EXPECT_GT(hot, cold * 1.05);

  // the same as that of a system made at 500 K
  std::vector<System> made;
  made.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 500.0));
  Preparation::prepareSystems(made);
  EXPECT_EQ(made[0].components[0].idealGasRosenbluthWeight.value(), hot);
}

// The chunks of grows are independent of the thread that makes them: the weight is the same number on one thread
// and on four.
TEST(IDEAL_GAS_REFERENCES, weight_does_not_depend_on_the_number_of_threads)
{
  const ForceField forceField = makeForceField();
  const System system = makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0);
  const Preparation::IdealGasSettings settings{.minimumNumberOfGrows = 40000, .maximumNumberOfGrows = 40000};

  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  const ThreadPool::ThreadingType before = pool.getThreadingType();
  pool.init(1uz, ThreadPool::ThreadingType::Serial);
  const Preparation::IdealGasReference serial = Preparation::idealGasRosenbluthWeight(system, 0, settings);
  pool.init(4uz, ThreadPool::ThreadingType::ThreadPool);
  const Preparation::IdealGasReference threaded = Preparation::idealGasRosenbluthWeight(system, 0, settings);
  pool.init(1uz, before);

  EXPECT_EQ(serial.logRosenbluthWeight, threaded.logRosenbluthWeight);
  EXPECT_EQ(serial.relativeError, threaded.relativeError);
  EXPECT_EQ(serial.numberOfGrows, threaded.numberOfGrows);
}

// How a weight was computed is kept in a restart: the resumed run does not compute it again.
TEST(IDEAL_GAS_REFERENCES, computed_weight_is_kept_in_the_archive)
{
  const ForceField forceField = makeForceField();
  std::vector<System> systems;
  systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
  Preparation::prepareSystems(systems);
  const Component& component = systems[0].components[0];

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ideal-gas-weight-component.bin";
  {
    std::ofstream file(path, std::ios::binary);
    Archive<std::ofstream> archive(file);
    archive << component;
  }
  Component read;
  {
    std::ifstream file(path, std::ios::binary);
    Archive<std::ifstream> archive(file);
    archive >> read;
  }
  std::filesystem::remove(path);

  ASSERT_TRUE(read.idealGasRosenbluthWeight.has_value());
  ASSERT_TRUE(read.idealGasRosenbluthWeightEstimate.has_value());
  EXPECT_EQ(read.idealGasRosenbluthWeight.value(), component.idealGasRosenbluthWeight.value());
  EXPECT_EQ(read.idealGasRosenbluthWeightEstimate->temperature, 300.0);
  EXPECT_EQ(read.idealGasRosenbluthWeightEstimate->numberOfGrows,
            component.idealGasRosenbluthWeightEstimate->numberOfGrows);
  EXPECT_EQ(read.idealGasRosenbluthWeightEstimate->relativeError,
            component.idealGasRosenbluthWeightEstimate->relativeError);
  EXPECT_EQ(read.idealGasRosenbluthWeightEstimate->converged, component.idealGasRosenbluthWeightEstimate->converged);
}

// The tolerance of a component ('IdealGasRosenbluthWeightTolerance') comes before the default: a tighter one takes
// more grows, and one that the maximum number of grows does not reach is told.
TEST(IDEAL_GAS_REFERENCES, tolerance_of_the_component_sets_the_number_of_grows)
{
  const ForceField forceField = makeForceField();
  auto prepared = [&](std::optional<double> tolerance, const Preparation::IdealGasSettings& settings)
  {
    std::vector<System> systems;
    systems.push_back(makeSystem(forceField, "butane", molecule_fixtures::kButaneJson, 300.0));
    systems[0].components[0].idealGasRosenbluthWeightTolerance = tolerance;
    Preparation::prepareSystems(systems, settings);
    return systems[0].components[0].idealGasRosenbluthWeightEstimate.value();
  };

  const Component::IdealGasRosenbluthWeightEstimate loose = prepared(std::nullopt, {});
  const Component::IdealGasRosenbluthWeightEstimate tight = prepared(1.0e-3, {});
  EXPECT_TRUE(loose.converged);
  EXPECT_TRUE(tight.converged);
  EXPECT_LE(tight.relativeError, 1.0e-3);
  EXPECT_GT(tight.numberOfGrows, loose.numberOfGrows);

  const Component::IdealGasRosenbluthWeightEstimate capped = prepared(1.0e-6, {.maximumNumberOfGrows = 40000});
  EXPECT_FALSE(capped.converged);
  EXPECT_EQ(capped.numberOfGrows, 40000uz);
}
