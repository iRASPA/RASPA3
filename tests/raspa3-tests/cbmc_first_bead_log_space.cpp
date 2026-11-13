#include <gtest/gtest.h>

import std;

import double3;
import randomnumbers;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import cbmc;

// The first-bead Rosenbluth sums are formed in log space. An existing bead at a strongly repulsive site, well
// inside the overlap criterion, has a Boltzmann factor that is zero in double precision; with a single first-bead
// position (no other trials) the retrace summed that zero alone and took its logarithm, minus infinity, which the
// deletion move then turned into an acceptance probability of x/0.

namespace
{

// two Lennard-Jones atoms 1.66 A apart: U = 4 eps (sigma/r)^12 ~ 5e5 K, below the overlap criterion of 1e6 K,
// and -beta U ~ -1600 at 300 K, far below the -745 at which exp underflows
System makeRepulsivePair()
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  const Component x(forceField, "X", 190.0, 4.5e6, 0.01, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {},
                    {}, 5, 21, probabilities, std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {x}, {}, {2}, 5);
  system.spanOfMolecule(0, 0)[0].position = double3(10.0, 10.0, 10.0);
  system.spanOfMolecule(0, 1)[0].position = double3(11.66, 10.0, 10.0);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}

}  // namespace

TEST(CBMC_FIRST_BEAD_LOG_SPACE, a_strongly_repulsive_existing_bead_retraces_to_a_finite_log_weight)
{
  System system = makeRepulsivePair();
  const double pairEnergy = system.runningEnergies.potentialEnergy();
  ASSERT_GT(pairEnergy, 2.0e5);
  ASSERT_LT(pairEnergy, system.forceField.energyOverlapCriteria);
  const double logBoltzmann = -system.beta * pairEnergy;
  ASSERT_LT(logBoltzmann, -1000.0);  // exp(logBoltzmann) is zero in double precision

  // a single first-bead position: the existing bead is the only term of the Rosenbluth sum
  CBMC::GrowContext context = system.makeGrowContext();
  context.settings.numberOfFirstBeadPositions = 1;
  RandomNumber random(3);
  const CBMC::RetraceResult retrace =
      CBMC::retraceMolecule(random, context, system.components[0], system.spanOfMolecule(0, 0));
  EXPECT_NEAR(retrace.logRosenbluthWeight, logBoltzmann, 1e-6);

  // the reinsertion scheme with no retained partial weight: the same single term
  const CBMC::RetraceResult reinsertion =
      CBMC::retraceMolecule(random, context, system.components[0], system.spanOfMolecule(0, 0),
                            {.firstBead = CBMC::FirstBeadScheme::Reinsertion, .logStoredR = std::nullopt});
  EXPECT_NEAR(reinsertion.logRosenbluthWeight, logBoltzmann, 1e-6);
}

// The retained partial weight of the multiple-first-bead reinsertion (Esselink et al., r) is kept as its logarithm,
// summed over the other trials themselves: formed as exp(log-sum) (1 - share of the selected trial), it overflowed for
// a strongly attractive site and cancelled to zero when the selected trial dominated. The grow's weight is the
// selected trial's Boltzmann factor plus r; a single position leaves no r at all; and a retrace with an r of e^800
// (inf in linear space) stays finite.
TEST(CBMC_FIRST_BEAD_LOG_SPACE, reinsertion_retains_the_partial_weight_in_log_space)
{
  System system = makeRepulsivePair();
  CBMC::GrowContext context = system.makeGrowContext();
  RandomNumber random(11);

  context.settings.numberOfFirstBeadPositions = 8;
  const std::optional<CBMC::GrowResult> grown =
      CBMC::regrowMolecule(random, context, system.components[0], system.moleculeData[0], system.spanOfMolecule(0, 0));
  ASSERT_TRUE(grown.has_value());
  ASSERT_TRUE(grown->firstBeadLogStoredR.has_value());
  const double logSelected = -system.beta * grown->energies.potentialEnergy();
  const double logStored = grown->firstBeadLogStoredR.value();
  const double maximum = std::max(logSelected, logStored);
  const double logTotal = maximum + std::log(std::exp(logSelected - maximum) + std::exp(logStored - maximum));
  EXPECT_NEAR(logTotal - std::log(8.0), grown->logRosenbluthWeight, 1e-9 * std::max(1.0, std::abs(logTotal)));

  context.settings.numberOfFirstBeadPositions = 1;
  const std::optional<CBMC::GrowResult> single =
      CBMC::regrowMolecule(random, context, system.components[0], system.moleculeData[0], system.spanOfMolecule(0, 0));
  ASSERT_TRUE(single.has_value());
  EXPECT_FALSE(single->firstBeadLogStoredR.has_value());

  const CBMC::RetraceResult retrace =
      CBMC::retraceMolecule(random, context, system.components[0], system.spanOfMolecule(0, 0),
                            {.firstBead = CBMC::FirstBeadScheme::Reinsertion, .logStoredR = 800.0});
  EXPECT_NEAR(retrace.logRosenbluthWeight, 800.0, 1e-9);
}
