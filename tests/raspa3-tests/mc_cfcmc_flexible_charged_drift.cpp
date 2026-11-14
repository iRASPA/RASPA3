#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import atom;
import forcefield;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import randomnumbers;
import mc_moves;
import reaction;
import mc_moves_lambda_exchange;

// CFCMC moves rescale the charges of fractional molecules in place. For a flexible molecule with charged
// non-excluded intramolecular pairs (the 1-4 pair of butane here) the Ewald Fourier sum counts such a pair with
// the scaled charges, and the intramolecular Coulomb term carries the compensating, equally scaled, kernel: the
// two change by opposite amounts when the scaling changes, and a move that updates only the Fourier part drifts
// from a full recomputation in both terms.

namespace
{

ForceField makeChargedAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.25, 0.0, 6, false}, {"CH2", false, 14.03, -0.25, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    false, true);
}

Component makeChargedButane(const ForceField& forceField, const MCMoveProbabilities& probabilities)
{
  TemporaryFile file("butane.json", molecule_fixtures::kButaneJson);
  Component butane(Component::Type::Adsorbate, 0, forceField, "butane", file.stemPath().string(), 5, 21, probabilities,
                   std::nullopt, false);
  butane.lnPartitionFunction = 1.0;
  return butane;
}

void runAndCheckDrift(const MCMoveProbabilities& probabilities)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  std::vector<System> systems{System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {},
                                     {makeChargedButane(forceField, probabilities)}, {}, {6}, 5)};
  ASSERT_TRUE(systems[0].forceField.usesEwaldFourier());
  ASSERT_NE(systems[0].computeTotalEnergies().intraCoul, 0.0);  // the charged 1-4 pairs are exercised

  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 17uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    const RunningEnergy recomputed = s.computeTotalEnergies();
    const RunningEnergy drift = s.runningEnergies - recomputed;
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

}  // namespace

// The reference geometry of a flexible component is a topology with placeholder coordinates (the butane fixture
// has every bead at the origin). The ideal-gas conformation must be grown from it, never retraced: a retrace of
// coincident charged beads has no weight, and a chain started on it never leaves it, so every inserted molecule
// would be the collapsed template.
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, ideal_gas_conformation_is_grown_from_the_placeholder_geometry)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {},
                {makeChargedButane(forceField, MCMoveProbabilities())}, {}, {0}, 5);
  RandomNumber random(5);
  for (int draw = 0; draw < 5; ++draw)
  {
    const auto [molecule, atoms] = system.equilibratedIdealGasMoleculeRandomInBox(random, 0);
    ASSERT_EQ(atoms.size(), 4uz);
    for (std::size_t bond = 0; bond + 1 < atoms.size(); ++bond)
    {
      EXPECT_NEAR((atoms[bond + 1].position - atoms[bond].position).length(), 1.54, 1e-3)
          << "draw " << draw << " bond " << bond;
    }
    EXPECT_GT((atoms[3].position - atoms[0].position).length(), 2.0) << "draw " << draw;
  }
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, swap_cfcmc)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, swap_cbcfcmc)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);
  runAndCheckDrift(probabilities);
}

// The drift of a missing kernel term depends on the end states only (each hop's error is an exact difference),
// and the kernel only switches on above lambda = 1/2: the drift is checked after every step, not at the end of
// a run that may happen to finish where it started.
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, pair_swap_cfcmc_on_two_charged_butanes)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::PairSwapCFCMC, 1.0);
  Component a = makeChargedButane(forceField, probabilities);
  Component b = makeChargedButane(forceField, probabilities);
  b.name = "butane-b";
  a.pairComponentId = 1;
  a.maximumPairDistance = 6.0;
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b}, {}, {4, 4}, 5);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(5);
  std::size_t fractionalMoleculeSystem = 0;
  double largestDrift = 0.0;
  bool coulombSwitchedOn = false;
  for (std::size_t step = 0; step < 400; ++step)
  {
    MC_Moves::performRandomMoveProduction(random, system, system, 0, fractionalMoleculeSystem, 0);
    coulombSwitchedOn = coulombSwitchedOn || system.components[0].lambdaPairSwap.currentBin > 10;
    largestDrift =
        std::max(largestDrift, std::abs((system.runningEnergies - system.computeTotalEnergies()).potentialEnergy()));
  }
  ASSERT_TRUE(coulombSwitchedOn);
  EXPECT_LT(largestDrift, 1e-6);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, group_swap_cfcmc)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GroupSwapCFCMC, 1.0);
  Component a = makeChargedButane(forceField, probabilities);
  Component b = makeChargedButane(forceField, MCMoveProbabilities{});
  b.name = "butane-b";
  a.groupComponentIds = {1};
  a.maximumGroupDistance = 6.0;
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b}, {}, {4, 4}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 23uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_swap_cfcmc_between_two_boxes)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  const Component butane = makeChargedButane(forceField, probabilities);
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 29uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

namespace
{
void runGibbsDrift(Move::Types move, std::size_t seed)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(move, 1.0);
  const Component butane = makeChargedButane(forceField, probabilities);
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, seed, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}
}  // namespace

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_swap_cbcfcmc_between_two_boxes)
{
  runGibbsDrift(Move::Types::GibbsSwapCBCFCMC, 31);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_conventional_cfcmc_between_two_boxes)
{
  runGibbsDrift(Move::Types::GibbsConventionalCFCMC, 37);
}

namespace
{
// two components of charged butane, the first driving a pair or group move onto the second
void runTwoComponentDrift(Move::Types move, std::size_t seed)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(move, 1.0);
  const bool pair = move == Move::Types::PairSwapCBCFCMC;
  Component a = makeChargedButane(forceField, probabilities);
  Component b = makeChargedButane(forceField, pair ? probabilities : MCMoveProbabilities{});
  b.name = "butane-b";
  if (pair)
  {
    a.pairComponentId = 1;
    a.maximumPairDistance = 6.0;
  }
  else
  {
    a.groupComponentIds = {1};
    a.maximumGroupDistance = 6.0;
  }
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b}, {}, {4, 4}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, seed, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}
}  // namespace

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, pair_swap_cbcfcmc) { runTwoComponentDrift(Move::Types::PairSwapCBCFCMC, 41); }

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, group_swap_cbcfcmc) { runTwoComponentDrift(Move::Types::GroupSwapCBCFCMC, 43); }

// a serial Rx/CFC reaction between two charged butane components: the lambda changes rescale the fractional
// molecules of the active side in place
// The lambda moves of the serial and the parallel reactions rescale the fractional molecules in place, and the
// boundary crossings of the parallel reactions turn molecules whole or fractional in place: the intramolecular
// Coulomb kernel of a charged flexible molecule changes with its Fourier term in every one of them.
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, serial_and_parallel_reactions)
{
  for (Move::Types move : {Move::Types::ReactionCFCMC, Move::Types::ReactionCBCFCMC,
                           Move::Types::ReactionConventionalCFCMC, Move::Types::ReactionConventionalCBCFCMC})
  {
    const ForceField forceField = makeChargedAlkaneForceField();
    MCMoveProbabilities componentProbabilities;
    componentProbabilities.setProbability(Move::Types::Translation, 1.0);
    MCMoveProbabilities systemProbabilities;
    systemProbabilities.setProbability(move, 1.0);
    Component a = makeChargedButane(forceField, componentProbabilities);
    Component b = makeChargedButane(forceField, componentProbabilities);
    b.name = "butane-b";
    std::vector<System> systems{System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b},
                                       {}, {5, 5}, 5, systemProbabilities)};
    Reaction reaction(0, {1, 0}, {0, 1});
    reaction.reactionMove = move;
    systems[0].reactions.list.push_back(reaction);
    systems[0].createReactionFractionalMolecules();

    MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 47uz, 5, false);
    mc.run();
    for (System& s : mc.systems)
    {
      const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
      EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6) << Move::moveNames[std::to_underlying(move)];
      EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6) << Move::moveNames[std::to_underlying(move)];
    }
  }
}

// a lambda exchange between two replicas with the GC lambda pinned at different bins: positions stay, the
// fractional molecules are rescaled in place
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, lambda_exchange)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  auto makeReplica = [&](std::size_t bin)
  {
    Component butane = makeChargedButane(forceField, probabilities);
    butane.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::GC;
    butane.fixedLambdaBin = bin;
    System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5);
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };
  System systemA = makeReplica(4);
  System systemB = makeReplica(16);

  // one accepted exchange (the move updates the running energies itself; an even number of exchanges would
  // bring both replicas back to where they started)
  RandomNumber random(53);
  bool exchanged = false;
  for (std::size_t attempt = 0; attempt < 200 && !exchanged; ++attempt)
  {
    exchanged = MC_Moves::LambdaExchange(random, systemA, systemB).has_value();
  }
  ASSERT_TRUE(exchanged);
  for (System* s : {&systemA, &systemB})
  {
    const RunningEnergy drift = s->runningEnergies - s->computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

namespace
{
// butane with a net charge of +0.5: the net-charge Ewald correction then depends on the fractional scaling
ForceField makeNetChargedAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.5, 0.0, 6, false}, {"CH2", false, 14.03, -0.25, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    false, true);
}

double recomputedNetCharge(const System& system)
{
  double charge = system.netChargeFramework;
  for (const Atom& atom : system.spanOfMoleculeAtoms()) charge += atom.scalingCoulomb * atom.charge;
  return charge;
}
}  // namespace

// The fractional molecule of a conventional Gibbs CFCMC transfer is created switched off and then set to its lambda
// in place: the effective pseudo-atom counts of the tail corrections and the net charge follow its scaling (both stayed
// at their values of lambda zero, and the first energies of the run used them). The per-group charge of the tagged
// atoms outside a move is what the Ewald difference needs for dU/dlambda.
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_conventional_fractional_molecule_sets_net_charge_and_tail_counts)
{
  const ForceField forceField = makeNetChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsConventionalCFCMC, 1.0);
  Component butane = makeChargedButane(forceField, probabilities);
  // the dU/dlambda group the input reader assigns
  butane.lambdaGibbs.dUdlambdaGroupId = 1;
  butane.lambdaGibbs.computeDUdlambda = true;
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5);
  ASSERT_TRUE(system.usesGibbsConventionalCFCMC());

  // the fractional molecule as created, switched off (with consistent counts), then set to its lambda
  const std::size_t slot = system.indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(0);
  for (Atom& atom : system.spanOfMolecule(0, slot)) atom.setScalingToFractional(0.0, 1);
  system.computeTailCorrectionCounts();
  const std::vector<double> countsOff = system.effectiveNumberOfPseudoAtomsVDW;
  system.initializeGibbsConventionalFractionalMolecules(false);

  System recomputed = system;
  recomputed.computeTailCorrectionCounts();
  // the Van der Waals coupling is on at the starting lambda (the Coulomb coupling of lambda 1/2 is still zero)
  ASSERT_GT(std::abs(recomputed.effectiveNumberOfPseudoAtomsVDW[0] - countsOff[0]), 1e-6);
  EXPECT_NEAR(system.netCharge, recomputedNetCharge(system), 1e-12);
  ASSERT_EQ(system.effectiveNumberOfPseudoAtomsVDW.size(), recomputed.effectiveNumberOfPseudoAtomsVDW.size());
  for (std::size_t type = 0; type < system.effectiveNumberOfPseudoAtomsVDW.size(); ++type)
  {
    EXPECT_NEAR(system.effectiveNumberOfPseudoAtomsVDW[type], recomputed.effectiveNumberOfPseudoAtomsVDW[type], 1e-12);
  }

  // the tagged charge of the fractional molecule (unscaled), and none outside it
  const std::span<const Atom> fractionalMolecule = system.spanOfMolecule(0, slot);
  const std::uint8_t groupId = fractionalMolecule.front().groupId;
  ASSERT_GT(groupId, 0);
  EXPECT_NEAR(system.ewaldMoveChargeDerivative({})[groupId - 1], 0.5, 1e-12);
  EXPECT_NEAR(system.ewaldMoveChargeDerivative(fractionalMolecule)[groupId - 1], 0.0, 1e-12);
}

// Two charged ion components, each with its own CFCMC fractional molecule and dU/dlambda group: a swap of one
// changes the net charge while the other's fractional ion carries group-tagged charge outside the move. The Ewald
// net-charge correction then contributes 2 u_ion dQ Q'_g to the running dU/dlambda of the other group, which the
// swap sub-steps left out (no external group charge was passed), so the running dU/dlambda drifted from a full
// recomputation at every accepted swap.
TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, swap_cfcmc_keeps_the_dudlambda_of_another_fractional_ion)
{
  const ForceField forceField({{"Na", false, 22.99, 0.6, 0.0, 11, false}, {"Cl", false, 35.45, -0.4, 0.0, 17, false}},
                              {{15.0, 2.6}, {40.0, 3.4}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                              true, false, true);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  Component sodium(forceField, "sodium", 0.0, 0.0, 0.0, {Atom({0, 0, 0}, 0.6, 1.0, 0, 0, 0, false, false)}, {}, {}, 5,
                   21, probabilities, std::nullopt, false);
  Component chloride(forceField, "chloride", 0.0, 0.0, 0.0, {Atom({0, 0, 0}, -0.4, 1.0, 0, 1, 1, false, false)}, {}, {},
                     5, 21, probabilities, std::nullopt, false);
  sodium.lnPartitionFunction = 1.0;
  chloride.lnPartitionFunction = 1.0;
  // thermodynamic integration on both: each fractional ion carries its own dU/dlambda group
  sodium.lambdaGC.computeDUdlambda = true;
  chloride.lambdaGC.computeDUdlambda = true;
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {sodium, chloride}, {}, {6, 6},
                5);
  ASSERT_TRUE(system.forceField.usesEwaldFourier());
  ASSERT_TRUE(system.components[0].hasFractionalMolecule);
  ASSERT_TRUE(system.components[1].hasFractionalMolecule);
  ASSERT_NE(system.components[0].lambdaGC.dUdlambdaGroupId, 0);
  ASSERT_NE(system.components[1].lambdaGC.dUdlambdaGroupId, 0);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(11);
  std::size_t fractionalMoleculeSystem = 0;
  double largestDrift = 0.0;
  std::size_t numberOfExchanges = 0;
  for (std::size_t step = 0; step < 3000; ++step)
  {
    const std::size_t before = system.numberOfIntegerMoleculesPerComponent[step % 2];
    MC_Moves::performRandomMoveProduction(random, system, system, step % 2, fractionalMoleculeSystem, 0);
    if (system.numberOfIntegerMoleculesPerComponent[step % 2] != before) ++numberOfExchanges;
    const RunningEnergy recomputed = system.computeTotalEnergies();
    for (std::size_t g = 0; g < 2; ++g)
      largestDrift =
          std::max(largestDrift, std::abs(system.runningEnergies.dudlambdaEwald[g] - recomputed.dudlambdaEwald[g]));
  }
  ASSERT_GT(numberOfExchanges, 0uz);
  EXPECT_LT(largestDrift, 1e-6);
}
