#include <gtest/gtest.h>

import std;

import atom;
import component;
import double3;
import forcefield;
import mc_moves;
import mc_moves_move_types;
import mc_moves_pair_deletion_cbmc;
import mc_moves_pair_insertion_cbmc;
import mc_moves_probabilities;
import mc_moves_random_rotation;
import mc_moves_statistics;
import move_statistics;
import molecule;
import randomnumbers;
import simulationbox;
import system;
import transition_matrix;

namespace
{

ForceField makeNonInteractingForceField()
{
  return ForceField({{"A", false, 1.0, 0.0, 0.0, 1, false}, {"B", false, 1.0, 0.0, 0.0, 1, false}},
                    {{0.0, 1.0}, {0.0, 1.0}}, ForceField::MixingRule::Lorentz_Berthelot, 4.0, 4.0, 4.0, false, false,
                    false);
}

Component makeIonComponent(const ForceField& forceField, std::size_t componentId, std::string_view name,
                           std::size_t type, const MCMoveProbabilities& probabilities)
{
  Component component = Component::makeIon(forceField, componentId, name, type, 0.0);
  component.mc_moves_probabilities = probabilities;
  component.fugacityCoefficient = 1.0;
  component.idealGasRosenbluthWeight = 1.0;
  return component;
}

System makePairDeletionSystem(bool twoNeighbors, double pressure = 1.0e12)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::PairSwap, 1.0);

  Component componentA = makeIonComponent(forceField, 0, "A", 0, probabilities);
  Component componentB = makeIonComponent(forceField, 1, "B", 1, {});
  componentA.pairComponentId = 1;
  componentB.pairComponentId = 0;
  componentA.maximumPairDistance = 1.5;
  componentB.maximumPairDistance = 1.5;

  const std::vector<std::vector<double3>> positions{
      {double3(0.5, 5.0, 5.0)},
      {double3(9.5, 5.0, 5.0), twoNeighbors ? double3(1.5, 5.0, 5.0) : double3(5.0, 5.0, 5.0), double3(5.0, 8.0, 5.0)}};
  return System(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, pressure, 1.0, {}, {componentA, componentB},
                positions, {0, 0}, 5);
}

System makePairInsertionSystem(const std::vector<double3>& positionsB, double pressure)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::PairSwap, 1.0);

  Component componentA = makeIonComponent(forceField, 0, "A", 0, probabilities);
  Component componentB = makeIonComponent(forceField, 1, "B", 1, {});
  componentA.pairComponentId = 1;
  componentB.pairComponentId = 0;
  componentA.maximumPairDistance = 2.0;
  componentB.maximumPairDistance = 2.0;

  return System(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, pressure, 1.0, {}, {componentA, componentB},
                {{}, positionsB}, {0, 0}, 5);
}

double totalTrials(const System& system, Move::Types move)
{
  const auto& statistics = std::get<MoveStatistics<double3>>(system.components[0].mc_moves_statistics[move]);
  return statistics.totalCounts.x + statistics.totalCounts.y + statistics.totalCounts.z;
}

System makeInitializationRoutingSystem(Move::Types move)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilitiesA;
  probabilitiesA.setProbability(move, 1.0);

  Component componentA = makeIonComponent(forceField, 0, "A", 0, probabilitiesA);
  if (move == Move::Types::PairSwapCFCMC || move == Move::Types::PairSwapCBCFCMC)
  {
    Component componentB = makeIonComponent(forceField, 1, "B", 1, {});
    componentA.pairComponentId = 1;
    componentB.pairComponentId = 0;
    componentA.maximumPairDistance = 2.0;
    componentB.maximumPairDistance = 2.0;
    return System(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, 1.0e5, 1.0, {}, {componentA, componentB},
                  {}, {1, 1}, 5);
  }

  return System(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, 1.0e5, 1.0, {}, {componentA}, {}, {1}, 5);
}

}  // namespace

TEST(MC_COMPONENT_MOVES, pair_deletion_uses_all_minimum_image_neighbors)
{
  const auto checkNeighborFactor = [](auto deletionMove)
  {
    System oneNeighbor = makePairDeletionSystem(false);
    System twoNeighbors = makePairDeletionSystem(true);
    RandomNumber randomOne(91);
    RandomNumber randomTwo(91);

    const auto [energyOne, acceptanceOne] = deletionMove(randomOne, oneNeighbor);
    const auto [energyTwo, acceptanceTwo] = deletionMove(randomTwo, twoNeighbors);

    EXPECT_FALSE(energyOne.has_value());
    EXPECT_FALSE(energyTwo.has_value());
    ASSERT_GT(acceptanceOne.x, 0.0);
    EXPECT_NEAR(acceptanceTwo.x / acceptanceOne.x, 2.0, 1.0e-12);
  };

  checkNeighborFactor([](RandomNumber& random, System& system)
                      { return MC_Moves::pairDeletionMove(random, system, 0, 0); });
  checkNeighborFactor([](RandomNumber& random, System& system)
                      { return MC_Moves::pairDeletionMoveCBMC(random, system, 0, 0); });
}

TEST(MC_COMPONENT_MOVES, pair_deletion_selects_each_neighbor_across_seeds)
{
  const auto checkUniformChoices = [](auto deletionMove)
  {
    std::size_t selectedBoundaryNeighbor = 0;
    std::size_t selectedInteriorNeighbor = 0;

    for (std::size_t seed = 0; seed < 64; ++seed)
    {
      System system = makePairDeletionSystem(true, 1.0e-12);
      RandomNumber random(seed);
      const auto [energy, acceptance] = deletionMove(random, system);

      ASSERT_TRUE(energy.has_value());
      ASSERT_GT(acceptance.x, 1.0);
      ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 2uz);

      bool boundaryNeighborRemains = false;
      bool interiorNeighborRemains = false;
      for (std::size_t moleculeB = 0; moleculeB < system.numberOfMoleculesPerComponent[1]; ++moleculeB)
      {
        const double x = system.spanOfMolecule(1, moleculeB).front().position.x;
        boundaryNeighborRemains = boundaryNeighborRemains || std::abs(x - 9.5) < 1.0e-12;
        interiorNeighborRemains = interiorNeighborRemains || std::abs(x - 1.5) < 1.0e-12;
      }
      selectedBoundaryNeighbor += !boundaryNeighborRemains;
      selectedInteriorNeighbor += !interiorNeighborRemains;
    }

    EXPECT_EQ(selectedBoundaryNeighbor + selectedInteriorNeighbor, 64uz);
    EXPECT_GE(selectedBoundaryNeighbor, 20uz);
    EXPECT_LE(selectedBoundaryNeighbor, 44uz);
    EXPECT_GE(selectedInteriorNeighbor, 20uz);
    EXPECT_LE(selectedInteriorNeighbor, 44uz);
  };

  checkUniformChoices([](RandomNumber& random, System& system)
                      { return MC_Moves::pairDeletionMove(random, system, 0, 0); });
  checkUniformChoices([](RandomNumber& random, System& system)
                      { return MC_Moves::pairDeletionMoveCBMC(random, system, 0, 0); });
}

TEST(MC_COMPONENT_MOVES, pair_insertion_uses_reverse_state_neighbor_count)
{
  const auto checkReverseNeighborCount = [](auto insertionMove)
  {
    constexpr std::size_t seed = 918;
    System probe = makePairInsertionSystem({double3(4.0, 4.0, 4.0), double3(6.0, 6.0, 6.0)}, 1.0e12);
    RandomNumber probeRandom(seed);
    const auto [probeEnergy, probeAcceptance] = insertionMove(probeRandom, probe);
    ASSERT_TRUE(probeEnergy.has_value());
    const double3 trialPositionA = probe.spanOfMolecule(0, 0).front().position;

    System twoReversePartners = makePairInsertionSystem(
        {trialPositionA + double3(9.0, 0.0, 0.0), trialPositionA + double3(5.0, 0.0, 0.0)}, 1.0e-12);
    System threeReversePartners = makePairInsertionSystem(
        {trialPositionA + double3(9.0, 0.0, 0.0), trialPositionA + double3(1.0, 0.0, 0.0)}, 1.0e-12);
    RandomNumber randomTwo(seed);
    RandomNumber randomThree(seed);
    const auto [energyTwo, acceptanceTwo] = insertionMove(randomTwo, twoReversePartners);
    const auto [energyThree, acceptanceThree] = insertionMove(randomThree, threeReversePartners);

    EXPECT_FALSE(energyTwo.has_value());
    EXPECT_FALSE(energyThree.has_value());
    ASSERT_GT(acceptanceTwo.z, 0.0);
    EXPECT_NEAR(acceptanceThree.z / acceptanceTwo.z, 2.0 / 3.0, 1.0e-12);
  };

  checkReverseNeighborCount([](RandomNumber& random, System& system)
                            { return MC_Moves::pairInsertionMove(random, system, 0); });
  checkReverseNeighborCount([](RandomNumber& random, System& system)
                            { return MC_Moves::pairInsertionMoveCBMC(random, system, 0); });
}

TEST(MC_COMPONENT_MOVES, pair_insertion_deletion_proposals_are_reciprocal)
{
  constexpr std::size_t seed = 918;
  System probe = makePairInsertionSystem({double3(4.0, 4.0, 4.0), double3(6.0, 6.0, 6.0)}, 1.0e12);
  RandomNumber probeRandom(seed);
  ASSERT_TRUE(MC_Moves::pairInsertionMove(probeRandom, probe, 0).first.has_value());
  const double3 trialPositionA = probe.spanOfMolecule(0, 0).front().position;

  System system = makePairInsertionSystem(
      {trialPositionA + double3(9.0, 0.0, 0.0), trialPositionA + double3(1.0, 0.0, 0.0)}, 1.0e12);
  RandomNumber insertionRandom(seed);
  const auto [insertionEnergy, insertionAcceptance] = MC_Moves::pairInsertionMove(insertionRandom, system, 0);
  ASSERT_TRUE(insertionEnergy.has_value());
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 1uz);

  RandomNumber deletionRandom(41);
  const auto [deletionEnergy, deletionAcceptance] = MC_Moves::pairDeletionMove(deletionRandom, system, 0, 0);
  EXPECT_FALSE(deletionEnergy.has_value());
  EXPECT_NEAR(insertionAcceptance.z * deletionAcceptance.x, 1.0, 1.0e-11);
}

TEST(MC_COMPONENT_MOVES, random_rotation_uses_uniform_rotation_matrix)
{
  const ForceField forceField = makeNonInteractingForceField();
  Component component(
      forceField, "rigid", 1.0, 1.0, 0.0,
      {Atom({1.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false), Atom({0.0, 1.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
      {}, {}, 5, 21);
  System system(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, 1.0e5, 1.0, {}, {component}, {}, {1}, 5);

  const Molecule oldMolecule = system.moleculeData[0];
  std::vector<Atom> oldAtoms(system.spanOfMolecule(0, 0).begin(), system.spanOfMolecule(0, 0).end());
  RandomNumber expectedRandom(1234);
  const auto expected =
      system.components[0].rotate(oldMolecule, oldAtoms, expectedRandom.randomRotationMatrix().quaternion());

  RandomNumber moveRandom(1234);
  ASSERT_TRUE(MC_Moves::randomRotationMove(moveRandom, system, 0, 0).has_value());
  const std::span<const Atom> actualAtoms = system.spanOfMolecule(0, 0);
  for (std::size_t i = 0; i < actualAtoms.size(); ++i)
  {
    EXPECT_NEAR(actualAtoms[i].position.x, expected.second[i].position.x, 1.0e-12);
    EXPECT_NEAR(actualAtoms[i].position.y, expected.second[i].position.y, 1.0e-12);
    EXPECT_NEAR(actualAtoms[i].position.z, expected.second[i].position.z, 1.0e-12);
  }
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.x, oldMolecule.centerOfMassPosition.x, 1.0e-12);
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.y, oldMolecule.centerOfMassPosition.y, 1.0e-12);
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.z, oldMolecule.centerOfMassPosition.z, 1.0e-12);
}

TEST(MC_COMPONENT_MOVES, initialization_routes_fractional_swap_handlers)
{
  const std::array moves{Move::Types::SwapCFCMC, Move::Types::SwapCBCFCMC, Move::Types::PairSwapCFCMC,
                         Move::Types::PairSwapCBCFCMC};

  for (Move::Types move : moves)
  {
    SCOPED_TRACE(std::to_underlying(move));
    System system = makeInitializationRoutingSystem(move);
    RandomNumber random(37);
    std::size_t fractionalMoleculeSystem = 0;

    MC_Moves::performRandomMoveInitialization(random, system, system, 0, fractionalMoleculeSystem);

    EXPECT_GT(totalTrials(system, move), 0.0);
    EXPECT_EQ(totalTrials(system, Move::Types::SwapCBMC), 0.0);
    EXPECT_EQ(totalTrials(system, Move::Types::PairSwapCBMC), 0.0);
  }
}

TEST(MC_COMPONENT_MOVES, neutral_trial_updates_tmmc_diagonal)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Widom, 1.0);
  Component component = makeIonComponent(forceField, 0, "A", 0, probabilities);
  System system(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, 1.0e5, 1.0, {}, {component}, {}, {1}, 5);

  system.tmmc.doTMMC = true;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 2;
  system.tmmc.initialize();

  RandomNumber random(37);
  std::size_t fractionalMoleculeSystem = 0;
  EXPECT_EQ(MC_Moves::performRandomMoveInitialization(random, system, system, 0, fractionalMoleculeSystem),
            Move::Types::Widom);

  EXPECT_DOUBLE_EQ(system.tmmc.cmatrix[1].x, 0.0);
  EXPECT_DOUBLE_EQ(system.tmmc.cmatrix[1].y, 1.0);
  EXPECT_DOUBLE_EQ(system.tmmc.cmatrix[1].z, 0.0);
}
