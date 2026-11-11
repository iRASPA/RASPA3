#include <gtest/gtest.h>

import std;

import atom;
import component;
import double3;
import forcefield;
import mc_moves;
import mc_moves_swap_cfcmc_cbmc;
import mc_moves_group_swap;
import mc_moves_tethered_proton_hop;
import mc_moves_move_types;
import mc_moves_pair_deletion_cbmc;
import mc_moves_pair_insertion_cbmc;
import mc_moves_probabilities;
import mc_moves_random_rotation;
import mc_moves_rotation_smart_mc;
import mc_moves_translation_rotation_smart_mc;
import mc_moves_translation_smart_mc_all;
import mc_moves_rotation_smart_mc_all;
import mc_moves_statistics;
import move_statistics;
import molecule;
import randomnumbers;
import running_energy;
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

System makeGroupSystem(const std::vector<double3>& positionsA, const std::vector<double3>& positionsB,
                       double pressure)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::GroupSwap, 1.0);

  Component componentA = makeIonComponent(forceField, 0, "A", 0, probabilities);
  Component componentB = makeIonComponent(forceField, 1, "B", 1, {});
  componentA.groupComponentIds = {1, 1};
  componentA.maximumGroupDistance = 2.0;

  return System(forceField, SimulationBox(10.0, 10.0, 10.0), false, 300.0, pressure, 1.0, {},
                {componentA, componentB}, {positionsA, positionsB}, {0, 0}, 5);
}

System makeTetheredProtonHopSystem(const std::vector<double3>& cartesianSites)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::TetheredProtonHop, 1.0);

  Component proton = makeIonComponent(forceField, 0, "A", 0, probabilities);

  // The move reads the candidate sites as fractional coordinates in the simulation box; convert the
  // requested Cartesian positions with the inverse cell so the test can reason in Cartesian space.
  const SimulationBox box(10.0, 10.0, 10.0);
  std::vector<double3> fractionalSites;
  fractionalSites.reserve(cartesianSites.size());
  for (const double3& site : cartesianSites)
  {
    fractionalSites.push_back(box.inverseCell * site);
  }
  proton.tetheredProtonHopSiteGroups = {fractionalSites};

  // A single proton placed on the first candidate site.
  return System(forceField, box, false, 300.0, 1.0e5, 1.0, {}, {proton}, {{cartesianSites.front()}}, {0}, 5);
}

System makeAutoPlacedProtonSystem(const std::vector<std::vector<double3>>& cartesianGroups)
{
  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::TetheredProtonHop, 1.0);

  Component proton = makeIonComponent(forceField, 0, "A", 0, probabilities);

  const SimulationBox box(10.0, 10.0, 10.0);
  std::vector<std::vector<double3>> fractionalGroups;
  for (const std::vector<double3>& group : cartesianGroups)
  {
    std::vector<double3> fractionalGroup;
    for (const double3& site : group) fractionalGroup.push_back(box.inverseCell * site);
    fractionalGroups.push_back(std::move(fractionalGroup));
  }
  proton.tetheredProtonHopSiteGroups = fractionalGroups;

  // No explicit initial positions and a zero create-count: the protons must be placed from the
  // tether groups alone (one proton per group, on the group's first site).
  return System(forceField, box, false, 300.0, 1.0e5, 1.0, {}, {proton}, {}, {0}, 5);
}

double totalTrials(const System& system, Move::Types move)
{
  const auto& statistics = std::get<MoveStatistics<double3>>(system.components[0].mc_moves_statistics[move]);
  return statistics.totalCounts.x + statistics.totalCounts.y + statistics.totalCounts.z;
}

void expectAtomScalingState(const std::vector<Atom>& actual, const std::vector<Atom>& expected)
{
  ASSERT_EQ(actual.size(), expected.size());
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(actual[i].scalingVDW, expected[i].scalingVDW);
    EXPECT_DOUBLE_EQ(actual[i].scalingCoulomb, expected[i].scalingCoulomb);
    EXPECT_EQ(actual[i].groupId, expected[i].groupId);
    EXPECT_EQ(actual[i].isFractional, expected[i].isFractional);
  }
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

TEST(MC_COMPONENT_MOVES, pair_insertion_acceptance_matches_grand_canonical_prefactor)
{
  // Insertion into an empty, non-interacting system: all energies and Rosenbluth ratios are one,
  // and the only in-range reverse partner is the inserted molecule B itself (k = 1). The
  // acceptance must then equal the bare grand-canonical prefactor
  //   (beta f_A V / (N_A+1)) * (beta f_B V_s b / k)
  // with V_s = 4 pi R_max^3 / 3 the pair-sphere volume and b the radial-proposal bias
  // (b = 1 for the conventional move, b = 3 r^2 / R_max^2 for the CBMC move).
  constexpr double R_max = 2.0;
  const double sphereVolume = (4.0 / 3.0) * std::numbers::pi * R_max * R_max * R_max;

  {
    System system = makePairInsertionSystem({}, 1.0e-12);
    RandomNumber random(918);
    const auto [energy, acceptance] = MC_Moves::pairInsertionMove(random, system, 0);
    EXPECT_FALSE(energy.has_value());
    const double expected = (system.beta * system.pressure * system.simulationBox.volume) *
                            (system.beta * system.pressure * sphereVolume);
    EXPECT_NEAR(acceptance.z / expected, 1.0, 1.0e-10);
  }

  {
    // High pressure so the CBMC insertion is accepted and the drawn pair distance r can be
    // recovered from the inserted positions.
    System system = makePairInsertionSystem({}, 1.0e12);
    RandomNumber random(918);
    const auto [energy, acceptance] = MC_Moves::pairInsertionMoveCBMC(random, system, 0);
    ASSERT_TRUE(energy.has_value());
    const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(
        system.spanOfMolecule(0, 0).front().position - system.spanOfMolecule(1, 0).front().position);
    const double r = dr.length();
    const double expected = (system.beta * system.pressure * system.simulationBox.volume) *
                            (system.beta * system.pressure * sphereVolume * (3.0 * r * r / (R_max * R_max)));
    EXPECT_NEAR(acceptance.z / expected, 1.0, 1.0e-10);
  }
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

TEST(MC_COMPONENT_MOVES, group_deletion_uses_sequential_candidate_counts)
{
  // Central molecule at the box center with two satellite slots of the same component. With two
  // in-range satellites the sequential candidate counts are (2, 1); with three they are (3, 2).
  // The candidate counts enter the deletion acceptance in the numerator, so the acceptance must
  // differ by a factor (3*2)/(2*1) = 3.
  const auto checkCandidateCounts = [](auto deletionMove)
  {
    System twoInRange = makeGroupSystem(
        {double3(5.0, 5.0, 5.0)},
        {double3(6.5, 5.0, 5.0), double3(3.5, 5.0, 5.0), double3(5.0, 5.0, 9.0)}, 1.0e12);
    System threeInRange = makeGroupSystem(
        {double3(5.0, 5.0, 5.0)},
        {double3(6.5, 5.0, 5.0), double3(3.5, 5.0, 5.0), double3(5.0, 6.5, 5.0)}, 1.0e12);
    RandomNumber randomTwo(91);
    RandomNumber randomThree(91);

    const auto [energyTwo, acceptanceTwo] = deletionMove(randomTwo, twoInRange);
    const auto [energyThree, acceptanceThree] = deletionMove(randomThree, threeInRange);

    EXPECT_FALSE(energyTwo.has_value());
    EXPECT_FALSE(energyThree.has_value());
    ASSERT_GT(acceptanceTwo.x, 0.0);
    EXPECT_NEAR(acceptanceThree.x / acceptanceTwo.x, 3.0, 1.0e-9);
  };

  checkCandidateCounts([](RandomNumber& random, System& system)
                       { return MC_Moves::groupDeletionMove(random, system, 0, 0); });
  // All in-range satellites sit at the same distance, so the CBMC distance-bias factors are equal
  // in both systems and the ratio isolates the same candidate-count factor.
  checkCandidateCounts([](RandomNumber& random, System& system)
                       { return MC_Moves::groupDeletionMoveCBMC(random, system, 0, 0); });
}

TEST(MC_COMPONENT_MOVES, group_insertion_deletion_proposals_are_reciprocal)
{
  // Insert a full group into an empty system and immediately propose the reverse deletion of the
  // same molecules; the product of the two acceptance probabilities must be exactly one because
  // both moves are evaluated in the same configuration (the energies are zero for the
  // non-interacting force field).
  const auto checkReciprocity = [](auto insertionMove, auto deletionMove)
  {
    System system = makeGroupSystem({}, {}, 1.0e12);
    RandomNumber insertionRandom(918);
    const auto [insertionEnergy, insertionAcceptance] = insertionMove(insertionRandom, system);
    ASSERT_TRUE(insertionEnergy.has_value());
    ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 1uz);
    ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 2uz);

    RandomNumber deletionRandom(41);
    const auto [deletionEnergy, deletionAcceptance] = deletionMove(deletionRandom, system, 0, 0);
    EXPECT_FALSE(deletionEnergy.has_value());
    EXPECT_NEAR(insertionAcceptance.z * deletionAcceptance.x, 1.0, 1.0e-9);
  };

  checkReciprocity([](RandomNumber& random, System& system)
                   { return MC_Moves::groupInsertionMove(random, system, 0); },
                   [](RandomNumber& random, System& system, std::size_t component, std::size_t molecule)
                   { return MC_Moves::groupDeletionMove(random, system, component, molecule); });
  checkReciprocity([](RandomNumber& random, System& system)
                   { return MC_Moves::groupInsertionMoveCBMC(random, system, 0); },
                   [](RandomNumber& random, System& system, std::size_t component, std::size_t molecule)
                   { return MC_Moves::groupDeletionMoveCBMC(random, system, component, molecule); });
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_requires_at_least_two_sites)
{
  // A group with a single candidate site offers no alternative position, so the move is a no-op.
  System system = makeTetheredProtonHopSystem({double3(0.5, 5.0, 5.0)});
  RandomNumber random(91);
  EXPECT_FALSE(MC_Moves::tetheredProtonHopMove(random, system, 0, 0).has_value());
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_sites_place_one_proton_per_group)
{
  // The number of protons is derived from the tether groups (one per group), each placed on the
  // group's first site, independent of any create-count or initial-position input.
  const std::vector<std::vector<double3>> groups{
      {double3(0.5, 5.0, 5.0), double3(1.5, 5.0, 5.0)},
      {double3(5.0, 8.0, 5.0), double3(5.0, 9.0, 5.0), double3(5.0, 7.0, 5.0)}};

  System system = makeAutoPlacedProtonSystem(groups);

  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 2uz);
  EXPECT_LT((system.spanOfMolecule(0, 0).front().position - groups[0].front()).length(), 1.0e-9);
  EXPECT_LT((system.spanOfMolecule(0, 1).front().position - groups[1].front()).length(), 1.0e-9);

  // the degrees of freedom are those of the two inserted protons, counted once (by the insertion itself)
  EXPECT_EQ(system.translationalDegreesOfFreedom, 2uz * system.components[0].translationalDegreesOfFreedom);
  EXPECT_EQ(system.rotationalDegreesOfFreedom, 2uz * system.components[0].rotationalDegreesOfFreedom);
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_relocates_to_a_different_site)
{
  // Non-interacting force field: the energy difference is zero, so the symmetric proposal is always
  // accepted. The proton must land exactly on one of the other candidate sites, never the one it
  // started on.
  const std::vector<double3> sites{double3(0.5, 5.0, 5.0), double3(1.5, 5.0, 5.0), double3(2.5, 5.0, 5.0)};

  for (std::size_t seed = 0; seed < 16; ++seed)
  {
    System system = makeTetheredProtonHopSystem(sites);
    RandomNumber random(seed);

    ASSERT_TRUE(MC_Moves::tetheredProtonHopMove(random, system, 0, 0).has_value());

    const double3 position = system.spanOfMolecule(0, 0).front().position;
    EXPECT_GT((position - sites[0]).length(), 1.0e-9);
    const bool onSite1 = (position - sites[1]).length() < 1.0e-9;
    const bool onSite2 = (position - sites[2]).length() < 1.0e-9;
    EXPECT_TRUE(onSite1 || onSite2);
  }
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_selects_each_alternative_site_across_seeds)
{
  // With two alternative sites the uniform proposal must reach both across seeds.
  const std::vector<double3> sites{double3(0.5, 5.0, 5.0), double3(1.5, 5.0, 5.0), double3(2.5, 5.0, 5.0)};

  std::size_t selectedSite1 = 0;
  std::size_t selectedSite2 = 0;
  for (std::size_t seed = 0; seed < 64; ++seed)
  {
    System system = makeTetheredProtonHopSystem(sites);
    RandomNumber random(seed);
    ASSERT_TRUE(MC_Moves::tetheredProtonHopMove(random, system, 0, 0).has_value());

    const double3 position = system.spanOfMolecule(0, 0).front().position;
    selectedSite1 += (position - sites[1]).length() < 1.0e-9;
    selectedSite2 += (position - sites[2]).length() < 1.0e-9;
  }

  EXPECT_EQ(selectedSite1 + selectedSite2, 64uz);
  EXPECT_GE(selectedSite1, 20uz);
  EXPECT_GE(selectedSite2, 20uz);
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_preserves_the_thermal_offset_and_is_reversible)
{
  // Combined with translation moves the proton sits off-centre at site + delta. The hop must carry the
  // offset along (site-to-site displacement) rather than snap onto the target centre: a snapped position
  // could never be regenerated by the reverse hop, which breaks detailed balance.
  const std::vector<double3> sites{double3(0.5, 5.0, 5.0), double3(1.5, 5.0, 5.0), double3(2.5, 5.0, 5.0)};
  const double3 delta(0.12, -0.07, 0.21);

  bool sawEveryTarget[3] = {true, false, false};
  for (std::size_t seed = 0; seed < 32; ++seed)
  {
    System system = makeTetheredProtonHopSystem(sites);
    system.spanOfMolecule(0, 0).front().position = sites[0] + delta;
    system.moleculeData[0].centerOfMassPosition = sites[0] + delta;
    const double3 start = system.spanOfMolecule(0, 0).front().position;

    RandomNumber random(seed);
    ASSERT_TRUE(MC_Moves::tetheredProtonHopMove(random, system, 0, 0).has_value());

    const double3 hopped = system.spanOfMolecule(0, 0).front().position;
    std::size_t target = 0;
    for (std::size_t i = 1; i < sites.size(); ++i)
    {
      if ((hopped - (sites[i] + delta)).length() < 1.0e-12) target = i;
    }
    ASSERT_NE(target, 0uz) << "the hop must land at a different site with the same offset";
    sawEveryTarget[target] = true;
    EXPECT_GT((hopped - sites[target]).length(), 0.1);
    EXPECT_LT((system.moleculeData[0].centerOfMassPosition - hopped).length(), 1.0e-12);

    // the reverse hop (target -> original site) restores the starting position exactly
    bool reversed = false;
    for (std::size_t attempt = 0; attempt < 64 && !reversed; ++attempt)
    {
      ASSERT_TRUE(MC_Moves::tetheredProtonHopMove(random, system, 0, 0).has_value());
      const double3 position = system.spanOfMolecule(0, 0).front().position;
      if ((position - start).length() < 1.0e-12)
      {
        reversed = true;
      }
    }
    EXPECT_TRUE(reversed) << "the reverse hop never regenerated the starting position (seed " << seed << ")";
  }
  EXPECT_TRUE(sawEveryTarget[1]);
  EXPECT_TRUE(sawEveryTarget[2]);
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_rejects_when_the_offset_would_belong_to_another_site)
{
  // Sites at x = 0.5, 1.5 and 1.9: a proton at 0.5 + 0.3 belongs to site 0. Hopping to site 1 would put it
  // at 1.8, which is nearest to site 2, so the reverse hop would never be proposed from there: the move
  // must reject that target. Hopping to site 2 puts it at 2.2, which belongs to site 2 and is accepted.
  const std::vector<double3> sites{double3(0.5, 5.0, 5.0), double3(1.5, 5.0, 5.0), double3(1.9, 5.0, 5.0)};
  const double3 start(0.8, 5.0, 5.0);

  std::size_t accepted = 0;
  std::size_t rejected = 0;
  for (std::size_t seed = 0; seed < 64; ++seed)
  {
    System system = makeTetheredProtonHopSystem(sites);
    system.spanOfMolecule(0, 0).front().position = start;
    system.moleculeData[0].centerOfMassPosition = start;

    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::tetheredProtonHopMove(random, system, 0, 0);
    const double3 position = system.spanOfMolecule(0, 0).front().position;
    if (result.has_value())
    {
      ++accepted;
      EXPECT_LT((position - double3(2.2, 5.0, 5.0)).length(), 1.0e-12);
    }
    else
    {
      ++rejected;
      EXPECT_LT((position - start).length(), 1.0e-12);
    }
  }
  EXPECT_GE(accepted, 16uz);
  EXPECT_GE(rejected, 16uz);
}

TEST(MC_COMPONENT_MOVES, tethered_proton_hop_maintains_the_cell_list)
{
  // The hop updates the molecule in place like a translation: it reports the moved atoms to the cell list
  // and is listed as a maintaining move, so the driver does not rebuild the whole list after every attempt.
  EXPECT_TRUE(MC_Moves::maintainsCellList(Move::Types::TetheredProtonHop));

  const ForceField forceField = makeNonInteractingForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::TetheredProtonHop, 1.0);
  Component proton = makeIonComponent(forceField, 0, "A", 0, probabilities);

  // a box large enough for the list to be enabled (at least four cells of the cut-off per direction)
  const SimulationBox box(20.0, 20.0, 20.0);
  std::vector<std::vector<double3>> groups;
  for (std::size_t g = 0; g < 8; ++g)
  {
    const double3 base(2.0 + 2.0 * static_cast<double>(g), 3.0, 4.0);
    groups.push_back({box.inverseCell * base, box.inverseCell * (base + double3(0.0, 1.0, 0.0)),
                      box.inverseCell * (base + double3(0.0, 0.0, 1.0))});
  }
  proton.tetheredProtonHopSiteGroups = groups;
  System system(forceField, box, false, 300.0, 1.0e5, 1.0, {}, {proton}, {}, {0}, 5);
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 8uz);

  ASSERT_TRUE(system.cellList().enabled);
  RandomNumber random(7);
  for (std::size_t step = 0; step < 200; ++step)
  {
    const std::size_t selectedMolecule = random.uniform_integer(0, 7);
    MC_Moves::tetheredProtonHopMove(random, system, 0, selectedMolecule);
    EXPECT_TRUE(system.moleculeCellList.isCurrent(system.spanOfMoleculeAtoms(), system.simulationBox,
                                                  system.cellListCutOff()));
    ASSERT_TRUE(system.verifyCellList());
  }
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
      system.components[0].rotate(oldMolecule, oldAtoms, expectedRandom.randomRotationMatrix().transpose().quaternion());

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

TEST(MC_COMPONENT_MOVES, cfcmc_cbmc_restores_atoms_after_upper_tmmc_window_rejection)
{
  System system = makeInitializationRoutingSystem(Move::Types::SwapCBCFCMC);
  Component& component = system.components[0];
  const std::size_t fractionalIndex =
      system.indexOfFractionalMoleculeForMove(Move::Types::SwapCBCFCMC, 0);
  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[0];

  system.tmmc.doTMMC = true;
  system.tmmc.rejectOutOfBound = true;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = oldN;
  system.tmmc.numberOfLambdaBins = component.lambdaGC.numberOfSamplePoints;
  system.tmmc.initialize();

  for (Atom& atom : system.spanOfMolecule(0, fractionalIndex))
  {
    atom.setScalingToFractional(1.0, component.lambdaGC.dUdlambdaGroupId);
  }
  component.lambdaGC.setCurrentBin(component.lambdaGC.numberOfSamplePoints - 1uz);
  const std::vector<Atom> before = system.atomData;

  RandomNumber random(91);
  bool testedBoundary = false;
  for (std::size_t attempt = 0; attempt < 100uz && !testedBoundary; ++attempt)
  {
    system.atomData = before;
    component.lambdaGC.setCurrentBin(component.lambdaGC.numberOfSamplePoints - 1uz);
    const auto [energy, pacc] = MC_Moves::swapMove_CFCMC_CBMC(random, system, 0, fractionalIndex);
    if (pacc.z > 0.0)
    {
      EXPECT_FALSE(energy.has_value());
      expectAtomScalingState(system.atomData, before);
      testedBoundary = true;
    }
  }
  EXPECT_TRUE(testedBoundary);
}

TEST(MC_COMPONENT_MOVES, cfcmc_cbmc_restores_atoms_after_lower_tmmc_window_rejection)
{
  System system = makeInitializationRoutingSystem(Move::Types::SwapCBCFCMC);
  Component& component = system.components[0];
  const std::size_t fractionalIndex =
      system.indexOfFractionalMoleculeForMove(Move::Types::SwapCBCFCMC, 0);
  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[0];

  system.tmmc.doTMMC = true;
  system.tmmc.rejectOutOfBound = true;
  system.tmmc.minMacrostate = oldN;
  system.tmmc.maxMacrostate = oldN + 2uz;
  system.tmmc.numberOfLambdaBins = component.lambdaGC.numberOfSamplePoints;
  system.tmmc.initialize();

  for (Atom& atom : system.spanOfMolecule(0, fractionalIndex))
  {
    atom.setScalingToFractional(0.0, component.lambdaGC.dUdlambdaGroupId);
  }
  component.lambdaGC.setCurrentBin(0uz);
  const std::vector<Atom> before = system.atomData;

  RandomNumber random(37);
  bool testedBoundary = false;
  for (std::size_t attempt = 0; attempt < 100uz && !testedBoundary; ++attempt)
  {
    system.atomData = before;
    component.lambdaGC.setCurrentBin(0uz);
    const auto [energy, pacc] = MC_Moves::swapMove_CFCMC_CBMC(random, system, 0, fractionalIndex);
    if (pacc.x > 0.0)
    {
      EXPECT_FALSE(energy.has_value());
      expectAtomScalingState(system.atomData, before);
      testedBoundary = true;
    }
  }
  EXPECT_TRUE(testedBoundary);
}

TEST(MC_COMPONENT_MOVES, neutral_trial_leaves_tmmc_matrix_untouched)
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

  // A move that cannot change N is not a transition of the macrostate chain: no row is touched.
  for (const double3& row : system.tmmc.cmatrix)
  {
    EXPECT_DOUBLE_EQ(row.x, 0.0);
    EXPECT_DOUBLE_EQ(row.y, 0.0);
    EXPECT_DOUBLE_EQ(row.z, 0.0);
  }
}

// The torque-biased rotation draws a rotation vector from a Gaussian about the drift. Every vector along the same
// axis whose angle differs by a multiple of 2 pi is the same rotation, so the proposal density of the rotation is
// the sum over those aliases, each weighted by theta^2 (the Haar measure in exponential coordinates). For a small
// step the principal term is all there is and the ratio is that of the two Gaussians.
TEST(MC_COMPONENT_MOVES, smart_rotation_proposal_ratio_is_the_gaussian_ratio_for_a_small_step)
{
  const double3 step(0.05, 0.02, -0.01);
  const double3 driftOld(0.01, -0.02, 0.03);
  const double3 driftNew(-0.02, 0.01, 0.0);
  const double sigma = 0.1;
  const double gaussian =
      ((step - driftOld).length_squared() - (step + driftNew).length_squared()) / (2.0 * sigma * sigma);
  EXPECT_NEAR(MC_Moves::logRotationProposalRatio(step, driftOld, driftNew, sigma), gaussian, 1e-12);
}

// A step and its reverse: the ratio of the reverse step is the inverse of the ratio of the step, for any step size.
TEST(MC_COMPONENT_MOVES, smart_rotation_proposal_ratio_of_the_reverse_step_is_the_inverse)
{
  const double3 step(1.9, -1.2, 0.8);  // 2.4 radians: beyond pi with its aliases in reach
  const double3 driftOld(0.4, 0.1, -0.2);
  const double3 driftNew(-0.3, 0.5, 0.2);
  const double sigma = 1.5;
  EXPECT_NEAR(MC_Moves::logRotationProposalRatio(step, driftOld, driftNew, sigma),
              -MC_Moves::logRotationProposalRatio(-step, driftNew, driftOld, sigma), 1e-12);
}

// At the largest step size the aliased rotation vectors carry weight: the ratio differs from the principal
// Gaussian ratio, and equals the sum carried out over many more aliases.
TEST(MC_COMPONENT_MOVES, smart_rotation_proposal_ratio_counts_the_aliased_rotation_vectors)
{
  const double3 step(2.6, 0.0, 0.0);
  const double3 driftOld(0.4, 0.0, 0.0);
  const double3 driftNew(0.9, 0.3, 0.0);  // the reverse step's alias at -3.7 radians lies close to its drift
  const double sigma = 1.5;
  const double gaussian =
      ((step - driftOld).length_squared() - (step + driftNew).length_squared()) / (2.0 * sigma * sigma);
  const double ratio = MC_Moves::logRotationProposalRatio(step, driftOld, driftNew, sigma);
  EXPECT_GT(std::abs(ratio - gaussian), 0.05);

  double forward = 0.0, reverse = 0.0;
  for (int k = -40; k <= 40; ++k)
  {
    const double theta = 2.6 + 2.0 * std::numbers::pi * static_cast<double>(k);
    const double3 vector(theta, 0.0, 0.0);
    forward += theta * theta * std::exp(-(vector - driftOld).length_squared() / (2.0 * sigma * sigma));
    reverse += theta * theta * std::exp(-(-vector - driftNew).length_squared() / (2.0 * sigma * sigma));
  }
  EXPECT_NEAR(ratio, std::log(reverse / forward), 1e-12);
}

// The collective smart moves evaluate the gradients of the trial configuration in place. A refused trial must
// leave the system as it was: the gradients of the old configuration (the next consumer may read them), the
// velocities, and the positions. A box of rigid CO2 with a huge step is refused on overlaps.
TEST(MC_COMPONENT_MOVES, refused_collective_smart_moves_leave_the_gradients_of_the_old_configuration)
{
  for (const bool rotation : {false, true})
  {
    const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
    const Component co2 = Component::makeCO2(forceField, 0, false);
    System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {60}, 5,
                           MCMoveProbabilities());
    system.runningEnergies = system.computeTotalEnergies();
    const Move::Types move = rotation ? Move::Types::RotationSmartMCAll : Move::Types::TranslationSmartMCAll;
    system.mc_moves_statistics.setMaxChange(move, rotation ? 3.0 : 10.0);
    for (std::size_t i = 0; i < system.atomDynamics.size(); ++i)
    {
      system.atomDynamics[i].velocity = double3(static_cast<double>(i), 1.0, -1.0);
    }

    RandomNumber random(5);
    std::optional<RunningEnergy> result = std::nullopt;
    std::size_t refused = 0;
    for (std::size_t attempt = 0; attempt < 50 && refused == 0; ++attempt)
    {
      result = rotation ? MC_Moves::rotationSmartMCMoveAll(random, system)
                        : MC_Moves::translationSmartMCMoveAll(random, system);
      if (!result.has_value()) ++refused;
    }
    ASSERT_EQ(refused, 1uz) << (rotation ? "rotation" : "translation");

    // the gradients are those of the configuration the system is in
    System reference = system;
    reference.precomputeTotalGradients();
    for (std::size_t i = 0; i < system.atomDynamics.size(); ++i)
    {
      EXPECT_LT((system.atomDynamics[i].gradient - reference.atomDynamics[i].gradient).length(), 1e-9)
          << (rotation ? "rotation" : "translation") << " atom " << i;
      EXPECT_EQ(system.atomDynamics[i].velocity, double3(static_cast<double>(i), 1.0, -1.0));
    }
  }
}

// A monoatomic species has no rotational degree of freedom, and the smart rotations refuse it. That does not
// change during a simulation, so the attempts are not trials: counted as such, the moves would report a vanishing
// acceptance and the optimizer would shrink their step for nothing.
TEST(MC_COMPONENT_MOVES, smart_rotations_of_a_monoatomic_species_are_not_counted_as_trials)
{
  const ForceField forceField = makeNonInteractingForceField();
  const Component ion = makeIonComponent(forceField, 0, "A", 0, {});
  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {ion}, {}, {5}, 5,
                         MCMoveProbabilities());
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(11);
  for (std::size_t i = 0; i != 20; ++i)
  {
    EXPECT_FALSE(MC_Moves::rotationSmartMCMove(random, system, 0, 0).has_value());
    EXPECT_FALSE(MC_Moves::translationRotationSmartMCMove(random, system, 0, 0).has_value());
  }
  const MoveStatistics<double>& rotation =
      std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::RotationSmartMC]);
  EXPECT_EQ(rotation.counts, 0.0);
  const MoveStatistics<double3>& translationRotation = std::get<MoveStatistics<double3>>(
      system.components[0].mc_moves_statistics[Move::Types::TranslationRotationSmartMC]);
  EXPECT_EQ(translationRotation.counts.x, 0.0);
  EXPECT_EQ(translationRotation.counts.y, 0.0);
}

// The collective smart rotation turns only the rigid multi-atomic molecules. Without one it has nothing to propose
// and is refused before it counts as a trial, instead of a no-op trial with two gradient evaluations and a total
// energy that is accepted every time and drives the step size to its upper limit.
TEST(MC_COMPONENT_MOVES, collective_smart_rotation_without_a_rigid_polyatomic_molecule_is_not_a_trial)
{
  const ForceField forceField = makeNonInteractingForceField();
  const Component ion = makeIonComponent(forceField, 0, "A", 0, {});
  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {ion}, {}, {5}, 5,
                         MCMoveProbabilities());
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(13);
  for (std::size_t i = 0; i != 20; ++i)
  {
    EXPECT_FALSE(MC_Moves::rotationSmartMCMoveAll(random, system).has_value());
  }
  const MoveStatistics<double>& statistics =
      std::get<MoveStatistics<double>>(system.mc_moves_statistics[Move::Types::RotationSmartMCAll]);
  EXPECT_EQ(statistics.counts, 0.0);
}
