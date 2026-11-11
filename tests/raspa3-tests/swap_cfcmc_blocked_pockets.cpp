#include <gtest/gtest.h>

import std;

import double3;
import double4;
import int3;
import randomnumbers;
import atom;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_swap_cfcmc;

// The blocking-pocket test of the conventional CFCMC swap: the pocket radius scales with the VDW coupling of the
// atom tested, so the molecule that becomes fractional at the new lambda is what both crossings must test. A
// deletion crossing tested the molecule it had just decoupled (never inside any pocket) and not the molecule it
// demotes, so it could demote a molecule into a pocket that an insertion at the same lambda would refuse.

namespace
{

constexpr std::size_t methaneType = 3;  // CH4 in the zeolite force field

System makeMethaneInFAU(std::size_t molecules)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, methaneType, 0, false, false)}, {}, {}, 5, 21, probabilities,
                          std::nullopt, false);
  return System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {framework}, {methane}, {}, {molecules}, 5);
}

// parks the fractional molecule decoupled at the lowest lambda bin: every attempt from there is an upward hop
// or a deletion crossing
void parkFractionalMoleculeDecoupled(System& system)
{
  const std::size_t slot = system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, 0);
  const std::uint8_t groupId = system.components[0].lambdaGC.dUdlambdaGroupId;
  for (Atom& atom : system.spanOfMolecule(0, slot)) atom.setScalingToFractional(0.0, groupId);
  system.components[0].lambdaGC.setCurrentBin(0);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
}

// the lowest integer molecule count seen over the attempts: a deletion crossing lowers it by one
std::size_t lowestIntegerCountOverSwaps(System& system, std::size_t attempts, std::size_t seed)
{
  RandomNumber random(seed);
  std::size_t lowest = system.numberOfIntegerMoleculesPerComponent[0];
  for (std::size_t attempt = 0; attempt < attempts; ++attempt)
  {
    const auto [energy, acceptance] = MC_Moves::swapMove_CFCMC(random, system, 0, false, false);
    if (energy.has_value()) system.runningEnergies += energy.value();
    lowest = std::min(lowest, system.numberOfIntegerMoleculesPerComponent[0]);
  }
  return lowest;
}

}  // namespace

TEST(SWAP_CFCMC_BLOCKED_POCKETS, a_deletion_crossing_tests_the_demoted_molecule_at_its_new_coupling)
{
  constexpr std::size_t molecules = 20;

  // a pocket whose radius covers the whole FAU cell at (nearly) full coupling (the radius is scaled by the
  // square root of scalingVDW): the molecule a deletion crossing demotes to lambda = 1 - epsilon is inside it
  // wherever it is, so no deletion crossing may ever be accepted
  System blocked = makeMethaneInFAU(molecules);
  blocked.components[0].blockingPockets.push_back(double4(0.5, 0.5, 0.5, 100.0));
  parkFractionalMoleculeDecoupled(blocked);
  ASSERT_EQ(blocked.numberOfIntegerMoleculesPerComponent[0], molecules);
  EXPECT_EQ(lowestIntegerCountOverSwaps(blocked, 3000, 7), molecules);

  // without the pocket the same attempts do delete molecules (the deletion crossing itself works)
  System open = makeMethaneInFAU(molecules);
  parkFractionalMoleculeDecoupled(open);
  EXPECT_LT(lowestIntegerCountOverSwaps(open, 3000, 7), molecules);
}
