#include <gtest/gtest.h>

import std;

import int3;
import double3;
import randomnumbers;
import atom;
import forcefield;
import framework;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;

// 'OmitInterInteractions' with Ewald electrostatics: the energy is the cross term of the molecules with the
// framework alone, without the molecule-molecule Fourier term, the self energy and the exclusions. The
// incremental moves have to see the same energy as the full evaluation, or the running energy drifts away
// from it from the first accepted move (before, every move carried the full molecule-molecule Fourier
// difference, the self energy and the exclusions of the moved molecule).

namespace
{

System makeOmitSystem(const MCMoveProbabilities& probabilities)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.omitInterInteractions = true;
  forceField.omitEwaldFourier = false;
  Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities;
  return System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {co2}, {}, {12}, 5);
}

void runAndCheckDrift(const MCMoveProbabilities& probabilities, std::size_t seed)
{
  std::vector<System> systems{makeOmitSystem(probabilities)};
  ASSERT_TRUE(systems[0].forceField.usesEwaldFourier());

  MonteCarlo mc({10, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, seed, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    const RunningEnergy recomputed = s.computeTotalEnergies();
    const RunningEnergy drift = s.runningEnergies - recomputed;
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);

    // the fixture exercises the framework-molecule Fourier cross term, and nothing of the molecules alone
    EXPECT_NE(recomputed.ewald_fourier, 0.0);
    EXPECT_DOUBLE_EQ(recomputed.ewald_self, 0.0);
    EXPECT_DOUBLE_EQ(recomputed.ewald_exclusion, 0.0);
    EXPECT_DOUBLE_EQ(recomputed.moleculeMoleculeCharge, 0.0);
    EXPECT_DOUBLE_EQ(recomputed.moleculeMoleculeVDW, 0.0);
  }
}

}  // namespace

TEST(MC_OMIT_INTER_INTERACTIONS_DRIFT, translation_rotation_reinsertion_and_swap)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  runAndCheckDrift(probabilities, 7);
}

TEST(MC_OMIT_INTER_INTERACTIONS_DRIFT, cfcmc_swap_with_a_fractional_molecule)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  runAndCheckDrift(probabilities, 11);
}
