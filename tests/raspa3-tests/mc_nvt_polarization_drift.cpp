#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import randomnumbers;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

// These drift tests validate the incremental polarization bookkeeping of the translation and rotation moves.
// Polarization is enabled together with molecule-molecule polarization (omitInterPolarization == false), so the
// electric field on the moved molecule *and* the change of the field on every other molecule have to be tracked
// incrementally. After the run the running polarization energy must still agree with a full recomputation.

static System makePolarizableCO2System(MCMoveProbabilities probabilities, double pressure,
                                       std::size_t numberOfMolecules)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  c.mc_moves_probabilities = probabilities;

  return System(forceField, std::nullopt, false, 300.0, pressure, 1.0, {f}, {c}, {}, {numberOfMolecules}, 5);
}

static void runAndCheckDrift(MCMoveProbabilities probabilities, double pressure = 1e4,
                             std::size_t numberOfMolecules = 20)
{
  System system = makePolarizableCO2System(probabilities, pressure, numberOfMolecules);

  std::vector<System> systems{system};
  std::size_t numberOfProductionCycles{20};
  std::size_t numberOfInitializationCycles{5};
  std::size_t numberOfEquilibrationCycles{5};
  std::size_t printEvery{1000};
  std::size_t writeBinaryRestartEvery{10000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                 writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {}, numberOfBlocks,
                outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);

    // the fixture must actually exercise a non-trivial polarization energy
    EXPECT_LT(recomputedEnergies.polarization, -1e-8);
  }
}

TEST(MC_NVT_POLARIZATION_DRIFT, translation)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_NVT_POLARIZATION_DRIFT, rotation)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_NVT_POLARIZATION_DRIFT, translation_and_rotation)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 0.5);
  probabilities.setProbability(Move::Types::Rotation, 0.5);
  runAndCheckDrift(probabilities);
}

TEST(MC_MUVT_POLARIZATION_DRIFT, swap)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  runAndCheckDrift(probabilities, 1e5, 4);
}

TEST(MC_MUVT_POLARIZATION_DRIFT, swap_cbmc)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  runAndCheckDrift(probabilities, 1e5, 4);
}

TEST(MC_NVT_POLARIZATION_DRIFT, reinsertion_cbmc)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_MUVT_POLARIZATION_DRIFT, swap_cfcmc)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  runAndCheckDrift(probabilities, 1e5, 4);
}

TEST(MC_MUVT_POLARIZATION_DRIFT, swap_cbcfcmc)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);
  runAndCheckDrift(probabilities, 1e5, 4);
}
