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

// Net-charged species under 'OmitInterInteractions': the net-charge correction of a move is the cross term of the
// moved charge with the framework alone, so the charge handed from one Ewald sub-step to the next stays the
// framework's. The moves added the adsorbate charge between their sub-steps, which put a spurious
// adsorbate-adsorbate term into the Fourier energy. It shows where whole charged molecules are moved one after
// another in one move, as in a group swap of a sodium and a chloride (in the CFCMC crossings the molecules near
// lambda = 0 carry no charge under the staged schedule, so the spurious term vanishes there).
TEST(MC_OMIT_INTER_INTERACTIONS_DRIFT, charged_ion_group_swaps)
{
  for (const Move::Types move : {Move::Types::GroupSwap, Move::Types::GroupSwapCBMC})
  {
    ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
    forceField.omitInterInteractions = true;
    const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
    MCMoveProbabilities probabilities;
    probabilities.setProbability(Move::Types::Translation, 1.0);
    probabilities.setProbability(move, 1.0);
    Component sodium = Component::makeIon(forceField, 0, "sodium", forceField.findPseudoAtom("Na+").value(), 1.0);
    Component chloride = Component::makeIon(forceField, 1, "chloride", forceField.findPseudoAtom("Cl-").value(), -1.0);
    sodium.mc_moves_probabilities = probabilities;
    // The chloride moves too (a component without any move is no longer drawn for a step). The ions do not
    // interact and the pressure is high, so the group insertions pile up ions without bound; sharing the steps
    // with the translations of the chloride keeps the count, and with it the magnitude of the running energy
    // whose rounding the drift test must resolve, where it was.
    MCMoveProbabilities translationOnly;
    translationOnly.setProbability(Move::Types::Translation, 1.0);
    chloride.mc_moves_probabilities = translationOnly;
    sodium.groupComponentIds = {1};
    sodium.maximumGroupDistance = 6.0;
    std::vector<System> systems{
        System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {sodium, chloride}, {}, {4, 4}, 5)};
    MonteCarlo mc({10, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 13uz, 5, false);
    mc.run();
    for (System& s : mc.systems)
    {
      const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
      // the running energy is a sum of thousands of accepted differences of a total of order 1e7 K: its own
      // rounding is resolved relative to that total
      const double tolerance = 1e-6 + 1e-13 * std::abs(s.runningEnergies.potentialEnergy());
      EXPECT_NEAR(drift.ewald_fourier, 0.0, tolerance) << Move::moveNames[std::to_underlying(move)];
      EXPECT_NEAR(drift.potentialEnergy(), 0.0, tolerance) << Move::moveNames[std::to_underlying(move)];
    }
  }
}
