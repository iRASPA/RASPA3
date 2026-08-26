#include <gtest/gtest.h>

import std;

import int3;
import double3;
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
import move_statistics;
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

// The identity-switch move exchanges the identities of one molecule of component A and one
// molecule of component B in the same box (the MCCCS-MN "swatch"). Because a pair is exchanged,
// the number of molecules of every component is conserved and the acceptance rule contains no
// fugacities: only the Rosenbluth ratio of the two grows over the two retraces, corrected for
// the Fourier-space Ewald difference, the polarization difference and the direct interaction
// between the two exchanged molecules.

// Two components built from the same pseudo-atom are physically identical and differ only in
// their label, so every factor in the acceptance rule cancels exactly and the move is always
// accepted. The composition must stay exactly where it started: unlike the semi-grand identity
// change, this move cannot convert one component into the other.
TEST(MC_IDENTITY_SWITCH_DRIFT, identical_components_are_always_accepted)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);

  const std::optional<std::size_t> typeCH4 = forceField.findPseudoAtom("CH4");
  ASSERT_TRUE(typeCH4.has_value());

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::IdentitySwitchCBMC, 1.0);

  Component beadA = Component::makeIon(forceField, 0, "beadA", typeCH4.value(), 0.0);
  beadA.mc_moves_probabilities = probabilities;
  beadA.identitySwitches = {1};

  Component beadB = Component::makeIon(forceField, 1, "beadB", typeCH4.value(), 0.0);
  beadB.mc_moves_probabilities = probabilities;
  beadB.identitySwitches = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  const std::size_t initialA{30};
  const std::size_t initialB{10};

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {beadA, beadB}, {}, {initialA, initialB}, 5,
                         MCMoveProbabilities());

  std::vector<System> systems{system};
  std::size_t numberOfProductionCycles{200};
  std::size_t numberOfInitializationCycles{50};
  std::size_t numberOfEquilibrationCycles{50};
  std::size_t printEvery{1000};
  std::size_t writeBinaryRestartEvery{10000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles,
                              printEvery, writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery},
                             systems, 42uz, numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    // the move conserves the number of molecules of each component
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], initialA);
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[1], initialB);

    for (std::size_t componentId = 0; componentId != 2; ++componentId)
    {
      const MoveStatistics<double> &statistics = std::get<MoveStatistics<double>>(
          s.components[componentId].mc_moves_statistics[Move::Types::IdentitySwitchCBMC]);

      // the move must actually have run, otherwise the check below is vacuous
      EXPECT_GT(statistics.totalConstructed, 0.0);

      // identical components: all Rosenbluth weights and energy differences cancel, so every
      // constructed trial is accepted
      EXPECT_DOUBLE_EQ(statistics.totalAccepted, statistics.totalConstructed);
    }

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

// Runs a box of CO2 and water molecules that switch identities with each other and checks that the
// running energies still match a full recomputation.
static void expectNoDriftForChargedMultisiteSwitch(const ForceField &forceField)
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::IdentitySwitchCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities;
  co2.identitySwitches = {1};

  Component water = Component::makeWater(forceField, 1, true);
  water.mc_moves_probabilities = probabilities;
  water.identitySwitches = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  const std::size_t initialCO2{20};
  const std::size_t initialWater{20};

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2, water}, {}, {initialCO2, initialWater}, 5,
                         MCMoveProbabilities());

  std::vector<System> systems{system};
  std::size_t numberOfProductionCycles{200};
  std::size_t numberOfInitializationCycles{50};
  std::size_t numberOfEquilibrationCycles{50};
  std::size_t printEvery{1000};
  std::size_t writeBinaryRestartEvery{10000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles,
                              printEvery, writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery},
                             systems, 42uz, numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], initialCO2);
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[1], initialWater);

    const MoveStatistics<double> &statistics =
        std::get<MoveStatistics<double>>(s.components[0].mc_moves_statistics[Move::Types::IdentitySwitchCBMC]);
    EXPECT_GT(statistics.totalAccepted, 0.0);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);

    if (forceField.computePolarization)
    {
      // the fixture must actually exercise a non-trivial polarization energy
      EXPECT_LT(recomputedEnergies.polarization, -1e-8);
    }
  }
}

// Charged, multi-site molecules of different sizes. Both molecules change at the same time, so
// the Ewald difference is evaluated for a span covering two molecules at once: the intramolecular
// exclusion terms must be matched per molecule and the cross terms between the two exchanged
// molecules must not be counted as exclusions. This is the case that a single-molecule move never
// exercises; an error in it shows up as accumulating energy drift rather than as a wrong average.
TEST(MC_IDENTITY_SWITCH_DRIFT, charged_multisite_components_energy_drift)
{
  expectNoDriftForChargedMultisiteSwitch(ForceField::makeZeoliteForceField(12.0, true, false, true));
}

// Same system with the dual cut-off scheme. The two exchanged molecules are grown and retraced in
// a nested background, so their mutual interaction lives inside the second Rosenbluth weight on
// each side of the move and is carried from the inner to the full cut-off by that molecule's dual
// cut-off correction. An intra-pair term left behind at the inner cut-off, or a retrace that is
// not corrected while its matching grow is, shows up here as energy drift.
TEST(MC_IDENTITY_SWITCH_DRIFT, charged_multisite_components_dual_cut_off_energy_drift)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.useDualCutOff = true;
  forceField.dualCutOff = 6.0;

  expectNoDriftForChargedMultisiteSwitch(forceField);
}

// Same system with molecule-molecule polarization. Two molecules change at once, so both the field
// on each of them (which must come from the surviving molecules plus its new partner, not from the
// molecules being replaced) and the field change felt by every surrounding molecule have to be
// tracked incrementally. A missing or double-counted contribution shows up as polarization drift.
TEST(MC_IDENTITY_SWITCH_DRIFT, charged_multisite_components_polarization_energy_drift)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  expectNoDriftForChargedMultisiteSwitch(forceField);
}
