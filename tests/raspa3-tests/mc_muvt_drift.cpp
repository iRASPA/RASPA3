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
import mc_moves_move_types;

// print loading:
// std::print("{}", system.averageLoadings.writeAveragesStatistics(system.components, system.frameworkMass()));
//
// print drifts
// system.runningEnergies.printMCDiff(recomputedEnergies);

TEST(MC_MUVT_DRIFT, insertion)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::Swap, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::Swap, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(Move::Types::Translation, 1.0);
  probabilities_water.setProbability(Move::Types::Rotation, 1.0);
  probabilities_water.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, insertionCBMC)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::SwapCBMC, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::SwapCBMC, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(Move::Types::Translation, 1.0);
  probabilities_water.setProbability(Move::Types::Rotation, 1.0);
  probabilities_water.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, insertionCFCMC)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(Move::Types::Translation, 1.0);
  probabilities_water.setProbability(Move::Types::Rotation, 1.0);
  probabilities_water.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, insertionCBCFCMC)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(Move::Types::Translation, 1.0);
  probabilities_water.setProbability(Move::Types::Rotation, 1.0);
  probabilities_water.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, identity_change)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities_co2;
  co2.identityChanges = {1};
  co2.molFraction = 0.25;
  co2.swappable = true;
  co2.idealGasRosenbluthWeight = 1.0;

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component methane = Component::makeMethane(forceField, 1);
  methane.mc_moves_probabilities = probabilities_methane;
  methane.identityChanges = {0};
  methane.molFraction = 0.75;
  methane.swappable = true;
  methane.idealGasRosenbluthWeight = 1.0;

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane}, {}, {10, 10}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, 42uz,
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.numberOfMoleculesPerComponent[0] + s.numberOfMoleculesPerComponent[1], 20uz);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, identity_change_cfcmc)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::IdentityChangeCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities_co2;
  co2.identityChanges = {1};
  co2.molFraction = 0.25;
  co2.swappable = true;
  co2.idealGasRosenbluthWeight = 1.0;

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::IdentityChangeCBMC, 1.0);
  probabilities_methane.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component methane = Component::makeMethane(forceField, 1);
  methane.mc_moves_probabilities = probabilities_methane;
  methane.identityChanges = {0};
  methane.molFraction = 0.75;
  methane.swappable = true;
  methane.idealGasRosenbluthWeight = 1.0;

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane}, {}, {10, 10}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, 42uz,
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.numberOfFractionalMoleculesPerComponent[0], 1uz);
    EXPECT_EQ(s.numberOfFractionalMoleculesPerComponent[1], 1uz);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

namespace {

std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

ForceField makeZeoliteAlkaneForceField()
{
  return ForceField({{"-", false, 0.0, 0.0, 0.0, 0, false},
                     {"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"Al", true, 26.982, 2.05, 0.0, 13, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"Na+", false, 12.0, 1.0, 0.0, 6, false},
                     {"Cl-", false, 15.9994, -1.0, 0.0, 8, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"C_co2", false, 12.0, 0.6512, 0.0, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.0, 8, false},
                     {"Ow", false, 15.9996, 0.0, 0.0, 8, false},
                     {"Hw", false, 1.0008, 0.241, 0.0, 1, false},
                     {"Lw", false, 0.0, -0.241, 0.0, 0, false},
                     {"probe-He", false, 4.002602, 0.0, 0.0, 2, false},
                     {"probe-Ar", false, 39.948, 0.0, 0.0, 18, false},
                     {"probe-CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"probe-N2", false, 14.00674, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{1.0, 1.0},
                     {22.0, 2.30},
                     {22.0, 2.30},
                     {53.0, 3.30},
                     {15.0966, 2.65755},
                     {142.562, 3.51932},
                     {158.5, 3.72},
                     {29.933, 2.745},
                     {85.671, 3.017},
                     {89.633, 3.097},
                     {0.0, 1.0},
                     {0.0, 1.0},
                     {10.9, 2.64},
                     {124.070, 3.38},
                     {158.5, 3.72},
                     {91.5, 3.681},
                     {98.0, 3.75},
                     {46.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
}

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                     true, false);
}

Component makeAlkaneFromExample(const ForceField &forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities &probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), moleculePath.string(), 5,
                   21, probabilities, std::nullopt, false);
}

}  // namespace

TEST(MC_MUVT_DRIFT, identity_change_co2_propane_butane)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities_co2;
  co2.identityChanges = {1, 2};
  co2.molFraction = 0.25;
  co2.swappable = true;
  co2.idealGasRosenbluthWeight = 1.0;

  MCMoveProbabilities probabilities_propane = probabilities_co2;
  Component propane = makeAlkaneFromExample(forceField, 1, "propane", probabilities_propane);
  propane.mc_moves_probabilities = probabilities_propane;
  propane.identityChanges = {0, 2};
  propane.molFraction = 0.375;
  propane.swappable = true;
  propane.idealGasRosenbluthWeight = 1.0;

  MCMoveProbabilities probabilities_butane = probabilities_co2;
  Component butane = makeAlkaneFromExample(forceField, 2, "butane", probabilities_butane);
  butane.mc_moves_probabilities = probabilities_butane;
  butane.identityChanges = {0, 1};
  butane.molFraction = 0.375;
  butane.swappable = true;
  butane.idealGasRosenbluthWeight = 1.0;

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, propane, butane}, {}, {5, 5, 5}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, 42uz,
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.numberOfMoleculesPerComponent[0] + s.numberOfMoleculesPerComponent[1] +
                  s.numberOfMoleculesPerComponent[2],
              15uz);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, identity_change_propane_butane_tail_corrections)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
  propane.identityChanges = {1};
  butane.identityChanges = {0};
  propane.molFraction = 0.5;
  butane.molFraction = 0.5;
  propane.idealGasRosenbluthWeight = 1.0;
  butane.idealGasRosenbluthWeight = 1.0;

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {50, 50}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.numberOfMoleculesPerComponent[0] + s.numberOfMoleculesPerComponent[1], 100uz);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.externalFieldCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, pair_swap_conventional_na_cl)
{
  ForceField forceField =
      ForceField({{"Na+", false, 22.99, 1.0, 0.0, 11, false}, {"Cl-", false, 35.45, -1.0, 0.0, 17, false}},
                 {{15.0966, 2.65755}, {142.562, 3.51932}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                 12.0, true, false, true);

  std::optional<std::size_t> type_na = forceField.findPseudoAtom("Na+");
  std::optional<std::size_t> type_cl = forceField.findPseudoAtom("Cl-");
  ASSERT_TRUE(type_na.has_value());
  ASSERT_TRUE(type_cl.has_value());

  MCMoveProbabilities probabilities_na = MCMoveProbabilities();
  probabilities_na.setProbability(Move::Types::Translation, 1.0);
  probabilities_na.setProbability(Move::Types::PairSwap, 1.0);

  MCMoveProbabilities probabilities_cl = MCMoveProbabilities();
  probabilities_cl.setProbability(Move::Types::Translation, 1.0);

  Component na = Component::makeIon(forceField, 0, "Na+", type_na.value(), 1.0);
  na.mc_moves_probabilities = probabilities_na;
  na.pairComponentId = 1;
  na.maximumPairDistance = 8.0;
  na.fugacityCoefficient = 1.0;
  na.idealGasRosenbluthWeight = 1.0;

  Component cl = Component::makeIon(forceField, 1, "Cl-", type_cl.value(), -1.0);
  cl.mc_moves_probabilities = probabilities_cl;
  cl.pairComponentId = 0;
  cl.maximumPairDistance = 8.0;
  cl.fugacityCoefficient = 1.0;
  cl.idealGasRosenbluthWeight = 1.0;

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {na, cl}, {}, {4, 4}, 5);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], s.numberOfIntegerMoleculesPerComponent[1]);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-5);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, pair_swap_cbmc_na_cl)
{
  ForceField forceField =
      ForceField({{"Na+", false, 22.99, 1.0, 0.0, 11, false}, {"Cl-", false, 35.45, -1.0, 0.0, 17, false}},
                 {{15.0966, 2.65755}, {142.562, 3.51932}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                  12.0, true, false, true);

  std::optional<std::size_t> type_na = forceField.findPseudoAtom("Na+");
  std::optional<std::size_t> type_cl = forceField.findPseudoAtom("Cl-");
  ASSERT_TRUE(type_na.has_value());
  ASSERT_TRUE(type_cl.has_value());

  MCMoveProbabilities probabilities_na = MCMoveProbabilities();
  probabilities_na.setProbability(Move::Types::Translation, 1.0);
  probabilities_na.setProbability(Move::Types::PairSwapCBMC, 1.0);

  MCMoveProbabilities probabilities_cl = MCMoveProbabilities();
  probabilities_cl.setProbability(Move::Types::Translation, 1.0);

  Component na = Component::makeIon(forceField, 0, "Na+", type_na.value(), 1.0);
  na.mc_moves_probabilities = probabilities_na;
  na.pairComponentId = 1;
  na.maximumPairDistance = 8.0;
  na.fugacityCoefficient = 1.0;
  na.idealGasRosenbluthWeight = 1.0;

  Component cl = Component::makeIon(forceField, 1, "Cl-", type_cl.value(), -1.0);
  cl.mc_moves_probabilities = probabilities_cl;
  cl.pairComponentId = 0;
  cl.maximumPairDistance = 8.0;
  cl.fugacityCoefficient = 1.0;
  cl.idealGasRosenbluthWeight = 1.0;

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {na, cl}, {}, {4, 4}, 5);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], s.numberOfIntegerMoleculesPerComponent[1]);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-5);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, pair_swap_cfcmc_na_cl)
{
  ForceField forceField =
      ForceField({{"Na+", false, 22.99, 1.0, 0.0, 11, false}, {"Cl-", false, 35.45, -1.0, 0.0, 17, false}},
                 {{15.0966, 2.65755}, {142.562, 3.51932}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                  12.0, true, false, true);

  std::optional<std::size_t> type_na = forceField.findPseudoAtom("Na+");
  std::optional<std::size_t> type_cl = forceField.findPseudoAtom("Cl-");
  ASSERT_TRUE(type_na.has_value());
  ASSERT_TRUE(type_cl.has_value());

  MCMoveProbabilities probabilities_na = MCMoveProbabilities();
  probabilities_na.setProbability(Move::Types::Translation, 1.0);
  probabilities_na.setProbability(Move::Types::PairSwapCFCMC, 1.0);

  MCMoveProbabilities probabilities_cl = MCMoveProbabilities();
  probabilities_cl.setProbability(Move::Types::Translation, 1.0);

  Component na = Component::makeIon(forceField, 0, "Na+", type_na.value(), 1.0);
  na.mc_moves_probabilities = probabilities_na;
  na.pairComponentId = 1;
  na.maximumPairDistance = 8.0;
  na.fugacityCoefficient = 1.0;
  na.idealGasRosenbluthWeight = 1.0;

  Component cl = Component::makeIon(forceField, 1, "Cl-", type_cl.value(), -1.0);
  cl.mc_moves_probabilities = probabilities_cl;
  cl.pairComponentId = 0;
  cl.maximumPairDistance = 8.0;
  cl.fugacityCoefficient = 1.0;
  cl.idealGasRosenbluthWeight = 1.0;

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {na, cl}, {}, {4, 4}, 5);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], s.numberOfIntegerMoleculesPerComponent[1]);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-5);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
  }
}

TEST(MC_MUVT_DRIFT, pair_swap_cbcfcmc_na_cl)
{
  ForceField forceField =
      ForceField({{"Na+", false, 22.99, 1.0, 0.0, 11, false}, {"Cl-", false, 35.45, -1.0, 0.0, 17, false}},
                 {{15.0966, 2.65755}, {142.562, 3.51932}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                  12.0, true, false, true);

  std::optional<std::size_t> type_na = forceField.findPseudoAtom("Na+");
  std::optional<std::size_t> type_cl = forceField.findPseudoAtom("Cl-");
  ASSERT_TRUE(type_na.has_value());
  ASSERT_TRUE(type_cl.has_value());

  MCMoveProbabilities probabilities_na = MCMoveProbabilities();
  probabilities_na.setProbability(Move::Types::Translation, 1.0);
  probabilities_na.setProbability(Move::Types::PairSwapCBCFCMC, 1.0);

  MCMoveProbabilities probabilities_cl = MCMoveProbabilities();
  probabilities_cl.setProbability(Move::Types::Translation, 1.0);

  Component na = Component::makeIon(forceField, 0, "Na+", type_na.value(), 1.0);
  na.mc_moves_probabilities = probabilities_na;
  na.pairComponentId = 1;
  na.maximumPairDistance = 8.0;
  na.fugacityCoefficient = 1.0;
  na.idealGasRosenbluthWeight = 1.0;

  Component cl = Component::makeIon(forceField, 1, "Cl-", type_cl.value(), -1.0);
  cl.mc_moves_probabilities = probabilities_cl;
  cl.pairComponentId = 0;
  cl.maximumPairDistance = 8.0;
  cl.fugacityCoefficient = 1.0;
  cl.idealGasRosenbluthWeight = 1.0;

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {na, cl}, {}, {4, 4}, 5);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfIntegerMoleculesPerComponent[0], s.numberOfIntegerMoleculesPerComponent[1]);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-5);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-5);
    EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
  }
}
