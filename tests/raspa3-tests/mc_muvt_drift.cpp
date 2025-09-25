#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <iostream>
#include <optional>
#include <print>
#include <ranges>
#include <span>
#include <tuple>
#include <vector>

import int3;
import double3;
import double3x3;
import randomnumbers;
import factory;
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
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::Swap, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::Swap, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCBMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCBMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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

TEST(MC_MUVT_DRIFT, insertionCFCMC_thermodynamic_integration_co2)
{
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, true);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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

TEST(MC_MUVT_DRIFT, insertionCBCFCMC_thermodynamic_integration_co2)
{
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, true);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, false);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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

TEST(MC_MUVT_DRIFT, insertionCFCMC_thermodynamic_integration_methane)
{
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, true);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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

TEST(MC_MUVT_DRIFT, insertionCBCFCMC_thermodynamic_integration_methane)
{
  const ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Framework f = TestFactories::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(MoveTypes::Translation, 1.0);
  probabilities_co2.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_co2.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component co2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, false);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(MoveTypes::Translation, 1.0);
  probabilities_methane.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_methane.setProbability(MoveTypes::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(MoveTypes::SwapCBCFCMC, 1.0);

  Component methane = Component(1, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, true);

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(MoveTypes::Translation, 1.0);
  probabilities_water.setProbability(MoveTypes::Rotation, 1.0);
  probabilities_water.setProbability(MoveTypes::ReinsertionCBMC, 1.0);

  Component water = Component(
      2, forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  std::vector<System> systems{system};
  size_t numberOfCycles{1000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  RandomNumber randomSeed(std::nullopt);

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, randomSeed,
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
