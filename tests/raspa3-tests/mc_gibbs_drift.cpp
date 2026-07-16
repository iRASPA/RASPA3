#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import molecule;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import reaction;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_move_types;

import randomnumbers;
import mc_moves_gibbs_swap_cbcfcmc;
import mc_moves_gibbs_conventional_common;
import mc_moves_parallel_tempering_swap;

namespace
{

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
                                const MCMoveProbabilities &probabilities, bool thermodynamicIntegration = false)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  Component component = Component(Component::Type::Adsorbate, componentId, forceField, std::string(name),
                                  moleculePath.string(), 5, 21, probabilities, std::nullopt, thermodynamicIntegration);
  if (thermodynamicIntegration)
  {
    component.idealGasRosenbluthWeight = 1.0;
  }
  return component;
}

MCMoveProbabilities makeGibbsCFCMCProbabilities()
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCBCFCMC, 1.0);
  return probabilities;
}

MCMoveProbabilities makeGibbsConventionalCFCMCProbabilities()
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::GibbsConventionalCFCMC, 1.0);
  return probabilities;
}

MCMoveProbabilities makeGibbsConventionalCBCFCMCProbabilities()
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::GibbsConventionalCBCFCMC, 1.0);
  return probabilities;
}

void enableThermodynamicIntegration(Component &component)
{
  component.lambdaGC.computeDUdlambda = true;
}

void expectNoEnergyDrift(System &system)
{
  RunningEnergy recomputedEnergies = system.computeTotalEnergies();
  RunningEnergy drift = system.runningEnergies - recomputedEnergies;

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
  EXPECT_NEAR(drift.totalDudlambdaVDW(), 0.0, 1e-6);
  EXPECT_NEAR(drift.totalDudlambdaCharge(), 0.0, 1e-6);
  EXPECT_NEAR(drift.totalDudlambdaEwald(), 0.0, 1e-6);
}

// Serial Gibbs CFCMC: active fractional molecule in the first system only (see MonteCarlo::setup).
void setupSerialGibbsCFCMC(std::vector<System> &systems)
{
  MonteCarlo mc = MonteCarlo({0, 0, 0, 0, 1000, 10000, 5000, 5000}, systems, {}, 5, false);
  mc.setup();
  systems = mc.systems;
  for (System &system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
  }
}

void setGibbsLambda(System &system, std::size_t componentId, std::size_t bin)
{
  Component &component = system.components[componentId];
  component.lambdaGibbs.setCurrentBin(bin);
  const double lambda = component.lambdaGibbs.lambdaValue();
  const std::size_t fractionalIndex =
      system.indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, componentId);
  for (Atom &atom : system.spanOfMolecule(componentId, fractionalIndex))
  {
    atom.setScalingToFractional(lambda, component.lambdaGibbs.dUdlambdaGroupId);
  }
  system.runningEnergies = system.computeTotalEnergies();
  system.trialEik = system.storedEik;
}

void favorAllBinsExceptCurrent(Component &component)
{
  const std::size_t oldBin = component.lambdaGibbs.currentBin;
  std::fill(component.lambdaGibbs.biasFactor.begin(), component.lambdaGibbs.biasFactor.end(), 100.0);
  component.lambdaGibbs.biasFactor[oldBin] = 0.0;
}

void rejectAllBinsExceptCurrent(Component &component)
{
  const std::size_t oldBin = component.lambdaGibbs.currentBin;
  std::fill(component.lambdaGibbs.biasFactor.begin(), component.lambdaGibbs.biasFactor.end(), -100.0);
  component.lambdaGibbs.biasFactor[oldBin] = 0.0;
}

}  // namespace

TEST(MC_GIBBS_DRIFT, translation_rotation_reinsertion_volume)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities, std::nullopt, false);

  Component methane =
      Component(forceField, "methane", 190.564, 45599200, 0.01142,
                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                false);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();

  System systemVapor = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2, methane, water}, {}, {50, 55, 80},
                              5, systemProbabilities);
  System systemLiquid = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2, methane, water}, {}, {50, 55, 80},
                               5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, binary_mixture_propane_butane)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsIdentityChangeCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
  propane.gibbsIdentityChanges = {1};
  butane.gibbsIdentityChanges = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {50, 50}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {50, 50}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, gibbs_mixture_methane_ethane_co2)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();

  MCMoveProbabilities probabilities_rigid = MCMoveProbabilities();
  probabilities_rigid.setProbability(Move::Types::Translation, 1.0);
  probabilities_rigid.setProbability(Move::Types::Rotation, 1.0);
  probabilities_rigid.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_rigid.setProbability(Move::Types::GibbsSwapCBMC, 1.0);
  probabilities_rigid.setProbability(Move::Types::GibbsIdentityChangeCBMC, 1.0);

  MCMoveProbabilities probabilities_ethane = probabilities_rigid;

  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities_rigid;
  methane.gibbsIdentityChanges = {1, 2};

  Component ethane = makeAlkaneFromExample(forceField, 1, "ethane", probabilities_ethane);
  ethane.gibbsIdentityChanges = {0, 2};

  Component co2 = Component::makeCO2(forceField, 2, true);
  co2.mc_moves_probabilities = probabilities_rigid;
  co2.gibbsIdentityChanges = {0, 1};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane, ethane, co2}, {}, {50, 50, 50}, 5,
                              systemProbabilities);
  System systemLiquid = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane, ethane, co2}, {}, {50, 50, 50},
                               5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, binary_mixture_propane_butane_tail_corrections)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsIdentityChangeCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
  propane.gibbsIdentityChanges = {1};
  butane.gibbsIdentityChanges = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {50, 50}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {50, 50}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, binary_mixture_propane_butane_cfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsCFCMCProbabilities();

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities, true);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities, true);
  propane.gibbsIdentityChanges = {1};
  butane.gibbsIdentityChanges = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {20, 20}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {20, 20}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{20};
  size_t numberOfInitializationCycles{5};
  size_t numberOfEquilibrationCycles{5};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, binary_mixture_propane_butane_cfcbmc_gibbs_only_single_accept)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::GibbsSwapCBCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities, true);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities, true);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();

  std::vector<System> systems{
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {10, 10}, 5, systemProbabilities),
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {10, 10}, 5, systemProbabilities)};
  setupSerialGibbsCFCMC(systems);
  System &systemVapor = systems[0];
  System &systemLiquid = systems[1];
  expectNoEnergyDrift(systemVapor);
  expectNoEnergyDrift(systemLiquid);

  RandomNumber random(42);
  std::size_t fractionalMoleculeSystem = 0;
  for (int attempt = 0; attempt < 500; ++attempt)
  {
    for (std::size_t component = 0; component < 2; ++component)
    {
      System &systemA = systemVapor.containsTheFractionalMolecule ? systemVapor : systemLiquid;
      System &systemB = systemVapor.containsTheFractionalMolecule ? systemLiquid : systemVapor;

      RunningEnergy beforeA = systemA.computeTotalEnergies();
      RunningEnergy beforeB = systemB.computeTotalEnergies();

      std::optional<std::pair<RunningEnergy, RunningEnergy>> delta =
          MC_Moves::GibbsSwapMove_CBCFCMC(random, systemA, systemB, component, fractionalMoleculeSystem);
      if (!delta)
      {
        continue;
      }

      systemA.runningEnergies += delta->first;
      systemB.runningEnergies += delta->second;

      RunningEnergy afterA = systemA.computeTotalEnergies();
      RunningEnergy afterB = systemB.computeTotalEnergies();
      RunningEnergy actualDeltaA = afterA - beforeA;
      RunningEnergy actualDeltaB = afterB - beforeB;

      EXPECT_NEAR(delta->first.potentialEnergy(), actualDeltaA.potentialEnergy(), 1e-6)
          << "attempt=" << attempt << " component=" << component;
      EXPECT_NEAR(delta->first.moleculeMoleculeVDW, actualDeltaA.moleculeMoleculeVDW, 1e-6);
      EXPECT_NEAR(delta->first.totalDudlambdaVDW(), actualDeltaA.totalDudlambdaVDW(), 1e-6);
      EXPECT_NEAR(delta->second.potentialEnergy(), actualDeltaB.potentialEnergy(), 1e-6);
      expectNoEnergyDrift(systemA);
      expectNoEnergyDrift(systemB);
      return;
    }
  }
  FAIL() << "No accepted GibbsSwapCBCFCMC move in 500 attempts";
}

TEST(MC_GIBBS_DRIFT, binary_mixture_propane_butane_cfcbmc_gibbs_only_many_accepts)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::GibbsSwapCBCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities, true);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities, true);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();

  std::vector<System> systems{
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {10, 10}, 5, systemProbabilities),
      System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, {10, 10}, 5, systemProbabilities)};
  setupSerialGibbsCFCMC(systems);
  System &systemVapor = systems[0];
  System &systemLiquid = systems[1];

  RandomNumber random(42);
  std::size_t fractionalMoleculeSystem = 0;
  std::size_t acceptCount = 0;
  for (int attempt = 0; attempt < 2000; ++attempt)
  {
    System &systemA = systemVapor.containsTheFractionalMolecule ? systemVapor : systemLiquid;
    System &systemB = systemVapor.containsTheFractionalMolecule ? systemLiquid : systemVapor;
    const std::size_t component = systemA.randomComponent(random);

    RunningEnergy beforeA = systemA.computeTotalEnergies();
    RunningEnergy beforeB = systemB.computeTotalEnergies();

    std::optional<std::pair<RunningEnergy, RunningEnergy>> delta =
        MC_Moves::GibbsSwapMove_CBCFCMC(random, systemA, systemB, component, fractionalMoleculeSystem);
    if (!delta)
    {
      continue;
    }

    systemA.runningEnergies += delta->first;
    systemB.runningEnergies += delta->second;
    ++acceptCount;

    RunningEnergy afterA = systemA.computeTotalEnergies();
    RunningEnergy afterB = systemB.computeTotalEnergies();
    SCOPED_TRACE("accept=" + std::to_string(acceptCount) + " component=" + std::to_string(component));
    ASSERT_NEAR(delta->first.moleculeMoleculeVDW, (afterA - beforeA).moleculeMoleculeVDW, 1e-6);
    ASSERT_NEAR(delta->first.totalDudlambdaVDW(), (afterA - beforeA).totalDudlambdaVDW(), 1e-6);
    ASSERT_NEAR(delta->second.moleculeMoleculeVDW, (afterB - beforeB).moleculeMoleculeVDW, 1e-6);
    ASSERT_NEAR(delta->second.totalDudlambdaVDW(), (afterB - beforeB).totalDudlambdaVDW(), 1e-6);

    // the inactive box has all its fractional molecules switched off with groupId 0,
    // so it must have no dUdlambda contribution at all
    const RunningEnergy &inactiveAfter = systemA.containsTheFractionalMolecule ? afterB : afterA;
    ASSERT_NEAR(inactiveAfter.totalDudlambdaVDW(), 0.0, 1e-8);

    expectNoEnergyDrift(systemA);
    expectNoEnergyDrift(systemB);
  }
  EXPECT_GT(acceptCount, 0u);
}

TEST(MC_GIBBS_DRIFT, methane_gibbs_cfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsCFCMCProbabilities();

  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  enableThermodynamicIntegration(methane);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, methane_gibbs_conventional_cfcmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsConventionalCFCMCProbabilities();

  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  enableThermodynamicIntegration(methane);
  methane.lambdaGibbs.computeDUdlambda = true;

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);
  EXPECT_TRUE(systemVapor.usesGibbsConventionalCFCMC());

  System systemLiquid =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, methane_gibbs_conventional_cfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsConventionalCBCFCMCProbabilities();

  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  enableThermodynamicIntegration(methane);
  methane.lambdaGibbs.computeDUdlambda = true;

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);
  EXPECT_TRUE(systemVapor.usesGibbsConventionalCFCMC());

  System systemLiquid =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {40}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  size_t numberOfProductionCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
  size_t printEvery{1000};
  size_t writeBinaryRestartEvery{10000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t numberOfBlocks{5};
  bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery}, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
  }
}

TEST(MC_GIBBS_DRIFT, conventional_mixed_insert_and_lambda_boundary)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = makeGibbsConventionalCFCMCProbabilities();
  methane.lambdaGibbs.computeDUdlambda = true;

  const SimulationBox box(30.0, 30.0, 30.0);
  std::vector<System> systems{
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {12}, 5),
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {12}, 5)};
  setupSerialGibbsCFCMC(systems);

  setGibbsLambda(systems[0], 0, 19);
  setGibbsLambda(systems[1], 0, 10);
  favorAllBinsExceptCurrent(systems[0].components[0]);
  favorAllBinsExceptCurrent(systems[1].components[0]);

  const std::size_t integersA = systems[0].numberOfIntegerMoleculesPerComponent[0];
  const std::size_t integersB = systems[1].numberOfIntegerMoleculesPerComponent[0];
  const RunningEnergy beforeA = systems[0].computeTotalEnergies();
  const RunningEnergy beforeB = systems[1].computeTotalEnergies();

  RandomNumber random(42);
  const auto delta =
      MC_Moves::GibbsConventionalCommon::gibbsConventionalMove(random, systems[0], systems[1], 0, false);
  ASSERT_TRUE(delta.has_value());
  systems[0].runningEnergies += delta->first;
  systems[1].runningEnergies += delta->second;

  EXPECT_EQ(systems[0].numberOfIntegerMoleculesPerComponent[0], integersA + 1);
  EXPECT_EQ(systems[1].numberOfIntegerMoleculesPerComponent[0], integersB);
  EXPECT_NEAR(delta->first.potentialEnergy(),
              (systems[0].computeTotalEnergies() - beforeA).potentialEnergy(), 1e-6);
  EXPECT_NEAR(delta->second.potentialEnergy(),
              (systems[1].computeTotalEnergies() - beforeB).potentialEnergy(), 1e-6);
  expectNoEnergyDrift(systems[0]);
  expectNoEnergyDrift(systems[1]);
}

TEST(MC_GIBBS_DRIFT, cbmc_mixed_lambda_and_delete_boundary)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = makeGibbsConventionalCBCFCMCProbabilities();
  methane.lambdaGibbs.computeDUdlambda = true;
  methane.idealGasRosenbluthWeight = 1.0;

  const SimulationBox box(30.0, 30.0, 30.0);
  std::vector<System> systems{
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {12}, 5),
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane}, {}, {12}, 5)};
  setupSerialGibbsCFCMC(systems);

  setGibbsLambda(systems[0], 0, 10);
  setGibbsLambda(systems[1], 0, 2);
  favorAllBinsExceptCurrent(systems[0].components[0]);
  favorAllBinsExceptCurrent(systems[1].components[0]);

  const std::size_t integersA = systems[0].numberOfIntegerMoleculesPerComponent[0];
  const std::size_t integersB = systems[1].numberOfIntegerMoleculesPerComponent[0];
  const RunningEnergy beforeA = systems[0].computeTotalEnergies();
  const RunningEnergy beforeB = systems[1].computeTotalEnergies();

  RandomNumber random(42);
  const auto delta =
      MC_Moves::GibbsConventionalCommon::gibbsConventionalMove(random, systems[0], systems[1], 0, true);
  ASSERT_TRUE(delta.has_value());
  systems[0].runningEnergies += delta->first;
  systems[1].runningEnergies += delta->second;

  EXPECT_EQ(systems[0].numberOfIntegerMoleculesPerComponent[0], integersA);
  EXPECT_EQ(systems[1].numberOfIntegerMoleculesPerComponent[0], integersB - 1);
  EXPECT_NEAR(delta->first.potentialEnergy(),
              (systems[0].computeTotalEnergies() - beforeA).potentialEnergy(), 1e-6);
  EXPECT_NEAR(delta->second.potentialEnergy(),
              (systems[1].computeTotalEnergies() - beforeB).potentialEnergy(), 1e-6);
  expectNoEnergyDrift(systems[0]);
  expectNoEnergyDrift(systems[1]);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_swaps_complete_configuration_only)
{
  const ForceField forceField({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{0.0, 1.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particle(forceField, "particle", 1.0, 1.0, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                           {}, {}, 5, 21);
  const SimulationBox boxA(20.0, 20.0, 20.0);
  const SimulationBox boxB(25.0, 25.0, 25.0);
  System systemA(forceField, boxA, false, 250.0, 1e4, 1.0, {}, {particle}, {}, {2}, 5);
  System systemB(forceField, boxB, false, 450.0, 1e4, 1.0, {}, {particle}, {}, {4}, 7);
  systemA.mc_moves_statistics.setMaxChange(Move::Types::ParallelTempering, 0.123);
  systemB.mc_moves_statistics.setMaxChange(Move::Types::ParallelTempering, 0.456);
  systemA.runningEnergies = systemA.computeTotalEnergies();
  systemA.trialEik = systemA.storedEik;
  systemB.runningEnergies = systemB.computeTotalEnergies();
  systemB.trialEik = systemB.storedEik;

  const double temperatureA = systemA.temperature;
  const double betaA = systemA.beta;
  const double temperatureB = systemB.temperature;
  const double betaB = systemB.beta;
  const MCMoveProbabilities probabilitiesA = systemA.mc_moves_probabilities;
  const MCMoveProbabilities probabilitiesB = systemB.mc_moves_probabilities;
  const std::size_t energyBlocksA = systemA.averageEnergies.numberOfBlocks;
  const std::size_t energyBlocksB = systemB.averageEnergies.numberOfBlocks;
  const std::vector<Atom> atomsA = systemA.atomData;
  const std::vector<Atom> atomsB = systemB.atomData;
  const std::vector<Molecule> moleculesA = systemA.moleculeData;
  const std::vector<Molecule> moleculesB = systemB.moleculeData;

  RandomNumber random(7);
  const auto energies = MC_Moves::ParallelTemperingSwap(random, systemA, systemB);
  ASSERT_TRUE(energies.has_value());

  EXPECT_EQ(systemA.temperature, temperatureA);
  EXPECT_EQ(systemA.beta, betaA);
  EXPECT_EQ(systemB.temperature, temperatureB);
  EXPECT_EQ(systemB.beta, betaB);
  EXPECT_EQ(systemA.mc_moves_probabilities, probabilitiesA);
  EXPECT_EQ(systemB.mc_moves_probabilities, probabilitiesB);
  EXPECT_DOUBLE_EQ(systemA.mc_moves_statistics.getMaxChange(Move::Types::ParallelTempering), 0.123);
  EXPECT_DOUBLE_EQ(systemB.mc_moves_statistics.getMaxChange(Move::Types::ParallelTempering), 0.456);
  EXPECT_EQ(systemA.averageEnergies.numberOfBlocks, energyBlocksA);
  EXPECT_EQ(systemB.averageEnergies.numberOfBlocks, energyBlocksB);
  ASSERT_EQ(systemA.atomData.size(), atomsB.size());
  ASSERT_EQ(systemB.atomData.size(), atomsA.size());
  EXPECT_EQ(systemA.atomData.front().position, atomsB.front().position);
  EXPECT_EQ(systemB.atomData.front().position, atomsA.front().position);
  ASSERT_EQ(systemA.moleculeData.size(), moleculesB.size());
  ASSERT_EQ(systemB.moleculeData.size(), moleculesA.size());
  EXPECT_EQ(systemA.moleculeData.front().centerOfMassPosition, moleculesB.front().centerOfMassPosition);
  EXPECT_EQ(systemB.moleculeData.front().centerOfMassPosition, moleculesA.front().centerOfMassPosition);
  EXPECT_EQ(systemA.moleculeData.front().numberOfAtoms, moleculesB.front().numberOfAtoms);
  EXPECT_EQ(systemB.moleculeData.front().numberOfAtoms, moleculesA.front().numberOfAtoms);
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], 4u);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], 2u);
  EXPECT_EQ(systemA.simulationBox, boxB);
  EXPECT_EQ(systemB.simulationBox, boxA);
  EXPECT_EQ(systemA.electricField.size(), systemA.atomData.size());
  EXPECT_EQ(systemB.electricField.size(), systemB.atomData.size());
  expectNoEnergyDrift(systemA);
  expectNoEnergyDrift(systemB);

  const auto reverseEnergies = MC_Moves::ParallelTemperingSwap(random, systemA, systemB);
  ASSERT_TRUE(reverseEnergies.has_value());
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], 2u);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], 4u);
  EXPECT_EQ(systemA.simulationBox, boxA);
  EXPECT_EQ(systemB.simulationBox, boxB);
  EXPECT_EQ(systemA.temperature, temperatureA);
  EXPECT_EQ(systemB.temperature, temperatureB);
  EXPECT_DOUBLE_EQ(systemA.mc_moves_statistics.getMaxChange(Move::Types::ParallelTempering), 0.123);
  EXPECT_DOUBLE_EQ(systemB.mc_moves_statistics.getMaxChange(Move::Types::ParallelTempering), 0.456);
  expectNoEnergyDrift(systemA);
  expectNoEnergyDrift(systemB);
}

TEST(MC_GIBBS_DRIFT, charged_ewald_mixed_boundary_chains_and_rolls_back)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = makeGibbsConventionalCFCMCProbabilities();
  co2.lambdaGibbs.computeDUdlambda = true;

  const SimulationBox box(30.0, 30.0, 30.0);
  std::vector<System> systems{
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2}, {}, {8}, 5),
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2}, {}, {8}, 5)};
  setupSerialGibbsCFCMC(systems);
  setGibbsLambda(systems[0], 0, 19);
  setGibbsLambda(systems[1], 0, 10);
  favorAllBinsExceptCurrent(systems[0].components[0]);
  favorAllBinsExceptCurrent(systems[1].components[0]);

  const RunningEnergy beforeA = systems[0].computeTotalEnergies();
  systems[0].trialEik = systems[0].storedEik;
  const RunningEnergy beforeB = systems[1].computeTotalEnergies();
  systems[1].trialEik = systems[1].storedEik;
  RandomNumber acceptedRandom(42);
  const auto delta =
      MC_Moves::GibbsConventionalCommon::gibbsConventionalMove(acceptedRandom, systems[0], systems[1], 0, false);
  ASSERT_TRUE(delta.has_value());
  systems[0].runningEnergies += delta->first;
  systems[1].runningEnergies += delta->second;

  const RunningEnergy afterA = systems[0].computeTotalEnergies();
  const RunningEnergy afterB = systems[1].computeTotalEnergies();
  EXPECT_NEAR(delta->first.ewald_fourier, (afterA - beforeA).ewald_fourier, 1e-6);
  EXPECT_NEAR(delta->first.ewald_self, (afterA - beforeA).ewald_self, 1e-6);
  EXPECT_NEAR(delta->first.ewald_exclusion, (afterA - beforeA).ewald_exclusion, 1e-6);
  EXPECT_NEAR(delta->second.ewald_fourier, (afterB - beforeB).ewald_fourier, 1e-6);
  expectNoEnergyDrift(systems[0]);
  expectNoEnergyDrift(systems[1]);

  setGibbsLambda(systems[0], 0, 19);
  setGibbsLambda(systems[1], 0, 10);
  rejectAllBinsExceptCurrent(systems[0].components[0]);
  rejectAllBinsExceptCurrent(systems[1].components[0]);
  const auto storedA = systems[0].storedEik;
  const auto totalA = systems[0].trialEik;
  const auto storedB = systems[1].storedEik;
  const auto totalB = systems[1].trialEik;
  const std::size_t countA = systems[0].numberOfIntegerMoleculesPerComponent[0];
  const std::size_t countB = systems[1].numberOfIntegerMoleculesPerComponent[0];
  const double scalingA = systems[0].spanOfMolecule(
      0, systems[0].indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, 0))[0].scalingCoulomb;
  const double scalingB = systems[1].spanOfMolecule(
      0, systems[1].indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, 0))[0].scalingCoulomb;

  RandomNumber rejectedRandom(42);
  EXPECT_FALSE(MC_Moves::GibbsConventionalCommon::gibbsConventionalMove(
                   rejectedRandom, systems[0], systems[1], 0, false)
                   .has_value());
  EXPECT_EQ(systems[0].storedEik, storedA);
  EXPECT_EQ(systems[0].trialEik, totalA);
  EXPECT_EQ(systems[1].storedEik, storedB);
  EXPECT_EQ(systems[1].trialEik, totalB);
  EXPECT_EQ(systems[0].numberOfIntegerMoleculesPerComponent[0], countA);
  EXPECT_EQ(systems[1].numberOfIntegerMoleculesPerComponent[0], countB);
  EXPECT_EQ(systems[0].spanOfMolecule(
                0, systems[0].indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, 0))[0]
                .scalingCoulomb,
            scalingA);
  EXPECT_EQ(systems[1].spanOfMolecule(
                0, systems[1].indexOfFractionalMoleculeForMove(Move::Types::GibbsConventionalCFCMC, 0))[0]
                .scalingCoulomb,
            scalingB);
  expectNoEnergyDrift(systems[0]);
  expectNoEnergyDrift(systems[1]);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_rejection_preserves_configuration)
{
  const ForceField forceField({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{0.0, 1.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particle(forceField, "particle", 1.0, 1.0, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System systemA(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e9, 1.0, {}, {particle}, {}, {2}, 5);
  System systemB(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1.0, 1.0, {}, {particle}, {}, {4}, 5);
  systemA.runningEnergies = systemA.computeTotalEnergies();
  systemA.trialEik = systemA.storedEik;
  systemB.runningEnergies = systemB.computeTotalEnergies();
  systemB.trialEik = systemB.storedEik;
  const auto atomsA = systemA.atomData;
  const auto atomsB = systemB.atomData;
  const SimulationBox boxA = systemA.simulationBox;
  const SimulationBox boxB = systemB.simulationBox;

  RandomNumber random(7);
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  ASSERT_EQ(systemA.atomData.size(), atomsA.size());
  ASSERT_EQ(systemB.atomData.size(), atomsB.size());
  EXPECT_EQ(systemA.atomData.front().position, atomsA.front().position);
  EXPECT_EQ(systemB.atomData.front().position, atomsB.front().position);
  EXPECT_EQ(systemA.simulationBox, boxA);
  EXPECT_EQ(systemB.simulationBox, boxB);
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], 2u);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], 4u);
  expectNoEnergyDrift(systemA);
  expectNoEnergyDrift(systemB);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_rejects_different_hamiltonians_before_random_draw)
{
  const ForceField forceFieldA({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{0.0, 1.0}},
                               ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const ForceField forceFieldB({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{10.0, 1.0}},
                               ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particleA(forceFieldA, "particle", 1.0, 1.0, 0.0,
                            {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  const Component particleB(forceFieldB, "particle", 1.0, 1.0, 0.0,
                            {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System systemA(forceFieldA, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {particleA}, {}, {2}, 5);
  System systemB(forceFieldB, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {particleB}, {}, {2}, 5);
  const std::size_t countA = systemA.numberOfIntegerMoleculesPerComponent[0];
  const std::size_t countB = systemB.numberOfIntegerMoleculesPerComponent[0];

  RandomNumber random(9);
  const std::size_t drawsBefore = random.count;
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  EXPECT_EQ(random.count, drawsBefore);
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], countA);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], countB);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_keeps_framework_bundle_fixed)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  const Framework framework = Framework::makeFAU(forceField);
  Component methane = Component::makeMethane(forceField, 0);
  System systemA(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {methane}, {}, {2}, 5);
  System systemB(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {methane}, {}, {4}, 5);
  systemA.runningEnergies = systemA.computeTotalEnergies();
  systemA.trialEik = systemA.storedEik;
  systemB.runningEnergies = systemB.computeTotalEnergies();
  systemB.trialEik = systemB.storedEik;
  const std::vector<Atom> frameworkAtomsA(systemA.spanOfFrameworkAtoms().begin(), systemA.spanOfFrameworkAtoms().end());
  const std::vector<Atom> frameworkAtomsB(systemB.spanOfFrameworkAtoms().begin(), systemB.spanOfFrameworkAtoms().end());
  const double3 mobilePositionA = systemA.spanOfMolecule(0, 0)[0].position;
  const double3 mobilePositionB = systemB.spanOfMolecule(0, 0)[0].position;

  RandomNumber random(11);
  ASSERT_TRUE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  ASSERT_EQ(systemA.spanOfFrameworkAtoms().size(), frameworkAtomsA.size());
  ASSERT_EQ(systemB.spanOfFrameworkAtoms().size(), frameworkAtomsB.size());
  EXPECT_EQ(systemA.spanOfFrameworkAtoms()[0].position, frameworkAtomsA[0].position);
  EXPECT_EQ(systemB.spanOfFrameworkAtoms()[0].position, frameworkAtomsB[0].position);
  EXPECT_EQ(systemA.spanOfMolecule(0, 0)[0].position, mobilePositionB);
  EXPECT_EQ(systemB.spanOfMolecule(0, 0)[0].position, mobilePositionA);
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], 4u);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], 2u);
  expectNoEnergyDrift(systemA);
  expectNoEnergyDrift(systemB);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_rejects_flexible_components_before_random_draw)
{
  const ForceField forceField = makeAlkaneForceField();
  const Component propane = makeAlkaneFromExample(forceField, 0, "propane", MCMoveProbabilities{});
  ASSERT_FALSE(propane.rigid);
  System systemA(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {propane}, {}, {2}, 5);
  System systemB(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {propane}, {}, {2}, 5);

  RandomNumber random(15);
  const std::size_t drawsBefore = random.count;
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  EXPECT_EQ(random.count, drawsBefore);
  EXPECT_EQ(systemA.numberOfIntegerMoleculesPerComponent[0], 2u);
  EXPECT_EQ(systemB.numberOfIntegerMoleculesPerComponent[0], 2u);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_rejects_fractional_state_before_random_draw)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = makeGibbsConventionalCFCMCProbabilities();
  System systemA(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {4}, 5);
  System systemB(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {4}, 5);
  ASSERT_GT(systemA.numberOfFractionalMoleculesPerComponent[0], 0u);
  ASSERT_GT(systemB.numberOfFractionalMoleculesPerComponent[0], 0u);

  RandomNumber random(16);
  const std::size_t drawsBefore = random.count;
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  EXPECT_EQ(random.count, drawsBefore);
}

TEST(MC_GIBBS_DRIFT, parallel_tempering_rejects_reaction_definitions_before_random_draw)
{
  const ForceField forceField({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{0.0, 1.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particle(forceField, "particle", 1.0, 1.0, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System systemA(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {particle}, {}, {2}, 5);
  System systemB(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {particle}, {}, {2}, 5);
  systemA.reactions.list.emplace_back(0, std::vector<std::size_t>{1}, std::vector<std::size_t>{1});
  systemB.reactions.list.emplace_back(0, std::vector<std::size_t>{1}, std::vector<std::size_t>{1});

  RandomNumber random(17);
  const std::size_t drawsBefore = random.count;
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value());
  EXPECT_EQ(random.count, drawsBefore);
}
