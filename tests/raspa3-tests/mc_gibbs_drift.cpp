#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
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

import randomnumbers;
import mc_moves_gibbs_swap_cbcfcmc;

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
  EXPECT_NEAR(drift.dudlambdaVDW, 0.0, 1e-6);
  EXPECT_NEAR(drift.dudlambdaCharge, 0.0, 1e-6);
  EXPECT_NEAR(drift.dudlambdaEwald, 0.0, 1e-6);
}

// Serial Gibbs CFCMC: active fractional molecule in the first system only (see MonteCarlo::setup).
void setupSerialGibbsCFCMC(std::vector<System> &systems)
{
  MonteCarlo mc = MonteCarlo(0, 0, 0, 1000, 10000, 5000, 5000, systems, {}, 5, false);
  mc.setup();
  systems = mc.systems;
  for (System &system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
  }
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
      EXPECT_NEAR(delta->first.dudlambdaVDW, actualDeltaA.dudlambdaVDW, 1e-6);
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
    ASSERT_NEAR(delta->first.dudlambdaVDW, (afterA - beforeA).dudlambdaVDW, 1e-6);
    ASSERT_NEAR(delta->second.moleculeMoleculeVDW, (afterB - beforeB).moleculeMoleculeVDW, 1e-6);
    ASSERT_NEAR(delta->second.dudlambdaVDW, (afterB - beforeB).dudlambdaVDW, 1e-6);

    // the inactive box has all its fractional molecules switched off with groupId 0,
    // so it must have no dUdlambda contribution at all
    const RunningEnergy &inactiveAfter = systemA.containsTheFractionalMolecule ? afterB : afterA;
    ASSERT_NEAR(inactiveAfter.dudlambdaVDW, 0.0, 1e-8);

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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
  size_t numberOfCycles{2000};
  size_t numberOfInitializationCycles{500};
  size_t numberOfEquilibrationCycles{1000};
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
    expectNoEnergyDrift(s);
  }
}
