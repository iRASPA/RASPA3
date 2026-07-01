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

namespace
{

std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

void expectNoThermodynamicIntegrationDrift(System& system)
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

void runShortMonteCarlo(std::vector<System>& systems, size_t numberOfCycles = 20)
{
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

  for (System& s : mc.systems)
  {
    expectNoThermodynamicIntegrationDrift(s);
  }
}

void runGibbsMonteCarlo(std::vector<System>& systems)
{
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

  for (System& s : mc.systems)
  {
    expectNoThermodynamicIntegrationDrift(s);
  }
}

void enableThermodynamicIntegration(Component& component)
{
  component.lambdaGC.computeDUdlambda = true;
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

Component makeAlkaneFromExample(const ForceField& forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities& probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  Component component = Component(Component::Type::Adsorbate, componentId, forceField, std::string(name),
                                  moleculePath.string(), 5, 21, probabilities, std::nullopt, true);
  component.idealGasRosenbluthWeight = 1.0;
  enableThermodynamicIntegration(component);
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

MCMoveProbabilities makeGibbsConventionalCFCMCCBMCProbabilities()
{
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::GibbsConventionalCFCMCCBMC, 1.0);
  return probabilities;
}

std::tuple<Component, Component, Component> makeMuvtComponents(const ForceField& forceField, Move::Types swapMove,
                                                                bool co2Ti, bool methaneTi)
{
  MCMoveProbabilities probabilities_co2 = MCMoveProbabilities();
  probabilities_co2.setProbability(Move::Types::Translation, 1.0);
  probabilities_co2.setProbability(Move::Types::Rotation, 1.0);
  probabilities_co2.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_co2.setProbability(swapMove, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities_co2, std::nullopt, co2Ti);

  MCMoveProbabilities probabilities_methane = MCMoveProbabilities();
  probabilities_methane.setProbability(Move::Types::Translation, 1.0);
  probabilities_methane.setProbability(Move::Types::Rotation, 1.0);
  probabilities_methane.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities_methane.setProbability(swapMove, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21,
                                probabilities_methane, std::nullopt, methaneTi);

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

  return {co2, methane, water};
}

void runMuvtDrift(Move::Types swapMove, bool co2Ti, bool methaneTi)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  auto [co2, methane, water] = makeMuvtComponents(forceField, swapMove, co2Ti, methaneTi);
  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);
  std::vector<System> systems{system};
  runShortMonteCarlo(systems);
}

ForceField makeZeoliteForceFieldWithGrids()
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.gridPseudoAtomIndices = {3, 4, 5, 6, 7, 8, 9};
  return forceField;
}

void runMuvtGridsDrift(Move::Types swapMove, bool co2Ti, bool methaneTi)
{
  ForceField forceField = makeZeoliteForceFieldWithGrids();
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  auto [co2, methane, water] = makeMuvtComponents(forceField, swapMove, co2Ti, methaneTi);
  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);
  std::vector<System> systems{system};
  runShortMonteCarlo(systems);
}

}  // namespace

// --- muVT (FAU zeolite) ---

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cfcmc_co2)
{
  runMuvtDrift(Move::Types::SwapCFCMC, true, false);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cbcfcbmc_co2)
{
  runMuvtDrift(Move::Types::SwapCBCFCMC, true, false);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cfcmc_methane)
{
  runMuvtDrift(Move::Types::SwapCFCMC, false, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cbcfcbmc_methane)
{
  runMuvtDrift(Move::Types::SwapCBCFCMC, false, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cfcmc_both_components)
{
  runMuvtDrift(Move::Types::SwapCFCMC, true, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_swap_cbcfcbmc_both_components)
{
  runMuvtDrift(Move::Types::SwapCBCFCMC, true, true);
}

// --- muVT with interpolation grids ---

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cfcmc_co2)
{
  runMuvtGridsDrift(Move::Types::SwapCFCMC, true, false);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cbcfcbmc_co2)
{
  runMuvtGridsDrift(Move::Types::SwapCBCFCMC, true, false);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cfcmc_methane)
{
  runMuvtGridsDrift(Move::Types::SwapCFCMC, false, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cbcfcbmc_methane)
{
  runMuvtGridsDrift(Move::Types::SwapCBCFCMC, false, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cfcmc_both_components)
{
  runMuvtGridsDrift(Move::Types::SwapCFCMC, true, true);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, muvt_grids_swap_cbcfcbmc_both_components)
{
  runMuvtGridsDrift(Move::Types::SwapCBCFCMC, true, true);
}

// --- Gibbs ensemble ---

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, gibbs_methane_swap_cbcfcbmc)
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
  runGibbsMonteCarlo(systems);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, gibbs_methane_conventional_cfcmc)
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
  runGibbsMonteCarlo(systems);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, gibbs_methane_conventional_cbcfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsConventionalCFCMCCBMCProbabilities();

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
  runGibbsMonteCarlo(systems);
}

TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, gibbs_methane_ethane_swap_cbcfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsCFCMCProbabilities();

  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  enableThermodynamicIntegration(methane);
  methane.gibbsIdentityChanges = {1};

  Component ethane = makeAlkaneFromExample(forceField, 1, "ethane", probabilities);
  ethane.gibbsIdentityChanges = {0};

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane, ethane}, {}, {30, 30}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {methane, ethane}, {}, {30, 30}, 5, systemProbabilities);

  std::vector<System> systems{systemVapor, systemLiquid};
  runGibbsMonteCarlo(systems);
}

// Flexible multi-bead Gibbs swap with TI: known to drift until TI bookkeeping is fixed for flexible chains.
TEST(MC_THERMODYNAMIC_INTEGRATION_DRIFT, DISABLED_gibbs_binary_propane_butane_swap_cbcfcbmc)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  MCMoveProbabilities probabilities = makeGibbsCFCMCProbabilities();

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
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
  runGibbsMonteCarlo(systems);
}
