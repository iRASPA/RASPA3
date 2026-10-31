#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import archive;
import randomnumbers;
import atom;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import monte_carlo;
import monte_carlo_transition_matrix;
import simulation_schedule;
import mc_moves_move_types;
import mc_moves_probabilities;

// Driver-level bookkeeping of the Monte Carlo drivers: the per-step sampling with a single system, the
// restart archive, and the first-run preparation of the transition-matrix driver.

namespace
{

// A single Lennard-Jones component with CFCMC swaps (so that lambda sampling is active) in an empty box.
System makeMethaneCFCMCSystem(std::size_t numberOfMolecules)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  Component methane(forceField, "methane", 190.564, 45599200, 0.01142,
                    {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                    false);
  methane.lambdaGC.computeDUdlambda = true;
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                {numberOfMolecules}, 5);
}

SimulationSchedule makeSchedule(std::size_t production, std::size_t initialization, std::size_t equilibration)
{
  return {production, 0, initialization, equilibration, 1000, 10000, 5000, 5000};
}

}  // namespace

// With a single system the "second" system of the adjacent pair aliases the first: the lambda bookkeeping of
// that system must be sampled once per step (as it is for every system of a multi-system run), not twice.
TEST(MC_DRIVER_BOOKKEEPING, single_system_is_sampled_once_per_step)
{
  std::vector<System> systems{makeMethaneCFCMCSystem(20)};
  MonteCarlo mc(makeSchedule(1, 0, 1), systems, 7uz, 5, false);
  mc.setup();
  for (System& system : mc.systems) system.runningEnergies = system.computeTotalEnergies();

  // one component: max(number of molecules, 20) * 1 steps per cycle (the count is taken before the cycle,
  // the CFCMC swaps may change it)
  const std::size_t expectedSteps = std::max(mc.systems[0].numberOfMolecules(), 20uz);
  mc.simulationStage = MonteCarlo::SimulationStage::Equilibration;
  mc.performCycle();
  EXPECT_EQ(mc.systems[0].components[0].lambdaGC.occupancyTotal, expectedSteps);

  // the transition-matrix driver takes the same path
  std::vector<System> tmmcSystems{makeMethaneCFCMCSystem(20)};
  RandomNumber random(7);
  MonteCarloTransitionMatrix tmmc(makeSchedule(1, 0, 1), tmmcSystems, random, 5);
  for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
  const std::size_t expectedTMMCSteps = std::max(tmmc.systems[0].numberOfMolecules(), 20uz);
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
  tmmc.performCycle();
  EXPECT_EQ(tmmc.systems[0].components[0].lambdaGC.occupancyTotal, expectedTMMCSteps);
}

// The overall cycle number indexes the property-evolution samples; it has to survive a binary restart.
TEST(MC_DRIVER_BOOKKEEPING, archive_restores_the_absolute_cycle_number)
{
  std::vector<System> systems{makeMethaneCFCMCSystem(10)};
  MonteCarlo original(makeSchedule(4, 3, 2), systems, 11uz, 5, false);
  original.run();
  ASSERT_EQ(original.absoluteCurrentCycle, 9uz);

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_monte_carlo_bookkeeping.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  MonteCarlo restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.absoluteCurrentCycle, original.absoluteCurrentCycle);
  EXPECT_EQ(restored.currentCycle, original.currentCycle);
  EXPECT_EQ(restored.simulationStage, original.simulationStage);
}

// The transition-matrix driver prepares the systems like the plain Monte Carlo driver: the fractional molecule of a
// Gibbs CFCMC swap is active in the first box only (its copy in the other box must carry groupId 0, otherwise it
// contributes a spurious dU/dlambda), and the requested framework interpolation grids are built.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_setup_prepares_fractional_molecules_and_grids)
{
  // charges on, so that the Ewald-real grid is built as well; only the VDW grid size is given explicitly and the
  // Coulomb grid size must then follow from its spacing (not read from the unset optional)
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const std::size_t methaneType = forceField.findPseudoAtom("CH4").value();
  forceField.settings.gridPseudoAtomIndices = {methaneType};
  forceField.settings.numberOfVDWGridPoints = uint3(8, 8, 8);
  forceField.settings.spacingCoulombGrid = 1.5;
  forceField.settings.numberOfGridTestPoints = 10;

  Framework itq29 = Framework::makeITQ29(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  methane.lambdaGC.computeDUdlambda = true;

  std::vector<System> systems{
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {itq29}, {methane}, {}, {4}, 5),
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {itq29}, {methane}, {}, {6}, 5)};
  ASSERT_EQ(systems[0].numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[0], 1uz);

  // the driver writes its output files relative to the working directory
  const std::filesystem::path previous = std::filesystem::current_path();
  const std::filesystem::path scratch = std::filesystem::temp_directory_path() / "raspa3_tmmc_setup_test";
  std::filesystem::create_directories(scratch);
  std::filesystem::current_path(scratch);
  {
    RandomNumber random(3);
    MonteCarloTransitionMatrix tmmc(makeSchedule(0, 0, 0), systems, random, 5);
    tmmc.setup();

    EXPECT_TRUE(tmmc.systems[0].containsTheFractionalMolecule);
    EXPECT_FALSE(tmmc.systems[1].containsTheFractionalMolecule);
    for (std::size_t systemId : {0uz, 1uz})
    {
      System& system = tmmc.systems[systemId];
      const std::size_t fractional = system.indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(0);
      const std::uint8_t expectedGroupId = systemId == 0 ? system.components[0].lambdaGC.dUdlambdaGroupId : 0;
      for (const Atom& atom : system.spanOfMolecule(0, fractional))
      {
        EXPECT_EQ(atom.groupId, expectedGroupId) << "system " << systemId;
      }
      EXPECT_TRUE(system.interpolationGrids[methaneType].has_value()) << "system " << systemId;
      EXPECT_TRUE(system.interpolationGrids.back().has_value()) << "system " << systemId;  // Ewald-real grid
    }
    EXPECT_EQ(tmmc.streams.size(), 2uz);
  }
  std::filesystem::current_path(previous);
  std::filesystem::remove_all(scratch);
}
