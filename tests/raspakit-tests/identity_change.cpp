#include <gtest/gtest.h>

import std;

import int3;
import double3;
import randomnumbers;
import atom;
import forcefield;
import framework;
import component;
import molecule;
import system;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_identity_change;
import mc_moves_reinsertion;

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

Component makeAlkaneFromExample(const ForceField &forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities &probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), moleculePath.string(), 5,
                   21, probabilities, std::nullopt, false);
}

void configureSemiGrandIdentityChange(Component &component, std::span<const std::size_t> identityChanges,
                                      double molFraction, const MCMoveProbabilities &probabilities)
{
  component.mc_moves_probabilities = probabilities;
  component.identityChanges.assign(identityChanges.begin(), identityChanges.end());
  component.molFraction = molFraction;
  component.swappable = true;
  component.idealGasRosenbluthWeight = 1.0;
}

}  // namespace

TEST(identity_change, conserves_total_molecule_count_and_energy)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  Component methane = Component::makeMethane(forceField, 1);

  for (Component *component : {&co2, &methane})
  {
    component->mc_moves_probabilities = probabilities;
    component->identityChanges = {0, 1};
    component->molFraction = (component == &co2) ? 0.25 : 0.75;
    component->swappable = true;
    component->idealGasRosenbluthWeight = 1.0;
  }

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {f}, {co2, methane}, {}, {8, 8}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(42);
  const std::size_t totalMolecules =
      system.numberOfMoleculesPerComponent[0] + system.numberOfMoleculesPerComponent[1];
  std::size_t accepted = 0;

  for (std::size_t i = 0; i < 1000; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 1);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
    {
      continue;
    }

    const std::optional<RunningEnergy> energyDifference =
        MC_Moves::identityChangeMove(random, system, selectedComponent);
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
      ++accepted;
    }
  }

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0] + system.numberOfMoleculesPerComponent[1], totalMolecules);
  EXPECT_GT(accepted, 0uz);

  const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
  EXPECT_NEAR(drift.tail, 0.0, 1e-5);
}

TEST(identity_change, reinsertion_after_identity_change_conserves_energy)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  Component methane = Component::makeMethane(forceField, 1);

  for (Component *component : {&co2, &methane})
  {
    component->mc_moves_probabilities = probabilities;
    component->identityChanges = {0, 1};
    component->molFraction = (component == &co2) ? 0.25 : 0.75;
    component->swappable = true;
    component->idealGasRosenbluthWeight = 1.0;
  }

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {f}, {co2, methane}, {}, {8, 8}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(42);
  std::size_t identityAccepted = 0;

  for (std::size_t i = 0; i < 5000 && identityAccepted == 0; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 1);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
    {
      continue;
    }

    const std::optional<RunningEnergy> energyDifference =
        MC_Moves::identityChangeMove(random, system, selectedComponent);
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
      ++identityAccepted;
    }
  }

  ASSERT_GT(identityAccepted, 0uz);

  for (std::size_t i = 0; i < 1000; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 1);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
    {
      continue;
    }

    const std::size_t selectedMolecule = system.randomMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
    const std::size_t molecule_index = system.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
    Molecule &molecule = system.moleculeData[molecule_index];

    const std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
        random, system, selectedComponent, selectedMolecule, molecule, molecule_atoms);
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
    }
  }

  const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
}

TEST(identity_change, incremental_delta_matches_recompute_per_accept)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  Component methane = Component::makeMethane(forceField, 1);

  for (Component *component : {&co2, &methane})
  {
    component->mc_moves_probabilities = probabilities;
    component->identityChanges = {0, 1};
    component->molFraction = (component == &co2) ? 0.25 : 0.75;
    component->swappable = true;
    component->idealGasRosenbluthWeight = 1.0;
  }

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {f}, {co2, methane}, {}, {8, 8}, 5);

  RandomNumber random(42);
  std::size_t accepted = 0;

  for (std::size_t i = 0; i < 5000 && accepted < 50; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 1);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
    {
      continue;
    }

    const RunningEnergy energyBefore = system.computeTotalEnergies();
    const std::optional<RunningEnergy> energyDifference =
        MC_Moves::identityChangeMove(random, system, selectedComponent);
    if (!energyDifference)
    {
      continue;
    }

    const RunningEnergy energyAfter = system.computeTotalEnergies();
    const RunningEnergy recomputeDelta = energyAfter - energyBefore;
    const RunningEnergy deltaError = energyDifference.value() - recomputeDelta;

    EXPECT_NEAR(deltaError.potentialEnergy(), 0.0, 1e-5) << "accept " << accepted;
    EXPECT_NEAR(deltaError.moleculeMoleculeVDW, 0.0, 1e-5);
    EXPECT_NEAR(deltaError.moleculeMoleculeCharge, 0.0, 1e-5);

    system.runningEnergies += energyDifference.value();
    ++accepted;
  }

  EXPECT_GT(accepted, 0uz);
}

TEST(identity_change, co2_propane_butane_conserves_total_molecule_count_and_energy)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);

  constexpr std::array<std::size_t, 3> identityChanges = {0, 1, 2};

  Component co2 = Component::makeCO2(forceField, 0, true);
  configureSemiGrandIdentityChange(co2, identityChanges, 0.25, probabilities);

  Component propane = makeAlkaneFromExample(forceField, 1, "propane", probabilities);
  configureSemiGrandIdentityChange(propane, identityChanges, 0.375, probabilities);

  Component butane = makeAlkaneFromExample(forceField, 2, "butane", probabilities);
  configureSemiGrandIdentityChange(butane, identityChanges, 0.375, probabilities);

  ASSERT_EQ(propane.growType, Component::GrowType::Flexible);
  ASSERT_EQ(butane.growType, Component::GrowType::Flexible);

  System system = System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {f}, {co2, propane, butane}, {},
                         {4, 4, 4}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(42);
  const std::size_t totalMolecules = system.numberOfMoleculesPerComponent[0] +
                                     system.numberOfMoleculesPerComponent[1] +
                                     system.numberOfMoleculesPerComponent[2];
  std::size_t accepted = 0;

  for (std::size_t i = 0; i < 5000; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 2);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
    {
      continue;
    }

    const std::optional<RunningEnergy> energyDifference =
        MC_Moves::identityChangeMove(random, system, selectedComponent);
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
      ++accepted;
    }
  }

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0] + system.numberOfMoleculesPerComponent[1] +
                system.numberOfMoleculesPerComponent[2],
            totalMolecules);
  EXPECT_GT(accepted, 0uz);

  const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
  EXPECT_NEAR(drift.intraVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.intraCoul, 0.0, 1e-5);
  EXPECT_NEAR(drift.tail, 0.0, 1e-5);
}

TEST(identity_change, cfcmc_preserves_fractional_molecules_and_integer_energy)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component co2 = Component::makeCO2(forceField, 0, true);
  Component methane = Component::makeMethane(forceField, 1);

  for (Component *component : {&co2, &methane})
  {
    component->mc_moves_probabilities = probabilities;
    component->identityChanges = {0, 1};
    component->molFraction = (component == &co2) ? 0.25 : 0.75;
    component->swappable = true;
    component->idealGasRosenbluthWeight = 1.0;
  }

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e6, 1.0, {f}, {co2, methane}, {}, {8, 8}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  ASSERT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 1uz);
  ASSERT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 1uz);

  RandomNumber random(42);
  const std::size_t totalIntegerMolecules = system.numberOfIntegerMoleculesPerComponent[0] +
                                            system.numberOfIntegerMoleculesPerComponent[1];
  std::size_t accepted = 0;

  for (std::size_t i = 0; i < 5000; ++i)
  {
    const std::size_t selectedComponent = random.uniform_integer(0, 1);
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0 &&
        system.numberOfIntegerMoleculesPerComponent[1 - selectedComponent] == 0)
    {
      continue;
    }

    const std::optional<RunningEnergy> energyDifference =
        MC_Moves::identityChangeMove(random, system, selectedComponent);
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
      ++accepted;
    }
  }

  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0] + system.numberOfIntegerMoleculesPerComponent[1],
            totalIntegerMolecules);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 1uz);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 1uz);
  EXPECT_GT(accepted, 0uz);

  const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-5);
  EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-5);
  EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-5);
  EXPECT_NEAR(drift.tail, 0.0, 1e-5);
}
