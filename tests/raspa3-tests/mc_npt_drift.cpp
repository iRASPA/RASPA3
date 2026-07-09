#include <gtest/gtest.h>

import std;

import double3;
import units;
import atom;
import forcefield;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;

TEST(MC_NPT_DRIFT, anisotropic_volume_change)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::AnisotropicVolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2}, {}, {10}, 5, systemProbabilities);
  system.pressureTensorDiagonal = double3(1e4, 1.2e4, 0.8e4) / Units::PressureConversionFactor;

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, isotropic_volume_change)
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
                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::VolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2, methane}, {}, {10, 15}, 5,
                         systemProbabilities);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, isotropic_volume_change_flexible)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  // Flexible (bonded) two-site "ethane"-like molecule: a harmonic bond makes the component flexible, so it is
  // grown/moved with intramolecular potentials. Volume moves must scale it by its center of mass, leaving the
  // internal geometry (and thus the bond energy) unchanged.
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component ethane = Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                              {Atom({0.0, 0.0, 0.77}, 0.0, 1.0, 0, 6, 0, false, false),
                               Atom({0.0, 0.0, -0.77}, 0.0, 1.0, 0, 6, 0, false, false)},
                              connectivityTable, intraMolecularPotentials, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::VolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {ethane}, {}, {20}, 5, systemProbabilities);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.bond, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, anisotropic_volume_change_flexible)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  // Flexible (bonded) two-site "ethane"-like molecule (see isotropic_volume_change_flexible). Anisotropic volume
  // moves scale each molecule by its center of mass, so bond lengths (and the bond energy) must be preserved even
  // when the three box directions scale by different factors.
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component ethane = Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                              {Atom({0.0, 0.0, 0.77}, 0.0, 1.0, 0, 6, 0, false, false),
                               Atom({0.0, 0.0, -0.77}, 0.0, 1.0, 0, 6, 0, false, false)},
                              connectivityTable, intraMolecularPotentials, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::AnisotropicVolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {ethane}, {}, {20}, 5, systemProbabilities);
  system.pressureTensorDiagonal = double3(1e4, 1.2e4, 0.8e4) / Units::PressureConversionFactor;

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo(20, 5, 5, 1000, 10000, 5000, 5000, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.bond, 0.0, 1e-6);
  }
}
