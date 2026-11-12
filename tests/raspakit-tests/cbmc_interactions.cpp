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
import simulationbox;
import running_energy;
import cbmc_external_energy;
import cbmc_torsion_selection;
import cbmc_grow_step;
import chiral_center;
import randomnumbers;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(cbmc_interactions, framework_molecule_1)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeMethane(forceField, 0);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);

  std::optional<RunningEnergy> frameworkMoleculeEnergy =
      CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                           system.framework, frameworkAtoms, 12.0, 12.0, atomData);

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);
}

TEST(cbmc_interactions, framework_molecule_2)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeMethane(forceField, 0);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);
  atomData[1].position = double3(0.6 * 11.8671, 0.7 * 11.8671, 0.65 * 11.8671);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);

  std::optional<RunningEnergy> frameworkMoleculeEnergy =
      CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                           system.framework, frameworkAtoms, 12.0, 12.0, atomData);

  std::optional<RunningEnergy> interMoleculeEnergy1 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomData, 12.0, 12.0, {atomData.begin(), 1});

  std::optional<RunningEnergy> interMoleculeEnergy2 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomData, 12.0, 12.0, {atomData.begin() + 1, 1});

  std::optional<RunningEnergy> interMoleculeEnergy3 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomData, 12.0, 12.0, {atomData.begin(), 2});

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy1->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy2->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy3->moleculeMoleculeVDW * Units::EnergyToKelvin / 2.0, 2352.42793591, 1e-6);
}

// A declared chiral centre whose parity the spin decides, with its placed neighbour on the spin axis: the parity is
// then the same for every spin, and declared the other way round every trial is behind the wall. That spin is dead
// (weight zero), not a selection among trials of zero weight, whose shifted sum would be NaN. The guard cannot rest
// on a log factor of minus infinity: the project's fast-math flags fold the finiteness test away.
TEST(cbmc_interactions, torsion_selection_with_every_spin_behind_the_chiral_wall_is_a_dead_spin)
{
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  const SimulationBox box(30.0, 30.0, 30.0);
  auto atom = [](double3 position)
  { return Atom(position, 0.0, 1.0, 0, std::uint16_t{0}, std::uint8_t{0}, false, false); };
  // previous (0) to current (1) along z, the placed neighbour (4) further along the axis, the grown beads (2, 3)
  // off it
  std::vector<Atom> chain{atom({0.0, 0.0, 0.0}), atom({0.0, 0.0, 1.5}), atom({1.0, 0.0, 2.0}), atom({0.0, 1.0, 2.0}),
                          atom({0.0, 0.0, 3.0})};
  const std::vector<Atom> base{chain[2], chain[3]};

  CBMC::GrowStep step{};
  step.previousBead = 0;
  step.currentBead = 1;
  step.nextBeads = {2, 3};
  ChiralCenter centre{};
  centre.ids = {1, 2, 4, 3};
  // the signed volume of this geometry is negative for every spin: the declaration of the other hand
  centre.type = ChiralCenter::Chirality::R_Chiral;
  step.spin.determinedChiralCenters = {centre};
  const double3 axis(0.0, 0.0, 1.0);

  RandomNumber random(7);
  const CBMC::TorsionOrientation dead =
      CBMC::selectTorsionOrientation(random, 4, 1.0, forceField, box, chain, base, step, axis, false);
  // a NaN weight compares false against anything, which the fast-math flags let the compiler fold; its text does not
  EXPECT_EQ(std::to_string(dead.rosenbluthWeight), std::to_string(0.0));
  EXPECT_EQ(dead.positions.size(), 2uz);
  EXPECT_EQ(chain[2].position.x, 1.0);
  EXPECT_EQ(chain[3].position.y, 1.0);

  // the declaration of this hand passes every spin: a live selection among four equal trials, weight one
  centre.type = ChiralCenter::Chirality::S_Chiral;
  step.spin.determinedChiralCenters = {centre};
  const CBMC::TorsionOrientation live =
      CBMC::selectTorsionOrientation(random, 4, 1.0, forceField, box, chain, base, step, axis, false);
  EXPECT_NEAR(live.rosenbluthWeight, 1.0, 1e-12);
}
