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
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

TEST(MC, translation)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {5, 3}, 5);

  [[maybe_unused]] std::span<Atom> atomData = system.spanOfMoleculeAtoms();
}

// 'rebuildConfigurationDerivedState' (replica exchange, solute tempering) recounts the pseudo-atoms after the
// framework atoms have been placed in front of the molecule atoms: the per-type totals must stay the adsorbate
// totals and must not absorb the framework atoms.
TEST(MC, rebuild_configuration_derived_state_keeps_molecule_pseudo_atom_counts)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(1, 1, 1));
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {5, 3}, 5);
  ASSERT_GT(system.numberOfFrameworkAtoms, 0uz);

  std::vector<std::size_t> expected(forceField.pseudoAtoms.size(), 0uz);
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    expected[static_cast<std::size_t>(atom.type)] += 1;
  }
  EXPECT_EQ(system.totalNumberOfPseudoAtoms, expected);

  system.rebuildConfigurationDerivedState();

  EXPECT_EQ(system.totalNumberOfPseudoAtoms, expected);
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const std::size_t count = std::accumulate(system.numberOfPseudoAtoms[componentId].begin(),
                                              system.numberOfPseudoAtoms[componentId].end(), 0uz);
    EXPECT_EQ(count, system.numberOfMoleculesPerComponent[componentId] * system.components[componentId].atoms.size());
  }
}
