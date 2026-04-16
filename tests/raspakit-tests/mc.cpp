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
