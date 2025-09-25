#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <ranges>
#include <span>
#include <tuple>
#include <vector>

import int3;
import double3;
import double3x3;
import factory;
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
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {5, 3}, 5);

  [[maybe_unused]] std::span<Atom> atomData = system.spanOfMoleculeAtoms();
}
