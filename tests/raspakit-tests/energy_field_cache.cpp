#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_backend;

// The molecular field held between routes. It is handed back when the same landscape is asked for again,
// and only then: a changed cutoff, a changed force field, a changed probe or a changed framework under the
// same names is a different landscape and must be built afresh.

namespace
{

constexpr double edge = 12.0;
const uint3 gridSize{6, 6, 6};

PairInteractions probeAndAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "O"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0}};
  interactions.cutOffVDW = 5.0;
  interactions.cutOffCoulomb = 5.0;
  return interactions;
}

Crystal oneAtom()
{
  Crystal framework;
  framework.name = "held";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.5 * edge, 0.5 * edge), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

MolecularField build(const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe)
{
  return buildMolecularField(cpuEnergyBackend(), interactions, framework, probe, gridSize, 1, 300.0, false, 1.0e-6);
}

}  // namespace


TEST(energy_field_cache, the_same_landscape_is_handed_back_and_a_changed_one_is_not)
{
  forgetMolecularField();
  PairInteractions interactions = probeAndAtom();
  Crystal framework = oneAtom();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();

  MolecularField first = build(interactions, framework, probe);
  EXPECT_FALSE(first.reused);
  MolecularField again = build(interactions, framework, probe);
  EXPECT_TRUE(again.reused);

  // The same names, a shorter dispersion cutoff.
  PairInteractions shorter = interactions;
  shorter.cutOffVDW = 3.5;
  EXPECT_FALSE(build(shorter, framework, probe).reused);

  // The same names, a different strength for the pair.
  PairInteractions stronger = interactions;
  for (PairParameters &pair : stronger.parameters) pair.strengthParameter = 200.0;
  EXPECT_FALSE(build(stronger, framework, probe).reused);

  // The same name, count and volume, the atom moved.
  Crystal moved = framework;
  moved.atoms[0].position = double3(0.25 * edge, 0.5 * edge, 0.5 * edge);
  moved.fractionalPositions[0] = double3(0.25, 0.5, 0.5);
  EXPECT_FALSE(build(interactions, moved, probe).reused);

  // The same probe name, a site hung off the centre.
  LinearProbe longer = probe;
  longer.sites[0].offset = 0.5;
  EXPECT_FALSE(build(interactions, framework, longer).reused);

  forgetMolecularField();
}
