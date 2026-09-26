#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import atom;
import molecule;
import units;
import forcefield;
import component;
import framework;
import simulationbox;
import interpolation_energy_grid;
import randomnumbers;
import mc_moves_probabilities;
import cbmc;

// The contract of the CBMC entry points: which first-bead schemes each accepts, and that the same
// molecule regrown under its own id is excluded from its own background without any caller action.

namespace
{

ForceField makeForceField()
{
  return ForceField({{"CH3", false, 15.0, 0.0, 0.0, 6, false}, {"CH2", false, 14.0, 0.0, 0.0, 6, false}},
                    {{108.0, 3.76}, {56.0, 3.96}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

struct Fixture
{
  ForceField forceField{makeForceField()};
  TemporaryFile file{"cbmc-api-butane.json", molecule_fixtures::kButaneJson};
  Component butane{Component::Type::Adsorbate, 0, forceField, "butane", file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false};
  SimulationBox box{30.0, 30.0, 30.0};
  double beta{1.0 / (Units::KB * 300.0)};
  std::optional<Framework> noFramework{};
  std::vector<std::optional<InterpolationEnergyGrid>> noGrids{forceField.pseudoAtoms.size() + 1};
  std::optional<InterpolationEnergyGrid> noExternalFieldGrid{};

  Fixture() { butane.prepareGrowthPlans(beta); }

  CBMC::GrowContext context(std::span<const Atom> background) const
  {
    return CBMC::GrowContext(false, forceField, box, noGrids, noExternalFieldGrid, noFramework,
                             std::span<const Atom>{}, background, beta, CBMC::CutOffMode::Full);
  }
};

}  // namespace

TEST(CBMC_API, grow_new_molecule_rejects_regrow_only_schemes)
{
  Fixture f;
  RandomNumber random(1);
  const CBMC::GrowContext empty = f.context({});
  const CBMC::NewMoleculeIdentity identity{.componentId = 0, .moleculeId = 0};

  EXPECT_THROW((void)CBMC::growNewMolecule(random, empty, f.butane, identity,
                                            {.firstBead = CBMC::FirstBeadScheme::Reinsertion}),
               std::invalid_argument);
  const std::vector<std::size_t> placed{0};
  EXPECT_THROW((void)CBMC::growNewMolecule(random, empty, f.butane, identity,
                                            {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced,
                                             .beadsAlreadyPlaced = placed}),
               std::invalid_argument);
  // Pinned and Fixed need a position.
  EXPECT_THROW((void)CBMC::growNewMolecule(random, empty, f.butane, identity, {.firstBead = CBMC::FirstBeadScheme::Pinned}),
               std::invalid_argument);
}

TEST(CBMC_API, regrow_molecule_rejects_new_molecule_schemes)
{
  Fixture f;
  RandomNumber random(2);
  const CBMC::GrowContext empty = f.context({});

  std::optional<CBMC::GrowResult> grown;
  while (!grown) grown = CBMC::growNewMolecule(random, empty, f.butane, {.componentId = 0, .moleculeId = 0});

  for (CBMC::FirstBeadScheme scheme :
       {CBMC::FirstBeadScheme::MultipleFirstBead, CBMC::FirstBeadScheme::Pinned, CBMC::FirstBeadScheme::Fixed})
  {
    EXPECT_THROW((void)CBMC::regrowMolecule(random, empty, f.butane, grown->molecule, grown->atoms,
                                            {.firstBead = scheme, .firstBeadPosition = grown->atoms[0].position}),
                 std::invalid_argument);
  }
  // AlreadyPlaced needs a placed set.
  EXPECT_THROW((void)CBMC::regrowMolecule(random, empty, f.butane, grown->molecule, grown->atoms,
                                          {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced}),
               std::invalid_argument);
}

// A molecule regrown against a background that still holds its old copy must not feel that copy:
// the exclusion is by molecule id, so the retrace of the old configuration in a background that
// consists only of itself sees no inter-molecular energy at all.
TEST(CBMC_API, molecule_is_excluded_from_its_own_background_by_id)
{
  Fixture f;
  RandomNumber random(3);
  const CBMC::GrowContext empty = f.context({});

  std::optional<CBMC::GrowResult> grown;
  while (!grown) grown = CBMC::growNewMolecule(random, empty, f.butane, {.componentId = 0, .moleculeId = 7});
  const std::vector<Atom> molecule = grown->atoms;

  // The background is the molecule itself: identical to no background at all.
  RandomNumber randomA(11), randomB(11);
  const CBMC::RetraceResult alone = CBMC::retraceMolecule(randomA, empty, f.butane, molecule);
  const CBMC::RetraceResult selfBackground = CBMC::retraceMolecule(randomB, f.context(molecule), f.butane, molecule);
  EXPECT_DOUBLE_EQ(alone.logRosenbluthWeight, selfBackground.logRosenbluthWeight);
  EXPECT_DOUBLE_EQ(alone.energies.moleculeMoleculeVDW, 0.0);
  EXPECT_DOUBLE_EQ(selfBackground.energies.moleculeMoleculeVDW, 0.0);

  // A copy under another id at the same place is a real neighbour (a hard overlap here), unless the
  // context skips it.
  std::vector<Atom> otherId = molecule;
  for (Atom &atom : otherId) atom.moleculeId = 8;
  EXPECT_THROW((void)CBMC::retraceMolecule(randomB, f.context(otherId), f.butane, molecule), std::runtime_error);

  RandomNumber randomC(11);
  const CBMC::RetraceResult skipped =
      CBMC::retraceMolecule(randomC, f.context(otherId).withSkippedMolecule(8), f.butane, molecule);
  EXPECT_DOUBLE_EQ(alone.logRosenbluthWeight, skipped.logRosenbluthWeight);
}
