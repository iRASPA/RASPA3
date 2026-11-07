#include <gtest/gtest.h>

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import mc_moves_move_types;
import mc_moves_probabilities;

// Two CFCMC groups sharing a satellite component (two salts with the same counter-ion): the shared component
// holds one group fractional slot per group, laid out driver by driver. Each group has to address its own
// slot, and a pinned lambda of one group has to reach only that group's fractional molecules.

namespace
{

ForceField makeForceField()
{
  return ForceField({{"Na", false, 22.99, 0.0, 0.0, 11, false},
                     {"K", false, 39.10, 0.0, 0.0, 19, false},
                     {"Cl", false, 35.45, 0.0, 0.0, 17, false}},
                    {{50.0, 2.5}, {50.0, 3.0}, {100.0, 4.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                    12.0, true, false, false);
}

Component makeIon(const ForceField& forceField, std::string name, std::uint16_t type, bool drivesGroup)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  if (drivesGroup) probabilities.setProbability(Move::Types::GroupSwapCFCMC, 1.0);
  Component ion(forceField, std::move(name), 0.0, 0.0, 0.0, {Atom({0, 0, 0}, 0.0, 1.0, 0, type, 0, false, false)}, {},
                {}, 5, 21, probabilities, std::nullopt, false);
  if (drivesGroup)
  {
    ion.groupComponentIds = {2};  // the chloride
    ion.maximumGroupDistance = 5.0;
  }
  return ion;
}

System makeTwoSalts()
{
  const ForceField forceField = makeForceField();
  return System(
      forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
      {makeIon(forceField, "Na", 0, true), makeIon(forceField, "K", 1, true), makeIon(forceField, "Cl", 2, false)}, {},
      {3, 3, 6}, 5);
}

}  // namespace

TEST(GROUP_FRACTIONAL_SLOTS, a_shared_satellite_holds_one_sub_slot_per_group_in_driver_order)
{
  const System system = makeTwoSalts();
  constexpr std::size_t sodium = 0, potassium = 1, chloride = 2;
  const Move::Types move = Move::Types::GroupSwapCFCMC;

  ASSERT_EQ(system.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[chloride], 2uz);
  ASSERT_EQ(system.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[sodium], 1uz);
  ASSERT_EQ(system.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[potassium], 1uz);

  EXPECT_EQ(system.groupSwapFractionalSubSlotRange(chloride, sodium, move), (std::pair{0uz, 1uz}));
  EXPECT_EQ(system.groupSwapFractionalSubSlotRange(chloride, potassium, move), (std::pair{1uz, 1uz}));
  EXPECT_EQ(system.groupSwapFractionalSubSlotRange(sodium, sodium, move), (std::pair{0uz, 1uz}));
  EXPECT_EQ(system.groupSwapFractionalSubSlotRange(sodium, potassium, move).second, 0uz);

  EXPECT_EQ(system.groupSwapLambdaDriverOfSubSlot(chloride, 0, move), std::optional{sodium});
  EXPECT_EQ(system.groupSwapLambdaDriverOfSubSlot(chloride, 1, move), std::optional{potassium});
  EXPECT_FALSE(system.groupSwapLambdaDriverOfSubSlot(chloride, 2, move).has_value());
  EXPECT_EQ(system.groupSwapLambdaDriver(chloride, move), std::optional{sodium});
}

TEST(GROUP_FRACTIONAL_SLOTS, a_pinned_lambda_reaches_only_its_own_groups_fractional_molecules)
{
  const ForceField forceField = makeForceField();
  Component sodium = makeIon(forceField, "Na", 0, true);
  Component potassium = makeIon(forceField, "K", 1, true);
  sodium.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::GroupSwap;
  sodium.fixedLambdaBin = 4;
  potassium.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::GroupSwap;
  potassium.fixedLambdaBin = 16;
  const System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
                      {sodium, potassium, makeIon(forceField, "Cl", 2, false)}, {}, {3, 3, 6}, 5);
  constexpr std::size_t chloride = 2;

  const std::size_t chlorideBegin = system.indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(chloride);
  const double sodiumScaling =
      system.spanOfMolecule(0, system.indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(0))[0].scalingVDW;
  const double potassiumScaling =
      system.spanOfMolecule(1, system.indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(1))[0].scalingVDW;
  ASSERT_NE(sodiumScaling, potassiumScaling);

  // chloride sub-slot 0 follows sodium's pinned lambda, sub-slot 1 potassium's
  EXPECT_DOUBLE_EQ(system.spanOfMolecule(chloride, chlorideBegin)[0].scalingVDW, sodiumScaling);
  EXPECT_DOUBLE_EQ(system.spanOfMolecule(chloride, chlorideBegin + 1)[0].scalingVDW, potassiumScaling);
}
