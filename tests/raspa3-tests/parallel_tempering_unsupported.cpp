#include <gtest/gtest.h>

import std;

import int3;
import forcefield;
import vdwparameters;
import framework;
import component;
import system;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;

// The replica-exchange swap can not handle a framework with mobile atoms (it exchanges the adsorbates only),
// temperature-dependent pair potentials (the replicas then have different Hamiltonians), or solute tempering
// with polarization (no polarization term in the solute Hamiltonian change). Such systems are refused with a
// reason, by the drivers up front and by the swap itself, instead of silently refusing every swap or swapping
// a state the energies do not describe.

namespace
{

System makeCO2InFAU(double temperature)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  return System(forceField, std::nullopt, false, temperature, 1e5, 1.0, {framework}, {co2}, {}, {4}, 5);
}

}  // namespace

TEST(PARALLEL_TEMPERING_UNSUPPORTED, a_plain_system_is_supported)
{
  const System system = makeCO2InFAU(300.0);
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(system, false).has_value());
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(system, true).has_value());
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(system, makeCO2InFAU(350.0)).has_value());
}

TEST(PARALLEL_TEMPERING_UNSUPPORTED, a_framework_with_mobile_atoms_is_refused)
{
  System systemA = makeCO2InFAU(300.0);
  System systemB = makeCO2InFAU(350.0);
  systemA.framework->rigid = false;
  systemB.framework->rigid = false;
  ASSERT_TRUE(systemA.framework->hasMobileAtoms());
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(systemA, false).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(systemA, systemB).has_value());
}

TEST(PARALLEL_TEMPERING_UNSUPPORTED, a_temperature_dependent_potential_is_refused)
{
  System systemA = makeCO2InFAU(300.0);
  System systemB = makeCO2InFAU(350.0);
  systemA.forceField.data.front().type = VDWParameters::Type::FeynmannHibbs;
  systemB.forceField.data.front().type = VDWParameters::Type::FeynmannHibbs;
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(systemA, false).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(systemA, systemB).has_value());
}

TEST(PARALLEL_TEMPERING_UNSUPPORTED, solute_tempering_with_polarization_is_refused)
{
  System system = makeCO2InFAU(300.0);
  system.forceField.computePolarization = true;
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(system, true).has_value());
  // plain temperature exchange with polarization stays supported (the energies are those of one Hamiltonian)
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(system, false).has_value());
}
