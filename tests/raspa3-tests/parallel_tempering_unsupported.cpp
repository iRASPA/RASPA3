#include <gtest/gtest.h>

import std;

import int3;
import forcefield;
import simulationbox;
import atom;
import double3;
import vdwparameters;
import framework;
import component;
import system;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;
import property_msd;
import property_vacf;

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

namespace
{
// a box of single-site Lennard-Jones particles, the cut-offs and Ewald parameters derived from its box as the
// drivers do
System makeParticleBox(double boxLength, double temperature)
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{120.0, 3.5}},
                              ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false, false);
  const Component particle(forceField, "particle", 190.0, 4.6e6, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, SimulationBox(boxLength, boxLength, boxLength), false, temperature, 1.0e5, 1.0, {},
                {particle}, {}, {20}, 5);
  // the Coulomb cut-off follows the box, as with the default input ('CutOffCoulomb' auto)
  system.forceField.cutOffCoulombAutomatic = true;
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}
}  // namespace

// Replicas at different volumes (NPT replica exchange) share one Hamiltonian: the cut-offs the force field derives
// from each box automatically are caches of that box, like the automatic Ewald parameters, and must not make the
// swap refuse every attempt.
TEST(PARALLEL_TEMPERING_UNSUPPORTED, replicas_at_different_volumes_with_automatic_cutoffs_share_the_hamiltonian)
{
  const System a = makeParticleBox(22.0, 300.0);
  const System b = makeParticleBox(24.0, 330.0);
  ASSERT_TRUE(a.forceField.cutOffCoulombAutomatic);
  ASSERT_NE(a.forceField.cutOffCoulomb, b.forceField.cutOffCoulomb);
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(a, b).has_value());
}

// After the configurations (and boxes) of two replicas are exchanged, the automatic cut-offs of each force field
// follow the box that arrived, so the minimum image holds in the smaller box.
TEST(PARALLEL_TEMPERING_UNSUPPORTED, an_exchanged_box_reinitializes_the_automatic_cutoffs)
{
  System a = makeParticleBox(22.0, 300.0);
  System b = makeParticleBox(24.0, 330.0);
  a.exchangeConfiguration(b);
  a.rebuildConfigurationDerivedState();
  b.rebuildConfigurationDerivedState();
  EXPECT_DOUBLE_EQ(a.simulationBox.lengthA, 24.0);
  EXPECT_NEAR(a.forceField.cutOffCoulomb, 12.0, 1e-9);
  EXPECT_NEAR(b.forceField.cutOffCoulomb, 11.0, 1e-9);
}

// An ideal-gas Rosenbluth weight belongs to one temperature: copied to every replica of a ladder it normalises the
// insertions of the other temperatures with the wrong intramolecular partition function.
TEST(PARALLEL_TEMPERING_UNSUPPORTED, an_ideal_gas_rosenbluth_weight_is_refused)
{
  System system = makeCO2InFAU(300.0);
  system.components[0].idealGasRosenbluthWeight = 0.37;
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(system, false).has_value());
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(system, true).has_value());
  system.components[0].idealGasRosenbluthWeight = 1.0;
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(system, false).has_value());
}

// A time-correlation property belongs to the trajectory it samples: a replica-exchange swap carries its time
// origins and sums along with the configuration, so each correlation is taken over one continuous walker.
TEST(PARALLEL_TEMPERING_UNSUPPORTED, time_correlation_bookkeeping_travels_with_the_configuration)
{
  System systemA = makeCO2InFAU(300.0);
  System systemB = makeCO2InFAU(350.0);
  systemA.setPropertyMSD(PropertyMeanSquaredDisplacement(systemA.numberOfMoleculesPerComponent,
                                                         systemA.moleculeData.size(), 0.0005, 25, 1, 100));
  systemB.setPropertyMSD(PropertyMeanSquaredDisplacement(systemB.numberOfMoleculesPerComponent,
                                                         systemB.moleculeData.size(), 0.0005, 25, 1, 100));
  systemA.propertyMSD->countMSD = 3;
  systemB.propertyMSD->countMSD = 7;
  systemA.setPropertyVACF(PropertyVelocityAutoCorrelationFunction(systemA.numberOfMoleculesPerComponent,
                                                                  systemA.moleculeData.size(), 0.0005, 4, 10, 1, 100));
  systemA.propertyVACF->countAccumulatedVACF = 11;

  systemA.exchangeConfiguration(systemB);

  ASSERT_TRUE(systemA.propertyMSD.has_value());
  ASSERT_TRUE(systemB.propertyMSD.has_value());
  EXPECT_EQ(systemA.propertyMSD->countMSD, 7uz);
  EXPECT_EQ(systemB.propertyMSD->countMSD, 3uz);
  EXPECT_FALSE(systemA.propertyVACF.has_value());
  ASSERT_TRUE(systemB.propertyVACF.has_value());
  EXPECT_EQ(systemB.propertyVACF->countAccumulatedVACF, 11uz);
}
