#include <gtest/gtest.h>

import std;

import int3;
import double3;
import units;
import atom;
import atom_dynamics;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import interactions_ewald;
import potential_coulomb_real_space;

// These tests verify that the finite-cutoff charge methods (Wolf and its relatives) include the
// intramolecular exclusion / completion of the shifted pair sum, i.e. terms S62-S63 of the
// Brick-CFCMC formulation (Dubbeldam et al.). Per intramolecular pair inside the Coulomb cutoff the
// correction is q_i q_j (V(r) - 1/r), where V(r) is the method's shifted real-space potential.

namespace
{
// Analytic per-molecule intramolecular exclusion energy for a finite-cutoff charge method.
double referenceExclusionEnergy(const ForceField& forceField, const SimulationBox& box, std::span<const Atom> atoms)
{
  const double cutOffSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  double energy = 0.0;
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j != atoms.size(); ++j)
    {
      double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
      double rr = double3::dot(dr, dr);
      if (rr >= cutOffSquared) continue;
      double r = std::sqrt(rr);
      const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
      energy += atoms[i].scalingCoulomb * atoms[j].scalingCoulomb * Units::CoulombicConversionFactor *
                atoms[i].charge * atoms[j].charge * (factors.potential - 1.0 / r);
    }
  }
  return energy;
}
}  // namespace

TEST(wolf_exclusion, intramolecular_exclusion_matches_analytic_sum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;

  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {c}, {}, {2}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  atomData[0].position = double3(-5.0, 0.0, 1.149);
  atomData[1].position = double3(-5.0, 0.0, 0.0);
  atomData[2].position = double3(-5.0, 0.0, -1.149);
  atomData[3].position = double3(5.0, 0.0, 1.149);
  atomData[4].position = double3(5.0, 0.0, 0.0);
  atomData[5].position = double3(5.0, 0.0, -1.149);

  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  double expectedExclusion = referenceExclusionEnergy(system.forceField, system.simulationBox, atomData.subspan(0, 3)) +
                             referenceExclusionEnergy(system.forceField, system.simulationBox, atomData.subspan(3, 3));

  // The Wolf method has no reciprocal part.
  EXPECT_EQ(energy.ewald_fourier, 0.0);

  // The exclusion must be present and equal to the analytic per-molecule sum (the previous
  // implementation left this term at zero).
  EXPECT_NE(expectedExclusion, 0.0);
  EXPECT_NEAR(energy.ewald_exclusion, expectedExclusion, 1e-8);

  // The self term is the standard Wolf self energy over all atoms.
  double sumChargeSquared = 0.0;
  for (const Atom& atom : atomData)
  {
    sumChargeSquared += (atom.scalingCoulomb * atom.charge) * (atom.scalingCoulomb * atom.charge);
  }
  double expectedSelf =
      Units::CoulombicConversionFactor * Potentials::coulombSelfEnergyPrefactor(system.forceField) * sumChargeSquared;
  EXPECT_NEAR(energy.ewald_self, expectedSelf, 1e-8);
}

TEST(wolf_exclusion, rigid_translation_leaves_self_and_exclusion_unchanged)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;

  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  atomData[0].position = double3(0.0, 0.0, 1.149);
  atomData[1].position = double3(0.0, 0.0, 0.0);
  atomData[2].position = double3(0.0, 0.0, -1.149);

  std::vector<Atom> oldatoms(atomData.begin(), atomData.end());

  // Rigidly translate the molecule; internal distances are unchanged.
  std::vector<Atom> newatoms = oldatoms;
  for (Atom& atom : newatoms) atom.position += double3(3.7, -2.1, 0.9);

  RunningEnergy difference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, std::span<const Atom>(newatoms), std::span<const Atom>(oldatoms));

  // A rigid translation changes neither the self energy nor the intramolecular exclusion.
  EXPECT_NEAR(difference.ewald_self, 0.0, 1e-10);
  EXPECT_NEAR(difference.ewald_exclusion, 0.0, 1e-10);
}
