#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import molecule;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_factor;
import gradient_factor;
import running_energy;
import randomnumbers;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import mc_moves_translation;
import mc_moves_rotation;

TEST(electrostatic_polarization, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 1.0;

  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 1.0;

  RunningEnergy energy = system.computeTotalEnergies();

  // Total polarization energy for the framework + reciprocal + molecule-molecule (real-space) field model.
  // The intramolecular reciprocal-exclusion field is intentionally not part of the field (the reciprocal field
  // is built from the fixed framework structure factor only), and the molecule-molecule contribution is included
  // because omitInterPolarization == false.
  EXPECT_NEAR(energy.polarization * Units::EnergyToKelvin, -1.3034969457316241, 1e-6);
}

static double maxFieldDifference(std::span<const double3> a, std::span<const double3> b)
{
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    m = std::max(m, (a[i] - b[i]).length());
  }
  return m;
}

TEST(electrostatic_polarization, stored_field_consistency_translation)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {10}, 5);

  RandomNumber random(42);

  RunningEnergy running = system.computeTotalEnergies();

  for (std::size_t step = 0; step < 400; ++step)
  {
    std::size_t selectedMolecule = system.randomMoleculeOfComponent(random, 0);

    std::optional<RunningEnergy> e =
        MC_Moves::translationMove(random, system, 0, selectedMolecule);
    if (e.has_value()) running += e.value();
  }

  // snapshot the incrementally-maintained field, then rebuild it from scratch and compare
  std::span<double3> maintained = system.spanOfMoleculeElectricField();
  std::vector<double3> snapshot(maintained.begin(), maintained.end());

  system.computeTotalElectricField();
  std::span<double3> rebuilt = system.spanOfMoleculeElectricField();

  EXPECT_LT(maxFieldDifference(snapshot, rebuilt), 1e-8);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(running.polarization - recomputed.polarization, 0.0, 1e-6);
  EXPECT_NEAR(running.potentialEnergy() - recomputed.potentialEnergy(), 0.0, 1e-6);
}
