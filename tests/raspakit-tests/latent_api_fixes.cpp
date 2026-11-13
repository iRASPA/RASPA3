#include <gtest/gtest.h>

import std;

import int3;
import double3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import simulationbox;
import system;
import running_energy;
import interpolation_energy_grid;
import interactions_intermolecular;
import interactions_polarization;
import interactions_framework_molecule;
import integrators_cputime;
import json;

namespace
{
ForceField makeChargedForceField()
{
  return ForceField({PseudoAtom("A", false, 12.0, 0.4, 1.5, 6, false)}, {VDWParameters(120.0, 3.4)},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
}
}  // namespace

// With the intermolecular polarization switched off, the field-difference routines still return the van der Waals
// and Coulomb energy difference (they returned zero, while their contract is the energy of the move) and leave the
// field buffers alone.
TEST(latent_api_fixes, field_difference_energies_do_not_depend_on_the_polarization_switch)
{
  ForceField forceField = makeChargedForceField();
  forceField.omitInterPolarization = true;
  const SimulationBox box(30.0, 30.0, 30.0);
  const std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.4, 1.0, 0, 0, 0, 0, false),
                                Atom({9.0, 5.0, 5.0}, -0.4, 1.0, 1, 0, 0, 0, false)};
  const std::vector<Atom> newAtoms{Atom({8.4, 5.6, 5.2}, -0.4, 1.0, 1, 0, 0, 0, false)};
  const std::span<const Atom> oldAtoms(&atoms[1], 1);

  const std::optional<RunningEnergy> reference =
      Interactions::computeInterMolecularEnergyDifference(forceField, box, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(reference.has_value());
  ASSERT_NE(reference->moleculeMoleculeVDW, 0.0);
  ASSERT_NE(reference->moleculeMoleculeCharge, 0.0);

  std::vector<double3> fieldMolecules(atoms.size());
  std::vector<double3> fieldMolecule(newAtoms.size());
  const std::optional<RunningEnergy> withField = Interactions::computeInterMolecularElectricFieldDifference(
      forceField, box, fieldMolecules, fieldMolecule, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(withField.has_value());
  EXPECT_NEAR(withField->moleculeMoleculeVDW, reference->moleculeMoleculeVDW, 1e-10);
  EXPECT_NEAR(withField->moleculeMoleculeCharge, reference->moleculeMoleculeCharge, 1e-10);

  std::vector<double3> neighbourDelta(atoms.size());
  std::vector<double3> fieldNew(newAtoms.size());
  std::vector<double3> fieldOld(newAtoms.size());
  const std::optional<RunningEnergy> polarizationVariant =
      Interactions::computeInterMolecularPolarizationElectricFieldDifference(forceField, box, neighbourDelta, fieldNew,
                                                                             fieldOld, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(polarizationVariant.has_value());
  EXPECT_NEAR(polarizationVariant->moleculeMoleculeVDW, reference->moleculeMoleculeVDW, 1e-10);
  EXPECT_NEAR(polarizationVariant->moleculeMoleculeCharge, reference->moleculeMoleculeCharge, 1e-10);

  for (const std::vector<double3>* buffer : {&fieldMolecules, &fieldMolecule, &neighbourDelta, &fieldNew, &fieldOld})
  {
    for (const double3& field : *buffer) EXPECT_EQ(double3::dot(field, field), 0.0);
  }
}

// The three-argument polarization difference takes the new field first; an insertion passes no old field, which
// gives no old term (the old field was read past its end).
TEST(latent_api_fixes, polarization_difference_of_an_insertion_has_no_old_term)
{
  const ForceField forceField = makeChargedForceField();
  std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.4, 1.0, 0, 0, 0, 0, false)};
  std::vector<double3> newField{double3(0.3, -0.2, 0.1)};
  std::vector<double3> oldField{double3(0.1, 0.1, 0.1)};

  const double polarizability = forceField.pseudoAtoms[0].polarizability / Units::CoulombicConversionFactor;
  const double newTerm = -0.5 * polarizability * double3::dot(newField[0], newField[0]);
  const double oldTerm = 0.5 * polarizability * double3::dot(oldField[0], oldField[0]);

  EXPECT_NEAR(Interactions::computePolarizationEnergyDifference(forceField, newField, {}, atoms).polarization, newTerm,
              1e-14);
  EXPECT_NEAR(Interactions::computePolarizationEnergyDifference(forceField, newField, oldField, atoms).polarization,
              newTerm + oldTerm, 1e-14);
}

// A caller without grids (an empty grid vector) gets the pair sum: the per-type grid was indexed without a bounds
// check, and the Coulomb grid read without checking it exists.
TEST(latent_api_fixes, framework_molecule_gradient_without_grids_evaluates_the_pairs)
{
  ForceField forceField = makeChargedForceField();
  const SimulationBox box(24.0, 24.0, 24.0);
  const Atom frameworkAtom({0.25, 0.25, 0.25}, 0.4, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "framework", box, 1, {frameworkAtom}, {frameworkAtom}, int3(1, 1, 1));
  Component component(forceField, "site", 100.0, 1.0e6, 0.1, {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, 0, 0, false, false)},
                      ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {component}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = {7.0, 8.0, 9.0};
  system.spanOfMoleculeAtoms()[0].position = {10.6, 8.4, 8.7};

  const RunningEnergy reference = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  ASSERT_NE(reference.frameworkMoleculeVDW, 0.0);

  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids{};
  const RunningEnergy gradient = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.spanOfMoleculeDynamics(), noGrids, system.framework);
  EXPECT_NEAR(gradient.frameworkMoleculeVDW, reference.frameworkMoleculeVDW, 1e-10);
  EXPECT_NEAR(gradient.frameworkMoleculeCharge, reference.frameworkMoleculeCharge, 1e-10);
}

// The molecular-dynamics difference table is defined (it was declared only), and the MD outputs carry the
// thermobarostat and conserved energies; the timing total leaves out the rotations nested in the free rotor.
TEST(latent_api_fixes, md_energy_outputs_carry_every_term_and_the_timing_total_counts_once)
{
  RunningEnergy energy{};
  energy.moleculeMoleculeVDW = 1.0;
  energy.translationalKineticEnergy = 2.0;
  energy.rotationalKineticEnergy = 0.5;
  energy.NoseHooverEnergy = 0.25;
  energy.thermobarostatEnergy = 0.125;
  RunningEnergy recomputed = energy;
  recomputed.thermobarostatEnergy = 0.0;

  const std::string difference = energy.printMDDiff(recomputed);
  EXPECT_NE(difference.find("Thermobarostat enthalpy"), std::string::npos);
  EXPECT_NE(difference.find("Conserved energy"), std::string::npos);
  EXPECT_NE(energy.printMD().find("Thermobarostat enthalpy"), std::string::npos);
  EXPECT_NE(energy.repr().find("Thermobarostat enthalpy"), std::string::npos);

  const auto json = energy.jsonMD();
  EXPECT_TRUE(json.contains("Conserved energy [K]"));
  EXPECT_TRUE(json.contains("Baro/thermostat energy [K]"));
  EXPECT_TRUE(json.contains("framework-molecule Real [K]"));
  EXPECT_NEAR(json["Conserved energy [K]"].get<double>(), Units::EnergyToKelvin * energy.conservedEnergy(), 1e-9);

  IntegratorsCPUTime time{};
  time.noSquishFreeRotorOrderTwo = std::chrono::duration<double>(3.0);
  time.noSquishRotate = std::chrono::duration<double>(2.0);
  EXPECT_DOUBLE_EQ(time.total().count(), 3.0);
}
