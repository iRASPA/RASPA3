#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_status;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import integrators;
import integrators_update;
import integrators_compute;
import interpolation_energy_grid;
import cif_reader;
import mc_moves_probabilities;

namespace
{
std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

std::string readLocalFile(const std::filesystem::path &path)
{
  std::ifstream stream(path);
  if (!stream)
  {
    throw std::runtime_error(std::format("Could not read '{}'", path.string()));
  }
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

void printEnergyDecomposition(const RunningEnergy &energy, const std::string &label)
{
  const double conv = Units::EnergyToKelvin;
  const auto toKelvin = [&](double internal) { return conv * internal; };
  const auto toKJmol = [&](double internal) { return conv * internal / 1.2027242847; };
  std::print("=== {} ===\n", label);
  std::print("  Total:                 {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.potentialEnergy()),
             toKelvin(energy.potentialEnergy()));
  std::print("  Host/ads VDW:          {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.frameworkMoleculeVDW),
             toKelvin(energy.frameworkMoleculeVDW));
  std::print("  Host/ads Real:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.frameworkMoleculeCharge),
             toKelvin(energy.frameworkMoleculeCharge));
  std::print("  Ewald Fourier:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_fourier),
             toKelvin(energy.ewald_fourier));
  std::print("  Ewald self:            {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_self),
             toKelvin(energy.ewald_self));
  std::print("  Ewald exclusion:       {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_exclusion),
             toKelvin(energy.ewald_exclusion));
  std::print("  Total VDW:             {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.VanDerWaalsEnergy()),
             toKelvin(energy.VanDerWaalsEnergy()));
  std::print("  Total Coulomb:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.CoulombEnergy()),
             toKelvin(energy.CoulombEnergy()));
  std::print("\n");
}
}  // namespace

TEST(energy_decomposition, CO2_in_IRMOF1_OMS_geometry)
{
  const std::filesystem::path exampleDir =
      repositoryRoot() / "examples/basic/19_minimization_co2_in_irmof_1";

  ForceField forceField = ForceField::readForceField(exampleDir.string(), "force_field.json").value();
  forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
  forceField.useCharge = true;

  const std::string cifContent = readLocalFile(exampleDir / "IRMOF-1.cif");
  const auto cif = CIFReader::readCIFString(cifContent, forceField, CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(cif.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, fractionalAtomsUnitCell] = cif.value();
  Framework framework = Framework(forceField, "IRMOF-1", simulationBox, spaceGroupHallNumber, definedAtoms,
                                  fractionalAtomsUnitCell, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  Component co2(Component::Type::Adsorbate, 0, forceField, "CO2", (exampleDir / "CO2").string(), 5, 21, probabilities,
                std::nullopt, false);

  auto evaluateAtPositions = [&](const std::string &label, const std::array<double3, 3> &positions)
  {
    System system = System(forceField, std::nullopt, false, 300.0, 0.0, 0.81, framework, {co2}, {}, {1}, 5);
    std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
    for (std::size_t atom = 0; atom < positions.size(); ++atom)
    {
      moleculeAtoms[atom].position = positions[atom];
    }
    RunningEnergy energy = system.computeTotalEnergies();
    printEnergyDecomposition(energy, label);
    return energy;
  };

  // These positions were recorded against the old P1 IRMOF-1.cif. The current symmetry-based
  // CIF describes the same crystal in a setting shifted by c/2, so the z coordinates carry a
  // +12.916 Angstrom (half unit cell) translation relative to the original values.
  const std::array<double3, 3> raspa2Oms = {double3(10.136569193677, 10.136569193677, 20.335907533962),
                                            double3(9.417951022032, 9.417951022032, 20.871982945414),
                                            double3(8.699332850387, 8.699332850387, 21.408058356865)};
  const std::array<double3, 3> raspa3Pore = {double3(9.4276655408314216, 17.808185239505558, 3.488334459168581),
                                             double3(9.919740710490089, 18.722474610990751, 2.996259289509912),
                                             double3(10.411815880148756, 19.636763982475944, 2.504184119851244)};

  const RunningEnergy omsEnergy = evaluateAtPositions("RASPA3 at RASPA2 OMS geometry", raspa2Oms);
  evaluateAtPositions("RASPA3 at RASPA3 pore geometry (seed 12345)", raspa3Pore);
  EXPECT_NEAR(omsEnergy.potentialEnergy() * Units::EnergyToKelvin, -2701.54, 5.0);
}

TEST(energy_decomposition, CO2_Methane_in_Box)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane, co2}, {}, {15, 30}, 5);

  RunningEnergy energy = system.computeTotalEnergies();

  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.atomData, system.atomDynamics);
  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge, strainDerivative.first.totalEnergy.energy,
              1e-6);
}

TEST(energy_decomposition, CO2_Methane_in_Box_Ewald)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);
  Component co2_2 = Component::makeCO2(forceField, 2, true);

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2, co2_2}, {}, {15, 30}, 5);

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  RunningEnergy energy = system.computeTotalEnergies();

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.forceField, system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
      system.netChargeFramework, system.netChargePerComponent);

  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion,
              strainDerivative.first.totalEnergy.energy, 1e-6);
}

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3> &lhs,
                                                   const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

TEST(energy_decomposition, CO2_Methane_in_Framework)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {10, 15}, 5);

  system.precomputeTotalRigidEnergy();

  RunningEnergy energy = system.computeTotalEnergies();

  RunningEnergy energyForces = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox,
      system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);

  std::pair<EnergyStatus, double3x3> strainDerivative = system.computeMolecularPressure();

  EXPECT_NEAR(energy.potentialEnergy(), strainDerivative.first.totalEnergy.energy, 1e-6);
  EXPECT_NEAR(energy.potentialEnergy(), energyForces.potentialEnergy(), 1e-6);
}
