#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <print>
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
import energy_factor;
import energy_status;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_internal;
import connectivity_table;
import interactions_ewald;
import intra_molecular_potentials;
import bond_potential;

TEST(MC_intramolecular_strain_tensor, Test_20_ethane_25x25x25_harmonic_bond_potential)
{
  double delta = 1e-5;
  double tolerance = 1e-4;

  ForceField forceField = ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 8, false},
                                      {"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                                      {"CH2", false, 14.02658, 0.0, 0.0, 8, false},
                                      {"CH", false, 13.01864, 0.0, 0.0, 8, false},
                                      {"C", false, 12.0, 0.0, 0.0, 8, false}},
                                     {{158.5, 3.72}, {108.0, 3.76}, {56.0, 3.96}, {17.0, 4.67}, {0.8, 6.38}},
                                     ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);

  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component c = Component(
      0, forceField, "ethane", 305.33, 4871800.0, 0.0993,
      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false), Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
      connectivityTable, intraMolecularPotentials, 5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {}, {20}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<double, double3x3> pressureInfo = Interactions::computeIntraMolecularBondStrainDerivative(
      system.components[0].intraMolecularPotentials, system.moleculeData, moleculeAtomPositions);

  double3 gradient{};
  for (size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    double3 saved_position = moleculeAtomPositions[i].position;

    // finite difference x
    moleculeAtomPositions[i].position.x = saved_position.x + 0.5 * delta;
    x2 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);

    moleculeAtomPositions[i].position.x = saved_position.x - 0.5 * delta;
    x1 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);
    moleculeAtomPositions[i].position.x = saved_position.x;

    // finite difference y
    moleculeAtomPositions[i].position.y = saved_position.y + 0.5 * delta;
    y2 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);

    moleculeAtomPositions[i].position.y = saved_position.y - 0.5 * delta;
    y1 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);
    moleculeAtomPositions[i].position.y = saved_position.y;

    // finite difference z
    moleculeAtomPositions[i].position.z = saved_position.z + 0.5 * delta;
    z2 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);

    moleculeAtomPositions[i].position.z = saved_position.z - 0.5 * delta;
    z1 = Interactions::computeIntraMolecularBondEnergy(system.components[0].intraMolecularPotentials,
                                                       system.moleculeData, moleculeAtomPositions);
    moleculeAtomPositions[i].position.z = saved_position.z;

    gradient.x = (x2.bond - x1.bond) / delta;
    gradient.y = (y2.bond - y1.bond) / delta;
    gradient.z = (z2.bond - z1.bond) / delta;

    EXPECT_NEAR(moleculeAtomPositions[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(moleculeAtomPositions[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(moleculeAtomPositions[i].gradient.z, gradient.z, tolerance);
  }

  std::vector<std::pair<double3x3, double>> strains{
      std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ax},
      std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.bx},
      std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cx},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ay},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.by},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cy},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}},
                pressureInfo.second.az},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}},
                pressureInfo.second.bz},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}},
                pressureInfo.second.cz}};

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  for (const std::pair<double3x3, double> strain : strains)
  {
    SimulationBox strainBox_forward2 =
        SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward2),
                   [&strainBox_forward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward2 = Interactions::computeIntraMolecularBondEnergy(
        system.components[0].intraMolecularPotentials, system.moleculeData, moleculeAtomPositions_forward2);

    SimulationBox strainBox_forward1 =
        SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward1),
                   [&strainBox_forward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward1 = Interactions::computeIntraMolecularBondEnergy(
        system.components[0].intraMolecularPotentials, system.moleculeData, moleculeAtomPositions_forward1);

    SimulationBox strainBox_backward1 =
        SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward1),
                   [&strainBox_backward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward1 = Interactions::computeIntraMolecularBondEnergy(
        system.components[0].intraMolecularPotentials, system.moleculeData, moleculeAtomPositions_backward1);

    SimulationBox strainBox_backward2 =
        SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward2),
                   [&strainBox_backward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward2 = Interactions::computeIntraMolecularBondEnergy(
        system.components[0].intraMolecularPotentials, system.moleculeData, moleculeAtomPositions_backward2);

    double strainDerivative = (-EnergyForward2.potentialEnergy() + 8.0 * EnergyForward1.potentialEnergy() -
                               8.0 * EnergyBackward1.potentialEnergy() + EnergyBackward2.potentialEnergy()) /
                              (6.0 * delta);

    EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }
}
