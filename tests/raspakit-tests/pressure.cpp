#include <gtest/gtest.h>

#include <cstddef>
#include <vector>
#include <tuple>
#include <algorithm>
#include <span>
#include <complex>

import int3;
import double3;
import double3x3;

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
import interactions_ewald;

TEST(MC_strain_tensor, Test_20_CH4_25x25x25_LJ)
{
  double delta = 1e-7;
  double tolerance = 1e-4;

  ForceField forceField = ForceField(
    { PseudoAtom("CH4", 16.04246, 0.0, 6, false) },
    { VDWParameters(158.5, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { 
      // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 0, 0, 0) 
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { c }, { 20 }, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = 
    Interactions::computeInterMolecularEnergyStrainDerivative(system.forceField, system.components, system.simulationBox, moleculeAtomPositions);

  std::vector<std::pair<double3x3, double>> strains{
     std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.ax},
     std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.bx},
     std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.cx},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.ay},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.by},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0} }, pressureInfo.second.cy},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0} }, pressureInfo.second.az},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0} }, pressureInfo.second.bz},
     std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta} }, pressureInfo.second.cz}
  };

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{ double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0} };

  for(const std::pair<double3x3, double> strain: strains)
  {
    SimulationBox strainBox_forward2 = SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_forward2),
      [&strainBox_forward2, &inv](const Atom &m) { return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyForward2;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward2, moleculeAtomPositions_forward2, EnergyForward2);

    SimulationBox strainBox_forward1 = SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_forward1),
      [&strainBox_forward1, &inv](const Atom& m) { return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyForward1;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward1, moleculeAtomPositions_forward1, EnergyForward1);
    
    
    SimulationBox strainBox_backward1 = SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_backward1),
      [&strainBox_backward1, &inv](const Atom& m) { return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyBackward1;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward1, moleculeAtomPositions_backward1, EnergyBackward1);

    SimulationBox strainBox_backward2 = SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_backward2),
      [&strainBox_backward2, &inv](const Atom& m) { return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyBackward2;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward2, moleculeAtomPositions_backward2, EnergyBackward2);
    
    double strainDerivative = (-EnergyForward2.total() + 8.0 * EnergyForward1.total() - 8.0 * EnergyBackward1.total() + EnergyBackward2.total()) / (6.0 * delta);
    
    EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }
}

TEST(MC_strain_tensor, Test_20_Na_Cl_25x25x25_LJ_Real)
{
  double delta = 1e-7;
  double tolerance = 1e-4;

  ForceField forceField = ForceField(
    { PseudoAtom("Si",    28.0855,   2.05,  14, false),
      PseudoAtom("O",     15.999,   -1.025,  8, false),
      PseudoAtom("CH4",   16.04246,  0.0,    6, false),
      PseudoAtom("Na+",   12.0,      1.0,    6, false),
      PseudoAtom("Cl-",   15.9994,  -1.0,    8, false),
    },
    { VDWParameters(22.0, 2.30),
      VDWParameters(53.0, 3.3),
      VDWParameters(158.5, 3.72),
      VDWParameters(15.0966, 2.65755),
      VDWParameters(142.562, 3.51932)
    },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
    {
      // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0, 0.0), 1.0, 1.0, 0, 3, 0, 0),
    }, 5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
    {
      // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0, 0.0), -1.0, 1.0, 0, 4, 1, 0),
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { na, cl }, { 1, 1 }, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  std::pair<EnergyStatus, double3x3> pressureInfo = 
    Interactions::computeInterMolecularEnergyStrainDerivative(system.forceField, system.components, system.simulationBox, moleculeAtomPositions);
  pressureInfo.first.sumTotal();

  std::vector<std::pair<double3x3, double>> strains{
     std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0} },   pressureInfo.second.ax},
     std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0} },   pressureInfo.second.bx},
     std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0} },   pressureInfo.second.cx},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0} },   pressureInfo.second.ay},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0} },   pressureInfo.second.by},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0} },   pressureInfo.second.cy},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0},   double3{delta, 0.0, 0.0} }, pressureInfo.second.az},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0},   double3{0.0, delta, 0.0} }, pressureInfo.second.bz},
     std::pair{double3x3{double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, 0.0},   double3{0.0, 0.0, delta} }, pressureInfo.second.cz}
  };

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{ double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0} };

  for(const std::pair<double3x3, double> strain: strains)
  {
    SimulationBox strainBox_forward2 = SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_forward2),
      [&strainBox_forward2, &inv](const Atom &m) { return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyForward2;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward2, moleculeAtomPositions_forward2, EnergyForward2);

    SimulationBox strainBox_forward1 = SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_forward1),
      [&strainBox_forward1, &inv](const Atom& m) { return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyForward1;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward1, moleculeAtomPositions_forward1, EnergyForward1);
    
    
    SimulationBox strainBox_backward1 = SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_backward1),
      [&strainBox_backward1, &inv](const Atom& m) { return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyBackward1;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward1, moleculeAtomPositions_backward1, EnergyBackward1);

    SimulationBox strainBox_backward2 = SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), std::back_inserter(moleculeAtomPositions_backward2),
      [&strainBox_backward2, &inv](const Atom& m) { return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId, m.groupId); });
    RunningEnergy EnergyBackward2;
    Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward2, moleculeAtomPositions_backward2, EnergyBackward2);
    
    double strainDerivative = (-EnergyForward2.total() + 8.0 * EnergyForward1.total() - 8.0 * EnergyBackward1.total() + EnergyBackward2.total()) / (6.0 * delta);
    
    EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }
}
