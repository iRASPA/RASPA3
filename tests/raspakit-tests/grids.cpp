#include <gtest/gtest.h>

#include <cstddef>
#include <algorithm>
#include <complex>
#include <span>
#include <vector>
#include <numbers>
#include <print>
#if defined(__has_include) && __has_include(<mdspan>)
  #include <mdspan>
#endif

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
import gradient_factor;
import hessian_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_external_field;
import interactions_ewald;
import energy_status;
import interpolation_energy_grid;
#if !(defined(__has_include) && __has_include(<mdspan>))
  import mdspan;
#endif

TEST(grids, Test_CHA_grid)
{

  ForceField forceField = ForceField(
    {
      PseudoAtom("Zn", true, 65.37,    1.275, 0.0, 30, false),
      PseudoAtom("O",  true, 15.9994, -1.5,   0.0, 8,  false),
      PseudoAtom("O",  true, 15.9994, -0.6,   0.0, 8,  false),
      PseudoAtom("C1", true, 12.0107, -0.475, 0.0, 6,  false),
      PseudoAtom("C2", true, 12.0107, -0.125, 0.0, 6,  false),
      PseudoAtom("C3", true, 12.0107, -0.15,  0.0, 6,  false),
      PseudoAtom("H",  true, 1.00794,  0.15,  0.0, 1,  false),
      PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
      PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
      PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
    },
    {
      VDWParameters(0.42, 2.7), 
      VDWParameters(700.0, 2.98), 
      VDWParameters(70.5,  3.11), 
      VDWParameters(47.0,  3.74), 
      VDWParameters(47.86, 3.47), 
      VDWParameters(47.86, 3.47), 
      VDWParameters(7.65, 2.85), 
      VDWParameters(158.5, 3.72), 
      VDWParameters(29.933, 2.745),
      VDWParameters(85.671, 3.017)
    },
    ForceField::MixingRule::Jorgensen, 12.0, 12.0, 12.0, false, false, true);

  std::cout << forceField.printPseudoAtomStatus();
  std::cout <<  forceField.printForceFieldStatus();


    Framework f = Framework(
      0, forceField, "IRMOF-1", SimulationBox(25.832, 25.832, 25.832,
        90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0), 225,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.2934, 0.2066, 0.2066),  1.275,  1.0, 0, 0, 0, 0),
        Atom(double3(0.25,   0.25,   0.25),   -1.5,    1.0, 0, 1, 0, 0),
        Atom(double3(0.2819, 0.2181, 0.134),   0.6,    1.0, 0, 2, 0, 0),
        Atom(double3(0.25,   0.25,   0.1113), -0.475,  1.0, 0, 3, 0, 0),
        Atom(double3(0.25,   0.25,   0.0538), -0.1255, 1.0, 0, 4, 0, 0),
        Atom(double3(0.2829, 0.2171, 0.0269), -0.15,   1.0, 0, 5, 0, 0),
        Atom(double3(0.3049, 0.1951, 0.0448),  0.15,   1.0, 0, 6, 0, 0),
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {1}, 5);

  size_t typeB = forceField.pseudoAtoms.size() - 1;
  //size_t typeB = 7;
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  int3 numberOfGridPoints{259, 259, 259};

  InterpolationEnergyGrid grid = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  //InterpolationEnergyGrid grid = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  grid.makeInterpolationGrid(ForceField::InterpolationGridType::LennardJones, system.forceField, system.frameworkComponents.front(), 7);

  double3 pos = double3(8,7,14);
  double3 s = system.simulationBox.inverseCell * pos;

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = pos;
  RunningEnergy energy =
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, frameworkAtoms, atomPositions);

  double interpolated_value = grid.interpolateVDWGrid(s);
  std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPosition(ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, pos, 7, frameworkAtoms);

  std::print("energy: {}  grid analytical: {} interpolated: {}\n", energy.frameworkMoleculeVDW, analytical[0], interpolated_value);

  InterpolationEnergyGrid grid_repulsion = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  //InterpolationEnergyGrid grid_repulsion = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  grid_repulsion.makeInterpolationGrid(ForceField::InterpolationGridType::LennardJonesRepulsion, system.forceField, system.frameworkComponents.front(), typeB);

  InterpolationEnergyGrid grid_attraction = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  //InterpolationEnergyGrid grid_attraction = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  grid_attraction.makeInterpolationGrid(ForceField::InterpolationGridType::LennardJonesAttraction, system.forceField, system.frameworkComponents.front(), typeB);



  double interpolated_repulsion_value = grid_repulsion.interpolateVDWGrid(s);
  double interpolated_attraction_value = grid_attraction.interpolateVDWGrid(s);
  std::print("{} {} {}\n", interpolated_repulsion_value, interpolated_attraction_value, interpolated_repulsion_value - interpolated_attraction_value);
}
