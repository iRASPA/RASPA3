#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <numbers>
#include <print>
#include <span>
#include <sstream>
#include <vector>
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
import interactions_framework_molecule_grid;
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
          PseudoAtom("Zn", true, 65.37, 1.275, 0.0, 30, false),
          PseudoAtom("O1", true, 15.9994, -1.5, 0.0, 8, false),
          PseudoAtom("O2", true, 15.9994, -0.6, 0.0, 8, false),
          PseudoAtom("C1", true, 12.0107, -0.475, 0.0, 6, false),
          PseudoAtom("C2", true, 12.0107, -0.125, 0.0, 6, false),
          PseudoAtom("C3", true, 12.0107, -0.15, 0.0, 6, false),
          PseudoAtom("H1", true, 1.00794, 0.15, 0.0, 1, false),
          PseudoAtom("CH4", false, 16.04246, 1.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(0.42, 2.7), VDWParameters(700.0, 2.98), VDWParameters(70.5, 3.11), VDWParameters(47.0, 3.74),
       VDWParameters(47.86, 3.47), VDWParameters(47.86, 3.47), VDWParameters(7.65, 2.85),
       // VDWParameters(158.5, 3.72),
       // VDWParameters(29.933, 2.745),
       // VDWParameters(85.671, 3.017)
       VDWParameters(0.0, 3.72), VDWParameters(0.0, 2.745), VDWParameters(0.0, 3.017)},
      ForceField::MixingRule::Jorgensen, 12.0, 12.0, 12.0, false, false, true);

  std::cout << forceField.printPseudoAtomStatus();
  std::cout << forceField.printForceFieldStatus();

  Framework framework = Framework(0, forceField, "IRMOF-1",
                                  SimulationBox(25.832, 25.832, 25.832, 90.0 * std::numbers::pi / 180.0,
                                                120.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0),
                                  523,
                                  {
                                      // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t
                                      // type, uint8_t componentId, uint8_t groupId
                                      Atom(double3(0.2934, 0.2066, 0.2066), 1.275, 1.0, 0, 0, 0, 0),
                                      Atom(double3(0.25, 0.25, 0.25), -1.5, 1.0, 0, 1, 0, 0),
                                      Atom(double3(0.2819, 0.2181, 0.134), 0.6, 1.0, 0, 2, 0, 0),
                                      Atom(double3(0.25, 0.25, 0.1113), -0.475, 1.0, 0, 3, 0, 0),
                                      Atom(double3(0.25, 0.25, 0.0538), -0.1255, 1.0, 0, 4, 0, 0),
                                      Atom(double3(0.2829, 0.2171, 0.0269), -0.15, 1.0, 0, 5, 0, 0),
                                      Atom(double3(0.3049, 0.1951, 0.0448), 0.15, 1.0, 0, 6, 0, 0),
                                  },
                                  int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {framework}, {c}, {1}, 5);

  // size_t typeB = forceField.pseudoAtoms.size() - 1;
  // size_t typeB = 7;
  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  std::print("atoms: {} {}\n", frameworkAtoms.size(), framework.atoms.size());

  double3 pos = double3(11, 12, 13);

  std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
      ForceField::InterpolationGridType::EwaldReal, system.forceField, system.simulationBox, pos, 7,
      framework.simulationBox, frameworkAtoms);

  double3 cartesian_force =
      framework.simulationBox.inverseCell.transpose() * double3(analytical[1], analytical[2], analytical[3]);
  double3x3 analytical_hessian = double3x3(analytical[4], analytical[5], analytical[6], analytical[5], analytical[7],
                                           analytical[8], analytical[6], analytical[8], analytical[9]);

  double3x3 cartesian_hessian =
      framework.simulationBox.inverseCell.transpose() * analytical_hessian * framework.simulationBox.inverseCell;

  std::print("analytical: {} gradient: {} {} {}\n", analytical[0], cartesian_force.x, cartesian_force.y,
             cartesian_force.z);
  std::print("hessian:  {} {} {} {} {} {}\n", analytical[4], analytical[5], analytical[6], analytical[7], analytical[8],
             analytical[9]);
  std::print("CORRECT cartesian hessian:  {} {} {} {} {} {}\n", cartesian_hessian.ax, cartesian_hessian.ay,
             cartesian_hessian.az, cartesian_hessian.by, cartesian_hessian.bz, cartesian_hessian.cz);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = pos;
  atomPositions[0].charge = 1.0;
  RunningEnergy energy = Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                                        frameworkAtoms, atomPositions);
  auto [l1, l2, hessian] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, pos,
                                                                           1.0, frameworkAtoms);

  std::print("energy: {}  gradient: {} {} {}\n", energy.frameworkMoleculeCharge, atomPositions[0].gradient.x,
             atomPositions[0].gradient.y, atomPositions[0].gradient.z);
  std::print("CORRECT analytical hessian: {} {} {} {} {} {}\n", hessian.ax, hessian.ay, hessian.az, hessian.by,
             hessian.bz, hessian.cz);

  int3 numberOfGridPoints{32, 32, 32};
  InterpolationEnergyGrid grid =
      InterpolationEnergyGrid(framework.simulationBox, numberOfGridPoints, ForceField::InterpolationScheme::Triquintic);
  // InterpolationEnergyGrid grid = InterpolationEnergyGrid(numberOfGridPoints,
  // InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  std::stringstream stream;
  grid.makeInterpolationGrid(stream, ForceField::InterpolationGridType::EwaldReal, system.forceField,
                             system.framework.value(), 7);

  // RunningEnergy energy =
  //     Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox,
  //     system.interpolationGrids,
  //                                                  system.framework, frameworkAtoms, atomPositions);

  auto [interpolated_value, interpolated_gradient, interpolated_hessian] = grid.interpolateHessian(pos);

  std::print("energy: {}  grid analytical: {} interpolated: {} {} {} {}\n", energy.frameworkMoleculeVDW, analytical[0],
             interpolated_value, interpolated_gradient.x, interpolated_gradient.y, interpolated_gradient.z);
  std::print("hessian: {} {} {} {} {} {}\n", interpolated_hessian.ax, interpolated_hessian.ay, interpolated_hessian.az,
             interpolated_hessian.by, interpolated_hessian.cy, interpolated_hessian.cz);

  // InterpolationEnergyGrid grid_repulsion = InterpolationEnergyGrid(numberOfGridPoints,
  // InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  ////InterpolationEnergyGrid grid_repulsion = InterpolationEnergyGrid(numberOfGridPoints,
  /// InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  // grid_repulsion.makeInterpolationGrid(ForceField::InterpolationGridType::LennardJonesRepulsion, system.forceField,
  // system.framework.value(), typeB);

  // InterpolationEnergyGrid grid_attraction = InterpolationEnergyGrid(numberOfGridPoints,
  // InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  ////InterpolationEnergyGrid grid_attraction = InterpolationEnergyGrid(numberOfGridPoints,
  /// InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  // grid_attraction.makeInterpolationGrid(ForceField::InterpolationGridType::LennardJonesAttraction, system.forceField,
  // system.framework.value(), typeB);

  // double interpolated_repulsion_value = grid_repulsion.interpolateVDWGrid(s);
  // double interpolated_attraction_value = grid_attraction.interpolateVDWGrid(s);
  // std::print("{} {} {}\n", interpolated_repulsion_value, interpolated_attraction_value, interpolated_repulsion_value
  // - interpolated_attraction_value);
}
