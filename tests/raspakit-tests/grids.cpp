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
import interactions_ewald;
import energy_status;
import interpolation_energy_grid;
#if !(defined(__has_include) && __has_include(<mdspan>))
  import mdspan;
#endif

TEST(grids, Test_saddle_point)
{
  InterpolationEnergyGrid grid = InterpolationEnergyGrid(int3(1,1,1));

  std::vector<int3> corners = {
    int3(0, 0, 0), int3(1, 0, 0), int3(0, 1, 0), int3(1, 1, 0),
    int3(0, 0, 1), int3(1, 0, 1), int3(0, 1, 1), int3(1, 1, 1)
  };

  std::mdspan<double, std::dextents<size_t, 4>, std::layout_left> data_cell(grid.data.data(), 27, 2, 2, 2);
  for(size_t i = 0 ; i < 8; ++i)
  {
    data_cell[0, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) *
                                                             (corners[i].y-0.5) * (corners[i].y-0.5) *
                                                             (corners[i].z-0.5);

    data_cell[ 1, corners[i].x, corners[i].y, corners[i].z] = 3.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[ 2, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[ 3, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5);

    data_cell[ 4, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[ 5, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[ 6, corners[i].x, corners[i].y, corners[i].z] = 3.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5);
    data_cell[ 7, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].z-0.5);
    data_cell[ 8, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5);
    data_cell[ 9, corners[i].x, corners[i].y, corners[i].z] = 0.0;

    data_cell[10, corners[i].x, corners[i].y, corners[i].z] = 12.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[11, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5);
    data_cell[12, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].z-0.5);
    data_cell[13, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5);
    data_cell[14, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5);
    data_cell[15, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[16, corners[i].x, corners[i].y, corners[i].z] = 0.0;

    data_cell[17, corners[i].x, corners[i].y, corners[i].z] = 12.0 * (corners[i].x-0.5) * (corners[i].z-0.5);
    data_cell[18, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[19, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[20, corners[i].x, corners[i].y, corners[i].z] = 12.0 * (corners[i].x-0.5) * (corners[i].y-0.5);
    data_cell[21, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5);
    data_cell[22, corners[i].x, corners[i].y, corners[i].z] = 0.0;

    data_cell[23, corners[i].x, corners[i].y, corners[i].z] = 12.0 * (corners[i].x-0.5);
    data_cell[24, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[25, corners[i].x, corners[i].y, corners[i].z] = 0.0;

    data_cell[26, corners[i].x, corners[i].y, corners[i].z] = 0.0;

    /*
    data_cell[0, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) *
                                                             (corners[i].y-0.5) * (corners[i].y-0.5) *
                                                             (corners[i].z-0.5);

    data_cell[1, corners[i].x, corners[i].y, corners[i].z] = 3.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[2, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[3, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5);

    data_cell[4, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[5, corners[i].x, corners[i].y, corners[i].z] = 3.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5);
    data_cell[6, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5);

    data_cell[7, corners[i].x, corners[i].y, corners[i].z] = 6.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5);
    */

    /*
    data_cell[0, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) *
                                                             (corners[i].y-0.5) * (corners[i].y-0.5) *
                                                             (corners[i].z-0.5) * (corners[i].z-0.5);

    data_cell[1, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5) * (corners[i].z-0.5);
    data_cell[2, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5) * (corners[i].z-0.5);
    data_cell[3, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);

    data_cell[4, corners[i].x, corners[i].y, corners[i].z] = 4.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5) * (corners[i].z-0.5);
    data_cell[5, corners[i].x, corners[i].y, corners[i].z] = 4.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[6, corners[i].x, corners[i].y, corners[i].z] = 4.0 * (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);

    data_cell[7, corners[i].x, corners[i].y, corners[i].z] = 8.0 * (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);
    */
    /*
    data_cell[0, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].x-0.5) * (corners[i].x-0.5) +
                                                             (corners[i].y-0.5) * (corners[i].y-0.5) *
                                                             (corners[i].z-0.5);

    data_cell[1, corners[i].x, corners[i].y, corners[i].z] = 3.0 * (corners[i].x-0.5) * (corners[i].x-0.5);
    data_cell[2, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[3, corners[i].x, corners[i].y, corners[i].z] = (corners[i].y-0.5) * (corners[i].y-0.5);

    data_cell[4, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[5, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    data_cell[6, corners[i].x, corners[i].y, corners[i].z] = 2.0 * (corners[i].y-0.5);

    data_cell[7, corners[i].x, corners[i].y, corners[i].z] = 0.0;
    */

    /*
    data_cell[0, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].y-0.5) * (corners[i].z-0.5);

    data_cell[1, corners[i].x, corners[i].y, corners[i].z] = (corners[i].y-0.5) * (corners[i].z-0.5);
    data_cell[2, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].z-0.5);
    data_cell[3, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5) * (corners[i].y-0.5);

    data_cell[4, corners[i].x, corners[i].y, corners[i].z] = (corners[i].z-0.5);
    data_cell[5, corners[i].x, corners[i].y, corners[i].z] = (corners[i].y-0.5);
    data_cell[6, corners[i].x, corners[i].y, corners[i].z] = (corners[i].x-0.5);

    data_cell[7, corners[i].x, corners[i].y, corners[i].z] = 1.0;
    */
  }

  //for(size_t i = 0; i < 100; ++i)
  //{
  //  double3 probe_point = double3(drand48(), drand48(), drand48());
  //  //double value = (probe_point.x-0.5) * (probe_point.x-0.5) * (probe_point.x-0.5) + (probe_point.y-0.5) * (probe_point.y-0.5) * (probe_point.z-0.5);
  //  //double value = (probe_point.x-0.5) * (probe_point.x-0.5) * (probe_point.y-0.5) * (probe_point.y-0.5) * (probe_point.z-0.5) * (probe_point.z-0.5);
  //  double value = (probe_point.x-0.5) * (probe_point.x-0.5) * (probe_point.x-0.5) * (probe_point.y-0.5) * (probe_point.y-0.5) * (probe_point.z-0.5);
  //  double inter_value = grid.interpolateVDWGrid(probe_point);
  //  std::print("point: {} {} {} energy: {} interpolated value: {}\n", probe_point.x, probe_point.y, probe_point.z, value, inter_value);
  //}
}

/*
TEST(grids, Test_grid_1x1x1)
{
  ForceField forceField = ForceField(
      {
        PseudoAtom("H", true, 1.0, 1.0, 0.0, 1, false),
      },
      {VDWParameters(1.0, 2.0, VDWParameters::Type::RepulsiveHarmonic)},
      ForceField::MixingRule::Lorentz_Berthelot, 2.0, 2.0, 2.0, false, false, true);

  Framework f = Framework(
      0, forceField, "Box", SimulationBox(1.0, 1.0, 1.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
       Atom(double3(0.5, 0.5, 0.5), 0.0, 1.0, 0, 0, 0, 0)},
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 0.0, 0.0, 0.0,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 0, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);


  InterpolationEnergyGrid grid = InterpolationEnergyGrid(int3(1,1,1));
  grid.makeVDWGrid(forceField, system.frameworkComponents.front(), 0);

  for(size_t i = 0; i < grid.data.size(); ++i)
  {
    std::print("output {} {}\n", i, grid.data[i] * Units::EnergyToKelvin);
  }
  std::print("\n");

  double value = grid.interpolateVDWGrid(double3(0.0, 0.0, 0.0));
  std::print("interpolated value: {}\n", value * Units::EnergyToKelvin);
  std::print("interpolated value 0.5, 0.5, 0.5: {}\n", grid.interpolateVDWGrid(double3(0.5, 0.5, 0.5))  * Units::EnergyToKelvin);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  //std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  //atomPositions[0].position = double3(0.15, 0.05, 0.25);

  for(size_t i=0;i<=20;i++)
  {
    double sx = static_cast<double>(i) / 20.0;
    double sy = 0.5;
    double sz = 0.5;

  auto [energy, _, _, _, _, _, _] =
          Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, double3(sx, sy, sz), 0, frameworkAtoms);

  double inter_value = grid.interpolateVDWGrid(double3(sx, sy, sz));

  std::print("{} {} {}\n", sx, energy *  Units::EnergyToKelvin, inter_value * Units::EnergyToKelvin);
  }
}
*/
