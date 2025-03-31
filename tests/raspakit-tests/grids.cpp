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
        PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
        PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
        PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
        PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
        PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
    },
    {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
     VDWParameters(85.671, 3.017)},
    ForceField::MixingRule::Lorentz_Berthelot, 21.8, 21.8, 21.8, false, false, true);


   forceField.hasExternalField = true;
   //forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction;
   forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction;
   forceField.potentialEnergySurfaceOrigin = double3(-3.0, -3.0, -3.0);

  /*
    Framework f = Framework(
      //0, forceField, "CHA", SimulationBox(9.459, 9.459, 9.459,
      0, forceField, "CHA", SimulationBox(3.0, 3.0, 3.0,
        90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, 0)
      },
      int3(1, 1, 1));
      */
    Framework f = Framework(
      0, forceField, "CHA", SimulationBox(6.0, 6.0, 6.0,
        90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  size_t typeB = 2;
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  int3 numberOfGridPoints{8, 8, 8};

  //InterpolationEnergyGrid grid = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Tricubic);
  InterpolationEnergyGrid grid = InterpolationEnergyGrid(numberOfGridPoints, InterpolationEnergyGrid::InterpolationOrder::Triquintic);
  grid.makeVDWGrid(system.forceField, system.frameworkComponents.front(), typeB);

  double highest_error = 0.0;
  double value = 0.0;
  double3 pos_value{};
  int3 N{101, 101, 101};
  double3 delta = double3(1.0 / static_cast<double>(N.x - 1),
                          1.0 / static_cast<double>(N.y - 1),
                          1.0 / static_cast<double>(N.z - 1));
  for(size_t i = 0; i < static_cast<size_t>(N.z); ++i)
  {
    for(size_t j = 0; j < static_cast<size_t>(N.y); ++j)
    {
      for(size_t k = 0; k < static_cast<size_t>(N.z); ++k)
      {
        double3 s = double3(static_cast<double>(i) * delta.x,
                            static_cast<double>(j) * delta.y,
                            static_cast<double>(k) * delta.z);
        double3 posB = system.simulationBox.cell * s + forceField.potentialEnergySurfaceOrigin;

        //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
        //std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
        //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
        std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
        double interpolated_value = grid.interpolateVDWGrid(s);

        double error = std::fabs(analytical[0] - interpolated_value);
        if(error > highest_error) 
        {
          highest_error = error;
          value = analytical[0];
          pos_value = posB + forceField.potentialEnergySurfaceOrigin;
        }
      }
    }
  }
  std::print("highest error: {} value: {} at: {} {} {}\n", highest_error, value, pos_value.x, pos_value.y, pos_value.z);
/*
  //double3 posB = double3(5.1, 4.25, 5.4);
  double3 posB = double3(-2.0, -2.0, -1.0);
  double3 s = (system.simulationBox.inverseCell * (posB - forceField.potentialEnergySurfaceOrigin)).fract();
  //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
  //std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
  //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
  std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
  double interpolated_value = grid.interpolateVDWGrid(s);
  std::print("compare: {} to interpolated value {}\n", analytical[0], interpolated_value);

  for(size_t i = 0; i < static_cast<size_t>(N.z); ++i)
  {
    for(size_t j = 0; j < static_cast<size_t>(N.y); ++j)
    {
      double3 s = double3(static_cast<double>(i) * delta.x,
                          static_cast<double>(j) * delta.y,
                          static_cast<double>(50) * delta.z);
      double3 posB = system.simulationBox.cell * s + forceField.potentialEnergySurfaceOrigin;

      //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
      //std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
      //std::array<double, 8> analytical = Interactions::calculateTricubicFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
      std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPositionExternalField(system.forceField, system.simulationBox, posB);
      double interpolated_value = grid.interpolateVDWGrid(s);

      std::print("{} {} {}\n", posB.x, posB.y, interpolated_value);
    }
    std::print("\n");
  }
*/
}
