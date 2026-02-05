module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <string>
#include <vector>
#endif

module mc_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import atom;
import randomnumbers;
import skspacegroupdatabase;
import framework;
import forcefield;
import units;

void MC_SurfaceArea::run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
                         std::string probePseudoAtom, std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps) const
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_iterations = numberOfIterations.value_or(1000);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);

  double accumulated_surface_area{};
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    for (std::make_signed_t<std::size_t> atom_index = 0; const Atom &atom : framework.unitCellAtoms)
    {
      std::size_t atomType = static_cast<std::size_t>(atom.type);
      double size_parameter = forceField(probeType.value(), atomType).sizeParameter();
      double equilibrium_distance = wellDepthFactor * size_parameter;

      double total_trials{};
      double counted{};
      for (std::size_t j = 0; j < number_of_inner_steps; ++j)
      {
        double3 vec = random.randomVectorOnUnitSphere();

        double3 position = atom.position + equilibrium_distance * vec;
        if (!framework.computeOverlap(forceField, position, wellDepthFactor, probeType.value(), atom_index))
        {
          counted += 1.0;
        }

        total_trials += 1.0;
      }

      double temp = (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
      accumulated_surface_area += temp;

      ++atom_index;
    }
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.sa.cpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
  std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", probePseudoAtom, wellDepthFactor, forceField[probeType.value()].sizeParameter());
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps: {}\n", number_of_inner_steps);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area / static_cast<double>(number_of_iterations) << " [A^2]" << std::endl;
  myfile << (accumulated_surface_area / static_cast<double>(number_of_iterations)) * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.unitCellMass
         << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * (accumulated_surface_area / static_cast<double>(number_of_iterations)) / framework.simulationBox.volume << " [m^2/cm^3]" << std::endl;

  myfile.close();
}
