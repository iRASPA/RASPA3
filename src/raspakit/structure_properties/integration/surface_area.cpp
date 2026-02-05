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

module integration_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import double3x3;
import atom;
import randomnumbers;
import framework;
import forcefield;
import units;

void Integration_SurfaceArea::run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
                         std::string probePseudoAtom, std::optional<std::size_t> numberOfSlices) const
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("Integration_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_slices = numberOfSlices.value_or(1024);

  double accumulated_surface_area{};
  for (std::make_signed_t<std::size_t> atom_index = 0; const Atom &atom : framework.unitCellAtoms)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    double size_parameter = forceField(probeType.value(), atomType).sizeParameter();
    double equilibrium_distance = wellDepthFactor * size_parameter;

    double total_trials{};
    double counted{};
    for(std::size_t stackNumber = 0; stackNumber <= number_of_slices; ++stackNumber)
    {
      for(std::size_t sliceNumber = 0; sliceNumber < number_of_slices; ++sliceNumber)
      {
        double theta = static_cast<double>(sliceNumber) * 2.0 * std::numbers::pi / static_cast<double>(number_of_slices);

        double u = static_cast<double>(stackNumber) / static_cast<double>(number_of_slices);
        double phi = std::acos(2.0 * u - 1.0);

        double sinTheta = std::sin(theta);
        double sinPhi = std::sin(phi);
        double cosTheta = std::cos(theta);
        double cosPhi = std::cos(phi);
        double3 unit_vector{sinPhi * cosTheta, sinPhi * sinTheta, cosPhi};

        double3 position = atom.position + equilibrium_distance * unit_vector;
        if (!framework.computeOverlap(forceField, position, wellDepthFactor, probeType.value(), atom_index))
        {
          counted += 1.0;
        }

        total_trials += 1.0;
      }
    }

    double temp = (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
    accumulated_surface_area += temp;

    ++atom_index;
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".integration.sa.cpu.txt");
  std::print(myfile, "# Surface area using integration method\n");
  std::print(myfile, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", probePseudoAtom, wellDepthFactor, forceField[probeType.value()].sizeParameter());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Number of integration points per atom: {}\n", (number_of_slices + 1) * number_of_slices);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area << " [A^2]" << std::endl;
  myfile << accumulated_surface_area * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.unitCellMass
         << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * accumulated_surface_area / framework.simulationBox.volume << " [m^2/cm^3]" << std::endl;

  myfile.close();
}
