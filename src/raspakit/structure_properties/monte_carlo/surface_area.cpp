module;

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

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import atom;
import randomnumbers;
import framework;
import forcefield;
import units;

void MC_SurfaceArea::run(const ForceField &forceField, const Framework &framework, double well_depth_factor,
                         std::string probe_pseudo_atom, std::size_t number_of_iterations) const
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probe_pseudo_atom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  double accumulated_surface_area{};
  for (std::make_signed_t<std::size_t> atom_index = 0; const Atom &atom : framework.unitCellAtoms)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    double size_parameter = forceField(probeType.value(), atomType).sizeParameter();
    double equilibrium_distance = well_depth_factor * size_parameter;

    double total_trials{};
    double counted{};
    for (std::size_t i = 0; i < number_of_iterations; ++i)
    {
      double3 vec = random.randomVectorOnUnitSphere();

      double3 position = atom.position + equilibrium_distance * vec;
      if (!framework.computeOverlap(forceField, position, well_depth_factor, probeType.value(), atom_index))
      {
        counted += 1.0;
      }

      total_trials += 1.0;
    }

    double temp = (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
    accumulated_surface_area += temp;

    ++atom_index;
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.sa.cpu.txt");
  std::print(myfile, "# Surface area using Mont Carlo-based method\n");
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area << " [A^2]" << std::endl;
  myfile << accumulated_surface_area * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.unitCellMass
         << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * accumulated_surface_area / framework.simulationBox.volume << " [m^2/cm^3]" << std::endl;

  myfile.close();
}
