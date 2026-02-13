module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#include "mdspanwrapper.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <istream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <optional>
#include "mdspanwrapper.h"
#endif

module energy_void_fraction;

#ifdef USE_STD_IMPORT
import std;
#endif

import opencl;
import float4;
import double2;
import double4;
import randomnumbers;
import skspacegroupdatabase;
import forcefield;
import framework;
import simulationbox;
import atom;
import potential_energy_vdw;
import energy_factor;
import units;
#if !(defined(__has_include) && __has_include(<mdspan>))
//import mdspan;
#endif

void EnergyVoidFraction::run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom, 
                             std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnersteps)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  std::optional<std::size_t> probe_type = forceField.findPseudoAtom(probePseudoAtom);
  if (!probe_type.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_iterations = numberOfIterations.value_or(1000);
  std::size_t number_of_inner_steps = numberOfInnersteps.value_or(1000);

  time_begin = std::chrono::system_clock::now();

  double cutoff = forceField.cutOffFrameworkVDW;
  int3 numberOfReplicas = framework.simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  SimulationBox simulation_box = framework.simulationBox.scaled(numberOfReplicas);
  double3x3 simulation_cell = simulation_box.cell;

  double sumBoltzmannWeight{};
  std::size_t count{};
  for (std::size_t l = 0; l < number_of_iterations; ++l)
  {
    for (std::size_t m = 0; m < number_of_inner_steps; ++m)
    {
      double3 random_position = {random.uniform(), random.uniform(), random.uniform()};

      double energy = 0.0;
      for (const Atom &atom: framework.fractionalUnitCellAtoms)
      {
        std::size_t typeB = atom.type;
        for (std::size_t i = 0; i < static_cast<std::size_t>(numberOfReplicas.x); i++)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfReplicas.y); j++)
          {
            for (std::size_t k = 0; k < static_cast<std::size_t>(numberOfReplicas.z); k++)
            {
              double3 translation(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
              double3 ds = (random_position - atom.position + translation) / numberOfReplicas;

              ds.x -= std::rint(ds.x);
              ds.y -= std::rint(ds.y);
              ds.z -= std::rint(ds.z);

              double3 dr = simulation_cell * ds;

              double rr = double3::dot(dr, dr);
              if (rr < cutoff * cutoff)
              {
                Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
                  forceField, 0, 0, 1.0, 1.0, rr, probe_type.value(), typeB);

                energy += energyFactor.energy;
              }
            }
          }
        }
      }

      sumBoltzmannWeight += std::exp(-(1.0 / (Units::KB * 298.0) * energy));
      count++;
    }
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".energy.vf.cpu.txt");
  std::print(myfile, "# Void-fraction using energy-based method\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
  std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Framework volume: {} [Å³]\n", framework.simulationBox.volume);
  std::print(myfile, "# Framework mass: {} [g/mol]\n", framework.unitCellMass);
  std::print(myfile, "# Framework density: {} [kg/m³]\n", 1e-3 * framework.unitCellMass /
      (framework.simulationBox.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant));
  std::print(myfile, "# Probe atom: {} strength-value/kʙ: {} [K] size-value: {} [Å]\n", probePseudoAtom, 
      forceField[probe_type.value()].strengthParameter() * Units::EnergyToKelvin, 
      forceField[probe_type.value()].sizeParameter());
  std::print(myfile, "# Cutoff: {} Å\n", cutoff);
  std::print(myfile, "# Number of unit cells: {}x{}x{}\n", numberOfReplicas.x, numberOfReplicas.y, numberOfReplicas.z);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << sumBoltzmannWeight / static_cast<double>(count) << std::endl;
  myfile.close();
}
