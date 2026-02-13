module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <string>
#include <tuple>
#include <vector>
#endif

module mc_void_fraction;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import double3x3;
import randomnumbers;
import skspacegroupdatabase;
import atom;
import framework;
import forcefield;
import component;
import system;
import mc_moves_widom;
import units;

void MC_VoidFraction::run(const ForceField &forceField, const Framework &framework, double wellDepthFactor, std::string probePseudoAtom, std::optional<std::size_t> numberOfIterations)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_VoidFraction: Unknown probe-atom type\n"));
  }

  std::size_t number_of_iterations = numberOfIterations.value_or(100000);

  double no_overlap{};
  double count{};
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double3 s = double3(random.uniform(), random.uniform(), random.uniform());
    double3 cartesian_position = framework.simulationBox.cell * s;

    if (!framework.computeOverlap(forceField, cartesian_position, wellDepthFactor, probeType.value(), -1))
    //if (!framework.computeVanDerWaalsRadiusOverlap(forceField, cartesian_position))
    {
      no_overlap += 1.0;
    }
    count += 1.0;
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.vf.cpu.txt");
  std::print(myfile, "# Void-fraction using Mont Carlo-based method\n");
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
  std::print(myfile, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", probePseudoAtom, wellDepthFactor, forceField[probeType.value()].sizeParameter());
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << no_overlap / count << std::endl;
  myfile.close();
}
