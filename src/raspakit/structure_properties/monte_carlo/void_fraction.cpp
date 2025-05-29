module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#include <optional>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <tuple>
#endif

module mc_void_fraction;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <print>;
import <string>;
import <vector>;
import <optional>;
import <limits>;
import <algorithm>;
import <iostream>;
import <fstream>;
import <tuple>;
#endif

import double3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import component;
import system;
import mc_moves_widom;



void MC_VoidFraction::run(const ForceField &forceField, const Framework &framework, size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  std::optional<size_t> probeType = forceField.findPseudoAtom("He");

  if(!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  Component helium = Component(0, forceField, "helium", 5.2, 228000.0, -0.39,
                               {Atom({0, 0, 0}, 0.0, 1.0, 0, static_cast<uint16_t>(probeType.value()), 0, 0)}, 5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {framework}, {helium}, {0}, 5);

  double no_overlap{};
  double count{};
  for(size_t i = 0; i < 20 * number_of_iterations; ++i)
  {
    double3 s = double3(random.uniform(), random.uniform(), random.uniform());
    double3 pos = framework.simulationBox.cell * s;

    if(!framework.computeVanDerWaalsRadiusOverlap(forceField, pos))
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
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << no_overlap / count << std::endl;
  myfile.close();
}
