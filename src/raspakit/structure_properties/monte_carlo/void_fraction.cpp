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

  std::optional<size_t> probeType = forceField.findPseudoAtom("He");

  if(!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  Component helium = Component(0, forceField, "helium", 5.2, 228000.0, -0.39,
                               {Atom({0, 0, 0}, 0.0, 1.0, 0, static_cast<uint16_t>(probeType.value()), 0, 0)}, 5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {framework}, {helium}, {0}, 5);

  double average_Rosenbluth_weight{};
  double count{};
  for(size_t i = 0; i < number_of_iterations; ++i)
  {
    std::pair<double, double> result = MC_Moves::WidomMove(random, system, 0);
    average_Rosenbluth_weight += result.first;
    count += 1.0;

  }
  std::ofstream myfile;
  myfile.open(framework.name + ".mc.vf.txt");
  myfile << average_Rosenbluth_weight / count << std::endl;
  myfile.close();
}
