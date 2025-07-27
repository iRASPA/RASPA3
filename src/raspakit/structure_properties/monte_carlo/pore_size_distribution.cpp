module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <random>
#include <string>
#include <vector>
#endif

module mc_pore_size_distribution;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import framework;
import forcefield;
import atom;

void MC_PoreSizeDistribution::run(const ForceField &forceField, const Framework &framework, double well_depth_factor,
                                  std::size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  double delta_r = 10.0 / static_cast<double>(data.size());

  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double3 sA = double3(random.uniform(), random.uniform(), random.uniform());
    double3 posA = framework.simulationBox.cell * sA;

    if (!framework.computeLargestNonOverlappingFreeRadius(forceField, posA, well_depth_factor).has_value()) continue;

    double largest_radius = std::numeric_limits<double>::lowest();
    for (std::size_t j = 0; j < number_of_iterations; ++j)
    {
      double3 sB = double3(random.uniform(), random.uniform(), random.uniform());
      double3 posB = framework.simulationBox.cell * sB;

      std::optional<double> radius =
          framework.computeLargestNonOverlappingFreeRadius(forceField, posB, well_depth_factor);

      if (!radius.has_value()) continue;

      double3 dr = posA - posB;
      dr = framework.simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      if (r > radius.value()) continue;

      largest_radius = std::max(largest_radius, radius.value());
    }

    std::size_t index = static_cast<std::size_t>(largest_radius / delta_r);
    if (index < data.size())
    {
      for (std::size_t k = 0; k <= index; ++k)
      {
        data[k] += 1.0;
      }
    }
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.psd.cpu.txt");
  std::print(myfile, "# Pore-size distribution using Mont Carlo-based method\n");
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  double normalization = 1.0 / static_cast<double>(10000);
  myfile << "# column 1: diameter d [A]\n";
  myfile << "# column 2: cumulative pore volume\n";
  myfile << "# column 3: PSD\n";
  myfile << "# value at d=0 is related to the void-fraction\n";

  for (std::size_t i = 0; i < data.size(); ++i)
  {
    if (data[i] > 0)
    {
      double derivative{};
      if ((i > 2) && (i < data.size() - 2))
      {
        derivative = (-data[i + 2] + 8.0 * data[i + 1] - 8.0 * data[i - 1] + data[i - 2]) / (12.0 * delta_r);
      }

      myfile << 2.0 * delta_r * (static_cast<double>(i) + 0.5) << " " << normalization * data[i] << " "
             << normalization * -derivative << "\n";
    }
  }
  myfile.close();
}
