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
#include <vector>
#pragma push_macro("__SSE3__")
#undef __SSE3__
#include <random>
#pragma pop_macro("__SSE3__")
#endif

module mc_pore_size_distribution;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import double3x3;
import randomnumbers;
import framework;
import forcefield;
import atom;

void MC_PoreSizeDistribution::run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
                                  std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps,
                                  std::optional<double> maximumRange)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin, time_end;

  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);

  time_begin = std::chrono::system_clock::now();

  double delta_r = maximumRange.value_or(10.0) / static_cast<double>(numberOfBins);

  std::size_t number_of_iterations = numberOfIterations.value_or(10000);

  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double3 sA = double3(random.uniform(), random.uniform(), random.uniform());
    double3 posA = framework.simulationBox.cell * sA;

    if (!framework.computeLargestNonOverlappingFreeRadius(forceField, posA, wellDepthFactor).has_value()) continue;

    double largest_radius = std::numeric_limits<double>::lowest();
    for (std::size_t j = 0; j < number_of_inner_steps; ++j)
    {
      double3 sB = double3(random.uniform(), random.uniform(), random.uniform());
      double3 posB = framework.simulationBox.cell * sB;

      std::optional<double> radius =
          framework.computeLargestNonOverlappingFreeRadius(forceField, posB, wellDepthFactor);

      if (!radius.has_value()) continue;

      double3 dr = posA - posB;
      dr = framework.simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);

      if (rr > radius.value() * radius.value()) continue;

      largest_radius = std::max(largest_radius, radius.value());
    }

    if(largest_radius >= 0.0)
    {
      std::size_t index = static_cast<std::size_t>(largest_radius / delta_r);
      if(index >= 0 && index < numberOfBins)
      {
        histogram[index]++;
        for(std::size_t k = 0; k <= index; ++k)
        {
           histogram_cummulative[k]++;
        }
      }
    }
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.psd.cpu.txt");
  std::print(myfile, "# Pore-size distribution using Mont Carlo-based method\n");
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << "# column 1: diameter d [A]\n";
  myfile << "# column 2: PSD\n";
  myfile << "# column 3: cumulative pore volume\n";
  myfile << "# value at d=0 is related to the void-fraction\n";

  double normalization = 1.0 / static_cast<double>(number_of_iterations);
  for(std::size_t index = 0; index < numberOfBins; ++index)
  {
    std::print(myfile, "{} {} {}\n", 2.0 * delta_r * (static_cast<double>(index) + 0.5),
               histogram[index] * normalization,
               histogram_cummulative[index] * normalization);
  }

  myfile.close();
}
