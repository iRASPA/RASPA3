module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <format>
#include <iostream>
#include <tuple>
#include <vector>
// #if defined(__has_include) && __has_include(<stacktrace>)
// #include <stacktrace>
// #endif
#include <print>
#endif

module cbmc_util;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import double3x3;
import double3;
import simd_quatd;
import randomnumbers;
import stringutils;

std::vector<Atom> CBMC::rotateRandomlyAround(RandomNumber &random, std::vector<Atom> atoms, std::size_t startingBead)
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> randomlyRotatedAtoms{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    Atom b = atoms[i];
    b.position = atoms[startingBead].position + randomRotationMatrix * (b.position - atoms[startingBead].position);
    randomlyRotatedAtoms.push_back(b);
  }
  return randomlyRotatedAtoms;
}

std::vector<Atom> CBMC::rotateRandomlyAround(simd_quatd &q, std::vector<Atom> atoms, std::size_t startingBead)
{
  double3x3 randomRotationMatrix = double3x3::buildRotationMatrixInverse(q);
  std::vector<Atom> randomlyRotatedAtoms{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    Atom b = atoms[i];
    b.position = atoms[startingBead].position + randomRotationMatrix * (b.position - atoms[startingBead].position);
    randomlyRotatedAtoms.push_back(b);
  }
  return randomlyRotatedAtoms;
}

// LogBoltzmannFactors are (-Beta U)
std::size_t CBMC::selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors)
{
  std::vector<double> ShiftedBoltzmannFactors(LogBoltzmannFactors.size());

  // Energies are always bounded from below [-U_max, infinity>
  // Find the lowest energy value, i.e. the largest value of (-Beta U)
  std::vector<double>::iterator match = std::max_element(LogBoltzmannFactors.begin(), LogBoltzmannFactors.end());
  
  if (match == LogBoltzmannFactors.end())
  {
    throw std::runtime_error("[cbmc-utils]: no maximum value found\n");
  }
  double largest_value = *match;

  // Standard trick: shift the Boltzmann factors down to avoid numerical problems
  // The largest value of 'ShiftedBoltzmannFactors' will be 1 (which corresponds to the lowest energy).
  double SumShiftedBoltzmannFactors = 0.0;
  for (std::size_t i = 0; i < LogBoltzmannFactors.size(); ++i)
  {
    ShiftedBoltzmannFactors[i] = std::exp(LogBoltzmannFactors[i] - largest_value);
    SumShiftedBoltzmannFactors += ShiftedBoltzmannFactors[i];
  }

  // select the Boltzmann factor
  std::size_t selected = 0;
  double cumw = ShiftedBoltzmannFactors[0];
  double ws = random.uniform() * SumShiftedBoltzmannFactors;
  while (cumw < ws) cumw += ShiftedBoltzmannFactors[++selected];

  return selected;
}
