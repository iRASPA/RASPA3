module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <tuple>
#include <vector>
#endif

export module cbmc_util;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import double3x3;
import double3;
import simd_quatd;
import randomnumbers;

export namespace CBMC
{
std::vector<Atom> rotateRandomlyAround(RandomNumber &random, std::span<Atom> atoms, std::size_t startingBead);
std::vector<Atom> rotateRandomlyAround(simd_quatd &q, std::span<Atom> atoms, std::size_t startingBead);

// LogBoltzmannFactors are (-Beta U)
std::size_t selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors);
}  // namespace CBMC
