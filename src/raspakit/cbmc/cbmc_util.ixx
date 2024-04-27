module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#endif

export module cbmc_util;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <cmath>;
import <algorithm>;
#endif

import atom;
import double3x3;
import double3;
import simd_quatd;
import randomnumbers;


export namespace CBMC
{
  std::vector<Atom> rotateRandomlyAround(RandomNumber &random, std::vector<Atom> atoms, size_t startingBead);
  std::vector<Atom> rotateRandomlyAround(simd_quatd &q, std::vector<Atom> atoms, size_t startingBead);

  // LogBoltzmannFactors are (-Beta U)
  size_t selectTrialPosition(RandomNumber &random, std::vector <double> LogBoltzmannFactors);
}
