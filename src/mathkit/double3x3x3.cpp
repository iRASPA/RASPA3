module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <map>
#include <numbers>
#include <ostream>
#include <vector>
#endif

module double3x3x3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3x3;
import double3x3;
import simd_quatd;
import double3;
import archive;
