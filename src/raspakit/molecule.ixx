module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
#endif

export module molecule;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <cstddef>;
import <istream>;
import <ostream>;
import <fstream>;
import <sstream>;
import <type_traits>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

// C++17 and higher: std::vector<T> is automatically properly aligned based on type T
export struct Molecule
{
  double3 centerOfMassPosition;
  simd_quatd orientation;    

  Molecule(double3 centerOfMassPosition, simd_quatd orientation) :
    centerOfMassPosition(centerOfMassPosition), orientation(orientation)
  {
  };

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Molecule &molecule);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Molecule &molecule);

  inline std::string repr() const
  {
    std::ostringstream stream;
  
    return stream.str();
  }
};

// should be 2 times double4 = 2x(8x4) = 2x32 = 64 bytes
static_assert(sizeof(Molecule) == 64, "struct Molecule size is not 64");
