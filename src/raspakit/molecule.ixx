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

// Note: C++17 and higher: std::vector<T> is automatically properly aligned based on type T
export struct Molecule
{
  double3 centerOfMassPosition;
  double3 velocity;
  double3 gradient;
  simd_quatd orientation;    
  simd_quatd orientationMomentum;
  simd_quatd orientationGradient;

  Molecule() noexcept = default;
  Molecule(double3 centerOfMassPosition, simd_quatd orientation) :
    centerOfMassPosition(centerOfMassPosition), velocity(0.0, 0.0, 0.0), gradient(0.0, 0.0, 0.0),
    orientation(orientation),orientationMomentum(0.0, 0.0, 0.0, 0.0), orientationGradient(0.0, 0.0, 0.0, 0.0)
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

// should be 6 times double4 = 6x(8x4) = 6x32 = 192 bytes
static_assert(sizeof(Molecule) == 192, "struct Molecule size is not 192");
