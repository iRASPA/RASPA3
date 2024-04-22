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

export module atom;

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
import stringutils;

import scaling;


// Note: C++17 and higher: std::vector<T> is automatically properly aligned based on type T
export struct Atom
{
  double3 position;
  double3 velocity{};    
  double3 gradient{};    
  double charge;
  double scalingVDW{ 1.0 };
  double scalingCoulomb{ 1.0 };
  uint32_t moleculeId{ 0 };
  uint16_t type{ 0 };
  uint8_t componentId{ 0 };
  uint8_t groupId{ 0 };   // defaults to false

  Atom() noexcept = default;
  Atom(double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId) :
    position(position), charge(charge), moleculeId(moleculeId), type(type), componentId(componentId), groupId(groupId)
  {
    scalingVDW = Scaling::scalingVDW(lambda);
    scalingCoulomb = Scaling::scalingCoulomb(lambda);
  };

  // scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
  void setScaling(double lambda)
  {
    scalingVDW = Scaling::scalingVDW(lambda);
    scalingCoulomb = Scaling::scalingCoulomb(lambda);
  }

  void setScalingFullyOn()
  {
    scalingVDW = 1.0;
    scalingCoulomb = 1.0;
  }

  void setScalingFullyOff()
  {
    scalingVDW = 0.0;
    scalingCoulomb = 0.0;
  }

  void setScalingToInteger()
  {
    scalingVDW = 1.0;
    scalingCoulomb = 1.0;
    groupId = uint8_t{ 0 };
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Atom &atom);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Atom &atom);

  inline std::string repr() const
  {
    std::ostringstream stream;
  
    std::print(stream, "({}, {}, {}, [{}, {}, {}, {}])\n", position.x, position.y, position.z,
        moleculeId, type, componentId, groupId);
  
    return stream.str();
  }
};

// should be 4 times double4 = 4x(8x4) = 4x32 = 128 bytes
static_assert(sizeof(Atom) == 128, "struct Atom size is not 128");
