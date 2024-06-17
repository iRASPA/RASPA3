module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <ostream>
#include <print>
#include <sstream>
#include <type_traits>
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
import <print>;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;
import json;

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
  Molecule(double3 centerOfMassPosition, simd_quatd orientation)
      : centerOfMassPosition(centerOfMassPosition),
        velocity(0.0, 0.0, 0.0),
        gradient(0.0, 0.0, 0.0),
        orientation(orientation),
        orientationMomentum(0.0, 0.0, 0.0, 0.0),
        orientationGradient(0.0, 0.0, 0.0, 0.0) {};

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Molecule &molecule);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Molecule &molecule);

  friend void to_json(nlohmann::json &, const Molecule &);
  friend void from_json(const nlohmann::json &, Molecule &);

  inline std::string repr() const
  {
    std::ostringstream stream;

    return stream.str();
  }
};

// should be 6 times double4 = 6x(8x4) = 6x32 = 192 bytes
static_assert(sizeof(Molecule) == 192, "struct Molecule size is not 192");

void to_json(nlohmann::json &j, const Molecule &a)
{
  j = nlohmann::json{{"centerOfMassPosition", a.centerOfMassPosition},
                     {"velocity", a.velocity},
                     {"gradient", a.gradient},
                     {"orientation", a.orientation},
                     {"orientationMomentum", a.orientationMomentum},
                     {"orientationGradient", a.orientationGradient}};
}

void from_json(const nlohmann::json &j, Molecule &a)
{
  j.at("centerOfMassPosition").get_to(a.centerOfMassPosition);
  j.at("velocity").get_to(a.velocity);
  j.at("gradient").get_to(a.gradient);
  j.at("orientation").get_to(a.orientation);
  j.at("orientationMomentum").get_to(a.orientationMomentum);
  j.at("orientationGradient").get_to(a.orientationGradient);
}
