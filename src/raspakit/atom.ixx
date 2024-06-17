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

export module atom;

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
import stringutils;
import json;

import scaling;

// Note: C++17 and higher: std::vector<T> is automatically properly aligned based on type T
export struct Atom
{
  double3 position;
  double3 velocity{};
  double3 gradient{};
  double charge;
  double scalingVDW{1.0};
  double scalingCoulomb{1.0};
  uint32_t moleculeId{0};
  uint16_t type{0};
  uint8_t componentId{0};
  uint8_t groupId{0};  // defaults to false

  Atom() noexcept = default;
  Atom(double3 position, double charge, double scalingVDW, double scalingCoulomb, uint32_t moleculeId, uint16_t type,
       uint8_t componentId, uint8_t groupId)
      : position(position),
        charge(charge),
        scalingVDW(scalingVDW),
        scalingCoulomb(scalingCoulomb),
        moleculeId(moleculeId),
        type(type),
        componentId(componentId),
        groupId(groupId)
  {
  }
  Atom(double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       uint8_t groupId)
      : position(position),
        charge(charge),
        moleculeId(moleculeId),
        type(type),
        componentId(componentId),
        groupId(groupId)
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
    groupId = uint8_t{0};
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Atom &atom);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Atom &atom);

  friend void to_json(nlohmann::json &, const Atom &);
  friend void from_json(const nlohmann::json &, Atom &);

  inline std::string repr() const
  {
    std::ostringstream stream;

    std::print(stream, "({}, {}, {}, [{}, {}, {}, {}])\n", position.x, position.y, position.z, moleculeId, type,
               componentId, groupId);

    return stream.str();
  }
};

// should be 4 times double4 = 4x(8x4) = 4x32 = 128 bytes
static_assert(sizeof(Atom) == 128, "struct Atom size is not 128");

void to_json(nlohmann::json &j, const Atom &a)
{
  j = nlohmann::json{{"position", a.position},       {"velocity", a.velocity},
                     {"gradient", a.gradient},       {"charge", a.charge},
                     {"scalingVDW", a.scalingVDW},   {"scalingCoulomb", a.scalingCoulomb},
                     {"moleculeId", a.moleculeId},   {"type", a.type},
                     {"componentId", a.componentId}, {"groupId", a.groupId}};
}

void from_json(const nlohmann::json &j, Atom &a)
{
  j.at("position").get_to(a.position);
  j.at("velocity").get_to(a.velocity);
  j.at("gradient").get_to(a.gradient);
  j.at("charge").get_to(a.charge);
  j.at("scalingVDW").get_to(a.scalingVDW);
  j.at("scalingCoulomb").get_to(a.scalingCoulomb);
  j.at("moleculeId").get_to(a.moleculeId);
  j.at("type").get_to(a.type);
  j.at("componentId").get_to(a.componentId);
  j.at("groupId").get_to(a.groupId);
}
