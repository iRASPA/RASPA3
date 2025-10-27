module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <istream>
#include <ostream>
#include <print>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#endif

export module atom;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import double3;
import stringutils;
import json;
import scaling;

/**
 * \brief Represents an atom in the simulation system.
 *
 * The Atom struct encapsulates the properties and behaviors of an individual atom
 * within the simulation. It includes positional data, velocity, gradient, charge,
 * scaling factors for van der Waals and Coulomb interactions, and identifiers for
 * molecule, type, component, and group associations. The struct provides constructors
 * for initializing atoms and methods to adjust scaling parameters dynamically.
 * Atom class type should have exactly a size of 128 bytes (4 times double4).
 */
export struct Atom
{
  double3 position;               ///< The position of the atom in 3D space.
  double3 velocity{};             ///< The velocity of the atom.
  double3 gradient{};             ///< The gradient acting on the atom.
  double charge;                  ///< The electric charge of the atom.
  double scalingVDW{1.0};         ///< Scaling factor for van der Waals interactions.
  double scalingCoulomb{1.0};     ///< Scaling factor for Coulomb interactions.
  std::uint32_t moleculeId{0};    ///< Identifier for the molecule this atom belongs to.
  std::uint16_t type{0};          ///< Type identifier of the atom.
  std::uint8_t componentId{0};    ///< Component identifier within the system.
  std::uint8_t groupId : 4;       ///< Group identifier, defaults to false.
  std::uint8_t isFractional : 4;  ///< Fractional or not, defaults to false.

  /**
   * \brief Default constructor for the Atom struct.
   *
   * Initializes an Atom object with default values.
   */
  Atom() noexcept = default;

  /**
   * \brief Constructs an Atom with specified parameters.
   *
   * Initializes an Atom with the provided position, charge, scaling factors,
   * molecule ID, type, component ID, and group ID.
   *
   * \param position The initial position of the atom.
   * \param charge The electric charge of the atom.
   * \param scalingVDW Scaling factor for van der Waals interactions.
   * \param scalingCoulomb Scaling factor for Coulomb interactions.
   * \param moleculeId Identifier for the molecule.
   * \param type Type identifier of the atom.
   * \param componentId Component identifier within the system.
   * \param groupId Group identifier, defaults to false.
   */
  Atom(double3 position, double charge, double scalingVDW, double scalingCoulomb, std::uint32_t moleculeId,
       std::uint16_t type, std::uint8_t componentId, std::uint8_t groupId, std::uint8_t isFractional)
      : position(position),
        charge(charge),
        scalingVDW(scalingVDW),
        scalingCoulomb(scalingCoulomb),
        moleculeId(moleculeId),
        type(type),
        componentId(componentId),
        groupId(groupId),
        isFractional(isFractional)
  {
  }

  /**
   * \brief Constructs an Atom with lambda-dependent scaling.
   *
   * Initializes an Atom with position, charge, and scaling factors determined
   * by the provided lambda value. Also sets molecule ID, type, component ID,
   * and group ID.
   *
   * \param position The initial position of the atom.
   * \param charge The electric charge of the atom.
   * \param lambda The scaling parameter for interactions.
   * \param moleculeId Identifier for the molecule.
   * \param type Type identifier of the atom.
   * \param componentId Component identifier within the system.
   * \param groupId Group identifier, defaults to false.
   */
  Atom(double3 position, double charge, double lambda, std::uint32_t moleculeId, std::uint16_t type,
       std::uint8_t componentId, std::uint8_t groupId, std::uint8_t isFractional)
      : position(position),
        charge(charge),
        moleculeId(moleculeId),
        type(type),
        componentId(componentId),
        groupId(groupId),
        isFractional(isFractional)
  {
    scalingVDW = Scaling::scalingVDW(lambda);
    scalingCoulomb = Scaling::scalingCoulomb(lambda);
  };

  /**
   * \brief Sets the scaling factors based on lambda.
   *
   * Adjusts the scaling factors for van der Waals and Coulomb interactions
   * linearly based on the provided lambda value. Scaling for van der Waals
   * is active from 0 to 0.5, and Coulomb scaling is active from 0.5 to 1.0.
   *
   * \param lambda The scaling parameter.
   */
  void setScaling(double lambda)
  {
    scalingVDW = Scaling::scalingVDW(lambda);
    scalingCoulomb = Scaling::scalingCoulomb(lambda);
  }

  /**
   * \brief Fully activates scaling factors.
   *
   * Sets both van der Waals and Coulomb scaling factors to 1.0, fully
   * activating the interactions.
   */
  void setScalingFullyOn()
  {
    scalingVDW = 1.0;
    scalingCoulomb = 1.0;
  }

  /**
   * \brief Fully deactivates scaling factors.
   *
   * Sets both van der Waals and Coulomb scaling factors to 0.0, fully
   * deactivating the interactions.
   */
  void setScalingFullyOff()
  {
    scalingVDW = 0.0;
    scalingCoulomb = 0.0;
  }

  /**
   * \brief Sets scaling factors to integer values.
   *
   * Sets both van der Waals and Coulomb scaling factors to 1.0 and resets
   * the group ID to 0.
   */
  void setScalingToInteger()
  {
    scalingVDW = 1.0;
    scalingCoulomb = 1.0;
    groupId = false;
    isFractional = false;
  }

  // scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Atom &atom);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Atom &atom);
  friend void to_json(nlohmann::json &, const Atom &);
  friend void from_json(const nlohmann::json &, Atom &);

  /**
   * \brief Returns a string representation of the Atom.
   *
   * Generates a string that includes the position coordinates and identifiers
   * associated with the atom.
   *
   * \return A string representing the Atom.
   */
  inline std::string repr() const
  {
    std::ostringstream stream;

    std::print(stream, "({}, {}, {}, [{}, {}, {}, {}])\n", position.x, position.y, position.z, moleculeId, type,
               componentId, groupId);

    return stream.str();
  }
};

// does not compile on llvm-18
/*
export template <>
struct std::formatter<Atom>: std::formatter<string_view>
{
  auto format(const Atom& atom, std::format_context& ctx) const
  {
    std::string temp{};
    std::format_to(std::back_inserter(temp), "(position: {}, charge: {}, scalings: {} {}, molecule-id: {}, type: {},
component-id: {}, group: {}, fractional: {})", atom.position, atom.charge, atom. scalingVDW, atom.scalingCoulomb,
atom.moleculeId, atom.type, atom.componentId, atom.groupId, atom.isFractional); return
std::formatter<string_view>::format(temp, ctx);
  }
};
*/

export void to_json(nlohmann::json &j, const Atom &a)
{
  j = nlohmann::json{{"position", a.position},
                     {"velocity", a.velocity},
                     {"gradient", a.gradient},
                     {"charge", a.charge},
                     {"scalingVDW", a.scalingVDW},
                     {"scalingCoulomb", a.scalingCoulomb},
                     {"moleculeId", a.moleculeId},
                     {"type", a.type},
                     {"componentId", a.componentId},
                     {"groupId", static_cast<std::uint8_t>(a.groupId)},
                     {"isFractional", static_cast<std::uint8_t>(a.isFractional)}};
}

export void from_json(const nlohmann::json &j, Atom &a)
{
  std::uint8_t groupIdBitField;
  std::uint8_t isFractionalBitField;
  j.at("position").get_to(a.position);
  j.at("velocity").get_to(a.velocity);
  j.at("gradient").get_to(a.gradient);
  j.at("charge").get_to(a.charge);
  j.at("scalingVDW").get_to(a.scalingVDW);
  j.at("scalingCoulomb").get_to(a.scalingCoulomb);
  j.at("moleculeId").get_to(a.moleculeId);
  j.at("type").get_to(a.type);
  j.at("componentId").get_to(a.componentId);
  j.at("groupId").get_to(groupIdBitField);
  a.groupId = groupIdBitField;
  j.at("isFractional").get_to(isFractionalBitField);
  a.isFractional = isFractionalBitField;
}
