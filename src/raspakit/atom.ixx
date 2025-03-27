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
  double3 position;            ///< The position of the atom in 3D space.
  double3 velocity{};          ///< The velocity of the atom.
  double3 gradient{};          ///< The gradient acting on the atom.
  double charge;               ///< The electric charge of the atom.
  double scalingVDW{1.0};      ///< Scaling factor for van der Waals interactions.
  double scalingCoulomb{1.0};  ///< Scaling factor for Coulomb interactions.
  uint32_t moleculeId{0};      ///< Identifier for the molecule this atom belongs to.
  uint16_t type{0};            ///< Type identifier of the atom.
  uint8_t componentId{0};      ///< Component identifier within the system.
  uint8_t groupId{0};          ///< Group identifier, defaults to false.

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
    groupId = uint8_t{0};
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

export void to_json(nlohmann::json &j, const Atom &a)
{
  j = nlohmann::json{{"position", a.position},       {"velocity", a.velocity},
                     {"gradient", a.gradient},       {"charge", a.charge},
                     {"scalingVDW", a.scalingVDW},   {"scalingCoulomb", a.scalingCoulomb},
                     {"moleculeId", a.moleculeId},   {"type", a.type},
                     {"componentId", a.componentId}, {"groupId", a.groupId}};
}

export void from_json(const nlohmann::json &j, Atom &a)
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
