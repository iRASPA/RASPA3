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

/**
 * \brief Represents a molecule in the simulation system.
 *
 * The Molecule struct encapsulates the properties and behaviors of a molecule
 * within the simulation. It includes the center of mass position, velocity,
 * gradient, orientation (as a quaternion), orientation momentum, and orientation
 * gradient. The struct provides constructors for initializing molecules and
 * methods for serialization and deserialization. Molecule class type should
 * have exactly a size of 192 bytes (6 times double4).
 */
export struct Molecule
{
  double3 centerOfMassPosition;    ///< The center of mass position of the molecule in 3D space.
  double3 velocity;                ///< The velocity of the molecule.
  double3 gradient;                ///< The gradient (force) acting on the molecule.
  simd_quatd orientation;          ///< The orientation of the molecule represented as a quaternion.
  simd_quatd orientationMomentum;  ///< The angular momentum of the molecule's orientation.
  simd_quatd orientationGradient;  ///< The gradient (torque) acting on the molecule's orientation.

  /**
   * \brief Default constructor for the Molecule struct.
   *
   * Initializes a Molecule object with default values.
   */
  Molecule() noexcept = default;

  /**
   * \brief Constructs a Molecule with specified center of mass position and orientation.
   *
   * Initializes a Molecule with the provided center of mass position and orientation.
   * Other members are initialized to zero.
   *
   * \param centerOfMassPosition The initial center of mass position of the molecule.
   * \param orientation The initial orientation of the molecule represented as a quaternion.
   */
  Molecule(double3 centerOfMassPosition, simd_quatd orientation)
      : centerOfMassPosition(centerOfMassPosition),
        velocity(0.0, 0.0, 0.0),
        gradient(0.0, 0.0, 0.0),
        orientation(orientation),
        orientationMomentum(0.0, 0.0, 0.0, 0.0),
        orientationGradient(0.0, 0.0, 0.0, 0.0) {};

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Molecule &molecule);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Molecule &molecule);
  friend void to_json(nlohmann::json &j, const Molecule &a);
  friend void from_json(const nlohmann::json &j, Molecule &a);

  /**
   * \brief Returns a string representation of the Molecule.
   *
   * Generates a string that includes the properties of the molecule.
   *
   * \return A string representing the Molecule.
   */
  inline std::string repr() const
  {
    std::ostringstream stream;

    return stream.str();
  }
};

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
