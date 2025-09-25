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
import std;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;
import json;

// Note: C++17 and higher: std::vector<T> is automatically properly aligned based on type T
export struct Molecule
{
  double3 centerOfMassPosition;    ///< The center of mass position of the molecule in 3D space.
  double3 velocity;                ///< The velocity of the molecule.
  double3 gradient;                ///< The gradient (force) acting on the molecule.
  simd_quatd orientation;          ///< The orientation of the molecule represented as a quaternion.
  simd_quatd orientationMomentum;  ///< The angular momentum of the molecule's orientation.
  simd_quatd orientationGradient;  ///< The gradient (torque) acting on the molecule's orientation.
  double mass;                     ///< Molecular mass
  double invMass;                  ///< 1/mass to save on computing (and set correct size)
  double padding1;
  double padding2;
  std::size_t atomIndex;      ///< Pointing to the index in the list of atoms
  std::size_t numberOfAtoms;  ///< Number of subatoms in this molecule
  std::size_t componentId;    ///< Pointing to the index in the components list
  std::size_t padding3;

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

  Molecule(double3 centerOfMassPosition, simd_quatd orientation, double mass, std::size_t componentId,
           std::size_t numberOfAtoms)
      : centerOfMassPosition(centerOfMassPosition),
        velocity(0.0, 0.0, 0.0),
        gradient(0.0, 0.0, 0.0),
        orientation(orientation),
        orientationMomentum(0.0, 0.0, 0.0, 0.0),
        orientationGradient(0.0, 0.0, 0.0, 0.0),
        mass(mass),
        invMass(1.0 / mass),
        numberOfAtoms(numberOfAtoms),
        componentId(componentId) {};

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

// should be 8 times double4 = 8x(8x4) = 8x32 = 256 bytes
static_assert(sizeof(Molecule) == 256, "struct Molecule size is not 256");

export void to_json(nlohmann::json &j, const Molecule &a)
{
  j = nlohmann::json{{"centerOfMassPosition", a.centerOfMassPosition},
                     {"velocity", a.velocity},
                     {"gradient", a.gradient},
                     {"orientation", a.orientation},
                     {"orientationMomentum", a.orientationMomentum},
                     {"orientationGradient", a.orientationGradient},
                     {"mass", a.mass},
                     {"atomIndex", a.atomIndex},
                     {"numberOfAtoms", a.numberOfAtoms},
                     {"componentId", a.componentId}};
}

export void from_json(const nlohmann::json &j, Molecule &a)
{
  j.at("centerOfMassPosition").get_to(a.centerOfMassPosition);
  j.at("velocity").get_to(a.velocity);
  j.at("gradient").get_to(a.gradient);
  j.at("orientation").get_to(a.orientation);
  j.at("orientationMomentum").get_to(a.orientationMomentum);
  j.at("orientationGradient").get_to(a.orientationGradient);
  j.at("mass").get_to(a.mass);
  j.at("atomIndex").get_to(a.atomIndex);
  j.at("numberOfAtoms").get_to(a.numberOfAtoms);
  j.at("componentId").get_to(a.componentId);
}
