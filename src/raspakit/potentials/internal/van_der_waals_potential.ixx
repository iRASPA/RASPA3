module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <type_traits>
#include <vector>
#endif

export module van_der_waals_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for vanDerWaals potentials.
 *
 * Defines the maximum number of parameters that can be associated with a vanDerWaals potential.
 */
export const std::size_t maximumNumberOfVanDerWaalsParameters{4};

/**
 * \brief Enumeration of different vanDerWaals types.
 *
 * Specifies the type of vanDerWaals potential to be used in simulations.
 */
export enum class VanDerWaalsType : std::size_t { LennardJones = 0 };

/**
 * \brief Represents a vanDerWaals potential between two particles.
 *
 * The vanDerWaalsPotential struct encapsulates the type of vanDerWaals and associated parameters between two particles.
 * It includes versioning for serialization, vanDerWaals type, identifiers of vanDerWaalsed particles, and vanDerWaals
 * parameters.
 */
export struct VanDerWaalsPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 2> identifiers;  ///< Identifiers of the two particles forming the vanDerWaals.
  VanDerWaalsType type;                    ///< The type of vanDerWaals potential.
  double scaling;
  std::array<double, maximumNumberOfVanDerWaalsParameters>
      parameters;  ///< Parameters associated with the vanDerWaals potential.

  /**
   * \brief Default constructor for vanDerWaalsPotential.
   *
   * Initializes a vanDerWaalsPotential object with Undefined vanDerWaals type and zeroed vanDerWaals IDs.
   */
  VanDerWaalsPotential() : identifiers({0, 0}), type(VanDerWaalsType::LennardJones), scaling(1.0) {}

  VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, VanDerWaalsType type,
                       std::vector<double> vector_parameters, double scaling);

  /**
   * \brief Constructs a vanDerWaalsPotential with specified type and vanDerWaals IDs.
   *
   * \param type The type of vanDerWaals potential.
   * \param identifiers A pair of particle identifiers forming the vanDerWaals.
   */
  VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, const VanDerWaalsType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(VanDerWaalsPotential const &) const = default;

  /**
   * \brief Generates a string representation of the vanDerWaals potential.
   *
   * Provides a formatted string containing vanDerWaals type, particle IDs, and parameters.
   *
   * \return A string describing the vanDerWaals potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each vanDerWaals type.
   *
   * A static vector indicating the number of parameters needed for each vanDerWaals type.
   */
  static inline std::array<std::size_t, 2> numberOfVanDerWaalsParameters{2};

  /**
   * \brief Mapping of vanDerWaals type strings to type enums.
   *
   * A static map that associates vanDerWaals type names with their corresponding type enumeration values.
   */
  static inline std::map<std::string, VanDerWaalsType, caseInsensitiveComparator> definitionForString{
      {"LENNARD_JONES", VanDerWaalsType::LennardJones}};

  double calculateEnergy(const double3 &posA, const double3 &posB) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VanDerWaalsPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VanDerWaalsPotential &b);
};
