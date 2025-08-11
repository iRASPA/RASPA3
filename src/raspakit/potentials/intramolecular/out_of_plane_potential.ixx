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
#include <tuple>
#include <type_traits>
#include <vector>
#endif

export module out_of_plane_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for out_of_plane_bend potentials.
 *
 * Defines the maximum number of parameters that can be associated with a out_of_plane_bend potential.
 */
export const std::size_t maximumNumberOfOutOfPlaneBendParameters{4};

/**
 * \brief Enumeration of different out_of_plane_bend types.
 *
 * Specifies the type of out_of_plane_bend potential to be used in simulations.
 */
export enum class OutOfPlaneBendType : std::size_t { Harmonic = 0 };

/**
 * \brief Represents a out_of_plane_bend potential between two particles.
 *
 * The OutOfPlaneBendPotential struct encapsulates the type of out_of_plane_bend and associated parameters between two
 * particles. It includes versioning for serialization, out_of_plane_bend type, identifiers of out_of_plane_bended
 * particles, and out_of_plane_bend parameters.
 */
export struct OutOfPlaneBendPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the out_of_plane_bend.
  OutOfPlaneBendType type;                 ///< The type of out_of_plane_bend potential.
  std::array<double, maximumNumberOfOutOfPlaneBendParameters>
      parameters;  ///< Parameters associated with the out_of_plane_bend potential.

  /**
   * \brief Default constructor for OutOfPlaneBendPotential.
   *
   * Initializes a OutOfPlaneBendPotential object with Undefined out_of_plane_bend type and zeroed out_of_plane_bend
   * IDs.
   */
  OutOfPlaneBendPotential() : identifiers({0, 0, 0, 0}), type(OutOfPlaneBendType::Harmonic) {}

  OutOfPlaneBendPotential(std::array<std::size_t, 4> identifiers, OutOfPlaneBendType type,
                          std::vector<double> vector_parameters);

  /**
   * \brief Constructs a OutOfPlaneBendPotential with specified type and out_of_plane_bend IDs.
   *
   * \param type The type of out_of_plane_bend potential.
   * \param identifiers A pair of particle identifiers forming the out_of_plane_bend.
   */
  OutOfPlaneBendPotential(std::array<std::size_t, 4> identifiers, const OutOfPlaneBendType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(OutOfPlaneBendPotential const &) const = default;

  /**
   * \brief Generates a string representation of the out_of_plane_bend potential.
   *
   * Provides a formatted string containing out_of_plane_bend type, particle IDs, and parameters.
   *
   * \return A string describing the out_of_plane_bend potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each out_of_plane_bend type.
   *
   * A static vector indicating the number of parameters needed for each out_of_plane_bend type.
   */
  static inline std::array<std::size_t, 1> numberOfOutOfPlaneBendParameters{2};

  /**
   * \brief Mapping of out_of_plane_bend type strings to OutOfPlaneBendType enums.
   *
   * A static map that associates out_of_plane_bend type names with their corresponding OutOfPlaneBendType enumeration
   * values.
   */
  static inline std::map<std::string, OutOfPlaneBendType, caseInsensitiveComparator> definitionForString{
      {"HARMONIC", OutOfPlaneBendType::Harmonic}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const OutOfPlaneBendPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, OutOfPlaneBendPotential &b);
};
