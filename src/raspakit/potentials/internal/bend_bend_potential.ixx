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

export module bend_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bend_bend potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bend_bend potential.
 */
export const std::size_t maximumNumberOfBendBendParameters{3};

/**
 * \brief Enumeration of different bend_bend types.
 *
 * Specifies the type of bend_bend potential to be used in simulations.
 */
export enum class BendBendType : std::size_t { CVFF = 0, CFF = 1, MM3 = 2 };

/**
 * \brief Represents a bend_bend potential between two particles.
 *
 * The BendBendPotential struct encapsulates the type of bend_bend and associated parameters between two particles.
 * It includes versioning for serialization, bend_bend type, identifiers of bend_bended particles, and bend_bend
 * parameters.
 */
export struct BendBendPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the bend_bend.
  BendBendType type;                       ///< The type of bend_bend potential.
  std::array<double, maximumNumberOfBendBendParameters>
      parameters;  ///< Parameters associated with the bend_bend potential.

  /**
   * \brief Default constructor for BendBendPotential.
   *
   * Initializes a BendBendPotential object with Undefined bend_bend type and zeroed bend_bend IDs.
   */
  BendBendPotential() : identifiers({0, 0, 0, 0}), type(BendBendType::CVFF) {}

  BendBendPotential(std::array<std::size_t, 4> identifiers, BendBendType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BendBendPotential with specified type and bend_bend IDs.
   *
   * \param type The type of bend_bend potential.
   * \param identifiers A pair of particle identifiers forming the bend_bend.
   */
  BendBendPotential(std::array<std::size_t, 4> identifiers, const BendBendType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(BendBendPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bend_bend potential.
   *
   * Provides a formatted string containing bend_bend type, particle IDs, and parameters.
   *
   * \return A string describing the bend_bend potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bend_bend type.
   *
   * A static vector indicating the number of parameters needed for each bend_bend type.
   */
  static inline std::array<std::size_t, 3> numberOfBendBendParameters{3, 3, 3};

  /**
   * \brief Mapping of bend_bend type strings to BendBendType enums.
   *
   * A static map that associates bend_bend type names with their corresponding BendBendType enumeration values.
   */
  static inline std::map<std::string, BendBendType, caseInsensitiveComparator> definitionForString{
      {"CVFF", BendBendType::CVFF}, {"CFF", BendBendType::CFF}, {"MM3", BendBendType::MM3}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendBendPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendBendPotential &b);
};
