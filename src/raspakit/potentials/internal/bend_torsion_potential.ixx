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

export module bend_torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bend_torsion potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bend_torsion potential.
 */
export const std::size_t maximumNumberOfBendTorsionParameters{4};

/**
 * \brief Enumeration of different bend_torsion types.
 *
 * Specifies the type of bend_torsion potential to be used in simulations.
 */
export enum class BendTorsionType : std::size_t {
  Smoothed = 0,
  SmoothedThreeCosine = 1,
  Nicholas = 2,
  CFF = 3,
  SmoothedCFF = 4,
  SmoothedCFF2 = 5,
  SmoothedCFF3 = 6,
  CVFF = 7
};

/**
 * \brief Represents a bend_torsion potential between two particles.
 *
 * The BendTorsionPotential struct encapsulates the type of bend_torsion and associated parameters between two
 * particles. It includes versioning for serialization, bend_torsion type, identifiers of bend_torsioned particles, and
 * bend_torsion parameters.
 */
export struct BendTorsionPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the bend_torsion.
  BendTorsionType type;                    ///< The type of bend_torsion potential.
  std::array<double, maximumNumberOfBendTorsionParameters>
      parameters;  ///< Parameters associated with the bend_torsion potential.

  /**
   * \brief Default constructor for BendTorsionPotential.
   *
   * Initializes a BendTorsionPotential object with Undefined bend_torsion type and zeroed bend_torsion IDs.
   */
  BendTorsionPotential() : identifiers({0, 0, 0, 0}), type(BendTorsionType::Smoothed) {}

  BendTorsionPotential(std::array<std::size_t, 4> identifiers, BendTorsionType type,
                       std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BendTorsionPotential with specified type and bend_torsion IDs.
   *
   * \param type The type of bend_torsion potential.
   * \param identifiers A pair of particle identifiers forming the bend_torsion.
   */
  BendTorsionPotential(std::array<std::size_t, 4> identifiers, const BendTorsionType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(BendTorsionPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bend_torsion potential.
   *
   * Provides a formatted string containing bend_torsion type, particle IDs, and parameters.
   *
   * \return A string describing the bend_torsion potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bend_torsion type.
   *
   * A static vector indicating the number of parameters needed for each bend_torsion type.
   */
  static inline std::array<std::size_t, 8> numberOfBendTorsionParameters{3, 3, 3, 3, 3, 3, 3, 3};

  /**
   * \brief Mapping of bend_torsion type strings to BendTorsionType enums.
   *
   * A static map that associates bend_torsion type names with their corresponding BendTorsionType enumeration values.
   */
  static inline std::map<std::string, BendTorsionType, caseInsensitiveComparator> definitionForString{
      {"SMOOTHED", BendTorsionType::Smoothed},          {"SMOOTHED_THREE_COSINE", BendTorsionType::SmoothedThreeCosine},
      {"NICHOLAS", BendTorsionType::Nicholas},          {"CFF", BendTorsionType::CFF},
      {"SMOOTHED_CFF", BendTorsionType::SmoothedCFF},   {"SMOOTHED_CFF2", BendTorsionType::SmoothedCFF2},
      {"SMOOTHED_CFF2", BendTorsionType::SmoothedCFF3}, {"CVFF", BendTorsionType::CVFF}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendTorsionPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendTorsionPotential &b);
};
