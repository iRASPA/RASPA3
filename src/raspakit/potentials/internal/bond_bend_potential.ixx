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

export module bond_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bond_bend potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bond_bend potential.
 */
export const std::size_t maximumNumberOfBondBendParameters{5};

/**
 * \brief Enumeration of different bond_bend types.
 *
 * Specifies the type of bond_bend potential to be used in simulations.
 */
export enum class BondBendType : std::size_t {
  CVFF = 0,
  CFF = 1,
  MM3 = 2,
  TruncatedHarmonic = 3,
  ScreenedHarmonic = 4,
  ScreenedVessal = 5,
  TruncatedVessal = 7
};

/**
 * \brief Represents a bond_bend potential between two particles.
 *
 * The BondBendPotential struct encapsulates the type of bond_bend and associated parameters between two particles.
 * It includes versioning for serialization, bond_bend type, identifiers of bond_bended particles, and bond_bend
 * parameters.
 */
export struct BondBendPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the bond_bend.
  BondBendType type;                       ///< The type of bond_bend potential.
  std::array<double, maximumNumberOfBondBendParameters>
      parameters;  ///< Parameters associated with the bond_bend potential.

  /**
   * \brief Default constructor for BondBendPotential.
   *
   * Initializes a BondBendPotential object with Undefined bond_bend type and zeroed bond_bend IDs.
   */
  BondBendPotential() : identifiers({0, 0, 0, 0}), type(BondBendType::CVFF) {}

  BondBendPotential(std::array<std::size_t, 4> identifiers, BondBendType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BondBendPotential with specified type and bond_bend IDs.
   *
   * \param type The type of bond_bend potential.
   * \param identifiers A pair of particle identifiers forming the bond_bend.
   */
  BondBendPotential(std::array<std::size_t, 4> identifiers, const BondBendType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(BondBendPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bond_bend potential.
   *
   * Provides a formatted string containing bond_bend type, particle IDs, and parameters.
   *
   * \return A string describing the bond_bend potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bond_bend type.
   *
   * A static vector indicating the number of parameters needed for each bond_bend type.
   */
  static inline std::array<std::size_t, 7> numberOfBondBendParameters{5, 5, 4, 3, 4, 4, 4};

  /**
   * \brief Mapping of bond_bend type strings to BondBendType enums.
   *
   * A static map that associates bond_bend type names with their corresponding BondBendType enumeration values.
   */
  static inline std::map<std::string, BondBendType, caseInsensitiveComparator> definitionForString{
      {"CVFF", BondBendType::CVFF},
      {"CFF", BondBendType::CFF},
      {"MM3", BondBendType::MM3},
      {"TRUNCATED_HARMONIC", BondBendType::TruncatedHarmonic},
      {"SCREENED_HARMONIC", BondBendType::ScreenedHarmonic},
      {"SCREENED_VESSAL", BondBendType::ScreenedVessal},
      {"TRUNCATED_VESSAL", BondBendType::TruncatedVessal}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondBendPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondBendPotential &b);
};
