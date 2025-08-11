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

export module bond_bond_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bond_bond potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bond_bond potential.
 */
export const std::size_t maximumNumberOfBondBondParameters{3};

/**
 * \brief Enumeration of different bond_bond types.
 *
 * Specifies the type of bond_bond potential to be used in simulations.
 */
export enum class BondBondType : std::size_t { CVFF = 0, CFF = 1 };

/**
 * \brief Represents a bond_bond potential between two particles.
 *
 * The BondBondPotential struct encapsulates the type of bond_bond and associated parameters between two particles.
 * It includes versioning for serialization, bond_bond type, identifiers of bond_bonded particles, and bond_bond
 * parameters.
 */
export struct BondBondPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 3> identifiers;  ///< Identifiers of the two particles forming the bond_bond.
  BondBondType type;                       ///< The type of bond_bond potential.
  std::array<double, maximumNumberOfBondBondParameters>
      parameters;  ///< Parameters associated with the bond_bond potential.

  /**
   * \brief Default constructor for BondBondPotential.
   *
   * Initializes a BondBondPotential object with Undefined bond_bond type and zeroed bond_bond IDs.
   */
  BondBondPotential() : identifiers({0, 0, 0}), type(BondBondType::CVFF) {}

  BondBondPotential(std::array<std::size_t, 3> identifiers, BondBondType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BondBondPotential with specified type and bond_bond IDs.
   *
   * \param type The type of bond_bond potential.
   * \param identifiers A pair of particle identifiers forming the bond_bond.
   */
  BondBondPotential(std::array<std::size_t, 3> identifiers, const BondBondType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(BondBondPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bond_bond potential.
   *
   * Provides a formatted string containing bond_bond type, particle IDs, and parameters.
   *
   * \return A string describing the bond_bond potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bond_bond type.
   *
   * A static vector indicating the number of parameters needed for each bond_bond type.
   */
  static inline std::array<std::size_t, 2> numberOfBondBondParameters{3, 3};

  /**
   * \brief Mapping of bond_bond type strings to BondBondType enums.
   *
   * A static map that associates bond_bond type names with their corresponding BondBondType enumeration values.
   */
  static inline std::map<std::string, BondBondType, caseInsensitiveComparator> definitionForString{
      {"CVFF", BondBondType::CVFF}, {"CFF", BondBondType::CFF}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondBondPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondBondPotential &b);
};
