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

export module bond_torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bond_torsion potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bond_torsion potential.
 */
export const std::size_t maximumNumberOfBondTorsionParameters{4};

/**
 * \brief Enumeration of different bond_torsion types.
 *
 * Specifies the type of bond_torsion potential to be used in simulations.
 */
export enum class BondTorsionType : std::size_t { MM3 = 0 };

/**
 * \brief Represents a bond_torsion potential between two particles.
 *
 * The BondTorsionPotential struct encapsulates the type of bond_torsion and associated parameters between two
 * particles. It includes versioning for serialization, bond_torsion type, identifiers of bond_torsioned particles, and
 * bond_torsion parameters.
 */
export struct BondTorsionPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the bond_torsion.
  BondTorsionType type;                    ///< The type of bond_torsion potential.
  std::array<double, maximumNumberOfBondTorsionParameters>
      parameters;  ///< Parameters associated with the bond_torsion potential.

  /**
   * \brief Default constructor for BondTorsionPotential.
   *
   * Initializes a BondTorsionPotential object with Undefined bond_torsion type and zeroed bond_torsion IDs.
   */
  BondTorsionPotential() : identifiers({0, 0, 0, 0}), type(BondTorsionType::MM3) {}

  BondTorsionPotential(std::array<std::size_t, 4> identifiers, BondTorsionType type,
                       std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BondTorsionPotential with specified type and bond_torsion IDs.
   *
   * \param type The type of bond_torsion potential.
   * \param identifiers A pair of particle identifiers forming the bond_torsion.
   */
  BondTorsionPotential(std::array<std::size_t, 4> identifiers, const BondTorsionType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(BondTorsionPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bond_torsion potential.
   *
   * Provides a formatted string containing bond_torsion type, particle IDs, and parameters.
   *
   * \return A string describing the bond_torsion potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bond_torsion type.
   *
   * A static vector indicating the number of parameters needed for each bond_torsion type.
   */
  static inline std::array<std::size_t, 1> numberOfBondTorsionParameters{4};

  /**
   * \brief Mapping of bond_torsion type strings to BondTorsionType enums.
   *
   * A static map that associates bond_torsion type names with their corresponding BondTorsionType enumeration values.
   */
  static inline std::map<std::string, BondTorsionType, caseInsensitiveComparator> definitionForString{
      {"MM3", BondTorsionType::MM3}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondTorsionPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondTorsionPotential &b);
};
