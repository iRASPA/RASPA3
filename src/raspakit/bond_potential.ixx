module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <type_traits>
#include <vector>
#endif

export module bond_potential;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <map>;
import <vector>;
import <array>;
import <fstream>;
import <type_traits>;
import <print>;
#endif

import stringutils;
import archive;

/**
 * \brief Maximum number of parameters allowed for bond potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bond potential.
 */
export const size_t maximumNumberOfBondParameters{4};

/**
 * \brief Enumeration of different bond types.
 *
 * Specifies the type of bond potential to be used in simulations.
 */
export enum class BondType : size_t { Undefined = 0, Fixed = 1, Rigid = 2, Harmonic = 3 };

/**
 * \brief Represents a bond potential between two particles.
 *
 * The BondPotential struct encapsulates the type of bond and associated parameters between two particles.
 * It includes versioning for serialization, bond type, identifiers of bonded particles, and bond parameters.
 */
export struct BondPotential
{
  uint64_t versionNumber{1};  ///< Version number for serialization.

  BondType bondType;                                             ///< The type of bond potential.
  std::pair<size_t, size_t> bondIds;                             ///< Identifiers of the two particles forming the bond.
  std::array<double, maximumNumberOfBondParameters> parameters;  ///< Parameters associated with the bond potential.

  /**
   * \brief Default constructor for BondPotential.
   *
   * Initializes a BondPotential object with Undefined bond type and zeroed bond IDs.
   */
  BondPotential() : bondType(BondType::Undefined), bondIds(0, 0) {}

  /**
   * \brief Constructs a BondPotential with specified type and bond IDs.
   *
   * \param bondType The type of bond potential.
   * \param bondIds A pair of particle identifiers forming the bond.
   */
  BondPotential(const BondType bondType, std::pair<size_t, size_t> bondIds) : bondType(bondType), bondIds(bondIds) {}

  bool operator==(BondPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bond potential.
   *
   * Provides a formatted string containing bond type, particle IDs, and parameters.
   *
   * \return A string describing the bond potential.
   */
  std::string print() const
  {
    switch (bondType)
    {
      case BondType::Fixed:
        return std::format("FIXED_BOND ({}-{})\n", bondIds.first, bondIds.second);
      case BondType::Rigid:
        return std::format("RIGID_BOND ({}-{})\n", bondIds.first, bondIds.second);
      case BondType::Harmonic:
        return std::format("HARMONIC_BOND ({}-{}): p_0/k_B={} [K/A^2], p_1={} [A]\n", bondIds.first, bondIds.second,
                           parameters[0], parameters[1]);
      default:
        return "Unknown potential";
    }
  }

  /**
   * \brief Number of parameters required for each bond type.
   *
   * A static vector indicating the number of parameters needed for each bond type.
   */
  static inline std::vector<size_t> numberOfBondParameters{0, 0, 2};

  /**
   * \brief Comparator struct for case-insensitive string comparison.
   *
   * Used for comparing bond type strings in a case-insensitive manner.
   */
  struct comp
  {
    /**
     * \brief Compares two strings case-insensitively.
     *
     * \param lhs Left-hand side string.
     * \param rhs Right-hand side string.
     * \return True if lhs is less than rhs.
     */
    bool operator()(const std::string &lhs, const std::string &rhs) const
    {
#if defined(_WIN32)
      return _stricmp(lhs.c_str(), rhs.c_str()) < 0;
#else
      return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
#endif
    }
  };

  /**
   * \brief Mapping of bond type strings to BondType enums.
   *
   * A static map that associates bond type names with their corresponding BondType enumeration values.
   */
  static inline std::map<std::string, BondType, comp> bondDefinitionForString{
      {"FIXED_BOND", BondType::Fixed}, {"RIGID_BOND", BondType::Rigid}, {"HARMONIC_BOND", BondType::Harmonic}};

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondPotential &b);
};
