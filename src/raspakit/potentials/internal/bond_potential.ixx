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

export module bond_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bond potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bond potential.
 */
export const std::size_t maximumNumberOfBondParameters{4};

/**
 * \brief Enumeration of different bond types.
 *
 * Specifies the type of bond potential to be used in simulations.
 */
export enum class BondType : std::size_t {
  Fixed = 0,
  Harmonic = 1,
  CoreShellSpring = 2,
  Morse = 3,
  LJ_12_6 = 4,
  LennardJones = 5,
  Buckingham = 6,
  RestrainedHarmonic = 7,
  Quartic = 8,
  CFF_Quartic = 9,
  MM3 = 10
};

/**
 * \brief Represents a bond potential between two particles.
 *
 * The BondPotential struct encapsulates the type of bond and associated parameters between two particles.
 * It includes versioning for serialization, bond type, identifiers of bonded particles, and bond parameters.
 */
export struct BondPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 2> identifiers;                        ///< Identifiers of the two particles forming the bond.
  BondType type;                                                 ///< The type of bond potential.
  std::array<double, maximumNumberOfBondParameters> parameters;  ///< Parameters associated with the bond potential.

  /**
   * \brief Default constructor for BondPotential.
   *
   * Initializes a BondPotential object with Undefined bond type and zeroed bond IDs.
   */
  BondPotential() : identifiers({0, 0}), type(BondType::Harmonic) {}

  BondPotential(std::array<std::size_t, 2> identifiers, BondType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BondPotential with specified type and bond IDs.
   *
   * \param type The type of bond potential.
   * \param identifiers A pair of particle identifiers forming the bond.
   */
  BondPotential(std::array<std::size_t, 2> identifiers, const BondType type) : identifiers(identifiers), type(type) {}

  bool operator==(BondPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bond potential.
   *
   * Provides a formatted string containing bond type, particle IDs, and parameters.
   *
   * \return A string describing the bond potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bond type.
   *
   * A static vector indicating the number of parameters needed for each bond type.
   */
  static inline std::array<std::size_t, 11> numberOfBondParameters{1, 2, 1, 3, 2, 2, 3, 3, 4, 4, 2};

  /**
   * \brief Mapping of bond type strings to BondType enums.
   *
   * A static map that associates bond type names with their corresponding BondType enumeration values.
   */
  static inline std::map<std::string, BondType> definitionForString{
      {"FIXED", BondType::Fixed},
      {"HARMONIC", BondType::Harmonic},
      {"CORE_SHELL_SPRING", BondType::CoreShellSpring},
      {"MORSE", BondType::Morse},
      {"LJ_12_6", BondType::LJ_12_6},
      {"LENNARD_JONES", BondType::LennardJones},
      {"BUCKINGHAM", BondType::Buckingham},
      {"RESTRAINED_HARMONIC", BondType::RestrainedHarmonic},
      {"QUARTIC", BondType::Quartic},
      {"CFF_QUARTIC", BondType::CFF_Quartic},
      {"MM3", BondType::MM3}};

  double generateBondLength(RandomNumber &random, double beta) const;

  double calculateEnergy(const double3 &posA, const double3 &posB) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BondPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BondPotential &b);
};
