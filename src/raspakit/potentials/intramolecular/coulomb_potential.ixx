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

export module coulomb_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for Coulomb potentials.
 *
 * Defines the maximum number of parameters that can be associated with a Coulomb potential.
 */
export const std::size_t maximumNumberOfCoulombParameters{4};

/**
 * \brief Enumeration of different Coulomb types.
 *
 * Specifies the type of Coulomb potential to be used in simulations.
 */
export enum class CoulombType : std::size_t { Coulomb = 0 };

/**
 * \brief Represents a Coulomb potential between two particles.
 *
 * The CoulombPotential struct encapsulates the type of Coulomb and associated parameters between two particles.
 * It includes versioning for serialization, Coulomb type, identifiers of Coulombed particles, and Coulomb parameters.
 */
export struct CoulombPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 2> identifiers;  ///< Identifiers of the two particles forming the Coulomb.
  CoulombType type;                        ///< The type of Coulomb potential.
  double scaling;
  std::array<double, maximumNumberOfCoulombParameters>
      parameters;  ///< Parameters associated with the Coulomb potential.

  /**
   * \brief Default constructor for CoulombPotential.
   *
   * Initializes a CoulombPotential object with Undefined Coulomb type and zeroed Coulomb IDs.
   */
  CoulombPotential() : identifiers({0, 0}), type(CoulombType::Coulomb), scaling(1.0) {}

  CoulombPotential(std::array<std::size_t, 2> identifiers, CoulombType type, std::vector<double> vector_parameters,
                   double scaling);

  /**
   * \brief Constructs a CoulombPotential with specified type and Coulomb IDs.
   *
   * \param CoulombType The type of Coulomb potential.
   * \param CoulombIds A pair of particle identifiers forming the Coulomb.
   */
  CoulombPotential(std::array<std::size_t, 2> identifiers, const CoulombType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(CoulombPotential const &) const = default;

  /**
   * \brief Generates a string representation of the Coulomb potential.
   *
   * Provides a formatted string containing Coulomb type, particle IDs, and parameters.
   *
   * \return A string describing the Coulomb potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each Coulomb type.
   *
   * A static vector indicating the number of parameters needed for each Coulomb type.
   */
  static inline std::array<std::size_t, 1> numberOfCoulombParameters{1};

  /**
   * \brief Mapping of Coulomb type strings to CoulombType enums.
   *
   * A static map that associates Coulomb type names with their corresponding CoulombType enumeration values.
   */
  static inline std::map<std::string, CoulombType, caseInsensitiveComparator> definitionForString{
      {"COULOMB", CoulombType::Coulomb}};

  double calculateEnergy(const double3 &posA, const double3 &posB) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CoulombPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CoulombPotential &b);
};
