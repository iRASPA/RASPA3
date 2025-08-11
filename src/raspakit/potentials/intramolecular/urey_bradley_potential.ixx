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

export module urey_bradley_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for ureyBradley potentials.
 *
 * Defines the maximum number of parameters that can be associated with a ureyBradley potential.
 */
export const std::size_t maximumNumberOfUreyBradleyParameters{4};

/**
 * \brief Enumeration of different ureyBradley types.
 *
 * Specifies the type of ureyBradley potential to be used in simulations.
 */
export enum class UreyBradleyType : std::size_t {
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
 * \brief Represents a ureyBradley potential between two particles.
 *
 * The UreyBradleyPotential struct encapsulates the type of ureyBradley and associated parameters between two particles.
 * It includes versioning for serialization, ureyBradley type, identifiers of ureyBradleyed particles, and ureyBradley
 * parameters.
 */
export struct UreyBradleyPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 2> identifiers;  ///< Identifiers of the two particles forming the ureyBradley.
  UreyBradleyType type;                    ///< The type of ureyBradley potential.
  std::array<double, maximumNumberOfUreyBradleyParameters>
      parameters;  ///< Parameters associated with the ureyBradley potential.

  /**
   * \brief Default constructor for UreyBradleyPotential.
   *
   * Initializes a UreyBradleyPotential object with Undefined ureyBradley type and zeroed ureyBradley IDs.
   */
  UreyBradleyPotential() : identifiers({0, 0}), type(UreyBradleyType::Harmonic) {}

  UreyBradleyPotential(std::array<std::size_t, 2> identifiers, UreyBradleyType type,
                       std::vector<double> vector_parameters);

  /**
   * \brief Constructs a UreyBradleyPotential with specified type and ureyBradley IDs.
   *
   * \param type The type of ureyBradley potential.
   * \param identifiers A pair of particle identifiers forming the ureyBradley.
   */
  UreyBradleyPotential(std::array<std::size_t, 2> identifiers, const UreyBradleyType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(UreyBradleyPotential const &) const = default;

  /**
   * \brief Generates a string representation of the ureyBradley potential.
   *
   * Provides a formatted string containing ureyBradley type, particle IDs, and parameters.
   *
   * \return A string describing the ureyBradley potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each ureyBradley type.
   *
   * A static vector indicating the number of parameters needed for each ureyBradley type.
   */
  static inline std::array<std::size_t, 11> numberOfUreyBradleyParameters{1, 2, 1, 3, 2, 2, 3, 3, 4, 4, 2};

  /**
   * \brief Mapping of ureyBradley type strings to UreyBradleyType enums.
   *
   * A static map that associates ureyBradley type names with their corresponding UreyBradleyType enumeration values.
   */
  static inline std::map<std::string, UreyBradleyType, caseInsensitiveComparator> definitionForString{
      {"FIXED", UreyBradleyType::Fixed},
      {"HARMONIC", UreyBradleyType::Harmonic},
      {"CORE_SHELL_SPRING", UreyBradleyType::CoreShellSpring},
      {"MORSE", UreyBradleyType::Morse},
      {"LJ_12_6", UreyBradleyType::LJ_12_6},
      {"LENNARD_JONES", UreyBradleyType::LennardJones},
      {"BUCKINGHAM", UreyBradleyType::Buckingham},
      {"RESTRAINED_HARMONIC", UreyBradleyType::RestrainedHarmonic},
      {"QUARTIC", UreyBradleyType::Quartic},
      {"CFF_QUARTIC", UreyBradleyType::CFF_Quartic},
      {"MM3", UreyBradleyType::MM3}};

  double generateUreyBradleyLength(RandomNumber &random, double beta) const;

  double calculateEnergy(const double3 &posA, const double3 &posB) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const UreyBradleyPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, UreyBradleyPotential &b);
};
