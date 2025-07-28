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

export module inversion_bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for inversion_bend potentials.
 *
 * Defines the maximum number of parameters that can be associated with a inversion_bend potential.
 */
export const std::size_t maximumNumberOfInversionBendParameters{4};

/**
 * \brief Enumeration of different inversion_bend types.
 *
 * Specifies the type of inversion_bend potential to be used in simulations.
 */
export enum class InversionBendType : std::size_t {
  Harmonic = 0,
  HarmonicCosine = 1,
  Planar = 2,
  Harmonic2 = 3,
  HarmonicCosine2 = 4,
  Planar2 = 5,
  MM3 = 6
};

/**
 * \brief Represents a inversion_bend potential between two particles.
 *
 * The InversionBendPotential struct encapsulates the type of inversion_bend and associated parameters between two
 * particles. It includes versioning for serialization, inversion_bend type, identifiers of inversion_bended particles,
 * and inversion_bend parameters.
 */
export struct InversionBendPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the inversion_bend.
  InversionBendType type;                  ///< The type of inversion_bend potential.
  std::array<double, maximumNumberOfInversionBendParameters>
      parameters;  ///< Parameters associated with the inversion_bend potential.

  /**
   * \brief Default constructor for InversionBendPotential.
   *
   * Initializes a InversionBendPotential object with Undefined inversion_bend type and zeroed inversion_bend IDs.
   */
  InversionBendPotential() : identifiers({0, 0, 0, 0}), type(InversionBendType::Harmonic) {}

  InversionBendPotential(std::array<std::size_t, 4> identifiers, InversionBendType type,
                         std::vector<double> vector_parameters);

  /**
   * \brief Constructs a InversionBendPotential with specified type and inversion_bend IDs.
   *
   * \param type The type of inversion_bend potential.
   * \param identifiers A pair of particle identifiers forming the inversion_bend.
   */
  InversionBendPotential(std::array<std::size_t, 4> identifiers, const InversionBendType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(InversionBendPotential const &) const = default;

  /**
   * \brief Generates a string representation of the inversion_bend potential.
   *
   * Provides a formatted string containing inversion_bend type, particle IDs, and parameters.
   *
   * \return A string describing the inversion_bend potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each inversion_bend type.
   *
   * A static vector indicating the number of parameters needed for each inversion_bend type.
   */
  static inline std::array<std::size_t, 7> numberOfInversionBendParameters{2, 2, 1, 2, 2, 1, 2};

  /**
   * \brief Mapping of inversion_bend type strings to InversionBendType enums.
   *
   * A static map that associates inversion_bend type names with their corresponding InversionBendType enumeration
   * values.
   */
  static inline std::map<std::string, InversionBendType, caseInsensitiveComparator> definitionForString{
      {"HARMONIC", InversionBendType::Harmonic},
      {"HARMONIC_COSINE", InversionBendType::HarmonicCosine},
      {"PLANAR", InversionBendType::Planar},
      {"HARMONIC2", InversionBendType::Harmonic2},
      {"HARMONIC_COSINE2", InversionBendType::HarmonicCosine2},
      {"PLANAR2", InversionBendType::Planar2},
      {"MM3", InversionBendType::MM3}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InversionBendPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InversionBendPotential &b);
};
