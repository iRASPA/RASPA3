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

export module torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for torsion potentials.
 *
 * Defines the maximum number of parameters that can be associated with a torsion potential.
 */
export const std::size_t maximumNumberOfTorsionParameters{6};

/**
 * \brief Enumeration of different torsion types.
 *
 * Specifies the type of torsion potential to be used in simulations.
 */
export enum class TorsionType : std::size_t {
  Fixed = 0,
  Harmonic = 1,
  HarmonicCosine = 2,
  ThreeCosine = 3,
  RyckaertBellemans = 4,
  TraPPE = 5,
  TraPPE_Extended = 6,
  ModifiedTraPPE = 7,
  CVFF = 8,
  CFF = 9,
  CFF2 = 10,
  OPLS = 11,
  MM3 = 12,
  FourierSeries = 13,
  FourierSeries2 = 14
};

/**
 * \brief Represents a torsion potential between two particles.
 *
 * The TorsionPotential struct encapsulates the type of torsion and associated parameters between two particles.
 * It includes versioning for serialization, torsion type, identifiers of torsioned particles, and torsion parameters.
 */
export struct TorsionPotential
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<std::size_t, 4> identifiers;  ///< Identifiers of the two particles forming the torsion.
  TorsionType type;                        ///< The type of torsion potential.
  std::array<double, maximumNumberOfTorsionParameters>
      parameters;  ///< Parameters associated with the torsion potential.

  /**
   * \brief Default constructor for TorsionPotential.
   *
   * Initializes a TorsionPotential object with Undefined torsion type and zeroed torsion IDs.
   */
  TorsionPotential() : identifiers({0, 0, 0, 0}), type(TorsionType::Harmonic) {}

  TorsionPotential(std::array<std::size_t, 4> identifiers, TorsionType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a TorsionPotential with specified type and torsion IDs.
   *
   * \param type The type of torsion potential.
   * \param identifiers A pair of particle identifiers forming the torsion.
   */
  TorsionPotential(std::array<std::size_t, 4> identifiers, const TorsionType type)
      : identifiers(identifiers), type(type)
  {
  }

  bool operator==(TorsionPotential const &) const = default;

  /**
   * \brief Generates a string representation of the torsion potential.
   *
   * Provides a formatted string containing torsion type, particle IDs, and parameters.
   *
   * \return A string describing the torsion potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each torsion type.
   *
   * A static vector indicating the number of parameters needed for each torsion type.
   */
  static inline std::array<std::size_t, 8> numberOfTorsionParameters{0, 2, 2, 1, 2, 2, 2, 1};

  /**
   * \brief Mapping of torsion type strings to TorsionType enums.
   *
   * A static map that associates torsion type names with their corresponding TorsionType enumeration values.
   */
  static inline std::map<std::string, TorsionType, caseInsensitiveComparator> definitionForString{
      {"Fixed", TorsionType::Fixed},
      {"HARMONIC", TorsionType::Harmonic},
      {"HARMONIC_COSINE", TorsionType::HarmonicCosine},
      {"THREE_COSINE", TorsionType::ThreeCosine},
      {"SIX_COSINE", TorsionType::RyckaertBellemans},
      {"TRAPPE", TorsionType::TraPPE},
      {"TRAPPE_EXTENDED", TorsionType::TraPPE_Extended},
      {"CVFF", TorsionType::CVFF},
      {"CFF", TorsionType::CFF},
      {"CFF2", TorsionType::CFF2},
      {"OPLS", TorsionType::OPLS},
      {"MM3", TorsionType::MM3},
      {"FOURIER_SERIES", TorsionType::FourierSeries},
      {"FOURIER_SERIES2", TorsionType::FourierSeries2}};

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, const double3 &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TorsionPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TorsionPotential &b);
};
