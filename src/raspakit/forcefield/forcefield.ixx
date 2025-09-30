module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <set>
#include <string>
#include <vector>
#endif

#ifndef USE_LEGACY_HEADERS
#include <string>
#endif

export module forcefield;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double4;
import double3;
import int3;
import pseudo_atom;
import vdwparameters;
import json;
import simulationbox;
/**
 * \brief Represents the force field used in simulations.
 *
 * The ForceField struct contains all parameters and methods related to force field calculations,
 * including van der Waals interactions, electrostatics, Ewald summation parameters,
 * mixing rules, and methods to initialize and compute these parameters.
 */
export struct ForceField
{
  /**
   * \brief Enumeration of methods used for electrostatic interactions.
   */
  enum class ChargeMethod : int
  {
    Ewald = 0,        ///< Ewald summation method.
    Coulomb = 1,      ///< Direct Coulomb interactions.
    Wolf = 2,         ///< Wolf summation method.
    ModifiedWolf = 3  ///< Modified Wolf method.
  };

  /**
   * \brief Enumeration of mixing rules for cross interactions.
   */
  enum class MixingRule : int
  {
    Lorentz_Berthelot = 0,  ///< Lorentz-Berthelot mixing rule.
    Jorgensen = 1
  };

  enum class PotentialEnergySurfaceType : std::size_t
  {
    None = 0,
    SecondOrderPolynomialTestFunction = 1,
    ThirdOrderPolynomialTestFunction = 2,
    FourthOrderPolynomialTestFunction = 3,
    FifthOrderPolynomialTestFunction = 4,
    SixthOrderPolynomialTestFunction = 5,
    ExponentialNonPolynomialTestFunction = 6,
    MullerBrown = 7,
    Eckhardt = 8,
    GonzalezSchlegel = 9  // https://sci-hub.se/https://doi.org/10.1063/1.465995
  };

  enum class InterpolationGridType : std::size_t
  {
    LennardJones = 0,
    LennardJonesRepulsion = 1,
    LennardJonesAttraction = 2,
    EwaldReal = 3
  };

  enum class InterpolationScheme : std::size_t
  {
    Polynomial = 1,
    Tricubic = 8,
    Triquintic = 27
  };

  std::uint64_t versionNumber{1};  ///< Version number of the force field format.

  std::vector<VDWParameters>
      data{};  ///< Interaction parameters between pseudo-atoms; size is numberOfPseudoAtoms squared.
  std::vector<bool> shiftPotentials{};  ///< Indicates if potential shift is applied between pairs of atoms.
  std::vector<bool> tailCorrections{};  ///< Indicates if tail corrections are applied between pairs of atoms.
  MixingRule mixingRule{MixingRule::Lorentz_Berthelot};  ///< Mixing rule used for cross interactions.
  bool cutOffFrameworkVDWAutomatic{false};
  double cutOffFrameworkVDW{12.0};  ///< Cut-off distance for VDW interactions between framework and molecules.
  bool cutOffMoleculeVDWAutomatic{false};
  double cutOffMoleculeVDW{12.0};  ///< Cut-off distance for VDW interactions between molecules.
  bool cutOffCoulombAutomatic{true};
  double cutOffCoulomb{12.0};  ///< Cut-off distance for Coulomb interactions.
  double dualCutOff{6.0};      ///< Inner cut-off distance when using dual cut-off scheme.

  std::size_t numberOfPseudoAtoms{0};     ///< Number of pseudo-atoms defined in the force field.
  std::vector<PseudoAtom> pseudoAtoms{};  ///< List of pseudo-atoms in the force field.

  ChargeMethod chargeMethod{ChargeMethod::Ewald};  ///< Method used for calculating electrostatic interactions.

  double EwaldPrecision{1e-6};        ///< Desired precision for Ewald summation.
  double EwaldAlpha{0.265058};        ///< Ewald convergence parameter alpha.
  int3 numberOfWaveVectors{8, 8, 8};  ///< Number of wave vectors in each direction for Ewald summation.
  std::size_t reciprocalIntegerCutOffSquared{
      std::numeric_limits<std::size_t>::max()};  ///< Squared integer cut-off in reciprocal space.
  double reciprocalCutOffSquared{
      std::numeric_limits<double>::max()};  ///< Squared cut-off distance in reciprocal space.
  bool automaticEwald{true};                ///< Indicates if Ewald parameters are computed automatically.

  bool useCharge{true};          ///< Indicates if charges are used in calculations.
  bool omitEwaldFourier{false};  ///< If true, omits the Fourier component in Ewald summation.

  double energyOverlapCriteria{1e6};  ///< Energy criteria for considering overlaps.

  std::size_t numberOfTrialDirections{ 10 };
  std::size_t numberOfTorsionTrialDirections{ 100 };
  std::size_t numberOfFirstBeadPositions{ 10 };
  std::size_t numberOfTrialMovesPerOpenBead{ 150 };
  double minimumRosenbluthFactor{ 1e-150 };  ///< Minimum allowed Rosenbluth factor.

  bool useDualCutOff{false};          ///< Indicates if dual cut-off scheme is used.
  bool omitInterInteractions{false};  ///< If true, omits interactions between molecules.

  bool computePolarization{false};   ///< Indicates if polarization effects are computed.
  bool omitInterPolarization{true};  ///< If true, omits polarization between molecules.

  bool hasExternalField{false};
  PotentialEnergySurfaceType potentialEnergySurfaceType{PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction};
  double3 potentialEnergySurfaceOrigin{0.0, 0.0, 0.0};

  std::vector<std::size_t> gridPseudoAtomIndices;
  double spacingVDWGrid{0.15};
  double spacingCoulombGrid{0.15};
  std::optional<int3> numberOfVDWGridPoints{};
  std::optional<int3> numberOfCoulombGridPoints{};
  std::size_t numberOfGridTestPoints{100000};
  bool interpolationSchemeAuto{true};
  InterpolationScheme interpolationScheme{InterpolationScheme::Polynomial};

  /**
   * \brief Default constructor for the ForceField struct.
   */
  ForceField() noexcept = default;

  /**
   * \brief Constructs a ForceField with specified parameters.
   *
   * Initializes the force field using the provided pseudo-atoms, van der Waals parameters,
   * mixing rule, cut-off distances, and flags for shifting potentials and tail corrections.
   *
   * \param pseudoAtoms Vector of pseudo-atoms.
   * \param parameters Vector of van der Waals self-interaction parameters.
   * \param mixingRule The mixing rule to use for cross interactions.
   * \param cutOffFrameworkVDW Cut-off distance for VDW interactions between framework and molecules.
   * \param cutOffMoleculeVDW Cut-off distance for VDW interactions between molecules.
   * \param cutOffCoulomb Cut-off distance for Coulomb interactions.
   * \param shifted If true, applies potential shift to interactions.
   * \param tailCorrections If true, applies tail corrections to interactions.
   * \param useCharge If true, includes electrostatic interactions.
   */
  ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> parameters, MixingRule mixingRule,
             double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb, bool shifted,
             bool tailCorrections, bool useCharge = true) noexcept(false);

  /**
   * \brief Constructs a ForceField by reading parameters from a file.
   *
   * Initializes the force field by parsing a force field file specified by the file path.
   *
   * \param filePath Path to the force field file.
   */
  ForceField(std::string filePath) noexcept(false);

  VDWParameters &operator()(std::size_t row, std::size_t col) { return data[row * numberOfPseudoAtoms + col]; }
  const VDWParameters &operator[](std::size_t row) const { return data[row * numberOfPseudoAtoms + row]; }
  const VDWParameters &operator()(std::size_t row, std::size_t col) const
  {
    return data[row * numberOfPseudoAtoms + col];
  }
  bool operator==(const ForceField &other) const;

  /**
   * \brief Applies the mixing rule to compute cross-interaction parameters.
   *
   * Calculates the interaction parameters between different pseudo-atoms
   * using the specified mixing rule.
   */
  void applyMixingRule();

  /**
   * \brief Returns the cut-off distance for van der Waals interactions between two pseudo-atoms.
   *
   * \param i Index of the first pseudo-atom.
   * \param j Index of the second pseudo-atom.
   * \return The cut-off distance for VDW interactions.
   */
  double cutOffVDW(std::size_t i, std::size_t j) const;

  /**
   * \brief Pre-computes the potential shift for interactions.
   *
   * Calculates the potential shift for each pair of pseudo-atoms if shifting is enabled.
   */
  void preComputePotentialShift();

  /**
   * \brief Pre-computes the tail corrections for interactions.
   *
   * Calculates the tail correction energy for each pair of pseudo-atoms if tail corrections are enabled.
   */
  void preComputeTailCorrection();

  /**
   * \brief Reads a ForceField from a file.
   *
   * Attempts to read the force field parameters from the specified file.
   *
   * \param directoryName Optional directory name where the file is located.
   * \param forceFieldFileName Name of the force field file.
   * \return An optional ForceField object if successful.
   */
  static std::optional<ForceField> readForceField(std::optional<std::string> directoryName,
                                                  std::string forceFieldFileName) noexcept(false);

  /**
   * \brief Returns a string representation of the pseudo-atom status.
   *
   * Generates a detailed string containing information about the pseudo-atoms.
   *
   * \return A string representing the pseudo-atom status.
   */
  std::string printPseudoAtomStatus() const;

  std::string printCutOffAutoStatus() const;

  /**
   * \brief Returns a string representation of the force field status.
   *
   * Generates a detailed string containing information about the force field parameters.
   *
   * \return A string representing the force field status.
   */
  std::string printForceFieldStatus() const;

  /**
   * \brief Returns a JSON representation of the pseudo-atom status.
   *
   * Generates a JSON array containing information about the pseudo-atoms.
   *
   * \return A vector of JSON objects representing the pseudo-atoms.
   */
  std::vector<nlohmann::json> jsonPseudoAtomStatus() const;

  /**
   * \brief Returns a JSON representation of the force field status.
   *
   * Generates a JSON object containing information about the force field parameters.
   *
   * \return A JSON object representing the force field status.
   */
  nlohmann::json jsonForceFieldStatus() const;

  /**
   * \brief Finds the index of a pseudo-atom by name.
   *
   * Searches for a pseudo-atom with the given name and returns its index if found.
   *
   * \param name Name of the pseudo-atom.
   * \return Optional index of the pseudo-atom.
   */
  std::optional<std::size_t> findPseudoAtom(const std::string &name) const;

  /**
   * \brief Finds the index of a pseudo-atom by name in a given list.
   *
   * Searches for a pseudo-atom with the given name in the provided vector and returns its index if found.
   *
   * \param pseudoAtoms Vector of pseudo-atoms to search.
   * \param name Name of the pseudo-atom.
   * \return Optional index of the pseudo-atom.
   */
  static std::optional<std::size_t> findPseudoAtom(const std::vector<PseudoAtom> pseudoAtoms, const std::string &name);

  /**
   * \brief Initializes the Ewald parameters based on the simulation box.
   *
   * Calculates the Ewald alpha parameter and the number of wave vectors required for the desired precision.
   *
   * \param simulationBox The simulation box to use for initialization.
   */
  void initializeEwaldParameters(const SimulationBox &simulationBox);

  void initializeAutomaticCutOff(const SimulationBox &simulationBox);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceField &f);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceField &f);

  /**
   * \brief Returns a string representation of the ForceField.
   *
   * Combines the pseudo-atom status and force field status into a single string.
   *
   * \return A string representing the ForceField.
   */
  std::string repr() const { return printPseudoAtomStatus() + "\n" + printForceFieldStatus(); }

  /**
   * \struct InsensitiveCompare
   * \brief Comparator for case-insensitive string comparison.
   *
   * This structure provides a functor to compare two strings without considering their case.
   * It is used to enable case-insensitive lookups within sets of strings.
   */
  struct InsensitiveCompare
  {
    /**
     * \brief Compares two strings in a case-insensitive manner.
     *
     * This operator overload allows for the comparison of two strings without regard to
     * their case, facilitating case-insensitive sorting and lookup.
     *
     * \param a The first string to compare.
     * \param b The second string to compare.
     * \return `true` if `a` is lexicographically less than `b` (case-insensitive), `false` otherwise.
     */
    bool operator()(const std::string &a, const std::string &b) const
    {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
      return _stricmp(a.c_str(), b.c_str()) < 0;
#else
      return strcasecmp(a.c_str(), b.c_str()) < 0;
#endif
    }
  };

  void validateInput(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  // Static Member Variables

  /**
   * \brief Set of general option keys accepted in the input data.
   *
   * This set contains all the general configuration keys that are recognized
   * at the top level of the input JSON. It is used to validate the presence
   * of only known keys.
   */
  static const std::set<std::string, InsensitiveCompare> options;
};
