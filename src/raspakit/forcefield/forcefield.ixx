module;

#include <string>

export module forcefield;

import std;

import archive;
import double4;
import double3;
import int3;
import pseudo_atom;
import vdwparameters;
import cmap_potential;
import json;
import simulationbox;
import potential_ewald_real_space_table;
export import forcefield_settings;

/**
 * \brief Represents the force field used in simulations.
 *
 * The ForceField struct holds the parameters that define the Hamiltonian: the pseudo-atoms, the pair potentials
 * between them (with mixing rule, truncation and tail corrections), the cut-offs, the charge method with its
 * Ewald parameters, and the external field. These are the only fields the energy, gradient and Hessian kernels
 * read. Everything that only steers how a simulation samples or evaluates this Hamiltonian (the CBMC sampling
 * parameters, the dual cut-off scheme, the interpolation grids) is kept apart in 'settings' (ForceFieldSettings)
 * and is read from the same force-field file.
 */
export struct ForceField
{
  /**
   * \brief Enumeration of methods used for electrostatic interactions.
   */
  enum class ChargeMethod : int
  {
    Ewald = 0,                 ///< Ewald summation method.
    Coulomb = 1,               ///< Direct Coulomb interactions.
    Wolf = 2,                  ///< Original damped shifted-potential Wolf summation.
    DampedShiftedForce = 3,    ///< Fennell-Gezelter damped shifted-force method.
    ModifiedShiftedForce = 4,  ///< Waibel-Feinler-Gross modified shifted-force method.
    ZeroDipole = 5             ///< Fukuda zero-dipole summation method.
  };

  /**
   * \brief Enumeration of mixing rules for cross interactions.
   */
  enum class MixingRule : int
  {
    Lorentz_Berthelot = 0,  ///< Lorentz-Berthelot mixing rule.
    Jorgensen = 1,          ///< Jorgensen mixing rule.
    SixthPower = 2          ///< Sixth-power (Waldman-Hagler) mixing rule for Class-II/CFF.
  };

  /**
   * \brief The global van der Waals truncation scheme of the force-field file ("TruncationMethod").
   *
   * 'Shifted' sets the per-pair 'shiftPotentials' flags; 'Switched' and 'ForceSwitched' convert every
   * Lennard-Jones pair of the tables to VDWParameters::Type::LennardJonesSwitched (the quintic potential switch of
   * OpenMM and GROMACS) or LennardJonesForceSwitched (the CHARMM force switch) in applySwitching(). The switching
   * starts at 'switchingDistance' ("SwitchingDistance", default: 2 Angstrom below the cutoff).
   */
  enum class TruncationMethod : int
  {
    Truncated = 0,
    Shifted = 1,
    Switched = 2,
    ForceSwitched = 3
  };

  enum class PotentialEnergySurfaceType : std::size_t
  {
    None = 0,
    GridFile = 1,
    SecondOrderPolynomialTestFunction = 2,
    ThirdOrderPolynomialTestFunction = 3,
    FourthOrderPolynomialTestFunction = 4,
    FifthOrderPolynomialTestFunction = 5,
    SixthOrderPolynomialTestFunction = 6,
    ExponentialNonPolynomialTestFunction = 7,
    MullerBrown = 8,
    Eckhardt = 9,
    GonzalezSchlegel = 10,  // https://sci-hub.se/https://doi.org/10.1063/1.465995
    CylinderX = 11,
    CylinderY = 12,
    CylinderZ = 13,
    RectangleX = 14,
    RectangleY = 15,
    RectangleZ = 16
  };

  std::uint64_t versionNumber{4};  ///< Version number of the force field format.

  std::vector<VDWParameters>
      data{};  ///< Interaction parameters between pseudo-atoms; size is numberOfPseudoAtoms squared.
  /// The pair parameters of the intramolecular 1-4 pairs (CHARMM-style force fields give a pseudo-atom separate
  /// Lennard-Jones parameters for its 1-4 interactions). Empty when the 1-4 pairs use the regular table 'data';
  /// otherwise a full table of numberOfPseudoAtoms squared, mixed, shifted and solute-scaled like 'data' (no tail
  /// corrections: the 1-4 pairs are counted explicitly). Read from "parameters14" of the self and binary
  /// interactions of the force-field file.
  std::vector<VDWParameters> data14{};
  std::vector<bool> shiftPotentials{};  ///< Indicates if potential shift is applied between pairs of atoms.
  std::vector<bool> tailCorrections{};  ///< Indicates if tail corrections are applied between pairs of atoms.
  MixingRule mixingRule{MixingRule::Lorentz_Berthelot};  ///< Mixing rule used for cross interactions.
  TruncationMethod truncationMethod{TruncationMethod::Truncated};  ///< The global truncation scheme.
  /// Where the switching function of the switched truncation methods starts [Angstrom]; 0 selects the default of
  /// 2 Angstrom below the (pair) cutoff.
  double switchingDistance{0.0};

  /// The CMAP correction maps defined in the force field ('CMAPs'); components refer to them by name.
  std::vector<CMAPMap> cmapMaps{};

  bool cutOffFrameworkVDWAutomatic{false};
  double cutOffFrameworkVDW{12.0};  ///< Cut-off distance for VDW interactions between framework and molecules.
  bool cutOffMoleculeVDWAutomatic{false};
  double cutOffMoleculeVDW{12.0};  ///< Cut-off distance for VDW interactions between molecules.
  bool cutOffCoulombAutomatic{true};
  double cutOffCoulomb{12.0};  ///< Cut-off distance for Coulomb interactions.

  double temperature{300.0};  ///< External temperature, used by temperature-dependent potentials (Feynman-Hibbs).

  std::size_t numberOfPseudoAtoms{0};     ///< Number of pseudo-atoms defined in the force field.
  std::vector<PseudoAtom> pseudoAtoms{};  ///< List of pseudo-atoms in the force field.

  ChargeMethod chargeMethod{ChargeMethod::Ewald};  ///< Method used for calculating electrostatic interactions.

  double EwaldPrecision{1e-6};        ///< Desired precision for Ewald summation.
  double EwaldAlpha{0.265058};        ///< Ewald convergence parameter alpha.
  double modifiedShiftedForceBeta{0.3};  ///< Exponential switching parameter beta (inverse length).
  int3 numberOfWaveVectors{8, 8, 8};  ///< Number of wave vectors in each direction for Ewald summation.
  std::size_t reciprocalIntegerCutOffSquared{
      std::numeric_limits<std::size_t>::max()};  ///< Squared integer cut-off in reciprocal space.
  double reciprocalCutOffSquared{
      std::numeric_limits<double>::max()};  ///< Squared cut-off distance in reciprocal space.
  bool automaticEwald{true};                ///< Indicates if Ewald parameters are computed automatically.

  /// Tabulated erfc(alpha r)/r for the Ewald real-space energy and gradient of fully coupled pairs (see
  /// Potentials::potentialCoulomb). Derived from 'EwaldAlpha' and 'cutOffCoulomb': rebuilt by
  /// updateEwaldRealSpaceTable, which every routine that changes either of them calls. Not serialized and not part
  /// of operator==. A stale table (alpha differs) is never used; the pair potential then falls back to the exact
  /// library functions.
  EwaldRealSpaceTable ewaldRealSpaceTable{};
  /// Set to false (before the table is (re)built) to evaluate the Ewald real-space term with the library erfc
  /// everywhere. The tabulated energy and gradient are accurate to 1e-9 relative, but their second derivative is
  /// only accurate to about 1e-6, so tests that finite-difference the energy against a closed-form Hessian through
  /// heavy cancellation need the exact evaluation. Not serialized.
  bool useEwaldRealSpaceTable{true};

  bool useCharge{true};          ///< Indicates if charges are used in calculations.
  bool omitEwaldFourier{false};  ///< If true, omits the Fourier component in Ewald summation.

  [[nodiscard]] bool usesEwaldFourier() const
  {
    return useCharge && chargeMethod == ChargeMethod::Ewald && !omitEwaldFourier;
  }

  /// Finite-cutoff shifted-potential charge methods (Wolf, damped-shifted-force, modified-shifted-force,
  /// zero-dipole) require the real-space per-atom self-energy and the intra-molecular exclusion / completion
  /// of the shifted pair sum. Plain Coulomb and the Ewald method (including the real-space-only debugging mode
  /// with omitted Fourier part) must NOT receive these corrections.
  [[nodiscard]] bool usesRealSpaceChargeCorrections() const
  {
    return useCharge && (chargeMethod == ChargeMethod::Wolf || chargeMethod == ChargeMethod::DampedShiftedForce ||
                         chargeMethod == ChargeMethod::ModifiedShiftedForce ||
                         chargeMethod == ChargeMethod::ZeroDipole);
  }

  double energyOverlapCriteria{1e6};  ///< Energy criteria for considering overlaps.

  bool omitInterInteractions{false};  ///< If true, omits interactions between molecules.

  bool computePolarization{false};   ///< Indicates if polarization effects are computed.
  bool omitInterPolarization{true};  ///< If true, omits polarization between molecules.

  // The external field: an analytic test surface or a potential read from a cube file. Which one and where it
  // sits are parameters of the Hamiltonian; whether it is evaluated through an interpolation grid is a setting.
  PotentialEnergySurfaceType potentialEnergySurfaceType{PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction};
  double3 potentialEnergySurfaceOrigin{0.0, 0.0, 0.0};
  std::string externalFieldGridFileName{ "external_field.cube" };
  double4 externalFieldGeometryParameters{5.0, 5.0, 0.0, 0.0};

  /// The sampling and numerical settings read from the force-field file (CBMC trial counts, recoil growth,
  /// dual cut-off, interpolation grids). No energy kernel reads these; see ForceFieldSettings.
  ForceFieldSettings settings{};

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

  /// Whether the 1-4 pairs have their own pair parameters ('data14').
  [[nodiscard]] bool hasPair14Parameters() const { return !data14.empty(); }

  /// The pair parameters of a 1-4 pair of the pseudo-atom types 'row' and 'col': the 1-4 table when the force
  /// field has one, the regular table otherwise.
  [[nodiscard]] const VDWParameters &pair14(std::size_t row, std::size_t col) const
  {
    return data14.empty() ? data[row * numberOfPseudoAtoms + col] : data14[row * numberOfPseudoAtoms + col];
  }

  /// The pair parameters of the pair of types (row, col): the 1-4 table for a 1-4 pair, the regular table
  /// otherwise.
  [[nodiscard]] const VDWParameters &pair(std::size_t row, std::size_t col, bool is14) const
  {
    return is14 ? pair14(row, col) : data[row * numberOfPseudoAtoms + col];
  }

  /// Creates the 1-4 table as a copy of the regular table when there is none yet (the entries are then
  /// overwritten by the 1-4 self interactions, and mixed by applyMixingRule).
  void ensurePair14Table();

  /// Sets the 1-4 self interaction of pseudo-atom 'type' (creating the 1-4 table when needed). The cross terms
  /// follow from the mixing rule: call applyMixingRule (and the pre-computations) afterwards.
  void setPair14SelfInteraction(std::size_t type, const VDWParameters &parameters);

  /**
   * \brief Applies the mixing rule to compute cross-interaction parameters.
   *
   * Calculates the interaction parameters between different pseudo-atoms
   * using the specified mixing rule.
   */
  void applyMixingRule();

  /// Applies the mixing rule to the cross terms of one pair table (the regular table or the 1-4 table).
  void mixTable(std::vector<VDWParameters> &table) const;

  /// Converts the Lennard-Jones pairs of both tables to the switched form of 'truncationMethod' (Switched,
  /// ForceSwitched) and clears their potential shift; the other truncation methods convert the switched forms back
  /// to plain Lennard-Jones. Called after mixing; the pre-computations must follow.
  void applySwitching();

  /// Sets the truncation method and the switching distance (0: default) and redoes the pre-computations.
  void setTruncationMethod(TruncationMethod method, double switchingDistanceVDW = 0.0);

  static std::string truncationMethodName(TruncationMethod method);

  /**
   * \brief Returns the cut-off distance for van der Waals interactions between two pseudo-atoms.
   *
   * \param i Index of the first pseudo-atom.
   * \param j Index of the second pseudo-atom.
   * \return The cut-off distance for VDW interactions.
   */
  double cutOffVDW(std::size_t i, std::size_t j) const;

  /**
   * \brief Pre-computes derived constants for each pair interaction.
   *
   * Computes per-pair derived constants: the Feynman-Hibbs temperature pre-factor, the
   * shifted-force cutoff constants, and the soft-core reference diameter used in the
   * continuous-fractional lambda-scaling. Must be called before preComputePotentialShift()
   * and preComputeTailCorrection(), and re-called when the temperature changes.
   */
  void preComputeDerivedParameters();

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
   * \brief Scales the interactions of a set of pseudo-atom types (solute tempering, REST2).
   *
   * Pair interactions between two solute types are scaled by lambda, between a solute type and
   * any other type by sqrt(lambda); the partial charges of the solute types by sqrt(lambda), so
   * that the Coulomb interactions follow the same pattern. The derived constants, potential
   * shifts and tail corrections are recomputed. 'soluteType' has one entry per pseudo-atom.
   */
  void scaleSoluteInteractions(const std::vector<bool> &soluteType, double lambda);

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

  /**
   * \brief Rebuilds 'ewaldRealSpaceTable' for the current 'EwaldAlpha' and 'cutOffCoulomb' when it is out of date.
   *
   * Must be called after any direct assignment to 'EwaldAlpha' or 'cutOffCoulomb' (the initialize routines and
   * the JSON reader do so themselves). Cheap when the table already matches.
   */
  void updateEwaldRealSpaceTable();

  static ChargeMethod chargeMethodFromString(const std::string &value);
  static std::string_view chargeMethodName(ChargeMethod method);

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
inline int SGA_stricmp(const char *a, const char *b) const {
  int ca, cb;
  do {
     ca = *reinterpret_cast<const unsigned char *>(a);
     cb = *reinterpret_cast<const unsigned char *>(b);
     ca = std::tolower(std::toupper(ca));
     cb = std::tolower(std::toupper(cb));
     a++;
     b++;
   } while (ca == cb && ca != '\0');
   return ca - cb;
}

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
      return SGA_stricmp(a.c_str(), b.c_str()) < 0;
#endif
    }
  };

  void validateInput(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  // Static Member Variables

  /**
   * \brief Set of general option keys accepted in the input data.
   *
   * This set contains the parameter keys that are recognized at the top level of the
   * force-field JSON; the keys of the sampling and numerical settings are listed in
   * 'ForceFieldSettings::options'. 'validateInput' accepts the union of both.
   */
  static const std::set<std::string, InsensitiveCompare> options;

  static ForceField makeZeoliteForceField(double rc = 12.0, bool shifted = true, bool tailCorrections = false, bool useEwald = false);
  static ForceField makeMetalOrganicFrameworkForceField(double rc = 12.0, bool shifted = true, bool tailCorrections = false, bool useEwald = false);
  static ForceField makeZeoPlusPlusForceField(double rc = 12.0, bool shifted = true, bool tailCorrections = false, bool useEwald = false);
};
