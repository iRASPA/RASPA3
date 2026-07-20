module;

export module component;

import std;

import stringutils;
import archive;
import randomnumbers;
import int3;
import double3;
import double4;
import simd_quatd;
import averages;
import atom;
import molecule;
import forcefield;
import property_lambda_probability_histogram;
import simulationbox;
import property_widom;
import property_gibbs_widom;
import move_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import mc_moves_cputime;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import intra_molecular_potentials;
import chiral_center;
import connectivity_table;
export import fragment;
export import fragment_graph;
export import cbmc_growth_plan;
import json;
import cbmc_move_statistics;

/**
 * \brief Represents a component within the simulation system.
 *
 * The Component struct encapsulates all properties and behaviors associated with a
 * simulation component, such as adsorbates or cations. It includes various physical
 * properties, molecular structures, and simulation parameters necessary for conducting
 * simulations. The struct provides constructors for initializing components from files
 * or programmatically and includes methods for computing physical properties and
 * manipulating molecular positions.
 */
export struct Component
{
  /**
   * \brief Enumeration of component types.
   */
  enum class Type : std::size_t
  {
    Adsorbate = 0,  ///< Represents an adsorbate component.
    Cation = 1      ///< Represents a cation component.
  };

  /**
   * \brief Enumeration of molecular shapes.
   */
  enum class Shape : std::size_t
  {
    NonLinear = 0,  ///< Non-linear molecular shape.
    Linear = 1,     ///< Linear molecular shape.
    Point = 2       ///< Point-particle molecular shape.
  };

  enum class Chirality : std::size_t
  {
    NoChirality = 0,
    S_Chiral = 1,
    R_Chiral = 2
  };

  /**
   * \brief Default constructor for the Component struct.
   *
   * Initializes a Component object with default values.
   */
  Component();

  /**
   * \brief Constructs a Component from a file.
   *
   * Initializes a Component by reading data from a specified file, using the provided
   * force field and simulation parameters.
   *
   * \param type The type of the component (Adsorbate or Cation).
   * \param currentComponent The current component ID.
   * \param forceField The force field used for parsing atom types and interactions.
   * \param componentName The name of the component.
   * \param fileName Optional filename containing component data.
   * \param numberOfBlocks The number of simulation blocks.
   * \param numberOfLambdaBins The number of lambda bins for scaling.
   * \param systemProbabilities Move probabilities for the Monte Carlo simulation.
   * \param fugacityCoefficient Optional fugacity coefficient.
   * \param thermodynamicIntegration Flag indicating if thermodynamic integration is used.
   *
   * \throws std::runtime_error If the component file cannot be read or parsed.
   */
  Component(Component::Type type, std::size_t componentId, const ForceField &forceField,
            const std::string &componentName, std::optional<const std::string> fileName, std::size_t numberOfBlocks,
            std::size_t numberOfLambdaBins, const MCMoveProbabilities &systemProbabilities = MCMoveProbabilities(),
            std::optional<double> fugacityCoefficient = std::nullopt,
            bool thermodynamicIntegration = false) noexcept(false);

  /**
   * \brief Constructs a Component programmatically.
   *
   * Initializes a Component with specified physical properties and molecular structure.
   *
   * \param forceField The force field used for defining atom properties.
   * \param componentName The name of the component.
   * \param T_c The critical temperature of the component.
   * \param P_c The critical pressure of the component.
   * \param w The acentric factor of the component.
   * \param atomList A list of atoms defining the component's molecular structure.
   * \param numberOfBlocks The number of simulation blocks.
   * \param numberOfLambdaBins The number of lambda bins for scaling.
   * \param systemProbabilities Move probabilities for the Monte Carlo simulation.
   * \param fugacityCoefficient Optional fugacity coefficient.
   * \param thermodynamicIntegration Flag indicating if thermodynamic integration is used.
   *
   * \throws std::runtime_error If pseudo-atoms are not recognized or data is invalid.
   */
  Component(const ForceField &forceField, std::string componentName, double T_c, double P_c,
            double w, std::vector<Atom> definedAtoms, const ConnectivityTable &connectivityTable,
            const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::size_t numberOfBlocks,
            std::size_t numberOfLambdaBins, const MCMoveProbabilities &particleProbabilities = MCMoveProbabilities(),
            std::optional<double> fugacityCoefficient = std::nullopt,
            bool thermodynamicIntegration = false, std::vector<double4> blockingPockets = {}) noexcept(false);

  std::uint64_t versionNumber{4};  ///< Version number for serialization.

  Type type{0};  ///< Type of the component (Adsorbate or Cation).

  std::string name{};                         ///< Name of the component.
  std::optional<std::string> filenameData{};  ///< Optional filename containing component data.
  std::string filename{};                     ///< Filename associated with the component.


  bool rigid{true};                             ///< Flag indicating if the component is rigid.
  std::size_t translationalDegreesOfFreedom{};  ///< Number of translational degrees of freedom.
  std::size_t rotationalDegreesOfFreedom{};     ///< Number of rotational degrees of freedom.

  double criticalTemperature{0.0};  ///< Critical temperature of the component [K].
  double criticalPressure{0.0};     ///< Critical pressure of the component [Pa].
  double acentricFactor{0.0};       ///< Acentric factor of the component [-].
  double molFraction{1.0};          ///< Mole fraction of the component [-].
  bool swappable{false};            ///< Flag indicating if the component is swappable.
  double partialPressure{0.0};      ///< Partial pressure of the component [Pa].

  double totalMass{0.0};                        ///< Total mass of the component [kg].
  std::optional<double> fugacityCoefficient{};  ///< Optional fugacity coefficient [-].
  double amountOfExcessMolecules{0.0};          ///< Amount of excess molecules [-].
  double bulkFluidDensity{0.0};                 ///< Bulk fluid density [kg/m³].
  double compressibility{0.0};                  ///< Compressibility of the component [-].

  std::optional<double> idealGasRosenbluthWeight{};  ///< Optional ideal gas Rosenbluth weight [-].
  std::optional<double> idealGasEnergy{};            ///< Optional ideal gas energy [J].

  /// Partner component index for distance-biased ion-pair GCMC (Orkoulas & Panagiotopoulos).
  std::optional<std::size_t> pairComponentId{};
  /// Maximum ion-pair separation R_max for distance-biased pair insertion/deletion [m].
  std::optional<double> maximumPairDistance{};

  double netCharge{0.0};                                ///< Net charge of the component [e].
  std::size_t startingBead{0};                          ///< Starting bead index for simulations.
  std::vector<std::pair<Atom, double>> definedAtoms{};  ///< List of defined atoms and their masses.

  double3 inertiaVector{};         ///< Inertia vector of the component.
  double3 inverseInertiaVector{};  ///< Inverse of the inertia vector.
  Shape shapeType;                 ///< Shape type of the molecule.
  std::vector<Atom> atoms{};       ///< List of atoms in the component.

  ConnectivityTable connectivityTable{};                            ///< Connectivity table for the component.
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};  ///< List of internal potentials.

  /// The fragment decomposition of the molecule and its deterministic CBMC growth structure. Every
  /// molecule has one: a fully flexible molecule is all single-atom fragments, a fully rigid molecule
  /// is a single fragment covering all atoms, and a semi-flexible molecule mixes rigid-body fragments
  /// with single-atom (flexible) fragments. Flexible rings (simple, fused, or bridged) appear as the
  /// closure bonds of the graph. Replaces RASPA2's 'MoleculeGroup'/'atomGroupIds'.
  FragmentGraph fragmentGraph{};
  // Warm-start reference geometry for the CBMC internal Monte-Carlo: the most recently grown,
  // thermalized, non-overlapping conformation of this component. The internal MC that samples a
  // step's base conformation only sees the bonded energy, so it cannot relax a self-overlapping
  // cold seed; reusing the previous accepted conformation (rigidly re-oriented per step) gives a
  // good non-overlapping starting point. Persistent per-component scratch: written by the grow entry
  // points (on a non-const Component) and read by the operators; not serialized (a restart simply
  // cold-starts the first grow) and not thread-shared (one grow of a given component at a time).
  std::vector<Atom> warmStartConformation{};
  // Persistent scratch conformation of a flexible molecule grown in isolation (ideal-gas), i.e. drawn
  // from exp(-beta * U_intra). Kept between calls so the CBMC reinsertion Markov chain that produces
  // equilibrated ideal-gas conformations (System::equilibratedIdealGasConformation) stays warm.
  mutable std::vector<Atom> grownIdealGasAtoms{};
  std::vector<std::vector<std::size_t>> partialReinsertionFixedAtoms{};
  // Cache of CBMC growth plans keyed by the set of already-placed beads. A plan is deterministic
  // (it depends only on the molecule's topology and the placed set), and building one filters the
  // intramolecular potentials per step, so the plans for the common placed sets (the starting bead
  // and the partial-reinsertion fixed sets) are built once and reused for every grow and retrace.
  // Derived data: rebuilt on demand, never serialized. Marked 'mutable' so retraces on a
  // 'const Component&' can populate it.
  mutable std::map<std::vector<std::size_t>, std::vector<CBMC::GrowStep>> growthPlanCache{};
  std::vector<std::size_t> identityChanges{};
  std::vector<std::size_t> gibbsIdentityChanges{};

  std::size_t initialNumberOfMolecules{0};  ///< Initial number of molecules in the component.

  PropertyLambdaProbabilityHistogram lambdaGC;     ///< Lambda probability histogram for Grand-Canonical simulations.
  PropertyLambdaProbabilityHistogram lambdaGibbs;  ///< Lambda probability histogram for Gibbs simulations.
  PropertyLambdaProbabilityHistogram lambdaPairSwap;    ///< Lambda probability histogram for ion-pair CFCMC swaps.
  PropertyLambdaProbabilityHistogram lambdaPairSwapCB;  ///< Lambda probability histogram for ion-pair CB/CFCMC swaps.
  bool hasFractionalMolecule{false};               ///< Flag indicating if the component has fractional molecules.

  MCMoveProbabilities mc_moves_probabilities;  ///< Move probabilities for Monte Carlo simulations.
  MCMoveStatistics mc_moves_statistics;
  MCMoveCpuTime mc_moves_cputime;  ///< CPU time statistics for Monte Carlo moves.

  // CBMC internal-move statistics counters. Marked 'mutable' so that the trial-orientation
  // Monte-Carlo scheme, which only updates these counters, can run on a 'const Component&'.
  mutable std::vector<CBMCMoveStatistics> cbmc_moves_statistics;

  PropertyWidom averageRosenbluthWeights;            ///< Average Rosenbluth weights for Widom insertion.
  PropertyGibbsWidom averageGibbsRosenbluthWeights;  ///< Average Rosenbluth weights for Widom insertion.
  
  std::vector<double4> blockingPockets{};  ///< List of blocking pockets defined by position and radius.

  double lnPartitionFunction{0};  ///< Natural logarithm of the partition function [-].

  //MultiSiteIsotherm isotherm{};            ///< Isotherm information for the component.
  //double massTransferCoefficient{0.0};     ///< Mass transfer coefficient [1/s].
  //double axialDispersionCoefficient{0.0};  ///< Axial dispersion coefficient [m²/s].
  //bool isCarrierGas{false};                ///< Flag indicating if the component is a carrier gas.

  //std::size_t columnPressure{0};  ///< Column index for pressure data.
  //std::size_t columnLoading{1};   ///< Column index for loading data.
  //std::size_t columnError{2};     ///< Column index for error data.

  /**
   * \brief Enumeration of pressure scaling types.
   */
  enum class PressureScale
  {
    Log = 0,    ///< Logarithmic pressure scaling.
    Normal = 1  ///< Normal pressure scaling.
  };

  PressureScale pressureScale{PressureScale::Log};  ///< Pressure scaling type.

  /**
   * \brief Reads component data from a file.
   *
   * Parses component information from the specified file using the provided force field.
   *
   * \param forceField The force field used for parsing atom types and interactions.
   * \param fileName The name of the file containing component data.
   *
   * \throws std::runtime_error If the file cannot be found, opened, or parsed correctly.
   */
  void readComponent(std::size_t componentId, const ForceField &forceField, const std::string &fileName);

  /**
   * \brief Generates a string representing the component's current status.
   *
   * Compiles various properties and parameters of the component into a formatted string.
   *
   * \param forceField The force field used for interpreting atom types.
   * \return A string detailing the component's status.
   */
  std::string printStatus(std::size_t componentId, const ForceField &forceField, double inputPressure) const;


  /**
   * \brief Serializes the component's status to JSON.
   *
   * Converts the component's current status and properties into a JSON object.
   *
   * \return A JSON object representing the component's status.
   */
  nlohmann::json jsonStatus() const;

  /**
   * \brief Computes rigid body properties of the component.
   *
   * Calculates properties such as the center of mass, inertia tensor, and degrees of freedom
   * based on the defined atoms.
   */
  void computeRigidProperties();

  /**
   * \brief Computes the center of mass for a given list of atoms.
   *
   * Calculates the center of mass position based on the provided atom list.
   *
   * \param atom_list The list of atoms to compute the center of mass for.
   * \return The center of mass position as a double3 vector.
   */
  double3 computeCenterOfMass(std::span<Atom> atom_list) const;

  /**
   * \brief Builds the 'Molecule' record for freshly grown atom positions.
   *
   * Computes the mass-weighted center of mass. For a fully rigid component the orientation
   * quaternion is recovered by an orthogonal Procrustes fit of the body-frame reference geometry
   * ('atoms', COM-centered in the principal frame) onto the laboratory positions, following the
   * rigid-molecule convention p = com + buildRotationMatrixInverse(q) * atoms[i].position. For a
   * flexible or semi-flexible component the orientation is the identity (per-fragment rigid-body
   * states are tracked separately).
   *
   * \param moleculeAtoms The grown laboratory-frame atom positions of one molecule.
   * \return The molecule record (center of mass, orientation, mass, component id, atom count).
   */
  Molecule createMoleculeRecord(std::span<const Atom> moleculeAtoms) const;

  /**
   * \brief Rotates the positions of all atoms using a given quaternion.
   *
   * Applies a rotation to the atom positions based on the provided quaternion.
   *
   * \param q The quaternion representing the rotation.
   * \return A vector of atoms with updated positions after rotation.
   */
  std::vector<Atom> rotatePositions(const simd_quatd &q) const;

  /**
   * \brief Creates a copy of the component's atoms based on a molecule span.
   *
   * Copies and adjusts the positions of atoms relative to a given molecule.
   *
   * \param molecule A span representing the molecule to copy atoms from.
   * \return A vector of copied atoms with adjusted positions.
   */
  std::vector<Atom> copiedAtoms(std::span<Atom> molecule) const;

  /**
   * \brief Translates a molecule by a specified displacement.
   *
   * Moves the entire molecule by the given displacement vector.
   *
   * \param molecule The molecule to translate.
   * \param molecule_atoms A span representing the atoms of the molecule.
   * \param displacement The displacement vector to apply.
   * \return A pair containing the translated molecule and its corresponding atoms.
   */
  std::pair<Molecule, std::vector<Atom>> translate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                   double3 displacement) const;

  /**
   * \brief Rotates a molecule using a specified rotation quaternion.
   *
   * Applies a rotation to the entire molecule based on the provided quaternion.
   *
   * \param molecule The molecule to rotate.
   * \param molecule_atoms A span representing the atoms of the molecule.
   * \param rotation The quaternion representing the rotation to apply.
   * \return A pair containing the rotated molecule and its corresponding atoms.
   */
  std::pair<Molecule, std::vector<Atom>> rotate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                simd_quatd rotation) const;

  ConnectivityTable readConnectivityTable(std::size_t size,
                                          const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  std::vector<BondPotential> readBondPotentials(const ForceField &forceField,
                                                const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BendPotential> readBendPotentials(const ForceField &forceField,
                                                const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<TorsionPotential> readTorsionPotentials(const ForceField &forceField,
                                                      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  std::vector<UreyBradleyPotential> readUreyBradleyPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<InversionBendPotential> readInversionBendPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<OutOfPlaneBendPotential> readOutOfPlaneBendPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<TorsionPotential> readImproperTorsionPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BondBondPotential> readBondBondPotentials(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BondBendPotential> readBondBendPotentials(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BondTorsionPotential> readBondTorsionPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BendBendPotential> readBendBendPotentials(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<BendTorsionPotential> readBendTorsionPotentials(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  std::vector<VanDerWaalsPotential> readVanDerWaalsPotentials(
      const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<CoulombPotential> readCoulombPotentials(
      const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  std::vector<std::vector<std::size_t>> readPartialReinsertionFixedAtoms(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Reads the optional 'ChiralCenters' of a molecule: a list of atom-index quadruples
   *        [center, neighbor1, neighbor2, neighbor3].
   *
   * The handedness of each center is taken from the reference geometry: the sign of the signed
   * volume of the tetrahedron (neighbor1 - center, neighbor2 - center, neighbor3 - center). CBMC
   * growth preserves this parity: branch swaps and ring-conformation moves that would invert a
   * declared center are rejected. Throws when an index is out of range, the four indices are not
   * distinct, or the reference geometry of the center is (near-)planar so its parity is undefined.
   */
  std::vector<ChiralCenter> readChiralCenters(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Reads the optional 'RigidBodies' of a molecule: a list of atom-index lists, each a rigid
   *        unit. Every atom not listed becomes a single-atom (flexible) fragment.
   *
   * Returns the rigid-body atom lists; the fragment graph is built from them by 'buildFragmentGraph'.
   * Throws when a rigid body references an out-of-range atom, when an atom appears twice, or when a
   * rigid body is not a connected subgraph of the connectivity table.
   */
  std::vector<std::vector<std::size_t>> readRigidBodies(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Builds 'fragmentGraph' from the connectivity table, the given rigid bodies, the starting
   *        bead, and the component reference geometry ('atoms').
   *
   * Computes each fragment's rigid-body reference data (from the reference geometry so the body-fixed
   * positions match the offsets the CBMC growth places), the spanning tree, growth order, and closure
   * bonds, and updates the component translational/rotational degrees of freedom.
   */
  void buildFragmentGraph(const std::vector<std::vector<std::size_t>> &rigidBodies);

  /**
   * \brief Returns the (cached) deterministic CBMC growth plan starting from 'beadsAlreadyPlaced'.
   *
   * The plan is built once per distinct placed set and reused afterwards; grow and retrace of the
   * same move obtain the identical plan, as required for detailed balance. The returned reference
   * stays valid for the lifetime of the component (the cache never erases entries).
   */
  const std::vector<CBMC::GrowStep> &growthPlan(const std::vector<std::size_t> &beadsAlreadyPlaced) const;

  /**
   * \brief Returns whether all given atom indices lie inside one and the same rigid-body fragment.
   *
   * Used to skip bond/bend/torsion potentials whose atoms all belong to a single rigid fragment:
   * the geometry of a rigid fragment is fixed and must not be sampled.
   */
  bool isInsideRigidFragment(std::span<const std::size_t> ids) const;

  /**
   * \brief Returns the index of the rigid-body fragment that 'bead' belongs to, or nullopt when the
   *        bead is a single-atom (flexible) fragment.
   */
  std::optional<std::size_t> rigidFragmentContaining(std::size_t bead) const;

  /**
   * \brief Returns whether this component is semi-flexible: it has at least one rigid body and more
   *        than one fragment (so it is neither fully rigid nor fully flexible).
   */
  bool isSemiFlexible() const;

  /// Number of rigid-body fragments. Zero for fully flexible molecules; one for fully rigid molecules
  /// (but such a molecule is not semi-flexible).
  std::size_t numberOfRigidFragments() const;

  /**
   * \brief Regenerates the laboratory positions of the atoms of fragment 'fragmentIndex' from its
   *        rigid-body state, writing into 'moleculeAtoms' (local molecule atom order).
   */
  void regenerateFragmentAtoms(const GroupState &state, std::size_t fragmentIndex,
                               std::span<Atom> moleculeAtoms) const;

  /**
   * \brief Recovers the rigid-body state (center of mass and orientation) of fragment
   *        'fragmentIndex' from the current laboratory atom positions.
   *
   * Exact for perfectly rigid geometry (orthogonal Procrustes fit); used to initialize GroupState
   * from atom positions after CBMC growth, on restart, or in any path that only maintains atoms.
   */
  GroupState deriveFragmentState(std::size_t fragmentIndex, std::span<const Atom> moleculeAtoms) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Component &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Component &c);

  /**
   * \brief Returns a string representation of the Component.
   *
   * Generates a simple string indicating a test representation of the Component.
   *
   * \return A string representing the Component.
   */
  std::string repr() const;

  static Component makeMethane(const ForceField &forceField, std::size_t id = 0);
  static Component makeCO2(const ForceField &forceField, std::size_t id = 0, bool useCharges = true);
  static Component makeWater(const ForceField &forceField, std::size_t id = 0, bool useCharges = true);
  static Component makeIon(const ForceField &forceField, std::size_t id, std::string_view name, std::size_t type, double q);
};

/**
 * \brief Parses a list of parameters from a string.
 *
 * Template function that extracts a list of parameters of type T from a given string,
 * stopping at comment indicators or invalid entries.
 *
 * \tparam T The type of parameters to parse.
 * \param arguments The string containing the list of parameters.
 * \param lineNumber The line number for error reporting.
 * \return A vector containing the parsed parameters of type T.
 *
 * \throws std::runtime_error If no valid values are found before encountering a comment or invalid entry.
 */
template <typename T>
std::vector<T> parseListOfParameters(const std::string &arguments, std::size_t lineNumber)
{
  std::vector<T> list{};

  std::string str;
  std::istringstream ss(arguments);

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
    T value;
    std::istringstream s(str);
    if (s >> value)
    {
      list.push_back(value);
    }
    else
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
  };

  return list;
}
