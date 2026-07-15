module;

export module framework;

import std;

import stringutils;
import archive;
import randomnumbers;
import int3;
import double3;
import double2;
import averages;
import atom;
import forcefield;
import simulationbox;
import property_widom;
import connectivity_table;
import intra_molecular_potentials;
import json;

export struct FrameworkIntraMolecularImageShifts
{
  std::vector<std::array<int3, 2>> bonds{};
  std::vector<std::array<int3, 3>> bends{};
  std::vector<std::array<int3, 4>> torsions{};
  std::vector<std::array<int3, 4>> improperTorsions{};
  std::vector<std::array<int3, 2>> vanDerWaals{};
  std::vector<std::array<int3, 2>> coulombs{};
};

/**
 * \brief A type-based intramolecular potential definition parsed from a flexible-framework definition.
 *
 * Stores the pseudo-atom types the potential applies to, the potential-type keyword (e.g. "HARMONIC"),
 * and its numeric parameters. These definitions are retained on the Framework so the concrete
 * intramolecular potentials can be re-derived on a different atom set (e.g. a reduced primitive cell)
 * without re-reading the definition file. \tparam N Number of atoms the potential couples (2 bond, 3 bend,
 * 4 torsion/improper).
 */
export template <std::size_t N>
struct FrameworkPotentialDefinition
{
  std::array<std::size_t, N> atomTypes{};
  std::string potentialType{};
  std::vector<double> parameters{};
};

/**
 * \brief Represents a framework in the simulation system.
 *
 * The Framework struct encapsulates the properties and behaviors of a framework
 * (e.g., a solid lattice or crystalline structure) within the simulation. It includes
 * information about the simulation box, atoms within the framework, and various
 * structural and force field parameters. The struct provides constructors for initializing
 * frameworks from files or programmatically, and methods to process and output framework data.
 */
export struct Framework
{
  /**
   * \brief Default constructor for the Framework struct.
   *
   * Initializes an empty Framework object with default values.
   */
  Framework();

  /**
   * \brief Constructs a Framework programmatically with specified parameters.
   *
   * Initializes a Framework using provided simulation box, space group number, defined atoms,
   * and unit cell information, allowing for programmatic creation of frameworks without file input.
   *
   * \param componentId Identifier for the framework component.
   * \param forceField Reference to the force field containing pseudo-atom definitions.
   * \param componentName Name of the framework component.
   * \param simulationBox Simulation box defining the unit cell dimensions.
   * \param spaceGroupHallNumber Space group number according to the Hall notation.
   * \param definedAtoms Vector of atoms defining the fractional positions and types within the asymmetric unit cell.
   * \param fractionalAtoms Vector of atoms defining the fractional positions and types within the unit cell.
   * \param numberOfUnitCells Number of unit cells in each dimension to construct the supercell.
   */
  Framework(const ForceField& forceField, std::string componentName, SimulationBox simulationBox,
            std::size_t spaceGroupHallNumber, const std::vector<Atom>& definedAtoms,
            const std::vector<Atom>& fractionalAtoms, int3 numberOfUnitCells) noexcept(false);

  /**
   * \brief Constructs a Framework programmatically with specified parameters.
   *
   * Initializes a Framework using provided simulation box, space group number, defined atoms,
   * and unit cell information, allowing for programmatic creation of frameworks without file input.
   *
   * \param componentId Identifier for the framework component.
   * \param forceField Reference to the force field containing pseudo-atom definitions.
   * \param componentName Name of the framework component.
   * \param simulationBox Simulation box defining the unit cell dimensions.
   * \param spaceGroupHallNumber Space group number according to the Hall notation.
   * \param definedAtoms Vector of atoms defining the fractional positions and types within the asymmetric unit cell.
   * \param numberOfUnitCells Number of unit cells in each dimension to construct the supercell.
   */
  Framework(const ForceField& forceField, std::string componentName, SimulationBox simulationBox,
            std::size_t spaceGroupHallNumber, const std::vector<Atom>& definedAtoms,
            int3 numberOfUnitCells) noexcept(false);

  std::uint64_t versionNumber{3};  ///< Version number for serialization purposes.

  SimulationBox simulationBox;          ///< Simulation box defining the unit cell dimensions.
  std::size_t spaceGroupHallNumber{1};  ///< Space group number according to the Hall notation.
  int3 numberOfUnitCells{1, 1, 1};      ///< Number of unit cells in each dimension for the supercell.

  std::string name{};      ///< Name of the framework component.
  std::string filename{};  ///< File name of the framework.
  std::size_t numberOfComponents{1};

  bool rigid{true};  ///< Flag indicating if the framework is rigid.

  double mass{0.0};          ///< Total mass of the framework.
  double unitCellMass{0.0};  ///< Mass of the unit cell.

  double netCharge{0.0};       ///< Net charge of the framework.
  double smallestCharge{0.0};  ///< Smallest atomic charge in the framework.
  double largestCharge{0.0};   ///< Largest atomic charge in the framework.

  std::vector<Atom> definedAtoms{};           ///< Fractional Atoms defining the unit cell before symmetry operations.
  std::vector<Atom> fractionalUnitCellAtoms;  ///< Fractional atoms in the unit cell after applying symmetry operations.
  std::vector<Atom> unitCellAtoms;            ///< Cartesian atoms in the unit cell after applying symmetry operations.
  std::vector<Atom> atoms{};  ///< All Cartesian atoms in the framework after constructing the supercell.
  std::unordered_set<std::size_t> uniqueAtomTypes{};

  ConnectivityTable connectivityTable{};                            ///< Periodic framework bond connectivity.
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};  ///< Indexed internal framework potentials.
  FrameworkIntraMolecularImageShifts intraMolecularImageShifts{};   ///< Periodic images for internal potentials.

  ///@{
  /// Retained nonbonded intra-framework generation options. Stored so the van der Waals pair list can
  /// be regenerated for a finalized simulation box and cutoff (e.g. with periodic-image replicas when a
  /// cell is smaller than twice the cutoff).
  bool excludeIntra12Interactions{true};
  bool excludeIntra13Interactions{true};
  bool excludeIntraBondInteractions{false};
  bool excludeIntraBendInteractions{false};
  double intra14VanDerWaalsScaling{0.0};
  double intra14ChargeChargeScaling{0.0};
  /// Name/path of the flexible-framework intramolecular definition (empty for rigid or programmatic
  /// frameworks). Retained so the intramolecular topology can be re-derived on a reduced primitive cell.
  std::string frameworkDefinitionName{};
  /// Parsed type-based potential definitions retained from the framework definition. Together with the
  /// exclusion/scaling options above these are the complete input to generateIntraMolecularPotentials(),
  /// so the intramolecular topology can be re-derived on any atom set without re-reading the file.
  std::vector<FrameworkPotentialDefinition<2>> bondDefinitions{};
  std::vector<FrameworkPotentialDefinition<3>> bendDefinitions{};
  std::vector<FrameworkPotentialDefinition<4>> torsionDefinitions{};
  std::vector<FrameworkPotentialDefinition<4>> improperTorsionDefinitions{};
  ///@}

  void determineUniqueAtomTypes();

  /**
   * \brief Derives connectivity and all intramolecular potentials from the stored type-based definitions.
   *
   * Detects the framework bond connectivity from covalent radii, matches every bond/bend/torsion/improper
   * against the retained potential definitions, and generates the nonbonded (van der Waals / Coulomb) pair
   * list honoring the retained exclusion and 1-4 scaling options. The periodic-image shifts for every
   * potential are recomputed for the current simulation box. Operates on the current \c atoms using the
   * stored \c bondDefinitions / \c bendDefinitions / \c torsionDefinitions / \c improperTorsionDefinitions,
   * so it can be re-run after reducing the framework to a primitive cell without touching the file system.
   *
   * \param forceField Force field providing pseudo-atom radii, van der Waals parameters, and charge usage.
   */
  void generateIntraMolecularPotentials(const ForceField& forceField);

  /**
   * \brief Rebuilds the intra-framework van der Waals pair list with periodic-image replicas if needed.
   *
   * When the smallest perpendicular width of \p box is at least twice \p cutOffFrameworkVDW a single
   * minimum image is sufficient and the list is left unchanged. Otherwise the van der Waals potentials
   * and their periodic-image shifts are regenerated so that every atom pair (including an atom with its
   * own periodic images) contributes one entry per replica cell that falls within the cutoff. The 1-2/1-3
   * exclusions and 1-4 scaling apply only to the covalently bonded (minimum) image; all other images are
   * treated as full van der Waals neighbors. The real-space Coulomb list is unaffected, since Ewald always
   * uses a half-box cutoff for which the minimum image is exact.
   *
   * \param forceField Force field providing the van der Waals parameters.
   * \param box Finalized simulation box.
   * \param cutOffFrameworkVDW Finalized framework van der Waals cutoff distance.
   */
  void regenerateVanDerWaalsImageList(const ForceField& forceField, const SimulationBox& box,
                                      double cutOffFrameworkVDW);

  /**
   * Reads type-based intramolecular definitions and expands them over the framework supercell.
   */
  void readFrameworkDefinition(const ForceField& forceField, const std::string& definitionName);

  /**
   * \brief Constructs the supercell by replicating the unit cell atoms.
   *
   * Generates the full framework structure by replicating the unit cell atoms according
   * to the specified number of unit cells in each dimension.
   */
  void makeSuperCell();

  std::vector<Atom> makeSuperCell(int3 numberOfCells) const;

  std::vector<double3> fractionalAtomPositionsUnitCell() const;
  std::vector<double3> cartesianAtomPositionsUnitCell() const;
  std::vector<double2> atomUnitCellLennardJonesPotentialParameters(const ForceField& forceField) const;

  std::optional<double> computeLargestNonOverlappingFreeRadius(const ForceField& forceField, double3 probe_position,
                                                               double well_depth_factor) const;
  bool computeVanDerWaalsRadiusOverlap(const ForceField& forceField, double3 probe_position) const;
  bool computeOverlap(const ForceField& forceField, double3 probe_position, double well_depth_factor,
                      std::size_t probe_type, std::make_signed_t<std::size_t> skip) const;

  /**
   * \brief Generates a string representation of the framework status.
   *
   * Creates a formatted string containing detailed information about the framework,
   * including atom types, positions, charges, and bond information.
   *
   * \param forceField Reference to the force field containing pseudo-atom definitions.
   * \return A string representing the framework status.
   */
  std::string printStatus(const ForceField& forceField) const;

  /**
   * \brief Generates a string representation of the breakthrough status.
   *
   * \return A string representing the breakthrough status.
   */
  std::string printBreakthroughStatus() const;

  /**
   * \brief Generates a JSON representation of the framework status.
   *
   * Creates a JSON object containing detailed information about the framework, suitable
   * for serialization or logging.
   *
   * \return A JSON object representing the framework status.
   */
  nlohmann::json jsonStatus() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Framework& c);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Framework& c);

  /**
   * \brief Returns a string representation of the Framework.
   *
   * Generates a string that includes the framework's ID, name, number of atoms, net charge,
   * mass, and other relevant information.
   *
   * \return A string representing the Framework.
   */
  std::string repr() const;

  static Framework makeFAU(const ForceField& forceField, int3 replicate = {1, 1, 1});
  static Framework makeITQ29(const ForceField& forceField, int3 replicate = {1, 1, 1});
  static Framework makeMFI(const ForceField& forceField, int3 replicate = {1, 1, 1});
  static Framework makeCHA(const ForceField& forceField, int3 replicate = {1, 1, 1});
};

/**
 * \brief Parses a string of parameters into a vector of values.
 *
 * Parses a whitespace-separated string of parameters and converts them into a vector of
 * values of type T. Parsing stops if a comment is encountered or if a value cannot be read.
 *
 * \tparam T The type of the values to parse (e.g., int, double).
 * \param arguments The string containing the parameters.
 * \param lineNumber The line number in the input file for error reporting.
 * \return A vector of parsed values of type T.
 *
 * \throws std::runtime_error If no values could be read.
 */
template <typename T>
std::vector<T> parseListOfParameters(const std::string& arguments, std::size_t lineNumber)
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
