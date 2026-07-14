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
  Framework(const ForceField &forceField, std::string componentName,
            SimulationBox simulationBox, std::size_t spaceGroupHallNumber, const std::vector<Atom> &definedAtoms,
            const std::vector<Atom> &fractionalAtoms, int3 numberOfUnitCells) noexcept(false);

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
  Framework(const ForceField &forceField, std::string componentName,
            SimulationBox simulationBox, std::size_t spaceGroupHallNumber, 
            const std::vector<Atom> &definedAtoms, int3 numberOfUnitCells) noexcept(false);

  std::uint64_t versionNumber{2};  ///< Version number for serialization purposes.

  SimulationBox simulationBox;          ///< Simulation box defining the unit cell dimensions.
  std::size_t spaceGroupHallNumber{1};  ///< Space group number according to the Hall notation.
  int3 numberOfUnitCells{1, 1, 1};      ///< Number of unit cells in each dimension for the supercell.

  std::string name{};                         ///< Name of the framework component.
  std::string filename{};                     ///< File name of the framework.
  std::size_t numberOfComponents{1};

  bool rigid{true};  ///< Flag indicating if the framework is rigid.

  double mass{0.0};          ///< Total mass of the framework.
  double unitCellMass{0.0};  ///< Mass of the unit cell.

  double netCharge{0.0};                                       ///< Net charge of the framework.
  double smallestCharge{0.0};                                  ///< Smallest atomic charge in the framework.
  double largestCharge{0.0};                                   ///< Largest atomic charge in the framework.

  std::vector<Atom> definedAtoms{};            ///< Fractional Atoms defining the unit cell before symmetry operations.
  std::vector<Atom> fractionalUnitCellAtoms;   ///< Fractional atoms in the unit cell after applying symmetry operations.
  std::vector<Atom> unitCellAtoms;             ///< Cartesian atoms in the unit cell after applying symmetry operations.
  std::vector<Atom> atoms{};                   ///< All Cartesian atoms in the framework after constructing the supercell.
  std::unordered_set<std::size_t> uniqueAtomTypes{};

  ConnectivityTable connectivityTable{};                            ///< Periodic framework bond connectivity.
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};  ///< Indexed internal framework potentials.

  void determineUniqueAtomTypes();

  /**
   * Reads type-based intramolecular definitions and expands them over the framework supercell.
   */
  void readFrameworkDefinition(const ForceField &forceField, const std::string &definitionName);

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
  std::vector<double2> atomUnitCellLennardJonesPotentialParameters(const ForceField &forceField) const;

  std::optional<double> computeLargestNonOverlappingFreeRadius(const ForceField &forceField, double3 probe_position,
                                                               double well_depth_factor) const;
  bool computeVanDerWaalsRadiusOverlap(const ForceField &forceField, double3 probe_position) const;
  bool computeOverlap(const ForceField &forceField, double3 probe_position, double well_depth_factor,
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
  std::string printStatus(const ForceField &forceField) const;

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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Framework &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Framework &c);

  /**
   * \brief Returns a string representation of the Framework.
   *
   * Generates a string that includes the framework's ID, name, number of atoms, net charge,
   * mass, and other relevant information.
   *
   * \return A string representing the Framework.
   */
  std::string repr() const;

  static Framework makeFAU(const ForceField &forceField, int3 replicate = {1, 1, 1});
  static Framework makeITQ29(const ForceField &forceField, int3 replicate = {1, 1, 1});
  static Framework makeMFI(const ForceField &forceField, int3 replicate = {1, 1, 1});
  static Framework makeCHA(const ForceField &forceField, int3 replicate = {1, 1, 1});
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
