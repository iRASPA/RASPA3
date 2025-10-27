module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <map>
#include <optional>
#include <ostream>
#include <print>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#endif

export module component;

#ifdef USE_STD_IMPORT
import std;
#endif

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
import isotherm;
import multi_site_isotherm;
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
import connectivity_table;
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
   * \brief Enumeration of growth types for the component.
   */
  enum class GrowType : std::size_t
  {
    Rigid = 0,     ///< Rigid growth type.
    Flexible = 1,  ///< Flexible growth type.
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
    NoCharility = 0,
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
  Component(Component::Type type, std::size_t currentComponent, const ForceField &forceField,
            const std::string &componentName, std::optional<const std::string> fileName, std::size_t numberOfBlocks,
            std::size_t numberOfLambdaBins, const MCMoveProbabilities &systemProbabilities = MCMoveProbabilities(),
            std::optional<double> fugacityCoefficient = std::nullopt,
            bool thermodynamicIntegration = false) noexcept(false);

  /**
   * \brief Constructs a Component programmatically.
   *
   * Initializes a Component with specified physical properties and molecular structure.
   *
   * \param componentId The unique identifier for the component.
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
  Component(std::size_t componentId, const ForceField &forceField, std::string componentName, double T_c, double P_c,
            double w, std::vector<Atom> definedAtoms, const ConnectivityTable &connectivityTable,
            const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::size_t numberOfBlocks,
            std::size_t numberOfLambdaBins, const MCMoveProbabilities &systemProbabilities = MCMoveProbabilities(),
            std::optional<double> fugacityCoefficient = std::nullopt,
            bool thermodynamicIntegration = false) noexcept(false);

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  Type type{0};          ///< Type of the component (Adsorbate or Cation).
  GrowType growType{0};  ///< Growth type of the component.

  std::size_t componentId{0};                 ///< Unique identifier for the component.
  std::string name{};                         ///< Name of the component.
  std::optional<std::string> filenameData{};  ///< Optional filename containing component data.
  std::string filename{};                     ///< Filename associated with the component.

  std::vector<double4> blockingPockets{};  ///< List of blocking pockets defined by position and radius.

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

  double netCharge{0.0};                                ///< Net charge of the component [e].
  std::size_t startingBead{0};                          ///< Starting bead index for simulations.
  std::vector<std::pair<Atom, double>> definedAtoms{};  ///< List of defined atoms and their masses.

  double3 inertiaVector{};         ///< Inertia vector of the component.
  double3 inverseInertiaVector{};  ///< Inverse of the inertia vector.
  Shape shapeType;                 ///< Shape type of the molecule.
  std::vector<Atom> atoms{};       ///< List of atoms in the component.

  ConnectivityTable connectivityTable{};                            ///< Connectivity table for the component.
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};  ///< List of internal potentials.
  std::vector<Atom> grownAtoms{};
  std::vector<std::vector<std::size_t>> partialReinsertionFixedAtoms{};

  std::size_t initialNumberOfMolecules{0};  ///< Initial number of molecules in the component.

  PropertyLambdaProbabilityHistogram lambdaGC;     ///< Lambda probability histogram for Grand-Canonical simulations.
  PropertyLambdaProbabilityHistogram lambdaGibbs;  ///< Lambda probability histogram for Gibbs simulations.
  bool hasFractionalMolecule{false};               ///< Flag indicating if the component has fractional molecules.

  MCMoveProbabilities mc_moves_probabilities;  ///< Move probabilities for Monte Carlo simulations.
  MCMoveStatistics mc_moves_statistics;
  MCMoveCpuTime mc_moves_cputime;  ///< CPU time statistics for Monte Carlo moves.

  std::vector<CBMCMoveStatistics> cbmc_moves_statistics;

  PropertyWidom averageRosenbluthWeights;  ///< Average Rosenbluth weights for Widom insertion.

  double lnPartitionFunction{0};  ///< Natural logarithm of the partition function [-].

  MultiSiteIsotherm isotherm{};            ///< Isotherm information for the component.
  double massTransferCoefficient{0.0};     ///< Mass transfer coefficient [1/s].
  double axialDispersionCoefficient{0.0};  ///< Axial dispersion coefficient [m²/s].
  bool isCarrierGas{false};                ///< Flag indicating if the component is a carrier gas.

  std::size_t columnPressure{0};  ///< Column index for pressure data.
  std::size_t columnLoading{1};   ///< Column index for loading data.
  std::size_t columnError{2};     ///< Column index for error data.

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
  void readComponent(const ForceField &forceField, const std::string &fileName);

  /**
   * \brief Generates a string representing the component's current status.
   *
   * Compiles various properties and parameters of the component into a formatted string.
   *
   * \param forceField The force field used for interpreting atom types.
   * \return A string detailing the component's status.
   */
  std::string printStatus(const ForceField &forceField, double inputPressure) const;

  /**
   * \brief Generates a string representing the breakthrough status of the component.
   *
   * Compiles breakthrough-related properties and parameters of the component into a formatted string.
   *
   * \return A string detailing the breakthrough status of the component.
   */
  std::string printBreakthroughStatus() const;

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
   * \brief Generates an equilibrated molecule randomly placed within a simulation box.
   *
   * Rotates and translates the molecule randomly to ensure equilibration within the simulation box.
   *
   * \param random A random number generator instance.
   * \param simulationBox The simulation box within which to place the molecule.
   * \return A pair containing the equilibrated molecule and its corresponding atoms.
   */
  std::pair<Molecule, std::vector<Atom>> equilibratedMoleculeRandomInBox(RandomNumber &random,
                                                                         const SimulationBox &simulationBox) const;

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

  std::vector<VanDerWaalsPotential> readVanDerWaalsPotentials(
      const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  std::vector<CoulombPotential> readCoulombPotentials(
      const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  std::vector<std::vector<std::size_t>> readPartialReinsertionFixedAtoms(
      const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

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
