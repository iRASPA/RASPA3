module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
#include <tuple>
#include <expected>
#endif

export module cif_reader;

#ifdef USE_STD_IMPORT
import std;
#endif

import scanner;
import characterset;
import atom;
import simulationbox;
import forcefield;

/**
 * \brief Handles the parsing of CIF (Crystallographic Information File) data.
 *
 * The CIFReader struct is responsible for reading and interpreting CIF content,
 * extracting relevant information such as atomic positions, simulation box parameters,
 * and symmetry details. It utilizes a Scanner to navigate through the CIF content
 * and populates structures like fractionalAtoms and simulationBox based on the parsed data.
 */
export struct CIFReader
{
  /**
   * \brief Enumeration for specifying the source of atomic charges.
   *
   * UseChargesFrom defines the source from which atomic charges should be obtained.
   * Options include using charges from pseudo-atoms, the CIF file, or by performing
   * charge equilibration calculations.
   */
  enum class UseChargesFrom : std::size_t
  {
    PseudoAtoms = 0,         ///< Use charges from pseudo-atoms defined in the force field.
    CIF_File = 1,            ///< Use charges specified in the CIF file.
    ChargeEquilibration = 2  ///< Compute charges using charge equilibration methods.
  };

  enum class ParseError : std::size_t
  {
    invalidInput = 0,
    invalidForceField = 1
  };

  /**
   * \brief Constructs a CIFReader with the given CIF content and force field.
   *
   * Initializes the CIFReader by setting up the scanner with the provided content
   * and processing the CIF data according to the specified force field.
   *
   * \param content The CIF file content as a string.
   * \param forceField The force field to be used for parsing atom types and interactions.
   */
  CIFReader(const std::string& content);

  static auto readCIFString(const std::string& content, const ForceField& forceField, CIFReader::UseChargesFrom useChargesFrom) ->
    std::expected<std::tuple<SimulationBox, std::size_t, std::vector<Atom>, std::vector<Atom>>, CIFReader::ParseError>;


  /**
   * \brief Parses a generic line from the CIF content.
   *
   * Processes a single line of CIF data, extracting relevant information based on the context.
   *
   * \param string The line of CIF content to parse.
   */
  void parseLine(std::string& string);

  /**
   * \brief Parses the audit section of the CIF content.
   *
   * Extracts audit-related information from the CIF data.
   *
   * \param string The audit line from the CIF content.
   */
  void parseAudit(std::string& string);

  /**
   * \brief Parses the chemical section of the CIF content.
   *
   * Extracts chemical composition and related information from the CIF data.
   *
   * \param string The chemical line from the CIF content.
   */
  void parseChemical(std::string& string);

  /**
   * \brief Parses the cell parameters from the CIF content.
   *
   * Extracts cell dimensions and angles to define the simulation box.
   *
   * \param string The cell parameter line from the CIF content.
   */
  void parseCell(std::string& string);

  /**
   * \brief Parses symmetry information from the CIF content.
   *
   * Extracts symmetry-related data, including space group information.
   *
   * \param string The symmetry line from the CIF content.
   */
  void parseSymmetry(std::string& string);

  /**
   * \brief Parses the name section of the CIF content.
   *
   * Extracts naming information from the CIF data.
   *
   * \param string The name line from the CIF content.
   */
  void parseName(std::string& string);

  /**
   * \brief Parses loop constructs within the CIF content.
   *
   * Handles loops in the CIF data, particularly those related to atomic sites.
   *
   * \param string The loop line from the CIF content.
   * \param forceField The force field used for interpreting atom types.
   */
  std::expected<void, ParseError> parseLoop([[maybe_unused]] std::string& string, const ForceField& forceField);

  /**
   * \brief Skips over comment lines in the CIF content.
   *
   * Advances the scanner past any comment lines to continue parsing relevant data.
   */
  void skipComment();

  /**
   * \brief Parses a value from the CIF content.
   *
   * Extracts a single value, ensuring it does not start a new data block or loop.
   *
   * \return An optional string containing the parsed value, or std::nullopt if none.
   */
  std::optional<std::string> parseValue();

  /**
   * \brief Scans and retrieves an integer value from the CIF content.
   *
   * Extracts an integer value, handling any non-alphanumeric characters appropriately.
   *
   * \return The scanned integer value.
   */
  std::size_t scanInt();

  /**
   * \brief Scans and retrieves a double value from the CIF content.
   *
   * Extracts a double-precision floating-point value from the CIF data.
   *
   * \return The scanned double value.
   */
  double scanDouble();

  /**
   * \brief Scans and retrieves a double value from a given string.
   *
   * Converts a string to a double, ignoring any non-digit characters.
   *
   * \param tempString The string to convert to a double.
   * \return The converted double value.
   */
  double scanDouble(std::string tempString);

  /**
   * \brief Scans and retrieves a string value from the CIF content.
   *
   * Extracts a string value up to the next newline character.
   *
   * \return An optional string containing the scanned value, or std::nullopt if none.
   */
  std::optional<std::string> scanString();

  static std::vector<Atom> expandDefinedAtomsToUnitCell(const SimulationBox &simulation_box,
                  std::size_t spaceGroupHallNumber, const std::vector<Atom> &definedAtoms);

  Scanner scanner;                                  ///< Scanner object for navigating through the CIF content.
  

  std::vector<Atom> fractionalAtoms;                ///< List of atoms with fractional coordinates.
  std::optional<std::size_t> spaceGroupHallNumber;  ///< Optional space group Hall number.
  double a;                                         ///< Cell length a.
  double b;                                         ///< Cell length b.
  double c;                                         ///< Cell length c.
  double alpha;                                     ///< Cell angle alpha in degrees.
  double beta;                                      ///< Cell angle beta in degrees.
  double gamma;                                     ///< Cell angle gamma in degrees.
};
