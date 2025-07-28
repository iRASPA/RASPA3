module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
#endif

export module cif_reader;

#ifndef USE_LEGACY_HEADERS
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
   * \brief Constructs a CIFReader with the given CIF content and force field.
   *
   * Initializes the CIFReader by setting up the scanner with the provided content
   * and processing the CIF data according to the specified force field.
   *
   * \param content The CIF file content as a string.
   * \param forceField The force field to be used for parsing atom types and interactions.
   */
  CIFReader(const std::string& content, const ForceField& forceField);

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
  void parseLoop([[maybe_unused]] std::string& string, const ForceField& forceField);

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

  Scanner _scanner;                                   ///< Scanner object for navigating through the CIF content.
  std::string::const_iterator _previousScanLocation;  ///< Iterator pointing to the previous scan location.

  std::vector<Atom> fractionalAtoms;                 ///< List of atoms with fractional coordinates.
  SimulationBox simulationBox;                       ///< The simulation box defined by cell parameters.
  std::optional<std::size_t> _spaceGroupHallNumber;  ///< Optional space group Hall number.
  double _a;                                         ///< Cell length a.
  double _b;                                         ///< Cell length b.
  double _c;                                         ///< Cell length c.
  double _alpha;                                     ///< Cell angle alpha in degrees.
  double _beta;                                      ///< Cell angle beta in degrees.
  double _gamma;                                     ///< Cell angle gamma in degrees.
};
