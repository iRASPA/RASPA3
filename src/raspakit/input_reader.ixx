module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cctype>
#include <complex>
#include <cstddef>
#include <cstring>
#include <format>
#include <fstream>
#include <istream>
#include <locale>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#endif

#ifdef USE_STD_IMPORT
#include <string.h>
#endif

export module input_reader;

#ifdef USE_STD_IMPORT
import std;
#endif

import stringutils;

import int3;
import double3;
import threadpool;
import json;

import system;
import atom;
import framework;
import component;
import simulationbox;
import forcefield;
import loadings;
import enthalpy_of_adsorption;
import energy_status;
import averages;

/**
 * \struct InputDataSystem
 * \brief Structure to hold system-specific input data for simulations.
 *
 * This structure contains configuration parameters for sampling various statistical
 * properties such as radial distribution functions (RDF) and density grids within
 * a simulation system.
 */
struct InputDataSystem
{
  // Sampling the conventional radial distribution function (RDF)

  bool computeConventionalRadialDistributionFunction{false};  ///< Flag to enable computation of conventional RDF.
  std::size_t sampleConventionalRadialDistributionFunctionEvery{10};   ///< Interval for sampling conventional RDF.
  std::size_t writeConventionalRadialDistributionFunctionEvery{5000};  ///< Interval for writing conventional RDF data.
  std::size_t conventionalRadialDistributionFunctionHistogramSize{
      128};                                                ///< Number of bins in the conventional RDF histogram.
  double conventionalRadialDistributionFunctionRange{12};  ///< Range (in appropriate units) for the conventional RDF.

  // Sampling the radial distribution function (RDF)

  bool computeRadialDistributionFunction{false};             ///< Flag to enable computation of RDF.
  std::size_t sampleRadialDistributionFunctionEvery{10};     ///< Interval for sampling RDF.
  std::size_t writeRadialDistributionFunctionEvery{5000};    ///< Interval for writing RDF data.
  std::size_t radialDistributionFunctionHistogramSize{128};  ///< Number of bins in the RDF histogram.
  double radialDistributionFunctionRange{15.0};              ///< Range (in appropriate units) for the RDF.

  // Sampling the 3D density grid

  bool computeDensityGrid{false};           ///< Flag to enable computation of the density grid.
  std::size_t sampleDensityGridEvery{1};    ///< Interval for sampling the density grid.
  std::size_t writeDensityGridEvery{5000};  ///< Interval for writing density grid data.
  int3 densityGridSize{128, 128, 128};      ///< Dimensions of the density grid (x, y, z).
};

/**
 * \class InputReader
 * \brief Class responsible for reading and parsing simulation input files.
 *
 * The InputReader class handles the parsing of simulation parameters from input files,
 * validation of the input data, and initialization of simulation systems based on the
 * parsed configurations. It supports various simulation types and manages related
 * configurations such as threading, random seeds, and force fields.
 */
export struct InputReader
{
  /**
   * \enum SimulationType
   * \brief Enumeration of supported simulation types.
   */
  enum class SimulationType : int
  {
    MonteCarlo = 0,                  ///< Standard Monte Carlo simulation.
    MonteCarloTransitionMatrix = 1,  ///< Monte Carlo simulation using a transition matrix approach.
    MolecularDynamics = 2,           ///< Molecular Dynamics simulation.
    Minimization = 3,                ///< Energy minimization simulation.
    Test = 4,                        ///< Test simulation type for debugging or development.
    Breakthrough = 5,                ///< Breakthrough simulation for adsorption studies.
    MixturePrediction = 6,           ///< Simulation for predicting mixtures.
    Fitting = 7,                     ///< Simulation type for fitting parameters.
    ParallelTempering = 8            ///< Parallel Tempering simulation for enhanced sampling.
  };

  /**
   * \brief Constructs an InputReader object by loading and parsing the specified input file.
   *
   * This constructor reads the simulation configuration from the provided input file,
   * parses the JSON data, validates the input, and initializes simulation systems based
   * on the parsed configurations.
   *
   * \param inputFile The path to the input file containing simulation parameters.
   *
   * \throws std::runtime_error If the input file does not exist, cannot be parsed, or contains invalid data.
   */
  InputReader(const std::string inputFile);

  // Member Variables

  std::string inputString;  ///< Raw input string loaded from the input file.

  SimulationType simulationType{SimulationType::MonteCarlo};  ///< Type of simulation to be executed.

  std::size_t numberOfBlocks{5};  ///< Number of simulation blocks.

  std::size_t numberOfCycles{0};                ///< Total number of simulation cycles.
  std::size_t numberOfInitializationCycles{0};  ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles{0};   ///< Number of equilibration cycles.
  std::size_t printEvery{5000};                 ///< Interval for printing simulation progress.
  std::size_t writeBinaryRestartEvery{5000};    ///< Interval for writing binary restart files.
  std::size_t rescaleWangLandauEvery{5000};     ///< Interval for rescaling in Wang-Landau sampling.
  std::size_t optimizeMCMovesEvery{5000};       ///< Interval for optimizing Monte Carlo moves.
  std::size_t writeEvery{100};                  ///< Interval for writing simulation data.
  bool restartFromBinary{false};  ///< Flag to indicate if the simulation should restart from a binary file.
  std::string restartFromBinaryFileName{"restart_data.bin"};  ///< Filename of the binary restart-file

  std::optional<unsigned long long> randomSeed{std::nullopt};  ///< Optional random seed for reproducibility.

  std::size_t numberOfThreads{1};  ///< Number of threads to be used in the simulation.
  ThreadPool::ThreadingType threadingType{ThreadPool::ThreadingType::Serial};  ///< Type of threading to be used.

  ForceField forceField;          ///< Force field used for defining interactions in the simulation.
  std::vector<System> systems{};  ///< Vector of simulation systems configured for the simulation.

  std::string displayName{"Column"};  ///< Name used for display purposes.
  double temperature{-1.0};           ///< Simulation temperature.

  // Member Functions

  /**
   * \brief Parses molecular simulation configurations from the provided JSON data.
   *
   * This function extracts and initializes parameters specific to molecular simulations,
   * such as system configurations, components, and simulation-specific settings.
   *
   * \param parsed_data The JSON object containing the parsed input data.
   */
  void parseMolecularSimulations(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Parses fitting simulation configurations from the provided JSON data.
   *
   * This function extracts and initializes parameters specific to fitting simulations,
   * allowing for the adjustment of simulation parameters to fit experimental data.
   *
   * \param parsed_data The JSON object containing the parsed input data.
   */
  void parseFitting(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Parses mixture prediction simulation configurations from the provided JSON data.
   *
   * This function extracts and initializes parameters specific to mixture prediction
   * simulations, enabling the study of mixtures within the simulation environment.
   *
   * \param parsed_data The JSON object containing the parsed input data.
   */
  void parseMixturePrediction(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Parses breakthrough simulation configurations from the provided JSON data.
   *
   * This function extracts and initializes parameters specific to breakthrough simulations,
   * which are used to study adsorption phenomena in porous materials.
   *
   * \param parsed_data The JSON object containing the parsed input data.
   */
  void parseBreakthrough(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  /**
   * \brief Validates the parsed input data to ensure correctness and completeness.
   *
   * This function checks for the presence of required keys and the validity of their
   * corresponding values within the provided JSON data. It throws exceptions if any
   * discrepancies or missing information are detected.
   *
   * \param parsed_data The JSON object containing the parsed input data.
   *
   * \throws std::runtime_error If unknown input keys are found or required keys are missing.
   */
  void validateInput(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  void parseUnits(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

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

  // Static Member Variables

  /**
   * \brief Set of general option keys accepted in the input data.
   *
   * This set contains all the general configuration keys that are recognized
   * at the top level of the input JSON. It is used to validate the presence
   * of only known keys.
   */
  static const std::set<std::string, InsensitiveCompare> generalOptions;

  /**
   * \brief Set of system-specific option keys accepted in the input data.
   *
   * This set contains all the configuration keys that are recognized within
   * each system's configuration in the input JSON. It is used to validate
   * the presence of only known system-specific keys.
   */
  static const std::set<std::string, InsensitiveCompare> systemOptions;

  /**
   * \brief Set of component-specific option keys accepted in the input data.
   *
   * This set contains all the configuration keys that are recognized within
   * each component's configuration in the input JSON. It is used to validate
   * the presence of only known component-specific keys.
   */
  static const std::set<std::string, InsensitiveCompare> componentOptions;

  static const int3 parseExternalFieldGridDimensions(const std::string& filename);

  /**
   * \brief Parses external field grid based on the provided CUBE file.
   *
   * This function extracts parameters specific to external field grids,
   * allowing to use this data in construction of external field interpolation grid in simulations.
   *
   * \param filename The name of the CUBE file containing the external field data.
   * \return A pair containing the 3D grid of external field values and the grid dimensions.
   */
  static const std::vector<double> parseExternalFieldGridCube(const std::string& filename);
};