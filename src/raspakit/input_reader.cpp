module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <set>
#include <sstream>
#include <streambuf>
#include <vector>
#endif

module input_reader;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import stringutils;
import json;
import system;
import atom;
import pseudo_atom;
import framework;
import component;
import simulationbox;
import forcefield;
import double3;
import double4;
import units;
import sample_movies;
import threadpool;
import isotherm;
import multi_site_isotherm;
import pressure_range;
import mc_moves_probabilities;
import mc_moves_move_types;
import reaction;
import reactions;
import transition_matrix;
import property_conventional_rdf;
import property_rdf;
import property_density_grid;
import property_energy_histogram;
import property_number_of_molecules_histogram;
import property_msd;
import property_vacf;
import write_lammps_data;
import thermostat;

int3 parseInt3(const std::string& item, auto json)
{
  if (json.is_array())
  {
    if (json.size() != 3)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
    int3 value{};
    try
    {
      value.x = json[0].template get<std::int32_t>();
      value.y = json[1].template get<std::int32_t>();
      value.z = json[2].template get<std::int32_t>();
      return value;
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
  }
  throw std::runtime_error(
      std::format("[Input reader]: key '{}', value {} should be array of 3 integer  numbers\n", item, json.dump()));
}

double3 parseDouble3(const std::string& item, auto json)
{
  if (json.is_array())
  {
    if (json.size() != 3)
    {
      throw std::runtime_error(std::format(
          "[Input reader]: key '{}', value {} should be array of 3 floatng point numbers\n", item, json.dump()));
    }
    double3 value{};
    try
    {
      value.x = json[0].template get<double>();
      value.y = json[1].template get<double>();
      value.z = json[2].template get<double>();
      return value;
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(std::format(
          "[Input reader]: key '{}', value {} should be array of 3 floating point numbers\n", item, json.dump()));
    }
  }
  throw std::runtime_error(std::format(
      "[Input reader]: key '{}', value {} should be array of 3 floating point numbers\n", item, json.dump()));
}

template <typename T>
std::vector<T> parseList(std::size_t size, const std::string& item, auto json)
{
  if (json.is_array())
  {
    std::vector<T> values{};
    try
    {
      values = json.template get<std::vector<T>>();
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(std::format(
          "[Input reader (parseList)]: key '{}', value {} should be array of numbers\n", item, json.dump()));
    }

    // resize to 'size' using the last value to fill the new ones
    values.resize(size, values.back());

    return values;
  }
  throw std::runtime_error(
      std::format("[Input reader (parseList)]: key '{}', value {} should be array of numbers\n", item, json.dump()));
}

InputReader::InputReader(const std::string inputFile)  // : inputStream(inputFile)
{
  if (!std::filesystem::exists(inputFile))
  {
    throw std::runtime_error(std::format("[Input reader]: File '{}' not found\n", inputFile));
  }

  std::ifstream input("simulation.json");

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  try
  {
    parsed_data = nlohmann::json::parse(input);
  }
  catch (nlohmann::json::parse_error& ex)
  {
    std::cerr << "parse error at byte " << ex.byte << std::endl;
  }

  validateInput(parsed_data);

  parseUnits(parsed_data);

  if (parsed_data.contains("SimulationType") && parsed_data["SimulationType"].is_string())
  {
    std::string simulationTypeString = parsed_data["SimulationType"].get<std::string>();
    if (caseInSensStringCompare(simulationTypeString, "MonteCarlo"))
    {
      simulationType = SimulationType::MonteCarlo;
      parseMolecularSimulations(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "MonteCarloTransitionMatrix"))
    {
      simulationType = SimulationType::MonteCarloTransitionMatrix;
      parseMolecularSimulations(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "MolecularDynamics"))
    {
      simulationType = SimulationType::MolecularDynamics;
      parseMolecularSimulations(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "Minimization"))
    {
      simulationType = SimulationType::Minimization;
      parseMolecularSimulations(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "Breakthrough"))
    {
      simulationType = SimulationType::Breakthrough;
      parseBreakthrough(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "MixturePrediction"))
    {
      simulationType = SimulationType::MixturePrediction;
      parseMixturePrediction(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "Fitting"))
    {
      simulationType = SimulationType::Fitting;
      parseFitting(parsed_data);
    }
    else if (caseInSensStringCompare(simulationTypeString, "ParallelTempering"))
    {
      simulationType = SimulationType::ParallelTempering;
      parseMolecularSimulations(parsed_data);
    }
    else
    {
      throw std::runtime_error(
          std::format("[Input reader]: {} not a valid simulation type", simulationTypeString, parsed_data.dump()));
    }
  }
}

void InputReader::parseFitting([[maybe_unused]] const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data) {}

void InputReader::parseMixturePrediction([[maybe_unused]] const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
}

void InputReader::parseBreakthrough(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  // count number of systems
  if (!parsed_data.contains("Systems"))
  {
    throw std::runtime_error(
        std::format("[Input reader]: no system defined with keyword 'Systems' and value of array-type\n"));
  }
  std::size_t jsonNumberOfSystems = parsed_data["Systems"].size();
  if (jsonNumberOfSystems == 0)
  {
    throw std::runtime_error(std::format("[Input reader]: keyword 'Systems' has empty value of array-type\n"));
  }

  systems = std::vector<System>(jsonNumberOfSystems);

  // count number of components
  std::size_t jsonNumberOfComponents{};
  if (parsed_data.contains("Components"))
  {
    jsonNumberOfComponents = parsed_data["Components"].size();
  }

  // pre-allocate the vector of vector (jsonNumberOfSystems x jsonNumberOfComponents), i.e. for each system a list of
  // components
  std::vector<std::vector<Component>> jsonComponents(jsonNumberOfSystems,
                                                     std::vector<Component>(jsonNumberOfComponents));

  // Parse component options
  for (std::size_t componentId = 0; auto& [_, item] : parsed_data["Components"].items())
  {
    Component component{};

    if (!item.contains("Name"))
    {
      throw std::runtime_error(
          std::format("[Input reader]: component must have a key 'Name' with a value of string-type'\n"));
    }
    component.name = item["Name"].get<std::string>();

    // construct Component
    for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
    {
      jsonComponents[i][componentId] = component;
    }

    componentId++;
  }

  for (std::size_t systemId = 0; auto& [key, value] : parsed_data["Systems"].items())
  {
    if (!value.contains("Type"))
    {
      throw std::runtime_error(std::format("[Input reader]: system must have a key 'Type' with value 'Framework'\n"));
    }
    std::string typeString = value["Type"].get<std::string>();

    if (caseInSensStringCompare(typeString, "Framework"))
    {
      Framework framework{};

      // Parse framework options
      if (!value.contains("Name"))
      {
        throw std::runtime_error(
            std::format("[Input reader]: framework must have a key 'Name' with a value of string-type'\n"));
      }

      framework.name = value["Name"].get<std::string>();

      double heliumVoidFraction{1.0};
      if (value.contains("HeliumVoidFraction"))
      {
        heliumVoidFraction = value["HeliumVoidFraction"].get<double>();
      }

      if (!value.contains("ExternalTemperature"))
      {
        throw std::runtime_error(
            std::format("[Input reader]: framework must have a key 'ExternalTemperature' with a value of "
                        "floating-point-type'\n"));
      }
      double T = value["ExternalTemperature"].get<double>();

      std::optional<double> P{};
      if (value.contains("ExternalPressure"))
      {
        P = value["ExternalPressure"].get<double>();
      }
      if (value.contains("ChemicalPotential"))
      {
        P = std::exp(value["ChemicalPotential"].get<double>() / (Units::KB * T));
      }

      // create system
      systems[systemId] = System(systemId, T, P, heliumVoidFraction, {framework}, jsonComponents[systemId]);
    }

    systemId++;
  }
}

void InputReader::parseMolecularSimulations(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  std::size_t jsonNumberOfBlocks{5};

  // count number of systems
  if (!parsed_data.contains("Systems"))
  {
    throw std::runtime_error(
        std::format("[Input reader]: no system defined with keyword 'Systems' and value of array-type\n"));
  }
  std::size_t jsonNumberOfSystems = parsed_data["Systems"].size();
  if (jsonNumberOfSystems == 0)
  {
    throw std::runtime_error(std::format("[Input reader]: keyword 'Systems' has empty value of array-type\n"));
  }

  systems = std::vector<System>(jsonNumberOfSystems);

  // Read the local 'force_field.json' if present. This file will be used if no 'ForceField' keyword is specified per
  // system
  std::optional<std::string> directoryName{};
  if (parsed_data.contains("ForceField") && parsed_data["ForceField"].is_string())
  {
    directoryName = parsed_data["ForceField"].get<std::string>();
  }
  std::vector<std::optional<ForceField>> forceFields = std::vector<std::optional<ForceField>>(jsonNumberOfSystems);
  const std::optional<ForceField> standard = ForceField::readForceField(directoryName, "force_field.json");
  for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
  {
    forceFields[i] = standard;
  }

  std::size_t jsonNumberOfLambdaBins{41};
  if (parsed_data.contains("NumberOfLambdaBins") && parsed_data["NumberOfLambdaBins"].is_number_unsigned())
  {
    jsonNumberOfLambdaBins = parsed_data["NumberOfLambdaBins"].get<std::size_t>();
  }

  if (parsed_data.contains("RestartFromBinaryFile") && parsed_data["RestartFromBinaryFile"].is_boolean())
  {
    restartFromBinary = parsed_data["RestartFromBinaryFile"].get<bool>();
  }

  if (parsed_data.contains("RandomSeed") && parsed_data["RandomSeed"].is_number_unsigned())
  {
    randomSeed = parsed_data["RandomSeed"].get<unsigned long long>();
  }

  if (parsed_data.contains("NumberOfCycles") && parsed_data["NumberOfCycles"].is_number_unsigned())
  {
    numberOfCycles = parsed_data["NumberOfCycles"].get<std::size_t>();
  }

  if (parsed_data.contains("NumberOfInitializationCycles") &&
      parsed_data["NumberOfInitializationCycles"].is_number_unsigned())
  {
    numberOfInitializationCycles = parsed_data["NumberOfInitializationCycles"].get<std::size_t>();
  }

  if (parsed_data.contains("NumberOfEquilibrationCycles") &&
      parsed_data["NumberOfEquilibrationCycles"].is_number_unsigned())
  {
    numberOfEquilibrationCycles = parsed_data["NumberOfEquilibrationCycles"].get<std::size_t>();
  }

  if (parsed_data.contains("PrintEvery") && parsed_data["PrintEvery"].is_number_unsigned())
  {
    printEvery = parsed_data["PrintEvery"].get<std::size_t>();
  }

  if (parsed_data.contains("WriteBinaryRestartEvery") && parsed_data["WriteBinaryRestartEvery"].is_number_unsigned())
  {
    writeBinaryRestartEvery = parsed_data["WriteBinaryRestartEvery"].get<std::size_t>();
  }

  if (parsed_data.contains("RescaleWangLandauEvery") && parsed_data["RescaleWangLandauEvery"].is_number_unsigned())
  {
    rescaleWangLandauEvery = parsed_data["RescaleWangLandauEvery"].get<std::size_t>();
  }

  if (parsed_data.contains("OptimizeMCMovesEvery") && parsed_data["OptimizeMCMovesEvery"].is_number_unsigned())
  {
    optimizeMCMovesEvery = parsed_data["OptimizeMCMovesEvery"].get<std::size_t>();
  }

  if (parsed_data.contains("ThreadingType") && parsed_data["ThreadingType"].is_string())
  {
    std::string threadingTypeString = parsed_data["ThreadingType"].get<std::string>();
    if (caseInSensStringCompare(threadingTypeString, "Serial"))
    {
      threadingType = ThreadPool::ThreadingType::Serial;
    }
    if (caseInSensStringCompare(threadingTypeString, "ThreadPool"))
    {
      threadingType = ThreadPool::ThreadingType::ThreadPool;
    }
    if (caseInSensStringCompare(threadingTypeString, "OpenMP"))
    {
      threadingType = ThreadPool::ThreadingType::OpenMP;
    }
    if (caseInSensStringCompare(threadingTypeString, "GPU-Offload"))
    {
      threadingType = ThreadPool::ThreadingType::GPU_Offload;
    }
  }

  if (parsed_data.contains("NumberOfThreads") && parsed_data["NumberOfThreads"].is_number_unsigned())
  {
    numberOfThreads = parsed_data["NumberOfThreads"].get<std::size_t>();
    if (numberOfThreads > 1)
      threadingType = ThreadPool::ThreadingType::ThreadPool;
    else
      threadingType = ThreadPool::ThreadingType::Serial;
  }

  // count number of components
  std::size_t jsonNumberOfComponents{};
  if (parsed_data.contains("Components"))
  {
    jsonNumberOfComponents = parsed_data["Components"].size();
  }

  // pre-allocate the vector of vector (jsonNumberOfSystems x jsonNumberOfComponents), i.e. for each system a list of
  // components
  std::vector<std::vector<Component>> jsonComponents(jsonNumberOfSystems,
                                                     std::vector<Component>(jsonNumberOfComponents));

  // for each system, a list of how many molecules to create for each component
  std::vector<std::vector<std::size_t>> jsonCreateNumberOfMolecules(jsonNumberOfSystems,
                                                                    std::vector<std::size_t>(jsonNumberOfComponents));

  // Parse component options
  if (parsed_data.contains("Components"))
  {
    for (std::size_t componentId = 0; auto& [_, item] : parsed_data["Components"].items())
    {
      std::vector<MCMoveProbabilities> move_probabilities(jsonNumberOfSystems);

      if (!item.contains("Name"))
      {
        throw std::runtime_error(
            std::format("[Input reader]: component must have a key 'Name' with a value of string-type'\n"));
      }
      std::string jsonComponentName = item["Name"].get<std::string>();

      Component::Type componentType = Component::Type::Adsorbate;
      if (item.contains("Type") && item["Type"].is_string())
      {
        std::string typeString = item["Type"].get<std::string>();
        if (caseInSensStringCompare(typeString, "Adsorbate"))
        {
          componentType = Component::Type::Adsorbate;
        }
        else if (caseInSensStringCompare(typeString, "Cation"))
        {
          componentType = Component::Type::Cation;
        }
      }

      // Convenience notation listing the properties as a single value. These will then be taken for all systems.
      // ========================================================================================================

      if (item.contains("TranslationProbability") && item["TranslationProbability"].is_number_float())
      {
        double translationProbability = item["TranslationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::Translation, translationProbability);
        }
      }

      if (item.contains("RandomTranslationProbability") && item["RandomTranslationProbability"].is_number_float())
      {
        double randomTranslationProbability = item["RandomTranslationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::RandomTranslation, randomTranslationProbability);
        }
      }

      if (item.contains("RotationProbability") && item["RotationProbability"].is_number_float())
      {
        double rotationProbability = item["RotationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::Rotation, rotationProbability);
        }
      }

      if (item.contains("RandomRotationProbability") && item["RandomRotationProbability"].is_number_float())
      {
        double randomRotationProbability = item["RandomRotationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::RandomRotation, randomRotationProbability);
        }
      }

      if (item.contains("ReinsertionProbability") && item["ReinsertionProbability"].is_number_float())
      {
        double reinsertionCBMCProbability = item["ReinsertionProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::ReinsertionCBMC, reinsertionCBMCProbability);
        }
      }

      if (item.contains("SwapConventionalProbability") && item["SwapConventionalProbability"].is_number_float())
      {
        double swapProbability = item["SwapConventionalProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::Swap, swapProbability);
        }
      }

      if (item.contains("SwapProbability") && item["SwapProbability"].is_number_float())
      {
        double swapCBMCProbability = item["SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::SwapCBMC, swapCBMCProbability);
        }
      }

      if (item.contains("CFCMC_SwapProbability") && item["CFCMC_SwapProbability"].is_number_float())
      {
        double swapCFCMCProbability = item["CFCMC_SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::SwapCFCMC, swapCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_CBMC_SwapProbability") && item["CFCMC_CBMC_SwapProbability"].is_number_float())
      {
        double swapCBCFCMCProbability = item["CFCMC_CBMC_SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::SwapCBCFCMC, swapCBCFCMCProbability);
        }
      }

      if (item.contains("Gibbs_CFCMC_SwapProbability") && item["Gibbs_CFCMC_SwapProbability"].is_number_float())
      {
        double gibbsSwapCFCMCProbability = item["Gibbs_CFCMC_SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::GibbsSwapCFCMC, gibbsSwapCFCMCProbability);
        }
      }

      if (item.contains("GibbsSwapProbability") && item["GibbsSwapProbability"].is_number_float())
      {
        double gibbsSwapCBMCProbability = item["GibbsSwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::GibbsSwapCBMC, gibbsSwapCBMCProbability);
        }
      }

      if (item.contains("WidomProbability") && item["WidomProbability"].is_number_float())
      {
        double widomProbability = item["WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::Widom, widomProbability);
        }
      }

      if (item.contains("CFCMC_WidomProbability") && item["CFCMC_WidomProbability"].is_number_float())
      {
        double widomCFCMCProbability = item["CFCMC_WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::WidomCFCMC, widomCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_CBMC_WidomProbability") && item["CFCMC_CBMC_WidomProbability"].is_number_float())
      {
        double widomCBCFCMCProbability = item["CFCMC_CBMC_WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(MoveTypes::WidomCBCFCMC, widomCBCFCMCProbability);
        }
      }

      if (item.contains("CreateNumberOfMolecules") && item["CreateNumberOfMolecules"].is_number_integer())
      {
        std::size_t n = item["CreateNumberOfMolecules"].get<std::size_t>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonCreateNumberOfMolecules[i][componentId] = n;
        }
      }

      // construct Component
      for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        if (!forceFields[i].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }

        jsonComponents[i][componentId] =
            Component(componentType, componentId, forceFields[i].value(), jsonComponentName, jsonComponentName,
                      jsonNumberOfBlocks, jsonNumberOfLambdaBins, move_probabilities[i]);
      }

      if (item.contains("StartingBead") && item["StartingBead"].is_number_integer())
      {
        std::size_t n = item["StartingBead"].get<std::size_t>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          if (n >= jsonComponents[i][componentId].definedAtoms.size())
          {
            throw std::runtime_error(std::format("[Input reader]: starting bead larger than the molecule size'\n"));
          }

          jsonComponents[i][componentId].startingBead = n;
        }
      }

      if (item.contains("FugacityCoefficient") && item["FugacityCoefficient"].is_number_float())
      {
        double fugacity_coefficient = item["FugacityCoefficient"].get<double>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].fugacityCoefficient = fugacity_coefficient;
        }
      }

      if (item.contains("IdealGasRosenbluthWeight") && item["IdealGasRosenbluthWeight"].is_number_float())
      {
        double ideal_gas_rosenbluth_weight = item["IdealGasRosenbluthWeight"].get<double>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].idealGasRosenbluthWeight = ideal_gas_rosenbluth_weight;
        }
      }

      if (item.contains("MolFraction") && item["MolFraction"].is_number_float())
      {
        double mol_fraction = item["MolFraction"].get<double>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].molFraction = mol_fraction;
        }
      }

      if (item.contains("ThermodynamicIntegration") && item["ThermodynamicIntegration"].is_boolean())
      {
        bool thermodynamic_integration = item["ThermodynamicIntegration"].get<bool>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].lambdaGC.computeDUdlambda = thermodynamic_integration;
        }
      }

      if (item.contains("LambdaBiasFileName") && item["LambdaBiasFileName"].is_string())
      {
        std::string lambda_bias_file_name = item["LambdaBiasFileName"].get<std::string>();

        try
        {
          std::filesystem::path path(lambda_bias_file_name);
          for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
          {
            jsonComponents[i][componentId].lambdaGC.readBiasingFile(path);
          }
        }
        catch (nlohmann::json::parse_error& ex)
        {
          std::cerr << "parse error at byte " << ex.byte << std::endl;
        }
      }

      if (item.contains("BlockingPockets") && item.contains("BlockingPockets"))
      {
        for (auto& [_, block_pockets_item] : item["BlockingPockets"].items())
        {
          if (!block_pockets_item.is_array())
          {
            throw std::runtime_error(
                std::format("[Component reader]: item {} must be an array\n", block_pockets_item.dump()));
          }

          if (block_pockets_item.size() != 4)
          {
            throw std::runtime_error(
                std::format("[Component reader]: item {} must be an array with four elements, "
                            "an array with the x,y,z positions, and a radius\n",
                            block_pockets_item.dump()));
          }

          std::vector<double> data =
              block_pockets_item.is_array() ? block_pockets_item.get<std::vector<double>>() : std::vector<double>{};
          for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
          {
            jsonComponents[i][componentId].blockingPockets.push_back(double4(data[0], data[1], data[2], data[3]));
          }
        }
      }

      // Explicit notation listing the properties as an array of the values for the particular systems
      // ========================================================================================================

      if (item.contains("FugacityCoefficient") && item["FugacityCoefficient"].is_array())
      {
        std::vector<double> fugacity_coefficients =
            parseList<double>(jsonNumberOfSystems, "FugacityCoefficient", item["FugacityCoefficient"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].fugacityCoefficient = fugacity_coefficients[i];
        }
      }

      if (item.contains("IdealGasRosenbluthWeight") && item["IdealGasRosenbluthWeight"].is_array())
      {
        std::vector<double> ideal_gas_rosenbluth_weight =
            parseList<double>(jsonNumberOfSystems, "IdealGasRosenbluthWeight", item["IdealGasRosenbluthWeight"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].idealGasRosenbluthWeight = ideal_gas_rosenbluth_weight[i];
        }
      }

      if (item.contains("MolFraction") && item["MolFraction"].is_array())
      {
        std::vector<double> mol_fractions = parseList<double>(jsonNumberOfSystems, "MolFraction", item["MolFraction"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].molFraction = mol_fractions[i];
        }
      }

      if (item.contains("CreateNumberOfMolecules") && item["CreateNumberOfMolecules"].is_array())
      {
        std::vector<std::size_t> initialNumberOfMolecule =
            parseList<std::size_t>(jsonNumberOfSystems, "CreateNumberOfMolecules", item["CreateNumberOfMolecules"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonCreateNumberOfMolecules[i][componentId] = initialNumberOfMolecule[i];
        }
      }

      if (item.contains("LambdaBiasFileName") && item["LambdaBiasFileName"].is_array())
      {
        std::vector<std::string> lambda_bias_file_names =
            parseList<std::string>(jsonNumberOfSystems, "LambdaBiasFileName", item["LambdaBiasFileName"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          try
          {
            std::filesystem::path path(lambda_bias_file_names[i]);
            jsonComponents[i][componentId].lambdaGC.readBiasingFile(path);
          }
          catch (nlohmann::json::parse_error& ex)
          {
            std::cerr << "parse error at byte " << ex.byte << std::endl;
          }
        }
      }

      componentId++;
    }
  }

  if (parsed_data.contains("Systems"))
  {
    for (std::size_t systemId = 0; auto& [key, value] : parsed_data["Systems"].items())
    {
      MCMoveProbabilities mc_moves_probabilities{};

      if (value.contains("ForceField") && value["ForceField"].is_string())
      {
        std::string name = parsed_data["ForceField"].get<std::string>();
        forceFields[systemId] = ForceField::readForceField(name, "force_field.json");
      }

      if (value.contains("CutOff") && value["CutOff"].is_number_float())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->cutOffFrameworkVDW = value["CutOff"].get<double>();
        forceFields[systemId]->cutOffMoleculeVDW = value["CutOff"].get<double>();
        forceFields[systemId]->preComputePotentialShift();
        forceFields[systemId]->preComputeTailCorrection();
      }

      if (value.contains("CutOffVDW") && value["CutOffVDW"].is_number_float())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->cutOffFrameworkVDW = value["CutOffVDW"].get<double>();
        forceFields[systemId]->cutOffMoleculeVDW = value["CutOffVDW"].get<double>();
        forceFields[systemId]->preComputePotentialShift();
        forceFields[systemId]->preComputeTailCorrection();
      }

      if (value.contains("CutOffFrameworkVDW") && value["CutOffFrameworkVDW"].is_number_float())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->cutOffFrameworkVDW = value["CutOffFrameworkVDW"].get<double>();
        forceFields[systemId]->preComputePotentialShift();
        forceFields[systemId]->preComputeTailCorrection();
      }

      if (value.contains("CutOffMoleculeVDW") && value["CutOffMoleculeVDW"].is_number_float())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->cutOffMoleculeVDW = value["CutOffMoleculeVDW"].get<double>();
        forceFields[systemId]->preComputePotentialShift();
        forceFields[systemId]->preComputeTailCorrection();
      }

      if (value.contains("CutOffCoulomb") && value["CutOffCoulomb"].is_number_float())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->cutOffCoulomb = value["CutOffCoulomb"].get<double>();
      }

      if (value.contains("OmitEwaldFourier") && value["OmitEwaldFourier"].is_boolean())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->omitEwaldFourier = value["OmitEwaldFourier"].get<bool>();
      }

      if (value.contains("ComputePolarization") && value["ComputePolarization"].is_boolean())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->computePolarization = value["ComputePolarization"].get<bool>();
      }

      if (value.contains("ChargeMethod") && value["ChargeMethod"].is_string())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }

        std::string chargeMethodString = value["ChargeMethod"].get<std::string>();

        if (caseInSensStringCompare(chargeMethodString, "Ewald"))
        {
          forceFields[systemId]->chargeMethod = ForceField::ChargeMethod::Ewald;
          forceFields[systemId]->useCharge = true;
        }
        if (caseInSensStringCompare(chargeMethodString, "None"))
        {
          forceFields[systemId]->chargeMethod = ForceField::ChargeMethod::Ewald;
          forceFields[systemId]->useCharge = false;
        }
      }

      Framework::UseChargesFrom useChargesFrom{Framework::UseChargesFrom::PseudoAtoms};
      if (value.contains("UseChargesFrom") && value["UseChargesFrom"].is_string())
      {
        std::string useChargesFromString = value["UseChargesFrom"].get<std::string>();

        if (caseInSensStringCompare(useChargesFromString, "PsuedoAtoms"))
        {
          useChargesFrom = Framework::UseChargesFrom::PseudoAtoms;
        }
        if (caseInSensStringCompare(useChargesFromString, "CIF_File"))
        {
          useChargesFrom = Framework::UseChargesFrom::CIF_File;
        }
        if (caseInSensStringCompare(useChargesFromString, "ChargeEquilibration"))
        {
          useChargesFrom = Framework::UseChargesFrom::ChargeEquilibration;
        }
      }

      if (value.contains("VolumeMoveProbability") && value["VolumeMoveProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(MoveTypes::VolumeChange, value["VolumeMoveProbability"].get<double>());
      }

      if (value.contains("GibbsVolumeMoveProbability") && value["GibbsVolumeMoveProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(MoveTypes::GibbsVolume,
                                              value["GibbsVolumeMoveProbability"].get<double>());
      }

      if (value.contains("ParallelTemperingSwapProbability") &&
          value["ParallelTemperingSwapProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(MoveTypes::ParallelTempering,
                                              value["ParallelTemperingSwapProbability"].get<double>());
      }
      if (value.contains("HybridMCProbability") && value["HybridMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(MoveTypes::HybridMC, value["HybridMCProbability"].get<double>());
      }

      if (!value.contains("Type"))
      {
        throw std::runtime_error(
            std::format("[Input reader]: system must have a key 'Type' with value 'Box' or 'Framework'\n"));
      }
      std::string typeString = value["Type"].get<std::string>();

      if (!value.contains("ExternalTemperature"))
      {
        throw std::runtime_error(
            std::format("[Input reader]: framework must have a key 'ExternalTemperature' with a value of "
                        "floating-point-type'\n"));
      }
      double T = 300.0;
      if (value.contains("ExternalTemperature"))
      {
        T = value["ExternalTemperature"].get<double>();
      }

      std::optional<double> P{};
      if (value.contains("ExternalPressure"))
      {
        P = value["ExternalPressure"].get<double>();
      }
      if (value.contains("ChemicalPotential"))
      {
        P = std::exp(value["ChemicalPotential"].get<double>() / (Units::KB * T));
      }

      if (caseInSensStringCompare(typeString, "Framework"))
      {
        // Parse framework options
        if (!value.contains("Name"))
        {
          throw std::runtime_error(
              std::format("[Input reader]: framework must have a key 'Name' with a value of string-type'\n"));
        }
        std::string frameworkNameString = value["Name"].get<std::string>();

        int3 jsonNumberOfUnitCells{1, 1, 1};
        if (value.contains("NumberOfUnitCells"))
        {
          jsonNumberOfUnitCells = parseInt3("NumberOfUnitCells", value["NumberOfUnitCells"]);
        }

        double heliumVoidFraction{1.0};
        if (value.contains("HeliumVoidFraction"))
        {
          heliumVoidFraction = value["HeliumVoidFraction"].get<double>();
        }

        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }

        std::optional<Framework> jsonFrameworkComponents{Framework(0, forceFields[systemId].value(),
                                                                   frameworkNameString, frameworkNameString,
                                                                   jsonNumberOfUnitCells, useChargesFrom)};

        // create system
        systems[systemId] = System(systemId, forceFields[systemId].value(), std::nullopt, T, P, heliumVoidFraction,
                                   jsonFrameworkComponents, jsonComponents[systemId],
                                   jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
      }
      else if (caseInSensStringCompare(typeString, "Box"))
      {
        // Parse box options

        double3 boxLengths{25.0, 25.0, 25.0};
        if (value.contains("BoxLengths"))
        {
          boxLengths = parseDouble3("BoxLengths", value["BoxLengths"]);
        }

        double3 boxAngles{90.0, 90.0, 90.0};
        if (value.contains("BoxAngles"))
        {
          boxAngles = parseDouble3("BoxAngles", value["BoxAngles"]);
        }
        boxAngles = boxAngles * (std::numbers::pi / 180.0);

        // create system
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        SimulationBox simulationBox{boxLengths.x, boxLengths.y, boxLengths.z, boxAngles.x, boxAngles.y, boxAngles.z};
        systems[systemId] =
            System(systemId, forceFields[systemId].value(), simulationBox, T, P, 1.0, {}, jsonComponents[systemId],
                   jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
      }
      else
      {
        throw std::runtime_error(
            std::format("[Input reader]: system key 'Type' must have value 'Box' or 'Framework'\n"));
      }

      if (value.contains("MacroStateUseBias") && value["MacroStateUseBias"].is_boolean())
      {
        systems[systemId].tmmc.useBias = value["MacroStateUseBias"].get<bool>();
      }

      if (value.contains("MacroStateMinimumNumberOfMolecules") &&
          value["MacroStateMinimumNumberOfMolecules"].is_number_unsigned())
      {
        systems[systemId].tmmc.minMacrostate = value["MacroStateMinimumNumberOfMolecules"].get<std::size_t>();
      }

      if (value.contains("MacroStateMaximumNumberOfMolecules") &&
          value["MacroStateMaximumNumberOfMolecules"].is_number_unsigned())
      {
        systems[systemId].tmmc.maxMacrostate = value["MacroStateMaximumNumberOfMolecules"].get<std::size_t>();
      }

      if (value.contains("ExternalField") && value["ExternalField"].is_boolean())
      {
        systems[systemId].hasExternalField = value["ExternalField"].get<bool>();
      }

      if (value.contains("ComputeEnergyHistogram") && value["ComputeEnergyHistogram"].is_boolean())
      {
        if (value["ComputeEnergyHistogram"].get<bool>())
        {
          std::size_t sampleEnergyHistogramEvery{1};
          if (value.contains("SampleEnergyHistogramEvery") && value["SampleEnergyHistogramEvery"].is_number_unsigned())
          {
            sampleEnergyHistogramEvery = value["SampleEnergyHistogramEvery"].get<std::size_t>();
          }

          std::size_t writeEnergyHistogramEvery{5000};
          if (value.contains("WriteEnergyHistogramEvery") && value["WriteEnergyHistogramEvery"].is_number_unsigned())
          {
            writeEnergyHistogramEvery = value["WriteEnergyHistogramEvery"].get<std::size_t>();
          }

          std::size_t numberOfBinsEnergyHistogram{128};
          if (value.contains("NumberOfBinsEnergyHistogram") &&
              value["NumberOfBinsEnergyHistogram"].is_number_unsigned())
          {
            numberOfBinsEnergyHistogram = value["NumberOfBinsEnergyHistogram"].get<std::size_t>();
          }

          double lowerLimitEnergyHistogram{-5000.0};
          if (value.contains("LowerLimitEnergyHistogram") && value["LowerLimitEnergyHistogram"].is_number_float())
          {
            lowerLimitEnergyHistogram = value["LowerLimitEnergyHistogram"].get<double>();
          }

          double upperLimitEnergyHistogram{1000.0};
          if (value.contains("UpperLimitEnergyHistogram") && value["UpperLimitEnergyHistogram"].is_number_float())
          {
            upperLimitEnergyHistogram = value["UpperLimitEnergyHistogram"].get<double>();
          }

          systems[systemId].averageEnergyHistogram = PropertyEnergyHistogram(
              jsonNumberOfBlocks, numberOfBinsEnergyHistogram, {lowerLimitEnergyHistogram, upperLimitEnergyHistogram},
              sampleEnergyHistogramEvery, writeEnergyHistogramEvery);
        }
      }

      if (value.contains("ComputeNumberOfMoleculesHistogram") &&
          value["ComputeNumberOfMoleculesHistogram"].is_boolean())
      {
        if (value["ComputeNumberOfMoleculesHistogram"].get<bool>())
        {
          std::size_t sampleNumberOfMoleculesHistogramEvery{1};
          if (value.contains("SampleNumberOfMoleculesHistogramEvery") &&
              value["SampleNumberOfMoleculesHistogramEvery"].is_number_unsigned())
          {
            sampleNumberOfMoleculesHistogramEvery = value["SampleNumberOfMoleculesHistogramEvery"].get<std::size_t>();
          }

          std::size_t writeNumberOfMoleculesHistogramEvery{5000};
          if (value.contains("WriteNumberOfMoleculesHistogramEvery") &&
              value["WriteNumberOfMoleculesHistogramEvery"].is_number_unsigned())
          {
            writeNumberOfMoleculesHistogramEvery = value["WriteNumberOfMoleculesHistogramEvery"].get<std::size_t>();
          }

          std::size_t minimumRangeNumberOfMoleculesHistogram{0};
          if (value.contains("LowerLimitNumberOfMoleculesHistogram") &&
              value["LowerLimitNumberOfMoleculesHistogram"].is_number_unsigned())
          {
            minimumRangeNumberOfMoleculesHistogram = value["LowerLimitNumberOfMoleculesHistogram"].get<std::size_t>();
          }

          std::size_t maximumRangeNumberOfMoleculesHistogram{200};
          if (value.contains("UpperLimitNumberOfMoleculesHistogram") &&
              value["UpperLimitNumberOfMoleculesHistogram"].is_number_unsigned())
          {
            maximumRangeNumberOfMoleculesHistogram = value["UpperLimitNumberOfMoleculesHistogram"].get<std::size_t>();
          }

          systems[systemId].averageNumberOfMoleculesHistogram = PropertyNumberOfMoleculesHistogram(
              jsonNumberOfBlocks, {minimumRangeNumberOfMoleculesHistogram, maximumRangeNumberOfMoleculesHistogram},
              systems[systemId].components.size(), sampleNumberOfMoleculesHistogramEvery,
              writeNumberOfMoleculesHistogramEvery);
        }
      }

      if (value.contains("ComputeRDF") && value["ComputeRDF"].is_boolean())
      {
        if (value["ComputeRDF"].get<bool>())
        {
          std::size_t sampleRDFEvery{10};
          if (value.contains("SampleRDFEvery") && value["SampleRDFEvery"].is_number_unsigned())
          {
            sampleRDFEvery = value["SampleRDFEvery"].get<std::size_t>();
          }

          std::size_t writeRDFEvery{5000};
          if (value.contains("WriteRDFEvery") && value["WriteRDFEvery"].is_number_unsigned())
          {
            writeRDFEvery = value["WriteRDFEvery"].get<std::size_t>();
          }

          std::size_t numberOfBinsRDF{128};
          if (value.contains("NumberOfBinsRDF") && value["NumberOfBinsRDF"].is_number_unsigned())
          {
            numberOfBinsRDF = value["NumberOfBinsRDF"].get<std::size_t>();
          }

          double rangeRDF{15.0};
          if (value.contains("UpperLimitRDF") && value["UpperLimitRDF"].is_number_float())
          {
            rangeRDF = value["UpperLimitRDF"].get<double>();
          }

          systems[systemId].propertyRadialDistributionFunction =
              PropertyRadialDistributionFunction(jsonNumberOfBlocks, systems[systemId].forceField.pseudoAtoms.size(),
                                                 numberOfBinsRDF, rangeRDF, sampleRDFEvery, writeRDFEvery);
        }
      }

      if (value.contains("ComputeConventionalRDF") && value["ComputeConventionalRDF"].is_boolean())
      {
        if (value["ComputeConventionalRDF"].get<bool>())
        {
          std::size_t sampleConventionalRDFEvery{10};
          if (value.contains("SampleConventionalRDFEvery") && value["SampleConventionalRDFEvery"].is_number_unsigned())
          {
            sampleConventionalRDFEvery = value["SampleConventionalRDFEvery"].get<std::size_t>();
          }

          std::size_t writeConventionalRDFEvery{5000};
          if (value.contains("WriteConventionalRDFEvery") && value["WriteConventionalRDFEvery"].is_number_unsigned())
          {
            writeConventionalRDFEvery = value["WriteConventionalRDFEvery"].get<std::size_t>();
          }

          std::size_t numberOfBinsConventionalRDF{128};
          if (value.contains("NumberOfBinsConventionalRDF") &&
              value["NumberOfBinsConventionalRDF"].is_number_unsigned())
          {
            numberOfBinsConventionalRDF = value["NumberOfBinsConventionalRDF"].get<std::size_t>();
          }

          double rangeConventionalRDF{15.0};
          if (value.contains("RangeConventionalRDF") && value["RangeConventionalRDF"].is_number_float())
          {
            rangeConventionalRDF = value["RangeConventionalRDF"].get<double>();
          }

          systems[systemId].propertyConventionalRadialDistributionFunction =
              PropertyConventionalRadialDistributionFunction(
                  jsonNumberOfBlocks, systems[systemId].forceField.pseudoAtoms.size(), numberOfBinsConventionalRDF,
                  rangeConventionalRDF, sampleConventionalRDFEvery, writeConventionalRDFEvery);
        }
      }

      if (value.contains("ComputeMSD") && value["ComputeMSD"].is_boolean())
      {
        if (value["ComputeMSD"].get<bool>())
        {
          std::size_t sampleMSDEvery{10};
          if (value.contains("SampleMSDEvery") && value["SampleMSDEvery"].is_number_unsigned())
          {
            sampleMSDEvery = value["SampleMSDEvery"].get<std::size_t>();
          }

          std::size_t writeMSDEvery{5000};
          if (value.contains("WriteMSDEvery") && value["WriteMSDEvery"].is_number_unsigned())
          {
            writeMSDEvery = value["WriteMSDEvery"].get<std::size_t>();
          }

          std::size_t numberOfBlockElementsMSD{25};
          if (value.contains("NumberOfBlockElementsMSD") && value["NumberOfBlockElementsMSD"].is_number_unsigned())
          {
            numberOfBlockElementsMSD = value["NumberOfBlockElementsMSD"].get<std::size_t>();
          }

          systems[systemId].propertyMSD = PropertyMeanSquaredDisplacement(
              systems[systemId].components.size(), systems[systemId].moleculePositions.size(), sampleMSDEvery,
              writeMSDEvery, numberOfBlockElementsMSD);
        }
      }

      if (value.contains("ComputeVACF") && value["ComputeVACF"].is_boolean())
      {
        if (value["ComputeVACF"].get<bool>())
        {
          std::size_t sampleVACFEvery{10};
          if (value.contains("SampleVACFEvery") && value["SampleVACFEvery"].is_number_unsigned())
          {
            sampleVACFEvery = value["SampleVACFEvery"].get<std::size_t>();
          }

          std::size_t writeVACFEvery{5000};
          if (value.contains("WriteVACFEvery") && value["WriteVACFEvery"].is_number_unsigned())
          {
            writeVACFEvery = value["WriteVACFEvery"].get<std::size_t>();
          }

          std::size_t numberOfBuffersVACF{20};
          if (value.contains("NumberOfBuffersVACF") && value["NumberOfBuffersVACF"].is_number_unsigned())
          {
            numberOfBuffersVACF = value["NumberOfBuffersVACF"].get<std::size_t>();
          }

          std::size_t bufferLengthVACF{1000};
          if (value.contains("BufferLengthVACF") && value["BufferLengthVACF"].is_number_unsigned())
          {
            bufferLengthVACF = value["BufferLengthVACF"].get<std::size_t>();
          }

          systems[systemId].propertyVACF = PropertyVelocityAutoCorrelationFunction(
              systems[systemId].components.size(), systems[systemId].moleculePositions.size(), numberOfBuffersVACF,
              bufferLengthVACF, sampleVACFEvery, writeVACFEvery);
        }
      }

      if (value.contains("ComputeDensityGrid") && value["ComputeDensityGrid"].is_boolean())
      {
        if (value["ComputeDensityGrid"].get<bool>())
        {
          std::size_t sampleDensityGridEvery{10};
          if (value.contains("SampleDensityGridEvery") && value["SampleDensityGridEvery"].is_number_unsigned())
          {
            sampleDensityGridEvery = value["SampleDensityGridEvery"].get<std::size_t>();
          }

          std::size_t writeDensityGridEvery{5000};
          if (value.contains("WriteDensityGridEvery") && value["WriteDensityGridEvery"].is_number_unsigned())
          {
            writeDensityGridEvery = value["WriteDensityGridEvery"].get<std::size_t>();
          }

          int3 densityGridSize{128, 128, 128};
          if (value.contains("DensityGridSize") && value["DensityGridSize"].is_array())
          {
            densityGridSize = parseInt3("DensityGridSize", value["DensityGridSize"]);
          }

          std::vector<std::size_t> densityGridPseudoAtomsList{};
          if (value.contains("DensityGridPseudoAtomsList") && value["DensityGridPseudoAtomsList"].is_array())
          {
            std::vector<std::string> string_list = value["DensityGridPseudoAtomsList"].get<std::vector<std::string>>();
            for (std::string string : string_list)
            {
              std::optional<std::size_t> atomType = systems[systemId].forceField.findPseudoAtom(string);
              if (atomType.has_value())
              {
                densityGridPseudoAtomsList.push_back(atomType.value());
              }
            }
          }

          PropertyDensityGrid::Normalization norm = PropertyDensityGrid::Normalization::Max;
          if (value.contains("DensityGridNormalization") && value["DensityGridNormalization"].is_string())
          {
            std::string normString = value["DensityGridNormalization"].get<std::string>();
            if (caseInSensStringCompare(normString, "NumberDensity"))
            {
              norm = PropertyDensityGrid::Normalization::NumberDensity;
            }
          }

          systems[systemId].propertyDensityGrid = PropertyDensityGrid(
              systems[systemId].framework ? 1 : 0, systems[systemId].components.size(), densityGridSize,
              sampleDensityGridEvery, writeDensityGridEvery, densityGridPseudoAtomsList, norm);
        }
      }

      if (value.contains("OutputPDBMovie") && value["OutputPDBMovie"].is_boolean())
      {
        if (value["OutputPDBMovie"].get<bool>())
        {
          std::size_t sampleMovieEvery{1};
          if (value.contains("SampleMovieEvery") && value["SampleMovieEvery"].is_number_unsigned())
          {
            sampleMovieEvery = value["SampleMovieEvery"].get<std::size_t>();
          }

          systems[systemId].samplePDBMovie = SampleMovie(systemId, sampleMovieEvery);
        }
      }

      if (value.contains("WriteLammpsData") && value["WriteLammpsData"].is_boolean())
      {
        if (value["WriteLammpsData"].get<bool>())
        {
          std::size_t writeLammpsDataEvery = 10;
          if (value.contains("WriteLammpsDataEvery") && value["WriteLammpsDataEvery"].is_number_unsigned())
          {
            writeLammpsDataEvery = value["WriteLammpsDataEvery"].get<std::size_t>();
          }
          systems[systemId].writeLammpsData = WriteLammpsData(systemId, writeLammpsDataEvery);
        }
      }

      if (value.contains("Ensemble") && value["Ensemble"].is_string())
      {
        std::size_t thermostatChainLength{5};
        std::size_t numberOfYoshidaSuzukiSteps{5};

        std::string ensembleString = value["Ensemble"].get<std::string>();
        if (caseInSensStringCompare(ensembleString, "NVT"))
        {
          systems[systemId].thermostat =
              Thermostat(systems[systemId].temperature, thermostatChainLength, numberOfYoshidaSuzukiSteps,
                         systems[systemId].timeStep, systems[systemId].translationalDegreesOfFreedom,
                         systems[systemId].rotationalDegreesOfFreedom);
        }
      }

      if (value.contains("TimeStep") && value["TimeStep"].is_number_float())
      {
        systems[systemId].timeStep = value["TimeStep"].get<double>();
      }
      if (value.contains("HybridMCMoveNumberOfSteps") && value["HybridMCMoveNumberOfSteps"].is_number_unsigned())
      {
        systems[systemId].numberOfHybridMCSteps = value["HybridMCMoveNumberOfSteps"].get<std::size_t>();
        if (value.contains("TimeStep") && value["TimeStep"].is_number_float())
        {
          systems[systemId].mc_moves_statistics.setMaxChange(MoveTypes::HybridMC, value["Timestep"].get<double>());
        }
      }

      systemId++;
    }
  }

  // Post-compute
  // ========================================================

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    systems[i].maxIsothermTerms = 0uz;
    if (!systems[i].components.empty())
    {
      std::vector<Component>::iterator maxIsothermTermsIterator = std::max_element(
          systems[i].components.begin(), systems[i].components.end(),
          [](Component& lhs, Component& rhs) { return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites; });
      systems[i].maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
    }
  }

  if (simulationType == SimulationType::MonteCarloTransitionMatrix)
  {
    for (std::size_t i = 0uz; i < systems.size(); ++i)
    {
      systems[i].tmmc.doTMMC = true;
      systems[i].tmmc.useBias = true;
      systems[i].tmmc.useTMBias = true;
    }
  }

  // Checks
  // ========================================================

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    for (std::size_t reactionId = 0uz; const Reaction& reaction : systems[i].reactions.list)
    {
      if (reaction.productStoichiometry.size() != systems[i].numerOfAdsorbateComponents() ||
          (reaction.productStoichiometry.size() != systems[i].numerOfAdsorbateComponents()))
      {
        throw std::runtime_error(
            std::format("Error [Reaction {}]: mismatch Stoichiometry ({} given not equal"
                        "to twice the number of components {})\n",
                        reactionId, reaction.productStoichiometry.size() + reaction.reactantStoichiometry.size(),
                        2uz * systems[i].numerOfAdsorbateComponents()));
      }

      ++reactionId;
    }
  }

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    std::size_t numberOfDUDlambda{0uz};
    for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].lambdaGC.computeDUdlambda)
      {
        ++numberOfDUDlambda;
      }
    }
    if (numberOfDUDlambda > 1)
    {
      throw std::runtime_error(
          std::format("Error [System {}]: multiple thermodynamic integrations present "
                      "(there can be only one)\n",
                      i));
    }
  }

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    double sum = 0.0;
    for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].type != Component::Type::Cation)
      {
        sum += systems[i].components[j].molFraction;
      }
    }
    if (std::abs(sum - 1.0) > 1e-15)
    {
      for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if (systems[i].components[j].type != Component::Type::Cation)
        {
          systems[i].components[j].molFraction /= sum;
        }
      }
    }
  }

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    systems[i].numberOfCarrierGases = 0uz;
    systems[i].carrierGasComponent = 0uz;
    for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].isCarrierGas)
      {
        systems[i].carrierGasComponent = j;
        std::vector<double> values{1.0, 0.0};
        const Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
        systems[i].components[systems[i].carrierGasComponent].isotherm.add(isotherm);
        systems[i].components[systems[i].carrierGasComponent].isotherm.numberOfSites = 1;

        systems[i].numberOfCarrierGases++;
      }
    }

    if (simulationType == SimulationType::Breakthrough)
    {
      if (systems[i].numberOfCarrierGases == 0uz)
      {
        throw std::runtime_error("Error [Breakthrough]: no carrier gas component present\n");
      }
      if (systems[i].numberOfCarrierGases > 1)
      {
        throw std::runtime_error(
            "Error [Breakthrough]: multiple carrier gas component present (there can be only one)\n");
      }
    }
  }

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    if (systems[i].tmmc.doTMMC)
    {
      if (systems[i].numerOfAdsorbateComponents() > 1)
      {
        throw std::runtime_error("Error: Multiple components for TMMC not yet implemented.\n");
      }

      // check initial number of molecules is in the range of the TMMC macrostates
      for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if (systems[i].components[j].type == Component::Type::Adsorbate)
        {
          std::size_t numberOfMolecules = systems[i].initialNumberOfMolecules[j];
          if (numberOfMolecules < systems[i].tmmc.minMacrostate || numberOfMolecules > systems[i].tmmc.maxMacrostate)
          {
            throw std::runtime_error(
                std::format("Error: Molecules created ({}) need to fit into the TMMC macrostate "
                            "range ({}-{})\n",
                            numberOfMolecules, systems[i].tmmc.minMacrostate, systems[i].tmmc.maxMacrostate));
          }
        }
      }
    }
  }
}

const std::set<std::string, InputReader::InsensitiveCompare> InputReader::generalOptions = {
    "SimulationType",
    "Units",
    "ForceField",
    "NumberOfLambdaBins",
    "RestartFromBinaryFile",
    "RandomSeed",
    "NumberOfCycles",
    "NumberOfInitializationCycles",
    "NumberOfEquilibrationCycles",
    "PrintEvery",
    "WriteBinaryRestartEvery",
    "RescaleWangLandauEvery",
    "OptimizeMCMovesEvery",
    "ThreadingType",
    "NumberOfThreads",
    "Components",
    "Systems"};

const std::set<std::string, InputReader::InsensitiveCompare> InputReader::systemOptions = {
    "ForceField",
    "CutOff",
    "CutOffVDW",
    "CutOffCoulomb",
    "OmitEwaldFourier",
    "ComputePolarization",
    "ChargeMethod",
    "VolumeMoveProbability",
    "GibbsVolumeMoveProbability",
    "ParallelTemperingSwapProbability",
    "HybridMCProbability",
    "HybridMCMoveNumberOfSteps",
    "Type",
    "ExternalTemperature",
    "ExternalPressure",
    "ChemicalPotential",
    "UseChargesFrom",
    "Framework",
    "Name",
    "NumberOfUnitCells",
    "HeliumVoidFraction",
    "BoxLengths",
    "BoxAngles",
    "ExternalField",
    "ComputeEnergyHistogram",
    "SampleEnergyHistogramEvery",
    "WriteEnergyHistogramEvery",
    "NumberOfBinsEnergyHistogram",
    "LowerLimitEnergyHistogram",
    "UpperLimitEnergyHistogram",
    "ComputeNumberOfMoleculesHistogram",
    "SampleNumberOfMoleculesHistogramEvery",
    "WriteNumberOfMoleculesHistogramEvery",
    "LowerLimitNumberOfMoleculesHistogram",
    "UpperLimitNumberOfMoleculesHistogram",
    "ComputeRDF",
    "SampleRDFEvery",
    "WriteRDFEvery",
    "NumberOfBinsRDF",
    "UpperLimitRDF",
    "ComputeConventionalRDF",
    "SampleConventionalRDFEvery",
    "WriteConventionalRDFEvery",
    "NumberOfBinsConventionalRDF",
    "RangeConventionalRDF",
    "ComputeMSD",
    "SampleMSDEvery",
    "WriteMSDEvery",
    "NumberOfBlockElementsMSD",
    "ComputeVACF",
    "SampleVACFEvery",
    "WriteVACFEvery",
    "NumberOfBuffersVACF",
    "BufferLengthVACF",
    "ComputeDensityGrid",
    "SampleDensityGridEvery",
    "WriteDensityGridEvery",
    "DensityGridSize",
    "DensityGridPseudoAtomsList",
    "DensityGridNormalization",
    "OutputPDBMovie",
    "SampleMovieEvery",
    "WriteLammpsData",
    "WriteLammpsDataEvery",
    "Ensemble",
    "TimeStep",
    "MacroStateUseBias",
    "MacroStateMinimumNumberOfMolecules",
    "MacroStateMaximumNumberOfMolecules"};

const std::set<std::string, InputReader::InsensitiveCompare> InputReader::componentOptions = {
    "Name",
    "Type",
    "MoleculeDefinition",
    "TranslationProbability",
    "RandomTranslationProbability",
    "RotationProbability",
    "RandomRotationProbability",
    "ReinsertionProbability",
    "SwapConventionalProbability",
    "SwapProbability",
    "CFCMC_SwapProbability",
    "CFCMC_CBMC_SwapProbability",
    "GibbsSwapProbability",
    "Gibbs_CFCMC_SwapProbability",
    "WidomProbability",
    "CFCMC_WidomProbability",
    "CFCMC_CBMC_WidomProbability",
    "CreateNumberOfMolecules",
    "StartingBead",
    "FugacityCoefficient",
    "IdealGasRosenbluthWeight",
    "MolFraction",
    "ThermodynamicIntegration",
    "LambdaBiasFileName",
    "BlockingPockets"};

void InputReader::validateInput(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  for (auto& [key, _] : parsed_data.items())
  {
    if (!generalOptions.contains(key))
    {
      throw std::runtime_error(std::format("Error: Unknown input '{}'\n", key));
    }
  }

  if (parsed_data.contains("Systems"))
  {
    for (auto& [_, value] : parsed_data["Systems"].items())
    {
      for (auto& [key, _] : value.items())
      {
        if (!systemOptions.contains(key))
        {
          throw std::runtime_error(std::format("Error: Unknown system input '{}'\n", key));
        }
      }
    }
  }

  if (parsed_data.contains("Components"))
  {
    for (auto& [_, value] : parsed_data["Components"].items())
    {
      for (auto& [key, _] : value.items())
      {
        if (!componentOptions.contains(key))
        {
          throw std::runtime_error(std::format("Error: Unknown component input '{}'\n", key));
        }
      }
    }
  }
}

void InputReader::parseUnits(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  if (parsed_data.contains("Units") && parsed_data["Units"].is_string())
  {
    std::string unitsString = parsed_data["Units"].get<std::string>();

    if (caseInSensStringCompare(unitsString, "Reduced"))
    {
      Units::setUnits(Units::System::ReducedUnits);
    }
  }
}
