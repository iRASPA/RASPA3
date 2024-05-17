module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
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
#include <sstream>
#include <streambuf>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module input_reader;

#ifndef USE_LEGACY_HEADERS
import <filesystem>;
import <fstream>;
import <streambuf>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <exception>;
import <numbers>;
import <vector>;
import <array>;
import <complex>;
import <ios>;
import <optional>;
import <algorithm>;
import <map>;
import <iterator>;
import <functional>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
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
import units;
import sample_movies;
import threadpool;
import threading;
import isotherm;
import multi_site_isotherm;
import pressure_range;
import mc_moves_probabilities_system;
import mc_moves_probabilities_particles;
import reaction;
import reactions;
import transition_matrix;
import property_conventional_rdf;
import property_rdf;
import property_density_grid;

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
      value.x = json[0].template get<int32_t>();
      value.y = json[1].template get<int32_t>();
      value.z = json[2].template get<int32_t>();
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
std::vector<T> parseList(size_t size, const std::string& item, auto json)
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

InputReader::InputReader(const std::string inputFile) : inputStream(inputFile)
{
  std::ifstream input("simulation.json");

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  size_t jsonNumberOfBlocks{5};

  try
  {
    parsed_data = nlohmann::json::parse(input);
  }
  catch (nlohmann::json::parse_error& ex)
  {
    std::cerr << "parse error at byte " << ex.byte << std::endl;
  }

  // count number of systems
  if (!parsed_data.contains("Systems"))
  {
    throw std::runtime_error(
        std::format("[Input reader]: no system defined with keyword 'Systems' and value of array-type\n"));
  }
  size_t jsonNumberOfSystems = parsed_data["Systems"].size();
  if (jsonNumberOfSystems == 0)
  {
    throw std::runtime_error(std::format("[Input reader]: keyword 'Systems' has empty value of array-type\n"));
  }

  systems = std::vector<System>(jsonNumberOfSystems);

  // Read the local 'force_field.json' if present. This file will be used if no 'ForceField' keyword is specified per
  // system
  std::optional<std::string> directoryName{};
  if (parsed_data["ForceField"].is_string())
  {
    directoryName = parsed_data["ForceField"].get<std::string>();
  }
  std::vector<std::optional<ForceField>> forceFields = std::vector<std::optional<ForceField>>(jsonNumberOfSystems);
  const std::optional<ForceField> standard = ForceField::readForceField(directoryName, "force_field.json");
  for (size_t i = 0; i != jsonNumberOfSystems; ++i)
  {
    forceFields[i] = standard;
  }

  if (parsed_data["NumberOfCycles"].is_number_unsigned())
  {
    numberOfCycles = parsed_data["NumberOfCycles"].get<size_t>();
  }

  if (parsed_data["NumberOfInitializationCycles"].is_number_unsigned())
  {
    numberOfInitializationCycles = parsed_data["NumberOfInitializationCycles"].get<size_t>();
  }

  if (parsed_data["NumberOfEquilibrationCycles"].is_number_unsigned())
  {
    numberOfEquilibrationCycles = parsed_data["NumberOfEquilibrationCycles"].get<size_t>();
  }

  if (parsed_data["PrintEvery"].is_number_unsigned())
  {
    printEvery = parsed_data["PrintEvery"].get<size_t>();
  }

  if (parsed_data["WriteBinaryRestartEvery"].is_number_unsigned())
  {
    writeBinaryRestartEvery = parsed_data["WriteBinaryRestartEvery"].get<size_t>();
  }

  if (parsed_data["RescaleWangLandauEvery"].is_number_unsigned())
  {
    rescaleWangLandauEvery = parsed_data["RescaleWangLandauEvery"].get<size_t>();
  }

  if (parsed_data["OptimizeMCMovesEvery"].is_number_unsigned())
  {
    optimizeMCMovesEvery = parsed_data["OptimizeMCMovesEvery"].get<size_t>();
  }


  if (parsed_data["SimulationType"].is_string())
  {
    std::string simulationTypeString = parsed_data["SimulationType"].get<std::string>();
    if (caseInSensStringCompare(simulationTypeString, "MonteCarlo"))
    {
      simulationType = SimulationType::MonteCarlo;
    }
    else if (caseInSensStringCompare(simulationTypeString, "MonteCarloTransitionMatrix"))
    {
      simulationType = SimulationType::MonteCarloTransitionMatrix;
    }
    else if (caseInSensStringCompare(simulationTypeString, "MolecularDynamics"))
    {
      simulationType = SimulationType::MolecularDynamics;
    }
    else if (caseInSensStringCompare(simulationTypeString, "Minimization"))
    {
      simulationType = SimulationType::Minimization;
    }
    else if (caseInSensStringCompare(simulationTypeString, "Test"))
    {
      simulationType = SimulationType::Test;
    }
    else if (caseInSensStringCompare(simulationTypeString, "Breakthrough"))
    {
      simulationType = SimulationType::Breakthrough;
    }
    else if (caseInSensStringCompare(simulationTypeString, "MixturePrediction"))
    {
      simulationType = SimulationType::MixturePrediction;
    }
    else if (caseInSensStringCompare(simulationTypeString, "Fitting"))
    {
      simulationType = SimulationType::Fitting;
    }
    else if (caseInSensStringCompare(simulationTypeString, "ParallelTempering"))
    {
      simulationType = SimulationType::ParallelTempering;
    }
    else
    {
      throw std::runtime_error(
          std::format("[Input reader]: {} not a valid simulation type", simulationTypeString, parsed_data.dump()));
    }
  }

  if (parsed_data["ThreadingType"].is_string())
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

  if (parsed_data["NumberOfThreads"].is_number_unsigned())
  {
    numberOfThreads = parsed_data["NumberOfThreads"].get<size_t>();
    if (numberOfThreads > 1) threadingType = ThreadPool::ThreadingType::ThreadPool;
    else threadingType = ThreadPool::ThreadingType::Serial;
  }

  // count number of components
  size_t jsonNumberOfComponents{};
  if (parsed_data.contains("Components"))
  {
    jsonNumberOfComponents = parsed_data["Components"].size();
  }

  // pre-allocate the vector of vector (jsonNumberOfSystems x jsonNumberOfComponents), i.e. for each system a list of
  // components
  std::vector<std::vector<Component>> jsonComponents(jsonNumberOfSystems,
                                                     std::vector<Component>(jsonNumberOfComponents));

  // for each system, a list of how many molecules to create for each component
  std::vector<std::vector<size_t>> jsonCreateNumberOfMolecules(jsonNumberOfSystems,
                                                               std::vector<size_t>(jsonNumberOfComponents));

  // Parse component options
  for (size_t componentId = 0; auto& [_, item] : parsed_data["Components"].items())
  {
    std::vector<MCMoveProbabilitiesParticles> move_probabilities(jsonNumberOfSystems);

    if (!item.contains("Name"))
    {
      throw std::runtime_error(
          std::format("[Input reader]: component must have a key 'Name' with a value of string-type'\n"));
    }
    std::string jsonComponentName = item["Name"].get<std::string>();

    // Convenience notation listing the properties as a single value. These will then be taken for all systems.
    // ========================================================================================================

    if (item["TranslationProbability"].is_number_float())
    {
      double probabilityTranslationMove = item["TranslationProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityTranslationMove = probabilityTranslationMove;
      }
    }

    if (item["RandomTranslationProbability"].is_number_float())
    {
      double probabilityRandomTranslationMove = item["RandomTranslationProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityRandomTranslationMove = probabilityRandomTranslationMove;
      }
    }

    if (item["RotationProbability"].is_number_float())
    {
      double probabilityRotationMove = item["RotationProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityRotationMove = probabilityRotationMove;
      }
    }

    if (item["RandomRotationProbability"].is_number_float())
    {
      double probabilityRandomTranslationMove = item["RandomRotationProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityRandomTranslationMove = probabilityRandomTranslationMove;
      }
    }

    if (item["ReinsertionProbability"].is_number_float())
    {
      double probabilityReinsertionMove_CBMC = item["ReinsertionProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityReinsertionMove_CBMC = probabilityReinsertionMove_CBMC;
      }
    }

    if (item["SwapConventionalProbability"].is_number_float())
    {
      double probabilitySwapMove = item["SwapConventionalProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilitySwapMove = probabilitySwapMove;
      }
    }

    if (item["SwapProbability"].is_number_float())
    {
      double probabilitySwapMove_CBMC = item["SwapProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilitySwapMove_CBMC = probabilitySwapMove_CBMC;
      }
    }

    if (item["CFCMC_SwapProbability"].is_number_float())
    {
      double probabilitySwapMove_CFCMC = item["CFCMC_SwapProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilitySwapMove_CFCMC = probabilitySwapMove_CFCMC;
      }
    }

    if (item["CFCMC_CBMC_SwapProbability"].is_number_float())
    {
      double probabilitySwapMove_CFCMC_CBMC = item["CFCMC_CBMC_SwapProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilitySwapMove_CFCMC_CBMC = probabilitySwapMove_CFCMC_CBMC;
      }
    }

    if (item["GibbsCFCMCSwapProbability"].is_number_float())
    {
      double probabilityGibbsSwapMove_CFCMC = item["GibbsCFCMCSwapProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityGibbsSwapMove_CFCMC = probabilityGibbsSwapMove_CFCMC;
      }
    }

    if (item["GibbsSwapProbability"].is_number_float())
    {
      double probabilityGibbsSwapMove_CBMC = item["GibbsSwapProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityGibbsSwapMove_CBMC = probabilityGibbsSwapMove_CBMC;
      }
    }

    if (item["WidomProbability"].is_number_float())
    {
      double probabilityWidomMove = item["WidomProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityWidomMove = probabilityWidomMove;
      }
    }

    if (item["CFCMC_WidomProbability"].is_number_float())
    {
      double probabilityWidomMove_CFCMC = item["CFCMC_WidomProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityWidomMove_CFCMC = probabilityWidomMove_CFCMC;
      }
    }

    if (item["CFCMC_CBMC_WidomProbability"].is_number_float())
    {
      double probabilityWidomMove_CFCMC_CBMC = item["CFCMC_CBMC_WidomProbability"].get<double>();
      for (size_t i = 0; i < move_probabilities.size(); ++i)
      {
        move_probabilities[i].probabilityWidomMove_CFCMC_CBMC = probabilityWidomMove_CFCMC_CBMC;
      }
    }


    if (item["CreateNumberOfMolecules"].is_number_integer())
    {
      size_t n = item["CreateNumberOfMolecules"].get<size_t>();
      for (size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        jsonCreateNumberOfMolecules[i][componentId] = n;
      }
    }

     size_t jsonNumberOfLambdaBins{31};
     if (parsed_data["NumberOfLambdaBins"].is_number_unsigned())
     {
       jsonNumberOfLambdaBins = parsed_data["NumberOfLambdaBins"].get<size_t>();
     }

    // Explicit notation listing the properties as an array of the values for the particular systems
    // ========================================================================================================

    if (item["CreateNumberOfMolecules"].is_array())
    {
      std::vector<size_t> initialNumberOfMolecule =
          parseList<size_t>(jsonNumberOfSystems, "CreateNumberOfMolecules", item["CreateNumberOfMolecules"]);
      for (size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        jsonCreateNumberOfMolecules[i][componentId] = initialNumberOfMolecule[i];
      }
    }

    // construct Component
    for (size_t i = 0; i != jsonNumberOfSystems; ++i)
    {
      if (!forceFields[i].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }

      jsonComponents[i][componentId] =
          Component(Component::Type::Adsorbate, componentId, forceFields[i].value(), jsonComponentName,
                    jsonComponentName, jsonNumberOfBlocks, jsonNumberOfLambdaBins, move_probabilities[i]);
    }

    if (item["StartingBead"].is_number_integer())
    {
      size_t n = item["StartingBead"].get<size_t>();
      for (size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        if (n >= jsonComponents[i][componentId].definedAtoms.size())
        {
          throw std::runtime_error(std::format("[Input reader]: starting bead larger than the molecule size'\n"));
        }

        jsonComponents[i][componentId].startingBead = n;
      }
    }

    if (item["FugacityCoefficient"].is_number_float())
    {
      double fugacity_coefficient = item["FugacityCoefficient"].get<double>();
      for (size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        jsonComponents[i][componentId].fugacityCoefficient = fugacity_coefficient;
      }
    }


    if (item["ThermodynamicIntegration"].is_boolean())
    {
      bool thermodynamic_integration = item["ThermodynamicIntegration"].get<bool>();
      for (size_t i = 0; i != jsonNumberOfSystems; ++i)
      {
        jsonComponents[i][componentId].lambdaGC.computeDUdlambda = thermodynamic_integration;
      }
    }

    componentId++;
  }

  for (size_t systemId = 0; auto& [key, value] : parsed_data["Systems"].items())
  {
    MCMoveProbabilitiesSystem mc_moves_probabilities{};

    if (value["ForceField"].is_string())
    {
      std::string name = parsed_data["ForceField"].get<std::string>();
      forceFields[systemId] = ForceField::readForceField(name, "force_field.json");
    }

    if (value["CutOffVDW"].is_number_float())
    {
      if (!forceFields[systemId].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }
      forceFields[systemId]->cutOffVDW = value["CutOffVDW"].get<double>();
      forceFields[systemId]->preComputePotentialShift();
    }

    if (value["CutOffCoulomb"].is_number_float())
    {
      if (!forceFields[systemId].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }
      forceFields[systemId]->cutOffCoulomb = value["CutOffCoulomb"].get<double>();
    }

    if (value["OmitEwaldFourier"].is_boolean())
    {
      if (!forceFields[systemId].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }
      forceFields[systemId]->omitEwaldFourier = value["OmitEwaldFourier"].get<bool>();
    }

    if (value["ChargeMethod"].is_string())
    {
      if (!forceFields[systemId].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }

      std::string chargeMethodString = value["ChargeMethod"].get<std::string>();

      if (caseInSensStringCompare(chargeMethodString, "Ewald"))
      {
        forceFields[systemId]->chargeMethod = ForceField::ChargeMethod::Ewald;
        forceFields[systemId]->noCharges = false;
      }
      if (caseInSensStringCompare(chargeMethodString, "None"))
      {
        forceFields[systemId]->chargeMethod = ForceField::ChargeMethod::Ewald;
        forceFields[systemId]->noCharges = true;
      }
    }

    if (value["VolumeMoveProbability"].is_number_float())
    {
      mc_moves_probabilities.probabilityVolumeMove = value["VolumeMoveProbability"].get<double>();
    }

    if (value["GibbsVolumeMoveProbability"].is_number_float())
    {
      mc_moves_probabilities.probabilityGibbsVolumeMove = value["GibbsVolumeMoveProbability"].get<double>();
    }

    if (value["ParallelTemperingSwapProbability"].is_number_float())
    {
      mc_moves_probabilities.probabilityParallelTemperingSwap = value["ParallelTemperingSwapProbability"].get<double>();
    }

    if (!value.contains("Type"))
    {
      throw std::runtime_error(
          std::format("[Input reader]: system must have a key 'Type' with value 'Box' or 'Framework'\n"));
    }
    std::string typeString = value["Type"].get<std::string>();

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

      if (!value.contains("ExternalTemperature"))
      {
        throw std::runtime_error(std::format(
            "[Input reader]: framework must have a key 'ExternalTemperature' with a value of floating-point-type'\n"));
      }
      double T = value["ExternalTemperature"].get<double>();

      std::optional<double> P{};
      if (value.contains("ExternalPressure"))
      {
        P = value["ExternalPressure"].get<double>();
      }

      if (!forceFields[systemId].has_value())
      {
        throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
      }

      std::vector<Framework> jsonFrameworkComponents{
          Framework(0, forceFields[systemId].value(), frameworkNameString, frameworkNameString, jsonNumberOfUnitCells)};

      // create system
      systems[systemId] =
          System(systemId, std::nullopt, T, P, forceFields[systemId].value(), jsonFrameworkComponents,
                 jsonComponents[systemId], jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
    }
    else if (caseInSensStringCompare(typeString, "Box"))
    {
      // Parse box options

      if (!value.contains("ExternalTemperature"))
      {
        throw std::runtime_error(std::format(
            "[Input reader]: framework must have a key 'ExternalTemperature' with a value of floating-point-type'\n"));
      }
      [[maybe_unused]] double T = value["ExternalTemperature"].get<double>();

      std::optional<double> P{};
      if (value.contains("ExternalPressure"))
      {
        P = value["ExternalPressure"].get<double>();
      }
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
          System(systemId, simulationBox, T, P, forceFields[systemId].value(), {}, jsonComponents[systemId],
                 jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
    }
    else
    {
      throw std::runtime_error(std::format("[Input reader]: system key 'Type' must have value 'Box' or 'Framework'\n"));
    }

    if (value["ComputeConventionalRDF"].is_boolean())
    {
      if(value["ComputeConventionalRDF"].get<bool>())
      {
        size_t numberOfBinsConventionalRDF{ 128 };
        if (value["NumberOfBinsConventionalRDF"].is_number_unsigned())
        {
          numberOfBinsConventionalRDF = value["NumberOfBinsConventionalRDF"].get<size_t>();
        }

        double rangeConventionalRDF{ 15.0 };
        if (value["RangeConventionalRDF"].is_number_float())
        {
          rangeConventionalRDF = value["RangeConventionalRDF"].get<double>();
        }

        size_t sampleConventionalRDFEvery{ 10 };
        if (value["SampleConventionalRDFEvery"].is_number_unsigned())
        {
          sampleConventionalRDFEvery = value["SampleConventionalRDFEvery"].get<size_t>();
        }

        size_t writeConventionalRDFEvery{ 5000 };
        if (value["WriteConventionalRDFEvery"].is_number_unsigned())
        {
          writeConventionalRDFEvery = value["WriteConventionalRDFEvery"].get<size_t>();
        }

        systems[systemId].propertyConventionalRadialDistributionFunction = 
          PropertyConventionalRadialDistributionFunction(jsonNumberOfBlocks, systems[systemId].forceField.pseudoAtoms.size(),
                            numberOfBinsConventionalRDF, rangeConventionalRDF, sampleConventionalRDFEvery, writeConventionalRDFEvery);
      }
    }

    if (value["ComputeRDF"].is_boolean())
    {
      if(value["ComputeRDF"].get<bool>())
      {
        size_t numberOfBinsRDF{ 128 };
        if (value["NumberOfBinsRDF"].is_number_unsigned())
        {
          numberOfBinsRDF = value["NumberOfBinsRDF"].get<size_t>();
        }

        double rangeRDF{ 15.0 };
        if (value["RangeRDF"].is_number_float())
        {
          rangeRDF = value["RangeRDF"].get<double>();
        }

        size_t sampleRDFEvery{ 10 };
        if (value["SampleRDFEvery"].is_number_unsigned())
        {
          sampleRDFEvery = value["SampleRDFEvery"].get<size_t>();
        }

        size_t writeRDFEvery{ 5000 };
        if (value["WriteRDFEvery"].is_number_unsigned())
        {
          writeRDFEvery = value["WriteRDFEvery"].get<size_t>();
        }

        systems[systemId].propertyRadialDistributionFunction = 
          PropertyRadialDistributionFunction(jsonNumberOfBlocks, systems[systemId].forceField.pseudoAtoms.size(),
                            numberOfBinsRDF, rangeRDF, sampleRDFEvery, writeRDFEvery);
      }
    }

    if (value["ComputeDensityGrid"].is_boolean())
    {
      if(value["ComputeDensityGrid"].get<bool>())
      {
        int3 densityGridSize{ 128, 128, 128 };
        if (value["DensityGridSize"].is_array())
        {
          densityGridSize = parseInt3("DensityGridSize", value["DensityGridSize"]);
        }

        size_t sampleDensityGridEvery{ 10 };
        if (value["SampleDensityGridEvery"].is_number_unsigned())
        {
          sampleDensityGridEvery = value["SampleDensityGridEvery"].get<size_t>();
        }

        size_t writeDensityGridEvery{ 5000 };
        if (value["WriteDensityGridEvery"].is_number_unsigned())
        {
          writeDensityGridEvery = value["WriteDensityGridEvery"].get<size_t>();
        }

        std::vector<size_t> densityGridPseudoAtomsList{};
        if (value["DensityGridPseudoAtomsList"].is_array())
        {
          std::vector<std::string> string_list = value["DensityGridPseudoAtomsList"].get<std::vector<std::string>>();
          for(std::string string: string_list)
          {
            std::optional<size_t> atomType = systems[systemId].forceField.findPseudoAtom(string);
            if(atomType.has_value())
            {
              densityGridPseudoAtomsList.push_back(atomType.value());
            }
          }
        }


        systems[systemId].propertyDensityGrid = PropertyDensityGrid(systems[systemId].frameworkComponents.size(), systems[systemId].components.size(),
                                      densityGridSize, sampleDensityGridEvery, writeDensityGridEvery, densityGridPseudoAtomsList);
      }
    }

    if (value["OutputPDBMovie"].is_boolean())
    {
      if(value["OutputPDBMovie"].get<bool>())
      {
        size_t sampleMovieEvery{ 1 };
        if (value["SampleMovieEvery"].is_number_unsigned())
        {
          sampleMovieEvery = value["SampleMovieEvery"].get<size_t>();
        }

        systems[systemId].samplePDBMovie = SampleMovie(systemId, sampleMovieEvery);
      }

    }


    systemId++;
  }

  // Post-compute
  // ========================================================

  for (size_t i = 0uz; i < systems.size(); ++i)
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
    for (size_t i = 0uz; i < systems.size(); ++i)
    {
      systems[i].tmmc.doTMMC = true;
    }
  }

  // Checks
  // ========================================================

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    for (size_t reactionId = 0uz; const Reaction& reaction : systems[i].reactions.list)
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

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    size_t numberOfDUDlambda{0uz};
    for (size_t j = 0uz; j < systems[i].components.size(); ++j)
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

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    double sum = 0.0;
    for (size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].type != Component::Type::Framework)
      {
        sum += systems[i].components[j].molFraction;
      }
    }
    if (std::abs(sum - 1.0) > 1e-15)
    {
      for (size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if (systems[i].components[j].type != Component::Type::Framework)
        {
          systems[i].components[j].molFraction /= sum;
        }
      }
    }
  }

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    systems[i].numberOfCarrierGases = 0uz;
    systems[i].carrierGasComponent = 0uz;
    for (size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].type != Component::Type::Framework)
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

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    if (systems[i].tmmc.doTMMC)
    {
      if (systems[i].numerOfAdsorbateComponents() > 1)
      {
        throw std::runtime_error("Error: Multiple components for TMMC not yet implemented.\n");
      }

      // check initial number of molecules is in the range of the TMMC macrostates
      for (size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if (systems[i].components[j].type == Component::Type::Adsorbate)
        {
          size_t numberOfMolecules = systems[i].initialNumberOfMolecules[j];
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
