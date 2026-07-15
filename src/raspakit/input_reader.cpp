module;

module input_reader;

import std;

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
import mc_moves_probabilities;
import mc_moves_move_types;
import reaction;
import reactions;
import partition_function;
import transition_matrix;
import property_conventional_rdf;
import property_rdf;
import property_density_grid;
import property_energy_histogram;
import property_number_of_molecules_histogram;
import property_molecule_properties;
import property_volume_evolution;
import property_number_of_molecules_evolution;
import property_msd;
import property_vacf;
import property_elastic_constants_fluctuation;
import write_lammps_data;
import thermostat;
import thermobarostat;
import cif_reader;
import running_energy;
import minimization_cell_layout;

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

InputReader::InputReader(const std::string inputFile)
{
  if (!std::filesystem::exists(inputFile))
  {
    throw std::runtime_error(std::format("[Input reader]: File '{}' not found\n", inputFile));
  }

  std::ifstream input(inputFile);

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

void InputReader::parseBreakthrough([[maybe_unused]] const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data) {}

void InputReader::parseMolecularSimulations(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  std::size_t jsonNumberOfBlocks{5};
  if (parsed_data.contains("NumberOfBlocks"))
  {
    if (!parsed_data["NumberOfBlocks"].is_number_unsigned())
      throw std::runtime_error("[Input reader]: NumberOfBlocks must be an integer of at least three");
    jsonNumberOfBlocks = parsed_data["NumberOfBlocks"].get<std::size_t>();
    if (jsonNumberOfBlocks < 3)
      throw std::runtime_error("[Input reader]: NumberOfBlocks must be at least three");
  }
  numberOfBlocks = jsonNumberOfBlocks;

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

  if (parsed_data.contains("BinaryRestartFileName") && parsed_data["BinaryRestartFileName"].is_string())
  {
    restartFromBinaryFileName = parsed_data["BinaryRestartFileName"].get<std::string>();
  }

  if (parsed_data.contains("RandomSeed") && parsed_data["RandomSeed"].is_number_unsigned())
  {
    randomSeed = parsed_data["RandomSeed"].get<unsigned long long>();
  }

  if (parsed_data.contains("NumberOfProductionCycles") && parsed_data["NumberOfProductionCycles"].is_number_unsigned())
  {
    numberOfProductionCycles = parsed_data["NumberOfProductionCycles"].get<std::size_t>();
  }

  if (parsed_data.contains("NumberOfPreInitializationCycles") &&
      parsed_data["NumberOfPreInitializationCycles"].is_number_unsigned())
  {
    numberOfPreInitializationCycles = parsed_data["NumberOfPreInitializationCycles"].get<std::size_t>();
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
  minimizationOptions.printEvery = printEvery;

  if (parsed_data.contains("MaximumNumberOfMinimizationSteps"))
  {
    if (!parsed_data["MaximumNumberOfMinimizationSteps"].is_number_unsigned())
    {
      throw std::runtime_error("[Input reader]: MaximumNumberOfMinimizationSteps must be a positive unsigned integer");
    }
    minimizationOptions.maximumNumberOfSteps = parsed_data["MaximumNumberOfMinimizationSteps"].get<std::size_t>();
  }
  auto parseMinimizationNumber = [&](const char* key, double& value)
  {
    if (!parsed_data.contains(key))
    {
      return;
    }
    if (!parsed_data[key].is_number())
    {
      throw std::runtime_error(std::format("[Input reader]: {} must be a finite number", key));
    }
    value = parsed_data[key].get<double>();
    if (!std::isfinite(value))
    {
      throw std::runtime_error(std::format("[Input reader]: {} must be finite", key));
    }
  };
  parseMinimizationNumber("MaximumStepLength", minimizationOptions.maximumStepLength);
  parseMinimizationNumber("MaximumCellStepLength", minimizationOptions.maximumCellStepLength);
  parseMinimizationNumber("RMSGradientTolerance", minimizationOptions.rmsGradientTolerance);
  parseMinimizationNumber("MaxGradientTolerance", minimizationOptions.maxGradientTolerance);
  parseMinimizationNumber("MinimizationConvergenceFactor", minimizationOptions.convergenceFactor);
  parseMinimizationNumber("MinimumEigenvalue", minimizationOptions.minimumEigenvalue);
  parseMinimizationNumber("ElasticEigenvalueTolerance", minimizationOptions.elasticEigenvalueTolerance);

  if (parsed_data.contains("ComputeElasticConstants"))
  {
    if (!parsed_data["ComputeElasticConstants"].is_boolean())
    {
      throw std::runtime_error("[Input reader]: ComputeElasticConstants must be a boolean");
    }
    minimizationOptions.computeElasticConstants = parsed_data["ComputeElasticConstants"].get<bool>();
  }

  if (parsed_data.contains("ComputeNormalModes"))
  {
    if (!parsed_data["ComputeNormalModes"].is_boolean())
    {
      throw std::runtime_error("[Input reader]: ComputeNormalModes must be a boolean");
    }
    minimizationOptions.computeNormalModes = parsed_data["ComputeNormalModes"].get<bool>();
  }

  if (parsed_data.contains("NormalModeMovies"))
  {
    if (!parsed_data["NormalModeMovies"].is_boolean())
    {
      throw std::runtime_error("[Input reader]: NormalModeMovies must be a boolean");
    }
    minimizationOptions.normalModeMovies = parsed_data["NormalModeMovies"].get<bool>();
  }

  if (parsed_data.contains("NormalModeMoviePeriods"))
  {
    if (!parsed_data["NormalModeMoviePeriods"].is_number_unsigned())
    {
      throw std::runtime_error("[Input reader]: NormalModeMoviePeriods must be a positive unsigned integer");
    }
    minimizationOptions.normalModeMoviePeriods = parsed_data["NormalModeMoviePeriods"].get<std::size_t>();
  }

  if (parsed_data.contains("NormalModeMoviePointsPerPeriod"))
  {
    if (!parsed_data["NormalModeMoviePointsPerPeriod"].is_number_unsigned())
    {
      throw std::runtime_error("[Input reader]: NormalModeMoviePointsPerPeriod must be a positive unsigned integer");
    }
    minimizationOptions.normalModeMoviePointsPerPeriod =
        parsed_data["NormalModeMoviePointsPerPeriod"].get<std::size_t>();
  }

  parseMinimizationNumber("NormalModeMovieAmplitude", minimizationOptions.normalModeMovieAmplitude);

  if (minimizationOptions.normalModeMovies &&
      (minimizationOptions.normalModeMoviePeriods == 0 ||
       minimizationOptions.normalModeMoviePointsPerPeriod == 0 || minimizationOptions.normalModeMovieAmplitude <= 0.0))
  {
    throw std::runtime_error(
        "[Input reader]: NormalModeMoviePeriods and NormalModeMoviePointsPerPeriod must be positive and "
        "NormalModeMovieAmplitude must be positive");
  }

  if (minimizationOptions.maximumNumberOfSteps == 0 || minimizationOptions.maximumStepLength <= 0.0 ||
      minimizationOptions.maximumCellStepLength <= 0.0 || minimizationOptions.rmsGradientTolerance <= 0.0 ||
      minimizationOptions.maxGradientTolerance <= 0.0 || minimizationOptions.minimumEigenvalue <= 0.0 ||
      minimizationOptions.elasticEigenvalueTolerance <= 0.0 || minimizationOptions.convergenceFactor < 0.0)
  {
    throw std::runtime_error(
        "[Input reader]: minimization limits and tolerances must be positive and convergence factor non-negative");
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

  std::vector<std::vector<std::vector<double3>>> jsonRestartFilePositions(
      jsonNumberOfSystems, std::vector<std::vector<double3>>(jsonNumberOfComponents, std::vector<double3>()));

  std::vector<std::vector<Reaction>> jsonReactions(jsonNumberOfSystems);

  // species names of components whose 'LnPartitionFunction' is computed from the embedded
  // thermochemical databases; resolved per system once 'ExternalTemperature' is known
  std::vector<std::optional<std::string>> jsonLnPartitionFunctionSpecies(jsonNumberOfComponents);

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
      //

      if (item.contains("TranslationProbability") && item["TranslationProbability"].is_number_float())
      {
        double translationProbability = item["TranslationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::Translation, translationProbability);
        }
      }

      if (item.contains("RandomTranslationProbability") && item["RandomTranslationProbability"].is_number_float())
      {
        double randomTranslationProbability = item["RandomTranslationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::RandomTranslation, randomTranslationProbability);
        }
      }

      if (item.contains("ForceBiasTranslationProbability") && item["ForceBiasTranslationProbability"].is_number_float())
      {
        double forceBiasTranslationProbability = item["ForceBiasTranslationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::ForceBiasTranslation, forceBiasTranslationProbability);
        }
      }

      if (item.contains("RotationProbability") && item["RotationProbability"].is_number_float())
      {
        double rotationProbability = item["RotationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::Rotation, rotationProbability);
        }
      }

      if (item.contains("RandomRotationProbability") && item["RandomRotationProbability"].is_number_float())
      {
        double randomRotationProbability = item["RandomRotationProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::RandomRotation, randomRotationProbability);
        }
      }

      if (item.contains("ReinsertionProbability") && item["ReinsertionProbability"].is_number_float())
      {
        double reinsertionCBMCProbability = item["ReinsertionProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::ReinsertionCBMC, reinsertionCBMCProbability);
        }
      }

      if (item.contains("PartialReinsertionProbability") && item["PartialReinsertionProbability"].is_number_float())
      {
        double partial_reinsertion_CBMC_probability = item["PartialReinsertionProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::PartialReinsertionCBMC,
                                               partial_reinsertion_CBMC_probability);
        }
      }

      if (item.contains("IdentityChangeProbability") && item["IdentityChangeProbability"].is_number_float())
      {
        double identity_change_CBMC_probability = item["IdentityChangeProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::IdentityChangeCBMC, identity_change_CBMC_probability);
        }
      }

      if (item.contains("SwapConventionalProbability") && item["SwapConventionalProbability"].is_number_float())
      {
        double swapProbability = item["SwapConventionalProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::Swap, swapProbability);
        }
      }

      if (item.contains("SwapProbability") && item["SwapProbability"].is_number_float())
      {
        double swapCBMCProbability = item["SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::SwapCBMC, swapCBMCProbability);
        }
      }

      if (item.contains("PairSwapConventionalProbability") && item["PairSwapConventionalProbability"].is_number_float())
      {
        double pairSwapProbability = item["PairSwapConventionalProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::PairSwap, pairSwapProbability);
        }
      }

      if (item.contains("PairSwapProbability") && item["PairSwapProbability"].is_number_float())
      {
        double pairSwapCBMCProbability = item["PairSwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::PairSwapCBMC, pairSwapCBMCProbability);
        }
      }

      if (item.contains("CFCMC_PairSwapProbability") && item["CFCMC_PairSwapProbability"].is_number_float())
      {
        double pairSwapCFCMCProbability = item["CFCMC_PairSwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::PairSwapCFCMC, pairSwapCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_CBMC_PairSwapProbability") && item["CFCMC_CBMC_PairSwapProbability"].is_number_float())
      {
        double pairSwapCBCFCMCProbability = item["CFCMC_CBMC_PairSwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::PairSwapCBCFCMC, pairSwapCBCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_SwapProbability") && item["CFCMC_SwapProbability"].is_number_float())
      {
        double swapCFCMCProbability = item["CFCMC_SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::SwapCFCMC, swapCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_CBMC_SwapProbability") && item["CFCMC_CBMC_SwapProbability"].is_number_float())
      {
        double swapCBCFCMCProbability = item["CFCMC_CBMC_SwapProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::SwapCBCFCMC, swapCBCFCMCProbability);
        }
      }

      if (item.contains("GibbsSwapCFCMCProbability") && item["GibbsSwapCFCMCProbability"].is_number_float())
      {
        double gibbsSwapCFCMCProbability = item["GibbsSwapCFCMCProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsSwapCFCMC, gibbsSwapCFCMCProbability);
        }
      }

      if (item.contains("GibbsSwapCBCFCMCProbability") && item["GibbsSwapCBCFCMCProbability"].is_number_float())
      {
        double gibbsSwapCBCFCMCProbability = item["GibbsSwapCBCFCMCProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsSwapCBCFCMC, gibbsSwapCBCFCMCProbability);
        }
      }

      if (item.contains("GibbsConventionalCBCFCMCProbability") &&
          item["GibbsConventionalCBCFCMCProbability"].is_number_float())
      {
        double gibbsConventionalCBCFCMCProbability = item["GibbsConventionalCBCFCMCProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsConventionalCBCFCMC,
                                               gibbsConventionalCBCFCMCProbability);
        }
      }

      if (item.contains("GibbsConventionalCFCMCProbability") &&
          item["GibbsConventionalCFCMCProbability"].is_number_float())
      {
        double gibbsConventionalCFCMCProbability = item["GibbsConventionalCFCMCProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsConventionalCFCMC, gibbsConventionalCFCMCProbability);
        }
      }

      if (item.contains("GibbsSwapCBMCProbability") && item["GibbsSwapCBMCProbability"].is_number_float())
      {
        double gibbsSwapCBMCProbability = item["GibbsSwapCBMCProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsSwapCBMC, gibbsSwapCBMCProbability);
        }
      }

      if (item.contains("GibbsIdentityChangeProbability") && item["GibbsIdentityChangeProbability"].is_number_float())
      {
        double gibbsIdentityChangeCBMCProbability = item["GibbsIdentityChangeProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::GibbsIdentityChangeCBMC,
                                               gibbsIdentityChangeCBMCProbability);
        }
      }

      if (item.contains("WidomProbability") && item["WidomProbability"].is_number_float())
      {
        double widomProbability = item["WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::Widom, widomProbability);
        }
      }

      if (item.contains("CFCMC_WidomProbability") && item["CFCMC_WidomProbability"].is_number_float())
      {
        double widomCFCMCProbability = item["CFCMC_WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::WidomCFCMC, widomCFCMCProbability);
        }
      }

      if (item.contains("CFCMC_CBMC_WidomProbability") && item["CFCMC_CBMC_WidomProbability"].is_number_float())
      {
        double widomCBCFCMCProbability = item["CFCMC_CBMC_WidomProbability"].get<double>();
        for (std::size_t i = 0; i < move_probabilities.size(); ++i)
        {
          move_probabilities[i].setProbability(Move::Types::WidomCBCFCMC, widomCBCFCMCProbability);
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

      if (item.contains("PairComponent") && item["PairComponent"].is_number_integer())
      {
        std::size_t pair_component = item["PairComponent"].get<std::size_t>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].pairComponentId = pair_component;
        }
      }

      if (item.contains("MaximumPairDistance") && item["MaximumPairDistance"].is_number_float())
      {
        double maximum_pair_distance = item["MaximumPairDistance"].get<double>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].maximumPairDistance = maximum_pair_distance;
        }
      }

      if (item.contains("IdentityChanges") && item["IdentityChanges"].is_array())
      {
        std::vector<std::size_t> identity_changes = item["IdentityChanges"].get<std::vector<std::size_t>>();
        if (std::ranges::contains(identity_changes, componentId))
        {
          throw std::runtime_error(
              std::format("[Input reader]: component '{}' (id {}) cannot list itself in 'IdentityChanges'\n",
                          jsonComponentName, componentId));
        }
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].identityChanges = identity_changes;
        }
      }

      if (item.contains("GibbsIdentityChanges") && item["GibbsIdentityChanges"].is_array())
      {
        std::vector<std::size_t> gibbs_identity_changes = item["GibbsIdentityChanges"].get<std::vector<std::size_t>>();
        if (std::ranges::contains(gibbs_identity_changes, componentId))
        {
          throw std::runtime_error(
              std::format("[Input reader]: component '{}' (id {}) cannot list itself in 'GibbsIdentityChanges'\n",
                          jsonComponentName, componentId));
        }
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].gibbsIdentityChanges = gibbs_identity_changes;
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

      // selects which lambda the dU/dlambda group-tagging (thermodynamic integration) applies to:
      // "CFCMC" (default lambdaGC, also used by CB/CFCMC swap), "CFCMC_PairSwap", or "CFCMC_CBMC_PairSwap"
      if (item.contains("ThermodynamicIntegration") && item["ThermodynamicIntegration"].is_string())
      {
        std::string thermodynamicIntegrationString = item["ThermodynamicIntegration"].get<std::string>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          if (caseInSensStringCompare(thermodynamicIntegrationString, "CFCMC") ||
              caseInSensStringCompare(thermodynamicIntegrationString, "CFCMC_Swap") ||
              caseInSensStringCompare(thermodynamicIntegrationString, "CFCMC_CBMC_Swap"))
          {
            jsonComponents[i][componentId].lambdaGC.computeDUdlambda = true;
          }
          else if (caseInSensStringCompare(thermodynamicIntegrationString, "CFCMC_PairSwap"))
          {
            jsonComponents[i][componentId].lambdaPairSwap.computeDUdlambda = true;
          }
          else if (caseInSensStringCompare(thermodynamicIntegrationString, "CFCMC_CBMC_PairSwap"))
          {
            jsonComponents[i][componentId].lambdaPairSwapCB.computeDUdlambda = true;
          }
          else
          {
            throw std::runtime_error(
                std::format("Error [Component {}]: unknown 'ThermodynamicIntegration' value '{}' "
                            "(allowed: 'CFCMC', 'CFCMC_PairSwap', 'CFCMC_CBMC_PairSwap')\n",
                            componentId, thermodynamicIntegrationString));
          }
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

      if (item.contains("LnPartitionFunction") && item["LnPartitionFunction"].is_array())
      {
        std::vector<double> ln_partition_functions =
            parseList<double>(jsonNumberOfSystems, "LnPartitionFunction", item["LnPartitionFunction"]);
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].lnPartitionFunction = ln_partition_functions[i];
        }
      }
      else if (item.contains("LnPartitionFunction") && item["LnPartitionFunction"].is_number())
      {
        const double lnPartitionFunction = item["LnPartitionFunction"].get<double>();
        for (std::size_t i = 0; i != jsonNumberOfSystems; ++i)
        {
          jsonComponents[i][componentId].lnPartitionFunction = lnPartitionFunction;
        }
      }
      // "auto" computes ln(q/V) from the embedded thermochemical databases using the component
      // name; any other string is used as the species name for the lookup. The value is evaluated
      // at the 'ExternalTemperature' of each system. The database is selected per system with
      // 'ThermochemicalDatabase' (default: "NASA", NASA polynomials with JANAF fallback).
      else if (item.contains("LnPartitionFunction") && item["LnPartitionFunction"].is_string())
      {
        std::string speciesName = item["LnPartitionFunction"].get<std::string>();
        if (caseInSensStringCompare(speciesName, "auto"))
        {
          speciesName = jsonComponentName;
        }
        jsonLnPartitionFunctionSpecies[componentId] = speciesName;
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
        forceFields[systemId]->preComputeDerivedParameters();
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
        forceFields[systemId]->preComputeDerivedParameters();
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
        forceFields[systemId]->preComputeDerivedParameters();
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
        forceFields[systemId]->preComputeDerivedParameters();
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

        if (caseInSensStringCompare(chargeMethodString, "None"))
        {
          forceFields[systemId]->chargeMethod = ForceField::ChargeMethod::Ewald;
          forceFields[systemId]->useCharge = false;
        }
        else
        {
          forceFields[systemId]->chargeMethod = ForceField::chargeMethodFromString(chargeMethodString);
          forceFields[systemId]->useCharge = true;
        }
      }

      if (value.contains("ModifiedShiftedForceBeta") && value["ModifiedShiftedForceBeta"].is_number())
      {
        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        forceFields[systemId]->modifiedShiftedForceBeta = value["ModifiedShiftedForceBeta"].get<double>();
        if (forceFields[systemId]->modifiedShiftedForceBeta <= 0.0)
        {
          throw std::runtime_error("[Input reader]: ModifiedShiftedForceBeta must be positive\n");
        }
      }

      CIFReader::UseChargesFrom useChargesFrom{CIFReader::UseChargesFrom::PseudoAtoms};
      if (value.contains("UseChargesFrom") && value["UseChargesFrom"].is_string())
      {
        std::string useChargesFromString = value["UseChargesFrom"].get<std::string>();

        if (caseInSensStringCompare(useChargesFromString, "PsuedoAtoms"))
        {
          useChargesFrom = CIFReader::UseChargesFrom::PseudoAtoms;
        }
        if (caseInSensStringCompare(useChargesFromString, "CIF_File"))
        {
          useChargesFrom = CIFReader::UseChargesFrom::CIF_File;
        }
        if (caseInSensStringCompare(useChargesFromString, "ChargeEquilibration"))
        {
          useChargesFrom = CIFReader::UseChargesFrom::ChargeEquilibration;
        }
      }

      if (value.contains("VolumeMoveProbability") && value["VolumeMoveProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::VolumeChange, value["VolumeMoveProbability"].get<double>());
      }

      if (value.contains("AnisotropicVolumeMoveProbability") &&
          value["AnisotropicVolumeMoveProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::AnisotropicVolumeChange,
                                              value["AnisotropicVolumeMoveProbability"].get<double>());
      }

      if (value.contains("GibbsVolumeMoveProbability") && value["GibbsVolumeMoveProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::GibbsVolume,
                                              value["GibbsVolumeMoveProbability"].get<double>());
      }

      if (value.contains("ParallelTemperingSwapProbability") &&
          value["ParallelTemperingSwapProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ParallelTempering,
                                              value["ParallelTemperingSwapProbability"].get<double>());
      }
      if (value.contains("HybridMCProbability") && value["HybridMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::HybridMC, value["HybridMCProbability"].get<double>());
      }

      if (value.contains("ForceBiasTranslationAllProbability") &&
          value["ForceBiasTranslationAllProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ForceBiasTranslationAll,
                                              value["ForceBiasTranslationAllProbability"].get<double>());
      }

      if (value.contains("ReactionCBMCProbability") && value["ReactionCBMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ReactionCBMC,
                                              value["ReactionCBMCProbability"].get<double>());
      }

      if (value.contains("ReactionConventionalCFCMCProbability") &&
          value["ReactionConventionalCFCMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ReactionConventionalCFCMC,
                                              value["ReactionConventionalCFCMCProbability"].get<double>());
      }

      if (value.contains("ReactionConventionalCBCFCMCProbability") &&
          value["ReactionConventionalCBCFCMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ReactionConventionalCBCFCMC,
                                              value["ReactionConventionalCBCFCMCProbability"].get<double>());
      }

      if (value.contains("ReactionCFCMCProbability") && value["ReactionCFCMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ReactionCFCMC,
                                              value["ReactionCFCMCProbability"].get<double>());
      }

      if (value.contains("ReactionCBCFCMCProbability") && value["ReactionCBCFCMCProbability"].is_number_float())
      {
        mc_moves_probabilities.setProbability(Move::Types::ReactionCBCFCMC,
                                              value["ReactionCBCFCMCProbability"].get<double>());
      }

      if (value.contains("Reactions") && value["Reactions"].is_array())
      {
        double defaultMaximumLambdaChange = 0.3;
        if (value.contains("MaximumReactionLambdaChange") && value["MaximumReactionLambdaChange"].is_number())
        {
          defaultMaximumLambdaChange = value["MaximumReactionLambdaChange"].get<double>();
        }

        double defaultMaximumLambdaChangeProducts = 0.3;
        if (value.contains("MaximumReactionLambdaChangeProducts") &&
            value["MaximumReactionLambdaChangeProducts"].is_number())
        {
          defaultMaximumLambdaChangeProducts = value["MaximumReactionLambdaChangeProducts"].get<double>();
        }

        double defaultLambdaSwitchPoint = 0.3;
        if (value.contains("LambdaSwitchPoint") && value["LambdaSwitchPoint"].is_number())
        {
          defaultLambdaSwitchPoint = value["LambdaSwitchPoint"].get<double>();
        }

        // Each reaction is driven by one reaction move, declared per reaction with the "Move" key.
        // The serial moves are ReactionCFCMC and ReactionCBCFCMC, the parallel moves the 'conventional'
        // variants. When the enabled reaction moves are unambiguous the "Move" key may be omitted.
        const bool serialRxCFCEnabled = mc_moves_probabilities.getProbability(Move::Types::ReactionCFCMC) > 0.0 ||
                                        mc_moves_probabilities.getProbability(Move::Types::ReactionCBCFCMC) > 0.0;
        const bool parallelRxCFCEnabled =
            mc_moves_probabilities.getProbability(Move::Types::ReactionConventionalCFCMC) > 0.0 ||
            mc_moves_probabilities.getProbability(Move::Types::ReactionConventionalCBCFCMC) > 0.0;
        const std::array reactionMoves{Move::Types::ReactionCBMC, Move::Types::ReactionCFCMC,
                                       Move::Types::ReactionCBCFCMC, Move::Types::ReactionConventionalCFCMC,
                                       Move::Types::ReactionConventionalCBCFCMC};
        std::vector<Move::Types> enabledReactionMoves;
        for (const Move::Types reactionMove : reactionMoves)
        {
          if (mc_moves_probabilities.getProbability(reactionMove) > 0.0)
          {
            enabledReactionMoves.push_back(reactionMove);
          }
        }

        const auto parseReactionMove = [](const std::string& moveName) -> std::optional<Move::Types>
        {
          if (moveName == "ReactionCBMC") return Move::Types::ReactionCBMC;
          if (moveName == "ReactionCFCMC") return Move::Types::ReactionCFCMC;
          if (moveName == "ReactionCBCFCMC") return Move::Types::ReactionCBCFCMC;
          if (moveName == "ReactionConventionalCFCMC") return Move::Types::ReactionConventionalCFCMC;
          if (moveName == "ReactionConventionalCBCFCMC") return Move::Types::ReactionConventionalCBCFCMC;
          return std::nullopt;
        };

        std::size_t reactionId = 0;
        for (const auto& reactionItem : value["Reactions"])
        {
          if (!reactionItem.contains("Reactants") || !reactionItem.contains("Products"))
          {
            throw std::runtime_error(
                std::format("[Input reader]: each reaction must contain 'Reactants' and 'Products'\n"));
          }
          Reaction reaction(reactionId, reactionItem["Reactants"].get<std::vector<std::size_t>>(),
                            reactionItem["Products"].get<std::vector<std::size_t>>());
          reaction.maximumLambdaChange = defaultMaximumLambdaChange;
          reaction.maximumLambdaChangeProducts = defaultMaximumLambdaChangeProducts;
          reaction.lambdaSwitchPoint = defaultLambdaSwitchPoint;

          if (reactionItem.contains("MaximumReactionLambdaChange") &&
              reactionItem["MaximumReactionLambdaChange"].is_number())
          {
            reaction.maximumLambdaChange = reactionItem["MaximumReactionLambdaChange"].get<double>();
          }
          if (reactionItem.contains("MaximumReactionLambdaChangeProducts") &&
              reactionItem["MaximumReactionLambdaChangeProducts"].is_number())
          {
            reaction.maximumLambdaChangeProducts = reactionItem["MaximumReactionLambdaChangeProducts"].get<double>();
          }
          if (reactionItem.contains("LambdaSwitchPoint") && reactionItem["LambdaSwitchPoint"].is_number())
          {
            reaction.lambdaSwitchPoint = reactionItem["LambdaSwitchPoint"].get<double>();
          }

          if (reactionItem.contains("Move"))
          {
            if (!reactionItem["Move"].is_string())
            {
              throw std::runtime_error(
                  std::format("[Input reader]: reaction {} key 'Move' must be a string\n", reactionId));
            }
            const std::string moveName = reactionItem["Move"].get<std::string>();
            const std::optional<Move::Types> parsedMove = parseReactionMove(moveName);
            if (!parsedMove.has_value())
            {
              throw std::runtime_error(
                  std::format("[Input reader]: reaction {} has unknown 'Move' value '{}'; valid values are "
                              "'ReactionCBMC', 'ReactionCFCMC', 'ReactionCBCFCMC', 'ReactionConventionalCFCMC', "
                              "'ReactionConventionalCBCFCMC'\n",
                              reactionId, moveName));
            }
            reaction.reactionMove = parsedMove.value();
          }
          else if (enabledReactionMoves.size() > 1)
          {
            throw std::runtime_error(
                std::format("[Input reader]: multiple reaction moves are enabled; reaction {} must declare "
                            "which exact move drives it with the 'Move' key\n",
                            reactionId));
          }
          else if (enabledReactionMoves.size() == 1)
          {
            reaction.reactionMove = enabledReactionMoves.front();
          }
          else
          {
            reaction.reactionMove = Move::Types::ReactionCBMC;
          }

          if (reaction.isSerialRxCFC() && !serialRxCFCEnabled)
          {
            throw std::runtime_error(
                std::format("[Input reader]: reaction {} declares serial move '{}' but neither "
                            "'ReactionCFCMCProbability' nor 'ReactionCBCFCMCProbability' is set\n",
                            reactionId, Move::moveNames[std::to_underlying(reaction.reactionMove)]));
          }
          if (reaction.isParallelRxCFC() && !parallelRxCFCEnabled)
          {
            throw std::runtime_error(
                std::format("[Input reader]: reaction {} declares parallel move '{}' but neither "
                            "'ReactionConventionalCFCMCProbability' nor 'ReactionConventionalCBCFCMCProbability' "
                            "is set\n",
                            reactionId, Move::moveNames[std::to_underlying(reaction.reactionMove)]));
          }
          if (mc_moves_probabilities.getProbability(reaction.reactionMove) <= 0.0 && !enabledReactionMoves.empty())
          {
            throw std::runtime_error(
                std::format("[Input reader]: reaction {} declares move '{}' but its probability is not enabled\n",
                            reactionId, Move::moveNames[std::to_underlying(reaction.reactionMove)]));
          }

          jsonReactions[systemId].push_back(reaction);
          ++reactionId;
        }
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

      // resolve 'LnPartitionFunction' entries given as species names at this system's temperature
      PartitionFunction::Source thermochemicalDatabase = PartitionFunction::Source::NASA;
      if (value.contains("ThermochemicalDatabase") && value["ThermochemicalDatabase"].is_string())
      {
        const std::string databaseString = value["ThermochemicalDatabase"].get<std::string>();
        if (caseInSensStringCompare(databaseString, "NASA"))
        {
          thermochemicalDatabase = PartitionFunction::Source::NASA;
        }
        else if (caseInSensStringCompare(databaseString, "JANAF"))
        {
          thermochemicalDatabase = PartitionFunction::Source::JANAF;
        }
        else
        {
          throw std::runtime_error(
              std::format("[Input reader]: system {}: unknown 'ThermochemicalDatabase' '{}'; use \"NASA\" or "
                          "\"JANAF\"\n",
                          systemId, databaseString));
        }
      }

      for (std::size_t i = 0; i != jsonNumberOfComponents; ++i)
      {
        if (jsonLnPartitionFunctionSpecies[i].has_value())
        {
          const std::string& speciesName = jsonLnPartitionFunctionSpecies[i].value();
          if (!PartitionFunction::contains(speciesName, thermochemicalDatabase))
          {
            throw std::runtime_error(
                std::format("[Input reader]: system {} component {}: no thermochemical data found for species "
                            "'{}' in the {} database; use a different species name, database, or specify "
                            "'LnPartitionFunction' as a number\n",
                            systemId, i, speciesName,
                            thermochemicalDatabase == PartitionFunction::Source::NASA ? "NASA/JANAF" : "JANAF"));
          }
          jsonComponents[systemId][i].lnPartitionFunction =
              PartitionFunction::logPartitionFunction(speciesName, T, thermochemicalDatabase);
        }
      }

      bool hasExternalField = false;
      if (value.contains("ExternalField") && value["ExternalField"].is_boolean())
      {
        hasExternalField = value["ExternalField"].get<bool>();
      }

      std::optional<SimulationBox> restart_simulation_box{};
      if (value.contains("RestartFileName") && value["RestartFileName"].is_string())
      {
        std::string restartFileName = value["RestartFileName"].get<std::string>();

        if (!std::filesystem::exists(restartFileName))
        {
          restartFileName += ".json";

          if (!std::filesystem::exists(restartFileName))
          {
            throw std::runtime_error(std::format("[Input reader]: File '{}' not found\n", restartFileName));
          }
        }

        std::ifstream input_restart_file(restartFileName);

        nlohmann::basic_json<nlohmann::raspa_map> restart_data{};

        try
        {
          restart_data = nlohmann::json::parse(input_restart_file);
        }
        catch (nlohmann::json::parse_error& ex)
        {
          std::cerr << "parse error at byte " << ex.byte << std::endl;
        }

        for (std::size_t i = 0; i < jsonComponents[systemId].size(); ++i)
        {
          std::string component_name = jsonComponents[systemId][i].name;
          if (restart_data.contains(component_name))
          {
            std::vector<double3> restart_positions = restart_data[component_name].get<std::vector<double3>>();
            jsonRestartFilePositions[systemId][i] = restart_positions;
          }
        }

        if (restart_data.contains("SimulationBox"))
        {
          restart_simulation_box = restart_data["SimulationBox"].get<SimulationBox>();
        }
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

        const std::string file_content = readFileContent(frameworkNameString, ".cif");

        if (const auto cif = CIFReader::readCIFString(file_content, forceFields[systemId].value(), useChargesFrom);
            cif.has_value())
        {
          auto [simulation_box, space_group_hall_symbol, defined_atoms, fractional_atoms_unit_cell] = cif.value();
          Framework framework =
              Framework(forceFields[systemId].value(), frameworkNameString, simulation_box, space_group_hall_symbol,
                        defined_atoms, fractional_atoms_unit_cell, jsonNumberOfUnitCells);

          if (value.contains("FrameworkDefinition"))
          {
            if (!value["FrameworkDefinition"].is_string())
            {
              throw std::runtime_error("[Input reader]: 'FrameworkDefinition' must have a value of string type\n");
            }
            framework.readFrameworkDefinition(forceFields[systemId].value(),
                                              value["FrameworkDefinition"].get<std::string>());
          }

          if (value.contains("FrameworkType"))
          {
            if (!value["FrameworkType"].is_string())
            {
              throw std::runtime_error("[Input reader]: 'FrameworkType' must have a value of string type\n");
            }
            const std::string frameworkType = value["FrameworkType"].get<std::string>();
            if (caseInSensStringCompare(frameworkType, "Flexible"))
            {
              framework.rigid = false;
            }
            else if (caseInSensStringCompare(frameworkType, "Rigid"))
            {
              framework.rigid = true;
            }
            else
            {
              throw std::runtime_error(std::format(
                  "[Input reader]: unknown 'FrameworkType' '{}'; expected 'Rigid' or 'Flexible'\n", frameworkType));
            }
          }

          std::optional<Framework> jsonFrameworkComponents{framework};

          // create system
          systems[systemId] =
              System(forceFields[systemId].value(), std::nullopt, hasExternalField, T, P, heliumVoidFraction,
                     jsonFrameworkComponents, jsonComponents[systemId], jsonRestartFilePositions[systemId],
                     jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
        }
        else if (cif.error() == CIFReader::ParseError::invalidInput)
        {
          std::print("Invalid input\n");
          std::exit(-1);
        }
        else if (cif.error() == CIFReader::ParseError::invalidForceField)
        {
          std::print("not all atoms defined in CIF-file\n");
          std::exit(-1);
        }
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

        if (!forceFields[systemId].has_value())
        {
          throw std::runtime_error(std::format("[Input reader]: No forcefield specified or found'\n"));
        }
        SimulationBox simulationBox{boxLengths.x, boxLengths.y, boxLengths.z, boxAngles.x, boxAngles.y, boxAngles.z};

        if (restart_simulation_box.has_value())
        {
          simulationBox = restart_simulation_box.value();
        }

        systems[systemId] = System(forceFields[systemId].value(), simulationBox, hasExternalField, T, P, 1.0, {},
                                   jsonComponents[systemId], jsonRestartFilePositions[systemId],
                                   jsonCreateNumberOfMolecules[systemId], jsonNumberOfBlocks, mc_moves_probabilities);
      }
      else
      {
        throw std::runtime_error(
            std::format("[Input reader]: system key 'Type' must have value 'Box' or 'Framework'\n"));
      }

      {
        if (value.contains("CellType"))
        {
          if (!value["CellType"].is_string())
          {
            throw std::runtime_error("[Input reader]: 'CellType' must have a value of string type\n");
          }
          const std::string cellType = value["CellType"].get<std::string>();
          const std::optional<CellMinimizationType> parsedCellType = cellMinimizationTypeFromString(cellType);
          if (!parsedCellType.has_value())
          {
            throw std::runtime_error(
                std::format("[Input reader]: unknown 'CellType' '{}'; expected 'Fixed', 'Isotropic', 'Anisotropic', "
                            "'Monoclinic', 'Regular', 'RegularUpperTriangle', or 'MonoclinicUpperTriangle'\n",
                            cellType));
          }
          systems[systemId].cellMinimizationType = parsedCellType.value();
        }
        if (value.contains("MonoclinicAngleType"))
        {
          if (!value["MonoclinicAngleType"].is_string())
          {
            throw std::runtime_error("[Input reader]: 'MonoclinicAngleType' must have a value of string type\n");
          }
          const std::string angleType = value["MonoclinicAngleType"].get<std::string>();
          const std::optional<MonoclinicAngleType> parsedAngleType = monoclinicAngleTypeFromString(angleType);
          if (!parsedAngleType.has_value())
          {
            throw std::runtime_error(std::format(
                "[Input reader]: unknown 'MonoclinicAngleType' '{}'; expected 'Alpha', 'Beta', or 'Gamma'\n",
                angleType));
          }
          systems[systemId].monoclinicAngleType = parsedAngleType.value();
        }

        double3 inputPressureDiagonal(systems[systemId].input_pressure, systems[systemId].input_pressure,
                                      systems[systemId].input_pressure);
        if (value.contains("ExternalPressureX") && value["ExternalPressureX"].is_number())
        {
          inputPressureDiagonal.x = value["ExternalPressureX"].get<double>();
        }
        if (value.contains("ExternalPressureY") && value["ExternalPressureY"].is_number())
        {
          inputPressureDiagonal.y = value["ExternalPressureY"].get<double>();
        }
        if (value.contains("ExternalPressureZ") && value["ExternalPressureZ"].is_number())
        {
          inputPressureDiagonal.z = value["ExternalPressureZ"].get<double>();
        }
        systems[systemId].input_pressureTensorDiagonal = inputPressureDiagonal;
        systems[systemId].pressureTensorDiagonal = inputPressureDiagonal / Units::PressureConversionFactor;
        if (systems[systemId].cellMinimizationType != CellMinimizationType::Fixed &&
            (value.contains("ExternalPressureX") || value.contains("ExternalPressureY") ||
             value.contains("ExternalPressureZ")))
        {
          throw std::runtime_error(
              "[Input reader]: variable-cell minimization supports hydrostatic 'ExternalPressure' only");
        }
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

      systems[systemId].reactions.list = std::move(jsonReactions[systemId]);
      for (const Reaction& reaction : systems[systemId].reactions.list)
      {
        const bool conventional = reaction.reactionMove == Move::Types::ReactionConventionalCFCMC ||
                                  reaction.reactionMove == Move::Types::ReactionCFCMC;
        if (!conventional)
        {
          continue;
        }

        const std::size_t componentCount =
            std::min({systems[systemId].components.size(), reaction.reactantStoichiometry.size(),
                      reaction.productStoichiometry.size()});
        for (std::size_t componentId = 0; componentId < componentCount; ++componentId)
        {
          const bool participates =
              reaction.reactantStoichiometry[componentId] > 0 || reaction.productStoichiometry[componentId] > 0;
          if (participates && systems[systemId].components[componentId].growType != Component::GrowType::Rigid)
          {
            throw std::runtime_error(
                std::format("[Input reader]: reaction {} uses conventional move '{}' with flexible component "
                            "'{}'; conventional reaction moves require every participating component to have "
                            "fixed internal conformation (component Type 'Rigid'). Use '{}' for flexible "
                            "components\n",
                            reaction.id, Move::moveNames[std::to_underlying(reaction.reactionMove)],
                            systems[systemId].components[componentId].name,
                            reaction.reactionMove == Move::Types::ReactionCFCMC ? "ReactionCBCFCMC"
                                                                                : "ReactionConventionalCBCFCMC"));
          }
        }
      }
      if (!systems[systemId].reactions.list.empty() && systems[systemId].usesReactionConventionalCFCMC())
      {
        systems[systemId].createReactionFractionalMolecules();
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
              jsonNumberOfBlocks, systems[systemId].components.size(),
              {minimumRangeNumberOfMoleculesHistogram, maximumRangeNumberOfMoleculesHistogram},
              sampleNumberOfMoleculesHistogramEvery, writeNumberOfMoleculesHistogramEvery);
        }
      }

      if (value.contains("ComputeMoleculeProperties") && value["ComputeMoleculeProperties"].is_boolean())
      {
        if (value["ComputeMoleculeProperties"].get<bool>())
        {
          std::size_t sampleMoleculePropertiesEvery{10};
          if (value.contains("SampleMoleculePropertiesEvery") &&
              value["SampleMoleculePropertiesEvery"].is_number_unsigned())
          {
            sampleMoleculePropertiesEvery = value["SampleMoleculePropertiesEvery"].get<std::size_t>();
          }

          std::size_t writeMoleculePropertiesEvery{5000};
          if (value.contains("WriteMoleculePropertiesEvery") &&
              value["WriteMoleculePropertiesEvery"].is_number_unsigned())
          {
            writeMoleculePropertiesEvery = value["WriteMoleculePropertiesEvery"].get<std::size_t>();
          }

          std::size_t numberOfBinsMoleculeProperties{128};
          if (value.contains("NumberOfBinsMoleculeProperties") &&
              value["NumberOfBinsMoleculeProperties"].is_number_unsigned())
          {
            numberOfBinsMoleculeProperties = value["NumberOfBinsMoleculeProperties"].get<std::size_t>();
          }

          double bondRangeMoleculeProperties{4.0};
          if (value.contains("BondRangeMoleculeProperties") && value["BondRangeMoleculeProperties"].is_number_float())
          {
            bondRangeMoleculeProperties = value["BondRangeMoleculeProperties"].get<double>();
          }

          systems[systemId].propertyMoleculeProperties = PropertyMoleculeProperties(
              jsonNumberOfBlocks, systems[systemId].components, numberOfBinsMoleculeProperties,
              bondRangeMoleculeProperties, sampleMoleculePropertiesEvery, writeMoleculePropertiesEvery);
        }
      }

      if (value.contains("ComputeNumberOfMoleculesEvolution") &&
          value["ComputeNumberOfMoleculesEvolution"].is_boolean())
      {
        if (value["ComputeNumberOfMoleculesEvolution"].get<bool>())
        {
          std::size_t sample_number_of_molecules_every{1};
          if (value.contains("SampleNumberOfMoleculesEvolutionEvery") &&
              value["SampleNumberOfMoleculesEvolutionEvery"].is_number_unsigned())
          {
            sample_number_of_molecules_every = value["SampleNumberOfMoleculesEvolutionEvery"].get<std::size_t>();
          }

          std::size_t write_number_of_molecules_evolution_every{5000};
          if (value.contains("WriteNumberOfMoleculesEvolutionEvery") &&
              value["WriteNumberOfMoleculesEvolutionEvery"].is_number_unsigned())
          {
            write_number_of_molecules_evolution_every =
                value["WriteNumberOfMoleculesEvolutionEvery"].get<std::size_t>();
          }

          systems[systemId].propertyNumberOfMoleculesEvolution = PropertyNumberOfMoleculesEvolution(
              numberOfProductionCycles + numberOfInitializationCycles + numberOfEquilibrationCycles,
              systems[systemId].components.size(), sample_number_of_molecules_every,
              write_number_of_molecules_evolution_every);
        }
      }

      if (value.contains("ComputeVolumeEvolution") && value["ComputeVolumeEvolution"].is_boolean())
      {
        if (value["ComputeVolumeEvolution"].get<bool>())
        {
          std::size_t sample_volume_evolution_every{1};
          if (value.contains("SampleVolumeEvolutionEvery") && value["SampleVolumeEvolutionEvery"].is_number_unsigned())
          {
            sample_volume_evolution_every = value["SampleVolumeEvolutionEvery"].get<std::size_t>();
          }

          std::size_t write_volume_evolution_every{5000};
          if (value.contains("writeVolumeEvolutionEvery") && value["writeVolumeEvolutionEvery"].is_number_unsigned())
          {
            write_volume_evolution_every = value["writeVolumeEvolutionEvery"].get<std::size_t>();
          }

          systems[systemId].propertyVolumeEvolution = PropertyVolumeEvolution(
              numberOfProductionCycles + numberOfInitializationCycles + numberOfEquilibrationCycles,
              sample_volume_evolution_every, write_volume_evolution_every);
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

      // parse the time step before it is captured by the MSD/VACF properties and the thermostat below
      if (value.contains("TimeStep") && value["TimeStep"].is_number_float())
      {
        systems[systemId].timeStep = value["TimeStep"].get<double>();
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
              systems[systemId].numberOfMoleculesPerComponent, systems[systemId].moleculeData.size(),
              systems[systemId].timeStep, numberOfBlockElementsMSD, sampleMSDEvery, writeMSDEvery);
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
              systems[systemId].numberOfMoleculesPerComponent, systems[systemId].moleculeData.size(),
              systems[systemId].timeStep, numberOfBuffersVACF, bufferLengthVACF, sampleVACFEvery, writeVACFEvery);
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

          PropertyDensityGrid::Binning binning = PropertyDensityGrid::Binning::Standard;
          if (value.contains("DensityGridBinning") && value["DensityGridBinning"].is_string())
          {
            std::string binningString = value["DensityGridBinning"].get<std::string>();
            if (caseInSensStringCompare(binningString, "Equitable"))
            {
              binning = PropertyDensityGrid::Binning::Equitable;
            }
            else if (caseInSensStringCompare(binningString, "Standard"))
            {
              binning = PropertyDensityGrid::Binning::Standard;
            }
            else
            {
              throw std::runtime_error(
                  std::format("Error: DensityGridBinning must be 'Standard' or 'Equitable', got '{}'", binningString));
            }
          }

          systems[systemId].propertyDensityGrid = PropertyDensityGrid(
              systems[systemId].framework ? 1 : 0, systems[systemId].components.size(), densityGridSize,
              sampleDensityGridEvery, writeDensityGridEvery, densityGridPseudoAtomsList, norm, binning);
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

          bool restrict_to_box{true};
          if (value.contains("RestrictMoviePositionsToBox") && value["RestrictMoviePositionsToBox"].is_boolean())
          {
            restrict_to_box = value["RestrictMoviePositionsToBox"].get<bool>();
          }

          systems[systemId].samplePDBMovie = SampleMovie(systemId, sampleMovieEvery, restrict_to_box, std::nullopt);
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
        std::size_t numberOfRespaSteps{5};
        std::size_t numberOfYoshidaSuzukiSteps{5};
        double timeScaleParameterThermostat{0.15};
        std::size_t barostatChainLength{5};
        double timeScaleParameterBarostat{1.0};

        if (value.contains("ThermostatChainLength") && value["ThermostatChainLength"].is_number_unsigned())
          thermostatChainLength = value["ThermostatChainLength"].get<std::size_t>();
        if (value.contains("BarostatChainLength") && value["BarostatChainLength"].is_number_unsigned())
          barostatChainLength = value["BarostatChainLength"].get<std::size_t>();
        if (value.contains("NumberOfRespaSteps") && value["NumberOfRespaSteps"].is_number_unsigned())
          numberOfRespaSteps = value["NumberOfRespaSteps"].get<std::size_t>();
        if (value.contains("NumberOfYoshidaSuzukiSteps") && value["NumberOfYoshidaSuzukiSteps"].is_number_unsigned())
          numberOfYoshidaSuzukiSteps = value["NumberOfYoshidaSuzukiSteps"].get<std::size_t>();
        if (value.contains("TimeScaleParameterThermostat") && value["TimeScaleParameterThermostat"].is_number())
          timeScaleParameterThermostat = value["TimeScaleParameterThermostat"].get<double>();
        if (value.contains("TimeScaleParameterBarostat") && value["TimeScaleParameterBarostat"].is_number())
          timeScaleParameterBarostat = value["TimeScaleParameterBarostat"].get<double>();
        if (thermostatChainLength == 0 || barostatChainLength == 0 || numberOfRespaSteps == 0)
          throw std::runtime_error(
              "[Input reader]: thermostat/barostat chain lengths and NumberOfRespaSteps must be positive");

        std::string ensembleString = value["Ensemble"].get<std::string>();
        const std::optional<MolecularDynamicsEnsemble> ensemble = molecularDynamicsEnsembleFromString(ensembleString);
        if (!ensemble.has_value())
          throw std::runtime_error(std::format(
              "[Input reader]: unknown MD 'Ensemble' '{}'; expected NVE, NVT, NPT, NPTPR, MuVT, MuPT, or MuPTPR",
              ensembleString));
        systems[systemId].molecularDynamicsEnsemble = ensemble.value();
        if (molecularDynamicsUsesThermostat(ensemble.value()))
        {
          systems[systemId].thermostat =
              Thermostat(systems[systemId].temperature, systems[systemId].timeStep,
                         systems[systemId].translationalDegreesOfFreedom, systems[systemId].rotationalDegreesOfFreedom,
                         thermostatChainLength, numberOfYoshidaSuzukiSteps, timeScaleParameterThermostat);
          systems[systemId].thermostat->numberOfRespaSteps = numberOfRespaSteps;
        }
        if (molecularDynamicsHasParticleExchange(ensemble.value()))
        {
          if (!value.contains("ExternalPressure") || !value["ExternalPressure"].is_number())
            throw std::runtime_error("[Input reader]: MuVT/MuPT/MuPTPR requires a numeric ExternalPressure");
          const bool hasSwapComponent = std::ranges::any_of(systems[systemId].components, [](const Component& component) {
            return component.mc_moves_probabilities.getProbability(Move::Types::SwapCBMC) > 0.0;
          });
          if (!hasSwapComponent)
            throw std::runtime_error(
                "[Input reader]: MuVT/MuPT/MuPTPR requires SwapProbability > 0 for at least one component");
        }
        if (molecularDynamicsUsesIsotropicBarostat(ensemble.value()) ||
            molecularDynamicsUsesFlexibleBarostat(ensemble.value()))
        {
          if (!value.contains("ExternalPressure") || !value["ExternalPressure"].is_number())
            throw std::runtime_error("[Input reader]: NPT/NPTPR/MuPT/MuPTPR requires a numeric ExternalPressure");
          if (value.contains("ExternalPressureX") || value.contains("ExternalPressureY") ||
              value.contains("ExternalPressureZ"))
            throw std::runtime_error("[Input reader]: NPT/NPTPR/MuPT/MuPTPR supports hydrostatic ExternalPressure only");
          CellMinimizationType cellType = systems[systemId].cellMinimizationType;
          if (molecularDynamicsUsesIsotropicBarostat(ensemble.value())) cellType = CellMinimizationType::Isotropic;
          if (molecularDynamicsUsesFlexibleBarostat(ensemble.value()) && cellType == CellMinimizationType::Fixed)
            cellType = CellMinimizationType::Regular;
          systems[systemId].thermobarostat = Thermobarostat(
              ensemble.value(), cellType, systems[systemId].monoclinicAngleType, systems[systemId].temperature,
              systems[systemId].pressure, systems[systemId].timeStep, systems[systemId].translationalDegreesOfFreedom,
              barostatChainLength, numberOfYoshidaSuzukiSteps, timeScaleParameterBarostat);
          systems[systemId].thermobarostat->numberOfRespaSteps = numberOfRespaSteps;
        }
      }

      if (value.contains("ComputeElasticConstantsFromFluctuations"))
      {
        if (!value["ComputeElasticConstantsFromFluctuations"].is_boolean())
          throw std::runtime_error("[Input reader]: ComputeElasticConstantsFromFluctuations must be a boolean");
        if (value["ComputeElasticConstantsFromFluctuations"].get<bool>())
        {
          std::size_t sampleEvery = 100;
          if (value.contains("ElasticConstantsSampleEvery"))
          {
            if (!value["ElasticConstantsSampleEvery"].is_number_unsigned())
              throw std::runtime_error("[Input reader]: ElasticConstantsSampleEvery must be a positive integer");
            sampleEvery = value["ElasticConstantsSampleEvery"].get<std::size_t>();
          }
          if (sampleEvery == 0)
            throw std::runtime_error("[Input reader]: ElasticConstantsSampleEvery must be positive");
          if (simulationType != SimulationType::MolecularDynamics ||
              systems[systemId].molecularDynamicsEnsemble != MolecularDynamicsEnsemble::NVT)
            throw std::runtime_error(
                "[Input reader]: stress-fluctuation elastic constants currently require NVT molecular dynamics");
          if (systems[systemId].hasExternalField)
            throw std::runtime_error(
                "[Input reader]: stress-fluctuation elastic constants do not support external fields");
          if (systems[systemId].forceField.computePolarization)
            throw std::runtime_error(
                "[Input reader]: stress-fluctuation elastic constants do not support polarization");
          systems[systemId].propertyElasticConstantsFluctuation =
              PropertyElasticConstantsFluctuation(jsonNumberOfBlocks);
          systems[systemId].elasticConstantsSampleEvery = sampleEvery;
        }
      }

      if (value.contains("HybridMCMoveNumberOfSteps") && value["HybridMCMoveNumberOfSteps"].is_number_unsigned())
      {
        systems[systemId].numberOfHybridMCSteps = value["HybridMCMoveNumberOfSteps"].get<std::size_t>();
        if (value.contains("TimeStep") && value["TimeStep"].is_number_float())
        {
          systems[systemId].mc_moves_statistics.setMaxChange(Move::Types::HybridMC, value["TimeStep"].get<double>());
        }
      }

      systemId++;
    }
  }

  // Post-compute
  // ========================================================

  // for (std::size_t i = 0uz; i < systems.size(); ++i)
  //{
  //   systems[i].maxIsothermTerms = 0uz;
  //   if (!systems[i].components.empty())
  //   {
  //     std::vector<Component>::iterator maxIsothermTermsIterator = std::max_element(
  //         systems[i].components.begin(), systems[i].components.end(),
  //         [](Component& lhs, Component& rhs) { return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites; });
  //     systems[i].maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
  //   }
  // }

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
      if (reaction.reactantStoichiometry.size() != systems[i].numerOfAdsorbateComponents() ||
          reaction.productStoichiometry.size() != systems[i].numerOfAdsorbateComponents())
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
      if (systems[i].components[j].lambdaGibbs.computeDUdlambda)
      {
        ++numberOfDUDlambda;
      }
      if (systems[i].components[j].lambdaPairSwap.computeDUdlambda)
      {
        ++numberOfDUDlambda;
      }
      if (systems[i].components[j].lambdaPairSwapCB.computeDUdlambda)
      {
        ++numberOfDUDlambda;
      }
    }
    // serial Rx/CFC reactions consume one dU/dlambda group; parallel Rx/CFC reactions consume two
    // (reactants coupled at 1-lambda and products at lambda need separate accumulators)
    for (const Reaction& reaction : systems[i].reactions.list)
    {
      if (reaction.isSerialRxCFC())
      {
        numberOfDUDlambda += 1;
      }
      else if (reaction.isParallelRxCFC())
      {
        numberOfDUDlambda += 2;
      }
    }
    if (numberOfDUDlambda > maximumNumberOfDUDlambdaGroups)
    {
      throw std::runtime_error(
          std::format("Error [System {}]: too many thermodynamic integrations present "
                      "({} requested, at most {} lambda coordinates can be followed simultaneously)\n",
                      i, numberOfDUDlambda, maximumNumberOfDUDlambdaGroups));
    }
  }

  // for (std::size_t i = 0uz; i < systems.size(); ++i)
  //{
  //   double sum = 0.0;
  //   for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
  //   {
  //     if (systems[i].components[j].type != Component::Type::Cation)
  //     {
  //       sum += systems[i].components[j].molFraction;
  //     }
  //   }
  //   if (std::abs(sum - 1.0) > 1e-15)
  //   {
  //     for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
  //     {
  //       if (systems[i].components[j].type != Component::Type::Cation)
  //       {
  //         systems[i].components[j].molFraction /= sum;
  //       }
  //     }
  //   }
  // }

  // for (std::size_t i = 0uz; i < systems.size(); ++i)
  //{
  //   systems[i].numberOfCarrierGases = 0uz;
  //   systems[i].carrierGasComponent = 0uz;
  //   for (std::size_t j = 0uz; j < systems[i].components.size(); ++j)
  //   {
  //     if (systems[i].components[j].isCarrierGas)
  //     {
  //       systems[i].carrierGasComponent = j;
  //       std::vector<double> values{1.0, 0.0};
  //       const Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
  //       systems[i].components[systems[i].carrierGasComponent].isotherm.add(isotherm);
  //       systems[i].components[systems[i].carrierGasComponent].isotherm.numberOfSites = 1;

  //      systems[i].numberOfCarrierGases++;
  //    }
  //  }

  //  if (simulationType == SimulationType::Breakthrough)
  //  {
  //    if (systems[i].numberOfCarrierGases == 0uz)
  //    {
  //      throw std::runtime_error("Error [Breakthrough]: no carrier gas component present\n");
  //    }
  //    if (systems[i].numberOfCarrierGases > 1)
  //    {
  //      throw std::runtime_error(
  //          "Error [Breakthrough]: multiple carrier gas component present (there can be only one)\n");
  //    }
  //  }
  //}

  for (std::size_t i = 0uz; i < systems.size(); ++i)
  {
    if (systems[i].tmmc.doTMMC)
    {
      if (systems[i].numerOfAdsorbateComponents() > 1)
      {
        throw std::runtime_error("Error: Multiple components for TMMC not yet implemented.\n");
      }

      const std::array unsupportedSystemMoves{
          std::pair{Move::Types::ReactionCBMC, "ReactionCBMCProbability"},
          std::pair{Move::Types::ReactionConventionalCFCMC, "ReactionConventionalCFCMCProbability"},
          std::pair{Move::Types::ReactionConventionalCBCFCMC, "ReactionConventionalCBCFCMCProbability"},
          std::pair{Move::Types::ReactionCFCMC, "ReactionCFCMCProbability"},
          std::pair{Move::Types::ReactionCBCFCMC, "ReactionCBCFCMCProbability"},
          std::pair{Move::Types::ParallelTempering, "ParallelTemperingSwapProbability"}};
      for (const auto& [move, keyword] : unsupportedSystemMoves)
      {
        if (systems[i].mc_moves_probabilities.getProbability(move) != 0.0)
        {
          throw std::runtime_error(
              std::format("Error [TMMC, system {}]: '{}' must be zero; move '{}' is incompatible with the "
                          "current one-dimensional nearest-neighbor TMMC macrostate (only N -> N and N -> N +/- 1 "
                          "transitions are supported).\n",
                          i, keyword, Move::moveNames[std::to_underlying(move)]));
        }
      }

      const std::array unsupportedComponentMoves{
          std::pair{Move::Types::IdentityChangeCBMC, "IdentityChangeProbability"},
          std::pair{Move::Types::GibbsIdentityChangeCBMC, "GibbsIdentityChangeProbability"}};
      for (std::size_t componentId = 0uz; componentId < systems[i].components.size(); ++componentId)
      {
        const Component& component = systems[i].components[componentId];
        for (const auto& [move, keyword] : unsupportedComponentMoves)
        {
          if (component.mc_moves_probabilities.getProbability(move) != 0.0)
          {
            throw std::runtime_error(
                std::format("Error [TMMC, system {}, component {} ('{}')]: '{}' must be zero; move '{}' is "
                            "incompatible with the current one-dimensional nearest-neighbor TMMC macrostate "
                            "(only N -> N and N -> N +/- 1 transitions are supported).\n",
                            i, componentId, component.name, keyword, Move::moveNames[std::to_underlying(move)]));
          }
        }
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
    "NumberOfProductionCycles",
    "NumberOfBlocks",
    "NumberOfPreInitializationCycles",
    "NumberOfInitializationCycles",
    "NumberOfEquilibrationCycles",
    "PrintEvery",
    "MaximumNumberOfMinimizationSteps",
    "MaximumStepLength",
    "MaximumCellStepLength",
    "RMSGradientTolerance",
    "MaxGradientTolerance",
    "MinimizationConvergenceFactor",
    "MinimumEigenvalue",
    "ComputeElasticConstants",
    "ElasticEigenvalueTolerance",
    "ComputeNormalModes",
    "NormalModeMovies",
    "NormalModeMoviePeriods",
    "NormalModeMoviePointsPerPeriod",
    "NormalModeMovieAmplitude",
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
    "ModifiedShiftedForceBeta",
    "VolumeMoveProbability",
    "AnisotropicVolumeMoveProbability",
    "GibbsVolumeMoveProbability",
    "ParallelTemperingSwapProbability",
    "HybridMCProbability",
    "HybridMCMoveNumberOfSteps",
    "ForceBiasTranslationAllProbability",
    "Type",
    "ExternalTemperature",
    "ExternalPressure",
    "ExternalPressureX",
    "ExternalPressureY",
    "ExternalPressureZ",
    "ChemicalPotential",
    "UseChargesFrom",
    "Framework",
    "FrameworkDefinition",
    "FrameworkType",
    "CellType",
    "MonoclinicAngleType",
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
    "ComputeMoleculeProperties",
    "SampleMoleculePropertiesEvery",
    "WriteMoleculePropertiesEvery",
    "NumberOfBinsMoleculeProperties",
    "BondRangeMoleculeProperties",
    "ComputeElasticConstantsFromFluctuations",
    "ElasticConstantsSampleEvery",
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
    "DensityGridBinning",
    "OutputPDBMovie",
    "SampleMovieEvery",
    "RestrictMoviePositionsToBox",
    "ComputeNumberOfMoleculesEvolution",
    "ComputeVolumeEvolution",
    "WriteLammpsData",
    "WriteLammpsDataEvery",
    "Ensemble",
    "ThermostatChainLength",
    "BarostatChainLength",
    "NumberOfRespaSteps",
    "NumberOfYoshidaSuzukiSteps",
    "TimeScaleParameterThermostat",
    "TimeScaleParameterBarostat",
    "TimeStep",
    "MacroStateUseBias",
    "MacroStateMinimumNumberOfMolecules",
    "MacroStateMaximumNumberOfMolecules",
    "RestartFileName",
    "ReactionCBMCProbability",
    "ReactionConventionalCFCMCProbability",
    "ReactionConventionalCBCFCMCProbability",
    "ReactionCFCMCProbability",
    "ReactionCBCFCMCProbability",
    "MaximumReactionLambdaChange",
    "MaximumReactionLambdaChangeProducts",
    "LambdaSwitchPoint",
    "Reactions",
    "ThermochemicalDatabase"};

const std::set<std::string, InputReader::InsensitiveCompare> InputReader::componentOptions = {
    "Name",
    "Type",
    "MoleculeDefinition",
    "TranslationProbability",
    "RandomTranslationProbability",
    "ForceBiasTranslationProbability",
    "RotationProbability",
    "RandomRotationProbability",
    "ReinsertionProbability",
    "PartialReinsertionProbability",
    "IdentityChangeProbability",
    "SwapConventionalProbability",
    "SwapProbability",
    "PairSwapConventionalProbability",
    "PairSwapProbability",
    "CFCMC_PairSwapProbability",
    "CFCMC_CBMC_PairSwapProbability",
    "PairComponent",
    "MaximumPairDistance",
    "CFCMC_SwapProbability",
    "CFCMC_CBMC_SwapProbability",
    "GibbsSwapCBMCProbability",
    "GibbsSwapCFCMCProbability",
    "GibbsSwapCBCFCMCProbability",
    "GibbsConventionalCFCMCProbability",
    "GibbsConventionalCBCFCMCProbability",
    "GibbsIdentityChangeProbability",
    "WidomProbability",
    "CFCMC_WidomProbability",
    "CFCMC_CBMC_WidomProbability",
    "CreateNumberOfMolecules",
    "StartingBead",
    "FugacityCoefficient",
    "IdealGasRosenbluthWeight",
    "MolFraction",
    "IdentityChanges",
    "GibbsIdentityChanges",
    "ThermodynamicIntegration",
    "LambdaBiasFileName",
    "BlockingPockets",
    "LnPartitionFunction"};

const std::set<std::string, InputReader::InsensitiveCompare> InputReader::reactionOptions = {
    "Reactants",        "Products", "Move", "MaximumReactionLambdaChange", "MaximumReactionLambdaChangeProducts",
    "LambdaSwitchPoint"};

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

      if (value.contains("Reactions") && value["Reactions"].is_array())
      {
        for (const auto& reactionItem : value["Reactions"])
        {
          for (auto& [key, _] : reactionItem.items())
          {
            if (!reactionOptions.contains(key))
            {
              throw std::runtime_error(std::format("Error: Unknown reaction input '{}'\n", key));
            }
          }
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
