module;

module input_reader;

import <filesystem>;
import <fstream>;
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
import <print>;

import int3;
import stringutils;
import system;
import atom;
import component;
import simulationbox;
import forcefield;
import double3;
import units;
import sample_movies;
import threadpool;
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


template<class T>
T parse(const std::string& arguments, [[maybe_unused]] const std::string& keyword, [[maybe_unused]] size_t lineNumber)
{
  T value;

  std::string str;
  std::istringstream ss(arguments);

  ss >> value;

  return value;
}

template<typename T>
std::vector<T> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
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
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n", 
                                             keyword, lineNumber));
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
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n", 
                                             keyword, lineNumber));
      }
      return list;
    }

  };

  if (list.empty()) 
  {
    throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}", 
                                         keyword, lineNumber));
  }
  return list;
}

int parseBoolean(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  bool value{};

  std::istringstream ss(arguments);

  if (ss >> std::boolalpha >> value)
  {
    return value;
  };

  std::string str;
  std::istringstream ss2(arguments);
  if (ss2 >> str)
  {
    if (caseInSensStringCompare(str, "yes")) return true;
    if (caseInSensStringCompare(str, "no")) return false;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

template<>
std::vector<bool> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  std::vector<bool> list{};

  std::string str, boolString;
  std::istringstream ss(arguments);

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n", 
                                             keyword, lineNumber));
      }
      return list;
    }
    bool value = false;
    std::istringstream s(str);
    std::istringstream s2(str);
    if (s >> std::boolalpha >>  value)
    {
      list.push_back(value);
    }
    else if (s2 >> boolString)
    {
      if (caseInSensStringCompare(boolString, "yes")) 
      {
        list.push_back(true);
      }
      else if (caseInSensStringCompare(boolString, "no")) 
      {
        list.push_back(false);
      } 
      else
      {
        if (list.empty()) 
        {
          throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n", 
                                               keyword, lineNumber));
        }
        return list;
      }
    }
  };

  if (list.empty()) 
  {
    throw std::runtime_error( std::format("No values could be read for keyword '{}' at line: {}\n", keyword, lineNumber));
  }
  return list;
}


size_t parseInteger(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  size_t value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value)
  {
    return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

double parseDouble(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  double value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value)
  {
      return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

double3 parseDouble3(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  double3 value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value.x >> value.y >> value.z)
  {
    return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

int3 parseInt3(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  int3 value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value.x >> value.y >> value.z)
  {
    return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

std::string parseString(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  std::string value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value)
  {
    return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}


InputReader::InputReader()
{
  const std::string simulationSettingsFileName{"simulation.input"};

  std::cout << "Path: " << std::filesystem::current_path() << std::endl;

  std::filesystem::path pathfile = std::filesystem::path(simulationSettingsFileName);
  if (!std::filesystem::exists(pathfile)) 
  {
    throw std::runtime_error(std::format("Required file '{}' not found\n", simulationSettingsFileName));
  }

  std::ifstream fileInput{ pathfile };
  if (!fileInput) 
  {
    throw std::runtime_error(std::format("File '{}' exists, but error opening file\n", simulationSettingsFileName));
  }

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  std::size_t numberOfSystems{ 0uz };
  size_t lineNumber{ 0uz };

  while (std::getline(fileInput, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);           

      if (caseInSensStringCompare(keyword, "RestartFromBinary"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        restartFromBinary = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "RandomSeed"))
      {
        unsigned long long value = parse<unsigned long long>(arguments, keyword, lineNumber);
        randomSeed = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "SimulationType"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "MonteCarlo"))
          {
            simulationType = SimulationType::MonteCarlo;
            continue;
          }
          if (caseInSensStringCompare(str, "MonteCarloTransitionMatrix"))
          {
            simulationType = SimulationType::MonteCarloTransitionMatrix;
            continue;
          }
          if (caseInSensStringCompare(str, "MolecularDynamics"))
          {
            simulationType = SimulationType::MolecularDynamics;
            continue;
          }
          if (caseInSensStringCompare(str, "Breakthrough"))
          {
            simulationType = SimulationType::Breakthrough;
            continue;
          }
          if (caseInSensStringCompare(str, "Minimization"))
          {
            simulationType = SimulationType::Minimization;
            continue;
          }
          if (caseInSensStringCompare(str, "MixturePrediction"))
          {
            simulationType = SimulationType::MixturePrediction;
            continue;
          }
          if (caseInSensStringCompare(str, "Fitting"))
          {
            simulationType = SimulationType::Fitting;
            continue;
          }
          if (caseInSensStringCompare(str, "Test"))
          {
            simulationType = SimulationType::Test;
            continue;
          }
        };
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfCycles")))
      {
        numberOfCycles = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfInitializationCycles")))
      {
        numberOfInitializationCycles = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfEquilibrationCycles")))
      {
        numberOfEquilibrationCycles = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("PrintEvery")))
      {
        printEvery = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("WriteBinaryRestartEvery")))
      {
        writeBinaryRestartEvery = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("RescaleWangLandauEvery")))
      {
        rescaleWangLandauEvery = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("OptimizeMCMovesEvery")))
      {
        optimizeMCMovesEvery = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfBlocks")))
      {
        numberOfBlocks = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfLambdaBins")))
      {
        numberOfLambdaBins = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("NumberOfThreads")))
      {
        numberOfThreads = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("ThreadingType")))
      {
        std::string str;
        std::istringstream ss2(arguments);
        if (ss2 >> str)
        {
          if (caseInSensStringCompare(str, "Serial")) 
          {
            threadingType = ThreadPool::ThreadingType::Serial;
            continue;
          }
          if (caseInSensStringCompare(str, "ThreadPool")) 
          {
            threadingType = ThreadPool::ThreadingType::ThreadPool;
            continue;
          }
          if (caseInSensStringCompare(str, "OpenMP")) 
          {
            threadingType = ThreadPool::ThreadingType::OpenMP;
            continue;
          }
          if (caseInSensStringCompare(str, "GPU-Offload"))
          {
            threadingType = ThreadPool::ThreadingType::GPU_Offload;
            continue;
          }
        };
      }

      if (caseInSensStringCompare(keyword, "ExternalField"))
      {
        requireExistingSystem(keyword, lineNumber);
        bool value = parseBoolean(arguments, keyword, lineNumber);
        systems.back().hasExternalField = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        ForceField currentForceField = ForceField(numberOfSystems);
        numberOfSystems += 1;
        systems.emplace_back(numberOfSystems - 1, currentForceField, std::vector<Component>{}, 
                             std::vector<size_t>{}, numberOfBlocks);
        continue;
      }

      if (caseInSensStringCompare(keyword, "BoxLengths"))
      {
        requireExistingSystem(keyword, lineNumber);
        double3 value = parseDouble3(arguments, keyword, lineNumber);
        systems.back().simulationBox.setBoxLengths(value);
        continue;
      }

      if (caseInSensStringCompare(keyword, "BoxAngles"))
      {
        requireExistingSystem(keyword, lineNumber);
        double3 value = parseDouble3(arguments, keyword, lineNumber);
        systems.back().simulationBox.setBoxAngles((std::numbers::pi / 180.0) * value);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        ForceField currentForceField = ForceField(numberOfSystems);
        numberOfSystems += 1;
        systems.emplace_back(numberOfSystems - 1, currentForceField, std::vector<Component>{}, 
                             std::vector<size_t>{}, numberOfBlocks);
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("FrameworkName")))
      {
        std::istringstream ss(arguments);
        std::string frameworkName;
        ss >> frameworkName;

        switch(simulationType)
        {
          case SimulationType::MonteCarlo:
          case SimulationType::MonteCarloTransitionMatrix:
          case SimulationType::MolecularDynamics:
          case SimulationType::Minimization:
          case SimulationType::Test:
          {

            systems.back().addComponent(Component(Component::Type::Framework, systems.back().components.size(), 
                                        frameworkName, frameworkName, numberOfBlocks, numberOfLambdaBins));
            continue;
            break;
          }
          case SimulationType::Breakthrough:
          case SimulationType::MixturePrediction:
          case SimulationType::Fitting:
          {
            systems.back().addComponent(Component(Component::Type::Framework, systems.back().components.size(), 
                                        frameworkName, std::nullopt, numberOfBlocks, numberOfLambdaBins));
            continue;
            break;
          }
        }
      }


      if (caseInSensStringCompare(keyword, "NumberOfUnitCells"))
      {
          requireExistingSystemAndComponent(keyword, lineNumber);
          int3 value = parseInt3(arguments, keyword, lineNumber);
          systems.back().components.back().numberOfUnitCells = value;
          continue;
      }

      if (caseInSensStringCompare(keyword, "ExternalTemperature"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().temperature = value;
        systems.back().beta = 1.0 / (Units::KB * value);
        continue;
      }

      if (caseInSensStringCompare(keyword, "ExternalPressure"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().input_pressure = value;
        systems.back().pressure = value / Units::PressureConversionFactor;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ChargeMethod"))
      {
        requireExistingSystem(keyword, lineNumber);
        std::string value = parseString(arguments, keyword, lineNumber);
        if (caseInSensStringCompare(value, "None"))
        {
          systems.back().forceField.noCharges = true;
          continue;
        }
        if (caseInSensStringCompare(value, "Ewald")) 
        {
          systems.back().forceField.noCharges = false;
          continue;
        }
      }

      if (caseInSensStringCompare(keyword, "EwaldPrecision"))
      {
          requireExistingSystem(keyword, lineNumber);
          double value = parseDouble(arguments, keyword, lineNumber);
          systems.back().forceField.automaticEwald = true;
          systems.back().forceField.EwaldPrecision = value;
          continue;
      }

      if (caseInSensStringCompare(keyword, "EwaldParameters"))
      {
          requireExistingSystem(keyword, lineNumber);
          systems.back().forceField.automaticEwald = false;

          std::istringstream iss1(arguments);
          std::string alpha, kvectors;
          iss1 >> alpha;
          std::getline(iss1, kvectors);

          double value = parseDouble(alpha, keyword, lineNumber);
          systems.back().forceField.EwaldAlpha = value;
          int3 values = parseInt3(kvectors, keyword, lineNumber);
          systems.back().forceField.numberOfWaveVectors = values;
          continue;
      }

      if (caseInSensStringCompare(keyword, "OmitEwaldFourier"))
      {
        requireExistingSystem(keyword, lineNumber);
        bool value = parseBoolean(arguments, keyword, lineNumber);
        systems.back().forceField.omitEwaldFourier = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "VolumeMoveProbability"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().mc_moves_probabilities.probabilityVolumeMove = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsVolumeMoveProbability"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().mc_moves_probabilities.probabilityGibbsVolumeMove = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "UseBiasOnMacrostate"))
      {
        requireExistingSystem(keyword, lineNumber);
        bool value = parseBoolean(arguments, keyword, lineNumber);
        systems.back().tmmc.useBias = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "TMMCMin"))
      {
        requireExistingSystem(keyword, lineNumber);
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        systems.back().tmmc.minMacrostate = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "TMMCMax"))
      {
        requireExistingSystem(keyword, lineNumber);
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        systems.back().tmmc.maxMacrostate = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "Reaction"))
      {
        requireExistingSystem(keyword, lineNumber);
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        systems.back().reactions.list.emplace_back(Reaction(systems.back().reactions.list.size(), 
                                                     std::vector(values.begin(), values.begin() + values.size() / 2),
                                                     std::vector(values.begin() + values.size() / 2, values.end())));
        continue;
      }

      if (caseInSensStringCompare(keyword, "Movies"))
      {
        //requireExistingSystem(keyword, lineNumber);
        //bool value = parseBoolean(arguments, keyword, lineNumber);
        //systems.back().sampleMovie.sample = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "WriteMoviesEvery"))
      {
        //requireExistingSystem(keyword, lineNumber);
        //size_t value = parseInteger(arguments, keyword, lineNumber);
        //systems.back().sampleMovie.writeEvery = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "PressureStart"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().pressure_range.pressureStart = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureEnd"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().pressure_range.pressureEnd = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "NumberOfPressurePoints"))
      {
        requireExistingSystem(keyword, lineNumber);
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        systems.back().pressure_range.numberOfPoints = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureScale"))
      {
        requireExistingSystem(keyword, lineNumber);
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "Log"))
          {
            systems.back().pressure_range.scale = PressureRange::Scale::Log;
            continue;
          }
          if (caseInSensStringCompare(str, "Linear"))
          {
            systems.back().pressure_range.scale = PressureRange::Scale::Linear;
            continue;
          }
        }
      }

      if (caseInSensStringCompare(keyword, "ColumnVoidFraction"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnVoidFraction = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ParticleDensity"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnParticleDensity = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "TotalPressure"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnTotalPressure = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureGradient"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnPressureGradient = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnEntranceVelocity"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnEntranceVelocity = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "TimeStep"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnTimeStep = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "NumberOfTimeSteps"))
      {
        requireExistingSystem(keyword, lineNumber);
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "auto"))
          {
            systems.back().columnAutoNumberOfTimeSteps = true;
          }
          else
          {
            size_t value = parse<size_t>(arguments, keyword, lineNumber);
            systems.back().columnNumberOfTimeSteps = value;
            systems.back().columnAutoNumberOfTimeSteps = false;
          }
          continue;
        }
      }
      if (caseInSensStringCompare(keyword, "WriteEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->writeEvery = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnLength"))
      {
        requireExistingSystem(keyword, lineNumber);
        double value = parseDouble(arguments, keyword, lineNumber);
        systems.back().columnLength = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "NumberOfGridPoints"))
      {
        requireExistingSystem(keyword, lineNumber);
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        systems.back().columnNumberOfGridPoints = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "MixturePredictionMethod"))
      {
        requireExistingSystem(keyword, lineNumber);
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "IAST")) 
          {
            systems.back().mixturePredictionMethod = MultiSiteIsotherm::PredictionMethod::IAST;
            continue;
          }
          if (caseInSensStringCompare(str, "SIAST")) 
          {
            systems.back().mixturePredictionMethod = MultiSiteIsotherm::PredictionMethod::SIAST;
            continue;
          }
          if (caseInSensStringCompare(str, "EI")) 
          {
            systems.back().mixturePredictionMethod = MultiSiteIsotherm::PredictionMethod::EI;
            continue;
          }
          if (caseInSensStringCompare(str, "SEI")) 
          {
            systems.back().mixturePredictionMethod = MultiSiteIsotherm::PredictionMethod::SEI;
            continue;
          }
        };
      }

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        requireExistingSystem(keyword, lineNumber);
        
        std::istringstream ss(arguments);
        std::string c, moleculeNameKeyword,remainder;
        ss >> c >> moleculeNameKeyword;
        std::getline(ss, remainder);

        std::vector<std::string> values = parseListOfSystemValues<std::string>(remainder, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          switch(simulationType)
          {
            case SimulationType::MonteCarlo:
            case SimulationType::MonteCarloTransitionMatrix:
            case SimulationType::MolecularDynamics:
            case SimulationType::Minimization:
            case SimulationType::Test:
              systems[i].addComponent(Component(Component::Type::Adsorbate, systems[i].components.size(),
                    values[i], values[i], numberOfBlocks, numberOfLambdaBins));
              break;
            case SimulationType::Breakthrough:
            case SimulationType::MixturePrediction:
            case SimulationType::Fitting:
              systems[i].addComponent(Component(Component::Type::Adsorbate, systems[i].components.size(),
                    values[i], std::nullopt, numberOfBlocks, numberOfLambdaBins));
              break;
          }
        }
        continue;
      }

      
      if (caseInSensStringCompare(keyword, "CreateNumberOfMolecules"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);

        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().initialNumberOfMolecules = values[i];
          systems[i].initialNumberOfMolecules[systems[i].components.size() - 1] = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "MolFraction"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().molFraction = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "ThermodynamicIntegration"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<bool> values = parseListOfSystemValues<bool>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().lambdaGC.computeDUdlambda = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "LnPartitionFunction"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().lnPartitionFunction = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "FileName"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          for (size_t i = 0uz; i < systems.size(); ++i)
          {
            systems[i].components.back().filename = str;
          }
          continue;
        }
      }
      if (caseInSensStringCompare(keyword, "CarrierGas"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<bool> values = parseListOfSystemValues<bool>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isCarrierGas = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "MassTransferCoefficient"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().massTransferCoefficient = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "AxialDispersionCoefficient"))
      {
        requireExistingSystem(keyword, lineNumber);

        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().axialDispersionCoefficient = values[i];
        }
        continue;
      }


      if (caseInSensStringCompare(keyword, "TranslationProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityTranslationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RandomTranslationProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityRandomTranslationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RotationProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityRotationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RandomRotationProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityRandomRotationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "ReinsertionProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityReinsertionMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "SwapProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilitySwapMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsSwapProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityGibbsSwapMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsCFCMCSwapProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CFCMC_SwapProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilitySwapMove_CFCMC = values[i];
        }
        continue;
      }
      
      if (caseInSensStringCompare(keyword, "CFCMC_CBMC_SwapProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "WidomProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityWidomMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CFCMCWidomProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityWidomMove_CFCMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CBCFCMCWidomProbability"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "MassTransferCoefficient"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().massTransferCoefficient = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "AxialDispersionCoefficient"))
      {
        requireExistingSystemAndComponent(keyword, lineNumber);
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().axialDispersionCoefficient = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "NumberOfIsothermSites"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.numberOfSites = value;
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Langmuir"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2uz)
        {
          throw std::runtime_error("Error: Langmuir requires two parameters\n");
        }
        values.resize(2uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Anti-Langmuir"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2uz)
        {
          throw std::runtime_error("Error: Anti-Langmuir requires two parameters\n");
        }
        values.resize(2uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Anti_Langmuir, values, 2);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "BET"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: BET requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::BET, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Henry"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 1uz)
        {
          throw std::runtime_error("Error: Henry requires one parameter\n");
        }
        values.resize(1uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Henry, values, 1);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Freundlich"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2uz)
        {
          throw std::runtime_error("Error: Freundlich requires two parameters\n");
        }
        values.resize(2uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Freundlich, values, 2);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Sips"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Sips requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Sips, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Langmuir-Freundlich"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Langmuir requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir_Freundlich, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Redlich-Peterson"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Redlich-Peterson requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Redlich_Peterson, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Toth"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Toth requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Toth, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Unilan"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Unilan requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Unilan, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "O'Brien&Myers"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: O'Brien&Myers requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::OBrien_Myers, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Quadratic"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Quadratic requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Quadratic, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "Temkin"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3uz)
        {
          throw std::runtime_error("Error: Temkin requires three parameters\n");
        }
        values.resize(3uz);
        const Isotherm isotherm = Isotherm(Isotherm::Type::Temkin, values, 3);
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().isotherm.add(isotherm);
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnPressure"))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().columnPressure = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnLoading"))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().columnLoading = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "ColumnError"))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components.back().columnError = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "computeConventionalRDF"))
      {
        std::vector<bool> values = parseListOfSystemValues<bool>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].computeConventionalRadialDistributionFunction = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "writeConventionalRDFEvery"))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].writeConventionalRadialDistributionFunctionEvery = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "conventionalRDFistogramSize"))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].conventionalRadialDistributionFunctionHistogramSize = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "conventionalRDFRange"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].conventionalRadialDistributionFunctionRange = values[i];
        }
        continue;
      }

      if (!(keyword.starts_with("//") || (keyword.starts_with("#"))))
      {
        throw std::runtime_error(std::format("Error [Input]: unrecognized keyword '{}' at line: {}\n", 
                                             keyword, lineNumber));
      }

    }

  }

  // Post-initialize
  // ========================================================
  
  // read and initialize the components from file
  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    systems[i].initializeComponents();
  }

  for(System &system: systems)
  {
    system.conventionalRadialDistributionFunction = 
      PropertyConventionalRadialDistributionFunction(numberOfBlocks, system.forceField.pseudoAtoms.size(), 
                                                     system.conventionalRadialDistributionFunctionHistogramSize,
                                                     system.conventionalRadialDistributionFunctionRange);
  }
  

  // Post-compute
  // ========================================================
  

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    systems[i].maxIsothermTerms = 0uz;
    if(!systems[i].components.empty())
    {
      std::vector<Component>::iterator maxIsothermTermsIterator = 
            std::max_element(systems[i].components.begin(), systems[i].components.end(),
              [] (Component& lhs, Component& rhs) {
                  return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites;
              });
      systems[i].maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
    }
  }

  if(simulationType == SimulationType::MonteCarloTransitionMatrix)
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
        throw std::runtime_error(std::format("Error [Reaction {}]: mismatch Stoichiometry ({} given not equal" 
                                             "to twice the number of components {})\n", 
                                             reactionId, reaction.productStoichiometry.size() + 
                                             reaction.reactantStoichiometry.size(), 
                                             2uz * systems[i].numerOfAdsorbateComponents()));
      }
    
      ++reactionId;
    }
  }

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    size_t numberOfDUDlambda{ 0uz };
    for (size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if (systems[i].components[j].lambdaGC.computeDUdlambda) 
      {
        ++numberOfDUDlambda;
      }
    }
    if(numberOfDUDlambda > 1)
    {
      throw std::runtime_error(std::format("Error [System {}]: multiple thermodynamic integrations present " 
                                           "(there can be only one)\n", i));
    }
  }

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    double sum = 0.0;
    for(size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if(systems[i].components[j].type != Component::Type::Framework)
      {
        sum += systems[i].components[j].molFraction;
      }
    }
    if(std::abs(sum - 1.0) > 1e-15)
    {
      std::cout << "Normalizing: Gas-phase molfractions did not sum exactly to unity!\n\n";
      for(size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if(systems[i].components[j].type != Component::Type::Framework)
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
    for(size_t j = 0uz; j < systems[i].components.size(); ++j)
    {
      if(systems[i].components[j].type != Component::Type::Framework)
      {
        if(systems[i].components[j].isCarrierGas)
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

    if(simulationType == SimulationType::Breakthrough)
    {
      if(systems[i].numberOfCarrierGases == 0uz)
      {
        throw std::runtime_error("Error [Breakthrough]: no carrier gas component present\n");
      }
      if(systems[i].numberOfCarrierGases > 1)
      {
        throw std::runtime_error("Error [Breakthrough]: multiple carrier gas component present (there can be only one)\n");
      }
    }
  }

  for (size_t i = 0uz; i < systems.size(); ++i)
  {
    if(systems[i].tmmc.doTMMC)
    {
      if(systems[i].numerOfAdsorbateComponents() > 1)
      {
        throw std::runtime_error("Error: Multiple components for TMMC not yet implemented.\n");
      }

      // check initial number of molecules is in the range of the TMMC macrostates
      for(size_t j = 0uz; j < systems[i].components.size(); ++j)
      {
        if(systems[i].components[j].type == Component::Type::Adsorbate)
        {
          size_t numberOfMolecules = systems[i].initialNumberOfMolecules[j];
          if(numberOfMolecules < systems[i].tmmc.minMacrostate || numberOfMolecules > systems[i].tmmc.maxMacrostate)
          {
            throw std::runtime_error(std::format("Error: Molecules created ({}) need to fit into the TMMC macrostate "
                                                 "range ({}-{})\n", numberOfMolecules, systems[i].tmmc.minMacrostate,
                                                 systems[i].tmmc.maxMacrostate));
          }
        }
      }
    }
  }
}

void InputReader::requireExistingSystem(const std::string& keyword, size_t lineNumber)
{
  if (systems.empty()) 
  {
    throw std::runtime_error(std::format("No system (Framework or Box) defined yet at keyword '{}' at line: {}\n", 
                                         keyword, lineNumber));
  }
}

void InputReader::requireExistingSystemAndComponent(const std:: string &keyword, size_t lineNumber)
{
  if (systems.empty()) 
  {
    throw std::runtime_error(
       std::format("No system (Framework or Box) defined yet at keyword '{}' at line: {}\n", keyword, lineNumber));
  }
  if (systems[0uz].components.empty()) 
  {
    throw std::runtime_error(
      std::format("No component defined yet at keyword '{}' at line: {}\n", keyword, lineNumber));
  }
}
