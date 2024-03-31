module;

module input_reader;

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
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#else
  import print;
#endif


import int3;
import stringutils;
import system;
import atom;
import framework;
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
import property_density_grid;



InputReader::InputReader(const std::string inputFile):
  inputStream(inputFile),
  inputString((std::istreambuf_iterator<char>(inputStream)), std::istreambuf_iterator<char>()),
  simulationType(scanSimulationType()),
  numberOfBlocks(scanNumberOfBlocks()),
  totalNumberOfSystems(scanNumberOfSystems()),
  totalNumberOfComponents(scanNumberOfComponents()),
  forceFields(scanForceFields()),
  inputDataSystem(totalNumberOfSystems),
  systems(totalNumberOfSystems)
{
  // pre-scan all the input
  scanForceFieldParameters(forceFields);

  std::vector<std::optional<SimulationBox>> boxes = scanBoxes();
  std::vector<std::vector<Framework>> frameworkComponents = scanFrameworkComponents();

  std::vector<std::vector<Component>> components = scanComponents();
  std::vector<std::vector<size_t>> initialNumberOfMolecules = scanInitialNumberOfMolecules();

  std::vector<double> systemTemperatures = scanSystemTemperatures();
  std::vector<std::optional<double>> systemPressures = scanSystemPressures();

  // construct the systems
  for(size_t i = 0; i < totalNumberOfSystems; ++i)
  {
    systems[i] = System(i, boxes[i], systemTemperatures[i], systemPressures[i], forceFields[i], 
                        frameworkComponents[i], components[i], initialNumberOfMolecules[i], numberOfBlocks);
  }

  // read in settings
  scanGeneralSettings();
  scanSystemSettings();
  scanComponentSettings();

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  size_t lineNumber{ 0uz };

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);           

      if(allKeywords.contains(keyword))
      {
        continue;
      }

      if (!(keyword.starts_with("//") || (keyword.starts_with("#"))))
      {
        throw std::runtime_error(std::format("Error [Input]: unrecognized keyword '{}' at line: {}\n", 
                                             keyword, lineNumber));
      }

    }

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

InputReader::SimulationType InputReader::scanSimulationType()
{
  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  InputReader::SimulationType scannedSimulationType{SimulationType::MonteCarlo};
  while (std::getline(s, line))
  {
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, "SimulationType"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "MonteCarlo"))
          {
            scannedSimulationType = SimulationType::MonteCarlo;
            continue;
          }
          if (caseInSensStringCompare(str, "MonteCarloTransitionMatrix"))
          {
            scannedSimulationType = SimulationType::MonteCarloTransitionMatrix;
            continue;
          }
          if (caseInSensStringCompare(str, "MolecularDynamics"))
          {
            scannedSimulationType = SimulationType::MolecularDynamics;
            continue;
          }
          if (caseInSensStringCompare(str, "Breakthrough"))
          {
            scannedSimulationType = SimulationType::Breakthrough;
            continue;
          }
          if (caseInSensStringCompare(str, "Minimization"))
          {
            scannedSimulationType = SimulationType::Minimization;
            continue;
          }
          if (caseInSensStringCompare(str, "MixturePrediction"))
          {
            scannedSimulationType = SimulationType::MixturePrediction;
            continue;
          }
          if (caseInSensStringCompare(str, "Fitting"))
          {
            scannedSimulationType = SimulationType::Fitting;
            continue;
          }
          if (caseInSensStringCompare(str, "Test"))
          {
            scannedSimulationType = SimulationType::Test;
            continue;
          }
        }
      }
    }
  }
  return scannedSimulationType;
}


size_t InputReader::scanNumberOfSystems()
{
  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  size_t scannedNumberOfSystems{};
  while (std::getline(s, line))
  {
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")) ||
          caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems++;
        continue;
      }
    }
  }
  return scannedNumberOfSystems;
}

size_t InputReader::scanNumberOfComponents()
{
  std::stringstream s(inputString);


  std::string line{};
  std::string keyword{};
  std::string arguments{};

  size_t scannedNumberOfComponents{};
  while (std::getline(s, line))
  {
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        scannedNumberOfComponents++;
        continue;
      }
    }
  }
  return scannedNumberOfComponents;
}


size_t InputReader::scanNumberOfBlocks()
{
  std::stringstream s(inputString);

  size_t lineNumber{};
  std::string line{};
  std::string keyword{};
  std::string arguments{};

  size_t scannedNumberOfBlocks{ 5 };
  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("NumberOfBlocks")))
      {
        scannedNumberOfBlocks = parse<size_t>(arguments, keyword, lineNumber);
        continue;
      }
    }
  }
  return scannedNumberOfBlocks;
}

// pre-read all the force fields
std::vector<ForceField> InputReader::scanForceFields()
{
  std::vector<ForceField> scannedForceFields(totalNumberOfSystems);

  size_t scannedNumberOfSystems{};

  // read global forcefield, and pre-set it to all the systems
  ForceField forcefield(0);
  for(size_t i = 0; i != forceFields.size(); ++i)
  {
    scannedForceFields[i] = forcefield;
  }

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        scannedNumberOfSystems++;
        continue;
      }
      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems++;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("ForceField")))
      {
        // read the forcefield defined for this particular system
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          // overwrite default forcefield
          scannedForceFields[scannedNumberOfSystems - 1] = ForceField(scannedNumberOfSystems - 1);
        }
      }
    }
  }
  return scannedForceFields;
}

void InputReader::scanForceFieldParameters(std::vector<ForceField> &scannedForceFields)
{
  size_t scannedNumberOfSystems{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        scannedNumberOfSystems++;
        continue;
      }
      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems++;
        continue;
      }


      if (caseInSensStringCompare(keyword, "ChargeMethod"))
      {
        std::string value = parseString(arguments, keyword, lineNumber);

        if (caseInSensStringCompare(value, "None"))
        {
          scannedForceFields[scannedNumberOfSystems - 1].noCharges = true;
          continue;
        }

        if (caseInSensStringCompare(value, "Ewald")) 
        {
          scannedForceFields[scannedNumberOfSystems - 1].noCharges = false;
          continue;
        }
      }

      if (caseInSensStringCompare(keyword, "EwaldPrecision"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        scannedForceFields[scannedNumberOfSystems - 1].automaticEwald = true;
        scannedForceFields[scannedNumberOfSystems - 1].EwaldPrecision = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "EwaldParameters"))
      {
        scannedForceFields[scannedNumberOfSystems - 1].automaticEwald = false;

        std::istringstream iss1(arguments);
        std::string alpha, kvectors;
        iss1 >> alpha;
        std::getline(iss1, kvectors);

        double value = parseDouble(alpha, keyword, lineNumber);
        scannedForceFields[scannedNumberOfSystems - 1].EwaldAlpha = value;
        int3 values = parseInt3(kvectors, keyword, lineNumber);
        scannedForceFields[scannedNumberOfSystems - 1].numberOfWaveVectors = values;
        continue;
      }

      if (caseInSensStringCompare(keyword, "OmitEwaldFourier"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        scannedForceFields[scannedNumberOfSystems - 1].omitEwaldFourier = value;
        continue;
      }
    }
  }
}

std::vector<std::vector<size_t>> InputReader::scanInitialNumberOfMolecules()
{
  std::vector<std::vector<size_t>> scannedInitialNumberOfMolecules(totalNumberOfSystems, std::vector<size_t>(totalNumberOfComponents, 0uz));

  size_t scannedNumberOfComponents{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        scannedNumberOfComponents++;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("CreateNumberOfMolecules")))
      {
        std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
        values.resize(totalNumberOfSystems, values.back());
        for (size_t i = 0uz; i < totalNumberOfSystems; ++i)
        {
          scannedInitialNumberOfMolecules[i][scannedNumberOfComponents - 1] = values[i];
        }
        continue;
      }
    }
  }

  return scannedInitialNumberOfMolecules;
}


std::vector<std::optional<SimulationBox>> InputReader::scanBoxes()
{
  std::vector<std::optional<SimulationBox>> scannedBoxes(totalNumberOfSystems);

  std::vector<double3> scannedBoxLengths(totalNumberOfSystems, double3(30.0, 30.0, 30.0));
  std::vector<double3> scannedBoxAngles(totalNumberOfSystems, double3(90.0, 90.0, 90.0));
  std::vector<bool> isBox(totalNumberOfSystems, false);

  size_t scannedNumberOfSystems{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        isBox[scannedNumberOfSystems] = true;
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        isBox[scannedNumberOfSystems] = false;
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, "BoxLengths"))
      {
        double3 value = parseDouble3(arguments, keyword, lineNumber);
        scannedBoxLengths[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "BoxAngles"))
      {
        requireExistingSystem(keyword, lineNumber);
        double3 value = parseDouble3(arguments, keyword, lineNumber);
        scannedBoxAngles[scannedNumberOfSystems - 1] = value;
        continue;
      }

    }
  }

  for(size_t i = 0; i != totalNumberOfSystems; ++i)
  {
    if(isBox[i])
    {
      double3 lengths = scannedBoxLengths[i];
      double3 angles = (std::numbers::pi / 180.0) * scannedBoxAngles[i];
      scannedBoxes[i] = SimulationBox(lengths.x, lengths.y, lengths.z, angles.x, angles.y, angles.z);
    }
  }

  return scannedBoxes;
}

std::vector<std::vector<Framework>> InputReader::scanFrameworkComponents()
{
  std::vector<std::vector<Framework>> scannedFrameworkComponents(totalNumberOfSystems, std::vector<Framework>());
  std::vector<std::vector<std::string>> scannedFrameworkNames(totalNumberOfSystems, std::vector<std::string>());
  std::vector<std::vector<int3>> scannedNumberOfUnitCells(totalNumberOfSystems, std::vector<int3>());
  std::vector<size_t> scannedNumberOfFrameworkComponents(totalNumberOfSystems, 0uz);

  size_t scannedNumberOfSystems{};
  size_t scannedNumberOfFrameworks{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        scannedNumberOfSystems += 1uz;
        scannedNumberOfFrameworkComponents[scannedNumberOfSystems - 1] = 0uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems += 1uz;
        scannedNumberOfFrameworkComponents[scannedNumberOfSystems - 1] = 0uz;
        scannedNumberOfFrameworks = 0uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("FrameworkName")))
      {
        std::istringstream ss(arguments);
        std::string frameworkName;
        ss >> frameworkName;

        scannedFrameworkNames[scannedNumberOfSystems - 1].push_back(frameworkName);
        scannedNumberOfUnitCells[scannedNumberOfSystems - 1].push_back(int3(1,1,1));

        scannedNumberOfFrameworkComponents[scannedNumberOfSystems - 1]++;

        scannedNumberOfFrameworks++;
      }

      if (caseInSensStringCompare(keyword, "NumberOfUnitCells"))
      {
        int3 value = parseInt3(arguments, keyword, lineNumber);

        // check needed when 'NumberOfUnitCells' is defined before 'FrameworkName'
        if(scannedNumberOfUnitCells[scannedNumberOfSystems - 1].size() == scannedNumberOfFrameworks)
        {
          scannedNumberOfUnitCells[scannedNumberOfSystems - 1][scannedNumberOfFrameworks - 1] = value;
        }
        else
        {
          scannedNumberOfUnitCells[scannedNumberOfSystems - 1].push_back(value);
        }
        continue;
      }

    }
  }

  for(size_t i = 0; i != totalNumberOfSystems; ++i)
  {
    for(size_t j = 0; j != scannedNumberOfFrameworkComponents[i]; ++j)
    {
      Framework f = Framework(j, forceFields[i], scannedFrameworkNames[i][j],
                                                   scannedFrameworkNames[i][j], scannedNumberOfUnitCells[i][j]);
      scannedFrameworkComponents[i].push_back(f);
    }
  }

  return scannedFrameworkComponents;
}


std::vector<std::vector<Component>> InputReader::scanComponents()
{
  std::vector<std::vector<Component>> components(totalNumberOfSystems, std::vector<Component>(totalNumberOfComponents));
  std::vector<std::string> componentNames(totalNumberOfComponents);
  std::vector<size_t> numberOfLambdaBins(totalNumberOfComponents, 21uz);

  size_t scannedNumberOfComponents{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        scannedNumberOfComponents++;
      }
      if (caseInSensStringCompare(keyword, std::string("MoleculeName")))
      {
        std::string name = parseString(arguments, keyword, lineNumber);
        componentNames[scannedNumberOfComponents - 1] = name;
      }
      if (caseInSensStringCompare(keyword, std::string("NumberOfLambdaBins")))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        numberOfLambdaBins[scannedNumberOfComponents - 1] = value;
      }
    }
  }

  for(size_t i = 0; i != totalNumberOfSystems; ++i)
  {
    for(size_t j = 0; j != totalNumberOfComponents; ++j)
    {
      components[i][j] = Component(Component::Type::Adsorbate, j, forceFields[i], componentNames[j],
                                   componentNames[j], numberOfBlocks, numberOfLambdaBins[j]);
    }
  }

  return components;
}

std::vector<double> InputReader::scanSystemTemperatures()
{
  std::vector<double> scannedSystemTemperatures(totalNumberOfSystems, 300.0);

  size_t scannedNumberOfSystems{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")) ||
          caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ExternalTemperature"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        scannedSystemTemperatures[scannedNumberOfSystems - 1] = value;
        continue;
      }
    }
  }

  return scannedSystemTemperatures;
}

std::vector<std::optional<double>> InputReader::scanSystemPressures()
{
  std::vector<std::optional<double>> scannedSystemPressures(totalNumberOfSystems);

  size_t scannedNumberOfSystems{};
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);

      if (caseInSensStringCompare(keyword, std::string("Box")) ||
          caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ExternalPressure"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        scannedSystemPressures[scannedNumberOfSystems - 1] = value;
        continue;
      }
    }
  }

  return scannedSystemPressures;
}


void InputReader::scanGeneralSettings()
{
  size_t lineNumber{};

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  while (std::getline(s, line))
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
        continue;
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


      if (caseInSensStringCompare(keyword, std::string("NumberOfLambdaBins")))
      {
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
        }
      }
    }
  }
}

void InputReader::scanSystemSettings()
{
  std::vector<bool> computeConventionalRadialDistributionFunction(totalNumberOfSystems);
  std::vector<size_t> sampleConventionalRadialDistributionFunctionEvery(totalNumberOfSystems);
  std::vector<size_t> writeConventionalRadialDistributionFunctionEvery(totalNumberOfSystems);
  std::vector<size_t> conventionalRadialDistributionFunctionHistogramSize(totalNumberOfSystems);
  std::vector<double> conventionalRadialDistributionFunctionRange(totalNumberOfSystems);

  std::vector<bool> computeRadialDistributionFunction(totalNumberOfSystems);
  std::vector<size_t> sampleRadialDistributionFunctionEvery(totalNumberOfSystems);
  std::vector<size_t> writeRadialDistributionFunctionEvery(totalNumberOfSystems);
  std::vector<size_t> radialDistributionFunctionHistogramSize(totalNumberOfSystems);
  std::vector<double> radialDistributionFunctionRange(totalNumberOfSystems);

  std::vector<bool> computeDensityGrid(totalNumberOfSystems);
  std::vector<size_t> sampleDensityGridEvery(totalNumberOfSystems);
  std::vector<size_t> writeDensityGridEvery(totalNumberOfSystems);
  std::vector<int3> densityGridSize(totalNumberOfSystems);

  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  std::size_t scannedNumberOfSystems{ 0uz };
  size_t lineNumber{ 0uz };

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);           

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ComputeConventionalRDF"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        computeConventionalRadialDistributionFunction[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "SampleConventionalRDFEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        sampleConventionalRadialDistributionFunctionEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "SampleConventionalRDFEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        sampleConventionalRadialDistributionFunctionEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "WriteConventionalRDFEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        writeConventionalRadialDistributionFunctionEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ConventionalRDFistogramSize"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        conventionalRadialDistributionFunctionHistogramSize[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ConventionalRDFRange"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        conventionalRadialDistributionFunctionRange[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ComputeRDF"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        computeRadialDistributionFunction[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "WriteRDFEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        writeRadialDistributionFunctionEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "RDFistogramSize"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        radialDistributionFunctionHistogramSize[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "RDFRange"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        radialDistributionFunctionRange[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ComputeDensityGrid"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        computeDensityGrid[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "SampleDensityGridEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        sampleDensityGridEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "WriteDensityGridEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        writeDensityGridEvery[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "DensityGridSize"))
      {
        int3 value = parseInt3(arguments, keyword, lineNumber);
        densityGridSize[scannedNumberOfSystems - 1] = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ExternalField"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        systems[scannedNumberOfSystems -1].hasExternalField = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "VolumeMoveProbability"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        systems[scannedNumberOfSystems - 1].mc_moves_probabilities.probabilityVolumeMove = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsVolumeMoveProbability"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        systems[scannedNumberOfSystems - 1].mc_moves_probabilities.probabilityGibbsVolumeMove = value;
        continue;
      }
/*
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
    */
    }
  }

  for(size_t i = 0; i < totalNumberOfSystems; ++i)
  {
    if(computeConventionalRadialDistributionFunction[i])
    {
      systems[i].propertyConventionalRadialDistributionFunction =
        PropertyConventionalRadialDistributionFunction(numberOfBlocks,
                                                       systems[i].forceField.pseudoAtoms.size(),
                                                       conventionalRadialDistributionFunctionHistogramSize[i],
                                                       conventionalRadialDistributionFunctionRange[i],
                                                       sampleConventionalRadialDistributionFunctionEvery[i],
                                                       writeConventionalRadialDistributionFunctionEvery[i]);
    }

    if(computeRadialDistributionFunction[i])
    {
      systems[i].propertyRadialDistributionFunction =
        PropertyRadialDistributionFunction(numberOfBlocks,
                                           systems[i].forceField.pseudoAtoms.size(),
                                           radialDistributionFunctionHistogramSize[i],
                                           radialDistributionFunctionRange[i],
                                           sampleRadialDistributionFunctionEvery[i],
                                           writeRadialDistributionFunctionEvery[i]);
    }

    if(computeDensityGrid[i])
    {
      systems[i].propertyDensityGrid = PropertyDensityGrid(systems[i].frameworkComponents.size(),
                                                           systems[i].components.size(),
                                                           densityGridSize[i],
                                                           sampleDensityGridEvery[i],
                                                           writeDensityGridEvery[i]);
    }
  }

}

void InputReader::scanComponentSettings()
{
  std::stringstream s(inputString);

  std::string line{};
  std::string keyword{};
  std::string arguments{};

  [[maybe_unused]]std::size_t scannedNumberOfSystems{ 0uz };
  [[maybe_unused]]std::size_t scannedNumberOfComponents{ 0uz };
  [[maybe_unused]]size_t lineNumber{ 0uz };

  while (std::getline(s, line))
  {
    lineNumber += 1uz;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);           

      if (caseInSensStringCompare(keyword, std::string("Box")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Framework")))
      {
        scannedNumberOfSystems += 1uz;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        scannedNumberOfComponents++;
        continue;
      }

      if (caseInSensStringCompare(keyword, "MolFraction"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].molFraction = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "FugacityCoefficient"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].fugacityCoefficient = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "ThermodynamicIntegration"))
      {
        std::vector<bool> values = parseListOfSystemValues<bool>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].lambdaGC.computeDUdlambda = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "LnPartitionFunction"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].lnPartitionFunction = values[i];
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
        std::vector<bool> values = parseListOfSystemValues<bool>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].isCarrierGas = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "MassTransferCoefficient"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].massTransferCoefficient = values[i];
        }
        continue;
      }
      if (caseInSensStringCompare(keyword, "AxialDispersionCoefficient"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].axialDispersionCoefficient = values[i];
        }
        continue;
      }


      if (caseInSensStringCompare(keyword, "TranslationProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityTranslationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RandomTranslationProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
    
        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityRandomTranslationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RotationProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityRotationMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "RandomRotationProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityRandomRotationMove = values[i];
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
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityReinsertionMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "SwapConventionalProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilitySwapMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "SwapProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilitySwapMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsSwapProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityGibbsSwapMove_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "GibbsCFCMCSwapProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CFCMC_SwapProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilitySwapMove_CFCMC = values[i];
        }
        continue;
      }
      
      if (caseInSensStringCompare(keyword, "CFCMC_CBMC_SwapProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "WidomProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityWidomMove = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CFCMCWidomProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityWidomMove_CFCMC = values[i];
        }
        continue;
      }

      if (caseInSensStringCompare(keyword, "CBCFCMCWidomProbability"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

        values.resize(systems.size(), values.back());
        for (size_t i = 0uz; i < systems.size(); ++i)
        {
          systems[i].components[scannedNumberOfComponents - 1].mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC = values[i];
        }
        continue;
      }
/*
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
      */
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
