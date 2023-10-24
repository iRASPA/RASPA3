module;

module monte_carlo;

import <iostream>;
import <algorithm>;
import <numeric>;
import <chrono>;
import <vector>;
import <span>;
import <string>;
import <optional>;
import <fstream>;
import <sstream>;
import <filesystem>;
import <tuple>;
import <ios>;
import <print>;
import <complex>;
import <exception>;
import <source_location>;

import stringutils;
import archive;
import system;
import randomnumbers;
import mc_moves;
import input_reader;
import component;
import averages;
import loadings;
import units;
import enthalpy_of_adsorption;
import simulationbox;
import forcefield;
import sample_movies;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import atom;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import property_widom;
import property_simulationbox;
import property_energy;
import property_loading;
import property_enthalpy;
import mc_moves_probabilities_particles;
import mc_moves_cputime;
import mc_moves_count;
import property_pressure;
import transition_matrix;

MonteCarlo::MonteCarlo(InputReader& reader) noexcept : 
    numberOfCycles(reader.numberOfCycles),
    numberOfInitializationCycles(reader.numberOfInitializationCycles),
    numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
    printEvery(reader.printEvery),
    systems(std::move(reader.systems)),
    random(reader.randomSeed),
    estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
    
}

System& MonteCarlo::randomSystem()
{
  return systems[size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MonteCarlo::run()
{
  switch(simulationStage)
  {
    case SimulationStage::Initialization:
      goto continueInitializationStage;
    case SimulationStage::Equilibration:
      goto continueEquilibrationStage;
    case SimulationStage::Production:
      goto continueProductionStage;
    default:
      break;
  }

  continueInitializationStage: initialize();
  continueEquilibrationStage: equilibrate();
  t1 = std::chrono::system_clock::now();
  continueProductionStage: production();
  t2 = std::chrono::system_clock::now();

  output();
}

void MonteCarlo::createOutputFiles()
{
  for (System &system: systems)
  {
    std::string directoryNameString = std::format("Output/System_{}/", system.systemId);
    std::filesystem::path directoryName{ directoryNameString };
    std::filesystem::create_directories(directoryName);

    std::string fileNameString = std::format("Output/System_{}/output_{}_{}.data",
        system.systemId, system.temperature, system.input_pressure);
    streams.emplace_back(fileNameString, std::ios::out );
  }
}


void MonteCarlo::initialize()
{
  size_t totalNumberOfMolecules{ 0 };
  size_t totalNumberOfComponents{ 0 };
  size_t numberOfStepsPerCycle{ 0 };

  if(simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  for (System &system: systems)
  {
    system.registerEwaldFourierEnergySingleIon(double3(0.0, 0.0, 0.0), 1.0);
    system.removeRedundantMoves();
    system.determineSwapableComponents();
    system.determineFractionalComponents();
    system.rescaleMoveProbabilities();
    system.rescaleMolarFractions();
    system.computeComponentFluidProperties();
    system.computeFrameworkDensity();
    system.computeNumberOfPseudoAtoms();

    double3 perpendicularWidths = system.simulationBox.perpendicularWidths();
    system.forceField.initializeEwaldParameters(perpendicularWidths);

    system.createInitialMolecules(random);

    system.averageEnthalpiesOfAdsorption.resize(system.swapableComponents.size());
  }

  createOutputFiles();

  for(System & system : systems)
  {
    // switch the fractional molecule on in the first system, and off in all others
    if (system.systemId == 0) system.containsTheFractionalMolecule = true;
    else system.containsTheFractionalMolecule = false;
  }

  for (const System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}", system.writeOutputHeaderHardware());
    std::print(stream, "{}", Units::printStatus());
    std::print(stream, "{}", system.simulationBox.printParameters());
    std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(stream, "{}", system.forceField.printForceFieldStatus());
    std::print(stream, "{}", system.writeComponentStatus());
    std::print(stream, "{}", system.reactions.printStatus());
  }

  for (System& system : systems)
  {

    system.precomputeTotalRigidEnergy();
    system.recomputeTotalEnergies();

    std::ostream stream(streams[system.systemId].rdbuf());
    system.runningEnergies.print(stream, "Recomputed from scratch");
  };
  
  for (currentCycle = 0; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
    totalNumberOfMolecules = std::transform_reduce(systems.begin(), systems.end(), size_t(0),
        [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules();});
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, size_t(20)) * totalNumberOfComponents;

    for (size_t j = 0; j != numberOfStepsPerCycle; j++)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);
      particleMoves.performRandomMove(random, selectedSystem, selectSecondSystem, selectedComponent, fractionalMoleculeSystem);

      for(System &system : systems)
      {
        for(Component &component: system.components)
        {
          component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
        }
      }
    }

    if (currentCycle % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::flush(stream);
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % printEvery == 0)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if(ofile) 
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

    continueInitializationStage: ;
  }
}

void MonteCarlo::equilibrate()
{
  size_t totalNumberOfMolecules{ 0 };
  size_t totalNumberOfComponents{ 0 };
  size_t numberOfStepsPerCycle{ 0 };

  if(simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.recomputeTotalEnergies();
    system.runningEnergies.print(stream, "Recomputed from scratch");

    for(Component &component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize, system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  for (currentCycle = 0; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    totalNumberOfMolecules = std::transform_reduce(systems.begin(), systems.end(), size_t(0),
        [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules();});
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, size_t(20)) * totalNumberOfComponents;

    for (size_t j = 0; j != numberOfStepsPerCycle; j++)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);
      particleMoves.performRandomMove(random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);

      selectedSystem.components[selectedComponent].lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, 
                  selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, 
                  selectedSecondSystem.containsTheFractionalMolecule);

      selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSecondSystem.containsTheFractionalMolecule);
    }

    if (currentCycle % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::print(stream, "{}", system.writeEquilibrationStatusReport(currentCycle, numberOfEquilibrationCycles));
        std::flush(stream);
      }
    }

    if (currentCycle % 5000 == 0)
    {
      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors, system.containsTheFractionalMolecule);
        }
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % printEvery == 0)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if(ofile) 
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }
    continueEquilibrationStage: ;

  }
}

void MonteCarlo::production()
{
  size_t totalNumberOfMolecules{ 0 };
  size_t totalNumberOfComponents{ 0 };
  size_t numberOfStepsPerCycle{ 0 };
  double minBias{ 0.0 };

  if(simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.recomputeTotalEnergies(); 
    system.runningEnergies.print(stream, "Recomputed from scratch");

    system.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();
    system.mc_moves_count.clearCountStatistics();

    for(Component &component : system.components)
    {
      component.mc_moves_probabilities.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();
      component.mc_moves_count.clearCountStatistics();
      
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize, system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  minBias = std::numeric_limits<double>::max();
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      double currentMinBias = *std::min_element(component.lambdaGC.biasFactor.cbegin(), component.lambdaGC.biasFactor.cend());
      minBias = currentMinBias < minBias ? currentMinBias : minBias;
    }
  }
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      component.lambdaGC.normalize(minBias);
    }
  }

  numberOfSteps = 0;
  for (currentCycle = 0; currentCycle != numberOfCycles; ++currentCycle)
  {

    estimation.setCurrentSample(currentCycle);

    totalNumberOfMolecules = std::transform_reduce(systems.begin(), systems.end(), size_t(0),
        [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules();});
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, size_t(20)) * totalNumberOfComponents;

    for (size_t j = 0; j != numberOfStepsPerCycle; j++)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);
      
      particleMoves.performRandomMoveProduction(random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem, estimation.currentBin);

      selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSecondSystem.containsTheFractionalMolecule);

      ++numberOfSteps;
    }


    for (System& system : systems)
    {
      system.sampleProperties(estimation.currentBin);

      std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
      std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
      system.currentEnergyStatus = molecularPressure.first;
      system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
      std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
      system.mc_moves_cputime.energyPressureComputation += (time2 - time1);

      // add the sample energy to the averages
      if (currentCycle % 10 == 0 || currentCycle % printEvery == 0)
      {
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }
    }

    if (currentCycle % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        std::print(stream, "{}", system.writeProductionStatusReport(currentCycle, numberOfCycles));
        std::flush(stream);
      }
    }



    if (currentCycle % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % printEvery == 0)
    {
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if(ofile) 
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }


    }
    //if (currentCycle % printEvery == 0)
    //{
    //  std::ofstream ofile("test.bin", std::ios::binary);
    //  Archive<std::ofstream> archive(ofile);
    //  archive << *this;
    //  ofile.close();

    //  std::ifstream ifile("test.bin", std::ios::binary);
    //  if(!ifile.is_open())
    //  {
    //      std::cout << "ERROR file doesn't exist.." << std::endl;
    //      return;
    //  }

    //  Archive<std::ifstream> archive2(ifile);
    //  MonteCarlo b{};
    //  archive2 >> b;
    //  ifile.close();

    //  std::cout << std::boolalpha << (*this == b) << std::endl;
    //}
    continueProductionStage: ;

  }
}

void MonteCarlo::output()
{
  MCMoveCpuTime total;
  MCMoveCount countTotal;
  for(const System &system: systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_count;
  }

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies.print(stream, "Running energies");
    
    RunningEnergy recomputedEnergies = system.computeTotalEnergies();
    recomputedEnergies.print(stream, "Recomputed from scratch");
    
    RunningEnergy drift = system.runningEnergies - recomputedEnergies;
    drift.print(stream, "Monte-Carlo energy drift");

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");
    
    std::print(stream, "{}", system.writeMCMoveStatistics());
    //std::print(stream, system.lambda.writeAveragesStatistics(system.beta));

    std::print(stream, "Production run counting of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for(const Component &component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_count.writeMCMoveCountComponentStatistics(numberOfSteps, component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_count.writeMCMoveCountSystemStatistics(numberOfSteps));

    std::print(stream, "Production run counting of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", countTotal.writeMCMoveCountAllSystemStatistics(numberOfSteps));

    std::print(stream, "\n\n");

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for(const Component &component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());

    std::print(stream, "Production run CPU timings of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::chrono::duration<double> totalSimulationTime = (t2 - t1);
    std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalSimulationTime));

    std::print(stream, "\n\n");
    //std::print(stream, "\nProduction simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count());

    std::print(stream, "{}", system.averageEnergies.writeAveragesStatistics(system.components));
    std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());
    std::print(stream, "{}", system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swapableComponents, system.components));
    std::print(stream, "{}", system.averageLoadings.writeAveragesStatistics(system.components, system.frameworkMass));
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MonteCarlo& mc)
{
  archive << mc.versionNumber;

  archive << mc.numberOfCycles;
  archive << mc.numberOfSteps;
  archive << mc.numberOfInitializationCycles;
  archive << mc.numberOfEquilibrationCycles;
  archive << mc.optimizeMCMovesEvery;
  archive << mc.printEvery;

  archive << mc.currentCycle;
  archive << mc.simulationStage;

  archive << mc.systems;
  archive << mc.random;

  archive << mc.fractionalMoleculeSystem;

  archive << mc.estimation;

  archive << mc.particleMoves;

  archive << static_cast<uint64_t>(0x6f6b6179); // magic number 'okay' in hex
 
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MonteCarlo& mc)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MonteCarlo' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }


  archive >> mc.numberOfCycles;
  archive >> mc.numberOfSteps;
  archive >> mc.numberOfInitializationCycles;
  archive >> mc.numberOfEquilibrationCycles;
  archive >> mc.optimizeMCMovesEvery;
  archive >> mc.printEvery;

  archive >> mc.currentCycle;
  archive >> mc.simulationStage;

  archive >> mc.systems;
  archive >> mc.random;

  archive >> mc.fractionalMoleculeSystem;

  archive >> mc.estimation;

  archive >> mc.particleMoves;


  uint64_t magicNumber;
  archive >> magicNumber;
  if(magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
  }
  std::print("Magic number read correctly: {} vs {}\n", magicNumber, static_cast<uint64_t>(0x6f6b6179));
  return archive;
}
