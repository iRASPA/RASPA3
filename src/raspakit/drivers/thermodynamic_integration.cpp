module;

module thermodynamic_integration;

import std;

import stringutils;
import hardware_info;
import system;
import framework;
import randomnumbers;
import input_reader;
import component;
import averages;
import property_loading;
import units;
import simulationbox;
import forcefield;
import energy_status;
import running_energy;
import atom;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import mc_moves;
import mc_moves_cputime;
import mc_moves_statistics;
import json;

ThermodynamicIntegration::ThermodynamicIntegration(InputReader& reader) noexcept
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      systems(std::move(reader.systems)),
      outputJsons(systems.size()),
      estimation(reader.numberOfBlocks, reader.numberOfProductionCycles)
{
}

void ThermodynamicIntegration::run()
{
  setup();
  initialize();
  equilibrate();
  production();
  output();
}

void ThermodynamicIntegration::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    // every system carries its own fractional molecule, pinned at the fixed lambda
    system.containsTheFractionalMolecule = true;
  }

  std::filesystem::create_directories("output");
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.txt", system.temperature, system.input_pressure, system_id);
    streams.emplace_back(fileNameString, std::ios::out);
    fileNameString = std::format("output/output_{}_{}.s{}.json", system.temperature, system.input_pressure, system_id);
    outputJsonFileNames.emplace_back(fileNameString);

    ++system_id;
  }

  for (std::size_t system_id{0}; const System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());
    std::print(stream, "{}", system.writeSystemStatus());
    std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(stream, "{}", system.forceField.printForceFieldStatus());
    std::print(stream, "{}", system.writeComponentStatus());
    std::print(stream, "{}", system.writeNumberOfPseudoAtoms());

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
    outputJsons[system_id]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif

    outputJsons[system_id]["seed"] = random.seed;
    outputJsons[system_id]["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
    outputJsons[system_id]["initialization"]["units"] = Units::jsonStatus();
    outputJsons[system_id]["initialization"]["initialConditions"] = system.jsonSystemStatus();
    outputJsons[system_id]["initialization"]["forceField"] = system.forceField.jsonForceFieldStatus();
    outputJsons[system_id]["initialization"]["forceField"]["pseudoAtoms"] = system.forceField.jsonPseudoAtomStatus();
    outputJsons[system_id]["initialization"]["components"] = system.jsonComponentStatus();

    std::ofstream json(outputJsonFileNames[system_id]);
    json << outputJsons[system_id].dump(4);

    ++system_id;
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());
    system.createExternalFieldInterpolationGrid(stream, system_id);
    system.createFrameworkInterpolationGrids(stream);

    ++system_id;
  }
}

void ThermodynamicIntegration::performCycle()
{
  std::size_t totalNumberOfMolecules{0uz};
  std::size_t totalNumberOfComponents{0uz};
  std::size_t numberOfStepsPerCycle{0uz};

  totalNumberOfMolecules = std::transform_reduce(
      systems.begin(), systems.end(), 0uz, [](const std::size_t& acc, const std::size_t& b) { return acc + b; },
      [](const System& system) { return system.numberOfMolecules(); });
  totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

  numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

  for (std::size_t j = 0uz; j != numberOfStepsPerCycle; j++)
  {
    std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
    System& selectedSystem = systems[selectedSystemPair.first];
    System& selectedSecondSystem = systems[selectedSystemPair.second];

    std::size_t selectedComponent = selectedSystem.randomComponent(random);

    // the fixed lambda-bin never changes: no lambda-changing moves are enabled and no
    // Wang-Landau biasing is needed
    switch (simulationStage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                  fractionalMoleculeSystem);
        break;
      case SimulationStage::Equilibration:
        MC_Moves::performRandomMoveEquilibration(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                 fractionalMoleculeSystem);
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                              fractionalMoleculeSystem, estimation.currentBin);
        numberOfSteps++;
        break;
    }

    selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSystem.containsTheFractionalMolecule);
    selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
        selectedSecondSystem.containsTheFractionalMolecule);
    // pinned pair- and group-swap lambda histograms
    selectedSystem.pairSwapLambdaSampleOccupancy();
    selectedSecondSystem.pairSwapLambdaSampleOccupancy();
  }
}

void ThermodynamicIntegration::initialize()
{
  std::chrono::steady_clock::time_point t1, t2;

  simulationStage = SimulationStage::Initialization;

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    std::ostream stream(streams[system_id].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    system.writeRestartFile(system_id);

    ++system_id;
  }

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; ++currentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::ostream stream(streams[system_id].rdbuf());
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::flush(stream);

        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.writeRestartFile(system_id);

        ++system_id;
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalInitializationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);
  }
}

void ThermodynamicIntegration::equilibrate()
{
  std::chrono::steady_clock::time_point t1, t2;

  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
  }

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::ostream stream(streams[system_id].rdbuf());
        std::print(stream, "{}", system.writeEquilibrationStatusReportMC(currentCycle, numberOfEquilibrationCycles));
        std::flush(stream);

        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.writeRestartFile(system_id);

        ++system_id;
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalEquilibrationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);
  }
}

void ThermodynamicIntegration::production()
{
  std::chrono::steady_clock::time_point t1, t2;

  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();

    system.mc_moves_statistics.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();

    for (Component& component : system.components)
    {
      component.mc_moves_statistics.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();

      // start production with fresh dU/dlambda accumulators; the bin stays fixed
      component.lambdaGC.clear();
    }
    // pinned pair- and group-swap lambda histograms
    system.pairSwapLambdaClearBookkeeping();
  }

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfProductionCycles; ++currentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    estimation.setCurrentSample(currentCycle);

    performCycle();

    for (std::size_t system_id{0}; System& system : systems)
    {
      system.sampleProperties(system_id, estimation.currentBin, currentCycle);

      if (currentCycle % 10uz == 0uz || currentCycle % printEvery == 0uz)
      {
        std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
        std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
        system.currentEnergyStatus = molecularPressure.first;
        system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
        std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();

        system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }

      ++system_id;
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::ostream stream(streams[system_id].rdbuf());
        std::string status_line{std::format("Current cycle: {} out of {}\n", currentCycle, numberOfProductionCycles)};
        std::print(stream, "{}", system.writeProductionStatusReportMC(status_line));
        std::flush(stream);

        ++system_id;
      }
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.writeRestartFile(system_id);

        ++system_id;
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalProductionSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);
  }

  // write final state after the simulation has finished
  for (std::size_t system_id{0}; System& system : systems)
  {
    system.writeRestartFile(system_id);

    std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
    system.currentEnergyStatus = molecularPressure.first;
    system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;

    system.loadings =
        LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

    std::ostream stream(streams[system_id].rdbuf());

    std::print(stream, "\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "                             Simulation finished!\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "\n");

    std::string status_line{std::format("Final state after {} cycles\n", numberOfProductionCycles)};
    std::print(stream, "{}", system.writeProductionStatusReportMC(status_line));
    std::flush(stream);

    ++system_id;
  }
}

void ThermodynamicIntegration::output()
{
  MCMoveCpuTime total;
  MCMoveStatistics countTotal;
  for (const System& system : systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_statistics;
    for (const Component& component : system.components)
    {
      countTotal += component.mc_moves_statistics;
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    RunningEnergy recomputedEnergies = system.computeTotalEnergies();
    RunningEnergy drift = system.runningEnergies - recomputedEnergies;
    stream << system.runningEnergies.printMCDiff(recomputedEnergies);

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", system.writeMCMoveStatistics());

    std::print(stream, "Production run counting of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));

    std::print(stream, "\n\n");

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));

      ++componentId;
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());

    std::print(stream, "Production run CPU timings of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
    std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
    std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
    std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
    std::print(stream, "Total simulation time:          {:14f} [s]\n", totalSimulationTime.count());
    std::print(stream, "\n\n");

    std::print(
        stream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));

    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());
    }

    std::print(stream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));

    std::print(stream, "{}", writeThermodynamicIntegrationPoint(system));

    std::flush(stream);

    // json statistics
    outputJsons[system_id]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    outputJsons[system_id]["output"]["recomputedEnergies"] = recomputedEnergies.jsonMC();
    outputJsons[system_id]["output"]["drift"] = drift.jsonMC();

    outputJsons[system_id]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();

    outputJsons[system_id]["output"]["cpuTimings"]["summedSystemsAndComponents"] =
        total.jsonOverallMCMoveCPUTimeStatistics(totalProductionSimulationTime);
    outputJsons[system_id]["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["total"] = totalSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();

    outputJsons[system_id]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    outputJsons[system_id]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();
    outputJsons[system_id]["properties"]["thermodynamicIntegration"] = jsonThermodynamicIntegrationPoint(system);

    for (const Component& component : system.components)
    {
      outputJsons[system_id]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }

    std::ofstream json(outputJsonFileNames[system_id]);
    json << outputJsons[system_id].dump(4);

    ++system_id;
  }
}

std::string ThermodynamicIntegration::writeThermodynamicIntegrationPoint(const System& system) const
{
  std::ostringstream stream;

  std::print(stream, "Thermodynamic integration at fixed lambda\n");
  std::print(stream, "===============================================================================\n\n");

  for (std::size_t componentId{0}; const Component& component : system.components)
  {
    if (component.fixedLambdaBin.has_value())
    {
      const PropertyLambdaProbabilityHistogram& lambda = component.fixedLambdaHistogram();
      const std::size_t binIndex = component.fixedLambdaBin.value();
      const double lambdaValue = static_cast<double>(binIndex) * lambda.delta;

      std::pair<std::vector<double>, std::vector<double>> dudlambda = lambda.averageDuDlambda();

      std::print(stream, "component {} ({})\n", componentId, component.name);
      std::print(stream, "    lambda-bin index: {} (out of 0-{})\n", binIndex, lambda.numberOfSamplePoints - 1);
      std::print(stream, "    lambda:           {:.6f}\n", lambdaValue);
      for (std::size_t blockIndex = 0; blockIndex < lambda.numberOfBlocks; ++blockIndex)
      {
        std::print(stream, "        Block[ {:2d}] <dU/dlambda>: {: .6e} [K]\n", blockIndex,
                   Units::EnergyToKelvin * lambda.averagedDUdlambda(blockIndex)[binIndex]);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    <dU/dlambda>: {: .6e} +/- {: .6e} [K]\n",
                 Units::EnergyToKelvin * dudlambda.first[binIndex], Units::EnergyToKelvin * dudlambda.second[binIndex]);
      std::print(stream, "    <dU/dlambda>: {: .6e} +/- {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * dudlambda.first[binIndex],
                 Units::EnergyToKJPerMol * dudlambda.second[binIndex]);
      std::print(stream, "\n");
    }
    ++componentId;
  }
  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json ThermodynamicIntegration::jsonThermodynamicIntegrationPoint(const System& system) const
{
  nlohmann::json status;

  for (const Component& component : system.components)
  {
    if (component.fixedLambdaBin.has_value())
    {
      const PropertyLambdaProbabilityHistogram& lambda = component.fixedLambdaHistogram();
      const std::size_t binIndex = component.fixedLambdaBin.value();

      std::pair<std::vector<double>, std::vector<double>> dudlambda = lambda.averageDuDlambda();

      std::vector<double> blockAverages(lambda.numberOfBlocks);
      for (std::size_t blockIndex = 0; blockIndex < lambda.numberOfBlocks; ++blockIndex)
      {
        blockAverages[blockIndex] = Units::EnergyToKelvin * lambda.averagedDUdlambda(blockIndex)[binIndex];
      }

      status[component.name]["lambdaBinIndex"] = binIndex;
      status[component.name]["numberOfLambdaBins"] = lambda.numberOfSamplePoints;
      status[component.name]["lambda"] = static_cast<double>(binIndex) * lambda.delta;
      status[component.name]["blockAverageDUdlambda[K]"] = blockAverages;
      status[component.name]["averageDUdlambda[K]"] = Units::EnergyToKelvin * dudlambda.first[binIndex];
      status[component.name]["confidenceDUdlambda[K]"] = Units::EnergyToKelvin * dudlambda.second[binIndex];
      status[component.name]["averageDUdlambda[kJ/mol]"] = Units::EnergyToKJPerMol * dudlambda.first[binIndex];
      status[component.name]["confidenceDUdlambda[kJ/mol]"] = Units::EnergyToKJPerMol * dudlambda.second[binIndex];
    }
  }

  return status;
}
