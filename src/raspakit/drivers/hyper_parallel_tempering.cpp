module;

module hyper_parallel_tempering;

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
import equation_of_states;
import energy_status;
import running_energy;
import atom;
import int3;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import mc_moves;
import mc_moves_cputime;
import mc_moves_statistics;
import mc_moves_parallel_tempering_swap;
import json;

// The analysis-property writers (RDFs, density grid, histograms, molecule properties) gate
// themselves on their own 'writeEvery'; a cycle argument of 0 forces the write (used for the
// final flush at the end of the run). The replica id keys the output filenames, so every replica
// writes its own set of files.
static void writeReplicaAnalysisOutputs(System& system, std::size_t replicaId, std::size_t cycle)
{
  if (system.propertyConventionalRadialDistributionFunction.has_value())
  {
    system.propertyConventionalRadialDistributionFunction->writeOutput(
        system.forceField, replicaId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyRadialDistributionFunction.has_value())
  {
    system.propertyRadialDistributionFunction->writeOutput(system.forceField, replicaId, system.simulationBox.volume,
                                                           system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyDensityGrid.has_value())
  {
    system.propertyDensityGrid->writeOutput(replicaId, system.simulationBox, system.forceField, system.framework,
                                            system.components, cycle);
  }
  if (system.averageEnergyHistogram.has_value())
  {
    system.averageEnergyHistogram->writeOutput(replicaId, cycle);
  }
  if (system.averageNumberOfMoleculesHistogram.has_value())
  {
    system.averageNumberOfMoleculesHistogram->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyMoleculeProperties.has_value())
  {
    system.propertyMoleculeProperties->writeOutput(replicaId, system.components, cycle);
  }
}

HyperParallelTempering::HyperParallelTempering(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      numberOfBlocks(reader.numberOfBlocks),
      parallelTemperingSwapEvery(reader.parallelTemperingSwapEvery),
      temperatures(reader.parallelTemperingTemperatures),
      pressures(reader.parallelTemperingPressures),
      numberOfTemperatures(temperatures.size()),
      numberOfPressures(pressures.size()),
      numberOfReplicas(numberOfTemperatures * numberOfPressures)
{
  // the single declared system is replicated into one replica per (temperature, pressure) grid point
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  systems.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  // replica (t, p) is pinned at temperature T_t and pressure P_p, with its own random-number stream
  randoms.reserve(numberOfReplicas);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      const std::size_t replicaId = replicaIndex(temperatureIndex, pressureIndex);
      System& system = systems[replicaId];
      const double T = temperatures[temperatureIndex];
      const double P = pressures[pressureIndex];

      system.temperature = T;
      system.beta = 1.0 / (Units::KB * T);

      system.input_pressure = P;
      system.pressure = P / Units::PressureConversionFactor;
      system.input_pressureTensorDiagonal = double3(system.input_pressure, system.input_pressure,
                                                    system.input_pressure);
      system.pressureTensorDiagonal = double3(system.pressure, system.pressure, system.pressure);

      // temperature-dependent potentials (Feynman-Hibbs) derive pair coefficients, shifts and
      // tail-corrections from the temperature
      if (system.forceField.temperature != T)
      {
        system.forceField.temperature = T;
        system.forceField.preComputeDerivedParameters();
        system.forceField.preComputePotentialShift();
        system.forceField.preComputeTailCorrection();
      }

      // convert the pressure to per-component fugacities at this grid point: the fugacity
      // coefficients are recomputed with the Peng-Robinson equation of state at (T_t, P_p)
      // (an explicitly given 'FugacityCoefficient' would only be valid at one state point)
      for (Component& component : system.components)
      {
        component.fugacityCoefficient = std::nullopt;
      }
      system.equationOfState = EquationOfState(EquationOfState::Type::PengRobinson,
                                               EquationOfState::MixingRules::VanDerWaals, T, P, system.simulationBox,
                                               system.heliumVoidFraction, system.components);

      // the CBMC ideal-gas conformation reservoirs are Boltzmann samples at the system temperature
      system.buildConformationReservoirs();

      randoms.emplace_back(random.seed + replicaId + 1);
    }
  }

  stepsPerReplica.assign(numberOfReplicas, 0uz);
  swapAttemptsPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAcceptedPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAttemptsPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
  swapAcceptedPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
}

void HyperParallelTempering::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void HyperParallelTempering::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  stream.open("output/output.hyper_parallel_tempering.txt", std::ios::out);
  outputJsonFileName = "output/output.hyper_parallel_tempering.json";

  const System& front = systems.front();
  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  std::print(stream, "Hyper-parallel tempering\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Number of temperatures:                      {}\n", numberOfTemperatures);
  std::print(stream, "Number of pressures:                         {}\n", numberOfPressures);
  std::print(stream, "Number of replicas / threads:                {}\n", numberOfReplicas);
  std::print(stream, "Temperature ladder:                         ");
  for (double T : temperatures)
  {
    std::print(stream, " {}", T);
  }
  std::print(stream, " [K]\n");
  std::print(stream, "Pressure ladder:                            ");
  for (double P : pressures)
  {
    std::print(stream, " {}", P);
  }
  std::print(stream, " [Pa]\n");
  if (parallelTemperingSwapEvery == 0uz)
  {
    std::print(stream, "Configuration swaps:                         disabled\n\n");
  }
  else
  {
    std::print(stream, "Configuration-swap sweep every:              {} cycles\n", parallelTemperingSwapEvery);
    std::print(stream, "  (sweeps alternate between the temperature and the pressure direction of the grid)\n\n");
  }

  std::print(stream, "Replica grid: replica (t, p) = t * {} + p, fugacity coefficients from Peng-Robinson\n",
             numberOfPressures);
  std::print(stream, "    replica    temperature [K]    pressure [Pa]\n");
  std::print(stream, "    ------------------------------------------------\n");
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      const std::size_t replicaId = replicaIndex(temperatureIndex, pressureIndex);
      std::print(stream, "    {:7d}    {:15.4f}    {:13.5e}\n", replicaId, temperatures[temperatureIndex],
                 pressures[pressureIndex]);
    }
  }
  std::print(stream, "\n");

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["temperatures"] = temperatures;
  outputJson["initialization"]["pressures"] = pressures;
  outputJson["initialization"]["parallelTemperingSwapEvery"] = parallelTemperingSwapEvery;

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-replica output files: each worker thread writes exclusively to its own stream
  replicaStreams.reserve(numberOfReplicas);
  replicaJsonFileNames.reserve(numberOfReplicas);
  replicaJsons.resize(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(std::format("output/output_{}_{}.hyper_parallel_tempering.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                std::ios::out);
    replicaJsonFileNames.emplace_back(std::format("output/output_{}_{}.hyper_parallel_tempering.r{}.json",
                                                  system.temperature, system.input_pressure, replicaId));

    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "Hyper-parallel tempering: replica {} of {} (temperature {} [K], pressure {} [Pa])\n",
               replicaId, numberOfReplicas, system.temperature, system.input_pressure);
    std::print(replicaStream, "Random seed of this replica: {}\n\n", randoms[replicaId].seed);
    std::print(replicaStream, "{}\n", HardwareInfo::writeInfo());
    std::print(replicaStream, "{}", Units::printStatus());
    std::print(replicaStream, "{}", system.writeSystemStatus());
    std::print(replicaStream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(replicaStream, "{}", system.forceField.printForceFieldStatus());
    std::print(replicaStream, "{}", system.writeComponentStatus());
    std::print(replicaStream, "{}", system.writeNumberOfPseudoAtoms());

#ifdef VERSION
    replicaJsons[replicaId]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
    replicaJsons[replicaId]["seed"] = randoms[replicaId].seed;
    replicaJsons[replicaId]["replicaId"] = replicaId;
    replicaJsons[replicaId]["temperature"] = system.temperature;
    replicaJsons[replicaId]["pressure"] = system.input_pressure;
    replicaJsons[replicaId]["initialization"]["initialConditions"] = system.jsonSystemStatus();
    replicaJsons[replicaId]["initialization"]["components"] = system.jsonComponentStatus();

    std::ofstream replicaJson(replicaJsonFileNames[replicaId]);
    replicaJson << replicaJsons[replicaId].dump(4);
  }

  // interpolation grids are computed once and shared (copied) between the replicas
  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t replicaId = 1; replicaId < systems.size(); ++replicaId)
  {
    systems[replicaId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[replicaId].interpolationGrids = systems.front().interpolationGrids;
  }
}

void HyperParallelTempering::performReplicaCycle(std::size_t replicaId, SimulationStage stage,
                                                 std::size_t currentBlock)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];

  // every replica is self-contained; the Gibbs-style moves that need a partner system are not
  // supported by this driver
  std::size_t fractionalMoleculeSystem = 0uz;

  const std::size_t numberOfStepsPerCycle =
      std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();

  for (std::size_t j = 0uz; j != numberOfStepsPerCycle; ++j)
  {
    std::size_t selectedComponent = system.randomComponent(rng);

    switch (stage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::PreInitialization:
        MC_Moves::performRandomMovePreInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Equilibration:
        MC_Moves::performRandomMoveEquilibration(rng, system, system, selectedComponent, fractionalMoleculeSystem);

        // Wang-Landau biasing of the CFCMC lambda moves (all state is owned by this replica)
        system.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, system.containsTheFractionalMolecule);
        system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem,
                                              currentBlock);
        ++stepsPerReplica[replicaId];
        break;
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
    system.pairSwapLambdaSampleOccupancy();
    system.reactionLambdaSampleOccupancy();
  }
}

void HyperParallelTempering::performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept
{
  // the replicas keep their (temperature, pressure) state points; the configurations migrate
  // through the grid. Sweeps alternate between the temperature direction and the pressure
  // direction; within a direction the pairing offset alternates between visits, so a configuration
  // can traverse the whole grid.
  bool temperatureDirection = (swapSweeps % 2uz == 0uz);
  if (temperatureDirection && numberOfTemperatures < 2uz) temperatureDirection = false;
  if (!temperatureDirection && numberOfPressures < 2uz) temperatureDirection = true;
  const std::size_t offset = (swapSweeps / 2uz) % 2uz;

  if (temperatureDirection)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      for (std::size_t temperatureIndex = offset; temperatureIndex + 1 < numberOfTemperatures; temperatureIndex += 2uz)
      {
        ++swapAttempts;
        ++swapAttemptsPerTemperaturePair[temperatureIndex];
        if (MC_Moves::ParallelTemperingSwap(random, systems[replicaIndex(temperatureIndex, pressureIndex)],
                                            systems[replicaIndex(temperatureIndex + 1, pressureIndex)])
                .has_value())
        {
          ++swapAccepted;
          ++swapAcceptedPerTemperaturePair[temperatureIndex];
        }
      }
    }
  }
  else
  {
    for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
    {
      for (std::size_t pressureIndex = offset; pressureIndex + 1 < numberOfPressures; pressureIndex += 2uz)
      {
        ++swapAttempts;
        ++swapAttemptsPerPressurePair[pressureIndex];
        if (MC_Moves::ParallelTemperingSwap(random, systems[replicaIndex(temperatureIndex, pressureIndex)],
                                            systems[replicaIndex(temperatureIndex, pressureIndex + 1)])
                .has_value())
        {
          ++swapAccepted;
          ++swapAcceptedPerPressurePair[pressureIndex];
        }
      }
    }
  }
  ++swapSweeps;
  ++sweepsThisStage;

  // progress report: all worker threads are parked on the barrier, so this is race-free
  const std::size_t printSweeps = std::max(1uz, printEvery / std::max(1uz, parallelTemperingSwapEvery));
  if (sweepsThisStage % printSweeps == 0uz)
  {
    const std::string_view stageName = (stage == SimulationStage::PreInitialization)  ? "pre-initialization"
                                       : (stage == SimulationStage::Initialization)   ? "initialization"
                                       : (stage == SimulationStage::Equilibration)    ? "equilibration"
                                                                                      : "production";
    const std::size_t cycle = std::min(sweepsThisStage * parallelTemperingSwapEvery, numberOfCycles);
    {
      std::scoped_lock lock(outputMutex);
      std::print(stream, "Hyper-parallel-tempering sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%)\n",
                 swapSweeps, stageName, cycle, numberOfCycles, swapAccepted, swapAttempts,
                 100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));
      std::flush(stream);
    }

    // running snapshot of the assembled adsorption isotherms: all worker threads are parked on
    // the barrier, so reading the per-replica loading averages is race-free
    if (stage == SimulationStage::Production)
    {
      writeIsothermSnapshot();
    }
  }
}

void HyperParallelTempering::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  simulationStage = stage;

  // serial per-stage preparation
  if (stage == SimulationStage::Equilibration)
  {
    for (System& system : systems)
    {
      for (Component& component : system.components)
      {
        component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                               system.containsTheFractionalMolecule);
        component.lambdaGC.clear();
      }
      system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
      system.pairSwapLambdaClearBookkeeping();
      system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
      system.reactionLambdaClearBookkeeping();
    }
  }
  if (stage == SimulationStage::Production)
  {
    for (System& system : systems)
    {
      system.mc_moves_statistics.clearMoveStatistics();
      system.mc_moves_cputime.clearTimingStatistics();

      for (Component& component : system.components)
      {
        component.mc_moves_statistics.clearMoveStatistics();
        component.mc_moves_cputime.clearTimingStatistics();

        component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                               system.containsTheFractionalMolecule);
        component.lambdaGC.clear();
      }
      system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize);
      system.pairSwapLambdaClearBookkeeping();
      system.reactionLambdaFinalize();
      system.reactionLambdaClearBookkeeping();
    }
    std::fill(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);
  }
  sweepsThisStage = 0uz;

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::PreInitialization)  ? "Pre-initialization"
                                       : (stage == SimulationStage::Initialization)   ? "Initialization"
                                       : (stage == SimulationStage::Equilibration)    ? "Equilibration"
                                                                                      : "Production";
    std::print(stream, "\n{} stage: {} cycles on {} replicas/threads\n", stageName, numberOfCycles, numberOfReplicas);
    std::flush(stream);
  }

  // one worker thread per replica; the barrier completion performs the swap sweeps
  auto onAllArrived = [this, stage, numberOfCycles]() noexcept { performSwapSweep(stage, numberOfCycles); };
  std::barrier synchronizationPoint(static_cast<std::ptrdiff_t>(numberOfReplicas), onAllArrived);

  const std::size_t stageCycleOffset = absoluteCycleOffset;

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, stageCycleOffset, &synchronizationPoint]()
          {
            System& system = systems[replicaId];

            // each thread computes the total energies of its own replica
            if (stage == SimulationStage::PreInitialization || stage == SimulationStage::Initialization)
            {
              system.precomputeTotalRigidEnergy();
            }
            system.runningEnergies = system.computeTotalEnergies();

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfProductionCycles));

            for (std::size_t cycle = 0uz; cycle != numberOfCycles; ++cycle)
            {
              if (stage == SimulationStage::Production)
              {
                estimation.setCurrentSample(cycle);
              }

              performReplicaCycle(replicaId, stage, estimation.currentBin);

              // time-evolution properties (number of molecules, volume): sampled over all stages,
              // indexed by the absolute cycle number; the writers gate on their own 'writeEvery'
              const std::size_t absoluteCycle = stageCycleOffset + cycle;
              system.samplePropertiesEvolution(absoluteCycle);
              if (system.propertyNumberOfMoleculesEvolution.has_value())
              {
                system.propertyNumberOfMoleculesEvolution->writeOutput(replicaId, absoluteCycle);
              }
              if (system.propertyVolumeEvolution.has_value())
              {
                system.propertyVolumeEvolution->writeOutput(replicaId, absoluteCycle);
              }

              if (stage == SimulationStage::Production)
              {
                system.sampleProperties(replicaId, estimation.currentBin, cycle);

                // analysis-property files (RDFs, density grid, histograms, molecule properties);
                // the writers gate on their own 'writeEvery'
                writeReplicaAnalysisOutputs(system, replicaId, cycle);

                // energy/pressure averages for the per-replica final report
                if (cycle % 10uz == 0uz || cycle % printEvery == 0uz)
                {
                  std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
                  std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
                  system.currentEnergyStatus = molecularPressure.first;
                  system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
                  std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();

                  system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
                  system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
                }
              }

              if (cycle % optimizeMCMovesEvery == 0uz)
              {
                system.optimizeMCMoves();
              }

              // Wang-Landau biasing-factor adjustment (all state is owned by this replica)
              if (stage == SimulationStage::Equilibration && cycle % rescaleWangLandauEvery == 0uz)
              {
                for (Component& component : system.components)
                {
                  component.lambdaGC.WangLandauIteration(
                      PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
                      system.containsTheFractionalMolecule);
                }
                system.pairSwapLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
                system.reactionLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
              }

              if (cycle % printEvery == 0uz)
              {
                // each thread writes exclusively to its own replica stream: no locking needed
                system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                              system.simulationBox);

                std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
                switch (stage)
                {
                  case SimulationStage::PreInitialization:
                    std::print(replicaStream, "{}", system.writePreInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Initialization:
                    std::print(replicaStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(replicaStream, "{}", system.writeEquilibrationStatusReportMC(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                  {
                    std::string status_line = std::format("Current cycle: {} out of {}\n", cycle, numberOfCycles);
                    std::print(replicaStream, "{}", system.writeProductionStatusReportMC(status_line));
                    break;
                  }
                  default:
                    break;
                }
                std::flush(replicaStream);
              }

              // the only synchronization point between the threads: the configuration-swap sweep
              if (parallelTemperingSwapEvery != 0uz && (cycle + 1uz) % parallelTemperingSwapEvery == 0uz)
              {
                synchronizationPoint.arrive_and_wait();
              }
            }
          });
    }
    // jthreads join on scope exit
  }

  absoluteCycleOffset += numberOfCycles;

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  switch (stage)
  {
    case SimulationStage::PreInitialization:
      totalPreInitializationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Initialization:
      totalInitializationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Equilibration:
      totalEquilibrationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Production:
      totalProductionSimulationTime += (t2 - t1);
      break;
    default:
      break;
  }
  totalSimulationTime += (t2 - t1);
}

void HyperParallelTempering::output()
{
  std::size_t numberOfSteps = std::accumulate(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);

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

  std::print(stream, "\n");
  std::print(stream, "===============================================================================\n");
  std::print(stream, "                             Simulation finished!\n");
  std::print(stream, "===============================================================================\n");
  std::print(stream, "\n");

  // energy drift check of every replica (energies recomputed in parallel, one thread per replica);
  // the final state used by the per-replica reports is refreshed in the same pass
  std::vector<RunningEnergy> recomputed(systems.size());
  {
    std::vector<std::jthread> threads;
    threads.reserve(systems.size());
    for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, &recomputed]()
          {
            System& system = systems[replicaId];
            recomputed[replicaId] = system.computeTotalEnergies();

            std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
            system.currentEnergyStatus = molecularPressure.first;
            system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
            system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                          system.simulationBox);
          });
    }
  }

  writeReplicaFinalReports(recomputed);

  // final assembled adsorption isotherms (also covers runs with the swaps disabled)
  writeIsothermSnapshot();

  std::print(stream, "Energy drift per replica\n");
  std::print(stream, "===============================================================================\n\n");
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    const RunningEnergy drift = systems[replicaId].runningEnergies - recomputed[replicaId];
    std::print(stream, "    replica {:4d} (temperature {:10.4f} [K], pressure {:13.5e} [Pa]): drift {: .6e} [K]\n",
               replicaId, systems[replicaId].temperature, systems[replicaId].input_pressure,
               Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  std::print(stream, "Hyper-parallel-tempering swap statistics\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    sweeps:    {}\n", swapSweeps);
  std::print(stream, "    attempts:  {}\n", swapAttempts);
  std::print(stream, "    accepted:  {} ({:.4f} %)\n\n", swapAccepted,
             100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));

  // low acceptance for a particular pair marks a bottleneck in the grid (configurations cannot
  // migrate past it); consider a denser ladder around such a pair
  if (numberOfTemperatures > 1uz)
  {
    std::print(stream, "    temperature direction (aggregated over the pressures)\n");
    std::print(stream, "    pair (indices)     temperature [K]           attempts    accepted    acceptance\n");
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    for (std::size_t temperatureIndex = 0; temperatureIndex + 1 < numberOfTemperatures; ++temperatureIndex)
    {
      std::print(stream, "    {:4d} - {:<4d}   {:10.4f} - {:<10.4f}   {:9d}   {:9d}    {:8.4f} %\n", temperatureIndex,
                 temperatureIndex + 1, temperatures[temperatureIndex], temperatures[temperatureIndex + 1],
                 swapAttemptsPerTemperaturePair[temperatureIndex], swapAcceptedPerTemperaturePair[temperatureIndex],
                 100.0 * static_cast<double>(swapAcceptedPerTemperaturePair[temperatureIndex]) /
                     static_cast<double>(std::max(1uz, swapAttemptsPerTemperaturePair[temperatureIndex])));
    }
    std::print(stream, "\n");
  }
  if (numberOfPressures > 1uz)
  {
    std::print(stream, "    pressure direction (aggregated over the temperatures)\n");
    std::print(stream, "    pair (indices)     pressure [Pa]                 attempts    accepted    acceptance\n");
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    for (std::size_t pressureIndex = 0; pressureIndex + 1 < numberOfPressures; ++pressureIndex)
    {
      std::print(stream, "    {:4d} - {:<4d}   {:12.5e} - {:<12.5e}   {:9d}   {:9d}    {:8.4f} %\n", pressureIndex,
                 pressureIndex + 1, pressures[pressureIndex], pressures[pressureIndex + 1],
                 swapAttemptsPerPressurePair[pressureIndex], swapAcceptedPerPressurePair[pressureIndex],
                 100.0 * static_cast<double>(swapAcceptedPerPressurePair[pressureIndex]) /
                     static_cast<double>(std::max(1uz, swapAttemptsPerPressurePair[pressureIndex])));
    }
    std::print(stream, "\n");
  }
  std::print(stream, "\n");

  std::print(stream, "Production run CPU timings of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n", totalPreInitializationSimulationTime.count());
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n", totalSimulationTime.count());
  std::print(stream, "\n\n");
  std::flush(stream);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["hyperParallelTempering"]["sweeps"] = swapSweeps;
  outputJson["output"]["hyperParallelTempering"]["attempts"] = swapAttempts;
  outputJson["output"]["hyperParallelTempering"]["accepted"] = swapAccepted;
  outputJson["output"]["hyperParallelTempering"]["attemptsPerTemperaturePair"] = swapAttemptsPerTemperaturePair;
  outputJson["output"]["hyperParallelTempering"]["acceptedPerTemperaturePair"] = swapAcceptedPerTemperaturePair;
  outputJson["output"]["hyperParallelTempering"]["attemptsPerPressurePair"] = swapAttemptsPerPressurePair;
  outputJson["output"]["hyperParallelTempering"]["acceptedPerPressurePair"] = swapAcceptedPerPressurePair;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["total"] = totalSimulationTime.count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void HyperParallelTempering::writeIsothermSnapshot() const
{
  // one isotherm file per temperature: every replica along the pressure ladder contributes one
  // point, so the whole isotherm is assembled from the per-replica block averages
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    std::ofstream isotherm(std::format("output/isotherm_{}.hyper_parallel_tempering.txt",
                                       temperatures[temperatureIndex]),
                           std::ios::trunc);

    std::print(isotherm, "# Hyper-parallel tempering: adsorption isotherm at {} [K]\n",
               temperatures[temperatureIndex]);
    std::print(isotherm, "# one block per component (gnuplot 'index'), rows ordered by increasing pressure\n");
    std::print(isotherm, "# errors are the block confidence-interval errors (zero until three blocks are filled)\n");
    std::print(isotherm, "# column 1: fugacity [Pa]\n");
    std::print(isotherm, "# column 2, 3: absolute loading, error [molecules/cell]\n");
    std::print(isotherm, "# column 4, 5: absolute loading, error [molecules/unit-cell]\n");
    std::print(isotherm, "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
    std::print(isotherm, "# column 8, 9: absolute loading, error [mg/g-framework]\n");
    std::print(isotherm, "# column 10: pressure [Pa]\n");

    const System& front = systems[replicaIndex(temperatureIndex, 0)];
    for (std::size_t componentId = 0; componentId < front.components.size(); ++componentId)
    {
      std::print(isotherm, "\n\n# component {} ({})\n", componentId, front.components[componentId].name);

      for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
      {
        const System& system = systems[replicaIndex(temperatureIndex, pressureIndex)];
        const Component& component = system.components[componentId];

        const double fugacity =
            component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.input_pressure;

        const auto [loadingAverage, loadingError] = system.averageLoadings.average();
        const double moleculesPerCell = loadingAverage.numberOfMolecules[componentId];
        const double moleculesPerCellError = loadingError.numberOfMolecules[componentId];

        // without a framework the per-unit-cell and per-framework-mass units are not defined
        double toMoleculesPerUnitCell = 1.0;
        double toMolePerKg = 0.0;
        double toMgPerG = 0.0;
        if (system.framework.has_value())
        {
          const int3 numberOfUnitCells = system.framework->numberOfUnitCells;
          toMoleculesPerUnitCell =
              1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
          const double frameworkMass = system.frameworkMass().value();
          toMolePerKg = 1000.0 / frameworkMass;
          toMgPerG = 1000.0 * component.totalMass / frameworkMass;
        }

        std::print(isotherm,
                   "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}\n",
                   fugacity, moleculesPerCell, moleculesPerCellError, toMoleculesPerUnitCell * moleculesPerCell,
                   toMoleculesPerUnitCell * moleculesPerCellError, toMolePerKg * moleculesPerCell,
                   toMolePerKg * moleculesPerCellError, toMgPerG * moleculesPerCell,
                   toMgPerG * moleculesPerCellError, system.input_pressure);
      }
    }
  }
}

void HyperParallelTempering::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
{
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    System& system = systems[replicaId];
    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());

    std::print(replicaStream, "\n");
    std::print(replicaStream, "===============================================================================\n");
    std::print(replicaStream, "                             Simulation finished!\n");
    std::print(replicaStream, "===============================================================================\n");
    std::print(replicaStream, "\n");

    std::string status_line = std::format("Final state after {} cycles\n", numberOfProductionCycles);
    std::print(replicaStream, "{}", system.writeProductionStatusReportMC(status_line));

    const RunningEnergy drift = system.runningEnergies - recomputed[replicaId];
    replicaStream << system.runningEnergies.printMCDiff(recomputed[replicaId]);
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "Monte-Carlo moves statistics\n");
    std::print(replicaStream, "===============================================================================\n\n");
    std::print(replicaStream, "{}", system.writeMCMoveStatistics());

    std::print(replicaStream, "Production run CPU timings of the MC moves of this replica\n");
    std::print(replicaStream, "===============================================================================\n\n");
    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(replicaStream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(replicaStream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "{}",
               system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework,
                                                              system.components));
    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(replicaStream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    std::print(replicaStream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    std::flush(replicaStream);

    // final flush of the analysis-property files (cycle 0 bypasses the 'writeEvery' gate)
    writeReplicaAnalysisOutputs(system, replicaId, 0uz);

    // final flush of the time-evolution files: the total cycle count is rounded up to a multiple
    // of the writer's own 'writeEvery' so its gate passes and all collected samples are written
    if (system.propertyNumberOfMoleculesEvolution.has_value() &&
        system.propertyNumberOfMoleculesEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyNumberOfMoleculesEvolution->writeEvery.value();
      system.propertyNumberOfMoleculesEvolution->writeOutput(
          replicaId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }
    if (system.propertyVolumeEvolution.has_value() && system.propertyVolumeEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyVolumeEvolution->writeEvery.value();
      system.propertyVolumeEvolution->writeOutput(
          replicaId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }

    // per-replica json statistics
    replicaJsons[replicaId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    replicaJsons[replicaId]["output"]["recomputedEnergies"] = recomputed[replicaId].jsonMC();
    replicaJsons[replicaId]["output"]["drift"] = drift.jsonMC();
    replicaJsons[replicaId]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    replicaJsons[replicaId]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();
    for (const Component& component : system.components)
    {
      replicaJsons[replicaId]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }
    replicaJsons[replicaId]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    replicaJsons[replicaId]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();

    std::ofstream json(replicaJsonFileNames[replicaId]);
    json << replicaJsons[replicaId].dump(4);
  }
}
