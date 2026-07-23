module;

module parallel_tempering;

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
import mc_moves_parallel_tempering_swap;
import json;

ParallelTempering::ParallelTempering(InputReader& reader)
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
      numberOfReplicas(temperatures.size())
{
  // the single declared system is replicated into one replica per temperature of the ladder
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  systems.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  // replica k is pinned at temperature T_k, with its own random-number stream
  randoms.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    System& system = systems[replicaId];
    const double T = temperatures[replicaId];

    system.temperature = T;
    system.beta = 1.0 / (Units::KB * T);

    // temperature-dependent potentials (Feynman-Hibbs) derive pair coefficients, shifts and
    // tail-corrections from the temperature
    if (system.forceField.temperature != T)
    {
      system.forceField.temperature = T;
      system.forceField.preComputeDerivedParameters();
      system.forceField.preComputePotentialShift();
      system.forceField.preComputeTailCorrection();
    }

    // the CBMC ideal-gas conformation reservoirs are Boltzmann samples at the system temperature
    system.buildConformationReservoirs();

    randoms.emplace_back(random.seed + replicaId + 1);
  }

  stepsPerReplica.assign(numberOfReplicas, 0uz);
  swapAttemptsPerPair.assign(numberOfReplicas - 1uz, 0uz);
  swapAcceptedPerPair.assign(numberOfReplicas - 1uz, 0uz);
}

void ParallelTempering::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ParallelTempering::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  stream.open("output/output.parallel_tempering.txt", std::ios::out);
  outputJsonFileName = "output/output.parallel_tempering.json";

  const System& front = systems.front();
  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  std::print(stream, "Parallel tempering\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Number of temperatures / replicas / threads: {}\n", numberOfReplicas);
  std::print(stream, "Temperature ladder:                         ");
  for (double T : temperatures)
  {
    std::print(stream, " {}", T);
  }
  std::print(stream, " [K]\n");
  if (parallelTemperingSwapEvery == 0uz)
  {
    std::print(stream, "Configuration swaps:                         disabled\n\n");
  }
  else
  {
    std::print(stream, "Configuration-swap sweep every:              {} cycles\n\n", parallelTemperingSwapEvery);
  }

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["temperatures"] = temperatures;
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
    replicaStreams.emplace_back(std::format("output/output_{}_{}.parallel_tempering.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                std::ios::out);
    replicaJsonFileNames.emplace_back(std::format("output/output_{}_{}.parallel_tempering.r{}.json",
                                                  system.temperature, system.input_pressure, replicaId));

    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "Parallel tempering: replica {} of {} (temperature {} [K])\n", replicaId,
               numberOfReplicas, system.temperature);
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

void ParallelTempering::performReplicaCycle(std::size_t replicaId, SimulationStage stage, std::size_t currentBlock)
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

void ParallelTempering::performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept
{
  // the replicas keep their temperatures; the configurations migrate through the ladder.
  // alternate the pairing offset between sweeps so configurations can traverse the whole ladder
  const std::size_t offset = swapSweeps % 2uz;
  for (std::size_t replicaId = offset; replicaId + 1 < numberOfReplicas; replicaId += 2uz)
  {
    ++swapAttempts;
    ++swapAttemptsPerPair[replicaId];
    if (MC_Moves::ParallelTemperingSwap(random, systems[replicaId], systems[replicaId + 1]).has_value())
    {
      ++swapAccepted;
      ++swapAcceptedPerPair[replicaId];
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
    std::scoped_lock lock(outputMutex);
    std::print(stream, "Parallel-tempering sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%)\n", swapSweeps,
               stageName, cycle, numberOfCycles, swapAccepted, swapAttempts,
               100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));
    std::flush(stream);
  }
}

void ParallelTempering::runStage(SimulationStage stage, std::size_t numberOfCycles)
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

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, &synchronizationPoint]()
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

              if (stage == SimulationStage::Production)
              {
                system.sampleProperties(replicaId, estimation.currentBin, cycle);

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

void ParallelTempering::output()
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

  std::print(stream, "Energy drift per replica\n");
  std::print(stream, "===============================================================================\n\n");
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    const RunningEnergy drift = systems[replicaId].runningEnergies - recomputed[replicaId];
    std::print(stream, "    replica {:4d} (temperature {:10.4f} [K]): drift {: .6e} [K]\n", replicaId,
               systems[replicaId].temperature, Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  std::print(stream, "Parallel-tempering swap statistics\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    sweeps:    {}\n", swapSweeps);
  std::print(stream, "    attempts:  {}\n", swapAttempts);
  std::print(stream, "    accepted:  {} ({:.4f} %)\n\n", swapAccepted,
             100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));

  // low acceptance for a particular pair marks a bottleneck in the temperature ladder
  // (configurations cannot migrate past it); consider a denser ladder around such a pair
  std::print(stream, "    pair (replicas)    temperature [K]           attempts    accepted    acceptance\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    std::print(stream, "    {:4d} - {:<4d}   {:10.4f} - {:<10.4f}   {:9d}   {:9d}    {:8.4f} %\n", replicaId,
               replicaId + 1, temperatures[replicaId], temperatures[replicaId + 1], swapAttemptsPerPair[replicaId],
               swapAcceptedPerPair[replicaId],
               100.0 * static_cast<double>(swapAcceptedPerPair[replicaId]) /
                   static_cast<double>(std::max(1uz, swapAttemptsPerPair[replicaId])));
  }
  std::print(stream, "\n\n");

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
  outputJson["output"]["parallelTempering"]["sweeps"] = swapSweeps;
  outputJson["output"]["parallelTempering"]["attempts"] = swapAttempts;
  outputJson["output"]["parallelTempering"]["accepted"] = swapAccepted;
  outputJson["output"]["parallelTempering"]["attemptsPerPair"] = swapAttemptsPerPair;
  outputJson["output"]["parallelTempering"]["acceptedPerPair"] = swapAcceptedPerPair;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["total"] = totalSimulationTime.count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ParallelTempering::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
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
