module;

module parallel_thermodynamic_integration;

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
import mc_moves_lambda_exchange;
import json;

ParallelThermodynamicIntegration::ParallelThermodynamicIntegration(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      numberOfBlocks(reader.numberOfBlocks),
      numberOfLambdaBins(reader.numberOfLambdaBins),
      lambdaExchangeEvery(reader.lambdaExchangeEvery)
{
  if (numberOfLambdaBins < 2)
  {
    throw std::runtime_error("[ParallelThermodynamicIntegration]: 'NumberOfLambdaBins' must be at least 2\n");
  }

  // the single declared system is replicated into one replica per lambda-bin
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  for (std::size_t componentId = 0; componentId < templateSystem.components.size(); ++componentId)
  {
    if (templateSystem.components[componentId].fixedLambdaBin.has_value())
    {
      tiComponentId = componentId;
    }
  }

  systems.reserve(numberOfLambdaBins);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfLambdaBins; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  // replica k starts pinned at lambda-bin k, with its own random-number stream
  randoms.reserve(numberOfLambdaBins);
  for (std::size_t replicaId = 0; replicaId < numberOfLambdaBins; ++replicaId)
  {
    systems[replicaId].components[tiComponentId].fixedLambdaBin = replicaId;
    systems[replicaId].initializeFixedLambdaFractionalMolecules();
    // re-pinning changed the scaling of the fractional molecule(s); refresh the effective
    // pseudo-atom counts used by the aggregated tail-correction accounting
    systems[replicaId].computeTailCorrectionCounts();

    randoms.emplace_back(random.seed + replicaId + 1);
  }

  stepsPerReplica.assign(numberOfLambdaBins, 0uz);
  exchangeAttemptsPerPair.assign(numberOfLambdaBins - 1uz, 0uz);
  exchangeAcceptedPerPair.assign(numberOfLambdaBins - 1uz, 0uz);
}

void ParallelThermodynamicIntegration::run()
{
  setup();
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ParallelThermodynamicIntegration::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    // every replica carries its own fractional molecule(s), pinned at its lambda
    system.containsTheFractionalMolecule = true;
  }

  std::filesystem::create_directories("output");
  const System& front = systems.front();
  std::string fileNameString =
      std::format("output/output_{}_{}.parallel_ti.txt", front.temperature, front.input_pressure);
  stream.open(fileNameString, std::ios::out);
  outputJsonFileName = std::format("output/output_{}_{}.parallel_ti.json", front.temperature, front.input_pressure);
  dudlambdaFileName = std::format("output/dudlambda_{}_{}.parallel_ti.txt", front.temperature, front.input_pressure);

  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  std::print(stream, "Parallel thermodynamic integration\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Number of lambda-bins / replicas / threads: {}\n", numberOfLambdaBins);
  std::print(stream, "Pinned component:                           {} ({})\n", tiComponentId,
             front.components[tiComponentId].name);
  if (lambdaExchangeEvery == 0uz)
  {
    std::print(stream, "Lambda-exchange:                            disabled\n\n");
  }
  else
  {
    std::print(stream, "Lambda-exchange sweep every:                {} cycles\n\n", lambdaExchangeEvery);
  }

  std::print(stream, "{}", front.writeSystemStatus());
  std::print(stream, "{}", front.forceField.printPseudoAtomStatus());
  std::print(stream, "{}", front.forceField.printForceFieldStatus());
  std::print(stream, "{}", front.writeComponentStatus());
  std::print(stream, "{}", front.writeNumberOfPseudoAtoms());

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["initialConditions"] = front.jsonSystemStatus();
  outputJson["initialization"]["forceField"] = front.forceField.jsonForceFieldStatus();
  outputJson["initialization"]["forceField"]["pseudoAtoms"] = front.forceField.jsonPseudoAtomStatus();
  outputJson["initialization"]["components"] = front.jsonComponentStatus();
  outputJson["initialization"]["numberOfLambdaBins"] = numberOfLambdaBins;
  outputJson["initialization"]["lambdaExchangeEvery"] = lambdaExchangeEvery;

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-replica output files: each worker thread writes exclusively to its own stream
  replicaStreams.reserve(numberOfLambdaBins);
  replicaJsonFileNames.reserve(numberOfLambdaBins);
  replicaJsons.resize(numberOfLambdaBins);
  for (std::size_t replicaId = 0; replicaId < numberOfLambdaBins; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(
        std::format("output/output_{}_{}.parallel_ti.r{}.txt", system.temperature, system.input_pressure, replicaId),
        std::ios::out);
    replicaJsonFileNames.emplace_back(
        std::format("output/output_{}_{}.parallel_ti.r{}.json", system.temperature, system.input_pressure, replicaId));

    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "Parallel thermodynamic integration: replica {} of {} (starting at lambda-bin {})\n",
               replicaId, numberOfLambdaBins, replicaId);
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
    replicaJsons[replicaId]["startingLambdaBin"] = replicaId;
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

void ParallelThermodynamicIntegration::performReplicaCycle(std::size_t replicaId, SimulationStage stage,
                                                           std::size_t currentBlock)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];

  // every replica owns its own fractional molecule(s); the Gibbs-style moves that could relocate
  // the fractional molecule between systems are not enabled in this driver
  std::size_t fractionalMoleculeSystem = 0uz;

  const std::size_t numberOfStepsPerCycle =
      std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();

  for (std::size_t j = 0uz; j != numberOfStepsPerCycle; ++j)
  {
    std::size_t selectedComponent = system.randomComponent(rng);

    // the pinned lambda-bin only changes through the lambda-exchange between replicas: no
    // lambda-changing moves are enabled and no Wang-Landau biasing is needed
    switch (stage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Equilibration:
        MC_Moves::performRandomMoveEquilibration(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem,
                                              currentBlock);
        ++stepsPerReplica[replicaId];
        break;
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
    // pinned pair- and group-swap lambda histograms
    system.pairSwapLambdaSampleOccupancy();
  }
}

void ParallelThermodynamicIntegration::performLambdaExchangeSweep(SimulationStage stage,
                                                                  std::size_t numberOfCycles) noexcept
{
  // exchanges migrate the lambda values over the replicas: pair the replicas that currently hold
  // neighboring lambda-bins, alternating the pairing offset between sweeps
  std::vector<std::size_t> replicaAtBin(numberOfLambdaBins);
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    replicaAtBin[systems[replicaId].components[tiComponentId].fixedLambdaBin.value()] = replicaId;
  }

  const std::size_t offset = exchangeSweeps % 2uz;
  for (std::size_t bin = offset; bin + 1 < numberOfLambdaBins; bin += 2uz)
  {
    ++exchangeAttempts;
    ++exchangeAttemptsPerPair[bin];
    if (MC_Moves::LambdaExchange(random, systems[replicaAtBin[bin]], systems[replicaAtBin[bin + 1]]).has_value())
    {
      ++exchangeAccepted;
      ++exchangeAcceptedPerPair[bin];
    }
  }
  ++exchangeSweeps;
  ++sweepsThisStage;

  // progress report: all worker threads are parked on the barrier, so this is race-free
  const std::size_t printSweeps = std::max(1uz, printEvery / std::max(1uz, lambdaExchangeEvery));
  if (sweepsThisStage % printSweeps == 0uz)
  {
    const std::string_view stageName = (stage == SimulationStage::Initialization) ? "initialization"
                                       : (stage == SimulationStage::Equilibration) ? "equilibration"
                                                                                   : "production";
    const std::size_t cycle = std::min(sweepsThisStage * lambdaExchangeEvery, numberOfCycles);
    {
      std::scoped_lock lock(outputMutex);
      std::print(stream, "Lambda-exchange sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%)\n", exchangeSweeps,
                 stageName, cycle, numberOfCycles, exchangeAccepted, exchangeAttempts,
                 100.0 * static_cast<double>(exchangeAccepted) / static_cast<double>(std::max(1uz, exchangeAttempts)));
      std::flush(stream);
    }

    if (stage == SimulationStage::Production)
    {
      writeDUdlambdaSnapshot(cycle, numberOfCycles);
    }
  }
}

void ParallelThermodynamicIntegration::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  simulationStage = stage;

  // serial per-stage preparation
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

        // start production with fresh dU/dlambda accumulators; the occupied bin is kept
        component.lambdaGC.clear();
      }
      // pinned pair- and group-swap lambda histograms
      system.pairSwapLambdaClearBookkeeping();
    }
    std::fill(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);
  }
  sweepsThisStage = 0uz;

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::Initialization) ? "Initialization"
                                       : (stage == SimulationStage::Equilibration) ? "Equilibration"
                                                                                   : "Production";
    std::print(stream, "\n{} stage: {} cycles on {} replicas/threads\n", stageName, numberOfCycles,
               numberOfLambdaBins);
    std::flush(stream);
  }

  // one worker thread per replica; the barrier completion performs the lambda-exchange sweeps
  auto onAllArrived = [this, stage, numberOfCycles]() noexcept { performLambdaExchangeSweep(stage, numberOfCycles); };
  std::barrier synchronizationPoint(static_cast<std::ptrdiff_t>(numberOfLambdaBins), onAllArrived);

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfLambdaBins);
    for (std::size_t replicaId = 0; replicaId < numberOfLambdaBins; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, &synchronizationPoint]()
          {
            System& system = systems[replicaId];

            // each thread computes the total energies of its own replica
            if (stage == SimulationStage::Initialization)
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

              if (cycle % printEvery == 0uz)
              {
                // each thread writes exclusively to its own replica stream: no locking needed
                system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                              system.simulationBox);

                std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
                switch (stage)
                {
                  case SimulationStage::Initialization:
                    std::print(replicaStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(replicaStream, "{}", system.writeEquilibrationStatusReportMC(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                  {
                    std::string status_line =
                        std::format("Current cycle: {} out of {} (currently at lambda-bin {})\n", cycle,
                                    numberOfCycles, system.components[tiComponentId].fixedLambdaBin.value());
                    std::print(replicaStream, "{}", system.writeProductionStatusReportMC(status_line));
                    break;
                  }
                  default:
                    break;
                }
                std::flush(replicaStream);
              }

              // the only synchronization point between the threads: the lambda-exchange sweep
              if (lambdaExchangeEvery != 0uz && (cycle + 1uz) % lambdaExchangeEvery == 0uz)
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

void ParallelThermodynamicIntegration::writeDUdlambdaSnapshot(std::size_t cycle, std::size_t numberOfCycles) const
{
  const PropertyLambdaProbabilityHistogram stitched = stitchedHistogram();
  const std::vector<double> average = stitched.averagedDUdlambda();

  // Mid-run only the blocks reached so far contain samples; the confidence intervals are computed
  // over the filled blocks only (empty blocks would enter the block analysis as spurious zeros).
  // With fewer than three filled blocks the error estimate is not yet defined and reported as zero.
  std::vector<std::vector<double>> filledBlockAverages;
  std::vector<double> filledBlockIntegrals;
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    double blockSampleCount = 0.0;
    for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
    {
      blockSampleCount += stitched.bookKeepingDUdlambda[blockIndex][bin].second;
    }
    if (blockSampleCount > 0.0)
    {
      filledBlockAverages.push_back(stitched.averagedDUdlambda(blockIndex));
      filledBlockIntegrals.push_back(cubicSplineIntegral(filledBlockAverages.back(), stitched.delta));
    }
  }

  const std::vector<double> errors = blockErrorEstimate(filledBlockAverages, average);
  const double splineIntegral = cubicSplineIntegral(average, stitched.delta);
  const double splineError = blockErrorEstimate<double>(filledBlockIntegrals, splineIntegral);
  const double simpsonIntegral = stitched.simpsonIntegral(average);
  const double trapezoidalIntegral = trapezoidIntegral(average, stitched.delta);

  std::ofstream snapshot(dudlambdaFileName, std::ios::trunc);
  std::print(snapshot, "# Parallel thermodynamic integration: stitched <dU/dlambda>(lambda) curve\n");
  std::print(snapshot, "# component {} ({})\n", tiComponentId, systems.front().components[tiComponentId].name);
  std::print(snapshot, "# production cycle {} of {}\n", cycle, numberOfCycles);
  std::print(snapshot, "# blocks filled: {} of {} (confidence intervals use the filled blocks only)\n",
             filledBlockAverages.size(), numberOfBlocks);
  std::print(snapshot, "# excess chemical potential (spline):    {: .6e} +/- {: .6e} [K]\n",
             Units::EnergyToKelvin * splineIntegral, Units::EnergyToKelvin * splineError);
  std::print(snapshot, "# excess chemical potential (spline):    {: .6e} +/- {: .6e} [kJ/mol]\n",
             Units::EnergyToKJPerMol * splineIntegral, Units::EnergyToKJPerMol * splineError);
  std::print(snapshot, "# excess chemical potential (Simpson):   {: .6e} [K]\n",
             Units::EnergyToKelvin * simpsonIntegral);
  std::print(snapshot, "# excess chemical potential (trapezoid): {: .6e} [K]\n",
             Units::EnergyToKelvin * trapezoidalIntegral);
  std::print(snapshot, "# column 1: lambda\n");
  std::print(snapshot, "# column 2: <dU/dlambda> [K]\n");
  std::print(snapshot, "# column 3: confidence-interval error [K]\n");
  for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
  {
    std::print(snapshot, "{:.6f} {: .8e} {: .8e}\n", static_cast<double>(bin) * stitched.delta,
               Units::EnergyToKelvin * average[bin], Units::EnergyToKelvin * errors[bin]);
  }
}

void ParallelThermodynamicIntegration::output()
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

  stream << systems.front().runningEnergies.printMCDiff(recomputed.front());
  std::print(stream, "\n");
  std::print(stream, "Energy drift per replica\n");
  std::print(stream, "===============================================================================\n\n");
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    const RunningEnergy drift = systems[replicaId].runningEnergies - recomputed[replicaId];
    std::print(stream, "    replica {:4d} (lambda-bin {:4d}): drift {: .6e} [K]\n", replicaId,
               systems[replicaId].components[tiComponentId].fixedLambdaBin.value(),
               Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  std::print(stream, "Lambda-exchange statistics\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    sweeps:    {}\n", exchangeSweeps);
  std::print(stream, "    attempts:  {}\n", exchangeAttempts);
  std::print(stream, "    accepted:  {} ({:.4f} %)\n\n", exchangeAccepted,
             100.0 * static_cast<double>(exchangeAccepted) / static_cast<double>(std::max(1uz, exchangeAttempts)));

  // low acceptance for a particular pair marks a bottleneck in the lambda ladder (replicas cannot
  // migrate past it); consider more lambda-bins if a pair falls well below its neighbors
  const double deltaLambda = 1.0 / static_cast<double>(numberOfLambdaBins - 1);
  std::print(stream, "    pair (bins)      lambda            attempts    accepted    acceptance\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t bin = 0; bin + 1 < numberOfLambdaBins; ++bin)
  {
    std::print(stream, "    {:4d} - {:<4d}   {:.5f} - {:.5f}   {:9d}   {:9d}    {:8.4f} %\n", bin, bin + 1,
               static_cast<double>(bin) * deltaLambda, static_cast<double>(bin + 1) * deltaLambda,
               exchangeAttemptsPerPair[bin], exchangeAcceptedPerPair[bin],
               100.0 * static_cast<double>(exchangeAcceptedPerPair[bin]) /
                   static_cast<double>(std::max(1uz, exchangeAttemptsPerPair[bin])));
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run CPU timings of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n", totalSimulationTime.count());
  std::print(stream, "\n\n");

  std::print(stream, "{}", writeStitchedThermodynamicIntegration());
  std::flush(stream);

  // final snapshot of the curve (also covers runs with the lambda-exchange disabled)
  writeDUdlambdaSnapshot(numberOfProductionCycles, numberOfProductionCycles);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["lambdaExchange"]["sweeps"] = exchangeSweeps;
  outputJson["output"]["lambdaExchange"]["attempts"] = exchangeAttempts;
  outputJson["output"]["lambdaExchange"]["accepted"] = exchangeAccepted;
  outputJson["output"]["lambdaExchange"]["attemptsPerPair"] = exchangeAttemptsPerPair;
  outputJson["output"]["lambdaExchange"]["acceptedPerPair"] = exchangeAcceptedPerPair;
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["total"] = totalSimulationTime.count();
  outputJson["properties"]["thermodynamicIntegration"] = jsonStitchedThermodynamicIntegration();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ParallelThermodynamicIntegration::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
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

    std::string status_line =
        std::format("Final state after {} cycles (ended at lambda-bin {})\n", numberOfProductionCycles,
                    system.components[tiComponentId].fixedLambdaBin.value());
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

    std::print(replicaStream, "{}", writeReplicaThermodynamicIntegration(replicaId));
    std::flush(replicaStream);

    // per-replica json statistics
    replicaJsons[replicaId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    replicaJsons[replicaId]["output"]["recomputedEnergies"] = recomputed[replicaId].jsonMC();
    replicaJsons[replicaId]["output"]["drift"] = drift.jsonMC();
    replicaJsons[replicaId]["output"]["finalLambdaBin"] = system.components[tiComponentId].fixedLambdaBin.value();
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

    // the replica's own per-bin dU/dlambda book-keeping (samples over the bins it visited)
    const PropertyLambdaProbabilityHistogram& histogram = system.components[tiComponentId].fixedLambdaHistogram();
    std::pair<std::vector<double>, std::vector<double>> dudlambda = histogram.averageDuDlambda();
    std::vector<double> sampleCounts(numberOfLambdaBins, 0.0);
    for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
    {
      for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        sampleCounts[bin] += histogram.bookKeepingDUdlambda[blockIndex][bin].second;
      }
    }
    std::vector<double> averagesInKelvin(numberOfLambdaBins);
    for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
    {
      averagesInKelvin[bin] = Units::EnergyToKelvin * dudlambda.first[bin];
    }
    replicaJsons[replicaId]["properties"]["thermodynamicIntegration"]["samplesPerLambdaBin"] = sampleCounts;
    replicaJsons[replicaId]["properties"]["thermodynamicIntegration"]["averageDUdlambda[K]"] = averagesInKelvin;

    std::ofstream json(replicaJsonFileNames[replicaId]);
    json << replicaJsons[replicaId].dump(4);
  }
}

std::string ParallelThermodynamicIntegration::writeReplicaThermodynamicIntegration(std::size_t replicaId) const
{
  std::ostringstream outStream;

  const System& system = systems[replicaId];
  const Component& component = system.components[tiComponentId];
  const PropertyLambdaProbabilityHistogram& histogram = component.fixedLambdaHistogram();

  std::pair<std::vector<double>, std::vector<double>> dudlambda = histogram.averageDuDlambda();

  std::print(outStream, "Thermodynamic-integration book-keeping of this replica\n");
  std::print(outStream, "===============================================================================\n\n");
  std::print(outStream, "component {} ({})\n", tiComponentId, component.name);
  std::print(outStream, "current lambda-bin: {}\n\n", component.fixedLambdaBin.value());
  std::print(outStream, "    bin   lambda     samples        <dU/dlambda> [K]\n");
  std::print(outStream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
  {
    double sampleCount = 0.0;
    for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      sampleCount += histogram.bookKeepingDUdlambda[blockIndex][bin].second;
    }
    if (sampleCount > 0.0)
    {
      std::print(outStream, "    {:3d}   {:.5f}   {:12.0f}   {: .8e}\n", bin,
                 static_cast<double>(bin) * histogram.delta, sampleCount,
                 Units::EnergyToKelvin * dudlambda.first[bin]);
    }
  }
  std::print(outStream, "\n");
  std::print(outStream,
             "(samples in bins other than the starting bin {} were collected after accepted lambda-exchanges;\n"
             " the stitched curve over all replicas is reported in the combined output file)\n\n\n",
             replicaId);

  return outStream.str();
}

PropertyLambdaProbabilityHistogram ParallelThermodynamicIntegration::stitchedHistogram() const
{
  PropertyLambdaProbabilityHistogram stitched(numberOfBlocks, numberOfLambdaBins);

  for (const System& system : systems)
  {
    const PropertyLambdaProbabilityHistogram& histogram =
        system.components[tiComponentId].fixedLambdaHistogram();
    for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
      {
        stitched.bookKeepingDUdlambda[blockIndex][bin].first +=
            histogram.bookKeepingDUdlambda[blockIndex][bin].first;
        stitched.bookKeepingDUdlambda[blockIndex][bin].second +=
            histogram.bookKeepingDUdlambda[blockIndex][bin].second;
      }
    }
  }

  return stitched;
}

// Integral over the full data range of the natural cubic spline through equidistant points with
// spacing 'delta'. The spline second derivatives follow from the standard tridiagonal system
// (natural boundary conditions), solved with the Thomas algorithm; per interval the exact integral
// is delta/2 (y_i + y_{i+1}) - delta^3/24 (M_i + M_{i+1}).
double ParallelThermodynamicIntegration::cubicSplineIntegral(const std::vector<double>& data, double delta)
{
  const std::size_t n = data.size();
  if (n < 2uz) return 0.0;
  if (n == 2uz) return 0.5 * delta * (data[0] + data[1]);

  std::vector<double> secondDerivative(n, 0.0);
  {
    // interior equations: (1/6) M[i-1] + (2/3) M[i] + (1/6) M[i+1] = (y[i+1] - 2 y[i] + y[i-1]) / delta^2
    const std::size_t interior = n - 2uz;
    std::vector<double> superDiagonal(interior, 0.0);
    std::vector<double> rhs(interior, 0.0);

    double pivot = 2.0 / 3.0;
    rhs[0] = (data[2] - 2.0 * data[1] + data[0]) / (delta * delta) / pivot;
    for (std::size_t i = 1; i < interior; ++i)
    {
      superDiagonal[i - 1] = (1.0 / 6.0) / pivot;
      pivot = 2.0 / 3.0 - (1.0 / 6.0) * superDiagonal[i - 1];
      rhs[i] = ((data[i + 2] - 2.0 * data[i + 1] + data[i]) / (delta * delta) - (1.0 / 6.0) * rhs[i - 1]) / pivot;
    }
    secondDerivative[interior] = rhs[interior - 1];
    for (std::size_t i = interior - 1; i > 0; --i)
    {
      secondDerivative[i] = rhs[i - 1] - superDiagonal[i - 1] * secondDerivative[i + 1];
    }
  }

  double integral = 0.0;
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    integral += 0.5 * delta * (data[i] + data[i + 1]) -
                (delta * delta * delta / 24.0) * (secondDerivative[i] + secondDerivative[i + 1]);
  }
  return integral;
}

double ParallelThermodynamicIntegration::trapezoidIntegral(const std::vector<double>& data, double delta)
{
  const std::size_t n = data.size();
  if (n < 2uz) return 0.0;

  double integral = 0.0;
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    integral += 0.5 * delta * (data[i] + data[i + 1]);
  }
  return integral;
}

std::string ParallelThermodynamicIntegration::writeStitchedThermodynamicIntegration() const
{
  std::ostringstream outStream;

  const PropertyLambdaProbabilityHistogram stitched = stitchedHistogram();
  const Component& component = systems.front().components[tiComponentId];

  std::print(outStream, "Parallel thermodynamic integration: stitched <dU/dlambda>(lambda) curve\n");
  std::print(outStream, "===============================================================================\n\n");

  std::print(outStream, "component {} ({})\n\n", tiComponentId, component.name);

  std::pair<std::vector<double>, std::vector<double>> dudlambda = stitched.averageDuDlambda();

  std::print(outStream, "    bin   lambda     <dU/dlambda> [K]       +/-    [K]\n");
  std::print(outStream, "    ---------------------------------------------------------------------------\n");
  for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
  {
    std::print(outStream, "    {:3d}   {:.5f}   {: .8e}   {: .8e}\n", bin,
               static_cast<double>(bin) * stitched.delta, Units::EnergyToKelvin * dudlambda.first[bin],
               Units::EnergyToKelvin * dudlambda.second[bin]);
  }
  std::print(outStream, "\n");

  // block-wise integrals give the confidence interval of each quadrature rule. The data live on
  // fixed equidistant lambda-bins, so the applicable rules are the natural-spline integral and the
  // Newton-Cotes rules (Gaussian quadrature would require samples at the non-equidistant Gauss
  // nodes). The three estimates agreeing within their error bars signals that the grid resolves
  // the curvature of the curve; the spline value is the recommended one.
  std::vector<double> blockSplineIntegrals(numberOfBlocks);
  std::vector<double> blockSimpsonIntegrals(numberOfBlocks);
  std::vector<double> blockTrapezoidIntegrals(numberOfBlocks);
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    const std::vector<double> blockAverage = stitched.averagedDUdlambda(blockIndex);
    blockSplineIntegrals[blockIndex] = cubicSplineIntegral(blockAverage, stitched.delta);
    blockSimpsonIntegrals[blockIndex] = stitched.simpsonIntegral(blockAverage);
    blockTrapezoidIntegrals[blockIndex] = trapezoidIntegral(blockAverage, stitched.delta);
    std::print(outStream, "        Block[ {:2d}] integral (spline): {: .6e} [K]\n", blockIndex,
               Units::EnergyToKelvin * blockSplineIntegrals[blockIndex]);
  }
  const double splineIntegral = cubicSplineIntegral(dudlambda.first, stitched.delta);
  const double splineError = blockErrorEstimate<double>(blockSplineIntegrals, splineIntegral);
  const double simpsonIntegral = stitched.simpsonIntegral(dudlambda.first);
  const double simpsonError = blockErrorEstimate<double>(blockSimpsonIntegrals, simpsonIntegral);
  const double trapezoidalIntegral = trapezoidIntegral(dudlambda.first, stitched.delta);
  const double trapezoidError = blockErrorEstimate<double>(blockTrapezoidIntegrals, trapezoidalIntegral);

  std::print(outStream, "    ---------------------------------------------------------------------------\n");
  std::print(outStream, "    excess chemical potential\n");
  std::print(outStream, "        spline:     {: .6e} +/- {: .6e} [K]\n", Units::EnergyToKelvin * splineIntegral,
             Units::EnergyToKelvin * splineError);
  std::print(outStream, "        Simpson:    {: .6e} +/- {: .6e} [K]\n", Units::EnergyToKelvin * simpsonIntegral,
             Units::EnergyToKelvin * simpsonError);
  std::print(outStream, "        trapezoid:  {: .6e} +/- {: .6e} [K]\n", Units::EnergyToKelvin * trapezoidalIntegral,
             Units::EnergyToKelvin * trapezoidError);
  std::print(outStream, "        spline:     {: .6e} +/- {: .6e} [kJ/mol]\n", Units::EnergyToKJPerMol * splineIntegral,
             Units::EnergyToKJPerMol * splineError);
  std::print(outStream, "        Simpson:    {: .6e} +/- {: .6e} [kJ/mol]\n",
             Units::EnergyToKJPerMol * simpsonIntegral, Units::EnergyToKJPerMol * simpsonError);
  std::print(outStream, "        trapezoid:  {: .6e} +/- {: .6e} [kJ/mol]\n",
             Units::EnergyToKJPerMol * trapezoidalIntegral, Units::EnergyToKJPerMol * trapezoidError);
  std::print(outStream,
             "    (quadrature spread far below the sampling error confirms the lambda-grid is fine enough)\n");
  std::print(outStream, "\n\n");

  return outStream.str();
}

nlohmann::json ParallelThermodynamicIntegration::jsonStitchedThermodynamicIntegration() const
{
  nlohmann::json status;

  const PropertyLambdaProbabilityHistogram stitched = stitchedHistogram();
  const Component& component = systems.front().components[tiComponentId];

  std::pair<std::vector<double>, std::vector<double>> dudlambda = stitched.averageDuDlambda();

  std::vector<double> lambdas(numberOfLambdaBins);
  std::vector<double> averagesInKelvin(numberOfLambdaBins);
  std::vector<double> errorsInKelvin(numberOfLambdaBins);
  for (std::size_t bin = 0; bin < numberOfLambdaBins; ++bin)
  {
    lambdas[bin] = static_cast<double>(bin) * stitched.delta;
    averagesInKelvin[bin] = Units::EnergyToKelvin * dudlambda.first[bin];
    errorsInKelvin[bin] = Units::EnergyToKelvin * dudlambda.second[bin];
  }

  std::vector<double> blockSplineIntegrals(numberOfBlocks);
  std::vector<double> blockSimpsonIntegrals(numberOfBlocks);
  std::vector<double> blockTrapezoidIntegrals(numberOfBlocks);
  std::vector<double> blockSplineIntegralsInKelvin(numberOfBlocks);
  for (std::size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    const std::vector<double> blockAverage = stitched.averagedDUdlambda(blockIndex);
    blockSplineIntegrals[blockIndex] = cubicSplineIntegral(blockAverage, stitched.delta);
    blockSimpsonIntegrals[blockIndex] = stitched.simpsonIntegral(blockAverage);
    blockTrapezoidIntegrals[blockIndex] = trapezoidIntegral(blockAverage, stitched.delta);
    blockSplineIntegralsInKelvin[blockIndex] = Units::EnergyToKelvin * blockSplineIntegrals[blockIndex];
  }
  const double splineIntegral = cubicSplineIntegral(dudlambda.first, stitched.delta);
  const double splineError = blockErrorEstimate<double>(blockSplineIntegrals, splineIntegral);
  const double simpsonIntegral = stitched.simpsonIntegral(dudlambda.first);
  const double simpsonError = blockErrorEstimate<double>(blockSimpsonIntegrals, simpsonIntegral);
  const double trapezoidalIntegral = trapezoidIntegral(dudlambda.first, stitched.delta);
  const double trapezoidError = blockErrorEstimate<double>(blockTrapezoidIntegrals, trapezoidalIntegral);

  status[component.name]["numberOfLambdaBins"] = numberOfLambdaBins;
  status[component.name]["lambda"] = lambdas;
  status[component.name]["averageDUdlambda[K]"] = averagesInKelvin;
  status[component.name]["confidenceDUdlambda[K]"] = errorsInKelvin;
  status[component.name]["blockIntegralsSpline[K]"] = blockSplineIntegralsInKelvin;
  status[component.name]["excessChemicalPotentialSpline[K]"] = Units::EnergyToKelvin * splineIntegral;
  status[component.name]["confidenceExcessChemicalPotentialSpline[K]"] = Units::EnergyToKelvin * splineError;
  status[component.name]["excessChemicalPotentialSpline[kJ/mol]"] = Units::EnergyToKJPerMol * splineIntegral;
  status[component.name]["confidenceExcessChemicalPotentialSpline[kJ/mol]"] = Units::EnergyToKJPerMol * splineError;
  status[component.name]["excessChemicalPotentialSimpson[K]"] = Units::EnergyToKelvin * simpsonIntegral;
  status[component.name]["confidenceExcessChemicalPotentialSimpson[K]"] = Units::EnergyToKelvin * simpsonError;
  status[component.name]["excessChemicalPotentialSimpson[kJ/mol]"] = Units::EnergyToKJPerMol * simpsonIntegral;
  status[component.name]["confidenceExcessChemicalPotentialSimpson[kJ/mol]"] = Units::EnergyToKJPerMol * simpsonError;
  status[component.name]["excessChemicalPotentialTrapezoid[K]"] = Units::EnergyToKelvin * trapezoidalIntegral;
  status[component.name]["confidenceExcessChemicalPotentialTrapezoid[K]"] = Units::EnergyToKelvin * trapezoidError;
  status[component.name]["excessChemicalPotentialTrapezoid[kJ/mol]"] = Units::EnergyToKJPerMol * trapezoidalIntegral;
  status[component.name]["confidenceExcessChemicalPotentialTrapezoid[kJ/mol]"] =
      Units::EnergyToKJPerMol * trapezoidError;

  return status;
}
