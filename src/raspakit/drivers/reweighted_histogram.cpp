module;

module reweighted_histogram;

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
import vdwparameters;
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

ReweightedHistogram::ReweightedHistogram(InputReader& reader)
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
      sampleReweightingEvery(std::max(1uz, reader.sampleReweightingEvery)),
      reweightingTemperatures(reader.reweightingTemperatures),
      reweightingNumberOfPressures(reader.reweightingNumberOfPressures),
      temperatures(reader.parallelTemperingTemperatures),
      pressures(reader.parallelTemperingPressures),
      numberOfTemperatures(temperatures.size()),
      numberOfPressures(pressures.size()),
      numberOfReplicas(numberOfTemperatures * numberOfPressures)
{
  // the reweighting evaluates one sampled potential energy U(x) at all the temperatures, so the
  // force field must not depend on the temperature
  for (const VDWParameters& parameter : reader.systems.front().forceField.data)
  {
    if (parameter.type == VDWParameters::Type::FeynmannHibbs)
    {
      throw std::runtime_error(
          std::format("[ReweightedHistogram]: temperature-dependent potentials (Feynman-Hibbs) are not "
                      "compatible with histogram reweighting (a sampled energy U(x) is reused at all "
                      "temperatures)\n"));
    }
  }

  // the reweighted isotherms default to the simulated temperatures and the simulated pressure span
  if (reweightingTemperatures.empty())
  {
    reweightingTemperatures = temperatures;
  }
  reweightingPressureRange =
      reader.reweightingPressureRange.value_or(std::make_pair(pressures.front(), pressures.back()));

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

      // pair coefficients, shifts and tail-corrections are temperature-independent here (checked
      // above), but the stored force-field temperature is kept consistent with the replica
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

  reweightingSamples.resize(numberOfReplicas);
  for (std::vector<Sample>& samples : reweightingSamples)
  {
    samples.reserve(numberOfProductionCycles / sampleReweightingEvery + 1uz);
  }

  stepsPerReplica.assign(numberOfReplicas, 0uz);
  swapAttemptsPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAcceptedPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAttemptsPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
  swapAcceptedPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
}

void ReweightedHistogram::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ReweightedHistogram::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  stream.open("output/output.reweighted_histogram.txt", std::ios::out);
  outputJsonFileName = "output/output.reweighted_histogram.json";

  const System& front = systems.front();
  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  std::print(stream, "Reweighted histogram\n");
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
    std::print(stream, "Configuration swaps:                         disabled\n");
  }
  else
  {
    std::print(stream, "Configuration-swap sweep every:              {} cycles\n", parallelTemperingSwapEvery);
    std::print(stream, "  (sweeps alternate between the temperature and the pressure direction of the grid)\n");
  }
  std::print(stream, "Reweighting (N, U) sample every:             {} cycles\n", sampleReweightingEvery);
  std::print(stream, "Reweighted-isotherm temperatures:           ");
  for (double T : reweightingTemperatures)
  {
    std::print(stream, " {}", T);
  }
  std::print(stream, " [K]\n");
  std::print(stream, "Reweighted-isotherm pressure range:          {:.5e} - {:.5e} [Pa], {} log-spaced points\n\n",
             reweightingPressureRange.first, reweightingPressureRange.second, reweightingNumberOfPressures);

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
  outputJson["initialization"]["sampleReweightingEvery"] = sampleReweightingEvery;
  outputJson["initialization"]["reweightingTemperatures"] = reweightingTemperatures;
  outputJson["initialization"]["reweightingPressureRange"] =
      std::vector<double>{reweightingPressureRange.first, reweightingPressureRange.second};
  outputJson["initialization"]["reweightingNumberOfPressures"] = reweightingNumberOfPressures;

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-replica output files: each worker thread writes exclusively to its own stream
  replicaStreams.reserve(numberOfReplicas);
  replicaJsonFileNames.reserve(numberOfReplicas);
  replicaJsons.resize(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(std::format("output/output_{}_{}.reweighted_histogram.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                std::ios::out);
    replicaJsonFileNames.emplace_back(std::format("output/output_{}_{}.reweighted_histogram.r{}.json",
                                                  system.temperature, system.input_pressure, replicaId));

    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "Reweighted histogram: replica {} of {} (temperature {} [K], pressure {} [Pa])\n",
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

void ReweightedHistogram::performReplicaCycle(std::size_t replicaId, SimulationStage stage, std::size_t currentBlock)
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

void ReweightedHistogram::performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept
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
      std::print(stream, "Reweighted-histogram sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%)\n",
                 swapSweeps, stageName, cycle, numberOfCycles, swapAccepted, swapAttempts,
                 100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));
      std::flush(stream);
    }

    // running snapshot of the directly-measured adsorption isotherms: all worker threads are
    // parked on the barrier, so reading the per-replica loading averages is race-free
    if (stage == SimulationStage::Production)
    {
      writeIsothermSnapshot();
    }
  }
}

void ReweightedHistogram::runStage(SimulationStage stage, std::size_t numberOfCycles)
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

                // the raw (N, U) reweighting sample of this replica: the molecule count of the
                // single adsorbate component and the potential energy, tagged with the block
                if (cycle % sampleReweightingEvery == 0uz)
                {
                  reweightingSamples[replicaId].push_back(
                      Sample{.energy = system.runningEnergies.potentialEnergy(),
                             .numberOfMolecules =
                                 static_cast<std::uint32_t>(system.numberOfIntegerMoleculesPerComponent[0]),
                             .block = static_cast<std::uint32_t>(estimation.currentBin)});
                }

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

void ReweightedHistogram::output()
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

  // final directly-measured adsorption isotherms (also covers runs with the swaps disabled)
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

  std::print(stream, "Replica-exchange swap statistics\n");
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

  // the multiple-histogram reweighting analysis: combines the pooled (N, U) samples of all
  // replicas and writes the reweighted isotherms and free energies
  performReweightingAnalysis();

  std::print(stream, "Production run CPU timings of the MC moves summed over replicas and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n", totalPreInitializationSimulationTime.count());
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Reweighting analysis time:      {:14f} [s]\n", totalReweightingAnalysisTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n",
             (totalSimulationTime + totalReweightingAnalysisTime).count());
  std::print(stream, "\n\n");
  std::flush(stream);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["replicaExchange"]["sweeps"] = swapSweeps;
  outputJson["output"]["replicaExchange"]["attempts"] = swapAttempts;
  outputJson["output"]["replicaExchange"]["accepted"] = swapAccepted;
  outputJson["output"]["replicaExchange"]["attemptsPerTemperaturePair"] = swapAttemptsPerTemperaturePair;
  outputJson["output"]["replicaExchange"]["acceptedPerTemperaturePair"] = swapAcceptedPerTemperaturePair;
  outputJson["output"]["replicaExchange"]["attemptsPerPressurePair"] = swapAttemptsPerPressurePair;
  outputJson["output"]["replicaExchange"]["acceptedPerPressurePair"] = swapAcceptedPerPressurePair;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["reweightingAnalysis"] = totalReweightingAnalysisTime.count();
  outputJson["output"]["cpuTimings"]["total"] = (totalSimulationTime + totalReweightingAnalysisTime).count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ReweightedHistogram::writeIsothermSnapshot() const
{
  // one isotherm file per temperature: every replica along the pressure ladder contributes one
  // directly-measured point, assembled from the per-replica block averages
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    std::ofstream isotherm(std::format("output/isotherm_{}.reweighted_histogram.txt",
                                       temperatures[temperatureIndex]),
                           std::ios::trunc);

    std::print(isotherm, "# Reweighted histogram: directly-measured adsorption isotherm at {} [K]\n",
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

void ReweightedHistogram::performReweightingAnalysis()
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  std::print(stream, "Multiple-histogram reweighting (WHAM) analysis\n");
  std::print(stream, "===============================================================================\n\n");

  std::size_t totalNumberOfSamples = 0uz;
  for (const std::vector<Sample>& samples : reweightingSamples)
  {
    totalNumberOfSamples += samples.size();
  }
  if (totalNumberOfSamples == 0uz)
  {
    std::print(stream, "    no (N, U) samples were collected (no production cycles); analysis skipped\n\n\n");
    return;
  }

  // per-state (replica) thermodynamic fields, all in internal units: the reduced potential of a
  // sample (N, U) at state i is u_i = beta_i U - N ln(beta_i f_i)
  std::vector<double> stateBeta(numberOfReplicas);
  std::vector<double> stateLogBetaFugacity(numberOfReplicas);
  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    const System& system = systems[stateId];
    const Component& component = system.components.front();
    const double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    stateBeta[stateId] = system.beta;
    stateLogBetaFugacity[stateId] = std::log(system.beta * fugacity);
  }

  // bin the pooled samples over (N, energy bin); the mean energy of the samples inside a bin is
  // used in the reduced potentials, which removes the leading-order binning error
  double energyMinimum = 1e300;
  double energyMaximum = -1e300;
  for (const std::vector<Sample>& samples : reweightingSamples)
  {
    for (const Sample& sample : samples)
    {
      energyMinimum = std::min(energyMinimum, sample.energy);
      energyMaximum = std::max(energyMaximum, sample.energy);
    }
  }
  constexpr std::size_t numberOfEnergyBins = 512uz;
  const double deltaEnergy =
      std::max((energyMaximum - energyMinimum) / static_cast<double>(numberOfEnergyBins), 1e-10);

  std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> binIds;
  std::vector<std::uint32_t> binMoleculeCount;      // N of the bin
  std::vector<double> binEnergySum;                 // sum of the sampled energies in the bin
  std::vector<std::vector<double>> binBlockCounts;  // per-bin sample count per block
  std::vector<std::vector<double>> stateBlockCounts(numberOfReplicas,
                                                    std::vector<double>(numberOfBlocks, 0.0));

  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    for (const Sample& sample : reweightingSamples[stateId])
    {
      const std::uint32_t energyBin = static_cast<std::uint32_t>(std::min(
          numberOfEnergyBins - 1uz,
          static_cast<std::size_t>(std::max(0.0, (sample.energy - energyMinimum) / deltaEnergy))));
      const auto [it, inserted] =
          binIds.try_emplace(std::make_pair(sample.numberOfMolecules, energyBin), binIds.size());
      if (inserted)
      {
        binMoleculeCount.push_back(sample.numberOfMolecules);
        binEnergySum.push_back(0.0);
        binBlockCounts.emplace_back(numberOfBlocks, 0.0);
      }
      const std::size_t binId = it->second;
      binEnergySum[binId] += sample.energy;
      binBlockCounts[binId][sample.block] += 1.0;
      stateBlockCounts[stateId][sample.block] += 1.0;
    }
  }
  const std::size_t numberOfBins = binIds.size();

  std::vector<double> binTotalCount(numberOfBins, 0.0);
  std::vector<double> binMeanEnergy(numberOfBins);
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    binTotalCount[binId] = std::accumulate(binBlockCounts[binId].begin(), binBlockCounts[binId].end(), 0.0);
    binMeanEnergy[binId] = binEnergySum[binId] / binTotalCount[binId];
  }

  // precomputed reduced potentials u_i(b) = beta_i U_b - N_b ln(beta_i f_i)
  std::vector<double> reducedPotential(numberOfBins * numberOfReplicas);
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      reducedPotential[binId * numberOfReplicas + stateId] =
          stateBeta[stateId] * binMeanEnergy[binId] -
          static_cast<double>(binMoleculeCount[binId]) * stateLogBetaFugacity[stateId];
    }
  }

  // the code base is compiled with -ffast-math, so infinities must not occur: empty bins/states
  // carry this finite 'log of zero' sentinel instead and drop out of the sums
  constexpr double logZero = -1e300;
  constexpr double logZeroThreshold = -1e299;

  // Solves the WHAM self-consistent equations in log space by direct iteration:
  //     ln Omega_b = ln M_b - ln sum_i n_i exp(g_i - u_i(b))
  //     g_i        = -ln sum_b Omega_b exp(-u_i(b))
  // (gauge-fixed to g_0 = 0). Returns the number of iterations performed and the final residual.
  auto solveWham = [&](const std::vector<double>& logStateCounts, const std::vector<double>& logBinCounts,
                       std::vector<double>& freeEnergies, std::vector<double>& logDensityOfStates,
                       std::size_t maximumNumberOfIterations) -> std::pair<std::size_t, double>
  {
    constexpr double tolerance = 1e-8;
    std::vector<double> updatedFreeEnergies(numberOfReplicas);
    logDensityOfStates.assign(numberOfBins, logZero);

    std::size_t iteration = 0uz;
    double residual = 1e30;
    for (; iteration < maximumNumberOfIterations && residual > tolerance; ++iteration)
    {
      for (std::size_t binId = 0; binId < numberOfBins; ++binId)
      {
        if (logBinCounts[binId] <= logZeroThreshold)
        {
          logDensityOfStates[binId] = logZero;
          continue;
        }
        double largest = logZero;
        for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
        {
          if (logStateCounts[stateId] <= logZeroThreshold) continue;
          largest = std::max(largest, logStateCounts[stateId] + freeEnergies[stateId] -
                                          reducedPotential[binId * numberOfReplicas + stateId]);
        }
        double sum = 0.0;
        for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
        {
          if (logStateCounts[stateId] <= logZeroThreshold) continue;
          sum += std::exp(logStateCounts[stateId] + freeEnergies[stateId] -
                          reducedPotential[binId * numberOfReplicas + stateId] - largest);
        }
        logDensityOfStates[binId] = logBinCounts[binId] - largest - std::log(sum);
      }

      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        double largest = logZero;
        for (std::size_t binId = 0; binId < numberOfBins; ++binId)
        {
          if (logDensityOfStates[binId] <= logZeroThreshold) continue;
          largest = std::max(largest,
                             logDensityOfStates[binId] - reducedPotential[binId * numberOfReplicas + stateId]);
        }
        if (largest <= logZeroThreshold)
        {
          updatedFreeEnergies[stateId] = 0.0;
          continue;
        }
        double sum = 0.0;
        for (std::size_t binId = 0; binId < numberOfBins; ++binId)
        {
          if (logDensityOfStates[binId] <= logZeroThreshold) continue;
          sum += std::exp(logDensityOfStates[binId] - reducedPotential[binId * numberOfReplicas + stateId] - largest);
        }
        updatedFreeEnergies[stateId] = -(largest + std::log(sum));
      }

      const double gauge = updatedFreeEnergies[0];
      residual = 0.0;
      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        updatedFreeEnergies[stateId] -= gauge;
        residual = std::max(residual, std::abs(updatedFreeEnergies[stateId] - freeEnergies[stateId]));
      }
      freeEnergies = updatedFreeEnergies;
    }
    return {iteration, residual};
  };

  // Reweights the density of states to an arbitrary state point (beta, ln(beta f)): the weight of
  // bin b is w_b ~ Omega_b exp(-beta U_b + N_b ln(beta f)). Returns the average molecule count
  // and the effective number of independent samples supporting the estimate (an overlap
  // diagnostic: small values flag extrapolation beyond the sampled (N, U) region).
  auto reweight = [&](const std::vector<double>& logDensityOfStates, double beta,
                      double logBetaFugacity) -> std::pair<double, double>
  {
    double largest = logZero;
    std::vector<double> logWeights(numberOfBins, logZero);
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logDensityOfStates[binId] <= logZeroThreshold) continue;
      logWeights[binId] = logDensityOfStates[binId] - beta * binMeanEnergy[binId] +
                          static_cast<double>(binMoleculeCount[binId]) * logBetaFugacity;
      largest = std::max(largest, logWeights[binId]);
    }
    if (largest <= logZeroThreshold)
    {
      return {0.0, 0.0};
    }
    double partitionSum = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      partitionSum += std::exp(logWeights[binId] - largest);
    }
    const double logPartition = largest + std::log(partitionSum);

    double averageMolecules = 0.0;
    double inverseEffectiveSamples = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      const double weight = std::exp(logWeights[binId] - logPartition);
      averageMolecules += weight * static_cast<double>(binMoleculeCount[binId]);
      inverseEffectiveSamples += weight * weight / binTotalCount[binId];
    }
    return {averageMolecules, 1.0 / std::max(inverseEffectiveSamples, 1e-300)};
  };

  // full-data solve
  std::vector<double> logStateCounts(numberOfReplicas);
  std::vector<double> logBinCounts(numberOfBins);
  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    const double count = std::accumulate(stateBlockCounts[stateId].begin(), stateBlockCounts[stateId].end(), 0.0);
    logStateCounts[stateId] = count > 0.0 ? std::log(count) : logZero;
  }
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    logBinCounts[binId] = std::log(binTotalCount[binId]);
  }

  std::vector<double> freeEnergies(numberOfReplicas, 0.0);
  std::vector<double> logDensityOfStates;
  const std::pair<std::size_t, double> convergence =
      solveWham(logStateCounts, logBinCounts, freeEnergies, logDensityOfStates, 100000uz);
  const std::size_t iterations = convergence.first;
  const double residual = convergence.second;

  std::print(stream, "    pooled samples:          {}\n", totalNumberOfSamples);
  std::print(stream, "    occupied (N, U) bins:    {} ({} energy bins, deltaU = {:.6e} [K])\n", numberOfBins,
             numberOfEnergyBins, Units::EnergyToKelvin * deltaEnergy);
  std::print(stream, "    WHAM iterations:         {} (residual {:.3e})\n\n", iterations, residual);
  if (residual > 1e-8)
  {
    std::print(stream, "    WARNING: the WHAM equations did not converge; the reweighted results are unreliable.\n");
    std::print(stream, "             Weak overlap between neighboring state points is the usual cause; use denser\n");
    std::print(stream, "             temperature/pressure ladders or longer production runs.\n\n");
  }

  // per-block solves (warm-started from the full solution) for the error bars
  std::vector<std::vector<double>> blockLogDensityOfStates(numberOfBlocks);
  std::vector<bool> blockIsValid(numberOfBlocks, false);
  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    std::vector<double> blockLogStateCounts(numberOfReplicas);
    std::vector<double> blockLogBinCounts(numberOfBins);
    double blockTotal = 0.0;
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      const double count = stateBlockCounts[stateId][block];
      blockTotal += count;
      blockLogStateCounts[stateId] = count > 0.0 ? std::log(count) : logZero;
    }
    if (blockTotal == 0.0) continue;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      const double count = binBlockCounts[binId][block];
      blockLogBinCounts[binId] = count > 0.0 ? std::log(count) : logZero;
    }

    std::vector<double> blockFreeEnergies = freeEnergies;
    solveWham(blockLogStateCounts, blockLogBinCounts, blockFreeEnergies, blockLogDensityOfStates[block], 100000uz);
    blockIsValid[block] = true;
  }

  // mean and confidence-interval error of the reweighted molecule count at (beta, ln(beta f));
  // the error is estimated from the spread of the per-block reweighting solutions
  auto reweightWithError = [&](double beta, double logBetaFugacity) -> std::tuple<double, double, double>
  {
    const auto [mean, effectiveSamples] = reweight(logDensityOfStates, beta, logBetaFugacity);
    std::vector<double> blockValues;
    blockValues.reserve(numberOfBlocks);
    for (std::size_t block = 0; block < numberOfBlocks; ++block)
    {
      if (!blockIsValid[block]) continue;
      blockValues.push_back(reweight(blockLogDensityOfStates[block], beta, logBetaFugacity).first);
    }
    const double error = blockErrorEstimate(blockValues, mean);
    return {mean, error, effectiveSamples};
  };

  // unit conversions shared by all replicas (identical framework)
  const System& front = systems.front();
  const Component& frontComponent = front.components.front();
  double toMoleculesPerUnitCell = 1.0;
  double toMolePerKg = 0.0;
  double toMgPerG = 0.0;
  if (front.framework.has_value())
  {
    const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
    toMoleculesPerUnitCell =
        1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
    const double frameworkMass = front.frameworkMass().value();
    toMolePerKg = 1000.0 / frameworkMass;
    toMgPerG = 1000.0 * frontComponent.totalMass / frameworkMass;
  }

  // free energies per state point, with the self-consistency check: the reweighted loading at
  // every simulated state point must agree with the directly-measured average loading
  {
    std::ofstream freeEnergyFile("output/reweighted_free_energies.reweighted_histogram.txt", std::ios::trunc);
    std::print(freeEnergyFile, "# Reweighted histogram: WHAM free energies per simulated state point\n");
    std::print(freeEnergyFile, "# g_i = -ln Xi_i (dimensionless, gauge g_0 = 0)\n");
    std::print(freeEnergyFile,
               "# the direct column is the block-averaged loading measured by the replica itself;\n");
    std::print(freeEnergyFile,
               "# the reweighted column is the WHAM estimate at the same state point (self-consistency check)\n");
    std::print(freeEnergyFile, "# column 1: replica id\n");
    std::print(freeEnergyFile, "# column 2: temperature [K]\n");
    std::print(freeEnergyFile, "# column 3: pressure [Pa]\n");
    std::print(freeEnergyFile, "# column 4: fugacity [Pa]\n");
    std::print(freeEnergyFile, "# column 5: free energy g_i [-]\n");
    std::print(freeEnergyFile, "# column 6: samples [-]\n");
    std::print(freeEnergyFile, "# column 7, 8: direct loading, error [molecules/cell]\n");
    std::print(freeEnergyFile, "# column 9, 10: reweighted loading, error [molecules/cell]\n");

    std::print(stream, "    self-consistency check (loading in molecules/cell)\n");
    std::print(stream, "    replica    temperature [K]    pressure [Pa]        direct              reweighted\n");
    std::print(stream, "    ---------------------------------------------------------------------------------\n");
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      const System& system = systems[stateId];
      const Component& component = system.components.front();
      const double fugacityPa =
          component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.input_pressure;

      const auto [loadingAverage, loadingError] = system.averageLoadings.average();
      const double directLoading = loadingAverage.numberOfMolecules[0];
      const double directError = loadingError.numberOfMolecules[0];

      const auto [reweightedLoading, reweightedError, effectiveSamples] =
          reweightWithError(stateBeta[stateId], stateLogBetaFugacity[stateId]);
      const double samples = std::accumulate(stateBlockCounts[stateId].begin(), stateBlockCounts[stateId].end(), 0.0);

      std::print(freeEnergyFile, "{:6d}   {:10.4f}   {: .6e}   {: .6e}   {: .10e}   {:10.0f}   "
                                 "{: .6e} {: .6e}   {: .6e} {: .6e}\n",
                 stateId, system.temperature, system.input_pressure, fugacityPa, freeEnergies[stateId], samples,
                 directLoading, directError, reweightedLoading, reweightedError);

      std::print(stream, "    {:7d}    {:15.4f}    {:13.5e}    {:9.4f} ± {:7.4f}   {:9.4f} ± {:7.4f}\n", stateId,
                 system.temperature, system.input_pressure, directLoading, directError, reweightedLoading,
                 reweightedError);

      outputJson["output"]["reweighting"]["states"][stateId] = {
          {"temperature", system.temperature},   {"pressure", system.input_pressure},
          {"fugacity", fugacityPa},              {"freeEnergy", freeEnergies[stateId]},
          {"samples", samples},                  {"directLoading", directLoading},
          {"directLoadingError", directError},   {"reweightedLoading", reweightedLoading},
          {"reweightedLoadingError", reweightedError}, {"effectiveSamples", effectiveSamples}};
    }
    std::print(stream, "\n");
  }

  // reweighted isotherms: one file per requested temperature, evaluated on a fine log-spaced
  // pressure grid (fugacity coefficients from the Peng-Robinson equation of state at every point)
  const std::vector<EquationOfState::FluidInput> fluidInputs = {
      {frontComponent.criticalTemperature, frontComponent.criticalPressure, frontComponent.acentricFactor, 1.0, true}};
  const double logPressureMinimum = std::log(reweightingPressureRange.first);
  const double logPressureMaximum = std::log(reweightingPressureRange.second);
  const std::size_t numberOfTargetPressures =
      reweightingPressureRange.second > reweightingPressureRange.first ? reweightingNumberOfPressures : 1uz;

  for (const double targetTemperature : reweightingTemperatures)
  {
    std::ofstream isotherm(std::format("output/reweighted_isotherm_{}.reweighted_histogram.txt", targetTemperature),
                           std::ios::trunc);
    std::print(isotherm, "# Reweighted histogram: reweighted adsorption isotherm at {} [K] ({})\n", targetTemperature,
               frontComponent.name);
    std::print(isotherm, "# errors from re-solving the WHAM equations per block\n");
    std::print(isotherm,
               "# the effective sample size (last column) diagnoses the overlap: small values flag\n"
               "# extrapolation beyond the sampled (N, U) region\n");
    std::print(isotherm, "# column 1: fugacity [Pa]\n");
    std::print(isotherm, "# column 2, 3: absolute loading, error [molecules/cell]\n");
    std::print(isotherm, "# column 4, 5: absolute loading, error [molecules/unit-cell]\n");
    std::print(isotherm, "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
    std::print(isotherm, "# column 8, 9: absolute loading, error [mg/g-framework]\n");
    std::print(isotherm, "# column 10: pressure [Pa]\n");
    std::print(isotherm, "# column 11: effective sample size [-]\n\n");

    const double targetBeta = 1.0 / (Units::KB * targetTemperature);
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfTargetPressures; ++pressureIndex)
    {
      const double targetPressure =
          numberOfTargetPressures == 1uz
              ? reweightingPressureRange.first
              : std::exp(logPressureMinimum + static_cast<double>(pressureIndex) *
                                                  (logPressureMaximum - logPressureMinimum) /
                                                  static_cast<double>(numberOfTargetPressures - 1uz));

      const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
          targetTemperature, targetPressure, fluidInputs, EquationOfState::Type::PengRobinson,
          EquationOfState::MixingRules::VanDerWaals);
      const double fugacityCoefficient = fluidResults.front().fugacityCoefficient.value_or(1.0);
      const double fugacityPa = fugacityCoefficient * targetPressure;
      const double logBetaFugacity = std::log(targetBeta * fugacityPa / Units::PressureConversionFactor);

      const auto [loading, loadingError, effectiveSamples] = reweightWithError(targetBeta, logBetaFugacity);

      std::print(isotherm,
                 "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}   "
                 "{: .6e}\n",
                 fugacityPa, loading, loadingError, toMoleculesPerUnitCell * loading,
                 toMoleculesPerUnitCell * loadingError, toMolePerKg * loading, toMolePerKg * loadingError,
                 toMgPerG * loading, toMgPerG * loadingError, targetPressure, effectiveSamples);
    }
  }

  // vapor-liquid equilibrium (bulk boxes only): at every requested subcritical temperature the
  // coexistence fugacity is located by the equal-weight criterion (Wilding): the reweighted
  // molecule-number distribution P(N) is bimodal, and the fugacity is bisected until the vapor
  // peak (N < N_cut) and the liquid peak (N >= N_cut) carry equal probability weight. The cut is
  // placed at the deepest valley of ln P(N) between the two peaks.
  if (!front.framework.has_value())
  {
    const double volume = front.simulationBox.volume;
    const std::uint32_t maximumMoleculeCount = *std::ranges::max_element(binMoleculeCount);

    std::vector<std::vector<std::size_t>> binsByMoleculeCount(maximumMoleculeCount + 1uz);
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      binsByMoleculeCount[binMoleculeCount[binId]].push_back(binId);
    }

    // ln P(N) (unnormalized) at (beta, ln(beta f)): the energy dimension is summed out
    auto moleculeDistribution = [&](const std::vector<double>& logDOS, double beta,
                                    double logBetaFugacity) -> std::vector<double>
    {
      std::vector<double> logProbability(maximumMoleculeCount + 1uz, logZero);
      for (std::size_t n = 0; n <= maximumMoleculeCount; ++n)
      {
        double largest = logZero;
        for (const std::size_t binId : binsByMoleculeCount[n])
        {
          if (logDOS[binId] <= logZeroThreshold) continue;
          largest = std::max(largest, logDOS[binId] - beta * binMeanEnergy[binId] +
                                          static_cast<double>(n) * logBetaFugacity);
        }
        if (largest <= logZeroThreshold) continue;
        double sum = 0.0;
        for (const std::size_t binId : binsByMoleculeCount[n])
        {
          if (logDOS[binId] <= logZeroThreshold) continue;
          sum += std::exp(logDOS[binId] - beta * binMeanEnergy[binId] +
                          static_cast<double>(n) * logBetaFugacity - largest);
        }
        logProbability[n] = largest + std::log(sum);
      }
      return logProbability;
    };

    // Places the phase cut at the deepest valley of ln P(N): the cut maximizing
    // min(peak left, peak right) - valley. Returns (cut, depth) or nothing when the distribution
    // is not bimodal (depth below threshold). Unsampled N are treated as deep-valley entries.
    constexpr double bimodalDepthThreshold = 0.2;
    auto findPhaseSplit = [&](const std::vector<double>& logProbability)
        -> std::optional<std::pair<std::size_t, double>>
    {
      const std::size_t size = logProbability.size();
      double lowestSampled = 1e300;
      for (const double value : logProbability)
      {
        if (value > logZeroThreshold) lowestSampled = std::min(lowestSampled, value);
      }
      std::vector<double> filled(size);
      for (std::size_t n = 0; n < size; ++n)
      {
        filled[n] = logProbability[n] > logZeroThreshold ? logProbability[n] : lowestSampled - 20.0;
      }

      std::vector<double> prefixMaximum(size), suffixMaximum(size);
      prefixMaximum[0] = filled[0];
      for (std::size_t n = 1; n < size; ++n) prefixMaximum[n] = std::max(prefixMaximum[n - 1], filled[n]);
      suffixMaximum[size - 1] = filled[size - 1];
      for (std::size_t n = size - 1; n > 0; --n) suffixMaximum[n - 1] = std::max(suffixMaximum[n], filled[n - 1]);

      std::size_t bestCut = 0uz;
      double bestDepth = -1e300;
      for (std::size_t cut = 1; cut + 1 < size; ++cut)
      {
        const double depth = std::min(prefixMaximum[cut - 1], suffixMaximum[cut + 1]) - filled[cut];
        if (depth > bestDepth)
        {
          bestDepth = depth;
          bestCut = cut;
        }
      }
      if (bestCut == 0uz || bestDepth < bimodalDepthThreshold)
      {
        return std::nullopt;
      }
      return std::make_pair(bestCut, bestDepth);
    };

    auto logSumExpRange = [&](const std::vector<double>& logProbability, std::size_t begin,
                              std::size_t end) -> double
    {
      double largest = logZero;
      for (std::size_t n = begin; n < end; ++n)
      {
        if (logProbability[n] > logZeroThreshold) largest = std::max(largest, logProbability[n]);
      }
      if (largest <= logZeroThreshold) return logZero;
      double sum = 0.0;
      for (std::size_t n = begin; n < end; ++n)
      {
        if (logProbability[n] > logZeroThreshold) sum += std::exp(logProbability[n] - largest);
      }
      return largest + std::log(sum);
    };

    struct CoexistencePoint
    {
      double logBetaFugacity;
      double fugacityPa;
      double saturationPressurePa;
      bool saturationPressureFromEmptyBox;  // exact normalization (false: ideal-gas reference)
      double vaporDensity;                  // [kg/m^3]
      double liquidDensity;                 // [kg/m^3]
      double vaporMolecules;
      double liquidMolecules;
      double valleyDepth;
      std::size_t cut;
      std::vector<double> logProbability;
    };

    // Scans the pressure range for a bracket where the liquid-vapor weight difference changes
    // sign (both endpoints bimodal), then bisects the fugacity to the equal-weight point.
    auto solveCoexistence = [&](const std::vector<double>& logDOS,
                                double temperature) -> std::optional<CoexistencePoint>
    {
      const double beta = 1.0 / (Units::KB * temperature);

      auto weightDifference = [&](double logBetaFugacity) -> std::optional<std::pair<double, std::size_t>>
      {
        const std::vector<double> logProbability = moleculeDistribution(logDOS, beta, logBetaFugacity);
        const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
        if (!split.has_value()) return std::nullopt;
        const double logVaporWeight = logSumExpRange(logProbability, 0uz, split->first);
        const double logLiquidWeight = logSumExpRange(logProbability, split->first, logProbability.size());
        return std::make_pair(logLiquidWeight - logVaporWeight, split->first);
      };

      auto logBetaFugacityOfPressure = [&](double pressure) -> double
      {
        const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
            temperature, pressure, fluidInputs, EquationOfState::Type::PengRobinson,
            EquationOfState::MixingRules::VanDerWaals);
        const double fugacityPa = fluidResults.front().fugacityCoefficient.value_or(1.0) * pressure;
        return std::log(beta * fugacityPa / Units::PressureConversionFactor);
      };

      // bracket search over the (log-spaced) scanned pressure range
      constexpr std::size_t numberOfScanPoints = 200uz;
      double lowerLogBetaFugacity = 0.0;
      double upperLogBetaFugacity = 0.0;
      double lowerDifference = 0.0;
      bool havePrevious = false;
      bool haveBracket = false;
      for (std::size_t scanIndex = 0; scanIndex < numberOfScanPoints; ++scanIndex)
      {
        const double pressure =
            std::exp(logPressureMinimum + static_cast<double>(scanIndex) *
                                              (logPressureMaximum - logPressureMinimum) /
                                              static_cast<double>(numberOfScanPoints - 1uz));
        const double logBetaFugacity = logBetaFugacityOfPressure(pressure);
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(logBetaFugacity);
        if (!difference.has_value())
        {
          havePrevious = false;
          continue;
        }
        if (havePrevious && lowerDifference < 0.0 && difference->first >= 0.0)
        {
          upperLogBetaFugacity = logBetaFugacity;
          haveBracket = true;
          break;
        }
        lowerLogBetaFugacity = logBetaFugacity;
        lowerDifference = difference->first;
        havePrevious = true;
      }
      if (!haveBracket) return std::nullopt;

      // bisection on ln(beta f) to the equal-weight point
      for (std::size_t iterationIndex = 0; iterationIndex < 100uz; ++iterationIndex)
      {
        const double midpoint = 0.5 * (lowerLogBetaFugacity + upperLogBetaFugacity);
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(midpoint);
        if (!difference.has_value()) break;  // lost bimodality mid-bracket; use the current interval
        if (difference->first < 0.0)
        {
          lowerLogBetaFugacity = midpoint;
        }
        else
        {
          upperLogBetaFugacity = midpoint;
        }
        if (std::abs(difference->first) < 1e-10) break;
      }

      const double logBetaFugacity = 0.5 * (lowerLogBetaFugacity + upperLogBetaFugacity);
      const std::vector<double> logProbability = moleculeDistribution(logDOS, beta, logBetaFugacity);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value()) return std::nullopt;
      const std::size_t cut = split->first;

      const double logVaporWeight = logSumExpRange(logProbability, 0uz, cut);
      const double logLiquidWeight = logSumExpRange(logProbability, cut, logProbability.size());
      double vaporMolecules = 0.0;
      double liquidMolecules = 0.0;
      for (std::size_t n = 0; n < logProbability.size(); ++n)
      {
        if (logProbability[n] <= logZeroThreshold) continue;
        if (n < cut)
        {
          vaporMolecules += static_cast<double>(n) * std::exp(logProbability[n] - logVaporWeight);
        }
        else
        {
          liquidMolecules += static_cast<double>(n) * std::exp(logProbability[n] - logLiquidWeight);
        }
      }

      // Saturation pressure from beta p V = ln Xi. The absolute normalization of Xi is fixed by
      // the empty-box state (N = 0, U = 0, weight exactly one) when it was sampled; otherwise it
      // is fixed approximately at the dilute end of the scanned pressure range, where the vapor
      // is nearly ideal and beta p V ~ beta f V.
      double saturationPressurePa;
      const double logPartitionAtCoexistence = logSumExpRange(logProbability, 0uz, logProbability.size());
      const bool saturationPressureFromEmptyBox = logProbability[0] > logZeroThreshold;
      if (saturationPressureFromEmptyBox)
      {
        saturationPressurePa =
            ((logPartitionAtCoexistence - logProbability[0]) / (beta * volume)) * Units::PressureConversionFactor;
      }
      else
      {
        const double referenceLogBetaFugacity = logBetaFugacityOfPressure(std::exp(logPressureMinimum));
        const std::vector<double> referenceLogProbability =
            moleculeDistribution(logDOS, beta, referenceLogBetaFugacity);
        const double logPartitionAtReference =
            logSumExpRange(referenceLogProbability, 0uz, referenceLogProbability.size());
        const double referenceBetaPressureVolume = std::exp(referenceLogBetaFugacity) * volume;
        saturationPressurePa =
            ((logPartitionAtCoexistence - logPartitionAtReference + referenceBetaPressureVolume) / (beta * volume)) *
            Units::PressureConversionFactor;
      }

      const double toKgPerCubicMeter = frontComponent.totalMass * Units::DensityConversionFactor / volume;
      return CoexistencePoint{
          .logBetaFugacity = logBetaFugacity,
          .fugacityPa = (std::exp(logBetaFugacity) / beta) * Units::PressureConversionFactor,
          .saturationPressurePa = saturationPressurePa,
          .saturationPressureFromEmptyBox = saturationPressureFromEmptyBox,
          .vaporDensity = vaporMolecules * toKgPerCubicMeter,
          .liquidDensity = liquidMolecules * toKgPerCubicMeter,
          .vaporMolecules = vaporMolecules,
          .liquidMolecules = liquidMolecules,
          .valleyDepth = split->second,
          .cut = cut,
          .logProbability = logProbability};
    };

    std::ofstream coexistence("output/vle_coexistence.reweighted_histogram.txt", std::ios::trunc);
    std::print(coexistence, "# Reweighted histogram: vapor-liquid coexistence of {} (equal-weight criterion)\n",
               frontComponent.name);
    std::print(coexistence, "# box volume {:.4f} [A^3]; errors from re-solving the WHAM equations per block\n",
               volume);
    std::print(coexistence,
               "# the saturation pressure follows from beta p V = ln Xi, normalized by the empty-box\n"
               "# state when sampled (column 13 = 0), otherwise approximately by an ideal-gas reference\n"
               "# at the dilute end of the scanned pressure range (column 13 = 1)\n");
    std::print(coexistence, "# column 1: temperature [K]\n");
    std::print(coexistence, "# column 2, 3: coexistence fugacity, error [Pa]\n");
    std::print(coexistence, "# column 4, 5: saturation pressure, error [Pa]\n");
    std::print(coexistence, "# column 6, 7: vapor density, error [kg/m^3]\n");
    std::print(coexistence, "# column 8, 9: liquid density, error [kg/m^3]\n");
    std::print(coexistence, "# column 10, 11: vapor, liquid peak [molecules]\n");
    std::print(coexistence, "# column 12: ln P(N) valley depth [-] (vanishes towards the critical point)\n");
    std::print(coexistence, "# column 13: saturation-pressure normalization (0 exact, 1 ideal-gas reference)\n\n");

    std::print(stream, "    vapor-liquid coexistence (equal-weight criterion)\n");
    std::print(stream, "    temperature [K]    fugacity [Pa]        P_sat [Pa]           "
                       "rho_vap [kg/m^3]     rho_liq [kg/m^3]\n");
    std::print(stream, "    -----------------------------------------------------------------"
                       "----------------------------------\n");

    for (const double targetTemperature : reweightingTemperatures)
    {
      const std::optional<CoexistencePoint> point = solveCoexistence(logDensityOfStates, targetTemperature);
      if (!point.has_value())
      {
        std::print(coexistence,
                   "# {} [K]: no coexistence found (supercritical, pressure range does not bracket the\n"
                   "#   transition, or the (N, U) sampling does not connect the two phases)\n",
                   targetTemperature);
        std::print(stream, "    {:15.4f}    no coexistence found in the scanned pressure range\n",
                   targetTemperature);
        continue;
      }

      // per-block coexistence solves for the error bars
      std::vector<double> blockFugacities, blockPressures, blockVaporDensities, blockLiquidDensities;
      for (std::size_t block = 0; block < numberOfBlocks; ++block)
      {
        if (!blockIsValid[block]) continue;
        const std::optional<CoexistencePoint> blockPoint =
            solveCoexistence(blockLogDensityOfStates[block], targetTemperature);
        if (!blockPoint.has_value()) continue;
        blockFugacities.push_back(blockPoint->fugacityPa);
        blockVaporDensities.push_back(blockPoint->vaporDensity);
        blockLiquidDensities.push_back(blockPoint->liquidDensity);
        blockPressures.push_back(blockPoint->saturationPressurePa);
      }
      const double fugacityError = blockErrorEstimate(blockFugacities, point->fugacityPa);
      const double pressureError = blockErrorEstimate(blockPressures, point->saturationPressurePa);
      const double vaporDensityError = blockErrorEstimate(blockVaporDensities, point->vaporDensity);
      const double liquidDensityError = blockErrorEstimate(blockLiquidDensities, point->liquidDensity);

      std::print(coexistence,
                 "{:10.4f}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   "
                 "{: .6e} {: .6e}   {: .6e}   {}\n",
                 targetTemperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
                 point->vaporDensity, vaporDensityError, point->liquidDensity, liquidDensityError,
                 point->vaporMolecules, point->liquidMolecules, point->valleyDepth,
                 point->saturationPressureFromEmptyBox ? 0 : 1);

      std::print(stream, "    {:15.4f}    {: .6e} ± {:.2e}   {: .6e} ± {:.2e}{}  {:9.4f} ± {:7.4f}   "
                         "{:9.4f} ± {:7.4f}\n",
                 targetTemperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
                 point->saturationPressureFromEmptyBox ? ' ' : '*', point->vaporDensity, vaporDensityError,
                 point->liquidDensity, liquidDensityError);

      // the coexistence molecule-number distribution (for inspection and finite-size scaling)
      std::ofstream distribution(std::format("output/vle_distribution_{}.reweighted_histogram.txt",
                                             targetTemperature),
                                 std::ios::trunc);
      std::print(distribution, "# Reweighted histogram: P(N) at coexistence, {} at {} [K]\n", frontComponent.name,
                 targetTemperature);
      std::print(distribution, "# coexistence fugacity {:.6e} [Pa], phase cut at N = {}\n", point->fugacityPa,
                 point->cut);
      std::print(distribution, "# column 1: N [molecules]\n");
      std::print(distribution, "# column 2: P(N) (normalized)\n");
      std::print(distribution, "# column 3: ln P(N)\n\n");
      const double logNormalization =
          logSumExpRange(point->logProbability, 0uz, point->logProbability.size());
      for (std::size_t n = 0; n < point->logProbability.size(); ++n)
      {
        if (point->logProbability[n] <= logZeroThreshold) continue;
        const double logNormalized = point->logProbability[n] - logNormalization;
        std::print(distribution, "{:6d}   {: .6e}   {: .6e}\n", n, std::exp(logNormalized), logNormalized);
      }

      nlohmann::json vleEntry;
      vleEntry["temperature"] = targetTemperature;
      vleEntry["fugacity"] = point->fugacityPa;
      vleEntry["fugacityError"] = fugacityError;
      vleEntry["saturationPressure"] = point->saturationPressurePa;
      vleEntry["saturationPressureError"] = pressureError;
      vleEntry["vaporDensity"] = point->vaporDensity;
      vleEntry["vaporDensityError"] = vaporDensityError;
      vleEntry["liquidDensity"] = point->liquidDensity;
      vleEntry["liquidDensityError"] = liquidDensityError;
      vleEntry["valleyDepth"] = point->valleyDepth;
      vleEntry["saturationPressureNormalization"] =
          point->saturationPressureFromEmptyBox ? "empty-box" : "ideal-gas reference";
      outputJson["output"]["reweighting"]["vaporLiquidEquilibrium"].push_back(vleEntry);
    }
    std::print(stream, "\n");
    std::print(stream, "    (* saturation pressure normalized approximately by an ideal-gas reference at the\n");
    std::print(stream, "       dilute end of the pressure range; the empty-box state was never sampled)\n");
    std::print(stream, "    note: near coexistence the direct and reweighted loadings in the self-consistency\n");
    std::print(stream, "          table may legitimately differ; a single replica stays in one (possibly\n");
    std::print(stream, "          metastable) phase while the reweighted average covers both phases\n\n");
    std::print(stream, "    coexistence written to output/vle_coexistence.reweighted_histogram.txt\n");
    std::print(stream, "    coexistence P(N) written to output/vle_distribution_{{T}}.reweighted_histogram.txt\n\n");
  }

  outputJson["output"]["reweighting"]["pooledSamples"] = totalNumberOfSamples;
  outputJson["output"]["reweighting"]["occupiedBins"] = numberOfBins;
  outputJson["output"]["reweighting"]["iterations"] = iterations;
  outputJson["output"]["reweighting"]["residual"] = residual;
  outputJson["output"]["reweighting"]["converged"] = residual <= 1e-8;
  outputJson["output"]["reweighting"]["freeEnergies"] = freeEnergies;

  std::print(stream, "    reweighted isotherms written to output/reweighted_isotherm_{{T}}.reweighted_histogram.txt\n");
  std::print(stream, "    free energies written to output/reweighted_free_energies.reweighted_histogram.txt\n\n\n");
  std::flush(stream);

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalReweightingAnalysisTime += (t2 - t1);
}

void ReweightedHistogram::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
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
