module;

module parallel_tmmc;

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
import transition_matrix;
import mc_moves;
import mc_moves_cputime;
import mc_moves_statistics;
import cbmc;
import cbmc_chain_data;
import json;

// The analysis-property writers (RDFs, density grid, histograms, molecule properties) gate
// themselves on their own 'writeEvery'; a cycle argument of 0 forces the write (used for the
// final flush at the end of the run). The walker id keys the output filenames, so every walker
// writes its own set of files.
static void writeWalkerAnalysisOutputs(System& system, std::size_t walkerId, std::size_t cycle)
{
  if (system.propertyConventionalRadialDistributionFunction.has_value())
  {
    system.propertyConventionalRadialDistributionFunction->writeOutput(
        system.forceField, walkerId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyRadialDistributionFunction.has_value())
  {
    system.propertyRadialDistributionFunction->writeOutput(system.forceField, walkerId, system.simulationBox.volume,
                                                           system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyDensityGrid.has_value())
  {
    system.propertyDensityGrid->writeOutput(walkerId, system.simulationBox, system.forceField, system.framework,
                                            system.components, cycle);
  }
  if (system.averageEnergyHistogram.has_value())
  {
    system.averageEnergyHistogram->writeOutput(walkerId, cycle);
  }
  if (system.averageNumberOfMoleculesHistogram.has_value())
  {
    system.averageNumberOfMoleculesHistogram->writeOutput(walkerId, system.components, cycle);
  }
  if (system.propertyMoleculeProperties.has_value())
  {
    system.propertyMoleculeProperties->writeOutput(walkerId, system.components, cycle);
  }
}

ParallelTMMC::ParallelTMMC(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      numberOfBlocks(reader.numberOfBlocks),
      reweightingNumberOfPressures(reader.reweightingNumberOfPressures),
      temperatures(reader.parallelTemperingTemperatures),
      numberOfWindows(reader.tmmcNumberOfWindows)
{
  // a single 'ExternalTemperature' gives a one-temperature run
  if (temperatures.empty())
  {
    temperatures.push_back(reader.systems.front().temperature);
  }
  numberOfTemperatures = temperatures.size();
  numberOfWalkers = numberOfTemperatures * numberOfWindows;

  referencePressure = reader.systems.front().input_pressure;

  // the isotherm/coexistence scan range defaults to four decades around the reference pressure
  reweightingPressureRange =
      reader.reweightingPressureRange.value_or(std::make_pair(0.01 * referencePressure, 100.0 * referencePressure));

  // the macrostate windows: window w spans [windowBoundaries[w], windowBoundaries[w+1]];
  // neighboring windows share their endpoint macrostate, which stitches the collection matrices
  minMacrostate = reader.systems.front().tmmc.minMacrostate;
  maxMacrostate = reader.systems.front().tmmc.maxMacrostate;
  windowBoundaries.resize(numberOfWindows + 1uz);
  for (std::size_t windowIndex = 0; windowIndex <= numberOfWindows; ++windowIndex)
  {
    windowBoundaries[windowIndex] = minMacrostate + (windowIndex * (maxMacrostate - minMacrostate)) / numberOfWindows;
  }

  // the single declared system is replicated into one walker per (temperature, window) pair
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  systems.reserve(numberOfWalkers);
  for (std::size_t walkerId = 0; walkerId + 1 < numberOfWalkers; ++walkerId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  // walker (t, w) is pinned at temperature T_t and macrostate window w, with its own
  // random-number stream and its own transition-matrix statistics file
  randoms.reserve(numberOfWalkers);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      System& system = systems[walkerId];
      const double T = temperatures[temperatureIndex];

      system.temperature = T;
      system.beta = 1.0 / (Units::KB * T);

      if (system.forceField.temperature != T)
      {
        system.forceField.temperature = T;
        system.forceField.preComputeDerivedParameters();
        system.forceField.preComputePotentialShift();
        system.forceField.preComputeTailCorrection();
      }

      // convert the reference pressure to the per-temperature reference fugacity: the fugacity
      // coefficient is recomputed with the Peng-Robinson equation of state at (T_t, P_ref)
      // (an explicitly given 'FugacityCoefficient' would only be valid at one temperature)
      for (Component& component : system.components)
      {
        component.fugacityCoefficient = std::nullopt;
      }
      system.equationOfState =
          EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals, T,
                          referencePressure, system.simulationBox, system.heliumVoidFraction, system.components);

      // the CBMC ideal-gas conformation reservoirs are Boltzmann samples at the system temperature
      system.buildConformationReservoirs();

      // the walker is confined to its window; the bias is cleared once at the start of the
      // equilibration (statistics of the relaxing initial configurations are dropped)
      system.tmmc.minMacrostate = windowBoundaries[windowIndex];
      system.tmmc.maxMacrostate = windowBoundaries[windowIndex + 1uz];
      system.tmmc.rezeroAfterInitialization = true;
      system.tmmc.statisticsFileName = std::format("tmmc/tmmc_statistics_{}_w{}.parallel_tmmc.txt", T, windowIndex);
      system.tmmc.initialize();

      randoms.emplace_back(random.seed + walkerId + 1);
    }
  }

  blockCollectionMatrices.resize(numberOfWalkers);
  stepsPerWalker.assign(numberOfWalkers, 0uz);
}

void ParallelTMMC::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ParallelTMMC::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  stream.open("output/output.parallel_tmmc.txt", std::ios::out);
  outputJsonFileName = "output/output.parallel_tmmc.json";

  const System& front = systems.front();
  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  std::print(stream, "Parallel transition-matrix Monte Carlo (TMMC)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Number of temperatures:                      {}\n", numberOfTemperatures);
  std::print(stream, "Number of macrostate windows:                {}\n", numberOfWindows);
  std::print(stream, "Number of walkers / threads:                 {}\n", numberOfWalkers);
  std::print(stream, "Temperature ladder:                         ");
  for (double T : temperatures)
  {
    std::print(stream, " {}", T);
  }
  std::print(stream, " [K]\n");
  std::print(stream, "Reference pressure:                          {:.5e} [Pa]\n", referencePressure);
  std::print(stream, "Macrostate range:                            [{}, {}] molecules\n", minMacrostate, maxMacrostate);
  std::print(stream, "Bias update every:                           {} steps\n", front.tmmc.updateTMEvery);
  std::print(stream, "Isotherm/coexistence pressure scan:          {:.5e} - {:.5e} [Pa], {} log-spaced points\n\n",
             reweightingPressureRange.first, reweightingPressureRange.second, reweightingNumberOfPressures);

  std::print(stream, "Walker grid: walker (t, w) = t * {} + w\n", numberOfWindows);
  std::print(stream, "    walker    temperature [K]    window [molecules]\n");
  std::print(stream, "    ----------------------------------------------\n");
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      std::print(stream, "    {:6d}    {:15.4f}    [{}, {}]\n", walkerId, temperatures[temperatureIndex],
                 windowBoundaries[windowIndex], windowBoundaries[windowIndex + 1uz]);
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
  outputJson["initialization"]["referencePressure"] = referencePressure;
  outputJson["initialization"]["macrostateRange"] = std::vector<std::size_t>{minMacrostate, maxMacrostate};
  outputJson["initialization"]["windowBoundaries"] = windowBoundaries;
  outputJson["initialization"]["reweightingPressureRange"] =
      std::vector<double>{reweightingPressureRange.first, reweightingPressureRange.second};
  outputJson["initialization"]["reweightingNumberOfPressures"] = reweightingNumberOfPressures;

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-walker output files: each worker thread writes exclusively to its own stream
  walkerStreams.reserve(numberOfWalkers);
  walkerJsonFileNames.reserve(numberOfWalkers);
  walkerJsons.resize(numberOfWalkers);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      const System& system = systems[walkerId];
      walkerStreams.emplace_back(
          std::format("output/output_{}_w{}.parallel_tmmc.r{}.txt", system.temperature, windowIndex, walkerId),
          std::ios::out);
      walkerJsonFileNames.emplace_back(
          std::format("output/output_{}_w{}.parallel_tmmc.r{}.json", system.temperature, windowIndex, walkerId));

      std::ostream walkerStream(walkerStreams[walkerId].rdbuf());
      std::print(walkerStream, "{}", system.writeOutputHeader());
      std::print(walkerStream, "Parallel TMMC: walker {} of {} (temperature {} [K], window [{}, {}] molecules)\n",
                 walkerId, numberOfWalkers, system.temperature, system.tmmc.minMacrostate, system.tmmc.maxMacrostate);
      std::print(walkerStream, "Random seed of this walker: {}\n\n", randoms[walkerId].seed);
      std::print(walkerStream, "{}\n", HardwareInfo::writeInfo());
      std::print(walkerStream, "{}", Units::printStatus());
      std::print(walkerStream, "{}", system.writeSystemStatus());
      std::print(walkerStream, "{}", system.forceField.printPseudoAtomStatus());
      std::print(walkerStream, "{}", system.forceField.printForceFieldStatus());
      std::print(walkerStream, "{}", system.writeComponentStatus());
      std::print(walkerStream, "{}", system.writeNumberOfPseudoAtoms());

#ifdef VERSION
      walkerJsons[walkerId]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
      walkerJsons[walkerId]["seed"] = randoms[walkerId].seed;
      walkerJsons[walkerId]["walkerId"] = walkerId;
      walkerJsons[walkerId]["temperature"] = system.temperature;
      walkerJsons[walkerId]["window"] = std::vector<std::size_t>{system.tmmc.minMacrostate, system.tmmc.maxMacrostate};
      walkerJsons[walkerId]["initialization"]["initialConditions"] = system.jsonSystemStatus();
      walkerJsons[walkerId]["initialization"]["components"] = system.jsonComponentStatus();

      std::ofstream walkerJson(walkerJsonFileNames[walkerId]);
      walkerJson << walkerJsons[walkerId].dump(4);
    }
  }

  // interpolation grids are computed once and shared (copied) between the walkers
  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t walkerId = 1; walkerId < systems.size(); ++walkerId)
  {
    systems[walkerId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[walkerId].interpolationGrids = systems.front().interpolationGrids;
  }

  // every walker starts inside its window: molecules are grown with CBMC up to the lower window
  // boundary (the walk then explores the window under the flattening bias), in parallel
  {
    std::print(stream, "Growing the initial configurations into their windows (CBMC)\n");
    std::flush(stream);
    std::vector<std::jthread> threads;
    threads.reserve(numberOfWalkers);
    for (std::size_t walkerId = 0; walkerId < numberOfWalkers; ++walkerId)
    {
      threads.emplace_back(
          [this, walkerId]()
          {
            System& system = systems[walkerId];
            RandomNumber& rng = randoms[walkerId];
            const std::size_t componentId = 0uz;
            const std::size_t target = system.tmmc.minMacrostate;

            while (system.numberOfIntegerMoleculesPerComponent[componentId] < target)
            {
              std::optional<ChainGrowData> growData = std::nullopt;
              bool insideBlockedPocket{false};
              do
              {
                do
                {
                  growData = CBMC::growMoleculeSwapInsertion(
                      rng,
                      CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox,
                                        system.interpolationGrids, system.externalFieldInterpolationGrid,
                                        system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                        system.beta, system.forceField.cutOffFrameworkVDW,
                                        system.forceField.cutOffMoleculeVDW, system.forceField.cutOffCoulomb},
                      system.components[componentId], componentId, system.numberOfMolecules(), 1.0, false, false);
                } while (!growData || growData->energies.potentialEnergy() > system.forceField.energyOverlapCriteria);

                std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());
                insideBlockedPocket = system.insideBlockedPockets(system.components[componentId], newMolecule);
              } while (insideBlockedPocket);

              system.insertMolecule(componentId, growData->molecule, growData->atoms);
            }
          });
    }
    // jthreads join on scope exit
  }
  std::print(stream, "Initial configurations ready\n\n");
  std::flush(stream);
}

void ParallelTMMC::performWalkerCycle(std::size_t walkerId, SimulationStage stage, std::size_t currentBlock)
{
  System& system = systems[walkerId];
  RandomNumber& rng = randoms[walkerId];

  // every walker is self-contained; the Gibbs-style moves that need a partner system are not
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

        // Wang-Landau biasing of the CFCMC lambda moves (all state is owned by this walker)
        system.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, system.containsTheFractionalMolecule);
        system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem,
                                              currentBlock);
        ++stepsPerWalker[walkerId];
        break;
    }

    // the TMMC state sampling: the moves already recorded the unbiased acceptance probabilities
    // into the collection matrix; here the visit histogram is updated and the flattening bias is
    // re-derived from the collection matrix every 'TMMCUpdateEvery' steps
    system.tmmc.updateHistogram(system.numberOfIntegerMoleculesPerComponent[0]);
    system.tmmc.numberOfSteps++;
    if (stage == SimulationStage::Equilibration || stage == SimulationStage::Production)
    {
      system.tmmc.adjustBias();
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
    system.pairSwapLambdaSampleOccupancy();
    system.reactionLambdaSampleOccupancy();
  }
}

std::vector<double3> ParallelTMMC::productionCollectionMatrix(std::size_t walkerId) const
{
  std::vector<double3> matrix = systems[walkerId].tmmc.cmatrix;
  const std::vector<double3>& start = productionStartCollectionMatrices[walkerId];
  for (std::size_t index = 0; index < matrix.size(); ++index)
  {
    matrix[index] -= start[index];
  }
  return matrix;
}

void ParallelTMMC::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  simulationStage = stage;

  // serial per-stage preparation
  if (stage == SimulationStage::Equilibration)
  {
    for (System& system : systems)
    {
      // drop the collection-matrix statistics of the relaxing initial configurations
      system.tmmc.clearCMatrix();

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
    productionStartCollectionMatrices.resize(systems.size());
    productionStartHistograms.resize(systems.size());
    for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
    {
      System& system = systems[walkerId];

      // the collection matrix and visit histogram are NOT reset: the recorded acceptance
      // probabilities are unbiased regardless of the applied bias, so the equilibration
      // statistics remain valid and keep accumulating. Resetting them would make the periodic
      // bias re-derivation (which normalizes over the visited range only) discontinuous at the
      // edge of the recently revisited region, trapping the walker behind an artificial bias
      // wall. The production-start snapshots below give the production-only block increments
      // and coverage diagnostics.
      productionStartCollectionMatrices[walkerId] = system.tmmc.cmatrix;
      productionStartHistograms[walkerId] = system.tmmc.histogram;

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
    std::fill(stepsPerWalker.begin(), stepsPerWalker.end(), 0uz);
  }

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::PreInitialization) ? "Pre-initialization"
                                       : (stage == SimulationStage::Initialization)  ? "Initialization"
                                       : (stage == SimulationStage::Equilibration)   ? "Equilibration"
                                                                                     : "Production";
    std::print(stream, "\n{} stage: {} cycles on {} walkers/threads\n", stageName, numberOfCycles, numberOfWalkers);
    std::flush(stream);
  }

  const std::size_t stageCycleOffset = absoluteCycleOffset;

  // one worker thread per walker; the walkers are fully independent within a stage
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfWalkers);
    for (std::size_t walkerId = 0; walkerId < numberOfWalkers; ++walkerId)
    {
      threads.emplace_back(
          [this, walkerId, stage, numberOfCycles, stageCycleOffset]()
          {
            System& system = systems[walkerId];

            // each thread computes the total energies of its own walker
            if (stage == SimulationStage::PreInitialization || stage == SimulationStage::Initialization)
            {
              system.precomputeTotalRigidEnergy();
            }
            system.runningEnergies = system.computeTotalEnergies();

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfProductionCycles));
            std::size_t currentBlock = 0uz;

            for (std::size_t cycle = 0uz; cycle != numberOfCycles; ++cycle)
            {
              if (stage == SimulationStage::Production)
              {
                estimation.setCurrentSample(cycle);

                // cumulative production-only collection-matrix snapshot at every block boundary
                // (the per-block increments give the error bars of the analysis)
                if (estimation.currentBin != currentBlock)
                {
                  blockCollectionMatrices[walkerId].push_back(productionCollectionMatrix(walkerId));
                  currentBlock = estimation.currentBin;
                }
              }

              performWalkerCycle(walkerId, stage, estimation.currentBin);

              // time-evolution properties (number of molecules, volume): sampled over all stages,
              // indexed by the absolute cycle number; the writers gate on their own 'writeEvery'
              const std::size_t absoluteCycle = stageCycleOffset + cycle;
              system.samplePropertiesEvolution(absoluteCycle);
              if (system.propertyNumberOfMoleculesEvolution.has_value())
              {
                system.propertyNumberOfMoleculesEvolution->writeOutput(walkerId, absoluteCycle);
              }
              if (system.propertyVolumeEvolution.has_value())
              {
                system.propertyVolumeEvolution->writeOutput(walkerId, absoluteCycle);
              }

              if (stage == SimulationStage::Production)
              {
                system.sampleProperties(walkerId, estimation.currentBin, cycle);

                // analysis-property files (RDFs, density grid, histograms, molecule properties);
                // the writers gate on their own 'writeEvery'
                writeWalkerAnalysisOutputs(system, walkerId, cycle);

                // energy/pressure averages for the per-walker final report
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

              // Wang-Landau biasing-factor adjustment (all state is owned by this walker)
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
                // each thread writes exclusively to its own walker stream: no locking needed
                system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                              system.simulationBox);

                std::ostream walkerStream(walkerStreams[walkerId].rdbuf());
                switch (stage)
                {
                  case SimulationStage::PreInitialization:
                    std::print(walkerStream, "{}", system.writePreInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Initialization:
                    std::print(walkerStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(walkerStream, "{}", system.writeEquilibrationStatusReportMC(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                  {
                    std::string status_line = std::format("Current cycle: {} out of {}\n", cycle, numberOfCycles);
                    std::print(walkerStream, "{}", system.writeProductionStatusReportMC(status_line));
                    break;
                  }
                  default:
                    break;
                }
                std::flush(walkerStream);

                // one combined progress line, from the first walker
                if (walkerId == 0uz)
                {
                  std::scoped_lock lock(outputMutex);
                  std::print(stream, "Parallel-TMMC cycle {} of {} (walker 0)\n", cycle, numberOfCycles);
                  std::flush(stream);
                }
              }
            }

            // final snapshot: the cumulative production-only collection matrix at the end
            if (stage == SimulationStage::Production)
            {
              blockCollectionMatrices[walkerId].push_back(productionCollectionMatrix(walkerId));
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

void ParallelTMMC::output()
{
  std::size_t numberOfSteps = std::accumulate(stepsPerWalker.begin(), stepsPerWalker.end(), 0uz);

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

  // energy drift check of every walker (energies recomputed in parallel, one thread per walker);
  // the final state used by the per-walker reports is refreshed in the same pass
  std::vector<RunningEnergy> recomputed(systems.size());
  {
    std::vector<std::jthread> threads;
    threads.reserve(systems.size());
    for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
    {
      threads.emplace_back(
          [this, walkerId, &recomputed]()
          {
            System& system = systems[walkerId];
            recomputed[walkerId] = system.computeTotalEnergies();

            std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
            system.currentEnergyStatus = molecularPressure.first;
            system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
            system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                          system.simulationBox);

            // final per-walker transition-matrix statistics file
            system.tmmc.writeStatistics();
          });
    }
  }

  writeWalkerFinalReports(recomputed);

  std::print(stream, "Energy drift per walker\n");
  std::print(stream, "===============================================================================\n\n");
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    const RunningEnergy drift = systems[walkerId].runningEnergies - recomputed[walkerId];
    std::print(stream, "    walker {:4d} (temperature {:10.4f} [K], window [{}, {}]): drift {: .6e} [K]\n", walkerId,
               systems[walkerId].temperature, systems[walkerId].tmmc.minMacrostate,
               systems[walkerId].tmmc.maxMacrostate, Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over walkers and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  // macrostate coverage per walker: every state of the window must be visited for the stitched
  // ln Pi(N) to be reliable; the min/max visit counts diagnose the flatness of the biased walk
  std::print(stream, "Macrostate coverage per walker (production)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    walker    temperature [K]    window            visited    min visits    max visits\n");
  std::print(stream, "    -------------------------------------------------------------------------------\n");
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    const System& system = systems[walkerId];

    // production-only visits: the cumulative histogram minus its production-start snapshot
    std::vector<std::size_t> histogram = system.tmmc.histogram;
    if (walkerId < productionStartHistograms.size() &&
        productionStartHistograms[walkerId].size() == histogram.size())
    {
      for (std::size_t index = 0; index < histogram.size(); ++index)
      {
        histogram[index] -= productionStartHistograms[walkerId][index];
      }
    }
    const std::size_t visited =
        static_cast<std::size_t>(std::ranges::count_if(histogram, [](std::size_t count) { return count > 0uz; }));
    const std::size_t minVisits = histogram.empty() ? 0uz : *std::ranges::min_element(histogram);
    const std::size_t maxVisits = histogram.empty() ? 0uz : *std::ranges::max_element(histogram);
    std::print(stream, "    {:6d}    {:15.4f}    [{:5d}, {:5d}]    {:4d}/{:<4d}   {:10d}    {:10d}\n", walkerId,
               system.temperature, system.tmmc.minMacrostate, system.tmmc.maxMacrostate, visited, histogram.size(),
               minVisits, maxVisits);

    outputJson["output"]["tmmc"]["coverage"][walkerId] = {
        {"temperature", system.temperature},
        {"window", std::vector<std::size_t>{system.tmmc.minMacrostate, system.tmmc.maxMacrostate}},
        {"visited", visited},
        {"states", histogram.size()},
        {"minVisits", minVisits},
        {"maxVisits", maxVisits}};
  }
  std::print(stream, "\n\n");

  // the transition-matrix analysis: combines the collection matrices of the windows per
  // temperature and writes ln Pi(N), the reweighted isotherms and the vapor-liquid coexistence
  performTransitionMatrixAnalysis();

  std::print(stream, "Production run CPU timings of the MC moves summed over walkers and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n", totalPreInitializationSimulationTime.count());
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Analysis time:                  {:14f} [s]\n", totalAnalysisTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n", (totalSimulationTime + totalAnalysisTime).count());
  std::print(stream, "\n\n");
  std::flush(stream);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["analysis"] = totalAnalysisTime.count();
  outputJson["output"]["cpuTimings"]["total"] = (totalSimulationTime + totalAnalysisTime).count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ParallelTMMC::performTransitionMatrixAnalysis()
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  std::print(stream, "Transition-matrix analysis\n");
  std::print(stream, "===============================================================================\n\n");

  const std::size_t numberOfMacrostates = maxMacrostate - minMacrostate + 1uz;

  // the code base is compiled with -ffast-math, so infinities must not occur: unsampled states
  // carry this finite 'log of zero' sentinel instead and drop out of the sums
  constexpr double logZero = -1e300;
  constexpr double logZeroThreshold = -1e299;

  // ln Pi(N) from a collection matrix over the full macrostate range by the detailed-balance
  // recursion ln Pi(N+1) = ln Pi(N) + ln P(N -> N+1) - ln P(N+1 -> N); states outside the visited
  // range stay at the 'log of zero' sentinel
  auto computeLogPi = [&](const std::vector<double3>& collectionMatrix) -> std::vector<double>
  {
    std::vector<double> logPi(numberOfMacrostates, logZero);

    auto rowTotal = [&](std::size_t index) -> double
    { return collectionMatrix[index].x + collectionMatrix[index].y + collectionMatrix[index].z; };

    std::size_t first = numberOfMacrostates;
    std::size_t last = 0uz;
    for (std::size_t index = 0; index < numberOfMacrostates; ++index)
    {
      if (rowTotal(index) > 0.0)
      {
        first = std::min(first, index);
        last = std::max(last, index);
      }
    }
    if (first >= numberOfMacrostates) return logPi;

    logPi[first] = 0.0;
    for (std::size_t index = first; index < last; ++index)
    {
      double value = logPi[index];
      if (collectionMatrix[index].z > 0.0)
      {
        value += std::log(collectionMatrix[index].z) - std::log(rowTotal(index));
      }
      if (collectionMatrix[index + 1uz].x > 0.0)
      {
        value += std::log(rowTotal(index + 1uz)) - std::log(collectionMatrix[index + 1uz].x);
      }
      logPi[index + 1uz] = value;
    }
    return logPi;
  };

  auto logSumExpRange = [&](const std::vector<double>& values, std::size_t begin, std::size_t end) -> double
  {
    double largest = logZero;
    for (std::size_t n = begin; n < end; ++n)
    {
      if (values[n] > logZeroThreshold) largest = std::max(largest, values[n]);
    }
    if (largest <= logZeroThreshold) return logZero;
    double sum = 0.0;
    for (std::size_t n = begin; n < end; ++n)
    {
      if (values[n] > logZeroThreshold) sum += std::exp(values[n] - largest);
    }
    return largest + std::log(sum);
  };

  // ln Pi(N; f) = ln Pi(N; f_ref) + N ln(f / f_ref) (exact); N is the global molecule count
  auto reweightedDistribution = [&](const std::vector<double>& logPi, double deltaLogFugacity) -> std::vector<double>
  {
    std::vector<double> logProbability(numberOfMacrostates, logZero);
    for (std::size_t index = 0; index < numberOfMacrostates; ++index)
    {
      if (logPi[index] <= logZeroThreshold) continue;
      logProbability[index] = logPi[index] + static_cast<double>(minMacrostate + index) * deltaLogFugacity;
    }
    return logProbability;
  };

  // conditional average of N over the macrostate subrange [begin, end)
  auto conditionalAverageMolecules = [&](const std::vector<double>& logProbability, std::size_t begin,
                                         std::size_t end) -> double
  {
    const double logPartition = logSumExpRange(logProbability, begin, end);
    if (logPartition <= logZeroThreshold) return 0.0;
    double average = 0.0;
    for (std::size_t index = begin; index < end; ++index)
    {
      if (logProbability[index] <= logZeroThreshold) continue;
      average += static_cast<double>(minMacrostate + index) * std::exp(logProbability[index] - logPartition);
    }
    return average;
  };

  auto averageMolecules = [&](const std::vector<double>& logProbability) -> double
  { return conditionalAverageMolecules(logProbability, 0uz, numberOfMacrostates); };

  // Places the phase cut at the deepest valley of ln P(N): the cut maximizing
  // min(peak left, peak right) - valley. Returns (cut, depth) or nothing when the distribution
  // is not bimodal (depth below threshold). Unsampled N are treated as deep-valley entries.
  constexpr double bimodalDepthThreshold = 0.2;
  auto findPhaseSplit = [&](const std::vector<double>& logProbability) -> std::optional<std::pair<std::size_t, double>>
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

  // unit conversions shared by all walkers (identical framework/box)
  const System& front = systems.front();
  const Component& frontComponent = front.components.front();
  const double volume = front.simulationBox.volume;
  double toMoleculesPerUnitCell = 1.0;
  double toMolePerKg = 0.0;
  double toMgPerG = 0.0;
  if (front.framework.has_value())
  {
    const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
    toMoleculesPerUnitCell = 1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
    const double frameworkMass = front.frameworkMass().value();
    toMolePerKg = 1000.0 / frameworkMass;
    toMgPerG = 1000.0 * frontComponent.totalMass / frameworkMass;
  }

  const std::vector<EquationOfState::FluidInput> fluidInputs = {
      {frontComponent.criticalTemperature, frontComponent.criticalPressure, frontComponent.acentricFactor, 1.0, true}};
  const double logPressureMinimum = std::log(reweightingPressureRange.first);
  const double logPressureMaximum = std::log(reweightingPressureRange.second);

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

  std::ofstream coexistence;
  const bool doVLE = !front.framework.has_value();
  if (doVLE)
  {
    coexistence.open("output/vle_coexistence.parallel_tmmc.txt", std::ios::trunc);
    std::print(coexistence, "# Parallel TMMC: vapor-liquid coexistence of {} (equal-weight criterion)\n",
               frontComponent.name);
    std::print(coexistence, "# box volume {:.4f} [A^3]; errors from the per-block collection-matrix increments\n",
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
  }

  bool anyApproximateNormalization = false;

  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    const double temperature = temperatures[temperatureIndex];
    const double beta = 1.0 / (Units::KB * temperature);

    // the reference fugacity of this temperature (all windows share it)
    const System& windowFront = systems[walkerIndex(temperatureIndex, 0uz)];
    const Component& component = windowFront.components.front();
    const double referenceFugacity =
        component.molFraction * component.fugacityCoefficient.value_or(1.0) * windowFront.pressure;
    const double referenceLogFugacity = std::log(referenceFugacity);

    auto logBetaFugacityToDelta = [&](double fugacityInternal) -> double
    { return std::log(fugacityInternal) - referenceLogFugacity; };

    auto fugacityOfPressure = [&](double pressurePa) -> double
    {
      const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
          temperature, pressurePa, fluidInputs, EquationOfState::Type::PengRobinson,
          EquationOfState::MixingRules::VanDerWaals);
      return fluidResults.front().fugacityCoefficient.value_or(1.0) * pressurePa / Units::PressureConversionFactor;
    };

    // combine the collection matrices of the windows: the final (total) matrix and the per-block
    // increments; entries at the shared boundary macrostates simply add (the recorded acceptance
    // probabilities are unbiased estimates of the same transition probabilities)
    std::vector<double3> totalCollectionMatrix(numberOfMacrostates, double3(0.0, 0.0, 0.0));
    std::vector<std::vector<double3>> blockCollectionMatrix(
        numberOfBlocks, std::vector<double3>(numberOfMacrostates, double3(0.0, 0.0, 0.0)));

    std::size_t availableBlocks = numberOfBlocks;
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      const System& system = systems[walkerId];
      const std::vector<std::vector<double3>>& snapshots = blockCollectionMatrices[walkerId];
      availableBlocks = std::min(availableBlocks, snapshots.size());

      const std::size_t offset = system.tmmc.minMacrostate - minMacrostate;

      // the total uses the full cumulative collection matrix (equilibration plus production):
      // the recorded acceptance probabilities are unbiased, so all statistics count
      for (std::size_t index = 0; index < system.tmmc.cmatrix.size(); ++index)
      {
        totalCollectionMatrix[offset + index] += system.tmmc.cmatrix[index];
      }

      // the error bars come from the production-only per-block increments
      if (snapshots.empty()) continue;
      for (std::size_t block = 0; block < std::min(numberOfBlocks, snapshots.size()); ++block)
      {
        for (std::size_t index = 0; index < snapshots[block].size(); ++index)
        {
          double3 increment = snapshots[block][index];
          if (block > 0uz)
          {
            increment -= snapshots[block - 1uz][index];
          }
          blockCollectionMatrix[block][offset + index] += increment;
        }
      }
    }

    const std::vector<double> logPi = computeLogPi(totalCollectionMatrix);

    std::vector<std::vector<double>> blockLogPi;
    blockLogPi.reserve(availableBlocks);
    for (std::size_t block = 0; block < availableBlocks; ++block)
    {
      blockLogPi.push_back(computeLogPi(blockCollectionMatrix[block]));
    }

    // total visit histogram over the windows (a pure diagnostic)
    std::vector<std::size_t> totalHistogram(numberOfMacrostates, 0uz);
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const System& system = systems[walkerIndex(temperatureIndex, windowIndex)];
      const std::size_t offset = system.tmmc.minMacrostate - minMacrostate;
      for (std::size_t index = 0; index < system.tmmc.histogram.size(); ++index)
      {
        totalHistogram[offset + index] += system.tmmc.histogram[index];
      }
    }

    // ln Pi(N) at the reference fugacity, normalized to unit total probability, with the error
    // from the per-block solutions (normalized the same way, which removes the arbitrary gauge)
    {
      std::ofstream lnpiFile(std::format("output/lnpi_{}.parallel_tmmc.txt", temperature), std::ios::trunc);
      std::print(lnpiFile, "# Parallel TMMC: macrostate probability distribution of {} at {} [K]\n",
                 frontComponent.name, temperature);
      std::print(lnpiFile, "# at the reference fugacity {:.6e} [Pa] (reference pressure {:.6e} [Pa])\n",
                 referenceFugacity * Units::PressureConversionFactor, referencePressure);
      std::print(lnpiFile, "# reweight exactly with ln Pi(N; f) = ln Pi(N) + N ln(f / f_ref)\n");
      std::print(lnpiFile, "# column 1: N [molecules]\n");
      std::print(lnpiFile, "# column 2, 3: ln Pi(N), error (normalized to unit total probability)\n");
      std::print(lnpiFile, "# column 4: visits (summed over the windows)\n\n");

      const double logNormalization = logSumExpRange(logPi, 0uz, numberOfMacrostates);
      std::vector<double> blockLogNormalization(availableBlocks);
      for (std::size_t block = 0; block < availableBlocks; ++block)
      {
        blockLogNormalization[block] = logSumExpRange(blockLogPi[block], 0uz, numberOfMacrostates);
      }

      for (std::size_t index = 0; index < numberOfMacrostates; ++index)
      {
        if (logPi[index] <= logZeroThreshold) continue;
        const double normalized = logPi[index] - logNormalization;

        std::vector<double> blockValues;
        blockValues.reserve(availableBlocks);
        for (std::size_t block = 0; block < availableBlocks; ++block)
        {
          if (blockLogPi[block][index] > logZeroThreshold && blockLogNormalization[block] > logZeroThreshold)
          {
            blockValues.push_back(blockLogPi[block][index] - blockLogNormalization[block]);
          }
        }
        const double error = blockErrorEstimate(blockValues, normalized);

        std::print(lnpiFile, "{:6d}   {: .10e} {: .6e}   {}\n", minMacrostate + index, normalized, error,
                   totalHistogram[index]);
      }
    }

    // the equilibrium, adsorption-branch and desorption-branch loadings at a given fugacity:
    // when the reweighted Pi(N) is bimodal the adsorption branch is the conditional average over
    // the low-density basin (N below the deepest valley; the metastable states followed on the
    // way up), the desorption branch the conditional average over the high-density basin, and
    // the equilibrium loading averages over both. Where Pi(N) is unimodal all three coincide.
    // branch order: 0 = equilibrium, 1 = adsorption, 2 = desorption
    auto branchLoadings = [&](const std::vector<double>& logPiInput,
                              double delta) -> std::pair<std::array<double, 3>, bool>
    {
      const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
      const double equilibrium = averageMolecules(logProbability);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value())
      {
        return {{equilibrium, equilibrium, equilibrium}, false};
      }
      return {{equilibrium, conditionalAverageMolecules(logProbability, 0uz, split->first),
               conditionalAverageMolecules(logProbability, split->first, numberOfMacrostates)},
              true};
    };

    // reweighted isotherms <N>(f) on a log-spaced pressure grid (exact in the fugacity; the
    // fugacity coefficients from the Peng-Robinson equation of state at every point): the
    // equilibrium isotherm plus the adsorption and desorption branches (the hysteresis loop)
    {
      const std::array<std::string, 3> branchFileNames = {
          std::format("output/reweighted_isotherm_{}.parallel_tmmc.txt", temperature),
          std::format("output/adsorption_isotherm_{}.parallel_tmmc.txt", temperature),
          std::format("output/desorption_isotherm_{}.parallel_tmmc.txt", temperature)};
      const std::array<std::string_view, 3> branchDescriptions = {
          "equilibrium reweighted adsorption isotherm (averaged over both basins of Pi(N))",
          "adsorption branch (conditional average over the low-density basin of Pi(N);\n"
          "# the metastable states followed on adsorption - equals the equilibrium isotherm\n"
          "# where Pi(N) is unimodal)",
          "desorption branch (conditional average over the high-density basin of Pi(N);\n"
          "# the metastable states followed on desorption - equals the equilibrium isotherm\n"
          "# where Pi(N) is unimodal)"};

      std::array<std::ofstream, 3> branchFiles;
      for (std::size_t branch = 0; branch < 3uz; ++branch)
      {
        branchFiles[branch].open(branchFileNames[branch], std::ios::trunc);
        std::print(branchFiles[branch], "# Parallel TMMC: {} at {} [K] ({})\n", branchDescriptions[branch], temperature,
                   frontComponent.name);
        std::print(branchFiles[branch], "# errors from the per-block collection-matrix increments\n");
        std::print(branchFiles[branch],
                   "# the isotherm saturates artificially near the upper macrostate bound N = {}; points with\n"
                   "# <N> approaching the bound are truncated by the finite macrostate range\n",
                   maxMacrostate);
        std::print(branchFiles[branch], "# column 1: fugacity [Pa]\n");
        std::print(branchFiles[branch], "# column 2, 3: absolute loading, error [molecules/cell]\n");
        std::print(branchFiles[branch], "# column 4, 5: absolute loading, error [molecules/unit-cell]\n");
        std::print(branchFiles[branch], "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
        std::print(branchFiles[branch], "# column 8, 9: absolute loading, error [mg/g-framework]\n");
        std::print(branchFiles[branch], "# column 10: pressure [Pa]\n");
        std::print(branchFiles[branch], "# column 11: bimodal (1 when Pi(N) has two basins at this pressure)\n\n");
      }

      for (std::size_t pressureIndex = 0; pressureIndex < reweightingNumberOfPressures; ++pressureIndex)
      {
        const double targetPressure =
            reweightingNumberOfPressures == 1uz
                ? reweightingPressureRange.first
                : std::exp(logPressureMinimum + static_cast<double>(pressureIndex) *
                                                    (logPressureMaximum - logPressureMinimum) /
                                                    static_cast<double>(reweightingNumberOfPressures - 1uz));
        const double fugacityInternal = fugacityOfPressure(targetPressure);
        const double delta = logBetaFugacityToDelta(fugacityInternal);

        const auto [loadings, bimodal] = branchLoadings(logPi, delta);

        std::array<std::vector<double>, 3> blockValues;
        for (std::size_t branch = 0; branch < 3uz; ++branch)
        {
          blockValues[branch].reserve(availableBlocks);
        }
        for (std::size_t block = 0; block < availableBlocks; ++block)
        {
          const std::array<double, 3> blockLoadings = branchLoadings(blockLogPi[block], delta).first;
          for (std::size_t branch = 0; branch < 3uz; ++branch)
          {
            blockValues[branch].push_back(blockLoadings[branch]);
          }
        }

        for (std::size_t branch = 0; branch < 3uz; ++branch)
        {
          const double loading = loadings[branch];
          const double loadingError = blockErrorEstimate(blockValues[branch], loading);
          std::print(branchFiles[branch],
                     "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}   "
                     "{}\n",
                     fugacityInternal * Units::PressureConversionFactor, loading, loadingError,
                     toMoleculesPerUnitCell * loading, toMoleculesPerUnitCell * loadingError, toMolePerKg * loading,
                     toMolePerKg * loadingError, toMgPerG * loading, toMgPerG * loadingError, targetPressure,
                     bimodal ? 1 : 0);
        }
      }
    }

    outputJson["output"]["tmmc"]["temperatures"].push_back(temperature);

    // vapor-liquid coexistence (bulk boxes only) by the equal-weight criterion (Wilding): the
    // fugacity is bisected until the vapor peak (N < N_cut) and the liquid peak (N >= N_cut) of
    // the reweighted Pi(N) carry equal probability weight; the cut is placed at the deepest
    // valley of ln Pi(N) between the two peaks
    if (!doVLE) continue;

    auto solveCoexistence = [&](const std::vector<double>& logPiInput) -> std::optional<CoexistencePoint>
    {
      auto weightDifference = [&](double delta) -> std::optional<std::pair<double, std::size_t>>
      {
        const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
        const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
        if (!split.has_value()) return std::nullopt;
        const double logVaporWeight = logSumExpRange(logProbability, 0uz, split->first);
        const double logLiquidWeight = logSumExpRange(logProbability, split->first, logProbability.size());
        return std::make_pair(logLiquidWeight - logVaporWeight, split->first);
      };

      // bracket search over the (log-spaced) scanned pressure range
      constexpr std::size_t numberOfScanPoints = 200uz;
      double lowerDelta = 0.0;
      double upperDelta = 0.0;
      double lowerDifference = 0.0;
      bool havePrevious = false;
      bool haveBracket = false;
      for (std::size_t scanIndex = 0; scanIndex < numberOfScanPoints; ++scanIndex)
      {
        const double pressure =
            std::exp(logPressureMinimum + static_cast<double>(scanIndex) * (logPressureMaximum - logPressureMinimum) /
                                              static_cast<double>(numberOfScanPoints - 1uz));
        const double delta = logBetaFugacityToDelta(fugacityOfPressure(pressure));
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(delta);
        if (!difference.has_value())
        {
          havePrevious = false;
          continue;
        }
        if (havePrevious && lowerDifference < 0.0 && difference->first >= 0.0)
        {
          upperDelta = delta;
          haveBracket = true;
          break;
        }
        lowerDelta = delta;
        lowerDifference = difference->first;
        havePrevious = true;
      }
      if (!haveBracket) return std::nullopt;

      // bisection on ln f to the equal-weight point
      for (std::size_t iterationIndex = 0; iterationIndex < 100uz; ++iterationIndex)
      {
        const double midpoint = 0.5 * (lowerDelta + upperDelta);
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(midpoint);
        if (!difference.has_value()) break;  // lost bimodality mid-bracket; use the current interval
        if (difference->first < 0.0)
        {
          lowerDelta = midpoint;
        }
        else
        {
          upperDelta = midpoint;
        }
        if (std::abs(difference->first) < 1e-10) break;
      }

      const double delta = 0.5 * (lowerDelta + upperDelta);
      const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value()) return std::nullopt;
      const std::size_t cut = split->first;

      const double logVaporWeight = logSumExpRange(logProbability, 0uz, cut);
      const double logLiquidWeight = logSumExpRange(logProbability, cut, logProbability.size());
      double vaporMolecules = 0.0;
      double liquidMolecules = 0.0;
      for (std::size_t index = 0; index < logProbability.size(); ++index)
      {
        if (logProbability[index] <= logZeroThreshold) continue;
        const double n = static_cast<double>(minMacrostate + index);
        if (index < cut)
        {
          vaporMolecules += n * std::exp(logProbability[index] - logVaporWeight);
        }
        else
        {
          liquidMolecules += n * std::exp(logProbability[index] - logLiquidWeight);
        }
      }

      const double coexistenceFugacityInternal = std::exp(referenceLogFugacity + delta);

      // Saturation pressure from beta p V = ln Xi. The absolute normalization of Xi is fixed by
      // the empty-box state (N = 0, weight exactly one) when the macrostate range starts at zero
      // and the state was visited; otherwise it is fixed approximately at the dilute end of the
      // scanned pressure range, where the vapor is nearly ideal and beta p V ~ beta f V.
      double saturationPressurePa;
      const double logPartitionAtCoexistence = logSumExpRange(logProbability, 0uz, logProbability.size());
      const bool saturationPressureFromEmptyBox = (minMacrostate == 0uz) && (logProbability[0] > logZeroThreshold);
      if (saturationPressureFromEmptyBox)
      {
        saturationPressurePa =
            ((logPartitionAtCoexistence - logProbability[0]) / (beta * volume)) * Units::PressureConversionFactor;
      }
      else
      {
        const double referenceFugacityDilute = fugacityOfPressure(std::exp(logPressureMinimum));
        const std::vector<double> referenceLogProbability =
            reweightedDistribution(logPiInput, logBetaFugacityToDelta(referenceFugacityDilute));
        const double logPartitionAtReference =
            logSumExpRange(referenceLogProbability, 0uz, referenceLogProbability.size());
        const double referenceBetaPressureVolume = beta * referenceFugacityDilute * volume;
        saturationPressurePa =
            ((logPartitionAtCoexistence - logPartitionAtReference + referenceBetaPressureVolume) / (beta * volume)) *
            Units::PressureConversionFactor;
      }

      const double toKgPerCubicMeter = frontComponent.totalMass * Units::DensityConversionFactor / volume;
      return CoexistencePoint{.logBetaFugacity = std::log(beta * coexistenceFugacityInternal),
                              .fugacityPa = coexistenceFugacityInternal * Units::PressureConversionFactor,
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

    if (temperatureIndex == 0uz)
    {
      std::print(stream, "    vapor-liquid coexistence (equal-weight criterion)\n");
      std::print(stream,
                 "    temperature [K]    fugacity [Pa]        P_sat [Pa]           "
                 "rho_vap [kg/m^3]     rho_liq [kg/m^3]\n");
      std::print(stream,
                 "    -----------------------------------------------------------------"
                 "----------------------------------\n");
    }

    const std::optional<CoexistencePoint> point = solveCoexistence(logPi);
    if (!point.has_value())
    {
      std::print(coexistence,
                 "# {} [K]: no coexistence found (supercritical, pressure range does not bracket the\n"
                 "#   transition, or the macrostate range does not span both phases)\n",
                 temperature);
      std::print(stream, "    {:15.4f}    no coexistence found in the scanned pressure range\n", temperature);
      continue;
    }

    // per-block coexistence solves for the error bars
    std::vector<double> blockFugacities, blockPressures, blockVaporDensities, blockLiquidDensities;
    for (std::size_t block = 0; block < availableBlocks; ++block)
    {
      const std::optional<CoexistencePoint> blockPoint = solveCoexistence(blockLogPi[block]);
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

    anyApproximateNormalization = anyApproximateNormalization || !point->saturationPressureFromEmptyBox;

    std::print(coexistence,
               "{:10.4f}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   "
               "{: .6e} {: .6e}   {: .6e}   {}\n",
               temperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
               point->vaporDensity, vaporDensityError, point->liquidDensity, liquidDensityError, point->vaporMolecules,
               point->liquidMolecules, point->valleyDepth, point->saturationPressureFromEmptyBox ? 0 : 1);

    std::print(stream,
               "    {:15.4f}    {: .6e} ± {:.2e}   {: .6e} ± {:.2e}{}  {:9.4f} ± {:7.4f}   "
               "{:9.4f} ± {:7.4f}\n",
               temperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
               point->saturationPressureFromEmptyBox ? ' ' : '*', point->vaporDensity, vaporDensityError,
               point->liquidDensity, liquidDensityError);

    // the coexistence molecule-number distribution (for inspection and finite-size scaling)
    std::ofstream distribution(std::format("output/vle_distribution_{}.parallel_tmmc.txt", temperature),
                               std::ios::trunc);
    std::print(distribution, "# Parallel TMMC: P(N) at coexistence, {} at {} [K]\n", frontComponent.name, temperature);
    std::print(distribution, "# coexistence fugacity {:.6e} [Pa], phase cut at N = {}\n", point->fugacityPa,
               minMacrostate + point->cut);
    std::print(distribution, "# column 1: N [molecules]\n");
    std::print(distribution, "# column 2: P(N) (normalized)\n");
    std::print(distribution, "# column 3: ln P(N)\n\n");
    const double logNormalization = logSumExpRange(point->logProbability, 0uz, point->logProbability.size());
    for (std::size_t index = 0; index < point->logProbability.size(); ++index)
    {
      if (point->logProbability[index] <= logZeroThreshold) continue;
      const double logNormalized = point->logProbability[index] - logNormalization;
      std::print(distribution, "{:6d}   {: .6e}   {: .6e}\n", minMacrostate + index, std::exp(logNormalized),
                 logNormalized);
    }

    nlohmann::json vleEntry;
    vleEntry["temperature"] = temperature;
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
    outputJson["output"]["tmmc"]["vaporLiquidEquilibrium"].push_back(vleEntry);
  }

  if (doVLE)
  {
    std::print(stream, "\n");
    if (anyApproximateNormalization)
    {
      std::print(stream, "    (* saturation pressure normalized approximately by an ideal-gas reference at the\n");
      std::print(stream, "       dilute end of the pressure range; the empty-box state N = 0 was never sampled -\n");
      std::print(stream, "       set 'MacroStateMinimumNumberOfMolecules' to 0 for the exact normalization)\n");
    }
    std::print(stream, "    coexistence written to output/vle_coexistence.parallel_tmmc.txt\n");
    std::print(stream, "    coexistence P(N) written to output/vle_distribution_{{T}}.parallel_tmmc.txt\n\n");
  }

  std::print(stream, "    ln Pi(N) written to output/lnpi_{{T}}.parallel_tmmc.txt\n");
  std::print(stream, "    equilibrium isotherms written to output/reweighted_isotherm_{{T}}.parallel_tmmc.txt\n");
  std::print(stream, "    adsorption/desorption branches (hysteresis loop) written to\n");
  std::print(stream,
             "    output/adsorption_isotherm_{{T}}.parallel_tmmc.txt and "
             "output/desorption_isotherm_{{T}}.parallel_tmmc.txt\n\n\n");
  std::flush(stream);

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalAnalysisTime += (t2 - t1);
}

void ParallelTMMC::writeWalkerFinalReports(std::vector<RunningEnergy>& recomputed)
{
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    System& system = systems[walkerId];
    std::ostream walkerStream(walkerStreams[walkerId].rdbuf());

    std::print(walkerStream, "\n");
    std::print(walkerStream, "===============================================================================\n");
    std::print(walkerStream, "                             Simulation finished!\n");
    std::print(walkerStream, "===============================================================================\n");
    std::print(walkerStream, "\n");

    std::string status_line = std::format("Final state after {} cycles\n", numberOfProductionCycles);
    std::print(walkerStream, "{}", system.writeProductionStatusReportMC(status_line));

    const RunningEnergy drift = system.runningEnergies - recomputed[walkerId];
    walkerStream << system.runningEnergies.printMCDiff(recomputed[walkerId]);
    std::print(walkerStream, "\n\n");

    std::print(walkerStream, "Monte-Carlo moves statistics\n");
    std::print(walkerStream, "===============================================================================\n\n");
    std::print(walkerStream, "{}", system.writeMCMoveStatistics());

    std::print(walkerStream, "Production run CPU timings of the MC moves of this walker\n");
    std::print(walkerStream, "===============================================================================\n\n");
    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(walkerStream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(walkerStream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(walkerStream, "\n\n");

    // under the flattening transition-matrix bias the direct averages are flat-histogram
    // averages, not grand-canonical ones; they are reported as diagnostics only
    std::print(walkerStream, "NOTE: the walker samples under the flattening transition-matrix bias; the\n");
    std::print(walkerStream, "      averages below are biased flat-histogram averages (diagnostics only).\n");
    std::print(walkerStream, "      The physical results are in the combined analysis files.\n\n");
    std::print(
        walkerStream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));
    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(walkerStream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    std::print(walkerStream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    std::flush(walkerStream);

    // final flush of the analysis-property files (cycle 0 bypasses the 'writeEvery' gate)
    writeWalkerAnalysisOutputs(system, walkerId, 0uz);

    // final flush of the time-evolution files: the total cycle count is rounded up to a multiple
    // of the writer's own 'writeEvery' so its gate passes and all collected samples are written
    if (system.propertyNumberOfMoleculesEvolution.has_value() &&
        system.propertyNumberOfMoleculesEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyNumberOfMoleculesEvolution->writeEvery.value();
      system.propertyNumberOfMoleculesEvolution->writeOutput(
          walkerId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }
    if (system.propertyVolumeEvolution.has_value() && system.propertyVolumeEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyVolumeEvolution->writeEvery.value();
      system.propertyVolumeEvolution->writeOutput(walkerId,
                                                  ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }

    // per-walker json statistics
    walkerJsons[walkerId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    walkerJsons[walkerId]["output"]["recomputedEnergies"] = recomputed[walkerId].jsonMC();
    walkerJsons[walkerId]["output"]["drift"] = drift.jsonMC();
    walkerJsons[walkerId]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    walkerJsons[walkerId]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();
    for (const Component& component : system.components)
    {
      walkerJsons[walkerId]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }
    walkerJsons[walkerId]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    walkerJsons[walkerId]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();

    std::ofstream json(walkerJsonFileNames[walkerId]);
    json << walkerJsons[walkerId].dump(4);
  }
}
