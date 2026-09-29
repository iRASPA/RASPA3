module;

module parallel_tempering_molecular_dynamics;

import std;

import stringutils;
import hardware_info;
import archive;
import graceful_shutdown;
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
import integrators;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import molecular_dynamics;
import thermobarostat;
import json;

// The analysis-property writers (RDFs, density grid, MSD, VACF, histograms, molecule properties)
// gate themselves on their own 'writeEvery'; a cycle argument of 0 forces the write (used for the
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
  if (system.propertyMSD.has_value())
  {
    system.propertyMSD->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyVACF.has_value())
  {
    system.propertyVACF->writeOutput(replicaId, system.components, cycle);
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
  if (system.propertyPolymerShape.has_value())
  {
    system.propertyPolymerShape->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyPolymerBackbone.has_value())
  {
    system.propertyPolymerBackbone->writeOutput(replicaId, system.components, cycle);
  }
}

// Kinetic and extended-system contributions of the current state, added to the potential-energy
// terms already in 'runningEnergies' (which the caller has just recomputed from the gradients).
static void refreshKineticAndExtendedEnergies(System& system)
{
  system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  system.runningEnergies.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                  system.framework, system.spanOfFrameworkGroupData());
  if (system.thermostat.has_value())
  {
    system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
  }
  if (system.thermobarostat.has_value())
  {
    system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);
  }
}

ParallelTemperingMolecularDynamics::ParallelTemperingMolecularDynamics(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
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

    // the heat bath of the replica: the thermostat masses are derived from the temperature when
    // the chain is initialized at the start of the equilibration stage
    if (system.thermostat.has_value())
    {
      system.thermostat->temperature = T;
    }
    if (system.thermobarostat.has_value())
    {
      system.thermobarostat->temperature = T;
    }

    // temperature-dependent potentials (Feynman-Hibbs) derive pair coefficients, shifts and
    // tail-corrections from the temperature
    if (system.forceField.temperature != T)
    {
      system.forceField.temperature = T;
      system.forceField.preComputeDerivedParameters();
      system.forceField.preComputePotentialShift();
      system.forceField.preComputeTailCorrection();
    }

    // the CBMC ideal-gas conformation reservoirs (used by the MC stages) are Boltzmann samples at
    // the system temperature
    system.buildConformationReservoirs();
    system.buildRecoilReferenceConformations();

    randoms.emplace_back(random.seed + replicaId + 1);
  }

  stepsPerReplica.assign(numberOfReplicas, 0uz);
  integratorsCPUTimePerReplica.assign(numberOfReplicas, IntegratorsCPUTime{});
  swapAttemptsPerPair.assign(numberOfReplicas - 1uz, 0uz);
  swapAcceptedPerPair.assign(numberOfReplicas - 1uz, 0uz);
  roundTrips.initialize(numberOfReplicas);
}

void ParallelTemperingMolecularDynamics::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ParallelTemperingMolecularDynamics::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    // the integrator needs gradients: a polynomial interpolation grid does not provide them
    if (system.forceField.interpolationScheme == ForceField::InterpolationScheme::Polynomial)
    {
      system.forceField.interpolationScheme = ForceField::InterpolationScheme::Tricubic;
    }
  }

  std::filesystem::create_directories("output");

  // on a binary-restart resume append to the existing output files (and skip re-printing the
  // headers) so each log continues where the interrupted run left off
  const bool resumedFromBinaryRestart = simulationStage != SimulationStage::Uninitialized;
  stream.open("output/output.parallel_tempering_md.txt", resumedFromBinaryRestart ? std::ios::app : std::ios::out);
  outputJsonFileName = "output/output.parallel_tempering_md.json";

  const System& front = systems.front();
  if (!resumedFromBinaryRestart)
  {
    std::print(stream, "{}", front.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());

    std::print(stream, "Replica-exchange molecular dynamics (parallel tempering)\n");
    std::print(stream, "===============================================================================\n\n");
    std::print(stream, "Number of temperatures / replicas / threads: {}\n", numberOfReplicas);
    std::print(stream, "Temperature ladder:                         ");
    for (double T : temperatures)
    {
      std::print(stream, " {}", T);
    }
    std::print(stream, " [K]\n");
    std::print(stream, "Ensemble:                                    {}\n",
               molecularDynamicsEnsembleName(front.molecularDynamicsEnsemble));
    std::print(stream, "Time step:                                   {} [ps]\n", front.timeStep);
    if (parallelTemperingSwapEvery == 0uz)
    {
      std::print(stream, "Configuration swaps:                         disabled\n\n");
    }
    else
    {
      std::print(stream, "Configuration-swap sweep every:              {} cycles (MD steps)\n\n",
                 parallelTemperingSwapEvery);
    }
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
  outputJson["initialization"]["ensemble"] = molecularDynamicsEnsembleName(front.molecularDynamicsEnsemble);
  outputJson["initialization"]["timeStep"] = front.timeStep;

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-replica output files: each worker thread writes exclusively to its own stream
  replicaStreams.reserve(numberOfReplicas);
  replicaJsonFileNames.reserve(numberOfReplicas);
  replicaJsons.resize(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(std::format("output/output_{}_{}.parallel_tempering_md.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                resumedFromBinaryRestart ? std::ios::app : std::ios::out);
    replicaJsonFileNames.emplace_back(std::format("output/output_{}_{}.parallel_tempering_md.r{}.json",
                                                  system.temperature, system.input_pressure, replicaId));

    if (!resumedFromBinaryRestart)
    {
      std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
      std::print(replicaStream, "{}", system.writeOutputHeader());
      std::print(replicaStream, "Replica-exchange molecular dynamics: replica {} of {} (temperature {} [K])\n",
                 replicaId, numberOfReplicas, system.temperature);
      std::print(replicaStream, "Random seed of this replica: {}\n\n", randoms[replicaId].seed);
      std::print(replicaStream, "{}\n", HardwareInfo::writeInfo());
      std::print(replicaStream, "{}", Units::printStatus());
      std::print(replicaStream, "{}", system.writeSystemStatus());
      std::print(replicaStream, "{}", system.forceField.printPseudoAtomStatus());
      std::print(replicaStream, "{}", system.forceField.printForceFieldStatus());
      std::print(replicaStream, "{}", system.writeComponentStatus());
      std::print(replicaStream, "{}", system.crossLinks.printStatus());
      std::print(replicaStream, "{}", system.writeNumberOfPseudoAtoms());
    }

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

void ParallelTemperingMolecularDynamics::performReplicaMonteCarloCycle(std::size_t replicaId, SimulationStage stage)
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
      case SimulationStage::PreInitialization:
        MC_Moves::performRandomMovePreInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      default:
        break;
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
  }
}

void ParallelTemperingMolecularDynamics::prepareReplicaMolecularDynamicsStage(std::size_t replicaId,
                                                                             SimulationStage stage)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];

  // the rigid-body state of the (semi-)rigid molecules is derived from the MC-generated positions
  system.initializeGroupData();
  system.initializeFrameworkGroupData();
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                        system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                        system.spanOfFrameworkGroupData());

  if (stage == SimulationStage::Equilibration)
  {
    // Maxwell-Boltzmann velocities at the replica temperature, drawn from the replica's own stream
    Integrators::initializeVelocities(rng, system.moleculeData, system.spanOfMoleculeAtoms(),
                                      system.spanOfMoleculeDynamics(), system.components, system.temperature,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(),
                                      &system.forceField, system.spanOfGroupData(),
                                      system.spanOfFrameworkGroupData());
    Integrators::removeCenterOfMassVelocityDrift(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());

    if (system.thermostat.has_value())
    {
      const bool flexibleFrameworkConstraint =
          system.framework && system.framework->hasMobileAtoms() &&
          system.numberOfFrameworkAtoms + system.spanOfMoleculeAtoms().size() > 1uz;
      if (flexibleFrameworkConstraint || (!system.framework.has_value() && system.numberOfMolecules() > 1uz))
      {
        system.translationalCenterOfMassConstraint = 3;
        system.thermostat->translationalCenterOfMassConstraint = 3;
      }
      system.thermostat->initialize(rng);
    }
  }

  system.precomputeTotalGradients();
  refreshKineticAndExtendedEnergies(system);
  system.conservedEnergy = system.runningEnergies.conservedEnergy();
  system.referenceEnergy = system.conservedEnergy;

  std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
  replicaStream << system.runningEnergies.printMD("Recomputed from scratch", system.referenceEnergy);
  std::print(replicaStream, "\n\n\n\n");
  std::flush(replicaStream);
}

void ParallelTemperingMolecularDynamics::performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept
{
  // in the MD stages the momenta travel with the configuration and are rescaled to the replica
  // temperature; in the MC stages there are no meaningful momenta yet (they are drawn at the start
  // of the equilibration stage)
  const bool molecularDynamicsStage =
      (stage == SimulationStage::Equilibration) || (stage == SimulationStage::Production);

  // the replicas keep their temperatures; the configurations migrate through the ladder.
  // alternate the pairing offset between sweeps so configurations can traverse the whole ladder
  const std::size_t offset = swapSweeps % 2uz;
  for (std::size_t replicaId = offset; replicaId + 1 < numberOfReplicas; replicaId += 2uz)
  {
    ++swapAttempts;
    ++swapAttemptsPerPair[replicaId];
    const bool accepted =
        molecularDynamicsStage
            ? MC_Moves::ParallelTemperingSwapMolecularDynamics(random, systems[replicaId], systems[replicaId + 1])
                  .has_value()
            : MC_Moves::ParallelTemperingSwap(random, systems[replicaId], systems[replicaId + 1]).has_value();
    if (accepted)
    {
      ++swapAccepted;
      ++swapAcceptedPerPair[replicaId];
      roundTrips.recordAcceptedSwap(replicaId, replicaId + 1);
    }
  }
  roundTrips.endOfSweep();
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
    std::print(stream, "Replica-exchange sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%), round trips {}\n",
               swapSweeps, stageName, cycle, numberOfCycles, swapAccepted, swapAttempts,
               100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)),
               roundTrips.roundTrips);
    std::flush(stream);
  }
}

void ParallelTemperingMolecularDynamics::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  // binary restart: stages the restart file was written after are skipped entirely; the stage it
  // was written in resumes from the checkpointed cycle (its per-stage preparation already ran
  // before the checkpoint and must not be repeated)
  if (stage < simulationStage)
  {
    return;
  }
  const std::size_t startCycle = (stage == simulationStage) ? cyclesCompletedThisStage : 0uz;

  simulationStage = stage;

  const bool molecularDynamicsStage =
      (stage == SimulationStage::Equilibration) || (stage == SimulationStage::Production);

  if (startCycle == 0uz)
  {
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
      }
    }
    if (stage == SimulationStage::Production)
    {
      for (System& system : systems)
      {
        system.mc_moves_statistics.clearMoveStatistics();
        system.mc_moves_cputime.clearTimingStatistics();
        system.accumulatedDrift = 0.0;

        for (Component& component : system.components)
        {
          component.mc_moves_statistics.clearMoveStatistics();
          component.mc_moves_cputime.clearTimingStatistics();

          component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                                 system.containsTheFractionalMolecule);
          component.lambdaGC.clear();
        }
      }
      std::fill(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);
      std::fill(integratorsCPUTimePerReplica.begin(), integratorsCPUTimePerReplica.end(), IntegratorsCPUTime{});
    }
    sweepsThisStage = 0uz;
  }

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::PreInitialization)  ? "Pre-initialization"
                                       : (stage == SimulationStage::Initialization)   ? "Initialization"
                                       : (stage == SimulationStage::Equilibration)    ? "Equilibration"
                                                                                      : "Production";
    const std::string_view method = molecularDynamicsStage ? "MD steps" : "MC cycles";
    if (startCycle == 0uz)
    {
      std::print(stream, "\n{} stage: {} {} on {} replicas/threads\n", stageName, numberOfCycles, method,
                 numberOfReplicas);
    }
    else
    {
      std::print(stream, "\n{} stage: resumed from binary restart at cycle {} of {} ({}) on {} replicas/threads\n",
                 stageName, startCycle, numberOfCycles, method, numberOfReplicas);
    }
    std::flush(stream);
  }

  // one worker thread per replica; the barrier completion performs the swap sweeps and writes the
  // periodic binary restart file (all worker threads are parked, so the state is consistent)
  auto onAllArrived = [this, stage, numberOfCycles]() noexcept
  {
    const std::size_t completedCycles = checkpointCycle.load(std::memory_order_relaxed);
    if (parallelTemperingSwapEvery != 0uz && completedCycles % parallelTemperingSwapEvery == 0uz)
    {
      performSwapSweep(stage, numberOfCycles);
    }
    const bool stopRequested = GracefulShutdown::requested();
    const bool binaryRestartDue =
        writeBinaryRestartEvery != 0uz &&
        (completedCycles % writeBinaryRestartEvery == 0uz || completedCycles == numberOfCycles);
    if (binaryRestartDue || stopRequested)
    {
      writeBinaryRestartFile(completedCycles);
    }
    // all worker threads are parked on this barrier, so the checkpoint just written is a
    // consistent snapshot: safe to exit here on a shutdown signal
    if (stopRequested)
    {
      // std::exit skips stack unwinding: flush the text output streams explicitly
      std::flush(stream);
      for (std::ofstream& replicaStream : replicaStreams) std::flush(replicaStream);
      GracefulShutdown::exitAfterCheckpoint();
    }
  };
  std::barrier synchronizationPoint(static_cast<std::ptrdiff_t>(numberOfReplicas), onAllArrived);

  const std::size_t stageCycleOffset = absoluteCycleOffset;

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, stageCycleOffset, startCycle, molecularDynamicsStage,
           &synchronizationPoint]()
          {
            System& system = systems[replicaId];

            // each thread prepares its own replica: total energies for the MC stages, the full
            // dynamical state (velocities, thermostat, gradients) for the MD stages
            if (molecularDynamicsStage)
            {
              if (startCycle == 0uz)
              {
                prepareReplicaMolecularDynamicsStage(replicaId, stage);
              }
              else
              {
                // resumed mid-stage: the checkpoint holds positions, momenta and thermostat state;
                // only the gradients (not serialized) have to be rebuilt
                system.precomputeTotalRigidEnergy();
                system.precomputeTotalGradients();
                Integrators::updateCenterOfMassAndQuaternionGradients(
                    system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                    system.components, system.spanOfGroupData(), system.framework, system.spanOfFrameworkDynamics(),
                    system.spanOfFrameworkGroupData());
                refreshKineticAndExtendedEnergies(system);
              }
            }
            else
            {
              system.precomputeTotalRigidEnergy();
              system.runningEnergies = system.computeTotalEnergies();
            }

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfProductionCycles));

            for (std::size_t cycle = startCycle; cycle != numberOfCycles; ++cycle)
            {
              if (stage == SimulationStage::Production)
              {
                estimation.setCurrentSample(cycle);

                // energy/pressure averages for the per-replica final report (sampled before the
                // step, as in the plain MD driver)
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

              if (molecularDynamicsStage)
              {
                // one time step per cycle
                system.runningEnergies = molecularDynamicsStep(system);
                system.conservedEnergy = system.runningEnergies.conservedEnergy();
                system.accumulatedDrift +=
                    std::abs((system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
                if (stage == SimulationStage::Production)
                {
                  ++stepsPerReplica[replicaId];
                }
              }
              else
              {
                performReplicaMonteCarloCycle(replicaId, stage);
              }

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
                // the gradients of the integrator step are current; reuse them for the force-based RDF
                if (system.forceBasedRDFSampleDue(cycle))
                {
                  system.sampleForceBasedRDFFromCurrentGradients(cycle, estimation.currentBin);
                }

                // analysis-property files (RDFs, density grid, MSD, VACF, histograms, molecule
                // properties); the writers gate on their own 'writeEvery'
                writeReplicaAnalysisOutputs(system, replicaId, cycle);
              }

              if (!molecularDynamicsStage && cycle % optimizeMCMovesEvery == 0uz)
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
                  case SimulationStage::PreInitialization:
                    std::print(replicaStream, "{}", system.writePreInitializationStatusReport(cycle, numberOfCycles));
                    std::print(replicaStream, "{}\n\n\n\n", system.runningEnergies.printMC(""));
                    break;
                  case SimulationStage::Initialization:
                    std::print(replicaStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    std::print(replicaStream, "{}\n\n\n\n", system.runningEnergies.printMC(""));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(replicaStream, "{}", system.writeEquilibrationStatusReportMD(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                    std::print(replicaStream, "{}", system.writeProductionStatusReportMD(cycle, numberOfCycles));
                    break;
                  default:
                    break;
                }
                std::flush(replicaStream);
              }

              // the only synchronization point between the threads: the configuration-swap sweep
              // and the periodic binary-restart checkpoint, both performed by the barrier completion
              const bool swapDue =
                  parallelTemperingSwapEvery != 0uz && (cycle + 1uz) % parallelTemperingSwapEvery == 0uz;
              const bool binaryRestartDue =
                  writeBinaryRestartEvery != 0uz &&
                  ((cycle + 1uz) % writeBinaryRestartEvery == 0uz || cycle + 1uz == numberOfCycles);
              if (swapDue || binaryRestartDue)
              {
                if (replicaId == 0uz)
                {
                  checkpointCycle.store(cycle + 1uz, std::memory_order_relaxed);
                }
                synchronizationPoint.arrive_and_wait();
              }
            }

            // the integrator timings are accumulated thread-locally; gather them per replica
            // (production only: the accumulators are reset at the start of that stage)
            if (stage == SimulationStage::Production)
            {
              integratorsCPUTimePerReplica[replicaId] += integratorsCPUTime;
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

void ParallelTemperingMolecularDynamics::output()
{
  std::size_t numberOfSteps = std::accumulate(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);

  MCMoveCpuTime total;
  IntegratorsCPUTime integratorsTotal;
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    total += systems[replicaId].mc_moves_cputime;
    integratorsTotal += integratorsCPUTimePerReplica[replicaId];
  }

  std::print(stream, "\n");
  std::print(stream, "===============================================================================\n");
  std::print(stream, "                             Simulation finished!\n");
  std::print(stream, "===============================================================================\n");
  std::print(stream, "\n");

  // potential-energy drift check of every replica (energies recomputed in parallel, one thread per
  // replica); the final state used by the per-replica reports is refreshed in the same pass
  std::vector<RunningEnergy> recomputed(systems.size());
  std::vector<double> potentialDrift(systems.size());
  {
    std::vector<std::jthread> threads;
    threads.reserve(systems.size());
    for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, &recomputed, &potentialDrift]()
          {
            System& system = systems[replicaId];

            // the running energies come from the integrator energy path (which, like the MD driver, omits the
            // tail corrections), so the drift has to be measured against that same path and not against the
            // Monte Carlo-style computeTotalEnergies()
            const RunningEnergy running = system.runningEnergies;
            system.precomputeTotalGradients();
            potentialDrift[replicaId] =
                Units::EnergyToKelvin * (running.potentialEnergy() - system.runningEnergies.potentialEnergy());
            system.runningEnergies = running;

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
  std::print(stream, "    (potential: running minus recomputed; conserved: accumulated relative drift of the\n");
  std::print(stream, "     conserved energy over the production steps, reset at every accepted swap)\n\n");
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    std::print(stream,
               "    replica {:4d} (temperature {:10.4f} [K]): potential {: .6e} [K]   conserved {: .6e}\n",
               replicaId, systems[replicaId].temperature, potentialDrift[replicaId],
               systems[replicaId].accumulatedDrift /
                   static_cast<double>(std::max(1uz, stepsPerReplica[replicaId])));
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run: {} MD steps summed over the replicas\n", numberOfSteps);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "Replica-exchange swap statistics\n");
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

  std::print(stream, "{}", roundTrips.writeStatistics(temperatures, parallelTemperingSwapEvery));

  std::print(stream, "Production run CPU timings summed over replicas\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics());
  std::print(stream, "{}", integratorsTotal.writeIntegratorsCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "\n");
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
  outputJson["output"]["parallelTempering"]["roundTrips"] = roundTrips.jsonStatistics();
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["total"] = totalSimulationTime.count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ParallelTemperingMolecularDynamics::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
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

    std::print(replicaStream, "{}",
               system.writeProductionStatusReportMD(numberOfProductionCycles, numberOfProductionCycles));

    const RunningEnergy drift = system.runningEnergies - recomputed[replicaId];
    std::print(replicaStream, "Potential energy drift (running minus recomputed): {: .6e} [K]\n",
               Units::EnergyToKelvin * drift.potentialEnergy());
    std::print(replicaStream, "Accumulated relative drift of the conserved energy: {: .6e} (per step, reset at swaps)\n",
               system.accumulatedDrift / static_cast<double>(std::max(1uz, stepsPerReplica[replicaId])));
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "Production run CPU timings of the MD simulation of this replica\n");
    std::print(replicaStream, "===============================================================================\n\n");
    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(replicaStream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(replicaStream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(replicaStream, "{}",
               integratorsCPUTimePerReplica[replicaId].writeIntegratorsCPUTimeStatistics(
                   totalProductionSimulationTime));
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "{}",
               system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework,
                                                              system.components));

    std::print(replicaStream, "Temperature averages and statistics:\n");
    std::print(replicaStream, "===============================================================================\n\n");
    std::print(replicaStream, "{}", system.averageTemperature.writeAveragesStatistics("Total"));
    std::print(replicaStream, "{}", system.averageTranslationalTemperature.writeAveragesStatistics("Translational"));
    std::print(replicaStream, "{}", system.averageRotationalTemperature.writeAveragesStatistics("Rotational"));

    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(replicaStream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    std::print(replicaStream, "{}",
               system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swappableComponents,
                                                                            system.components));
    std::print(replicaStream, "{}",
               system.averagePartialMolarProperties.writeAveragesStatistics(system.swappableComponents,
                                                                            system.components));
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
    replicaJsons[replicaId]["output"]["accumulatedConservedEnergyDrift"] = system.accumulatedDrift;
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

void ParallelTemperingMolecularDynamics::writeBinaryRestartFile(std::size_t cyclesCompleted) noexcept
{
  cyclesCompletedThisStage = cyclesCompleted;

  ::writeBinaryRestartFile(*this);
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ParallelTemperingMolecularDynamics& pt)
{
  archive << pt.versionNumber;

  archive << pt.random;

  archive << pt.numberOfProductionCycles;
  archive << pt.numberOfPreInitializationCycles;
  archive << pt.numberOfInitializationCycles;
  archive << pt.numberOfEquilibrationCycles;

  archive << pt.printEvery;
  archive << pt.optimizeMCMovesEvery;
  archive << pt.writeBinaryRestartEvery;

  archive << pt.numberOfBlocks;
  archive << pt.parallelTemperingSwapEvery;

  archive << pt.simulationStage;
  archive << pt.cyclesCompletedThisStage;

  archive << pt.temperatures;
  archive << pt.numberOfReplicas;
  archive << pt.systems;
  archive << pt.randoms;

  archive << pt.stepsPerReplica;
  archive << pt.integratorsCPUTimePerReplica;
  archive << pt.absoluteCycleOffset;

  archive << pt.swapSweeps;
  archive << pt.sweepsThisStage;
  archive << pt.swapAttempts;
  archive << pt.swapAccepted;
  archive << pt.swapAttemptsPerPair;
  archive << pt.swapAcceptedPerPair;
  archive << pt.roundTrips;

  archive << pt.totalPreInitializationSimulationTime;
  archive << pt.totalInitializationSimulationTime;
  archive << pt.totalEquilibrationSimulationTime;
  archive << pt.totalProductionSimulationTime;
  archive << pt.totalSimulationTime;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ParallelTemperingMolecularDynamics& pt)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > pt.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'ParallelTemperingMolecularDynamics' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> pt.random;

  archive >> pt.numberOfProductionCycles;
  archive >> pt.numberOfPreInitializationCycles;
  archive >> pt.numberOfInitializationCycles;
  archive >> pt.numberOfEquilibrationCycles;

  archive >> pt.printEvery;
  archive >> pt.optimizeMCMovesEvery;
  archive >> pt.writeBinaryRestartEvery;

  archive >> pt.numberOfBlocks;
  archive >> pt.parallelTemperingSwapEvery;

  archive >> pt.simulationStage;
  archive >> pt.cyclesCompletedThisStage;

  archive >> pt.temperatures;
  archive >> pt.numberOfReplicas;
  archive >> pt.systems;
  archive >> pt.randoms;

  archive >> pt.stepsPerReplica;
  archive >> pt.integratorsCPUTimePerReplica;
  archive >> pt.absoluteCycleOffset;

  archive >> pt.swapSweeps;
  archive >> pt.sweepsThisStage;
  archive >> pt.swapAttempts;
  archive >> pt.swapAccepted;
  archive >> pt.swapAttemptsPerPair;
  archive >> pt.swapAcceptedPerPair;
  archive >> pt.roundTrips;

  archive >> pt.totalPreInitializationSimulationTime;
  archive >> pt.totalInitializationSimulationTime;
  archive >> pt.totalEquilibrationSimulationTime;
  archive >> pt.totalProductionSimulationTime;
  archive >> pt.totalSimulationTime;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error("ParallelTemperingMolecularDynamics: error in binary restart\n");
  }

  return archive;
}
