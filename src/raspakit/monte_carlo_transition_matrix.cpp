module;

module monte_carlo_transition_matrix;

import std;

import stringutils;
import hardware_info;
import archive;
import system;
import randomnumbers;
import input_reader;
import component;
import averages;
import property_loading;
import units;
import property_enthalpy;
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
import mc_moves;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_statistics;
import property_pressure;
import transition_matrix;
import interactions_ewald;
import equation_of_states;
import interpolation_energy_grid;
import simulation_schedule;

namespace
{

bool isTMMCCrossSystemMove(Move::Types moveType)
{
  switch (moveType)
  {
    case Move::Types::GibbsVolume:
    case Move::Types::GibbsSwapCBMC:
    case Move::Types::GibbsSwapCFCMC:
    case Move::Types::GibbsSwapCBCFCMC:
    case Move::Types::GibbsConventionalCFCMC:
    case Move::Types::GibbsConventionalCBCFCMC:
      return true;
    default:
      return false;
  }
}

void sampleTMMCState(System& system, std::size_t componentId, bool adjustBias)
{
  const std::size_t N = system.numberOfIntegerMoleculesPerComponent[componentId];
  system.tmmc.updateHistogram(N);
  system.tmmc.numberOfSteps++;
  if (adjustBias)
  {
    system.tmmc.adjustBias();
  }
}

}  // namespace

MonteCarloTransitionMatrix::MonteCarloTransitionMatrix() : random(std::nullopt) {};

MonteCarloTransitionMatrix::MonteCarloTransitionMatrix(InputReader& reader) noexcept
    : numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      systems(std::move(reader.systems)),
      random(reader.randomSeed),
      outputJsons(systems.size()),
      estimation(reader.numberOfBlocks, reader.numberOfProductionCycles)
{
}

MonteCarloTransitionMatrix::MonteCarloTransitionMatrix(const SimulationSchedule& schedule, std::vector<System>& systems,
                                                       RandomNumber& randomSeed, std::size_t numberOfBlocks)
    : numberOfProductionCycles(schedule.numberOfProductionCycles),
      numberOfPreInitializationCycles(schedule.numberOfPreInitializationCycles),
      numberOfInitializationCycles(schedule.numberOfInitializationCycles),
      numberOfEquilibrationCycles(schedule.numberOfEquilibrationCycles),
      printEvery(schedule.printEvery),
      writeBinaryRestartEvery(schedule.writeBinaryRestartEvery),
      rescaleWangLandauEvery(schedule.rescaleWangLandauEvery),
      optimizeMCMovesEvery(schedule.optimizeMCMovesEvery),
      systems(systems),
      random(randomSeed),
      outputJsons(systems.size()),
      estimation(numberOfBlocks, schedule.numberOfProductionCycles)
{
}

System& MonteCarloTransitionMatrix::randomSystem()
{
  return systems[std::size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MonteCarloTransitionMatrix::run()
{
  switch (simulationStage)
  {
    case SimulationStage::PreInitialization:
      goto continuePreInitializationStage;
    case SimulationStage::Initialization:
      goto continueInitializationStage;
    case SimulationStage::Equilibration:
      goto continueEquilibrationStage;
    case SimulationStage::Production:
      goto continueProductionStage;
    default:
      break;
  }

continuePreInitializationStage:
  preInitialize();
continueInitializationStage:
  initialize();
continueEquilibrationStage:
  equilibrate();
continueProductionStage:
  production();

  output();
}

void MonteCarloTransitionMatrix::createOutputFiles()
{
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
}

void MonteCarloTransitionMatrix::performCycle()
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
    // move to 'slide' when implemented in llvm
    // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
    // std::ranges::views::slide(s, 2uz);

    std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
    System& selectedSystem = systems[selectedSystemPair.first];
    System& selectedSecondSystem = systems[selectedSystemPair.second];

    std::size_t selectedComponent = selectedSystem.randomComponent(random);
    Move::Types performedMove = Move::Types::Count;

    switch (simulationStage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::PreInitialization:
        performedMove = MC_Moves::performRandomMovePreInitialization(random, selectedSystem, selectedSecondSystem,
                                                                     selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Initialization:
        performedMove = MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectedSecondSystem,
                                                                  selectedComponent, fractionalMoleculeSystem);

        sampleTMMCState(selectedSystem, selectedComponent, false);
        if (isTMMCCrossSystemMove(performedMove) && &selectedSecondSystem != &selectedSystem)
        {
          sampleTMMCState(selectedSecondSystem, selectedComponent, false);
        }
        break;
      case SimulationStage::Equilibration:
        performedMove = MC_Moves::performRandomMoveEquilibration(random, selectedSystem, selectedSecondSystem,
                                                                 selectedComponent, fractionalMoleculeSystem);

        sampleTMMCState(selectedSystem, selectedComponent, true);
        if (isTMMCCrossSystemMove(performedMove) && &selectedSecondSystem != &selectedSystem)
        {
          sampleTMMCState(selectedSecondSystem, selectedComponent, true);
        }

        selectedSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, selectedSystem.containsTheFractionalMolecule);
        selectedSecondSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
            selectedSecondSystem.containsTheFractionalMolecule);

        selectedSystem.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        selectedSecondSystem.pairSwapLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);

        selectedSystem.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        selectedSecondSystem.reactionLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        break;
      case SimulationStage::Production:
        performedMove =
            MC_Moves::performRandomMoveProduction(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                  fractionalMoleculeSystem, estimation.currentBin);

        sampleTMMCState(selectedSystem, selectedComponent, true);
        if (isTMMCCrossSystemMove(performedMove) && &selectedSecondSystem != &selectedSystem)
        {
          sampleTMMCState(selectedSecondSystem, selectedComponent, true);
        }

        numberOfSteps++;
        break;
    }

    selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSystem.containsTheFractionalMolecule);
    selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
        selectedSecondSystem.containsTheFractionalMolecule);
    selectedSystem.pairSwapLambdaSampleOccupancy();
    selectedSecondSystem.pairSwapLambdaSampleOccupancy();
    selectedSystem.reactionLambdaSampleOccupancy();
    selectedSecondSystem.reactionLambdaSampleOccupancy();
  }
}

void MonteCarloTransitionMatrix::preInitialize()
{
  std::chrono::steady_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::PreInitialization) goto continuePreInitializationStage;
  simulationStage = SimulationStage::PreInitialization;

  if (streams.empty())
  {
    createOutputFiles();

    for (std::size_t system_id{0}; System& system : systems)
    {
      // switch the fractional molecule on in the first system, and off in all others
      if (system_id == 0uz)
        system.containsTheFractionalMolecule = true;
      else
        system.containsTheFractionalMolecule = false;

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
      std::print(stream, "{}", system.reactions.printStatus());

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
      outputJsons[system_id]["initialization"]["reactions"] = system.reactions.jsonStatus();

      std::ofstream json(outputJsonFileNames[system_id]);
      json << outputJsons[system_id].dump(4);

      ++system_id;
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    std::ostream stream(streams[system_id].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    system.writeRestartFile(system_id);

    ++system_id;
  };

  for (currentCycle = 0uz; currentCycle != numberOfPreInitializationCycles; currentCycle++)
  {
    t1 = std::chrono::steady_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::ostream stream(streams[system_id].rdbuf());
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::print(stream, "{}",
                   system.writePreInitializationStatusReport(currentCycle, numberOfPreInitializationCycles));
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

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalPreInitializationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continuePreInitializationStage:;
  }
}

void MonteCarloTransitionMatrix::initialize()
{
  std::chrono::steady_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  if (streams.empty())
  {
    createOutputFiles();

    for (std::size_t system_id{0}; System& system : systems)
    {
      // switch the fractional molecule on in the first system, and off in all others
      if (system_id == 0uz)
        system.containsTheFractionalMolecule = true;
      else
        system.containsTheFractionalMolecule = false;

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
      std::print(stream, "{}", system.reactions.printStatus());

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
      outputJsons[system_id]["initialization"]["reactions"] = system.reactions.jsonStatus();

      std::ofstream json(outputJsonFileNames[system_id]);
      json << outputJsons[system_id].dump(4);

      ++system_id;
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.tmmc.initialize();

    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    std::ostream stream(streams[system_id].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    system.writeRestartFile(system_id);

    ++system_id;
  };

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
    t1 = std::chrono::steady_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::ostream stream(streams[system_id].rdbuf());
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

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

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalInitializationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continueInitializationStage:;
  }
}

void MonteCarloTransitionMatrix::equilibrate()
{
  std::chrono::steady_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();

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

    ++system_id;
  };

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::ostream stream(streams[system_id].rdbuf());
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

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

    if (currentCycle % rescaleWangLandauEvery == 0uz)
    {
      for (System& system : systems)
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
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalEquilibrationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continueEquilibrationStage:;
  }
}

void MonteCarloTransitionMatrix::production()
{
  double minBias{0uz};
  std::chrono::steady_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();

    system.tmmc.numberOfSteps = 0;

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

    ++system_id;
  }

  minBias = std::numeric_limits<double>::max();
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      double currentMinBias =
          *std::min_element(component.lambdaGC.biasFactor.cbegin(), component.lambdaGC.biasFactor.cend());
      minBias = currentMinBias < minBias ? currentMinBias : minBias;
    }

    const double pairSwapMinBias = system.pairSwapLambdaMinBias();
    if (pairSwapMinBias < minBias)
    {
      minBias = pairSwapMinBias;
    }

    if (system.usesReactionConventionalCFCMC())
    {
      const double reactionMinBias = system.reactionLambdaMinBias();
      if (reactionMinBias < minBias)
      {
        minBias = reactionMinBias;
      }
    }
  }
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      component.lambdaGC.normalize(minBias);
    }
    system.pairSwapLambdaNormalize(minBias);
    system.reactionLambdaNormalize(minBias);
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
      if (system.forceBasedRDFSampleDue(currentCycle))
      {
        system.sampleForceBasedRDFWithFullGradients(currentCycle, estimation.currentBin);
      }

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
        std::ostream stream(streams[system_id].rdbuf());
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::string status_line{std::format("Current cycle: {} out of {}\n", currentCycle, numberOfProductionCycles)};
        std::print(stream, "{}", system.writeProductionStatusReportMC(status_line));
        std::flush(stream);

        ++system_id;
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }

    t2 = std::chrono::steady_clock::now();

    totalProductionSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continueProductionStage:;
  }

  // Write the collection matrix
  for (System& system : systems)
  {
    system.tmmc.writeStatistics();
  }
}

void MonteCarloTransitionMatrix::output()
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
    std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n",
               totalPreInitializationSimulationTime.count());
    std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
    std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
    std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
    std::print(stream, "Total simulation time:          {:14f} [s]\n", totalSimulationTime.count());
    std::print(stream, "\n\n");

    // json statistics
    outputJsons[system_id]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    outputJsons[system_id]["output"]["recomputedEnergies"] = recomputedEnergies.jsonMC();
    outputJsons[system_id]["output"]["drift"] = drift.jsonMC();

    outputJsons[system_id]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    // outputJsons[system.systemId]["output"]["MCMoveStatistics"]["summedOverAllSystems"] =
    //    countTotal.jsonAllSystemStatistics(numberOfSteps);

    outputJsons[system_id]["output"]["cpuTimings"]["summedSystemsAndComponents"] =
        total.jsonOverallMCMoveCPUTimeStatistics(totalProductionSimulationTime);
    outputJsons[system_id]["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["total"] = totalSimulationTime.count();
    outputJsons[system_id]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();

    for (const Component& component : system.components)
    {
      // outputJsons[system.systemId]["output"]["MCMoveStatistics"][component.name]["percentage"] =
      //     component.mc_moves_statistics.jsonMCMoveStatistics(numberOfSteps);
      outputJsons[system_id]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }

    std::ofstream json(outputJsonFileNames[system_id]);
    json << outputJsons[system_id].dump(4);

    ++system_id;
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MonteCarloTransitionMatrix& mc)
{
  archive << mc.versionNumber;

  archive << mc.numberOfProductionCycles;
  archive << mc.numberOfPreInitializationCycles;
  archive << mc.numberOfInitializationCycles;
  archive << mc.numberOfEquilibrationCycles;
  archive << mc.numberOfSteps;

  archive << mc.printEvery;
  archive << mc.writeBinaryRestartEvery;
  archive << mc.rescaleWangLandauEvery;
  archive << mc.optimizeMCMovesEvery;

  archive << mc.currentCycle;
  archive << mc.simulationStage;

  archive << mc.systems;
  archive << mc.random;

  archive << mc.fractionalMoleculeSystem;

  archive << mc.estimation;

  archive << mc.totalPreInitializationSimulationTime;
  archive << mc.totalInitializationSimulationTime;
  archive << mc.totalEquilibrationSimulationTime;
  archive << mc.totalProductionSimulationTime;
  archive << mc.totalSimulationTime;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MonteCarloTransitionMatrix& mc)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MonteCarloTransitionMatrix' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> mc.numberOfProductionCycles;
  archive >> mc.numberOfPreInitializationCycles;
  archive >> mc.numberOfInitializationCycles;
  archive >> mc.numberOfEquilibrationCycles;
  archive >> mc.numberOfSteps;

  archive >> mc.printEvery;
  archive >> mc.writeBinaryRestartEvery;
  archive >> mc.rescaleWangLandauEvery;
  archive >> mc.optimizeMCMovesEvery;

  archive >> mc.currentCycle;
  archive >> mc.simulationStage;

  archive >> mc.systems;
  archive >> mc.random;

  archive >> mc.fractionalMoleculeSystem;

  archive >> mc.estimation;

  archive >> mc.totalPreInitializationSimulationTime;
  archive >> mc.totalInitializationSimulationTime;
  archive >> mc.totalEquilibrationSimulationTime;
  archive >> mc.totalProductionSimulationTime;
  archive >> mc.totalSimulationTime;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
  }
  std::cout << std::format("Magic number read correctly: {} vs {}\n", magicNumber,
                           static_cast<std::uint64_t>(0x6f6b6179));
  return archive;
}
