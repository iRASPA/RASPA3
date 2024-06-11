module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module monte_carlo;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <algorithm>;
import <numeric>;
import <ranges>;
import <chrono>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <span>;
import <string>;
import <optional>;
import <fstream>;
import <sstream>;
import <filesystem>;
import <tuple>;
import <ios>;
import <complex>;
import <exception>;
import <source_location>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#if defined(__has_include) && __has_include(<mdspan>)
import <mdspan>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

import stringutils;
import hardware_info;
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
import interactions_ewald;
import equation_of_states;
import hdf5;

MonteCarlo::MonteCarlo() : random(std::nullopt){};

MonteCarlo::MonteCarlo(InputReader& reader) noexcept
    : numberOfCycles(reader.numberOfCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      systems(std::move(reader.systems)),
      random(reader.randomSeed),
      estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
}

System& MonteCarlo::randomSystem() { return systems[size_t(random.uniform() * static_cast<double>(systems.size()))]; }

void MonteCarlo::run()
{
  switch (simulationStage)
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

continueInitializationStage:
  initialize();
continueEquilibrationStage:
  equilibrate();
continueProductionStage:
  production();

  output();
  log();
}

void MonteCarlo::createOutputFiles()
{
  std::filesystem::create_directories("output");
  for (System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.data", system.temperature, system.input_pressure, system.systemId);
    streams.emplace_back(fileNameString, std::ios::out);

    std::string fileNameStringHDF5 =
        std::format("output/output_{}_{}.s{}.h5", system.temperature, system.input_pressure, system.systemId);
    logs.emplace_back(HDF5Writer(fileNameStringHDF5));
  }
}

void MonteCarlo::initialize()
{
  std::chrono::system_clock::time_point t1, t2;
  size_t totalNumberOfMolecules{ 0uz };
  size_t totalNumberOfComponents{ 0uz };
  size_t numberOfStepsPerCycle{ 0uz };

  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  createOutputFiles();

  for (System& system : systems)
  {
    // switch the fractional molecule on in the first system, and off in all others
    if (system.systemId == 0uz)
      system.containsTheFractionalMolecule = true;
    else
      system.containsTheFractionalMolecule = false;
  }

  for (const System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());
    std::print(stream, "{}", system.writeSystemStatus());
    std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(stream, "{}", system.forceField.printForceFieldStatus());
    std::print(stream, "{}", system.writeComponentStatus());
    std::print(stream, "{}", system.reactions.printStatus());

    system.logMetaData(logs[system.systemId]);
    logs[system.systemId].writeMetaInfo("/", "seed", std::to_string(random.seed));
    HardwareInfo::logInfo(logs[system.systemId]);
    Units::logStatus(logs[system.systemId]);
    system.logSystemStatus(logs[system.systemId]);
    system.forceField.logPseudoAtomStatus(logs[system.systemId]);
    system.forceField.logForceFieldStatus(logs[system.systemId]);
    system.logComponentStatus(logs[system.systemId]);
  }

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    std::ostream stream(streams[system.systemId].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");
  };

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
    t1 = std::chrono::system_clock::now();
    totalNumberOfMolecules = std::transform_reduce(systems.begin(), systems.end(), 0uz,
        [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules();});
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

    for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
    {
      // move to 'slide' when implemented in llvm
      // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
      // std::ranges::views::slide(s, 2uz);

      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMove(random, selectedSystem, selectSecondSystem, selectedComponent,
                                  fractionalMoleculeSystem);

      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
        }
      }
    }

    //for (System& system : systems)
    //{
    //  system.checkCartesianPositions();
    //  system.checkMoleculeIds();
    //}

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::print(stream, "\n");
        std::flush(stream);
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
    t2 = std::chrono::system_clock::now();
    totalInitializationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continueInitializationStage:;
  }
}

void MonteCarlo::equilibrate()
{
  std::chrono::system_clock::time_point t1, t2;
  size_t totalNumberOfMolecules{ 0uz };
  size_t totalNumberOfComponents{ 0uz };
  size_t numberOfStepsPerCycle{ 0uz };

  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();

    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    t1 = std::chrono::system_clock::now();
    totalNumberOfMolecules = std::transform_reduce(systems.begin(), systems.end(), 0uz,
        [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules();});
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

    for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMove(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                  fractionalMoleculeSystem);

      selectedSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
          PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
          PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
          selectedSecondSystem.containsTheFractionalMolecule);

      selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
          selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
          selectedSecondSystem.containsTheFractionalMolecule);
    }

    //for (System& system : systems)
    //{
    //  system.checkCartesianPositions();
    //  system.checkMoleculeIds();
    //}

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::print(stream, "{}", system.writeEquilibrationStatusReportMC(currentCycle, numberOfEquilibrationCycles));
        std::flush(stream);
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
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % printEvery == 0uz)
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
    t2 = std::chrono::system_clock::now();
    totalEquilibrationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

    continueEquilibrationStage: ;

  }
}

void MonteCarlo::production()
{
  std::chrono::system_clock::time_point t1, t2;
  size_t totalNumberOfMolecules{0uz};
  size_t totalNumberOfComponents{0uz};
  size_t numberOfStepsPerCycle{0uz};
  double minBias{0.0};

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();

    system.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();
    system.mc_moves_count.clearCountStatistics();

    for (Component& component : system.components)
    {
      component.mc_moves_statistics.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();
      component.mc_moves_count.clearCountStatistics();

      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  minBias = std::numeric_limits<double>::max();
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      double currentMinBias =
          *std::min_element(component.lambdaGC.biasFactor.cbegin(), component.lambdaGC.biasFactor.cend());
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

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfCycles; ++currentCycle)
  {
    t1 = std::chrono::system_clock::now();
    estimation.setCurrentSample(currentCycle);

    totalNumberOfMolecules = std::transform_reduce(
        systems.begin(), systems.end(), 0uz, [](const size_t& acc, const size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules(); });
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

    for (size_t j = 0uz; j != numberOfStepsPerCycle; j++)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectedSecondSystem = systems[selectedSystemPair.second];

      size_t selectedComponent = selectedSystem.randomComponent(random);

      MC_Moves::performRandomMoveProduction(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                            fractionalMoleculeSystem, estimation.currentBin);

      // sample the occupancy within the innerloop
      selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
          selectedSystem.containsTheFractionalMolecule);
      selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
          selectedSecondSystem.containsTheFractionalMolecule);

      ++numberOfSteps;
    }

    //for (System& system : systems)
    //{
    //  system.checkCartesianPositions();
    //  system.checkMoleculeIds();
    //}

    // sample properties
    for (System& system : systems)
    {
      system.sampleProperties(estimation.currentBin, currentCycle);
    }

    for (System& system : systems)
    {
      // add the sample energy to the averages
      if (currentCycle % 10uz == 0uz || currentCycle % printEvery == 0uz)
      {
        std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
        std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
        system.currentEnergyStatus = molecularPressure.first;
        system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
        std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
        system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        std::print(stream, "{}", system.writeProductionStatusReportMC(currentCycle, numberOfCycles));
        std::print(stream, "\n");
        std::flush(stream);
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    // output properties to files
    for (System& system : systems)
    {
      if (system.propertyConventionalRadialDistributionFunction.has_value())
      {
        system.propertyConventionalRadialDistributionFunction->writeOutput(
            system.forceField, system.systemId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms,
            currentCycle);
      }

      if (system.propertyRadialDistributionFunction.has_value())
      {
        system.propertyRadialDistributionFunction->writeOutput(system.forceField, system.systemId,
                                                               system.simulationBox.volume,
                                                               system.totalNumberOfPseudoAtoms, currentCycle);
      }
      if (system.propertyDensityGrid.has_value())
      {
        system.propertyDensityGrid->writeOutput(system.systemId, system.simulationBox, system.forceField,
                                                system.frameworkComponents, system.components, currentCycle);
      }
    }

    // write binary-restart file
    if (currentCycle % printEvery == 0uz)
    {
      std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
      Archive<std::ofstream> archive(ofile);
      archive << *this;
      ofile.close();
      if (ofile)
      {
        std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
      }
    }
    t2 = std::chrono::system_clock::now();
    totalProductionSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

    continueProductionStage: ;

  }
}

void MonteCarlo::output()
{
  MCMoveCpuTime total;
  MCMoveCount countTotal;
  for (const System& system : systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_count;
  }

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    stream << system.runningEnergies.printMC("Running energies");
    std::print(stream, "\n\n\n\n");

    RunningEnergy recomputedEnergies = system.computeTotalEnergies();
    stream << recomputedEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    RunningEnergy drift = system.runningEnergies - recomputedEnergies;
    stream << drift.printMC("Monte-Carlo energy drift");
    std::print(stream, "\n\n\n\n");

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", system.writeMCMoveStatistics());

    std::print(stream, "Production run counting of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for (const Component& component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_statistics.writeMCMoveStatistics(numberOfSteps, 
                                                             component.componentId, component.name));
    }
    //std::print(stream, "{}", system.mc_moves_statistics.writeSystemStatistics(numberOfSteps));

    std::print(stream, "Production run counting of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", countTotal.writeAllSystemStatistics(numberOfSteps));

    std::print(stream, "\n\n");

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");

    for (const Component& component : system.components)
    {
      std::print(stream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(component.componentId, component.name));
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

    std::print(stream, "{}",
               system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.frameworkComponents,
                                                              system.components));
    std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());
    std::print(
        stream, "{}",
        system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swapableComponents, system.components));
    std::print(stream, "{}", system.averageLoadings.writeAveragesStatistics(system.components, system.frameworkMass));
  }
}

void MonteCarlo::log()
{
  for (System& system : systems)
  {
    // recomputedEnergies.printMC
    // drift.printMC
    system.runningEnergies.logMC(logs[system.systemId], "running_energies");

    RunningEnergy recomputedEnergies = system.computeTotalEnergies();
    recomputedEnergies.logMC(logs[system.systemId], "recomputed_energies");

    RunningEnergy drift = system.runningEnergies - recomputedEnergies;
    drift.logMC(logs[system.systemId], "drift");

    // system.writeMCMoveStatistics
    // for component: copmonent.mc_moves_count.writeComponentStatistics
    // system.mc_moves_count.writeSystemStatistics
    // countTotal.writeAllSystemStatistics
    // for componnet: component.mc_moves_cputime.writeMCMoveCPUTimeStatistics
    // system.mc_moves_cputime.writeMCMoveCPUTimeStatistics
    // total.writeMCMoveCPUTimeStatistics
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MonteCarlo& mc)
{
  archive << mc.versionNumber;

  archive << mc.numberOfCycles;
  archive << mc.numberOfSteps;
  archive << mc.numberOfInitializationCycles;
  archive << mc.numberOfEquilibrationCycles;

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

  archive << mc.totalInitializationSimulationTime;
  archive << mc.totalEquilibrationSimulationTime;
  archive << mc.totalProductionSimulationTime;
  archive << mc.totalSimulationTime;

  archive << static_cast<uint64_t>(0x6f6b6179); // magic number 'okay' in hex
 
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MonteCarlo& mc)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MonteCarlo' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> mc.numberOfCycles;
  archive >> mc.numberOfSteps;
  archive >> mc.numberOfInitializationCycles;
  archive >> mc.numberOfEquilibrationCycles;

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

  archive >> mc.totalInitializationSimulationTime;
  archive >> mc.totalEquilibrationSimulationTime;
  archive >> mc.totalProductionSimulationTime;
  archive >> mc.totalSimulationTime;

  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
  }
  std::cout << std::format("Magic number read correctly: {} vs {}\n", magicNumber, static_cast<uint64_t>(0x6f6b6179));
  return archive;
}
