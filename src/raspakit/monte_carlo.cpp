module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <print>
#include <ranges>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#endif

module monte_carlo;

#ifdef USE_STD_IMPORT
import std;
#endif

import stringutils;
import hardware_info;
import archive;
import system;
import framework;
import randomnumbers;
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
import mc_moves;
import mc_moves_probabilities;
import mc_moves_cputime;
import mc_moves_statistics;
import mc_moves_move_types;
import molecule;
import property_pressure;
import transition_matrix;
import interactions_ewald;
import equation_of_states;
import interpolation_energy_grid;

MonteCarlo::MonteCarlo() : outputToFiles(false), random(std::nullopt) {};

MonteCarlo::MonteCarlo(InputReader& reader) noexcept
    : outputToFiles(true),
      random(reader.randomSeed),
      numberOfCycles(reader.numberOfCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      writeRestartEvery(5000),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      systems(std::move(reader.systems)),
      outputJsons(systems.size()),
      estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
}

MonteCarlo::MonteCarlo(std::size_t numberOfCycles, std::size_t numberOfInitializationCycles,
                       std::size_t numberOfEquilibrationCycles, std::size_t printEvery,
                       std::size_t writeBinaryRestartEvery, std::size_t rescaleWangLandauEvery,
                       std::size_t optimizeMCMovesEvery, std::vector<System>& systems, RandomNumber& randomSeed,
                       std::size_t numberOfBlocks, bool outputToFiles)
    : outputToFiles(outputToFiles),
      random(randomSeed),
      numberOfCycles(numberOfCycles),
      numberOfInitializationCycles(numberOfInitializationCycles),
      numberOfEquilibrationCycles(numberOfEquilibrationCycles),
      printEvery(printEvery),
      writeRestartEvery(5000),
      writeBinaryRestartEvery(writeBinaryRestartEvery),
      rescaleWangLandauEvery(rescaleWangLandauEvery),
      optimizeMCMovesEvery(optimizeMCMovesEvery),
      systems(systems),
      outputJsons(systems.size()),
      estimation(numberOfBlocks, numberOfCycles)
{
}

System& MonteCarlo::randomSystem()
{
  return systems[std::size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MonteCarlo::run()
{
  switch (simulationStage)
  {
    case SimulationStage::Uninitialized:
      // this case only happens at first run, not when using a binart-restart file
      for (System& system : systems)
      {
        system.forceField.initializeAutomaticCutOff(system.simulationBox);
        system.forceField.initializeEwaldParameters(system.simulationBox);

        // switch the fractional molecule on in the first system, and off in all others
        if (system.systemId == 0uz)
          system.containsTheFractionalMolecule = true;
        else
          system.containsTheFractionalMolecule = false;

        // if the MC/MD hybrid move is on, make sure that interpolation-method include gradients
        if (system.mc_moves_probabilities.getProbability(MoveTypes::HybridMC) > 0.0 &&
            system.forceField.interpolationScheme == ForceField::InterpolationScheme::Polynomial)
        {
          system.forceField.interpolationScheme = ForceField::InterpolationScheme::Tricubic;
        }
      }
      if (outputToFiles)
      {
        createOutputFiles();
        writeOutputHeader();
      }

      createInterpolationGrids();
      break;
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

  if (outputToFiles)
  {
    output();
  }
}

void MonteCarlo::createOutputFiles()
{
  std::filesystem::create_directories("output");
  for (System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.txt", system.temperature, system.input_pressure, system.systemId);
    streams.emplace_back(fileNameString, std::ios::out);
    fileNameString =
        std::format("output/output_{}_{}.s{}.json", system.temperature, system.input_pressure, system.systemId);
    outputJsonFileNames.emplace_back(fileNameString);
  }
}

void MonteCarlo::writeOutputHeader()
{
  if (outputToFiles)
  {
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
      std::print(stream, "{}", system.writeNumberOfPseudoAtoms());
      std::print(stream, "{}", system.reactions.printStatus());

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
      outputJsons[system.systemId]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif

      outputJsons[system.systemId]["seed"] = random.seed;
      outputJsons[system.systemId]["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
      outputJsons[system.systemId]["initialization"]["units"] = Units::jsonStatus();
      outputJsons[system.systemId]["initialization"]["initialConditions"] = system.jsonSystemStatus();
      outputJsons[system.systemId]["initialization"]["forceField"] = system.forceField.jsonForceFieldStatus();
      outputJsons[system.systemId]["initialization"]["forceField"]["pseudoAtoms"] =
          system.forceField.jsonPseudoAtomStatus();
      outputJsons[system.systemId]["initialization"]["components"] = system.jsonComponentStatus();
      outputJsons[system.systemId]["initialization"]["reactions"] = system.reactions.jsonStatus();

      std::ofstream json(outputJsonFileNames[system.systemId]);
      json << outputJsons[system.systemId].dump(4);
    }
  }
}

void MonteCarlo::createInterpolationGrids()
{
  for (System& system : systems)
  {
    if (outputToFiles)
    {
      std::ostream stream(streams[system.systemId].rdbuf());
      system.createInterpolationGrids(stream);
    }
    else
    {
      std::ostringstream local;
      system.createInterpolationGrids(local);
    }
  }
}

void MonteCarlo::performCycle()
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

        selectedSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, selectedSystem.containsTheFractionalMolecule);

        selectedSecondSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
            selectedSecondSystem.containsTheFractionalMolecule);
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
  }
}

void MonteCarlo::initialize()
{
  std::chrono::system_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    if (outputToFiles)
    {
      std::ostream stream(streams[system.systemId].rdbuf());
      stream << system.runningEnergies.printMC("Recomputed from scratch");
      std::print(stream, "\n\n\n\n");

      system.writeRestartFile();
    }
  };

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
    t1 = std::chrono::system_clock::now();

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system.systemId].rdbuf());
          std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
          std::flush(stream);
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

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
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
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        for (System& system : systems)
        {
          system.writeRestartFile();
        }
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

  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
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

    performCycle();

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system.systemId].rdbuf());
          std::print(stream, "{}", system.writeEquilibrationStatusReportMC(currentCycle, numberOfEquilibrationCycles));
          std::flush(stream);
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

        if (outputToFiles)
        {
          std::filesystem::create_directories("bias_factors");
          for (Component& component : system.components)
          {
            component.lambdaGC.writeBiasingFile(
                std::format("bias_factors/lambda_bias_{}.s{}.json", component.name, system.systemId));
          }
        }
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
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
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        for (System& system : systems)
        {
          system.writeRestartFile();
        }
      }
    }

    t2 = std::chrono::system_clock::now();

    totalEquilibrationSimulationTime += (t2 - t1);
    totalSimulationTime += (t2 - t1);

  continueEquilibrationStage:;
  }
}

void MonteCarlo::production()
{
  double minBias{0uz};
  std::chrono::system_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
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

  if (outputToFiles)
  {
    std::filesystem::create_directories("bias_factors");
    for (System& system : systems)
    {
      for (Component& component : system.components)
      {
        component.lambdaGC.writeBiasingFile(
            std::format("bias_factors/lambda_bias_{}.s{}.json", component.name, system.systemId));
      }
    }
  }

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfCycles; ++currentCycle)
  {
    t1 = std::chrono::system_clock::now();

    estimation.setCurrentSample(currentCycle);

    performCycle();

    for (System& system : systems)
    {
      system.sampleProperties(estimation.currentBin, currentCycle);
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
        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system.systemId].rdbuf());
          std::string status_line{std::format("Current cycle: {} out of {}\n", currentCycle, numberOfCycles)};
          std::print(stream, "{}", system.writeProductionStatusReportMC(status_line));
          std::flush(stream);
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

    if (outputToFiles)
    {
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
                                                  system.framework, system.components, currentCycle);
        }
        if (system.averageEnergyHistogram.has_value())
        {
          system.averageEnergyHistogram->writeOutput(system.systemId, currentCycle);
        }
        if (system.averageNumberOfMoleculesHistogram.has_value())
        {
          system.averageNumberOfMoleculesHistogram->writeOutput(system.systemId, system.components, currentCycle);
        }
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
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
    }

    if (currentCycle % writeRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        for (System& system : systems)
        {
          system.writeRestartFile();
        }
      }
    }

    t2 = std::chrono::system_clock::now();
    switch (simulationStage)
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
        throw std::runtime_error("Unexpected simulation stage");
    }
    totalSimulationTime += (t2 - t1);

  continueProductionStage:;
  }

  // Carry out last write of properties after simulation has finished
  if (outputToFiles)
  {
    // Write out last status
    for (System& system : systems)
    {
      system.writeRestartFile();

      system.forceField.initializeAutomaticCutOff(system.simulationBox);

      std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
      system.currentEnergyStatus = molecularPressure.first;
      system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;

      system.loadings =
          Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

      std::ostream stream(streams[system.systemId].rdbuf());

      std::print(stream, "\n");
      std::print(stream, "===============================================================================\n");
      std::print(stream, "                             Simulation finished!\n");
      std::print(stream, "===============================================================================\n");
      std::print(stream, "\n");

      std::string status_line{std::format("Final state after {} cycles\n", numberOfCycles)};
      std::print(stream, "{}", system.writeProductionStatusReportMC(status_line));
      std::flush(stream);

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
                                                system.framework, system.components, currentCycle);
      }
      if (system.averageEnergyHistogram.has_value())
      {
        system.averageEnergyHistogram->writeOutput(system.systemId, currentCycle);
      }
      if (system.averageNumberOfMoleculesHistogram.has_value())
      {
        system.averageNumberOfMoleculesHistogram->writeOutput(system.systemId, system.components, currentCycle);
      }
    }
  }
}

void MonteCarlo::output()
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

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

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

    for (const Component& component : system.components)
    {
      std::print(stream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());

    std::print(stream, "Production run CPU timings of the MC moves summed over systems and components\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
    if (totalGridCreationTime.count() > 1e-4)
    {
      std::print(stream, "Grid creation time:             {:14f} [s]\n", totalGridCreationTime.count());
    }
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

    std::print(
        stream, "{}",
        system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swappableComponents, system.components));
    std::print(stream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));

    // json statistics
    outputJsons[system.systemId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    outputJsons[system.systemId]["output"]["recomputedEnergies"] = recomputedEnergies.jsonMC();
    outputJsons[system.systemId]["output"]["drift"] = drift.jsonMC();

    outputJsons[system.systemId]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    // outputJsons[system.systemId]["output"]["MCMoveStatistics"]["summedOverAllSystems"] =
    //    countTotal.jsonAllSystemStatistics(numberOfSteps);

    outputJsons[system.systemId]["output"]["cpuTimings"]["summedSystemsAndComponents"] =
        total.jsonOverallMCMoveCPUTimeStatistics(totalProductionSimulationTime);
    outputJsons[system.systemId]["output"]["cpuTimings"]["gridCreation"] = totalGridCreationTime.count();
    outputJsons[system.systemId]["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
    outputJsons[system.systemId]["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
    outputJsons[system.systemId]["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
    outputJsons[system.systemId]["output"]["cpuTimings"]["total"] = totalSimulationTime.count();
    outputJsons[system.systemId]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();

    outputJsons[system.systemId]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    outputJsons[system.systemId]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();
    outputJsons[system.systemId]["properties"]["averageEnthalpy"] =
        system.averageEnthalpiesOfAdsorption.jsonAveragesStatistics(system.swappableComponents, system.components);

    for (const Component& component : system.components)
    {
      // outputJsons[system.systemId]["output"]["MCMoveStatistics"][component.name]["percentage"] =
      //     component.mc_moves_statistics.jsonMCMoveStatistics(numberOfSteps);
      outputJsons[system.systemId]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }

    std::ofstream json(outputJsonFileNames[system.systemId]);
    json << outputJsons[system.systemId].dump(4);
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MonteCarlo& mc)
{
  archive << mc.versionNumber;

  archive << mc.outputToFiles;
  archive << mc.random;

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
  archive << mc.fractionalMoleculeSystem;

  archive << mc.estimation;

  // not written: totalGridCreationTime

  archive << mc.totalInitializationSimulationTime;
  archive << mc.totalEquilibrationSimulationTime;
  archive << mc.totalProductionSimulationTime;
  archive << mc.totalSimulationTime;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MonteCarlo& mc)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MonteCarlo' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> mc.outputToFiles;
  archive >> mc.random;

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
  archive >> mc.fractionalMoleculeSystem;

  archive >> mc.estimation;

  // not written: totalGridCreationTime

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
