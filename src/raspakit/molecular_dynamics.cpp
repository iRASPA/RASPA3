module;

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

module molecular_dynamics;

#ifndef USE_LEGACY_HEADERS
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
import property_msd;
import property_vacf;
import mc_moves;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_cputime;
import mc_moves_statistics;
import property_pressure;
import transition_matrix;
import interactions_ewald;
import equation_of_states;
import integrators;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import interpolation_energy_grid;

MolecularDynamics::MolecularDynamics() : random(std::nullopt) {};

MolecularDynamics::MolecularDynamics(InputReader& reader) noexcept
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

System& MolecularDynamics::randomSystem()
{
  return systems[std::size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MolecularDynamics::run()
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
}

void MolecularDynamics::createOutputFiles()
{
  std::filesystem::create_directories("output");
  for (System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.txt", system.temperature, system.input_pressure, system.systemId);
    streams.emplace_back(fileNameString, std::ios::out);
  }
}

void MolecularDynamics::createInterpolationGrids()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.createInterpolationGrids(stream);
  }
}

void MolecularDynamics::initialize()
{
  std::size_t totalNumberOfMolecules{0uz};
  std::size_t totalNumberOfComponents{0uz};
  std::size_t numberOfStepsPerCycle{0uz};

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

    // if the MC/MD hybrid move is on, make sure that interpolation-method include gradients
    if (system.forceField.interpolationScheme == ForceField::InterpolationScheme::Polynomial)
    {
      system.forceField.interpolationScheme = ForceField::InterpolationScheme::Tricubic;
    }
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
  }

  createInterpolationGrids();

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy =
        Integrators::computeTranslationalKineticEnergy(system.moleculeData);
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components);

    std::ostream stream(streams[system.systemId].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");
  };

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; currentCycle++)
  {
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
      System& selectSecondSystem = systems[selectedSystemPair.second];

      std::size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem);

      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
        }
      }
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::print(stream, "{}\n\n\n\n", system.runningEnergies.printMC(""));
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

  continueInitializationStage:;
  }
}

void MolecularDynamics::equilibrate()
{
  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
    Integrators::initializeVelocities(random, system.moleculeData, system.components, system.temperature);

    Integrators::removeCenterOfMassVelocityDrift(system.moleculeData);
    if (system.thermostat.has_value())
    {
      if (!system.framework.has_value() && system.numberOfMolecules() > 1uz)
      {
        system.translationalCenterOfMassConstraint = 3;
        system.thermostat->translationalCenterOfMassConstraint = 3;
      }
      system.thermostat->initialize(random);
    }

    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy =
        Integrators::computeTranslationalKineticEnergy(system.moleculeData);
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components);
    if (system.thermostat.has_value())
    {
      system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
    }
    system.referenceEnergy = system.runningEnergies.conservedEnergy();

    stream << system.runningEnergies.printMD("Recomputed from scratch", system.referenceEnergy);
    std::print(stream, "\n\n\n\n");

    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }
  };

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle)
  {
    for (System& system : systems)
    {
      system.runningEnergies = Integrators::velocityVerlet(
          system.moleculeData, system.spanOfMoleculeAtoms(), system.components, system.timeStep, system.thermostat,
          system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, system.eik_x, system.eik_y,
          system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
          system.numberOfMoleculesPerComponent);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift +=
          std::abs(Units::EnergyToKelvin * (system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::print(stream, "{}", system.writeEquilibrationStatusReportMD(currentCycle, numberOfEquilibrationCycles));
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
  continueEquilibrationStage:;
  }
}

void MolecularDynamics::production()
{
  std::chrono::system_clock::time_point t1, t2;
  double minBias{0.0};

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy =
        Integrators::computeTranslationalKineticEnergy(system.moleculeData);
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components);
    if (system.thermostat.has_value())
    {
      system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
    }
    system.referenceEnergy = system.runningEnergies.conservedEnergy();

    stream << system.runningEnergies.printMD("Recomputed from scratch", system.referenceEnergy);
    std::print(stream, "\n");

    system.mc_moves_statistics.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();
    integratorsCPUTime.clearTimingStatistics();

    system.accumulatedDrift = 0.0;

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

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfCycles; ++currentCycle)
  {
    t1 = std::chrono::system_clock::now();

    estimation.setCurrentSample(currentCycle);

    for ([[maybe_unused]] System& system : systems)
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

    for (System& system : systems)
    {
      system.runningEnergies = Integrators::velocityVerlet(
          system.moleculeData, system.spanOfMoleculeAtoms(), system.components, system.timeStep, system.thermostat,
          system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, system.eik_x, system.eik_y,
          system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
          system.numberOfMoleculesPerComponent);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift +=
          std::abs(Units::EnergyToKelvin * (system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
    }

    // sample properties
    for (System& system : systems)
    {
      system.sampleProperties(estimation.currentBin, currentCycle);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        std::print(stream, "{}", system.writeProductionStatusReportMD(currentCycle, numberOfCycles));
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
                                                system.framework, system.components, currentCycle);
      }

      if (system.propertyMSD.has_value())
      {
        system.propertyMSD->writeOutput(system.systemId, system.components, system.numberOfIntegerMoleculesPerComponent,
                                        system.timeStep, currentCycle);
      }

      if (system.propertyVACF.has_value())
      {
        system.propertyVACF->writeOutput(system.systemId, system.components,
                                         system.numberOfIntegerMoleculesPerComponent, system.timeStep, currentCycle);
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
    totalSimulationTime += (t2 - t1);
  continueProductionStage:;
  }

  // output properties to files
  for (System& system : systems)
  {
    if (system.propertyConventionalRadialDistributionFunction.has_value())
    {
      system.propertyConventionalRadialDistributionFunction->writeOutput(system.forceField, system.systemId,
                                                                         system.simulationBox.volume,
                                                                         system.totalNumberOfPseudoAtoms, currentCycle);
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

    if (system.propertyMSD.has_value())
    {
      system.propertyMSD->writeOutput(system.systemId, system.components, system.numberOfIntegerMoleculesPerComponent,
                                      system.timeStep, currentCycle);
    }

    if (system.propertyVACF.has_value())
    {
      system.propertyVACF->writeOutput(system.systemId, system.components, system.numberOfIntegerMoleculesPerComponent,
                                       system.timeStep, currentCycle);
    }
  }
}

void MolecularDynamics::output()
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

    std::print(stream, "\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "                             Simulation finished!\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "\n");

    std::string status_line{std::format("Final state after {} cycles\n\n", numberOfCycles)};

    std::print(stream, "Production run CPU timings of the MD simulation\n");
    std::print(stream, "===============================================================================\n\n");

    for (const Component& component : system.components)
    {
      std::print(stream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(component.componentId, component.name));
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(stream, "{}", integratorsCPUTime.writeIntegratorsCPUTimeStatistics(totalSimulationTime));
    std::print(stream, "\n\n");

    std::print(
        stream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));

    std::print(stream, "Temperature averages and statistics:\n");
    std::print(stream, "===============================================================================\n\n");
    std::print(stream, "{}", system.averageTemperature.writeAveragesStatistics("Total"));
    std::print(stream, "{}", system.averageTranslationalTemperature.writeAveragesStatistics("Translational"));
    std::print(stream, "{}", system.averageRotationalTemperature.writeAveragesStatistics("Rotational"));

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
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MolecularDynamics& mc)
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

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MolecularDynamics& mc)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MolecularDynamics' at line {} in file {}\n",
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

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
  }
  std::cout << std::format("Magic number read correctly: {} vs {}\n", magicNumber,
                           static_cast<std::uint64_t>(0x6f6b6179));
  return archive;
}
