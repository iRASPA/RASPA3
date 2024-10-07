module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <numeric>
#include <optional>
#include <print>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

module monte_carlo_transition_matrix;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <algorithm>;
import <numeric>;
import <chrono>;
import <vector>;
import <span>;
import <string>;
import <optional>;
import <fstream>;
import <filesystem>;
import <tuple>;
import <ios>;
import <complex>;
import <print>;
#endif

import stringutils;
import hardware_info;
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
import property_pressure;
import transition_matrix;
import interactions_ewald;

MonteCarloTransitionMatrix::MonteCarloTransitionMatrix(InputReader& reader) noexcept
    : numberOfCycles(reader.numberOfCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      systems(std::move(reader.systems)),
      random(reader.randomSeed),
      estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
}

System& MonteCarloTransitionMatrix::randomSystem()
{
  return systems[size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MonteCarloTransitionMatrix::run()
{
  initialize();
  equilibrate();
  t1 = std::chrono::system_clock::now();
  production();
  t2 = std::chrono::system_clock::now();

  output();
  cleanup();
}

void MonteCarloTransitionMatrix::initialize()
{
  for (System& system : systems)
  {
    system.tmmc.initialize();

    std::string directoryNameString = std::format("output/system_{}/", system.systemId);
    std::filesystem::path directoryName{directoryNameString};
    std::filesystem::create_directories(directoryName);

    std::string fileNameString =
        std::format("output/system_{}/output_{}_{}.txt", system.systemId, system.temperature, system.input_pressure);
    streams.emplace_back(fileNameString, std::ios::out);
  }

  for (System& system : systems)
  {
    // switch the fractional molecule on in the first system, and off in all others
    if (system.systemId == 0)
      system.containsTheFractionalMolecule = true;
    else
      system.containsTheFractionalMolecule = false;
  }

  for (const System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    std::print(stream, "{}", system.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());
    std::print(stream, "{}", system.writeSystemStatus());
    std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(stream, "{}", system.forceField.printForceFieldStatus());
    std::print(stream, "{}", system.writeComponentStatus());
    std::print(stream, "{}", system.reactions.printStatus());
  }

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();

    std::ostream stream(streams[system.systemId].rdbuf());
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");
  };

  for (size_t i = 0; i != numberOfInitializationCycles; i++)
  {
    for (size_t j = 0; j != systems.size(); ++j)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      size_t numberOfSteps =
          std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent(random);
        MC_Moves::performRandomMove(random, selectedSystem, selectSecondSystem, selectedComponent,
                                    fractionalMoleculeSystem);

        size_t N = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
        selectedSystem.tmmc.updateHistogram(N);
        selectedSystem.tmmc.numberOfSteps++;
      }
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        std::print(stream, "{}", system.writeInitializationStatusReport(i, numberOfInitializationCycles));
      }
    }

    if (i % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }
  }
}

void MonteCarloTransitionMatrix::equilibrate()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
    }
  };

  for (size_t j = 0; j != systems.size(); ++j)
  {
    systems[j].tmmc.numberOfSteps = 0;
  }

  for (size_t i = 0; i != numberOfEquilibrationCycles; i++)
  {
    for (size_t j = 0; j != systems.size(); ++j)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      size_t numberOfSteps =
          std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent(random);
        MC_Moves::performRandomMove(random, selectedSystem, selectSecondSystem, selectedComponent,
                                    fractionalMoleculeSystem);

        size_t N = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
        selectedSystem.tmmc.updateHistogram(N);
        selectedSystem.tmmc.numberOfSteps++;
        selectedSystem.tmmc.adjustBias();
      }
    }

    for (System& system : systems)
    {
      for (Component& component : system.components)
      {
        if (component.hasFractionalMolecule)
        {
          if (system.containsTheFractionalMolecule)
          {
            component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
                                                   system.containsTheFractionalMolecule);
          }
        }
      }
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings =
            Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        std::print(stream, "{}", system.writeEquilibrationStatusReportMC(i, numberOfEquilibrationCycles));
        for (Component& component : system.components)
        {
          if (component.hasFractionalMolecule)
          {
            component.lambdaGC.WangLandauIteration(
                PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
                system.containsTheFractionalMolecule);
          }
        }
      }
    }

    if (i % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }
  }
}

void MonteCarloTransitionMatrix::production()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.runningEnergies = system.computeTotalEnergies();
    stream << system.runningEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");
    // system.sampleMovie.initialize();

    system.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();

    for (Component& component : system.components)
    {
      component.mc_moves_statistics.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();
      if (component.hasFractionalMolecule)
      {
        component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                               system.containsTheFractionalMolecule);
      }
    }
  };

  for (size_t j = 0; j != systems.size(); ++j)
  {
    systems[j].tmmc.numberOfSteps = 0;
  }

  for (size_t i = 0; i != numberOfCycles; i++)
  {
    estimation.setCurrentSample(i);

    for (size_t j = 0; j != systems.size(); ++j)
    {
      std::pair<size_t, size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      size_t numberOfSteps =
          std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent(random);
        MC_Moves::performRandomMoveProduction(random, selectedSystem, selectSecondSystem, selectedComponent,
                                              fractionalMoleculeSystem, estimation.currentBin);

        size_t N = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
        selectedSystem.tmmc.updateHistogram(N);
        selectedSystem.tmmc.numberOfSteps++;
        selectedSystem.tmmc.adjustBias();
      }
    }

    for (System& system : systems)
    {
      std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
      std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
      system.currentEnergyStatus = molecularPressure.first;
      system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
      std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
      system.mc_moves_cputime.energyPressureComputation += (time2 - time1);

      // add the sample energy to the averages
      if (i % 10 == 0 || i % printEvery == 0)
      {
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }

      system.sampleProperties(estimation.currentBin, estimation.currentBin);
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        std::print(stream, "{}", system.writeProductionStatusReportMC(i, numberOfCycles));
      }
    }

    if (i % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    // for (System& system : systems)
    //{
    // system.sampleMovie.update(i);
    //}
  }

  // Write the collection matrix
  for (System& system : systems)
  {
    system.tmmc.writeStatistics();
  }
}

void MonteCarloTransitionMatrix::output()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    RunningEnergy runningEnergies = system.runningEnergies;
    stream << runningEnergies.printMC("Running energies");
    std::print(stream, "\n\n\n\n");

    system.runningEnergies = system.computeTotalEnergies();
    RunningEnergy recomputedEnergies = system.runningEnergies;
    stream << recomputedEnergies.printMC("Recomputed from scratch");
    std::print(stream, "\n\n\n\n");

    RunningEnergy drift = runningEnergies - recomputedEnergies;
    stream << drift.printMC("Monte-Carlo energy drift");
    std::print(stream, "\n\n\n\n");

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");

    std::print(stream, "{}", system.writeMCMoveStatistics());
    // std::print(stream, system.lambda.writeAveragesStatistics(system.beta));

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");
    system.writeCPUTimeStatistics(stream);
    std::chrono::duration<double> totalSimulationTime = (t2 - t1);
    std::print(stream, "\nProduction simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count());

    // std::print(stream, system.averageEnergies.writeAveragesStatistics(system.components));
    // std::print(stream, system.averagePressure.writeAveragesStatistics());
    // std::print(stream, system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swappableComponents,
    // system.components)); std::print(stream, system.averageLoadings.writeAveragesStatistics(system.components,
    // system.frameworkMass));
  }
}
void MonteCarloTransitionMatrix::cleanup()
{
  // for (System& system : systems)
  //{
  //   //system.sampleMovie.closeOutputFile();
  // }
}
