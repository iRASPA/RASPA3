module;

module monte_carlo;

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
import print;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import atom;
import double3;
import double3x3;
import lambda;
import property_widom;
import property_dudlambda;
import property_simulationbox;
import property_energy;
import property_loading;
import property_enthalpy;

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

MonteCarlo::MonteCarlo(InputReader& reader) noexcept : 
    numberOfCycles(reader.numberOfCycles),
    numberOfInitializationCycles(reader.numberOfInitializationCycles),
    numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
    printEvery(reader.printEvery),
    systems(std::move(reader.systems)),
    estimation(reader.numberOfBlocks, reader.numberOfCycles)
{
    
}

System& MonteCarlo::randomSystem()
{
  return systems[size_t(RandomNumber::Uniform() * static_cast<double>(systems.size()))];
}

void MonteCarlo::run()
{
  initialize();
  equilibrate();
  t1 = std::chrono::system_clock::now();
  production();
  t2 = std::chrono::system_clock::now();

  output();
  cleanup();
}

void MonteCarlo::initialize()
{
  for (System &system: systems)
  {
    system.registerEwaldFourierEnergySingleIon(double3(0.0, 0.0, 0.0), 0.0);
    system.removeRedundantMoves();
    system.determineSwapableComponents();
    system.determineFractionalComponents();
    system.rescaleMoveProbabilities();
    system.rescaleMolarFractions();
    system.computeComponentFluidProperties();
    system.computeFrameworkDensity();
    system.computeNumberOfPseudoAtoms();

    double3 perpendicularWidths = system.simulationBox.perpendicularWidths();
    system.forceField.initializeEwaldParameters(perpendicularWidths);

    system.createInitialMolecules();

    system.averageEnthalpiesOfAdsorption.resize(system.swapableComponents.size());

    std::string directoryNameString = std::print("Output/System_{}/", system.systemId);
    std::filesystem::path directoryName{ directoryNameString };
    std::filesystem::create_directories(directoryName);

    std::string fileNameString = std::print("Output/System_{}/output_{}_{}.data",
        system.systemId, system.temperature, system.input_pressure);
    streams.emplace_back(fileNameString, std::ios::out );
  }

  for (const System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.writeOutputHeader(stream);
    system.writeOutputHeaderHardware(stream);
    Units::printStatus(stream);
    system.simulationBox.printParameters(stream);
    system.forceField.printPseudoAtomStatus(stream);
    system.forceField.printForceFieldStatus(stream);
    system.writeComponentStatus(stream);
  }

  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.precomputeTotalRigidEnergy();
    system.computeTotalEnergies();
    system.runningEnergies.print(stream, "Recomputed from scratch");
  };
  
  for (size_t i = 0; i != numberOfInitializationCycles; i++)
  {
    for (size_t j = 0; j != systems.size(); ++j)
    {
      System& selectedSystem = randomSystem();
      System& selectSecondSystem = systems[(selectedSystem.systemId + 1) % systems.size()];

      size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent();
        particleMoves.performRandomMove(selectedSystem, selectSecondSystem, selectedComponent);
      }
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
        system.writeInitializationStatusReport(stream, i, numberOfInitializationCycles);
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

void MonteCarlo::equilibrate()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.computeTotalEnergies();
    system.runningEnergies.print(stream, "Recomputed from scratch");

    for(Component &component : system.components)
    {
      component.lambda.WangLandauIteration(Lambda::WangLandauPhase::Initialize);
    }
  };

  for (size_t i = 0; i != numberOfEquilibrationCycles; i++)
  {
    for (size_t j = 0; j != systems.size(); ++j)
    {
      System& selectedSystem = randomSystem();
      System& selectSecondSystem = systems[(selectedSystem.systemId + 1) % systems.size()];

      size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent();
        particleMoves.performRandomMove(selectedSystem, selectSecondSystem, selectedComponent);
      }
    }

    for (System& system : systems)
    {
      for(Component &component : system.components)
      {
        if(component.hasFractionalMolecule)
        {
          component.lambda.WangLandauIteration(Lambda::WangLandauPhase::Sample);
        }
      }
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());

        system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        system.writeEquilibrationStatusReport(stream, i, numberOfEquilibrationCycles);
        for(Component &component : system.components)
        {
          if(component.hasFractionalMolecule)
          {
            component.lambda.WangLandauIteration(Lambda::WangLandauPhase::AdjustBiasingFactors);
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

void MonteCarlo::production()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    system.computeTotalEnergies(); 
    system.runningEnergies.print(stream, "Recomputed from scratch");
    //system.sampleMovie.initialize();

    system.clearMoveStatistics();
    system.clearTimingStatistics();

    for(Component &component : system.components)
    {
      component.mc_moves_probabilities.clearMoveStatistics();
      component.mc_moves_probabilities.clearTimingStatistics();
      if(component.hasFractionalMolecule)
      {
        component.lambda.WangLandauIteration(Lambda::WangLandauPhase::Finalize);
      }
    }
  };

  
  for (size_t i = 0; i != numberOfCycles; i++)
  {
    estimation.setCurrentSample(i);

    for (size_t j = 0; j != systems.size(); ++j)
    {
      System& selectedSystem = randomSystem();
      System& selectSecondSystem = systems[(selectedSystem.systemId + 1) % systems.size()];
      
      size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
      for (size_t k = 0; k != numberOfSteps; k++)
      {
        size_t selectedComponent = selectedSystem.randomComponent();
        particleMoves.performRandomMoveProduction(selectedSystem, selectSecondSystem, selectedComponent, estimation.currentBin);
      }
    }

    for (System& system : systems)
    {
      std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
      std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
      system.currentEnergyStatus = molecularPressure.first;
      system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
      std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
      system.cpuTime_Pressure += (time2 - time1);

      // add the sample energy to the averages
      if (i % 10 == 0 || i % printEvery == 0)
      {
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }

      system.sampleProperties(estimation.currentBin);
    }

    if (i % printEvery == 0)
    {
      for (System& system : systems)
      {
        std::ostream stream(streams[system.systemId].rdbuf());
        system.writeProductionStatusReport(stream, i, numberOfCycles);
      }
    }

    if (i % optimizeMCMovesEvery == 0)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    //for (System& system : systems)
    //{
        //system.sampleMovie.update(i);
    //}
  }

}

void MonteCarlo::output()
{
  for (System& system : systems)
  {
    std::ostream stream(streams[system.systemId].rdbuf());

    RunningEnergy runningEnergies = system.runningEnergies;
    runningEnergies.print(stream, "Running energies");
    
    system.computeTotalEnergies();
    RunningEnergy recomputedEnergies = system.runningEnergies;
    recomputedEnergies.print(stream, "Recomputed from scratch");
    
    RunningEnergy drift = runningEnergies - recomputedEnergies;
    drift.print(stream, "Monte-Carlo energy drift");

    std::print(stream, "\n\n");

    std::print(stream, "Monte-Carlo moves statistics\n");
    std::print(stream, "===============================================================================\n\n");
    
    system.writeMCMoveStatistics(stream);

    std::print(stream, "Production run CPU timings of the MC moves\n");
    std::print(stream, "===============================================================================\n\n");
    system.writeCPUTimeStatistics(stream);
    std::chrono::duration<double> totalSimulationTime = (t2 - t1);
    std::print(stream, "\nProduction simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count());

    system.writeEnergyAveragesStatistics(stream);
    system.writePressureAveragesStatistics(stream);
    system.writeEnthalpyOfAdsorption(stream);
  }
}
void MonteCarlo::cleanup()
{
  //for (System& system : systems)
  //{
  //  //system.sampleMovie.closeOutputFile();
  //}
}
