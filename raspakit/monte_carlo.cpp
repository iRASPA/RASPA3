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

        system.createInitialMolecules();

        system.averageEnthalpiesOfAdsorption.resize(system.swapableComponents.size());
    }
    for (System& system : systems)
    {
        system.createOutputFile();
        system.writeOutputHeader();
        system.writeOutputHeaderHardware();
        system.writeToOutputFile(Units::printStatus());
        system.writeToOutputFile(system.simulationBox.printParameters());
        system.writeToOutputFile(system.forceField.printPseudoAtomStatus());
        system.writeToOutputFile(system.forceField.printForceFieldStatus());
        system.writeComponentStatus();
    }

    for (System& system : systems)
    {
        system.computeTotalEnergies();
        system.writeToOutputFile(system.runningEnergies.print("Recomputed from scratch"));
    };

    for (size_t i = 0; i != numberOfInitializationCycles; i++)
    {
        for (size_t j = 0; j != systems.size(); ++j)
        {
            System& selectedSystem = randomSystem();

            size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
            for (size_t k = 0; k != numberOfSteps; k++)
            {
                size_t selectedComponent = selectedSystem.randomComponent();
                particleMoves.performRandomMove(selectedSystem, selectedComponent);
            }
        }

        if (i % printEvery == 0)
        {
            for (System& system : systems)
            {
                system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

                system.writeInitializationStatusReport(i, numberOfInitializationCycles);
            }

            std::cout << "Init iteration: " << i << std::endl;

        }
    }
}

void MonteCarlo::equilibrate()
{
    for (System& system : systems)
    {
        system.computeTotalEnergies();
        system.writeToOutputFile(system.runningEnergies.print("Recomputed from scratch"));

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

            size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
            for (size_t k = 0; k != numberOfSteps; k++)
            {
                size_t selectedComponent = selectedSystem.randomComponent();
                particleMoves.performRandomMove(selectedSystem, selectedComponent);
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
                system.loadings = Loadings(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

                system.writeEquilibrationStatusReport(i, numberOfEquilibrationCycles);
                for(Component &component : system.components)
                {
                   if(component.hasFractionalMolecule)
                   {
                     component.lambda.WangLandauIteration(Lambda::WangLandauPhase::AdjustBiasingFactors);
                   }
                }
            }

            std::cout << "Equilibration iteration: " << i << std::endl;

        }
    }
}

void MonteCarlo::production()
{
    for (System& system : systems)
    {
        system.computeTotalEnergies(); 
        system.writeToOutputFile(system.runningEnergies.print("Recomputed from scratch"));
        system.sampleMovie.initialize();

        system.clearMoveStatistics();
        system.clearTimingStatistics();

        for(Component &component : system.components)
        {
          component.clearMoveStatistics();
          component.clearTimingStatistics();
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
            
            size_t numberOfSteps = std::max(selectedSystem.numberOfMolecules(), size_t(20)) * selectedSystem.numerOfAdsorbateComponents();
            for (size_t k = 0; k != numberOfSteps; k++)
            {
                size_t selectedComponent = selectedSystem.randomComponent();
                particleMoves.performRandomMoveProduction(selectedSystem, selectedComponent, estimation.currentBin);
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
                system.writeProductionStatusReport(i, numberOfCycles);
            }

            std::cout << "iteration: " << i << std::endl;
        }

        for (System& system : systems)
        {
            system.sampleMovie.update(i);
        }
    }

}

void MonteCarlo::output()
{
    for (System& system : systems)
    {
        RunningEnergy runningEnergies = system.runningEnergies;
        system.writeToOutputFile(runningEnergies.print("Running energies"));
        
        system.computeTotalEnergies();
        RunningEnergy recomputedEnergies = system.runningEnergies;
        system.writeToOutputFile(recomputedEnergies.print("Recomputed from scratch"));
        
        RunningEnergy drift = runningEnergies - recomputedEnergies;
        system.writeToOutputFile(drift.print("Monte-Carlo energy drift"));

        system.writeToOutputFile("\n\n");
        
        system.writeToOutputFile("Monte-Carlo moves statistics\n");
        system.writeToOutputFile("===============================================================================\n\n");
        
        system.writeMCMoveStatistics();

        system.writeToOutputFile("Production run CPU timings of the MC moves\n");
        system.writeToOutputFile("===============================================================================\n\n");
         system.writeToOutputFile(system.writeCPUTimeStatistics());
        std::chrono::duration<double> totalSimulationTime = (t2 - t1);
        system.writeToOutputFile(std::print("\nProduction simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count()));

        system.writeToOutputFile(system.writeEnergyAveragesStatistics());
        system.writeToOutputFile(system.writePressureAveragesStatistics());
        system.writeToOutputFile(system.writeEnthalpyOfAdsorption());
    }
}
void MonteCarlo::cleanup()
{
    for (System& system : systems)
    {
        system.closeOutputFile();
        system.sampleMovie.closeOutputFile();
    }
}
