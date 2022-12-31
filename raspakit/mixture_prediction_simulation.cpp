module;

module mixture_prediction_simulation;

import <vector>;
import <span>;
import <cmath>;
import <string>;
import <iostream>;
import <fstream>;
import <limits>;
import <filesystem>;
import <algorithm>;
import <numeric>;
import <sstream>;
import <chrono>;

import print;
import input_reader;
import component;
import system;
import simulationbox;
import multi_site_isotherm;
import mixture_prediction;


MixturePredictionSimulation::MixturePredictionSimulation(InputReader &inputReader):
    systems(std::move(inputReader.systems))
{
  for(System &system: systems)
  {
    std::string directoryNameString = std::print("Output/System_{}/", system.systemId);
    std::filesystem::path directoryName{ directoryNameString };
    std::filesystem::create_directories(directoryName);
  }
}

void MixturePredictionSimulation::run()
{
  for(System &system: systems)
  {
    MixturePrediction mixture(system);

    std::string fileNameString = std::print("Output/System_{}/output_{}_{}.data",
        system.systemId, system.simulationBox.temperature, system.simulationBox.input_pressure);
    std::ofstream fstream(fileNameString, std::ios::out );
    std::ostream stream(fstream.rdbuf());
    //std::ostream stream(std::cout.rdbuf());

    system.writeOutputHeader(stream);
    system.writeOutputHeaderHardware(stream);
    mixture.writeHeader(stream);

    mixture.createPureComponentsPlotScript();
    mixture.createMixturePlotScript();
    mixture.createMixtureAdsorbedMolFractionPlotScript();
    mixture.createPlotScript();

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

    mixture.run(stream);

    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

    std::chrono::duration<double> totalSimulationTime = (t2 - t1);
    std::print(stream, "\nMixture prediction simulation time: {:14f} [s]\n\n\n", totalSimulationTime.count());
  }
}

