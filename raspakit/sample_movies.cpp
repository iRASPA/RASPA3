module;

module sample_movies;

import <string>;
import <vector>;
import <iostream>;
import <fstream>;
import <streambuf>;
import <filesystem>;
import <numbers>;
import <print>;

import double3;
import stringutils;
import atom;
import simulationbox;
import forcefield;
import units;
import skelement;

SampleMovie::SampleMovie(size_t systemId, const ForceField& forceField, const SimulationBox& simulationBox, const std::vector<Atom>& atomPositions) :
    systemId(systemId),
    forceField(forceField),
    simulationBox(simulationBox),
    atomPositions(atomPositions)
{ 
}


SampleMovie::SampleMovie(SampleMovie&& s) noexcept:
    systemId(s.systemId),
    forceField(s.forceField),
    simulationBox(std::move(s.simulationBox)),
    atomPositions(std::move(s.atomPositions))
    //outputFile(std::move(s.outputFile))
{

}

void SampleMovie::initialize()
{
    if (sample)
    {
        std::filesystem::path cwd = std::filesystem::current_path();
        std::filesystem::path directoryName = cwd / std::format("Movies/System_{}/", systemId);
        //std::filesystem::path fileName = cwd / std::print("Movies/System_{}/movie_all_{}_{}.pdb",
        //    systemId, temperature, pressure * Units::PressureConversionFactor);
        std::filesystem::create_directories(directoryName);

        //outputFile = std::ofstream(fileName, std::ios::out);
    }
}

void SampleMovie::update(size_t cycle)
{
    if (sample)
    {
        if (cycle % writeEvery == 0)
        {
            //std::print(outputFile, "MODEL {}\n", modelNumber);
            //std::print(outputFile, "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n", simulationBox.lengthA, simulationBox.lengthB, simulationBox.lengthC,
            //    simulationBox.angleAlpha * 180.0 / std::numbers::pi, simulationBox.angleBeta * 180.0 / std::numbers::pi, simulationBox.angleGamma * 180.0 / std::numbers::pi);

            //for (int index = 1; const Atom & atom : atomPositions)
            //{
            //    size_t atomicNumber = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].atomicNumber;
            //    std::string chemicalElement = PredefinedElements::predefinedElements[atomicNumber]._chemicalSymbol;
            //    std::print(outputFile, "ATOM  {:>5} {:<4}{:1}{:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4}{:>2}\n",
            //               index, "name", ' ', " ", ' ', 0, ' ', atom.position.x, atom.position.y, atom.position.z, 1.0, 0.0, ' ', chemicalElement);
            //    ++index;
            //}
            //outputFile << "ENDMDL\n";
            ++modelNumber;
        }
    }
}

void SampleMovie::closeOutputFile()
{
    //if (sample)
    //{
    //    if (outputFile.is_open()) outputFile.close();
    //}
}

