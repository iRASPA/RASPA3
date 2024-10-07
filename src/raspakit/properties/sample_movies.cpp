module;

#ifdef USE_LEGACY_HEADERS
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <span>
#include <streambuf>
#include <string>
#include <vector>
#endif

module sample_movies;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <span>;
import <iostream>;
import <fstream>;
import <streambuf>;
import <filesystem>;
import <numbers>;
import <print>;
#endif

import double3;
import stringutils;
import atom;
import simulationbox;
import forcefield;
import units;
import skelement;

SampleMovie::SampleMovie(size_t systemId, size_t sampleEvery) : sampleEvery(sampleEvery)
{
  std::filesystem::create_directory("movies");
  std::ofstream stream(std::format("movies/movie.s{}.pdb", systemId));
}

void SampleMovie::update(const ForceField &forceField, size_t systemId, const SimulationBox simulationBox,
                         std::span<Atom> atomPositions, size_t currentCycle)
{
  if (currentCycle % sampleEvery == 0)
  {
    std::filesystem::create_directory("movies");
    std::ofstream stream(std::format("movies/movie.s{}.pdb", systemId), std::ios_base::app);

    std::print(stream, "MODEL {:>4}\n", modelNumber);
    std::print(stream, "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n", simulationBox.lengthA,
               simulationBox.lengthB, simulationBox.lengthC, simulationBox.angleAlpha * 180.0 / std::numbers::pi,
               simulationBox.angleBeta * 180.0 / std::numbers::pi, simulationBox.angleGamma * 180.0 / std::numbers::pi);

    for (int index = 1; const Atom &atom : atomPositions)
    {
      size_t atomicNumber = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].atomicNumber;
      std::string name = std::format("{:<4}", forceField.pseudoAtoms[static_cast<size_t>(atom.type)].name);
      std::string chemicalElement = PredefinedElements::predefinedElements[atomicNumber]._chemicalSymbol;
      std::print(stream,
                 "ATOM  {:>5} {:4}{:1}{:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4}{:>2}\n",
                 index, name.substr(0, 4), ' ', " ", ' ', 0, ' ', atom.position.x, atom.position.y, atom.position.z,
                 1.0, 0.0, ' ', chemicalElement);
      ++index;
    }
    stream << "ENDMDL\n";
    ++modelNumber;
  }
}
