module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <source_location>
#include <span>
#include <streambuf>
#include <string>
#include <vector>
#endif

module sample_movies;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3;
import stringutils;
import atom;
import simulationbox;
import forcefield;
import units;
import skelement;

SampleMovie::SampleMovie(std::size_t systemId, std::size_t sampleEvery) : sampleEvery(sampleEvery)
{
  std::filesystem::create_directory("movies");
  std::ofstream stream(std::format("movies/movie.s{}.pdb", systemId));
}

void SampleMovie::update(const ForceField &forceField, std::size_t systemId, const SimulationBox simulationBox,
                         std::span<Atom> atomData, std::size_t currentCycle)
{
  if (currentCycle % sampleEvery == 0)
  {
    std::filesystem::create_directory("movies");
    std::ofstream stream(std::format("movies/movie.s{}.pdb", systemId), std::ios_base::app);

    std::print(stream, "MODEL {:>4}\n", modelNumber);
    std::print(stream, "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n", simulationBox.lengthA,
               simulationBox.lengthB, simulationBox.lengthC, simulationBox.angleAlpha * 180.0 / std::numbers::pi,
               simulationBox.angleBeta * 180.0 / std::numbers::pi, simulationBox.angleGamma * 180.0 / std::numbers::pi);

    for (int index = 1; const Atom &atom : atomData)
    {
      std::size_t atomicNumber = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].atomicNumber;
      std::string name = std::format("{:<4}", forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].name);
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const SampleMovie &m)
{
  archive << m.versionNumber;

  archive << m.sampleEvery;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, SampleMovie &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'SampleMovie' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> m.sampleEvery;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("SampleMovie: Error in binary restart\n"));
  }
#endif

  return archive;
}
