module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#endif

module property_msd;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <string>;
import <iostream>;
import <fstream>;
import <sstream>;
import <tuple>;
import <vector>;
import <algorithm>;
import <complex>;
import <format>;
import <numbers>;
import <span>;
import <array>;
import <cmath>;
import <exception>;
import <source_location>;
import <print>;
import <filesystem>;
#endif

import archive;
import double3;
import double4;
import atom;
import simulationbox;
import forcefield;
import component;
import averages;

void PropertyMeanSquaredDisplacement::addSample(size_t currentCycle, const std::vector<Component> &components,
                                                const std::vector<size_t> &numberOfMoleculesPerComponent,
                                                std::vector<Molecule> &moleculePositions)
{
  if (currentCycle % sampleEvery != 0uz) return;

  // determine current number of blocks
  numberOfBlocksMSD = 1;
  size_t p = countMSD / numberOfBlockElementsMSD;
  while (p != 0)
  {
    ++numberOfBlocksMSD;
    p /= numberOfBlockElementsMSD;
  }

  if (numberOfBlocksMSD > maxNumberOfBlocksMSD)
  {
    blockLengthMSD.resize(numberOfBlocksMSD);
    maxNumberOfBlocksMSD = numberOfBlocksMSD;

    msdSelfCount.resize(numberOfBlocksMSD, std::vector<std::vector<size_t>>(
                                               numberOfComponents, std::vector<size_t>(numberOfBlockElementsMSD, 0uz)));
    blockDataMSDSelf.resize(numberOfBlocksMSD,
                            std::vector<std::vector<double3>>(
                                numberOfParticles, std::vector<double3>(numberOfBlockElementsMSD, double3())));
    msdSelf.resize(numberOfBlocksMSD,
                   std::vector<std::vector<double4>>(numberOfComponents,
                                                     std::vector<double4>(numberOfBlockElementsMSD, double4())));

    msdOnsagerCount.resize(
        numberOfBlocksMSD,
        std::vector<std::vector<size_t>>(numberOfComponents, std::vector<size_t>(numberOfBlockElementsMSD, 0uz)));
    blockDataMSDOnsager.resize(numberOfBlocksMSD,
                               std::vector<std::vector<double3>>(
                                   numberOfComponents, std::vector<double3>(numberOfBlockElementsMSD, double3())));
    msdOnsager.resize(
        numberOfBlocksMSD,
        std::vector<std::vector<std::vector<double4>>>(
            numberOfComponents, std::vector<std::vector<double4>>(
                                    numberOfComponents, std::vector<double4>(numberOfBlockElementsMSD, double4()))));
  }

  for (size_t currentBlock = 0; currentBlock < numberOfBlocksMSD; currentBlock++)
  {
    // test for blocking operation: CountMSD is a multiple of NumberOfBlockElementsMSD^CurrentBlock
    if (countMSD % static_cast<size_t>(std::pow(numberOfBlockElementsMSD, currentBlock)) == 0)
    {
      // increase the current block-length
      blockLengthMSD[currentBlock]++;

      // limit length to numberOfBlockElementsMSD
      size_t currentBlocklength = std::min(blockLengthMSD[currentBlock], numberOfBlockElementsMSD);

      size_t molecule_index{0};
      for (size_t i = 0; i != components.size(); ++i)
      {
        // self diffusion
        for (size_t m = 0; m != numberOfMoleculesPerComponent[i]; ++m)
        {
          double3 value = moleculePositions[molecule_index].centerOfMassPosition;

          std::shift_right(begin(blockDataMSDSelf[currentBlock][molecule_index]),
                           end(blockDataMSDSelf[currentBlock][molecule_index]), 1);

          blockDataMSDSelf[currentBlock][molecule_index][0] = value;

          for (size_t k = 0; k < currentBlocklength; ++k)
          {
            // msd for each component
            ++msdSelfCount[currentBlock][i][k];

            double dr_x = blockDataMSDSelf[currentBlock][molecule_index][k].x - value.x;
            double dr_y = blockDataMSDSelf[currentBlock][molecule_index][k].y - value.y;
            double dr_z = blockDataMSDSelf[currentBlock][molecule_index][k].z - value.z;

            double msd_x = dr_x * dr_x;
            double msd_y = dr_y * dr_y;
            double msd_z = dr_z * dr_z;

            msdSelf[currentBlock][i][k].x += msd_x;
            msdSelf[currentBlock][i][k].y += msd_y;
            msdSelf[currentBlock][i][k].z += msd_z;
            msdSelf[currentBlock][i][k].w += msd_x + msd_y + msd_z;
          }

          ++molecule_index;
        }
      }

      molecule_index = 0;
      std::vector<double3> value_onsager(numberOfComponents);
      for (size_t i = 0; i != components.size(); ++i)
      {
        // self diffusion
        for (size_t m = 0; m != numberOfMoleculesPerComponent[i]; ++m)
        {
          double3 value = moleculePositions[molecule_index].centerOfMassPosition;
          value_onsager[i] += value;
          ++molecule_index;
        }
      }

      for (size_t i = 0; i != components.size(); ++i)
      {
        std::shift_right(begin(blockDataMSDOnsager[currentBlock][i]), end(blockDataMSDOnsager[currentBlock][i]), 1);
        blockDataMSDOnsager[currentBlock][i][0] = value_onsager[i];
      }

      for (size_t k = 0; k < currentBlocklength; k++)
      {
        for (size_t i = 0; i != components.size(); ++i)
        {
          ++msdOnsagerCount[currentBlock][i][k];
          for (size_t j = 0; j != components.size(); ++j)
          {
            double msd_x = (blockDataMSDOnsager[currentBlock][i][k].x - value_onsager[i].x) *
                           (blockDataMSDOnsager[currentBlock][j][k].x - value_onsager[j].x);
            double msd_y = (blockDataMSDOnsager[currentBlock][i][k].y - value_onsager[i].y) *
                           (blockDataMSDOnsager[currentBlock][j][k].y - value_onsager[j].y);
            double msd_z = (blockDataMSDOnsager[currentBlock][i][k].z - value_onsager[i].z) *
                           (blockDataMSDOnsager[currentBlock][j][k].z - value_onsager[j].z);

            msdOnsager[currentBlock][i][j][k].x += msd_x;
            msdOnsager[currentBlock][i][j][k].y += msd_y;
            msdOnsager[currentBlock][i][j][k].z += msd_z;
            msdOnsager[currentBlock][i][j][k].w += msd_x + msd_y + msd_z;
          }
        }
      }
    }
  }

  ++countMSD;
}

void PropertyMeanSquaredDisplacement::writeOutput(size_t systemId, const std::vector<Component> &components,
                                                  const std::vector<size_t> &numberOfMoleculesPerComponent,
                                                  double deltaT, size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("msd");

  for (size_t i = 0; i < components.size(); ++i)
  {
    std::ofstream stream_msd_self_output(std::format("msd/msd_self_{}.s{}.txt", components[i].name, systemId));

    stream_msd_self_output << std::format("# msd, number of counts: {}\n", countMSD);
    stream_msd_self_output << "# column 1: time [ps]\n";
    stream_msd_self_output << "# column 2: msd xyz [A^2]\n";
    stream_msd_self_output << "# column 3: msd x [A^2]\n";
    stream_msd_self_output << "# column 4: msd y [A^2]\n";
    stream_msd_self_output << "# column 5: msd z [A^2]\n";
    stream_msd_self_output << "# column 6: number of samples [-]\n";

    for (size_t currentBlock = 0; currentBlock < numberOfBlocksMSD; ++currentBlock)
    {
      size_t currentBlocklength = std::min(blockLengthMSD[currentBlock], numberOfBlockElementsMSD);
      double dt = static_cast<double>(sampleEvery) * deltaT * std::pow(numberOfBlockElementsMSD, currentBlock);
      for (size_t k = 1; k < currentBlocklength; ++k)
      {
        if (msdSelfCount[currentBlock][i][k] > 0)
        {
          double fac = 1.0 / static_cast<double>(msdSelfCount[currentBlock][i][k]);

          stream_msd_self_output << std::format(
              "{} {} {} {} {} (count: {})\n", static_cast<double>(k) * dt, fac * msdSelf[currentBlock][i][k].w,
              fac * msdSelf[currentBlock][i][k].x, fac * msdSelf[currentBlock][i][k].y,
              fac * msdSelf[currentBlock][i][k].z, msdSelfCount[currentBlock][i][k]);
        }
      }
    }
  }

  for (size_t i = 0; i < components.size(); ++i)
  {
    for (size_t j = 0; j < components.size(); ++j)
    {
      std::ofstream stream_msd_collective_output(
          std::format("msd/msd_onsager_{}_{}.s{}.txt", components[i].name, components[j].name, systemId));

      stream_msd_collective_output << std::format("# msd, number of counts: {}\n", countMSD);
      stream_msd_collective_output << "# column 1: time [ps]\n";
      stream_msd_collective_output << "# column 2: msd xyz [A^2]\n";
      stream_msd_collective_output << "# column 3: msd x [A^2]\n";
      stream_msd_collective_output << "# column 4: msd y [A^2]\n";
      stream_msd_collective_output << "# column 5: msd z [A^2]\n";
      stream_msd_collective_output << "# column 6: number of samples [-]\n";

      for (size_t currentBlock = 0; currentBlock < numberOfBlocksMSD; ++currentBlock)
      {
        size_t currentBlocklength = std::min(blockLengthMSD[currentBlock], numberOfBlockElementsMSD);
        double dt = static_cast<double>(sampleEvery) * deltaT * std::pow(numberOfBlockElementsMSD, currentBlock);
        for (size_t k = 1; k < currentBlocklength; ++k)
        {
          if (msdOnsagerCount[currentBlock][i][k] > 0)
          {
            double fac =
                1.0 / static_cast<double>(numberOfMoleculesPerComponent[i] * msdOnsagerCount[currentBlock][i][k]);

            stream_msd_collective_output << std::format(
                "{} {} {} {} {} (count: {})\n", static_cast<double>(k) * dt, fac * msdOnsager[currentBlock][i][j][k].w,
                fac * msdOnsager[currentBlock][i][j][k].x, fac * msdOnsager[currentBlock][i][j][k].y,
                fac * msdOnsager[currentBlock][i][j][k].z, msdOnsagerCount[currentBlock][i][k]);
          }
        }
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyMeanSquaredDisplacement &msd)
{
  archive << msd.versionNumber;

  archive << msd.numberOfBlocks;
  archive << msd.sampleEvery;
  archive << msd.writeEvery;
  archive << msd.countMSD;
  archive << msd.numberOfBlocksMSD;
  archive << msd.maxNumberOfBlocksMSD;
  archive << msd.numberOfBlockElementsMSD;
  archive << msd.blockLengthMSD;
  archive << msd.msdSelfCount;
  archive << msd.blockDataMSDSelf;
  archive << msd.msdSelf;
  archive << msd.msdOnsagerCount;
  archive << msd.blockDataMSDOnsager;
  archive << msd.msdOnsager;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMeanSquaredDisplacement &msd)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > msd.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyRadialDistributionFunction' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> msd.numberOfBlocks;
  archive >> msd.sampleEvery;
  archive >> msd.writeEvery;
  archive >> msd.countMSD;
  archive >> msd.numberOfBlocksMSD;
  archive >> msd.maxNumberOfBlocksMSD;
  archive >> msd.numberOfBlockElementsMSD;
  archive >> msd.blockLengthMSD;
  archive >> msd.msdSelfCount;
  archive >> msd.blockDataMSDSelf;
  archive >> msd.msdSelf;
  archive >> msd.msdOnsagerCount;
  archive >> msd.blockDataMSDOnsager;
  archive >> msd.msdOnsager;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyMeanSquaredDisplacement: Error in binary restart\n"));
  }
#endif

  return archive;
}
