module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <ranges>
#include <source_location>
#include <sstream>
#include <vector>
#endif

module property_energy_histogram;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;
import archive;
import units;

void PropertyEnergyHistogram::addSample(std::size_t blockIndex, std::size_t currentCycle, double4 energy,
                                        const double &weight)
{
  std::size_t bin;

  if (currentCycle % sampleEvery != 0uz) return;

  bin = static_cast<std::size_t>((energy.x * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin >= 0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].x += weight;
  }
  bin = static_cast<std::size_t>((energy.y * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin >= 0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].y += weight;
  }
  bin = static_cast<std::size_t>((energy.z * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin >= 0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].z += weight;
  }
  bin = static_cast<std::size_t>((energy.w * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin >= 0 && bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].w += weight;
  }

  numberOfCounts[blockIndex] += weight;
  totalNumberOfCounts += weight;
}

std::vector<double4> PropertyEnergyHistogram::averagedProbabilityHistogram(std::size_t blockIndex) const
{
  std::vector<double4> averagedData(numberOfBins);
  std::transform(bookKeepingEnergyHistogram[blockIndex].begin(), bookKeepingEnergyHistogram[blockIndex].end(),
                 averagedData.begin(), [&](const double4 &sample) { return sample / numberOfCounts[blockIndex]; });
  return averagedData;
}

std::vector<double4> PropertyEnergyHistogram::averagedProbabilityHistogram() const
{
  std::vector<double4> summedBlocks(numberOfBins);
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingEnergyHistogram[blockIndex].begin(),
                   summedBlocks.begin(), [](const double4 &a, const double4 &b) { return a + b; });
  }
  std::vector<double4> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double4 &sample) { return sample / totalNumberOfCounts; });

  return average;
}

std::pair<std::vector<double4>, std::vector<double4>> PropertyEnergyHistogram::averageProbabilityHistogram() const
{
  std::size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  std::vector<double4> average = averagedProbabilityHistogram();

  std::vector<double4> sumOfSquares(numberOfBins);
  std::size_t numberOfSamples = 0;
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<double4> blockAverage = averagedProbabilityHistogram(blockIndex);

    if (numberOfCounts[blockIndex] > 0.0)
    {
      for (std::size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        double4 value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
      ++numberOfSamples;
    }
  }
  std::vector<double4> confidenceIntervalError(numberOfBins);
  if (numberOfSamples >= 3)
  {
    std::vector<double4> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const double4 &sumofsquares)
                   { return sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

    std::vector<double4> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const double4 &sigma) { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const double4 &error) { return intermediateStandardNormalDeviate * error; });
  }

  return std::make_pair(average, confidenceIntervalError);
}

void PropertyEnergyHistogram::writeOutput(std::size_t systemId, std::size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("energy_histogram");

  std::ofstream stream_output(std::format("energy_histogram/energy_histogram.s{}.txt", systemId));

  stream_output << std::format("# energy_histogram, number of counts: {}\n", totalNumberOfCounts);
  stream_output << "# column 1: energy [K]\n";
  stream_output << "# column 2: total energy histogram [-]\n";
  stream_output << "# column 3: total energy histogram error [-]\n";
  stream_output << "# column 4: VDW energy histogram [-]\n";
  stream_output << "# column 5: VDW energy histogram error [-]\n";
  stream_output << "# column 6: Coulombic energy histogram [-]\n";
  stream_output << "# column 7: Coulombic energy histogram error [-]\n";
  stream_output << "# column 8: Polarization energy histogram [-]\n";
  stream_output << "# column 9: Polarization energy histogram error [-]\n";

  auto [average, error] = averageProbabilityHistogram();

  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    double energy =
        static_cast<double>(bin) * std::fabs(range.second - range.first) / static_cast<double>(numberOfBins) +
        range.first;
    if (average[bin].x > 0.0)
    {
      stream_output << std::format("{} {} {} {} {} {} {} {} {}\n", energy, average[bin].x, error[bin].x, average[bin].y,
                                   error[bin].y, average[bin].z, error[bin].z, average[bin].w, error[bin].w);
    }
  }
}

std::string PropertyEnergyHistogram::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "Energy histogram:\n");
  std::print(stream, "    sample every: {}\n", sampleEvery);
  std::print(stream, "    write every: {}\n", writeEvery);
  std::print(stream, "    number of bins: {}\n", numberOfBins);
  std::print(stream, "    range: ({}) - ({})\n", range.first, range.second);
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist)
{
  archive << hist.versionNumber;

  archive << hist.numberOfBlocks;
  archive << hist.numberOfBins;
  archive << hist.range;
  archive << hist.sampleEvery;
  archive << hist.writeEvery;
  archive << hist.bookKeepingEnergyHistogram;
  archive << hist.numberOfCounts;
  archive << hist.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &hist)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > hist.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyHistogram' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> hist.numberOfBlocks;
  archive >> hist.numberOfBins;
  archive >> hist.range;
  archive >> hist.sampleEvery;
  archive >> hist.writeEvery;
  archive >> hist.bookKeepingEnergyHistogram;
  archive >> hist.numberOfCounts;
  archive >> hist.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyEnergyHistogram: Error in binary restart\n"));
  }
#endif

  return archive;
}
