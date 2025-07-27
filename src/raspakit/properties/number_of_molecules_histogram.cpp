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
#include <iostream>
#include <map>
#include <print>
#include <ranges>
#include <source_location>
#include <vector>
#endif

module property_number_of_molecules_histogram;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import units;
import component;

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<double>());
  return result;
}

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<double>());
  return result;
}

std::vector<double> operator*(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::multiplies<double>());
  return result;
}

std::vector<double> operator*(const double &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(b.size());

  std::transform(b.begin(), b.end(), std::back_inserter(result), [&a](double v) { return a * v; });

  return result;
}

std::vector<double> operator/(const std::vector<double> &a, const double &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), std::back_inserter(result), [&b](double v) { return v / b; });
  return result;
}

std::vector<double> sqrt(const std::vector<double> &a)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), std::back_inserter(result), [](double v) { return std::sqrt(v); });
  return result;
}

void PropertyNumberOfMoleculesHistogram::addSample(std::size_t blockIndex, std::size_t currentCycle,
                                                   std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                                   const double &weight)
{
  std::size_t bin;

  if (currentCycle % sampleEvery != 0uz) return;

  for (std::size_t i = 0; i < numberOfIntegerMoleculesPerComponent.size(); ++i)
  {
    bin = numberOfIntegerMoleculesPerComponent[i] - range.first;
    if (bin >= 0 && bin < numberOfBins)
    {
      bookKeepingEnergyHistogram[blockIndex][bin][i] += weight;
    }
  }

  numberOfCounts[blockIndex] += weight;
  totalNumberOfCounts += weight;
}

std::vector<std::vector<double>> PropertyNumberOfMoleculesHistogram::averagedProbabilityHistogram(
    std::size_t blockIndex) const
{
  std::vector<std::vector<double>> averagedData(numberOfBins, std::vector<double>(size));
  std::transform(bookKeepingEnergyHistogram[blockIndex].begin(), bookKeepingEnergyHistogram[blockIndex].end(),
                 averagedData.begin(),
                 [&](const std::vector<double> &sample) { return sample / numberOfCounts[blockIndex]; });
  return averagedData;
}

std::vector<std::vector<double>> PropertyNumberOfMoleculesHistogram::averagedProbabilityHistogram() const
{
  std::vector<std::vector<double>> summedBlocks(numberOfBins, std::vector<double>(size));
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingEnergyHistogram[blockIndex].begin(),
                   summedBlocks.begin(),
                   [](const std::vector<double> &a, const std::vector<double> &b) { return a + b; });
  }

  std::vector<std::vector<double>> average(numberOfBins, std::vector<double>(size));
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const std::vector<double> &sample) { return sample / totalNumberOfCounts; });

  return average;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
PropertyNumberOfMoleculesHistogram::averageProbabilityHistogram() const
{
  std::vector<std::vector<double>> average = averagedProbabilityHistogram();

  std::vector<std::vector<double>> sumOfSquares(numberOfBins, std::vector<double>(size));
  std::size_t numberOfSamples = 0;
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<std::vector<double>> blockAverage = averagedProbabilityHistogram(blockIndex);

    if (numberOfCounts[blockIndex] > 0.0)
    {
      for (std::size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        std::vector<double> value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] = sumOfSquares[binIndex] + value * value;
      }
      ++numberOfSamples;
    }
  }
  std::vector<std::vector<double>> confidenceIntervalError(numberOfBins, std::vector<double>(size));
  if (numberOfSamples >= 3)
  {
    std::size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<std::vector<double>> standardDeviation(numberOfBins, std::vector<double>(size));
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const std::vector<double> &sumofsquares)
                   { return sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

    std::vector<std::vector<double>> standardError(numberOfBins, std::vector<double>(size));
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const std::vector<double> &sigma)
                   { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const std::vector<double> &error) { return intermediateStandardNormalDeviate * error; });
  }

  return std::make_pair(average, confidenceIntervalError);
}

void PropertyNumberOfMoleculesHistogram::writeOutput(std::size_t systemId, std::vector<Component> &components,
                                                     std::size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("number_of_molecules_histogram");

  std::ofstream stream_output(
      std::format("number_of_molecules_histogram/number_of_molecules_histogram_histogram.s{}.txt", systemId));

  stream_output << std::format("# number_of_molecules_histogram, number of counts: {}\n", totalNumberOfCounts);
  stream_output << "# column 1: number of molecules [-]\n";
  for (std::size_t i = 0; i < size; ++i)
  {
    stream_output << std::format("# column {}: number of molecules component {} [-]\n", 2 * i + 2, components[i].name);
    stream_output << std::format("# column {}: number of molecules component {} error [-]\n", 2 * i + 3,
                                 components[i].name);
  }

  auto [average, error] = averageProbabilityHistogram();

  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    stream_output << std::format("{}", bin + range.first);

    for (std::size_t i = 0; i < size; ++i)
    {
      stream_output << std::format(" {} {}", average[bin][i], error[bin][i]);
    }
    stream_output << "\n";
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyNumberOfMoleculesHistogram &hist)
{
  archive << hist.versionNumber;

  archive << hist.numberOfBlocks;
  archive << hist.numberOfBins;
  archive << hist.range;
  archive << hist.size;
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesHistogram &hist)
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
  archive >> hist.size;
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
    throw std::runtime_error(std::format("PropertyNumberOfMoleculesHistogram: Error in binary restart\n"));
  }
#endif

  return archive;
}
