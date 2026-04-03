module;

module property_energy_histogram;

import std;

import archive;
import units;
import average_energy_type;

void PropertyEnergyHistogram::addSample(std::size_t blockIndex, std::size_t currentCycle, AverageEnergyType energy,
                                        const double &weight)
{
  std::size_t bin;

  if (currentCycle % sampleEvery != 0uz) return;

  bin = static_cast<std::size_t>((energy.totalEnergy * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].totalEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.VanDerWaalsEnergy * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].VanDerWaalsEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.CoulombEnergy * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].CoulombEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.polarizationEnergy * Units::EnergyToKelvin - range.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(range.second - range.first));
  if (bin < numberOfBins)
  {
    bookKeepingEnergyHistogram[blockIndex][bin].polarizationEnergy += weight;
  }

  numberOfCounts[blockIndex] += weight;
  totalNumberOfCounts += weight;
}

std::vector<AverageEnergyType> PropertyEnergyHistogram::averagedProbabilityHistogram(std::size_t blockIndex) const
{
  std::vector<AverageEnergyType> averagedData(numberOfBins);
  std::transform(bookKeepingEnergyHistogram[blockIndex].begin(), bookKeepingEnergyHistogram[blockIndex].end(),
                 averagedData.begin(), [&](const AverageEnergyType &sample) { return sample / numberOfCounts[blockIndex]; });
  return averagedData;
}

std::vector<AverageEnergyType> PropertyEnergyHistogram::averagedProbabilityHistogram() const
{
  std::vector<AverageEnergyType> summedBlocks(numberOfBins);
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingEnergyHistogram[blockIndex].begin(),
                   summedBlocks.begin(), [](const AverageEnergyType &a, const AverageEnergyType &b) { return a + b; });
  }
  std::vector<AverageEnergyType> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const AverageEnergyType &sample) { return sample / totalNumberOfCounts; });

  return average;
}

std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>> PropertyEnergyHistogram::result() const
{
  std::size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];

  std::vector<double> bins(numberOfBins);
  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    bins[bin] = static_cast<double>(bin) * std::fabs(range.second - range.first) / static_cast<double>(numberOfBins) +
        range.first;
  }


  std::vector<AverageEnergyType> average = averagedProbabilityHistogram();

  std::vector<AverageEnergyType> sumOfSquares(numberOfBins);
  std::size_t numberOfSamples = 0;
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<AverageEnergyType> blockAverage = averagedProbabilityHistogram(blockIndex);

    if (numberOfCounts[blockIndex] > 0.0)
    {
      for (std::size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        AverageEnergyType value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
      ++numberOfSamples;
    }
  }
  std::vector<AverageEnergyType> confidenceIntervalError(numberOfBins);
  if (numberOfSamples >= 3)
  {
    std::vector<AverageEnergyType> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const AverageEnergyType &sumofsquares)
                   { return sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

    std::vector<AverageEnergyType> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const AverageEnergyType &sigma) { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const AverageEnergyType &error) { return intermediateStandardNormalDeviate * error; });
  }

  return {bins, average, confidenceIntervalError};
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

  auto [energies, average, error] = result();

  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    if (average[bin].totalEnergy > 0.0)
    {
      stream_output << std::format("{} {} {} {} {} {} {} {} {}\n", energies[bin], 
          average[bin].totalEnergy, error[bin].totalEnergy, average[bin].VanDerWaalsEnergy, error[bin].VanDerWaalsEnergy,
          average[bin].CoulombEnergy, error[bin].CoulombEnergy, average[bin].polarizationEnergy, error[bin].polarizationEnergy);
    }
  }
}


std::vector<double> PropertyEnergyHistogram::bins() const
{
  std::vector<double> bins(numberOfBins);
  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    bins[bin] = static_cast<double>(bin) * std::fabs(range.second - range.first) / static_cast<double>(numberOfBins) +
        range.first;
  }
  return bins;
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
