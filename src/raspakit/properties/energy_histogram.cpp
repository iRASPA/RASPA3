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

  bin = static_cast<std::size_t>((energy.totalEnergy - valueRange.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(valueRange.second - valueRange.first));
  if (bin < numberOfBins)
  {
    histogram(blockIndex, 0, bin).totalEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.VanDerWaalsEnergy - valueRange.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(valueRange.second - valueRange.first));
  if (bin < numberOfBins)
  {
    histogram(blockIndex, 0, bin).VanDerWaalsEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.CoulombEnergy - valueRange.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(valueRange.second - valueRange.first));
  if (bin < numberOfBins)
  {
    histogram(blockIndex, 0, bin).CoulombEnergy += weight;
  }
  bin = static_cast<std::size_t>((energy.polarizationEnergy - valueRange.first) * static_cast<double>(numberOfBins) /
                                 std::fabs(valueRange.second - valueRange.first));
  if (bin < numberOfBins)
  {
    histogram(blockIndex, 0, bin).polarizationEnergy += weight;
  }

  histogram.addCount(blockIndex, weight);
}

std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>>
PropertyEnergyHistogram::result() const
{
  auto [average, confidenceIntervalError] = histogram.average(0);
  return {bins(), average, confidenceIntervalError};
}

void PropertyEnergyHistogram::writeOutput(std::size_t systemId, std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  std::filesystem::create_directory("energy_histogram");

  std::ofstream stream_output(std::format("energy_histogram/energy_histogram.s{}.txt", systemId));

  stream_output << std::format("# energy_histogram, number of counts: {}\n", histogram.totalNumberOfCounts);
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
      stream_output << std::format("{} {} {} {} {} {} {} {} {}\n", energies[bin] * Units::EnergyToKelvin, 
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
    bins[bin] = static_cast<double>(bin) * std::fabs(valueRange.second - valueRange.first) / static_cast<double>(numberOfBins) +
        valueRange.first;
  }
  return bins;
}

std::string PropertyEnergyHistogram::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "Energy histogram:\n");
  std::print(stream, "    sample every: {}\n", sampleEvery);
  if(writeEvery.has_value())
  {
    std::print(stream, "    write every: {}\n", writeEvery.value());
  }
  std::print(stream, "    number of bins: {}\n", numberOfBins);
  std::print(stream, "    valueRange: ({}) - ({})\n", valueRange.first, valueRange.second);
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist)
{
  archive << hist.versionNumber;

  archive << hist.numberOfBins;
  archive << hist.valueRange;
  archive << hist.sampleEvery;
  archive << hist.writeEvery;
  archive << hist.histogram;

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

  archive >> hist.numberOfBins;
  archive >> hist.valueRange;
  archive >> hist.sampleEvery;
  archive >> hist.writeEvery;
  archive >> hist.histogram;

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
