module;

module property_number_of_molecules_histogram;

import std;

import archive;
import units;
import component;

void PropertyNumberOfMoleculesHistogram::addSample(std::size_t blockIndex, std::size_t currentCycle,
                                                   std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                                   const double &weight)
{
  std::make_signed_t<std::size_t> bin;

  if (currentCycle % sampleEvery != 0uz) return;

  for (std::size_t i = 0; i < numberOfIntegerMoleculesPerComponent.size(); ++i)
  {
    bin = std::make_signed_t<std::size_t>(numberOfIntegerMoleculesPerComponent[i]) - std::make_signed_t<std::size_t>(range.first);
    if (bin > 0 && bin < std::make_signed_t<std::size_t>(numberOfBins))
    {
      histogram(blockIndex, i, static_cast<std::size_t>(bin)) += weight;
    }
  }

  histogram.addCount(blockIndex, weight);
}

std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>
PropertyNumberOfMoleculesHistogram::result() const
{
  std::size_t numberOfComponents = histogram.numberOfChannels;
  std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>  
    result_data(numberOfComponents,{std::vector<double>(numberOfBins), std::vector<double>(numberOfBins), 
                                    std::vector<double>(numberOfBins)});

  std::vector<double> bin_data(numberOfBins);
  for (std::size_t i = 0; i != numberOfBins; ++i)
  {
    bin_data[i] = static_cast<double>(i) + static_cast<double>(range.first);
  }

  for(std::size_t i = 0; i < numberOfComponents; ++i)
  {
    auto [average, confidenceIntervalError] = histogram.average(i);

    result_data[i] = {bin_data, average, confidenceIntervalError};
  }

  return result_data;
}

void PropertyNumberOfMoleculesHistogram::writeOutput(std::size_t systemId, std::vector<Component> &components,
                                                     std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  std::filesystem::create_directory("number_of_molecules_histogram");

  std::ofstream stream_output(
      std::format("number_of_molecules_histogram/number_of_molecules_histogram_histogram.s{}.txt", systemId));

  std::size_t numberOfComponents = histogram.numberOfChannels;

  stream_output << std::format("# number_of_molecules_histogram, number of counts: {}\n",
                               histogram.totalNumberOfCounts);
  stream_output << "# column 1: number of molecules [-]\n";
  for (std::size_t i = 0; i < numberOfComponents; ++i)
  {
    stream_output << std::format("# column {}: number of molecules component {} [-]\n", 2 * i + 2, components[i].name);
    stream_output << std::format("# column {}: number of molecules component {} error [-]\n", 2 * i + 3,
                                 components[i].name);
  }

  std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> results_data = result();

  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    for (std::size_t i = 0; i < numberOfComponents; ++i)
    {
      stream_output << std::format("{}", bin + range.first);

      auto [bins, average, error] = results_data[i];

      stream_output << std::format(" {} {}", average[bin], error[bin]);
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
  archive << hist.sampleEvery;
  archive << hist.writeEvery;
  archive << hist.histogram;

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
    throw std::runtime_error(std::format("Invalid version reading 'PropertyNumberOfMoleculesHistogram' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> hist.numberOfBlocks;
  archive >> hist.numberOfBins;
  archive >> hist.range;
  archive >> hist.sampleEvery;
  archive >> hist.writeEvery;
  archive >> hist.histogram;

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
