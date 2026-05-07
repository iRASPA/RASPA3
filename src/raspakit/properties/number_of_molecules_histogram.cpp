module;

module property_number_of_molecules_histogram;

import std;

import archive;
import units;
import component;

inline std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<double>());
  return result;
}

inline std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<double>());
  return result;
}

inline std::vector<double> operator*(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::multiplies<double>());
  return result;
}

inline std::vector<double> operator*(const double &a, const std::vector<double> &b)
{
  std::vector<double> result;
  result.reserve(b.size());

  std::transform(b.begin(), b.end(), std::back_inserter(result), [&a](double v) { return a * v; });

  return result;
}

inline std::vector<double> operator/(const std::vector<double> &a, const double &b)
{
  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), std::back_inserter(result), [&b](double v) { return v / b; });
  return result;
}

inline std::vector<double> sqrt(const std::vector<double> &a)
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
  std::make_signed_t<std::size_t> bin;

  if (currentCycle % sampleEvery != 0uz) return;

  for (std::size_t i = 0; i < numberOfIntegerMoleculesPerComponent.size(); ++i)
  {
    bin = std::make_signed_t<std::size_t>(numberOfIntegerMoleculesPerComponent[i]) - std::make_signed_t<std::size_t>(range.first);
    if (bin > 0 && bin < std::make_signed_t<std::size_t>(numberOfBins))
    {
      bookKeepingEnergyHistogram[blockIndex][i][static_cast<std::size_t>(bin)] += weight;
    }
  }

  numberOfCounts[blockIndex] += weight;
  totalNumberOfCounts += weight;
}

std::vector<double> PropertyNumberOfMoleculesHistogram::averagedProbabilityHistogram(std::size_t blockIndex, std::size_t component_id) const
{
  std::vector<double> averagedData(numberOfBins);

  // loop over components, and operate on the vector of bins
  std::transform(bookKeepingEnergyHistogram[blockIndex][component_id].begin(), bookKeepingEnergyHistogram[blockIndex][component_id].end(),
                 averagedData.begin(),
                 [&](const double &sample) { return sample / numberOfCounts[blockIndex]; });
  return averagedData;
}

std::vector<double> PropertyNumberOfMoleculesHistogram::averagedProbabilityHistogram( std::size_t component_id) const
{
  std::vector<double> summedBlocks(numberOfBins);
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingEnergyHistogram[blockIndex][component_id].begin(),
                   summedBlocks.begin(),
                   [](const double &a, const double &b) { return a + b; });
  }

  std::vector<double> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double &sample) { return sample / totalNumberOfCounts; });

  return average;
}

std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>
PropertyNumberOfMoleculesHistogram::result() const
{
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
    std::vector<double> average = averagedProbabilityHistogram(i);

    std::vector<double> sumOfSquares(numberOfBins);

    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex, i);

      if (numberOfCounts[blockIndex] > 0.0)
      {
        for (std::size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
        {
          double value = blockAverage[binIndex] - average[binIndex];
          sumOfSquares[binIndex] = sumOfSquares[binIndex] + value * value;
        }
        ++numberOfSamples;
      }
    }

    std::vector<double> confidenceIntervalError(numberOfBins);
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfBlocks - 1;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      std::vector<double> standardDeviation(numberOfBins);
      std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                     [&](const double &sumofsquares)
                     { return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

      std::vector<double> standardError(numberOfBins);
      std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                     [&](const double &sigma)
                     { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

      std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                     [&](const double &error) { return intermediateStandardNormalDeviate * error; });
    }

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

  stream_output << std::format("# number_of_molecules_histogram, number of counts: {}\n", totalNumberOfCounts);
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
  archive << hist.numberOfComponents;
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
  archive >> hist.numberOfComponents;
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
    throw std::runtime_error(std::format("PropertyNumberOfMoleculesHistogram: Error in binary restart\n"));
  }
#endif

  return archive;
}
