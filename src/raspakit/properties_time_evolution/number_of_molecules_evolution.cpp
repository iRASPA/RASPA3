module;

module property_number_of_molecules_evolution;

import std;

import archive;
import units;
import component;

void PropertyNumberOfMoleculesEvolution::addSample(std::size_t absoluteCurrentCycle,
                                                   const std::vector<std::size_t> &numberOfIntegerMoleculesPerComponent)
{
  if(absoluteCurrentCycle % sampleEvery == 0uz)
  {
    for(std::size_t i = 0; i != numberOfIntegerMoleculesPerComponent.size(); ++i)
    {
      std::size_t index = absoluteCurrentCycle / sampleEvery;
      if(index < result[i].size())
      {
        result[i][index] = numberOfIntegerMoleculesPerComponent[i];
      }
    }
  }
}

void PropertyNumberOfMoleculesEvolution::writeOutput(std::size_t systemId, std::size_t absoluteCurrentCycle)
{
  if(!writeEvery.has_value()) return;

  if (absoluteCurrentCycle % writeEvery.value() != 0uz) return;

  std::filesystem::create_directory("number_of_molecules_evolution");

  std::ofstream stream_output(std::format("number_of_molecules_evolution/number_of_molecules_evolution.s{}.txt", systemId));

  std::size_t currentIndex = absoluteCurrentCycle / sampleEvery;
  for(std::size_t index = 0; index < std::min(currentIndex, result[0].size()); ++index)
  {
    for(std::size_t i = 0; i < numberOfComponents; ++i)
    {
      stream_output << std::format("{} ", result[i][index]);
    }
    stream_output << "\n";
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyNumberOfMoleculesEvolution &evolution)
{
  archive << evolution.versionNumber;

  archive << evolution.numberOfComponents;
  archive << evolution.sampleEvery;
  archive << evolution.writeEvery;
  archive << evolution.totalSamples;
  archive << evolution.result;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesEvolution &evolution)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > evolution.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyNumberOfMoleculesEvolution' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> evolution.numberOfComponents;
  archive >> evolution.sampleEvery;
  archive >> evolution.writeEvery;
  archive >> evolution.totalSamples;
  archive >> evolution.result;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyNumberOfMoleculesEvolution: Error in binary restart\n"));
  }
#endif

  return archive;
}
