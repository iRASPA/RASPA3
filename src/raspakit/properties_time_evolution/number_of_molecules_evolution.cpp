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
      result.at(i).at(absoluteCurrentCycle / sampleEvery) = numberOfIntegerMoleculesPerComponent[i];
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyNumberOfMoleculesEvolution &evolution)
{
  archive << evolution.versionNumber;

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
