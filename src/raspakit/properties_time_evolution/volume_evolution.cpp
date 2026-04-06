module;

module property_volume_evolution;

import std;

import archive;
import units;
import component;

void PropertyVolumeEvolution::addSample(std::size_t absoluteCurrentCycle, double currentVolume)
{
  if(absoluteCurrentCycle % sampleEvery == 0uz)
  {
    result.at(absoluteCurrentCycle / sampleEvery) = currentVolume;
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyVolumeEvolution &evolution)
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

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyVolumeEvolution &evolution)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > evolution.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyVolumeEvolution' at line {} in file {}\n",
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
    throw std::runtime_error(std::format("PropertyVolumeEvolution: Error in binary restart\n"));
  }
#endif

  return archive;
}
