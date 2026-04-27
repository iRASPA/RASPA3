module;

module property_conserved_energy_evolution;

import std;

import archive;
import units;
import running_energy;

void PropertyConservedEnergyEvolution::addSample(std::size_t absoluteCurrentCycle,
                                         const RunningEnergy &running_energy)
{
  if(absoluteCurrentCycle % sampleEvery == 0uz)
  {
    std::size_t index = absoluteCurrentCycle / sampleEvery;
    if(index < data.size())
    {
      data[index] = running_energy;
    }
  }
}

void PropertyConservedEnergyEvolution::writeOutput(std::size_t systemId, std::size_t absoluteCurrentCycle)
{
  if(!writeEvery.has_value()) return;

  if (absoluteCurrentCycle % writeEvery.value() != 0uz) return;

  std::filesystem::create_directory("conserved_energy");

  std::ofstream stream_output(std::format("conserved_energy/conserved_energy.s{}.txt", systemId));

  std::size_t currentIndex = absoluteCurrentCycle / sampleEvery;
  for(std::size_t index = 0; index < std::min(currentIndex, data.size()); ++index)
  {
    stream_output << std::format("{} {} {} {} {}\n", index, data[index].conservedEnergy(),
       data[index].potentialEnergy(), data[index].kineticEnergy(), data[index].NoseHooverEnergy);
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyConservedEnergyEvolution &evolution)
{
  archive << evolution.versionNumber;

  archive << evolution.sampleEvery;
  archive << evolution.writeEvery;
  archive << evolution.totalSamples;
  archive << evolution.data;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyConservedEnergyEvolution &evolution)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > evolution.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyConservedEnergyEvolution' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> evolution.sampleEvery;
  archive >> evolution.writeEvery;
  archive >> evolution.totalSamples;
  archive >> evolution.data;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyConservedEnergyEvolution: Error in binary restart\n"));
  }
#endif

  return archive;
}
