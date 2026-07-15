module;

module energy_dudlambda;

import std;

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyDuDlambda &e)
{
  archive << e.energy;
  archive << e.dUdlambda;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyDuDlambda &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EnergyDuDlambda: Error in binary restart\n"));
  }
#endif

  return archive;
}
