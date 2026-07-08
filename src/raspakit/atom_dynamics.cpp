module;

module atom_dynamics;

import std;

import archive;
import double3;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const AtomDynamics &atom)
{
  archive << atom.velocity;
  archive << atom.gradient;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, AtomDynamics &atom)
{
  archive >> atom.velocity;
  archive >> atom.gradient;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("AtomDynamics: Error in binary restart\n"));
  }
#endif

  return archive;
}
