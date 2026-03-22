module;

module chiral_center;

import std;

import archive;
import randomnumbers;
import double3;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ChiralCenter &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.ids;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ChiralCenter &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ChiralCenter' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.ids;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("ChiralCenter: Error in binary restart\n"));
  }
#endif

  return archive;
}
