module;

module widom_data;

import std;

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const WidomData &l)
{
  archive << l.versionNumber;

  archive << l.total;
  archive << l.excess;
  archive << l.idealGas;

#if DEBUG_ARCHIVE
  archive << static_caststd::<int64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, WidomData &l)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > l.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'WidomData' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> l.total;
  archive >> l.excess;
  archive >> l.idealGas;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("WidomData: Error in binary restart\n"));
  }
#endif

  return archive;
}
