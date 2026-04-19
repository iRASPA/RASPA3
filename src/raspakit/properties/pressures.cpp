module;

module pressures;

import std;

import archive;
import int3;
import stringutils;
import component;
import units;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Pressures &l)
{
  archive << l.versionNumber;

  archive << l.totalPressureTensor;
  archive << l.excessPressureTensor;
  archive << l.idealGasPressureTensor;
  archive << l.totalPressure;
  archive << l.excessPressure;
  archive << l.idealGasPressure;

#if DEBUG_ARCHIVE
  archive << static_caststd::<int64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Pressures &l)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > l.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Pressures' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> l.totalPressureTensor;
  archive >> l.excessPressureTensor;
  archive >> l.idealGasPressureTensor;
  archive >> l.totalPressure;
  archive >> l.excessPressure;
  archive >> l.idealGasPressure;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Pressures: Error in binary restart\n"));
  }
#endif

  return archive;
}
