module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <istream>
#include <map>
#include <ostream>
#include <source_location>
#include <utility>
#include <vector>
#endif

module averages;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BlockErrorEstimation &blockerror)
{
  archive << blockerror.versionNumber;

  archive << blockerror.numberOfBins;
  archive << blockerror.currentSample;
  archive << blockerror.numberOfSamples;
  archive << blockerror.currentBin;
  archive << blockerror.binSize;
  archive << blockerror.nextBin;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BlockErrorEstimation &blockerror)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > blockerror.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BlockErrorEstimation' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> blockerror.numberOfBins;
  archive >> blockerror.currentSample;
  archive >> blockerror.numberOfSamples;
  archive >> blockerror.currentBin;
  archive >> blockerror.binSize;
  archive >> blockerror.nextBin;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("BlockErrorEstimation: Error in binary restart\n"));
  }
#endif

  return archive;
}
