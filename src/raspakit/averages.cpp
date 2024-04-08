module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <vector>
#include <complex>
#endif

module averages;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <istream>;
import <ostream>;
import <fstream>;
import <vector>;
import <complex>;
#endif

import archive;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BlockErrorEstimation &blockerror)
{
  archive << blockerror.numberOfBins;
  archive << blockerror.currentSample;
  archive << blockerror.numberOfSamples;
  archive << blockerror.currentBin;
  archive << blockerror.binSize;
  archive << blockerror.nextBin;

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BlockErrorEstimation &blockerror)
{
  archive >> blockerror.numberOfBins;
  archive >> blockerror.currentSample;
  archive >> blockerror.numberOfSamples;
  archive >> blockerror.currentBin;
  archive >> blockerror.binSize;
  archive >> blockerror.nextBin;

  return archive;
}
