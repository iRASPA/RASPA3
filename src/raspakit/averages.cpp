module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <fstream>
#include <iostream>
#include <istream>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#endif

module averages;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <istream>;
import <ostream>;
import <fstream>;
import <complex>;
import <vector>;
import <map>;
import <array>;
import <utility>;
import <algorithm>;
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
