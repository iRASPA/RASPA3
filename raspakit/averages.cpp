module;

module averages;

import <iostream>;
import <istream>;
import <ostream>;
import <fstream>;
import <vector>;
import <complex>;

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
