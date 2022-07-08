module;

module lambda;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <cmath>;
import <numbers>;

import units;
import print;
import averages;

void Lambda::WangLandauIteration(Lambda::WangLandauPhase phase)
{
  switch(phase)
  {
    case Lambda::WangLandauPhase::Initialize:
      WangLandauScalingFactor = 1.0;
      std::fill(histogram.begin(), histogram.end(), 0.0);
      std::fill(biasFactor.begin(), biasFactor.end(), 0.0);
      break;
    case Lambda::WangLandauPhase::Sample:
      biasFactor[currentBin] -= WangLandauScalingFactor;
      histogram[currentBin] += 1.0;
      break;
    case Lambda::WangLandauPhase::AdjustBiasingFactors:
      {
        std::vector<double>::iterator minValueIterator = std::min_element(histogram.begin(), histogram.end());
        double minimumValue = *minValueIterator;
        if(minimumValue>0.01)
        {
          double sumOfHistogram = std::accumulate(histogram.begin(), histogram.end(), 0.0);
          if(minimumValue / sumOfHistogram > 0.01)
          {
             WangLandauScalingFactor *= 0.5;
          }
        }
        std::fill(histogram.begin(), histogram.end(), 0.0);
      }
      break;
    case Lambda::WangLandauPhase::Finalize:
      std::fill(histogram.begin(), histogram.end(), 0.0);

      double normalize = biasFactor[0];
      for(double &bias : biasFactor)
      {
        bias -= normalize;
      }
      break;
  }
}
