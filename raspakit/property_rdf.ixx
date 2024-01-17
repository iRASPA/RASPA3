export module property_rdf;

import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;


export struct PropertyRadialDistributionFunction
{
  PropertyRadialDistributionFunction(size_t numberOfBlocks, size_t numberOfBins) :
    numberOfBlocks(numberOfBlocks),
    numberOfBins(numberOfBins),
    sumProperty(std::vector<std::vector<double>>(numberOfBlocks, std::vector<double>(numberOfBins)))
  {
  }

  size_t numberOfBlocks;
  size_t numberOfBins;
  std::vector<std::vector<double>> sumProperty;
};
