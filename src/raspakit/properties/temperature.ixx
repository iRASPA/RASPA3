module;

export module property_temperature;

import std;

import archive;
import averages;
export import property_block_average;

export struct PropertyTemperature : BlockAverage<double>
{
  PropertyTemperature() = default;

  PropertyTemperature(std::size_t numberOfBlocks) : BlockAverage<double>(numberOfBlocks) {}

  std::string writeAveragesStatistics(const std::string tag);
};
