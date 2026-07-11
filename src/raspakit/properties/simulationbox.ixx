module;

export module property_simulationbox;

import std;

import archive;
import averages;
import simulationbox;
export import property_block_average;

export struct PropertySimulationBox : BlockAverage<SimulationBox>
{
  PropertySimulationBox() = default;

  PropertySimulationBox(std::size_t numberOfBlocks) : BlockAverage<SimulationBox>(numberOfBlocks) {}

  bool operator==(PropertySimulationBox const &) const = default;
};
