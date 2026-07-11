module;

export module property_energy;

import std;

import archive;
import averages;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import framework;
import component;
import json;
export import property_block_average;

export struct PropertyEnergy : BlockAverage<EnergyStatus>
{
  PropertyEnergy() = default;

  PropertyEnergy(std::size_t numberOfBlocks, std::size_t numberOfExternalFields, std::size_t numberOfFrameworks,
                 std::size_t numberOfComponents)
      : BlockAverage<EnergyStatus>(numberOfBlocks,
                                   EnergyStatus(numberOfExternalFields, numberOfFrameworks, numberOfComponents))
  {
  }

  std::vector<EnergyStatus> blockEnergy() const
  {
    std::vector<EnergyStatus> blockEnergies(numberOfBlocks);
    std::transform(bookKeeping.begin(), bookKeeping.end(), blockEnergies.begin(),
                   [](std::pair<EnergyStatus, double> block) { return block.first / block.second; });
    return blockEnergies;
  }

  std::string writeAveragesStatistics(bool externalField, std::optional<Framework> &framework,
                                      std::vector<Component> &components) const;
  nlohmann::json jsonAveragesStatistics(bool externalField, std::optional<Framework> &framework,
                                        std::vector<Component> &components) const;

  std::string repr() const;
};
