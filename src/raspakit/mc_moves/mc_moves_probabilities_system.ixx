export module mc_moves_probabilities_system;

import <string>;
import <chrono>;
import <fstream>;

import archive;
import double3;

export struct MCMoveProbabilitiesSystem
{
  uint64_t versionNumber{ 1 };

  bool operator==(MCMoveProbabilitiesSystem const&) const = default;

  double probabilityVolumeMove{ 0.0 };
  double probabilityGibbsVolumeMove{ 0.0 };

  void optimizeAcceptance();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesSystem &p);
};
