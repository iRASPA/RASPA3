module;

module energy_factor;

import <fstream>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import archive;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyFactor &e)
{
  archive << e.energy;
  archive << e.dUdlambda;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyFactor &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;

  return archive;
}
