module;

module average_energy_type;

import std;

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const AverageEnergyType &e)
{
  archive << e.totalEnergy << e.VanDerWaalsEnergy << e.CoulombEnergy << e.polarizationEnergy;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, AverageEnergyType &e)
{
  archive >> e.totalEnergy >> e.VanDerWaalsEnergy >> e.CoulombEnergy >> e.polarizationEnergy;
  return archive;
}
