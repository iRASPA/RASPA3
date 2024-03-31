module;

module enthalpy_of_adsorption;

import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#else
  import print;
#endif


import archive;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorption &p)
{
  archive << p.size;
  archive << p.values;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorption &p)
{
  archive >> p.size;
  archive >> p.values;

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnthalpyOfAdsorptionTerms &p)
{
  archive << p.size;
  archive << p.swapableComponents;
  archive << p.totalEnergyTimesNumberOfMolecules;
  archive << p.numberOfMoleculesSquared;
  archive << p.numberOfMolecules;
  archive << p.temperature;
  archive << p.totalEnergy;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnthalpyOfAdsorptionTerms &p)
{
  archive >> p.size;
  archive >> p.swapableComponents;
  archive >> p.totalEnergyTimesNumberOfMolecules;
  archive >> p.numberOfMoleculesSquared;
  archive >> p.numberOfMolecules;
  archive >> p.temperature;
  archive >> p.totalEnergy;

  return archive;
}



