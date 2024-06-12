module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <print>
#endif

module enthalpy_of_adsorption;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
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



