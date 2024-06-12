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
  
module energy_status_intra;

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


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyIntra &e)
{
  archive << e.versionNumber;

  archive << e.bond;
  archive << e.bend;
  archive << e.inversionBend;
  archive << e.ureyBradley;
  archive << e.torsion;
  archive << e.improperTorsion;
  archive << e.bondBond;
  archive << e.bondBend;
  archive << e.bondTorsion;
  archive << e.bendBend;
  archive << e.bendTorsion;
  archive << e.intraVDW;
  archive << e.intraChargeCharge;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyIntra &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyIntra' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.bond;
  archive >> e.bend;
  archive >> e.inversionBend;
  archive >> e.ureyBradley;
  archive >> e.torsion;
  archive >> e.improperTorsion;
  archive >> e.bondBond;
  archive >> e.bondBend;
  archive >> e.bondTorsion;
  archive >> e.bendBend;
  archive >> e.bendTorsion;
  archive >> e.intraVDW;
  archive >> e.intraChargeCharge;

  return archive;
}

