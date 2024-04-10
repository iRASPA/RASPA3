module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif
  
module energy_status_inter;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif
  
import archive;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyInter &e)
{
  archive << e.versionNumber;

  archive << e.VanDerWaals;
  archive << e.VanDerWaalsTailCorrection;
  archive << e.CoulombicReal;
  archive << e.CoulombicFourier;
  archive << e.totalInter;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyInter &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyInter' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.VanDerWaals;
  archive >> e.VanDerWaalsTailCorrection;
  archive >> e.CoulombicReal;
  archive >> e.CoulombicFourier;
  archive >> e.totalInter;

  return archive;
}
