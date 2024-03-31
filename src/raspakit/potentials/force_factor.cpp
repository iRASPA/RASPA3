module;

module force_factor;

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


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceFactor &e)
{
  archive << e.energy;
  archive << e.forceFactor;
  archive << e.dUdlambda;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceFactor &e)
{
  archive >> e.energy;
  archive >> e.forceFactor;
  archive >> e.dUdlambda;

  return archive;
}
