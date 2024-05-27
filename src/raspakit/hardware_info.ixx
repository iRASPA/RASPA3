module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#endif

export module hardware_info;

#ifndef USE_LEGACY_HEADERS
import <string>;
#endif

import hdf5;

export namespace HardwareInfo
{
  std::string writeInfo();
  void logInfo(HDF5Writer& hdf5);
}
