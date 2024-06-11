module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#endif

export module hardware_info;

#ifndef USE_LEGACY_HEADERS
import <string>;
#endif

import json;

export namespace HardwareInfo
{
  std::string writeInfo();
  nlohmann::json jsonInfo();
}
