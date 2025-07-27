module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <vector>
#endif

export module skspacegroupdatabase;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import skspacegroupsetting;

export struct SKSpaceGroupDataBase
{
  SKSpaceGroupDataBase();

  static const std::array<SKSpaceGroupSetting, 531> spaceGroupData;
  static const std::vector<std::vector<std::size_t>> spaceGroupHallData;
};
