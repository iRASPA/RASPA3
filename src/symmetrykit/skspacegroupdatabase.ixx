module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <vector>
#endif

export module skspacegroupdatabase;

#ifdef USE_STD_IMPORT
import std;
#endif

import skspacegroupsetting;

export struct SKSpaceGroupDataBase
{
  SKSpaceGroupDataBase();

  static const std::array<SKSpaceGroupSetting, 531> spaceGroupData;
  static const std::vector<std::vector<std::size_t>> spaceGroupHallData;
};
