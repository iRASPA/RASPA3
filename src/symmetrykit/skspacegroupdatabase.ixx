module;

export module skspacegroupdatabase;

import std;

import skspacegroupsetting;

export struct SKSpaceGroupDataBase
{
  SKSpaceGroupDataBase();

  static const std::array<SKSpaceGroupSetting, 531> spaceGroupData;
  static const std::vector<std::vector<std::size_t>> spaceGroupHallData;
};
