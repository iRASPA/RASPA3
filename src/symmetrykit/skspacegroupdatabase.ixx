module;

export module skspacegroupdatabase;

import std;

import skspacegroupsetting;

export struct SKSpaceGroupDataBase
{
  SKSpaceGroupDataBase();

  static const std::array<SKSpaceGroupSetting, 531> spaceGroupData;

  /// The setting of a Hall number, checked: a number past the table (a structure given one by hand) is refused
  /// instead of read past the end of it.
  static const SKSpaceGroupSetting& setting(std::size_t hallNumber)
  {
    if (hallNumber >= spaceGroupData.size())
    {
      throw std::out_of_range(std::format("space-group Hall number {} is not in 0..{}", hallNumber,
                                          spaceGroupData.size() - 1));
    }
    return spaceGroupData[hallNumber];
  }
  static const std::vector<std::vector<std::size_t>> spaceGroupHallData;
};
