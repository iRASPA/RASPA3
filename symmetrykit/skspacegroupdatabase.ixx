export module skspacegroupdatabase;

import <array>;
import <vector>;

import skspacegroupsetting;

export struct SKSpaceGroupDataBase
{
	SKSpaceGroupDataBase();

	static const std::array<SKSpaceGroupSetting, 531> spaceGroupData;
	static const std::vector<std::vector<int>> spaceGroupHallData;
};