module;

module uint3;

import std;

import ring;
import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const uint3 &vec)
{
  archive << vec.x << vec.y << vec.z;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, uint3 &vec)
{
  archive >> vec.x >> vec.y >> vec.z;
  return archive;
}
