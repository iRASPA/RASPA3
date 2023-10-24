module;

module int3;

import <fstream>;
import <complex>;

import ring;
import archive;

int3 int3::greatestCommonDivisor(int3 a, int b)
{
  return int3(Ring::greatestCommonDivisor(a.x, b), Ring::greatestCommonDivisor(a.y, b), Ring::greatestCommonDivisor(a.z, b));
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const int3 &vec)
{
  archive << vec.x << vec.y << vec.z;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, int3 &vec)
{
  archive >> vec.x >> vec.y >> vec.z;
  return archive;
}

