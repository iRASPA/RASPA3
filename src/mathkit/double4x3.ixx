module;

export module double4x3;

import std;

import double3;
import double4;

// Four column vectors. A double3 holds four doubles (x, y, z and the padding w), so a flat array of twelve doubles
// (or named elements) would not overlay the columns: only the columns are accessible.
export union double4x3
{
  double3 v[4];

  double4x3(double3 v1, double3 v2, double3 v3, double3 v4) : v{v1, v2, v3, v4} {}

  inline double3& operator[](std::size_t i) { return v[i]; }
  inline const double3& operator[](std::size_t i) const { return v[i]; }
};
