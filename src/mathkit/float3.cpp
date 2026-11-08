module;

module float3;

import std;

float3 float3::normalise()
{
  float magnitude = std::sqrt((x * x) + (y * y) + (z * z));

  if (magnitude != 0.0f)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }
  return *this;
}

float3 float3::fract()
{
  // the same form as double3::fract: x - floor(x), with 1.0f (the rounding of a tiny negative x plus one)
  // wrapped to 0.0f
  float3 s(x - std::floor(x), y - std::floor(y), z - std::floor(z));
  if (s.x >= 1.0f) s.x = 0.0f;
  if (s.y >= 1.0f) s.y = 0.0f;
  if (s.z >= 1.0f) s.z = 0.0f;
  return s;
}

std::ostream& operator<<(std::ostream& out, const float3& vec)  // output
{
  out << vec.x;
  out << vec.y;
  out << vec.z;
  return out;
}
