module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

module float3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

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
  float3 s = float3(x, y, z);
  s.x -= std::rint(x);
  s.y -= std::rint(y);
  s.z -= std::rint(z);

  if (s.x < 0.0f)
  {
    s.x += 1.0f;
  }
  if (s.x > 1.0f)
  {
    s.x -= 1.0f;
  }

  if (s.y < 0.0f)
  {
    s.y += 1.0f;
  }
  if (s.y > 1.0f)
  {
    s.y -= 1.0f;
  }

  if (s.z < 0.0f)
  {
    s.z += 1.0f;
  }
  if (s.z > 1.0f)
  {
    s.z -= 1.0f;
  }
  return s;
}

std::ostream& operator<<(std::ostream& out, const float3& vec)  // output
{
  out << vec.x;
  out << vec.y;
  out << vec.z;
  return out;
}
