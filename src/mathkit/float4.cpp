module;

module float4;

import std;

void float4::normalise()
{
  float magnitude = std::sqrt((x * x) + (y * y) + (z * z) * (w * w));

  if (magnitude != 0)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
    w /= magnitude;
  }
}
