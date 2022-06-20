module;

module int3;

import ring;

int3 int3::greatestCommonDivisor(int3 a, int b)
{
  return int3(Ring::greatestCommonDivisor(a.x, b), Ring::greatestCommonDivisor(a.y, b), Ring::greatestCommonDivisor(a.z, b));
}

