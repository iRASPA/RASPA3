module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numeric>
#include <vector>
#endif

module int3x3;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import ring;
import int3;

int3x3::int3x3() {}

int3x3::int3x3(int value)
{
  m11 = value;
  m22 = value;
  m33 = value;
}

int3x3::int3x3(int3 v1, int3 v2, int3 v3)
{
  m11 = v1.x;
  m21 = v1.y;
  m31 = v1.z;
  m12 = v2.x;
  m22 = v2.y;
  m32 = v2.z;
  m13 = v3.x;
  m23 = v3.y;
  m33 = v3.z;
}

int int3x3::determinant(void) const
{
  int determinant = +m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  return determinant;
}

int3x3 int3x3::adjugate() const
{
  int3x3 result{};
  result.m11 = this->m22 * this->m33 - this->m32 * this->m23;
  result.m12 = this->m13 * this->m32 - this->m12 * this->m33;
  result.m13 = this->m12 * this->m23 - this->m13 * this->m22;
  result.m21 = this->m23 * this->m31 - this->m21 * this->m33;
  result.m22 = this->m11 * this->m33 - this->m13 * this->m31;
  result.m23 = this->m21 * this->m13 - this->m11 * this->m23;
  result.m31 = this->m21 * this->m32 - this->m31 * this->m22;
  result.m32 = this->m31 * this->m12 - this->m11 * this->m32;
  result.m33 = this->m11 * this->m22 - this->m21 * this->m12;

  return result;
}

int int3x3::greatestCommonDivisor()
{
  std::vector<int> v1{m11, m12, m13, m21, m22, m23, m31, m32, m33};
  return std::accumulate(v1.begin(), v1.end(), 0, [&](int a, int b) { return Ring::greatestCommonDivisor(a, b); });
}

int int3x3::trace() const { return m11 + m22 + m33; }
