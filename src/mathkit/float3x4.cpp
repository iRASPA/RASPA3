module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>
#endif

module float3x4;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
import <ostream>;
#endif

float3x4::float3x4(const float3x4& a)  // copy constructor
{
  m11 = a.m11;
  m21 = a.m21;
  m31 = a.m31;
  m41 = a.m41;
  m12 = a.m12;
  m22 = a.m22;
  m32 = a.m32;
  m42 = a.m42;
  m13 = a.m13;
  m23 = a.m23;
  m33 = a.m33;
  m43 = a.m43;
}

float3x4& float3x4::operator=(const float3x4& a)  // assignment operator
{
  if (this != &a)  // protect against self-assingment
  {
    m11 = a.m11;
    m21 = a.m21;
    m31 = a.m31;
    m41 = a.m41;
    m12 = a.m12;
    m22 = a.m22;
    m32 = a.m32;
    m42 = a.m42;
    m13 = a.m13;
    m23 = a.m23;
    m33 = a.m33;
    m43 = a.m43;
  }
  return *this;
}
