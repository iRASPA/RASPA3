module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <ostream>
#endif

export module double3x3x3;

// #if defined(WIN32)
//     import <intrin.h>;
// #elif defined(__AVX__)
//     import <immintrin.h>;
// #endif

#ifndef USE_LEGACY_HEADERS
import <istream>;
import <ostream>;
import <fstream>;
import <cmath>;
#endif

import archive;
import double3;
import simd_quatd;
import double3x3;
import json;

export union double3x3x3
{
  double m[64];
  double mm[4][4][4];
  double3x3 v[4];
  struct
  {
    double m111, m121, m131, m141, m112, m122, m132, m142, m113, m123, m133, m143, m114, m124, m134, m144, m211, m221,
        m231, m241, m212, m222, m232, m242, m213, m223, m233, m243, m214, m224, m234, m244, m311, m321, m331, m341,
        m312, m322, m332, m342, m313, m323, m333, m343, m314, m324, m334, m344, m411, m421, m431, m441, m412, m422,
        m432, m442, m413, m423, m433, m443, m414, m424, m434, m444;
  };

  double3x3x3()
      : m111(0.0),
        m121(0.0),
        m131(0.0),
        m141(0.0),
        m112(0.0),
        m122(0.0),
        m132(0.0),
        m142(0.0),
        m113(0.0),
        m123(0.0),
        m133(0.0),
        m143(0.0),
        m114(0.0),
        m124(0.0),
        m134(0.0),
        m144(0.0),
        m211(0.0),
        m221(0.0),
        m231(0.0),
        m241(0.0),
        m212(0.0),
        m222(0.0),
        m232(0.0),
        m242(0.0),
        m213(0.0),
        m223(0.0),
        m233(0.0),
        m243(0.0),
        m214(0.0),
        m224(0.0),
        m234(0.0),
        m244(0.0),
        m311(0.0),
        m321(0.0),
        m331(0.0),
        m341(0.0),
        m312(0.0),
        m322(0.0),
        m332(0.0),
        m342(0.0),
        m313(0.0),
        m323(0.0),
        m333(0.0),
        m343(0.0),
        m314(0.0),
        m324(0.0),
        m334(0.0),
        m344(0.0),
        m411(0.0),
        m421(0.0),
        m431(0.0),
        m441(0.0),
        m412(0.0),
        m422(0.0),
        m432(0.0),
        m442(0.0),
        m413(0.0),
        m423(0.0),
        m433(0.0),
        m443(0.0),
        m414(0.0),
        m424(0.0),
        m434(0.0),
        m444(0.0) {};

  inline double3x3x3& operator+=(const double3x3x3& b)
  {
    this->m111 += b.m111;
    this->m121 += b.m121;
    this->m131 += b.m131;
    this->m112 += b.m112;
    this->m122 += b.m122;
    this->m132 += b.m132;
    this->m113 += b.m113;
    this->m123 += b.m123;
    this->m133 += b.m133;

    this->m211 += b.m211;
    this->m221 += b.m221;
    this->m231 += b.m231;
    this->m212 += b.m212;
    this->m222 += b.m222;
    this->m232 += b.m232;
    this->m213 += b.m213;
    this->m223 += b.m223;
    this->m233 += b.m233;

    this->m311 += b.m311;
    this->m321 += b.m321;
    this->m331 += b.m331;
    this->m312 += b.m312;
    this->m322 += b.m322;
    this->m332 += b.m332;
    this->m313 += b.m313;
    this->m323 += b.m323;
    this->m333 += b.m333;

    return *this;
  }
};
