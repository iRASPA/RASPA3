module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdint>
#include <string>
#endif

module skonethirdseitzmatrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

SKOneThirdSeitzMatrix::SKOneThirdSeitzMatrix(std::string text, std::uint8_t encoding, std::int8_t r1, 
                                             std::int8_t r2, std::int8_t r3, std::int8_t t)
{
  _text = text;
  _encoding = encoding;
  _r1 = r1;
  _r2 = r2;
  _r3 = r3;
  _t = t;
}
