module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#endif

module skonethirdseitzmatrix;

#ifndef USE_LEGACY_HEADERS
import <string>;
#endif

SKOneThirdSeitzMatrix::SKOneThirdSeitzMatrix(std::string text, uint8_t encoding, int8_t r1, int8_t r2, int8_t r3, int8_t t)
{
    _text = text;
    _encoding = encoding;
    _r1 = r1;
    _r2 = r2;
    _r3 = r3;
    _t = t;
}
