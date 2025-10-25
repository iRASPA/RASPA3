module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <sstream>
#endif

module float2;

#ifdef USE_STD_IMPORT
import std;
#endif

std::ostream& operator<<(std::ostream& stream, const float2& vec) { return stream << vec.x << vec.y; }

std::istream& operator>>(std::istream& stream, float2& vec) { return stream >> vec.x >> vec.y; }
