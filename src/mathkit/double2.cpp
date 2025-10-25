module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <sstream>
#endif

module double2;

#ifdef USE_STD_IMPORT
import std;
#endif

std::ostream& operator<<(std::ostream& stream, const double2& vec) { return stream << vec.x << vec.y; }

std::istream& operator>>(std::istream& stream, double2& vec) { return stream >> vec.x >> vec.y; }
