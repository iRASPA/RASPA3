module;

#ifdef USE_LEGACY_HEADERS
#include <sstream>
#endif

module double2;

#ifndef USE_LEGACY_HEADERS
import <sstream>;
#endif

std::ostream& operator<<(std::ostream& stream, const double2& vec)
{
  return stream << vec.x << vec.y;
}


std::istream& operator>>(std::istream& stream, double2& vec)
{
  return stream >> vec.x >> vec.y;
}
