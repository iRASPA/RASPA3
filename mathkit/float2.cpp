module;

#include <sstream>

module float2;


std::ostream& operator<<(std::ostream& stream, const float2& vec)
{
    return stream << vec.x << vec.y;
}


std::istream& operator>>(std::istream& stream, float2& vec)
{
    return stream >> vec.x >> vec.y;
}