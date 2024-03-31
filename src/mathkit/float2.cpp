module;

module float2;

import <sstream>;

std::ostream& operator<<(std::ostream& stream, const float2& vec)
{
    return stream << vec.x << vec.y;
}


std::istream& operator>>(std::istream& stream, float2& vec)
{
    return stream >> vec.x >> vec.y;
}
