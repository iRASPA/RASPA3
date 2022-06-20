module;

#include <sstream>

module double2;



std::ostream& operator<<(std::ostream& stream, const double2& vec)
{
	return stream << vec.x << vec.y;
}


std::istream& operator>>(std::istream& stream, double2& vec)
{
	return stream >> vec.x >> vec.y;
}