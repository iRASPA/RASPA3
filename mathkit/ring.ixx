export module ring;

import <tuple>;

export struct Ring
{
	Ring();
	static int floorDivision(int a, int b);
	static int modulo(int a, int b);
	static int greatestCommonDivisor(int a, int b);
	static std::tuple<int, int, int> extendedGreatestCommonDivisor(int a, int b);
	static std::pair<int, int> divisionModulo(int a, int b);
};