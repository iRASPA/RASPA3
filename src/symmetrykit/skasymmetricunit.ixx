module;

export module skasymmetricunit;

import std;

import int3;
import int3x3;
import double3;
import double3x3;

export struct SKAsymmetricUnit
{
  std::pair<int, int> a;
  std::pair<int, int> b;
  std::pair<int, int> c;

  SKAsymmetricUnit();
  SKAsymmetricUnit(std::pair<int, int> p1, std::pair<int, int> p2, std::pair<int, int> p3) : a(p1), b(p2), c(p3) {};

  SKAsymmetricUnit(const SKAsymmetricUnit& t)
  {
    this->a = t.a;
    this->b = t.b;
    this->c = t.c;
  }

  SKAsymmetricUnit& operator=(const SKAsymmetricUnit& other)
  {
    this->a = other.a;
    this->b = other.b;
    this->c = other.c;
    return *this;
  }

  // The image of 'point' (shifted by whole cells along each axis) that lies in the unit, if there is one. The unit
  // may extend outside [0, 1) (141 settings do, e.g. z in [-1/2, 0]), so this image, not the point wrapped into the
  // cell, is the representative. Closed bounds are widened by 'tolerance' so that round-off at a special position
  // (-1e-17 for 0) does not put a whole orbit outside; open bounds stay exact.
  std::optional<double3> imageInside(double3 point, double tolerance = 1e-8) const;
  bool contains(double3 point, double tolerance = 1e-8) const { return imageInside(point, tolerance).has_value(); }

  // Note: the unit is stored as an axis-aligned box in 1/48ths. For settings whose asymmetric unit is bounded by
  // diagonal planes this box encloses the unit, so a point inside the box need not be the unique representative.

 private:
  std::optional<double> coordinateInside(double point, std::pair<int, int> range, double tolerance) const;
};
