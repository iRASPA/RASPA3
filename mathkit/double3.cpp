module;

module double3;

import <iostream>;
import <fstream>;
import <cmath>;
import <complex>;

import archive;

double3 double3::normalise()
{
  double magnitude = sqrt((x * x) + (y * y) + (z * z));

  if (magnitude != 0)
  {
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }
  return *this;
}

double3 double3::fract() const
{
  double3 s = double3(x, y, z);
  s.x -= std::rint(x);
  s.y -= std::rint(y);
  s.z -= std::rint(z);

  if (s.x < 0.0)
  {
    s.x += 1.0;
  }
  if (s.x > 1.0)
  {
    s.x -= 1.0;
  }

  if (s.y < 0.0)
  {
    s.y += 1.0;
  }
  if (s.y > 1.0)
  {
    s.y -= 1.0;
  }

  if (s.z < 0.0)
  {
    s.z += 1.0;
  }
  if (s.z > 1.0)
  {
    s.z -= 1.0;
  }
  return s;
}

std::ostream& operator<<(std::ostream& out, const double3& vec) // output
{
  out << vec.x;
  out << ",";
  out << vec.y;
  out << ",";
  out << vec.z;
  return out;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const double3 &vec)
{
  archive << vec.x;
  archive << vec.y;
  archive << vec.z;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, double3 &vec)
{
  archive >> vec.x;
  archive >> vec.y;
  archive >> vec.z;
  return archive;
}

