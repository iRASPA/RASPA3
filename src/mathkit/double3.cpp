module;

module double3;

import std;

import archive;

double3 double3::normalized() const
{
  double magnitude = std::sqrt((x * x) + (y * y) + (z * z));

  if (magnitude != 0)
  {
    return double3(x / magnitude, y / magnitude, z / magnitude);
  }
  return *this;
}

double3 double3::fract() const
{
  // one implementation: x - floor(x), with 1.0 (the rounding of a tiny negative x plus one) wrapped to 0.0
  return double3::fract(*this);
}

std::string double3::to_string()
{
  return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")";
}

std::ostream &operator<<(std::ostream &out, const double3 &vec)  // output
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
  // the four words of the memory image (the padding word among them), so that a vector of these streams in bulk
  archive << vec.x;
  archive << vec.y;
  archive << vec.z;
  archive << vec.w;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, double3 &vec)
{
  archive >> vec.x;
  archive >> vec.y;
  archive >> vec.z;
  archive >> vec.w;
  return archive;
}
