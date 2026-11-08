export module bool3;

import std;

// A union like int3 and double3: b[0..2] are x, y and z.
export union bool3
{
  bool b[3];
  struct
  {
    bool x, y, z;
  };

  bool3() : x(false), y(false), z(false) {}
  bool3(bool x, bool y, bool z) : x(x), y(y), z(z) {}
  inline bool& operator[](std::size_t i) { return b[i]; }
  inline const bool& operator[](std::size_t i) const { return b[i]; }
};