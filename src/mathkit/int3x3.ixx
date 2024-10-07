export module int3x3;

import int3;
import double3;
import hashcombine;

export union int3x3
{
  int m[9];
  int mm[3][3];
  int3 v[3];  // columns
  struct
  {
    int m11, m21, m31,  // 1st column
        m12, m22, m32,  // 2nd column
        m13, m23, m33;  // 3rd column
  };
  struct
  {
    int ax, ay, az,  // 1st column
        bx, by, bz,  // 2nd column
        cx, cy, cz;  // 3rd column
  };

  int3x3();
  int3x3(int value);
  int3x3(int3 v1, int3 v2, int3 v3);

  inline int3& operator[](int i) { return v[i]; }
  inline const int3& operator[](int i) const { return v[i]; }

  inline union int3x3& operator+=(const union int3x3& b)
  {
    this->m11 += b.m11;
    this->m21 += b.m21;
    this->m31 += b.m31;
    this->m12 += b.m12;
    this->m22 += b.m22;
    this->m32 += b.m32;
    this->m13 += b.m13;
    this->m23 += b.m23;
    this->m33 += b.m33;
    return *this;
  }
  inline union int3x3& operator-=(const union int3x3& b)
  {
    this->m11 -= b.m11;
    this->m21 -= b.m21;
    this->m31 -= b.m31;
    this->m12 -= b.m12;
    this->m22 -= b.m22;
    this->m32 -= b.m32;
    this->m13 -= b.m13;
    this->m23 -= b.m23;
    this->m33 -= b.m33;
    return *this;
  }

  union int3x3 operator-() const
  {
    union int3x3 r;
    r.m11 = -m11;
    r.m21 = -m21;
    r.m31 = -m31;
    r.m12 = -m12;
    r.m22 = -m22;
    r.m32 = -m32;
    r.m13 = -m13;
    r.m23 = -m23;
    r.m33 = -m33;
    return r;
  }

  inline bool operator==(const union int3x3& b) const
  {
    return (this->m11 == b.m11) && (this->m21 == b.m21) && (this->m31 == b.m31) && (this->m12 == b.m12) &&
           (this->m22 == b.m22) && (this->m32 == b.m32) && (this->m13 == b.m13) && (this->m23 == b.m23) &&
           (this->m33 == b.m33);
  }

  int greatestCommonDivisor();
  int determinant(void) const;
  int3x3 adjugate() const;
  int trace(void) const;
};

export inline union int3x3 operator+(const union int3x3& a, const union int3x3& b)
{
  union int3x3 r;

  r.m11 = a.m11 + b.m11;
  r.m21 = a.m21 + b.m21;
  r.m31 = a.m31 + b.m31;
  r.m12 = a.m12 + b.m12;
  r.m22 = a.m22 + b.m22;
  r.m32 = a.m32 + b.m32;
  r.m13 = a.m13 + b.m13;
  r.m23 = a.m23 + b.m23;
  r.m33 = a.m33 + b.m33;

  return r;
}

export inline union int3x3 operator-(const union int3x3& a, const union int3x3& b)
{
  union int3x3 r;

  r.m11 = a.m11 - b.m11;
  r.m21 = a.m21 - b.m21;
  r.m31 = a.m31 - b.m31;
  r.m12 = a.m12 - b.m12;
  r.m22 = a.m22 - b.m22;
  r.m32 = a.m32 - b.m32;
  r.m13 = a.m13 - b.m13;
  r.m23 = a.m23 - b.m23;
  r.m33 = a.m33 - b.m33;

  return r;
}

export inline int3x3 operator*(const union int3x3& a, const union int3x3& b)
{
  union int3x3 r;

  r.m11 = a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31;
  r.m21 = a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31;
  r.m31 = a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31;

  r.m12 = a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32;
  r.m22 = a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32;
  r.m32 = a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32;

  r.m13 = a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33;
  r.m23 = a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33;
  r.m33 = a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33;

  return r;
}

export inline int3 operator*(const union int3x3& a, const int3& b)
{
  int3 r(0, 0, 0);

  r.x = a.m11 * b.x + a.m12 * b.y + a.m13 * b.z;
  r.y = a.m21 * b.x + a.m22 * b.y + a.m23 * b.z;
  r.z = a.m31 * b.x + a.m32 * b.y + a.m33 * b.z;

  return r;
}

export inline double3 operator*(const union int3x3& a, const double3& b)
{
  double3 r(0, 0, 0);

  r.x = a.m11 * b.x + a.m12 * b.y + a.m13 * b.z;
  r.y = a.m21 * b.x + a.m22 * b.y + a.m23 * b.z;
  r.z = a.m31 * b.x + a.m32 * b.y + a.m33 * b.z;

  return r;
}

export inline union int3x3 operator/(const union int3x3& a, const int& b)
{
  union int3x3 r;

  r.m11 = a.m11 / b;
  r.m21 = a.m21 / b;
  r.m31 = a.m31 / b;
  r.m12 = a.m12 / b;
  r.m22 = a.m22 / b;
  r.m32 = a.m32 / b;
  r.m13 = a.m13 / b;
  r.m23 = a.m23 / b;
  r.m33 = a.m33 / b;

  return r;
}

/*
namespace std
{
    template <> struct hash<union int3x3>
    {
        size_t operator()(const union int3x3& k) const
        {
            std::size_t h = 0;
            hash_combine(h, k.m11);
            hash_combine(h, k.m12);
            hash_combine(h, k.m13);
            hash_combine(h, k.m21);
            hash_combine(h, k.m22);
            hash_combine(h, k.m23);
            hash_combine(h, k.m31);
            hash_combine(h, k.m32);
            hash_combine(h, k.m33);
            return h;
        }
    };
}*/
