module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#endif

export module double4;

//#if defined(WIN32)
//    import <intrin.h>;
//#elif defined(__AVX__)
//    import <immintrin.h>;
//#endif

#ifndef USE_LEGACY_HEADERS
import <fstream>;
#endif

import archive;


export union double4
{
//    #ifdef __AVX__
//      __m256d value;
//    #endif
    double v[4];
    struct { double x, y, z, w; };

    double4(double x = 0, double y = 0, double z = 0, double w = 0) :x(x), y(y), z(z), w(w) {}

    bool operator==(double4 const& rhs) const
    {
      return (x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (w == rhs.w);
    }

    inline double& operator [] (int i) { return v[i]; }
    inline const double& operator [] (int i) const { return v[i]; }

    void normalise();
    inline static double dot(const double4& v1, const double4& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w; }

    // http://www.gamedev.net/topic/456301-cross-product-vector-4d/
    inline static double4 cross(const double4& v1, const double4& v2) { return double4(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x, v1.x * v2.y - v1.y * v2.x, 0.0); }
    inline static double4 cross(const double4& v1, const double4& v2, const double4& v3)
    {
        return double4(v1.y * (v2.z * v3.w - v3.z * v2.w) - v1.z * (v2.y * v3.w - v3.y * v2.w) + v1.w * (v2.y * v3.z - v3.y * v2.z),
            -v1.x * (v2.z * v3.w - v3.z * v2.w) + v1.z * (v2.x * v3.w - v3.x * v2.w) - v1.w * (v2.x * v3.z - v3.x * v2.z),
            v1.x * (v2.y * v3.w - v3.y * v2.w) - v1.y * (v2.x * v3.w - v3.x * v2.w) + v1.w * (v2.x * v3.y - v3.x * v2.y),
            -v1.x * (v2.y * v3.z - v3.y * v2.z) + v1.y * (v2.x * v3.z - v3.x * v2.z) - v1.z * (v2.x * v3.y - v3.x * v2.y));
    }

    double4 operator-() const { return double4(-this->x, -this->y, -this->z, -this->w); }
    double4& operator+=(const double4& b) { this->x += b.x, this->y += b.y, this->z += b.z; this->w += b.w; return *this; }
    double4& operator-=(const double4& b) { this->x -= b.x, this->y -= b.y, this->z -= b.z; this->w += b.w; return *this; }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const double4 &vec);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, double4 &vec);
};

export inline double4 operator+(const double4& a, const double4& b)
{
    return double4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

export inline double4 operator*(const double4& a, const double4& b)
{
    return double4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

export inline double4 operator*(const double4& a, const double& b)
{
    return double4(a.x * b, a.y * b, a.z * b, a.w * b);
}

export inline double4 operator*(const double& a, const double4& b)
{
    return double4(a * b.x, a * b.y, a * b.z, a * b.w);
}
