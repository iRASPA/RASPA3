module;

#ifdef USE_LEGACY_HEADERS
#include <ostream>
#include <fstream>
#endif 


export module simd_quatd;

#ifndef USE_LEGACY_HEADERS
import <ostream>;
import <fstream>;
#endif


import double3;
import archive;

export union simd_quatd
{
  double v[4];
  struct { double ix, iy, iz, r; };

  simd_quatd(double ix, double iy, double iz, double r);
  simd_quatd(double real, double3 imag);
  simd_quatd(double3 EulerAngles);
  static simd_quatd fromAxisAngle(double angle, double3 axis);
  double3 EulerAngles();
  simd_quatd normalized();
  static simd_quatd yaw(double angle);
  static simd_quatd pitch(double angle);
  static simd_quatd roll(double angle);

  static const simd_quatd data120[120];
  static const simd_quatd data60[60];
  static const simd_quatd data600[600];
  static const simd_quatd data300[300];
  static const simd_quatd data1992[1992];
  static const double weights1992[1992];
  static const simd_quatd data360[360];
  static const double weights360[360];

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const simd_quatd &q);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, simd_quatd &q);
};

export inline simd_quatd operator+(const simd_quatd& a, const simd_quatd& b)
{
  return simd_quatd(a.ix + b.ix, a.iy + b.iy, a.iz + b.iz, a.r + b.r);
}

export inline simd_quatd operator/(const simd_quatd& a, const double& b)
{
  return simd_quatd(a.ix / b, a.iy / b, a.iz / b, a.r / b);
}

export inline simd_quatd operator*(const simd_quatd& a, const simd_quatd& b)
{
  return simd_quatd(a.r * b.r - a.ix * b.ix - a.iy * b.iy - a.iz * b.iz,
        double3(a.r * b.ix + a.ix * b.r + a.iy * b.iz - a.iz * b.iy,
            a.r * b.iy - a.ix * b.iz + a.iy * b.r + a.iz * b.ix,
            a.r * b.iz + a.ix * b.iy - a.iy * b.ix + a.iz * b.r));
}
