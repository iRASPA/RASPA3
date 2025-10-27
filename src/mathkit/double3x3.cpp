module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <map>
#include <numbers>
#include <ostream>
#include <print>
#include <span>
#include <tuple>
#include <vector>
#endif

#define sqr(x) ((x) * (x))
#define SIGN(a, b) ((b) >= 0.0 ? std::fabs(a) : -std::fabs(a))

#ifdef BLAS_ILP64
typedef long long blas_int;
#else
typedef int blas_int;
#endif

/* DSYEV prototype */
extern "C"
{
  void dsyev_(char* jobz, char* uplo, blas_int* n, double* a, blas_int* lda, double* w, double* work, blas_int* lwork,
              blas_int* info);

  void dgesvd_(char* jobu, char* jobvt, blas_int* m, blas_int* n, double* a, blas_int* lda, double* s, double* u,
               blas_int* ldu, double* vt, blas_int* ldvt, double* work, blas_int* lwork, blas_int* info);
}

module double3x3;

#ifdef USE_STD_IMPORT
import std;
#endif

import int3x3;
import simd_quatd;
import double3;
import archive;

// dlambda_limit, below which two lambdas are relatively equal
double dlambda_limit = 1.0E-3;
double iszero_limit = 1.0E-20;

double3x3::double3x3(simd_quatd q)
{
  double sqw = q.r * q.r;
  double sqx = q.ix * q.ix;
  double sqy = q.iy * q.iy;
  double sqz = q.iz * q.iz;

  // invs (inverse square length) is only required if quaternion is not already normalised
  double invs = 1.0 / (sqx + sqy + sqz + sqw);
  m11 = (sqx - sqy - sqz + sqw) * invs;  // since sqw + sqx + sqy + sqz =1/invs*invs
  m22 = (-sqx + sqy - sqz + sqw) * invs;
  m33 = (-sqx - sqy + sqz + sqw) * invs;

  double tmp1 = q.ix * q.iy;
  double tmp2 = q.iz * q.r;
  m21 = 2.0 * (tmp1 + tmp2) * invs;
  m12 = 2.0 * (tmp1 - tmp2) * invs;

  tmp1 = q.ix * q.iz;
  tmp2 = q.iy * q.r;
  m31 = 2.0 * (tmp1 - tmp2) * invs;
  m13 = 2.0 * (tmp1 + tmp2) * invs;
  tmp1 = q.iy * q.iz;
  tmp2 = q.ix * q.r;
  m32 = 2.0 * (tmp1 + tmp2) * invs;
  m23 = 2.0 * (tmp1 - tmp2) * invs;
}

double3x3::double3x3(const int3x3& m)
{
  m11 = m.m11;
  m21 = m.m21;
  m31 = m.m31;
  m12 = m.m12;
  m22 = m.m22;
  m32 = m.m32;
  m13 = m.m13;
  m23 = m.m23;
  m33 = m.m33;
}
double double3x3::determinant(void) const
{
  double determinant = m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  return determinant;
}

double3x3::double3x3(double lattice[3][3])
{
  m11 = lattice[0][0];
  m21 = lattice[1][0];
  m31 = lattice[2][0];
  m12 = lattice[0][1];
  m22 = lattice[1][1];
  m32 = lattice[2][1];
  m13 = lattice[0][2];
  m23 = lattice[1][2];
  m33 = lattice[2][2];
}

double3x3 double3x3::identity()
{
  return double3x3(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0));
}

double double3x3::trace(void) const { return m11 + m22 + m33; }

const double3x3 double3x3::inverse() const
{
  double determinant = m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  double3x3 inverse;
  inverse.m11 = (m22 * m33 - m32 * m23) / determinant;
  inverse.m21 = -(m21 * m33 - m31 * m23) / determinant;
  inverse.m31 = (m21 * m32 - m31 * m22) / determinant;
  inverse.m12 = -(m12 * m33 - m32 * m13) / determinant;
  inverse.m22 = (m11 * m33 - m31 * m13) / determinant;
  inverse.m32 = -(m11 * m32 - m31 * m12) / determinant;
  inverse.m13 = (m12 * m23 - m22 * m13) / determinant;
  inverse.m23 = -(m11 * m23 - m21 * m13) / determinant;
  inverse.m33 = (m11 * m22 - m21 * m12) / determinant;

  return inverse;
}

double3x3 double3x3::inverse(const double3x3& a)
{
  double determinant = a.m11 * (a.m22 * a.m33 - a.m23 * a.m32) - a.m12 * (a.m21 * a.m33 - a.m23 * a.m31) +
                       a.m13 * (a.m21 * a.m32 - a.m22 * a.m31);

  double3x3 inverse;
  inverse.m11 = +(a.m22 * a.m33 - a.m32 * a.m23) / determinant;
  inverse.m21 = -(a.m21 * a.m33 - a.m31 * a.m23) / determinant;
  inverse.m31 = +(a.m21 * a.m32 - a.m31 * a.m22) / determinant;
  inverse.m12 = -(a.m12 * a.m33 - a.m32 * a.m13) / determinant;
  inverse.m22 = +(a.m11 * a.m33 - a.m31 * a.m13) / determinant;
  inverse.m32 = -(a.m11 * a.m32 - a.m31 * a.m12) / determinant;
  inverse.m13 = +(a.m12 * a.m23 - a.m22 * a.m13) / determinant;
  inverse.m23 = -(a.m11 * a.m23 - a.m21 * a.m13) / determinant;
  inverse.m33 = +(a.m11 * a.m22 - a.m21 * a.m12) / determinant;

  return inverse;
}

double3x3 double3x3::transpose(const double3x3& right)
{
  double3x3 res;

  res.m11 = right.m11;
  res.m12 = right.m21;
  res.m13 = right.m31;
  res.m21 = right.m12;
  res.m22 = right.m22;
  res.m23 = right.m32;
  res.m31 = right.m13;
  res.m32 = right.m23;
  res.m33 = right.m33;

  return res;
}

double3x3 const double3x3::transpose(void) const
{
  double3x3 res;

  res.m11 = m11;
  res.m12 = m21;
  res.m13 = m31;
  res.m21 = m12;
  res.m22 = m22;
  res.m23 = m32;
  res.m31 = m13;
  res.m32 = m23;
  res.m33 = m33;

  return res;
}

double3x3 double3x3::inversetranpose()
{
  double determinant = +m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31);

  double3x3 inverse;
  inverse.m11 = +(m22 * m33 - m32 * m23) / determinant;
  inverse.m12 = -(m21 * m33 - m31 * m23) / determinant;
  inverse.m13 = +(m21 * m32 - m31 * m22) / determinant;
  inverse.m21 = -(m12 * m33 - m32 * m13) / determinant;
  inverse.m22 = +(m11 * m33 - m31 * m13) / determinant;
  inverse.m23 = -(m11 * m32 - m31 * m12) / determinant;
  inverse.m31 = +(m12 * m23 - m22 * m13) / determinant;
  inverse.m32 = -(m11 * m23 - m21 * m13) / determinant;
  inverse.m33 = +(m11 * m22 - m21 * m12) / determinant;

  return inverse;
}

double trunc_sqrt(double x) { return (x <= 0.0 ? 0.0 : std::sqrt(x)); }

double trunc_acos(double x)
{
  if (x >= 1.0) return 0.0;
  if (x <= -1.0) return std::numbers::pi;
  return std::acos(x);
}

static double sign(double x) { return (x < 0.0 ? -1.0 : 1.0); }

double angle(double x, double y)
{
  if (x == 0.0) return (y == 0.0 ? 0.0 : 0.5 * std::numbers::pi * sign(y));
  return (x < 0.0 ? std::atan(y / x) + std::numbers::pi * sign(y) : std::atan(y / x));
}

void double3x3::EigenSystemSymmetric(double3& eigenvalues, double3x3& eigenvectors)
{
  char decompositionJobV = 'V';
  char upload = 'U';
  std::vector<double> matrix = std::vector<double>{ax, ay, az, bx, by, bz, cx, cy, cz};
  std::vector<double> work(9 * 3);
  blas_int lwork = 9 * 3;
  std::vector<double> e = std::vector<double>(3);
  blas_int error = 0;
  blas_int N = 3;
  blas_int M = 3;

  dsyev_(&decompositionJobV, &upload, &M, matrix.data(), &N, e.data(), work.data(), &lwork, &error);

  eigenvalues = double3(e[2], e[1], e[0]);
  double3 v1, v2, v3;
  v1 = double3(matrix[0], matrix[1], matrix[2]);
  v2 = double3(matrix[3], matrix[4], matrix[5]);
  v3 = double3(matrix[6], matrix[7], matrix[8]);
  eigenvectors = double3x3(v1, v2, v3);
  if (eigenvectors.determinant() < 0)
  {
    eigenvectors = double3x3(v1, v3, v2);
  }
}

// The space-fixed coordinate system is chosen in an inertial reference frame in which Newton's equations of motion are
// valid. The body-fixed coordinate system is an origin and three orthogonal axes (unit vectors) that are fixed to the
// body and rotate, tumble, spin and twist along with it.

// The unit quaternion q is introduced in order to generate a minimal, nonsingular, representation of the rotation
// matrix from a space-fixed denoted ‘s’ to a body-fixed coordinate system denoted ‘b’. The Quaternion Rotation Operator
// is L_q(v)=q* v q; this operator represents a rotation through an angle alpha about a vector q as its axis. We apply
// the quaternion rotation operator to a 3D vector v (a pure quaternion defined in the space-fixed frame), and express
// it as w in the body-fixed frame w_b=R(q) v_s
double3x3 double3x3::buildRotationMatrix(const simd_quatd& q)
{
  double3x3 R{};

  R.ax = 2.0 * (q.r * q.r + q.ix * q.ix) - 1.0;
  R.bx = 2.0 * (q.ix * q.iy + q.r * q.iz);
  R.cx = 2.0 * (q.ix * q.iz - q.r * q.iy);
  R.ay = 2.0 * (q.ix * q.iy - q.r * q.iz);
  R.by = 2.0 * (q.iy * q.iy + q.r * q.r) - 1.0;
  R.cy = 2.0 * (q.iy * q.iz + q.r * q.ix);
  R.az = 2.0 * (q.ix * q.iz + q.r * q.iy);
  R.bz = 2.0 * (q.iy * q.iz - q.r * q.ix);
  R.cz = 2.0 * (q.r * q.r + q.iz * q.iz) - 1.0;

  return R;
}

// To do the opposite, namely apply the quaternion rotation operator to a 3D vector v (a pure quaternion defined in the
// body-fixed frame), and express it as w in the space-fixed frame, we have w_s=R(q)^-1 v_b=R(q)^T v_b
double3x3 double3x3::buildRotationMatrixInverse(const simd_quatd& q)
{
  double3x3 R{};
  R.ax = 2.0 * (q.r * q.r + q.ix * q.ix) - 1.0;
  R.bx = 2.0 * (q.ix * q.iy - q.r * q.iz);
  R.cx = 2.0 * (q.ix * q.iz + q.r * q.iy);
  R.ay = 2.0 * (q.ix * q.iy + q.r * q.iz);
  R.by = 2.0 * (q.iy * q.iy + q.r * q.r) - 1.0;
  R.cy = 2.0 * (q.iy * q.iz - q.r * q.ix);
  R.az = 2.0 * (q.ix * q.iz - q.r * q.iy);
  R.bz = 2.0 * (q.iy * q.iz + q.r * q.ix);
  R.cz = 2.0 * (q.r * q.r + q.iz * q.iz) - 1.0;

  return R;
}

// Function   | Transform a rotation matrix to a quaternion
// Parameters | The 3x3 rotation matrix R
// Note       | The algorithm returns a single quaternion q, but -q represents the same rotation
// Ref.       | Page 305 "3-D Computer Graphics A Mathematical Introduction with OpenGL", Samuel R. Buss

simd_quatd double3x3::quaternion()
{
  simd_quatd q;

  double trace = this->ax + this->by + this->cz;

  std::size_t Switch{};
  double max = trace;
  if (this->ax > max)
  {
    max = this->ax;
    Switch = 1;
  }
  if (this->by > max)
  {
    max = this->by;
    Switch = 2;
  }
  if (this->cz > max)
  {
    max = this->cz;
    Switch = 3;
  }

  switch (Switch)
  {
    case 0:
      q.r = 0.5 * std::sqrt(trace + 1.0);
      q.ix = (this->cy - this->bz) / (4.0 * q.r);
      q.iy = (this->az - this->cx) / (4.0 * q.r);
      q.iz = (this->bx - this->ay) / (4.0 * q.r);
      break;
    case 1:
      q.ix = 0.5 * std::sqrt(2.0 * this->ax - trace + 1.0);
      q.r = (this->cy - this->bz) / (4.0 * q.ix);
      q.iy = (this->bx + this->ay) / (4.0 * q.ix);
      q.iz = (this->az + this->cx) / (4.0 * q.ix);
      break;
    case 2:
      q.iy = 0.5 * std::sqrt(2.0 * this->by - trace + 1.0);
      q.r = (this->az - this->cx) / (4.0 * q.iy);
      q.ix = (this->bx + this->ay) / (4.0 * q.iy);
      q.iz = (this->cy + this->bz) / (4.0 * q.iy);
      break;
    case 3:
      q.iz = 0.5 * std::sqrt(2.0 * this->cz - trace + 1.0);
      q.r = (this->bx - this->ay) / (4.0 * q.iz);
      q.ix = (this->az + this->cx) / (4.0 * q.iz);
      q.iy = (this->cy + this->bz) / (4.0 * q.iz);
      break;
  }
  return q;
}

std::tuple<double3x3, double3, double3x3> double3x3::singularValueDecomposition() const
{
  blas_int m1 = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;

  char jobU = 'A';
  char jobVT = 'A';
  double wkopt;

  double s[3], u[9], vt[9];
  std::vector<double> matrix =
      std::vector<double>{this->ax, this->ay, this->az, this->bx, this->by, this->bz, this->cx, this->cy, this->cz};

  lwork = -1;
  dgesvd_(&jobU, &jobVT, &m1, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);

  lwork = static_cast<blas_int>(wkopt);
  std::vector<double> work(static_cast<std::size_t>(lwork));
  dgesvd_(&jobU, &jobVT, &m1, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, work.data(), &lwork, &info);

  if (info > 0)
  {
    std::print("The algorithm computing SVD failed to converge.\n");
  }

  double3 sigma = double3(s[0], s[1], s[2]);
  double3x3 matrix_u = double3x3(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
  double3x3 matrix_vt = double3x3(vt[0], vt[1], vt[2], vt[3], vt[4], vt[5], vt[6], vt[7], vt[8]);

  return {matrix_u, sigma, matrix_vt};
}

// https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
// compute the rotation matrix to map set 'A' onto set 'B'
double3x3 double3x3::computeRotationMatrix(double3 center_of_mass_A, std::span<double3> positions_A,
                                           double3 center_of_mass_B, std::span<double3> positions_B)
{
  double3x3 H{};

  for (std::size_t i = 0; i < positions_A.size(); ++i)
  {
    double3 vec_i = positions_A[i] - center_of_mass_A;
    double3 vec_j = positions_B[i] - center_of_mass_B;

    H.ax += vec_i.x * vec_j.x;
    H.ay += vec_i.y * vec_j.x;
    H.az += vec_i.z * vec_j.x;

    H.bx += vec_i.x * vec_j.y;
    H.by += vec_i.y * vec_j.y;
    H.bz += vec_i.z * vec_j.y;

    H.cx += vec_i.x * vec_j.z;
    H.cy += vec_i.y * vec_j.z;
    H.cz += vec_i.z * vec_j.z;
  }

  blas_int m1 = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;

  char jobU = 'A';
  char jobVT = 'A';
  double wkopt;

  double s[3], u[9], vt[9];
  std::vector<double> matrix = std::vector<double>{H.ax, H.ay, H.az, H.bx, H.by, H.bz, H.cx, H.cy, H.cz};

  lwork = -1;
  dgesvd_(&jobU, &jobVT, &m1, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);

  lwork = static_cast<blas_int>(wkopt);
  std::vector<double> work(static_cast<std::size_t>(lwork));
  dgesvd_(&jobU, &jobVT, &m1, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, work.data(), &lwork, &info);

  if (info > 0)
  {
    std::print("The algorithm computing SVD failed to converge.\n");
  }

  double3x3 sigma = double3x3(s[0], 0.0, 0.0, 0.0, s[1], 0.0, 0.0, 0.0, s[2]);
  double3x3 matrix_u = double3x3(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
  double3x3 matrix_vt = double3x3(vt[0], vt[1], vt[2], vt[3], vt[4], vt[5], vt[6], vt[7], vt[8]);

  double3x3 result = matrix_vt.transpose() * matrix_u.transpose();

  double determinant = result.determinant();
  if (determinant < 0.0)
  {
    return matrix_vt.transpose() * double3x3(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, -1.0)) *
           matrix_u.transpose();
  }

  return result;
}

// compute the rotation matrix to map 'vec_i' onto 'vec_j'
double3x3 double3x3::computeRotationMatrix(double3 vec_i, double3 vec_j)
{
  double3x3 H{};

  H.ax = vec_i.x * vec_j.x;
  H.ay = vec_i.y * vec_j.x;
  H.az = vec_i.z * vec_j.x;

  H.bx = vec_i.x * vec_j.y;
  H.by = vec_i.y * vec_j.y;
  H.bz = vec_i.z * vec_j.y;

  H.cx = vec_i.x * vec_j.z;
  H.cy = vec_i.y * vec_j.z;
  H.cz = vec_i.z * vec_j.z;

  blas_int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;

  char jobU = 'A';
  char jobVT = 'A';
  double wkopt;

  double s[3], u[9], vt[9];
  std::vector<double> matrix = std::vector<double>{H.ax, H.ay, H.az, H.bx, H.by, H.bz, H.cx, H.cy, H.cz};

  lwork = -1;
  dgesvd_(&jobU, &jobVT, &m, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);

  lwork = static_cast<blas_int>(wkopt);
  std::vector<double> work(static_cast<std::size_t>(lwork));
  dgesvd_(&jobU, &jobVT, &m, &n, matrix.data(), &lda, s, u, &ldu, vt, &ldvt, work.data(), &lwork, &info);

  if (info > 0)
  {
    std::print("The algorithm computing SVD failed to converge.\n");
  }

  double3x3 sigma = double3x3(s[0], 0.0, 0.0, 0.0, s[1], 0.0, 0.0, 0.0, s[2]);
  double3x3 matrix_u = double3x3(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
  double3x3 matrix_vt = double3x3(vt[0], vt[1], vt[2], vt[3], vt[4], vt[5], vt[6], vt[7], vt[8]);

  double3x3 result = matrix_vt.transpose() * matrix_u.transpose();

  double determinant = result.determinant();
  if (determinant < 0.0)
  {
    return matrix_vt.transpose() * double3x3(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, -1.0)) *
           matrix_u.transpose();
  }

  return result;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const double3x3& vec)
{
  archive << vec.ax << vec.ay << vec.az;
  archive << vec.bx << vec.by << vec.bz;
  archive << vec.cx << vec.cy << vec.cz;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, double3x3& vec)
{
  archive >> vec.ax >> vec.ay >> vec.az;
  archive >> vec.bx >> vec.by >> vec.bz;
  archive >> vec.cx >> vec.cy >> vec.cz;
  return archive;
}
