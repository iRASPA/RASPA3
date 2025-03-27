module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#endif

module cubic;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <numeric>;
import <algorithm>;
#endif

int signR(double Z)
{
  int ret;

  if (Z > 0.0)
    ret = 1;
  else if (Z < 0.0)
    ret = -1;
  else
    ret = 0;

  return ret;
}

double CBRT(double Z)
{
  double ret;
  const double THIRD = 1. / 3.;
  // define cubic root as statement function
  // SIGN has different meanings in both C and Fortran
  //  Was unable to use the sign command of C, so wrote my own
  //  that why a new variable needs to be introduced that keeps track of the sign of
  //  SIGN is supposed to return a 1, -1 or 0 depending on what the sign of the argument is
  ret = std::abs(pow(std::abs(Z), THIRD)) * static_cast<double>(signR(Z));
  return ret;
}

int cubic(double A[4], double X[3], int* L)
{
  const double PI = 3.1415926535897932;
  const double THIRD = 1. / 3.;
  double U[3], W, P, Q, DIS, PHI;
  int i;

  // define cubic root as statement function
  //  In C, the function is defined outside of the cubic fct

  // ====determine the degree of the polynomial ====

  if (A[3] != 0.0)
  {
    // cubic problem
    W = A[2] / A[3] * THIRD;
    P = pow((A[1] / A[3] * THIRD - pow(W, 2)), 3);
    Q = -.5 * (2.0 * pow(W, 3) - (A[1] * W - A[0]) / A[3]);
    DIS = pow(Q, 2) + P;
    if (DIS < 0.0)
    {
      // three real solutions!
      // Confine the argument of ACOS to the interval [-1;1]!
      PHI = std::acos(std::min(1.0, std::max(-1.0, Q / std::sqrt(-P))));
      P = 2.0 * std::pow((-P), (5.e-1 * THIRD));
      for (i = 0; i < 3; i++) U[i] = P * cos((PHI + 2.0 * (static_cast<double>(i)) * PI) * THIRD) - W;
      X[0] = std::min(U[0], std::min(U[1], U[2]));
      X[1] = std::max(std::min(U[0], U[1]), std::max(std::min(U[0], U[2]), std::min(U[1], U[2])));
      X[2] = std::max(U[0], std::max(U[1], U[2]));
      *L = 3;
    }
    else
    {
      // only one real solution!
      DIS = std::sqrt(DIS);
      X[0] = CBRT(Q + DIS) + CBRT(Q - DIS) - W;
      *L = 1;
    }
  }
  else if (A[2] != 0.0)
  {
    // quadratic problem
    P = 0.5 * A[1] / A[2];
    DIS = std::pow(P, 2) - A[0] / A[2];
    if (DIS > 0.0)
    {
      // 2 real solutions
      X[0] = -P - std::sqrt(DIS);
      X[1] = -P + std::sqrt(DIS);
      *L = 2;
    }
    else
    {
      // no real solution
      *L = 0;
    }
  }
  else if (A[1] != 0.0)
  {
    // linear equation
    X[0] = A[0] / A[1];
    *L = 1;
  }
  else
  {
    // no equation
    *L = 0;
  }
  /*
   *     ==== perform one step of a newton iteration in order to minimize
   *          round-off errors ====
   */
  for (i = 0; i < *L; i++)
  {
    X[i] =
        X[i] - (A[0] + X[i] * (A[1] + X[i] * (A[2] + X[i] * A[3]))) / (A[1] + X[i] * (2.0 * A[2] + X[i] * 3.0 * A[3]));
    //  printf("\n X inside cubic %.15e\n", X[i]);
  }

  return 0;
}
