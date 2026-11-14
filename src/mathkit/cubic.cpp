module;

module cubic;

import std;

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
  ret = std::abs(std::pow(std::abs(Z), THIRD)) * static_cast<double>(signR(Z));
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

  // relative rounding tolerance of the quantities below; a repeated root makes p, Q or DIS vanish exactly
  // only when the coefficients are dyadic, otherwise they come out a few ulps away from zero
  const double epsilon = 16.0 * std::numeric_limits<double>::epsilon();

  // ====determine the degree of the polynomial ====

  if (A[3] != 0.0)
  {
    // cubic problem: x = t - W turns it into the depressed cubic t^3 + 3 p t - 2 Q = 0
    W = A[2] / A[3] * THIRD;
    const double p = A[1] / A[3] * THIRD - W * W;
    P = p * p * p;
    Q = -.5 * (2.0 * W * W * W - (A[1] * W - A[0]) / A[3]);
    DIS = Q * Q + P;
    // p and Q are differences of terms of these magnitudes, so they carry rounding errors of that size
    const double pScale = std::abs(A[1] / A[3]) * THIRD + W * W;
    const double qScale = 2.0 * std::abs(W * W * W) + (std::abs(A[1] * W) + std::abs(A[0])) / std::abs(A[3]);
    if (std::abs(p) <= epsilon * pScale && std::abs(Q) <= epsilon * qScale)
    {
      // triple root
      X[0] = X[1] = X[2] = -W;
      *L = 3;
    }
    else if (p < 0.0 && std::abs(DIS) <= epsilon * (3.0 * p * p * pScale + 2.0 * std::abs(Q) * qScale))
    {
      // vanishing discriminant: the depressed cubic has the double root d = -cbrt(Q) = -sign(Q) sqrt(-p) and
      // the simple root -2 d (a small p > 0 with a small Q leaves a single real root, DIS = Q^2 + p^3 > 0, which
      // the branch below handles); the trigonometric form
      // below would split the double root by sqrt(|DIS|) through the ill-conditioned acos near +-1
      const double d = std::copysign(std::sqrt(-p), -Q);
      X[0] = std::min(-2.0 * d, d) - W;
      X[1] = d - W;
      X[2] = std::max(-2.0 * d, d) - W;
      *L = 3;
    }
    else if (DIS < 0.0)
    {
      // three real solutions!
      // Confine the argument of ACOS to the interval [-1;1]!
      PHI = std::acos(std::min(1.0, std::max(-1.0, Q / std::sqrt(-P))));
      P = 2.0 * std::sqrt(-p);
      for (i = 0; i < 3; i++) U[i] = P * std::cos((PHI + 2.0 * (static_cast<double>(i)) * PI) * THIRD) - W;
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
    if (DIS >= 0.0)
    {
      // 2 real solutions (equal ones for DIS == 0, a repeated root)
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
    // linear equation A[0] + A[1] x = 0
    X[0] = -A[0] / A[1];
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
    // at a repeated root the derivative is rounding noise (it vanishes exactly only for dyadic coefficients)
    // and the step would be O(1) garbage: skip it there, and keep a step only when it lowers the residual
    const double derivative = A[1] + X[i] * (2.0 * A[2] + X[i] * 3.0 * A[3]);
    const double derivativeScale =
        std::abs(A[1]) + std::abs(X[i]) * (2.0 * std::abs(A[2]) + std::abs(X[i]) * 3.0 * std::abs(A[3]));
    if (std::abs(derivative) <= epsilon * derivativeScale) continue;
    const double residual = A[0] + X[i] * (A[1] + X[i] * (A[2] + X[i] * A[3]));
    const double refined = X[i] - residual / derivative;
    if (std::abs(A[0] + refined * (A[1] + refined * (A[2] + refined * A[3]))) < std::abs(residual)) X[i] = refined;
    //  printf("\n X inside cubic %.15e\n", X[i]);
  }

  return 0;
}
