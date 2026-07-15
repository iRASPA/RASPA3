module;

export module potential_pair_derivatives;

import std;

export namespace Potentials
{
/**
 * \brief Radial derivatives of a pair potential up to a compile-time order.
 *
 * PairDerivatives<Order> is the unified return type of the pair-potential evaluators
 * potentialVDW<Order> and potentialCoulomb<Order>. The struct only carries the fields
 * required for the requested derivative order:
 *   - Order 0: energy and dUdlambda
 *   - Order 1: adds firstDerivativeFactor
 *   - Order 2: adds secondDerivativeFactor
 *
 * Conventions:
 *   - energy is the scaled pair energy scalingA * scalingB * u.
 *   - dUdlambda holds the symmetric derivative factor X such that
 *       dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *   - firstDerivativeFactor is (1/r) dU/dr, so the Cartesian force on atom A is
 *       firstDerivativeFactor * dr with dr = posA - posB.
 *   - secondDerivativeFactor is (U'' - U'/r) / r^2, so the Cartesian Hessian block is
 *       secondDerivativeFactor * (dr ⊗ dr) + firstDerivativeFactor * I.
 */
template <std::size_t Order>
struct PairDerivatives;

template <>
struct PairDerivatives<0>
{
  double energy;
  double dUdlambda;
};

template <>
struct PairDerivatives<1>
{
  double energy;
  double dUdlambda;
  double firstDerivativeFactor;
};

template <>
struct PairDerivatives<2>
{
  double energy;
  double dUdlambda;
  double firstDerivativeFactor;
  double secondDerivativeFactor;
};
}  // namespace Potentials
