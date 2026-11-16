module;

export module mc_moves_acceptance;

import std;

export namespace MC_Moves
{
/// The acceptance probability of a move with and without the bias of its sampling.
struct BiasedAcceptance
{
  double physical;  ///< exp(ln P): the unbiased probability, what the transition matrix records.
  double sampling;  ///< exp(ln P + bias): what the Metropolis test uses.
};

/**
 * \brief The physical and the biased acceptance probability from the logarithm of the first.
 *
 * The bias of a lambda histogram (or of a fugacity ladder) is there to compensate the physical term: the two
 * have opposite signs and can each be far beyond what a double holds as an exponential (|ln P| > 745) while
 * their sum is of order one. Formed as exp(ln P) * exp(bias) the product is 0 * exp(bias) = 0 (a move that
 * is never accepted) or inf * 0 (not a number, which fails the comparison too). The sum is exponentiated
 * instead, and both exponents are held to the range of a double.
 */
inline BiasedAcceptance biasedAcceptance(double logPhysical, double biasTerm)
{
  constexpr double limit = 700.0;
  return {std::exp(std::clamp(logPhysical, -limit, limit)),
          std::exp(std::clamp(logPhysical + biasTerm, -limit, limit))};
}

/// The same for a probability preFactor * exp(logRest); a prefactor of zero is a move that can not be accepted.
inline BiasedAcceptance biasedAcceptance(double preFactor, double logRest, double biasTerm)
{
  if (!(preFactor > 0.0)) return {0.0, 0.0};
  return biasedAcceptance(std::log(preFactor) + logRest, biasTerm);
}

/**
 * \brief The Metropolis test of a probability under a macrostate bias factor: uniform < biasFactor * probability.
 *
 * A bias factor of zero is a refused move (a step out of the window, past the ceiling of an (N, lambda) chain, or
 * an integer swap on such a chain), whatever the probability. Formed as a product, a zero factor times an
 * overflowed probability is not a number, and with -ffast-math the comparison against it is not defined to fail:
 * the refused move could be taken. The factor is looked at first, an overflowed probability is an accepted
 * proposal (the factor is never below exp(-500)), and the product is only formed where it is finite.
 */
inline bool acceptedWithBias(double uniform, double biasFactor, double probability)
{
  if (!(biasFactor > 0.0)) return false;
  if (!(probability > 0.0 && probability < 1.0e300)) return probability >= 1.0e300;
  if (biasFactor >= 1.0 && probability >= 1.0) return true;
  return uniform < biasFactor * probability;
}
}  // namespace MC_Moves
