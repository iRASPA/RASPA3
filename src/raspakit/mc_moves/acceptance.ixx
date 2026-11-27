module;

export module mc_moves_acceptance;

import std;

export namespace MC_Moves
{
/**
 * \brief exp(x) with the exponent held to the range of a double.
 *
 * The acceptance rules compare a uniform number with exp(ln P). Beyond |ln P| of about 709 the exponential is
 * infinite or zero; the code is built with -ffast-math, under which a comparison against an infinity is not
 * defined, and a product with a second factor (inf * 0) is not a number. The value is held finite: above 700 it
 * is exp(700), which is accepted whatever the uniform number, and below -700 it is zero, which is rejected.
 */
inline double boundedExp(double exponent)
{
  constexpr double limit = 700.0;
  if (!(exponent > -limit)) return 0.0;
  return std::exp(std::min(exponent, limit));
}

/**
 * \brief The logarithm of the prefactor of a swap move: beta f V / (N + 1) of an insertion from a reservoir of
 * fugacity f, N / (beta f V) of a deletion into it.
 *
 * A reservoir of no fugacity (no external pressure, or a mole fraction of zero: a desorption into vacuum, a component
 * loaded beforehand) gives no insertion and accepts every deletion. Its logarithm is held far beyond the range of the
 * exponential ('boundedExp'), where the logarithm of 0 and the division by 0 gave -inf, +inf or NaN, which the
 * fast-math flags of the build assume never occur (and NaN fails every comparison of the Metropolis test).
 */
inline double logSwapPrefactor(double numerator, double denominator)
{
  constexpr double beyondRange = 1.0e300;
  if (!(numerator > 0.0)) return -beyondRange;
  if (!(denominator > 0.0)) return beyondRange;
  return std::log(numerator / denominator);
}

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
 *
 * Held by 'boundedExp', so that what is below the range is zero and not exp(-700). The physical probability
 * goes into the collection matrix of a transition-matrix run, which takes a link for supported once anything
 * above zero has been recorded on both sides of it: a move that can not happen then counted as one that had,
 * and ln Pi took a step of 700 across the link where it is to be carried flat.
 */
inline BiasedAcceptance biasedAcceptance(double logPhysical, double biasTerm)
{
  return {boundedExp(logPhysical), boundedExp(logPhysical + biasTerm)};
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
