module;

export module dctst_wham;

import std;

/**
 * \brief Binned multistate reweighting (WHAM / binned MBAR) for a one-dimensional reaction coordinate.
 *
 * K states sample the same coordinate under reduced (dimensionless, in units of kT) bias potentials
 * u_k(b) that are constant within every bin b of a fine grid. With M_b the total number of samples that
 * fell into bin b (summed over the states) and n_k the number of samples of state k, the unbiased
 * density Omega_b and the dimensionless free energies g_k of the states obey the self-consistent
 * equations
 *
 *     Omega_b = M_b / sum_k n_k exp(g_k - u_k(b)),      g_k = -ln sum_b Omega_b exp(-u_k(b)),
 *
 * which are solved by direct iteration in log space (gauge g_0 = 0). When the bias depends on the
 * coordinate only, these are exactly the MBAR equations in the limit of vanishing bin width; with the
 * bias evaluated at the bin that is sampled (piecewise constant), they are exact at any bin width.
 *
 * Empty bins and empty states are handled with a finite sentinel instead of -infinity because the
 * library is compiled with -ffast-math.
 */
export struct DcTSTWham
{
  /// Sentinel for ln(0) and the threshold below which a value counts as 'zero'.
  static constexpr double logZero = -1e300;
  static constexpr double logZeroThreshold = -1e299;

  struct Result
  {
    std::vector<double> logDensity;    ///< ln Omega_b per bin (unnormalized); 'logZero' for empty bins.
    std::vector<double> freeEnergies;  ///< g_k per state (gauge g_0 = 0).
    std::size_t iterations{0};
    double residual{0.0};
    bool converged{false};
  };

  /**
   * \param counts           counts[k][b]: samples of state k in bin b.
   * \param reducedPotential reducedPotential[k][b]: bias u_k(b) in units of kT.
   * \param tolerance        convergence criterion on max_k |g_k^{new} - g_k|.
   * \param maximumIterations
   * \param initialFreeEnergies optional starting values of g_k (e.g. from a previous block).
   */
  [[nodiscard]] static Result solve(const std::vector<std::vector<double>>& counts,
                                    const std::vector<std::vector<double>>& reducedPotential, double tolerance,
                                    std::size_t maximumIterations,
                                    const std::vector<double>& initialFreeEnergies = {});
};
