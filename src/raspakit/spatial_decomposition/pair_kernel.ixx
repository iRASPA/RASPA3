module;

export module spatial_decomposition_pair_kernel;

import std;

/**
 * \brief Tabulated Ewald real-space Coulomb term u(r^2) = erfc(alpha r) / r for the fast pair kernel.
 *
 * The table is a cubic Hermite spline on a uniform grid in r^2 (not r), so a charged pair needs neither a square
 * root nor a division: the gradient factor dU/dr / r equals 2 dU/d(r^2), the derivative of the interpolant. The
 * interpolant is C1 and the force is the exact derivative of the interpolated energy, so energy is conserved to
 * rounding regardless of the table error. With 16384 intervals over [1, cutoff^2] Angstrom^2 the interpolation
 * error is below 1e-9 relative at r = 1 Angstrom and falls quickly with r (the fourth derivative of u with
 * respect to r^2 scales as r^-9). Pairs closer than 1 Angstrom (which do not occur between molecules in MD) fall
 * back to the library erfc.
 */
export struct EwaldRealSpaceTable
{
  double alpha{0.0};
  double rrMin{1.0};
  double rrMax{0.0};
  double spacing{1.0};
  double inverseSpacing{1.0};
  std::vector<double> value{};       ///< u at the nodes.
  std::vector<double> derivative{};  ///< du/d(r^2) at the nodes.

  void build(double alphaValue, double rMaxValue, std::size_t intervals = 16384)
  {
    alpha = alphaValue;
    rrMin = 1.0;
    rrMax = std::max(rMaxValue * rMaxValue * 1.002, rrMin + 1.0);
    spacing = (rrMax - rrMin) / static_cast<double>(intervals);
    inverseSpacing = 1.0 / spacing;
    value.resize(intervals + 2);
    derivative.resize(intervals + 2);
    for (std::size_t k = 0; k < intervals + 2; ++k)
    {
      const double rr = rrMin + static_cast<double>(k) * spacing;
      exact(rr, value[k], derivative[k]);
    }
  }

  /// u(r^2) and du/d(r^2) from the library functions.
  void exact(double rr, double& u, double& dudrr) const
  {
    const double r = std::sqrt(rr);
    const double inverseR = 1.0 / r;
    const double erfcTerm = std::erfc(alpha * r);
    const double gaussian = std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double>;
    u = erfcTerm * inverseR;
    // du/dr = -(erfc(alpha r) / r^2 + 2 alpha exp(-alpha^2 r^2) / (sqrt(pi) r)), du/d(r^2) = du/dr / (2 r)
    dudrr = -0.5 * (erfcTerm * inverseR * inverseR + 2.0 * alpha * gaussian * inverseR) * inverseR;
  }

  /// u(r^2) and du/d(r^2) by cubic Hermite interpolation (exact fallback below rrMin; rr must be below rrMax).
  [[clang::always_inline]] inline void evaluate(double rr, double& u, double& dudrr) const
  {
    if (rr < rrMin) [[unlikely]]
    {
      exact(rr, u, dudrr);
      return;
    }
    const double position = (rr - rrMin) * inverseSpacing;
    const std::size_t k = static_cast<std::size_t>(position);
    const double t = position - static_cast<double>(k);
    const double p0 = value[k];
    const double p1 = value[k + 1];
    const double m0 = derivative[k] * spacing;
    const double m1 = derivative[k + 1] * spacing;
    const double t2 = t * t;
    const double t3 = t2 * t;
    // Hermite basis: h00 = 2t^3 - 3t^2 + 1, h10 = t^3 - 2t^2 + t, h01 = -2t^3 + 3t^2, h11 = t^3 - t^2
    u = (2.0 * t3 - 3.0 * t2 + 1.0) * p0 + (t3 - 2.0 * t2 + t) * m0 + (-2.0 * t3 + 3.0 * t2) * p1 + (t3 - t2) * m1;
    // d/dt of the basis, divided by the spacing
    dudrr = ((6.0 * t2 - 6.0 * t) * p0 + (3.0 * t2 - 4.0 * t + 1.0) * m0 + (-6.0 * t2 + 6.0 * t) * p1 +
             (3.0 * t2 - 2.0 * t) * m1) *
            inverseSpacing;
  }
};

/// Lennard-Jones parameters of one pair of pseudo-atom types, laid out for the fast kernel.
export struct LennardJonesPair
{
  double epsilon4{0.0};       ///< 4 epsilon [energy units].
  double inverseSigma2{1.0};  ///< 1 / sigma^2.
  double shift{0.0};          ///< Energy shift at the cutoff (0 for truncated potentials).
};
