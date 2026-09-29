module;

export module potential_ewald_real_space_table;

import std;

/**
 * \brief Tabulated Ewald real-space Coulomb term u(r^2) = erfc(alpha r) / r.
 *
 * The table is a cubic Hermite spline on a uniform grid in r^2 (not r), so a fully coupled charged pair needs
 * neither erfc, exp, a square root, nor a division: the gradient factor dU/dr / r equals 2 dU/d(r^2), the
 * derivative of the interpolant. The interpolant is C1 and the force is the exact derivative of the interpolated
 * energy, so energy is conserved to rounding regardless of the table error.
 *
 * The node spacing is fixed ('defaultSpacing', 0.00875 Angstrom^2, i.e. 16384 intervals for a 12 Angstrom
 * cutoff): the interpolation error is below 1e-9 relative at r = 1 Angstrom and falls quickly with r (the fourth
 * derivative of u with respect to r^2 scales as r^-9). The range is [1 Angstrom^2, cutoff^2] but is capped at
 * 'defaultMaximumIntervals' nodes (4 MB, about 48 Angstrom), so an unusually large cutoff does not blow up the table
 * or degrade its resolution; pairs beyond the covered range and pairs closer than 1 Angstrom (which do not occur
 * between molecules except in overlapping Monte Carlo trial configurations) fall back to the library erfc.
 *
 * The table is owned by the ForceField (see ForceField::updateEwaldRealSpaceTable) and used by
 * Potentials::potentialCoulomb for the energy and gradient of fully coupled pairs, and by the spatial-decomposition
 * molecular-dynamics pair kernel.
 */
export struct EwaldRealSpaceTable
{
  static constexpr double defaultSpacing{0.00875};               ///< Node spacing in r^2 [Angstrom^2].
  static constexpr std::size_t defaultMaximumIntervals{262144};  ///< Cap on the number of intervals.

  double alpha{0.0};
  double rrMin{1.0};
  double rrMax{0.0};
  double spacing{1.0};
  double inverseSpacing{1.0};
  bool rangeCapped{false};           ///< The requested range exceeded the cap; rrMax is below the cutoff squared.
  std::vector<double> value{};       ///< u at the nodes.
  std::vector<double> derivative{};  ///< du/d(r^2) at the nodes.

  void build(double alphaValue, double rMaxValue, double spacingValue = defaultSpacing,
             std::size_t maximumIntervals = defaultMaximumIntervals)
  {
    alpha = alphaValue;
    rrMin = 1.0;
    spacing = spacingValue;
    inverseSpacing = 1.0 / spacing;
    const double requestedRRMax = std::max(rMaxValue * rMaxValue * 1.002, rrMin + spacing);
    std::size_t intervals = static_cast<std::size_t>(std::ceil((requestedRRMax - rrMin) * inverseSpacing));
    rangeCapped = intervals > maximumIntervals;
    if (rangeCapped) intervals = maximumIntervals;
    rrMax = rrMin + static_cast<double>(intervals) * spacing;
    value.resize(intervals + 2);
    derivative.resize(intervals + 2);
    for (std::size_t k = 0; k < intervals + 2; ++k)
    {
      const double rr = rrMin + static_cast<double>(k) * spacing;
      exact(rr, value[k], derivative[k]);
    }
  }

  /// True when the interpolant covers every distance up to \p rMaxValue (no exact fallback beyond rrMin needed).
  [[nodiscard]] bool spans(double rMaxValue) const { return !value.empty() && rMaxValue * rMaxValue * 1.001 <= rrMax; }

  /// True when the table was built for \p alphaValue and needs no rebuild for the distance \p rMaxValue: it either
  /// spans it or is already at the range cap.
  [[nodiscard]] bool matches(double alphaValue, double rMaxValue) const
  {
    return !value.empty() && alpha == alphaValue && (spans(rMaxValue) || rangeCapped);
  }

  /// True when the table was built for \p alphaValue and the interpolant is valid at the squared distance \p rr.
  [[clang::always_inline]] [[nodiscard]] inline bool covers(double alphaValue, double rr) const
  {
    return alpha == alphaValue && rr >= rrMin && rr < rrMax;
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

  /// u(r^2) and du/d(r^2) by cubic Hermite interpolation; rr must satisfy rrMin <= rr < rrMax (see covers).
  [[clang::always_inline]] inline void interpolate(double rr, double& u, double& dudrr) const
  {
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

  /// u(r^2) and du/d(r^2) by cubic Hermite interpolation (exact fallback below rrMin; rr must be below rrMax).
  [[clang::always_inline]] inline void evaluate(double rr, double& u, double& dudrr) const
  {
    if (rr < rrMin) [[unlikely]]
    {
      exact(rr, u, dudrr);
      return;
    }
    interpolate(rr, u, dudrr);
  }
};
