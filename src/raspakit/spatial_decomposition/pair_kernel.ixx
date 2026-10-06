module;

export module spatial_decomposition_pair_kernel;

import std;

// The tabulated Ewald real-space term (EwaldRealSpaceTable) lives in the potentials layer, shared with
// Potentials::potentialCoulomb; re-exported here for the fast kernel.
export import potential_ewald_real_space_table;

/// Lennard-Jones parameters of one pair of pseudo-atom types, laid out for the fast kernel.
export struct LennardJonesPair
{
  double epsilon4{0.0};       ///< 4 epsilon [energy units].
  double inverseSigma2{1.0};  ///< 1 / sigma^2.
  double shift{0.0};          ///< Energy shift at the cutoff (0 for truncated potentials; for the force-switched
                              ///< form the constant offset of the energy below the switching distance).
};

/// The global switching of the Lennard-Jones pairs in the fast kernels (ForceField::TruncationMethod Switched:
/// the quintic potential switch; ForceSwitched: the CHARMM force switch), applied on [distance, cutoff]. The
/// constants derive from the switching distance r_s and the cutoff rc.
export struct LennardJonesSwitching
{
  std::uint32_t mode{0};       ///< 0: no switching, 1: potential switch, 2: force switch.
  double distance{0.0};        ///< r_s.
  double distanceSquared{0.0};  ///< r_s^2.
  double inverseWidth{0.0};    ///< 1 / (rc - r_s) (potential switch).
  double inverseCutOff3{0.0};  ///< rc^-3 (force switch).
  double a12{1.0};             ///< rc^6 / (rc^6 - r_s^6) (force switch).
  double a6{1.0};              ///< rc^3 / (rc^3 - r_s^3) (force switch).

  static LennardJonesSwitching make(std::uint32_t mode, double distance, double cutOff)
  {
    LennardJonesSwitching s{};
    s.mode = mode;
    if (mode == 0) return s;
    s.distance = distance;
    s.distanceSquared = distance * distance;
    s.inverseWidth = 1.0 / (cutOff - distance);
    s.inverseCutOff3 = 1.0 / (cutOff * cutOff * cutOff);
    const double q = (distance / cutOff) * (distance / cutOff) * (distance / cutOff);
    s.a12 = 1.0 / (1.0 - q * q);
    s.a6 = 1.0 / (1.0 - q);
    return s;
  }
};

/// Replaces the plain Lennard-Jones energy 'u' (unshifted) and gradient factor 'f' ((1/r) dU/dr) of a pair in
/// the switching region (rr > distanceSquared, below the cutoff) by the switched values.
export inline void applyLennardJonesSwitch(const LennardJonesSwitching& s, double epsilon4, double sigma6, double rr,
                                           double& u, double& f)
{
  const double r = std::sqrt(rr);
  if (s.mode == 1)
  {
    const double x = (r - s.distance) * s.inverseWidth;
    const double x2 = x * x;
    const double oneMinusX = 1.0 - x;
    const double sw = 1.0 + x2 * x * (-10.0 + x * (15.0 - 6.0 * x));
    const double dsw = -30.0 * x2 * oneMinusX * oneMinusX * s.inverseWidth;
    f = f * sw + u * dsw / r;
    u *= sw;
    return;
  }
  const double c6 = epsilon4 * sigma6;
  const double c12 = c6 * sigma6;
  const double invR3 = 1.0 / (r * rr);
  const double invR4 = invR3 / r;
  const double d3 = invR3 - s.inverseCutOff3;
  const double d6 = invR3 * invR3 - s.inverseCutOff3 * s.inverseCutOff3;
  u = c12 * s.a12 * d6 * d6 - c6 * s.a6 * d3 * d3;
  f = (-12.0 * c12 * s.a12 * d6 * invR3 * invR4 + 6.0 * c6 * s.a6 * d3 * invR4) / r;
}

/// Force on one atom of the compact local (owned + ghost image) array of a sub-domain; always double.
export struct LocalForce
{
  double x{0.0};
  double y{0.0};
  double z{0.0};
};
