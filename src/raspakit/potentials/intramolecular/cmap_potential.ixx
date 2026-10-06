module;

export module cmap_potential;

import std;

import archive;
import double3;
import double3x3;

/**
 * \brief A CMAP correction map: a periodic two-dimensional grid of energies over two dihedral angles,
 * interpolated with a bicubic spline (the CHARMM / AMBER ff19SB / OpenMM CMAPTorsionForce scheme).
 *
 * The grid has `resolution` points per dimension with spacing 360/resolution degrees. The energy of grid point
 * (i, j) is `energies[i * resolution + j]` and belongs to the angles
 *   phi_i = -180 + i * 360 / resolution,   psi_j = -180 + j * 360 / resolution   [degrees]
 * (the AMBER prmtop CMAP_PARAMETER layout: the first angle is the slow index). The energies are in Kelvin.
 *
 * The map is made C^2-smooth along the grid lines with periodic cubic splines, from which the gradients and the
 * cross derivative at the grid points follow; every cell then carries the 16 coefficients of its bicubic patch,
 * which reproduces the values, gradients and cross derivatives at the corners and is C^1 across cells (the
 * construction of OpenMM's CMAPTorsionForceImpl).
 */
export struct CMAPMap
{
  std::uint64_t versionNumber{1};

  std::string name{};
  std::size_t resolution{0};
  std::vector<double> energies{};  ///< resolution * resolution grid energies [K], first angle is the slow index

  /// Per cell (i, j), stored at i * resolution + j: the coefficients c[k * 4 + l] of
  /// sum_{k,l} c_{kl} t^k u^l with t, u in [0, 1) the fractional cell coordinates of the two angles.
  std::vector<std::array<double, 16>> coefficients{};

  struct Evaluation
  {
    double energy{};
    double dPhi{};     ///< dU/dphi [K/rad]
    double dPsi{};     ///< dU/dpsi [K/rad]
    double dPhiPhi{};  ///< d2U/dphi2
    double dPhiPsi{};  ///< d2U/dphi dpsi
    double dPsiPsi{};  ///< d2U/dpsi2
  };

  CMAPMap() = default;
  CMAPMap(std::string name, std::size_t resolution, std::vector<double> energies);

  bool operator==(const CMAPMap &other) const
  {
    return name == other.name && resolution == other.resolution && energies == other.energies;
  }

  /// (Re)computes the bicubic patch coefficients from `energies`.
  void buildSpline();

  /// The interpolated energy and its derivatives at the angles phi, psi [radians, any range].
  Evaluation evaluate(double phi, double psi) const;

  /// Multiplies the map energies (and the spline coefficients) by `factor`.
  void scaleEnergy(double factor);

  std::string print() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CMAPMap &m);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CMAPMap &m);
};

/**
 * \brief A CMAP correction term over five consecutive atoms A-B-C-D-E: the map energy at the dihedral angles
 * phi = (A, B, C, D) and psi = (B, C, D, E).
 *
 * The dihedral angles follow the IUPAC (AMBER / CHARMM / OpenMM) sign convention. The map is referenced by
 * its index in the `cmapMaps` list of the owning `IntraMolecularPotentials`.
 */
export struct CMAPPotential
{
  std::uint64_t versionNumber{1};

  std::array<std::size_t, 5> identifiers{};  ///< the atoms A, B, C, D, E of the term
  std::size_t mapIndex{0};                   ///< the map of the term

  CMAPPotential() = default;
  CMAPPotential(std::array<std::size_t, 5> identifiers, std::size_t mapIndex)
      : identifiers(identifiers), mapIndex(mapIndex)
  {
  }

  bool operator==(const CMAPPotential &) const = default;

  std::string print() const;

  /**
   * \brief The IUPAC dihedral angle of A-B-C-D and its Cartesian gradient.
   *
   * phi = atan2(|b2| b1.(b2 x b3), (b1 x b2).(b2 x b3)) with b1 = B - A, b2 = C - B, b3 = D - C.
   */
  static std::pair<double, std::array<double3, 4>> dihedralAngleAndGradient(const double3 &posA, const double3 &posB,
                                                                              const double3 &posC,
                                                                              const double3 &posD);

  /// The IUPAC dihedral angle of A-B-C-D.
  static double dihedralAngle(const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD);

  double calculateEnergy(const CMAPMap &map, const double3 &posA, const double3 &posB, const double3 &posC,
                         const double3 &posD, const double3 &posE) const;

  /// Energy, gradient on the five atoms and the strain derivative (sum_i (r_i - r_B) (x) g_i).
  std::tuple<double, std::array<double3, 5>, double3x3> potentialEnergyGradientStrain(
      const CMAPMap &map, const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD,
      const double3 &posE) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CMAPPotential &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CMAPPotential &p);
};
