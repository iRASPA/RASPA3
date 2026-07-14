module;

export module interactions_ewald_kvector;

import std;

import double3;
import double3x3;
import atom;

export namespace Interactions::Ewald
{
/**
 * \brief Builds the per-atom phase-factor tables exp(i k.r) used by all Ewald Fourier routines.
 *
 * Fills eik_x/eik_y/eik_z with exp(i 2 pi k s) for k = 0..kMax along each axis, where s is the
 * fractional coordinate of the atom. The tables are laid out with atom-major stride:
 * eik_x[i + kx * numberOfAtoms]. Values for k = 0 and k = 1 are computed explicitly; higher k
 * follow by the multiplicative recurrence. The eik_xy scratch vector is sized for later use by
 * fillEikXYRow. All vectors are grown when too small and never shrunk (callers reuse them as
 * persistent scratch space).
 *
 * \param position Callable returning the Cartesian position (double3) of atom i.
 */
template <typename PositionGetter>
inline void buildEikTables(std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
                           std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
                           std::size_t numberOfAtoms, std::size_t kxMax, std::size_t kyMax, std::size_t kzMax,
                           const double3x3& inverseCell, PositionGetter&& position)
{
  if (numberOfAtoms * (kxMax + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kxMax + 1));
  if (numberOfAtoms * (kyMax + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (kyMax + 1));
  if (numberOfAtoms * (kzMax + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kzMax + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inverseCell * position(i));
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kxMax; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= kyMax; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kzMax; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }
}

/**
 * \brief Convenience overload of buildEikTables for a span of atoms.
 */
inline void buildEikTables(std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
                           std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
                           std::span<const Atom> atoms, std::size_t kxMax, std::size_t kyMax, std::size_t kzMax,
                           const double3x3& inverseCell)
{
  buildEikTables(eik_x, eik_y, eik_z, eik_xy, atoms.size(), kxMax, kyMax, kzMax, inverseCell,
                 [atoms](std::size_t i) { return atoms[i].position; });
}

/**
 * \brief Convenience overload of buildEikTables for the concatenation of two atom spans.
 *
 * Used by the Monte Carlo difference routines: atoms 0..old-1 are the old positions and atoms
 * old..old+new-1 are the new positions, matching the index convention of the k-loop bodies.
 */
inline void buildEikTables(std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
                           std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
                           std::span<const Atom> oldatoms, std::span<const Atom> newatoms, std::size_t kxMax,
                           std::size_t kyMax, std::size_t kzMax, const double3x3& inverseCell)
{
  buildEikTables(eik_x, eik_y, eik_z, eik_xy, oldatoms.size() + newatoms.size(), kxMax, kyMax, kzMax, inverseCell,
                 [oldatoms, newatoms](std::size_t i)
                 { return i < oldatoms.size() ? oldatoms[i].position : newatoms[i - oldatoms.size()].position; });
}

/**
 * \brief Precomputes eik_x * eik_y for the current (kx, ky) into eik_xy, outside the kz-loop.
 *
 * Negative ky is handled by conjugating the stored positive-ky phase factor.
 */
inline void fillEikXYRow(std::vector<std::complex<double>>& eik_xy, const std::vector<std::complex<double>>& eik_x,
                         const std::vector<std::complex<double>>& eik_y, std::size_t numberOfAtoms,
                         std::make_signed_t<std::size_t> kx, std::make_signed_t<std::size_t> ky)
{
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
    eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
    eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
  }
}

/**
 * \brief Returns the full phase factor exp(i k.r_i) for atom i and the current wave vector.
 *
 * Combines the (kx, ky) row stored in eik_xy by fillEikXYRow with the kz phase; negative kz is
 * handled by conjugating the stored positive-kz factor.
 */
[[clang::always_inline]] inline std::complex<double> eikPhase(const std::vector<std::complex<double>>& eik_z,
                                                              const std::vector<std::complex<double>>& eik_xy,
                                                              std::size_t numberOfAtoms, std::size_t i,
                                                              std::make_signed_t<std::size_t> kz)
{
  std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
  eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
  return eik_xy[i] * eikz_temp;
}
}  // namespace Interactions::Ewald
