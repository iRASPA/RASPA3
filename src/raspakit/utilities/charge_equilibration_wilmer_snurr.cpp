module;

#if !defined(_WIN32)
#include <assert.h>
#endif

module charge_equilibration_wilmer_snurr;

import std;

import double3;
import double3x3;
import skelement;
import atom;
import simulationbox;
import interactions_ewald;

static void solveMatrix(std::span<double> matrix, std::size_t rows, std::size_t cols, std::span<double> b,
                        std::span<double> x)
{
  for (std::size_t i = 0; i < rows; ++i)
  {
    x[i] = b[i];
  }

  std::vector<double> d(cols);

  auto element = [matrix, cols](std::size_t row, std::size_t col) -> double& {
    return matrix[row * cols + col];
  };

  for (std::size_t i = 0; i < cols; ++i)
  {
    double r = 0.0;
    for (std::size_t k = i; k < rows; ++k)
    {
      r += element(k, i) * element(k, i);
    }

    if (r == 0.0)
    {
      throw std::runtime_error("[charge equilibration]: matrix is rank deficient\n");
    }

    const double alpha = (element(i, i) < 0.0) ? -std::sqrt(r) : std::sqrt(r);
    const double ak = 1.0 / (r + alpha * element(i, i));

    element(i, i) += alpha;
    d[i] = -alpha;

    for (std::size_t k = i + 1; k < cols; ++k)
    {
      double f = 0.0;
      for (std::size_t j = i; j < rows; ++j)
      {
        f += element(j, k) * element(j, i);
      }
      f *= ak;

      for (std::size_t j = i; j < rows; ++j)
      {
        element(j, k) -= f * element(j, i);
      }
    }

    if (std::abs(alpha) < 1.0e-5)
    {
      throw std::runtime_error("[charge equilibration]: apparent singularity in matrix\n");
    }

    double f = 0.0;
    for (std::size_t j = i; j < rows; ++j)
    {
      f += x[j] * element(j, i);
    }
    f *= ak;

    for (std::size_t j = i; j < rows; ++j)
    {
      x[j] -= f * element(j, i);
    }
  }

  for (std::size_t i = cols; i-- > 0;)
  {
    double sum = 0.0;
    for (std::size_t k = i + 1; k < cols; ++k)
    {
      sum += element(i, k) * x[k];
    }
    x[i] = (x[i] - sum) / d[i];
  }
}

static double getJ(const SimulationBox& simulationBox, std::span<Atom> frameworkAtoms, const std::vector<double>& J,
                   std::size_t i, std::size_t j, ChargeEquilibration::Type type)
{
  double k = 14.4;      // [Angstroms * electron volts]
  double lambda = 1.2;  // Global hardness scaling parameter

  switch (type)
  {
    case ChargeEquilibration::Type::NonPeriodic:
    {
      {
        if (i == j)
        {
          return J[i];  // Return the hardness/idempotential
        }
        else
        {
          double3 dr = frameworkAtoms[i].position - frameworkAtoms[j].position;

          double r2 = double3::dot(dr, dr);
          double r = std::sqrt(r2);

          double a = std::sqrt(J[i] * J[j]) / k;
          double orbital_overlap_term = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

          double Jab = lambda * (k / 2.0) * ((1.0 / r) + orbital_overlap_term);

          return Jab;
        }
      }
    }
    break;
    case ChargeEquilibration::Type::Periodic:
    {
      if (i == j)
      {
        return J[i];  // Return the hardness/idempotential
      }
      else
      {
        double3 dr = frameworkAtoms[i].position - frameworkAtoms[j].position;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);

        double r2 = double3::dot(dr, dr);
        double r = std::sqrt(r2);

        double a = std::sqrt(J[i] * J[j]) / k;
        double orbital_overlap_term = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

        double Jab = lambda * (k / 2.0) * ((1.0 / r) + orbital_overlap_term);

        return Jab;
      }
    }
    break;
    case ChargeEquilibration::Type::PeriodicDirectSum:
    {
      double minCellLength = 250;
      double3x3 box = simulationBox.cell;
      // the replicas along each lattice vector follow from the perpendicular width of the cell (the diagonal of a
      // triclinic cell is not), and at least the neighboring cells are included so a large cell still sees the
      // periodic images
      const double3 widths = simulationBox.perpendicularWidths();
      long long aVnum = std::max(1LL, static_cast<long long>(std::ceil(minCellLength / (2.0 * widths.x))) - 1);
      long long bVnum = std::max(1LL, static_cast<long long>(std::ceil(minCellLength / (2.0 * widths.y))) - 1);
      long long cVnum = std::max(1LL, static_cast<long long>(std::ceil(minCellLength / (2.0 * widths.z))) - 1);
      double3 dr;

      if (i == j)
      {
        double sigmaStar = 0;
        for (long long u = -aVnum; u <= aVnum; u++)
        {
          for (long long v = -bVnum; v <= bVnum; v++)
          {
            for (long long w = -cVnum; w <= cVnum; w++)
            {
              if (!((u == 0) && (v == 0) && (w == 0)))
              {
                dr.x =
                    static_cast<double>(u) * box.ax + static_cast<double>(v) * box.bx + static_cast<double>(w) * box.cx;
                dr.y =
                    static_cast<double>(u) * box.ay + static_cast<double>(v) * box.by + static_cast<double>(w) * box.cy;
                dr.z =
                    static_cast<double>(u) * box.az + static_cast<double>(v) * box.bz + static_cast<double>(w) * box.cz;
                double r2 = double3::dot(dr, dr);
                double r = std::sqrt(r2);

                double a = std::sqrt(J[i] * J[j]) / k;
                double orbitalOverlapTerm = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

                sigmaStar += (1.0 / r) + orbitalOverlapTerm;
              }
            }
          }
        }
        return J[i] + lambda * (k / 2.0) * sigmaStar;
      }
      else
      {
        // the window of images is centred on the minimum image of the pair: centred on the raw separation it shifted
        // by whole cells with the side of the cell boundary the two atoms happened to lie on, so a pair at +0.04a and
        // its symmetry-equivalent at -0.96a summed different images (the truncated sum is not shift-invariant) and
        // got different J_ij, and the charges changed when the structure was merely translated
        const double3 separation =
            simulationBox.applyPeriodicBoundaryConditions(frameworkAtoms[i].position - frameworkAtoms[j].position);
        double sigma = 0;
        for (long long u = -aVnum; u <= aVnum; u++)
        {
          for (long long v = -bVnum; v <= bVnum; v++)
          {
            for (long long w = -cVnum; w <= cVnum; w++)
            {
              dr.x = separation.x + static_cast<double>(u) * box.ax + static_cast<double>(v) * box.bx +
                     static_cast<double>(w) * box.cx;
              dr.y = separation.y + static_cast<double>(u) * box.ay + static_cast<double>(v) * box.by +
                     static_cast<double>(w) * box.cy;
              dr.z = separation.z + static_cast<double>(u) * box.az + static_cast<double>(v) * box.bz +
                     static_cast<double>(w) * box.cz;
              double r2 = double3::dot(dr, dr);
              double r = std::sqrt(r2);

              double a = std::sqrt(J[i] * J[j]) / k;
              double orbitalOverlapTerm = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

              sigma += (1.0 / r) + orbitalOverlapTerm;
            }
          }
        }
        return lambda * (k / 2.0) * sigma;
      }
    }
    break;
    // Type::PeriodicEwaldSum is handled in computeChargeEquilibration using the efficient
    // Ewald routine Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix
    default:
      break;
  }

  return 0.0;
}

// The orbital-overlap correction exp(-a^2 r^2)(2a - a^2 r - 1/r) summed over the lattice images of a separation, as
// Wilmer and Snurr add it to the Coulomb term of every image (for a site and its own images the origin is left out).
// It decays as a Gaussian, so only the images within 6/a, where the Gaussian is below 1e-15, are visited; a cell of a
// few Angstrom has several of them. The refactoring to the Ewald routine had kept the minimum image alone for a pair
// and no image at all for a site's own, which is what EQeq's charges for small cells rest on.
static double orbitalOverlapLatticeSum(const SimulationBox& simulationBox, const double3& separation, double a,
                                       bool skipOrigin)
{
  const double cutoff = 6.0 / a;
  const double3 widths = simulationBox.perpendicularWidths();
  const double reach = cutoff + separation.length();
  const long long nu = static_cast<long long>(std::ceil(reach / widths.x));
  const long long nv = static_cast<long long>(std::ceil(reach / widths.y));
  const long long nw = static_cast<long long>(std::ceil(reach / widths.z));

  double sum = 0.0;
  for (long long u = -nu; u <= nu; ++u)
  {
    for (long long v = -nv; v <= nv; ++v)
    {
      for (long long w = -nw; w <= nw; ++w)
      {
        if (skipOrigin && u == 0 && v == 0 && w == 0) continue;
        const double3 dr = separation + simulationBox.cell * double3(static_cast<double>(u), static_cast<double>(v),
                                                                     static_cast<double>(w));
        const double r2 = double3::dot(dr, dr);
        if (r2 > cutoff * cutoff) continue;
        const double r = std::sqrt(r2);
        sum += std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);
      }
    }
  }
  return sum;
}

void ChargeEquilibration::computeChargeEquilibration(const ForceField& forceField, const SimulationBox& simulationBox,
                                                     std::span<Atom> frameworkAtoms, ChargeEquilibration::Type type)
{
  // nothing to equilibrate: the assembly below writes the total-charge row into b[0] and runs its loops from 1 up to
  // size with !=, which for no atoms is a write past an empty vector and a loop towards SIZE_MAX
  if (frameworkAtoms.empty()) return;

  const double k = 14.4;      // [Angstroms * electron volts]
  const double gamma2 = 0.5;  // Global atomic radii scaling parameter
  std::vector<std::size_t> hydrogen_list{};

  std::size_t size = frameworkAtoms.size();
  std::vector<double> J(size);
  std::vector<double> X(size);
  std::vector<double> Xc(size);
  std::vector<double> R(size);

  for (std::size_t index = 0; const Atom& atom : frameworkAtoms)
  {
    std::size_t pseudo_atom_type = static_cast<std::size_t>(atom.type);
    std::size_t element_index = forceField.pseudoAtoms[pseudo_atom_type].atomicNumber;

    // the reference table holds ionization energies, so it has nothing for an anion (the sign used to be dropped,
    // which read a cation's data for it)
    const std::int64_t signedOxidationState = forceField.pseudoAtoms[pseudo_atom_type].oxidationState;
    if (signedOxidationState < 0)
    {
      throw std::runtime_error(std::format(
          "[charge equilibration]: negative oxidation state {} for pseudo-atom '{}' (atomic number {}); the "
          "reference table holds ionization energies only\n",
          signedOxidationState, forceField.pseudoAtoms[pseudo_atom_type].name, element_index));
    }
    const std::size_t oxidationState = static_cast<std::size_t>(signedOxidationState);

    // the reference table holds the elements up to polonium (Z = 84), with a blank entry at index 0, and room for
    // ten ionization energies per element of which only the measured ones are filled in; the rest are zero. The
    // check on the array size alone let an unmeasured (zero) energy through, which made J negative and the charges
    // NaN (e.g. Ba2+, Hf4+, Ta5+, W6+, Pt2+, Au+, Ir3+).
    if (element_index == 0 || element_index >= chargeEquilibrationElementData.size() ||
        oxidationState >= chargeEquilibrationElementData[element_index].ionizationEnergy.size() ||
        chargeEquilibrationElementData[element_index].ionizationEnergy[oxidationState] == 0.0)
    {
      throw std::runtime_error(std::format(
          "[charge equilibration]: no reference data for pseudo-atom '{}' (atomic number {}, oxidation state {})\n",
          forceField.pseudoAtoms[pseudo_atom_type].name, element_index, oxidationState));
    }

    if (element_index == 1)
    {
      hydrogen_list.push_back(index);
    }

    J[index] = referenceTableJ(forceField, pseudo_atom_type, element_index);
    X[index] = referenceTableX(forceField, pseudo_atom_type, element_index);
    Xc[index] = referenceTableXc(forceField, pseudo_atom_type, element_index);

    R[index] = gamma2 * k * (1.0 / J[index]);

    ++index;
  }

  for (std::size_t i : hydrogen_list)
  {
    X[i] = 0.5 * (13.598 - 2.0);
    J[i] = 13.598 + 2.0;
  }

  // two atoms on one site (or on each other's periodic images) have no Coulomb interaction to speak of: every branch
  // divides by their distance, and the infinities and NaNs that gave went through solveMatrix unnoticed (its checks
  // compare against zero, which a NaN never equals, and under -ffast-math no test for NaN survives), so every charge
  // came out NaN without a word. The distances are checked here, before any is inverted.
  for (std::size_t i = 0; i < size; ++i)
  {
    for (std::size_t j = i + 1; j < size; ++j)
    {
      double3 dr = frameworkAtoms[i].position - frameworkAtoms[j].position;
      if (type != ChargeEquilibration::Type::NonPeriodic) dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      if (double3::dot(dr, dr) < 1e-3 * 1e-3)
      {
        throw std::runtime_error(std::format(
            "[charge equilibration]: atoms {} ('{}') and {} ('{}') are less than 0.001 A apart at ({}, {}, {}); "
            "coincident atoms have no finite interaction\n",
            i, forceField.pseudoAtoms[static_cast<std::size_t>(frameworkAtoms[i].type)].name, j,
            forceField.pseudoAtoms[static_cast<std::size_t>(frameworkAtoms[j].type)].name, frameworkAtoms[i].position.x,
            frameworkAtoms[i].position.y, frameworkAtoms[i].position.z));
      }
    }
  }

  // Build the full matrix of hardness/interaction elements J_ab
  std::vector<double> Jab(size * size);

  if (type == ChargeEquilibration::Type::PeriodicEwaldSum)
  {
    // Efficient Ewald summation: compute the periodic Coulomb potential matrix
    // V_ij = sum over lattice images of 1/|r_i - r_j + L| in a single pass
    std::vector<std::complex<double>> eik_x;
    std::vector<std::complex<double>> eik_y;
    std::vector<std::complex<double>> eik_z;
    std::vector<std::complex<double>> eik_xy;
    Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix(eik_x, eik_y, eik_z, eik_xy, simulationBox,
                                                                        frameworkAtoms, Jab);

    const double lambda = 1.2;  // Global hardness scaling parameter
    for (std::size_t i = 0; i != size; ++i)
    {
      for (std::size_t j = i; j != size; ++j)
      {
        // the Ewald routine gives the Coulomb lattice sum; the orbital-overlap correction of every image is added
        // here (see orbitalOverlapLatticeSum), for a site's own images as well as for a pair's
        if (i == j)
        {
          const double orbital = orbitalOverlapLatticeSum(simulationBox, double3(0.0, 0.0, 0.0), J[i] / k, true);
          Jab[i * size + i] = J[i] + lambda * (k / 2.0) * (Jab[i * size + i] + orbital);
        }
        else
        {
          const double3 dr =
              simulationBox.applyPeriodicBoundaryConditions(frameworkAtoms[i].position - frameworkAtoms[j].position);
          const double orbital = orbitalOverlapLatticeSum(simulationBox, dr, std::sqrt(J[i] * J[j]) / k, false);

          double value = lambda * (k / 2.0) * (Jab[i * size + j] + orbital);
          Jab[i * size + j] = value;
          Jab[j * size + i] = value;
        }
      }
    }
  }
  else
  {
    // J_ij = J_ji, so each pair is summed once (the direct sum visits a quarter of a million images per pair)
    for (std::size_t i = 0; i != size; ++i)
    {
      for (std::size_t j = i; j != size; ++j)
      {
        const double value = getJ(simulationBox, frameworkAtoms, J, i, j, type);
        Jab[i * size + j] = value;
        Jab[j * size + i] = value;
      }
    }
  }

  std::vector<double> A(size * size);
  std::vector<double> b(size);
  std::vector<double> charges(size);

  auto element = [&A, size](std::size_t row, std::size_t col) -> double& { return A[row * size + col]; };

  // First row of A is all ones
  for (std::size_t i = 0; i < size; ++i)
  {
    element(0, i) = 1.0;
  }

  // First element in b is the total charge (FIX!)
  // Can be non-zero for non-periodic structures
  b[0] = 0.0;

  // Rest of elements in b are the differences in electronegativity
  // Use lattice potential. If it hasn't been calculated yet, it will just default to zero.
  for (std::size_t i = 1; i < size; ++i)
  {
    b[i] = (X[i] - Xc[i]) - (X[i - 1] - Xc[i - 1]);
  }

  // Fill in 2nd to Nth rows of A
  for (std::size_t i = 1; i < size; ++i)
  {
    for (std::size_t j = 0; j < size; ++j)
    {
      element(i, j) = Jab[(i - 1) * size + j] - Jab[i * size + j];
    }
  }

  solveMatrix(A, size, size, b, charges);

  for (std::size_t i = 0; i < size; ++i)
  {
    frameworkAtoms[i].charge = charges[i];
  }
}
