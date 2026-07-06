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
      long long aVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.ax))) - 1;
      long long bVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.by))) - 1;
      long long cVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.cz))) - 1;
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
        double sigma = 0;
        for (long long u = -aVnum; u <= aVnum; u++)
        {
          for (long long v = -bVnum; v <= bVnum; v++)
          {
            for (long long w = -cVnum; w <= cVnum; w++)
            {
              dr.x = (frameworkAtoms[i].position.x - frameworkAtoms[j].position.x) + static_cast<double>(u) * box.ax +
                     static_cast<double>(v) * box.bx + static_cast<double>(w) * box.cx;
              dr.y = (frameworkAtoms[i].position.y - frameworkAtoms[j].position.y) + static_cast<double>(u) * box.ay +
                     static_cast<double>(v) * box.by + static_cast<double>(w) * box.cy;
              dr.z = (frameworkAtoms[i].position.z - frameworkAtoms[j].position.z) + static_cast<double>(u) * box.az +
                     static_cast<double>(v) * box.bz + static_cast<double>(w) * box.cz;
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

void ChargeEquilibration::computeChargeEquilibration(const ForceField& forceField, const SimulationBox& simulationBox,
                                                     std::span<Atom> frameworkAtoms, ChargeEquilibration::Type type)
{
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
        if (i == j)
        {
          Jab[i * size + i] = J[i] + lambda * (k / 2.0) * Jab[i * size + i];
        }
        else
        {
          // Short-ranged orbital-overlap correction; it decays as exp(-a^2 r^2),
          // so the minimum image is the only image that contributes
          double3 dr = frameworkAtoms[i].position - frameworkAtoms[j].position;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          double r2 = double3::dot(dr, dr);
          double r = std::sqrt(r2);

          double a = std::sqrt(J[i] * J[j]) / k;
          double orbital_overlap_term = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

          double value = lambda * (k / 2.0) * (Jab[i * size + j] + orbital_overlap_term);
          Jab[i * size + j] = value;
          Jab[j * size + i] = value;
        }
      }
    }
  }
  else
  {
    for (std::size_t i = 0; i != size; ++i)
    {
      for (std::size_t j = 0; j != size; ++j)
      {
        Jab[i * size + j] = getJ(simulationBox, frameworkAtoms, J, i, j, type);
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
  for (std::size_t i = 1; i != size; ++i)
  {
    b[i] = (X[i] - Xc[i]) - (X[i - 1] - Xc[i - 1]);
  }

  // Fill in 2nd to Nth rows of A
  for (std::size_t i = 1; i != size; ++i)
  {
    for (std::size_t j = 0; j != size; ++j)
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
