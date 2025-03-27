module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <numbers>
#include <span>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

#if !defined(_WIN32)
#include <assert.h>
#endif

module charge_equilibration_wilmer_snurr;

#ifndef USE_LEGACY_HEADERS
import <cstdint>;
import <vector>;
import <span>;
import <mdspan>;
import <cmath>;
import <iostream>;
import <numbers>;
import <exception>;
#endif

import double3;
import double3x3;
import skelement;
import atom;
import simulationbox;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

extern "C"
{
  void dgesv_(const long long* n, const long long* nrhs, double* a, const long long* lda, long long* ipiv, double* b,
              const long long* ldb, long long* info);
}

static double getJ(const SimulationBox& simulationBox, std::span<Atom> frameworkAtoms, const std::vector<double>& J,
                   size_t i, size_t j, ChargeEquilibration::Type type)
{
  double k = 14.4;      // [Angstroms * electron volts]
  double lambda = 1.2;  // Global hardness scaling parameter
  double eta = 50.0;

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
    case ChargeEquilibration::Type::PeriodicEwaldSum:
    {
      double minCellLength = 250;
      double3x3 box = simulationBox.cell;
      double3x3 inverse_box = simulationBox.inverseCell;
      double volume = simulationBox.volume;
      long long aVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.ax))) - 1;
      long long bVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.by))) - 1;
      long long cVnum = static_cast<long long>(std::ceil(minCellLength / (2.0 * box.cz))) - 1;
      double3 dr;

      if (i == j)
      {
        double orbital = 0;
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
                double orbital_overlap_term = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

                orbital += orbital_overlap_term;
              }
            }
          }
        }

        double alphaStar = 0.0;
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

                alphaStar += std::erfc(r / eta) / r;
              }
            }
          }
        }

        double betaStar = 0;
        double h = 0.0;
        double b = 0.0;
        for (long long u = -aVnum; u <= aVnum; u++)
        {
          for (long long v = -bVnum; v <= bVnum; v++)
          {
            for (long long w = -cVnum; w <= cVnum; w++)
            {
              if (!((u == 0) && (v == 0) && (w == 0)))
              {
                double3 kv;
                kv.x = static_cast<double>(u) * inverse_box.ax + static_cast<double>(v) * inverse_box.bx +
                       static_cast<double>(w) * inverse_box.cx;
                kv.y = static_cast<double>(u) * inverse_box.ay + static_cast<double>(v) * inverse_box.by +
                       static_cast<double>(w) * inverse_box.cy;
                kv.z = static_cast<double>(u) * inverse_box.az + static_cast<double>(v) * inverse_box.bz +
                       static_cast<double>(w) * inverse_box.cz;

                kv.x *= 2.0 * std::numbers::pi;
                kv.y *= 2.0 * std::numbers::pi;
                kv.z *= 2.0 * std::numbers::pi;

                h = std::sqrt(kv.x * kv.x + kv.y * kv.y + kv.z * kv.z);
                b = 0.5 * h * eta;

                betaStar += 1.0 / (h * h) * exp(-b * b);
              }
            }
          }
        }
        betaStar *= 4.0 * std::numbers::pi / volume;

        return J[i] + lambda * (k / 2.0) * (alphaStar + betaStar + orbital - 2.0 / (eta * sqrt(std::numbers::pi)));
      }
      else
      {
        double orbital = 0.0;
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
              double orbital_overlap_term = std::exp(-(a * a * r2)) * (2.0 * a - a * a * r - 1.0 / r);

              orbital += orbital_overlap_term;
            }
          }
        }

        double alpha = 0.0;
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

              alpha += std::erfc(r / eta) / r;
            }
          }
        }

        double beta = 0.0;
        for (long long u = -aVnum; u <= aVnum; u++)
        {
          for (long long v = -bVnum; v <= bVnum; v++)
          {
            for (long long w = -cVnum; w <= cVnum; w++)
            {
              if (!((u == 0) && (v == 0) && (w == 0)))
              {
                double3 kv;
                kv.x = static_cast<double>(u) * inverse_box.ax + static_cast<double>(v) * inverse_box.bx +
                       static_cast<double>(w) * inverse_box.cx;
                kv.y = static_cast<double>(u) * inverse_box.ay + static_cast<double>(v) * inverse_box.by +
                       static_cast<double>(w) * inverse_box.cy;
                kv.z = static_cast<double>(u) * inverse_box.az + static_cast<double>(v) * inverse_box.bz +
                       static_cast<double>(w) * inverse_box.cz;

                kv.x *= 2.0 * std::numbers::pi;
                kv.y *= 2.0 * std::numbers::pi;
                kv.z *= 2.0 * std::numbers::pi;

                double h = std::sqrt(kv.x * kv.x + kv.y * kv.y + kv.z * kv.z);
                double b = 0.5 * h * eta;

                dr.x = (frameworkAtoms[i].position.x - frameworkAtoms[j].position.x);
                dr.y = (frameworkAtoms[i].position.y - frameworkAtoms[j].position.y);
                dr.z = (frameworkAtoms[i].position.z - frameworkAtoms[j].position.z);

                beta += std::cos(kv.x * dr.x + kv.y * dr.y + kv.z * dr.z) / (h * h) * std::exp(-b * b);
              }
            }
          }
        }
        beta *= 4.0 * std::numbers::pi / volume;

        return lambda * (k / 2.0) * (alpha + beta + orbital);
      }
    }
    break;
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
  std::vector<size_t> hydrogen_list{};

  size_t size = frameworkAtoms.size();
  std::vector<double> J(size);
  std::vector<double> X(size);
  std::vector<double> Xc(size);
  std::vector<double> R(size);
  std::vector<double> Q(size);

  for (size_t index = 0; const Atom& atom : frameworkAtoms)
  {
    size_t pseudo_atom_type = static_cast<size_t>(atom.type);
    size_t element_index = forceField.pseudoAtoms[pseudo_atom_type].atomicNumber;

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

  for (size_t i : hydrogen_list)
  {
    X[i] = 0.5 * (13.598 - 2.0);
    J[i] = 13.598 + 2.0;
  }

  std::vector<double> A(size * size);
  std::mdspan<double, std::dextents<size_t, 2>, std::layout_left> As(A.data(), size, size);
  std::vector<double> b(size);
  std::vector<double> x0(size);

  // First row of A is all ones
  for (size_t i = 0; i < size; ++i)
  {
    As[0, i] = 1.0;
  }

  // First element in b is the total charge (FIX!)
  // Can be non-zero for non-periodic structures
  b[0] = 0.0;

  // Rest of elements in b are the differences in electronegativity
  // Use lattice potential. If it hasn't been calculated yet, it will just default to zero.
  for (size_t i = 1; i != size; ++i)
  {
    b[i] = (X[i] - Xc[i]) - (X[i - 1] - Xc[i - 1]);
  }

  // Fill in 2nd to Nth rows of A
  for (size_t i = 1; i != size; ++i)
  {
    for (size_t j = 0; j != size; ++j)
    {
      As[i, j] =
          getJ(simulationBox, frameworkAtoms, J, i - 1, j, type) - getJ(simulationBox, frameworkAtoms, J, i, j, type);
    }
  }

  long long N = static_cast<long long>(size);
  long long nrhs = 1;
  long long lda = static_cast<long long>(size);
  std::vector<long long> ipiv(size);
  long long ldb = static_cast<long long>(size);
  long long info;

  dgesv_(&N, &nrhs, A.data(), &lda, ipiv.data(), b.data(), &ldb, &info);

  if (info > 0)
  {
    throw std::runtime_error(std::format("[charge equilibration]: no solution found'\n"));
  }

  for (size_t i = 0; i < size; ++i)
  {
    frameworkAtoms[i].charge = b[i];
  }
}
