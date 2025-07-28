module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <print>
#include <span>
#include <type_traits>
#include <vector>
#endif

module interactions_ewald;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import gradient_factor;
import running_energy;
import framework;
import component;
import forcefield;

// TODO:
// An Exact Ewald Summation Method in Theory and Practice
// S. Stenberg and B. Stenqvist
// J. Phys. Chem. A 2020, 124, 3943−3946; https://doi.org/10.1021/acs.jpca.0c01684
//
// Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections
// to the Ewald potential
// Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks
// J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
//
double Interactions::computeEwaldFourierEnergySingleIon(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy, const ForceField &forceField,
    const SimulationBox &simulationBox, double3 position, double charge)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (eik_x.size() < (kx_max_unsigned + 1)) eik_x.resize(kx_max_unsigned + 1);
  if (eik_y.size() < (kx_max_unsigned + 1)) eik_y.resize(ky_max_unsigned + 1);
  if (eik_z.size() < (kx_max_unsigned + 1)) eik_z.resize(kz_max_unsigned + 1);
  if (eik_xy.size() < 1) eik_xy.resize(1);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  eik_x[0] = std::complex<double>(1.0, 0.0);
  eik_y[0] = std::complex<double>(1.0, 0.0);
  eik_z[0] = std::complex<double>(1.0, 0.0);
  double3 s = 2.0 * std::numbers::pi * (inv_box * position);
  eik_x[1] = std::complex<double>(std::cos(s.x), std::sin(s.x));
  eik_y[1] = std::complex<double>(std::cos(s.y), std::sin(s.y));
  eik_z[1] = std::complex<double>(std::cos(s.z), std::sin(s.z));

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    eik_x[kx] = eik_x[kx - 1] * eik_x[1];
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    eik_y[ky] = eik_y[ky - 1] * eik_y[1];
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    eik_z[kz] = eik_z[kz - 1] * eik_z[1];
  }

  double energy_sum = 0.0;
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? 1.0 : 2.0;

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      std::complex<double> eiky_temp = eik_y[static_cast<std::size_t>(std::abs(ky))];
      eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
      eik_xy[0] = eik_x[static_cast<std::size_t>(kx)] * eiky_temp;

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double rksq = (kvec_x + kvec_y + kvec_z).length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          std::complex<double> cksum(0.0, 0.0);
          std::complex<double> eikz_temp = eik_z[static_cast<std::size_t>(std::abs(kz))];
          eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
          cksum += charge * (eik_xy[0] * eikz_temp);

          energy_sum += factor * std::norm(cksum) * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
        }
      }
    }
  }

  return -Units::CoulombicConversionFactor *
         ((2.0 * std::numbers::pi / simulationBox.volume) * energy_sum - alpha / std::sqrt(std::numbers::pi));
}

void Interactions::precomputeEwaldFourierRigid(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> rigidFrameworkAtoms)
{
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  std::size_t numberOfAtoms = rigidFrameworkAtoms.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors)
  {
    fixedFrameworkStoredEik.resize(numberOfWaveVectors);
  }

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * rigidFrameworkAtoms[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double rksq = (kvec_x + kvec_y + kvec_z).length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          std::pair<std::complex<double>, std::complex<double>> cksum(0.0, 0.0);
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = rigidFrameworkAtoms[i].charge;
            double scaling = rigidFrameworkAtoms[i].scalingCoulomb;
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
          }

          fixedFrameworkStoredEik[nvec] = cksum;
          ++nvec;
        }
      }
    }
  }
}

// Energy, called with 'storedEik'
// Volume-move, called with 'totalEik' for 'storedEik'
RunningEnergy Interactions::computeEwaldFourierEnergy(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtomPositions)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  bool omitInterInteractions = forceField.omitInterInteractions;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  RunningEnergy energySum{};

  if (!forceField.useCharge) return energySum;
  if (forceField.omitEwaldFourier) return energySum;

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * moleculeAtomPositions[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double rksq = (kvec_x + kvec_y + kvec_z).length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          std::pair<std::complex<double>, std::complex<double>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(moleculeAtomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::complex<double>> total;
          total.first = rigid.first + cksum.first;
          total.second = rigid.second + cksum.second;

          double rigidEnergy =
              temp * (rigid.first.real() * rigid.first.real() + rigid.first.imag() * rigid.first.imag());

          energySum.ewald_fourier +=
              temp * (total.first.real() * total.first.real() + total.first.imag() * total.first.imag()) - rigidEnergy;

          if (omitInterInteractions)
          {
            energySum.ewald_fourier -=
                temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          }

          energySum.dudlambdaEwald +=
              2.0 * temp * (total.first.real() * total.second.real() + total.first.imag() * total.second.imag());
          energySum.dudlambdaEwald -=
              2.0 * temp * (rigid.first.real() * rigid.second.real() + rigid.first.imag() * rigid.second.imag());

          storedEik[nvec] = total;
          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy
    double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (std::size_t i = 0; i != moleculeAtomPositions.size(); ++i)
    {
      double charge = moleculeAtomPositions[i].charge;
      double scaling = moleculeAtomPositions[i].scalingCoulomb;
      bool groupIdA = static_cast<bool>(moleculeAtomPositions[i].groupId);
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      energySum.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<const Atom> span = std::span(&moleculeAtomPositions[index], size);
        for (std::size_t i = 0; i != span.size() - 1; i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            bool groupIdB = static_cast<bool>(span[j].groupId);
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double r = std::sqrt(double3::dot(dr, dr));

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energySum.ewald_exclusion -= scalingA * scalingB * temp;
            energySum.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
          }
        }
        index += size;
      }
    }
  }

  return energySum;
}

// compute gradient
RunningEnergy Interactions::computeEwaldFourierGradient(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik,
    const std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<Atom> atomPositions)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  bool omitInterInteractions = forceField.omitInterInteractions;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  RunningEnergy energySum{};

  if (!forceField.useCharge) return energySum;
  if (forceField.omitEwaldFourier) return energySum;

  std::size_t numberOfAtoms = atomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * atomPositions[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          std::pair<std::complex<double>, std::complex<double>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::complex<double>> total;
          total.first = rigid.first + cksum.first;
          total.second = rigid.second + cksum.second;

          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          double rigidEnergy =
              temp * (rigid.first.real() * rigid.first.real() + rigid.first.imag() * rigid.first.imag());

          energySum.ewald_fourier +=
              temp * (total.first.real() * total.first.real() + total.first.imag() * total.first.imag()) - rigidEnergy;

          if (forceField.omitInterInteractions)
          {
            energySum.ewald_fourier -=
                temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          }

          energySum.dudlambdaEwald +=
              2.0 * temp * (total.first.real() * total.second.real() + total.first.imag() * total.second.imag());
          energySum.dudlambdaEwald -=
              2.0 * temp * (rigid.first.real() * rigid.second.real() + rigid.first.imag() * rigid.second.imag());

          if (forceField.omitInterInteractions)
          {
            total.first -= cksum.first;
            total.second -= cksum.second;
          }

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            atomPositions[i].gradient -= scaling * charge * 2.0 * temp *
                                         (cki.imag() * total.first.real() - cki.real() * total.first.imag()) * rk;
          }

          totalEik[nvec] = total;
          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy
    double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      double charge = atomPositions[i].charge;
      double scaling = atomPositions[i].scalingCoulomb;
      bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      energySum.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = std::span(&atomPositions[index], size);
        for (std::size_t i = 0; i != span.size() - 1; i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            bool groupIdB = static_cast<bool>(span[j].groupId);
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double rr = double3::dot(dr, dr);
            double r = std::sqrt(rr);

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energySum.ewald_exclusion -= scalingA * scalingB * temp;
            energySum.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);

            temp = Units::CoulombicConversionFactor * (2.0 * std::numbers::inv_sqrtpi) * alpha *
                   std::exp(-(alpha * alpha * r * r)) / rr;
            double Bt0 = -Units::CoulombicConversionFactor * std::erf(alpha * r) / r;
            double Bt1 = temp + Bt0 / rr;
            temp = chargeA * chargeB * Bt1;
            span[i].gradient -= temp * dr;
            span[j].gradient += temp * dr;
          }
        }
        index += size;
      }
    }
  }

  // Handle net-charges
  // for (std::size_t i = 0; i != components.size(); ++i)
  //{
  //  for (std::size_t j = 0; j != components.size(); ++j)
  //  {
  //    //energy.energy += CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j];
  //  }
  //}

  return energySum;
}

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  RunningEnergy energy;

  if (!forceField.useCharge) return energy;
  if (forceField.omitEwaldFourier) return energy;

  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  std::size_t numberOfAtoms = newatoms.size() + oldatoms.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * oldatoms[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }
  for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * newatoms[i - oldatoms.size()].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  std::pair<std::complex<double>, std::complex<double>> cksum_old;
  std::pair<std::complex<double>, std::complex<double>> cksum_new;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double rksq = (kvec_x + kvec_y + kvec_z).length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
            cksum_old.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_old.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            bool groupIdA = static_cast<bool>(newatoms[i - oldatoms.size()].groupId);
            cksum_new.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_new.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          energy.ewald_fourier += temp * std::norm(storedEik[nvec].first + cksum_new.first - cksum_old.first);
          energy.ewald_fourier -= temp * std::norm(storedEik[nvec].first);

          energy.dudlambdaEwald += 2.0 * temp *
                                   ((storedEik[nvec].first + cksum_new.first - cksum_old.first).real() *
                                        (storedEik[nvec].second + cksum_new.second - cksum_old.second).real() +
                                    (storedEik[nvec].first + cksum_new.first - cksum_old.first).imag() *
                                        (storedEik[nvec].second + cksum_new.second - cksum_old.second).imag());
          energy.dudlambdaEwald -= 2.0 * temp *
                                   ((storedEik[nvec].first).real() * (storedEik[nvec].second).real() +
                                    (storedEik[nvec].first).imag() * (storedEik[nvec].second).imag());

          totalEik[nvec].first = storedEik[nvec].first + cksum_new.first - cksum_old.first;
          totalEik[nvec].second = storedEik[nvec].second + cksum_new.second - cksum_old.second;

          ++nvec;
        }
      }
    }
  }

  for (std::size_t i = 0; i != oldatoms.size(); i++)
  {
    double chargeA = oldatoms[i].charge;
    double scalingA = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    double3 posA = oldatoms[i].position;
    for (std::size_t j = i + 1; j != oldatoms.size(); j++)
    {
      double chargeB = oldatoms[j].charge;
      double scalingB = oldatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(oldatoms[j].groupId);
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion += scalingA * scalingB * temp;
      energy.dudlambdaEwald += (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

  for (std::size_t i = 0; i != newatoms.size(); i++)
  {
    double chargeA = newatoms[i].charge;
    double scalingA = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    double3 posA = newatoms[i].position;
    for (std::size_t j = i + 1; j != newatoms.size(); j++)
    {
      double chargeB = newatoms[j].charge;
      double scalingB = newatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(newatoms[j].groupId);
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion -= scalingA * scalingB * temp;
      energy.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    energy.ewald_self += prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald += groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }
  for (std::size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    energy.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  return energy;
}

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  RunningEnergy energy;

  if (!forceField.useCharge) return energy;
  if (forceField.omitEwaldFourier) return energy;

  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  // bool omitInterInteractions = forceField.omitInterInteractions;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  std::size_t numberOfAtoms = newatoms.size() + oldatoms.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * oldatoms[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }
  for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * newatoms[i - oldatoms.size()].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  std::pair<std::complex<double>, std::complex<double>> cksum_old;
  std::pair<std::complex<double>, std::complex<double>> cksum_new;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
            cksum_old.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_old.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
            electricFieldOld[i] -=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            bool groupIdA = static_cast<bool>(newatoms[i - oldatoms.size()].groupId);
            cksum_new.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_new.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
            electricFieldNew[i - oldatoms.size()] +=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          energy.ewald_fourier += temp * std::norm(storedEik[nvec].first + cksum_new.first - cksum_old.first);
          energy.ewald_fourier -= temp * std::norm(storedEik[nvec].first);

          energy.dudlambdaEwald += 2.0 * temp *
                                   ((storedEik[nvec].first + cksum_new.first - cksum_old.first).real() *
                                        (storedEik[nvec].second + cksum_new.second - cksum_old.second).real() +
                                    (storedEik[nvec].first + cksum_new.first - cksum_old.first).imag() *
                                        (storedEik[nvec].second + cksum_new.second - cksum_old.second).imag());
          energy.dudlambdaEwald -= 2.0 * temp *
                                   ((storedEik[nvec].first).real() * (storedEik[nvec].second).real() +
                                    (storedEik[nvec].first).imag() * (storedEik[nvec].second).imag());

          totalEik[nvec].first = storedEik[nvec].first + cksum_new.first - cksum_old.first;
          totalEik[nvec].second = storedEik[nvec].second + cksum_new.second - cksum_old.second;

          ++nvec;
        }
      }
    }
  }

  for (std::size_t i = 0; i != oldatoms.size(); i++)
  {
    double chargeA = oldatoms[i].charge;
    double scalingA = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    double3 posA = oldatoms[i].position;
    for (std::size_t j = i + 1; j != oldatoms.size(); j++)
    {
      double chargeB = oldatoms[j].charge;
      double scalingB = oldatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(oldatoms[j].groupId);
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      double r = std::sqrt(rr);

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion += scalingA * scalingB * temp;
      energy.dudlambdaEwald += (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

  for (std::size_t i = 0; i != newatoms.size(); i++)
  {
    double chargeA = newatoms[i].charge;
    double scalingA = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    double3 posA = newatoms[i].position;
    for (std::size_t j = i + 1; j != newatoms.size(); j++)
    {
      double chargeB = newatoms[j].charge;
      double scalingB = newatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(newatoms[j].groupId);
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      double r = std::sqrt(rr);

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion -= scalingA * scalingB * temp;
      energy.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    energy.ewald_self += prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald += groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  for (std::size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    energy.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  return energy;
}

// Used to compute the difference in electricField for a grown or retraced state
// Used in insertion_CBCMC and deletion_CBCMC

void Interactions::computeEwaldFourierElectricFieldDifference(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  // bool omitInterInteractions = forceField.omitInterInteractions;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  std::size_t numberOfAtoms = newatoms.size() + oldatoms.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * oldatoms[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }
  for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * newatoms[i - oldatoms.size()].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  std::pair<std::complex<double>, std::complex<double>> cksum_old;
  std::pair<std::complex<double>, std::complex<double>> cksum_new;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            electricFieldOld[i] -=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            electricFieldNew[i - oldatoms.size()] +=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          ++nvec;
        }
      }
    }
  }
}

void Interactions::acceptEwaldMove(const ForceField &forceField,
                                   std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
                                   std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik)
{
  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  storedEik = totalEik;
}

std::pair<EnergyStatus, double3x3> Interactions::computeEwaldFourierEnergyStrainDerivative(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::optional<Framework> &framework,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> atomPositions, double UIon, double netChargeFramework,
    std::vector<double> netChargePerComponent) noexcept
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  EnergyStatus energy(1, framework.has_value() ? 1uz : 0uz, components.size());
  double3x3 strainDerivative;

  if (!forceField.useCharge || forceField.omitEwaldFourier) return std::make_pair(energy, strainDerivative);

  std::size_t numberOfAtoms = atomPositions.size();
  std::size_t numberOfComponents = components.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * atomPositions[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  std::vector<std::complex<double>> cksum(numberOfComponents, std::complex<double>(0.0, 0.0));
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          std::complex<double> test{0.0, 0.0};
          std::fill(cksum.begin(), cksum.end(), std::complex<double>(0.0, 0.0));
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::size_t comp = static_cast<std::size_t>(atomPositions[i].componentId);
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            cksum[comp] += scaling * charge * (eik_xy[i] * eikz_temp);
            test += scaling * charge * (eik_xy[i] * eikz_temp);
          }

          test += fixedFrameworkStoredEik[nvec].first;

          for (std::size_t i = 0; i != numberOfComponents; ++i)
          {
            energy.frameworkComponentEnergy(0, i).CoulombicFourier +=
                Potentials::EnergyFactor(2.0 * temp *
                                             (fixedFrameworkStoredEik[nvec].first.real() * cksum[i].real() +
                                              fixedFrameworkStoredEik[nvec].first.imag() * cksum[i].imag()),
                                         0.0);
            for (std::size_t j = 0; j != numberOfComponents; ++j)
            {
              energy.componentEnergy(i, j).CoulombicFourier += Potentials::EnergyFactor(
                  temp * (cksum[i].real() * cksum[j].real() + cksum[i].imag() * cksum[j].imag()), 0.0);
            }
          }

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;

            atomPositions[i].gradient -= scaling * charge * 2.0 * temp *
                                         (cki.imag() * test.real() - cki.real() * test.imag()) *
                                         (kvec_x + kvec_y + kvec_z);
          }

          double currentEnergy = temp * (test.real() * test.real() + test.imag() * test.imag());
          double fac = 2.0 * (1.0 / rksq + 0.25 / (alpha * alpha)) * currentEnergy;
          strainDerivative.ax -= currentEnergy - fac * rk.x * rk.x;
          strainDerivative.bx -= -fac * rk.x * rk.y;
          strainDerivative.cx -= -fac * rk.x * rk.z;

          strainDerivative.ay -= -fac * rk.y * rk.x;
          strainDerivative.by -= currentEnergy - fac * rk.y * rk.y;
          strainDerivative.cy -= -fac * rk.y * rk.z;

          strainDerivative.az -= -fac * rk.z * rk.x;
          strainDerivative.bz -= -fac * rk.z * rk.y;
          strainDerivative.cz -= currentEnergy - fac * rk.z * rk.z;

          ++nvec;
        }
      }
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (std::size_t i = 0; i != atomPositions.size(); ++i)
  {
    double charge = atomPositions[i].charge;
    double scaling = atomPositions[i].scalingCoulomb;
    std::size_t comp = static_cast<std::size_t>(atomPositions[i].componentId);
    energy.componentEnergy(comp, comp).CoulombicFourier -=
        Potentials::EnergyFactor(prefactor_self * scaling * charge * scaling * charge, 0.0);
  }

  // Subtract exclusion-energy
  std::size_t index{0};
  for (std::size_t l = 0; l != components.size(); ++l)
  {
    std::size_t size = components[l].atoms.size();
    for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
    {
      std::span<Atom> span = std::span(&atomPositions[index], size);
      for (std::size_t i = 0; i != span.size() - 1; i++)
      {
        double chargeA = span[i].charge;
        double scalingA = span[i].scalingCoulomb;
        // bool groupIdA = static_cast<bool>(span[i].groupId);
        double3 posA = span[i].position;
        for (std::size_t j = i + 1; j != span.size(); j++)
        {
          double chargeB = span[j].charge;
          double scalingB = span[j].scalingCoulomb;
          // bool groupIdB = static_cast<bool>(span[j].groupId);
          double3 posB = span[j].position;

          double3 dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          double rr = double3::dot(dr, dr);
          double r = std::sqrt(rr);

          energy.componentEnergy(l, l).CoulombicFourier -= Potentials::EnergyFactor(
              Units::CoulombicConversionFactor * scalingA * chargeA * scalingB * chargeB * std::erf(alpha * r) / r,
              0.0);

          double temp = Units::CoulombicConversionFactor * (2.0 * std::numbers::inv_sqrtpi) * alpha *
                        std::exp(-(alpha * alpha * r * r)) / rr;
          double Bt0 = -Units::CoulombicConversionFactor * std::erf(alpha * r) / r;
          double Bt1 = temp + Bt0 / rr;
          temp = chargeA * chargeB * Bt1;
          span[i].gradient -= temp * dr;
          span[j].gradient += temp * dr;

          strainDerivative.ax -= temp * dr.x * dr.x;
          strainDerivative.bx -= temp * dr.y * dr.x;
          strainDerivative.cx -= temp * dr.z * dr.x;

          strainDerivative.ay -= temp * dr.x * dr.y;
          strainDerivative.by -= temp * dr.y * dr.y;
          strainDerivative.cy -= temp * dr.z * dr.y;

          strainDerivative.az -= temp * dr.x * dr.z;
          strainDerivative.bz -= temp * dr.y * dr.z;
          strainDerivative.cz -= temp * dr.z * dr.z;
        }
      }
      index += size;
    }
  }

  // Handle net-charges
  for (std::size_t i = 0; i != components.size(); ++i)
  {
    energy.frameworkComponentEnergy(0, i).CoulombicFourier +=
        Potentials::EnergyFactor(2.0 * UIon * netChargeFramework * netChargePerComponent[i], 0.0);
  }

  for (std::size_t i = 0; i != components.size(); ++i)
  {
    for (std::size_t j = 0; j != components.size(); ++j)
    {
      energy.componentEnergy(i, j).CoulombicFourier +=
          Potentials::EnergyFactor(UIon * netChargePerComponent[i] * netChargePerComponent[j], 0.0);
    }
  }

  return std::make_pair(energy, strainDerivative);
}

void Interactions::computeEwaldFourierElectrostaticPotential(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::span<double> electricPotentialMolecules, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> moleculeAtomPositions)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  bool omitInterInteractions = forceField.omitInterInteractions;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  RunningEnergy energySum{};

  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  // std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * moleculeAtomPositions[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          std::pair<std::complex<double>, std::complex<double>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(moleculeAtomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::complex<double>> total;
          total.first = rigid.first + cksum.first;
          total.second = rigid.second + cksum.second;

          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          if (forceField.omitInterInteractions)
          {
            total.first -= cksum.first;
            total.second -= cksum.second;
          }

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            electricPotentialMolecules[i] +=
                2.0 * temp * (cki.real() * total.first.real() + cki.imag() * total.first.imag());
          }

          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy
    double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (std::size_t i = 0; i != moleculeAtomPositions.size(); ++i)
    {
      double charge = moleculeAtomPositions[i].charge;
      double scaling = moleculeAtomPositions[i].scalingCoulomb;
      electricPotentialMolecules[i] -= 2.0 * prefactor_self * scaling * charge;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<const Atom> span = std::span(&moleculeAtomPositions[index], size);
        std::span<double> electricPotential = std::span(&electricPotentialMolecules[index], size);
        for (std::size_t i = 0; i != span.size(); i++)
        {
          double3 posA = span[i].position;
          for (std::size_t j = 0; j != span.size(); j++)
          {
            if (i != j)
            {
              double chargeB = span[j].charge;
              double scalingB = span[j].scalingCoulomb;
              double3 posB = span[j].position;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double rr = double3::dot(dr, dr);
              double r = std::sqrt(rr);

              electricPotential[i] -= Units::CoulombicConversionFactor * scalingB * chargeB * std::erf(alpha * r) / r;
            }
          }
        }
        index += size;
      }
    }
  }
}

RunningEnergy Interactions::computeEwaldFourierElectricField(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldMolecules,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> moleculeAtomPositions)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  bool omitInterInteractions = forceField.omitInterInteractions;
  bool omitInterPolarization = forceField.omitInterPolarization;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  RunningEnergy energySum{};

  if (!forceField.useCharge) return energySum;
  if (forceField.omitEwaldFourier) return energySum;

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * moleculeAtomPositions[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }

  // Calculate remaining positive kx, ky and kz by recurrence
  for (std::size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (std::size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  std::size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        double3 rk = kvec_x + kvec_y + kvec_z;
        double rksq = rk.length_squared();

        // Ommit kvec==0
        std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          std::pair<std::complex<double>, std::complex<double>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(moleculeAtomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          std::pair<std::complex<double>, std::complex<double>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::complex<double>> total = rigid;
          // if (!omitInterInteractions || !omitInterPolarization)
          //{
          total.first += cksum.first;
          total.second += cksum.second;
          //}

          energySum.ewald_fourier +=
              temp * (total.first.real() * total.first.real() + total.first.imag() * total.first.imag());

          energySum.ewald_fourier -=
              temp * (rigid.first.real() * rigid.first.real() + rigid.first.imag() * rigid.first.imag());

          if (omitInterInteractions)
          {
            energySum.ewald_fourier -=
                temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          }

          energySum.dudlambdaEwald +=
              2.0 * temp * (total.first.real() * total.second.real() + total.first.imag() * total.second.imag());
          energySum.dudlambdaEwald -=
              2.0 * temp * (rigid.first.real() * rigid.second.real() + rigid.first.imag() * rigid.second.imag());

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;

            electricFieldMolecules[i] +=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          storedEik[nvec] = total;
          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy
    double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (std::size_t i = 0; i != moleculeAtomPositions.size(); ++i)
    {
      double charge = moleculeAtomPositions[i].charge;
      double scaling = moleculeAtomPositions[i].scalingCoulomb;
      bool groupIdA = static_cast<bool>(moleculeAtomPositions[i].groupId);
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      energySum.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = std::span(&moleculeAtomPositions[index], size);
        std::span<double3> electricField = std::span(&electricFieldMolecules[index], size);
        for (std::size_t i = 0; i != span.size(); i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            if (i != j)
            {
              double chargeB = span[j].charge;
              double scalingB = span[j].scalingCoulomb;
              bool groupIdB = static_cast<bool>(span[j].groupId);
              double3 posB = span[j].position;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double rr = double3::dot(dr, dr);
              double r = std::sqrt(rr);

              if (!omitInterInteractions)
              {
                double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
                energySum.ewald_exclusion -= scalingA * scalingB * temp;
                energySum.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
              }

              double temp = Units::CoulombicConversionFactor * (2.0 * std::numbers::inv_sqrtpi) * alpha *
                            std::exp(-(alpha * alpha * r * r)) / rr;
              double Bt0 = -Units::CoulombicConversionFactor * std::erf(alpha * r) / r;
              double Bt1 = temp + Bt0 / rr;
              if (!omitInterPolarization)
              {
                electricField[i] += scalingB * chargeB * Bt1 * dr;
                electricField[j] -= scalingA * chargeA * Bt1 * dr;
              }
            }
          }
        }
        index += size;
      }
    }
  }

  return energySum;
}
