module;

module system;

import <complex>;
import <span>;
import <numbers>;
import <cmath>;
import <vector>;
import <iostream>;
import <algorithm>;

import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;

std::pair<EnergyStatus, double3x3> System::computeEwaldFourierEnergyStrainDerivative() noexcept
{
  double alpha = simulationBox.alpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  EnergyStatus energy(components.size());
  double3x3 strainDerivative;

  if(noCharges || omitEwaldFourier) return std::make_pair(energy, strainDerivative);

  size_t numberOfAtoms = atomPositions.size();
  size_t numberOfComponents = components.size();

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  //size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for(size_t i = 0; i != numberOfAtoms; ++i)
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
  for(size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for(size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for(size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for(size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for(size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for(size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  std::vector<std::complex<double>> cksum(numberOfComponents, std::complex<double>(0.0, 0.0));
  for(std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for(std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for(size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<size_t>(std::abs(ky))];
        eiky_temp.imag(ky>=0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<size_t>(kx)] * eiky_temp;
      }

      for(std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        // Ommit kvec==0
        if((kx * kx + ky * ky + kz * kz) != 0)
        {
          double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;

          std::fill(cksum.begin(), cksum.end(), std::complex<double>(0.0, 0.0));
          for(size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            size_t comp = static_cast<size_t>(atomPositions[i].componentId);
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            cksum[comp] += scaling * charge * (eik_xy[i] * eikz_temp);
          }


          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          for(size_t i = 0; i != numberOfComponents; ++i)
          {
            for(size_t j = 0; j != numberOfComponents; ++j)
            {
              energy(i,j).CoulombicFourier += EnergyFactor(temp * (cksum[i].real() * cksum[j].real() + cksum[i].imag() * cksum[j].imag()), 0.0);
            }
          }

        }
      }
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * simulationBox.alpha / std::sqrt(std::numbers::pi);
  for(size_t i = 0; i != numberOfAtoms; ++i)
  {
    double charge = atomPositions[i].charge;
    double scaling = atomPositions[i].scalingCoulomb;
    size_t comp = static_cast<size_t>(atomPositions[i].componentId);
    energy(comp,comp).CoulombicFourier -= EnergyFactor(prefactor_self * scaling * charge * scaling * charge, 0.0);
  }

  // Subtract exclusion-energy
  for(size_t l = 0; l != components.size(); ++l)
  {
      for(size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
          std::span<Atom> span = spanOfMolecule(l, m);
          for(size_t i = 0; i != span.size() - 1; i++)
          {
            double factorA = span[i].scalingCoulomb * span[i].charge;
            double3 posA = span[i].position;
            for(size_t j = i + 1; j != span.size(); j++)
            {
              double factorB = span[j].scalingCoulomb * span[j].charge;
              double3 posB = span[j].position;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double r = std::sqrt(double3::dot(dr, dr));

              energy(l,l).CoulombicFourier -= EnergyFactor(Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
            }
          }
      }
  }


  // Handle net-charges
  for(size_t i = 0; i != numberOfComponents; ++i)
  {
    for(size_t j = 0; j != numberOfComponents; ++j)
    {
      energy(i,j).CoulombicFourier += EnergyFactor(CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j], 0.0);
    }
  }

  return std::make_pair(energy, strainDerivative);
}
