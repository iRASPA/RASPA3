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

// TODO:
// An Exact Ewald Summation Method in Theory and Practice
// S. Stenberg and B. Stenqvist
// J. Phys. Chem. A 2020, 124, 3943−3946
//
// Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections to the Ewald potential
// Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks
// J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
//
void System::registerEwaldFourierEnergySingleIon(double3 position, double charge)
{
  double alpha_squared = simulationBox.alpha * simulationBox.alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if(eik_x.size() < (kx_max_unsigned + 1)) eik_x.resize(kx_max_unsigned + 1);
  if(eik_y.size() < (kx_max_unsigned + 1)) eik_y.resize(ky_max_unsigned + 1);
  if(eik_z.size() < (kx_max_unsigned + 1)) eik_z.resize(kz_max_unsigned + 1);
  if(eik_xy.size() < 1 ) eik_xy.resize(1);


  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  eik_x[0] = std::complex<double>(1.0, 0.0);
  eik_y[0] = std::complex<double>(1.0, 0.0);
  eik_z[0] = std::complex<double>(1.0, 0.0);
  double3 s = 2.0 * std::numbers::pi * (inv_box * position);
  eik_x[1] = std::complex<double>(std::cos(s.x), std::sin(s.x));
  eik_y[1] = std::complex<double>(std::cos(s.y), std::sin(s.y));
  eik_z[1] = std::complex<double>(std::cos(s.z), std::sin(s.z));

  // Calculate remaining positive kx, ky and kz by recurrence
  for(size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    eik_x[kx] = eik_x[kx - 1] * eik_x[1];
  }
  for(size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    eik_y[ky] = eik_y[ky - 1] * eik_y[1];
  }
  for(size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    eik_z[kz] = eik_z[kz - 1] * eik_z[1];
  }

  double energy_sum = 0.0;
  for(std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? 1.0 : 2.0;

    for(std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      std::complex<double> eiky_temp = eik_y[static_cast<size_t>(std::abs(ky))];
      eiky_temp.imag(ky>=0 ? eiky_temp.imag() : -eiky_temp.imag());
      eik_xy[0] = eik_x[static_cast<size_t>(kx)] * eiky_temp;

      for(std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        // Ommit kvec==0
        if((kx * kx + ky * ky + kz * kz) != 0)
        {
          double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;

          std::complex<double> cksum(0.0, 0.0);
          std::complex<double> eikz_temp = eik_z[static_cast<size_t>(std::abs(kz))];
          eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
          cksum += charge * (eik_xy[0] * eikz_temp);

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          energy_sum += factor*std::norm(cksum) * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
        }
      }
    }
  }

  CoulombicFourierEnergySingleIon = (2.0 * std::numbers::pi / simulationBox.volume) * energy_sum; 
}

double System::energy(std::span<Atom> atoms, const SimulationBox &box)
{
  double alpha_squared = box.alpha * box.alpha;
  double3x3 inv_box = box.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  size_t numberOfAtoms = atoms.size();

  if(noCharges) return 0.0;
  if(omitEwaldFourier) return 0.0;

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * (ky_max_unsigned + 1) * (kz_max_unsigned + 1);
  if(storedEik.size() < components.size() * numberOfWaveVectors) storedEik.resize(components.size() * numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for(size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * atoms[i].position);
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

  double energy_sum = 0.0;
  for(std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? 1.0 : 2.0;

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

          std::complex<double> cksum(0.0, 0.0);
          for(size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            cksum += atoms[i].charge * (eik_xy[i] * eikz_temp);
          }

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          energy_sum += factor*std::norm(cksum) * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
        }
      }
    }
  }

  return (2.0 * std::numbers::pi / box.volume) * energy_sum; 
}

void System::computeEwaldFourierEnergy()
{
  double alpha = simulationBox.alpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if(noCharges) return;
  if(omitEwaldFourier) return;

  size_t numberOfAtoms = atomPositions.size();
  size_t numberOfComponents = components.size();

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if(storedEik.size() < components.size() * numberOfWaveVectors) storedEik.resize(components.size() * numberOfWaveVectors);

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

  size_t nvec = 0;
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
              runningEnergies(i,j).CoulombicFourier += EnergyFactor(temp * (cksum[i].real() * cksum[j].real() + cksum[i].imag() * cksum[j].imag()), 0.0);
            }
          }

          for(size_t i = 0; i != numberOfComponents; ++i)
          {
            storedEik[i + nvec * numberOfComponents] = cksum[i];
          }
          ++nvec;
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
    runningEnergies(comp,comp).CoulombicFourier -= EnergyFactor(prefactor_self * scaling * charge * scaling * charge, 0.0);
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

              runningEnergies(l,l).CoulombicFourier -= EnergyFactor(Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
            }
          }
      }
  }


  // Handle net-charges
  for(size_t i = 0; i != numberOfComponents; ++i)
  {
    for(size_t j = 0; j != numberOfComponents; ++j)
    {
      runningEnergies(i,j).CoulombicFourier += EnergyFactor(CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j], 0.0);
    }
  }
}



EnergyStatus System::energyDifferenceEwaldFourier(std::vector<std::complex<double>> &storedWavevectors, std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  if(noCharges) return EnergyStatus(components.size());
  if(omitEwaldFourier) return EnergyStatus(components.size());

  double alpha = simulationBox.alpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  size_t numberOfAtoms = newatoms.size() + oldatoms.size();

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if(storedEik.size() < components.size() * numberOfWaveVectors) storedEik.resize(components.size() * numberOfWaveVectors);
  if(totalEik.size() < components.size() * numberOfWaveVectors) totalEik.resize(components.size() * numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for(size_t i = 0; i != oldatoms.size(); ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    double3 s = 2.0 * std::numbers::pi * (inv_box * oldatoms[i].position);
    eik_x[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.x), std::sin(s.x));
    eik_y[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.y), std::sin(s.y));
    eik_z[i + 1 * numberOfAtoms] = std::complex<double>(std::cos(s.z), std::sin(s.z));
  }
  for(size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
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

  size_t nvec = 0;
  size_t numberOfComponents = components.size();
  EnergyStatus energies{numberOfComponents};
  std::vector<std::complex<double>> cksum_old(numberOfComponents, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> cksum_new(numberOfComponents, std::complex<double>(0.0, 0.0));
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
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

          std::fill(cksum_old.begin(), cksum_old.end(), std::complex<double>(0.0, 0.0));
          for(size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            size_t comp = static_cast<size_t>(oldatoms[i].componentId);
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            cksum_old[comp] += scaling * charge * (eik_xy[i] * eikz_temp);
          }


          std::fill(cksum_new.begin(), cksum_new.end(), std::complex<double>(0.0, 0.0));
          for(size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            size_t comp = static_cast<size_t>(newatoms[i - oldatoms.size()].componentId);
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            cksum_new[comp] += scaling * charge * (eik_xy[i] * eikz_temp);
          }

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          for(size_t i = 0; i != numberOfComponents; ++i)
          {
            for(size_t j = 0; j != numberOfComponents; ++j)
            {
              energies(i,j).CoulombicFourier += EnergyFactor(temp *
                  ((storedWavevectors[i + nvec * numberOfComponents] + cksum_new[i] - cksum_old[i]).real() *
                   (storedWavevectors[j + nvec * numberOfComponents] + cksum_new[j] - cksum_old[j]).real() +
                   (storedWavevectors[i + nvec * numberOfComponents] + cksum_new[i] - cksum_old[i]).imag() *
                   (storedWavevectors[j + nvec * numberOfComponents] + cksum_new[j] - cksum_old[j]).imag()), 0.0);
              energies(i,j).CoulombicFourier -= EnergyFactor(temp *
                  ((storedWavevectors[i + nvec * numberOfComponents]).real() *
                   (storedWavevectors[j + nvec * numberOfComponents]).real() +
                   (storedWavevectors[i + nvec * numberOfComponents]).imag() *
                   (storedWavevectors[j + nvec * numberOfComponents]).imag()), 0.0);
            }
          }

          for(size_t i = 0; i != numberOfComponents; ++i)
          {
            totalEik[i + nvec * numberOfComponents] = storedWavevectors[i + nvec * numberOfComponents] + cksum_new[i] - cksum_old[i];
          }

          ++nvec;
        }
      }
    }
  }

  for(size_t i = 0; i != oldatoms.size(); i++)
  {
    size_t compA = static_cast<size_t>(oldatoms[i].componentId);
    double factorA = oldatoms[i].scalingCoulomb * oldatoms[i].charge;
    double3 posA = oldatoms[i].position;
    for(size_t j = i + 1; j != oldatoms.size(); j++)
    {
      size_t compB = static_cast<size_t>(oldatoms[j].componentId);
      double factorB = oldatoms[j].scalingCoulomb * oldatoms[j].charge;
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      energies(compA,compB).CoulombicFourier += EnergyFactor(0.5 * Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
      energies(compB,compA).CoulombicFourier += EnergyFactor(0.5 * Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
    }
  }

  for(size_t i = 0; i != newatoms.size(); i++)
  {
    size_t compA = static_cast<size_t>(newatoms[i].componentId);
    double factorA = newatoms[i].scalingCoulomb * newatoms[i].charge;
    double3 posA = newatoms[i].position;
    for(size_t j = i + 1; j != newatoms.size(); j++)
    {
      size_t compB = static_cast<size_t>(newatoms[j].componentId);
      double factorB = newatoms[j].scalingCoulomb * newatoms[j].charge;
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      energies(compA,compB).CoulombicFourier -= EnergyFactor(0.5 * Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
      energies(compB,compA).CoulombicFourier -= EnergyFactor(0.5 * Units::CoulombicConversionFactor * factorA * factorB * std::erf(alpha * r) / r, 0.0);
    }
  }

   // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * simulationBox.alpha / std::sqrt(std::numbers::pi);
  for(size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    size_t comp = static_cast<size_t>(oldatoms[i].componentId);
    energies(comp,comp).CoulombicFourier += EnergyFactor(prefactor_self * scaling * scaling * charge * charge, 0.0);
  }
  for(size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    size_t comp = static_cast<size_t>(newatoms[i].componentId);
    energies(comp,comp).CoulombicFourier -= EnergyFactor(prefactor_self * scaling * scaling * charge * charge, 0.0);
  }

  energies.sumTotal();
  return energies;
}

void System::acceptEwaldMove()
{
    if(noCharges) return;
    if(omitEwaldFourier) return;

    storedEik = totalEik;
}
