module;

module system;

import <complex>;
import <span>;
import <numbers>;
import <cmath>;
import <vector>;
import <iostream>;
import <algorithm>;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import force_factor;
import running_energy;
import component;
import forcefield;

// TODO:
// An Exact Ewald Summation Method in Theory and Practice
// S. Stenberg and B. Stenqvist
// J. Phys. Chem. A 2020, 124, 3943−3946; https://doi.org/10.1021/acs.jpca.0c01684
//
// Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections to the Ewald potential
// Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks
// J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
//
void System::registerEwaldFourierEnergySingleIon(double3 position, double charge)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
  size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
  size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

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

  CoulombicFourierEnergySingleIon = -Units::CoulombicConversionFactor * ((2.0 * std::numbers::pi / simulationBox.volume) * energy_sum
                       - alpha/std::sqrt(std::numbers::pi));
}


void System::computeEwaldFourierRigidEnergy(const SimulationBox& box, RunningEnergy& energyStatus)
{
    double alpha = forceField.EwaldAlpha;
    double alpha_squared = alpha * alpha;
    double3x3 inv_box = box.inverseUnitCell;
    double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
    double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
    double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

    if (noCharges) return;
    if (omitEwaldFourier) return;

    std::span<const Atom> atoms = spanOfRigidFrameworkAtoms();
    size_t numberOfAtoms = atoms.size();

    size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
    size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
    size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

    std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
    std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
    std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

    if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
    if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
    if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
    if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

    size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
    if (fixedFrameworkStoredEik.size() < numberOfWaveVectors)
    {
      fixedFrameworkStoredEik.resize(numberOfWaveVectors);
    }

    // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
    for (size_t i = 0; i != numberOfAtoms; ++i)
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
    for (size_t kx = 2; kx <= kx_max_unsigned; ++kx)
    {
        for (size_t i = 0; i != numberOfAtoms; ++i)
        {
            eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
        }
    }
    for (size_t ky = 2; ky <= ky_max_unsigned; ++ky)
    {
        for (size_t i = 0; i != numberOfAtoms; ++i)
        {
            eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
        }
    }
    for (size_t kz = 2; kz <= kz_max_unsigned; ++kz)
    {
        for (size_t i = 0; i != numberOfAtoms; ++i)
        {
            eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
        }
    }

    size_t nvec = 0;
    double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / box.volume);
    for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
    {
        double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

        // Only positive kx are used, the negative kx are taken into account by the factor of two
        double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

        for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
        {
            double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

            // Precompute and store eik_x * eik_y outside the kz-loop
            for (size_t i = 0; i != numberOfAtoms; ++i)
            {
                std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<size_t>(std::abs(ky))];
                eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
                eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<size_t>(kx)] * eiky_temp;
            }

            for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
            {
                // Ommit kvec==0
                if ((kx * kx + ky * ky + kz * kz) != 0)
                {
                    double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;

                    std::pair<std::complex<double>, std::complex<double>> cksum(0.0, 0.0);
                    for (size_t i = 0; i != numberOfAtoms; ++i)
                    {
                        std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
                        eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
                        double charge = atoms[i].charge;
                        double scaling = atoms[i].scalingCoulomb;
                        cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
                        cksum.second += charge * (eik_xy[i] * eikz_temp);
                    }

                    double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
                    double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
                    energyStatus.ewald += temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());

                    fixedFrameworkStoredEik[nvec] = cksum;
                    ++nvec;
                }
            }
        }
    }

    // Subtract self-energy
    double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
        double charge = atoms[i].charge;
        double scaling = atoms[i].scalingCoulomb;
        energyStatus.ewald -= prefactor_self * scaling * charge * scaling * charge;
    }
}

void System::computeEwaldFourierEnergy(const SimulationBox &box, RunningEnergy &energyStatus)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = box.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if(noCharges) return;
  if(omitEwaldFourier) return;

  size_t numberOfAtoms = atomPositions.size();

  size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
  size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
  size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if(storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);

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
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / box.volume);
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

          std::pair<std::complex<double>, std::complex<double>> cksum;
          for(size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          energyStatus.ewald += temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          energyStatus.dudlambdaEwald += 2.0 * temp * (cksum.first.real() * cksum.second.real() + cksum.first.imag() * cksum.second.imag());

          storedEik[nvec] = cksum;
          ++nvec;
        }
      }
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for(size_t i = 0; i != numberOfAtoms; ++i)
  {
    double charge = atomPositions[i].charge;
    double scaling = atomPositions[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
    energyStatus.ewald -= prefactor_self * scaling * charge * scaling * charge;
    energyStatus.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  // Subtract exclusion-energy
  for(size_t l = 0; l != components.size(); ++l)
  {
    if(!(components[l].type == Component::Type::Framework && components[l].rigid))
    {
      for(size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = spanOfMolecule(l, m);
        for(size_t i = 0; i != span.size() - 1; i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double3 posA = span[i].position;
          for(size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            bool groupIdB = static_cast<bool>(span[j].groupId);
            double3 posB = span[j].position;
      
            double3 dr = posA - posB;
            dr = box.applyPeriodicBoundaryConditions(dr);
            double r = std::sqrt(double3::dot(dr, dr));
      
            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energyStatus.ewald -= scalingA * scalingB * temp;
            energyStatus.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
          }
        }
      }
    }
  }


  // Handle net-charges
  for(size_t i = 0; i != components.size(); ++i)
  {
    for(size_t j = 0; j != components.size(); ++j)
    {
      //energyStatus.ewald += CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j];
    }
  }
}

// compute gradient
ForceFactor System::computeEwaldFourierGradient()
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  ForceFactor energy{ 0.0, 0.0, 0.0 };

  if (noCharges) return energy;
  if (omitEwaldFourier) return energy;

  size_t numberOfAtoms = atomPositions.size();

  size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
  size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
  size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (size_t i = 0; i != numberOfAtoms; ++i)
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
  for (size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  

  //size_t nvec = 0;
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
      for (size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        // Ommit kvec==0
        if ((kx * kx + ky * ky + kz * kz) != 0)
        {
          double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
          double3 rk = kvec_x + kvec_y + kvec_z;
          double rksq = rk.length_squared();

          std::complex<double> cksum(0.0, 0.0);
          std::complex<double> cksum2(0.0, 0.0);
          for (size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
            cksum += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum2 += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          energy.energy += temp * (cksum.real() * cksum.real() + cksum.imag() * cksum.imag());
          energy.dUdlambda += 2.0 * temp * (cksum.real() * cksum2.real() + cksum.imag() * cksum2.imag());

          for (size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = atomPositions[i].charge;
            double scaling = atomPositions[i].scalingCoulomb;
            atomPositions[i].gradient -= scaling * charge * 2.0 * temp * (cki.imag() * cksum.real() - cki.real() * cksum.imag()) * rk;
          }

          //storedEik[nvec] = cksum;
//          ++nvec;
        }
      }
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (size_t i = 0; i != numberOfAtoms; ++i)
  {
    double charge = atomPositions[i].charge;
    double scaling = atomPositions[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
    energy.energy -= prefactor_self * scaling * charge * scaling * charge;
    energy.dUdlambda -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  // Subtract exclusion-energy
  for (size_t l = 0; l != components.size(); ++l)
  {
    if (!(components[l].type == Component::Type::Framework && components[l].rigid))
    {
      for (size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = spanOfMolecule(l, m);
        for (size_t i = 0; i != span.size() - 1; i++)
        {
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double chargeA = span[i].charge;
          double3 posA = span[i].position;
          for (size_t j = i + 1; j != span.size(); j++)
          {
            double scalingB = span[j].scalingCoulomb;
            bool groupIdB = static_cast<bool>(span[j].groupId);
            double chargeB = span[j].charge;
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double r = std::sqrt(double3::dot(dr, dr));

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energy.energy -= scalingA * scalingB * temp;
            energy.dUdlambda -=  (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
          }
        }
      }
    }
  }


  // Handle net-charges
  for (size_t i = 0; i != components.size(); ++i)
  {
    for (size_t j = 0; j != components.size(); ++j)
    {
      //energy.energy += CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j];
    }
  }

  return energy;
}

// used for volume move
void System::computeEwaldFourierEnergy(const SimulationBox& box, std::span<const Atom> moleculeAtomPositions, RunningEnergy& energyStatus)
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = box.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if (noCharges) return;
  if (omitEwaldFourier) return;

  size_t numberOfAtoms = moleculeAtomPositions.size();

  size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
  size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
  size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (size_t i = 0; i != numberOfAtoms; ++i)
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
  for (size_t kx = 2; kx <= kx_max_unsigned; ++kx)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_x[i + kx * numberOfAtoms] = eik_x[i + (kx - 1) * numberOfAtoms] * eik_x[i + 1 * numberOfAtoms];
    }
  }
  for (size_t ky = 2; ky <= ky_max_unsigned; ++ky)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_y[i + ky * numberOfAtoms] = eik_y[i + (ky - 1) * numberOfAtoms] * eik_y[i + 1 * numberOfAtoms];
    }
  }
  for (size_t kz = 2; kz <= kz_max_unsigned; ++kz)
  {
    for (size_t i = 0; i != numberOfAtoms; ++i)
    {
      eik_z[i + kz * numberOfAtoms] = eik_z[i + (kz - 1) * numberOfAtoms] * eik_z[i + 1 * numberOfAtoms];
    }
  }

  size_t nvec = 0;
  double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / box.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        // Ommit kvec==0
        if ((kx * kx + ky * ky + kz * kz) != 0)
        {
          double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;

          std::pair<std::complex<double>, std::complex<double>> cksum(0.0, 0.0);
          for (size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          energyStatus.ewald += temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          energyStatus.dudlambdaEwald += 2.0 * temp * (cksum.first.real() * cksum.second.real() + cksum.first.imag() * cksum.second.imag());

          totalEik[nvec].first = fixedFrameworkStoredEik[nvec].first + cksum.first;
          totalEik[nvec].second = fixedFrameworkStoredEik[nvec].second + cksum.second;
          ++nvec;
        }
      }
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (size_t i = 0; i != numberOfAtoms; ++i)
  {
    double charge = moleculeAtomPositions[i].charge;
    double scaling = atomPositions[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(atomPositions[i].groupId);
    energyStatus.ewald -= prefactor_self * scaling * charge * scaling * charge;
    energyStatus.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  // Subtract exclusion-energy
  for (size_t l = 0; l != components.size(); ++l)
  {
    if (!(components[l].type == Component::Type::Framework && components[l].rigid))
    {
      for (size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = spanOfMolecule(l, m);
        for (size_t i = 0; i != span.size() - 1; i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          bool groupIdA = static_cast<bool>(span[i].groupId);
          double3 posA = span[i].position;
          for (size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            bool groupIdB = static_cast<bool>(span[j].groupId);
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = box.applyPeriodicBoundaryConditions(dr);
            double r = std::sqrt(double3::dot(dr, dr));

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energyStatus.ewald -= scalingA * scalingB * temp;
            energyStatus.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
          }
        }
      }
    }
  }


  // Handle net-charges
  for (size_t i = 0; i != components.size(); ++i)
  {
    for (size_t j = 0; j != components.size(); ++j)
    {
      //energyStatus.ewald += CoulombicFourierEnergySingleIon * netCharge[i] * netCharge[j];
    }
  }
}


RunningEnergy System::energyDifferenceEwaldFourier(std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedWavevectors, std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  RunningEnergy energy;

  if(noCharges) return energy;
  if(omitEwaldFourier) return energy;

  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double3x3 inv_box = simulationBox.inverseUnitCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);
  size_t numberOfAtoms = newatoms.size() + oldatoms.size();

  size_t kx_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.x);
  size_t ky_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.y);
  size_t kz_max_unsigned = static_cast<size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  if(numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if(numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if(numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if(numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if(storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if(totalEik.size() < numberOfWaveVectors) totalEik.resize(numberOfWaveVectors);

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
  std::pair<std::complex<double>, std::complex<double>> cksum_old;
  std::pair<std::complex<double>, std::complex<double>> cksum_new;
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

          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for(size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
            cksum_old.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_old.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
          for(size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<size_t>(std::abs(kz))];
            eikz_temp.imag(kz>=0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            bool groupIdA = static_cast<bool>(newatoms[i - oldatoms.size()].groupId);
            cksum_new.first += scaling * charge * (eik_xy[i] * eikz_temp);
            cksum_new.second += groupIdA ? charge * eik_xy[i] * eikz_temp : 0.0;
          }

          double rksq = (kvec_x + kvec_y + kvec_z).length_squared();
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          energy.ewald += temp * std::norm(storedWavevectors[nvec].first + cksum_new.first - cksum_old.first);
          energy.ewald -= temp * std::norm(storedWavevectors[nvec].first);

          energy.dudlambdaEwald += 2.0 * temp * ((storedWavevectors[nvec].first + cksum_new.first - cksum_old.first).real() * (storedWavevectors[nvec].second + cksum_new.second - cksum_old.second).real() +
                                            (storedWavevectors[nvec].first + cksum_new.first - cksum_old.first).imag() * (storedWavevectors[nvec].second + cksum_new.second - cksum_old.second).imag());
          energy.dudlambdaEwald -= 2.0 * temp * ((storedWavevectors[nvec].first).real() * (storedWavevectors[nvec].second).real() +
                                            (storedWavevectors[nvec].first).imag() * (storedWavevectors[nvec].second).imag());

          totalEik[nvec].first = storedWavevectors[nvec].first + cksum_new.first - cksum_old.first;
          totalEik[nvec].second = storedWavevectors[nvec].second + cksum_new.second - cksum_old.second;

          ++nvec;
        }
      }
    }
  }

  for(size_t i = 0; i != oldatoms.size(); i++)
  {
    double chargeA = oldatoms[i].charge;
    double scalingA = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    double3 posA = oldatoms[i].position;
    for(size_t j = i + 1; j != oldatoms.size(); j++)
    {
      double chargeB = oldatoms[j].charge;
      double scalingB = oldatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(oldatoms[j].groupId);
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald += scalingA * scalingB * temp;
      energy.dudlambdaEwald += (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

  for(size_t i = 0; i != newatoms.size(); i++)
  {
    double chargeA = newatoms[i].charge;
    double scalingA = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    double3 posA = newatoms[i].position;
    for(size_t j = i + 1; j != newatoms.size(); j++)
    {
      double chargeB = newatoms[j].charge;
      double scalingB = newatoms[j].scalingCoulomb;
      bool groupIdB = static_cast<bool>(newatoms[j].groupId);
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald -= scalingA * scalingB * temp;
      energy.dudlambdaEwald -= (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0);
    }
  }

   // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for(size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(oldatoms[i].groupId);
    energy.ewald += prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald += groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }
  for(size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    bool groupIdA = static_cast<bool>(newatoms[i].groupId);
    energy.ewald -= prefactor_self * scaling * charge * scaling * charge;
    energy.dudlambdaEwald -= groupIdA ? 2.0 * prefactor_self * scaling * charge * charge : 0.0;
  }

  return energy;
}

void System::acceptEwaldMove()
{
  if(noCharges) return;
  if(omitEwaldFourier) return;

  storedEik = totalEik;
}
