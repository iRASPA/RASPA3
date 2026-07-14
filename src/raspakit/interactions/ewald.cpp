module;

module interactions_ewald;

import std;

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
import coulomb_potential;
import forcefield;

// Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections
// to the Ewald potential
// Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks
// J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
//

// Per-group structure-factor derivative sums: element g holds sum_i charge_i * exp(ik.r_i) over the
// atoms tagged with dU/dlambda group id g+1 (Atom::groupId, 0 means untracked).
using GroupComplexSums = std::array<std::complex<double>, maximumNumberOfDUDlambdaGroups>;

static inline GroupComplexSums operator+(const GroupComplexSums &a, const GroupComplexSums &b)
{
  GroupComplexSums result;
  for (std::size_t g = 0; g != a.size(); ++g) result[g] = a[g] + b[g];
  return result;
}

static inline GroupComplexSums operator-(const GroupComplexSums &a, const GroupComplexSums &b)
{
  GroupComplexSums result;
  for (std::size_t g = 0; g != a.size(); ++g) result[g] = a[g] - b[g];
  return result;
}

static inline GroupComplexSums &operator+=(GroupComplexSums &a, const GroupComplexSums &b)
{
  for (std::size_t g = 0; g != a.size(); ++g) a[g] += b[g];
  return a;
}

static inline GroupComplexSums &operator-=(GroupComplexSums &a, const GroupComplexSums &b)
{
  for (std::size_t g = 0; g != a.size(); ++g) a[g] -= b[g];
  return a;
}

// Accumulates the Fourier-space dU/dlambda contribution factor * Re(sk * conj(dsk[g])) for each group g.
static inline void addFourierDUdlambda(RunningEnergy &energy, double factor, const std::complex<double> &sk,
                                       const GroupComplexSums &dsk)
{
  for (std::size_t g = 0; g != dsk.size(); ++g)
  {
    energy.dudlambdaEwald[g] += factor * (sk.real() * dsk[g].real() + sk.imag() * dsk[g].imag());
  }
}

// Net-charge correction difference for a Monte Carlo move that replaces 'oldatoms' by 'newatoms';
// see Bogusz et al., J. Chem. Phys. 108, 7070 (1998). 'netCharge' is the total net charge of the
// system (framework plus adsorbates) before the move, 'netChargeDerivativeExternal' is the per-group
// charge of group-tagged atoms outside 'oldatoms'/'newatoms' (nonzero when other dU/dlambda-tagged
// molecules exist, e.g. the partner molecule in chained pair moves), and 'singleIonFourierSum' is
// the Fourier sum of a single unit charge for the current box and wave vectors.
static void addNetChargeCorrectionDifference(
    RunningEnergy &energy, double singleIonFourierSum, double alpha, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups> &netChargeDerivativeExternal,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  double deltaCharge = 0.0;
  std::array<double, maximumNumberOfDUDlambdaGroups> chargeDerivativeNew = netChargeDerivativeExternal;
  std::array<double, maximumNumberOfDUDlambdaGroups> chargeDerivativeOld = netChargeDerivativeExternal;
  for (const Atom &atom : oldatoms)
  {
    deltaCharge -= atom.scalingCoulomb * atom.charge;
    if (atom.groupId != 0) chargeDerivativeOld[atom.groupId - 1] += atom.charge;
  }
  for (const Atom &atom : newatoms)
  {
    deltaCharge += atom.scalingCoulomb * atom.charge;
    if (atom.groupId != 0) chargeDerivativeNew[atom.groupId - 1] += atom.charge;
  }

  double uIon = -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
  double netChargeNew = netCharge + deltaCharge;
  energy.ewald_fourier += uIon * (netChargeNew * netChargeNew - netCharge * netCharge);
  for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
  {
    energy.dudlambdaEwald[g] += 2.0 * uIon * (netChargeNew * chargeDerivativeNew[g] - netCharge * chargeDerivativeOld[g]);
  }
}

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
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
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
          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum{};
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
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtomPositions,
    double netChargeFramework)
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
  double singleIonFourierSum = 0.0;

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
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

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
          singleIonFourierSum += temp;

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            std::uint8_t groupIdA = moleculeAtomPositions[i].groupId;
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> total;
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

          addFourierDUdlambda(energySum, 2.0 * temp, total.first, total.second);
          addFourierDUdlambda(energySum, -2.0 * temp, rigid.first, rigid.second);

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
      std::uint8_t groupIdA = moleculeAtomPositions[i].groupId;
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      if (groupIdA != 0) energySum.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
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
          std::uint8_t groupIdA = span[i].groupId;
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            std::uint8_t groupIdB = span[j].groupId;
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double r = std::sqrt(double3::dot(dr, dr));

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energySum.ewald_exclusion -= scalingA * scalingB * temp;
            energySum.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, -temp);
          }
        }
        index += size;
      }
    }
  }

  // Net-charge correction: neutralizing background plus removal of the spurious interaction of the
  // net charge with its own periodic images; see Bogusz et al., J. Chem. Phys. 108, 7070 (1998).
  // The framework-framework part is excluded, consistent with the subtraction of the rigid-framework
  // Fourier energy above.
  {
    double netChargeAdsorbates = 0.0;
    std::array<double, maximumNumberOfDUDlambdaGroups> netChargeDerivative{};
    for (std::size_t i = 0; i != moleculeAtomPositions.size(); ++i)
    {
      netChargeAdsorbates += moleculeAtomPositions[i].scalingCoulomb * moleculeAtomPositions[i].charge;
      if (moleculeAtomPositions[i].groupId != 0) netChargeDerivative[moleculeAtomPositions[i].groupId - 1] += moleculeAtomPositions[i].charge;
    }
    double uIon =
        -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
    if (omitInterInteractions)
    {
      energySum.ewald_fourier += 2.0 * uIon * netChargeFramework * netChargeAdsorbates;
      for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
      {
        energySum.dudlambdaEwald[g] += 2.0 * uIon * netChargeFramework * netChargeDerivative[g];
      }
    }
    else
    {
      energySum.ewald_fourier += uIon * (2.0 * netChargeFramework + netChargeAdsorbates) * netChargeAdsorbates;
      for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
      {
        energySum.dudlambdaEwald[g] +=
            2.0 * uIon * (netChargeFramework + netChargeAdsorbates) * netChargeDerivative[g];
      }
    }
  }

  return energySum;
}

// compute gradient
RunningEnergy Interactions::computeEwaldFourierGradient(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &totalEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> atomData,
    std::span<AtomDynamics> atomDynamics, double netChargeFramework, const std::optional<Framework>& framework,
    std::span<const Atom> frameworkAtoms, std::span<AtomDynamics> frameworkDynamics)
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
  double singleIonFourierSum = 0.0;

  if (!forceField.useCharge) return energySum;
  if (forceField.omitEwaldFourier) return energySum;

  const bool flexibleFramework =
      framework && !framework->rigid && frameworkDynamics.size() == frameworkAtoms.size();
  const std::size_t frameworkOffset = flexibleFramework ? frameworkAtoms.size() : 0;
  std::vector<Atom> liveAtoms;
  if (flexibleFramework)
  {
    liveAtoms.reserve(frameworkAtoms.size() + atomData.size());
    liveAtoms.insert(liveAtoms.end(), frameworkAtoms.begin(), frameworkAtoms.end());
    liveAtoms.insert(liveAtoms.end(), atomData.begin(), atomData.end());
  }
  const std::span<const Atom> atoms = flexibleFramework ? std::span<const Atom>(liveAtoms) : atomData;
  const std::size_t numberOfAtoms = atoms.size();
  const auto addGradient = [&](std::size_t atom, const double3& gradient)
  {
    if (flexibleFramework && atom < frameworkOffset)
    {
      frameworkDynamics[atom].gradient += gradient;
    }
    else
    {
      atomDynamics[atom - frameworkOffset].gradient += gradient;
    }
  };

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
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
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
          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = atoms[i].charge;
            double scaling = atoms[i].scalingCoulomb;
            std::uint8_t groupIdA = atoms[i].groupId;
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid{};
          if (!flexibleFramework) rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> total;
          total.first = rigid.first + cksum.first;
          total.second = rigid.second + cksum.second;

          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          singleIonFourierSum += temp;

          double rigidEnergy =
              temp * (rigid.first.real() * rigid.first.real() + rigid.first.imag() * rigid.first.imag());

          energySum.ewald_fourier +=
              temp * (total.first.real() * total.first.real() + total.first.imag() * total.first.imag()) - rigidEnergy;

          if (forceField.omitInterInteractions)
          {
            energySum.ewald_fourier -=
                temp * (cksum.first.real() * cksum.first.real() + cksum.first.imag() * cksum.first.imag());
          }

          addFourierDUdlambda(energySum, 2.0 * temp, total.first, total.second);
          addFourierDUdlambda(energySum, -2.0 * temp, rigid.first, rigid.second);

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
            double charge = atoms[i].charge;
            double scaling = atoms[i].scalingCoulomb;
            addGradient(i, -scaling * charge * 2.0 * temp *
                               (cki.imag() * total.first.real() - cki.real() * total.first.imag()) * rk);
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
      double charge = atoms[i].charge;
      double scaling = atoms[i].scalingCoulomb;
      std::uint8_t groupIdA = atoms[i].groupId;
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      if (groupIdA != 0) energySum.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<const Atom> span = std::span(&atomData[index], size);
        std::span<AtomDynamics> spanDynamics = std::span(&atomDynamics[index], size);
        for (std::size_t i = 0; i != span.size() - 1; i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          std::uint8_t groupIdA = span[i].groupId;
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            double chargeB = span[j].charge;
            double scalingB = span[j].scalingCoulomb;
            std::uint8_t groupIdB = span[j].groupId;
            double3 posB = span[j].position;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double rr = double3::dot(dr, dr);
            double r = std::sqrt(rr);

            double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
            energySum.ewald_exclusion -= scalingA * scalingB * temp;
            energySum.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, -temp);

            temp = Units::CoulombicConversionFactor * (2.0 * std::numbers::inv_sqrtpi) * alpha *
                   std::exp(-(alpha * alpha * r * r)) / rr;
            double Bt0 = -Units::CoulombicConversionFactor * std::erf(alpha * r) / r;
            double Bt1 = temp + Bt0 / rr;
            temp = chargeA * chargeB * Bt1;
            spanDynamics[i].gradient -= temp * dr;
            spanDynamics[j].gradient += temp * dr;
          }
        }
        index += size;
      }
    }

    if (flexibleFramework)
    {
      std::set<std::array<std::size_t, 2>> excludedPairs;
      std::map<std::array<std::size_t, 2>, double> coulombScaling;
      for (const CoulombPotential& potential : framework->intraMolecularPotentials.coulombs)
      {
        coulombScaling[{std::min(potential.identifiers[0], potential.identifiers[1]),
                        std::max(potential.identifiers[0], potential.identifiers[1])}] = potential.scaling;
      }
      const auto excludeIfAbsentOrScaled = [&](const std::array<std::size_t, 2>& pair)
      {
        const auto scaling = coulombScaling.find(pair);
        if (scaling == coulombScaling.end() || scaling->second != 1.0) excludedPairs.insert(pair);
      };
      if (!framework->connectivityTable.table.empty())
      {
        for (const std::array<std::size_t, 2>& bond : framework->connectivityTable.findAllBonds())
        {
          excludeIfAbsentOrScaled({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
        }
        for (const std::array<std::size_t, 3>& bend : framework->connectivityTable.findAllBends())
        {
          excludeIfAbsentOrScaled({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
        }
        for (const std::array<std::size_t, 4>& torsion : framework->connectivityTable.findAllTorsions())
        {
          excludeIfAbsentOrScaled({std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])});
        }
      }
      for (const std::array<std::size_t, 2>& pair : excludedPairs)
      {
        const std::size_t i = pair[0];
        const std::size_t j = pair[1];
        const double3 dr =
            simulationBox.applyPeriodicBoundaryConditions(frameworkAtoms[i].position - frameworkAtoms[j].position);
        const double rr = double3::dot(dr, dr);
        const double r = std::sqrt(rr);
        const double chargeProduct = Units::CoulombicConversionFactor * frameworkAtoms[i].scalingCoulomb *
                                     frameworkAtoms[j].scalingCoulomb * frameworkAtoms[i].charge *
                                     frameworkAtoms[j].charge;
        const double erfTerm = std::erf(alpha * r);
        energySum.ewald_exclusion -= chargeProduct * erfTerm / r;
        const double gaussTerm = 2.0 * alpha * std::numbers::inv_sqrtpi * std::exp(-alpha_squared * rr);
        const double f1 = -chargeProduct * (gaussTerm / rr - erfTerm / (r * rr));
        const double3 gradient = f1 * dr;
        frameworkDynamics[i].gradient += gradient;
        frameworkDynamics[j].gradient -= gradient;
      }
    }
  }

  // Net-charge correction: neutralizing background plus removal of the spurious interaction of the
  // net charge with its own periodic images; see Bogusz et al., J. Chem. Phys. 108, 7070 (1998).
  // The correction is independent of the atom positions, so it contributes no gradient.
  {
    double netChargeAdsorbates = 0.0;
    std::array<double, maximumNumberOfDUDlambdaGroups> netChargeDerivative{};
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      netChargeAdsorbates += atoms[i].scalingCoulomb * atoms[i].charge;
      if (atoms[i].groupId != 0) netChargeDerivative[atoms[i].groupId - 1] += atoms[i].charge;
    }
    const double rigidFrameworkCharge = flexibleFramework ? 0.0 : netChargeFramework;
    double uIon =
        -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
    if (omitInterInteractions)
    {
      energySum.ewald_fourier += 2.0 * uIon * rigidFrameworkCharge * netChargeAdsorbates;
      for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
      {
        energySum.dudlambdaEwald[g] += 2.0 * uIon * rigidFrameworkCharge * netChargeDerivative[g];
      }
    }
    else
    {
      energySum.ewald_fourier +=
          uIon * (2.0 * rigidFrameworkCharge + netChargeAdsorbates) * netChargeAdsorbates;
      for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
      {
        energySum.dudlambdaEwald[g] +=
            2.0 * uIon * (rigidFrameworkCharge + netChargeAdsorbates) * netChargeDerivative[g];
      }
    }
  }

  return energySum;
}


// Used in smart-MC
void Interactions::computeEwaldFourierGradientSingleMolecule(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> atoms,
    std::span<AtomDynamics> atomDynamics)
{
  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  std::size_t numberOfAtoms = atoms.size();
  if (numberOfAtoms == 0) return;

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

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
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

  // Iterate over the exact same set/ordering of wave vectors used to build 'storedEik' so that the
  // nvec index selects the matching total structure factor S(k).
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
          std::complex<double> total = storedEik[nvec].first;
          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = atoms[i].charge;
            double scaling = atoms[i].scalingCoulomb;
            atomDynamics[i].gradient -=
                scaling * charge * 2.0 * temp * (cki.imag() * total.real() - cki.real() * total.imag()) * rk;
          }

          ++nvec;
        }
      }
    }
  }
}

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<const Atom> newatoms, std::span<const Atom> oldatoms,
    double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups> &netChargeDerivativeExternal)
{
  RunningEnergy energy;
  double singleIonFourierSum = 0.0;

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
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_old;
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_new;
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
          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            std::uint8_t groupIdA = oldatoms[i].groupId;
            cksum_old.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum_old.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
          for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            std::uint8_t groupIdA = newatoms[i - oldatoms.size()].groupId;
            cksum_new.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum_new.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          singleIonFourierSum += temp;

          energy.ewald_fourier += temp * std::norm(storedEik[nvec].first + cksum_new.first - cksum_old.first);
          energy.ewald_fourier -= temp * std::norm(storedEik[nvec].first);

          addFourierDUdlambda(energy, 2.0 * temp, storedEik[nvec].first + cksum_new.first - cksum_old.first,
                              storedEik[nvec].second + cksum_new.second - cksum_old.second);
          addFourierDUdlambda(energy, -2.0 * temp, storedEik[nvec].first, storedEik[nvec].second);

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
    std::uint8_t groupIdA = oldatoms[i].groupId;
    double3 posA = oldatoms[i].position;
    for (std::size_t j = i + 1; j != oldatoms.size(); j++)
    {
      double chargeB = oldatoms[j].charge;
      double scalingB = oldatoms[j].scalingCoulomb;
      std::uint8_t groupIdB = oldatoms[j].groupId;
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion += scalingA * scalingB * temp;
      energy.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, temp);
    }
  }

  for (std::size_t i = 0; i != newatoms.size(); i++)
  {
    double chargeA = newatoms[i].charge;
    double scalingA = newatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = newatoms[i].groupId;
    double3 posA = newatoms[i].position;
    for (std::size_t j = i + 1; j != newatoms.size(); j++)
    {
      double chargeB = newatoms[j].charge;
      double scalingB = newatoms[j].scalingCoulomb;
      std::uint8_t groupIdB = newatoms[j].groupId;
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion -= scalingA * scalingB * temp;
      energy.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, -temp);
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = oldatoms[i].groupId;
    energy.ewald_self += prefactor_self * scaling * charge * scaling * charge;
    if (groupIdA != 0) energy.dudlambdaEwald[groupIdA - 1] += 2.0 * prefactor_self * scaling * charge * charge;
  }
  for (std::size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = newatoms[i].groupId;
    energy.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
    if (groupIdA != 0) energy.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
  }

  addNetChargeCorrectionDifference(energy, singleIonFourierSum, alpha, netCharge, netChargeDerivativeExternal, newatoms, oldatoms);

  return energy;
}

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups> &netChargeDerivativeExternal)
{
  RunningEnergy energy;
  double singleIonFourierSum = 0.0;

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
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

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
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_old;
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_new;
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

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid = fixedFrameworkStoredEik[nvec];

          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = oldatoms[i].charge;
            double scaling = oldatoms[i].scalingCoulomb;
            std::uint8_t groupIdA = oldatoms[i].groupId;
            cksum_old.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum_old.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
            electricFieldOld[i] -=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
          for (std::size_t i = oldatoms.size(); i != oldatoms.size() + newatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            double charge = newatoms[i - oldatoms.size()].charge;
            double scaling = newatoms[i - oldatoms.size()].scalingCoulomb;
            std::uint8_t groupIdA = newatoms[i - oldatoms.size()].groupId;
            cksum_new.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum_new.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
            electricFieldNew[i - oldatoms.size()] +=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          singleIonFourierSum += temp;
          energy.ewald_fourier += temp * std::norm(storedEik[nvec].first + cksum_new.first - cksum_old.first);
          energy.ewald_fourier -= temp * std::norm(storedEik[nvec].first);

          addFourierDUdlambda(energy, 2.0 * temp, storedEik[nvec].first + cksum_new.first - cksum_old.first,
                              storedEik[nvec].second + cksum_new.second - cksum_old.second);
          addFourierDUdlambda(energy, -2.0 * temp, storedEik[nvec].first, storedEik[nvec].second);

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
    std::uint8_t groupIdA = oldatoms[i].groupId;
    double3 posA = oldatoms[i].position;
    for (std::size_t j = i + 1; j != oldatoms.size(); j++)
    {
      double chargeB = oldatoms[j].charge;
      double scalingB = oldatoms[j].scalingCoulomb;
      std::uint8_t groupIdB = oldatoms[j].groupId;
      double3 posB = oldatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      double r = std::sqrt(rr);

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion += scalingA * scalingB * temp;
      energy.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, temp);
    }
  }

  for (std::size_t i = 0; i != newatoms.size(); i++)
  {
    double chargeA = newatoms[i].charge;
    double scalingA = newatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = newatoms[i].groupId;
    double3 posA = newatoms[i].position;
    for (std::size_t j = i + 1; j != newatoms.size(); j++)
    {
      double chargeB = newatoms[j].charge;
      double scalingB = newatoms[j].scalingCoulomb;
      std::uint8_t groupIdB = newatoms[j].groupId;
      double3 posB = newatoms[j].position;

      double3 dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      double r = std::sqrt(rr);

      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
      energy.ewald_exclusion -= scalingA * scalingB * temp;
      energy.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, -temp);
    }
  }

  // Subtract self-energy
  double prefactor_self = Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
  for (std::size_t i = 0; i != oldatoms.size(); ++i)
  {
    double charge = oldatoms[i].charge;
    double scaling = oldatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = oldatoms[i].groupId;
    energy.ewald_self += prefactor_self * scaling * charge * scaling * charge;
    if (groupIdA != 0) energy.dudlambdaEwald[groupIdA - 1] += 2.0 * prefactor_self * scaling * charge * charge;
  }

  for (std::size_t i = 0; i != newatoms.size(); ++i)
  {
    double charge = newatoms[i].charge;
    double scaling = newatoms[i].scalingCoulomb;
    std::uint8_t groupIdA = newatoms[i].groupId;
    energy.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
    if (groupIdA != 0) energy.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
  }

  addNetChargeCorrectionDifference(energy, singleIonFourierSum, alpha, netCharge, netChargeDerivativeExternal, newatoms, oldatoms);

  return energy;
}

// Used to compute the difference in electricField for a grown or retraced state
// Used in insertion_CBCMC and deletion_CBCMC

void Interactions::computeEwaldFourierElectricFieldDifference(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &totalEik, const ForceField &forceField,
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
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

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
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_old;
  std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum_new;
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

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid = fixedFrameworkStoredEik[nvec];

          cksum_old = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
          for (std::size_t i = 0; i != oldatoms.size(); ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            std::complex<double> cki = eik_xy[i] * eikz_temp;
            electricFieldOld[i] -=
                2.0 * temp * (cki.imag() * rigid.first.real() - cki.real() * rigid.first.imag()) * rk;
          }

          cksum_new = std::make_pair(std::complex<double>(0.0, 0.0), GroupComplexSums{});
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
                                   std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
                                   std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &totalEik)
{
  if (!forceField.useCharge) return;
  if (forceField.omitEwaldFourier) return;

  storedEik = totalEik;
}

std::pair<EnergyStatus, double3x3> Interactions::computeEwaldFourierEnergyStrainDerivative(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::optional<Framework> &framework,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> atomData, std::span<AtomDynamics> atomDynamics, double netChargeFramework,
    std::vector<double> netChargePerComponent) noexcept
{
  double alpha = forceField.EwaldAlpha;
  double alpha_squared = alpha * alpha;
  double singleIonFourierSum = 0.0;

  // Total net charge, used for the net-charge correction to the strain derivative; the strain
  // derivative includes the rigid-framework contribution, so the framework charge is included.
  double netChargeTotal = netChargeFramework;
  for (double q : netChargePerComponent)
  {
    netChargeTotal += q;
  }
  double netChargeTotalSquared = netChargeTotal * netChargeTotal;
  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  EnergyStatus energy(1, framework.has_value() ? 1uz : 0uz, components.size());
  double3x3 strainDerivative;

  if (!forceField.useCharge || forceField.omitEwaldFourier) return std::make_pair(energy, strainDerivative);

  std::size_t numberOfAtoms = atomData.size();
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
    double3 s = 2.0 * std::numbers::pi * (inv_box * atomData[i].position);
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
            std::size_t comp = static_cast<std::size_t>(atomData[i].componentId);
            double charge = atomData[i].charge;
            double scaling = atomData[i].scalingCoulomb;
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
            double charge = atomData[i].charge;
            double scaling = atomData[i].scalingCoulomb;

            atomDynamics[i].gradient -= scaling * charge * 2.0 * temp *
                                        (cki.imag() * test.real() - cki.real() * test.imag()) *
                                        (kvec_x + kvec_y + kvec_z);
          }

          singleIonFourierSum += temp;

          // Include the net-charge correction in the strain derivative: per wave vector its energy
          // contribution is -temp * Q_total^2, with the same k- and volume-dependence as the
          // regular Fourier term (the correction's self part, alpha/sqrt(pi), is strain-independent).
          double currentEnergy =
              temp * (test.real() * test.real() + test.imag() * test.imag() - netChargeTotalSquared);
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
  for (std::size_t i = 0; i != atomData.size(); ++i)
  {
    double charge = atomData[i].charge;
    double scaling = atomData[i].scalingCoulomb;
    std::size_t comp = static_cast<std::size_t>(atomData[i].componentId);
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
      std::span<const Atom> span = std::span(&atomData[index], size);
      std::span<AtomDynamics> spanDynamics = std::span(&atomDynamics[index], size);
      for (std::size_t i = 0; i != span.size() - 1; i++)
      {
        double chargeA = span[i].charge;
        double scalingA = span[i].scalingCoulomb;
        // std::uint8_t groupIdA = span[i].groupId;
        double3 posA = span[i].position;
        for (std::size_t j = i + 1; j != span.size(); j++)
        {
          double chargeB = span[j].charge;
          double scalingB = span[j].scalingCoulomb;
          // std::uint8_t groupIdB = span[j].groupId;
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
          spanDynamics[i].gradient -= temp * dr;
          spanDynamics[j].gradient += temp * dr;

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

  // Handle net-charges: neutralizing background plus removal of the spurious interaction of the
  // net charge with its own periodic images; see Bogusz et al., J. Chem. Phys. 108, 7070 (1998).
  // The single-ion energy is computed internally from the wave-vector sum so that it is always
  // consistent with the current simulation box.
  double uIon =
      -(singleIonFourierSum - Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi));
  for (std::size_t i = 0; i != components.size(); ++i)
  {
    energy.frameworkComponentEnergy(0, i).CoulombicFourier +=
        Potentials::EnergyFactor(2.0 * uIon * netChargeFramework * netChargePerComponent[i], 0.0);
  }

  for (std::size_t i = 0; i != components.size(); ++i)
  {
    for (std::size_t j = 0; j != components.size(); ++j)
    {
      energy.componentEnergy(i, j).CoulombicFourier +=
          Potentials::EnergyFactor(uIon * netChargePerComponent[i] * netChargePerComponent[j], 0.0);
    }
  }

  return std::make_pair(energy, strainDerivative);
}

void Interactions::computeEwaldFourierElectrostaticPotential(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
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
          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            std::uint8_t groupIdA = moleculeAtomPositions[i].groupId;
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> total;
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
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldMolecules,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> moleculeAtomPositions)
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
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

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

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> cksum;
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            double charge = moleculeAtomPositions[i].charge;
            double scaling = moleculeAtomPositions[i].scalingCoulomb;
            std::uint8_t groupIdA = moleculeAtomPositions[i].groupId;
            cksum.first += scaling * charge * (eik_xy[i] * eikz_temp);
            if (groupIdA != 0) cksum.second[groupIdA - 1] += charge * eik_xy[i] * eikz_temp;
          }

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> rigid = fixedFrameworkStoredEik[nvec];

          std::pair<std::complex<double>, std::array<std::complex<double>, 4>> total = rigid;
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

          addFourierDUdlambda(energySum, 2.0 * temp, total.first, total.second);
          addFourierDUdlambda(energySum, -2.0 * temp, rigid.first, rigid.second);

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
      std::uint8_t groupIdA = moleculeAtomPositions[i].groupId;
      energySum.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      if (groupIdA != 0) energySum.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
    }

    // Subtract exclusion-energy
    std::size_t index{0};
    for (std::size_t l = 0; l != components.size(); ++l)
    {
      std::size_t size = components[l].atoms.size();
      for (std::size_t m = 0; m != numberOfMoleculesPerComponent[l]; ++m)
      {
        std::span<Atom> span = std::span(&moleculeAtomPositions[index], size);
        for (std::size_t i = 0; i != span.size(); i++)
        {
          double chargeA = span[i].charge;
          double scalingA = span[i].scalingCoulomb;
          std::uint8_t groupIdA = span[i].groupId;
          double3 posA = span[i].position;
          for (std::size_t j = i + 1; j != span.size(); j++)
          {
            if (i != j)
            {
              double chargeB = span[j].charge;
              double scalingB = span[j].scalingCoulomb;
              std::uint8_t groupIdB = span[j].groupId;
              double3 posB = span[j].position;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double rr = double3::dot(dr, dr);
              double r = std::sqrt(rr);

              if (!omitInterInteractions)
              {
                double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erf(alpha * r) / r;
                energySum.ewald_exclusion -= scalingA * scalingB * temp;
                energySum.addDudlambdaEwald(groupIdA, groupIdB, scalingA, scalingB, -temp);
              }

              // NOTE: the intra-molecular Ewald reciprocal-space exclusion does NOT contribute to the
              // polarization electric field in this model. The reciprocal field is built solely from the
              // (fixed) framework structure factor, so there is no intra-molecular reciprocal term to
              // exclude. Adsorbate-adsorbate polarization is handled entirely in real space by
              // computeInterMolecularElectricField / -Difference (different molecules only). Adding the
              // Bt1 term here would make the stored field inconsistent with the incremental Monte-Carlo
              // moves and introduce energy drift. The exclusion energy itself is still accounted for above.
            }
          }
        }
        index += size;
      }
    }
  }

  return energySum;
}

void Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    const SimulationBox &simulationBox, std::span<const Atom> atoms, std::span<double> potentialMatrix)
{
  const std::size_t numberOfAtoms = atoms.size();
  const double volume = simulationBox.volume;
  const double3x3 inv_box = simulationBox.inverseCell;
  const double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  const double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  const double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  // Derive the Ewald parameters from the box using the same formulas as
  // ForceField::initializeEwaldParameters. The force-field Ewald parameters can not be used here,
  // because charge equilibration runs during CIF-reading, before they are initialized.
  const double3 perpendicularWidths = simulationBox.perpendicularWidths();
  const double cutOff =
      0.5 * std::min({perpendicularWidths.x, perpendicularWidths.y, perpendicularWidths.z});
  const double eps = 1.0e-8;
  const double tol = std::sqrt(std::abs(std::log(eps * cutOff)));
  const double alpha = std::sqrt(std::abs(std::log(eps * cutOff * tol))) / cutOff;
  const double tol1 = std::sqrt(-std::log(eps * cutOff * (2.0 * tol * alpha) * (2.0 * tol * alpha)));

  const std::size_t kx_max_unsigned = static_cast<std::size_t>(
      std::rint(0.25 + perpendicularWidths.x * alpha * tol1 / std::numbers::pi));
  const std::size_t ky_max_unsigned = static_cast<std::size_t>(
      std::rint(0.25 + perpendicularWidths.y * alpha * tol1 / std::numbers::pi));
  const std::size_t kz_max_unsigned = static_cast<std::size_t>(
      std::rint(0.25 + perpendicularWidths.z * alpha * tol1 / std::numbers::pi));
  const std::size_t recip_integer_cutoff_squared =
      std::max({kx_max_unsigned, ky_max_unsigned, kz_max_unsigned}) *
      std::max({kx_max_unsigned, ky_max_unsigned, kz_max_unsigned});

  const std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  const std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  const std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::fill(potentialMatrix.begin(), potentialMatrix.end(), 0.0);

  // Real-space part: minimum-image erfc(alpha * r) / r; the cutoff equals half the smallest
  // perpendicular width, so the minimum image is the only image that contributes.
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    for (std::size_t j = i + 1; j != numberOfAtoms; ++j)
    {
      double3 dr = atoms[i].position - atoms[j].position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double r = std::sqrt(double3::dot(dr, dr));
      potentialMatrix[i * numberOfAtoms + j] += std::erfc(alpha * r) / r;
    }
  }

  if (numberOfAtoms * (kx_max_unsigned + 1) > eik_x.size()) eik_x.resize(numberOfAtoms * (kx_max_unsigned + 1));
  if (numberOfAtoms * (ky_max_unsigned + 1) > eik_y.size()) eik_y.resize(numberOfAtoms * (ky_max_unsigned + 1));
  if (numberOfAtoms * (kz_max_unsigned + 1) > eik_z.size()) eik_z.resize(numberOfAtoms * (kz_max_unsigned + 1));
  if (numberOfAtoms > eik_xy.size()) eik_xy.resize(numberOfAtoms);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
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

  // Fourier part: for every wave vector the contribution to the matrix,
  // prefactor * cos(k.(r_i - r_j)) = prefactor * Re(e^{ik.r_i} conj(e^{ik.r_j})),
  // is accumulated as a symmetric rank-2 update using the per-atom phases.
  std::vector<std::complex<double>> phase(numberOfAtoms);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    double factor = (kx == 0) ? 1.0 : 2.0;

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
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared))
        {
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            phase[i] = eik_xy[i] * eikz_temp;
          }

          double prefactor = factor * (4.0 * std::numbers::pi / volume) *
                             std::exp((-0.25 / (alpha * alpha)) * rksq) / rksq;

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            const double re_i = prefactor * phase[i].real();
            const double im_i = prefactor * phase[i].imag();
            double *row = &potentialMatrix[i * numberOfAtoms];
            for (std::size_t j = i; j != numberOfAtoms; ++j)
            {
              row[j] += re_i * phase[j].real() + im_i * phase[j].imag();
            }
          }
        }
      }
    }
  }

  // Self-energy of the Gaussian screening charge and the neutralizing-background correction,
  // which together make the matrix independent of the choice of alpha.
  const double background = std::numbers::pi / (alpha * alpha * volume);
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    for (std::size_t j = i; j != numberOfAtoms; ++j)
    {
      potentialMatrix[i * numberOfAtoms + j] -= background;
    }
    potentialMatrix[i * numberOfAtoms + i] -= 2.0 * alpha / std::sqrt(std::numbers::pi);
  }

  // Mirror the upper triangle into the lower triangle
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    for (std::size_t j = i + 1; j != numberOfAtoms; ++j)
    {
      potentialMatrix[j * numberOfAtoms + i] = potentialMatrix[i * numberOfAtoms + j];
    }
  }
}
