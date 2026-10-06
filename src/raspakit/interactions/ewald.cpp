module;

module interactions_ewald;

import std;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import potential_coulomb_real_space;
import energy_status_inter;
import units;
import energy_dudlambda;
import running_energy;
import framework;
import component;
import intra_molecular_exclusions;
import coulomb_potential;
import forcefield;
import interactions_ewald_kvector;

namespace
{
void addRealSpaceSelfEnergy(RunningEnergy& energy, const ForceField& forceField, std::span<const Atom> atoms,
                            double sign = 1.0)
{
  const double prefactor = sign * Units::CoulombicConversionFactor * Potentials::coulombSelfEnergyPrefactor(forceField);
  for (const Atom& atom : atoms)
  {
    const double scaledCharge = atom.scalingCoulomb * atom.charge;
    energy.ewald_self += prefactor * scaledCharge * scaledCharge;
    if (atom.groupId != 0)
    {
      energy.dudlambdaEwald[atom.groupId - 1] += 2.0 * prefactor * atom.scalingCoulomb * atom.charge * atom.charge;
    }
  }
}

RunningEnergy realSpaceSelfEnergyDifference(const ForceField& forceField, std::span<const Atom> newAtoms,
                                            std::span<const Atom> oldAtoms)
{
  RunningEnergy energy{};
  if (forceField.omitInterInteractions) return energy;
  addRealSpaceSelfEnergy(energy, forceField, newAtoms);
  addRealSpaceSelfEnergy(energy, forceField, oldAtoms, -1.0);
  return energy;
}

// Charge exclusion corrections of the excluded intramolecular pairs (IntraMolecularExclusions: the 1-2 and 1-3
// pairs and the pairs inside one rigid fragment; every pair for a rigid molecule). The non-excluded pairs are
// regular pairs of the force field and are evaluated by the intramolecular pair terms of the component
// (Potentials::intraMolecularCoulomb), which carry their own Ewald/shifted-potential completion.
//
// Ewald: the Fourier sum counts every intramolecular pair, so an excluded pair receives
// -q_i q_j erf(alpha r)/r (the 'exclusion' term), without cutoff.
//
// Finite-cutoff charge methods (Wolf, damped-shifted-force, modified-shifted-force, zero-dipole): for every
// excluded pair inside the Coulomb cutoff the correction is q_i q_j (V(r) - 1/r), where V(r) is the method's
// shifted real-space potential. This completes the shifted pair sum over all atoms inside the cutoff (so the
// per-atom self term is balanced) and removes the bare 1/r Coulomb, matching Eqs. (S61)-(S64) of Dubbeldam et al.
// (the Brick-CFCMC formulation, terms S62 and S63).
//
// The 'sign' argument is +1 to add the exclusion of a configuration and -1 to remove it, matching the old/new
// convention of the energy-difference routines. The gradient uses the RASPA factor convention f = (dU/dr)/r.
struct ExclusionPairTerm
{
  double energy{0.0};              ///< scaled energy (sign included)
  double dUdlambda{0.0};           ///< argument of RunningEnergy::addDudlambdaEwald (sign included)
  double firstDerivativeFactor{0.0};  ///< scaled (dU/dr)/r (sign included)
  double secondDerivativeFactor{0.0};  ///< scaled d/dr[(dU/dr)/r]/r (sign included)
  bool active{false};
};

ExclusionPairTerm chargeExclusionPairTerm(const ForceField& forceField, const Atom& atomA, const Atom& atomB,
                                          double rr, double sign)
{
  ExclusionPairTerm term{};
  const double scalingA = atomA.scalingCoulomb;
  const double scalingB = atomB.scalingCoulomb;
  const double prefactor = sign * Units::CoulombicConversionFactor * atomA.charge * atomB.charge;
  if (forceField.usesEwaldFourier())
  {
    const double r = std::sqrt(rr);
    const Potentials::EwaldExclusionFactors exclusion =
        Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, scalingA * scalingB, r);
    term.energy = -scalingA * scalingB * prefactor * exclusion.potential;
    term.dUdlambda = -prefactor * exclusion.dUdlambda;
    term.firstDerivativeFactor = -scalingA * scalingB * prefactor * exclusion.firstDerivativeFactor;
    term.secondDerivativeFactor = -scalingA * scalingB * prefactor * exclusion.secondDerivativeFactor;
    term.active = true;
  }
  else if (forceField.usesRealSpaceChargeCorrections())
  {
    if (rr >= forceField.cutOffCoulomb * forceField.cutOffCoulomb) return term;
    const double r = std::sqrt(rr);
    const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
    term.energy = scalingA * scalingB * prefactor * (factors.potential - 1.0 / r);
    term.dUdlambda = prefactor * (factors.potential - 1.0 / r);
    term.firstDerivativeFactor = scalingA * scalingB * prefactor * (factors.firstDerivativeFactor + 1.0 / (rr * r));
    term.secondDerivativeFactor =
        scalingA * scalingB * prefactor * (factors.secondDerivativeFactor - 3.0 / (rr * rr * r));
    term.active = true;
  }
  return term;
}

void addChargeExclusionPairEnergy(RunningEnergy& energy, const ForceField& forceField,
                                  const SimulationBox& simulationBox, const Atom& atomA, const Atom& atomB,
                                  double sign)
{
  const double3 dr = simulationBox.applyPeriodicBoundaryConditions(atomA.position - atomB.position);
  const ExclusionPairTerm term = chargeExclusionPairTerm(forceField, atomA, atomB, double3::dot(dr, dr), sign);
  if (!term.active) return;
  energy.ewald_exclusion += term.energy;
  energy.addDudlambdaEwald(atomA.groupId, atomB.groupId, atomA.scalingCoulomb, atomB.scalingCoulomb, term.dUdlambda);
}

// Exclusion energy of the molecules in 'atoms' (whole molecules, see forEachExcludedPair).
void addChargeExclusionEnergy(RunningEnergy& energy, const ForceField& forceField, const SimulationBox& simulationBox,
                              std::span<const Component> components, std::span<const Atom> atoms, double sign = 1.0)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier() && !forceField.usesRealSpaceChargeCorrections()) return;
  forEachExcludedPair(components, atoms,
                      [&](std::size_t i, std::size_t j)
                      { addChargeExclusionPairEnergy(energy, forceField, simulationBox, atoms[i], atoms[j], sign); });
}

// Gradient-aware variant of addChargeExclusionEnergy: also accumulates the pair contribution to the atomic
// gradients and, if requested, to the strain derivative.
void addChargeExclusionGradient(RunningEnergy& energy, const ForceField& forceField,
                                const SimulationBox& simulationBox, std::span<const Component> components,
                                std::span<const Atom> atoms, std::span<AtomDynamics> dynamics, double sign = 1.0,
                                double3x3* strainDerivative = nullptr)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier() && !forceField.usesRealSpaceChargeCorrections()) return;
  forEachExcludedPair(components, atoms,
                      [&](std::size_t i, std::size_t j)
                      {
                        const double3 dr =
                            simulationBox.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
                        const ExclusionPairTerm term =
                            chargeExclusionPairTerm(forceField, atoms[i], atoms[j], double3::dot(dr, dr), sign);
                        if (!term.active) return;
                        energy.ewald_exclusion += term.energy;
                        energy.addDudlambdaEwald(atoms[i].groupId, atoms[j].groupId, atoms[i].scalingCoulomb,
                                                 atoms[j].scalingCoulomb, term.dUdlambda);
                        const double3 f = term.firstDerivativeFactor * dr;
                        dynamics[i].gradient += f;
                        dynamics[j].gradient -= f;
                        if (strainDerivative)
                        {
                          strainDerivative->ax += f.x * dr.x;
                          strainDerivative->bx += f.y * dr.x;
                          strainDerivative->cx += f.z * dr.x;
                          strainDerivative->ay += f.x * dr.y;
                          strainDerivative->by += f.y * dr.y;
                          strainDerivative->cy += f.z * dr.y;
                          strainDerivative->az += f.x * dr.z;
                          strainDerivative->bz += f.y * dr.z;
                          strainDerivative->cz += f.z * dr.z;
                        }
                      });
}

// Exclusion corrections with the per-component energy bookkeeping of the strain-derivative routine: energy into
// the diagonal CoulombicFourier entry of the component, atomic gradients, and the strain derivative.
void addChargeExclusionStrainDerivative(EnergyStatus& energy, double3x3& strainDerivative,
                                        const ForceField& forceField, const SimulationBox& simulationBox,
                                        std::span<const Component> components, std::span<const Atom> atoms,
                                        std::span<AtomDynamics> dynamics)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier() && !forceField.usesRealSpaceChargeCorrections()) return;
  forEachExcludedPair(components, atoms,
                      [&](std::size_t i, std::size_t j)
                      {
                        const double3 dr =
                            simulationBox.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
                        const ExclusionPairTerm term =
                            chargeExclusionPairTerm(forceField, atoms[i], atoms[j], double3::dot(dr, dr), 1.0);
                        if (!term.active) return;
                        const std::size_t comp = static_cast<std::size_t>(atoms[i].componentId);
                        energy.componentEnergy(comp, comp).CoulombicFourier += EnergyDuDlambda(term.energy, 0.0);

                        const double3 f = term.firstDerivativeFactor * dr;
                        dynamics[i].gradient += f;
                        dynamics[j].gradient -= f;

                        strainDerivative.ax += f.x * dr.x;
                        strainDerivative.bx += f.y * dr.x;
                        strainDerivative.cx += f.z * dr.x;
                        strainDerivative.ay += f.x * dr.y;
                        strainDerivative.by += f.y * dr.y;
                        strainDerivative.cy += f.z * dr.y;
                        strainDerivative.az += f.x * dr.z;
                        strainDerivative.bz += f.y * dr.z;
                        strainDerivative.cz += f.z * dr.z;
                      });
}
}  // namespace

void Interactions::addChargeSelfEnergy(RunningEnergy& energy, const ForceField& forceField, std::span<const Atom> atoms)
{
  if (!forceField.useCharge || forceField.omitInterInteractions) return;
  if (forceField.usesEwaldFourier())
  {
    const double prefactor_self =
        Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi);
    for (const Atom& atom : atoms)
    {
      const double charge = atom.charge;
      const double scaling = atom.scalingCoulomb;
      const std::uint8_t groupIdA = atom.groupId;
      energy.ewald_self -= prefactor_self * scaling * charge * scaling * charge;
      if (groupIdA != 0) energy.dudlambdaEwald[groupIdA - 1] -= 2.0 * prefactor_self * scaling * charge * charge;
    }
  }
  else if (forceField.usesRealSpaceChargeCorrections())
  {
    addRealSpaceSelfEnergy(energy, forceField, atoms);
  }
}

void Interactions::addIntraMolecularChargeExclusionGradient(RunningEnergy& energy, const ForceField& forceField,
                                                            const SimulationBox& simulationBox,
                                                            std::span<const Component> components,
                                                            std::span<const Atom> moleculeAtoms,
                                                            std::span<AtomDynamics> moleculeDynamics,
                                                            double3x3* strainDerivative)
{
  addChargeExclusionGradient(energy, forceField, simulationBox, components, moleculeAtoms, moleculeDynamics, 1.0,
                             strainDerivative);
}

Interactions::ChargeExclusionPairTerm Interactions::chargeExclusionPairTerm(const ForceField& forceField,
                                                                            const Atom& atomA, const Atom& atomB,
                                                                            double rr)
{
  ChargeExclusionPairTerm result{};
  if (!forceField.useCharge) return result;
  if (!forceField.usesEwaldFourier() && !forceField.usesRealSpaceChargeCorrections()) return result;
  const ExclusionPairTerm term = ::chargeExclusionPairTerm(forceField, atomA, atomB, rr, 1.0);
  result.energy = term.energy;
  result.dUdlambda = term.dUdlambda;
  result.firstDerivativeFactor = term.firstDerivativeFactor;
  result.active = term.active;
  return result;
}

RunningEnergy Interactions::computeChargeSelfAndExclusionGradient(
    const ForceField& forceField, const SimulationBox& simulationBox, const std::vector<Component>& components,
    [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent, std::span<const Atom> atomData,
    std::span<AtomDynamics> atomDynamics)
{
  RunningEnergy energy{};
  if (!forceField.useCharge || forceField.omitInterInteractions) return energy;
  if (!forceField.usesEwaldFourier() && !forceField.usesRealSpaceChargeCorrections()) return energy;

  addChargeSelfEnergy(energy, forceField, atomData);
  addChargeExclusionGradient(energy, forceField, simulationBox, components, atomData, atomDynamics);
  return energy;
}

// Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections
// to the Ewald potential
// Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks
// J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
//

// Per-group structure-factor derivative sums: element g holds sum_i charge_i * exp(ik.r_i) over the
// atoms tagged with dU/dlambda group id g+1 (Atom::groupId, 0 means untracked).
using GroupComplexSums = std::array<std::complex<double>, maximumNumberOfDUDlambdaGroups>;

static inline GroupComplexSums operator+(const GroupComplexSums& a, const GroupComplexSums& b)
{
  GroupComplexSums result;
  for (std::size_t g = 0; g != a.size(); ++g) result[g] = a[g] + b[g];
  return result;
}

static inline GroupComplexSums operator-(const GroupComplexSums& a, const GroupComplexSums& b)
{
  GroupComplexSums result;
  for (std::size_t g = 0; g != a.size(); ++g) result[g] = a[g] - b[g];
  return result;
}

static inline GroupComplexSums& operator+=(GroupComplexSums& a, const GroupComplexSums& b)
{
  for (std::size_t g = 0; g != a.size(); ++g) a[g] += b[g];
  return a;
}

static inline GroupComplexSums& operator-=(GroupComplexSums& a, const GroupComplexSums& b)
{
  for (std::size_t g = 0; g != a.size(); ++g) a[g] -= b[g];
  return a;
}

// Accumulates the Fourier-space dU/dlambda contribution factor * Re(sk * conj(dsk[g])) for each group g.
static inline void addFourierDUdlambda(RunningEnergy& energy, double factor, const std::complex<double>& sk,
                                       const GroupComplexSums& dsk)
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
    RunningEnergy& energy, double singleIonFourierSum, double alpha, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  double deltaCharge = 0.0;
  std::array<double, maximumNumberOfDUDlambdaGroups> chargeDerivativeNew = netChargeDerivativeExternal;
  std::array<double, maximumNumberOfDUDlambdaGroups> chargeDerivativeOld = netChargeDerivativeExternal;
  for (const Atom& atom : oldatoms)
  {
    deltaCharge -= atom.scalingCoulomb * atom.charge;
    if (atom.groupId != 0) chargeDerivativeOld[atom.groupId - 1] += atom.charge;
  }
  for (const Atom& atom : newatoms)
  {
    deltaCharge += atom.scalingCoulomb * atom.charge;
    if (atom.groupId != 0) chargeDerivativeNew[atom.groupId - 1] += atom.charge;
  }

  double uIon = -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
  double netChargeNew = netCharge + deltaCharge;
  energy.ewald_fourier += uIon * (netChargeNew * netChargeNew - netCharge * netCharge);
  for (std::size_t g = 0; g != maximumNumberOfDUDlambdaGroups; ++g)
  {
    energy.dudlambdaEwald[g] +=
        2.0 * uIon * (netChargeNew * chargeDerivativeNew[g] - netCharge * chargeDerivativeOld[g]);
  }
}

double Interactions::computeEwaldFourierEnergySingleIon(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy, const ForceField& forceField,
    const SimulationBox& simulationBox, double3 position, double charge)
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

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, 1, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned, inv_box,
                        [&position](std::size_t) { return position; });

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
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> rigidFrameworkAtoms)
{
  double3x3 inv_box = simulationBox.inverseCell;
  double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier()) return;

  std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  std::size_t numberOfAtoms = rigidFrameworkAtoms.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors)
  {
    fixedFrameworkStoredEik.resize(numberOfWaveVectors);
  }

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, rigidFrameworkAtoms, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

  std::size_t nvec = 0;
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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
// Volume-move, called with 'trialEik' for 'storedEik'
RunningEnergy Interactions::computeEwaldFourierEnergy(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    const ForceField& forceField, const SimulationBox& simulationBox, const std::vector<Component>& components,
    [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtomPositions,
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
  if (!forceField.usesEwaldFourier())
  {
    if (forceField.usesRealSpaceChargeCorrections() && !forceField.omitInterInteractions)
    {
      addRealSpaceSelfEnergy(energySum, forceField, moleculeAtomPositions);
      // Intramolecular exclusion / completion of the shifted pair sum (see addChargeExclusionEnergy).
      addChargeExclusionEnergy(energySum, forceField, simulationBox, components, moleculeAtomPositions);
    }
    return energySum;
  }

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, moleculeAtomPositions, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

    // Subtract exclusion-energy (the excluded intramolecular pairs only)
    addChargeExclusionEnergy(energySum, forceField, simulationBox, components, moleculeAtomPositions);
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
      if (moleculeAtomPositions[i].groupId != 0)
        netChargeDerivative[moleculeAtomPositions[i].groupId - 1] += moleculeAtomPositions[i].charge;
    }
    double uIon = -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
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
        energySum.dudlambdaEwald[g] += 2.0 * uIon * (netChargeFramework + netChargeAdsorbates) * netChargeDerivative[g];
      }
    }
  }

  return energySum;
}

// compute gradient
RunningEnergy Interactions::computeEwaldFourierGradient(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    const ForceField& forceField, const SimulationBox& simulationBox, const std::vector<Component>& components,
    [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent, std::span<const Atom> atomData,
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
  if (!forceField.usesEwaldFourier())
  {
    if (forceField.usesRealSpaceChargeCorrections() && !forceField.omitInterInteractions)
    {
      addRealSpaceSelfEnergy(energySum, forceField, atomData);

      // Intramolecular exclusion / completion of the shifted pair sum, including atomic gradients.
      addChargeExclusionGradient(energySum, forceField, simulationBox, components, atomData, atomDynamics);

      // Flexible-framework counterpart: the framework charges carry a self term and the bonded (1-2, 1-3, 1-4)
      // framework pairs excluded from the real-space pair sum need the shifted completion q_i q_j (V(r) - 1/r),
      // with the resulting force on the framework atoms. Mirrors the erf-based Ewald framework exclusion block.
      const bool liveFramework =
          framework && framework->hasMobileAtoms() && frameworkDynamics.size() == frameworkAtoms.size();
      if (liveFramework)
      {
        addRealSpaceSelfEnergy(energySum, forceField, frameworkAtoms);

        const double cutOffSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
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
          double3 dr =
              simulationBox.applyPeriodicBoundaryConditions(frameworkAtoms[i].position - frameworkAtoms[j].position);
          const double rr = double3::dot(dr, dr);
          if (rr >= cutOffSquared) continue;
          const double r = std::sqrt(rr);

          const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
          const double chargeProduct = Units::CoulombicConversionFactor * frameworkAtoms[i].scalingCoulomb *
                                       frameworkAtoms[j].scalingCoulomb * frameworkAtoms[i].charge *
                                       frameworkAtoms[j].charge;

          energySum.ewald_exclusion += chargeProduct * (factors.potential - 1.0 / r);

          const double gradientFactor = chargeProduct * (factors.firstDerivativeFactor + 1.0 / (rr * r));
          const double3 gradient = gradientFactor * dr;
          frameworkDynamics[i].gradient += gradient;
          frameworkDynamics[j].gradient -= gradient;
        }
      }
    }
    return energySum;
  }

  // Live Fourier hosts are mobile framework atoms only; lab-fixed atoms contribute via fixedFrameworkStoredEik.
  // When no Framework is passed (legacy callers), treat the host as fixed and use the precomputed eik.
  // When the framework is fully flexible, fixedFrameworkAtomCount is zero and the eik is omitted.
  const bool liveFramework =
      framework && framework->hasMobileAtoms() && frameworkDynamics.size() == frameworkAtoms.size();
  const std::size_t fixedFrameworkAtomCount = framework ? framework->numberOfFixedAtoms() : 0;
  const std::size_t mobileFrameworkAtomCount = framework ? framework->numberOfMobileAtoms() : 0;
  const bool useFixedFrameworkEik = !liveFramework || fixedFrameworkAtomCount > 0;
  const std::size_t frameworkOffset = liveFramework ? mobileFrameworkAtomCount : 0;
  std::vector<Atom> liveAtoms;
  if (liveFramework)
  {
    liveAtoms.reserve(mobileFrameworkAtomCount + atomData.size());
    liveAtoms.insert(liveAtoms.end(),
                     frameworkAtoms.begin() + static_cast<std::ptrdiff_t>(fixedFrameworkAtomCount),
                     frameworkAtoms.end());
    liveAtoms.insert(liveAtoms.end(), atomData.begin(), atomData.end());
  }
  const std::span<const Atom> atoms = liveFramework ? std::span<const Atom>(liveAtoms) : atomData;
  const std::size_t numberOfAtoms = atoms.size();
  const auto addGradient = [&](std::size_t atom, const double3& gradient)
  {
    if (liveFramework && atom < frameworkOffset)
    {
      frameworkDynamics[fixedFrameworkAtomCount + atom].gradient += gradient;
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

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (trialEik.size() < numberOfWaveVectors) trialEik.resize(numberOfWaveVectors);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, atoms, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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
          if (useFixedFrameworkEik) rigid = fixedFrameworkStoredEik[nvec];

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

          trialEik[nvec] = total;
          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy (of the molecule atoms and, for a live framework, its mobile atoms)
    addChargeSelfEnergy(energySum, forceField, atoms);

    // Subtract exclusion-energy of every molecule (the excluded intramolecular pairs only)
    addChargeExclusionGradient(energySum, forceField, simulationBox, components, atomData, atomDynamics);

    if (liveFramework)
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
    double rigidFrameworkCharge = 0.0;
    if (useFixedFrameworkEik)
    {
      if (!liveFramework)
      {
        rigidFrameworkCharge = netChargeFramework;
      }
      else
      {
        for (std::size_t i = 0; i < fixedFrameworkAtomCount; ++i)
        {
          rigidFrameworkCharge += frameworkAtoms[i].scalingCoulomb * frameworkAtoms[i].charge;
        }
      }
    }
    double uIon = -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
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
      energySum.ewald_fourier += uIon * (2.0 * rigidFrameworkCharge + netChargeAdsorbates) * netChargeAdsorbates;
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
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> atoms,
    std::span<AtomDynamics> atomDynamics)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier()) return;

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

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, atoms, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

RunningEnergy fourierSelfNetChargeDifference(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal);

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Component> components,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal)
{
  RunningEnergy energy = fourierSelfNetChargeDifference(eik_x, eik_y, eik_z, eik_xy, storedEik, trialEik, forceField,
                                                        simulationBox, newatoms, oldatoms, netCharge,
                                                        netChargeDerivativeExternal);
  if (!forceField.useCharge) return energy;
  if (!forceField.usesEwaldFourier() && forceField.omitInterInteractions) return energy;
  // Intramolecular exclusion difference: add the new configuration, remove the old one.
  addChargeExclusionEnergy(energy, forceField, simulationBox, components, newatoms, 1.0);
  addChargeExclusionEnergy(energy, forceField, simulationBox, components, oldatoms, -1.0);
  return energy;
}

// The Fourier sum (and the structure-factor update into 'trialEik'), the self energy, and the net-charge
// correction of energyDifferenceEwaldFourier: everything but the intramolecular exclusion, so that it can be
// evaluated for a subset of the atoms of a molecule (energyDifferenceEwaldFourierMovedAtoms).
RunningEnergy fourierSelfNetChargeDifference(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal)
{
  RunningEnergy energy;
  double singleIonFourierSum = 0.0;

  if (!forceField.useCharge) return energy;
  if (!forceField.usesEwaldFourier())
  {
    if (!forceField.usesRealSpaceChargeCorrections()) return RunningEnergy{};
    return realSpaceSelfEnergyDifference(forceField, newatoms, oldatoms);
  }

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

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (trialEik.size() < numberOfWaveVectors) trialEik.resize(numberOfWaveVectors);

  Interactions::Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, oldatoms, newatoms, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Interactions::Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

          trialEik[nvec].first = storedEik[nvec].first + cksum_new.first - cksum_old.first;
          trialEik[nvec].second = storedEik[nvec].second + cksum_new.second - cksum_old.second;

          ++nvec;
        }
      }
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

  addNetChargeCorrectionDifference(energy, singleIonFourierSum, alpha, netCharge, netChargeDerivativeExternal, newatoms,
                                   oldatoms);

  return energy;
}

std::vector<std::size_t> Interactions::movedAtomIndices(std::span<const Atom> newMolecule,
                                                        std::span<const Atom> oldMolecule)
{
  std::vector<std::size_t> moved;
  const std::size_t common = std::min(newMolecule.size(), oldMolecule.size());
  for (std::size_t i = 0; i != common; ++i)
  {
    const double3& a = newMolecule[i].position;
    const double3& b = oldMolecule[i].position;
    if (a.x != b.x || a.y != b.y || a.z != b.z) moved.push_back(i);
  }
  for (std::size_t i = common; i != std::max(newMolecule.size(), oldMolecule.size()); ++i) moved.push_back(i);
  return moved;
}

RunningEnergy Interactions::energyDifferenceEwaldFourierMovedAtoms(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Component> components,
    std::span<const Atom> newMolecule, std::span<const Atom> oldMolecule, std::span<const std::size_t> movedIndices,
    double netCharge, const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal)
{
  RunningEnergy energy;
  if (!forceField.useCharge) return energy;

  // Nothing moved: no change (the structure factor stays as stored). The caller's 'acceptEwaldMove' copies
  // 'trialEik' into 'storedEik', so keep the two consistent.
  if (movedIndices.empty())
  {
    trialEik = storedEik;
    return energy;
  }

  // Everything moved, or no common indexing between the two spans: the general routine.
  if (newMolecule.size() != oldMolecule.size() || movedIndices.size() >= newMolecule.size())
  {
    return energyDifferenceEwaldFourier(eik_x, eik_y, eik_z, eik_xy, storedEik, trialEik, forceField, simulationBox,
                                        components, newMolecule, oldMolecule, netCharge, netChargeDerivativeExternal);
  }

  std::vector<Atom> movedNew;
  std::vector<Atom> movedOld;
  movedNew.reserve(movedIndices.size());
  movedOld.reserve(movedIndices.size());
  for (std::size_t index : movedIndices)
  {
    movedNew.push_back(newMolecule[index]);
    movedOld.push_back(oldMolecule[index]);
  }

  // Fourier sum (and the structure-factor update into 'trialEik'), self energy, and net-charge correction: all of
  // these only involve the moved atoms, because an unmoved atom contributes identically to the new and the old
  // configuration.
  energy = fourierSelfNetChargeDifference(eik_x, eik_y, eik_z, eik_xy, storedEik, trialEik, forceField, simulationBox,
                                          movedNew, movedOld, netCharge, netChargeDerivativeExternal);

  // The excluded intramolecular pairs with at least one moved atom: new configuration with sign +1, old
  // configuration with sign -1 (the convention of energyDifferenceEwaldFourier).
  if (!forceField.usesEwaldFourier() && (!forceField.usesRealSpaceChargeCorrections() || forceField.omitInterInteractions))
  {
    return energy;
  }
  std::vector<bool> isMoved(newMolecule.size(), false);
  for (std::size_t index : movedIndices) isMoved[index] = true;
  forEachExcludedPair(components, newMolecule,
                      [&](std::size_t i, std::size_t j)
                      {
                        if (!isMoved[i] && !isMoved[j]) return;
                        addChargeExclusionPairEnergy(energy, forceField, simulationBox, newMolecule[i], newMolecule[j],
                                                     1.0);
                        addChargeExclusionPairEnergy(energy, forceField, simulationBox, oldMolecule[i], oldMolecule[j],
                                                     -1.0);
                      });
  return energy;
}

RunningEnergy Interactions::energyDifferenceEwaldFourier(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Component> components,
    std::span<double3> electricFieldNew, std::span<double3> electricFieldOld, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms, double netCharge,
    const std::array<double, maximumNumberOfDUDlambdaGroups>& netChargeDerivativeExternal)
{
  RunningEnergy energy;
  double singleIonFourierSum = 0.0;

  if (!forceField.useCharge) return energy;
  if (!forceField.usesEwaldFourier())
  {
    if (!forceField.usesRealSpaceChargeCorrections()) return RunningEnergy{};
    RunningEnergy realSpaceEnergy = realSpaceSelfEnergyDifference(forceField, newatoms, oldatoms);
    if (!forceField.omitInterInteractions)
    {
      // Intramolecular exclusion difference: add the new configuration, remove the old one.
      addChargeExclusionEnergy(realSpaceEnergy, forceField, simulationBox, components, newatoms, 1.0);
      addChargeExclusionEnergy(realSpaceEnergy, forceField, simulationBox, components, oldatoms, -1.0);
    }
    return realSpaceEnergy;
  }

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

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (trialEik.size() < numberOfWaveVectors) trialEik.resize(numberOfWaveVectors);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, oldatoms, newatoms, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

          trialEik[nvec].first = storedEik[nvec].first + cksum_new.first - cksum_old.first;
          trialEik[nvec].second = storedEik[nvec].second + cksum_new.second - cksum_old.second;

          ++nvec;
        }
      }
    }
  }

  // Intramolecular exclusion difference (the excluded pairs only): add the new configuration, remove the old one.
  addChargeExclusionEnergy(energy, forceField, simulationBox, components, newatoms, 1.0);
  addChargeExclusionEnergy(energy, forceField, simulationBox, components, oldatoms, -1.0);

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

  addNetChargeCorrectionDifference(energy, singleIonFourierSum, alpha, netCharge, netChargeDerivativeExternal, newatoms,
                                   oldatoms);

  return energy;
}

// Used to compute the difference in electricField for a grown or retraced state
// Used in insertion_CBCMC and deletion_CBCMC

void Interactions::computeEwaldFourierElectricFieldDifference(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<double3> electricFieldNew,
    std::span<double3> electricFieldOld, std::span<const Atom> newatoms, std::span<const Atom> oldatoms)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier()) return;

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

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (trialEik.size() < numberOfWaveVectors) trialEik.resize(numberOfWaveVectors);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, oldatoms, newatoms, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

void Interactions::acceptEwaldMove(
    const ForceField& forceField,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik)
{
  if (!forceField.useCharge) return;
  if (!forceField.usesEwaldFourier()) return;

  storedEik = trialEik;
}

std::pair<EnergyStatus, double3x3> Interactions::computeEwaldFourierEnergyStrainDerivative(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    const ForceField& forceField, const SimulationBox& simulationBox, const std::optional<Framework>& framework,
    const std::vector<Component>& components, [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent,
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

  // Finite-cutoff charge methods (Wolf, damped-shifted-force, modified-shifted-force, zero-dipole) have no
  // reciprocal-space contribution. Their real-space inter-molecular pair virial is accumulated in the
  // inter-molecular and framework-molecule strain routines. The remaining electrostatic bookkeeping handled
  // here is (i) the per-atom self-energy, which depends only on the charges, the damping and the fixed cutoff
  // and is therefore strain-independent, and (ii) the intra-molecular exclusion / completion of the shifted
  // pair sum, q_i q_j (V(r) - 1/r), which IS distance dependent and thus contributes to both the atomic
  // gradient and the strain derivative (mirroring the erf-based Ewald exclusion block below). Populating the
  // energy decomposition and the strain tensor keeps computeEwaldFourierEnergyStrainDerivative consistent with
  // the strain response of computeEwaldFourierEnergy for these methods.
  if (!forceField.usesEwaldFourier())
  {
    if (!forceField.usesRealSpaceChargeCorrections() || forceField.omitInterInteractions)
      return std::make_pair(energy, strainDerivative);

    // Self-energy (strain-independent).
    double selfPrefactor = Units::CoulombicConversionFactor * Potentials::coulombSelfEnergyPrefactor(forceField);
    for (std::size_t i = 0; i != atomData.size(); ++i)
    {
      double scaledCharge = atomData[i].scalingCoulomb * atomData[i].charge;
      std::size_t comp = static_cast<std::size_t>(atomData[i].componentId);
      energy.componentEnergy(comp, comp).CoulombicFourier +=
          EnergyDuDlambda(selfPrefactor * scaledCharge * scaledCharge, 0.0);
    }

    // Intra-molecular exclusion / completion q_i q_j (V(r) - 1/r) of the excluded pairs: energy, gradient and
    // strain derivative.
    addChargeExclusionStrainDerivative(energy, strainDerivative, forceField, simulationBox, components, atomData,
                                       atomDynamics);

    return std::make_pair(energy, strainDerivative);
  }

  std::size_t numberOfAtoms = atomData.size();
  std::size_t numberOfComponents = components.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, atomData, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned,
                        inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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
                EnergyDuDlambda(2.0 * temp *
                                             (fixedFrameworkStoredEik[nvec].first.real() * cksum[i].real() +
                                              fixedFrameworkStoredEik[nvec].first.imag() * cksum[i].imag()),
                                         0.0);
            for (std::size_t j = 0; j != numberOfComponents; ++j)
            {
              energy.componentEnergy(i, j).CoulombicFourier += EnergyDuDlambda(
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
          double currentEnergy = temp * (test.real() * test.real() + test.imag() * test.imag() - netChargeTotalSquared);
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
        EnergyDuDlambda(prefactor_self * scaling * charge * scaling * charge, 0.0);
  }

  // Subtract exclusion-energy (the excluded intramolecular pairs only), with gradient and strain derivative
  addChargeExclusionStrainDerivative(energy, strainDerivative, forceField, simulationBox, components, atomData,
                                     atomDynamics);

  // Handle net-charges: neutralizing background plus removal of the spurious interaction of the
  // net charge with its own periodic images; see Bogusz et al., J. Chem. Phys. 108, 7070 (1998).
  // The single-ion energy is computed internally from the wave-vector sum so that it is always
  // consistent with the current simulation box.
  double uIon =
      -(singleIonFourierSum - Units::CoulombicConversionFactor * forceField.EwaldAlpha / std::sqrt(std::numbers::pi));
  for (std::size_t i = 0; i != components.size(); ++i)
  {
    energy.frameworkComponentEnergy(0, i).CoulombicFourier +=
        EnergyDuDlambda(2.0 * uIon * netChargeFramework * netChargePerComponent[i], 0.0);
  }

  for (std::size_t i = 0; i != components.size(); ++i)
  {
    for (std::size_t j = 0; j != components.size(); ++j)
    {
      energy.componentEnergy(i, j).CoulombicFourier +=
          EnergyDuDlambda(uIon * netChargePerComponent[i] * netChargePerComponent[j], 0.0);
    }
  }

  return std::make_pair(energy, strainDerivative);
}

void Interactions::computeEwaldFourierElectrostaticPotential(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    std::span<double> electricPotentialMolecules, const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<Component>& components, [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent,
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
  if (!forceField.usesEwaldFourier()) return;

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  // std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, moleculeAtomPositions, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

    // Subtract the exclusion potential of the excluded intramolecular pairs
    forEachExcludedPair(components, moleculeAtomPositions,
                        [&](std::size_t i, std::size_t j)
                        {
                          const double3 dr = simulationBox.applyPeriodicBoundaryConditions(
                              moleculeAtomPositions[i].position - moleculeAtomPositions[j].position);
                          const double r = std::sqrt(double3::dot(dr, dr));
                          const double erfTerm = Units::CoulombicConversionFactor * std::erf(alpha * r) / r;
                          electricPotentialMolecules[i] -=
                              moleculeAtomPositions[j].scalingCoulomb * moleculeAtomPositions[j].charge * erfTerm;
                          electricPotentialMolecules[j] -=
                              moleculeAtomPositions[i].scalingCoulomb * moleculeAtomPositions[i].charge * erfTerm;
                        });
  }
}

RunningEnergy Interactions::computeEwaldFourierElectricField(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& storedEik,
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<double3> electricFieldMolecules,
    const std::vector<Component>& components, [[maybe_unused]] const std::vector<std::size_t>& numberOfMoleculesPerComponent,
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
  if (!forceField.usesEwaldFourier())
  {
    if (forceField.usesRealSpaceChargeCorrections() && !forceField.omitInterInteractions)
      addRealSpaceSelfEnergy(energySum, forceField, moleculeAtomPositions);
    return energySum;
  }

  std::size_t numberOfAtoms = moleculeAtomPositions.size();

  std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::size_t numberOfWaveVectors = (kx_max_unsigned + 1) * 2 * (ky_max_unsigned + 1) * 2 * (kz_max_unsigned + 1);
  if (storedEik.size() < numberOfWaveVectors) storedEik.resize(numberOfWaveVectors);
  if (fixedFrameworkStoredEik.size() < numberOfWaveVectors) fixedFrameworkStoredEik.resize(numberOfWaveVectors);

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, moleculeAtomPositions, kx_max_unsigned, ky_max_unsigned,
                        kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

    // Subtract exclusion-energy (the excluded intramolecular pairs only).
    //
    // NOTE: the intra-molecular Ewald reciprocal-space exclusion does NOT contribute to the polarization electric
    // field in this model. The reciprocal field is built solely from the (fixed) framework structure factor, so
    // there is no intra-molecular reciprocal term to exclude. Adsorbate-adsorbate polarization is handled entirely
    // in real space by computeInterMolecularElectricField / -Difference (different molecules only). Adding the Bt1
    // term here would make the stored field inconsistent with the incremental Monte-Carlo moves and introduce
    // energy drift. The exclusion energy itself is still accounted for here.
    addChargeExclusionEnergy(energySum, forceField, simulationBox, components, moleculeAtomPositions);
  }

  return energySum;
}

void Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix(
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    const SimulationBox& simulationBox, std::span<const Atom> atoms, std::span<double> potentialMatrix)
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
  const double cutOff = 0.5 * std::min({perpendicularWidths.x, perpendicularWidths.y, perpendicularWidths.z});
  const double eps = 1.0e-8;
  const double tol = std::sqrt(std::abs(std::log(eps * cutOff)));
  const double alpha = std::sqrt(std::abs(std::log(eps * cutOff * tol))) / cutOff;
  const double tol1 = std::sqrt(-std::log(eps * cutOff * (2.0 * tol * alpha) * (2.0 * tol * alpha)));

  const std::size_t kx_max_unsigned =
      static_cast<std::size_t>(std::rint(0.25 + perpendicularWidths.x * alpha * tol1 / std::numbers::pi));
  const std::size_t ky_max_unsigned =
      static_cast<std::size_t>(std::rint(0.25 + perpendicularWidths.y * alpha * tol1 / std::numbers::pi));
  const std::size_t kz_max_unsigned =
      static_cast<std::size_t>(std::rint(0.25 + perpendicularWidths.z * alpha * tol1 / std::numbers::pi));
  const std::size_t recip_integer_cutoff_squared = std::max({kx_max_unsigned, ky_max_unsigned, kz_max_unsigned}) *
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

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, atoms, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned, inv_box);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

          double prefactor =
              factor * (4.0 * std::numbers::pi / volume) * std::exp((-0.25 / (alpha * alpha)) * rksq) / rksq;

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            const double re_i = prefactor * phase[i].real();
            const double im_i = prefactor * phase[i].imag();
            double* row = &potentialMatrix[i * numberOfAtoms];
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
