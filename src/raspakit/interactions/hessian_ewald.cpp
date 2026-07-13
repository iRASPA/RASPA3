module;

module interactions_hessian_ewald;

import std;

import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import molecule;
import component;
import forcefield;
import simulationbox;
import minimization_dof_layout;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;

namespace
{
struct EwaldSite
{
  std::size_t moleculeIndex{};
  std::size_t localAtom{};
  bool rigid{};
  double3 internalOffset{};  // pos - com for rigid molecules, zero otherwise
  std::optional<std::size_t> positionBase{};     // flexible atom xyz or rigid center-of-mass base DOF
  std::optional<std::size_t> orientationBase{};  // rigid orientation base DOF
  const Minimization::RigidAtomDerivatives *derivatives{};
};

/** DOF index and phase projection P = d(k.r)/dtheta for one site, at most 3 position + 3 orientation. */
struct SiteProjection
{
  std::array<std::size_t, 6> dof{};
  std::array<double, 6> projection{};
  std::size_t count{};
  std::size_t orientationStart{};  // entries [orientationStart, count) are orientation DOFs
};

SiteProjection buildSiteProjection(const EwaldSite &site, const double3 &rk)
{
  SiteProjection result{};
  if (site.positionBase)
  {
    result.dof[result.count] = *site.positionBase + 0;
    result.projection[result.count++] = rk.x;
    result.dof[result.count] = *site.positionBase + 1;
    result.projection[result.count++] = rk.y;
    result.dof[result.count] = *site.positionBase + 2;
    result.projection[result.count++] = rk.z;
  }
  result.orientationStart = result.count;
  if (site.orientationBase && site.derivatives != nullptr)
  {
    const std::array<double3, 3> dVec = {site.derivatives->dVecX, site.derivatives->dVecY, site.derivatives->dVecZ};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      result.dof[result.count] = *site.orientationBase + axis;
      result.projection[result.count++] = double3::dot(rk, dVec[axis]);
    }
  }
  return result;
}
}  // namespace

RunningEnergy Interactions::computeEwaldFourierHessian(const System &system, const MinimizationDofLayout &layout,
                                                       GeneralizedHessian &hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energySum{};

  const ForceField &forceField = system.forceField;
  if (!forceField.useCharge) return energySum;
  if (forceField.omitEwaldFourier) return energySum;

  const SimulationBox &simulationBox = system.simulationBox;
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::size_t numberOfAtoms = atoms.size();
  if (numberOfAtoms == 0) return energySum;

  const double alpha = forceField.EwaldAlpha;
  const double alpha_squared = alpha * alpha;
  const std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  const double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  const bool omitInterInteractions = forceField.omitInterInteractions;
  const bool computeStrain = (hessian.numStrain() == 1);

  const double3x3 inv_box = simulationBox.inverseCell;
  const double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  const double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  const double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, atoms);

  // Per-atom site metadata in global atom order.
  std::vector<EwaldSite> sites(numberOfAtoms);
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = system.moleculeData[moleculeIndex];
    const bool rigid = layout.molecules()[moleculeIndex].rigid;
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      EwaldSite &site = sites[molecule.atomIndex + localAtom];
      site.moleculeIndex = moleculeIndex;
      site.localAtom = localAtom;
      site.rigid = rigid;
      if (rigid)
      {
        site.internalOffset = atoms[molecule.atomIndex + localAtom].position - molecule.centerOfMassPosition;
        site.positionBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
        site.orientationBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
        site.derivatives = &rigidCache.atom(moleculeIndex, localAtom);
      }
      else
      {
        site.positionBase = layout.flexibleAtomDof(moleculeIndex, localAtom, MinimizationDofAxis::X);
      }
    }
  }

  const std::size_t kx_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.x);
  const std::size_t ky_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.y);
  const std::size_t kz_max_unsigned = static_cast<std::size_t>(forceField.numberOfWaveVectors.z);

  const std::make_signed_t<std::size_t> kx_max = static_cast<std::make_signed_t<std::size_t>>(kx_max_unsigned);
  const std::make_signed_t<std::size_t> ky_max = static_cast<std::make_signed_t<std::size_t>>(ky_max_unsigned);
  const std::make_signed_t<std::size_t> kz_max = static_cast<std::make_signed_t<std::size_t>>(kz_max_unsigned);

  std::vector<std::complex<double>> eik_x(numberOfAtoms * (kx_max_unsigned + 1));
  std::vector<std::complex<double>> eik_y(numberOfAtoms * (ky_max_unsigned + 1));
  std::vector<std::complex<double>> eik_z(numberOfAtoms * (kz_max_unsigned + 1));
  std::vector<std::complex<double>> eik_xy(numberOfAtoms);
  std::vector<std::complex<double>> eikr(numberOfAtoms);

  // Construct exp(ik.r) for atoms and k-vectors kx, ky, kz = 0, 1 explicitly
  for (std::size_t i = 0; i != numberOfAtoms; ++i)
  {
    eik_x[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_y[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    eik_z[i + 0 * numberOfAtoms] = std::complex<double>(1.0, 0.0);
    const double3 s = 2.0 * std::numbers::pi * (inv_box * atoms[i].position);
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

  const std::size_t numDofs = layout.numDofs();
  std::vector<std::complex<double>> dofPhase(numDofs);          // V_a = sum_i P_i^a e_i
  std::vector<double> positionStrainScratch(computeStrain ? numDofs : 0);

  std::size_t nvec = 0;
  double singleIonFourierSum = 0.0;
  double3x3 singleIonStrainGradient{};
  double singleIonStrainStrain = 0.0;
  const double prefactor = Units::CoulombicConversionFactor * (2.0 * std::numbers::pi / simulationBox.volume);
  for (std::make_signed_t<std::size_t> kx = 0; kx <= kx_max; ++kx)
  {
    const double3 kvec_x = 2.0 * std::numbers::pi * static_cast<double>(kx) * ax;

    // Only positive kx are used, the negative kx are taken into account by the factor of two
    const double factor = (kx == 0) ? (1.0 * prefactor) : (2.0 * prefactor);

    for (std::make_signed_t<std::size_t> ky = -ky_max; ky <= ky_max; ++ky)
    {
      const double3 kvec_y = 2.0 * std::numbers::pi * static_cast<double>(ky) * ay;

      // Precompute and store eik_x * eik_y outside the kz-loop
      for (std::size_t i = 0; i != numberOfAtoms; ++i)
      {
        std::complex<double> eiky_temp = eik_y[i + numberOfAtoms * static_cast<std::size_t>(std::abs(ky))];
        eiky_temp.imag(ky >= 0 ? eiky_temp.imag() : -eiky_temp.imag());
        eik_xy[i] = eik_x[i + numberOfAtoms * static_cast<std::size_t>(kx)] * eiky_temp;
      }

      for (std::make_signed_t<std::size_t> kz = -kz_max; kz <= kz_max; ++kz)
      {
        const double3 kvec_z = 2.0 * std::numbers::pi * static_cast<double>(kz) * az;
        const double3 rk = kvec_x + kvec_y + kvec_z;
        const double rksq = rk.length_squared();

        // Omit kvec==0
        const std::size_t ksq = static_cast<std::size_t>(kx * kx + ky * ky + kz * kz);
        if ((ksq != 0uz) && (ksq <= recip_integer_cutoff_squared) && (rksq < recip_cutoff_squared))
        {
          std::complex<double> cksum(0.0, 0.0);
          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            std::complex<double> eikz_temp = eik_z[i + numberOfAtoms * static_cast<std::size_t>(std::abs(kz))];
            eikz_temp.imag(kz >= 0 ? eikz_temp.imag() : -eikz_temp.imag());
            eikr[i] = atoms[i].scalingCoulomb * atoms[i].charge * (eik_xy[i] * eikz_temp);
            cksum += eikr[i];
          }

          const std::complex<double> rigidSF =
              (nvec < system.fixedFrameworkStoredEik.size()) ? system.fixedFrameworkStoredEik[nvec].first
                                                             : std::complex<double>(0.0, 0.0);

          std::complex<double> total = rigidSF + cksum;

          const double temp = factor * std::exp((-0.25 / alpha_squared) * rksq) / rksq;
          singleIonFourierSum += temp;

          const double rigidEnergy = temp * std::norm(rigidSF);
          double energyThisWaveVector = temp * std::norm(total) - rigidEnergy;
          if (omitInterInteractions)
          {
            energyThisWaveVector -= temp * std::norm(cksum);
            total -= cksum;
          }
          energySum.ewald_fourier += energyThisWaveVector;

          // Strain prefactor matrix Theta_ab = delta_ab - 2 k_a k_b lambda (RASPA2 convention).
          const double inverseLambdaSquared = 0.25 / alpha_squared + 1.0 / rksq;
          const double traceTheta = 1.0 - rksq / (2.0 * alpha_squared);
          double3x3 theta{};
          theta.ax = 1.0 - 2.0 * rk.x * rk.x * inverseLambdaSquared;
          theta.ay = -2.0 * rk.x * rk.y * inverseLambdaSquared;
          theta.az = -2.0 * rk.x * rk.z * inverseLambdaSquared;
          theta.bx = -2.0 * rk.y * rk.x * inverseLambdaSquared;
          theta.by = 1.0 - 2.0 * rk.y * rk.y * inverseLambdaSquared;
          theta.bz = -2.0 * rk.y * rk.z * inverseLambdaSquared;
          theta.cx = -2.0 * rk.z * rk.x * inverseLambdaSquared;
          theta.cy = -2.0 * rk.z * rk.y * inverseLambdaSquared;
          theta.cz = 1.0 - 2.0 * rk.z * rk.z * inverseLambdaSquared;

          hessian.strainGradient() -= energyThisWaveVector * theta;

          // uIon = alpha/sqrt(pi) - sum_k temp.  For a cell-strain component eta_ab,
          // d(temp)/deta_ab = -temp Theta_ab.  The isotropic exp(epsilon) curvature is
          // d2(temp)/dEpsilon2 = temp[(tr Theta)^2 - k^2/alpha^2].
          singleIonStrainGradient += temp * theta;
          singleIonStrainStrain -=
              temp * (traceTheta * traceTheta - rksq / alpha_squared);

          std::ranges::fill(dofPhase, std::complex<double>(0.0, 0.0));
          if (computeStrain)
          {
            std::ranges::fill(positionStrainScratch, 0.0);
          }

          std::complex<double> rigidPhaseSum(0.0, 0.0);  // W = sum_i (k.d_i) e_i, rigid sites only
          double strainStrainSites = 0.0;

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            const EwaldSite &site = sites[i];
            const std::complex<double> e = eikr[i];

            // f1 = 2 factor Re[i e conj(S)]; f2 = -2 factor Re[e conj(S)].
            const double f1 = 2.0 * temp * (e.real() * total.imag() - e.imag() * total.real());
            const double f2 = -2.0 * temp * (e.real() * total.real() + e.imag() * total.imag());

            dynamics[i].gradient += f1 * rk;

            const SiteProjection proj = buildSiteProjection(site, rk);
            for (std::size_t a = 0; a < proj.count; ++a)
            {
              dofPhase[proj.dof[a]] += proj.projection[a] * e;
            }

            // Same-site curvature: f2 P_a P_b, plus f1 (k . d2r/domega_a domega_b) for orientations.
            for (std::size_t a = 0; a < proj.count; ++a)
            {
              for (std::size_t b = 0; b < proj.count; ++b)
              {
                hessian.add(proj.dof[a], proj.dof[b], f2 * proj.projection[a] * proj.projection[b]);
              }
            }
            if (site.orientationBase && site.derivatives != nullptr)
            {
              const Minimization::RigidAtomDerivatives &derivatives = *site.derivatives;
              const std::array<std::array<double3, 3>, 3> ddVec = {
                  {{derivatives.ddVecAX, derivatives.ddVecAY, derivatives.ddVecAZ},
                   {derivatives.ddVecAY, derivatives.ddVecBY, derivatives.ddVecBZ},
                   {derivatives.ddVecAZ, derivatives.ddVecBZ, derivatives.ddVecCZ}}};
              for (std::size_t a = 0; a < 3; ++a)
              {
                for (std::size_t b = 0; b < 3; ++b)
                {
                  hessian.add(*site.orientationBase + a, *site.orientationBase + b,
                              f1 * double3::dot(rk, ddVec[a][b]));
                }
              }
            }

            if (site.rigid)
            {
              // Rigid strain first derivative: -f1 sym(d (x) k) (second term of Eq. 41,
              // Dubbeldam, Krishna, Snurr 2009).
              const double3 d = site.internalOffset;
              double3x3 correction{};
              correction.ax = f1 * d.x * rk.x;
              correction.by = f1 * d.y * rk.y;
              correction.cz = f1 * d.z * rk.z;
              correction.ay = correction.bx = 0.5 * f1 * (d.x * rk.y + d.y * rk.x);
              correction.az = correction.cx = 0.5 * f1 * (d.x * rk.z + d.z * rk.x);
              correction.bz = correction.cy = 0.5 * f1 * (d.y * rk.z + d.z * rk.y);
              hessian.strainGradient() -= correction;
            }

            if (computeStrain)
            {
              const double kd = site.rigid ? double3::dot(rk, site.internalOffset) : 0.0;
              if (kd != 0.0)
              {
                rigidPhaseSum += kd * e;
                strainStrainSites += (2.0 * traceTheta + 1.0) * kd * f1 + kd * kd * f2;
                for (std::size_t a = 0; a < proj.count; ++a)
                {
                  positionStrainScratch[proj.dof[a]] -= kd * f2 * proj.projection[a];
                }
              }
              // Orientation projections scale with exp(-epsilon): extra -f1 P term.
              for (std::size_t a = proj.orientationStart; a < proj.count; ++a)
              {
                positionStrainScratch[proj.dof[a]] -= f1 * proj.projection[a];
              }
            }
          }

          // Pair term: 2 factor Re[dS/dtheta_a conj(dS/dtheta_b)] = 2 factor Re[V_a conj(V_b)].
          for (std::size_t a = 0; a < numDofs; ++a)
          {
            const std::complex<double> va = dofPhase[a];
            if (va.real() == 0.0 && va.imag() == 0.0)
            {
              continue;
            }
            for (std::size_t b = 0; b < numDofs; ++b)
            {
              const std::complex<double> vb = dofPhase[b];
              hessian.add(a, b, 2.0 * temp * (va.real() * vb.real() + va.imag() * vb.imag()));
            }
          }

          if (computeStrain)
          {
            // Strain-strain: prefactor curvature + rigid-site corrections.
            double strainStrain = energyThisWaveVector * (traceTheta * traceTheta - rksq / alpha_squared);
            strainStrain += strainStrainSites;
            strainStrain += 2.0 * temp * std::norm(rigidPhaseSum);
            hessian.addStrainStrain(0, 0, strainStrain);

            // Position-strain: -traceTheta * gradient - pair coupling to rigid phases + site terms.
            for (std::size_t a = 0; a < numDofs; ++a)
            {
              const std::complex<double> va = dofPhase[a];
              const double gradientProjection = 2.0 * temp * (va.real() * total.imag() - va.imag() * total.real());
              const double rigidPairCoupling =
                  2.0 * temp * (va.real() * rigidPhaseSum.real() + va.imag() * rigidPhaseSum.imag());
              const double value = -traceTheta * gradientProjection - rigidPairCoupling + positionStrainScratch[a];
              if (value != 0.0)
              {
                hessian.addPositionStrain(a, 0, value);
              }
            }
          }

          ++nvec;
        }
      }
    }
  }

  if (!omitInterInteractions)
  {
    // Subtract self-energy (independent of all degrees of freedom since alpha is fixed).
    const double prefactor_self = Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi);
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      const double scaledCharge = atoms[i].scalingCoulomb * atoms[i].charge;
      energySum.ewald_self -= prefactor_self * scaledCharge * scaledCharge;
    }

    // Subtract exclusion-energy U = -C q_i q_j erf(alpha r)/r within each molecule. For rigid
    // molecules the internal distances are constant, so only the energy contributes.
    for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
    {
      const Molecule &molecule = system.moleculeData[moleculeIndex];
      if (molecule.numberOfAtoms < 2)
      {
        continue;
      }
      const bool rigid = layout.molecules()[moleculeIndex].rigid;
      std::span<const Atom> span = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      std::span<AtomDynamics> dynamicsSpan = dynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms);

      for (std::size_t i = 0; i != span.size() - 1; ++i)
      {
        for (std::size_t j = i + 1; j != span.size(); ++j)
        {
          double3 dr = span[i].position - span[j].position;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          const double rr = double3::dot(dr, dr);
          const double r = std::sqrt(rr);

          const double chargeProduct = Units::CoulombicConversionFactor * span[i].scalingCoulomb *
                                       span[j].scalingCoulomb * span[i].charge * span[j].charge;
          const double erfTerm = std::erf(alpha * r);
          const double gaussTerm = 2.0 * alpha * std::numbers::inv_sqrtpi * std::exp(-alpha_squared * rr);

          energySum.ewald_exclusion -= chargeProduct * erfTerm / r;

          if (rigid)
          {
            continue;
          }

          // U(r) = -chargeProduct erf(alpha r)/r; RASPA convention f1 = U'/r, f2 = (U'' - U'/r)/r^2.
          const double f1 = -chargeProduct * (gaussTerm / rr - erfTerm / (r * rr));
          const double f2 =
              chargeProduct * (2.0 * alpha_squared * gaussTerm / rr + 3.0 * gaussTerm / (rr * rr) -
                               3.0 * erfTerm / (r * rr * rr));

          const double3 gradientA = f1 * dr;
          dynamicsSpan[i].gradient += gradientA;
          dynamicsSpan[j].gradient -= gradientA;

          double3x3 strain{};
          strain.ax = dr.x * gradientA.x;
          strain.bx = dr.y * gradientA.x;
          strain.cx = dr.z * gradientA.x;
          strain.ay = dr.x * gradientA.y;
          strain.by = dr.y * gradientA.y;
          strain.cy = dr.z * gradientA.y;
          strain.az = dr.x * gradientA.z;
          strain.bz = dr.y * gradientA.z;
          strain.cz = dr.z * gradientA.z;
          hessian.strainGradient() += strain;

          Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, i, moleculeIndex, j, f1, f2, dr);
          if (computeStrain)
          {
            Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, i, moleculeIndex, j, f1,
                                                               f2, dr);
            Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, span[i].position, span[i].position,
                                                             span[j].position, span[j].position, false, false);
          }
        }
      }
    }
  }

  // Net-charge correction (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)). It is independent
  // of positions, so H_pp and H_pEpsilon vanish. Its cell dependence enters entirely through
  // uIon and contributes to the strain gradient and isotropic exp(epsilon) curvature.
  {
    double netChargeAdsorbates = 0.0;
    for (std::size_t i = 0; i != numberOfAtoms; ++i)
    {
      netChargeAdsorbates += atoms[i].scalingCoulomb * atoms[i].charge;
    }
    const double uIon = -(singleIonFourierSum - Units::CoulombicConversionFactor * alpha / std::sqrt(std::numbers::pi));
    double netChargeFactor = 0.0;
    if (omitInterInteractions)
    {
      netChargeFactor = 2.0 * system.netChargeFramework * netChargeAdsorbates;
    }
    else
    {
      netChargeFactor = (2.0 * system.netChargeFramework + netChargeAdsorbates) * netChargeAdsorbates;
    }
    energySum.ewald_fourier += uIon * netChargeFactor;
    hessian.strainGradient() += netChargeFactor * singleIonStrainGradient;
    if (computeStrain)
    {
      hessian.addStrainStrain(0, 0, netChargeFactor * singleIonStrainStrain);
    }
  }

  return energySum;
}
