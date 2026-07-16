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
import coulomb_potential;
import forcefield;
import simulationbox;
import minimization_dof_layout;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;
import minimization_cell_layout;
import potential_coulomb_real_space;
import interactions_ewald_kvector;

namespace
{
struct EwaldSite
{
  std::size_t moleculeIndex{};
  std::size_t localAtom{};
  bool rigid{};
  double3 internalOffset{};                      // pos - com for rigid molecules, zero otherwise
  std::optional<std::size_t> positionBase{};     // flexible atom xyz or rigid center-of-mass base DOF
  std::optional<std::size_t> orientationBase{};  // rigid orientation base DOF
  const Minimization::RigidAtomDerivatives* derivatives{};
};

/** DOF index and phase projection P = d(k.r)/dtheta for one site, at most 3 position + 3 orientation. */
struct SiteProjection
{
  std::array<std::size_t, 6> dof{};
  std::array<double, 6> projection{};
  std::array<double3, 6> direction{};
  std::size_t count{};
  std::size_t orientationStart{};  // entries [orientationStart, count) are orientation DOFs
};

SiteProjection buildSiteProjection(const EwaldSite& site, const double3& rk)
{
  SiteProjection result{};
  if (site.positionBase)
  {
    result.dof[result.count] = *site.positionBase + 0;
    result.projection[result.count] = rk.x;
    result.direction[result.count++] = double3(1.0, 0.0, 0.0);
    result.dof[result.count] = *site.positionBase + 1;
    result.projection[result.count] = rk.y;
    result.direction[result.count++] = double3(0.0, 1.0, 0.0);
    result.dof[result.count] = *site.positionBase + 2;
    result.projection[result.count] = rk.z;
    result.direction[result.count++] = double3(0.0, 0.0, 1.0);
  }
  result.orientationStart = result.count;
  if (site.orientationBase && site.derivatives != nullptr)
  {
    const std::array<double3, 3> dVec = {site.derivatives->dVecX, site.derivatives->dVecY, site.derivatives->dVecZ};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      result.dof[result.count] = *site.orientationBase + axis;
      result.projection[result.count] = double3::dot(rk, dVec[axis]);
      result.direction[result.count++] = dVec[axis];
    }
  }
  return result;
}

double realInner(const std::complex<double>& lhs, const std::complex<double>& rhs)
{
  return lhs.real() * rhs.real() + lhs.imag() * rhs.imag();
}

// Intra-molecular exclusion / completion Hessian for the finite-cutoff shifted charge methods (Wolf,
// damped-shifted-force, modified-shifted-force, zero-dipole). These have no reciprocal space, so this mirrors
// the erf-based Ewald exclusion block but uses the shifted real-space potential completion
// U(r) = q_i q_j (V(r) - 1/r). In RASPA's factor convention (f1 = U'/r, f2 = (U'' - U'/r)/r^2) the shifted
// potential contributes factors.firstDerivativeFactor / .secondDerivativeFactor, and the bare -1/r term adds
// +1/r^3 to f1 and -3/r^5 to f2. Only intra-molecular adsorbate pairs inside the Coulomb cutoff are treated,
// matching the energy/gradient/strain paths. Rigid molecules contribute energy only: the term is invariant
// under their center-of-mass and orientation degrees of freedom. The self-energy is position independent.
void addShiftedExclusionHessian(RunningEnergy& energySum, const ForceField& forceField,
                                const SimulationBox& simulationBox, const std::optional<Framework>& framework,
                                std::span<const Molecule> moleculeData, const MinimizationDofLayout& layout,
                                GeneralizedHessian& hessian, std::span<const Atom> frameworkAtoms,
                                std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
                                std::span<AtomDynamics> frameworkDynamics, const CellMinimizationLayout& cellLayout)
{
  if (forceField.omitInterInteractions) return;

  const bool flexibleFramework =
      framework.has_value() && !framework->rigid && layout.numberOfFrameworkAtoms() == frameworkAtoms.size();
  if (moleculeAtoms.empty() && !flexibleFramework) return;

  // Per-atom self-energy (position and strain independent). For a flexible framework the framework charges also
  // carry a self term, mirroring the Ewald self loop over the combined framework + molecule atom list.
  const double selfPrefactor = Units::CoulombicConversionFactor * Potentials::coulombSelfEnergyPrefactor(forceField);
  for (const Atom& atom : moleculeAtoms)
  {
    const double scaledCharge = atom.scalingCoulomb * atom.charge;
    energySum.ewald_self += selfPrefactor * scaledCharge * scaledCharge;
  }
  if (flexibleFramework)
  {
    for (const Atom& atom : frameworkAtoms)
    {
      const double scaledCharge = atom.scalingCoulomb * atom.charge;
      energySum.ewald_self += selfPrefactor * scaledCharge * scaledCharge;
    }
  }

  const bool computeStrain = (hessian.numStrain() == 1);
  const std::size_t numberOfCellDofs = layout.numberOfCellDofs();
  const double cutOffSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  const auto addFlexiblePairCellDerivatives =
      [&](std::size_t baseI, std::size_t baseJ, double f1, double f2, const double3& dr)
  {
    for (std::size_t a = 0; a < numberOfCellDofs; ++a)
    {
      const double3 displacementA = cellLayout.bases[a] * dr;
      const double3 mixedA =
          f1 * displacementA + f2 * double3::dot(dr, displacementA) * dr + cellLayout.bases[a] * (f1 * dr);
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double value = (&mixedA.x)[axis];
        hessian.add(baseI + axis, *layout.cellDof(a), value);
        hessian.add(*layout.cellDof(a), baseI + axis, value);
        hessian.add(baseJ + axis, *layout.cellDof(a), -value);
        hessian.add(*layout.cellDof(a), baseJ + axis, -value);
      }
      for (std::size_t b = 0; b < numberOfCellDofs; ++b)
      {
        const double3 displacementB = cellLayout.bases[b] * dr;
        const double3 secondDisplacement = cellStrainSecondDerivative(cellLayout, a, b) * dr;
        const double value = f1 * double3::dot(displacementA, displacementB) +
                             f2 * double3::dot(dr, displacementA) * double3::dot(dr, displacementB) +
                             f1 * double3::dot(dr, secondDisplacement);
        hessian.add(*layout.cellDof(a), *layout.cellDof(b), value);
      }
    }
  };

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (molecule.numberOfAtoms < 2) continue;
    const bool rigid = layout.molecules()[moleculeIndex].rigid;
    std::span<const Atom> span = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
    std::span<AtomDynamics> dynamicsSpan = moleculeDynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms);

    for (std::size_t i = 0; i != span.size() - 1; ++i)
    {
      for (std::size_t j = i + 1; j != span.size(); ++j)
      {
        double3 dr = simulationBox.applyPeriodicBoundaryConditions(span[i].position - span[j].position);
        const double rr = double3::dot(dr, dr);
        if (rr >= cutOffSquared) continue;
        const double r = std::sqrt(rr);

        const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
        const double chargeProduct = Units::CoulombicConversionFactor * span[i].scalingCoulomb *
                                     span[j].scalingCoulomb * span[i].charge * span[j].charge;

        energySum.ewald_exclusion += chargeProduct * (factors.potential - 1.0 / r);

        if (rigid) continue;

        const double f1 = chargeProduct * (factors.firstDerivativeFactor + 1.0 / (rr * r));
        const double f2 = chargeProduct * (factors.secondDerivativeFactor - 3.0 / (rr * rr * r));

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
        if (numberOfCellDofs != 0)
        {
          addFlexiblePairCellDerivatives(layout.flexibleAtomDofBase(moleculeIndex, i),
                                         layout.flexibleAtomDofBase(moleculeIndex, j), f1, f2, dr);
        }
        if (computeStrain)
        {
          Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, i, moleculeIndex, j, f1, f2,
                                                             dr);
          Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, span[i].position, span[i].position,
                                                           span[j].position, span[j].position, false, false);
        }
      }
    }
  }

  // Flexible-framework intra-molecular exclusion / completion (bonded 1-2, 1-3, 1-4 pairs excluded from the
  // real-space shifted pair sum). Only the framework atomic degrees of freedom couple; mirrors the erf-based
  // Ewald framework exclusion block but uses the shifted completion q_i q_j (V(r) - 1/r).
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

      const double f1 = chargeProduct * (factors.firstDerivativeFactor + 1.0 / (rr * r));
      const double f2 = chargeProduct * (factors.secondDerivativeFactor - 3.0 / (rr * rr * r));

      const double3 gradientI = f1 * dr;
      if (frameworkDynamics.size() == frameworkAtoms.size())
      {
        frameworkDynamics[i].gradient += gradientI;
        frameworkDynamics[j].gradient -= gradientI;
      }

      double3x3 strain{};
      strain.ax = dr.x * gradientI.x;
      strain.bx = dr.y * gradientI.x;
      strain.cx = dr.z * gradientI.x;
      strain.ay = dr.x * gradientI.y;
      strain.by = dr.y * gradientI.y;
      strain.cy = dr.z * gradientI.y;
      strain.az = dr.x * gradientI.z;
      strain.bz = dr.y * gradientI.z;
      strain.cz = dr.z * gradientI.z;
      hessian.strainGradient() += strain;

      Minimization::scatterAtomicPositionPositionByDof(hessian, *layout.frameworkAtomDof(i, MinimizationDofAxis::X),
                                                       *layout.frameworkAtomDof(j, MinimizationDofAxis::X), f1, f2, dr);
      if (numberOfCellDofs != 0)
      {
        addFlexiblePairCellDerivatives(*layout.frameworkAtomDof(i, MinimizationDofAxis::X),
                                       *layout.frameworkAtomDof(j, MinimizationDofAxis::X), f1, f2, dr);
      }
    }
  }
}
}  // namespace

RunningEnergy Interactions::computeEwaldFourierHessian(
    const ForceField& forceField, const SimulationBox& simulationBox, const std::optional<Framework>& framework,
    std::span<const std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> fixedFrameworkStoredEik,
    double netChargeFramework, std::span<const Molecule> moleculeData, std::span<const Component> components,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, const MinimizationDofLayout& layout,
    GeneralizedHessian& hessian, std::span<AtomDynamics> moleculeDynamics, std::span<AtomDynamics> frameworkDynamics,
    const CellMinimizationLayout& cellLayout)
{
  RunningEnergy energySum{};

  if (!forceField.useCharge) return energySum;
  if (!forceField.usesEwaldFourier())
  {
    // Finite-cutoff shifted methods (Wolf, DSF, mDSF, zero-dipole) have no reciprocal space. Add the
    // position-independent self-energy and the intra-molecular exclusion / completion Hessian so the
    // assembled second derivative is consistent with the gradient for multi-site molecules. Plain Coulomb
    // and the Ewald real-space-only debugging mode (omitted Fourier) require no such corrections.
    if (forceField.usesRealSpaceChargeCorrections())
    {
      addShiftedExclusionHessian(energySum, forceField, simulationBox, framework, moleculeData, layout, hessian,
                                 frameworkAtoms, moleculeAtoms, moleculeDynamics, frameworkDynamics, cellLayout);
    }
    return energySum;
  }

  const bool flexibleFramework =
      framework && !framework->rigid && layout.numberOfFrameworkAtoms() == frameworkAtoms.size();
  const std::size_t frameworkOffset = flexibleFramework ? frameworkAtoms.size() : 0;
  // For a flexible framework the Fourier sums run over the framework and molecule atoms as one
  // indexed array; assemble that view locally so the spans need not be contiguous in memory.
  std::vector<Atom> combinedAtoms{};
  if (flexibleFramework)
  {
    combinedAtoms.reserve(frameworkAtoms.size() + moleculeAtoms.size());
    combinedAtoms.insert(combinedAtoms.end(), frameworkAtoms.begin(), frameworkAtoms.end());
    combinedAtoms.insert(combinedAtoms.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  }
  std::span<const Atom> atoms = flexibleFramework ? std::span<const Atom>(combinedAtoms) : moleculeAtoms;
  const std::size_t numberOfAtoms = atoms.size();
  if (numberOfAtoms == 0) return energySum;

  const double alpha = forceField.EwaldAlpha;
  const double alpha_squared = alpha * alpha;
  const std::size_t recip_integer_cutoff_squared = forceField.reciprocalIntegerCutOffSquared;
  const double recip_cutoff_squared = forceField.reciprocalCutOffSquared;
  const bool omitInterInteractions = forceField.omitInterInteractions;
  const bool computeStrain = (hessian.numStrain() == 1);
  const std::size_t numberOfCellDofs = layout.numberOfCellDofs();
  if (numberOfCellDofs != cellLayout.size())
  {
    throw std::invalid_argument("computeEwaldFourierHessian: cell layout does not match minimization DOFs");
  }

  const double3x3 inv_box = simulationBox.inverseCell;
  const double3 ax = double3(inv_box.ax, inv_box.bx, inv_box.cx);
  const double3 ay = double3(inv_box.ay, inv_box.by, inv_box.cy);
  const double3 az = double3(inv_box.az, inv_box.bz, inv_box.cz);

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  // Per-atom site metadata in global atom order.
  std::vector<EwaldSite> sites(numberOfAtoms);
  if (flexibleFramework)
  {
    for (std::size_t atom = 0; atom < frameworkOffset; ++atom)
    {
      sites[atom].positionBase = layout.frameworkAtomDof(atom, MinimizationDofAxis::X);
    }
  }
  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    const bool rigid = layout.molecules()[moleculeIndex].rigid;
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      EwaldSite& site = sites[frameworkOffset + molecule.atomIndex + localAtom];
      site.moleculeIndex = moleculeIndex;
      site.localAtom = localAtom;
      site.rigid = rigid;
      if (rigid)
      {
        site.internalOffset =
            atoms[frameworkOffset + molecule.atomIndex + localAtom].position - molecule.centerOfMassPosition;
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
  const auto addGradient = [&](std::size_t atom, const double3& gradient)
  {
    if (flexibleFramework && atom < frameworkOffset)
    {
      if (frameworkDynamics.size() == frameworkOffset) frameworkDynamics[atom].gradient += gradient;
    }
    else
    {
      const std::size_t moleculeAtom = atom - frameworkOffset;
      if (moleculeDynamics.size() == moleculeAtoms.size())
      {
        moleculeDynamics[moleculeAtom].gradient += gradient;
      }
    }
  };
  const auto addFlexiblePairCellDerivatives =
      [&](std::size_t baseI, std::size_t baseJ, double f1, double f2, const double3& dr)
  {
    for (std::size_t a = 0; a < numberOfCellDofs; ++a)
    {
      const double3 displacementA = cellLayout.bases[a] * dr;
      const double3 mixedA =
          f1 * displacementA + f2 * double3::dot(dr, displacementA) * dr + cellLayout.bases[a] * (f1 * dr);
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double value = (&mixedA.x)[axis];
        hessian.add(baseI + axis, *layout.cellDof(a), value);
        hessian.add(*layout.cellDof(a), baseI + axis, value);
        hessian.add(baseJ + axis, *layout.cellDof(a), -value);
        hessian.add(*layout.cellDof(a), baseJ + axis, -value);
      }
      for (std::size_t b = 0; b < numberOfCellDofs; ++b)
      {
        const double3 displacementB = cellLayout.bases[b] * dr;
        const double3 secondDisplacement = cellStrainSecondDerivative(cellLayout, a, b) * dr;
        const double value = f1 * double3::dot(displacementA, displacementB) +
                             f2 * double3::dot(dr, displacementA) * double3::dot(dr, displacementB) +
                             f1 * double3::dot(dr, secondDisplacement);
        hessian.add(*layout.cellDof(a), *layout.cellDof(b), value);
      }
    }
  };

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

  Ewald::buildEikTables(eik_x, eik_y, eik_z, eik_xy, atoms, kx_max_unsigned, ky_max_unsigned, kz_max_unsigned, inv_box);

  const std::size_t numDofs = layout.numDofs();
  std::vector<std::complex<double>> dofPhase(numDofs);  // V_a = sum_i P_i^a e_i
  std::vector<double> positionStrainScratch(computeStrain ? numDofs : 0);
  std::vector<std::complex<double>> cellStructureFirst(numberOfCellDofs);
  std::vector<std::complex<double>> cellDynamicFirst(numberOfCellDofs);
  std::vector<std::complex<double>> cellStructureSecond(numberOfCellDofs * numberOfCellDofs);
  std::vector<std::complex<double>> cellDynamicSecond(numberOfCellDofs * numberOfCellDofs);
  std::vector<std::complex<double>> dofPhaseCell(numberOfCellDofs * numDofs);
  std::vector<double> singleIonCellCell(numberOfCellDofs * numberOfCellDofs);

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
      Ewald::fillEikXYRow(eik_xy, eik_x, eik_y, numberOfAtoms, kx, ky);

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

          const std::complex<double> rigidSF = (!flexibleFramework && nvec < fixedFrameworkStoredEik.size())
                                                   ? fixedFrameworkStoredEik[nvec].first
                                                   : std::complex<double>(0.0, 0.0);

          std::complex<double> total = rigidSF + cksum;
          const std::complex<double> fullTotal = total;

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
          singleIonStrainStrain -= temp * (traceTheta * traceTheta - rksq / alpha_squared);

          std::ranges::fill(dofPhase, std::complex<double>(0.0, 0.0));
          std::ranges::fill(cellStructureFirst, std::complex<double>(0.0, 0.0));
          std::ranges::fill(cellDynamicFirst, std::complex<double>(0.0, 0.0));
          std::ranges::fill(cellStructureSecond, std::complex<double>(0.0, 0.0));
          std::ranges::fill(cellDynamicSecond, std::complex<double>(0.0, 0.0));
          std::ranges::fill(dofPhaseCell, std::complex<double>(0.0, 0.0));
          if (computeStrain)
          {
            std::ranges::fill(positionStrainScratch, 0.0);
          }

          std::complex<double> rigidPhaseSum(0.0, 0.0);  // W = sum_i (k.d_i) e_i, rigid sites only
          double strainStrainSites = 0.0;

          for (std::size_t i = 0; i != numberOfAtoms; ++i)
          {
            const EwaldSite& site = sites[i];
            const std::complex<double> e = eikr[i];

            // f1 = 2 factor Re[i e conj(S)]; f2 = -2 factor Re[e conj(S)].
            const double f1 = 2.0 * temp * (e.real() * total.imag() - e.imag() * total.real());
            const double f2 = -2.0 * temp * (e.real() * total.real() + e.imag() * total.imag());

            addGradient(i, f1 * rk);

            const SiteProjection proj = buildSiteProjection(site, rk);
            for (std::size_t a = 0; a < proj.count; ++a)
            {
              dofPhase[proj.dof[a]] += proj.projection[a] * e;
            }

            if (numberOfCellDofs != 0)
            {
              std::array<double, 6> phaseFirst{};
              for (std::size_t a = 0; a < numberOfCellDofs; ++a)
              {
                const double3 transformedWaveVector = cellLayout.bases[a] * rk;
                phaseFirst[a] = site.rigid ? -double3::dot(transformedWaveVector, site.internalOffset) : 0.0;
                if (site.rigid)
                {
                  const std::complex<double> first = std::complex<double>(0.0, phaseFirst[a]) * e;
                  cellStructureFirst[a] += first;
                  cellDynamicFirst[a] += first;
                }

                for (std::size_t p = 0; p < proj.count; ++p)
                {
                  const double projectionFirst =
                      p < proj.orientationStart ? 0.0 : -double3::dot(transformedWaveVector, proj.direction[p]);
                  dofPhaseCell[a * numDofs + proj.dof[p]] +=
                      (std::complex<double>(0.0, projectionFirst) - proj.projection[p] * phaseFirst[a]) * e;
                }
              }
              if (site.rigid)
              {
                for (std::size_t a = 0; a < numberOfCellDofs; ++a)
                {
                  for (std::size_t b = 0; b < numberOfCellDofs; ++b)
                  {
                    const double phaseSecond =
                        double3::dot(cellStrainSecondDerivative(cellLayout, a, b) * rk, site.internalOffset);
                    const std::complex<double> second =
                        (std::complex<double>(0.0, phaseSecond) - phaseFirst[a] * phaseFirst[b]) * e;
                    cellStructureSecond[a * numberOfCellDofs + b] += second;
                    cellDynamicSecond[a * numberOfCellDofs + b] += second;
                  }
                }
              }
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
              const Minimization::RigidAtomDerivatives& derivatives = *site.derivatives;
              const std::array<std::array<double3, 3>, 3> ddVec = {
                  {{derivatives.ddVecAX, derivatives.ddVecAY, derivatives.ddVecAZ},
                   {derivatives.ddVecAY, derivatives.ddVecBY, derivatives.ddVecBZ},
                   {derivatives.ddVecAZ, derivatives.ddVecBZ, derivatives.ddVecCZ}}};
              for (std::size_t a = 0; a < 3; ++a)
              {
                for (std::size_t b = 0; b < 3; ++b)
                {
                  hessian.add(*site.orientationBase + a, *site.orientationBase + b, f1 * double3::dot(rk, ddVec[a][b]));
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
          // It cancels against the subtracted |C|^2 channel when intermolecular interactions are omitted.
          if (!omitInterInteractions)
          {
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
          }

          if (numberOfCellDofs != 0)
          {
            const double c = 0.25 / alpha_squared;
            std::array<double, 6> radiusFirst{};
            std::array<double, 6> logTempFirst{};
            std::array<double, 6> amplitudeFirst{};
            for (std::size_t a = 0; a < numberOfCellDofs; ++a)
            {
              radiusFirst[a] = -2.0 * double3::dot(rk, cellLayout.bases[a] * rk);
              logTempFirst[a] = -cellLayout.bases[a].trace() + (-c - 1.0 / rksq) * radiusFirst[a];
              amplitudeFirst[a] = 2.0 * realInner(cellStructureFirst[a], fullTotal) -
                                  (omitInterInteractions ? 2.0 * realInner(cellDynamicFirst[a], cksum) : 0.0);
            }

            const double amplitude =
                std::norm(fullTotal) - std::norm(rigidSF) - (omitInterInteractions ? std::norm(cksum) : 0.0);
            for (std::size_t a = 0; a < numberOfCellDofs; ++a)
            {
              for (std::size_t b = 0; b < numberOfCellDofs; ++b)
              {
                const double3x3 secondBasis = cellStrainSecondDerivative(cellLayout, a, b);
                const double radiusSecond = 4.0 * double3::dot(rk, secondBasis * rk);
                const double logTempSecond =
                    (-c - 1.0 / rksq) * radiusSecond + radiusFirst[a] * radiusFirst[b] / (rksq * rksq);
                const double tempFirstA = temp * logTempFirst[a];
                const double tempFirstB = temp * logTempFirst[b];
                const double tempSecond = temp * (logTempFirst[a] * logTempFirst[b] + logTempSecond);
                const std::size_t ab = a * numberOfCellDofs + b;
                const double amplitudeSecond =
                    2.0 * realInner(cellStructureSecond[ab], fullTotal) +
                    2.0 * realInner(cellStructureFirst[a], cellStructureFirst[b]) -
                    (omitInterInteractions ? 2.0 * realInner(cellDynamicSecond[ab], cksum) +
                                                 2.0 * realInner(cellDynamicFirst[a], cellDynamicFirst[b])
                                           : 0.0);
                const double cellCell = tempSecond * amplitude + tempFirstA * amplitudeFirst[b] +
                                        tempFirstB * amplitudeFirst[a] + temp * amplitudeSecond;
                hessian.add(*layout.cellDof(a), *layout.cellDof(b), cellCell);
                singleIonCellCell[ab] -= tempSecond;
              }
            }

            for (std::size_t p = 0; p < layout.numberOfPositionDofs(); ++p)
            {
              const std::complex<double> phasePosition = std::complex<double>(0.0, 1.0) * dofPhase[p];
              const double amplitudePosition = 2.0 * realInner(phasePosition, fullTotal) -
                                               (omitInterInteractions ? 2.0 * realInner(phasePosition, cksum) : 0.0);
              for (std::size_t a = 0; a < numberOfCellDofs; ++a)
              {
                const std::complex<double> phasePositionCell = dofPhaseCell[a * numDofs + p];
                const double amplitudePositionCell =
                    2.0 * realInner(phasePositionCell, fullTotal) +
                    2.0 * realInner(phasePosition, cellStructureFirst[a]) -
                    (omitInterInteractions ? 2.0 * realInner(phasePositionCell, cksum) +
                                                 2.0 * realInner(phasePosition, cellDynamicFirst[a])
                                           : 0.0);
                const double mixed = temp * (logTempFirst[a] * amplitudePosition + amplitudePositionCell);
                hessian.add(p, *layout.cellDof(a), mixed);
                hessian.add(*layout.cellDof(a), p, mixed);
              }
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
    for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
    {
      const Molecule& molecule = moleculeData[moleculeIndex];
      if (molecule.numberOfAtoms < 2)
      {
        continue;
      }
      const bool rigid = layout.molecules()[moleculeIndex].rigid;
      std::span<const Atom> span = atoms.subspan(frameworkOffset + molecule.atomIndex, molecule.numberOfAtoms);
      std::span<AtomDynamics> dynamicsSpan = moleculeDynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms);

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
          const double f2 = chargeProduct * (2.0 * alpha_squared * gaussTerm / rr + 3.0 * gaussTerm / (rr * rr) -
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
          if (numberOfCellDofs != 0)
          {
            addFlexiblePairCellDerivatives(layout.flexibleAtomDofBase(moleculeIndex, i),
                                           layout.flexibleAtomDofBase(moleculeIndex, j), f1, f2, dr);
          }
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

    if (flexibleFramework)
    {
      const std::span<const Atom> frameworkAtoms = atoms.first(frameworkOffset);
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
          const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
          excludeIfAbsentOrScaled(pair);
        }
      }
      for (const std::array<std::size_t, 2>& pair : excludedPairs)
      {
        const std::size_t i = pair[0];
        const std::size_t j = pair[1];
        double3 dr =
            simulationBox.applyPeriodicBoundaryConditions(frameworkAtoms[i].position - frameworkAtoms[j].position);
        const double rr = double3::dot(dr, dr);
        const double r = std::sqrt(rr);
        const double chargeProduct = Units::CoulombicConversionFactor * frameworkAtoms[i].scalingCoulomb *
                                     frameworkAtoms[j].scalingCoulomb * frameworkAtoms[i].charge *
                                     frameworkAtoms[j].charge;
        const double erfTerm = std::erf(alpha * r);
        const double gaussTerm = 2.0 * alpha * std::numbers::inv_sqrtpi * std::exp(-alpha_squared * rr);
        energySum.ewald_exclusion -= chargeProduct * erfTerm / r;
        const double f1 = -chargeProduct * (gaussTerm / rr - erfTerm / (r * rr));
        const double f2 = chargeProduct * (2.0 * alpha_squared * gaussTerm / rr + 3.0 * gaussTerm / (rr * rr) -
                                           3.0 * erfTerm / (r * rr * rr));
        const double3 gradientI = f1 * dr;
        addGradient(i, gradientI);
        addGradient(j, -gradientI);
        double3x3 strain{};
        strain.ax = dr.x * gradientI.x;
        strain.bx = dr.y * gradientI.x;
        strain.cx = dr.z * gradientI.x;
        strain.ay = dr.x * gradientI.y;
        strain.by = dr.y * gradientI.y;
        strain.cy = dr.z * gradientI.y;
        strain.az = dr.x * gradientI.z;
        strain.bz = dr.y * gradientI.z;
        strain.cz = dr.z * gradientI.z;
        hessian.strainGradient() += strain;
        Minimization::scatterAtomicPositionPositionByDof(hessian, *layout.frameworkAtomDof(i, MinimizationDofAxis::X),
                                                         *layout.frameworkAtomDof(j, MinimizationDofAxis::X), f1, f2,
                                                         dr);
        if (numberOfCellDofs != 0)
        {
          addFlexiblePairCellDerivatives(*layout.frameworkAtomDof(i, MinimizationDofAxis::X),
                                         *layout.frameworkAtomDof(j, MinimizationDofAxis::X), f1, f2, dr);
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
    const double rigidFrameworkCharge = flexibleFramework ? 0.0 : netChargeFramework;
    if (omitInterInteractions)
    {
      netChargeFactor = 2.0 * rigidFrameworkCharge * netChargeAdsorbates;
    }
    else
    {
      netChargeFactor = (2.0 * rigidFrameworkCharge + netChargeAdsorbates) * netChargeAdsorbates;
    }
    energySum.ewald_fourier += uIon * netChargeFactor;
    hessian.strainGradient() += netChargeFactor * singleIonStrainGradient;
    for (std::size_t a = 0; a < numberOfCellDofs; ++a)
    {
      for (std::size_t b = 0; b < numberOfCellDofs; ++b)
      {
        hessian.add(*layout.cellDof(a), *layout.cellDof(b),
                    netChargeFactor * singleIonCellCell[a * numberOfCellDofs + b]);
      }
    }
    if (computeStrain)
    {
      hessian.addStrainStrain(0, 0, netChargeFactor * singleIonStrainStrain);
    }
  }

  return energySum;
}
