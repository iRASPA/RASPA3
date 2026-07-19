module;

module interactions_hessian_framework_molecule;

import std;

import double3;
import double3x3;
import atom;
import molecule;
import framework;
import forcefield;
import simulationbox;
import interactions_pair_kernel;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;
import minimization_cell_layout;

namespace
{
void addFrameworkMoleculeCellCorrection(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                        const Minimization::RigidDerivativeCache& rigidCache,
                                        const CellMinimizationLayout& cellLayout, std::size_t frameworkAtom,
                                        bool flexibleFramework, std::size_t molecule, std::size_t localAtom, bool rigid,
                                        const double3& moleculePosition, const double3& moleculeCom,
                                        const double3& frameworkPosition, double f1, double f2, const double3& dr)
{
  const double3 offset = rigid ? moleculePosition - moleculeCom : double3{};
  const double3 exactBase = dr - offset;
  const double3 rawBase = moleculeCom - (flexibleFramework ? frameworkPosition : double3{});
  const double3 baseCorrection = exactBase - rawBase;
  const double3 gradient = f1 * dr;
  double3x3 strain{};
  strain.ax = exactBase.x * gradient.x;
  strain.bx = exactBase.y * gradient.x;
  strain.cx = exactBase.z * gradient.x;
  strain.ay = exactBase.x * gradient.y;
  strain.by = exactBase.y * gradient.y;
  strain.cy = exactBase.z * gradient.y;
  strain.az = exactBase.x * gradient.z;
  strain.bz = exactBase.y * gradient.z;
  strain.cz = exactBase.z * gradient.z;
  hessian.strainGradient() += strain;

  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const std::size_t cellA = *layout.cellDof(a);
    const double3 deltaV = cellLayout.bases[a] * baseCorrection;
    const double3 mixedCorrection = f1 * deltaV + f2 * double3::dot(dr, deltaV) * dr;
    // atomRigidComDof covers both whole rigid molecules and rigid groups of semi-flexible molecules.
    const std::size_t moleculeBase =
        rigid ? *layout.atomRigidComDof(molecule, localAtom) : layout.flexibleAtomDofBase(molecule, localAtom);
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      hessian.add(moleculeBase + axis, cellA, (&mixedCorrection.x)[axis]);
      hessian.add(cellA, moleculeBase + axis, (&mixedCorrection.x)[axis]);
      if (flexibleFramework)
      {
        const std::size_t frameworkBase = *layout.frameworkAtomDof(frameworkAtom, MinimizationDofAxis::X);
        hessian.add(frameworkBase + axis, cellA, -(&mixedCorrection.x)[axis]);
        hessian.add(cellA, frameworkBase + axis, -(&mixedCorrection.x)[axis]);
      }
    }
    if (rigid)
    {
      const std::size_t orientationBase = *layout.atomRigidOrientationDof(molecule, localAtom);
      const Minimization::RigidAtomDerivatives& derivatives = rigidCache.atom(molecule, localAtom);
      const std::array<double3, 3> directions = {derivatives.dVecX, derivatives.dVecY, derivatives.dVecZ};
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double value = double3::dot(directions[axis], mixedCorrection);
        hessian.add(orientationBase + axis, cellA, value);
        hessian.add(cellA, orientationBase + axis, value);
      }
    }

    const double3 exactV = cellLayout.bases[a] * exactBase;
    const double3 rawV = cellLayout.bases[a] * rawBase;
    for (std::size_t b = 0; b < cellLayout.size(); ++b)
    {
      const double3 exactW = cellLayout.bases[b] * exactBase;
      const double3 rawW = cellLayout.bases[b] * rawBase;
      const double exactValue = f1 * double3::dot(exactV, exactW) +
                                f2 * double3::dot(dr, exactV) * double3::dot(dr, exactW) +
                                double3::dot(gradient, cellStrainSecondDerivative(cellLayout, a, b) * exactBase);
      const double rawValue = f1 * double3::dot(rawV, rawW) + f2 * double3::dot(dr, rawV) * double3::dot(dr, rawW) +
                              double3::dot(gradient, cellStrainSecondDerivative(cellLayout, a, b) * rawBase);
      hessian.add(cellA, *layout.cellDof(b), exactValue - rawValue);
    }
  }
}
}  // namespace

RunningEnergy Interactions::computeFrameworkMoleculeHessian(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Molecule> moleculeData,
    std::span<const Component> components, std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> moleculeDynamics,
    std::span<AtomDynamics> frameworkDynamics, const CellMinimizationLayout& cellLayout, const Framework* framework)
{
  RunningEnergy energies{};

  if (frameworkAtoms.empty() || moleculeAtoms.empty())
  {
    return energies;
  }

  const bool mixed = framework != nullptr && framework->isMixed();
  const Minimization::RigidDerivativeCache rigidCache = Minimization::RigidDerivativeCache::build(
      moleculeData, components, moleculeAtoms, mixed ? framework : nullptr, frameworkAtoms);

  const SimulationBox& box = simulationBox;
  const bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const bool flexibleFramework = !mixed && layout.numberOfFrameworkAtoms() == frameworkAtoms.size();

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    const std::size_t moleculeIndex = static_cast<std::size_t>(it1->moleculeId);
    const Molecule& molecule = moleculeData[moleculeIndex];
    const std::size_t localAtom = static_cast<std::size_t>(it1 - moleculeAtoms.begin()) - molecule.atomIndex;
    const bool rigid = layout.atomRigidComDof(moleculeIndex, localAtom).has_value();
    // Framework atoms scale affinely with the cell under strain (consistent with the Ewald
    // Fourier treatment), so the strain derivative of dr only loses the rigid internal offset
    // pos - com of the molecule (or rigid group) side.
    const double3 comA = rigid ? rigidCache.bodyCenterOfMass(moleculeIndex, localAtom) : it1->position;

    for (std::span<const Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
    {
      const std::size_t frameworkAtom = static_cast<std::size_t>(it2 - frameworkAtoms.begin());

      // Shared gradient/Hessian scatter for both the VDW and Coulomb contributions of this pair.
      const auto applyPairHessian = [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
      {
        if (moleculeDynamics.size() == moleculeAtoms.size())
        {
          moleculeDynamics[static_cast<std::size_t>(it1 - moleculeAtoms.begin())].gradient +=
              factors.firstDerivativeFactor * dr;
        }
        if (frameworkDynamics.size() == frameworkAtoms.size())
        {
          frameworkDynamics[frameworkAtom].gradient -= factors.firstDerivativeFactor * dr;
        }
        if (mixed)
        {
          // Mixed frameworks: project both sides onto their generalized DOFs (fixed framework
          // atoms carry an empty site and drop out; rigid-group atoms couple through the group's
          // center-of-mass and orientation DOFs). Variable-cell terms are not supported here.
          const double3 gradientA = factors.firstDerivativeFactor * dr;
          Minimization::scatterRadialHessianSites(
              hessian,
              {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, localAtom),
               Minimization::makeFrameworkHessianSite(layout, rigidCache, frameworkAtom)},
              {gradientA, -gradientA}, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
          return;
        }
        if (flexibleFramework)
        {
          Minimization::scatterFlexibleFrameworkMoleculeHessian(
              hessian, layout, rigidCache, frameworkAtom, moleculeIndex, localAtom, rigid,
              factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        else
        {
          Minimization::scatterFrameworkMoleculeHessian(hessian, layout, rigidCache, moleculeIndex, localAtom, rigid,
                                                        factors.firstDerivativeFactor, factors.secondDerivativeFactor,
                                                        dr);
        }
        addFrameworkMoleculeCellCorrection(hessian, layout, rigidCache, cellLayout, frameworkAtom, flexibleFramework,
                                           moleculeIndex, localAtom, rigid, it1->position, comA, it2->position,
                                           factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        if (hessian.numStrain() == 1)
        {
          const double3 drStrainDerivative = dr - (it1->position - comA);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeIndex, localAtom, rigid,
                                                           1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterAtomicStrainStrainIsotropic(hessian, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, it1->position, comA,
                                                           it2->position, it2->position, rigid, false);
        }
      };

      Interactions::evaluatePair<2>(
          forceField, box, *it1, *it2, cutOffFrameworkVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
          {
            energies.frameworkMoleculeVDW += factors.energy;
            applyPairHessian(factors, dr);
          },
          [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
          {
            energies.frameworkMoleculeCharge += factors.energy;
            applyPairHessian(factors, dr);
          });
    }
  }

  return energies;
}
