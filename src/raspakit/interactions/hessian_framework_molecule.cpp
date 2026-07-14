module;

module interactions_hessian_framework_molecule;

import std;

import double3;
import double3x3;
import atom;
import molecule;
import forcefield;
import simulationbox;
import potential_hessian_vdw;
import potential_hessian_coulomb;
import hessian_factor;
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
    const std::size_t moleculeBase =
        rigid ? *layout.rigidMoleculeDof(molecule, RigidDof::ComX) : layout.flexibleAtomDofBase(molecule, localAtom);
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
      const std::size_t orientationBase = *layout.rigidMoleculeDof(molecule, RigidDof::OriX);
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

RunningEnergy Interactions::computeFrameworkMoleculeHessian(const System& system, const MinimizationDofLayout& layout,
                                                            GeneralizedHessian& hessian,
                                                            std::span<AtomDynamics> moleculeDynamics,
                                                            std::span<AtomDynamics> frameworkDynamics)
{
  RunningEnergy energies{};

  if (!system.framework.has_value())
  {
    return energies;
  }

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  if (frameworkAtoms.empty() || moleculeAtoms.empty())
  {
    return energies;
  }

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, moleculeAtoms);

  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  const bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const bool flexibleFramework = layout.numberOfFrameworkAtoms() == frameworkAtoms.size();

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    const std::size_t moleculeIndex = static_cast<std::size_t>(it1->moleculeId);
    const bool rigid = layout.molecules()[moleculeIndex].rigid;
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const std::size_t localAtom = static_cast<std::size_t>(it1 - moleculeAtoms.begin()) - molecule.atomIndex;
    // Framework atoms scale affinely with the cell under strain (consistent with the Ewald
    // Fourier treatment), so the strain derivative of dr only loses the rigid internal offset
    // pos - com of the molecule side.
    const double3 comA = rigid ? molecule.centerOfMassPosition : it1->position;
    const std::size_t typeB = static_cast<std::size_t>(it1->type);
    const double scalingVDWB = it1->scalingVDW;
    const double scalingCoulombB = it1->scalingCoulomb;
    const double chargeB = it1->charge;

    for (std::span<const Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
    {
      const std::size_t frameworkAtom = static_cast<std::size_t>(it2 - frameworkAtoms.begin());
      double3 dr = it1->position - it2->position;
      dr = box.applyPeriodicBoundaryConditions(dr);
      const double rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        const std::size_t typeA = static_cast<std::size_t>(it2->type);
        const double scalingVDWA = it2->scalingVDW;
        const Potentials::HessianFactor factors =
            Potentials::potentialVDWHessian(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        energies.frameworkMoleculeVDW += factors.energy;
        if (moleculeDynamics.size() == moleculeAtoms.size())
        {
          moleculeDynamics[static_cast<std::size_t>(it1 - moleculeAtoms.begin())].gradient +=
              factors.firstDerivativeFactor * dr;
        }
        if (frameworkDynamics.size() == frameworkAtoms.size())
        {
          frameworkDynamics[frameworkAtom].gradient -= factors.firstDerivativeFactor * dr;
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
      }

      if (useCharge && rr < cutOffChargeSquared)
      {
        const double r = std::sqrt(rr);
        const double scalingCoulombA = it2->scalingCoulomb;
        const double chargeA = it2->charge;
        const Potentials::HessianFactor factors =
            Potentials::potentialCoulombHessian(forceField, scalingCoulombA, scalingCoulombB, rr, r, chargeA, chargeB);
        energies.frameworkMoleculeCharge += factors.energy;
        if (moleculeDynamics.size() == moleculeAtoms.size())
        {
          moleculeDynamics[static_cast<std::size_t>(it1 - moleculeAtoms.begin())].gradient +=
              factors.firstDerivativeFactor * dr;
        }
        if (frameworkDynamics.size() == frameworkAtoms.size())
        {
          frameworkDynamics[frameworkAtom].gradient -= factors.firstDerivativeFactor * dr;
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
      }
    }
  }

  return energies;
}
