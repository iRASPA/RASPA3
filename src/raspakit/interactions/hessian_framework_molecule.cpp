module;

module interactions_hessian_framework_molecule;

import std;

import double3;
import atom;
import molecule;
import forcefield;
import simulationbox;
import potential_hessian_vdw;
import potential_hessian_coulomb;
import hessian_factor;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;

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
        if (layout.numberOfFrameworkAtoms() == frameworkAtoms.size())
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
        if (layout.numberOfFrameworkAtoms() == frameworkAtoms.size())
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
