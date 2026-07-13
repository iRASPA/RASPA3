module;

module interactions_hessian_intermolecular;

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

namespace
{
std::size_t localAtomIndex(std::span<const Atom> atoms, const Molecule &molecule, const Atom &atom)
{
  return static_cast<std::size_t>(&atom - &atoms[molecule.atomIndex]);
}
}  // namespace

RunningEnergy Interactions::computeInterMolecularHessian(const System &system, const MinimizationDofLayout &layout,
                                                         GeneralizedHessian &hessian,
                                                         std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  if (system.forceField.omitInterInteractions)
  {
    return energies;
  }

  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  if (atoms.size() < 2)
  {
    return energies;
  }

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, atoms);

  const ForceField &forceField = system.forceField;
  const SimulationBox &box = system.simulationBox;
  const bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = atoms.begin(); it1 != atoms.end() - 1; ++it1)
  {
    const std::size_t moleculeA = static_cast<std::size_t>(it1->moleculeId);
    const bool rigidA = layout.molecules()[moleculeA].rigid;
    const double3 comA =
        rigidA ? system.moleculeData[moleculeA].centerOfMassPosition : it1->position;
    const std::size_t typeA = static_cast<std::size_t>(it1->type);
    const double scalingVDWA = it1->scalingVDW;
    const double scalingCoulombA = it1->scalingCoulomb;
    const double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != atoms.end(); ++it2)
    {
      const std::size_t moleculeB = static_cast<std::size_t>(it2->moleculeId);
      if (moleculeA == moleculeB)
      {
        continue;
      }

      const bool rigidB = layout.molecules()[moleculeB].rigid;
      const double3 comB =
          rigidB ? system.moleculeData[moleculeB].centerOfMassPosition : it2->position;

      double3 dr = it1->position - it2->position;
      dr = box.applyPeriodicBoundaryConditions(dr);
      const double rr = double3::dot(dr, dr);

      const Molecule &molA = system.moleculeData[moleculeA];
      const Molecule &molB = system.moleculeData[moleculeB];
      const std::size_t localAtomA = localAtomIndex(atoms, molA, *it1);
      const std::size_t localAtomB = localAtomIndex(atoms, molB, *it2);

      if (rr < cutOffMoleculeVDWSquared)
      {
        const std::size_t typeB = static_cast<std::size_t>(it2->type);
        const double scalingVDWB = it2->scalingVDW;
        const Potentials::HessianFactor factors =
            Potentials::potentialVDWHessian(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        energies.moleculeMoleculeVDW += factors.energy;
        if (dynamics.size() == atoms.size())
        {
          const double3 gradient = factors.firstDerivativeFactor * dr;
          dynamics[static_cast<std::size_t>(it1 - atoms.begin())].gradient += gradient;
          dynamics[static_cast<std::size_t>(it2 - atoms.begin())].gradient -= gradient;
        }
        if (rigidA || rigidB)
        {
          Minimization::scatterInteractionHessian(hessian, layout, rigidCache, moleculeA, localAtomA, rigidA,
                                                  it1->position, comA, moleculeB, localAtomB, rigidB, it2->position,
                                                  comB, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        else
        {
          Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeA, localAtomA, moleculeB, localAtomB,
                                                      factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        if (hessian.numStrain() == 1)
        {
          const double3 drStrainDerivative = dr + (it2->position - comB) - (it1->position - comA);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeA, localAtomA, rigidA,
                                                           1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeB, localAtomB, rigidB,
                                                           -1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterAtomicStrainStrainIsotropic(
              hessian, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr, it1->position, comA,
              it2->position, comB, rigidA, rigidB);
        }
      }

      if (useCharge && rr < cutOffChargeSquared)
      {
        const double r = std::sqrt(rr);
        const double scalingCoulombB = it2->scalingCoulomb;
        const double chargeB = it2->charge;
        const Potentials::HessianFactor factors = Potentials::potentialCoulombHessian(
            forceField, scalingCoulombA, scalingCoulombB, rr, r, chargeA, chargeB);
        energies.moleculeMoleculeCharge += factors.energy;
        if (dynamics.size() == atoms.size())
        {
          const double3 gradient = factors.firstDerivativeFactor * dr;
          dynamics[static_cast<std::size_t>(it1 - atoms.begin())].gradient += gradient;
          dynamics[static_cast<std::size_t>(it2 - atoms.begin())].gradient -= gradient;
        }
        if (rigidA || rigidB)
        {
          Minimization::scatterInteractionHessian(hessian, layout, rigidCache, moleculeA, localAtomA, rigidA,
                                                  it1->position, comA, moleculeB, localAtomB, rigidB, it2->position,
                                                  comB, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        else
        {
          Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeA, localAtomA, moleculeB, localAtomB,
                                                      factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        if (hessian.numStrain() == 1)
        {
          const double3 drStrainDerivative = dr + (it2->position - comB) - (it1->position - comA);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeA, localAtomA, rigidA,
                                                           1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeB, localAtomB, rigidB,
                                                           -1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterAtomicStrainStrainIsotropic(
              hessian, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr, it1->position, comA,
              it2->position, comB, rigidA, rigidB);
        }
      }
    }
  }

  return energies;
}
