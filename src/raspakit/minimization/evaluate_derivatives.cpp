module;

module minimization_evaluate_derivatives;

import std;

import double3;
import running_energy;

void evaluateDerivatives(System& system, const MinimizationDofLayout& layout, DerivativeCapabilities capabilities,
                         DerivativeResults& results)
{
  results.energy = 0.0;

  const std::size_t nStrain = (capabilities.hessianPositionStrain || capabilities.hessianStrainStrain) ? 1 : 0;
  if (results.hessian.numDofs() != layout.numDofs() || results.hessian.numStrain() != nStrain)
  {
    results.hessian.resize(layout.numDofs(), nStrain);
  }

  if (capabilities.hessianPositionPosition || capabilities.hessianPositionStrain || capabilities.hessianStrainStrain)
  {
    results.hessian.zero();
  }

  if (!results.gradient.empty())
  {
    std::ranges::fill(results.gradient, 0.0);
  }

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<AtomDynamics> frameworkDynamics = system.spanOfFrameworkDynamics();

  for (AtomDynamics& dynamics : system.atomDynamics)
  {
    dynamics.gradient = double3(0.0, 0.0, 0.0);
  }

  if (capabilities.energy || capabilities.gradient || capabilities.hessianPositionPosition ||
      capabilities.hessianPositionStrain || capabilities.hessianStrainStrain)
  {
    RunningEnergy intraEnergy = Interactions::computeIntraMolecularBondHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += intraEnergy.bond;

    RunningEnergy bendEnergy = Interactions::computeIntraMolecularBendHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bendEnergy.bend;

    RunningEnergy ureyBradleyEnergy = Interactions::computeIntraMolecularUreyBradleyHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += ureyBradleyEnergy.ureyBradley;

    RunningEnergy torsionEnergy = Interactions::computeIntraMolecularTorsionHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += torsionEnergy.torsion + torsionEnergy.improperTorsion;

    RunningEnergy intraVDWEnergy = Interactions::computeIntraMolecularVanDerWaalsHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += intraVDWEnergy.intraVDW;

    RunningEnergy intraCoulombEnergy = Interactions::computeIntraMolecularCoulombHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += intraCoulombEnergy.intraCoul;

    RunningEnergy bondBondEnergy = Interactions::computeIntraMolecularBondBondHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bondBondEnergy.bondBond;

    RunningEnergy bondBendEnergy = Interactions::computeIntraMolecularBondBendHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bondBendEnergy.bondBend;

    RunningEnergy bendBendEnergy = Interactions::computeIntraMolecularBendBendHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bendBendEnergy.bendBend;

    RunningEnergy bondTorsionEnergy = Interactions::computeIntraMolecularBondTorsionHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bondTorsionEnergy.bondTorsion;

    RunningEnergy bendTorsionEnergy = Interactions::computeIntraMolecularBendTorsionHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += bendTorsionEnergy.bendTorsion;

    RunningEnergy inversionBendEnergy = Interactions::computeIntraMolecularInversionBendHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += inversionBendEnergy.inversionBend;

    if (system.framework && !system.framework->rigid)
    {
      const RunningEnergy frameworkIntraEnergy = Interactions::computeFrameworkIntraMolecularHessian(
          *system.framework, system.simulationBox, frameworkAtomPositions, layout, results.hessian, frameworkDynamics);
      results.energy += frameworkIntraEnergy.potentialEnergy();
    }

    if (capabilities.hessianPositionPosition || capabilities.hessianPositionStrain || capabilities.hessianStrainStrain)
    {
      RunningEnergy interEnergy =
          Interactions::computeInterMolecularHessian(system, layout, results.hessian, moleculeDynamics);
      results.energy += interEnergy.moleculeMoleculeVDW + interEnergy.moleculeMoleculeCharge;

      RunningEnergy frameworkEnergy = Interactions::computeFrameworkMoleculeHessian(
          system, layout, results.hessian, moleculeDynamics, frameworkDynamics);
      results.energy += frameworkEnergy.frameworkMoleculeVDW + frameworkEnergy.frameworkMoleculeCharge;

      RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierHessian(system, layout, results.hessian,
                                                                           moleculeDynamics, frameworkDynamics);
      results.energy += ewaldEnergy.ewald_fourier + ewaldEnergy.ewald_self + ewaldEnergy.ewald_exclusion;
    }

    if (capabilities.gradient && !results.gradient.empty())
    {
      for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
      {
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          if (const auto dof = layout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis)))
          {
            results.gradient[*dof] = (&frameworkDynamics[atom].gradient.x)[axis];
          }
        }
      }

      const Minimization::RigidDerivativeCache rigidCache =
          Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, moleculeAtomPositions);
      std::size_t moleculeIndex = 0;
      for (const Molecule& molecule : system.moleculeData)
      {
        const Component& component = system.components[molecule.componentId];
        if (!component.rigid)
        {
          for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
          {
            const AtomDynamics& dynamics = moleculeDynamics[molecule.atomIndex + localAtom];
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
              if (auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis)))
              {
                results.gradient[*dof] = (&dynamics.gradient.x)[axis];
              }
            }
          }
        }
        else
        {
          double3 centerOfMassGradient{};
          double3 orientationGradient{};
          for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
          {
            const AtomDynamics& dynamics = moleculeDynamics[molecule.atomIndex + localAtom];
            const Minimization::RigidAtomDerivatives& derivatives = rigidCache.atom(moleculeIndex, localAtom);
            centerOfMassGradient += dynamics.gradient;
            orientationGradient.x += double3::dot(dynamics.gradient, derivatives.dVecX);
            orientationGradient.y += double3::dot(dynamics.gradient, derivatives.dVecY);
            orientationGradient.z += double3::dot(dynamics.gradient, derivatives.dVecZ);
          }
          if (const auto base = layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX))
          {
            results.gradient[*base + 0] = centerOfMassGradient.x;
            results.gradient[*base + 1] = centerOfMassGradient.y;
            results.gradient[*base + 2] = centerOfMassGradient.z;
          }
          if (const auto base = layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX))
          {
            results.gradient[*base + 0] = orientationGradient.x;
            results.gradient[*base + 1] = orientationGradient.y;
            results.gradient[*base + 2] = orientationGradient.z;
          }
        }
        ++moleculeIndex;
      }
    }
  }
}
