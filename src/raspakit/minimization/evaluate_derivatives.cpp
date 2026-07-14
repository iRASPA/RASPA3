module;

module minimization_evaluate_derivatives;

import std;

import double3;
import double3x3;
import running_energy;
import minimization_cell_layout;
import minimization_rigid_kinematics;
import interactions_intermolecular;
import interactions_framework_molecule;

namespace
{
void assemblePositionGradient(const System& system, const MinimizationDofLayout& layout,
                              std::span<const Atom> moleculeAtomPositions,
                              std::span<const AtomDynamics> moleculeDynamics,
                              std::span<const AtomDynamics> frameworkDynamics, std::span<double> gradient)
{
  for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
  {
    const std::size_t base = *layout.frameworkAtomDof(atom, MinimizationDofAxis::X);
    gradient[base + 0] = frameworkDynamics[atom].gradient.x;
    gradient[base + 1] = frameworkDynamics[atom].gradient.y;
    gradient[base + 2] = frameworkDynamics[atom].gradient.z;
  }

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(system.moleculeData, system.components, moleculeAtomPositions);
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const Component& component = system.components[molecule.componentId];
    if (!component.rigid)
    {
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        const AtomDynamics& dynamics = moleculeDynamics[molecule.atomIndex + localAtom];
        const std::size_t base = layout.flexibleAtomDofBase(moleculeIndex, localAtom);
        gradient[base + 0] = dynamics.gradient.x;
        gradient[base + 1] = dynamics.gradient.y;
        gradient[base + 2] = dynamics.gradient.z;
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
      const std::size_t centerBase = *layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
      gradient[centerBase + 0] = centerOfMassGradient.x;
      gradient[centerBase + 1] = centerOfMassGradient.y;
      gradient[centerBase + 2] = centerOfMassGradient.z;
      const std::size_t orientationBase = *layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
      gradient[orientationBase + 0] = orientationGradient.x;
      gradient[orientationBase + 1] = orientationGradient.y;
      gradient[orientationBase + 2] = orientationGradient.z;
    }
  }
}

std::vector<double> affinePositionDirection(const System& system, const MinimizationDofLayout& layout,
                                            const double3x3& basis)
{
  std::vector<double> direction(layout.numberOfPositionDofs(), 0.0);
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
  {
    const double3 value = basis * frameworkAtoms[atom].position;
    const std::size_t base = *layout.frameworkAtomDof(atom, MinimizationDofAxis::X);
    direction[base + 0] = value.x;
    direction[base + 1] = value.y;
    direction[base + 2] = value.z;
  }

  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    if (layout.molecules()[moleculeIndex].rigid)
    {
      const double3 value = basis * molecule.centerOfMassPosition;
      const std::size_t base = *layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
      direction[base + 0] = value.x;
      direction[base + 1] = value.y;
      direction[base + 2] = value.z;
    }
    else
    {
      for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
      {
        const double3 value = basis * atoms[molecule.atomIndex + atom].position;
        const std::size_t base = layout.flexibleAtomDofBase(moleculeIndex, atom);
        direction[base + 0] = value.x;
        direction[base + 1] = value.y;
        direction[base + 2] = value.z;
      }
    }
  }
  return direction;
}

std::vector<double> affineGradientDirection(const System& system, const MinimizationDofLayout& layout,
                                            const double3x3& basis, std::span<const double> gradient)
{
  std::vector<double> direction(layout.numberOfPositionDofs(), 0.0);
  for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
  {
    const std::size_t base = *layout.frameworkAtomDof(atom, MinimizationDofAxis::X);
    const double3 value = basis * double3(gradient[base + 0], gradient[base + 1], gradient[base + 2]);
    direction[base + 0] = value.x;
    direction[base + 1] = value.y;
    direction[base + 2] = value.z;
  }
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    if (layout.molecules()[moleculeIndex].rigid)
    {
      const std::size_t base = *layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
      const double3 value = basis * double3(gradient[base + 0], gradient[base + 1], gradient[base + 2]);
      direction[base + 0] = value.x;
      direction[base + 1] = value.y;
      direction[base + 2] = value.z;
    }
    else
    {
      for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
      {
        const std::size_t base = layout.flexibleAtomDofBase(moleculeIndex, atom);
        const double3 value = basis * double3(gradient[base + 0], gradient[base + 1], gradient[base + 2]);
        direction[base + 0] = value.x;
        direction[base + 1] = value.y;
        direction[base + 2] = value.z;
      }
    }
  }
  return direction;
}

double matrixContraction(const double3x3& lhs, const double3x3& rhs)
{
  double result = 0.0;
  for (std::size_t column = 0; column < 3; ++column)
  {
    for (std::size_t row = 0; row < 3; ++row)
    {
      result += lhs.mm[column][row] * rhs.mm[column][row];
    }
  }
  return result;
}

void addAffineCellDerivatives(const System& system, const MinimizationDofLayout& layout,
                              const CellMinimizationLayout& cellLayout, std::span<double> gradient,
                              GeneralizedHessian& hessian)
{
  std::vector<std::vector<double>> firstDirections;
  std::vector<std::vector<double>> secondDirections(cellLayout.size() * cellLayout.size());
  firstDirections.reserve(cellLayout.size());
  for (const double3x3& basis : cellLayout.bases)
  {
    firstDirections.push_back(affinePositionDirection(system, layout, basis));
  }
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    for (std::size_t b = 0; b < cellLayout.size(); ++b)
    {
      secondDirections[a * cellLayout.size() + b] =
          affinePositionDirection(system, layout, cellStrainSecondDerivative(cellLayout, a, b));
    }
  }

  const std::size_t positionDofs = layout.numberOfPositionDofs();
  const std::size_t totalDofs = layout.numDofs();
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const std::size_t cellA = *layout.cellDof(a);
    const std::span<const double> positionGradient = gradient.first(positionDofs);
    const std::vector<double> geometricMixed =
        affineGradientDirection(system, layout, cellLayout.bases[a], positionGradient);
    for (std::size_t position = 0; position < positionDofs; ++position)
    {
      double mixed = 0.0;
      for (std::size_t q = 0; q < positionDofs; ++q)
      {
        mixed += hessian.positionPosition()[position * totalDofs + q] * firstDirections[a][q];
      }
      mixed += geometricMixed[position];
      hessian.add(position, cellA, mixed);
      hessian.add(cellA, position, mixed);
    }
    for (std::size_t b = 0; b < cellLayout.size(); ++b)
    {
      double value = std::inner_product(positionGradient.begin(), positionGradient.end(),
                                        secondDirections[a * cellLayout.size() + b].begin(), 0.0);
      for (std::size_t p = 0; p < positionDofs; ++p)
      {
        for (std::size_t q = 0; q < positionDofs; ++q)
        {
          value += firstDirections[a][p] * hessian.positionPosition()[p * totalDofs + q] * firstDirections[b][q];
        }
      }
      hessian.add(cellA, *layout.cellDof(b), value);
    }
  }
}
}  // namespace

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
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);

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

    RunningEnergy outOfPlaneBendEnergy = Interactions::computeIntraMolecularOutOfPlaneBendHessian(
        system.moleculeData, moleculeAtomPositions, system.components, layout, results.hessian, moleculeDynamics);
    results.energy += outOfPlaneBendEnergy.outOfPlaneBend;

    if (system.framework && !system.framework->rigid)
    {
      const RunningEnergy frameworkIntraEnergy = Interactions::computeFrameworkIntraMolecularHessian(
          system.forceField, *system.framework, system.simulationBox, frameworkAtomPositions, layout, results.hessian,
          frameworkDynamics, cellLayout);
      results.energy += frameworkIntraEnergy.potentialEnergy();
    }

    if (capabilities.hessianPositionPosition || capabilities.hessianPositionStrain || capabilities.hessianStrainStrain)
    {
      RunningEnergy interEnergy = Interactions::computeInterMolecularHessian(
          system.forceField, system.simulationBox, system.moleculeData, system.components, moleculeAtomPositions,
          layout, results.hessian, moleculeDynamics, cellLayout);
      results.energy += interEnergy.moleculeMoleculeVDW + interEnergy.moleculeMoleculeCharge;

      RunningEnergy frameworkEnergy = Interactions::computeFrameworkMoleculeHessian(
          system.forceField, system.simulationBox, system.moleculeData, system.components, frameworkAtomPositions,
          moleculeAtomPositions, layout, results.hessian, moleculeDynamics, frameworkDynamics, cellLayout);
      results.energy += frameworkEnergy.frameworkMoleculeVDW + frameworkEnergy.frameworkMoleculeCharge;

      if (layout.numberOfCellDofs() != cellLayout.size())
      {
        throw std::invalid_argument("evaluateDerivatives: cell layout does not match minimization DOFs");
      }
      if (layout.numberOfCellDofs() != 0)
      {
        if (!capabilities.gradient || results.gradient.size() != layout.numDofs())
        {
          throw std::invalid_argument("evaluateDerivatives: variable-cell derivatives require the full gradient");
        }
        assemblePositionGradient(system, layout, moleculeAtomPositions, moleculeDynamics, frameworkDynamics,
                                 results.gradient);
        addAffineCellDerivatives(system, layout, cellLayout, results.gradient, results.hessian);
        for (std::size_t a = 0; a < cellLayout.size(); ++a)
        {
          results.gradient[*layout.cellDof(a)] =
              matrixContraction(results.hessian.strainGradient(), cellLayout.bases[a]);
        }
      }

      const RunningEnergy intermolecularTail =
          Interactions::computeInterMolecularTailEnergy(system.forceField, system.simulationBox, moleculeAtomPositions);
      const RunningEnergy frameworkTail = Interactions::computeFrameworkMoleculeTailEnergy(
          system.forceField, system.simulationBox, frameworkAtomPositions, moleculeAtomPositions);
      const double tailEnergy = intermolecularTail.tail + frameworkTail.tail;
      results.energy += tailEnergy;
      for (std::size_t a = 0; a < cellLayout.size(); ++a)
      {
        const double traceA = cellLayout.bases[a].trace();
        results.gradient[*layout.cellDof(a)] -= tailEnergy * traceA;
        for (std::size_t b = 0; b < cellLayout.size(); ++b)
        {
          results.hessian.add(*layout.cellDof(a), *layout.cellDof(b),
                              tailEnergy * traceA * cellLayout.bases[b].trace());
        }
      }

      const double3x3 strainBeforeEwald = results.hessian.strainGradient();
      RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierHessian(
          system.forceField, system.simulationBox, system.framework, system.fixedFrameworkStoredEik,
          system.netChargeFramework, system.moleculeData, system.components, frameworkAtomPositions,
          moleculeAtomPositions, layout, results.hessian, moleculeDynamics, frameworkDynamics, cellLayout);
      results.energy += ewaldEnergy.ewald_fourier + ewaldEnergy.ewald_self + ewaldEnergy.ewald_exclusion;
      if (layout.numberOfCellDofs() != 0)
      {
        const double3x3 ewaldStrainGradient = results.hessian.strainGradient() - strainBeforeEwald;
        for (std::size_t a = 0; a < cellLayout.size(); ++a)
        {
          results.gradient[*layout.cellDof(a)] += matrixContraction(ewaldStrainGradient, cellLayout.bases[a]);
        }

        const double pressureVolume = system.pressure * system.simulationBox.volume;
        results.energy += pressureVolume;
        for (std::size_t a = 0; a < cellLayout.size(); ++a)
        {
          const double traceA = cellLayout.bases[a].trace();
          results.gradient[*layout.cellDof(a)] += pressureVolume * traceA;
          for (std::size_t b = 0; b < cellLayout.size(); ++b)
          {
            results.hessian.add(*layout.cellDof(a), *layout.cellDof(b),
                                pressureVolume * traceA * cellLayout.bases[b].trace());
          }
        }
      }
    }

    if (capabilities.gradient && !results.gradient.empty())
    {
      assemblePositionGradient(system, layout, moleculeAtomPositions, moleculeDynamics, frameworkDynamics,
                               results.gradient);
    }
  }
}
