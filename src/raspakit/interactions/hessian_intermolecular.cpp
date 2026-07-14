module;

module interactions_hessian_intermolecular;

import std;

import double3;
import double3x3;
import atom;
import molecule;
import forcefield;
import simulationbox;
import interactions_pair_kernel;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;
import minimization_cell_layout;

namespace
{
void addStrainGradient(GeneralizedHessian& hessian, const double3& strainDisplacement, const double3& gradient)
{
  double3x3 strain{};
  strain.ax = strainDisplacement.x * gradient.x;
  strain.bx = strainDisplacement.y * gradient.x;
  strain.cx = strainDisplacement.z * gradient.x;
  strain.ay = strainDisplacement.x * gradient.y;
  strain.by = strainDisplacement.y * gradient.y;
  strain.cy = strainDisplacement.z * gradient.y;
  strain.az = strainDisplacement.x * gradient.z;
  strain.bz = strainDisplacement.y * gradient.z;
  strain.cz = strainDisplacement.z * gradient.z;
  hessian.strainGradient() += strain;
}

void addPeriodicCellCorrection(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                               const Minimization::RigidDerivativeCache& rigidCache,
                               const CellMinimizationLayout& cellLayout, std::size_t moleculeA, std::size_t localAtomA,
                               bool rigidA, const double3& positionA, const double3& comA, std::size_t moleculeB,
                               std::size_t localAtomB, bool rigidB, const double3& positionB, const double3& comB,
                               double f1, double f2, const double3& dr)
{
  const double3 offsetA = rigidA ? positionA - comA : double3{};
  const double3 offsetB = rigidB ? positionB - comB : double3{};
  const double3 exactBase = dr + offsetB - offsetA;
  const double3 rawBase = comA - comB;
  const double3 baseCorrection = exactBase - rawBase;
  const double3 gradient = f1 * dr;
  addStrainGradient(hessian, exactBase, gradient);

  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const double3 deltaV = cellLayout.bases[a] * baseCorrection;
    const double3 mixedCorrection = f1 * deltaV + f2 * double3::dot(dr, deltaV) * dr;
    const std::size_t cellA = *layout.cellDof(a);
    const auto scatterSite = [&](std::size_t molecule, std::size_t localAtom, bool rigid, double sign)
    {
      const std::size_t positionBase =
          rigid ? *layout.rigidMoleculeDof(molecule, RigidDof::ComX) : layout.flexibleAtomDofBase(molecule, localAtom);
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        hessian.add(positionBase + axis, cellA, sign * (&mixedCorrection.x)[axis]);
        hessian.add(cellA, positionBase + axis, sign * (&mixedCorrection.x)[axis]);
      }
      if (rigid)
      {
        const std::size_t orientationBase = *layout.rigidMoleculeDof(molecule, RigidDof::OriX);
        const Minimization::RigidAtomDerivatives& derivatives = rigidCache.atom(molecule, localAtom);
        const std::array<double3, 3> directions = {derivatives.dVecX, derivatives.dVecY, derivatives.dVecZ};
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          const double value = sign * double3::dot(directions[axis], mixedCorrection);
          hessian.add(orientationBase + axis, cellA, value);
          hessian.add(cellA, orientationBase + axis, value);
        }
      }
    };
    scatterSite(moleculeA, localAtomA, rigidA, 1.0);
    scatterSite(moleculeB, localAtomB, rigidB, -1.0);

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

namespace
{
std::size_t localAtomIndex(std::span<const Atom> atoms, const Molecule& molecule, const Atom& atom)
{
  return static_cast<std::size_t>(&atom - &atoms[molecule.atomIndex]);
}
}  // namespace

RunningEnergy Interactions::computeInterMolecularHessian(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Molecule> moleculeData,
    std::span<const Component> components, std::span<const Atom> moleculeAtoms, const MinimizationDofLayout& layout,
    GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics, const CellMinimizationLayout& cellLayout)
{
  RunningEnergy energies{};

  if (forceField.omitInterInteractions)
  {
    return energies;
  }

  std::span<const Atom> atoms = moleculeAtoms;
  if (atoms.size() < 2)
  {
    return energies;
  }

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  const SimulationBox& box = simulationBox;
  const bool useCharge = forceField.useCharge;
  const double cutOffMoleculeVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = atoms.begin(); it1 != atoms.end() - 1; ++it1)
  {
    const std::size_t moleculeA = static_cast<std::size_t>(it1->moleculeId);
    const bool rigidA = layout.molecules()[moleculeA].rigid;
    const double3 comA = rigidA ? moleculeData[moleculeA].centerOfMassPosition : it1->position;

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != atoms.end(); ++it2)
    {
      const std::size_t moleculeB = static_cast<std::size_t>(it2->moleculeId);
      if (moleculeA == moleculeB)
      {
        continue;
      }

      const bool rigidB = layout.molecules()[moleculeB].rigid;
      const double3 comB = rigidB ? moleculeData[moleculeB].centerOfMassPosition : it2->position;

      const Molecule& molA = moleculeData[moleculeA];
      const Molecule& molB = moleculeData[moleculeB];
      const std::size_t localAtomA = localAtomIndex(atoms, molA, *it1);
      const std::size_t localAtomB = localAtomIndex(atoms, molB, *it2);

      // Shared gradient/Hessian scatter for both the VDW and Coulomb contributions of this pair.
      const auto applyPairHessian = [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
      {
        if (dynamics.size() == atoms.size())
        {
          const double3 gradient = factors.firstDerivativeFactor * dr;
          dynamics[static_cast<std::size_t>(it1 - atoms.begin())].gradient += gradient;
          dynamics[static_cast<std::size_t>(it2 - atoms.begin())].gradient -= gradient;
        }
        if (rigidA || rigidB)
        {
          Minimization::scatterInteractionHessian(
              hessian, layout, rigidCache, moleculeA, localAtomA, rigidA, it1->position, comA, moleculeB, localAtomB,
              rigidB, it2->position, comB, factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        }
        else
        {
          Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeA, localAtomA, moleculeB, localAtomB,
                                                      factors.firstDerivativeFactor, factors.secondDerivativeFactor,
                                                      dr);
        }
        addPeriodicCellCorrection(hessian, layout, rigidCache, cellLayout, moleculeA, localAtomA, rigidA, it1->position,
                                  comA, moleculeB, localAtomB, rigidB, it2->position, comB,
                                  factors.firstDerivativeFactor, factors.secondDerivativeFactor, dr);
        if (hessian.numStrain() == 1)
        {
          const double3 drStrainDerivative = dr + (it2->position - comB) - (it1->position - comA);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeA, localAtomA, rigidA,
                                                           1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterSitePositionStrainIsotropic(hessian, layout, rigidCache, moleculeB, localAtomB, rigidB,
                                                           -1.0, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, drStrainDerivative);
          Minimization::scatterAtomicStrainStrainIsotropic(hessian, factors.firstDerivativeFactor,
                                                           factors.secondDerivativeFactor, dr, it1->position, comA,
                                                           it2->position, comB, rigidA, rigidB);
        }
      };

      Interactions::evaluatePair<2>(
          forceField, box, *it1, *it2, cutOffMoleculeVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
          {
            energies.moleculeMoleculeVDW += factors.energy;
            applyPairHessian(factors, dr);
          },
          [&](const Potentials::PairDerivatives<2>& factors, const double3& dr)
          {
            energies.moleculeMoleculeCharge += factors.energy;
            applyPairHessian(factors, dr);
          });
    }
  }

  return energies;
}
