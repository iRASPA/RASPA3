module;

export module interactions_pair_kernel;

import std;

import double3;
import double3x3;
import atom;
import simulationbox;
import forcefield;
export import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;

export namespace Interactions
{
/**
 * \brief Evaluates the VDW and Coulomb interaction of one atom pair at the requested derivative order.
 *
 * Shared pair kernel for the energy (Order 0), gradient (Order 1), and Hessian (Order 2) routines:
 * computes the minimum-image separation, applies the cutoffs, dispatches to the unified pair
 * potentials, and hands the resulting factors to the caller-supplied sinks. Each sink is invoked
 * with (const Potentials::PairDerivatives<Order>&, const double3& dr), where dr = posA - posB
 * after periodic boundary conditions; the Cartesian force on atom A is firstDerivativeFactor * dr.
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
[[clang::always_inline]] inline void evaluatePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                                  const Atom& atomA, const Atom& atomB, double cutOffVDWSquared,
                                                  double cutOffChargeSquared, bool useCharge, VDWSink&& vdwSink,
                                                  CoulombSink&& coulombSink)
{
  double3 dr = atomA.position - atomB.position;
  dr = simulationBox.applyPeriodicBoundaryConditions(dr);
  const double rr = double3::dot(dr, dr);

  if (rr < cutOffVDWSquared)
  {
    const Potentials::PairDerivatives<Order> factors =
        Potentials::potentialVDW<Order>(forceField, atomA.scalingVDW, atomB.scalingVDW, rr,
                                        static_cast<std::size_t>(atomA.type), static_cast<std::size_t>(atomB.type));
    vdwSink(factors, dr);
  }
  if (useCharge && rr < cutOffChargeSquared)
  {
    const double r = std::sqrt(rr);
    const Potentials::PairDerivatives<Order> factors = Potentials::potentialCoulomb<Order>(
        forceField, atomA.scalingCoulomb, atomB.scalingCoulomb, r, atomA.charge, atomB.charge);
    coulombSink(factors, dr);
  }
}

/**
 * \brief Loops over all intermolecular atom pairs and evaluates them at the requested derivative order.
 *
 * Triangular double loop with the molecule self-interaction skip (atoms sharing a molecule id do
 * not interact). Uses the molecule-molecule VDW cutoff. Sinks are invoked with
 * (std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
 *  const Potentials::PairDerivatives<Order>&, const double3& dr).
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
inline void forEachMoleculeMoleculePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                        std::span<const Atom> moleculeAtoms, VDWSink&& vdwSink,
                                        CoulombSink&& coulombSink)
{
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return;

  for (std::size_t indexA = 0; indexA + 1 < moleculeAtoms.size(); ++indexA)
  {
    const Atom& atomA = moleculeAtoms[indexA];
    for (std::size_t indexB = indexA + 1; indexB < moleculeAtoms.size(); ++indexB)
    {
      const Atom& atomB = moleculeAtoms[indexB];

      // skip interactions within the same molecule
      if (atomA.moleculeId == atomB.moleculeId) continue;

      evaluatePair<Order>(
          forceField, simulationBox, atomA, atomB, cutOffVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
          { vdwSink(indexA, indexB, atomA, atomB, factors, dr); },
          [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
          { coulombSink(indexA, indexB, atomA, atomB, factors, dr); });
    }
  }
}

/**
 * \brief Loops over all framework atoms interacting with one molecule atom at the requested derivative order.
 *
 * Uses the framework-molecule VDW cutoff. dr = moleculeAtom.position - frameworkAtom.position, so
 * firstDerivativeFactor * dr is the force on the molecule atom. Sinks are invoked with
 * (std::size_t frameworkIndex, const Atom& frameworkAtom,
 *  const Potentials::PairDerivatives<Order>&, const double3& dr).
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
inline void forEachFrameworkMoleculePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                         const Atom& moleculeAtom, std::span<const Atom> frameworkAtoms,
                                         VDWSink&& vdwSink, CoulombSink&& coulombSink)
{
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::size_t index = 0; index < frameworkAtoms.size(); ++index)
  {
    const Atom& frameworkAtom = frameworkAtoms[index];
    evaluatePair<Order>(
        forceField, simulationBox, moleculeAtom, frameworkAtom, cutOffVDWSquared, cutOffChargeSquared, useCharge,
        [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
        { vdwSink(index, frameworkAtom, factors, dr); },
        [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
        { coulombSink(index, frameworkAtom, factors, dr); });
  }
}

/**
 * \brief Accumulates the strain-derivative contribution gradient (x) displacement into the tensor.
 *
 * Column layout matches the existing strain conventions: tensor.{a,b,c}{x,y,z} +=
 * gradient.{x,y,z} * displacement.{x,y,z}.
 */
inline void accumulateStrainDerivative(double3x3& tensor, const double3& gradient, const double3& displacement)
{
  tensor.ax += gradient.x * displacement.x;
  tensor.bx += gradient.y * displacement.x;
  tensor.cx += gradient.z * displacement.x;

  tensor.ay += gradient.x * displacement.y;
  tensor.by += gradient.y * displacement.y;
  tensor.cy += gradient.z * displacement.y;

  tensor.az += gradient.x * displacement.z;
  tensor.bz += gradient.y * displacement.z;
  tensor.cz += gradient.z * displacement.z;
}
}  // namespace Interactions
