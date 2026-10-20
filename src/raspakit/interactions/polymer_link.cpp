module;

module interactions_polymer_link;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import simulationbox;
import forcefield;
import component;
import polymer_links;
import units;
import bond_potential;
import bend_potential;
import torsion_potential;
import potential_pair_derivatives;
import potential_intra_pair;
import potential_coulomb_real_space;
import interactions_pair_kernel;

namespace
{
using PairState = PolymerLinkTable::PairState;
using PairScalingInfo = PolymerLinkTable::PairScalingInfo;

/// The energy (by slot) and the first-derivative factor of a pair in one of its three treatments; the gradient on
/// the first atom is factor * dr.
struct PairValue
{
  double vdw{0.0};
  double coulomb{0.0};
  double exclusion{0.0};
  double factor{0.0};
};

template <std::size_t Order>
PairValue pairValue(const ForceField &forceField, const Atom &a, const Atom &b, const double3 &dr, double rr,
                    PairState state, const PairScalingInfo &scaling)
{
  PairValue value{};
  const double r = std::sqrt(rr);
  switch (state)
  {
    case PairState::Full:
    {
      const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
      const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
      Interactions::evaluatePairAtSeparation<Order>(
          forceField, a, b, dr, rr, cutOffVDWSquared, cutOffChargeSquared, forceField.useCharge,
          [&](const Potentials::PairDerivatives<Order> &factors, const double3 &)
          {
            value.vdw += factors.energy;
            if constexpr (Order >= 1) value.factor += factors.firstDerivativeFactor;
          },
          [&](const Potentials::PairDerivatives<Order> &factors, const double3 &)
          {
            value.coulomb += factors.energy;
            if constexpr (Order >= 1) value.factor += factors.firstDerivativeFactor;
          });
      break;
    }
    case PairState::Scaled:
    {
      const Potentials::PairDerivatives<Order> vdw = Potentials::intraMolecularVDW<Order>(
          forceField, scaling.scalingVDW, rr, static_cast<std::size_t>(a.type), static_cast<std::size_t>(b.type),
          scaling.pair14);
      value.vdw += vdw.energy;
      if constexpr (Order >= 1) value.factor += vdw.firstDerivativeFactor;
      if (forceField.useCharge)
      {
        const Potentials::PairDerivatives<Order> coulomb = Potentials::intraMolecularCoulomb<Order>(
            forceField, scaling.scalingCoulomb, a.scalingCoulomb, b.scalingCoulomb, r, a.charge, b.charge);
        value.coulomb += coulomb.energy;
        if constexpr (Order >= 1) value.factor += coulomb.firstDerivativeFactor;
      }
      break;
    }
    case PairState::Excluded:
    {
      if (!forceField.useCharge || r <= 0.0) break;
      const double prefactor = Units::CoulombicConversionFactor * a.charge * b.charge;
      if (prefactor == 0.0) break;
      const double scalingTotal = a.scalingCoulomb * b.scalingCoulomb;
      if (forceField.usesEwaldFourier())
      {
        // the reciprocal sum counted the pair: remove it the way the intramolecular exclusion does (no cutoff)
        const Potentials::EwaldExclusionFactors exclusionFactors =
            Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, scalingTotal, r);
        value.exclusion -= scalingTotal * prefactor * exclusionFactors.potential;
        if constexpr (Order >= 1) value.factor -= scalingTotal * prefactor * exclusionFactors.firstDerivativeFactor;
      }
      else if (forceField.usesRealSpaceChargeCorrections() && rr < forceField.cutOffCoulomb * forceField.cutOffCoulomb)
      {
        // completion of the shifted pair sum of the finite-cutoff charge methods: q_i q_j (V(r) - 1/r)
        const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
        value.exclusion += scalingTotal * prefactor * (factors.potential - 1.0 / r);
        if constexpr (Order >= 1)
        {
          value.factor += scalingTotal * prefactor * (factors.firstDerivativeFactor + 1.0 / (rr * r));
        }
      }
      break;
    }
  }
  return value;
}

inline void addOuterProduct(double3x3 &tensor, const double3 &arm, const double3 &gradient)
{
  tensor.ax += arm.x * gradient.x;
  tensor.ay += arm.x * gradient.y;
  tensor.az += arm.x * gradient.z;
  tensor.bx += arm.y * gradient.x;
  tensor.by += arm.y * gradient.y;
  tensor.bz += arm.y * gradient.z;
  tensor.cx += arm.z * gradient.x;
  tensor.cy += arm.z * gradient.y;
  tensor.cz += arm.z * gradient.z;
}

/**
 * The evaluation of the overlay at derivative order 0 (energies) or 1 (gradients).
 *
 * 'linkGradients' receives the gradients of the terms that act between molecules (the polymer bonds and their
 * bends / torsions, the pair corrections); 'overrideGradients' those of the intramolecular overrides. 'strain'
 * and 'correction' (optional) receive the strain derivative and the atomic-to-molecular virial correction of the
 * inter-molecular terms only.
 */
template <std::size_t Order, typename LinkSink, typename OverrideSink>
RunningEnergy evaluate(const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> atoms,
                       const PolymerLinkTable &table, LinkSink &&linkGradients, OverrideSink &&overrideGradients,
                       double3x3 *strain, double3x3 *correction)
{
  RunningEnergy energy{};
  if (!table.derivedValid)
  {
    throw std::runtime_error("[Polymerization]: the derived topology of the polymer bonds is out of date\n");
  }
  const PolymerLinkTable::Layout &layout = table.layout;

  // the mass-weighted center of mass of the molecule of an atom (stored, unwrapped positions), on demand
  std::unordered_map<std::size_t, double3> centerOfMass{};
  auto comOf = [&](std::uint32_t atom) -> double3
  {
    const PolymerAtomRef r = layout.ref(atom);
    const std::size_t molecule = layout.moleculeOffset[r.componentId] + r.moleculeIndex;
    auto it = centerOfMass.find(molecule);
    if (it != centerOfMass.end()) return it->second;
    const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
    const std::size_t count = layout.atomsPerMolecule[r.componentId];
    double totalMass = 0.0;
    double3 com(0.0, 0.0, 0.0);
    for (std::size_t i = start; i < start + count; ++i)
    {
      const double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atoms[i].type)].mass;
      com += mass * atoms[i].position;
      totalMass += mass;
    }
    com = com / totalMass;
    centerOfMass.emplace(molecule, com);
    return com;
  };
  auto applyLink = [&](std::uint32_t atom, const double3 &gradient)
  {
    if constexpr (Order >= 1)
    {
      linkGradients(atom, gradient);
      if (correction) addOuterProduct(*correction, atoms[atom].position - comOf(atom), gradient);
    }
  };

  // positions of a chain of bonded atoms as one consistent geometry (the molecules are stored whole, so a bond
  // between two molecules needs the periodic image of the partner)
  auto chain = [&](std::span<const std::uint32_t> ids)
  {
    std::array<double3, 4> positions{};
    positions[0] = atoms[ids[0]].position;
    for (std::size_t i = 1; i < ids.size(); ++i)
    {
      positions[i] = positions[i - 1] +
                     simulationBox.applyPeriodicBoundaryConditions(atoms[ids[i]].position - positions[i - 1]);
    }
    return positions;
  };
  auto lambdaOfBond = [&](std::uint32_t bondId) { return table.rampFactor(table.bonds[bondId].formedAtStep); };

  // --- the bonded terms across the polymer bonds ------------------------------------------------------------------
  for (const PolymerLinkTable::LinkBondTerm &term : table.linkBonds)
  {
    const double lambda = lambdaOfBond(term.bondId);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.crossLink += lambda * term.potential.calculateEnergy(p[0], p[1]);
    }
    else
    {
      const auto [e, g, s] = term.potential.potentialEnergyGradientStrain(p[0], p[1]);
      energy.crossLink += lambda * e;
      applyLink(term.atoms[0], lambda * g[0]);
      applyLink(term.atoms[1], lambda * g[1]);
      if (strain) *strain += lambda * s;
    }
  }
  for (const PolymerLinkTable::LinkBendTerm &term : table.linkBends)
  {
    const double lambda = lambdaOfBond(term.bondId);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.crossLink += lambda * term.potential.calculateEnergy(p[0], p[1], p[2], std::nullopt);
    }
    else
    {
      const auto [e, g, s] = term.potential.potentialEnergyGradientStrain(p[0], p[1], p[2]);
      energy.crossLink += lambda * e;
      for (std::size_t i = 0; i < 3; ++i) applyLink(term.atoms[i], lambda * g[i]);
      if (strain) *strain += lambda * s;
    }
  }
  for (const PolymerLinkTable::LinkTorsionTerm &term : table.linkTorsions)
  {
    const double lambda = lambdaOfBond(term.bondId);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.crossLink += lambda * term.potential.calculateEnergy(p[0], p[1], p[2], p[3]);
    }
    else
    {
      const auto [e, g, s] = term.potential.potentialEnergyGradientStrain(p[0], p[1], p[2], p[3]);
      energy.crossLink += lambda * e;
      for (std::size_t i = 0; i < 4; ++i) applyLink(term.atoms[i], lambda * g[i]);
      if (strain) *strain += lambda * s;
    }
  }

  // --- the restraints of the stereocentres: 1/2 k (V - chi V_0)^2 on the normalised chiral volume --------------
  if (table.chiralForceConstant > 0.0)
  {
    for (const PolymerLinkTable::StereoCentreTerm &term : table.stereoCentres)
    {
      const double lambda = lambdaOfBond(term.bondId);
      if (lambda == 0.0) continue;
      const double3 centre = atoms[term.atoms[0]].position;
      const double3 bp = simulationBox.applyPeriodicBoundaryConditions(atoms[term.atoms[1]].position - centre);
      const double3 bn = simulationBox.applyPeriodicBoundaryConditions(atoms[term.atoms[2]].position - centre);
      const double3 bz = simulationBox.applyPeriodicBoundaryConditions(atoms[term.atoms[3]].position - centre);
      const ChiralVolume volume = ChiralVolume::compute(bp, bn, bz);
      const double deviation = volume.value - static_cast<double>(term.chirality) * table.chiralReferenceVolume;
      energy.crossLink += lambda * 0.5 * table.chiralForceConstant * deviation * deviation;
      if constexpr (Order >= 1)
      {
        const double factor = lambda * table.chiralForceConstant * deviation;
        const std::array<double3, 3> arms{bp, bn, bz};
        for (std::size_t i = 0; i < 4; ++i)
        {
          const double3 gradient = factor * volume.gradient[i];
          applyLink(term.atoms[i], gradient);
          if (strain && i > 0) addOuterProduct(*strain, arms[i - 1], gradient);
        }
      }
    }
  }

  // --- the pair corrections ---------------------------------------------------------------------------------------
  for (const PolymerLinkTable::PairCorrection &pair : table.pairCorrections)
  {
    const double lambda = lambdaOfBond(pair.bondId);
    if (lambda == 0.0) continue;
    const Atom &a = atoms[pair.atoms[0]];
    const Atom &b = atoms[pair.atoms[1]];
    const double3 dr = simulationBox.applyPeriodicBoundaryConditions(a.position - b.position);
    const double rr = double3::dot(dr, dr);
    const PairValue wanted = pairValue<Order>(forceField, a, b, dr, rr, pair.wanted, pair.wantedScaling);
    const PairValue counted = pairValue<Order>(forceField, a, b, dr, rr, pair.counted, pair.countedScaling);
    energy.moleculeMoleculeVDW += lambda * (wanted.vdw - counted.vdw);
    energy.moleculeMoleculeCharge += lambda * (wanted.coulomb - counted.coulomb);
    energy.ewald_exclusion += lambda * (wanted.exclusion - counted.exclusion);
    if constexpr (Order >= 1)
    {
      const double factor = lambda * (wanted.factor - counted.factor);
      if (factor == 0.0) continue;
      const double3 gradient = factor * dr;
      applyLink(pair.atoms[0], gradient);
      applyLink(pair.atoms[1], -gradient);
      if (strain) Interactions::accumulateStrainDerivative(*strain, gradient, dr);
    }
  }

  // --- the overrides of the retyped monomer terms (intramolecular: no strain, no virial correction) --------------
  for (const PolymerLinkTable::OverrideBondTerm &term : table.overrideBonds)
  {
    const double lambda = table.rampFactor(term.step);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.bond +=
          lambda * (term.replacement.calculateEnergy(p[0], p[1]) - term.original.calculateEnergy(p[0], p[1]));
    }
    else
    {
      const auto [eNew, gNew, sNew] = term.replacement.potentialEnergyGradientStrain(p[0], p[1]);
      const auto [eOld, gOld, sOld] = term.original.potentialEnergyGradientStrain(p[0], p[1]);
      energy.bond += lambda * (eNew - eOld);
      for (std::size_t i = 0; i < 2; ++i) overrideGradients(term.atoms[i], lambda * (gNew[i] - gOld[i]));
    }
  }
  for (const PolymerLinkTable::OverrideBendTerm &term : table.overrideBends)
  {
    const double lambda = table.rampFactor(term.step);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.bend += lambda * (term.replacement.calculateEnergy(p[0], p[1], p[2], std::nullopt) -
                               term.original.calculateEnergy(p[0], p[1], p[2], std::nullopt));
    }
    else
    {
      const auto [eNew, gNew, sNew] = term.replacement.potentialEnergyGradientStrain(p[0], p[1], p[2]);
      const auto [eOld, gOld, sOld] = term.original.potentialEnergyGradientStrain(p[0], p[1], p[2]);
      energy.bend += lambda * (eNew - eOld);
      for (std::size_t i = 0; i < 3; ++i) overrideGradients(term.atoms[i], lambda * (gNew[i] - gOld[i]));
    }
  }
  for (const PolymerLinkTable::OverrideTorsionTerm &term : table.overrideTorsions)
  {
    const double lambda = table.rampFactor(term.step);
    if (lambda == 0.0) continue;
    const std::array<double3, 4> p = chain(term.atoms);
    if constexpr (Order == 0)
    {
      energy.torsion += lambda * (term.replacement.calculateEnergy(p[0], p[1], p[2], p[3]) -
                                  term.original.calculateEnergy(p[0], p[1], p[2], p[3]));
    }
    else
    {
      const auto [eNew, gNew, sNew] = term.replacement.potentialEnergyGradientStrain(p[0], p[1], p[2], p[3]);
      const auto [eOld, gOld, sOld] = term.original.potentialEnergyGradientStrain(p[0], p[1], p[2], p[3]);
      energy.torsion += lambda * (eNew - eOld);
      for (std::size_t i = 0; i < 4; ++i) overrideGradients(term.atoms[i], lambda * (gNew[i] - gOld[i]));
    }
  }

  return energy;
}

void ensureDerived(const ForceField &forceField, const std::vector<Component> &components,
                   const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
                   PolymerLinkTable &table)
{
  if (!table.derivedValid) table.rebuild(forceField, components, numberOfMoleculesPerComponent, moleculeAtoms);
}
}  // namespace

RunningEnergy Interactions::computePolymerLinkEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                     const std::vector<Component> &components,
                                                     const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                                     std::span<const Atom> moleculeAtoms, PolymerLinkTable &table)
{
  if (table.empty()) return {};
  ensureDerived(forceField, components, numberOfMoleculesPerComponent, moleculeAtoms, table);
  auto noGradient = [](std::uint32_t, const double3 &) {};
  return evaluate<0>(forceField, simulationBox, moleculeAtoms, table, noGradient, noGradient, nullptr, nullptr);
}

RunningEnergy Interactions::computePolymerLinkGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                                       const std::vector<Component> &components,
                                                       const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                                       std::span<const Atom> moleculeAtoms,
                                                       std::span<AtomDynamics> moleculeDynamics, PolymerLinkTable &table,
                                                       double3x3 *strain, double3x3 *molecularCorrection)
{
  if (table.empty()) return {};
  ensureDerived(forceField, components, numberOfMoleculesPerComponent, moleculeAtoms, table);
  auto addGradient = [&](std::uint32_t atom, const double3 &g) { moleculeDynamics[atom].gradient += g; };
  return evaluate<1>(forceField, simulationBox, moleculeAtoms, table, addGradient, addGradient, strain,
                     molecularCorrection);
}

std::pair<RunningEnergy, double3x3> Interactions::computePolymerLinkEnergyStrainDerivative(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    std::span<AtomDynamics> pressureDynamics, PolymerLinkTable &table)
{
  if (table.empty()) return {};
  ensureDerived(forceField, components, numberOfMoleculesPerComponent, moleculeAtoms, table);
  double3x3 strain{};
  auto addGradient = [&](std::uint32_t atom, const double3 &g) { pressureDynamics[atom].gradient += g; };
  auto noGradient = [](std::uint32_t, const double3 &) {};
  const RunningEnergy energy =
      evaluate<1>(forceField, simulationBox, moleculeAtoms, table, addGradient, noGradient, &strain, nullptr);
  return {energy, strain};
}
