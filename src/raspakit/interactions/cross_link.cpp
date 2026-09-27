module;

module interactions_cross_link;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import simulationbox;
import forcefield;
import component;
import cross_links;
import units;
import bond_potential;
import bend_potential;
import connectivity_table;
import interactions_pair_kernel;
import potential_coulomb_real_space;

namespace
{

// Positions of the two molecules of a link as one consistent geometry: molecule B is shifted by the
// periodic image that brings its site next to the site of molecule A. Molecules are stored whole
// (unwrapped), so a single shift per molecule is exact.
struct LinkGeometry
{
  double3 posA;
  double3 posB;   // image of B's site next to A's site
  double3 shiftB; // posB - stored position of B's site
};

// Removes the non-bonded pair interaction of an excluded pair (x, y) that the inter-molecular sum has
// counted, and adds what an intramolecular pair receives instead (the Ewald exclusion or the
// shifted-potential completion). 'addGradient(onX, g)' receives the gradient on x (and -g on y).
template <std::size_t Order, typename GradientSink>
void excludePair(const ForceField &forceField, const SimulationBox &simulationBox, const Atom &x, const Atom &y,
                 RunningEnergy &energy, GradientSink &&addGradient, double3x3 *strain)
{
  const double3 dr = simulationBox.applyPeriodicBoundaryConditions(x.position - y.position);
  const double rr = double3::dot(dr, dr);
  const double r = std::sqrt(rr);
  const double scaling = x.scalingCoulomb * y.scalingCoulomb;

  auto accumulate = [&](double3 gradientOnX)
  {
    if constexpr (Order >= 1)
    {
      addGradient(gradientOnX);
      if (strain) Interactions::accumulateStrainDerivative(*strain, gradientOnX, dr);
    }
  };

  if (!forceField.omitInterInteractions)
  {
    const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    Interactions::evaluatePair<Order>(
        forceField, simulationBox, x, y, cutOffVDWSquared, cutOffChargeSquared, forceField.useCharge,
        [&](const Potentials::PairDerivatives<Order> &factors, const double3 &d)
        {
          energy.moleculeMoleculeVDW -= factors.energy;
          if constexpr (Order >= 1) accumulate(-factors.firstDerivativeFactor * d);
        },
        [&](const Potentials::PairDerivatives<Order> &factors, const double3 &d)
        {
          energy.moleculeMoleculeCharge -= factors.energy;
          if constexpr (Order >= 1) accumulate(-factors.firstDerivativeFactor * d);
        });
  }

  if (!forceField.useCharge || r <= 0.0) return;
  const double prefactor = Units::CoulombicConversionFactor * x.charge * y.charge;

  if (forceField.usesEwaldFourier())
  {
    // The reciprocal sum counted this pair: remove it the way the intramolecular exclusion does.
    const Potentials::EwaldExclusionFactors exclusion =
        Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, scaling, r);
    energy.ewald_exclusion -= scaling * prefactor * exclusion.potential;
    if constexpr (Order >= 1) accumulate(-scaling * prefactor * exclusion.firstDerivativeFactor * dr);
  }
  else if (forceField.usesRealSpaceChargeCorrections() && !forceField.omitInterInteractions &&
           rr < forceField.cutOffCoulomb * forceField.cutOffCoulomb)
  {
    // Completion of the shifted pair sum for the finite-cutoff charge methods: q_i q_j (V(r) - 1/r).
    const Potentials::CoulombRealSpaceFactors factors = Potentials::coulombRealSpaceFactors(forceField, r);
    energy.ewald_exclusion += scaling * prefactor * (factors.potential - 1.0 / r);
    if constexpr (Order >= 1)
    {
      accumulate(scaling * prefactor * (factors.firstDerivativeFactor + 1.0 / (rr * r)) * dr);
    }
  }
}

std::vector<std::size_t> neighboursOf(const Component &component, std::size_t atom)
{
  if (atom >= component.connectivityTable.numberOfBeads) return {};
  return component.connectivityTable.findAllNeighbors(atom);
}

// The terms of one link between site 'siteA' of molecule A and site 'siteB' of molecule B.
// 'atomA(local)' / 'atomB(local)' return the atoms of the two molecules by local index;
// 'addGradientA(local, g)' / 'addGradientB(local, g)' receive the gradients (Order >= 1 only).
// 'neighboursA' / 'neighboursB' are the intramolecular neighbours of the two sites.
template <std::size_t Order, typename AtomA, typename AtomB, typename GradientA, typename GradientB>
void evaluateLinkTerms(const ForceField &forceField, const SimulationBox &simulationBox,
                       const CrossLinkBondType &bondType, std::size_t siteA, std::size_t siteB,
                       std::span<const std::size_t> neighboursA, std::span<const std::size_t> neighboursB,
                       AtomA &&atomA, AtomB &&atomB, RunningEnergy &energy, GradientA &&addGradientA,
                       GradientB &&addGradientB, double3x3 *strain)
{
  const Atom &atomSiteA = atomA(siteA);
  const Atom &atomSiteB = atomB(siteB);

  LinkGeometry geometry;
  geometry.posA = atomSiteA.position;
  const double3 drAB = simulationBox.applyPeriodicBoundaryConditions(atomSiteA.position - atomSiteB.position);
  geometry.posB = atomSiteA.position - drAB;
  geometry.shiftB = geometry.posB - atomSiteB.position;

  // Bond (plus the constant formation energy).
  if constexpr (Order == 0)
  {
    energy.crossLink += bondType.bond.calculateEnergy(geometry.posA, geometry.posB) + bondType.formationEnergy;
  }
  else
  {
    const auto [bondEnergy, gradients, bondStrain] =
        bondType.bond.potentialEnergyGradientStrain(geometry.posA, geometry.posB);
    energy.crossLink += bondEnergy + bondType.formationEnergy;
    addGradientA(siteA, gradients[0]);
    addGradientB(siteB, gradients[1]);
    if (strain) *strain += bondStrain;
  }

  // Junction bends (n_A, A, B) and (A, B, n_B).
  if (bondType.junctionBend.has_value())
  {
    const BendPotential &bend = bondType.junctionBend.value();
    for (std::size_t nA : neighboursA)
    {
      const double3 posN = atomA(nA).position;
      if constexpr (Order == 0)
      {
        energy.crossLink += bend.calculateEnergy(posN, geometry.posA, geometry.posB, std::nullopt);
      }
      else
      {
        const auto [bendEnergy, gradients, bendStrain] =
            bend.potentialEnergyGradientStrain(posN, geometry.posA, geometry.posB);
        energy.crossLink += bendEnergy;
        addGradientA(nA, gradients[0]);
        addGradientA(siteA, gradients[1]);
        addGradientB(siteB, gradients[2]);
        if (strain) *strain += bendStrain;
      }
    }
    for (std::size_t nB : neighboursB)
    {
      const double3 posN = atomB(nB).position + geometry.shiftB;
      if constexpr (Order == 0)
      {
        energy.crossLink += bend.calculateEnergy(geometry.posA, geometry.posB, posN, std::nullopt);
      }
      else
      {
        const auto [bendEnergy, gradients, bendStrain] =
            bend.potentialEnergyGradientStrain(geometry.posA, geometry.posB, posN);
        energy.crossLink += bendEnergy;
        addGradientA(siteA, gradients[0]);
        addGradientB(siteB, gradients[1]);
        addGradientB(nB, gradients[2]);
        if (strain) *strain += bendStrain;
      }
    }
  }

  // Non-bonded exclusions across the link: the 1-2 pair and the 1-3 pairs.
  auto exclude = [&](std::size_t localA, std::size_t localB)
  {
    excludePair<Order>(
        forceField, simulationBox, atomA(localA), atomB(localB), energy,
        [&](const double3 &g)
        {
          addGradientA(localA, g);
          addGradientB(localB, -g);
        },
        strain);
  };
  exclude(siteA, siteB);
  for (std::size_t nA : neighboursA) exclude(nA, siteB);
  for (std::size_t nB : neighboursB) exclude(siteA, nB);
}

// Evaluates one link of the table. 'atomAt(site, localAtom)' returns the atom 'localAtom' of the
// molecule that 'site' belongs to; 'addGradient(site, localAtom, g)' receives the gradients (Order >= 1
// only).
template <std::size_t Order, typename AtomAccessor, typename GradientSink>
void evaluateLink(const ForceField &forceField, const SimulationBox &simulationBox,
                  const std::vector<Component> &components, const CrossLinkTable &table, const CrossLink &link,
                  AtomAccessor &&atomAt, RunningEnergy &energy, GradientSink &&addGradient, double3x3 *strain)
{
  const std::vector<std::size_t> neighboursA = neighboursOf(components[link.a.componentId], link.a.atomIndex);
  const std::vector<std::size_t> neighboursB = neighboursOf(components[link.b.componentId], link.b.atomIndex);

  evaluateLinkTerms<Order>(
      forceField, simulationBox, table.bondTypes[link.bondTypeId], link.a.atomIndex, link.b.atomIndex, neighboursA,
      neighboursB, [&](std::size_t local) -> const Atom & { return atomAt(link.a, local); },
      [&](std::size_t local) -> const Atom & { return atomAt(link.b, local); }, energy,
      [&](std::size_t local, const double3 &g) { addGradient(link.a, local, g); },
      [&](std::size_t local, const double3 &g) { addGradient(link.b, local, g); }, strain);
}

}  // namespace

Interactions::CrossLinkAtomLayout Interactions::CrossLinkAtomLayout::make(
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent)
{
  CrossLinkAtomLayout layout;
  layout.componentOffset.resize(components.size());
  layout.atomsPerMolecule.resize(components.size());
  std::size_t offset = 0;
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    layout.componentOffset[c] = offset;
    layout.atomsPerMolecule[c] = components[c].atoms.size();
    offset += numberOfMoleculesPerComponent[c] * components[c].atoms.size();
  }
  return layout;
}

RunningEnergy Interactions::computeCrossLinkEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                   const std::vector<Component> &components,
                                                   const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                                   std::span<const Atom> moleculeAtoms, const CrossLinkTable &table)
{
  return computeCrossLinkEnergyOfLinks(forceField, simulationBox, components, numberOfMoleculesPerComponent,
                                       moleculeAtoms, table, table.links);
}

RunningEnergy Interactions::computeCrossLinkEnergyOfLinks(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    const CrossLinkTable &table, std::span<const CrossLink> links)
{
  RunningEnergy energy{};
  if (links.empty()) return energy;

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site, std::size_t localAtom) -> const Atom &
  { return moleculeAtoms[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom]; };
  auto noGradient = [](const CrossLinkSite &, std::size_t, const double3 &) {};

  for (const CrossLink &link : links)
  {
    evaluateLink<0>(forceField, simulationBox, components, table, link, atomAt, energy, noGradient, nullptr);
  }
  return energy;
}

RunningEnergy Interactions::computeCrossLinkEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    const CrossLinkTable &table, std::size_t componentId, std::size_t moleculeIndex, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms)
{
  if (table.empty()) return {};
  const std::span<const std::size_t> linkIds = table.linksOfMolecule(componentId, moleculeIndex);
  if (linkIds.empty()) return {};

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto isMoved = [&](const CrossLinkSite &site)
  { return site.componentId == componentId && site.moleculeIndex == moleculeIndex; };
  auto atomNew = [&](const CrossLinkSite &site, std::size_t localAtom) -> const Atom &
  {
    if (isMoved(site)) return newAtoms[localAtom];
    return moleculeAtoms[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom];
  };
  auto atomOld = [&](const CrossLinkSite &site, std::size_t localAtom) -> const Atom &
  {
    if (isMoved(site)) return oldAtoms[localAtom];
    return moleculeAtoms[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom];
  };
  auto noGradient = [](const CrossLinkSite &, std::size_t, const double3 &) {};

  RunningEnergy energyNew{};
  RunningEnergy energyOld{};
  for (std::size_t id : linkIds)
  {
    const CrossLink &link = table.links[id];
    evaluateLink<0>(forceField, simulationBox, components, table, link, atomNew, energyNew, noGradient, nullptr);
    evaluateLink<0>(forceField, simulationBox, components, table, link, atomOld, energyOld, noGradient, nullptr);
  }
  return energyNew - energyOld;
}

std::vector<CrossLinkTether> Interactions::makeCrossLinkTethers(
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> moleculeAtoms, const CrossLinkTable &table, std::size_t componentId,
    std::size_t moleculeIndex)
{
  std::vector<CrossLinkTether> tethers{};
  if (table.empty()) return tethers;
  const std::span<const std::size_t> linkIds = table.linksOfMolecule(componentId, moleculeIndex);
  if (linkIds.empty()) return tethers;

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  tethers.reserve(linkIds.size());
  for (std::size_t id : linkIds)
  {
    const CrossLink &link = table.links[id];
    const bool aIsOwn = link.a.componentId == componentId && link.a.moleculeIndex == moleculeIndex;
    const CrossLinkSite &own = aIsOwn ? link.a : link.b;
    const CrossLinkSite &partner = aIsOwn ? link.b : link.a;

    CrossLinkTether tether{};
    tether.siteAtom = own.atomIndex;
    tether.bondType = &table.bondTypes[link.bondTypeId];
    const std::size_t partnerStart = layout.moleculeStart(partner.componentId, partner.moleculeIndex);
    tether.partnerSiteAtom = partner.atomIndex;
    tether.partnerSite = moleculeAtoms[partnerStart + partner.atomIndex];
    for (std::size_t n : neighboursOf(components[partner.componentId], partner.atomIndex))
    {
      tether.partnerNeighbours.emplace_back(n, moleculeAtoms[partnerStart + n]);
    }
    tethers.push_back(std::move(tether));
  }
  return tethers;
}

RunningEnergy Interactions::computeCrossLinkTetherEnergy(const ForceField &forceField,
                                                         const SimulationBox &simulationBox,
                                                         const Component &component,
                                                         std::span<const Atom> moleculeAtoms,
                                                         std::span<const CrossLinkTether> tethers,
                                                         std::span<const std::size_t> selected)
{
  RunningEnergy energy{};
  auto noGradient = [](std::size_t, const double3 &) {};

  for (std::size_t index : selected)
  {
    const CrossLinkTether &tether = tethers[index];
    const std::vector<std::size_t> neighboursA = neighboursOf(component, tether.siteAtom);
    std::vector<std::size_t> neighboursB{};
    neighboursB.reserve(tether.partnerNeighbours.size());
    for (const auto &[local, atom] : tether.partnerNeighbours) neighboursB.push_back(local);

    // The partner molecule is known only through its frozen site and neighbours.
    auto atomB = [&](std::size_t local) -> const Atom &
    {
      if (local == tether.partnerSiteAtom) return tether.partnerSite;
      const auto it = std::ranges::find_if(tether.partnerNeighbours,
                                           [local](const auto &entry) { return entry.first == local; });
      return it->second;
    };

    evaluateLinkTerms<0>(
        forceField, simulationBox, *tether.bondType, tether.siteAtom, tether.partnerSiteAtom, neighboursA,
        neighboursB, [&](std::size_t local) -> const Atom & { return moleculeAtoms[local]; }, atomB, energy,
        noGradient, noGradient, nullptr);
  }
  return energy;
}

RunningEnergy Interactions::computeCrossLinkGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                                     const std::vector<Component> &components,
                                                     const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                                     std::span<const Atom> moleculeAtoms,
                                                     std::span<AtomDynamics> moleculeDynamics,
                                                     const CrossLinkTable &table)
{
  RunningEnergy energy{};
  if (table.empty()) return energy;

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site, std::size_t localAtom) -> const Atom &
  { return moleculeAtoms[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom]; };
  auto addGradient = [&](const CrossLinkSite &site, std::size_t localAtom, const double3 &g)
  { moleculeDynamics[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom].gradient += g; };

  for (const CrossLink &link : table.links)
  {
    evaluateLink<1>(forceField, simulationBox, components, table, link, atomAt, energy, addGradient, nullptr);
  }
  return energy;
}

std::pair<RunningEnergy, double3x3> Interactions::computeCrossLinkEnergyStrainDerivative(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    std::span<AtomDynamics> moleculeDynamics, const CrossLinkTable &table)
{
  RunningEnergy energy{};
  double3x3 strain{};
  if (table.empty()) return {energy, strain};

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site, std::size_t localAtom) -> const Atom &
  { return moleculeAtoms[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom]; };
  auto addGradient = [&](const CrossLinkSite &site, std::size_t localAtom, const double3 &g)
  { moleculeDynamics[layout.moleculeStart(site.componentId, site.moleculeIndex) + localAtom].gradient += g; };

  for (const CrossLink &link : table.links)
  {
    evaluateLink<1>(forceField, simulationBox, components, table, link, atomAt, energy, addGradient, &strain);
  }
  return {energy, strain};
}
