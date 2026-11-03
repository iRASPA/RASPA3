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

// Which terms of an excluded pair a sum over links holds.
enum class ExclusionTerms
{
  None,       // the pair is held by another link (or tether)
  Full,       // the inter-molecular sum counted the pair: its pair interaction is removed, the exclusion term added
  EwaldOnly,  // the pair was left out of the inter-molecular sum (a CBMC regrowth, see CBMC::OverlapExemption)
};

// Removes the non-bonded pair interaction of an excluded pair (x, y) that the inter-molecular sum has
// counted (unless 'pairCounted' is false), and adds what an intramolecular pair receives instead (the Ewald
// exclusion or the shifted-potential completion). 'addGradient(onX, g)' receives the gradient on x (and -g on y).
template <std::size_t Order, typename GradientSink>
void excludePair(const ForceField &forceField, const SimulationBox &simulationBox, const Atom &x, const Atom &y,
                 RunningEnergy &energy, GradientSink &&addGradient, double3x3 *strain, bool pairCounted = true)
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

  if (!forceField.omitInterInteractions && pairCounted)
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
  else if (forceField.usesRealSpaceChargeCorrections() && !forceField.omitInterInteractions && pairCounted &&
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

// The bend between two links of the reactive site 'atom' of a component, if the site declares one.
const BendPotential *linkBendOf(const Component &component, std::size_t atom)
{
  for (const ReactiveSite &site : component.reactiveSites)
  {
    if (site.atom == atom) return site.linkBend.has_value() ? &site.linkBend.value() : nullptr;
  }
  return nullptr;
}

// The terms of one link between site 'siteA' of molecule A and site 'siteB' of molecule B.
// 'atomA(local)' / 'atomB(local)' return the atoms of the two molecules by local index;
// 'addGradientA(local, g)' / 'addGradientB(local, g)' receive the gradients (Order >= 1 only).
// 'neighboursA' / 'neighboursB' are the intramolecular neighbours of the two sites.
// 'takesExclusion(localA, localB)' tells which terms (ExclusionTerms) of the exclusion of a pair across the link
// these terms hold: several links can reach the same pair, which the inter-molecular sum counted once.
template <std::size_t Order, typename AtomA, typename AtomB, typename GradientA, typename GradientB,
          typename ExclusionFilter>
void evaluateLinkTerms(const ForceField &forceField, const SimulationBox &simulationBox,
                       const CrossLinkBondType &bondType, std::size_t siteA, std::size_t siteB,
                       std::span<const std::size_t> neighboursA, std::span<const std::size_t> neighboursB,
                       AtomA &&atomA, AtomB &&atomB, RunningEnergy &energy, GradientA &&addGradientA,
                       GradientB &&addGradientB, ExclusionFilter &&takesExclusion, double3x3 *strain)
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
    const ExclusionTerms terms = takesExclusion(localA, localB);
    if (terms == ExclusionTerms::None) return;
    excludePair<Order>(
        forceField, simulationBox, atomA(localA), atomB(localB), energy,
        [&](const double3 &g)
        {
          addGradientA(localA, g);
          addGradientB(localB, -g);
        },
        strain, terms == ExclusionTerms::Full);
  };
  exclude(siteA, siteB);
  for (std::size_t nA : neighboursA) exclude(nA, siteB);
  for (std::size_t nB : neighboursB) exclude(siteA, nB);
}

// Evaluates a set of links as if they were the only links: the terms of every link, and the terms that two links
// at one site form together (the bend of the site between its two partners, and the 1-3 exclusion of the two
// partners when they are on different molecules), with every excluded pair once however many links reach it.
// 'atomAt(site)' returns the atom that 'site' names (component, molecule, atom index); 'addGradient(site, g)'
// receives the gradients (Order >= 1 only).
template <std::size_t Order, typename AtomAccessor, typename GradientSink>
RunningEnergy evaluateLinkSet(const ForceField &forceField, const SimulationBox &simulationBox,
                              const std::vector<Component> &components, const CrossLinkTable &table,
                              std::span<const CrossLink> links, AtomAccessor &&atomAt, GradientSink &&addGradient,
                              double3x3 *strain)
{
  RunningEnergy energy{};
  if (links.empty()) return energy;

  auto atomOf = [](const CrossLinkSite &site, std::size_t local)
  { return CrossLinkSite{site.componentId, site.moleculeIndex, static_cast<std::uint32_t>(local)}; };

  std::set<std::pair<std::uint64_t, std::uint64_t>> excluded{};
  auto takesExclusion = [&](const CrossLinkSite &x, const CrossLinkSite &y)
  { return excluded.insert(std::minmax(x.siteKey(), y.siteKey())).second; };

  // the partners of every site of the set
  std::map<CrossLinkSite, std::vector<CrossLinkSite>> partnersOfSite{};

  for (const CrossLink &link : links)
  {
    const std::vector<std::size_t> neighboursA = neighboursOf(components[link.a.componentId], link.a.atomIndex);
    const std::vector<std::size_t> neighboursB = neighboursOf(components[link.b.componentId], link.b.atomIndex);
    evaluateLinkTerms<Order>(
        forceField, simulationBox, table.bondTypes[link.bondTypeId], link.a.atomIndex, link.b.atomIndex,
        neighboursA, neighboursB, [&](std::size_t local) -> const Atom & { return atomAt(atomOf(link.a, local)); },
        [&](std::size_t local) -> const Atom & { return atomAt(atomOf(link.b, local)); }, energy,
        [&](std::size_t local, const double3 &g) { addGradient(atomOf(link.a, local), g); },
        [&](std::size_t local, const double3 &g) { addGradient(atomOf(link.b, local), g); },
        [&](std::size_t localA, std::size_t localB)
        {
          return takesExclusion(atomOf(link.a, localA), atomOf(link.b, localB)) ? ExclusionTerms::Full
                                                                                : ExclusionTerms::None;
        },
        strain);
    partnersOfSite[link.a].push_back(link.b);
    partnersOfSite[link.b].push_back(link.a);
  }

  // The terms through two links of one site y: for the partners x and z, the bend (x, y, z) and the exclusion of
  // the pair (x, z), which the inter-molecular sum counted when x and z are on different molecules. When x and z
  // are on the same molecule the pair keeps the treatment of that molecule's own topology.
  for (const auto &[y, partners] : partnersOfSite)
  {
    if (partners.size() < 2) continue;
    const BendPotential *bend = linkBendOf(components[y.componentId], y.atomIndex);
    const double3 posY = atomAt(y).position;
    auto imageNextToY = [&](const CrossLinkSite &site)
    { return posY + simulationBox.applyPeriodicBoundaryConditions(atomAt(site).position - posY); };
    for (std::size_t i = 0; i < partners.size(); ++i)
    {
      for (std::size_t j = i + 1; j < partners.size(); ++j)
      {
        const CrossLinkSite &x = partners[i];
        const CrossLinkSite &z = partners[j];
        if (x == z) continue;
        if (bend != nullptr)
        {
          const double3 posX = imageNextToY(x);
          const double3 posZ = imageNextToY(z);
          if constexpr (Order == 0)
          {
            energy.crossLink += bend->calculateEnergy(posX, posY, posZ, std::nullopt);
          }
          else
          {
            const auto [bendEnergy, gradients, bendStrain] = bend->potentialEnergyGradientStrain(posX, posY, posZ);
            energy.crossLink += bendEnergy;
            addGradient(x, gradients[0]);
            addGradient(y, gradients[1]);
            addGradient(z, gradients[2]);
            if (strain) *strain += bendStrain;
          }
        }
        if (!x.sameMolecule(z) && takesExclusion(x, z))
        {
          excludePair<Order>(
              forceField, simulationBox, atomAt(x), atomAt(z), energy,
              [&](const double3 &g)
              {
                addGradient(x, g);
                addGradient(z, -g);
              },
              strain);
        }
      }
    }
  }
  return energy;
}

using MoleculeId = std::pair<std::uint32_t, std::uint32_t>;

MoleculeId moleculeOf(const CrossLinkSite &site) { return {site.componentId, site.moleculeIndex}; }

// The links whose terms can involve, or share an excluded pair with, a term of the links of the given molecules:
// the links of those molecules and the links of the molecules they are linked to (a pair reached through two links
// joins a molecule to the partner of its partner). Each link once.
std::vector<CrossLink> linksAround(const CrossLinkTable &table, const std::set<MoleculeId> &molecules)
{
  std::set<MoleculeId> around = molecules;
  for (const MoleculeId &molecule : molecules)
  {
    for (std::size_t id : table.linksOfMolecule(molecule.first, molecule.second))
    {
      const CrossLink &link = table.links[id];
      around.insert(moleculeOf(link.a));
      around.insert(moleculeOf(link.b));
    }
  }
  std::set<std::size_t> ids{};
  for (const MoleculeId &molecule : around)
  {
    for (std::size_t id : table.linksOfMolecule(molecule.first, molecule.second)) ids.insert(id);
  }
  std::vector<CrossLink> links{};
  links.reserve(ids.size());
  for (std::size_t id : ids) links.push_back(table.links[id]);
  return links;
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
  if (links.empty()) return {};

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site) -> const Atom & { return moleculeAtoms[layout.index(site)]; };
  auto noGradient = [](const CrossLinkSite &, const double3 &) {};
  return evaluateLinkSet<0>(forceField, simulationBox, components, table, links, atomAt, noGradient, nullptr);
}

RunningEnergy Interactions::computeCrossLinkTopologyChangeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    const CrossLinkTable &table, std::span<const CrossLink> removed, std::span<const CrossLink> added)
{
  // The links around the molecules of the changed links enter both sides; the terms they form among themselves
  // cancel, the terms they form with a changed link, or an excluded pair they share with one, do not.
  std::set<MoleculeId> molecules{};
  for (const CrossLink &link : removed) molecules.insert({moleculeOf(link.a), moleculeOf(link.b)});
  for (const CrossLink &link : added) molecules.insert({moleculeOf(link.a), moleculeOf(link.b)});

  std::vector<CrossLink> before = linksAround(table, molecules);
  std::erase_if(before, [&](const CrossLink &link)
                { return std::ranges::contains(removed, link) || std::ranges::contains(added, link); });
  std::vector<CrossLink> after = before;
  before.insert(before.end(), removed.begin(), removed.end());
  after.insert(after.end(), added.begin(), added.end());

  return computeCrossLinkEnergyOfLinks(forceField, simulationBox, components, numberOfMoleculesPerComponent,
                                       moleculeAtoms, table, after) -
         computeCrossLinkEnergyOfLinks(forceField, simulationBox, components, numberOfMoleculesPerComponent,
                                       moleculeAtoms, table, before);
}

RunningEnergy Interactions::computeCrossLinkEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    const CrossLinkTable &table, std::size_t componentId, std::size_t moleculeIndex, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms)
{
  if (table.empty()) return {};
  if (!table.moleculeIsLinked(componentId, moleculeIndex)) return {};

  // the links of the molecule and of its partners: the terms that involve an atom of the molecule, and the
  // excluded pairs those terms share with other links (the terms without an atom of the molecule cancel)
  const std::vector<CrossLink> links =
      linksAround(table, {MoleculeId{static_cast<std::uint32_t>(componentId), static_cast<std::uint32_t>(moleculeIndex)}});

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto isMoved = [&](const CrossLinkSite &site)
  { return site.componentId == componentId && site.moleculeIndex == moleculeIndex; };
  auto atomNew = [&](const CrossLinkSite &site) -> const Atom &
  { return isMoved(site) ? newAtoms[site.atomIndex] : moleculeAtoms[layout.index(site)]; };
  auto atomOld = [&](const CrossLinkSite &site) -> const Atom &
  { return isMoved(site) ? oldAtoms[site.atomIndex] : moleculeAtoms[layout.index(site)]; };
  auto noGradient = [](const CrossLinkSite &, const double3 &) {};

  return evaluateLinkSet<0>(forceField, simulationBox, components, table, links, atomNew, noGradient, nullptr) -
         evaluateLinkSet<0>(forceField, simulationBox, components, table, links, atomOld, noGradient, nullptr);
}

namespace
{
// Whether the van der Waals energy of a pair is above the overlap criterion: a CBMC regrowth leaves such an
// excluded pair out of its external energy (CBMC::OverlapExemption), and the tether out of the exclusion.
bool aboveOverlapCriterion(const ForceField &forceField, const SimulationBox &simulationBox, const Atom &x,
                           const Atom &y)
{
  const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  bool above = false;
  Interactions::evaluatePair<0>(
      forceField, simulationBox, x, y, cutOffVDWSquared, cutOffChargeSquared, false,
      [&](const Potentials::PairDerivatives<0> &factors, const double3 &)
      { above = factors.energy > forceField.energyOverlapCriteria; },
      [](const Potentials::PairDerivatives<0> &, const double3 &) {});
  return above;
}
}  // namespace

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
  // (atom of this molecule, partner component, partner molecule, atom of the partner): each pair held by one tether
  std::set<std::tuple<std::size_t, std::uint32_t, std::uint32_t, std::size_t>> excluded{};
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

    auto holdExclusion = [&](std::size_t local, std::size_t partnerLocal)
    {
      if (excluded.insert({local, partner.componentId, partner.moleculeIndex, partnerLocal}).second)
      {
        tether.exclusions.push_back({local, partnerLocal, partnerStart + partnerLocal});
      }
    };
    holdExclusion(own.atomIndex, partner.atomIndex);
    for (std::size_t n : neighboursOf(components[componentId], own.atomIndex)) holdExclusion(n, partner.atomIndex);
    for (const auto &[n, atom] : tether.partnerNeighbours) holdExclusion(own.atomIndex, n);

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
        noGradient, noGradient,
        [&](std::size_t local, std::size_t partnerLocal)
        {
          const bool held =
              std::ranges::any_of(tether.exclusions, [&](const CrossLinkTether::Exclusion &exclusion)
                                  { return exclusion.atom == local && exclusion.partnerAtom == partnerLocal; });
          if (!held) return ExclusionTerms::None;
          // the regrowth left a pair above the overlap criterion out of the external energy
          return aboveOverlapCriterion(forceField, simulationBox, moleculeAtoms[local], atomB(partnerLocal))
                     ? ExclusionTerms::EwaldOnly
                     : ExclusionTerms::Full;
        },
        nullptr);
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
  if (table.empty()) return {};

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site) -> const Atom & { return moleculeAtoms[layout.index(site)]; };
  auto addGradient = [&](const CrossLinkSite &site, const double3 &g)
  { moleculeDynamics[layout.index(site)].gradient += g; };
  return evaluateLinkSet<1>(forceField, simulationBox, components, table, table.links, atomAt, addGradient,
                            nullptr);
}

std::pair<RunningEnergy, double3x3> Interactions::computeCrossLinkEnergyStrainDerivative(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    std::span<AtomDynamics> moleculeDynamics, const CrossLinkTable &table)
{
  double3x3 strain{};
  if (table.empty()) return {RunningEnergy{}, strain};

  const CrossLinkAtomLayout layout = CrossLinkAtomLayout::make(components, numberOfMoleculesPerComponent);
  auto atomAt = [&](const CrossLinkSite &site) -> const Atom & { return moleculeAtoms[layout.index(site)]; };
  auto addGradient = [&](const CrossLinkSite &site, const double3 &g)
  { moleculeDynamics[layout.index(site)].gradient += g; };
  const RunningEnergy energy = evaluateLinkSet<1>(forceField, simulationBox, components, table, table.links, atomAt,
                                                  addGradient, &strain);
  return {energy, strain};
}
