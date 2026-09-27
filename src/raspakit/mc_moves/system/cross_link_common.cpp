module;

module mc_moves_cross_link_common;

import std;

import double3;
import atom;
import running_energy;
import system;
import component;
import cross_links;
import interactions_cross_link;

namespace
{
std::size_t reactiveSiteIndex(const Component &component, std::size_t atom)
{
  for (std::size_t i = 0; i < component.reactiveSites.size(); ++i)
  {
    if (component.reactiveSites[i].atom == atom) return i;
  }
  throw std::runtime_error(
      std::format("[Cross-link moves]: atom {} of component '{}' is not a reactive site\n", atom, component.name));
}
}  // namespace

const ReactiveSite &MC_Moves::CrossLinkCommon::reactiveSiteOf(const System &system, const CrossLinkSite &site)
{
  const Component &component = system.components[site.componentId];
  return component.reactiveSites[reactiveSiteIndex(component, site.atomIndex)];
}

double3 MC_Moves::CrossLinkCommon::positionOf(const System &system, const CrossLinkSite &site)
{
  return system.spanOfMolecule(site.componentId, site.moleculeIndex)[site.atomIndex].position;
}

bool MC_Moves::CrossLinkCommon::isFree(const System &system, const CrossLinkSite &site)
{
  return system.crossLinks.linkCount(site) < reactiveSiteOf(system, site).valence;
}

std::vector<CrossLinkSite> MC_Moves::CrossLinkCommon::freeSites(const System &system)
{
  std::vector<CrossLinkSite> sites;
  for (std::size_t c = 0; c < system.components.size(); ++c)
  {
    const Component &component = system.components[c];
    if (component.reactiveSites.empty()) continue;
    for (std::size_t m = 0; m < system.numberOfIntegerMoleculesPerComponent[c]; ++m)
    {
      for (const ReactiveSite &reactiveSite : component.reactiveSites)
      {
        CrossLinkSite site{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(m),
                           static_cast<std::uint32_t>(reactiveSite.atom)};
        if (system.crossLinks.linkCount(site) < reactiveSite.valence) sites.push_back(site);
      }
    }
  }
  return sites;
}

std::optional<std::size_t> MC_Moves::CrossLinkCommon::bondTypeFor(const System &system, const CrossLinkSite &a,
                                                                  const CrossLinkSite &b)
{
  return system.crossLinks.findBondType(reactiveSiteOf(system, a).siteType, reactiveSiteOf(system, b).siteType);
}

bool MC_Moves::CrossLinkCommon::withinCaptureRadius(const System &system, const CrossLinkSite &a,
                                                    const CrossLinkSite &b, std::size_t bondTypeId)
{
  const double captureRadius = system.crossLinks.bondTypes[bondTypeId].captureRadius;
  const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(positionOf(system, a) - positionOf(system, b));
  return double3::dot(dr, dr) <= captureRadius * captureRadius;
}

namespace
{
// The sites that could be linked to 'pivot' (different molecule, matching bond type, within the capture
// radius, not linked to it), restricted to free sites ('wantLinked' false: the partners of a new link)
// or to linked sites ('wantLinked' true: the partners of a bond exchange).
std::vector<MC_Moves::CrossLinkCommon::Candidate> candidatesOf(const System &system, const CrossLinkSite &pivot,
                                                               bool wantLinked)
{
  using MC_Moves::CrossLinkCommon::Candidate;
  std::vector<Candidate> candidates;
  const std::string &pivotType = MC_Moves::CrossLinkCommon::reactiveSiteOf(system, pivot).siteType;
  const double3 pivotPosition = MC_Moves::CrossLinkCommon::positionOf(system, pivot);

  for (std::size_t c = 0; c < system.components.size(); ++c)
  {
    const Component &component = system.components[c];
    if (component.reactiveSites.empty()) continue;

    // Bond types per reactive site of this component (independent of the molecule).
    std::vector<std::optional<std::size_t>> bondTypeOfSite(component.reactiveSites.size());
    bool anyBondType = false;
    for (std::size_t s = 0; s < component.reactiveSites.size(); ++s)
    {
      bondTypeOfSite[s] = system.crossLinks.findBondType(pivotType, component.reactiveSites[s].siteType);
      anyBondType = anyBondType || bondTypeOfSite[s].has_value();
    }
    if (!anyBondType) continue;

    for (std::size_t m = 0; m < system.numberOfIntegerMoleculesPerComponent[c]; ++m)
    {
      if (c == pivot.componentId && m == pivot.moleculeIndex) continue;  // no intramolecular links
      const std::span<const Atom> atoms = system.spanOfMolecule(c, m);
      for (std::size_t s = 0; s < component.reactiveSites.size(); ++s)
      {
        if (!bondTypeOfSite[s].has_value()) continue;
        const ReactiveSite &reactiveSite = component.reactiveSites[s];
        const CrossLinkSite site{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(m),
                                 static_cast<std::uint32_t>(reactiveSite.atom)};
        const std::uint32_t linkCount = system.crossLinks.linkCount(site);
        if (wantLinked ? (linkCount == 0) : (linkCount >= reactiveSite.valence)) continue;
        if (system.crossLinks.isLinked(pivot, site)) continue;

        const double captureRadius = system.crossLinks.bondTypes[bondTypeOfSite[s].value()].captureRadius;
        const double3 dr =
            system.simulationBox.applyPeriodicBoundaryConditions(pivotPosition - atoms[reactiveSite.atom].position);
        if (double3::dot(dr, dr) > captureRadius * captureRadius) continue;

        candidates.push_back(Candidate{site, bondTypeOfSite[s].value()});
      }
    }
  }
  return candidates;
}
}  // namespace

std::vector<MC_Moves::CrossLinkCommon::Candidate> MC_Moves::CrossLinkCommon::partnerCandidates(
    const System &system, const CrossLinkSite &pivot)
{
  return candidatesOf(system, pivot, false);
}

std::vector<MC_Moves::CrossLinkCommon::Candidate> MC_Moves::CrossLinkCommon::linkedCandidates(
    const System &system, const CrossLinkSite &pivot)
{
  return candidatesOf(system, pivot, true);
}

std::vector<std::size_t> MC_Moves::CrossLinkCommon::linksOfSite(const System &system, const CrossLinkSite &site)
{
  std::vector<std::size_t> ids;
  for (std::size_t id : system.crossLinks.linksOfMolecule(site.componentId, site.moleculeIndex))
  {
    const CrossLink &link = system.crossLinks.links[id];
    if (link.a == site || link.b == site) ids.push_back(id);
  }
  return ids;
}

RunningEnergy MC_Moves::CrossLinkCommon::linkEnergy(const System &system, const CrossLink &link)
{
  return Interactions::computeCrossLinkEnergyOfLinks(system.forceField, system.simulationBox, system.components,
                                                     system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(),
                                                     system.crossLinks, std::span<const CrossLink>(&link, 1));
}
