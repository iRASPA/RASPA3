module;

module cross_links;

import std;

import archive;
import json;
import units;
import bond_potential;
import bend_potential;

std::optional<std::size_t> CrossLinkTable::findBondType(const std::string &siteTypeA,
                                                        const std::string &siteTypeB) const
{
  for (std::size_t i = 0; i < bondTypes.size(); ++i)
  {
    if (bondTypes[i].matches(siteTypeA, siteTypeB)) return i;
  }
  return std::nullopt;
}

std::span<const std::size_t> CrossLinkTable::linksOfMolecule(std::size_t componentId, std::size_t moleculeIndex) const
{
  CrossLinkSite site{static_cast<std::uint32_t>(componentId), static_cast<std::uint32_t>(moleculeIndex), 0};
  auto it = linksPerMolecule.find(site.moleculeKey());
  if (it == linksPerMolecule.end()) return {};
  return it->second;
}

bool CrossLinkTable::moleculeIsLinked(std::size_t componentId, std::size_t moleculeIndex) const
{
  return !linksOfMolecule(componentId, moleculeIndex).empty();
}

std::uint32_t CrossLinkTable::linkCount(const CrossLinkSite &site) const
{
  auto it = linkCountPerSite.find(site.siteKey());
  return it == linkCountPerSite.end() ? 0u : it->second;
}

bool CrossLinkTable::isLinked(const CrossLinkSite &a, const CrossLinkSite &b) const
{
  for (std::size_t id : linksOfMolecule(a.componentId, a.moleculeIndex))
  {
    const CrossLink &link = links[id];
    if ((link.a == a && link.b == b) || (link.a == b && link.b == a)) return true;
  }
  return false;
}

std::size_t CrossLinkTable::addLink(const CrossLink &link)
{
  std::size_t id = links.size();
  links.push_back(link);
  linksPerMolecule[link.a.moleculeKey()].push_back(id);
  linksPerMolecule[link.b.moleculeKey()].push_back(id);
  ++linkCountPerSite[link.a.siteKey()];
  ++linkCountPerSite[link.b.siteKey()];
  return id;
}

void CrossLinkTable::removeLink(std::size_t id)
{
  const CrossLink removed = links[id];

  auto eraseFromMolecule = [&](const CrossLinkSite &site, std::size_t linkId)
  {
    auto it = linksPerMolecule.find(site.moleculeKey());
    if (it == linksPerMolecule.end()) return;
    std::erase(it->second, linkId);
    if (it->second.empty()) linksPerMolecule.erase(it);
  };
  auto decrementSite = [&](const CrossLinkSite &site)
  {
    auto it = linkCountPerSite.find(site.siteKey());
    if (it == linkCountPerSite.end()) return;
    if (--it->second == 0) linkCountPerSite.erase(it);
  };

  eraseFromMolecule(removed.a, id);
  eraseFromMolecule(removed.b, id);
  decrementSite(removed.a);
  decrementSite(removed.b);

  std::size_t last = links.size() - 1;
  if (id != last)
  {
    // The last link moves into the freed slot: repoint its molecule-index entries.
    const CrossLink moved = links[last];
    links[id] = moved;
    for (const CrossLinkSite &site : {moved.a, moved.b})
    {
      auto it = linksPerMolecule.find(site.moleculeKey());
      if (it == linksPerMolecule.end()) continue;
      for (std::size_t &entry : it->second)
      {
        if (entry == last) entry = id;
      }
    }
  }
  links.pop_back();
}

void CrossLinkTable::rebuildIndex()
{
  linksPerMolecule.clear();
  linkCountPerSite.clear();
  for (std::size_t id = 0; id < links.size(); ++id)
  {
    const CrossLink &link = links[id];
    linksPerMolecule[link.a.moleculeKey()].push_back(id);
    linksPerMolecule[link.b.moleculeKey()].push_back(id);
    ++linkCountPerSite[link.a.siteKey()];
    ++linkCountPerSite[link.b.siteKey()];
  }
}

void CrossLinkTable::moleculeDeleted(std::size_t componentId, std::size_t moleculeIndex)
{
  if (links.empty()) return;
  bool changed = false;
  for (CrossLink &link : links)
  {
    for (CrossLinkSite *site : {&link.a, &link.b})
    {
      if (site->componentId != componentId) continue;
      if (site->moleculeIndex == moleculeIndex)
      {
        throw std::runtime_error(
            std::format("[CrossLinkTable]: molecule {} of component {} was deleted while carrying a cross-link\n",
                        moleculeIndex, componentId));
      }
      if (site->moleculeIndex > moleculeIndex)
      {
        --site->moleculeIndex;
        changed = true;
      }
    }
  }
  if (changed) rebuildIndex();
}

void CrossLinkTable::moleculeInserted(std::size_t componentId, std::size_t moleculeIndex)
{
  if (links.empty()) return;
  bool changed = false;
  for (CrossLink &link : links)
  {
    for (CrossLinkSite *site : {&link.a, &link.b})
    {
      if (site->componentId != componentId) continue;
      if (site->moleculeIndex >= moleculeIndex)
      {
        ++site->moleculeIndex;
        changed = true;
      }
    }
  }
  if (changed) rebuildIndex();
}

void CrossLinkTable::swapTopology(CrossLinkTable &other)
{
  std::swap(links, other.links);
  std::swap(linksPerMolecule, other.linksPerMolecule);
  std::swap(linkCountPerSite, other.linkCountPerSite);
}

std::string CrossLinkTable::printStatus() const
{
  if (!enabled()) return {};

  std::ostringstream stream;
  std::print(stream, "Cross-links\n");
  std::print(stream, "========================================================================================================================\n\n");
  std::print(stream, "Number of cross-link bond types: {}\n", bondTypes.size());
  for (std::size_t i = 0; i < bondTypes.size(); ++i)
  {
    const CrossLinkBondType &type = bondTypes[i];
    std::print(stream, "  type {}: sites '{}'-'{}'\n", i, type.siteTypeA, type.siteTypeB);
    std::print(stream, "    bond:              {}", type.bond.print());
    if (type.junctionBend.has_value())
    {
      std::print(stream, "    junction bend:     {}", type.junctionBend->print());
    }
    std::print(stream, "    capture radius:    {:g} [Å]\n", type.captureRadius);
    std::print(stream, "    formation energy:  {:g} [K]\n", type.formationEnergy * Units::EnergyToKelvin);
  }
  std::print(stream, "Number of cross-links: {}\n\n", links.size());
  return stream.str();
}

nlohmann::json CrossLinkTable::jsonStatus() const
{
  nlohmann::json status;
  if (!enabled()) return status;

  nlohmann::json types = nlohmann::json::array();
  for (const CrossLinkBondType &type : bondTypes)
  {
    nlohmann::json item;
    item["sites"] = {type.siteTypeA, type.siteTypeB};
    item["captureRadius"] = type.captureRadius;
    item["formationEnergy"] = type.formationEnergy * Units::EnergyToKelvin;
    types.push_back(item);
  }
  status["bondTypes"] = types;
  status["numberOfLinks"] = links.size();

  nlohmann::json linkList = nlohmann::json::array();
  for (const CrossLink &link : links)
  {
    linkList.push_back({link.a.componentId, link.a.moleculeIndex, link.a.atomIndex, link.b.componentId,
                        link.b.moleculeIndex, link.b.atomIndex});
  }
  status["links"] = linkList;
  return status;
}

// ---------------------------------------------------------------------------------------------------
// Serialization
// ---------------------------------------------------------------------------------------------------

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactiveSite &s)
{
  archive << s.versionNumber;
  archive << s.atom;
  archive << s.siteType;
  archive << s.valence;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactiveSite &s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ReactiveSite' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> s.atom;
  archive >> s.siteType;
  archive >> s.valence;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkBondType &b)
{
  archive << b.versionNumber;
  archive << b.siteTypeA;
  archive << b.siteTypeB;
  archive << b.bond;
  archive << b.junctionBend;
  archive << b.captureRadius;
  archive << b.formationEnergy;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkBondType &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CrossLinkBondType' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> b.siteTypeA;
  archive >> b.siteTypeB;
  archive >> b.bond;
  archive >> b.junctionBend;
  archive >> b.captureRadius;
  archive >> b.formationEnergy;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkSite &s)
{
  archive << s.componentId;
  archive << s.moleculeIndex;
  archive << s.atomIndex;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkSite &s)
{
  archive >> s.componentId;
  archive >> s.moleculeIndex;
  archive >> s.atomIndex;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLink &l)
{
  archive << l.a;
  archive << l.b;
  archive << l.bondTypeId;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLink &l)
{
  archive >> l.a;
  archive >> l.b;
  archive >> l.bondTypeId;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkTable &t)
{
  archive << t.versionNumber;
  archive << t.bondTypes;
  archive << t.links;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkTable &t)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CrossLinkTable' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> t.bondTypes;
  archive >> t.links;
  t.rebuildIndex();

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CrossLinkTable: Error in binary restart\n"));
  }
#endif

  return archive;
}
