module;

module polymer_links;

import std;

import archive;
import json;
import double3;
import double3x3;
import simulationbox;
import atom;
import molecule;
import forcefield;
import component;
import reactive_group;
import polymer_reactions;
import bond_potential;
import bend_potential;
import torsion_potential;
import intra_molecular_exclusions;

// ---- PolymerAtomRef / PolymerBond ---------------------------------------------------------------------------------

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerAtomRef &r)
{
  archive << r.componentId;
  archive << r.moleculeIndex;
  archive << r.atomIndex;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerAtomRef &r)
{
  archive >> r.componentId;
  archive >> r.moleculeIndex;
  archive >> r.atomIndex;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerBond &b)
{
  archive << b.a;
  archive << b.b;
  archive << b.reaction;
  archive << b.formedAtStep;
  archive << b.formationDistance;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerBond &b)
{
  archive >> b.a;
  archive >> b.b;
  archive >> b.reaction;
  archive >> b.formedAtStep;
  archive >> b.formationDistance;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BrokenBond &b)
{
  archive << b.a;
  archive << b.b;
  archive << b.reaction;
  archive << b.brokenAtStep;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BrokenBond &b)
{
  archive >> b.a;
  archive >> b.b;
  archive >> b.reaction;
  archive >> b.brokenAtStep;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Retyping &r)
{
  archive << r.atom;
  archive << r.from;
  archive << r.to;
  archive << r.step;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Retyping &r)
{
  archive >> r.atom;
  archive >> r.from;
  archive >> r.to;
  archive >> r.step;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const StereoCentre &c)
{
  archive << c.atoms;
  archive << c.chirality;
  archive << c.bond;
  archive << c.group;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, StereoCentre &c)
{
  archive >> c.atoms;
  archive >> c.chirality;
  archive >> c.bond;
  archive >> c.group;
  return archive;
}

// ---- PolymerTemplates ---------------------------------------------------------------------------------------------

namespace
{
template <std::size_t N>
std::array<std::string, N> readTypeNames(const nlohmann::basic_json<nlohmann::raspa_map> &item,
                                         const ForceField &forceField, std::string_view what)
{
  if (!item.is_array() || item.size() != 3 || !item[0].is_array() || item[0].size() != N)
  {
    throw std::runtime_error(std::format(
        "[Input reader]: 'Polymerization' {} entry {} must be [[{} pseudo-atom names], \"TYPE\", [parameters]]\n",
        what, item.dump(), N));
  }
  std::array<std::string, N> names{};
  for (std::size_t i = 0; i < N; ++i)
  {
    if (!item[0][i].is_string())
    {
      throw std::runtime_error(std::format(
          "[Input reader]: 'Polymerization' {} entry {}: the atoms must be pseudo-atom names (the tables are "
          "type-based)\n",
          what, item.dump()));
    }
    names[i] = item[0][i].get<std::string>();
    if (!forceField.findPseudoAtom(names[i]).has_value())
    {
      throw std::runtime_error(std::format("[Input reader]: 'Polymerization' {} entry {}: unknown pseudo-atom '{}'\n",
                                           what, item.dump(), names[i]));
    }
  }
  return names;
}

std::vector<double> readParameters(const nlohmann::basic_json<nlohmann::raspa_map> &item)
{
  return item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};
}

template <std::size_t N>
std::array<std::uint16_t, N> reversed(const std::array<std::uint16_t, N> &key)
{
  std::array<std::uint16_t, N> result = key;
  std::ranges::reverse(result);
  return result;
}

template <std::size_t N>
std::array<std::uint32_t, N> canonical(std::array<std::uint32_t, N> key)
{
  std::array<std::uint32_t, N> other = key;
  std::ranges::reverse(other);
  return other < key ? other : key;
}
}  // namespace

void PolymerTemplates::read(const nlohmann::basic_json<nlohmann::raspa_map> &value, const ForceField &forceField)
{
  if (value.contains("Bonds"))
  {
    for (const auto &item : value["Bonds"])
    {
      try
      {
        BondTemplate t{};
        t.types = readTypeNames<2>(item, forceField, "'Bonds'");
        t.type = BondPotential::definitionForString.at(item[1].get<std::string>());
        t.parameters = readParameters(item);
        bonds.push_back(std::move(t));
      }
      catch (std::out_of_range const &)
      {
        throw std::runtime_error(
            std::format("[Input reader]: 'Polymerization' 'Bonds' entry {}: unknown bond potential\n", item.dump()));
      }
    }
  }
  if (value.contains("Bends"))
  {
    for (const auto &item : value["Bends"])
    {
      try
      {
        BendTemplate t{};
        t.types = readTypeNames<3>(item, forceField, "'Bends'");
        t.type = BendPotential::definitionForString.at(item[1].get<std::string>());
        t.parameters = readParameters(item);
        bends.push_back(std::move(t));
      }
      catch (std::out_of_range const &)
      {
        throw std::runtime_error(
            std::format("[Input reader]: 'Polymerization' 'Bends' entry {}: unknown bend potential\n", item.dump()));
      }
    }
  }
  if (value.contains("Torsions"))
  {
    for (const auto &item : value["Torsions"])
    {
      try
      {
        TorsionTemplate t{};
        t.types = readTypeNames<4>(item, forceField, "'Torsions'");
        t.type = TorsionPotential::definitionForString.at(item[1].get<std::string>());
        t.parameters = readParameters(item);
        torsions.push_back(std::move(t));
      }
      catch (std::out_of_range const &)
      {
        throw std::runtime_error(std::format(
            "[Input reader]: 'Polymerization' 'Torsions' entry {}: unknown torsion potential\n", item.dump()));
      }
    }
  }
  resolved = false;
}

void PolymerTemplates::resolve(const ForceField &forceField)
{
  bondIndex.clear();
  bendIndex.clear();
  torsionIndex.clear();
  auto typeOf = [&](const std::string &name) -> std::uint16_t
  {
    const std::optional<std::size_t> type = forceField.findPseudoAtom(name);
    if (!type.has_value())
    {
      throw std::runtime_error(
          std::format("[Polymerization]: unknown pseudo-atom '{}' in the bonded templates\n", name));
    }
    return static_cast<std::uint16_t>(type.value());
  };
  auto ambiguous = [](std::string_view what, std::span<const std::string> names)
  {
    std::string list{};
    for (const std::string &name : names) list += std::format("{} ", name);
    throw std::runtime_error(
        std::format("[Polymerization]: ambiguous '{}' template for the types {}(listed twice)\n", what, list));
  };
  for (std::size_t i = 0; i < bonds.size(); ++i)
  {
    const std::array<std::uint16_t, 2> key{typeOf(bonds[i].types[0]), typeOf(bonds[i].types[1])};
    if (!bondIndex.try_emplace(key, i).second && bondIndex.at(key) != i) ambiguous("Bonds", bonds[i].types);
    bondIndex.try_emplace(reversed(key), i);
  }
  for (std::size_t i = 0; i < bends.size(); ++i)
  {
    const std::array<std::uint16_t, 3> key{typeOf(bends[i].types[0]), typeOf(bends[i].types[1]),
                                           typeOf(bends[i].types[2])};
    if (!bendIndex.try_emplace(key, i).second && bendIndex.at(key) != i) ambiguous("Bends", bends[i].types);
    bendIndex.try_emplace(reversed(key), i);
  }
  for (std::size_t i = 0; i < torsions.size(); ++i)
  {
    const std::array<std::uint16_t, 4> key{typeOf(torsions[i].types[0]), typeOf(torsions[i].types[1]),
                                           typeOf(torsions[i].types[2]), typeOf(torsions[i].types[3])};
    if (!torsionIndex.try_emplace(key, i).second && torsionIndex.at(key) != i) ambiguous("Torsions", torsions[i].types);
    torsionIndex.try_emplace(reversed(key), i);
  }

  // the minima of the angular templates (the reference of the orientation strain): scanned numerically on a
  // synthetic geometry (the templates depend on the angle / dihedral only)
  bendMinimum.assign(bends.size(), 0.0);
  for (std::size_t i = 0; i < bends.size(); ++i)
  {
    const BendPotential potential({0, 1, 2}, bends[i].type, bends[i].parameters);
    double minimum = std::numeric_limits<double>::max();
    for (std::size_t k = 0; k <= 360; ++k)
    {
      const double theta = std::numbers::pi * static_cast<double>(k) / 360.0;
      minimum = std::min(minimum, potential.calculateEnergy(double3(std::cos(theta), std::sin(theta), 0.0),
                                                            double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0),
                                                            std::nullopt));
    }
    bendMinimum[i] = minimum;
  }
  torsionMinimum.assign(torsions.size(), 0.0);
  for (std::size_t i = 0; i < torsions.size(); ++i)
  {
    const TorsionPotential potential({0, 1, 2, 3}, torsions[i].type, torsions[i].parameters);
    double minimum = std::numeric_limits<double>::max();
    for (std::size_t k = 0; k < 720; ++k)
    {
      const double phi = 2.0 * std::numbers::pi * static_cast<double>(k) / 720.0;
      minimum = std::min(minimum, potential.calculateEnergy(double3(1.0, 0.0, 0.0), double3(0.0, 0.0, 0.0),
                                                            double3(0.0, 0.0, 1.0),
                                                            double3(std::cos(phi), std::sin(phi), 1.0)));
    }
    torsionMinimum[i] = minimum;
  }
  resolved = true;
}

std::optional<double> PolymerTemplates::bendStrain(const std::array<std::uint16_t, 3> &types, const double3 &a,
                                                   const double3 &b, const double3 &c) const
{
  const auto it = bendIndex.find(types);
  if (it == bendIndex.end()) return std::nullopt;
  const BendTemplate &t = bends[it->second];
  return BendPotential({0, 1, 2}, t.type, t.parameters).calculateEnergy(a, b, c, std::nullopt) -
         bendMinimum[it->second];
}

std::optional<double> PolymerTemplates::torsionStrain(const std::array<std::uint16_t, 4> &types, const double3 &a,
                                                      const double3 &b, const double3 &c, const double3 &d) const
{
  const auto it = torsionIndex.find(types);
  if (it == torsionIndex.end()) return std::nullopt;
  const TorsionTemplate &t = torsions[it->second];
  return TorsionPotential({0, 1, 2, 3}, t.type, t.parameters).calculateEnergy(a, b, c, d) -
         torsionMinimum[it->second];
}

std::optional<BondPotential> PolymerTemplates::bond(const std::array<std::uint16_t, 2> &types,
                                                    const std::array<std::size_t, 2> &identifiers) const
{
  const auto it = bondIndex.find(types);
  if (it == bondIndex.end()) return std::nullopt;
  const BondTemplate &t = bonds[it->second];
  return BondPotential(identifiers, t.type, t.parameters);
}

std::optional<BendPotential> PolymerTemplates::bend(const std::array<std::uint16_t, 3> &types,
                                                    const std::array<std::size_t, 3> &identifiers) const
{
  const auto it = bendIndex.find(types);
  if (it == bendIndex.end()) return std::nullopt;
  const BendTemplate &t = bends[it->second];
  return BendPotential(identifiers, t.type, t.parameters);
}

std::optional<TorsionPotential> PolymerTemplates::torsion(const std::array<std::uint16_t, 4> &types,
                                                          const std::array<std::size_t, 4> &identifiers) const
{
  const auto it = torsionIndex.find(types);
  if (it == torsionIndex.end()) return std::nullopt;
  const TorsionTemplate &t = torsions[it->second];
  return TorsionPotential(identifiers, t.type, t.parameters);
}

std::string PolymerTemplates::printStatus() const
{
  std::ostringstream stream;
  std::print(stream, "    product bonded templates:        {} bonds, {} bends, {} torsions\n", bonds.size(),
             bends.size(), torsions.size());
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerTemplates &t)
{
  archive << t.versionNumber;
  archive << static_cast<std::uint64_t>(t.bonds.size());
  for (const PolymerTemplates::BondTemplate &b : t.bonds)
  {
    archive << b.types;
    archive << b.type;
    archive << b.parameters;
  }
  archive << static_cast<std::uint64_t>(t.bends.size());
  for (const PolymerTemplates::BendTemplate &b : t.bends)
  {
    archive << b.types;
    archive << b.type;
    archive << b.parameters;
  }
  archive << static_cast<std::uint64_t>(t.torsions.size());
  for (const PolymerTemplates::TorsionTemplate &b : t.torsions)
  {
    archive << b.types;
    archive << b.type;
    archive << b.parameters;
  }
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerTemplates &t)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    throw std::runtime_error(std::format("Invalid version reading 'PolymerTemplates' at line {} in file {}\n",
                                         std::source_location::current().line(),
                                         std::source_location::current().file_name()));
  }
  std::uint64_t count{};
  archive >> count;
  t.bonds.resize(count);
  for (PolymerTemplates::BondTemplate &b : t.bonds)
  {
    archive >> b.types;
    archive >> b.type;
    archive >> b.parameters;
  }
  archive >> count;
  t.bends.resize(count);
  for (PolymerTemplates::BendTemplate &b : t.bends)
  {
    archive >> b.types;
    archive >> b.type;
    archive >> b.parameters;
  }
  archive >> count;
  t.torsions.resize(count);
  for (PolymerTemplates::TorsionTemplate &b : t.torsions)
  {
    archive >> b.types;
    archive >> b.type;
    archive >> b.parameters;
  }
  t.resolved = false;
  return archive;
}

// ---- ChiralVolume -------------------------------------------------------------------------------------------------

ChiralVolume ChiralVolume::compute(const double3 &bp, const double3 &bn, const double3 &bz)
{
  ChiralVolume result{};
  const double lp = std::sqrt(double3::dot(bp, bp));
  const double ln = std::sqrt(double3::dot(bn, bn));
  const double lz = std::sqrt(double3::dot(bz, bz));
  const double norm = lp * ln * lz;
  if (norm <= 0.0) return result;
  const double3 pn = double3::cross(bp, bn);
  const double3 nz = double3::cross(bn, bz);
  const double3 zp = double3::cross(bz, bp);
  const double value = double3::dot(pn, bz) / norm;
  result.value = value;
  // d(V/N)/db_i = (dV/db_i)/N - V b_i/|b_i|^2  with dV/db_p = b_n x b_z, dV/db_n = b_z x b_p, dV/db_z = b_p x b_n
  result.gradient[1] = nz / norm - (value / (lp * lp)) * bp;
  result.gradient[2] = zp / norm - (value / (ln * ln)) * bn;
  result.gradient[3] = pn / norm - (value / (lz * lz)) * bz;
  result.gradient[0] = -(result.gradient[1] + result.gradient[2] + result.gradient[3]);
  return result;
}

// ---- PolymerLinkTable: layout and groups --------------------------------------------------------------------------

PolymerAtomRef PolymerLinkTable::Layout::ref(std::size_t atom) const
{
  for (std::size_t c = componentOffset.size(); c-- > 0;)
  {
    if (atom >= componentOffset[c] && atomsPerMolecule[c] > 0)
    {
      const std::size_t local = atom - componentOffset[c];
      return PolymerAtomRef{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(local / atomsPerMolecule[c]),
                            static_cast<std::uint32_t>(local % atomsPerMolecule[c])};
    }
  }
  throw std::runtime_error("[Polymerization]: atom index outside of the molecule layout\n");
}

std::size_t PolymerLinkTable::Layout::moleculeOf(std::size_t atom) const
{
  const PolymerAtomRef r = ref(atom);
  return moleculeOffset[r.componentId] + r.moleculeIndex;
}

PolymerLinkTable::Layout PolymerLinkTable::makeLayout(
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent) const
{
  Layout result{};
  result.componentOffset.resize(components.size());
  result.atomsPerMolecule.resize(components.size());
  result.moleculeOffset.resize(components.size());
  std::size_t atomOffset = 0;
  std::size_t moleculeOffset = 0;
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    result.componentOffset[c] = atomOffset;
    result.moleculeOffset[c] = moleculeOffset;
    result.atomsPerMolecule[c] = components[c].atoms.size();
    atomOffset += numberOfMoleculesPerComponent[c] * components[c].atoms.size();
    moleculeOffset += numberOfMoleculesPerComponent[c];
  }
  result.numberOfAtoms = atomOffset;
  result.numberOfMolecules = moleculeOffset;
  return result;
}

std::size_t PolymerLinkTable::totalReactions() const
{
  return std::accumulate(reactionCounts.begin(), reactionCounts.end(), 0uz);
}

void PolymerLinkTable::activate(const std::vector<Component> &components,
                                const std::vector<std::size_t> &numberOfMoleculesPerComponent)
{
  if (!scheme.resolved)
  {
    throw std::runtime_error("[Polymerization]: the reaction scheme is not resolved (no 'Reactions' read)\n");
  }
  // the per-group type and initial state (derived; also after a restart)
  auto fillTypes = [&]()
  {
    groupType.assign(groupState.size(), 0);
    groupInitial.assign(groupState.size(), 0);
    for (std::size_t c = 0; c < components.size(); ++c)
    {
      for (std::size_t m = 0; m < moleculesPerComponent[c]; ++m)
      {
        for (std::size_t g = 0; g < groupsPerMolecule[c]; ++g)
        {
          const ReactiveGroup &definition = components[c].reactiveGroups[g];
          const std::uint16_t type = scheme.typeIndex(definition.type).value();
          const std::size_t index = groupIndex(c, m, g);
          groupType[index] = type;
          groupInitial[index] = scheme.stateIndex(type, definition.initialState).value();
        }
      }
    }
  };
  if (activated)
  {
    if (moleculesPerComponent != numberOfMoleculesPerComponent)
    {
      throw std::runtime_error(
          "[Polymerization]: the number of molecules changed after the polymerization was activated\n");
    }
    layout = makeLayout(components, numberOfMoleculesPerComponent);
    if (groupChirality.size() != groupState.size()) groupChirality.assign(groupState.size(), 0);
    if (previousGroup.size() != groupState.size()) previousGroup.assign(groupState.size(), noGroup);
    if (reactionCounts.size() != scheme.rules.size()) reactionCounts.resize(scheme.rules.size(), 0);
    fillTypes();
    return;
  }
  groupsPerMolecule.resize(components.size());
  groupOffset.resize(components.size());
  moleculesPerComponent = numberOfMoleculesPerComponent;
  std::size_t offset = 0;
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    groupsPerMolecule[c] = components[c].reactiveGroups.size();
    groupOffset[c] = offset;
    offset += groupsPerMolecule[c] * numberOfMoleculesPerComponent[c];
  }
  groupState.assign(offset, 0);
  groupStep.assign(offset, 0);
  groupChirality.assign(offset, 0);
  previousGroup.assign(offset, noGroup);
  reactionCounts.assign(scheme.rules.size(), 0);
  fillTypes();
  groupState = groupInitial;
  layout = makeLayout(components, numberOfMoleculesPerComponent);
  activated = true;
  derivedValid = false;
}

PolymerAtomRef PolymerLinkTable::groupRef(std::size_t index) const
{
  for (std::size_t c = groupOffset.size(); c-- > 0;)
  {
    if (index >= groupOffset[c] && groupsPerMolecule[c] > 0)
    {
      const std::size_t local = index - groupOffset[c];
      return PolymerAtomRef{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(local / groupsPerMolecule[c]),
                            static_cast<std::uint32_t>(local % groupsPerMolecule[c])};
    }
  }
  throw std::runtime_error("[Polymerization]: group index outside of the group layout\n");
}

const ReactiveGroup &PolymerLinkTable::definitionOf(std::size_t index, const std::vector<Component> &components) const
{
  const PolymerAtomRef r = groupRef(index);
  return components[r.componentId].reactiveGroups[r.atomIndex];
}

PolymerAtomRef PolymerLinkTable::atomRefOf(std::size_t index, std::uint16_t role,
                                           const std::vector<Component> &components) const
{
  PolymerAtomRef r = groupRef(index);
  const ReactiveGroup &definition = components[r.componentId].reactiveGroups[r.atomIndex];
  const std::optional<std::size_t> atom = definition.atomOf(scheme.types[groupType[index]].roles[role]);
  if (!atom.has_value())
  {
    throw std::runtime_error(std::format("[Polymerization]: group '{}' of component '{}' has no role '{}'\n",
                                         definition.type, components[r.componentId].name,
                                         scheme.types[groupType[index]].roles[role]));
  }
  r.atomIndex = static_cast<std::uint32_t>(atom.value());
  return r;
}

std::optional<std::size_t> PolymerLinkTable::groupOfAtom(std::uint32_t atom,
                                                         const std::vector<Component> &components) const
{
  const PolymerAtomRef r = layout.ref(atom);
  const std::vector<ReactiveGroup> &groups = components[r.componentId].reactiveGroups;
  for (std::size_t g = 0; g < groups.size(); ++g)
  {
    if (groups[g].contains(r.atomIndex)) return groupIndex(r.componentId, r.moleculeIndex, g);
  }
  return std::nullopt;
}

std::size_t PolymerLinkTable::countGroups(std::uint16_t type, std::uint8_t state) const
{
  std::size_t count = 0;
  for (std::size_t i = 0; i < groupState.size(); ++i)
  {
    if (groupType[i] == type && groupState[i] == state) ++count;
  }
  return count;
}

std::size_t PolymerLinkTable::countGroups(std::uint16_t type) const
{
  return static_cast<std::size_t>(std::ranges::count(groupType, type));
}

std::size_t PolymerLinkTable::countActive() const
{
  std::size_t count = 0;
  for (std::size_t i = 0; i < groupState.size(); ++i)
  {
    if (active(i)) ++count;
  }
  return count;
}

double PolymerLinkTable::conversion() const
{
  if (groupState.empty()) return 0.0;
  if (scheme.conversion.empty())
  {
    std::size_t changed = 0;
    for (std::size_t i = 0; i < groupState.size(); ++i) changed += reacted(i);
    return static_cast<double>(changed) / static_cast<double>(groupState.size());
  }
  std::size_t counted = 0, total = 0;
  for (std::size_t i = 0; i < groupState.size(); ++i)
  {
    for (const ReactionScheme::Conversion &entry : scheme.conversion)
    {
      if (entry.typeIndex != groupType[i]) continue;
      ++total;
      if (entry.stateMask[groupState[i]]) ++counted;
    }
  }
  return total > 0 ? static_cast<double>(counted) / static_cast<double>(total) : 0.0;
}

std::optional<std::size_t> PolymerLinkTable::centreOf(std::size_t group) const
{
  for (std::size_t k = centres.size(); k-- > 0;)
  {
    if (centres[k].group == group) return k;
  }
  return std::nullopt;
}

std::array<std::uint32_t, 4> PolymerLinkTable::stereoAtoms(const ReactionRule &rule, std::size_t k,
                                                           std::span<const std::size_t> groups,
                                                           const std::vector<Component> &components) const
{
  const ReactionRule::ResolvedStereo &stereo = rule.stereoAtoms[k];
  auto resolve = [&](const RoleIndex &ref) { return atomOf(groups[ref.reactant], ref.role, components); };
  return {resolve(stereo.centre), resolve(stereo.substituents[0]), resolve(stereo.substituents[1]),
          resolve(stereo.substituents[2])};
}

double PolymerLinkTable::chiralVolumeOf(const std::array<std::uint32_t, 4> &atoms, const SimulationBox &simulationBox,
                                        std::span<const Atom> systemAtoms) const
{
  const double3 centre = systemAtoms[atoms[0]].position;
  const double3 b1 = simulationBox.applyPeriodicBoundaryConditions(systemAtoms[atoms[1]].position - centre);
  const double3 b2 = simulationBox.applyPeriodicBoundaryConditions(systemAtoms[atoms[2]].position - centre);
  const double3 b3 = simulationBox.applyPeriodicBoundaryConditions(systemAtoms[atoms[3]].position - centre);
  return ChiralVolume::compute(b1, b2, b3).value;
}

std::int8_t PolymerLinkTable::geometricChirality(const std::array<std::uint32_t, 4> &atoms,
                                                 const SimulationBox &simulationBox,
                                                 std::span<const Atom> systemAtoms) const
{
  return chiralVolumeOf(atoms, simulationBox, systemAtoms) < 0.0 ? -1 : 1;
}

bool PolymerLinkTable::rampsActive() const
{
  if (relaxationSteps == 0) return false;
  for (const PolymerBond &bond : bonds)
  {
    if (currentStep < bond.formedAtStep + relaxationSteps) return true;
  }
  for (const BrokenBond &bond : brokenBonds)
  {
    if (currentStep < bond.brokenAtStep + relaxationSteps) return true;
  }
  for (std::size_t i = 0; i < groupState.size(); ++i)
  {
    if (reacted(i) && currentStep < groupStep[i] + relaxationSteps) return true;
  }
  return false;
}

// ---- reactions ----------------------------------------------------------------------------------------------------

void PolymerLinkTable::retype(std::size_t globalAtom, std::uint16_t type, std::uint64_t step,
                              const ForceField &forceField, std::span<Atom> atoms)
{
  Atom &atom = atoms[globalAtom];
  if (atom.type == type) return;
  retypings.push_back(Retyping{layout.ref(static_cast<std::uint32_t>(globalAtom)), atom.type, type, step});
  atom.type = type;
  atom.charge = forceField.pseudoAtoms[static_cast<std::size_t>(type)].charge;
}

void PolymerLinkTable::updateMoleculeMass(const PolymerAtomRef &ref, const ForceField &forceField,
                                          std::span<const Atom> atoms, std::span<Molecule> moleculeData) const
{
  const std::size_t start = layout.moleculeStart(ref.componentId, ref.moleculeIndex);
  const std::size_t count = layout.atomsPerMolecule[ref.componentId];
  double mass = 0.0;
  for (std::size_t i = start; i < start + count; ++i)
  {
    mass += forceField.pseudoAtoms[static_cast<std::size_t>(atoms[i].type)].mass;
  }
  Molecule &molecule = moleculeData[layout.moleculeOffset[ref.componentId] + ref.moleculeIndex];
  molecule.mass = mass;
  molecule.invMass = mass > 0.0 ? 1.0 / mass : 0.0;
}

std::vector<std::size_t> PolymerLinkTable::apply(std::size_t ruleIndex, std::span<const std::size_t> groups,
                                                 std::span<const std::int8_t> chiralities,
                                                 std::span<const std::int8_t> nextChiralities,
                                                 const ForceField &forceField, const std::vector<Component> &components,
                                                 std::span<Atom> atoms, std::span<Molecule> moleculeData,
                                                 const SimulationBox &simulationBox, std::uint64_t step)
{
  const ReactionRule &rule = scheme.rules[ruleIndex];
  if (groups.size() != rule.reactants.size())
  {
    throw std::runtime_error(std::format("[Polymerization]: reaction '{}' takes {} group(s), {} given\n", rule.name,
                                         rule.reactants.size(), groups.size()));
  }
  if (chiralities.size() != rule.stereoAtoms.size() || nextChiralities.size() != rule.stereoAtoms.size())
  {
    throw std::runtime_error(std::format("[Polymerization]: reaction '{}' declares {} stereocentre(s), {} configurations given\n",
                                         rule.name, rule.stereoAtoms.size(), chiralities.size()));
  }
  for (std::size_t i = 0; i < groups.size(); ++i)
  {
    if (groupType[groups[i]] != rule.reactantType[i] || groupState[groups[i]] != rule.reactantState[i])
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: reaction '{}': reactant {} is a '{}' group in state '{}', not '{}' in state '{}'\n",
          rule.name, rule.reactants[i].label, scheme.types[groupType[groups[i]]].name, stateName(groups[i]),
          rule.reactants[i].type, rule.reactants[i].state));
    }
  }
  if (groups.size() == 2 && groups[0] == groups[1])
  {
    throw std::runtime_error(std::format("[Polymerization]: reaction '{}' needs two different groups\n", rule.name));
  }

  // the atoms of the bonds and the stereocentres are resolved before the retyping (the roles do not change, but
  // the bookkeeping below is clearer with everything in hand)
  std::vector<std::array<std::uint32_t, 2>> bondAtoms{};
  for (const std::array<RoleIndex, 2> &bond : rule.bondAtoms)
  {
    bondAtoms.push_back({atomOf(groups[bond[0].reactant], bond[0].role, components),
                         atomOf(groups[bond[1].reactant], bond[1].role, components)});
  }
  std::vector<std::array<std::uint32_t, 4>> centreAtoms{};
  for (std::size_t k = 0; k < rule.stereoAtoms.size(); ++k) centreAtoms.push_back(stereoAtoms(rule, k, groups, components));

  // the bonds broken: bonds of the molecule topology that are still intact
  std::vector<std::array<std::uint32_t, 2>> breakPairs{};
  for (const std::array<RoleIndex, 2> &bond : rule.breakAtoms)
  {
    const std::array<std::uint32_t, 2> pair{atomOf(groups[bond[0].reactant], bond[0].role, components),
                                            atomOf(groups[bond[1].reactant], bond[1].role, components)};
    const PolymerAtomRef ra = layout.ref(pair[0]);
    const PolymerAtomRef rb = layout.ref(pair[1]);
    const bool alreadyBroken = std::ranges::any_of(
        brokenBonds, [&](const BrokenBond &b) { return (b.a == ra && b.b == rb) || (b.a == rb && b.b == ra); });
    if (alreadyBroken)
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: reaction '{}' breaks the bond {}-{} of molecule {} of component '{}', which an earlier "
          "reaction already broke\n",
          rule.name, ra.atomIndex, rb.atomIndex, ra.moleculeIndex, components[ra.componentId].name));
    }
    breakPairs.push_back(pair);
  }

  // the products: retyped atoms and new states
  for (std::size_t i = 0; i < groups.size(); ++i)
  {
    for (const auto &[role, type] : rule.retype[i])
    {
      retype(atomOf(groups[i], role, components), type, step, forceField, atoms);
    }
    groupState[groups[i]] = rule.productState[i];
    groupStep[groups[i]] = step;
  }
  for (std::size_t i = 0; i < groups.size(); ++i)
  {
    const PolymerAtomRef r = groupRef(groups[i]);
    if (i > 0 && r.sameMolecule(groupRef(groups[0]))) continue;
    updateMoleculeMass(r, forceField, atoms, moleculeData);
  }

  // the bonds
  std::vector<std::size_t> newBonds{};
  for (const std::array<std::uint32_t, 2> &pair : bondAtoms)
  {
    const double3 dr = simulationBox.applyPeriodicBoundaryConditions(atoms[pair[0]].position - atoms[pair[1]].position);
    bonds.push_back(PolymerBond{layout.ref(pair[0]), layout.ref(pair[1]), static_cast<std::uint16_t>(ruleIndex), step,
                                std::sqrt(double3::dot(dr, dr))});
    newBonds.push_back(bonds.size() - 1);
  }
  for (const std::array<std::uint32_t, 2> &pair : breakPairs)
  {
    brokenBonds.push_back(BrokenBond{layout.ref(pair[0]), layout.ref(pair[1]), static_cast<std::uint16_t>(ruleIndex), step});
  }

  // the stereocentres: the ramp of the restraint follows the (first) bond that involves the centre atom
  for (std::size_t k = 0; k < rule.stereoAtoms.size(); ++k)
  {
    StereoCentre centre{};
    centre.atoms = centreAtoms[k];
    centre.chirality = chiralities[k];
    centre.group = groups[rule.stereoAtoms[k].centre.reactant];
    centre.bond = static_cast<std::uint32_t>(newBonds.front());
    for (std::size_t b = 0; b < bondAtoms.size(); ++b)
    {
      if (bondAtoms[b][0] == centre.atoms[0] || bondAtoms[b][1] == centre.atoms[0])
      {
        centre.bond = static_cast<std::uint32_t>(newBonds[b]);
        break;
      }
    }
    centres.push_back(centre);
    groupChirality[centre.group] = chiralities[k];
    if (rule.stereoAtoms[k].next.has_value())
    {
      const std::size_t next = groups[rule.stereoAtoms[k].next.value()];
      groupChirality[next] = nextChiralities[k];
      previousGroup[next] = centre.group;
    }
  }

  ++reactionCounts[ruleIndex];
  derivedValid = false;
  return newBonds;
}

// ---- derived topology ---------------------------------------------------------------------------------------------

std::vector<std::uint32_t> PolymerLinkTable::neighboursOf(std::uint32_t atom,
                                                          const std::vector<Component> &components) const
{
  const PolymerAtomRef r = layout.ref(atom);
  const Component &component = components[r.componentId];
  const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
  std::vector<std::uint32_t> result{};
  if (r.atomIndex < component.connectivityTable.numberOfBeads)
  {
    for (std::size_t local : component.connectivityTable.findAllNeighbors(r.atomIndex))
    {
      const std::uint32_t partner = static_cast<std::uint32_t>(start + local);
      if (!brokenPairs.empty() && brokenPairs.contains(canonical<2>({atom, partner}))) continue;
      result.push_back(partner);
    }
  }
  const auto it = linkNeighbours.find(atom);
  if (it != linkNeighbours.end())
  {
    for (std::uint32_t partner : it->second)
    {
      if (!std::ranges::contains(result, partner)) result.push_back(partner);
    }
  }
  return result;
}

void PolymerLinkTable::rebuild(const ForceField &forceField, const std::vector<Component> &components,
                               const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                               std::span<const Atom> atoms)
{
  layout = makeLayout(components, numberOfMoleculesPerComponent);
  if (!templates.resolved) templates.resolve(forceField);

  linkBonds.clear();
  linkBends.clear();
  linkTorsions.clear();
  pairCorrections.clear();
  overrideBonds.clear();
  overrideBends.clear();
  overrideTorsions.clear();
  stereoCentres.clear();
  linkNeighbours.clear();
  brokenPairs.clear();
  nextGroup.assign(groupState.size(), noGroup);
  for (std::size_t g = 0; g < groupState.size(); ++g)
  {
    if (previousGroup[g] != noGroup) nextGroup[previousGroup[g]] = g;
  }

  for (std::uint32_t k = 0; k < bonds.size(); ++k)
  {
    const std::uint32_t a = static_cast<std::uint32_t>(layout.index(bonds[k].a));
    const std::uint32_t b = static_cast<std::uint32_t>(layout.index(bonds[k].b));
    linkNeighbours[a].push_back(b);
    linkNeighbours[b].push_back(a);
  }
  for (const BrokenBond &bond : brokenBonds)
  {
    const std::uint32_t a = static_cast<std::uint32_t>(layout.index(bond.a));
    const std::uint32_t b = static_cast<std::uint32_t>(layout.index(bond.b));
    auto [it, inserted] = brokenPairs.try_emplace(canonical<2>({a, b}), bond.brokenAtStep);
    if (!inserted) it->second = std::max(it->second, bond.brokenAtStep);
  }

  // the restraints of the stereocentres with a held configuration
  for (const StereoCentre &centre : centres)
  {
    if (centre.chirality == 0) continue;
    stereoCentres.push_back(StereoCentreTerm{centre.atoms, centre.chirality, centre.bond});
  }

  auto typeName = [&](std::uint32_t atom) -> const std::string &
  { return forceField.pseudoAtoms[static_cast<std::size_t>(atoms[atom].type)].name; };
  auto describe = [&](std::span<const std::uint32_t> term)
  {
    std::string text{};
    for (std::uint32_t atom : term)
    {
      const PolymerAtomRef r = layout.ref(atom);
      text += std::format("{}[{}:{}:{}] ", typeName(atom), components[r.componentId].name, r.moleculeIndex, r.atomIndex);
    }
    return text;
  };
  auto newer = [&](std::uint32_t k, std::uint32_t other) { return bonds[k].formedAtStep > bonds[other].formedAtStep; };

  // The bonded terms that contain a polymer bond, each once (a bend or torsion with two polymer bonds is reached
  // from both; the newest bond sets its ramp).
  std::map<std::array<std::uint32_t, 2>, std::uint32_t> bondTerms{};
  std::map<std::array<std::uint32_t, 3>, std::uint32_t> bendTerms{};
  std::map<std::array<std::uint32_t, 4>, std::uint32_t> torsionTerms{};
  auto keep = [&](auto &map, const auto &key, std::uint32_t k)
  {
    auto [it, inserted] = map.try_emplace(key, k);
    if (!inserted && newer(k, it->second)) it->second = k;
  };
  auto without = [](std::vector<std::uint32_t> list, std::uint32_t exclude)
  {
    std::erase(list, exclude);
    return list;
  };

  for (std::uint32_t k = 0; k < bonds.size(); ++k)
  {
    const std::uint32_t A = static_cast<std::uint32_t>(layout.index(bonds[k].a));
    const std::uint32_t B = static_cast<std::uint32_t>(layout.index(bonds[k].b));
    keep(bondTerms, canonical<2>({A, B}), k);
    const std::vector<std::uint32_t> nA = without(neighboursOf(A, components), B);
    const std::vector<std::uint32_t> nB = without(neighboursOf(B, components), A);
    for (std::uint32_t X : nA) keep(bendTerms, canonical<3>({X, A, B}), k);
    for (std::uint32_t Y : nB) keep(bendTerms, canonical<3>({A, B, Y}), k);
    for (std::uint32_t X : nA)
    {
      for (std::uint32_t Y : nB)
      {
        if (X != Y) keep(torsionTerms, canonical<4>({X, A, B, Y}), k);
      }
      for (std::uint32_t W : neighboursOf(X, components))
      {
        if (W != A && W != B) keep(torsionTerms, canonical<4>({W, X, A, B}), k);
      }
    }
    for (std::uint32_t Y : nB)
    {
      for (std::uint32_t Z : neighboursOf(Y, components))
      {
        if (Z != B && Z != A) keep(torsionTerms, canonical<4>({A, B, Y, Z}), k);
      }
    }
  }

  for (const auto &[key, k] : bondTerms)
  {
    const std::optional<BondPotential> potential =
        templates.bond({atoms[key[0]].type, atoms[key[1]].type}, {key[0], key[1]});
    if (!potential.has_value())
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: no 'Bonds' template for the polymer bond {}(add [[\"{}\", \"{}\"], ...] to the "
          "'Polymerization' block)\n",
          describe(key), typeName(key[0]), typeName(key[1])));
    }
    linkBonds.push_back(LinkBondTerm{key, potential.value(), k});
  }
  for (const auto &[key, k] : bendTerms)
  {
    const std::optional<BendPotential> potential = templates.bend(
        {atoms[key[0]].type, atoms[key[1]].type, atoms[key[2]].type}, {key[0], key[1], key[2]});
    if (!potential.has_value())
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: no 'Bends' template for the bend {}(add [[\"{}\", \"{}\", \"{}\"], ...] to the "
          "'Polymerization' block)\n",
          describe(key), typeName(key[0]), typeName(key[1]), typeName(key[2])));
    }
    linkBends.push_back(LinkBendTerm{key, potential.value(), k});
  }
  for (const auto &[key, k] : torsionTerms)
  {
    const std::optional<TorsionPotential> potential =
        templates.torsion({atoms[key[0]].type, atoms[key[1]].type, atoms[key[2]].type, atoms[key[3]].type},
                          {key[0], key[1], key[2], key[3]});
    if (!potential.has_value())
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: no 'Torsions' template for the torsion {}(add [[\"{}\", \"{}\", \"{}\", \"{}\"], ...] "
          "to the 'Polymerization' block)\n",
          describe(key), typeName(key[0]), typeName(key[1]), typeName(key[2]), typeName(key[3])));
    }
    linkTorsions.push_back(LinkTorsionTerm{key, potential.value(), k});
  }

  // The pairs whose non-bonded treatment changes: within three bonds of each other through a formed bond, or
  // through a broken bond. The events (formations and breakings) are replayed in the order of their steps on a
  // copy of the topology that grows with them, so that every pair is assigned the treatment the topology of the
  // moment requires and the ramp of the event that changed it. A pair changed by several events (the bond's own
  // pair is within three bonds of the next bond of its atom) keeps the entries of all of them, a chain: the entry
  // of an event takes the pair from the treatment the previous entry established ('counted') to the new one
  // ('wanted'), lambda (wanted - counted); an entry whose treatment the newer event leaves as it is is not
  // repeated (its ramp would restart, and with it the full non-bonded repulsion of two atoms already pulled to
  // the bond length).
  struct TopologyEvent
  {
    std::uint64_t step;
    bool breaking;
    std::array<std::uint32_t, 2> atoms;
  };
  std::vector<TopologyEvent> events{};
  events.reserve(bonds.size() + brokenBonds.size());
  for (const PolymerBond &bond : bonds)
  {
    events.push_back(TopologyEvent{bond.formedAtStep, false,
                                   {static_cast<std::uint32_t>(layout.index(bond.a)),
                                    static_cast<std::uint32_t>(layout.index(bond.b))}});
  }
  for (const BrokenBond &bond : brokenBonds)
  {
    events.push_back(TopologyEvent{bond.brokenAtStep, true,
                                   {static_cast<std::uint32_t>(layout.index(bond.a)),
                                    static_cast<std::uint32_t>(layout.index(bond.b))}});
  }
  std::ranges::stable_sort(events, {}, &TopologyEvent::step);

  // the topology as of the event being replayed
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> formedSoFar{};
  std::set<std::array<std::uint32_t, 2>> brokenSoFar{};
  auto neighboursAt = [&](std::uint32_t atom)
  {
    const PolymerAtomRef r = layout.ref(atom);
    const Component &component = components[r.componentId];
    const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
    std::vector<std::uint32_t> result{};
    if (r.atomIndex < component.connectivityTable.numberOfBeads)
    {
      for (std::size_t local : component.connectivityTable.findAllNeighbors(r.atomIndex))
      {
        const std::uint32_t partner = static_cast<std::uint32_t>(start + local);
        if (!brokenSoFar.empty() && brokenSoFar.contains(canonical<2>({atom, partner}))) continue;
        result.push_back(partner);
      }
    }
    const auto it = formedSoFar.find(atom);
    if (it != formedSoFar.end())
    {
      for (std::uint32_t partner : it->second)
      {
        if (!std::ranges::contains(result, partner)) result.push_back(partner);
      }
    }
    return result;
  };
  // the atoms within 'maxDepth' bonds of 'start' (with their distance); 'extra' is a bond added to the topology
  // (a broken bond, for the neighbourhood through it)
  auto breadthFirst = [&](std::uint32_t start, std::size_t maxDepth,
                          std::optional<std::array<std::uint32_t, 2>> extra = std::nullopt)
  {
    std::unordered_map<std::uint32_t, std::size_t> distance{{start, 0}};
    std::vector<std::uint32_t> frontier{start};
    for (std::size_t depth = 1; depth <= maxDepth && !frontier.empty(); ++depth)
    {
      std::vector<std::uint32_t> next{};
      for (std::uint32_t atom : frontier)
      {
        std::vector<std::uint32_t> neighbours = neighboursAt(atom);
        if (extra.has_value())
        {
          if (atom == (*extra)[0]) neighbours.push_back((*extra)[1]);
          if (atom == (*extra)[1]) neighbours.push_back((*extra)[0]);
        }
        for (std::uint32_t neighbour : neighbours)
        {
          if (distance.try_emplace(neighbour, depth).second) next.push_back(neighbour);
        }
      }
      frontier = std::move(next);
    }
    return distance;
  };
  std::unordered_map<std::uint32_t, std::unordered_map<std::uint32_t, std::size_t>> distanceCache{};
  auto graphDistance = [&](std::uint32_t u, std::uint32_t v) -> std::size_t
  {
    auto it = distanceCache.find(u);
    if (it == distanceCache.end()) it = distanceCache.emplace(u, breadthFirst(u, 3)).first;
    const auto found = it->second.find(v);
    return found == it->second.end() ? 4 : found->second;
  };
  const bool pair14Parameters = forceField.hasPair14Parameters();
  struct PairTreatment
  {
    PairState state{PairState::Full};
    PairScalingInfo scaling{};
    bool operator==(const PairTreatment &other) const
    {
      if (state != other.state) return false;
      if (state != PairState::Scaled) return true;
      return scaling.scalingVDW == other.scaling.scalingVDW && scaling.scalingCoulomb == other.scaling.scalingCoulomb &&
             scaling.pair14 == other.scaling.pair14;
    }
  };
  // what the pair sums count for the pair (the exclusions of the molecule as the component defines it)
  auto countedBySums = [&](const std::array<std::uint32_t, 2> &key)
  {
    PairTreatment counted{};
    const PolymerAtomRef ru = layout.ref(key[0]);
    const PolymerAtomRef rv = layout.ref(key[1]);
    if (!ru.sameMolecule(rv)) return counted;
    const IntraMolecularExclusions &exclusions = components[ru.componentId].intraMolecularPotentials.exclusions;
    if (exclusions.isExcluded(ru.atomIndex, rv.atomIndex))
    {
      counted.state = PairState::Excluded;
      return counted;
    }
    const IntraMolecularExclusions::ScaledPair scaled = exclusions.scalingOf(ru.atomIndex, rv.atomIndex);
    if (scaled.scalingVDW != 1.0 || scaled.scalingCoulomb != 1.0 || scaled.pair14)
    {
      counted.state = PairState::Scaled;
      counted.scaling = {scaled.scalingVDW, scaled.scalingCoulomb, scaled.pair14};
    }
    return counted;
  };
  // what the polymer topology requires of a pair at graph distance d
  auto requiredAt = [&](std::size_t d, const Component &component)
  {
    PairTreatment wanted{};
    if (d <= 2)
    {
      wanted.state = PairState::Excluded;
    }
    else if (d == 3)
    {
      wanted.state = PairState::Scaled;
      wanted.scaling = {component.intra14VanDerWaalsScaling, component.intra14ChargeChargeScaling, pair14Parameters};
      if (wanted.scaling.scalingVDW == 1.0 && wanted.scaling.scalingCoulomb == 1.0 && !wanted.scaling.pair14)
      {
        wanted.state = PairState::Full;
        wanted.scaling = {};
      }
    }
    return wanted;
  };
  struct PairChain
  {
    PairTreatment counted{};        ///< What the pair sums count.
    std::vector<std::size_t> entries{};  ///< Indices into 'chain' in the order of the events.
  };
  std::vector<PairCorrection> chain{};
  std::map<std::array<std::uint32_t, 2>, PairChain> chains{};
  // The pairs (u, v) within three bonds of each other through the bond A-B ('fromA', 'fromB': the neighbourhoods
  // of its two atoms in the graph the bond belongs to); 'step' is the step of the formation or breaking of A-B,
  // 'breaking' whether A-B is a broken bond (whose 1-2 and 1-3 pairs the new topology may release).
  auto correctPairs = [&](const std::unordered_map<std::uint32_t, std::size_t> &fromA,
                          const std::unordered_map<std::uint32_t, std::size_t> &fromB, std::uint64_t step,
                          bool breaking)
  {
    for (const auto &[u, du] : fromA)
    {
      for (const auto &[v, dv] : fromB)
      {
        if (u == v || du + 1 + dv > 3) continue;
        const std::array<std::uint32_t, 2> key = canonical<2>({u, v});
        const std::size_t d = graphDistance(key[0], key[1]);  // in the topology of the event (4: beyond three bonds)

        auto chainIt = chains.find(key);
        const PairTreatment counted = chainIt == chains.end() ? countedBySums(key) : chainIt->second.counted;
        const PairTreatment current =
            chainIt == chains.end() ? counted
                                    : PairTreatment{chain[chainIt->second.entries.back()].wanted,
                                                    chain[chainIt->second.entries.back()].wantedScaling};
        // an excluded pair of the molecule is left alone (a formed bond can only bring atoms closer; an exclusion
        // for another reason, a rigid fragment, stands), unless it was a 1-2 or 1-3 pair through a broken bond
        // that the new topology no longer excludes: that pair gets the required treatment below
        if (counted.state == PairState::Excluded && current.state == PairState::Excluded &&
            (!breaking || du + 1 + dv > 2 || d <= 2))
          continue;

        const PairTreatment wanted = requiredAt(d, components[layout.ref(key[0]).componentId]);
        if (wanted == current) continue;

        PairCorrection correction{};
        correction.atoms = key;
        correction.step = step;
        correction.counted = current.state;
        correction.countedScaling = current.scaling;
        correction.wanted = wanted.state;
        correction.wantedScaling = wanted.scaling;
        correction.listed = true;
        if (chainIt == chains.end()) chainIt = chains.emplace(key, PairChain{counted, {}}).first;
        chainIt->second.entries.push_back(chain.size());
        chain.push_back(correction);
      }
    }
  };
  for (std::size_t first = 0; first < events.size();)
  {
    std::size_t last = first;
    while (last < events.size() && events[last].step == events[first].step) ++last;
    // the topology after the events of this step
    for (std::size_t e = first; e < last; ++e)
    {
      const TopologyEvent &event = events[e];
      if (event.breaking)
      {
        brokenSoFar.insert(canonical<2>(event.atoms));
      }
      else
      {
        formedSoFar[event.atoms[0]].push_back(event.atoms[1]);
        formedSoFar[event.atoms[1]].push_back(event.atoms[0]);
      }
    }
    distanceCache.clear();
    for (std::size_t e = first; e < last; ++e)
    {
      const TopologyEvent &event = events[e];
      if (event.breaking)
      {
        // around a broken bond the affected pairs are the ones that were within three bonds through it: the
        // neighbourhoods in the topology with the bond still present
        const std::array<std::uint32_t, 2> bond = canonical<2>(event.atoms);
        correctPairs(breadthFirst(bond[0], 2, bond), breadthFirst(bond[1], 2, bond), event.step, true);
      }
      else
      {
        correctPairs(breadthFirst(event.atoms[0], 2), breadthFirst(event.atoms[1], 2), event.step, false);
      }
    }
    first = last;
  }
  // The pairs the pair lists hold at full strength that the topology excludes at some point of their chain leave
  // the pair lists (see PairCorrection): the first entry of the chain then supplies the whole
  // lambda value(wanted) + (1 - lambda) value(counted) in place of the lists, the later entries stay differences.
  pairCorrections.reserve(chain.size());
  std::vector<std::array<std::uint32_t, 2>> unlistedPairs{};
  for (auto &[key, pairChain] : chains)
  {
    const bool unlistedPair =
        pairChain.counted.state == PairState::Full &&
        std::ranges::any_of(pairChain.entries,
                            [&](std::size_t index) { return chain[index].wanted == PairState::Excluded; });
    if (unlistedPair)
    {
      chain[pairChain.entries.front()].listed = false;
      unlistedPairs.push_back(key);
    }
    for (std::size_t index : pairChain.entries) pairCorrections.push_back(chain[index]);
  }
  unlisted = PairListExclusions::fromPairs(layout.numberOfAtoms, unlistedPairs);

  // The bonded terms of the molecules that contain a retyped atom or run through a broken bond: overridden. A
  // term is taken through its history: from the parameters of the component to the template of the types after
  // each retyping that changes one of its atoms (the ramp of that reaction), and to zero at the breaking of a
  // bond it runs through (the ramp of the breaking; nothing after that). The overlay adds 'replacement - original'
  // per link of the chain.
  struct MoleculeHistory
  {
    std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint64_t, std::uint16_t>>> types{};  // local atom -> (step, type after)
    std::vector<std::pair<std::array<std::uint32_t, 2>, std::uint64_t>> broken{};                     // local pair, step
  };
  std::map<std::pair<std::uint32_t, std::uint32_t>, MoleculeHistory> histories{};
  for (const Retyping &retyping : retypings)
  {
    histories[{retyping.atom.componentId, retyping.atom.moleculeIndex}].types[retyping.atom.atomIndex].push_back(
        {retyping.step, retyping.to});
  }
  for (const auto &[pair, step] : brokenPairs)
  {
    const PolymerAtomRef ra = layout.ref(pair[0]);
    const PolymerAtomRef rb = layout.ref(pair[1]);
    histories[{ra.componentId, ra.moleculeIndex}].broken.push_back(
        {std::array<std::uint32_t, 2>{ra.atomIndex, rb.atomIndex}, step});
  }
  for (auto &[moleculeKey, history] : histories)
  {
    const Component &component = components[moleculeKey.first];
    const std::size_t start = layout.moleculeStart(moleculeKey.first, moleculeKey.second);
    for (auto &[local, changes] : history.types) std::ranges::sort(changes);

    // the type of a local atom after the reactions up to 'step'
    auto typeAt = [&](std::size_t local, std::uint64_t step) -> std::uint16_t
    {
      std::uint16_t type = component.atoms[local].type;
      const auto it = history.types.find(static_cast<std::uint32_t>(local));
      if (it == history.types.end()) return type;
      for (const auto &[changeStep, to] : it->second)
      {
        if (changeStep > step) break;
        type = to;
      }
      return type;
    };
    // the steps at which the types of a term change, and the step at which a bond it runs through breaks
    struct TermHistory
    {
      std::vector<std::uint64_t> retypeSteps{};
      std::optional<std::uint64_t> breakStep{};
    };
    auto termHistory = [&](std::span<const std::size_t> identifiers, bool adjacentOnly)
    {
      TermHistory result{};
      for (std::size_t id : identifiers)
      {
        const auto it = history.types.find(static_cast<std::uint32_t>(id));
        if (it == history.types.end()) continue;
        for (const auto &[step, to] : it->second) result.retypeSteps.push_back(step);
      }
      std::ranges::sort(result.retypeSteps);
      const auto duplicates = std::ranges::unique(result.retypeSteps);
      result.retypeSteps.erase(duplicates.begin(), duplicates.end());
      for (const auto &[pair, step] : history.broken)
      {
        bool through = false;
        if (adjacentOnly)
        {
          for (std::size_t i = 0; i + 1 < identifiers.size(); ++i)
          {
            if ((identifiers[i] == pair[0] && identifiers[i + 1] == pair[1]) ||
                (identifiers[i] == pair[1] && identifiers[i + 1] == pair[0]))
              through = true;
          }
        }
        else
        {
          through = std::ranges::contains(identifiers, pair[0]) && std::ranges::contains(identifiers, pair[1]);
        }
        if (through) result.breakStep = std::min(result.breakStep.value_or(step), step);
      }
      // the retypings from the breaking on do not matter for a term that is gone
      if (result.breakStep.has_value())
      {
        std::erase_if(result.retypeSteps, [&](std::uint64_t step) { return step >= *result.breakStep; });
      }
      return result;
    };
    auto missing = [&](const char *what, std::span<const std::uint32_t> global, std::uint64_t step)
    {
      std::string types{};
      std::string names{};
      for (std::uint32_t atom : global)
      {
        const std::string &name = forceField.pseudoAtoms[typeAt(layout.ref(atom).atomIndex, step)].name;
        types += std::format("{} ", name);
        names += std::format("{}\"{}\"", names.empty() ? "" : ", ", name);
      }
      return std::runtime_error(std::format(
          "[Polymerization]: no '{}' template for the retyped term {}of types {}at step {} (add [[{}], ...] to the "
          "'Polymerization' block)\n",
          what, describe(global), types, step, names));
    };

    for (const BondPotential &bond : component.intraMolecularPotentials.bonds)
    {
      const TermHistory terms = termHistory(bond.identifiers, false);
      if (terms.retypeSteps.empty() && !terms.breakStep.has_value()) continue;
      const std::array<std::uint32_t, 2> global{static_cast<std::uint32_t>(start + bond.identifiers[0]),
                                                static_cast<std::uint32_t>(start + bond.identifiers[1])};
      BondPotential previous = bond;
      previous.identifiers = {global[0], global[1]};
      for (std::uint64_t step : terms.retypeSteps)
      {
        const std::optional<BondPotential> replacement = templates.bond(
            {typeAt(bond.identifiers[0], step), typeAt(bond.identifiers[1], step)}, {global[0], global[1]});
        if (!replacement.has_value()) throw missing("Bonds", global, step);
        overrideBonds.push_back(OverrideBondTerm{global, previous, replacement.value(), step});
        previous = replacement.value();
      }
      if (terms.breakStep.has_value())
      {
        overrideBonds.push_back(OverrideBondTerm{global, previous, BondPotential({global[0], global[1]}, BondType::None, {}),
                                                 *terms.breakStep});
      }
    }
    for (const BendPotential &bend : component.intraMolecularPotentials.bends)
    {
      const TermHistory terms = termHistory(bend.identifiers, true);
      if (terms.retypeSteps.empty() && !terms.breakStep.has_value()) continue;
      const std::array<std::uint32_t, 3> global{static_cast<std::uint32_t>(start + bend.identifiers[0]),
                                                static_cast<std::uint32_t>(start + bend.identifiers[1]),
                                                static_cast<std::uint32_t>(start + bend.identifiers[2])};
      BendPotential previous = bend;
      previous.identifiers = {global[0], global[1], global[2]};
      for (std::uint64_t step : terms.retypeSteps)
      {
        const std::optional<BendPotential> replacement = templates.bend(
            {typeAt(bend.identifiers[0], step), typeAt(bend.identifiers[1], step), typeAt(bend.identifiers[2], step)},
            {global[0], global[1], global[2]});
        if (!replacement.has_value()) throw missing("Bends", global, step);
        overrideBends.push_back(OverrideBendTerm{global, previous, replacement.value(), step});
        previous = replacement.value();
      }
      if (terms.breakStep.has_value())
      {
        overrideBends.push_back(OverrideBendTerm{
            global, previous, BendPotential({global[0], global[1], global[2]}, BendType::Fixed, {}), *terms.breakStep});
      }
    }
    for (const TorsionPotential &torsion : component.intraMolecularPotentials.torsions)
    {
      const TermHistory terms = termHistory(torsion.identifiers, true);
      if (terms.retypeSteps.empty() && !terms.breakStep.has_value()) continue;
      const std::array<std::uint32_t, 4> global{static_cast<std::uint32_t>(start + torsion.identifiers[0]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[1]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[2]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[3])};
      TorsionPotential previous = torsion;
      previous.identifiers = {global[0], global[1], global[2], global[3]};
      for (std::uint64_t step : terms.retypeSteps)
      {
        const std::optional<TorsionPotential> replacement =
            templates.torsion({typeAt(torsion.identifiers[0], step), typeAt(torsion.identifiers[1], step),
                               typeAt(torsion.identifiers[2], step), typeAt(torsion.identifiers[3], step)},
                              {global[0], global[1], global[2], global[3]});
        if (!replacement.has_value()) throw missing("Torsions", global, step);
        overrideTorsions.push_back(OverrideTorsionTerm{global, previous, replacement.value(), step});
        previous = replacement.value();
      }
      if (terms.breakStep.has_value())
      {
        overrideTorsions.push_back(OverrideTorsionTerm{
            global, previous, TorsionPotential({global[0], global[1], global[2], global[3]}, TorsionType::Fixed, {}),
            *terms.breakStep});
      }
    }
    // the term types the overlay does not carry
    const auto &intra = component.intraMolecularPotentials;
    for (const auto &[pair, step] : history.broken)
    {
      auto hasBoth = [&](std::span<const std::size_t> identifiers)
      { return std::ranges::contains(identifiers, pair[0]) && std::ranges::contains(identifiers, pair[1]); };
      const bool unsupported =
          std::ranges::any_of(intra.ureyBradleys, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.inversionBends, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.outOfPlaneBends, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.improperTorsions, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.bondBonds, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.bondBends, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.bondTorsions, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.bendBends, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.bendTorsions, [&](const auto &t) { return hasBoth(t.identifiers); }) ||
          std::ranges::any_of(intra.cmaps, [&](const auto &t) { return hasBoth(t.identifiers); });
      if (unsupported)
      {
        throw std::runtime_error(std::format(
            "[Polymerization]: the broken bond {}-{} of component '{}' takes part in a Urey-Bradley, inversion, "
            "out-of-plane, improper, cross or CMAP term; only bonds, bends and torsions can be switched off\n",
            pair[0], pair[1], component.name));
      }
    }
  }

  derivedValid = true;
}

// ---- statistics and output ----------------------------------------------------------------------------------------

double PolymerLinkTable::angularStrain(std::uint32_t a, std::uint32_t b, const SimulationBox &simulationBox,
                                        const std::vector<Component> &components, std::span<const Atom> atoms,
                                        const std::map<std::uint32_t, std::uint16_t> &retyped) const
{
  auto typeOf = [&](std::uint32_t atom) -> std::uint16_t
  {
    const auto it = retyped.find(atom);
    return it != retyped.end() ? it->second : atoms[atom].type;
  };
  // positions imaged along the chain of the term (the first atom at its own position)
  auto next = [&](const double3 &from, std::uint32_t to)
  { return from + simulationBox.applyPeriodicBoundaryConditions(atoms[to].position - from); };
  auto without = [](std::vector<std::uint32_t> list, std::uint32_t exclude)
  {
    std::erase(list, exclude);
    return list;
  };

  const double3 pA = atoms[a].position;
  const double3 pB = next(pA, b);
  const std::vector<std::uint32_t> nA = without(neighboursOf(a, components), b);
  const std::vector<std::uint32_t> nB = without(neighboursOf(b, components), a);

  double strain = 0.0;
  auto addBend = [&](std::uint32_t x, std::uint32_t y, std::uint32_t z, const double3 &px, const double3 &py,
                     const double3 &pz)
  {
    const std::optional<double> value = templates.bendStrain({typeOf(x), typeOf(y), typeOf(z)}, px, py, pz);
    if (value.has_value()) strain += value.value();
  };
  auto addTorsion = [&](std::uint32_t w, std::uint32_t x, std::uint32_t y, std::uint32_t z, const double3 &pw,
                        const double3 &px, const double3 &py, const double3 &pz)
  {
    const std::optional<double> value =
        templates.torsionStrain({typeOf(w), typeOf(x), typeOf(y), typeOf(z)}, pw, px, py, pz);
    if (value.has_value()) strain += value.value();
  };

  for (std::uint32_t X : nA)
  {
    const double3 pX = next(pA, X);
    addBend(X, a, b, pX, pA, pB);
    for (std::uint32_t Y : nB)
    {
      if (X == Y) continue;
      addTorsion(X, a, b, Y, pX, pA, pB, next(pB, Y));
    }
    for (std::uint32_t W : neighboursOf(X, components))
    {
      if (W == a || W == b) continue;
      addTorsion(W, X, a, b, next(pX, W), pX, pA, pB);
    }
  }
  for (std::uint32_t Y : nB)
  {
    const double3 pY = next(pB, Y);
    addBend(a, b, Y, pA, pB, pY);
    for (std::uint32_t Z : neighboursOf(Y, components))
    {
      if (Z == b || Z == a) continue;
      addTorsion(a, b, Y, Z, pA, pB, pY, next(pY, Z));
    }
  }
  return strain;
}

PolymerLinkTable::NetworkStatistics PolymerLinkTable::networkStatistics(const SimulationBox &simulationBox,
                                                                         const std::vector<Component> &components,
                                                                         std::span<const Atom> atoms) const
{
  NetworkStatistics statistics{};
  if (bonds.empty()) return statistics;

  // the molecules carrying bonds (keys: component, molecule) are the nodes of the molecule graph
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> index{};
  auto idOf = [&](const PolymerAtomRef &r)
  {
    return index.try_emplace({r.componentId, r.moleculeIndex}, index.size()).first->second;
  };

  // unwrapped positions per molecule (the atoms of a molecule made contiguous by walking the connectivity with
  // minimum-image bond vectors from atom 0), so that the lattice vector of a bond is independent of how the
  // atoms happen to be wrapped
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<double3>> unwrapped{};
  auto unwrappedPosition = [&](const PolymerAtomRef &r) -> double3
  {
    auto [it, inserted] = unwrapped.try_emplace({r.componentId, r.moleculeIndex});
    if (inserted)
    {
      const Component &component = components[r.componentId];
      const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
      const std::size_t n = layout.atomsPerMolecule[r.componentId];
      std::vector<double3> &positions = it->second;
      positions.resize(n);
      std::vector<std::uint8_t> visited(n, 0);
      positions[0] = atoms[start].position;
      visited[0] = 1;
      std::vector<std::size_t> queue{0};
      for (std::size_t head = 0; head < queue.size(); ++head)
      {
        const std::size_t current = queue[head];
        if (current >= component.connectivityTable.numberOfBeads) continue;
        for (std::size_t neighbour : component.connectivityTable.findAllNeighbors(current))
        {
          if (neighbour >= n || visited[neighbour]) continue;
          positions[neighbour] = positions[current] + simulationBox.applyPeriodicBoundaryConditions(
                                                          atoms[start + neighbour].position - atoms[start + current].position);
          visited[neighbour] = 1;
          queue.push_back(neighbour);
        }
      }
      // atoms not reached through the connectivity (disconnected fragments): minimum image relative to atom 0
      for (std::size_t k = 0; k < n; ++k)
      {
        if (!visited[k])
        {
          positions[k] = positions[0] + simulationBox.applyPeriodicBoundaryConditions(atoms[start + k].position -
                                                                                      atoms[start].position);
        }
      }
    }
    return it->second[r.atomIndex];
  };

  using Lattice = std::array<std::int64_t, 3>;
  auto add = [](const Lattice &x, const Lattice &y) { return Lattice{x[0] + y[0], x[1] + y[1], x[2] + y[2]}; };
  auto subtract = [](const Lattice &x, const Lattice &y) { return Lattice{x[0] - y[0], x[1] - y[1], x[2] - y[2]}; };

  // the edges with their lattice vectors: position(b) - position(a) = (unwrapped difference) + cell * shift
  struct Edge
  {
    std::size_t a;
    std::size_t b;
    Lattice shift;
  };
  std::vector<Edge> edges{};
  edges.reserve(bonds.size());
  for (const PolymerBond &bond : bonds)
  {
    const std::size_t a = idOf(bond.a);
    const std::size_t b = idOf(bond.b);
    if (a == b) ++statistics.intramolecularBonds;
    const double3 fractional = simulationBox.inverseCell * (unwrappedPosition(bond.b) - unwrappedPosition(bond.a));
    const Lattice shift{-static_cast<std::int64_t>(std::round(fractional.x)),
                        -static_cast<std::int64_t>(std::round(fractional.y)),
                        -static_cast<std::int64_t>(std::round(fractional.z))};
    edges.push_back(Edge{a, b, shift});
  }

  // union-find with lattice offsets: offset[x] = position(x) - position(parent[x]) in lattice vectors
  std::vector<std::size_t> parent(index.size());
  std::iota(parent.begin(), parent.end(), 0uz);
  std::vector<Lattice> offset(index.size(), Lattice{0, 0, 0});
  auto find = [&](std::size_t x) -> std::pair<std::size_t, Lattice>
  {
    std::vector<std::size_t> path{};
    while (parent[x] != x)
    {
      path.push_back(x);
      x = parent[x];
    }
    const std::size_t root = x;
    Lattice accumulated{0, 0, 0};
    for (auto it = path.rbegin(); it != path.rend(); ++it)
    {
      accumulated = add(accumulated, offset[*it]);
      offset[*it] = accumulated;
      parent[*it] = root;
    }
    return {root, path.empty() ? Lattice{0, 0, 0} : offset[path.front()]};
  };
  for (const Edge &edge : edges)
  {
    const auto [rootA, offsetA] = find(edge.a);
    const auto [rootB, offsetB] = find(edge.b);
    if (rootA != rootB)
    {
      parent[rootB] = rootA;
      offset[rootB] = subtract(add(edge.shift, offsetA), offsetB);
    }
  }

  // the winding vectors of the cycles: an edge inside a cluster whose lattice vector differs from the one implied
  // by the spanning tree closes a cycle with winding = shift - (offset(b) - offset(a))
  std::map<std::size_t, std::vector<Lattice>> windings{};
  for (const Edge &edge : edges)
  {
    const auto [rootA, offsetA] = find(edge.a);
    const auto [rootB, offsetB] = find(edge.b);
    const Lattice winding = subtract(edge.shift, subtract(offsetB, offsetA));
    if (winding != Lattice{0, 0, 0})
    {
      windings[rootA].push_back(winding);
      for (std::size_t d = 0; d < 3; ++d) statistics.percolates[d] = statistics.percolates[d] || (winding[d] != 0);
    }
  }
  // the rank of a set of integer vectors (Gaussian elimination; the vectors are small integers)
  auto rankOf = [](const std::vector<Lattice> &vectors) -> std::size_t
  {
    std::vector<std::array<double, 3>> rows{};
    for (const Lattice &v : vectors)
    {
      rows.push_back({static_cast<double>(v[0]), static_cast<double>(v[1]), static_cast<double>(v[2])});
    }
    std::size_t rank = 0;
    for (std::size_t column = 0; column < 3 && rank < rows.size(); ++column)
    {
      std::size_t pivot = rank;
      for (std::size_t r = rank; r < rows.size(); ++r)
      {
        if (std::abs(rows[r][column]) > std::abs(rows[pivot][column])) pivot = r;
      }
      if (std::abs(rows[pivot][column]) < 1e-9) continue;
      std::swap(rows[pivot], rows[rank]);
      for (std::size_t r = 0; r < rows.size(); ++r)
      {
        if (r == rank) continue;
        const double factor = rows[r][column] / rows[rank][column];
        for (std::size_t c = 0; c < 3; ++c) rows[r][c] -= factor * rows[rank][c];
      }
      ++rank;
    }
    return rank;
  };

  std::map<std::size_t, std::size_t> clusterSize{};
  for (std::size_t i = 0; i < parent.size(); ++i) ++clusterSize[find(i).first];
  statistics.moleculesInNetwork = index.size();
  statistics.numberOfClusters = clusterSize.size();
  std::size_t &spanningMolecules = statistics.spanningMolecules;
  for (const auto &[root, size] : clusterSize)
  {
    statistics.largestCluster = std::max(statistics.largestCluster, size);
    const auto it = windings.find(root);
    if (it == windings.end()) continue;
    const std::size_t dimension = rankOf(it->second);
    ++statistics.spanningClusters;
    spanningMolecules += size;
    if (size > statistics.largestSpanningCluster)
    {
      statistics.largestSpanningCluster = size;
      statistics.percolationDimension = dimension;
    }
  }
  statistics.cycleRank = bonds.size() + statistics.numberOfClusters - statistics.moleculesInNetwork;
  statistics.effectiveCycleRank = statistics.cycleRank - statistics.spanningClusters;

  std::size_t moleculesWithSites = 0;
  for (std::size_t c = 0; c < groupsPerMolecule.size(); ++c)
  {
    if (groupsPerMolecule[c] > 0) moleculesWithSites += moleculesPerComponent[c];
  }
  if (moleculesWithSites > 0)
  {
    statistics.gelFraction = static_cast<double>(statistics.largestCluster) / static_cast<double>(moleculesWithSites);
    statistics.spanningFraction = static_cast<double>(spanningMolecules) / static_cast<double>(moleculesWithSites);
  }
  return statistics;
}

PolymerLinkTable::StereoStatistics PolymerLinkTable::stereoStatistics(const SimulationBox &simulationBox,
                                                                       const std::vector<Component> &,
                                                                       std::span<const Atom> atoms) const
{
  StereoStatistics statistics{};
  if (!derivedValid || nextGroup.size() != groupState.size()) return statistics;

  // the sign of the chiral volume of the (most recent) centre of every group (0: no centre)
  std::vector<std::int8_t> sign(groupState.size(), 0);
  for (const StereoCentre &centre : centres)
  {
    const std::int8_t geometric = geometricChirality(centre.atoms, simulationBox, atoms);
    sign[centre.group] = geometric;
    ++statistics.stereoCentres;
    if (centre.chirality != 0)
    {
      ++statistics.restrained;
      if (centre.chirality != geometric) ++statistics.inverted;
    }
  }
  // the dyads (previous group, this group) and the triads (previous, this, next)
  auto dyadIsMeso = [&](std::size_t i, std::size_t j) -> std::optional<bool>
  {
    if (i == noGroup || j == noGroup || sign[i] == 0 || sign[j] == 0) return std::nullopt;
    return sign[i] == sign[j];
  };
  for (std::size_t g = 0; g < groupState.size(); ++g)
  {
    const std::optional<bool> dyad = dyadIsMeso(previousGroup[g], g);
    if (!dyad.has_value()) continue;
    if (dyad.value()) ++statistics.mesoDyads; else ++statistics.racemoDyads;
    const std::optional<bool> next = dyadIsMeso(g, nextGroup[g]);
    if (!next.has_value()) continue;
    if (dyad.value() && next.value())
      ++statistics.mm;
    else if (!dyad.value() && !next.value())
      ++statistics.rr;
    else
      ++statistics.mr;
  }
  return statistics;
}

std::string PolymerLinkTable::printStatus(const SimulationBox &simulationBox, const std::vector<Component> &components,
                                          std::span<const Atom> atoms) const
{
  std::ostringstream stream;
  const NetworkStatistics network = networkStatistics(simulationBox, components, atoms);
  const StereoStatistics stereo = stereoStatistics(simulationBox, components, atoms);
  std::print(stream, "Polymerization status\n");
  std::print(stream, "------------------------------------------------------------------------------------------------------------------------\n");
  std::print(stream, "    groups:                          {} total, {} active (living chain ends)\n", numberOfGroups(),
             countActive());
  for (std::size_t t = 0; t < scheme.types.size(); ++t)
  {
    const ReactionScheme::GroupType &type = scheme.types[t];
    std::string counts{};
    for (std::size_t s = 0; s < type.states.size(); ++s)
    {
      counts += std::format("{} {}{}", countGroups(static_cast<std::uint16_t>(t), static_cast<std::uint8_t>(s)),
                            type.states[s], s + 1 < type.states.size() ? ", " : "");
    }
    std::print(stream, "    '{}' groups:{}{} total: {}\n", type.name,
               std::string(std::max<std::size_t>(1, 22 - std::min<std::size_t>(22, type.name.size())), ' '),
               countGroups(static_cast<std::uint16_t>(t)), counts);
  }
  std::print(stream, "    conversion:                      {:.4f}\n", conversion());
  std::string events{};
  for (std::size_t r = 0; r < scheme.rules.size(); ++r)
  {
    events += std::format("{} {}{}", r < reactionCounts.size() ? reactionCounts[r] : 0, scheme.rules[r].name,
                          r + 1 < scheme.rules.size() ? ", " : "");
  }
  std::print(stream, "    events:                          {}\n", events);
  std::print(stream, "    polymer bonds:                   {} ({} intramolecular / primary cycles)\n", bonds.size(),
             network.intramolecularBonds);
  if (!brokenBonds.empty() || scheme.hasBreaking())
  {
    std::print(stream, "    broken bonds:                    {}\n", brokenBonds.size());
  }
  std::print(stream, "    network:                         {} molecules in {} clusters, largest {} molecules (fraction {:.4f})\n",
             network.moleculesInNetwork, network.numberOfClusters, network.largestCluster, network.gelFraction);
  std::print(stream, "    cycles:                          rank {} of the molecule graph, {} independent loops per box in the unfolded network\n",
             network.cycleRank, network.effectiveCycleRank);
  if (network.spanningClusters > 0)
  {
    std::print(stream, "    percolation:                     {} spanning cluster(s) with {} molecules (fraction {:.4f}), largest spans {} direction(s) [{}{}{}]\n",
               network.spanningClusters, network.spanningMolecules, network.spanningFraction,
               network.percolationDimension, network.percolates[0] ? "a" : "-",
               network.percolates[1] ? "b" : "-", network.percolates[2] ? "c" : "-");
  }
  else
  {
    std::print(stream, "    percolation:                     no cluster connected to its periodic image (sol)\n");
  }
  for (std::size_t c = 0; c < groupsPerMolecule.size() && c < components.size(); ++c)
  {
    if (groupsPerMolecule[c] == 0) continue;
    std::size_t unreacted = 0, activeCount = 0, reactedCount = 0;
    for (std::size_t i = groupOffset[c]; i < groupOffset[c] + groupsPerMolecule[c] * moleculesPerComponent[c]; ++i)
    {
      if (!reacted(i))
        ++unreacted;
      else
        ++reactedCount;
      if (active(i)) ++activeCount;
    }
    std::print(stream, "    component {} ({}):{}{} groups per molecule, {} unreacted, {} reacted ({} active)\n", c,
               components[c].name,
               std::string(std::max<std::size_t>(1, 15 - std::min<std::size_t>(15, components[c].name.size())), ' '),
               groupsPerMolecule[c], unreacted, reactedCount, activeCount);
  }
  if (derivedValid)
  {
    std::print(stream,
               "    overlay terms:                   {} bonds, {} bends, {} torsions, {} pair corrections, {} "
               "retyped-term overrides, {} stereocentre restraints\n",
               linkBonds.size(), linkBends.size(), linkTorsions.size(), pairCorrections.size(),
               overrideBonds.size() + overrideBends.size() + overrideTorsions.size(), stereoCentres.size());
  }
  if (stereo.stereoCentres > 0)
  {
    std::print(stream, "    stereocentres:                   {} ({} restrained, {} inverted)\n", stereo.stereoCentres,
               stereo.restrained, stereo.inverted);
    std::print(stream, "    tacticity (dyads):               {} meso, {} racemo, meso fraction {:.4f}\n",
               stereo.mesoDyads, stereo.racemoDyads, stereo.mesoFraction());
    if (stereo.mm + stereo.mr + stereo.rr > 0)
    {
      const double total = static_cast<double>(stereo.mm + stereo.mr + stereo.rr);
      std::print(stream, "    tacticity (triads):              mm {:.4f}, mr {:.4f}, rr {:.4f}\n",
                 static_cast<double>(stereo.mm) / total, static_cast<double>(stereo.mr) / total,
                 static_cast<double>(stereo.rr) / total);
    }
  }
  std::print(stream, "    ramps in progress:               {}\n", rampsActive() ? "yes" : "no");
  std::print(stream, "\n");
  return stream.str();
}

nlohmann::json PolymerLinkTable::jsonGroupCounts() const
{
  nlohmann::json counts = nlohmann::json::object();
  for (std::size_t t = 0; t < scheme.types.size(); ++t)
  {
    nlohmann::json states = nlohmann::json::object();
    for (std::size_t s = 0; s < scheme.types[t].states.size(); ++s)
    {
      states[scheme.types[t].states[s]] = countGroups(static_cast<std::uint16_t>(t), static_cast<std::uint8_t>(s));
    }
    counts[scheme.types[t].name] = states;
  }
  return counts;
}

nlohmann::json PolymerLinkTable::jsonStatus(const SimulationBox &simulationBox,
                                            const std::vector<Component> &components,
                                            std::span<const Atom> atoms) const
{
  nlohmann::json status{};
  const NetworkStatistics network = networkStatistics(simulationBox, components, atoms);
  status["numberOfGroups"] = numberOfGroups();
  status["activeGroups"] = countActive();
  status["groups"] = jsonGroupCounts();
  status["conversion"] = conversion();
  nlohmann::json events = nlohmann::json::object();
  for (std::size_t r = 0; r < scheme.rules.size(); ++r)
  {
    events[scheme.rules[r].name] = r < reactionCounts.size() ? reactionCounts[r] : 0;
  }
  status["events"] = events;
  status["bonds"] = bonds.size();
  status["brokenBonds"] = brokenBonds.size();
  status["intramolecularBonds"] = network.intramolecularBonds;
  status["moleculesInNetwork"] = network.moleculesInNetwork;
  status["clusters"] = network.numberOfClusters;
  status["largestCluster"] = network.largestCluster;
  status["gelFraction"] = network.gelFraction;
  status["cycleRank"] = network.cycleRank;
  status["effectiveCycleRank"] = network.effectiveCycleRank;
  status["spanningClusters"] = network.spanningClusters;
  status["spanningMolecules"] = network.spanningMolecules;
  status["largestSpanningCluster"] = network.largestSpanningCluster;
  status["spanningFraction"] = network.spanningFraction;
  status["percolationDimension"] = network.percolationDimension;
  status["percolates"] = {network.percolates[0], network.percolates[1], network.percolates[2]};
  const StereoStatistics stereo = stereoStatistics(simulationBox, components, atoms);
  status["stereoCentres"] = stereo.stereoCentres;
  status["restrainedStereoCentres"] = stereo.restrained;
  status["invertedStereoCentres"] = stereo.inverted;
  status["mesoDyads"] = stereo.mesoDyads;
  status["racemoDyads"] = stereo.racemoDyads;
  status["mesoFraction"] = stereo.mesoFraction();
  status["triads"] = {{"mm", stereo.mm}, {"mr", stereo.mr}, {"rr", stereo.rr}};
  return status;
}

nlohmann::json PolymerLinkTable::jsonTopology() const
{
  nlohmann::json topology = nlohmann::json::array();
  for (const PolymerBond &bond : bonds)
  {
    nlohmann::json entry{};
    entry["atoms"] = {{bond.a.componentId, bond.a.moleculeIndex, bond.a.atomIndex},
                      {bond.b.componentId, bond.b.moleculeIndex, bond.b.atomIndex}};
    entry["reaction"] = bond.reaction < scheme.rules.size() ? scheme.rules[bond.reaction].name : std::string("?");
    entry["step"] = bond.formedAtStep;
    entry["formationDistance"] = bond.formationDistance;
    topology.push_back(entry);
  }
  return topology;
}

nlohmann::json PolymerLinkTable::jsonBrokenBonds() const
{
  nlohmann::json topology = nlohmann::json::array();
  for (const BrokenBond &bond : brokenBonds)
  {
    nlohmann::json entry{};
    entry["atoms"] = {{bond.a.componentId, bond.a.moleculeIndex, bond.a.atomIndex},
                      {bond.b.componentId, bond.b.moleculeIndex, bond.b.atomIndex}};
    entry["reaction"] = bond.reaction < scheme.rules.size() ? scheme.rules[bond.reaction].name : std::string("?");
    entry["step"] = bond.brokenAtStep;
    topology.push_back(entry);
  }
  return topology;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerLinkTable &t)
{
  archive << t.versionNumber;
  archive << t.templates;
  archive << t.scheme;
  archive << t.bonds;
  archive << t.brokenBonds;
  archive << t.retypings;
  archive << t.centres;
  archive << t.activated;
  archive << t.groupsPerMolecule;
  archive << t.groupOffset;
  archive << t.moleculesPerComponent;
  archive << t.groupState;
  archive << t.groupStep;
  archive << t.groupChirality;
  archive << t.previousGroup;
  archive << t.reactionCounts;
  archive << t.relaxationSteps;
  archive << t.currentStep;
  archive << t.chiralForceConstant;
  archive << t.chiralReferenceVolume;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerLinkTable &t)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    throw std::runtime_error(std::format("Invalid version reading 'PolymerLinkTable' at line {} in file {}\n",
                                         std::source_location::current().line(),
                                         std::source_location::current().file_name()));
  }
  archive >> t.templates;
  archive >> t.scheme;
  archive >> t.bonds;
  archive >> t.brokenBonds;
  archive >> t.retypings;
  archive >> t.centres;
  archive >> t.activated;
  archive >> t.groupsPerMolecule;
  archive >> t.groupOffset;
  archive >> t.moleculesPerComponent;
  archive >> t.groupState;
  archive >> t.groupStep;
  archive >> t.groupChirality;
  archive >> t.previousGroup;
  archive >> t.reactionCounts;
  archive >> t.relaxationSteps;
  archive >> t.currentStep;
  archive >> t.chiralForceConstant;
  archive >> t.chiralReferenceVolume;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PolymerLinkTable: Invalid magic number {} at the end of the restart data\n",
                                         magicNumber));
  }
#endif
  t.derivedValid = false;
  return archive;
}
