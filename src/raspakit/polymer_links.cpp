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
import polymerization_site;
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
  archive << b.kind;
  archive << b.formedAtStep;
  archive << b.formationDistance;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerBond &b)
{
  archive >> b.a;
  archive >> b.b;
  archive >> b.kind;
  archive >> b.formedAtStep;
  archive >> b.formationDistance;
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
  resolved = true;
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

// ---- PolymerLinkTable: layout and sites ---------------------------------------------------------------------------

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

void PolymerLinkTable::activate(const std::vector<Component> &components,
                                const std::vector<std::size_t> &numberOfMoleculesPerComponent)
{
  if (activated)
  {
    if (moleculesPerComponent != numberOfMoleculesPerComponent)
    {
      throw std::runtime_error(
          "[Polymerization]: the number of molecules changed after the polymerization was activated\n");
    }
    layout = makeLayout(components, numberOfMoleculesPerComponent);
    if (siteChirality.size() != siteState.size()) siteChirality.assign(siteState.size(), 0);
    return;
  }
  sitesPerMolecule.resize(components.size());
  siteOffset.resize(components.size());
  moleculesPerComponent = numberOfMoleculesPerComponent;
  std::size_t offset = 0;
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    sitesPerMolecule[c] = components[c].polymerizationSites.size();
    siteOffset[c] = offset;
    offset += sitesPerMolecule[c] * numberOfMoleculesPerComponent[c];
  }
  siteState.assign(offset, static_cast<std::uint8_t>(PolymerSiteState::Monomer));
  siteStep.assign(offset, 0);
  siteChirality.assign(offset, 0);
  layout = makeLayout(components, numberOfMoleculesPerComponent);
  activated = true;
  derivedValid = false;
}

PolymerAtomRef PolymerLinkTable::siteRef(std::size_t index) const
{
  for (std::size_t c = siteOffset.size(); c-- > 0;)
  {
    if (index >= siteOffset[c] && sitesPerMolecule[c] > 0)
    {
      const std::size_t local = index - siteOffset[c];
      return PolymerAtomRef{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(local / sitesPerMolecule[c]),
                            static_cast<std::uint32_t>(local % sitesPerMolecule[c])};
    }
  }
  throw std::runtime_error("[Polymerization]: site index outside of the site layout\n");
}

PolymerAtomRef PolymerLinkTable::headOf(std::size_t index, const std::vector<Component> &components) const
{
  PolymerAtomRef r = siteRef(index);
  r.atomIndex = static_cast<std::uint32_t>(components[r.componentId].polymerizationSites[r.atomIndex].head);
  return r;
}

PolymerAtomRef PolymerLinkTable::tailOf(std::size_t index, const std::vector<Component> &components) const
{
  PolymerAtomRef r = siteRef(index);
  r.atomIndex = static_cast<std::uint32_t>(components[r.componentId].polymerizationSites[r.atomIndex].tail);
  return r;
}

std::optional<PolymerAtomRef> PolymerLinkTable::sideGroupOf(std::size_t index,
                                                            const std::vector<Component> &components) const
{
  PolymerAtomRef r = siteRef(index);
  const std::optional<std::size_t> sideGroup = components[r.componentId].polymerizationSites[r.atomIndex].sideGroup;
  if (!sideGroup.has_value()) return std::nullopt;
  r.atomIndex = static_cast<std::uint32_t>(sideGroup.value());
  return r;
}

std::optional<std::pair<std::size_t, bool>> PolymerLinkTable::siteOfAtom(std::uint32_t atom,
                                                                         const std::vector<Component> &components) const
{
  const PolymerAtomRef r = layout.ref(atom);
  const std::vector<PolymerizationSite> &sites = components[r.componentId].polymerizationSites;
  for (std::size_t s = 0; s < sites.size(); ++s)
  {
    if (sites[s].tail == r.atomIndex) return std::pair{siteIndex(r.componentId, r.moleculeIndex, s), true};
    if (sites[s].head == r.atomIndex) return std::pair{siteIndex(r.componentId, r.moleculeIndex, s), false};
  }
  return std::nullopt;
}

std::size_t PolymerLinkTable::countSites(PolymerSiteState state) const
{
  return static_cast<std::size_t>(std::ranges::count(siteState, static_cast<std::uint8_t>(state)));
}

std::optional<double> PolymerLinkTable::chiralVolumeOf(std::size_t site, std::uint32_t partner,
                                                       const SimulationBox &simulationBox,
                                                       const std::vector<Component> &components,
                                                       std::span<const Atom> atoms) const
{
  const std::optional<PolymerAtomRef> sideGroup = sideGroupOf(site, components);
  if (!sideGroup.has_value()) return std::nullopt;
  const double3 centre = atoms[layout.index(tailOf(site, components))].position;
  const double3 head = atoms[layout.index(headOf(site, components))].position;
  const double3 side = atoms[layout.index(sideGroup.value())].position;
  const double3 bp = simulationBox.applyPeriodicBoundaryConditions(head - centre);
  const double3 bn = simulationBox.applyPeriodicBoundaryConditions(atoms[partner].position - centre);
  const double3 bz = simulationBox.applyPeriodicBoundaryConditions(side - centre);
  return ChiralVolume::compute(bp, bn, bz).value;
}

std::int8_t PolymerLinkTable::geometricChirality(std::size_t site, std::uint32_t partner,
                                                 const SimulationBox &simulationBox,
                                                 const std::vector<Component> &components,
                                                 std::span<const Atom> atoms) const
{
  const std::optional<double> volume = chiralVolumeOf(site, partner, simulationBox, components, atoms);
  if (!volume.has_value()) return 0;
  return volume.value() < 0.0 ? -1 : 1;
}

std::size_t PolymerLinkTable::predecessorOf(std::size_t site, const std::vector<Component> &components) const
{
  const PolymerAtomRef head = headOf(site, components);
  for (const PolymerBond &bond : bonds)
  {
    if (bond.kind == PolymerBondKind::Propagation && bond.b == head)
    {
      const std::uint32_t tail = static_cast<std::uint32_t>(layout.index(bond.a));
      const std::optional<std::pair<std::size_t, bool>> previous = siteOfAtom(tail, components);
      if (previous.has_value() && previous->second) return previous->first;
    }
  }
  return noSite;
}

double PolymerLinkTable::conversion() const
{
  if (siteState.empty()) return 0.0;
  return static_cast<double>(siteState.size() - countSites(PolymerSiteState::Monomer)) /
         static_cast<double>(siteState.size());
}

bool PolymerLinkTable::rampsActive() const
{
  if (relaxationSteps == 0) return false;
  for (const PolymerBond &bond : bonds)
  {
    if (currentStep < bond.formedAtStep + relaxationSteps) return true;
  }
  for (std::size_t i = 0; i < siteState.size(); ++i)
  {
    if (siteState[i] != static_cast<std::uint8_t>(PolymerSiteState::Monomer) &&
        currentStep < siteStep[i] + relaxationSteps)
      return true;
  }
  return false;
}

// ---- reactions ----------------------------------------------------------------------------------------------------

void PolymerLinkTable::retype(std::size_t globalAtom, const std::string &typeName, const ForceField &forceField,
                              std::span<Atom> atoms)
{
  const std::optional<std::size_t> type = forceField.findPseudoAtom(typeName);
  if (!type.has_value())
  {
    throw std::runtime_error(std::format("[Polymerization]: unknown product pseudo-atom '{}'\n", typeName));
  }
  Atom &atom = atoms[globalAtom];
  atom.type = static_cast<std::uint16_t>(type.value());
  atom.charge = forceField.pseudoAtoms[type.value()].charge;
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

void PolymerLinkTable::initiate(std::size_t site, const ForceField &forceField, const std::vector<Component> &components,
                                std::span<Atom> atoms, std::span<Molecule> moleculeData, std::uint64_t step)
{
  if (state(site) != PolymerSiteState::Monomer)
  {
    throw std::runtime_error("[Polymerization]: only a monomer site can be initiated\n");
  }
  const PolymerAtomRef r = siteRef(site);
  const PolymerizationSite &definition = components[r.componentId].polymerizationSites[r.atomIndex];
  const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
  retype(start + definition.head, definition.initiatedHeadType.value_or(definition.headProductType), forceField, atoms);
  retype(start + definition.tail, definition.tailProductType, forceField, atoms);
  updateMoleculeMass(r, forceField, atoms, moleculeData);
  siteState[site] = static_cast<std::uint8_t>(PolymerSiteState::Radical);
  siteStep[site] = step;
  siteChirality[site] = 0;
  ++numberOfInitiations;
  derivedValid = false;
}

void PolymerLinkTable::propagate(std::size_t radicalSite, std::size_t monomerSite, const ForceField &forceField,
                                 const std::vector<Component> &components, std::span<Atom> atoms,
                                 std::span<Molecule> moleculeData, std::uint64_t step, double distance,
                                 std::int8_t chiralityRadical, std::int8_t chiralityMonomer)
{
  if (state(radicalSite) != PolymerSiteState::Radical || state(monomerSite) != PolymerSiteState::Monomer)
  {
    throw std::runtime_error("[Polymerization]: propagation needs a radical and a monomer site\n");
  }
  const PolymerAtomRef m = siteRef(monomerSite);
  const PolymerizationSite &definition = components[m.componentId].polymerizationSites[m.atomIndex];
  const std::size_t start = layout.moleculeStart(m.componentId, m.moleculeIndex);
  retype(start + definition.head, definition.headProductType, forceField, atoms);
  retype(start + definition.tail, definition.tailProductType, forceField, atoms);
  updateMoleculeMass(m, forceField, atoms, moleculeData);

  bonds.push_back(PolymerBond{tailOf(radicalSite, components), headOf(monomerSite, components),
                              PolymerBondKind::Propagation, step, distance});
  siteState[radicalSite] = static_cast<std::uint8_t>(PolymerSiteState::Dead);
  siteStep[radicalSite] = step;
  siteChirality[radicalSite] = chiralityRadical;
  siteState[monomerSite] = static_cast<std::uint8_t>(PolymerSiteState::Radical);
  siteStep[monomerSite] = step;
  siteChirality[monomerSite] = chiralityMonomer;
  ++numberOfPropagations;
  derivedValid = false;
}

void PolymerLinkTable::terminate(std::size_t siteA, std::size_t siteB, const ForceField &,
                                 const std::vector<Component> &components, std::span<Atom>, std::span<Molecule>,
                                 std::uint64_t step, double distance, std::int8_t chiralityA, std::int8_t chiralityB)
{
  if (state(siteA) != PolymerSiteState::Radical || state(siteB) != PolymerSiteState::Radical || siteA == siteB)
  {
    throw std::runtime_error("[Polymerization]: termination needs two different radical sites\n");
  }
  bonds.push_back(PolymerBond{tailOf(siteA, components), tailOf(siteB, components), PolymerBondKind::Termination, step,
                              distance});
  siteState[siteA] = static_cast<std::uint8_t>(PolymerSiteState::Dead);
  siteStep[siteA] = step;
  siteChirality[siteA] = chiralityA;
  siteState[siteB] = static_cast<std::uint8_t>(PolymerSiteState::Dead);
  siteStep[siteB] = step;
  siteChirality[siteB] = chiralityB;
  ++numberOfTerminations;
  derivedValid = false;
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
      result.push_back(static_cast<std::uint32_t>(start + local));
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
  tailPartner.assign(siteState.size(), noAtom);
  tailBond.assign(siteState.size(), noAtom);
  previousSite.assign(siteState.size(), noSite);
  nextSite.assign(siteState.size(), noSite);

  for (std::uint32_t k = 0; k < bonds.size(); ++k)
  {
    const PolymerBond &bond = bonds[k];
    const std::uint32_t a = static_cast<std::uint32_t>(layout.index(bond.a));
    const std::uint32_t b = static_cast<std::uint32_t>(layout.index(bond.b));
    linkNeighbours[a].push_back(b);
    linkNeighbours[b].push_back(a);

    // the chain bookkeeping of the sites: 'a' is always the tail of a site; 'b' is the head of the next unit
    // (propagation) or the tail of the other chain (termination)
    const std::optional<std::pair<std::size_t, bool>> siteA = siteOfAtom(a, components);
    const std::optional<std::pair<std::size_t, bool>> siteB = siteOfAtom(b, components);
    if (siteA.has_value() && siteA->second)
    {
      tailPartner[siteA->first] = b;
      tailBond[siteA->first] = k;
    }
    if (siteB.has_value())
    {
      if (siteB->second)
      {
        tailPartner[siteB->first] = a;
        tailBond[siteB->first] = k;
      }
      else if (siteA.has_value() && siteA->second)
      {
        previousSite[siteB->first] = siteA->first;
        nextSite[siteA->first] = siteB->first;
      }
    }
  }

  // the restraints of the reacted stereocentres (a dead site whose tail has a partner and a side group)
  for (std::size_t s = 0; s < siteState.size(); ++s)
  {
    if (siteState[s] != static_cast<std::uint8_t>(PolymerSiteState::Dead) || siteChirality[s] == 0) continue;
    if (tailPartner[s] == noAtom) continue;
    const std::optional<PolymerAtomRef> sideGroup = sideGroupOf(s, components);
    if (!sideGroup.has_value()) continue;
    StereoCentreTerm term{};
    term.atoms = {static_cast<std::uint32_t>(layout.index(tailOf(s, components))),
                  static_cast<std::uint32_t>(layout.index(headOf(s, components))), tailPartner[s],
                  static_cast<std::uint32_t>(layout.index(sideGroup.value()))};
    term.chirality = siteChirality[s];
    term.bondId = tailBond[s];
    stereoCentres.push_back(term);
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

  // The pairs whose non-bonded treatment changes: within three bonds of each other through a polymer bond.
  auto breadthFirst = [&](std::uint32_t start, std::size_t maxDepth)
  {
    std::unordered_map<std::uint32_t, std::size_t> distance{{start, 0}};
    std::vector<std::uint32_t> frontier{start};
    for (std::size_t depth = 1; depth <= maxDepth && !frontier.empty(); ++depth)
    {
      std::vector<std::uint32_t> next{};
      for (std::uint32_t atom : frontier)
      {
        for (std::uint32_t neighbour : neighboursOf(atom, components))
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
  std::map<std::array<std::uint32_t, 2>, PairCorrection> corrections{};
  for (std::uint32_t k = 0; k < bonds.size(); ++k)
  {
    const std::uint32_t A = static_cast<std::uint32_t>(layout.index(bonds[k].a));
    const std::uint32_t B = static_cast<std::uint32_t>(layout.index(bonds[k].b));
    const auto fromA = breadthFirst(A, 2);
    const auto fromB = breadthFirst(B, 2);
    for (const auto &[u, du] : fromA)
    {
      for (const auto &[v, dv] : fromB)
      {
        if (u == v || du + 1 + dv > 3) continue;
        const std::array<std::uint32_t, 2> key = canonical<2>({u, v});
        const std::size_t d = graphDistance(key[0], key[1]);
        if (d > 3) continue;

        PairCorrection correction{};
        correction.atoms = key;
        correction.bondId = k;
        const PolymerAtomRef ru = layout.ref(key[0]);
        const PolymerAtomRef rv = layout.ref(key[1]);
        const Component &component = components[ru.componentId];

        // what the pair sums count for the pair now
        correction.counted = PairState::Full;
        if (ru.sameMolecule(rv))
        {
          const IntraMolecularExclusions &exclusions = component.intraMolecularPotentials.exclusions;
          if (exclusions.isExcluded(ru.atomIndex, rv.atomIndex))
          {
            correction.counted = PairState::Excluded;
          }
          else
          {
            const IntraMolecularExclusions::ScaledPair scaled = exclusions.scalingOf(ru.atomIndex, rv.atomIndex);
            if (scaled.scalingVDW != 1.0 || scaled.scalingCoulomb != 1.0 || scaled.pair14)
            {
              correction.counted = PairState::Scaled;
              correction.countedScaling = {scaled.scalingVDW, scaled.scalingCoulomb, scaled.pair14};
            }
          }
        }
        if (correction.counted == PairState::Excluded) continue;

        // what the polymer topology requires
        if (d <= 2)
        {
          correction.wanted = PairState::Excluded;
        }
        else
        {
          correction.wanted = PairState::Scaled;
          correction.wantedScaling = {component.intra14VanDerWaalsScaling, component.intra14ChargeChargeScaling,
                                      pair14Parameters};
          if (correction.wantedScaling.scalingVDW == 1.0 && correction.wantedScaling.scalingCoulomb == 1.0 &&
              !correction.wantedScaling.pair14)
          {
            correction.wanted = PairState::Full;
          }
        }
        if (correction.wanted == correction.counted &&
            (correction.wanted != PairState::Scaled ||
             (correction.wantedScaling.scalingVDW == correction.countedScaling.scalingVDW &&
              correction.wantedScaling.scalingCoulomb == correction.countedScaling.scalingCoulomb &&
              correction.wantedScaling.pair14 == correction.countedScaling.pair14)))
          continue;

        auto [it, inserted] = corrections.try_emplace(key, correction);
        if (!inserted && newer(k, it->second.bondId)) it->second = correction;
      }
    }
  }
  pairCorrections.reserve(corrections.size());
  for (const auto &[key, correction] : corrections) pairCorrections.push_back(correction);

  // The bonded terms of the molecules that contain a retyped atom.
  struct RetypedMolecule
  {
    std::unordered_map<std::uint32_t, std::uint64_t> atoms;  // local atom -> step of the retyping
  };
  std::map<std::pair<std::uint32_t, std::uint32_t>, RetypedMolecule> retyped{};
  for (std::size_t i = 0; i < siteState.size(); ++i)
  {
    if (siteState[i] == static_cast<std::uint8_t>(PolymerSiteState::Monomer)) continue;
    const PolymerAtomRef r = siteRef(i);
    const Component &component = components[r.componentId];
    const PolymerizationSite &definition = component.polymerizationSites[r.atomIndex];
    const std::size_t start = layout.moleculeStart(r.componentId, r.moleculeIndex);
    RetypedMolecule &molecule = retyped[{r.componentId, r.moleculeIndex}];
    for (std::size_t local : {definition.head, definition.tail})
    {
      if (atoms[start + local].type != component.atoms[local].type)
      {
        auto [it, inserted] = molecule.atoms.try_emplace(static_cast<std::uint32_t>(local), siteStep[i]);
        if (!inserted) it->second = std::max(it->second, siteStep[i]);
      }
    }
  }
  for (const auto &[moleculeKey, molecule] : retyped)
  {
    if (molecule.atoms.empty()) continue;
    const Component &component = components[moleculeKey.first];
    const std::size_t start = layout.moleculeStart(moleculeKey.first, moleculeKey.second);
    auto affected = [&](std::span<const std::size_t> identifiers, std::uint64_t &step)
    {
      bool any = false;
      step = 0;
      for (std::size_t id : identifiers)
      {
        const auto it = molecule.atoms.find(static_cast<std::uint32_t>(id));
        if (it != molecule.atoms.end())
        {
          any = true;
          step = std::max(step, it->second);
        }
      }
      return any;
    };
    for (const BondPotential &bond : component.intraMolecularPotentials.bonds)
    {
      std::uint64_t step{};
      if (!affected(bond.identifiers, step)) continue;
      const std::array<std::uint32_t, 2> global{static_cast<std::uint32_t>(start + bond.identifiers[0]),
                                                static_cast<std::uint32_t>(start + bond.identifiers[1])};
      const std::optional<BondPotential> replacement =
          templates.bond({atoms[global[0]].type, atoms[global[1]].type}, {global[0], global[1]});
      if (!replacement.has_value())
      {
        throw std::runtime_error(std::format(
            "[Polymerization]: no 'Bonds' template for the retyped bond {}(add [[\"{}\", \"{}\"], ...] to the "
            "'Polymerization' block)\n",
            describe(global), typeName(global[0]), typeName(global[1])));
      }
      BondPotential original = bond;
      original.identifiers = {global[0], global[1]};
      overrideBonds.push_back(OverrideBondTerm{global, original, replacement.value(), step});
    }
    for (const BendPotential &bend : component.intraMolecularPotentials.bends)
    {
      std::uint64_t step{};
      if (!affected(bend.identifiers, step)) continue;
      const std::array<std::uint32_t, 3> global{static_cast<std::uint32_t>(start + bend.identifiers[0]),
                                                static_cast<std::uint32_t>(start + bend.identifiers[1]),
                                                static_cast<std::uint32_t>(start + bend.identifiers[2])};
      const std::optional<BendPotential> replacement = templates.bend(
          {atoms[global[0]].type, atoms[global[1]].type, atoms[global[2]].type}, {global[0], global[1], global[2]});
      if (!replacement.has_value())
      {
        throw std::runtime_error(std::format(
            "[Polymerization]: no 'Bends' template for the retyped bend {}(add [[\"{}\", \"{}\", \"{}\"], ...] to "
            "the 'Polymerization' block)\n",
            describe(global), typeName(global[0]), typeName(global[1]), typeName(global[2])));
      }
      BendPotential original = bend;
      original.identifiers = {global[0], global[1], global[2]};
      overrideBends.push_back(OverrideBendTerm{global, original, replacement.value(), step});
    }
    for (const TorsionPotential &torsion : component.intraMolecularPotentials.torsions)
    {
      std::uint64_t step{};
      if (!affected(torsion.identifiers, step)) continue;
      const std::array<std::uint32_t, 4> global{static_cast<std::uint32_t>(start + torsion.identifiers[0]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[1]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[2]),
                                                static_cast<std::uint32_t>(start + torsion.identifiers[3])};
      const std::optional<TorsionPotential> replacement = templates.torsion(
          {atoms[global[0]].type, atoms[global[1]].type, atoms[global[2]].type, atoms[global[3]].type},
          {global[0], global[1], global[2], global[3]});
      if (!replacement.has_value())
      {
        throw std::runtime_error(std::format(
            "[Polymerization]: no 'Torsions' template for the retyped torsion {}(add [[\"{}\", \"{}\", \"{}\", "
            "\"{}\"], ...] to the 'Polymerization' block)\n",
            describe(global), typeName(global[0]), typeName(global[1]), typeName(global[2]), typeName(global[3])));
      }
      TorsionPotential original = torsion;
      original.identifiers = {global[0], global[1], global[2], global[3]};
      overrideTorsions.push_back(OverrideTorsionTerm{global, original, replacement.value(), step});
    }
  }

  derivedValid = true;
}

// ---- statistics and output ----------------------------------------------------------------------------------------

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
  for (std::size_t c = 0; c < sitesPerMolecule.size(); ++c)
  {
    if (sitesPerMolecule[c] > 0) moleculesWithSites += moleculesPerComponent[c];
  }
  if (moleculesWithSites > 0)
  {
    statistics.gelFraction = static_cast<double>(statistics.largestCluster) / static_cast<double>(moleculesWithSites);
    statistics.spanningFraction = static_cast<double>(spanningMolecules) / static_cast<double>(moleculesWithSites);
  }
  return statistics;
}

PolymerLinkTable::StereoStatistics PolymerLinkTable::stereoStatistics(const SimulationBox &simulationBox,
                                                                       const std::vector<Component> &components,
                                                                       std::span<const Atom> atoms) const
{
  StereoStatistics statistics{};
  if (!derivedValid || tailPartner.size() != siteState.size()) return statistics;

  // the sign of the chain-ordered chiral volume of every reacted stereocentre (0: not a stereocentre)
  std::vector<std::int8_t> sign(siteState.size(), 0);
  for (std::size_t s = 0; s < siteState.size(); ++s)
  {
    if (siteState[s] != static_cast<std::uint8_t>(PolymerSiteState::Dead) || tailPartner[s] == noAtom) continue;
    sign[s] = geometricChirality(s, tailPartner[s], simulationBox, components, atoms);
    if (sign[s] == 0) continue;
    ++statistics.stereoCentres;
    if (siteChirality[s] != 0)
    {
      ++statistics.restrained;
      if (siteChirality[s] != sign[s]) ++statistics.inverted;
    }
  }
  // the dyads (previous unit, this unit) and the triads (previous, this, next)
  auto dyadIsMeso = [&](std::size_t i, std::size_t j) -> std::optional<bool>
  {
    if (i == noSite || j == noSite || sign[i] == 0 || sign[j] == 0) return std::nullopt;
    return sign[i] == sign[j];
  };
  for (std::size_t s = 0; s < siteState.size(); ++s)
  {
    const std::optional<bool> dyad = dyadIsMeso(previousSite[s], s);
    if (!dyad.has_value()) continue;
    if (dyad.value()) ++statistics.mesoDyads; else ++statistics.racemoDyads;
    const std::optional<bool> next = dyadIsMeso(s, nextSite[s]);
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
  std::print(stream, "    sites:                           {} total, {} monomer, {} radical, {} dead\n",
             numberOfSites(), countSites(PolymerSiteState::Monomer), countSites(PolymerSiteState::Radical),
             countSites(PolymerSiteState::Dead));
  std::print(stream, "    conversion (reacted sites):      {:.4f}\n", conversion());
  std::print(stream, "    events:                          {} initiations, {} propagations, {} terminations\n",
             numberOfInitiations, numberOfPropagations, numberOfTerminations);
  std::print(stream, "    polymer bonds:                   {} ({} intramolecular / primary cycles)\n", bonds.size(),
             network.intramolecularBonds);
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
  for (std::size_t c = 0; c < sitesPerMolecule.size() && c < components.size(); ++c)
  {
    if (sitesPerMolecule[c] == 0) continue;
    std::size_t monomer = 0, radical = 0, dead = 0;
    for (std::size_t i = siteOffset[c]; i < siteOffset[c] + sitesPerMolecule[c] * moleculesPerComponent[c]; ++i)
    {
      switch (state(i))
      {
        case PolymerSiteState::Monomer:
          ++monomer;
          break;
        case PolymerSiteState::Radical:
          ++radical;
          break;
        case PolymerSiteState::Dead:
          ++dead;
          break;
      }
    }
    std::print(stream, "    component {} ({}):              {} sites per molecule, {} monomer, {} radical, {} dead\n", c,
               components[c].name, sitesPerMolecule[c], monomer, radical, dead);
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

nlohmann::json PolymerLinkTable::jsonStatus(const SimulationBox &simulationBox,
                                            const std::vector<Component> &components,
                                            std::span<const Atom> atoms) const
{
  nlohmann::json status{};
  const NetworkStatistics network = networkStatistics(simulationBox, components, atoms);
  status["numberOfSites"] = numberOfSites();
  status["monomerSites"] = countSites(PolymerSiteState::Monomer);
  status["radicalSites"] = countSites(PolymerSiteState::Radical);
  status["deadSites"] = countSites(PolymerSiteState::Dead);
  status["conversion"] = conversion();
  status["initiations"] = numberOfInitiations;
  status["propagations"] = numberOfPropagations;
  status["terminations"] = numberOfTerminations;
  status["bonds"] = bonds.size();
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
    entry["kind"] = bond.kind == PolymerBondKind::Propagation ? "propagation" : "termination";
    entry["step"] = bond.formedAtStep;
    entry["formationDistance"] = bond.formationDistance;
    topology.push_back(entry);
  }
  return topology;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerLinkTable &t)
{
  archive << t.versionNumber;
  archive << t.templates;
  archive << t.bonds;
  archive << t.activated;
  archive << t.sitesPerMolecule;
  archive << t.siteOffset;
  archive << t.moleculesPerComponent;
  archive << t.siteState;
  archive << t.siteStep;
  archive << t.relaxationSteps;
  archive << t.currentStep;
  archive << t.numberOfInitiations;
  archive << t.numberOfPropagations;
  archive << t.numberOfTerminations;
  archive << t.siteChirality;
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
  archive >> t.bonds;
  archive >> t.activated;
  archive >> t.sitesPerMolecule;
  archive >> t.siteOffset;
  archive >> t.moleculesPerComponent;
  archive >> t.siteState;
  archive >> t.siteStep;
  archive >> t.relaxationSteps;
  archive >> t.currentStep;
  archive >> t.numberOfInitiations;
  archive >> t.numberOfPropagations;
  archive >> t.numberOfTerminations;
  if (versionNumber >= 2)
  {
    archive >> t.siteChirality;
    archive >> t.chiralForceConstant;
    archive >> t.chiralReferenceVolume;
  }
  else
  {
    t.siteChirality.assign(t.siteState.size(), 0);
  }
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
