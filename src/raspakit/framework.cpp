module;

#if !defined(_WIN32)
#include <assert.h>
#endif

module framework;

import std;

import archive;
import int3;
import double3;
import double2;
import double4;
import double3x3;
import skposcarparser;
import skelement;
import characterset;
import stringutils;
import skparser;
import skposcarparser;
import skstructure;
import skatom;
import skcell;
import skspacegroup;
import forcefield;
import atom;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import connectivity_table;
import intra_molecular_potentials;
import units;
import json;
import molecule;
import simd_quatd;

namespace
{
std::filesystem::path frameworkDefinitionPath(const std::string& definitionName)
{
  std::filesystem::path directPath{definitionName};
  if (directPath.extension().empty()) directPath += ".json";
  if (std::filesystem::exists(directPath)) return directPath;

  std::filesystem::path relativePath = std::filesystem::path("frameworks") / definitionName / "framework.json";
  if (std::filesystem::exists(relativePath)) return relativePath;

  if (const char* raspaDirectory = std::getenv("RASPA_DIR"))
  {
    std::filesystem::path installedPath = std::filesystem::path(raspaDirectory) / relativePath;
    if (std::filesystem::exists(installedPath)) return installedPath;
  }

  throw std::runtime_error(std::format("[Framework reader]: definition '{}' not found (expected '{}')\n",
                                       definitionName, relativePath.string()));
}

template <std::size_t N>
std::vector<FrameworkPotentialDefinition<N>> readPotentialDefinitions(
    const ForceField& forceField, const nlohmann::basic_json<nlohmann::raspa_map>& data, const std::string& section)
{
  std::vector<FrameworkPotentialDefinition<N>> definitions;
  if (!data.contains(section)) return definitions;
  if (!data[section].is_array())
  {
    throw std::runtime_error(std::format("[Framework reader]: '{}' must be an array\n", section));
  }

  definitions.reserve(data[section].size());
  for (const auto& item : data[section])
  {
    if (!item.is_array() || item.size() != 3 || !item[0].is_array() || item[0].size() != N || !item[1].is_string() ||
        !item[2].is_array())
    {
      throw std::runtime_error(std::format(
          "[Framework reader]: item {} in '{}' must be [[{} atom-type strings], potential type, parameters]\n",
          item.dump(), section, N));
    }

    FrameworkPotentialDefinition<N> definition;
    for (std::size_t i = 0; i < N; ++i)
    {
      if (!item[0][i].is_string())
      {
        throw std::runtime_error(
            std::format("[Framework reader]: atom identifiers in '{}' must be strings\n", section));
      }
      const std::string atomType = item[0][i].template get<std::string>();
      const std::optional<std::size_t> type = forceField.findPseudoAtom(atomType);
      if (!type)
      {
        throw std::runtime_error(
            std::format("[Framework reader]: unknown pseudo-atom '{}' in '{}'\n", atomType, section));
      }
      definition.atomTypes[i] = *type;
    }
    definition.potentialType = item[1].template get<std::string>();
    definition.parameters = item[2].template get<std::vector<double>>();
    definitions.push_back(std::move(definition));
  }
  return definitions;
}

double readScalingFactor(const nlohmann::basic_json<nlohmann::raspa_map>& data, const std::string& key)
{
  if (!data.contains(key)) return 0.0;
  if (!data[key].is_number())
  {
    throw std::runtime_error(std::format("[Framework reader]: '{}' must be a number\n", key));
  }

  const double scaling = data[key].get<double>();
  if (!std::isfinite(scaling) || scaling < 0.0)
  {
    throw std::runtime_error(std::format("[Framework reader]: '{}' must be a finite non-negative number\n", key));
  }
  return scaling;
}

bool readExclusionOption(const nlohmann::basic_json<nlohmann::raspa_map>& data, const std::string& key,
                         bool defaultValue = true)
{
  if (!data.contains(key)) return defaultValue;
  if (!data[key].is_boolean())
  {
    throw std::runtime_error(std::format("[Framework reader]: '{}' must be a boolean\n", key));
  }
  return data[key].get<bool>();
}

int3 periodicImageShift(const SimulationBox& box, const double3& positionA, const double3& positionB)
{
  const double3 rawDisplacement = positionA - positionB;
  const double3 wrappedDisplacement = box.applyPeriodicBoundaryConditions(rawDisplacement);
  const double3 fractionalShift = box.inverseCell * (rawDisplacement - wrappedDisplacement);
  return int3(static_cast<int>(std::llround(fractionalShift.x)), static_cast<int>(std::llround(fractionalShift.y)),
              static_cast<int>(std::llround(fractionalShift.z)));
}

template <std::size_t N>
bool matchesTypes(const std::array<std::size_t, N>& indices, const std::array<std::size_t, N>& types,
                  const std::vector<Atom>& atoms)
{
  bool forward = true;
  bool reverse = true;
  for (std::size_t i = 0; i < N; ++i)
  {
    forward = forward && static_cast<std::size_t>(atoms[indices[i]].type) == types[i];
    reverse = reverse && static_cast<std::size_t>(atoms[indices[i]].type) == types[N - i - 1];
  }
  return forward || reverse;
}

template <std::size_t N>
std::optional<std::size_t> matchingDefinition(const std::array<std::size_t, N>& indices,
                                              const std::vector<FrameworkPotentialDefinition<N>>& definitions,
                                              const std::vector<Atom>& atoms, const std::string& section)
{
  std::optional<std::size_t> match;
  for (std::size_t i = 0; i < definitions.size(); ++i)
  {
    if (matchesTypes(indices, definitions[i].atomTypes, atoms))
    {
      if (match)
      {
        throw std::runtime_error(
            std::format("[Framework reader]: interaction in '{}' matches multiple definitions\n", section));
      }
      match = i;
    }
  }
  return match;
}

BondType bondType(const FrameworkPotentialDefinition<2>& definition)
{
  const auto iterator = BondPotential::definitionForString.find(definition.potentialType);
  if (iterator == BondPotential::definitionForString.end())
  {
    throw std::runtime_error(
        std::format("[Framework reader]: unknown bond potential '{}'\n", definition.potentialType));
  }
  if (definition.parameters.size() != BondPotential::numberOfBondParameters[std::to_underlying(iterator->second)])
  {
    throw std::runtime_error(std::format("[Framework reader]: wrong number of parameters for bond potential '{}'\n",
                                         definition.potentialType));
  }
  return iterator->second;
}

BendType bendType(const FrameworkPotentialDefinition<3>& definition)
{
  const auto iterator = BendPotential::definitionForString.find(definition.potentialType);
  if (iterator == BendPotential::definitionForString.end())
  {
    throw std::runtime_error(
        std::format("[Framework reader]: unknown bend potential '{}'\n", definition.potentialType));
  }
  if (definition.parameters.size() != BendPotential::numberOfBendParameters[std::to_underlying(iterator->second)])
  {
    throw std::runtime_error(std::format("[Framework reader]: wrong number of parameters for bend potential '{}'\n",
                                         definition.potentialType));
  }
  return iterator->second;
}

TorsionType torsionType(const FrameworkPotentialDefinition<4>& definition, const std::string& section)
{
  const auto iterator = TorsionPotential::definitionForString.find(definition.potentialType);
  if (iterator == TorsionPotential::definitionForString.end())
  {
    throw std::runtime_error(
        std::format("[Framework reader]: unknown {} potential '{}'\n", section, definition.potentialType));
  }
  if (definition.parameters.size() != TorsionPotential::numberOfTorsionParameters[std::to_underlying(iterator->second)])
  {
    throw std::runtime_error(std::format("[Framework reader]: wrong number of parameters for {} potential '{}'\n",
                                         section, definition.potentialType));
  }
  return iterator->second;
}

template <typename Potential>
std::string potentialSummary(const Potential& potential, const std::vector<Atom>& atoms, const ForceField& forceField,
                             bool improper = false)
{
  constexpr std::size_t numberOfAtoms = std::tuple_size_v<decltype(potential.identifiers)>;
  std::array<std::string, numberOfAtoms> atomTypes;
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    atomTypes[i] = forceField.pseudoAtoms[static_cast<std::size_t>(atoms[potential.identifiers[i]].type)].name;
  }

  if constexpr (numberOfAtoms == 4)
  {
    if (improper)
    {
      std::array<std::string, 3> outerTypes{atomTypes[0], atomTypes[2], atomTypes[3]};
      std::ranges::sort(outerTypes);
      atomTypes = {outerTypes[0], atomTypes[1], outerTypes[1], outerTypes[2]};
    }
    else
    {
      std::array<std::string, numberOfAtoms> reversedTypes = atomTypes;
      std::ranges::reverse(reversedTypes);
      atomTypes = std::min(atomTypes, reversedTypes);
    }
  }
  else
  {
    std::array<std::string, numberOfAtoms> reversedTypes = atomTypes;
    std::ranges::reverse(reversedTypes);
    atomTypes = std::min(atomTypes, reversedTypes);
  }

  std::string description = potential.print();
  const std::size_t separator = description.find(" : ");
  if (separator != std::string::npos) description.erase(0, separator + 3);
  if (!description.empty() && description.back() == '\n') description.pop_back();

  std::string typeDescription = atomTypes[0];
  for (std::size_t i = 1; i < numberOfAtoms; ++i)
  {
    typeDescription += std::format(" - {}", atomTypes[i]);
  }
  return std::format("{} : {}", typeDescription, description);
}

template <typename Potential>
void printPotentialSummaries(std::ostringstream& stream, const std::vector<Potential>& potentials,
                             const std::vector<Atom>& atoms, const ForceField& forceField, bool improper = false)
{
  std::map<std::string, std::size_t> summaries;
  for (const Potential& potential : potentials)
  {
    ++summaries[potentialSummary(potential, atoms, forceField, improper)];
  }
  for (const auto& [description, count] : summaries)
  {
    std::print(stream, "        {} ({} found)\n", description, count);
  }
}
}  // namespace

// default constructor, needed for binary restart-file
Framework::Framework() {}

Framework::Framework(const ForceField& forceField, std::string structureName, SimulationBox simulationBox,
                     std::size_t spaceGroupHallNumber, const std::vector<Atom>& definedAtoms,
                     const std::vector<Atom>& fractionalUnitCellAtoms, int3 numberOfUnitCells) noexcept(false)
    : simulationBox(simulationBox),
      spaceGroupHallNumber(spaceGroupHallNumber),
      numberOfUnitCells(numberOfUnitCells),
      name(structureName),
      definedAtoms(definedAtoms),
      fractionalUnitCellAtoms(fractionalUnitCellAtoms)
{
  for (std::size_t i = 0; i < definedAtoms.size(); ++i)
  {
    this->definedAtoms[i].moleculeId = 0;
  }

  unitCellAtoms = fractionalUnitCellAtoms;
  for (std::size_t i = 0; i < fractionalUnitCellAtoms.size(); ++i)
  {
    unitCellAtoms[i].position = simulationBox.cell * fractionalUnitCellAtoms[i].position;
  }

  makeSuperCell();

  unitCellMass = 0.0;
  for (const Atom& atom : unitCellAtoms)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    unitCellMass += forceField.pseudoAtoms[atomType].mass;
  }

  mass = 0.0;
  netCharge = 0.0;
  smallestCharge = 0.0;
  largestCharge = 0.0;
  if (!atoms.empty())
  {
    smallestCharge = std::numeric_limits<double>::max();
    largestCharge = std::numeric_limits<double>::lowest();
    for (const Atom& atom : atoms)
    {
      std::size_t atomType = static_cast<std::size_t>(atom.type);
      mass += forceField.pseudoAtoms[atomType].mass;
      netCharge += atom.charge;

      if (atom.charge > largestCharge) largestCharge = atom.charge;
      if (atom.charge < smallestCharge) smallestCharge = atom.charge;
    }
  }

  for (std::size_t i = 0; i < unitCellAtoms.size(); ++i)
  {
    unitCellAtoms[i].componentId = static_cast<std::uint8_t>(0);
    unitCellAtoms[i].moleculeId = 0;
  }

  determineUniqueAtomTypes();
}

Framework::Framework(const ForceField& forceField, std::string structureName, SimulationBox simulationBox,
                     std::size_t spaceGroupHallNumber, const std::vector<Atom>& definedAtoms,
                     int3 numberOfUnitCells) noexcept(false)
    : Framework(forceField, structureName, simulationBox, spaceGroupHallNumber, definedAtoms,
                CIFReader::expandDefinedAtomsToUnitCell(simulationBox, spaceGroupHallNumber, definedAtoms),
                numberOfUnitCells)
{
}

void Framework::determineUniqueAtomTypes()
{
  std::unordered_set<Atom, AtomTypeEqual, AtomTypeEqual> atom_types(atoms.begin(), atoms.end());

  uniqueAtomTypes.reserve(atom_types.size());
  std::transform(atom_types.begin(), atom_types.end(), std::inserter(uniqueAtomTypes, uniqueAtomTypes.begin()),
                 [](const Atom& atom) { return static_cast<std::size_t>(atom.type); });
}

void Framework::readFrameworkDefinition(const ForceField& forceField, const std::string& definitionName)
{
  frameworkDefinitionName = definitionName;
  const std::filesystem::path path = frameworkDefinitionPath(definitionName);
  std::ifstream input(path);
  if (!input)
  {
    throw std::runtime_error(std::format("[Framework reader]: unable to open '{}'\n", path.string()));
  }

  nlohmann::basic_json<nlohmann::raspa_map> data;
  try
  {
    data = nlohmann::json::parse(input);
  }
  catch (const nlohmann::json::parse_error& error)
  {
    throw std::runtime_error(
        std::format("[Framework reader]: JSON parse error in '{}' at byte {}\n", path.string(), error.byte));
  }

  rigid = true;
  if (data.contains("Type"))
  {
    if (!data["Type"].is_string())
    {
      throw std::runtime_error("[Framework reader]: 'Type' must be a string\n");
    }
    const std::string type = data["Type"].get<std::string>();
    if (caseInSensStringCompare(type, "Flexible"))
    {
      rigid = false;
    }
    else if (!caseInSensStringCompare(type, "Rigid"))
    {
      throw std::runtime_error(std::format("[Framework reader]: unknown framework type '{}'\n", type));
    }
  }

  // Retain the parsed type-based potential definitions and the nonbonded generation options so the
  // intramolecular topology can be re-derived from memory (e.g. on a reduced primitive cell) without
  // re-reading the definition file.
  bondDefinitions = readPotentialDefinitions<2>(forceField, data, "Bonds");
  bendDefinitions = readPotentialDefinitions<3>(forceField, data, "Bends");
  torsionDefinitions = readPotentialDefinitions<4>(forceField, data, "Torsions");
  improperTorsionDefinitions = readPotentialDefinitions<4>(forceField, data, "ImproperTorsions");

  excludeIntra12Interactions = readExclusionOption(data, "ExcludeIntra12Interactions");
  excludeIntra13Interactions = readExclusionOption(data, "ExcludeIntra13Interactions");
  excludeIntraBondInteractions = readExclusionOption(data, "ExcludeIntraBondInteractions", false);
  excludeIntraBendInteractions = readExclusionOption(data, "ExcludeIntraBendInteractions", false);
  intra14VanDerWaalsScaling = readScalingFactor(data, "Intra14VanDerWaalsScalingValue");
  intra14ChargeChargeScaling = readScalingFactor(data, "Intra14ChargeChargeScalingValue");

  readGroups(data);
  if (mixed && rigid)
  {
    throw std::runtime_error(
        "[Framework reader]: 'Groups' require Type 'Flexible'; cannot combine Groups with a rigid framework\n");
  }
  computeGroupRigidProperties(forceField);

  generateIntraMolecularPotentials(forceField);
}

void Framework::readGroups(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  groups.clear();
  atomGroupIds.clear();
  mixed = false;
  fixedAtomCount = 0;
  rigidGroupAtomCount = 0;
  flexibleAtomCount = 0;
  rigidGroupCount = 0;

  if (!parsed_data.contains("Groups"))
  {
    return;
  }

  if (rigid)
  {
    throw std::runtime_error(
        "[Framework reader]: 'Groups' require Type 'Flexible'; cannot combine Groups with Type 'Rigid'\n");
  }

  if (!parsed_data["Groups"].is_array())
  {
    throw std::runtime_error("[Framework reader]: 'Groups' must be an array\n");
  }

  const std::size_t numberOfAtoms = atoms.size();
  std::vector<std::make_signed_t<std::size_t>> assigned(numberOfAtoms, -1);

  for (auto &[_, item] : parsed_data["Groups"].items())
  {
    if (!item.is_object())
    {
      throw std::runtime_error(
          "[Framework reader]: each entry of 'Groups' must be an object with keys 'Type' and 'Atoms'\n");
    }

    if (!item.contains("Type") || !item["Type"].is_string())
    {
      throw std::runtime_error("[Framework reader]: each group must contain a string 'Type'\n");
    }
    const std::string typeString = item["Type"].get<std::string>();
    FrameworkGroupType groupType = FrameworkGroupType::Flexible;
    if (caseInSensStringCompare(typeString, "Fixed"))
    {
      groupType = FrameworkGroupType::Fixed;
    }
    else if (caseInSensStringCompare(typeString, "Rigid"))
    {
      groupType = FrameworkGroupType::Rigid;
    }
    else if (caseInSensStringCompare(typeString, "Flexible"))
    {
      groupType = FrameworkGroupType::Flexible;
    }
    else if (caseInSensStringCompare(typeString, "Cycle"))
    {
      throw std::runtime_error("[Framework reader]: group Type 'Cycle' is not supported for frameworks\n");
    }
    else
    {
      throw std::runtime_error(std::format(
          "[Framework reader]: group 'Type' must be 'Fixed', 'Rigid', or 'Flexible', got {}\n", typeString));
    }

    if (!item.contains("Atoms") || !item["Atoms"].is_array())
    {
      throw std::runtime_error("[Framework reader]: each group must contain an 'Atoms' array\n");
    }

    std::vector<std::size_t> groupAtoms = item["Atoms"].get<std::vector<std::size_t>>();
    if (groupAtoms.empty())
    {
      throw std::runtime_error("[Framework reader]: a group must contain at least one atom\n");
    }

    const std::size_t groupIndex = groups.size();
    for (std::size_t atom : groupAtoms)
    {
      if (atom >= numberOfAtoms)
      {
        throw std::runtime_error(std::format(
            "[Framework reader]: group atom index {} out of range (framework has {} atoms)\n", atom, numberOfAtoms));
      }
      if (assigned[atom] != -1)
      {
        throw std::runtime_error(
            std::format("[Framework reader]: atom {} is listed in more than one group\n", atom));
      }
      assigned[atom] = static_cast<std::make_signed_t<std::size_t>>(groupIndex);
    }

    groups.emplace_back(groupType, std::move(groupAtoms));
  }

  for (std::size_t atom = 0; atom != numberOfAtoms; ++atom)
  {
    if (assigned[atom] == -1)
    {
      throw std::runtime_error(std::format(
          "[Framework reader]: atom {} is not assigned to any group; 'Groups' must cover all atoms\n", atom));
    }
  }

  // Reorder atoms to [Fixed | Rigid-group atoms | Flexible] and remap group atom indices.
  std::vector<Atom> reordered;
  reordered.reserve(numberOfAtoms);

  const auto appendByType = [&](FrameworkGroupType type)
  {
    for (FrameworkGroup &group : groups)
    {
      if (group.type != type) continue;
      std::vector<std::size_t> remapped;
      remapped.reserve(group.atoms.size());
      for (std::size_t oldIndex : group.atoms)
      {
        remapped.push_back(reordered.size());
        reordered.push_back(atoms[oldIndex]);
      }
      group.atoms = std::move(remapped);
    }
  };
  appendByType(FrameworkGroupType::Fixed);
  appendByType(FrameworkGroupType::Rigid);
  appendByType(FrameworkGroupType::Flexible);
  atoms = std::move(reordered);

  atomGroupIds.assign(numberOfAtoms, 0);
  for (std::size_t groupIndex = 0; groupIndex < groups.size(); ++groupIndex)
  {
    for (std::size_t atom : groups[groupIndex].atoms)
    {
      atomGroupIds[atom] = groupIndex;
    }
  }

  finalizeGroupCache();
}

void Framework::finalizeGroupCache()
{
  fixedAtomCount = 0;
  rigidGroupAtomCount = 0;
  flexibleAtomCount = 0;
  rigidGroupCount = 0;
  for (const FrameworkGroup &group : groups)
  {
    if (group.isFixed())
    {
      fixedAtomCount += group.atoms.size();
    }
    else if (group.isRigidBody())
    {
      ++rigidGroupCount;
      rigidGroupAtomCount += group.atoms.size();
    }
    else
    {
      flexibleAtomCount += group.atoms.size();
    }
  }
  mixed = !groups.empty();
}

void Framework::computeGroupRigidProperties(const ForceField& forceField)
{
  finalizeGroupCache();

  std::vector<double3> referencePositions(atoms.size());
  std::vector<double> masses(atoms.size());
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    referencePositions[i] = atoms[i].position;
    masses[i] = forceField.pseudoAtoms[atoms[i].type].mass;
  }

  for (FrameworkGroup &group : groups)
  {
    if (!group.isRigidBody())
    {
      // Fixed/Flexible groups carry no rigid-body data; reset the Fragment part to its point-mass
      // defaults (keeping only the atom subset).
      static_cast<Fragment &>(group) = Fragment(std::move(group.atoms));
      continue;
    }
    group.computeRigidProperties(referencePositions, masses);
  }
}

std::size_t Framework::numberOfFixedAtoms() const
{
  if (rigid) return atoms.size();
  if (mixed) return fixedAtomCount;
  return 0;
}

std::size_t Framework::numberOfMobileAtoms() const
{
  if (rigid) return 0;
  if (mixed) return rigidGroupAtomCount + flexibleAtomCount;
  return atoms.size();
}

bool Framework::isInsideFixedOrRigidGroup(std::span<const std::size_t> ids) const
{
  if (atomGroupIds.empty() || ids.empty()) return false;

  const std::size_t g = atomGroupIds[ids[0]];
  if (!groups[g].isFixed() && !groups[g].isRigidBody()) return false;
  for (std::size_t id : ids)
  {
    if (atomGroupIds[id] != g) return false;
  }
  return true;
}

std::optional<std::size_t> Framework::rigidGroupContaining(std::size_t atom) const
{
  if (atomGroupIds.empty() || atom >= atomGroupIds.size()) return std::nullopt;
  const std::size_t g = atomGroupIds[atom];
  if (!groups[g].isRigidBody()) return std::nullopt;
  return g;
}

std::optional<std::size_t> Framework::fixedGroupContaining(std::size_t atom) const
{
  if (atomGroupIds.empty() || atom >= atomGroupIds.size()) return std::nullopt;
  const std::size_t g = atomGroupIds[atom];
  if (!groups[g].isFixed()) return std::nullopt;
  return g;
}

bool Framework::isFixedAtom(std::size_t atom) const
{
  if (rigid) return true;
  if (!mixed) return false;
  return fixedGroupContaining(atom).has_value();
}

bool Framework::isFlexibleAtom(std::size_t atom) const
{
  if (rigid) return false;
  if (!mixed) return true;
  if (atomGroupIds.empty() || atom >= atomGroupIds.size()) return false;
  return groups[atomGroupIds[atom]].isFlexible();
}

void Framework::regenerateGroupAtoms(const GroupState &state, std::size_t groupIndex,
                                     std::span<Atom> frameworkAtoms) const
{
  groups[groupIndex].regenerateAtoms(state, frameworkAtoms);
}

GroupState Framework::deriveGroupState(std::size_t groupIndex, std::span<const Atom> frameworkAtoms) const
{
  return groups[groupIndex].deriveState(frameworkAtoms);
}

void Framework::generateIntraMolecularPotentials(const ForceField& forceField)
{
  std::vector<BondType> bondTypes;
  std::vector<BendType> bendTypes;
  std::vector<TorsionType> torsionTypes;
  std::vector<TorsionType> improperTypes;
  bondTypes.reserve(bondDefinitions.size());
  bendTypes.reserve(bendDefinitions.size());
  torsionTypes.reserve(torsionDefinitions.size());
  improperTypes.reserve(improperTorsionDefinitions.size());
  std::ranges::transform(bondDefinitions, std::back_inserter(bondTypes),
                         [](const auto& definition) { return bondType(definition); });
  std::ranges::transform(bendDefinitions, std::back_inserter(bendTypes),
                         [](const auto& definition) { return bendType(definition); });
  std::ranges::transform(torsionDefinitions, std::back_inserter(torsionTypes),
                         [](const auto& definition) { return torsionType(definition, "torsion"); });
  std::ranges::transform(improperTorsionDefinitions, std::back_inserter(improperTypes),
                         [](const auto& definition) { return torsionType(definition, "improper torsion"); });

  ConnectivityTable connectivity(atoms.size());
  const SimulationBox superCellBox = simulationBox.scaled(numberOfUnitCells);
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    const std::size_t typeA = static_cast<std::size_t>(atoms[i].type);
    const std::size_t atomicNumberA = forceField.pseudoAtoms.at(typeA).atomicNumber;
    const double radiusA = PredefinedElements::predefinedElements.at(atomicNumberA)._covalentRadius;
    for (std::size_t j = i + 1; j < atoms.size(); ++j)
    {
      const std::size_t typeB = static_cast<std::size_t>(atoms[j].type);
      const std::size_t atomicNumberB = forceField.pseudoAtoms.at(typeB).atomicNumber;
      const double radiusB = PredefinedElements::predefinedElements.at(atomicNumberB)._covalentRadius;
      const double cutoff = radiusA + radiusB + 0.56;
      const double3 dr = superCellBox.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
      if (double3::dot(dr, dr) < cutoff * cutoff)
      {
        connectivity[i, j] = true;
        connectivity[j, i] = true;
      }
    }
  }

  // Rigid groups must form a connected subgraph so they can be integrated as a single rigid body.
  for (const FrameworkGroup& group : groups)
  {
    if (!group.isRigidBody()) continue;
    if (!connectivity.checkIsConnectedSubgraph(group.atoms))
    {
      std::stringstream result{};
      std::copy(group.atoms.begin(), group.atoms.end(), std::ostream_iterator<std::size_t>(result, " "));
      throw std::runtime_error(
          std::format("[Framework reader]: Rigid group ({}) is not a connected subgraph\n", result.str()));
    }
  }

  Potentials::IntraMolecularPotentials potentials;
  for (const std::array<std::size_t, 2>& indices : connectivity.findAllBonds())
  {
    if (isInsideFixedOrRigidGroup(indices)) continue;
    if (const std::optional<std::size_t> match = matchingDefinition(indices, bondDefinitions, atoms, "Bonds"))
    {
      potentials.bonds.emplace_back(indices, bondTypes[*match], bondDefinitions[*match].parameters);
    }
  }
  for (const std::array<std::size_t, 3>& indices : connectivity.findAllBends())
  {
    if (isInsideFixedOrRigidGroup(indices)) continue;
    if (const std::optional<std::size_t> match = matchingDefinition(indices, bendDefinitions, atoms, "Bends"))
    {
      potentials.bends.emplace_back(indices, bendTypes[*match], bendDefinitions[*match].parameters);
    }
  }
  for (const std::array<std::size_t, 4>& indices : connectivity.findAllTorsions())
  {
    if (isInsideFixedOrRigidGroup(indices)) continue;
    if (const std::optional<std::size_t> match = matchingDefinition(indices, torsionDefinitions, atoms, "Torsions"))
    {
      potentials.torsions.emplace_back(indices, torsionTypes[*match], torsionDefinitions[*match].parameters);
    }
  }

  std::map<std::array<std::size_t, 4>, std::size_t> improperOwners;
  for (std::size_t definitionIndex = 0; definitionIndex < improperTorsionDefinitions.size(); ++definitionIndex)
  {
    const auto& definition = improperTorsionDefinitions[definitionIndex];
    for (std::size_t center = 0; center < atoms.size(); ++center)
    {
      if (static_cast<std::size_t>(atoms[center].type) != definition.atomTypes[1]) continue;
      const std::vector<std::size_t> neighbors = connectivity.findAllNeighbors(center);
      for (const std::size_t atomA : neighbors)
      {
        if (static_cast<std::size_t>(atoms[atomA].type) != definition.atomTypes[0]) continue;
        for (const std::size_t atomC : neighbors)
        {
          if (atomC == atomA || static_cast<std::size_t>(atoms[atomC].type) != definition.atomTypes[2]) continue;
          for (const std::size_t atomD : neighbors)
          {
            if (atomD == atomA || atomD == atomC ||
                static_cast<std::size_t>(atoms[atomD].type) != definition.atomTypes[3])
            {
              continue;
            }

            std::array<std::size_t, 3> outerAtoms{atomA, atomC, atomD};
            std::ranges::sort(outerAtoms);
            const std::array<std::size_t, 4> key{center, outerAtoms[0], outerAtoms[1], outerAtoms[2]};
            if (const auto owner = improperOwners.find(key); owner != improperOwners.end())
            {
              if (owner->second != definitionIndex)
              {
                throw std::runtime_error("[Framework reader]: improper torsion matches multiple definitions\n");
              }
              continue;
            }

            const std::array<std::size_t, 4> improperIds{atomA, center, atomC, atomD};
            improperOwners.emplace(key, definitionIndex);
            if (isInsideFixedOrRigidGroup(improperIds)) continue;
            potentials.improperTorsions.emplace_back(improperIds, improperTypes[definitionIndex], definition.parameters);
          }
        }
      }
    }
  }

  const auto addVanDerWaals = [&](const std::array<std::size_t, 2>& indices, double scaling)
  {
    const std::size_t typeA = static_cast<std::size_t>(atoms[indices[0]].type);
    const std::size_t typeB = static_cast<std::size_t>(atoms[indices[1]].type);
    const double4 parameters = forceField(typeA, typeB).parameters;
    potentials.vanDerWaals.emplace_back(
        indices, VanDerWaalsType::LennardJones,
        std::vector<double>{parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, scaling);
  };

  const auto addCoulomb = [&](const std::array<std::size_t, 2>& indices, double scaling)
  {
    if (forceField.useCharge)
    {
      potentials.coulombs.emplace_back(indices, CoulombType::Coulomb, atoms[indices[0]].charge,
                                       atoms[indices[1]].charge, scaling);
    }
  };

  std::set<std::array<std::size_t, 2>> pairs12;
  for (const std::array<std::size_t, 2>& bond : connectivity.findAllBonds())
  {
    pairs12.insert({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
  }
  std::set<std::array<std::size_t, 2>> definedBondPairs;
  for (const BondPotential& bond : potentials.bonds)
  {
    definedBondPairs.insert(
        {std::min(bond.identifiers[0], bond.identifiers[1]), std::max(bond.identifiers[0], bond.identifiers[1])});
  }
  std::set<std::array<std::size_t, 2>> pairs13;
  for (const std::array<std::size_t, 3>& bend : connectivity.findAllBends())
  {
    const std::array<std::size_t, 2> pair{std::min(bend[0], bend[2]), std::max(bend[0], bend[2])};
    if (!pairs12.contains(pair)) pairs13.insert(pair);
  }
  std::set<std::array<std::size_t, 2>> definedBendPairs;
  for (const BendPotential& bend : potentials.bends)
  {
    const std::array<std::size_t, 2> pair{std::min(bend.identifiers[0], bend.identifiers[2]),
                                         std::max(bend.identifiers[0], bend.identifiers[2])};
    definedBendPairs.insert(pair);
  }
  std::set<std::array<std::size_t, 2>> pairs14;
  for (const std::array<std::size_t, 4>& torsion : connectivity.findAllTorsions())
  {
    const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
    if (!pairs12.contains(pair) && !pairs13.contains(pair)) pairs14.insert(pair);
  }

  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j < atoms.size(); ++j)
    {
      const std::array<std::size_t, 2> pair{i, j};
      if (isInsideFixedOrRigidGroup(pair)) continue;
      if ((excludeIntra12Interactions && pairs12.contains(pair)) ||
          (excludeIntraBondInteractions && definedBondPairs.contains(pair)) ||
          (excludeIntra13Interactions && pairs13.contains(pair)) ||
          (excludeIntraBendInteractions && definedBendPairs.contains(pair)))
      {
        continue;
      }
      if (pairs14.contains(pair))
      {
        if (intra14VanDerWaalsScaling > 0.0) addVanDerWaals(pair, intra14VanDerWaalsScaling);
        if (intra14ChargeChargeScaling > 0.0) addCoulomb(pair, intra14ChargeChargeScaling);
      }
      else
      {
        addVanDerWaals(pair, 1.0);
        addCoulomb(pair, 1.0);
      }
    }
  }

  std::map<std::array<std::size_t, 2>, int3> directedBondShifts;
  for (const std::array<std::size_t, 2>& bond : connectivity.findAllBonds())
  {
    const int3 shift = periodicImageShift(superCellBox, atoms[bond[0]].position, atoms[bond[1]].position);
    directedBondShifts[{bond[0], bond[1]}] = shift;
    directedBondShifts[{bond[1], bond[0]}] = -shift;
  }

  FrameworkIntraMolecularImageShifts imageShifts;
  imageShifts.bonds.reserve(potentials.bonds.size());
  for (const BondPotential& bond : potentials.bonds)
  {
    imageShifts.bonds.push_back({int3{}, directedBondShifts.at(bond.identifiers)});
  }
  imageShifts.bends.reserve(potentials.bends.size());
  for (const BendPotential& bend : potentials.bends)
  {
    const auto& ids = bend.identifiers;
    imageShifts.bends.push_back(
        {directedBondShifts.at({ids[1], ids[0]}), int3{}, directedBondShifts.at({ids[1], ids[2]})});
  }
  imageShifts.torsions.reserve(potentials.torsions.size());
  for (const TorsionPotential& torsion : potentials.torsions)
  {
    const auto& ids = torsion.identifiers;
    const int3 shiftC = directedBondShifts.at({ids[1], ids[2]});
    imageShifts.torsions.push_back(
        {directedBondShifts.at({ids[1], ids[0]}), int3{}, shiftC, shiftC + directedBondShifts.at({ids[2], ids[3]})});
  }
  imageShifts.improperTorsions.reserve(potentials.improperTorsions.size());
  for (const TorsionPotential& torsion : potentials.improperTorsions)
  {
    const auto& ids = torsion.identifiers;
    imageShifts.improperTorsions.push_back({directedBondShifts.at({ids[1], ids[0]}), int3{},
                                            directedBondShifts.at({ids[1], ids[2]}),
                                            directedBondShifts.at({ids[1], ids[3]})});
  }
  imageShifts.vanDerWaals.reserve(potentials.vanDerWaals.size());
  for (const VanDerWaalsPotential& potential : potentials.vanDerWaals)
  {
    const auto& ids = potential.identifiers;
    imageShifts.vanDerWaals.push_back(
        {int3{}, periodicImageShift(superCellBox, atoms[ids[0]].position, atoms[ids[1]].position)});
  }
  imageShifts.coulombs.reserve(potentials.coulombs.size());
  for (const CoulombPotential& potential : potentials.coulombs)
  {
    const auto& ids = potential.identifiers;
    imageShifts.coulombs.push_back(
        {int3{}, periodicImageShift(superCellBox, atoms[ids[0]].position, atoms[ids[1]].position)});
  }

  connectivityTable = std::move(connectivity);
  intraMolecularPotentials = std::move(potentials);
  intraMolecularImageShifts = std::move(imageShifts);
}

void Framework::regenerateVanDerWaalsImageList(const ForceField& forceField, const SimulationBox& box,
                                               double cutOffFrameworkVDW)
{
  if (rigid || cutOffFrameworkVDW <= 0.0 || atoms.empty()) return;

  // Only frameworks whose construction generated a nonbonded van der Waals pair list are replicated.
  // An empty list means the framework has no nonbonded van der Waals interactions (e.g. a purely bonded
  // model or a programmatically assembled framework), so there is nothing to extend over replica cells.
  if (intraMolecularPotentials.vanDerWaals.empty()) return;

  const double3 widths = box.perpendicularWidths();
  const double smallestWidth = std::min({widths.x, widths.y, widths.z});
  // A single minimum image already covers the cutoff; keep the construction-time list untouched.
  if (smallestWidth >= 2.0 * cutOffFrameworkVDW) return;

  const double cutOffSquared = cutOffFrameworkVDW * cutOffFrameworkVDW;

  // Topological classification (unordered pairs), reconstructed from the stored connectivity exactly as
  // in the construction-time nonbonded generation.
  std::set<std::array<std::size_t, 2>> pairs12;
  std::set<std::array<std::size_t, 2>> pairs13;
  std::set<std::array<std::size_t, 2>> pairs14;
  if (!connectivityTable.table.empty())
  {
    for (const std::array<std::size_t, 2>& bond : connectivityTable.findAllBonds())
    {
      pairs12.insert({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
    }
    for (const std::array<std::size_t, 3>& bend : connectivityTable.findAllBends())
    {
      const std::array<std::size_t, 2> pair{std::min(bend[0], bend[2]), std::max(bend[0], bend[2])};
      if (!pairs12.contains(pair)) pairs13.insert(pair);
    }
    for (const std::array<std::size_t, 4>& torsion : connectivityTable.findAllTorsions())
    {
      const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
      if (!pairs12.contains(pair) && !pairs13.contains(pair)) pairs14.insert(pair);
    }
  }
  std::set<std::array<std::size_t, 2>> definedBondPairs;
  for (const BondPotential& bond : intraMolecularPotentials.bonds)
  {
    definedBondPairs.insert(
        {std::min(bond.identifiers[0], bond.identifiers[1]), std::max(bond.identifiers[0], bond.identifiers[1])});
  }
  std::set<std::array<std::size_t, 2>> definedBendPairs;
  for (const BendPotential& bend : intraMolecularPotentials.bends)
  {
    definedBendPairs.insert(
        {std::min(bend.identifiers[0], bend.identifiers[2]), std::max(bend.identifiers[0], bend.identifiers[2])});
  }

  std::vector<VanDerWaalsPotential> vanDerWaals;
  std::vector<std::array<int3, 2>> vanDerWaalsImages;

  const auto emitPair = [&](std::size_t i, std::size_t j, const int3& shift, double scaling)
  {
    const std::size_t typeA = static_cast<std::size_t>(atoms[i].type);
    const std::size_t typeB = static_cast<std::size_t>(atoms[j].type);
    const double4 parameters = forceField(typeA, typeB).parameters;
    vanDerWaals.emplace_back(
        std::array<std::size_t, 2>{i, j}, VanDerWaalsType::LennardJones,
        std::vector<double>{parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, scaling);
    vanDerWaalsImages.push_back({int3{}, shift});
  };

  // Range of replica shells to scan; ceil(2 cutoff / width) is a safe upper bound per axis.
  const int3 shells = box.smallestNumberOfUnitCellsForMinimumImagesConvention(cutOffFrameworkVDW);

  // Distinct-atom pairs: atom i sits in the home cell, atom j roams over all periodic images. Fixing
  // i in the home cell and iterating i < j counts every unordered pair-image exactly once.
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j < atoms.size(); ++j)
    {
      const std::array<std::size_t, 2> pair{i, j};
      const bool excludedBonded = (excludeIntra12Interactions && pairs12.contains(pair)) ||
                                  (excludeIntraBondInteractions && definedBondPairs.contains(pair)) ||
                                  (excludeIntra13Interactions && pairs13.contains(pair)) ||
                                  (excludeIntraBendInteractions && definedBendPairs.contains(pair));
      const bool scaled14 = pairs14.contains(pair);
      // The covalently bonded neighbour is the minimum image; only that image is excluded or 1-4 scaled.
      const int3 bondedImage = periodicImageShift(box, atoms[i].position, atoms[j].position);

      for (int a = -shells.x; a <= shells.x; ++a)
      {
        for (int b = -shells.y; b <= shells.y; ++b)
        {
          for (int c = -shells.z; c <= shells.z; ++c)
          {
            const int3 shift(a, b, c);
            const bool isBondedImage = (shift == bondedImage);
            if (isBondedImage && excludedBonded) continue;
            if (isBondedImage && scaled14)
            {
              if (intra14VanDerWaalsScaling > 0.0) emitPair(i, j, shift, intra14VanDerWaalsScaling);
              continue;
            }
            const double3 dr = atoms[i].position - (atoms[j].position + box.cell * double3(static_cast<double>(a),
                                                                                           static_cast<double>(b),
                                                                                           static_cast<double>(c)));
            if (double3::dot(dr, dr) < cutOffSquared) emitPair(i, j, shift, 1.0);
          }
        }
      }
    }
  }

  // Self-images: an atom interacts with its own periodic images. n and -n describe the same unordered
  // pair, so only a positive half-shell is enumerated to avoid double counting. Self-images are never
  // covalently bonded in the supercell atom indexing, so no exclusion applies.
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    for (int a = -shells.x; a <= shells.x; ++a)
    {
      for (int b = -shells.y; b <= shells.y; ++b)
      {
        for (int c = -shells.z; c <= shells.z; ++c)
        {
          const bool positiveHalf = (a > 0) || (a == 0 && b > 0) || (a == 0 && b == 0 && c > 0);
          if (!positiveHalf) continue;
          const double3 dr =
              box.cell * double3(static_cast<double>(-a), static_cast<double>(-b), static_cast<double>(-c));
          if (double3::dot(dr, dr) < cutOffSquared) emitPair(i, i, int3(a, b, c), 1.0);
        }
      }
    }
  }

  intraMolecularPotentials.vanDerWaals = std::move(vanDerWaals);
  intraMolecularImageShifts.vanDerWaals = std::move(vanDerWaalsImages);
}

void Framework::makeSuperCell()
{
  for (std::int32_t i = 0; i < numberOfUnitCells.x; ++i)
  {
    for (std::int32_t j = 0; j < numberOfUnitCells.y; ++j)
    {
      for (std::int32_t k = 0; k < numberOfUnitCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position +=
              simulationBox.cell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          atoms.push_back(atomCopy);
        }
      }
    }
  }
}

std::vector<Atom> Framework::makeSuperCell(int3 numberOfCells) const
{
  std::vector<Atom> superCellAtoms{};

  for (std::int32_t i = 0; i < numberOfCells.x; ++i)
  {
    for (std::int32_t j = 0; j < numberOfCells.y; ++j)
    {
      for (std::int32_t k = 0; k < numberOfCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position +=
              simulationBox.cell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          superCellAtoms.push_back(atomCopy);
        }
      }
    }
  }

  return superCellAtoms;
}

std::optional<double> Framework::computeLargestNonOverlappingFreeRadius(const ForceField& forceField,
                                                                        double3 probe_position,
                                                                        double well_depth_factor) const
{
  double smallest_radius = std::numeric_limits<double>::max();

  // if inside blockingpocket, then return
  //    return std::nullopt;

  for (const Atom& atom : unitCellAtoms)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    double size_parameter = forceField[atomType].sizeParameter();
    double3 dr = probe_position - atom.position;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    double mixing_radius = 0.5 * well_depth_factor * size_parameter;
    if (rr < mixing_radius * mixing_radius)
    {
      return std::nullopt;
    }

    double radius = std::sqrt(rr) - mixing_radius;
    smallest_radius = std::min(smallest_radius, radius);
  }

  return smallest_radius;
}

bool Framework::computeVanDerWaalsRadiusOverlap(const ForceField& forceField, double3 probe_position) const
{
  for (const Atom& atom : unitCellAtoms)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    std::size_t atomicNumber = forceField.pseudoAtoms[atomType].atomicNumber;
    double radius = PredefinedElements::predefinedElements[atomicNumber]._VDWRadius;
    double3 dr = probe_position - atom.position;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < radius * radius)
    {
      return true;
    }
  }

  return false;
}
bool Framework::computeOverlap(const ForceField& forceField, double3 probe_position, double well_depth_factor,
                               std::size_t probe_type, std::make_signed_t<std::size_t> skip) const
{
  for (std::make_signed_t<std::size_t> atom_index = 0; const Atom& atom : unitCellAtoms)
  {
    if (atom_index != skip)
    {
      std::size_t atomType = static_cast<std::size_t>(atom.type);
      double size_parameter = forceField(probe_type, atomType).sizeParameter();
      double3 dr = probe_position - atom.position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);

      double radius = well_depth_factor * size_parameter;
      if (rr < radius * radius)
      {
        return true;
      }
    }

    ++atom_index;
  }

  return false;
}

std::vector<double3> Framework::fractionalAtomPositionsUnitCell() const
{
  std::vector<double3> positions;
  positions.reserve(unitCellAtoms.size());

  double3x3 inverseCell = simulationBox.inverseCell;
  for (const Atom& atom : unitCellAtoms)
  {
    double3 s = (inverseCell * atom.position).fract();
    positions.push_back(s);
  }

  return positions;
}

std::vector<double3> Framework::cartesianAtomPositionsUnitCell() const
{
  std::vector<double3> positions;
  positions.reserve(unitCellAtoms.size());

  for (const Atom& atom : unitCellAtoms)
  {
    positions.push_back(atom.position);
  }

  return positions;
}

std::vector<double2> Framework::atomUnitCellLennardJonesPotentialParameters(const ForceField& forceField) const
{
  std::vector<double2> parameters;
  parameters.reserve(unitCellAtoms.size());

  for (const Atom& atom : unitCellAtoms)
  {
    std::size_t type = atom.type;
    double2 parameter = double2(forceField[type].parameters.x, forceField[type].parameters.y);
    parameters.push_back(parameter);
  }

  return parameters;
}

std::string Framework::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", 0, name);

  std::print(stream, "    Box:     {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.ax, simulationBox.cell.bx,
             simulationBox.cell.cx);
  std::print(stream, "             {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.ay, simulationBox.cell.by,
             simulationBox.cell.cy);
  std::print(stream, "             {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.cell.az, simulationBox.cell.bz,
             simulationBox.cell.cz);
  std::print(stream, "    Lengths: {:9.5f} {:9.5f} {:9.5f}\n", simulationBox.lengthA, simulationBox.lengthB,
             simulationBox.lengthC);
  double conv = 180.0 / std::numbers::pi;
  std::print(stream, "    Angles:  {:9.5f} {:9.5f} {:9.5f}\n", conv * simulationBox.angleAlpha,
             conv * simulationBox.angleBeta, conv * simulationBox.angleGamma);
  double3 widths = simulationBox.perpendicularWidths();
  std::print(stream, "    Perpendicular widths:  {:9.5f} {:9.5f} {:9.5f}\n\n", widths.x, widths.y, widths.z);

  std::print(stream, "    number Of Atoms:          {:>12d} [-]\n", unitCellAtoms.size());
  std::print(stream, "    mass:                     {:>12.5f} [amu]\n", mass);

  for (std::size_t i = 0; i != definedAtoms.size(); ++i)
  {
    std::size_t atomType = static_cast<std::size_t>(definedAtoms[i].type);

    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomTypeString,
               definedAtoms[i].position.x, definedAtoms[i].position.y, definedAtoms[i].position.z,
               definedAtoms[i].charge);
  }
  // switch (useChargesFrom)
  //{
  //   case UseChargesFrom::PseudoAtoms:
  //     std::print(stream, "    use charge from: pseudo-atoms\n");
  //     break;
  //   case UseChargesFrom::CIF_File:
  //     std::print(stream, "    use charge from: CIF-file\n");
  //     break;
  //   case UseChargesFrom::ChargeEquilibration:
  //     std::print(stream, "    use charge from: charge-equilibration method\n");
  //     break;
  // }
  std::print(stream, "    net charge:                 {:>12.5f} [e]\n", netCharge);
  std::print(stream, "    smallest charge:            {:>12.5f} [e]\n", smallestCharge);
  std::print(stream, "    largest charge:             {:>12.5f} [e]\n\n", largestCharge);

  if (rigid)
  {
    std::print(stream, "    framework type: Rigid\n");
  }
  else if (mixed)
  {
    std::print(stream, "    framework type: Mixed (Fixed {}, Rigid groups {}, Flexible {})\n", fixedAtomCount,
               rigidGroupCount, flexibleAtomCount);
  }
  else
  {
    std::print(stream, "    framework type: Flexible\n");
  }
  std::print(stream, "    number of bonds: {}\n", intraMolecularPotentials.bonds.size());
  printPotentialSummaries(stream, intraMolecularPotentials.bonds, atoms, forceField);
  std::print(stream, "\n");

  std::print(stream, "    number of bend potentials: {} (out of {})\n", intraMolecularPotentials.bends.size(),
             connectivityTable.findAllBends().size());
  printPotentialSummaries(stream, intraMolecularPotentials.bends, atoms, forceField);
  std::print(stream, "\n");

  std::print(stream, "    number of torsion potentials: {} (out of {})\n", intraMolecularPotentials.torsions.size(),
             connectivityTable.findAllTorsions().size());
  printPotentialSummaries(stream, intraMolecularPotentials.torsions, atoms, forceField);
  std::print(stream, "\n");

  std::print(stream, "    number of improper-torsion potentials: {}\n",
             intraMolecularPotentials.improperTorsions.size());
  printPotentialSummaries(stream, intraMolecularPotentials.improperTorsions, atoms, forceField, true);
  std::print(stream, "\n");

  std::print(stream, "    number of intra-framework Van der Waals potentials: {}\n",
             intraMolecularPotentials.vanDerWaals.size());
  printPotentialSummaries(stream, intraMolecularPotentials.vanDerWaals, atoms, forceField);
  std::print(stream, "\n");

  std::print(stream, "    number of intra-framework Coulomb potentials: {}\n",
             intraMolecularPotentials.coulombs.size());
  printPotentialSummaries(stream, intraMolecularPotentials.coulombs, atoms, forceField);
  std::print(stream, "\n");

  return stream.str();
}

nlohmann::json Framework::jsonStatus() const
{
  nlohmann::json status;
  status["name"] = name;
  status["id"] = 0;
  status["mass"] = mass;

  // TODO I feel that the masses, positions and charges belong in the hdf5.
  status["n_bonds"] = intraMolecularPotentials.bonds.size();

  std::vector<std::string> bondTypes(intraMolecularPotentials.bonds.size());
  for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
  {
    bondTypes[i] = intraMolecularPotentials.bonds[i].print();
  }
  status["bondTypes"] = bondTypes;
  status["n_bends"] = intraMolecularPotentials.bends.size();
  status["n_torsions"] = intraMolecularPotentials.torsions.size();
  status["n_improper_torsions"] = intraMolecularPotentials.improperTorsions.size();
  status["n_intra_vdw"] = intraMolecularPotentials.vanDerWaals.size();
  status["n_intra_coulomb"] = intraMolecularPotentials.coulombs.size();

  return status;
}

// Version 2: the rigid-body reference data moved into the Fragment base class.
static constexpr std::uint64_t frameworkGroupVersionNumber{2};

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const FrameworkGroup& g)
{
  archive << frameworkGroupVersionNumber;
  archive << static_cast<std::uint8_t>(g.type);
  archive << static_cast<const Fragment&>(g);
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, FrameworkGroup& g)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber != frameworkGroupVersionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'FrameworkGroup' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  std::uint8_t typeValue{};
  archive >> typeValue;
  g.type = static_cast<FrameworkGroupType>(typeValue);
  archive >> static_cast<Fragment&>(g);
  return archive;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Framework& c)
{
  archive << c.versionNumber;

  archive << c.simulationBox;
  archive << c.spaceGroupHallNumber;
  archive << c.numberOfUnitCells;

  archive << c.name;
  archive << c.filename;

  archive << c.rigid;

  archive << c.mass;

  archive << c.netCharge;
  archive << c.smallestCharge;
  archive << c.largestCharge;

  archive << c.definedAtoms;
  archive << c.unitCellAtoms;
  archive << c.atoms;

  archive << c.connectivityTable;
  archive << c.intraMolecularPotentials;
  archive << c.intraMolecularImageShifts.bonds;
  archive << c.intraMolecularImageShifts.bends;
  archive << c.intraMolecularImageShifts.torsions;
  archive << c.intraMolecularImageShifts.improperTorsions;
  archive << c.intraMolecularImageShifts.vanDerWaals;
  archive << c.intraMolecularImageShifts.coulombs;

  archive << c.groups;
  archive << c.atomGroupIds;
  archive << c.fixedAtomCount;
  archive << c.rigidGroupAtomCount;
  archive << c.flexibleAtomCount;
  archive << c.rigidGroupCount;
  archive << c.mixed;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Framework& c)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > c.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Framework' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> c.simulationBox;
  archive >> c.spaceGroupHallNumber;
  archive >> c.numberOfUnitCells;

  archive >> c.name;
  archive >> c.filename;

  archive >> c.rigid;

  archive >> c.mass;

  archive >> c.netCharge;
  archive >> c.smallestCharge;
  archive >> c.largestCharge;

  archive >> c.definedAtoms;
  archive >> c.unitCellAtoms;
  archive >> c.atoms;

  if (versionNumber == 1)
  {
    [[maybe_unused]] std::vector<std::size_t> chiralCenters;
    std::vector<BondPotential> bonds;
    archive >> chiralCenters;
    archive >> bonds;
    c.intraMolecularPotentials.bonds = std::move(bonds);
  }
  else
  {
    archive >> c.connectivityTable;
    archive >> c.intraMolecularPotentials;
    if (versionNumber >= 3)
    {
      archive >> c.intraMolecularImageShifts.bonds;
      archive >> c.intraMolecularImageShifts.bends;
      archive >> c.intraMolecularImageShifts.torsions;
      archive >> c.intraMolecularImageShifts.improperTorsions;
      archive >> c.intraMolecularImageShifts.vanDerWaals;
      archive >> c.intraMolecularImageShifts.coulombs;
    }
  }

  if (versionNumber >= 4)
  {
    archive >> c.groups;
    archive >> c.atomGroupIds;
    archive >> c.fixedAtomCount;
    archive >> c.rigidGroupAtomCount;
    archive >> c.flexibleAtomCount;
    archive >> c.rigidGroupCount;
    archive >> c.mixed;
  }
  else
  {
    c.groups.clear();
    c.atomGroupIds.clear();
    c.fixedAtomCount = 0;
    c.rigidGroupAtomCount = 0;
    c.flexibleAtomCount = 0;
    c.rigidGroupCount = 0;
    c.mixed = false;
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Framework: Error in binary restart\n"));
  }
#endif

  return archive;
}

std::string Framework::repr() const
{
  std::ostringstream stream;

  std::print(stream, "Framework {} [{}]\n\n", 0, name);
  std::print(stream, "    number of unit cells: {}x{}x{}\n", numberOfUnitCells.x, numberOfUnitCells.y,
             numberOfUnitCells.z);

  std::print(stream, "    number Of Atoms:  {}\n", unitCellAtoms.size());
  std::print(stream, "    net charge:       {:12.5f} [e]\n", netCharge);
  std::print(stream, "    mass:             {:12.5f} [amu]\n", mass);

  for (std::size_t i = 0; i != definedAtoms.size(); ++i)
  {
    std::size_t atomType = static_cast<std::size_t>(definedAtoms[i].type);

    std::print(stream, "    {:3d}: type {:3d} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomType,
               definedAtoms[i].position.x, definedAtoms[i].position.y, definedAtoms[i].position.z,
               definedAtoms[i].charge);
  }

  std::print(stream, "    number of bonds: {}\n", intraMolecularPotentials.bonds.size());
  for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
  {
    std::print(stream, "        {}", intraMolecularPotentials.bonds[i].print());
  }
  std::print(stream, "\n");

  return stream.str();
}

Framework Framework::makeFAU(const ForceField& forceField, int3 replicate)
{
  std::optional<std::size_t> type_si = forceField.findPseudoAtom("Si");
  if (!type_si.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Si"));
  }

  std::optional<std::size_t> type_o = forceField.findPseudoAtom("O");
  if (!type_o.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "O"));
  }

  std::size_t space_group_hall_number = 526;
  SimulationBox simulation_box = SimulationBox(24.2576, 24.2576, 24.2576);
  std::vector<Atom> fractional_atoms_asymmetric_unit_cell{
      Atom({-0.05392, 0.1253, 0.03589}, 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false, false),
      Atom({0, -0.10623, 0.10623}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom({-0.00323, -0.00323, 0.14066}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom({0.0757, 0.0757, -0.03577}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom({0.07063, 0.07063, 0.32115}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false)};
  std::vector<Atom> fractional_atoms_unit_cell = CIFReader::expandDefinedAtomsToUnitCell(
      simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell);

  return Framework(forceField, "FAU", simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell,
                   fractional_atoms_unit_cell, replicate);
}

Framework Framework::makeITQ29(const ForceField& forceField, int3 replicate)
{
  std::optional<std::size_t> type_si = forceField.findPseudoAtom("Si");
  if (!type_si.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Si"));
  }

  std::optional<std::size_t> type_o = forceField.findPseudoAtom("O");
  if (!type_o.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "O"));
  }

  std::size_t space_group_hall_number = 517;
  SimulationBox simulation_box = SimulationBox(11.8671, 11.8671, 11.8671);
  std::vector<Atom> fractional_atoms_asymmetric_unit_cell{
      Atom({0.3683, 0.1847, 0}, 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false, false),
      Atom({0.5, 0.2179, 0}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom({0.2939, 0.2939, 0}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom({0.3429, 0.1098, 0.1098}, -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false)};
  std::vector<Atom> fractional_atoms_unit_cell = CIFReader::expandDefinedAtomsToUnitCell(
      simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell);

  return Framework(forceField, "ITQ-29", simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell,
                   fractional_atoms_unit_cell, replicate);
}

Framework Framework::makeMFI(const ForceField& forceField, int3 replicate)
{
  std::optional<std::size_t> type_si = forceField.findPseudoAtom("Si");
  if (!type_si.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Si"));
  }

  std::optional<std::size_t> type_o = forceField.findPseudoAtom("O");
  if (!type_o.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "O"));
  }

  std::size_t space_group_hall_number = 292;
  SimulationBox simulation_box = SimulationBox(20.022, 19.899, 13.383);
  std::vector<Atom> fractional_atoms_asymmetric_unit_cell{
      Atom(double3(0.42238, 0.0565, -0.33598), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.30716, 0.02772, -0.1893), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.27911, 0.06127, 0.0312), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.12215, 0.06298, 0.0267), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.07128, 0.02722, -0.18551), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.18641, 0.05896, -0.32818), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.42265, -0.1725, -0.32718), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.30778, -0.13016, -0.18548), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.27554, -0.17279, 0.03109), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.12058, -0.1731, 0.02979), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.07044, -0.13037, -0.182), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.18706, -0.17327, -0.31933), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0, false,
           false),
      Atom(double3(0.3726, 0.0534, -0.2442), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.3084, 0.0587, -0.0789), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.2007, 0.0592, 0.0289), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.0969, 0.0611, -0.0856), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.1149, 0.0541, -0.2763), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.2435, 0.0553, -0.246), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.3742, -0.1561, -0.2372), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.3085, -0.1552, -0.0728), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.198, -0.1554, 0.0288), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.091, -0.1614, -0.0777), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.1169, -0.1578, -0.2694), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.2448, -0.1594, -0.2422), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.3047, -0.051, -0.1866), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.0768, -0.0519, -0.1769), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.4161, 0.1276, -0.3896), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.4086, -0.0017, -0.4136), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.402, -0.1314, -0.4239), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.1886, 0.1298, -0.3836), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.194, 0.0007, -0.4082), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.1951, -0.1291, -0.419), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(-0.0037, 0.0502, -0.208), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(-0.004, -0.1528, -0.2078), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.4192, -0.25, -0.354), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom(double3(0.1884, -0.25, -0.3538), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false,
           false),
      Atom(double3(0.2883, -0.25, 0.0579), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
      Atom(double3(0.1085, -0.25, 0.0611), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0, false, false),
  };
  std::vector<Atom> fractional_atoms_unit_cell = CIFReader::expandDefinedAtomsToUnitCell(
      simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell);

  return Framework(forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), space_group_hall_number,
                   fractional_atoms_asymmetric_unit_cell, fractional_atoms_unit_cell, replicate);
}

Framework Framework::makeCHA(const ForceField& forceField, int3 replicate)
{
  std::optional<std::size_t> type_si = forceField.findPseudoAtom("Si");
  if (!type_si.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Si"));
  }

  std::optional<std::size_t> type_o = forceField.findPseudoAtom("O");
  if (!type_o.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "O"));
  }
  std::size_t space_group_hall_number = 1;
  SimulationBox simulation_box = SimulationBox(9.459, 9.459, 9.459, 94.07 * std::numbers::pi / 180.0,
                                               94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0);
  std::vector<Atom> fractional_atoms_asymmetric_unit_cell{
      Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, static_cast<std::uint16_t>(type_si.value()), 0,
           false, false),
      Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false),
      Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, static_cast<std::uint16_t>(type_o.value()), 0,
           false, false)};
  std::vector<Atom> fractional_atoms_unit_cell = CIFReader::expandDefinedAtomsToUnitCell(
      simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell);

  return Framework(forceField, "CHA", simulation_box, space_group_hall_number, fractional_atoms_asymmetric_unit_cell,
                   fractional_atoms_unit_cell, replicate);
}
