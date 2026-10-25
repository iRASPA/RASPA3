module;

module polymer_reactions;

import std;

import archive;
import json;
import forcefield;
import component;
import reactive_group;

// ---- ReactionRule -------------------------------------------------------------------------------------------------

namespace
{
Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RoleRef &r)
{
  archive << r.label;
  archive << r.role;
  return archive;
}
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RoleRef &r)
{
  archive >> r.label;
  archive >> r.role;
  return archive;
}
Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RoleIndex &r)
{
  archive << r.reactant;
  archive << r.role;
  return archive;
}
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RoleIndex &r)
{
  archive >> r.reactant;
  archive >> r.role;
  return archive;
}
}  // namespace

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactionRule &r)
{
  archive << r.name;
  archive << static_cast<std::uint64_t>(r.reactants.size());
  for (const ReactionRule::Reactant &reactant : r.reactants)
  {
    archive << reactant.label;
    archive << reactant.type;
    archive << reactant.state;
  }
  archive << static_cast<std::uint64_t>(r.products.size());
  for (const ReactionRule::Product &product : r.products)
  {
    archive << product.state;
    archive << product.retype;
  }
  archive << r.capture.has_value();
  if (r.capture.has_value())
  {
    archive << r.capture.value()[0];
    archive << r.capture.value()[1];
  }
  archive << r.captureRadius;
  archive << static_cast<std::uint64_t>(r.formBonds.size());
  for (const std::array<RoleRef, 2> &bond : r.formBonds)
  {
    archive << bond[0];
    archive << bond[1];
  }
  archive << static_cast<std::uint64_t>(r.breakBonds.size());
  for (const std::array<RoleRef, 2> &bond : r.breakBonds)
  {
    archive << bond[0];
    archive << bond[1];
  }
  archive << static_cast<std::uint64_t>(r.stereo.size());
  for (const ReactionRule::Stereo &stereo : r.stereo)
  {
    archive << stereo.centre;
    for (const RoleRef &substituent : stereo.substituents) archive << substituent;
    archive << stereo.next;
  }
  archive << r.probability;
  archive << r.rateConstant;
  archive << r.orientationTemperature;
  archive << r.allowIntramolecular;
  archive << r.initialEvents;

  archive << r.reactantType;
  archive << r.reactantState;
  archive << r.productState;
  archive << r.retype[0];
  archive << r.retype[1];
  archive << r.captureAtoms[0];
  archive << r.captureAtoms[1];
  archive << static_cast<std::uint64_t>(r.bondAtoms.size());
  for (const std::array<RoleIndex, 2> &bond : r.bondAtoms)
  {
    archive << bond[0];
    archive << bond[1];
  }
  archive << static_cast<std::uint64_t>(r.breakAtoms.size());
  for (const std::array<RoleIndex, 2> &bond : r.breakAtoms)
  {
    archive << bond[0];
    archive << bond[1];
  }
  archive << static_cast<std::uint64_t>(r.stereoAtoms.size());
  for (const ReactionRule::ResolvedStereo &stereo : r.stereoAtoms)
  {
    archive << stereo.centre;
    for (const RoleIndex &substituent : stereo.substituents) archive << substituent;
    archive << stereo.next;
  }
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactionRule &r)
{
  std::uint64_t count{};
  archive >> r.name;
  archive >> count;
  r.reactants.resize(count);
  for (ReactionRule::Reactant &reactant : r.reactants)
  {
    archive >> reactant.label;
    archive >> reactant.type;
    archive >> reactant.state;
  }
  archive >> count;
  r.products.resize(count);
  for (ReactionRule::Product &product : r.products)
  {
    archive >> product.state;
    archive >> product.retype;
  }
  bool hasCapture{};
  archive >> hasCapture;
  if (hasCapture)
  {
    std::array<RoleRef, 2> capture{};
    archive >> capture[0];
    archive >> capture[1];
    r.capture = capture;
  }
  archive >> r.captureRadius;
  archive >> count;
  r.formBonds.resize(count);
  for (std::array<RoleRef, 2> &bond : r.formBonds)
  {
    archive >> bond[0];
    archive >> bond[1];
  }
  archive >> count;
  r.breakBonds.resize(count);
  for (std::array<RoleRef, 2> &bond : r.breakBonds)
  {
    archive >> bond[0];
    archive >> bond[1];
  }
  archive >> count;
  r.stereo.resize(count);
  for (ReactionRule::Stereo &stereo : r.stereo)
  {
    archive >> stereo.centre;
    for (RoleRef &substituent : stereo.substituents) archive >> substituent;
    archive >> stereo.next;
  }
  archive >> r.probability;
  archive >> r.rateConstant;
  archive >> r.orientationTemperature;
  archive >> r.allowIntramolecular;
  archive >> r.initialEvents;

  archive >> r.reactantType;
  archive >> r.reactantState;
  archive >> r.productState;
  archive >> r.retype[0];
  archive >> r.retype[1];
  archive >> r.captureAtoms[0];
  archive >> r.captureAtoms[1];
  archive >> count;
  r.bondAtoms.resize(count);
  for (std::array<RoleIndex, 2> &bond : r.bondAtoms)
  {
    archive >> bond[0];
    archive >> bond[1];
  }
  archive >> count;
  r.breakAtoms.resize(count);
  for (std::array<RoleIndex, 2> &bond : r.breakAtoms)
  {
    archive >> bond[0];
    archive >> bond[1];
  }
  archive >> count;
  r.stereoAtoms.resize(count);
  for (ReactionRule::ResolvedStereo &stereo : r.stereoAtoms)
  {
    archive >> stereo.centre;
    for (RoleIndex &substituent : stereo.substituents) archive >> substituent;
    archive >> stereo.next;
  }
  return archive;
}

// ---- ReactionScheme -----------------------------------------------------------------------------------------------

std::optional<std::uint16_t> ReactionScheme::typeIndex(std::string_view name) const
{
  for (std::size_t i = 0; i < types.size(); ++i)
  {
    if (types[i].name == name) return static_cast<std::uint16_t>(i);
  }
  return std::nullopt;
}

std::optional<std::uint8_t> ReactionScheme::stateIndex(std::uint16_t type, std::string_view name) const
{
  const std::vector<std::string> &states = types[type].states;
  for (std::size_t i = 0; i < states.size(); ++i)
  {
    if (states[i] == name) return static_cast<std::uint8_t>(i);
  }
  return std::nullopt;
}

std::optional<std::uint16_t> ReactionScheme::roleIndex(std::uint16_t type, std::string_view name) const
{
  const std::vector<std::string> &roles = types[type].roles;
  for (std::size_t i = 0; i < roles.size(); ++i)
  {
    if (roles[i] == name) return static_cast<std::uint16_t>(i);
  }
  return std::nullopt;
}

bool ReactionScheme::orientationWeighted() const
{
  return std::ranges::any_of(rules, [](const ReactionRule &rule) { return rule.orientationTemperature.has_value(); });
}

bool ReactionScheme::hasStereo() const
{
  return std::ranges::any_of(rules, [](const ReactionRule &rule) { return !rule.stereo.empty(); });
}

bool ReactionScheme::hasBreaking() const
{
  return std::ranges::any_of(rules, [](const ReactionRule &rule) { return !rule.breakBonds.empty(); });
}

void ReactionScheme::resolve(const std::vector<Component> &components, const ForceField &forceField)
{
  types.clear();

  // the group types of the components: the roles from the first group of a type, the same set required of the others
  for (const Component &component : components)
  {
    for (const ReactiveGroup &group : component.reactiveGroups)
    {
      std::optional<std::uint16_t> index = typeIndex(group.type);
      if (!index.has_value())
      {
        GroupType type{};
        type.name = group.type;
        type.roles = group.roles();
        types.push_back(std::move(type));
        index = static_cast<std::uint16_t>(types.size() - 1);
      }
      GroupType &type = types[index.value()];
      std::vector<std::string> roles = group.roles();
      std::vector<std::string> expected = type.roles;
      std::ranges::sort(roles);
      std::ranges::sort(expected);
      if (roles != expected)
      {
        std::string have{}, want{};
        for (const std::string &role : roles) have += role + " ";
        for (const std::string &role : expected) want += role + " ";
        throw std::runtime_error(std::format(
            "[Polymerization]: the '{}' group of component '{}' declares the roles {}but the first group of this type "
            "declares {}(all groups of one type must have the same roles)\n",
            group.type, component.name, have, want));
      }
      if (!std::ranges::contains(type.states, group.initialState)) type.states.push_back(group.initialState);
    }
  }
  if (types.empty())
  {
    throw std::runtime_error("[Polymerization]: no component declares 'ReactiveGroups' in its molecule file\n");
  }

  // the states of the rules
  auto addState = [&](std::uint16_t type, const std::string &state)
  {
    if (!std::ranges::contains(types[type].states, state))
    {
      if (types[type].states.size() >= 255)
      {
        throw std::runtime_error(
            std::format("[Polymerization]: group type '{}' has more than 255 states\n", types[type].name));
      }
      types[type].states.push_back(state);
    }
  };
  std::set<std::string> names{};
  for (ReactionRule &rule : rules)
  {
    const std::string where = std::format("[Polymerization]: reaction '{}': ", rule.name);
    if (rule.name.empty()) throw std::runtime_error("[Polymerization]: a reaction needs a 'Name'\n");
    if (!names.insert(rule.name).second)
    {
      throw std::runtime_error(std::format("[Polymerization]: the reaction name '{}' is used twice\n", rule.name));
    }
    if (rule.reactants.empty() || rule.reactants.size() > 2)
    {
      throw std::runtime_error(where + "'Reactants' must list one or two groups\n");
    }
    if (rule.products.size() != rule.reactants.size())
    {
      throw std::runtime_error(where + "internal: the products are not aligned with the reactants\n");
    }
    if (rule.bimolecular() && rule.reactants[0].label == rule.reactants[1].label)
    {
      throw std::runtime_error(where + "the two reactants need different labels\n");
    }
    for (std::size_t i = 0; i < rule.reactants.size(); ++i)
    {
      const ReactionRule::Reactant &reactant = rule.reactants[i];
      const std::optional<std::uint16_t> type = typeIndex(reactant.type);
      if (!type.has_value())
      {
        std::string known{};
        for (const GroupType &t : types) known += std::format("'{}' ", t.name);
        throw std::runtime_error(std::format("{}reactant '{}' has the unknown group type '{}' (the components declare {})\n",
                                             where, reactant.label, reactant.type, known));
      }
      rule.reactantType[i] = type.value();
      addState(type.value(), reactant.state);
      rule.reactantState[i] = stateIndex(type.value(), reactant.state).value();
      const ReactionRule::Product &product = rule.products[i];
      if (product.state.has_value()) addState(type.value(), product.state.value());
      rule.productState[i] = stateIndex(type.value(), product.state.value_or(reactant.state)).value();
      rule.retype[i].clear();
      for (const auto &[role, typeName] : product.retype)
      {
        const std::optional<std::uint16_t> roleId = roleIndex(type.value(), role);
        if (!roleId.has_value())
        {
          throw std::runtime_error(std::format("{}'Retype' of '{}' names the role '{}' which the '{}' groups do not have\n",
                                               where, reactant.label, role, reactant.type));
        }
        const std::optional<std::size_t> pseudoAtom = forceField.findPseudoAtom(typeName);
        if (!pseudoAtom.has_value())
        {
          throw std::runtime_error(std::format("{}'Retype' of '{}' names the unknown pseudo-atom '{}' (define it in the "
                                               "force field)\n",
                                               where, reactant.label, typeName));
        }
        if (std::ranges::any_of(rule.retype[i], [&](const auto &entry) { return entry.first == roleId.value(); }))
        {
          throw std::runtime_error(std::format("{}'Retype' of '{}' lists the role '{}' twice\n", where, reactant.label, role));
        }
        rule.retype[i].emplace_back(roleId.value(), static_cast<std::uint16_t>(pseudoAtom.value()));
      }
    }

    auto resolveRef = [&](const RoleRef &ref, std::string_view what) -> RoleIndex
    {
      const std::optional<std::uint8_t> reactant = rule.reactantOf(ref.label);
      if (!reactant.has_value())
      {
        throw std::runtime_error(std::format("{}{} '{}' refers to the label '{}', which is not a reactant of this reaction\n",
                                             where, what, ref.text(), ref.label));
      }
      const std::optional<std::uint16_t> role = roleIndex(rule.reactantType[reactant.value()], ref.role);
      if (!role.has_value())
      {
        throw std::runtime_error(std::format("{}{} '{}': the '{}' groups have no role '{}'\n", where, what, ref.text(),
                                             rule.reactants[reactant.value()].type, ref.role));
      }
      return RoleIndex{reactant.value(), role.value()};
    };

    if (rule.bimolecular())
    {
      if (!rule.capture.has_value())
      {
        throw std::runtime_error(where + "a bimolecular reaction needs 'Capture' : {\"Atoms\" : [\"A.x\", \"B.y\"], "
                                         "\"Radius\" : r}\n");
      }
      rule.captureAtoms = {resolveRef(rule.capture.value()[0], "'Capture'"),
                           resolveRef(rule.capture.value()[1], "'Capture'")};
      if (rule.captureAtoms[0] == rule.captureAtoms[1])
      {
        throw std::runtime_error(where + "'Capture' names the same atom twice\n");
      }
      if (rule.captureRadius <= 0.0) throw std::runtime_error(where + "'Capture' 'Radius' must be positive\n");
    }
    else
    {
      if (rule.capture.has_value()) throw std::runtime_error(where + "a unimolecular reaction has no 'Capture'\n");
      if (!rule.formBonds.empty()) throw std::runtime_error(where + "a unimolecular reaction can not form bonds\n");
      if (!rule.stereo.empty()) throw std::runtime_error(where + "a unimolecular reaction creates no stereocentre\n");
      if (rule.orientationTemperature.has_value())
      {
        throw std::runtime_error(where + "'OrientationTemperature' needs bonds to evaluate (a bimolecular reaction)\n");
      }
    }
    if (rule.initialEvents > 0 && rule.bimolecular())
    {
      throw std::runtime_error(where + "'InitialEvents' applies to unimolecular reactions only\n");
    }
    rule.bondAtoms.clear();
    for (const std::array<RoleRef, 2> &bond : rule.formBonds)
    {
      const std::array<RoleIndex, 2> atoms{resolveRef(bond[0], "'FormBonds'"), resolveRef(bond[1], "'FormBonds'")};
      if (atoms[0] == atoms[1])
      {
        throw std::runtime_error(std::format("{}'FormBonds' entry [{}, {}] bonds an atom to itself\n", where,
                                             bond[0].text(), bond[1].text()));
      }
      if (atoms[0].reactant == atoms[1].reactant)
      {
        throw std::runtime_error(std::format("{}'FormBonds' entry [{}, {}] lies within one group (the bonds join the two "
                                             "reactants)\n",
                                             where, bond[0].text(), bond[1].text()));
      }
      rule.bondAtoms.push_back(atoms);
    }
    rule.breakAtoms.clear();
    for (const std::array<RoleRef, 2> &bond : rule.breakBonds)
    {
      const std::array<RoleIndex, 2> atoms{resolveRef(bond[0], "'BreakBonds'"), resolveRef(bond[1], "'BreakBonds'")};
      if (atoms[0] == atoms[1])
      {
        throw std::runtime_error(std::format("{}'BreakBonds' entry [{}, {}] names one atom twice\n", where,
                                             bond[0].text(), bond[1].text()));
      }
      if (atoms[0].reactant != atoms[1].reactant)
      {
        throw std::runtime_error(std::format("{}'BreakBonds' entry [{}, {}] spans two groups (only bonds of the molecule "
                                             "topology, within one group, can be broken)\n",
                                             where, bond[0].text(), bond[1].text()));
      }
      if (std::ranges::any_of(rule.bondAtoms, [&](const std::array<RoleIndex, 2> &formed)
                              { return (formed[0] == atoms[0] && formed[1] == atoms[1]) ||
                                       (formed[0] == atoms[1] && formed[1] == atoms[0]); }))
      {
        throw std::runtime_error(std::format("{}'BreakBonds' entry [{}, {}] is also formed by the reaction\n", where,
                                             bond[0].text(), bond[1].text()));
      }
      // the bond has to exist in every component that carries the group type
      const std::uint16_t type = rule.reactantType[atoms[0].reactant];
      const std::string &roleA = types[type].roles[atoms[0].role];
      const std::string &roleB = types[type].roles[atoms[1].role];
      for (const Component &component : components)
      {
        for (const ReactiveGroup &group : component.reactiveGroups)
        {
          if (group.type != types[type].name) continue;
          const std::size_t a = group.atomOf(roleA).value();
          const std::size_t b = group.atomOf(roleB).value();
          if (a >= component.connectivityTable.numberOfBeads || b >= component.connectivityTable.numberOfBeads ||
              !component.connectivityTable[a, b])
          {
            throw std::runtime_error(std::format(
                "{}'BreakBonds' entry [{}, {}]: atoms {} and {} of component '{}' are not bonded in the molecule "
                "topology\n",
                where, bond[0].text(), bond[1].text(), a, b, component.name));
          }
        }
      }
      rule.breakAtoms.push_back(atoms);
    }
    rule.stereoAtoms.clear();
    for (const ReactionRule::Stereo &stereo : rule.stereo)
    {
      if (rule.formBonds.empty())
      {
        throw std::runtime_error(where + "'Stereo' needs 'FormBonds' (the restraint of the centre is ramped in with the "
                                         "bond)\n");
      }
      ReactionRule::ResolvedStereo resolvedStereo{};
      resolvedStereo.centre = resolveRef(stereo.centre, "'Stereo' 'Centre'");
      for (std::size_t k = 0; k < 3; ++k)
      {
        resolvedStereo.substituents[k] = resolveRef(stereo.substituents[k], "'Stereo' 'Substituents'");
      }
      const std::array<RoleIndex, 4> all{resolvedStereo.centre, resolvedStereo.substituents[0],
                                         resolvedStereo.substituents[1], resolvedStereo.substituents[2]};
      for (std::size_t i = 0; i < 4; ++i)
      {
        for (std::size_t j = i + 1; j < 4; ++j)
        {
          if (all[i] == all[j])
          {
            throw std::runtime_error(std::format("{}'Stereo' centre '{}': the centre and its substituents must be four "
                                                 "different atoms\n",
                                                 where, stereo.centre.text()));
          }
        }
      }
      if (stereo.next.has_value())
      {
        const std::optional<std::uint8_t> next = rule.reactantOf(stereo.next.value());
        if (!next.has_value())
        {
          throw std::runtime_error(std::format("{}'Stereo' 'Next' '{}' is not a reactant label\n", where, stereo.next.value()));
        }
        if (next.value() == resolvedStereo.centre.reactant)
        {
          throw std::runtime_error(std::format("{}'Stereo' 'Next' must name the other reactant\n", where));
        }
        resolvedStereo.next = next;
      }
      rule.stereoAtoms.push_back(resolvedStereo);
    }
    if (rule.probability < 0.0 || rule.probability > 1.0)
    {
      throw std::runtime_error(where + "'Probability' must lie in [0, 1]\n");
    }
    if (rule.rateConstant.has_value() && rule.rateConstant.value() < 0.0)
    {
      throw std::runtime_error(where + "'RateConstant' must be non-negative\n");
    }
    if (rule.orientationTemperature.has_value() && rule.orientationTemperature.value() <= 0.0)
    {
      throw std::runtime_error(where + "'OrientationTemperature' must be positive\n");
    }
  }

  // the active states: reactant of a bimolecular rule and not an initial state of any group
  for (GroupType &type : types) type.active.assign(type.states.size(), 1);
  for (const Component &component : components)
  {
    for (const ReactiveGroup &group : component.reactiveGroups)
    {
      const std::uint16_t type = typeIndex(group.type).value();
      types[type].active[stateIndex(type, group.initialState).value()] = 0;
    }
  }
  for (GroupType &type : types)
  {
    for (std::size_t s = 0; s < type.states.size(); ++s)
    {
      if (type.active[s] == 0) continue;
      const bool reactant = std::ranges::any_of(rules,
                                                [&](const ReactionRule &rule)
                                                {
                                                  if (!rule.bimolecular()) return false;
                                                  for (std::size_t i = 0; i < 2; ++i)
                                                  {
                                                    if (types[rule.reactantType[i]].name == type.name &&
                                                        rule.reactantState[i] == s)
                                                      return true;
                                                  }
                                                  return false;
                                                });
      type.active[s] = reactant ? 1 : 0;
    }
  }

  // the conversion definition
  for (Conversion &entry : conversion)
  {
    const std::optional<std::uint16_t> type = typeIndex(entry.type);
    if (!type.has_value())
    {
      throw std::runtime_error(std::format("[Polymerization]: 'Conversion' names the unknown group type '{}'\n", entry.type));
    }
    entry.typeIndex = type.value();
    entry.stateMask.assign(types[type.value()].states.size(), 0);
    for (const std::string &state : entry.states)
    {
      const std::optional<std::uint8_t> index = stateIndex(type.value(), state);
      if (!index.has_value())
      {
        throw std::runtime_error(std::format("[Polymerization]: 'Conversion' names the state '{}' which no reaction or "
                                             "group of type '{}' uses\n",
                                             state, entry.type));
      }
      entry.stateMask[index.value()] = 1;
    }
  }
  resolved = true;
}

std::string ReactionScheme::describe(std::size_t r, double scanTime) const
{
  const ReactionRule &rule = rules[r];
  std::string text{};
  for (std::size_t i = 0; i < rule.reactants.size(); ++i)
  {
    if (i > 0) text += " + ";
    text += std::format("{}[{}:{}]", rule.reactants[i].label, rule.reactants[i].type, rule.reactants[i].state);
  }
  text += " -> ";
  for (std::size_t i = 0; i < rule.reactants.size(); ++i)
  {
    if (i > 0) text += " + ";
    text += std::format("{}:{}", rule.reactants[i].label, rule.products[i].state.value_or(rule.reactants[i].state));
    if (!rule.products[i].retype.empty())
    {
      text += "{";
      for (std::size_t k = 0; k < rule.products[i].retype.size(); ++k)
      {
        if (k > 0) text += ", ";
        text += std::format("{} -> {}", rule.products[i].retype[k].first, rule.products[i].retype[k].second);
      }
      text += "}";
    }
  }
  if (rule.bimolecular())
  {
    text += std::format("; capture {} - {} within {:g} A", rule.capture.value()[0].text(), rule.capture.value()[1].text(),
                        rule.captureRadius);
    for (const std::array<RoleRef, 2> &bond : rule.formBonds)
    {
      text += std::format(", bond {} - {}", bond[0].text(), bond[1].text());
    }
    text += rule.allowIntramolecular ? ", intramolecular allowed" : ", intramolecular forbidden";
  }
  for (const std::array<RoleRef, 2> &bond : rule.breakBonds)
  {
    text += std::format("; breaks {} - {}", bond[0].text(), bond[1].text());
  }
  text += std::format("; p = {:g} per scan", rule.probabilityPerScan(scanTime));
  if (rule.rateConstant.has_value()) text += std::format(" (k = {:g} 1/ps)", rule.rateConstant.value());
  if (rule.orientationTemperature.has_value()) text += std::format(", T* = {:g} K", rule.orientationTemperature.value());
  if (rule.initialEvents > 0) text += std::format(", {} initial events", rule.initialEvents);
  for (const ReactionRule::Stereo &stereo : rule.stereo)
  {
    text += std::format("; stereocentre {} ({}, {}, {})", stereo.centre.text(), stereo.substituents[0].text(),
                        stereo.substituents[1].text(), stereo.substituents[2].text());
    if (stereo.next.has_value()) text += std::format(" -> next {}", stereo.next.value());
  }
  return text;
}

std::string ReactionScheme::printStatus(double scanTime) const
{
  std::ostringstream stream;
  for (const GroupType &type : types)
  {
    std::string roles{}, states{};
    for (const std::string &role : type.roles) roles += (roles.empty() ? "" : ", ") + role;
    for (std::size_t s = 0; s < type.states.size(); ++s)
    {
      states += std::format("{}{}{}", s > 0 ? ", " : "", type.states[s], type.active[s] ? "*" : "");
    }
    std::print(stream, "    group type '{}':{}roles {}; states {} (* active)\n", type.name,
               std::string(std::max<std::size_t>(1, 20 - std::min<std::size_t>(20, type.name.size())), ' '), roles,
               states);
  }
  for (std::size_t r = 0; r < rules.size(); ++r)
  {
    std::print(stream, "    reaction '{}':{}{}\n", rules[r].name,
               std::string(std::max<std::size_t>(1, 22 - std::min<std::size_t>(22, rules[r].name.size())), ' '),
               describe(r, scanTime));
  }
  if (conversion.empty())
  {
    std::print(stream, "    conversion:                      fraction of the groups no longer in their initial state\n");
  }
  else
  {
    std::string text{};
    for (const Conversion &entry : conversion)
    {
      std::string states{};
      for (const std::string &state : entry.states) states += (states.empty() ? "" : ", ") + state;
      text += std::format("{}'{}' groups in the states {}", text.empty() ? "" : "; ", entry.type, states);
    }
    std::print(stream, "    conversion:                      fraction of the {}\n", text);
  }
  return stream.str();
}

nlohmann::json ReactionScheme::json() const
{
  nlohmann::json result{};
  nlohmann::json typesJson = nlohmann::json::array();
  for (const GroupType &type : types)
  {
    nlohmann::json active = nlohmann::json::array();
    for (std::size_t s = 0; s < type.states.size(); ++s)
    {
      if (type.active[s]) active.push_back(type.states[s]);
    }
    typesJson.push_back({{"name", type.name}, {"roles", type.roles}, {"states", type.states}, {"activeStates", active}});
  }
  result["groupTypes"] = typesJson;
  nlohmann::json rulesJson = nlohmann::json::array();
  for (const ReactionRule &rule : rules)
  {
    nlohmann::json entry{};
    entry["name"] = rule.name;
    nlohmann::json reactants = nlohmann::json::array();
    nlohmann::json products = nlohmann::json::array();
    for (std::size_t i = 0; i < rule.reactants.size(); ++i)
    {
      reactants.push_back(
          {{"label", rule.reactants[i].label}, {"type", rule.reactants[i].type}, {"state", rule.reactants[i].state}});
      nlohmann::json retype = nlohmann::json::object();
      for (const auto &[role, typeName] : rule.products[i].retype) retype[role] = typeName;
      products.push_back({{"label", rule.reactants[i].label},
                          {"state", rule.products[i].state.value_or(rule.reactants[i].state)},
                          {"retype", retype}});
    }
    entry["reactants"] = reactants;
    entry["products"] = products;
    if (rule.capture.has_value())
    {
      entry["capture"] = {{"atoms", {rule.capture.value()[0].text(), rule.capture.value()[1].text()}},
                          {"radius", rule.captureRadius}};
    }
    nlohmann::json bonds = nlohmann::json::array();
    for (const std::array<RoleRef, 2> &bond : rule.formBonds) bonds.push_back({bond[0].text(), bond[1].text()});
    entry["formBonds"] = bonds;
    nlohmann::json broken = nlohmann::json::array();
    for (const std::array<RoleRef, 2> &bond : rule.breakBonds) broken.push_back({bond[0].text(), bond[1].text()});
    entry["breakBonds"] = broken;
    nlohmann::json stereo = nlohmann::json::array();
    for (const ReactionRule::Stereo &s : rule.stereo)
    {
      nlohmann::json item{{"centre", s.centre.text()},
                          {"substituents", {s.substituents[0].text(), s.substituents[1].text(), s.substituents[2].text()}}};
      if (s.next.has_value()) item["next"] = s.next.value();
      stereo.push_back(item);
    }
    entry["stereo"] = stereo;
    entry["probability"] = rule.probability;
    if (rule.rateConstant.has_value()) entry["rateConstant"] = rule.rateConstant.value();
    if (rule.orientationTemperature.has_value()) entry["orientationTemperature"] = rule.orientationTemperature.value();
    entry["allowIntramolecular"] = rule.allowIntramolecular;
    entry["initialEvents"] = rule.initialEvents;
    rulesJson.push_back(entry);
  }
  result["reactions"] = rulesJson;
  nlohmann::json conversionJson = nlohmann::json::array();
  for (const Conversion &entry : conversion) conversionJson.push_back({{"type", entry.type}, {"states", entry.states}});
  result["conversion"] = conversionJson;
  return result;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactionScheme &s)
{
  archive << s.versionNumber;
  archive << s.rules;
  archive << static_cast<std::uint64_t>(s.conversion.size());
  for (const ReactionScheme::Conversion &entry : s.conversion)
  {
    archive << entry.type;
    archive << entry.states;
    archive << entry.typeIndex;
    archive << entry.stateMask;
  }
  archive << static_cast<std::uint64_t>(s.types.size());
  for (const ReactionScheme::GroupType &type : s.types)
  {
    archive << type.name;
    archive << type.roles;
    archive << type.states;
    archive << type.active;
  }
  archive << s.resolved;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactionScheme &s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    throw std::runtime_error(std::format("Invalid version reading 'ReactionScheme' at line {} in file {}\n",
                                         std::source_location::current().line(),
                                         std::source_location::current().file_name()));
  }
  archive >> s.rules;
  std::uint64_t count{};
  archive >> count;
  s.conversion.resize(count);
  for (ReactionScheme::Conversion &entry : s.conversion)
  {
    archive >> entry.type;
    archive >> entry.states;
    archive >> entry.typeIndex;
    archive >> entry.stateMask;
  }
  archive >> count;
  s.types.resize(count);
  for (ReactionScheme::GroupType &type : s.types)
  {
    archive >> type.name;
    archive >> type.roles;
    archive >> type.states;
    archive >> type.active;
  }
  archive >> s.resolved;
  return archive;
}
