module;

export module reactive_group;

import std;

import archive;

/**
 * \brief A reactive functional group of a molecule for the 'Polymerization' driver (see ReactionScheme).
 *
 * Declared per component in the molecule file, separately from the equilibrium 'ReactiveSites' of the
 * cross-link moves:
 *
 *   "ReactiveGroups" : [
 *     {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 0,  "tail" : 1,  "side" : 2}},
 *     {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 13, "tail" : 12, "side" : 11}}
 *   ]
 *
 * 'Type' names the kind of group; the reaction rules of the 'Polymerization' block refer to it. 'Atoms' maps
 * the role names of the group to atom indices of the molecule (into 'pseudoAtoms'); the rules address the atoms
 * by 'label.role'. All groups of one type (across the components of a system) must declare the same roles.
 * 'State' is the state the group starts in (default "monomer"); the rules are transitions between states and
 * may retype the atoms of the group (the double bond of a vinyl group becoming the saturated backbone, the
 * primary amine becoming a secondary one). A molecule file carries no product types: the chemistry is in the
 * rules of the system, the molecule file only names the atoms.
 */
export struct ReactiveGroup
{
  std::uint64_t versionNumber{1};

  std::string type{};
  std::vector<std::pair<std::string, std::size_t>> atoms{};  ///< role -> atom index of the molecule (in input order).
  std::string initialState{"monomer"};

  bool operator==(ReactiveGroup const &) const = default;

  /// The atom of role 'role' (none when the group has no such role).
  std::optional<std::size_t> atomOf(std::string_view role) const
  {
    for (const auto &[name, atom] : atoms)
    {
      if (name == role) return atom;
    }
    return std::nullopt;
  }
  /// Whether 'atom' belongs to the group.
  bool contains(std::size_t atom) const
  {
    return std::ranges::any_of(atoms, [&](const auto &entry) { return entry.second == atom; });
  }
  /// The role names, in input order.
  std::vector<std::string> roles() const
  {
    std::vector<std::string> result{};
    for (const auto &[name, atom] : atoms) result.push_back(name);
    return result;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactiveGroup &g)
  {
    archive << g.versionNumber;
    archive << g.type;
    archive << g.atoms;
    archive << g.initialState;
    return archive;
  }
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactiveGroup &g)
  {
    std::uint64_t versionNumber;
    archive >> versionNumber;
    if (versionNumber > g.versionNumber)
    {
      throw std::runtime_error(std::format("Invalid version reading 'ReactiveGroup' at line {} in file {}\n",
                                           std::source_location::current().line(),
                                           std::source_location::current().file_name()));
    }
    archive >> g.type;
    archive >> g.atoms;
    archive >> g.initialState;
    return archive;
  }
};
