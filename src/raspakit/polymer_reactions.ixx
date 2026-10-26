module;

export module polymer_reactions;

import std;

import archive;
import json;
import forcefield;
import component;
import reactive_group;

/**
 * \brief A reference to an atom of a reactant of a reaction rule by its role: 'label.role' in the input.
 */
export struct RoleRef
{
  std::string label{};
  std::string role{};

  bool operator==(RoleRef const &) const = default;
  std::string text() const { return std::format("{}.{}", label, role); }
};

/// A resolved RoleRef: the reactant (0 or 1) and the index of the role in the role list of its group type.
export struct RoleIndex
{
  std::uint8_t reactant{};
  std::uint16_t role{};

  bool operator==(RoleIndex const &) const = default;
};

/**
 * \brief One reaction rule of the polymerization: a transition of one group (unimolecular) or of a pair of
 * groups (bimolecular) between their states, with the atoms retyped, the bonds formed, and the stereocentres
 * created.
 *
 * In the input (the 'Reactions' list of the 'Polymerization' block):
 *
 *   {
 *     "Name" : "propagation",
 *     "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical"},
 *                    {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
 *     "Capture"   : {"Atoms" : ["R.tail", "M.head"], "Radius" : 4.0},
 *     "FormBonds" : [["R.tail", "M.head"]],
 *     "Products"  : [{"Label" : "R", "State" : "dead"},
 *                    {"Label" : "M", "State" : "radical", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
 *     "Stereo"    : [{"Centre" : "R.tail", "Substituents" : ["R.head", "M.head", "R.side"], "Next" : "M"}],
 *     "Probability" : 0.5, "AllowIntramolecular" : true, "OrientationTemperature" : 2000.0
 *   }
 *
 * A bimolecular rule is tried for every pair of groups in the reactant states whose capture atoms are within
 * the capture radius (closest pairs first, one reaction per group per scan), a unimolecular rule for every group
 * in its reactant state; each is accepted with 'Probability' per scan (or 1 - exp(-k dt_scan) from
 * 'RateConstant'), times the orientation weight exp(-dU / k_B T*) of the bonds it forms when
 * 'OrientationTemperature' is given. 'InitialEvents' of a unimolecular rule are applied to random eligible
 * groups when the reactions are activated (the initial radicals). The products give, per reactant, the new
 * state (default: unchanged) and the atoms of the group that change pseudo-atom type.
 *
 * 'BreakBonds' lists bonds of the molecule topology, within one group, that the reaction breaks (the ring
 * opening of an epoxide: "BreakBonds" : [["E.carbon", "E.oxygen"]]); the bond and the bends and torsions
 * through it are switched off and the non-bonded exclusions follow the new connectivity. The atoms stay in
 * their molecule; no atom is removed, so a transferred proton has to be an implicit hydrogen of a united-atom
 * bead (the retyping of the bead carries it).
 *
 * A stereo entry declares a stereocentre created by the reaction: the centre atom and its three substituents
 * (the fourth is the implicit hydrogen of a united-atom bead), whose configuration is the sign of the normalised
 * chiral volume (b1 x b2) . b3 of the bond vectors to the substituents in the listed order (see ChiralVolume).
 * For a chain the order (previous unit, next unit, side group) makes consecutive centres comparable: equal
 * signs are a meso dyad, opposite signs a racemo dyad. 'Next' names the reactant whose own centre is the next
 * one along the chain: it inherits the configuration it has to take from this centre and the tacticity
 * ('Tacticity' of the settings), and the dyad statistics pair the two.
 */
export struct ReactionRule
{
  struct Reactant
  {
    std::string label{};
    std::string type{};
    std::string state{};
    bool operator==(Reactant const &) const = default;
  };
  struct Product
  {
    std::optional<std::string> state{};                           ///< None: the state is unchanged.
    std::vector<std::pair<std::string, std::string>> retype{};  ///< role -> pseudo-atom name.
    bool operator==(Product const &) const = default;
  };
  struct Stereo
  {
    RoleRef centre{};
    std::array<RoleRef, 3> substituents{};
    std::optional<std::string> next{};  ///< Label of the reactant carrying the next centre of the chain.
    // declared here and defaulted out of line (polymer_reactions.cpp): clang 23 crashes generating the implicit
    // destructor (the std::array<RoleRef, 3> member) when it is first needed in an importing module unit. A
    // user-declared destructor keeps the struct an aggregate (designated initializers still work).
    ~Stereo();
    bool operator==(Stereo const &) const = default;
  };

  std::string name{};
  std::vector<Reactant> reactants{};  ///< One (unimolecular) or two (bimolecular).
  std::vector<Product> products{};    ///< Per reactant (same size as 'reactants').
  std::optional<std::array<RoleRef, 2>> capture{};  ///< The pair of atoms whose distance is the capture criterion.
  double captureRadius{4.0};                        ///< [Angstrom]
  std::vector<std::array<RoleRef, 2>> formBonds{};
  std::vector<std::array<RoleRef, 2>> breakBonds{};  ///< Bonds of the molecule topology (within one group) broken.
  std::vector<Stereo> stereo{};
  double probability{1.0};
  std::optional<double> rateConstant{};            ///< [1/ps]
  std::optional<double> orientationTemperature{};  ///< [K]
  bool allowIntramolecular{true};
  std::size_t initialEvents{0};

  // ---- resolved by ReactionScheme::resolve() -------------------------------------------------------------------
  std::array<std::uint16_t, 2> reactantType{};
  std::array<std::uint8_t, 2> reactantState{};
  std::array<std::uint8_t, 2> productState{};
  std::array<std::vector<std::pair<std::uint16_t, std::uint16_t>>, 2> retype{};  ///< Per reactant: role -> type.
  std::array<RoleIndex, 2> captureAtoms{};
  std::vector<std::array<RoleIndex, 2>> bondAtoms{};
  std::vector<std::array<RoleIndex, 2>> breakAtoms{};
  struct ResolvedStereo
  {
    RoleIndex centre{};
    std::array<RoleIndex, 3> substituents{};
    std::optional<std::uint8_t> next{};
  };
  std::vector<ResolvedStereo> stereoAtoms{};

  /// Equality of the input (the resolved tables are derived).
  bool operator==(ReactionRule const &other) const
  {
    return name == other.name && reactants == other.reactants && products == other.products &&
           capture == other.capture && captureRadius == other.captureRadius && formBonds == other.formBonds &&
           breakBonds == other.breakBonds && stereo == other.stereo && probability == other.probability &&
           rateConstant == other.rateConstant &&
           orientationTemperature == other.orientationTemperature &&
           allowIntramolecular == other.allowIntramolecular && initialEvents == other.initialEvents;
  }

  bool bimolecular() const { return reactants.size() == 2; }
  /// Both reactants have the same type and state (the pair is unordered; the roles are assigned at random).
  bool symmetric() const
  {
    return bimolecular() && reactants[0].type == reactants[1].type && reactants[0].state == reactants[1].state;
  }
  /// The index of the reactant with label 'label'.
  std::optional<std::uint8_t> reactantOf(std::string_view label) const
  {
    for (std::size_t i = 0; i < reactants.size(); ++i)
    {
      if (reactants[i].label == label) return static_cast<std::uint8_t>(i);
    }
    return std::nullopt;
  }
  /// The acceptance probability per scan: the probability, or 1 - exp(-k dt) from the rate constant.
  double probabilityPerScan(double scanTime) const
  {
    if (rateConstant.has_value()) return 1.0 - std::exp(-rateConstant.value() * scanTime);
    return probability;
  }
  /// The orientation weight exp(-strain / k_B T*) of a candidate with angular strain 'strainInKelvin' (1 without T*).
  double orientationWeight(double strainInKelvin) const
  {
    if (!orientationTemperature.has_value()) return 1.0;
    return std::exp(-std::max(0.0, strainInKelvin) / orientationTemperature.value());
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactionRule &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactionRule &r);
};

/**
 * \brief The reaction scheme of a system: the group types of its components with their roles and states, the
 * reaction rules between the states, and the definition of the conversion.
 *
 * The group types and their roles come from the 'ReactiveGroups' of the components, the states from the initial
 * states of the groups and from the rules. A state of a type is *active* when some bimolecular rule has it as a
 * reactant and no group starts in it: the groups in active states are the living chain ends (the radicals of a
 * chain growth, the partly reacted amines of a step growth) whose waiting times the PropertyRadicalWaitingTimes
 * collects.
 *
 * The conversion is, by default, the fraction of the groups that are no longer in their initial state;
 * 'Conversion' : [{"Type" : "acrylate", "States" : ["radical", "dead"]}] counts the groups of the listed types
 * in the listed states instead (over the groups of those types).
 */
export struct ReactionScheme
{
  std::uint64_t versionNumber{1};

  struct GroupType
  {
    std::string name{};
    std::vector<std::string> roles{};
    std::vector<std::string> states{};
    std::vector<std::uint8_t> active{};  ///< Per state.
    bool operator==(GroupType const &) const = default;
  };
  struct Conversion
  {
    std::string type{};
    std::vector<std::string> states{};
    std::uint16_t typeIndex{};               ///< Resolved.
    std::vector<std::uint8_t> stateMask{};  ///< Resolved, per state of the type.
    bool operator==(Conversion const &) const = default;
  };

  std::vector<ReactionRule> rules{};
  std::vector<Conversion> conversion{};
  std::vector<GroupType> types{};
  bool resolved{false};

  bool operator==(ReactionScheme const &other) const { return rules == other.rules && conversion == other.conversion; }
  bool empty() const { return rules.empty(); }

  std::optional<std::uint16_t> typeIndex(std::string_view name) const;
  std::optional<std::uint8_t> stateIndex(std::uint16_t type, std::string_view name) const;
  std::optional<std::uint16_t> roleIndex(std::uint16_t type, std::string_view name) const;
  const std::string &stateName(std::uint16_t type, std::uint8_t state) const { return types[type].states[state]; }
  bool isActive(std::uint16_t type, std::uint8_t state) const { return types[type].active[state] != 0; }
  /// Whether any rule carries an orientation temperature.
  bool orientationWeighted() const;
  /// Whether any rule declares a stereocentre.
  bool hasStereo() const;

  /**
   * \brief Collects the group types of the components and resolves the rules against them and the force field.
   *
   * Throws with a message naming the entry for: an unknown group type or role, inconsistent roles of one type
   * across the components, an unknown pseudo-atom in a 'Retype', a bimolecular rule without 'Capture', a
   * unimolecular rule with 'Capture', 'FormBonds' or 'Stereo', a bond or a stereocentre on coinciding atoms,
   * a 'Stereo' in a rule without 'FormBonds' (the ramp of the restraint follows the bond), and a 'BreakBonds'
   * entry that spans two groups, is also formed by the rule, or is not a bond of the molecule topology of every
   * component carrying the group type.
   */
  void resolve(const std::vector<Component> &components, const ForceField &forceField);
  /// Whether any rule breaks bonds.
  bool hasBreaking() const;

  /// A readable description of rule 'r' for the output file.
  std::string describe(std::size_t r, double scanTime) const;
  std::string printStatus(double scanTime) const;
  nlohmann::json json() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactionScheme &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactionScheme &s);
};
