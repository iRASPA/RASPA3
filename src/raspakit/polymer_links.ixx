module;

export module polymer_links;

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

/**
 * \brief Reference to one atom of the system: an atom of a molecule of a component.
 */
export struct PolymerAtomRef
{
  std::uint32_t componentId{};
  std::uint32_t moleculeIndex{};
  std::uint32_t atomIndex{};

  bool operator==(PolymerAtomRef const &) const = default;
  auto operator<=>(PolymerAtomRef const &) const = default;

  bool sameMolecule(const PolymerAtomRef &other) const
  {
    return componentId == other.componentId && moleculeIndex == other.moleculeIndex;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerAtomRef &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerAtomRef &r);
};

/**
 * \brief One bond formed by a reaction: two atoms of (in general) different molecules.
 */
export struct PolymerBond
{
  PolymerAtomRef a{};
  PolymerAtomRef b{};
  std::uint16_t reaction{0};      ///< Index of the reaction rule that formed it (see ReactionScheme).
  std::uint64_t formedAtStep{0};  ///< MD step at which the bond was formed (the start of its force ramp).
  double formationDistance{0.0};  ///< Distance of the two atoms when the bond was formed [Angstrom].

  bool operator==(PolymerBond const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerBond &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerBond &b);
};

/**
 * \brief One bond of the molecule topology broken by a reaction ('BreakBonds' of a ReactionRule): two atoms of
 * one molecule. The bond and the bends and torsions of the molecule through it are switched off (ramped out)
 * and the non-bonded exclusions of the pairs around it follow the new connectivity.
 */
export struct BrokenBond
{
  PolymerAtomRef a{};
  PolymerAtomRef b{};
  std::uint16_t reaction{0};      ///< Index of the reaction rule that broke it.
  std::uint64_t brokenAtStep{0};  ///< MD step at which the bond was broken (the start of its ramp).

  bool operator==(BrokenBond const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BrokenBond &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BrokenBond &b);
};

/**
 * \brief The change of pseudo-atom type of an atom by a reaction ('Retype' of a product).
 *
 * The log of the retypings gives the types of the atoms at any step; the bonded terms of a molecule that contain
 * a retyped atom are overridden step by step along it (a term whose atoms are retyped twice goes from the
 * parameters of the component to those of the first product types with the ramp of the first reaction, and from
 * there to those of the second with the ramp of the second).
 */
export struct Retyping
{
  PolymerAtomRef atom{};
  std::uint16_t from{0};  ///< The pseudo-atom type before.
  std::uint16_t to{0};    ///< The pseudo-atom type after.
  std::uint64_t step{0};  ///< MD step of the reaction (the start of the ramp of the overridden terms).

  bool operator==(Retyping const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Retyping &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Retyping &r);
};

/**
 * \brief A stereocentre created by a reaction ('Stereo' of a ReactionRule): the centre atom and its three
 * substituents (global atom indices, in the order of the rule), the configuration it is held in (the sign of its
 * normalised chiral volume; 0: free, not restrained), the bond whose ramp it follows, and the group it belongs
 * to (the group of the centre atom; the dyad statistics pair the centres of consecutive groups of a chain).
 */
export struct StereoCentre
{
  std::array<std::uint32_t, 4> atoms{};
  std::int8_t chirality{0};
  std::uint32_t bond{0};
  std::size_t group{0};

  bool operator==(StereoCentre const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const StereoCentre &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, StereoCentre &c);
};

/**
 * \brief The life of one radical (one living chain end): the identity that sits on an active group and is handed
 * from group to group by the reactions.
 *
 * A radical is born when a reaction puts a group into an active state without taking it from another reactant
 * (initiation), moves with the reaction that deactivates one reactant and activates another (propagation, chain
 * transfer) and dies when its group is deactivated without a successor (termination, both radicals of a
 * combination). The record holds the steps and rules of its birth and death, the group it was born on and the
 * group it sits on (or died on), the number of transfers ('propagations': the kinetic chain length of a chain
 * growth) and the number of reactions it took part in.
 */
export struct RadicalRecord
{
  std::uint64_t bornStep{0};
  std::uint16_t bornRule{0};
  std::size_t bornGroup{0};
  std::uint64_t diedStep{0};  ///< 0 (and 'alive'): still living.
  std::uint16_t diedRule{0};
  bool alive{true};
  std::size_t lastGroup{0};  ///< The group it sits on (living) or died on.
  std::size_t propagations{0};
  std::size_t events{0};

  bool operator==(RadicalRecord const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RadicalRecord &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RadicalRecord &r);
};

/**
 * \brief The type-based bonded parameter tables of the polymer ('Bonds', 'Bends', 'Torsions' of the
 * 'Polymerization' block, in the syntax of the molecule files).
 *
 * The parameters of every bonded term that is created or changed by a reaction (the new bond, the
 * bends and torsions across it, and the terms of the monomer that contain a retyped atom) are looked
 * up here by the pseudo-atom types of the atoms after the reaction. A missing entry is an error at
 * the time of the reaction.
 */
export struct PolymerTemplates
{
  std::uint64_t versionNumber{1};

  struct BondTemplate
  {
    std::array<std::string, 2> types{};
    BondType type{BondType::Harmonic};
    std::vector<double> parameters{};
    bool operator==(BondTemplate const &) const = default;
  };
  struct BendTemplate
  {
    std::array<std::string, 3> types{};
    BendType type{BendType::Harmonic};
    std::vector<double> parameters{};
    bool operator==(BendTemplate const &) const = default;
  };
  struct TorsionTemplate
  {
    std::array<std::string, 4> types{};
    TorsionType type{TorsionType::Harmonic};
    std::vector<double> parameters{};
    bool operator==(TorsionTemplate const &) const = default;
  };

  std::vector<BondTemplate> bonds{};
  std::vector<BendTemplate> bends{};
  std::vector<TorsionTemplate> torsions{};

  // Derived (resolve()): lookup by pseudo-atom type indices, both orientations.
  std::map<std::array<std::uint16_t, 2>, std::size_t> bondIndex{};
  std::map<std::array<std::uint16_t, 3>, std::size_t> bendIndex{};
  std::map<std::array<std::uint16_t, 4>, std::size_t> torsionIndex{};
  std::vector<double> bendMinimum{};     ///< Derived: the minimum of every bend template over the angle.
  std::vector<double> torsionMinimum{};  ///< Derived: the minimum of every torsion template over the dihedral.
  bool resolved{false};

  bool operator==(PolymerTemplates const &other) const
  {
    return bonds == other.bonds && bends == other.bends && torsions == other.torsions;
  }
  bool empty() const { return bonds.empty() && bends.empty() && torsions.empty(); }

  /// Reads the 'Bonds', 'Bends' and 'Torsions' arrays of 'value' (type-name based entries only) and appends them.
  void read(const nlohmann::basic_json<nlohmann::raspa_map> &value, const ForceField &forceField);

  /// Resolves the type names against the force field (throws for unknown pseudo-atoms).
  void resolve(const ForceField &forceField);

  std::optional<BondPotential> bond(const std::array<std::uint16_t, 2> &types,
                                    const std::array<std::size_t, 2> &identifiers) const;
  std::optional<BendPotential> bend(const std::array<std::uint16_t, 3> &types,
                                    const std::array<std::size_t, 3> &identifiers) const;
  std::optional<TorsionPotential> torsion(const std::array<std::uint16_t, 4> &types,
                                          const std::array<std::size_t, 4> &identifiers) const;

  /// The energy of the bend / torsion template of the given types at the given (imaged) positions, relative to
  /// the minimum of the template; nullopt when there is no template for the types.
  std::optional<double> bendStrain(const std::array<std::uint16_t, 3> &types, const double3 &a, const double3 &b,
                                   const double3 &c) const;
  std::optional<double> torsionStrain(const std::array<std::uint16_t, 4> &types, const double3 &a, const double3 &b,
                                      const double3 &c, const double3 &d) const;

  std::string printStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerTemplates &t);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerTemplates &t);
};

/**
 * \brief The normalised chiral volume of a centre with three substituents and its gradient.
 *
 * For the centre c with the substituents p, n and z (b_i = r_i - r_c):
 *
 *   V = (b_p x b_n) . b_z / (|b_p| |b_n| |b_z|),     V in [-1, 1],
 *
 * the signed volume of the parallelepiped of the unit bond vectors: zero for a planar centre, +-0.77 for a
 * tetrahedral one, +-0.61 for three substituents at 114 degrees; its sign is the configuration of the centre.
 * For a backbone stereocentre the substituents are ordered along the chain (p: the head, n: the polymer partner)
 * with z the side group, so that two consecutive centres of a chain have the same sign in a meso dyad and
 * opposite signs in a racemo dyad.
 */
export struct ChiralVolume
{
  double value{0.0};
  std::array<double3, 4> gradient{};  ///< d value / d r of the centre, p, n and z (in that order).

  static ChiralVolume compute(const double3 &bp, const double3 &bn, const double3 &bz);
};

/**
 * \brief The bonds formed by the reactions of a polymerization and everything their energy needs.
 *
 * Owned by System. The molecules keep their component-defined topology; the polymer bonds are an
 * overlay of inter-molecular (in the case of a cycle: intramolecular) bonded interactions, like the
 * cross-links of CrossLinkTable but irreversible, with the complete bonded topology of a real bond
 * and with the reacted atoms changing type:
 *
 *  - the bond itself, the bends across it and the torsions across it (the parameters by the product
 *    types from 'templates');
 *  - the pairs within three bonds across the new bond lose the non-bonded interaction the pair sums
 *    count for them (1-2, 1-3: excluded; 1-4: scaled with the 1-4 scaling of the component) -- the pair
 *    kernels are left untouched and the overlay subtracts exactly what they added;
 *  - the atoms of a reacted group change pseudo-atom type as the rule prescribes (the sp2 beads of the
 *    double bond become the sp3 beads of the backbone); the bonded terms of the molecule that contain a
 *    retyped atom are overridden by the product parameters (the overlay adds 'new - old'; 'retypings' logs the
 *    type changes, a term retyped by several reactions is taken from one parameter set to the next with the
 *    ramp of each);
 *  - a bond of the molecule topology broken by a rule ('brokenBonds'; the ring opening of an epoxide) is
 *    overridden to zero together with the bends and torsions of the molecule through it (the overlay adds
 *    '0 - old'), the connectivity the overlay follows loses the bond, and the pairs around it whose
 *    exclusion no longer holds (1-3 pairs that become 1-4 or farther) receive the non-bonded interaction the
 *    molecule's exclusions withhold from them.
 *
 * Every term of a bond is ramped in over 'relaxationSteps' MD steps after its formation (a linear
 * factor lambda on the term and on the pair corrections), so that a bond formed at the capture
 * distance pulls the atoms together smoothly instead of releasing the full bond energy at once. A later
 * reaction at an atom of an earlier bond does not restart the ramps of what it leaves unchanged: the pair
 * corrections and the overrides form chains in the order of the reactions, every link ramped by its own.
 *
 * The groups are the ReactiveGroup's of the molecules, enumerated per component, molecule and group definition
 * once the table is activated; every group has a type and a state of the ReactionScheme ('scheme'). A reaction
 * rule applied to one or two groups (apply()) retypes their atoms, changes their states, forms the bonds and
 * records the stereocentres of the rule. A stereocentre is held in its configuration by the restraint
 *
 *   U = 1/2 k (V - chi V_0)^2        (chi = +-1, k 'chiralForceConstant', V_0 'chiralReferenceVolume'),
 *
 * on its normalised chiral volume (see ChiralVolume), ramped in with the bond that created it; a centre with
 * chi = 0 (no tacticity control) is recorded for the statistics but not restrained. 'groupChirality' holds per
 * group the configuration its next centre has to take (set by the driver from the tacticity through the 'Next'
 * of a stereo entry; 0: free), and 'previousGroup' the group whose centre precedes it along the chain.
 *
 * The radicals (the living chain ends, see RadicalRecord) are numbered in the order of their birth; 'groupRadical'
 * holds per group the radical sitting on it and 'groupLastRadical' the last radical that passed through it, so a
 * dead group can be traced to its chain. apply() hands the radicals over by the state changes of the reactants.
 *
 * The authoritative state is 'bonds', 'brokenBonds', 'retypings', 'centres', the group states, chiralities and
 * radicals and the atom types of the system (the latter are stored with the atoms); the per-term caches are
 * derived and rebuilt by rebuild() after every reaction and after a restart.
 */
export struct PolymerLinkTable
{
  std::uint64_t versionNumber{2};

  PolymerTemplates templates{};
  ReactionScheme scheme{};
  std::vector<PolymerBond> bonds{};
  std::vector<BrokenBond> brokenBonds{};
  std::vector<Retyping> retypings{};  ///< In the order of the reactions.
  std::vector<StereoCentre> centres{};

  // Group bookkeeping. The groups are enumerated per component, molecule and group definition; the layout is
  // fixed when the table is activated (the start of the reactive stage; the molecule counts do not change afterwards).
  bool activated{false};
  std::vector<std::size_t> groupsPerMolecule{};     ///< Per component.
  std::vector<std::size_t> groupOffset{};           ///< Per component: index of the first group of the component.
  std::vector<std::size_t> moleculesPerComponent{}; ///< Per component: the molecule counts at activation.
  std::vector<std::uint8_t> groupState{};           ///< Per group: state index within its type (see 'scheme').
  std::vector<std::uint64_t> groupStep{};           ///< Per group: the step of its last state change (ramp).
  std::vector<std::int8_t> groupChirality{};        ///< Per group: the configuration its next centre has to take.
  std::vector<std::size_t> previousGroup{};         ///< Per group: the group preceding it along the chain (noGroup: none).
  std::vector<std::size_t> reactionCounts{};        ///< Per reaction rule: the number of times it was applied.
  std::vector<std::size_t> groupRadical{};          ///< Per group: the radical sitting on it (noRadical: none).
  std::vector<std::size_t> groupLastRadical{};      ///< Per group: the last radical that passed through it (noRadical: none).
  std::vector<RadicalRecord> radicals{};            ///< The radicals in the order of their birth (the index is the id).

  std::size_t relaxationSteps{0};  ///< Length of the force ramp of a new bond in MD steps (0: no ramp).
  std::uint64_t currentStep{0};    ///< The current MD step (set by the driver; the ramps refer to it).
  double chiralForceConstant{0.0};    ///< Of the restraint on the chiral volume [internal energy units].
  double chiralReferenceVolume{0.6};  ///< Magnitude of the normalised chiral volume the restraint holds.

  // ---- derived caches (activate() / rebuild()) ------------------------------------------------------------------

  /// Resolves (component, molecule, atom) to the index in the molecule-atom span, and molecules to global indices.
  struct Layout
  {
    std::vector<std::size_t> componentOffset{};
    std::vector<std::size_t> atomsPerMolecule{};
    std::vector<std::size_t> moleculeOffset{};
    std::size_t numberOfAtoms{0};
    std::size_t numberOfMolecules{0};

    std::size_t moleculeStart(std::size_t componentId, std::size_t moleculeIndex) const
    {
      return componentOffset[componentId] + moleculeIndex * atomsPerMolecule[componentId];
    }
    std::size_t index(const PolymerAtomRef &ref) const
    {
      return moleculeStart(ref.componentId, ref.moleculeIndex) + ref.atomIndex;
    }
    /// The component, the molecule (per component) and the local atom index of a global atom index.
    PolymerAtomRef ref(std::size_t atom) const;
    /// The global molecule index of a global atom index.
    std::size_t moleculeOf(std::size_t atom) const;
  };
  Layout layout{};
  std::vector<std::uint16_t> groupType{};     ///< Per group: the group type index of the scheme.
  std::vector<std::uint8_t> groupInitial{};   ///< Per group: the initial state index.

  struct LinkBondTerm
  {
    std::array<std::uint32_t, 2> atoms{};
    BondPotential potential{};
    std::uint32_t bondId{};
  };
  struct LinkBendTerm
  {
    std::array<std::uint32_t, 3> atoms{};
    BendPotential potential{};
    std::uint32_t bondId{};
  };
  struct LinkTorsionTerm
  {
    std::array<std::uint32_t, 4> atoms{};
    TorsionPotential potential{};
    std::uint32_t bondId{};
  };
  /// How a pair is treated by the non-bonded sums: at full strength, as a scaled (1-4) pair, or excluded.
  enum class PairState : std::uint8_t
  {
    Full = 0,
    Scaled = 1,
    Excluded = 2
  };
  struct PairScalingInfo
  {
    double scalingVDW{1.0};
    double scalingCoulomb{1.0};
    bool pair14{false};
  };
  /// A pair whose treatment before the bond ('counted': what the pair sums of the components count) differs from
  /// what the polymer topology requires ('wanted'). With the ramp factor lambda of the bond the pair has to be at
  /// lambda value(wanted) + (1 - lambda) value(counted). A 'listed' pair is still supplied by the regular sums
  /// (the pair lists at full strength, or the bonded work of its molecule for a scaled pair): the overlay adds
  /// lambda (value(wanted) - value(counted)). The pairs the topology excludes (1-2 and 1-3 through the new bond)
  /// that the pair lists held at full strength are left out of the lists altogether ('unlisted' below): the
  /// overlay then supplies the whole lambda value(wanted) + (1 - lambda) value(counted). (In the pair lists the
  /// full-strength Lennard-Jones of two bonded atoms at the bond length, 1e8 K/A and more, would be cancelled by
  /// the overlay afterwards; in the single precision of the device pair kernel that cancellation leaves a biased
  /// residual of the order of an ulp of the large term, a net force on the system.)
  struct PairCorrection
  {
    std::array<std::uint32_t, 2> atoms{};
    PairState counted{PairState::Full};
    PairState wanted{PairState::Excluded};
    PairScalingInfo countedScaling{};
    PairScalingInfo wantedScaling{};
    std::uint64_t step{};  ///< The step of the (newest) bond formation or breaking that changed the pair (ramp).
    bool listed{true};
  };
  /// A bonded term of a molecule whose parameters change because an atom of it was retyped, or that is switched
  /// off because it runs through a broken bond (the replacement is then the zero potential: BondType::None,
  /// BendType::Fixed, TorsionType::Fixed): the overlay adds 'replacement - original'. The ramp follows the step
  /// of the state change of the group that retyped the atom (of the reaction that broke the bond).
  struct OverrideBondTerm
  {
    std::array<std::uint32_t, 2> atoms{};
    BondPotential original{};
    BondPotential replacement{};
    std::uint64_t step{};
  };
  struct OverrideBendTerm
  {
    std::array<std::uint32_t, 3> atoms{};
    BendPotential original{};
    BendPotential replacement{};
    std::uint64_t step{};
  };
  struct OverrideTorsionTerm
  {
    std::array<std::uint32_t, 4> atoms{};
    TorsionPotential original{};
    TorsionPotential replacement{};
    std::uint64_t step{};
  };

  /// The restraint of a stereocentre: 'atoms' are the centre and its three substituents; 'chirality' the sign
  /// held; the ramp is the one of the bond that created the centre.
  struct StereoCentreTerm
  {
    std::array<std::uint32_t, 4> atoms{};
    std::int8_t chirality{0};
    std::uint32_t bondId{};
  };

  std::vector<LinkBondTerm> linkBonds{};
  std::vector<LinkBendTerm> linkBends{};
  std::vector<LinkTorsionTerm> linkTorsions{};
  std::vector<PairCorrection> pairCorrections{};
  std::vector<OverrideBondTerm> overrideBonds{};
  std::vector<OverrideBendTerm> overrideBends{};
  std::vector<OverrideTorsionTerm> overrideTorsions{};
  std::vector<StereoCentreTerm> stereoCentres{};
  /// The pairs the pair lists leave out (the pair corrections with listed == false), as a table over the atoms.
  PairListExclusions unlisted{};
  /// Link partners of the atoms that carry polymer bonds (global atom indices).
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> linkNeighbours{};
  /// The broken bonds as canonical pairs of global atom indices -> the step of the breaking.
  std::map<std::array<std::uint32_t, 2>, std::uint64_t> brokenPairs{};
  /// Per group: the group following it along the chain (the inverse of 'previousGroup'; noGroup: none).
  std::vector<std::size_t> nextGroup{};
  static constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
  static constexpr std::size_t noGroup = std::numeric_limits<std::size_t>::max();
  static constexpr std::size_t noRadical = std::numeric_limits<std::size_t>::max();
  bool derivedValid{false};

  bool operator==(PolymerLinkTable const &other) const
  {
    return templates == other.templates && scheme == other.scheme && bonds == other.bonds &&
           brokenBonds == other.brokenBonds && retypings == other.retypings && centres == other.centres &&
           groupState == other.groupState &&
           groupChirality == other.groupChirality && groupRadical == other.groupRadical &&
           radicals == other.radicals;
  }

  /// True when the overlay has nothing to evaluate (no bond formed or broken and no group reacted).
  bool empty() const { return bonds.empty() && brokenBonds.empty() && totalReactions() == 0; }
  std::size_t numberOfBonds() const { return bonds.size(); }
  std::size_t numberOfBrokenBonds() const { return brokenBonds.size(); }
  std::size_t totalReactions() const;

  /// Fixes the group layout for the current molecule counts; all groups start in their initial states.
  void activate(const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent);

  std::size_t numberOfGroups() const { return groupState.size(); }
  std::size_t groupIndex(std::size_t componentId, std::size_t moleculeIndex, std::size_t group) const
  {
    return groupOffset[componentId] + moleculeIndex * groupsPerMolecule[componentId] + group;
  }
  /// The component, molecule and group definition of group index 'index' (atomIndex holds the definition index).
  PolymerAtomRef groupRef(std::size_t index) const;
  const ReactiveGroup &definitionOf(std::size_t index, const std::vector<Component> &components) const;
  std::uint16_t typeOf(std::size_t index) const { return groupType[index]; }
  std::uint8_t state(std::size_t index) const { return groupState[index]; }
  const std::string &stateName(std::size_t index) const { return scheme.stateName(groupType[index], groupState[index]); }
  /// Whether the group has left its initial state.
  bool reacted(std::size_t index) const { return groupState[index] != groupInitial[index]; }
  /// Whether the group is in an active state (a living chain end, see ReactionScheme).
  bool active(std::size_t index) const { return scheme.isActive(groupType[index], groupState[index]); }
  std::int8_t chirality(std::size_t index) const { return groupChirality[index]; }
  /// The radical sitting on group 'index' (noRadical: none).
  std::size_t radicalOf(std::size_t index) const { return groupRadical[index]; }
  /// The last radical that passed through group 'index' (noRadical: none).
  std::size_t lastRadicalOf(std::size_t index) const { return groupLastRadical[index]; }
  /// The number of radicals created so far.
  std::size_t numberOfRadicals() const { return radicals.size(); }
  /// The number of radicals still living (equals the number of active groups).
  std::size_t countLivingRadicals() const;
  /// The atom of role 'role' (index into the roles of the group's type) of group 'index'.
  PolymerAtomRef atomRefOf(std::size_t index, std::uint16_t role, const std::vector<Component> &components) const;
  std::uint32_t atomOf(std::size_t index, std::uint16_t role, const std::vector<Component> &components) const
  {
    return static_cast<std::uint32_t>(layout.index(atomRefOf(index, role, components)));
  }
  /// The group that contains the global atom 'atom' (none when the atom belongs to no group).
  std::optional<std::size_t> groupOfAtom(std::uint32_t atom, const std::vector<Component> &components) const;
  std::size_t countGroups(std::uint16_t type, std::uint8_t state) const;
  std::size_t countGroups(std::uint16_t type) const;
  /// The groups in active states.
  std::size_t countActive() const;
  /// The conversion as the scheme defines it (default: the fraction of the groups no longer in their initial state).
  double conversion() const;
  /// The most recent stereocentre of a group (none: the group has no centre).
  std::optional<std::size_t> centreOf(std::size_t group) const;

  /// The atoms (centre, substituents) of stereo entry 'k' of rule 'rule' applied to 'groups'.
  std::array<std::uint32_t, 4> stereoAtoms(const ReactionRule &rule, std::size_t k, std::span<const std::size_t> groups,
                                           const std::vector<Component> &components) const;
  /// The normalised chiral volume of a centre with the given atoms (minimum-image geometry).
  double chiralVolumeOf(const std::array<std::uint32_t, 4> &atoms, const SimulationBox &simulationBox,
                        std::span<const Atom> systemAtoms) const;
  /// The sign of chiralVolumeOf(): +1 / -1.
  std::int8_t geometricChirality(const std::array<std::uint32_t, 4> &atoms, const SimulationBox &simulationBox,
                                 std::span<const Atom> systemAtoms) const;

  /// The ramp factor of a term created at 'step'.
  double rampFactor(std::uint64_t step) const
  {
    if (relaxationSteps == 0 || currentStep >= step + relaxationSteps) return 1.0;
    if (currentStep <= step) return 0.0;
    return static_cast<double>(currentStep - step) / static_cast<double>(relaxationSteps);
  }
  /// Whether any ramp is still in progress.
  bool rampsActive() const;

  // ---- reactions (the caller rebuilds the derived caches afterwards) --------------------------------------------

  /**
   * \brief Applies reaction rule 'ruleIndex' to 'groups' (one group for a unimolecular rule, two for a bimolecular
   * one, in the order of the reactants of the rule).
   *
   * Retypes the atoms of the products, changes the states of the groups, forms the bonds of the rule (at the
   * current distances of their atoms), breaks the bonds of its 'BreakBonds' (throws when one of them is already
   * broken) and records its stereocentres with the configurations 'chiralities' (one per stereo entry of the rule;
   * 0: free). A stereo entry with 'Next' sets the required configuration of the next group to 'nextChiralities[k]'
   * and makes the centre's group its predecessor.
   *
   * The radicals are handed over by the state changes of the reactants (active before / after the rule): a
   * reactant that becomes active inherits the radical of the one reactant that becomes inactive in the same
   * reaction (propagation, chain transfer), or receives a new radical when there is none (initiation); a
   * reactant that stays active keeps its radical; the radical of a reactant that becomes inactive dies when no
   * other reactant inherits it (termination).
   *
   * Returns the indices of the new bonds and, per reactant, the radical before and after the reaction and
   * whether it was born or died in it.
   */
  struct ApplyResult
  {
    std::vector<std::size_t> bonds{};                                    ///< Indices of the new bonds.
    std::array<std::size_t, 2> radicalBefore{noRadical, noRadical};  ///< Per reactant: the radical on it before.
    std::array<std::size_t, 2> radicalAfter{noRadical, noRadical};   ///< Per reactant: the radical on it after.
    std::array<bool, 2> born{false, false};                              ///< Per reactant: its radical was created here.
    std::array<bool, 2> died{false, false};                              ///< Per reactant: its radical died here.
  };
  ApplyResult apply(std::size_t ruleIndex, std::span<const std::size_t> groups, std::span<const std::int8_t> chiralities,
                    std::span<const std::int8_t> nextChiralities, const ForceField &forceField,
                    const std::vector<Component> &components, std::span<Atom> atoms, std::span<Molecule> moleculeData,
                    const SimulationBox &simulationBox, std::uint64_t step);

  /// Rebuilds the layout and the derived term lists for the current bonds, group states and atom types.
  void rebuild(const ForceField &forceField, const std::vector<Component> &components,
               const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> atoms);

  /// The neighbours of a global atom index in the current topology: the component connectivity without the
  /// broken bonds, plus the polymer bonds.
  std::vector<std::uint32_t> neighboursOf(std::uint32_t atom, const std::vector<Component> &components) const;

  /**
   * \brief The angular strain a bond between the atoms 'a' and 'b' would carry in the current geometry.
   *
   * The sum of the bend and torsion templates through the bond (the terms rebuild() generates for it), each
   * relative to its minimum, with the atoms typed as after the reaction ('retyped': global atom index -> product
   * type; the other atoms keep their current types). The stretch of the bond itself and the non-bonded
   * corrections are not included (the distance is the capture criterion), and a missing template is skipped
   * (the reaction reports it). Internal energy units. Used for the orientation-weighted acceptance of the
   * Polymerization driver: p = p0 exp(-strain / k_B T*).
   */
  double angularStrain(std::uint32_t a, std::uint32_t b, const SimulationBox &simulationBox,
                       const std::vector<Component> &components, std::span<const Atom> atoms,
                       const std::map<std::uint32_t, std::uint16_t> &retyped) const;

  /**
   * \brief Cluster and cycle analysis of the molecule graph (nodes: molecules, edges: polymer bonds).
   *
   * Every edge carries the lattice vector of the periodic image through which the bond is made (from the
   * minimum-image bond vector between the unwrapped molecules). A cycle of the molecule graph whose lattice
   * vectors sum to zero lifts to a closed loop of the unfolded network; a cycle with a non-zero sum is a *winding
   * loop*: in the unfolded (infinite) network it is a chain that connects the cluster to its own periodic image,
   * i.e. the cluster spans the box in that direction. The rank of the winding vectors of a cluster is the number
   * of directions in which it percolates (0 to 3). The number of independent loops per box of an infinite
   * periodic network is bonds - molecules (the '+1' of a finite connected graph does not apply), so the
   * effective cycle rank is the cycle rank of the molecule graph minus one per spanning cluster.
   */
  struct NetworkStatistics
  {
    std::size_t moleculesInNetwork{0};  ///< Molecules carrying at least one polymer bond.
    std::size_t numberOfClusters{0};    ///< Connected components of the bonded molecules.
    std::size_t largestCluster{0};      ///< Molecules in the largest cluster.
    std::size_t intramolecularBonds{0}; ///< Bonds closing a cycle within one molecule (primary cycles).
    std::size_t cycleRank{0};           ///< bonds - moleculesInNetwork + clusters (independent cycles of the molecule graph).
    double gelFraction{0.0};            ///< Molecules of the largest cluster over all molecules with groups.

    std::size_t effectiveCycleRank{0};    ///< cycleRank - spanningClusters: independent loops of the unfolded network per box.
    std::size_t spanningClusters{0};      ///< Clusters connected to one of their own periodic images.
    std::size_t spanningMolecules{0};     ///< Molecules in the spanning clusters.
    std::size_t largestSpanningCluster{0};///< Molecules in the largest spanning cluster (0: no percolation).
    std::size_t percolationDimension{0};  ///< Directions (0-3) in which the largest spanning cluster percolates.
    std::array<bool, 3> percolates{false, false, false};  ///< Some cluster winds along the a, b, c axis.
    double spanningFraction{0.0};         ///< Molecules of the spanning clusters over all molecules with groups.
  };
  NetworkStatistics networkStatistics(const SimulationBox &simulationBox, const std::vector<Component> &components,
                                      std::span<const Atom> atoms) const;

  /**
   * \brief Tacticity of the chains from the current geometry (needs the derived caches).
   *
   * A dyad is a pair of stereocentres of two consecutive groups of a chain ('previousGroup'); it is meso when
   * both have the same sign of the chiral volume and racemo otherwise. The triads (mm, mr, rr) are the pairs of
   * consecutive dyads. 'inverted' counts the restrained centres whose geometry currently has the sign opposite to
   * the one held.
   */
  struct StereoStatistics
  {
    std::size_t stereoCentres{0};  ///< Centres created by the reactions.
    std::size_t restrained{0};     ///< Centres with a configuration held by the restraint.
    std::size_t inverted{0};       ///< Restrained centres currently in the other configuration.
    std::size_t mesoDyads{0};
    std::size_t racemoDyads{0};
    std::size_t mm{0};
    std::size_t mr{0};
    std::size_t rr{0};

    std::size_t dyads() const { return mesoDyads + racemoDyads; }
    double mesoFraction() const
    {
      return dyads() > 0 ? static_cast<double>(mesoDyads) / static_cast<double>(dyads()) : 0.0;
    }
  };
  StereoStatistics stereoStatistics(const SimulationBox &simulationBox, const std::vector<Component> &components,
                                    std::span<const Atom> atoms) const;

  std::string printStatus(const SimulationBox &simulationBox, const std::vector<Component> &components,
                          std::span<const Atom> atoms) const;
  nlohmann::json jsonStatus(const SimulationBox &simulationBox, const std::vector<Component> &components,
                            std::span<const Atom> atoms) const;
  /// The bonds as [[c, m, a], [c, m, a]] with their reaction, for the topology output.
  nlohmann::json jsonTopology() const;
  /// The broken bonds in the same form.
  nlohmann::json jsonBrokenBonds() const;
  /// The group counts per type and state, for the status outputs.
  nlohmann::json jsonGroupCounts() const;
  /// The radicals (see RadicalRecord) as a table, for the topology output: per radical its id, birth and death
  /// (step and rule), the groups it was born on and sits on as [c, m, g] with the state of the latter (which
  /// locates the radical within the group when the scheme distinguishes, e.g., 'radical_tail' and
  /// 'radical_head'), its propagations and events.
  nlohmann::json jsonRadicals() const;
  /// Radical id or null, for the JSON outputs.
  static nlohmann::json jsonRadicalId(std::size_t radical);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerLinkTable &t);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerLinkTable &t);

 private:
  void retype(std::size_t globalAtom, std::uint16_t type, std::uint64_t step, const ForceField &forceField,
              std::span<Atom> atoms);
  void updateMoleculeMass(const PolymerAtomRef &ref, const ForceField &forceField, std::span<const Atom> atoms,
                          std::span<Molecule> moleculeData) const;
  Layout makeLayout(const std::vector<Component> &components,
                    const std::vector<std::size_t> &numberOfMoleculesPerComponent) const;
};
