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
import polymerization_site;
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

/// The state of a polymerization site (see PolymerizationSite).
export enum class PolymerSiteState : std::uint8_t
{
  Monomer = 0,  ///< Unreacted double bond.
  Radical = 1,  ///< Active chain end: the tail carries the unpaired electron.
  Dead = 2      ///< Fully reacted (both the head and the tail are bonded).
};

export enum class PolymerBondKind : std::uint8_t
{
  Propagation = 0,  ///< Radical tail + monomer head.
  Termination = 1   ///< Radical tail + radical tail (combination).
};

/**
 * \brief One bond formed by the polymerization: two atoms of (in general) different molecules.
 */
export struct PolymerBond
{
  PolymerAtomRef a{};  ///< The tail of the attacking radical.
  PolymerAtomRef b{};  ///< The head of the attacked monomer, or the tail of the other radical.
  PolymerBondKind kind{PolymerBondKind::Propagation};
  std::uint64_t formedAtStep{0};  ///< MD step at which the bond was formed (the start of its force ramp).
  double formationDistance{0.0};  ///< Distance of the two atoms when the bond was formed [Angstrom].

  bool operator==(PolymerBond const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerBond &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerBond &b);
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
 * \brief The bonds formed by a chain-growth polymerization and everything their energy needs.
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
 *  - the head and the tail of a reacted site change pseudo-atom type (the sp2 beads of the double bond
 *    become the sp3 beads of the backbone); the bonded terms of the monomer that contain a retyped atom
 *    are overridden by the product parameters (the overlay adds 'new - old').
 *
 * Every term of a bond is ramped in over 'relaxationSteps' MD steps after its formation (a linear
 * factor lambda on the term and on the pair corrections), so that a bond formed at the capture
 * distance pulls the atoms together smoothly instead of releasing the full bond energy at once.
 *
 * The tail of a reacted site is a stereocentre of the backbone (its substituents: the head, the side group of
 * the site definition and the polymer partner; the implicit hydrogen of the united-atom bead is the fourth).
 * 'siteChirality' holds per site the sign of its normalised chiral volume (see ChiralVolume): for a dead site
 * the configuration the centre is held in by the restraint
 *
 *   U = 1/2 k (V - chi V_0)^2        (chi = +-1, k 'chiralForceConstant', V_0 'chiralReferenceVolume'),
 *
 * ramped in with the bond that created the centre; for a radical the configuration its tail has to take when it
 * reacts (set by the driver from the tacticity); zero for a free centre (no tacticity control, or a site without
 * side group), which is neither restrained nor selected.
 *
 * The authoritative state is 'bonds', the site states and chiralities and the atom types of the system (the
 * latter are stored with the atoms); the per-term caches are derived and rebuilt by rebuild() after every
 * reaction and after a restart.
 */
export struct PolymerLinkTable
{
  std::uint64_t versionNumber{2};

  PolymerTemplates templates{};
  std::vector<PolymerBond> bonds{};

  // Site bookkeeping. The sites are enumerated per component, molecule and site definition; the layout is fixed
  // when the table is activated (the start of the reactive stage; the molecule counts do not change afterwards).
  bool activated{false};
  std::vector<std::size_t> sitesPerMolecule{};      ///< Per component.
  std::vector<std::size_t> siteOffset{};            ///< Per component: index of the first site of the component.
  std::vector<std::size_t> moleculesPerComponent{}; ///< Per component: the molecule counts at activation.
  std::vector<std::uint8_t> siteState{};            ///< Per site: PolymerSiteState.
  std::vector<std::uint64_t> siteStep{};            ///< Per site: the step of its last state change (ramp).
  std::vector<std::int8_t> siteChirality{};         ///< Per site: configuration (dead) or required one (radical).

  std::size_t relaxationSteps{0};  ///< Length of the force ramp of a new bond in MD steps (0: no ramp).
  std::uint64_t currentStep{0};    ///< The current MD step (set by the driver; the ramps refer to it).
  double chiralForceConstant{0.0};    ///< Of the restraint on the chiral volume [internal energy units].
  double chiralReferenceVolume{0.6};  ///< Magnitude of the normalised chiral volume the restraint holds.

  std::size_t numberOfInitiations{0};
  std::size_t numberOfPropagations{0};
  std::size_t numberOfTerminations{0};

  // ---- derived caches (rebuild()) ------------------------------------------------------------------------------

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
  /// A pair whose treatment by the pair sums ('counted') differs from what the polymer topology requires
  /// ('wanted'); the overlay adds value(wanted) - value(counted).
  struct PairCorrection
  {
    std::array<std::uint32_t, 2> atoms{};
    PairState counted{PairState::Full};
    PairState wanted{PairState::Excluded};
    PairScalingInfo countedScaling{};
    PairScalingInfo wantedScaling{};
    std::uint32_t bondId{};
  };
  /// A bonded term of a molecule whose parameters change because an atom of it was retyped: the overlay adds
  /// 'replacement - original'. The ramp follows the step of the state change of the site that retyped the atom.
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

  /// The restraint of a reacted stereocentre: 'atoms' are the centre (the tail), the head, the polymer partner and
  /// the side group; 'chirality' the sign held; the ramp is the one of the bond that created the centre.
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
  /// Link partners of the atoms that carry polymer bonds (global atom indices).
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> linkNeighbours{};
  /// Per site: the polymer partner of the tail (global atom index) and the bond that made it; none: no partner.
  std::vector<std::uint32_t> tailPartner{};
  std::vector<std::uint32_t> tailBond{};
  /// Per site: the site whose tail is bonded to the head (the previous unit of the chain), and the site whose
  /// head the tail is bonded to (the next unit); none: chain end or termination junction.
  std::vector<std::size_t> previousSite{};
  std::vector<std::size_t> nextSite{};
  static constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
  static constexpr std::size_t noSite = std::numeric_limits<std::size_t>::max();
  bool derivedValid{false};

  bool operator==(PolymerLinkTable const &other) const
  {
    return templates == other.templates && bonds == other.bonds && siteState == other.siteState &&
           siteChirality == other.siteChirality;
  }

  /// True when the overlay has nothing to evaluate (no bond formed and no site retyped).
  bool empty() const { return bonds.empty() && numberOfInitiations == 0; }
  std::size_t numberOfBonds() const { return bonds.size(); }

  /// Fixes the site layout for the current molecule counts; all sites start as monomers.
  void activate(const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent);

  std::size_t numberOfSites() const { return siteState.size(); }
  std::size_t siteIndex(std::size_t componentId, std::size_t moleculeIndex, std::size_t site) const
  {
    return siteOffset[componentId] + moleculeIndex * sitesPerMolecule[componentId] + site;
  }
  /// The component, molecule and site definition of site index 'index' (atomIndex holds the site definition).
  PolymerAtomRef siteRef(std::size_t index) const;
  PolymerSiteState state(std::size_t index) const { return static_cast<PolymerSiteState>(siteState[index]); }
  std::int8_t chirality(std::size_t index) const { return siteChirality[index]; }
  PolymerAtomRef headOf(std::size_t index, const std::vector<Component> &components) const;
  PolymerAtomRef tailOf(std::size_t index, const std::vector<Component> &components) const;
  /// The side group of the site definition (none: the tail is not a stereocentre).
  std::optional<PolymerAtomRef> sideGroupOf(std::size_t index, const std::vector<Component> &components) const;
  /// The site whose head or tail is the global atom 'atom', with whether the atom is the tail.
  std::optional<std::pair<std::size_t, bool>> siteOfAtom(std::uint32_t atom,
                                                         const std::vector<Component> &components) const;
  std::size_t countSites(PolymerSiteState state) const;
  /// Fraction of the sites that have reacted (radical or dead): the double-bond conversion.
  double conversion() const;

  /// The normalised chiral volume of the tail of 'site' with the global atom 'partner' as its polymer neighbour
  /// (minimum-image geometry); none when the site has no side group.
  std::optional<double> chiralVolumeOf(std::size_t site, std::uint32_t partner, const SimulationBox &simulationBox,
                                       const std::vector<Component> &components, std::span<const Atom> atoms) const;
  /// The sign of chiralVolumeOf(): +1 / -1, or 0 when the site has no side group.
  std::int8_t geometricChirality(std::size_t site, std::uint32_t partner, const SimulationBox &simulationBox,
                                 const std::vector<Component> &components, std::span<const Atom> atoms) const;
  /// The site whose tail is bonded to the head of 'site' through a propagation bond (from the bond list; valid
  /// at any time), or none.
  std::size_t predecessorOf(std::size_t site, const std::vector<Component> &components) const;

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

  /// Monomer -> radical by an initiator: the head and the tail take their product types. The centre the tail
  /// will become is free (no predecessor in the chain).
  void initiate(std::size_t site, const ForceField &forceField, const std::vector<Component> &components,
                std::span<Atom> atoms, std::span<Molecule> moleculeData, std::uint64_t step);
  /// Radical 'radicalSite' adds to the head of monomer 'monomerSite': bond tail(radical) - head(monomer); the
  /// radical site is dead afterwards with configuration 'chiralityRadical' and the monomer site is the new
  /// radical, which has to take 'chiralityMonomer' when it reacts (0: free).
  void propagate(std::size_t radicalSite, std::size_t monomerSite, const ForceField &forceField,
                 const std::vector<Component> &components, std::span<Atom> atoms, std::span<Molecule> moleculeData,
                 std::uint64_t step, double distance, std::int8_t chiralityRadical = 0,
                 std::int8_t chiralityMonomer = 0);
  /// Two radicals combine: bond tail - tail; both sites are dead afterwards with the given configurations.
  void terminate(std::size_t siteA, std::size_t siteB, const ForceField &forceField,
                 const std::vector<Component> &components, std::span<Atom> atoms, std::span<Molecule> moleculeData,
                 std::uint64_t step, double distance, std::int8_t chiralityA = 0, std::int8_t chiralityB = 0);

  /// Rebuilds the layout and the derived term lists for the current bonds, site states and atom types.
  void rebuild(const ForceField &forceField, const std::vector<Component> &components,
               const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> atoms);

  /// The neighbours (component connectivity and polymer bonds) of a global atom index.
  std::vector<std::uint32_t> neighboursOf(std::uint32_t atom, const std::vector<Component> &components) const;

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
    double gelFraction{0.0};            ///< Molecules of the largest cluster over all molecules with sites.

    std::size_t effectiveCycleRank{0};    ///< cycleRank - spanningClusters: independent loops of the unfolded network per box.
    std::size_t spanningClusters{0};      ///< Clusters connected to one of their own periodic images.
    std::size_t spanningMolecules{0};     ///< Molecules in the spanning clusters.
    std::size_t largestSpanningCluster{0};///< Molecules in the largest spanning cluster (0: no percolation).
    std::size_t percolationDimension{0};  ///< Directions (0-3) in which the largest spanning cluster percolates.
    std::array<bool, 3> percolates{false, false, false};  ///< Some cluster winds along the a, b, c axis.
    double spanningFraction{0.0};         ///< Molecules of the spanning clusters over all molecules with sites.
  };
  NetworkStatistics networkStatistics(const SimulationBox &simulationBox, const std::vector<Component> &components,
                                      std::span<const Atom> atoms) const;

  /**
   * \brief Tacticity of the chains from the current geometry (needs the derived caches).
   *
   * A dyad is a pair of consecutive stereocentres of a chain (the tails of two sites joined by a propagation
   * bond, both reacted); it is meso when both have the same sign of the chain-ordered chiral volume and racemo
   * otherwise. The triads (mm, mr, rr) are the pairs of consecutive dyads. 'inverted' counts the restrained
   * centres whose geometry currently has the sign opposite to the one held.
   */
  struct StereoStatistics
  {
    std::size_t stereoCentres{0};  ///< Reacted tails with a side group.
    std::size_t restrained{0};     ///< Stereocentres with a configuration held by the restraint.
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
  /// The bonds as [[c, m, a], [c, m, a]] with their kind, for the topology output.
  nlohmann::json jsonTopology() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerLinkTable &t);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerLinkTable &t);

 private:
  void retype(std::size_t globalAtom, const std::string &typeName, const ForceField &forceField, std::span<Atom> atoms);
  void updateMoleculeMass(const PolymerAtomRef &ref, const ForceField &forceField, std::span<const Atom> atoms,
                          std::span<Molecule> moleculeData) const;
  Layout makeLayout(const std::vector<Component> &components,
                    const std::vector<std::size_t> &numberOfMoleculesPerComponent) const;
};
