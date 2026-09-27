module;

export module cross_links;

import std;

import archive;
import json;
import atom;
import bond_potential;
import bend_potential;

/**
 * \brief A reactive site of a molecule: an atom that may take part in cross-links.
 *
 * Declared per component in the molecule file as
 *   "ReactiveSites" : [{"Atom" : 5, "Type" : "epoxide", "Valence" : 1}, ...]
 * The 'Type' names the site chemistry and selects the cross-link bond type (see CrossLinkBondType);
 * 'Valence' is the maximum number of cross-links the site may carry at once (default 1).
 */
export struct ReactiveSite
{
  std::uint64_t versionNumber{1};

  std::size_t atom{};        ///< Index of the atom within the molecule.
  std::string siteType{};    ///< Site type name, matched against the cross-link bond types.
  std::size_t valence{1};    ///< Maximum number of simultaneous cross-links of this site.

  bool operator==(ReactiveSite const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ReactiveSite &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ReactiveSite &s);
};

/**
 * \brief The interaction of one kind of cross-link: which two site types it joins and what it costs.
 *
 * A cross-link is an inter-molecular bonded interaction. Its energy is
 *   u_bond(r_AB) + formationEnergy + sum of junction bends (n_A, A, B) and (A, B, n_B)
 *   - (the non-bonded pair interactions of A-B, n_A-B and A-n_B that the inter-molecular sum counted),
 * so a linked pair of atoms and the 1-3 pairs across the link are treated exactly like the
 * corresponding intramolecular pairs of a single molecule. The 'formationEnergy' is a constant
 * added per link (a bond formation energy, negative when bonding is favourable); it drives the
 * bond-count in the formation/scission ensemble and cancels in the bond-swap move.
 *
 * 'captureRadius' bounds the proposal region of the topology moves: a link may only be formed (and
 * only be broken) between sites closer than this distance. It is a proposal parameter, not part of
 * the potential, and enters the acceptance rules only through the neighbour counts.
 */
export struct CrossLinkBondType
{
  std::uint64_t versionNumber{1};

  std::string siteTypeA{};
  std::string siteTypeB{};
  BondPotential bond{};                         ///< Bond potential between the two sites (identifiers unused).
  std::optional<BendPotential> junctionBend{};  ///< Optional bend applied to the 1-3 triples across the link.
  double captureRadius{2.0};                    ///< Proposal radius of the topology moves [Angstrom].
  double formationEnergy{0.0};                  ///< Constant energy per link (internal units).

  bool operator==(CrossLinkBondType const &) const = default;

  bool matches(const std::string &a, const std::string &b) const
  {
    return (siteTypeA == a && siteTypeB == b) || (siteTypeA == b && siteTypeB == a);
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkBondType &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkBondType &b);
};

/**
 * \brief Reference to one reactive site: an atom of a molecule of a component.
 *
 * Molecule indices are the per-component indices used throughout System (they shift when a molecule
 * of the same component is deleted or inserted in front; the table renumbers itself through
 * CrossLinkTable::moleculeDeleted / moleculeInserted).
 */
export struct CrossLinkSite
{
  std::uint32_t componentId{};
  std::uint32_t moleculeIndex{};
  std::uint32_t atomIndex{};

  bool operator==(CrossLinkSite const &) const = default;
  auto operator<=>(CrossLinkSite const &) const = default;

  std::uint64_t moleculeKey() const
  {
    return (static_cast<std::uint64_t>(componentId) << 32) | static_cast<std::uint64_t>(moleculeIndex);
  }
  std::uint64_t siteKey() const
  {
    return (static_cast<std::uint64_t>(componentId) << 56) | (static_cast<std::uint64_t>(moleculeIndex) << 24) |
           static_cast<std::uint64_t>(atomIndex);
  }
  bool sameMolecule(const CrossLinkSite &other) const
  {
    return componentId == other.componentId && moleculeIndex == other.moleculeIndex;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkSite &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkSite &s);
};

/**
 * \brief One cross-link: two sites and the bond type that joins them.
 */
export struct CrossLink
{
  CrossLinkSite a{};
  CrossLinkSite b{};
  std::size_t bondTypeId{};

  bool operator==(CrossLink const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLink &l);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLink &l);
};

/**
 * \brief One cross-link of a molecule that is being regrown by CBMC, seen from that molecule.
 *
 * The partner side is frozen: the partner's site and its intramolecular neighbours are copied as
 * atoms (position in the partner's stored image, charge, scaling), so the regrowth can evaluate the
 * link's terms against a trial conformation of the regrown molecule alone. Built by
 * 'Interactions::makeCrossLinkTethers' for the molecule of a reinsertion move and carried by the
 * CBMC grow context ('CBMC::GrowContext::crossLinkTethers').
 *
 * The regrowth keeps every linked site of the molecule in place (the linked sites are the placed set
 * of the partial regrowth), so of a tether's terms only those with a regrown neighbour of the site
 * -- the junction bends (n_A, A, B) and the 1-3 exclusions (n_A, B) -- vary; the others are constant
 * and cancel between grow and retrace.
 */
export struct CrossLinkTether
{
  std::size_t siteAtom{};                   ///< Local atom index of the site in the regrown molecule.
  const CrossLinkBondType *bondType{};      ///< Bond type of the link (owned by the system's table).
  std::size_t partnerSiteAtom{};            ///< Local atom index of the partner's site in its molecule.
  Atom partnerSite{};                       ///< The partner's site atom.
  std::vector<std::pair<std::size_t, Atom>> partnerNeighbours{};  ///< (local index, atom) of the partner's neighbours.
};

/**
 * \brief The mutable bond table of a system: every cross-link currently present.
 *
 * Owned by System. The molecules themselves keep their (component-defined, immutable) topology;
 * the cross-links are extra inter-molecular bonded interactions between atoms of different molecules,
 * so all component-level caches stay valid and a percolating network needs no giant molecule.
 *
 * 'bondTypes' and 'links' are the authoritative state (serialized); the per-molecule and per-site
 * indices are derived and rebuilt on demand.
 */
export struct CrossLinkTable
{
  std::uint64_t versionNumber{1};

  std::vector<CrossLinkBondType> bondTypes{};
  std::vector<CrossLink> links{};

  // Derived indices (not serialized; rebuilt by rebuildIndex()).
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> linksPerMolecule{};
  std::unordered_map<std::uint64_t, std::uint32_t> linkCountPerSite{};

  bool operator==(CrossLinkTable const &other) const
  {
    return bondTypes == other.bondTypes && links == other.links;
  }

  /// True when cross-linking is configured for this system (at least one bond type).
  bool enabled() const { return !bondTypes.empty(); }
  /// True when no link is present (fast path for every energy routine).
  bool empty() const { return links.empty(); }
  std::size_t numberOfLinks() const { return links.size(); }

  std::optional<std::size_t> findBondType(const std::string &siteTypeA, const std::string &siteTypeB) const;

  std::span<const std::size_t> linksOfMolecule(std::size_t componentId, std::size_t moleculeIndex) const;
  bool moleculeIsLinked(std::size_t componentId, std::size_t moleculeIndex) const;
  std::uint32_t linkCount(const CrossLinkSite &site) const;
  /// True when a link between exactly these two sites exists.
  bool isLinked(const CrossLinkSite &a, const CrossLinkSite &b) const;

  /// Appends a link and updates the indices; returns its id (index in 'links').
  std::size_t addLink(const CrossLink &link);
  /// Removes link 'id' (swap-with-last) and updates the indices.
  void removeLink(std::size_t id);
  void rebuildIndex();

  /// Renumbering after System deleted molecule 'moleculeIndex' of 'componentId' (the deleted molecule
  /// must not carry links: the moves guarantee that). Links to molecules with higher index shift down.
  void moleculeDeleted(std::size_t componentId, std::size_t moleculeIndex);
  /// Renumbering after System inserted a molecule at 'moleculeIndex' of 'componentId'.
  void moleculeInserted(std::size_t componentId, std::size_t moleculeIndex);

  /// Exchanges the links (not the bond types) with another table (parallel tempering).
  void swapTopology(CrossLinkTable &other);

  std::string printStatus() const;
  nlohmann::json jsonStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CrossLinkTable &t);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CrossLinkTable &t);
};
