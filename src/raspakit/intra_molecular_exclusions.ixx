module;

export module intra_molecular_exclusions;

import std;

import archive;
import atom;
import connectivity_table;

/**
 * \brief The non-bonded exclusion topology of a component.
 *
 * The intramolecular non-bonded interactions of a molecule follow the convention of the protein and polymer
 * force fields (and of LAMMPS, GROMACS, OpenMM): every pair of atoms interacts with the regular force-field pair
 * potential (truncated/shifted van der Waals inside the cutoff, Ewald or shifted Coulomb) EXCEPT the pairs in the
 * exclusion set E. E holds the 1-2 (bonded) and 1-3 (bend) pairs and the pairs inside one rigid fragment, whose
 * separation is fixed. The 1-4 pairs are not excluded but scaled by the component's 1-4 scaling factors (the
 * 'scaledPairs' below), which may be zero.
 *
 * For a rigid molecule every pair is excluded. The Ewald Fourier sum counts every intramolecular pair, so the
 * excluded pairs still receive the correction -q_i q_j erf(alpha r)/r (the 'exclusion' term); the non-excluded
 * pairs are evaluated with the real-space Ewald potential like any other pair.
 *
 * The excluded pairs are stored twice: as a sorted list of pairs (i < j) for the loops over the exclusion
 * corrections, and as a CSR table of the partners of every atom for the membership test of the pair loops.
 *
 * A generic pair loop (the cell lists of the spatial-decomposition engine) evaluates every pair of the system at
 * full strength: the same-molecule pairs it has to skip are the excluded pairs and the scaled pairs, which are
 * evaluated separately with their scaling. The union is kept as a second CSR table ('isExcludedFromPairList').
 *
 * The non-excluded pairs themselves are never stored: the memory of this structure is linear in the molecule size
 * (bonds, bends, rigid fragments, 1-4 pairs), and the pair loops that need them enumerate them with
 * 'forEachNonExcludedPair'. A rigid molecule is the exception: its excluded pairs are all pairs and are listed for
 * the Ewald exclusion corrections, which need them anyway.
 */
export struct IntraMolecularExclusions
{
  std::uint64_t versionNumber{2};  ///< Version number for serialization.

  /// A non-excluded pair whose van der Waals or Coulomb interaction is scaled (the 1-4 pairs).
  struct ScaledPair
  {
    std::uint32_t atomA{0};
    std::uint32_t atomB{0};
    double scalingVDW{1.0};
    double scalingCoulomb{1.0};

    bool operator==(const ScaledPair &) const = default;
  };

  std::size_t numberOfAtoms{0};
  bool allExcluded{false};                            ///< Every pair excluded (rigid molecule).
  std::vector<std::array<std::uint32_t, 2>> pairs{};  ///< The excluded pairs, atomA < atomB, sorted.
  std::vector<std::uint32_t> offsets{};               ///< CSR offsets of the excluded partners per atom.
  std::vector<std::uint32_t> partners{};              ///< Excluded partners of every atom (sorted per atom).
  std::vector<ScaledPair> scaledPairs{};              ///< The non-excluded pairs with a scaling other than one.
  std::vector<std::uint32_t> pairListOffsets{};       ///< CSR offsets of the excluded + scaled partners per atom.
  std::vector<std::uint32_t> pairListPartners{};      ///< Excluded + scaled partners of every atom (sorted).

  bool operator==(const IntraMolecularExclusions &) const = default;

  /// Whether the pair (i, j) is excluded from the non-bonded interactions.
  [[nodiscard]] bool isExcluded(std::size_t i, std::size_t j) const
  {
    if (allExcluded) return true;
    if (offsets.empty()) return false;
    const std::span<const std::uint32_t> list = partnersOf(i);
    return std::ranges::binary_search(list, static_cast<std::uint32_t>(j));
  }

  /// Whether the pair (i, j) has to be skipped by a generic full-strength pair loop: excluded, or scaled (the
  /// scaled pairs are evaluated separately with their scaling).
  [[nodiscard]] bool isExcludedFromPairList(std::size_t i, std::size_t j) const
  {
    if (allExcluded) return true;
    if (pairListOffsets.empty()) return false;
    const std::span<const std::uint32_t> list = pairListPartnersOf(i);
    return std::ranges::binary_search(list, static_cast<std::uint32_t>(j));
  }

  /// The excluded partners of atom 'i', sorted.
  [[nodiscard]] std::span<const std::uint32_t> partnersOf(std::size_t i) const
  {
    if (offsets.empty()) return {};
    return std::span<const std::uint32_t>(partners).subspan(offsets[i], offsets[i + 1] - offsets[i]);
  }

  /// The excluded and scaled partners of atom 'i', sorted (the same-molecule partners a generic pair loop skips).
  [[nodiscard]] std::span<const std::uint32_t> pairListPartnersOf(std::size_t i) const
  {
    if (pairListOffsets.empty()) return {};
    return std::span<const std::uint32_t>(pairListPartners)
        .subspan(pairListOffsets[i], pairListOffsets[i + 1] - pairListOffsets[i]);
  }

  /// Sets the scaled pairs (normalised to atomA < atomB and sorted; an excluded pair is rejected) and rebuilds the
  /// pair-list exclusion table (excluded + scaled partners).
  void setScaledPairs(std::vector<ScaledPair> scaled);

  /// The (van der Waals, Coulomb) scaling of the non-excluded pair (i, j): the factors of a scaled pair, (1, 1)
  /// otherwise.
  [[nodiscard]] std::pair<double, double> scalingOf(std::size_t i, std::size_t j) const
  {
    const std::array<std::uint32_t, 2> key{static_cast<std::uint32_t>(std::min(i, j)),
                                           static_cast<std::uint32_t>(std::max(i, j))};
    const auto it = std::ranges::lower_bound(scaledPairs, key, {},
                                             [](const ScaledPair &p) { return std::array<std::uint32_t, 2>{p.atomA, p.atomB}; });
    if (it != scaledPairs.end() && it->atomA == key[0] && it->atomB == key[1]) return {it->scalingVDW, it->scalingCoulomb};
    return {1.0, 1.0};
  }

  /// The number of pairs of the molecule that are not excluded.
  [[nodiscard]] std::size_t numberOfNonExcludedPairs() const
  {
    if (numberOfAtoms < 2) return 0;
    return numberOfAtoms * (numberOfAtoms - 1) / 2 - pairs.size();
  }

  /**
   * \brief Calls 'function(i, j, scalingVDW, scalingCoulomb)' for every non-excluded pair i < j of the molecule.
   *
   * This enumerates the implicit non-bonded pairs: the pairs are not stored (their number grows with the square of
   * the molecule size), they follow from the exclusions and the scaled pairs. The scaled pairs keep their factors
   * (also zero: with Ewald electrostatics a zero-scaled pair still carries the exclusion kernel), the others get
   * (1, 1). The excluded partners and the scaled pairs are consumed in a merge walk over their sorted lists, so the
   * cost is linear in the number of pairs of the molecule.
   */
  template <typename Function>
  void forEachNonExcludedPair(Function &&function) const
  {
    if (allExcluded) return;
    auto scaled = scaledPairs.begin();
    for (std::size_t i = 0; i + 1 < numberOfAtoms; ++i)
    {
      const std::span<const std::uint32_t> excluded = partnersOf(i);
      auto nextExcluded = std::ranges::lower_bound(excluded, static_cast<std::uint32_t>(i + 1));
      while (scaled != scaledPairs.end() && scaled->atomA < i) ++scaled;
      for (std::size_t j = i + 1; j < numberOfAtoms; ++j)
      {
        if (nextExcluded != excluded.end() && *nextExcluded == j)
        {
          ++nextExcluded;
          continue;
        }
        double scalingVDW = 1.0;
        double scalingCoulomb = 1.0;
        if (scaled != scaledPairs.end() && scaled->atomA == i && scaled->atomB == j)
        {
          scalingVDW = scaled->scalingVDW;
          scalingCoulomb = scaled->scalingCoulomb;
          ++scaled;
        }
        function(i, j, scalingVDW, scalingCoulomb);
      }
    }
  }

  /// Every pair excluded (rigid molecule).
  static IntraMolecularExclusions allPairs(std::size_t numberOfAtoms);

  /// Exclusions from an explicit pair list (duplicates and order are normalised).
  static IntraMolecularExclusions fromPairs(std::size_t numberOfAtoms, std::vector<std::array<std::size_t, 2>> excluded);

  /// The 1-2 and 1-3 pairs of a bond graph.
  static std::vector<std::array<std::size_t, 2>> topologicalPairs(const ConnectivityTable &connectivity);

  /// The 1-4 pairs of a bond graph that are not also 1-2 or 1-3 pairs (ring closures), normalised (a < b).
  static std::vector<std::array<std::size_t, 2>> pairs14(const ConnectivityTable &connectivity);

  std::string printStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const IntraMolecularExclusions &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, IntraMolecularExclusions &e);
};

/**
 * \brief Calls 'function(i, j)' for every excluded pair (i < j, indices into 'atoms') of the molecules in 'atoms'.
 *
 * 'atoms' holds whole molecules, each as a contiguous run of atoms with the same componentId and moleculeId (the
 * layout of the system's molecule atoms and of the per-molecule spans of the Monte Carlo moves; several molecules
 * may follow each other, as in the reaction moves). 'components' is indexed by componentId and provides the
 * exclusion topology through its 'intraMolecularPotentials.exclusions' member. A run whose length differs from the number of
 * atoms of its component is a partial molecule and is rejected: the exclusion corrections need the complete molecule.
 */
export template <typename Components, typename Function>
void forEachExcludedPair(const Components &components, std::span<const Atom> atoms, Function &&function)
{
  std::size_t start = 0;
  while (start < atoms.size())
  {
    const std::size_t componentId = static_cast<std::size_t>(atoms[start].componentId);
    const std::uint32_t moleculeId = atoms[start].moleculeId;
    std::size_t end = start + 1;
    while (end < atoms.size() && static_cast<std::size_t>(atoms[end].componentId) == componentId &&
           atoms[end].moleculeId == moleculeId)
    {
      ++end;
    }

    const IntraMolecularExclusions &exclusions = components[componentId].intraMolecularPotentials.exclusions;
    if (end - start != exclusions.numberOfAtoms)
    {
      throw std::runtime_error(std::format(
          "forEachExcludedPair: the span holds {} atoms of molecule {} of component {}, which has {} atoms; the "
          "exclusion corrections need whole molecules\n",
          end - start, moleculeId, componentId, exclusions.numberOfAtoms));
    }
    for (const std::array<std::uint32_t, 2> &pair : exclusions.pairs)
    {
      function(start + pair[0], start + pair[1]);
    }
    start = end;
  }
}
