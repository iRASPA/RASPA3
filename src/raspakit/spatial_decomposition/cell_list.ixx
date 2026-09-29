module;

export module spatial_decomposition_cell_list;

import std;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import spatial_decomposition_domain_decomposition;

/// Position and charge of an atom in the compact per-sub-domain arrays (array-of-structures for the pair kernel).
export struct LocalAtom
{
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double charge{0.0};
};

/**
 * \brief Cell list plus per-sub-domain Verlet neighbour lists for the spatial-decomposition MD engine.
 *
 * The atoms are binned on a grid of cells (a fixed fraction of `cutoff + skin` wide, in the perpendicular widths
 * of the cell, so triclinic boxes are handled in fractional coordinates) and stored in a structure-of-arrays
 * copy sorted by cell. Every sub-domain of the DomainDecomposition owns a balanced share of the atoms and builds
 * one neighbour list over them in which every pair of the system appears exactly once (each cell evaluates the
 * pairs with the "upper half" of its neighbour stencil, and the pairs inside a cell once with j > i).
 *
 * Each sub-domain works on a *compact local atom array*: first its owned atoms, then the *ghost images* it
 * interacts with. A ghost image is an atom (owned by this or another sub-domain) together with the periodic shift
 * that brings it next to the owned atoms, so the pair kernel needs no minimum-image operation: the image
 * positions are gathered once per step (gatherPositions) with the shift of the current box applied, and the
 * neighbour lists hold local indices. The force on a ghost image is accumulated in the sub-domain's private force
 * buffer; the owner of the atom adds it in a reduction phase after the pair phase (`imagesByOwner` lists, per
 * owner, the image slots that owner has to collect, so every ghost contribution is touched once and no two
 * threads write the same atom).
 *
 * The shift of a pair stays valid until the next rebuild: with `cutoff + skin` at most half the smallest
 * perpendicular width of the box (checked in `setup`) no other image of the neighbour can come within the cutoff
 * before an atom has moved more than `skin / 2`, which is the rebuild criterion (`refreshPositionsAndCheck`).
 */
export struct CellList
{
  struct DomainLists
  {
    std::vector<std::uint32_t> ownedAtoms{};  ///< Sorted-order indices of the owned atoms (local indices 0..).

    // ghost images, local index = ownedAtoms.size() + slot
    std::vector<std::uint32_t> imageAtom{};                   ///< Slot -> sorted-order index of the atom.
    std::vector<std::uint32_t> imageShift{};                  ///< Slot -> packed periodic shift (see shiftVector).
    std::vector<std::uint32_t> imageOwner{};                  ///< Slot -> owning sub-domain of the atom.
    std::vector<std::vector<std::uint32_t>> imagesByOwner{};  ///< Per owner (this one included): its slots.

    std::vector<std::uint32_t> neighbourStart{};  ///< CSR offsets into neighbourList per owned atom.
    std::vector<std::uint32_t> neighbourList{};   ///< Local indices of the neighbours (owned or image).

    // compact per-sub-domain atom data (owned atoms, then images); positions refreshed by gatherPositions
    std::vector<LocalAtom> positions{};
    std::vector<std::uint16_t> localType{};
    std::vector<double> localScalingVDW{}, localScalingCoulomb{};

    // build-time lookup tables (sorted-order index -> local index / first image slot), kept between builds
    std::vector<std::uint32_t> localIndexOfAtom{};
    std::vector<std::uint32_t> imageHead{};
    std::vector<std::uint32_t> imageNext{};

    double maximumDisplacementSquared{0.0};  ///< Result of the last refreshPositionsAndCheck.

    std::size_t numberOfLocalAtoms() const { return ownedAtoms.size() + imageAtom.size(); }
  };

  double cutoff{0.0};      ///< Interaction cutoff (largest of VDW and Coulomb) [Angstrom].
  double skin{0.0};        ///< Verlet skin [Angstrom].
  double listCutoff{0.0};  ///< cutoff + skin.
  int3 numberOfCells{1, 1, 1};
  int subdivision{3};          ///< Cells per list cutoff along an axis (cell width >= listCutoff / subdivision).
  int3 stencilRange{1, 1, 1};  ///< Neighbour-cell offsets needed to cover the list cutoff, per axis.
  double3x3 cellAtBuild{};     ///< Box matrix at the last build.
  DomainDecomposition decomposition{};
  std::size_t numberOfDomains{1};
  std::optional<int3> requestedGrid{};

  // Structure-of-arrays copy of the atoms, in cell-sorted order
  std::size_t numberOfAtoms{0};
  std::vector<std::uint32_t> sortedToOriginal{};
  std::vector<std::uint32_t> originalToSorted{};
  std::vector<double> x{}, y{}, z{};           ///< Current positions (refreshed every step by the engine).
  std::vector<double> refX{}, refY{}, refZ{};  ///< Positions at the last build.
  std::vector<double> charge{}, scalingVDW{}, scalingCoulomb{};
  std::vector<std::uint16_t> type{};
  std::vector<std::uint32_t> moleculeId{};
  std::vector<std::uint32_t> cellOfAtom{};
  std::vector<std::uint32_t> ownerOfAtom{};
  std::vector<double> wrappedX{}, wrappedY{}, wrappedZ{};  ///< Positions at the last build, wrapped into the box.
  std::vector<int3> wrap{};                                ///< Integer box translation removed by the wrapping.
  std::vector<std::uint32_t> cellStart{};     ///< CSR offsets of the cells into the sorted atoms (size cells + 1).
  std::vector<std::uint32_t> stencilStart{};  ///< CSR offsets of the evaluated neighbour cells per cell.

  /// A neighbour cell whose pairs a cell evaluates, with the box translation that places the neighbour cell next
  /// to the cell (the periodic wrap of the stencil offset). When the stencil wraps around a small axis so that the
  /// neighbour is reached with two different offsets, `ambiguous` is set and the pairs use the minimum image.
  struct StencilEntry
  {
    std::uint32_t cell{0};
    std::int8_t wrapX{0}, wrapY{0}, wrapZ{0};
    std::uint8_t ambiguous{0};
  };
  std::vector<StencilEntry> stencilList{};  ///< Neighbour cells whose pairs this cell evaluates (half stencil).

  std::vector<DomainLists> domains{};

  std::size_t numberOfBuilds{0};

  /**
   * \brief Fixes cutoff, skin and thread count and derives the cell grid and the sub-domain grid for the box.
   *
   * Throws when cutoff + skin exceeds half the smallest perpendicular width of the box, or when the threads
   * cannot be laid out on the resulting cell grid.
   */
  void setup(const SimulationBox& box, double cutoff, double skin, std::size_t numberOfThreads,
             std::optional<int3> requestedGrid);

  /// Recomputes the cell grid for the (possibly changed) box; the sub-domain grid is independent of it.
  void updateCellGrid(const SimulationBox& box);

  /**
   * \brief Bins all atoms into cells (counting sort), refreshes the sorted structure-of-arrays copy and rebalances
   * the sub-domains on the atom positions.
   *
   * Serial; O(N log N). Fills the sorted atom data, the cell offsets, the ownership and the per-domain owned-atom
   * lists, and records the box and the reference positions. Must be followed by buildLists for every domain.
   */
  void bin(const SimulationBox& box, std::span<const Atom> atoms);

  /**
   * \brief Builds the neighbour list, the ghost images and the compact atom data of one sub-domain from the current
   * binning (parallel per domain).
   *
   * Cell-pair based: for an owned atom the candidates of a neighbour cell are tested against the wrapped positions
   * with the box translation of the cell pair applied once per cell pair (three subtractions and a dot product per
   * candidate); the periodic shift of an accepted pair follows from that translation and the wrap counts of the
   * two atoms, so the minimum image is only computed for the rare ambiguous stencil entries.
   */
  void buildLists(std::size_t domain, const SimulationBox& box);

  /// Fills the compact position array of a domain from the shared positions, with the periodic shifts of the
  /// current box applied to the images (parallel per domain; after all owners refreshed their positions).
  void gatherPositions(std::size_t domain, const SimulationBox& box);

  /// Copies the current positions of the atoms owned by a domain into the x/y/z arrays and returns whether an
  /// owned atom has moved more than skin / 2 since the last build.
  bool refreshPositionsAndCheck(std::size_t domain, std::span<const Atom> atoms);

  bool boxChanged(const SimulationBox& box) const;

  std::size_t totalPairs() const;
  std::size_t totalImages() const;

  std::string status() const;

  int3 cellGridFor(const SimulationBox& box) const;
  int3 cellGridFor(const SimulationBox& box, int cellsPerCutoff) const;
  std::uint32_t cellIndexOf(const SimulationBox& box, const double3& position) const;

  /// Periodic shift vector (cell * n) of a packed shift code (n_i + 128 in bits 8 i .. 8 i + 7) for the box.
  static double3 shiftVector(const SimulationBox& box, std::uint32_t code);

 private:
  void buildStencil(const SimulationBox& box);
  static double3 wrappedFractional(const SimulationBox& box, const double3& position);
};
