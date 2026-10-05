module;

export module mc_cell_list;

import std;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;

/**
 * \brief Persistent cell list over the molecule atoms of a system for the Monte Carlo energy differences.
 *
 * The box is divided into a grid of cells in fractional coordinates whose perpendicular width is at least
 * the molecule-molecule cut-off, so every atom within the cut-off of a trial position lies in the cell of
 * that position or in one of its 26 neighbours. A trial energy then visits the atoms of 27 cells instead
 * of all N molecule atoms, which is what the brute-force 'computeInterMolecularEnergyDifference' and the
 * CBMC 'computeInterMolecularEnergy' do. Nothing about the 27-cell neighbourhood is stored: the cells are
 * enumerated at query time, so moving an atom from one cell to another is a removal from one bucket and an
 * insertion into another, independent of where the two cells are.
 *
 * Storage ("fixed-capacity buckets"): one array with 'capacity' slots per cell; cell c owns the slots
 * [c*capacity, c*capacity + counts[c]). A slot holds a copy of the position (so that the distance test
 * streams through contiguous 32-byte records), the molecule id (for the same-molecule exclusion) and the
 * index of the atom in the molecule-atom span (to fetch the type, charge and scaling factors of the few
 * atoms that pass the cut-off). Removal is a swap with the last slot of the bucket; both update paths are
 * O(1) per atom through the per-atom back-references 'cellOfAtom'/'slotOfAtom'. A bucket overflow (which
 * the 1.5x margin makes rare) invalidates the list; the owner rebuilds it on the next query.
 *
 * Grid selection: n_i = floor(w_i / cutOff) cells along each box direction (w_i the perpendicular width).
 * A direction with fewer than four cells gives no reduction (the 3-wide stencil would visit every cell),
 * so it is collapsed to a single cell; when all three directions collapse, the list is 'disabled' and the
 * callers use their brute-force loops. For a 75 A box with a 14 A cut-off this gives 5x5x5 cells, a query
 * visits 27/125 of the atoms.
 *
 * Validity: the list is tied to a particular atom array (its size and ordering), box and cut-off; the
 * owner (System) checks 'isCurrent' before every query and rebuilds when any of these changed. Position
 * changes of individual molecules are applied incrementally with 'updateAtoms' by the accepted moves;
 * moves that change positions in other ways invalidate the list.
 */
export struct MCCellList
{
  struct Record
  {
    double x, y, z;  ///< The position (three plain doubles: 'double3' carries a fourth, padding component).
    std::uint32_t atomIndex;
    std::uint32_t moleculeId;

    Record() = default;
    Record(const double3& position, std::uint32_t atomIndex, std::uint32_t moleculeId)
        : x(position.x), y(position.y), z(position.z), atomIndex(atomIndex), moleculeId(moleculeId)
    {
    }
    [[nodiscard]] inline double3 position() const { return double3(x, y, z); }
    inline void setPosition(const double3& position)
    {
      x = position.x;
      y = position.y;
      z = position.z;
    }
  };
  static_assert(sizeof(Record) == 32, "MCCellList::Record should be 32 bytes");

  // ---- state -------------------------------------------------------------------------------------------
  bool valid{false};            ///< The buckets describe 'numberOfAtoms' atoms in 'cellMatrix' at 'cutOff'.
  bool enabled{false};          ///< At least one direction has four or more cells (otherwise brute force).
  int3 numberOfCells{1, 1, 1};  ///< Cells along a, b and c (1 for a collapsed direction).
  std::size_t totalNumberOfCells{1};
  std::size_t capacity{0};  ///< Slots per cell.
  double cutOff{0.0};       ///< Cut-off the grid was sized for (the largest pair cut-off in use).
  double3x3 cellMatrix{};   ///< Box cell matrix at build time (a changed box triggers a rebuild).
  double3x3 inverseCellMatrix{};
  std::size_t numberOfAtoms{0};  ///< Size of the atom span at build time.

  std::vector<std::uint32_t> counts;      ///< Occupancy per cell.
  std::vector<Record> slots;              ///< totalNumberOfCells * capacity records.
  std::vector<std::uint32_t> cellOfAtom;  ///< Back-reference per atom index.
  std::vector<std::uint32_t> slotOfAtom;  ///< Slot within the bucket per atom index.

  // ---- statistics ---------------------------------------------------------------------------------------
  std::size_t numberOfBuilds{0};
  std::size_t numberOfAtomUpdates{0};

  /// The grid the list would use for this box and cut-off ('int3(1,1,1)' when it would be disabled).
  [[nodiscard]] static int3 gridFor(const SimulationBox& simulationBox, double cutOff);

  /// Whether the list would accelerate queries at all for this box and cut-off (any direction >= 4 cells).
  [[nodiscard]] static bool wouldBeEnabled(const SimulationBox& simulationBox, double cutOff);

  /// Rebuilds the buckets for 'atoms' in 'simulationBox' with the grid for 'cutOff'.
  void build(std::span<const Atom> atoms, const SimulationBox& simulationBox, double cutOff);

  /// Marks the list stale; the next query rebuilds it.
  void invalidate() { valid = false; }

  /// Whether the list describes exactly this atom array, box and cut-off.
  [[nodiscard]] bool isCurrent(std::span<const Atom> atoms, const SimulationBox& simulationBox, double cutOff) const;

  /**
   * \brief Applies the new positions of the atoms [first, first + count) of 'atoms' (the same span the list
   * was built from). Atoms that stayed in their cell get their stored position overwritten, the others are
   * moved to their new bucket. A bucket overflow invalidates the list.
   */
  void updateAtoms(std::span<const Atom> atoms, std::size_t first, std::size_t count);

  /// Cell index of a position (any position, also outside the box).
  [[nodiscard]] inline std::size_t cellIndexOf(const double3& position) const
  {
    double3 s = inverseCellMatrix * position;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    const int ix = std::min(static_cast<int>(s.x * static_cast<double>(numberOfCells.x)), numberOfCells.x - 1);
    const int iy = std::min(static_cast<int>(s.y * static_cast<double>(numberOfCells.y)), numberOfCells.y - 1);
    const int iz = std::min(static_cast<int>(s.z * static_cast<double>(numberOfCells.z)), numberOfCells.z - 1);
    return cellIndex(ix, iy, iz);
  }

  [[nodiscard]] inline std::size_t cellIndex(int ix, int iy, int iz) const
  {
    return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(numberOfCells.y) + static_cast<std::size_t>(iy)) *
               static_cast<std::size_t>(numberOfCells.z) +
           static_cast<std::size_t>(iz);
  }

  /**
   * \brief Visits every record in the 27-cell neighbourhood of 'position' (fewer along collapsed directions).
   *
   * The callback receives 'const Record &'. The caller does the distance test (it knows its cut-offs and the
   * minimum-image convention of its box); the neighbourhood is complete for any cut-off <= 'cutOff'.
   */
  template <typename Callback>
  inline void forEachNeighbourRecord(const double3& position, Callback&& callback) const
  {
    double3 s = inverseCellMatrix * position;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    const int cx = std::min(static_cast<int>(s.x * static_cast<double>(numberOfCells.x)), numberOfCells.x - 1);
    const int cy = std::min(static_cast<int>(s.y * static_cast<double>(numberOfCells.y)), numberOfCells.y - 1);
    const int cz = std::min(static_cast<int>(s.z * static_cast<double>(numberOfCells.z)), numberOfCells.z - 1);

    // a collapsed direction (one cell) has itself as its only neighbour
    const int dxMin = numberOfCells.x > 1 ? -1 : 0, dxMax = numberOfCells.x > 1 ? 1 : 0;
    const int dyMin = numberOfCells.y > 1 ? -1 : 0, dyMax = numberOfCells.y > 1 ? 1 : 0;
    const int dzMin = numberOfCells.z > 1 ? -1 : 0, dzMax = numberOfCells.z > 1 ? 1 : 0;

    for (int dx = dxMin; dx <= dxMax; ++dx)
    {
      const int ix = wrap(cx + dx, numberOfCells.x);
      for (int dy = dyMin; dy <= dyMax; ++dy)
      {
        const int iy = wrap(cy + dy, numberOfCells.y);
        for (int dz = dzMin; dz <= dzMax; ++dz)
        {
          const int iz = wrap(cz + dz, numberOfCells.z);
          const std::size_t cell = cellIndex(ix, iy, iz);
          const Record* bucket = slots.data() + cell * capacity;
          const std::uint32_t count = counts[cell];
          for (std::uint32_t k = 0; k != count; ++k)
          {
            callback(bucket[k]);
          }
        }
      }
    }
  }

  /**
   * \brief Visits every unordered pair of records that can lie within 'cutOff' of each other exactly once.
   *
   * The pairs within a cell (slot k < l) and the pairs between a cell and the 13 cells of its positive
   * half-shell (fewer along collapsed directions) are enumerated. The grid never has 2 or 3 cells along a
   * direction (such a direction is collapsed to one cell), so the 26 neighbours of a cell are distinct
   * cells and the half-shell enumerates every neighbouring cell pair once. The callback receives
   * '(const Record &a, const Record &b)' in no particular order; the caller does the distance test and
   * the same-molecule exclusion. This is the full-energy counterpart of 'forEachNeighbourRecord'.
   */
  template <typename Callback>
  inline void forEachPairOnce(Callback&& callback) const
  {
    // the 13 offsets of the positive half-shell: (dx > 0) or (dx == 0, dy > 0) or (dx == dy == 0, dz > 0)
    const int dxMax = numberOfCells.x > 1 ? 1 : 0;
    const int dyMax = numberOfCells.y > 1 ? 1 : 0;
    const int dzMax = numberOfCells.z > 1 ? 1 : 0;

    for (int cx = 0; cx < numberOfCells.x; ++cx)
    {
      for (int cy = 0; cy < numberOfCells.y; ++cy)
      {
        for (int cz = 0; cz < numberOfCells.z; ++cz)
        {
          const std::size_t cell = cellIndex(cx, cy, cz);
          const Record* bucket = slots.data() + cell * capacity;
          const std::uint32_t count = counts[cell];
          if (count == 0) continue;

          // pairs within the cell
          for (std::uint32_t k = 0; k + 1 < count; ++k)
          {
            for (std::uint32_t l = k + 1; l < count; ++l)
            {
              callback(bucket[k], bucket[l]);
            }
          }

          // pairs with the half-shell cells
          for (int dx = 0; dx <= dxMax; ++dx)
          {
            const int ix = wrap(cx + dx, numberOfCells.x);
            const int dyMin = dx > 0 ? -dyMax : 0;
            for (int dy = dyMin; dy <= dyMax; ++dy)
            {
              const int iy = wrap(cy + dy, numberOfCells.y);
              const int dzMin = (dx > 0 || dy > 0) ? -dzMax : 1;
              for (int dz = dzMin; dz <= dzMax; ++dz)
              {
                const int iz = wrap(cz + dz, numberOfCells.z);
                const std::size_t other = cellIndex(ix, iy, iz);
                const Record* otherBucket = slots.data() + other * capacity;
                const std::uint32_t otherCount = counts[other];
                for (std::uint32_t k = 0; k != count; ++k)
                {
                  for (std::uint32_t l = 0; l != otherCount; ++l)
                  {
                    callback(bucket[k], otherBucket[l]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /// Checks the buckets against a fresh assignment of 'atoms' (membership, stored positions, back-references).
  [[nodiscard]] bool verify(std::span<const Atom> atoms, const SimulationBox& simulationBox) const;

  /// Mean number of records a query visits (sum over cells of the occupancy of its neighbourhood / cells).
  [[nodiscard]] double averageNeighbourhoodSize() const;

 private:
  [[nodiscard]] static inline int wrap(int i, int n)
  {
    if (i < 0) return i + n;
    if (i >= n) return i - n;
    return i;
  }

  void removeAtom(std::size_t atomIndex);
  [[nodiscard]] bool insertAtom(const Atom& atom, std::size_t atomIndex);
};
