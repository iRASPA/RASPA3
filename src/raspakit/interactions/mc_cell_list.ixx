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
 * The box is divided into a grid of cells in fractional coordinates. The cells are a fraction of the
 * molecule-molecule cut-off wide ('subdivision' cells per cut-off along each axis, the geometry of the
 * spatial-decomposition CellList of the MD engine), so that the cells visited for a query hug the sphere of
 * the query radius instead of the 27 cut-off-sized cells of a classic cell list: a query of radius r_c visits
 * the (2 n + 1)^3 cells with n = ceil(r_c / w), about (7/3 r_c)^3 = 12.7 r_c^3 of space for w = r_c / 3
 * against 27 r_c^3, and the corner cells beyond the radius are skipped in an orthorhombic box.
 *
 * The queries are sphere gathers ('forEachRecordWithin'): every record whose nearest periodic image to a
 * centre lies within a radius is visited exactly once, together with that image position. The periodic
 * shift of a visited cell follows from its stencil offset, so no minimum-image operation is needed per
 * record; only when the stencil wraps around a short box axis (fewer cells along the axis than the stencil
 * is wide) are the records of that axis placed with the minimum image. This is what the gathered
 * neighbourhoods of the Monte Carlo trial energies (Interactions::Neighbourhood) are built from, and what
 * the full-energy pair enumeration ('forEachPairOnce') is the counterpart of.
 *
 * Storage ("fixed-capacity buckets", structure of arrays): 'capacity' slots per cell in every array; cell c
 * owns the slots [c*capacity, c*capacity + counts[c]). A slot holds a copy of everything a pair energy needs
 * -- the position wrapped into the box (so that the per-cell shift applies and the distance tests stream
 * through contiguous arrays), the charge, type, scaling factors, group id, the molecule id (for the
 * same-molecule exclusion) and the index of the atom in the atom span. A gather into an
 * Interactions::Neighbourhood is therefore a filtered block copy per cell ('forEachCellWithin'), without
 * touching the Atom array. Removal is a swap with the last slot of the bucket; both update paths are O(1)
 * per atom through the per-atom back-references 'cellOfAtom'/'slotOfAtom'. A bucket overflow (which the
 * 1.5x margin makes rare) invalidates the list; the owner rebuilds it on the next query.
 *
 * Grid selection: n_i = floor(subdivision w_i / cutOff) cells along each box direction (w_i the
 * perpendicular width). A direction along which a query of the cut-off radius visits every cell anyway
 * (n_i <= 2 ceil(cutOff n_i / w_i) + 1) gives no reduction and is collapsed to a single cell; when all
 * three directions collapse the list is 'disabled' and the callers gather from all atoms. With
 * subdivision 3 a direction is divided when it holds at least 8 cells, i.e. is at least 8/3 cut-offs wide.
 * The subdivision (1 to 3 cells per cut-off) is chosen per build from the atom density, see 'gridFor'.
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
    double x, y, z;  ///< The position wrapped into the box (three plain doubles: 'double3' carries a fourth).
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

  /// The finest subdivision considered: cells of a third of the cut-off (a 7^3 stencil encloses the cut-off
  /// sphere about three times tighter than cells of the full cut-off).
  static constexpr int maximumSubdivision = 3;
  /// Below this many atoms per cell the per-cell overhead of a query (stencil walk, shift, bucket set-up)
  /// outweighs the tighter stencil, and a coarser subdivision is chosen (see 'gridFor').
  static constexpr double minimumAtomsPerCell = 4.0;
  /// Cells per cut-off along an axis of the current grid (the cells are at least cutOff / subdivision wide).
  int subdivision{maximumSubdivision};

  // ---- state -------------------------------------------------------------------------------------------
  bool valid{false};            ///< The buckets describe 'numberOfAtoms' atoms in 'cellMatrix' at 'cutOff'.
  bool enabled{false};          ///< At least one direction is divided (otherwise the callers gather all atoms).
  int3 numberOfCells{1, 1, 1};  ///< Cells along a, b and c (1 for a collapsed direction).
  std::size_t totalNumberOfCells{1};
  std::size_t capacity{0};  ///< Slots per cell.
  double cutOff{0.0};       ///< Cut-off the grid was sized for (the largest pair cut-off in use).
  double3x3 cellMatrix{};   ///< Box cell matrix at build time (a changed box triggers a rebuild).
  double3x3 inverseCellMatrix{};
  double3 perpendicularWidths{};  ///< Perpendicular widths of the box at build time.
  bool orthorhombic{false};       ///< Diagonal cell matrix: the per-axis cell distances combine as a norm.
  std::size_t numberOfAtoms{0};   ///< Size of the atom span at build time.

  std::vector<std::uint32_t> counts;  ///< Occupancy per cell.
  /// Atoms with a scaling factor below one or a group id other than zero (a gather of a list without any copies
  /// the constants instead of the per-slot values).
  std::size_t numberOfPartiallyCoupledAtoms{0};
  // ---- the slots (totalNumberOfCells * capacity entries each) ----------------------------------------
  std::vector<double> slotX, slotY, slotZ;  ///< The position wrapped into the box.
  std::vector<double> slotCharge, slotScalingVDW, slotScalingCoulomb;
  std::vector<double> slotMoleculeIdValue;  ///< The molecule id as a double (exact), for the lane compare.
  std::vector<std::uint32_t> slotType, slotMoleculeId, slotAtomIndex;
  std::vector<std::uint8_t> slotGroupId;
  std::vector<std::uint32_t> cellOfAtom;  ///< Back-reference per atom index.
  std::vector<std::uint32_t> slotOfAtom;  ///< Slot within the bucket per atom index.

  /// The position, atom index and molecule id of a slot as a record.
  [[nodiscard]] inline Record record(std::size_t slot) const
  {
    return Record(double3(slotX[slot], slotY[slot], slotZ[slot]), slotAtomIndex[slot], slotMoleculeId[slot]);
  }

  // ---- statistics ---------------------------------------------------------------------------------------
  std::size_t numberOfBuilds{0};
  std::size_t numberOfAtomUpdates{0};

  /// The grid of 'subdivision' cells per cut-off for this box ('int3(1,1,1)' when it would be disabled).
  [[nodiscard]] static int3 gridForSubdivision(const SimulationBox& simulationBox, double cutOff, int subdivision);

  /// The grid the list would use for this box, cut-off and atom count, and its subdivision: the finest
  /// subdivision whose cells hold 'minimumAtomsPerCell' atoms on average (the per-cell overhead of a query is
  /// paid per stencil cell, the per-atom work per candidate: a dilute system is better served by fewer, larger
  /// cells), or else the coarsest subdivision that still divides a direction. 'int3(1,1,1)' when none does.
  [[nodiscard]] static std::pair<int3, int> gridFor(const SimulationBox& simulationBox, double cutOff,
                                                    std::size_t numberOfAtoms);

  /// Whether the list would accelerate queries at all for this box and cut-off (any direction divided).
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

  /// The position wrapped into the box (fractional coordinates in [0, 1)).
  [[nodiscard]] inline double3 wrappedPosition(const double3& position) const
  {
    double3 s = inverseCellMatrix * position;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    return cellMatrix * s;
  }

  /// Cell index of a position (any position, also outside the box).
  [[nodiscard]] inline std::size_t cellIndexOf(const double3& position) const
  {
    double3 s = inverseCellMatrix * position;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    return cellIndex(cellCoordinate(s.x, numberOfCells.x), cellCoordinate(s.y, numberOfCells.y),
                     cellCoordinate(s.z, numberOfCells.z));
  }

  [[nodiscard]] inline std::size_t cellIndex(int ix, int iy, int iz) const
  {
    return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(numberOfCells.y) + static_cast<std::size_t>(iy)) *
               static_cast<std::size_t>(numberOfCells.z) +
           static_cast<std::size_t>(iz);
  }

  /// Number of cells on either side of a cell that cover 'radius' along axis 'axis'.
  [[nodiscard]] inline int stencilRange(std::size_t axis, double radius) const
  {
    const int n = numberOfCells[static_cast<int>(axis)];
    if (n <= 1) return 0;
    const double cellWidth = perpendicularWidths[axis] / static_cast<double>(n);
    return static_cast<int>(std::ceil(radius / cellWidth));
  }

  /**
   * \brief Visits the periodic images of a point that lie within 'radius' of a centre, given the displacement
   * 'd' of the image nearest to the centre (the minimum image, |d| at most half the box in every direction).
   *
   * The callback receives the displacement of each image from the centre. Without 'allImages' only the nearest
   * image is visited (when within the radius). With 'allImages' every further lattice image within the radius
   * is visited as well -- a radius beyond half a box width reaches more than one image of a point (a radius of
   * a box length or more, as for a batch spanning the box, reaches images on both sides and further), and a
   * gather that hands out all of them lets a kernel skip the minimum image per pair as long as its own cut-off
   * stays within half the box (at most one image per pair then). Along axis a the image shifted by k lattice
   * vectors lies at a perpendicular distance |s_a + k| w_a from the centre (s the fractional displacement, w the
   * perpendicular 'widths'), so only |s_a + k| w_a <= radius needs testing.
   */
  template <typename Callback>
  static inline void forEachImageWithin(const SimulationBox& simulationBox, const double3& widths, const double3& d,
                                        double radius, bool allImages, Callback&& callback)
  {
    const double radiusSquared = radius * radius;
    if (double3::dot(d, d) <= radiusSquared) callback(d);
    if (!allImages) return;

    const double3 s = simulationBox.inverseCell * d;
    int low[3], high[3];
    bool any = false;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const double reach = radius / widths[axis];
      low[axis] = static_cast<int>(std::ceil(-reach - s[axis]));
      high[axis] = static_cast<int>(std::floor(reach - s[axis]));
      any = any || low[axis] < 0 || high[axis] > 0;
    }
    if (!any) return;

    for (int kx = low[0]; kx <= high[0]; ++kx)
    {
      for (int ky = low[1]; ky <= high[1]; ++ky)
      {
        for (int kz = low[2]; kz <= high[2]; ++kz)
        {
          if (kx == 0 && ky == 0 && kz == 0) continue;
          const double3 shift(static_cast<double>(kx), static_cast<double>(ky), static_cast<double>(kz));
          const double3 image = d + simulationBox.cell * shift;
          if (double3::dot(image, image) <= radiusSquared) callback(image);
        }
      }
    }
  }

  /**
   * \brief Visits every record whose nearest periodic image to 'centre' lies within 'radius' exactly once (or,
   * with 'allImages', every periodic image of every record within the radius, one call per image).
   *
   * The callback receives '(const Record &record, const double3 &image)' with 'image' the position of that
   * image (so 'image - centre' needs no minimum-image operation, and for any point p within 'radius - r_c' of
   * the centre, 'image - p' is the minimum image of every pair within r_c as long as 'radius' is below half the
   * smallest perpendicular width of the box -- or, with 'allImages', as long as r_c is). The cells visited are
   * those within 'radius' of the centre; along an axis the stencil wraps around (or that is collapsed) the
   * records are placed with the minimum image to the centre (see 'forEachImageWithin' for the additional
   * images). Any radius is allowed, also beyond the list's cut-off. A disabled or invalid list holds no records
   * (the callers gather from the atom array instead).
   */
  template <typename Callback>
  inline void forEachRecordWithin(const SimulationBox& simulationBox, const double3& centre, double radius,
                                  Callback&& callback, bool allImages = false) const
  {
    const std::array<bool, 3> whole = wholeAxes(radius);
    const bool ambiguous = whole[0] || whole[1] || whole[2];
    const double radiusSquared = radius * radius;
    forEachCellWithin(centre, radius,
                      [&](std::size_t cell, const double3& shift)
                      {
                        const std::size_t base = cell * capacity;
                        const std::uint32_t count = counts[cell];
                        if (ambiguous)
                        {
                          for (std::uint32_t k = 0; k != count; ++k)
                          {
                            const Record entry = record(base + k);
                            const double3 d = simulationBox.applyPeriodicBoundaryConditions(entry.position() - centre);
                            forEachImageWithin(simulationBox, perpendicularWidths, d, radius, allImages,
                                               [&](const double3& image) { callback(entry, centre + image); });
                          }
                          return;
                        }
                        for (std::uint32_t k = 0; k != count; ++k)
                        {
                          const Record entry = record(base + k);
                          const double3 image = entry.position() + shift;
                          const double3 d = image - centre;
                          if (double3::dot(d, d) <= radiusSquared) callback(entry, image);
                        }
                      });
  }

  /// The axes along which a query of 'radius' visits every cell (the stencil is as wide as the axis, or the
  /// axis is collapsed): the per-cell shift of 'forEachCellWithin' is zero along them and the records of a
  /// visited cell need the minimum image along them.
  [[nodiscard]] inline std::array<bool, 3> wholeAxes(double radius) const
  {
    std::array<bool, 3> whole{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const int n = numberOfCells[static_cast<int>(axis)];
      whole[axis] = (n <= 1) || (2 * stencilRange(axis, radius) + 1 >= n);
    }
    return whole;
  }

  /**
   * \brief Visits every non-empty cell that can hold a record within 'radius' of 'centre': the cells of the
   * stencil around the centre's cell whose lower-bound distance to the centre is within the radius.
   *
   * The callback receives '(std::size_t cell, const double3 &shift)': the records of the cell are the images
   * 'position + shift' along the axes the stencil does not wrap ('wholeAxes'); along a whole axis the shift is
   * zero and the record's minimum image along that axis has to be taken by the caller. This is the block
   * primitive the neighbourhood gathers are built on; 'forEachRecordWithin' is the per-record form.
   */
  template <typename Callback>
  inline void forEachCellWithin(const double3& centre, double radius, Callback&& callback) const
  {
    if (!valid || !enabled) return;
    const double radiusSquared = radius * radius;

    // the centre in fractional coordinates: its box (integer part) and its cell
    const double3 s = inverseCellMatrix * centre;
    const double3 box(std::floor(s.x), std::floor(s.y), std::floor(s.z));
    const double3 sw = s - box;
    int centreCell[3], ranges[3];
    bool whole[3];  // every cell of the axis is visited: the stencil is as wide as the axis (or wider)
    double centreCoordinate[3], cellWidth[3];
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const int n = numberOfCells[static_cast<int>(axis)];
      centreCell[axis] = cellCoordinate(sw[axis], n);
      ranges[axis] = stencilRange(axis, radius);
      whole[axis] = (n <= 1) || (2 * ranges[axis] + 1 >= n);
      cellWidth[axis] = perpendicularWidths[axis] / static_cast<double>(n);
      centreCoordinate[axis] = sw[axis] * perpendicularWidths[axis];
    }

    // along an axis the cells visited are either all of them, or the stencil around the centre cell
    const int beginX = whole[0] ? 0 : -ranges[0], endX = whole[0] ? numberOfCells.x - 1 : ranges[0];
    const int beginY = whole[1] ? 0 : -ranges[1], endY = whole[1] ? numberOfCells.y - 1 : ranges[1];
    const int beginZ = whole[2] ? 0 : -ranges[2], endZ = whole[2] ? numberOfCells.z - 1 : ranges[2];

    for (int dx = beginX; dx <= endX; ++dx)
    {
      int ix, tx;
      double gapX;
      stencilCell(0, dx, whole[0], centreCell[0], cellWidth[0], centreCoordinate[0], ix, tx, gapX);
      if (gapX * gapX > radiusSquared) continue;  // the whole slab is beyond the radius
      for (int dy = beginY; dy <= endY; ++dy)
      {
        int iy, ty;
        double gapY;
        stencilCell(1, dy, whole[1], centreCell[1], cellWidth[1], centreCoordinate[1], iy, ty, gapY);
        // the whole row is beyond the radius (the per-axis gaps combine as a norm in an orthorhombic box)
        if ((orthorhombic ? gapX * gapX + gapY * gapY : gapY * gapY) > radiusSquared) continue;
        for (int dz = beginZ; dz <= endZ; ++dz)
        {
          int iz, tz;
          double gapZ;
          stencilCell(2, dz, whole[2], centreCell[2], cellWidth[2], centreCoordinate[2], iz, tz, gapZ);

          const std::size_t cell = cellIndex(ix, iy, iz);
          if (counts[cell] == 0) continue;

          // a lower bound on the distance of the centre to the cell: the per-axis gaps (a projection on the
          // face normal each) combine as a norm in an orthorhombic box, and as their maximum otherwise
          const double lowerBoundSquared = orthorhombic ? gapX * gapX + gapY * gapY + gapZ * gapZ
                                                        : std::max({gapX * gapX, gapY * gapY, gapZ * gapZ});
          if (lowerBoundSquared > radiusSquared) continue;

          // the translation of the cell relative to the wrapped positions (zero along a whole axis)
          const double3 shift = cellMatrix * double3(whole[0] ? 0.0 : box.x + static_cast<double>(tx),
                                                     whole[1] ? 0.0 : box.y + static_cast<double>(ty),
                                                     whole[2] ? 0.0 : box.z + static_cast<double>(tz));
          callback(cell, shift);
        }
      }
    }
  }

  /**
   * \brief Visits every unordered pair of records that can lie within the list's cut-off of each other
   * exactly once.
   *
   * The pairs within a cell (slot k < l) and the pairs between a cell and every distinct cell of its
   * stencil with a larger index are enumerated (a cell the stencil reaches twice around a short axis is
   * taken once). The callback receives '(const Record &a, const Record &b)' in no particular order; the
   * caller does the distance test (with the minimum image) and the same-molecule exclusion. This is the
   * full-energy counterpart of 'forEachRecordWithin'.
   */
  template <typename Callback>
  inline void forEachPairOnce(Callback&& callback) const
  {
    if (!valid || !enabled) return;
    int ranges[3];
    bool whole[3];
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const int n = numberOfCells[static_cast<int>(axis)];
      ranges[axis] = stencilRange(axis, cutOff);
      whole[axis] = (n <= 1) || (2 * ranges[axis] + 1 >= n);
    }
    const int beginX = whole[0] ? 0 : -ranges[0], endX = whole[0] ? numberOfCells.x - 1 : ranges[0];
    const int beginY = whole[1] ? 0 : -ranges[1], endY = whole[1] ? numberOfCells.y - 1 : ranges[1];
    const int beginZ = whole[2] ? 0 : -ranges[2], endZ = whole[2] ? numberOfCells.z - 1 : ranges[2];

    // the stencil of a cell can reach the same cell by two offsets only around a short axis; such a cell is
    // taken once through a per-cell stamp
    std::vector<std::uint32_t> stamp(totalNumberOfCells, 0u);

    for (int cx = 0; cx < numberOfCells.x; ++cx)
    {
      for (int cy = 0; cy < numberOfCells.y; ++cy)
      {
        for (int cz = 0; cz < numberOfCells.z; ++cz)
        {
          const std::size_t cell = cellIndex(cx, cy, cz);
          const std::size_t base = cell * capacity;
          const std::uint32_t count = counts[cell];
          if (count == 0) continue;

          // pairs within the cell
          for (std::uint32_t k = 0; k + 1 < count; ++k)
          {
            const Record a = record(base + k);
            for (std::uint32_t l = k + 1; l < count; ++l)
            {
              const Record b = record(base + l);
              callback(a, b);
            }
          }

          // pairs with the distinct stencil cells of larger index (each unordered cell pair once)
          const std::uint32_t mark = static_cast<std::uint32_t>(cell) + 1u;
          for (int dx = beginX; dx <= endX; ++dx)
          {
            const int ix = whole[0] ? dx : wrap(cx + dx, numberOfCells.x);
            for (int dy = beginY; dy <= endY; ++dy)
            {
              const int iy = whole[1] ? dy : wrap(cy + dy, numberOfCells.y);
              for (int dz = beginZ; dz <= endZ; ++dz)
              {
                const int iz = whole[2] ? dz : wrap(cz + dz, numberOfCells.z);
                const std::size_t other = cellIndex(ix, iy, iz);
                if (other <= cell || stamp[other] == mark) continue;
                stamp[other] = mark;
                const std::size_t otherBase = other * capacity;
                const std::uint32_t otherCount = counts[other];
                for (std::uint32_t k = 0; k != count; ++k)
                {
                  const Record a = record(base + k);
                  for (std::uint32_t l = 0; l != otherCount; ++l)
                  {
                    const Record b = record(otherBase + l);
                    callback(a, b);
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

  /// Mean number of records in the stencil of a cell (the atoms a query of the cut-off radius visits before
  /// the corner cells are skipped; the full atom count for a disabled list).
  [[nodiscard]] double averageNeighbourhoodSize() const;

 private:
  [[nodiscard]] static inline int wrap(int i, int n)
  {
    if (i < 0) return i + n;
    if (i >= n) return i - n;
    return i;
  }

  [[nodiscard]] static inline int cellCoordinate(double fractional, int n)
  {
    return std::min(static_cast<int>(fractional * static_cast<double>(n)), n - 1);
  }

  /// The cell index, box translation and distance gap (to the centre coordinate, along the face normal) of
  /// stencil offset 'd' along 'axis'. For an axis visited whole 'd' is the cell itself (no translation, no gap).
  inline void stencilCell(std::size_t axis, int d, bool wholeAxis, int centreCell, double cellWidth,
                          double centreCoordinate, int& index, int& translation, double& gap) const
  {
    const int n = numberOfCells[static_cast<int>(axis)];
    if (wholeAxis)
    {
      index = d;
      translation = 0;
      gap = 0.0;
      return;
    }
    const int unwrapped = centreCell + d;
    index = wrap(unwrapped, n);
    translation = (unwrapped - index) / n;
    // the slab [low, high) of the (unwrapped) cell along the face normal, relative to the box origin
    const double low = static_cast<double>(unwrapped) * cellWidth;
    const double high = low + cellWidth;
    gap = std::max({0.0, low - centreCoordinate, centreCoordinate - high});
  }

  void removeAtom(std::size_t atomIndex);
  [[nodiscard]] bool insertAtom(const Atom& atom, std::size_t atomIndex);
  /// Fills slot 'slot' from 'atom' (the position wrapped into the box).
  void writeSlot(std::size_t slot, const Atom& atom, std::size_t atomIndex);
  void copySlot(std::size_t to, std::size_t from);
  void resizeSlots(std::size_t numberOfSlots);
  /// Whether the slot holds a scaled (fractional) atom or one of a non-zero group.
  [[nodiscard]] bool slotPartiallyCoupled(std::size_t slot) const;
};
