module;

export module spatial_decomposition_domain_decomposition;

import std;

import int3;
import double3;

/**
 * \brief Partition of the box into rectangular sub-domains (in fractional coordinates), one per worker thread.
 *
 * The sub-domains form a `grid` of staggered slabs: the box is cut along x into `grid.x` slabs, every slab along y
 * into `grid.y` columns and every column along z into `grid.z` blocks. The cut positions are chosen from the atom
 * positions at every neighbour-list rebuild so that all sub-domains own (nearly) the same number of atoms, which
 * balances the pair work over the threads independently of the cell grid of the neighbour search. Ownership is
 * per atom (an atom belongs to the sub-domain whose block contains its wrapped fractional position); the
 * neighbour lists of a sub-domain reference the atoms of other sub-domains as ghost images, so the thread owning
 * a sub-domain never writes another thread's atoms.
 */
export struct DomainDecomposition
{
  int3 grid{1, 1, 1};  ///< Sub-domains per axis; the product is the number of threads.

  std::vector<double> planesX{};                            ///< grid.x + 1 fractional cut positions along x.
  std::vector<std::vector<double>> planesY{};               ///< Per x slab: grid.y + 1 cut positions along y.
  std::vector<std::vector<std::vector<double>>> planesZ{};  ///< Per (x slab, y column): grid.z + 1 cuts along z.

  /**
   * \brief Chooses the sub-domain grid for a thread count.
   *
   * Among all factorizations nx*ny*nz = threads with n_i <= cells_i (a sub-domain thinner than a cell of the
   * neighbour search is pointless), the one with the smallest total boundary area (sum over axes of the domain
   * cross-section, measured with the box widths) is taken; that minimizes the number of ghost atoms per thread. A
   * requested grid is validated and used as is. Throws when no admissible grid exists.
   */
  static int3 chooseGrid(std::size_t threads, int3 cells, double3 widths, std::optional<int3> requested);

  /// Places the cut planes for `grid` so that the given wrapped fractional positions are spread evenly.
  void build(int3 grid, std::span<const double3> fractionalPositions);

  /// Sub-domain that owns a wrapped fractional position (components in [0, 1)).
  std::uint32_t ownerOf(const double3& fractional) const;

  std::size_t numberOfDomains() const
  {
    return static_cast<std::size_t>(grid.x) * static_cast<std::size_t>(grid.y) * static_cast<std::size_t>(grid.z);
  }

  std::string status() const;
};
