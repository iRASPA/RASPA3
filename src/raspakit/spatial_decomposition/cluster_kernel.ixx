module;

export module spatial_decomposition_cluster_kernel;

import std;

import double3;
import double3x3;
import spatial_decomposition_cell_list;
import spatial_decomposition_pair_kernel;

namespace ClusterKernelDetail
{
/// Width of the SIMD registers of the target in bytes.
#if defined(__AVX512F__) || defined(__AVX2__)
constexpr std::size_t simdBytes = 32;
#else
constexpr std::size_t simdBytes = 16;
#endif

/**
 * \brief A fixed-width vector of lanes with element-wise arithmetic, written as plain loops over a fixed count.
 *
 * The loops are trivially vectorizable (fixed trip count, no dependences), so the compiler maps them onto the
 * SIMD unit of the target (NEON, AVX2, AVX-512) without intrinsics; the type is portable to every compiler the
 * code builds with. Comparisons produce 0 / 1 masks of the same type, and masking is a multiplication (the
 * operands are kept finite, so no NaN can leak through a masked lane).
 */
template <typename Real, std::size_t N>
struct alignas(sizeof(Real) * N) Lanes
{
  Real v[N];

  static Lanes broadcast(Real a)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = a;
    return r;
  }
  static Lanes zero() { return broadcast(Real(0)); }
  static Lanes load(const Real* p)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = p[n];
    return r;
  }
  /// 1 in the lanes whose bit is set, 0 elsewhere.
  static Lanes fromBits(std::uint32_t bits)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = ((bits >> n) & 1u) ? Real(1) : Real(0);
    return r;
  }
  /// Bit n set where lane n is non-zero (the inverse of fromBits for 0 / 1 masks).
  std::uint32_t toBits() const
  {
    std::uint32_t bits = 0;
    for (std::size_t n = 0; n < N; ++n) bits |= (v[n] != Real(0) ? 1u : 0u) << n;
    return bits;
  }
  static Lanes gather(const Real* table, const std::uint32_t* index)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = table[index[n]];
    return r;
  }

  Real sum() const
  {
    Real s = Real(0);
    for (std::size_t n = 0; n < N; ++n) s += v[n];
    return s;
  }

  friend Lanes operator+(Lanes a, const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] += b.v[n];
    return a;
  }
  friend Lanes operator+(Real a, Lanes b)
  {
    for (std::size_t n = 0; n < N; ++n) b.v[n] += a;
    return b;
  }
  friend Lanes operator-(Lanes a, const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] -= b.v[n];
    return a;
  }
  friend Lanes operator*(Lanes a, const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] *= b.v[n];
    return a;
  }
  friend Lanes operator*(Real a, Lanes b)
  {
    for (std::size_t n = 0; n < N; ++n) b.v[n] *= a;
    return b;
  }
  friend Lanes operator*(Lanes a, Real b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] *= b;
    return a;
  }
  friend Lanes operator-(Real a, Lanes b)
  {
    for (std::size_t n = 0; n < N; ++n) b.v[n] = a - b.v[n];
    return b;
  }
  friend Lanes operator-(Lanes a, Real b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] -= b;
    return a;
  }
  friend Lanes operator-(Lanes a)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = -a.v[n];
    return a;
  }
  Lanes& operator+=(const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) v[n] += b.v[n];
    return *this;
  }
  Lanes& operator-=(const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) v[n] -= b.v[n];
    return *this;
  }
  friend Lanes reciprocal(Lanes a)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = Real(1) / a.v[n];
    return a;
  }
  friend Lanes sqrt(Lanes a)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = std::sqrt(a.v[n]);
    return a;
  }
  friend Lanes max(Lanes a, Real b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = a.v[n] > b ? a.v[n] : b;
    return a;
  }
  friend Lanes min(Lanes a, Real b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = a.v[n] < b ? a.v[n] : b;
    return a;
  }
  /// 1 where a < b, 0 elsewhere.
  friend Lanes less(const Lanes& a, Real b)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = a.v[n] < b ? Real(1) : Real(0);
    return r;
  }
  /// 1 where a != 0, 0 elsewhere.
  friend Lanes nonZero(const Lanes& a)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = a.v[n] != Real(0) ? Real(1) : Real(0);
    return r;
  }
  /// Integer truncation of non-negative lanes.
  void truncate(std::uint32_t* index) const
  {
    for (std::size_t n = 0; n < N; ++n) index[n] = static_cast<std::uint32_t>(v[n]);
  }
};

/// exp(a) for single-precision lanes with a <= 0 (arguments below -87 give the smallest normal numbers): round
/// to a multiple of ln 2, a degree-6 polynomial on the remainder (Cephes expf) and the power of two through the
/// exponent bits. Relative error about 2e-7; no NaN or infinity for any finite input.
template <std::size_t N>
Lanes<float, N> expNonPositive(Lanes<float, N> a)
{
  constexpr float log2e = 1.44269504088896341f;
  constexpr float ln2High = 0.693359375f;
  constexpr float ln2Low = -2.12194440e-4f;
  Lanes<float, N> result;
  for (std::size_t n = 0; n < N; ++n)
  {
    const float x = std::max(a.v[n], -87.0f);
    const float k = std::round(x * log2e);
    const float r = (x - k * ln2High) - k * ln2Low;
    float p = 1.9875691500e-4f;
    p = p * r + 1.3981999507e-3f;
    p = p * r + 8.3334519073e-3f;
    p = p * r + 4.1665795894e-2f;
    p = p * r + 1.6666665459e-1f;
    p = p * r + 5.0000001201e-1f;
    const float y = (p * r * r + r) + 1.0f;
    const std::int32_t exponent = static_cast<std::int32_t>(k) + 127;
    const float scale = std::bit_cast<float>(static_cast<std::uint32_t>(exponent) << 23);
    result.v[n] = y * scale;
  }
  return result;
}

/// Spreads the low 10 bits of `a` to every third bit (Morton interleaving).
inline std::uint32_t spreadBits(std::uint32_t a)
{
  a &= 0x3FFu;
  a = (a | (a << 16)) & 0x030000FFu;
  a = (a | (a << 8)) & 0x0300F00Fu;
  a = (a | (a << 4)) & 0x030C30C3u;
  a = (a | (a << 2)) & 0x09249249u;
  return a;
}

/// Appends the local indices [begin, end) to `order`, sorted along a Morton (Z-order) curve over a grid whose
/// cell holds about one atom, so that consecutive atoms, and so the clusters made of them, are spatially compact.
inline void appendSpatialOrder(const LocalAtom* positions, std::size_t begin, std::size_t end,
                               std::vector<std::pair<std::uint32_t, std::uint32_t>>& scratch,
                               std::vector<std::uint32_t>& order)
{
  if (begin >= end) return;
  double3 low(positions[begin].x, positions[begin].y, positions[begin].z);
  double3 high = low;
  for (std::size_t l = begin; l < end; ++l)
  {
    low.x = std::min(low.x, positions[l].x);
    low.y = std::min(low.y, positions[l].y);
    low.z = std::min(low.z, positions[l].z);
    high.x = std::max(high.x, positions[l].x);
    high.y = std::max(high.y, positions[l].y);
    high.z = std::max(high.z, positions[l].z);
  }
  const double3 extent = high - low;
  const double volume = std::max(extent.x, 1.0) * std::max(extent.y, 1.0) * std::max(extent.z, 1.0);
  const double cellSize = std::max(std::cbrt(volume / static_cast<double>(end - begin)), 1.0);
  const double inverseCellSize = 1.0 / cellSize;
  constexpr double maximumIndex = 1023.0;

  scratch.clear();
  scratch.reserve(end - begin);
  for (std::size_t l = begin; l < end; ++l)
  {
    const std::uint32_t ix =
        static_cast<std::uint32_t>(std::clamp((positions[l].x - low.x) * inverseCellSize, 0.0, maximumIndex));
    const std::uint32_t iy =
        static_cast<std::uint32_t>(std::clamp((positions[l].y - low.y) * inverseCellSize, 0.0, maximumIndex));
    const std::uint32_t iz =
        static_cast<std::uint32_t>(std::clamp((positions[l].z - low.z) * inverseCellSize, 0.0, maximumIndex));
    const std::uint32_t key = spreadBits(ix) | (spreadBits(iy) << 1) | (spreadBits(iz) << 2);
    scratch.emplace_back(key, static_cast<std::uint32_t>(l));
  }
  std::sort(scratch.begin(), scratch.end());
  for (const auto& [key, l] : scratch) order.push_back(l);
}
}  // namespace ClusterKernelDetail

/**
 * \brief Cluster pair kernel of the spatial-decomposition MD engine: plain 12-6 Lennard-Jones plus the Ewald
 * real-space term between fully coupled atoms, evaluated on clusters of atoms with SIMD-friendly arithmetic in the
 * precision `Real`.
 *
 * The owned atoms of a sub-domain are grouped into i-clusters of `clusterI` atoms and all local atoms (owned,
 * then ghost images) into j-clusters of `clusterJ` atoms, in a spatial (Morton) order of the atoms so that the
 * clusters are compact; the per-atom neighbour lists of the CellList are folded into a list of (j-cluster,
 * interaction mask) blocks per i-cluster, the mask holding one bit per (i, j) pair of the block that is in the
 * neighbour list. Every pair of the system appears in exactly one block (the CellList lists are half lists), so
 * Newton's third law applies: the kernel evaluates a block as `clusterI` rows of `clusterJ` lanes, accumulating
 * the i-forces in lane registers over all blocks of an i-cluster and the j-forces in lane registers over the rows
 * of a block, and touches the force array once per block. Pairs outside the cutoffs and the padding of the
 * clusters are masked.
 *
 * Dual list: the block list built at a rebuild spans the cutoff plus the full Verlet skin. From it a pruned list
 * is derived with a smaller skin (`pruneSkin`), keeping only the blocks and pair bits within cutoff + pruneSkin at
 * the current positions; the pruned list is valid, and reused, while no local atom (owned or image, so the check
 * is local to the sub-domain) has moved more than pruneSkin / 2 since the prune, and is then rebuilt from the
 * outer list. The kernel thus evaluates far fewer pairs outside the cutoff per step than the lists contain, for
 * the price of a distance pass over the outer list every few steps.
 *
 * Positions are stored relative to the center of the sub-domain's owned atoms, so single precision keeps the
 * distances to about 1e-6 Angstrom in a box of any size. The per-atom forces (LocalForce), the energies and the
 * strain derivative are accumulated in double at the end of every block or i-cluster, whichever applies.
 *
 * The Ewald real-space term erfc(alpha r) / r is tabulated (cubic Hermite in r^2, the same table as the scalar
 * kernel) for `Real` = double, which makes `ClusterPairKernel<double>` reproduce the scalar kernel to rounding,
 * and evaluated in closed form for `Real` = float (Abramowitz-Stegun 7.1.26 with a vectorized exponential, absolute
 * error 1.5e-7 on erfc; no table gathers).
 */
export template <typename Real>
class ClusterPairKernel
{
 public:
  static constexpr std::size_t clusterI = 4;
  /// One SIMD register of lanes: 4 float / 2 double on 128-bit units (NEON, SSE), 8 float / 4 double on AVX2 and
  /// above (capped so that a block fits the 32-bit interaction mask).
  static constexpr std::size_t clusterJ = std::min<std::size_t>(ClusterKernelDetail::simdBytes / sizeof(Real), 8);
  static_assert(clusterI * clusterJ <= 32, "the interaction mask of a block is 32 bits");
  static constexpr bool analyticEwald = std::is_same_v<Real, float>;
  using Vec = ClusterKernelDetail::Lanes<Real, clusterJ>;

  struct Domain
  {
    std::size_t owned{0};
    std::size_t local{0};
    std::size_t iClusters{0};
    std::size_t jClusters{0};
    std::size_t padded{0};  ///< Length of the atom arrays (covers the padded i- and j-clusters).
    double3 origin{};       ///< Positions are stored relative to this point.

    // kernel order: owned atoms first (a spatial permutation of the owned local atoms), then the images
    std::vector<std::uint32_t> order{};  ///< Kernel index -> local index.
    std::vector<std::uint32_t> rank{};   ///< Local index -> kernel index.

    std::vector<Real> x{}, y{}, z{}, charge{};
    std::vector<std::uint32_t> typeOffset{};  ///< type * numberOfTypes per atom (row of the parameter tables).
    std::vector<std::uint32_t> typeIndex{};   ///< type per atom.
    std::vector<LocalForce> force{};          ///< Scratch: forces in kernel order.

    // outer list (cutoff + Verlet skin), built at the rebuild
    std::vector<std::uint32_t> pairStart{};    ///< CSR offsets per i-cluster.
    std::vector<std::uint32_t> pairCluster{};  ///< j-cluster of the block.
    std::vector<std::uint32_t> pairMask{};     ///< Bit i * clusterJ + j set when pair (i, j) is in the list.

    // pruned list (cutoff + pruneSkin at the positions of the last prune)
    bool pruned{false};
    std::vector<std::uint32_t> innerStart{}, innerCluster{}, innerMask{};
    std::vector<double> pruneX{}, pruneY{}, pruneZ{};  ///< Local positions at the last prune (local order).

    std::size_t prunes{0};

    std::vector<std::uint32_t> lastEntryOfCluster{};  ///< Build scratch: j-cluster -> last block index.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> sortScratch{};
  };

  /// Fixes the pair-type tables, the cutoffs, the Ewald parameters (when `useCharge`; the table must span the
  /// Coulomb cutoff, checked by the caller) and the pruning skin (0 or at least the Verlet skin: no pruning).
  void setParameters(std::span<const LennardJonesPair> lennardJones, std::size_t types, bool charge, double cutOffVDW,
                     double cutOffCharge, double conversionFactor, const EwaldRealSpaceTable* table, double alpha,
                     double verletSkin, double pruneSkinValue)
  {
    numberOfTypes = types;
    epsilon4.resize(types * types);
    sigma6.resize(types * types);
    shift.resize(types * types);
    for (std::size_t k = 0; k < types * types; ++k)
    {
      const LennardJonesPair& pair = lennardJones[k];
      const double sigma2 = 1.0 / pair.inverseSigma2;
      epsilon4[k] = static_cast<Real>(pair.epsilon4);
      sigma6[k] = static_cast<Real>(sigma2 * sigma2 * sigma2);
      shift[k] = static_cast<Real>(pair.shift);
    }
    cutOffVDWSquared = static_cast<Real>(cutOffVDW * cutOffVDW);
    cutOffChargeSquared = static_cast<Real>(cutOffCharge * cutOffCharge);
    coulombFactor = static_cast<Real>(conversionFactor);
    useCharge = charge;
    ewaldAlpha = static_cast<Real>(alpha);
    ewaldAlphaExact = alpha;
    if (useCharge && table)
    {
      const std::size_t n = table->value.size();
      tableValue.resize(n);
      tableSlope.resize(n);
      for (std::size_t k = 0; k < n; ++k)
      {
        tableValue[k] = static_cast<Real>(table->value[k]);
        tableSlope[k] = static_cast<Real>(table->derivative[k] * table->spacing);
      }
      tableRRMin = static_cast<Real>(table->rrMin);
      tableInverseSpacing = static_cast<Real>(table->inverseSpacing);
      tableMaximumPosition = static_cast<Real>(n - 2);  // k <= n - 2 so that k + 1 is a valid node
    }

    const double cutoff = charge ? std::max(cutOffVDW, cutOffCharge) : cutOffVDW;
    pruneSkin = (pruneSkinValue > 0.0 && pruneSkinValue < verletSkin) ? pruneSkinValue : 0.0;
    // a hair above cutoff + skin: a pair rounded (in Real) to just beyond the threshold at the prune cannot come
    // within the cutoff before the next prune by more than that rounding
    const double inner = cutoff + pruneSkin + 1e-4;
    innerCutoffSquared = static_cast<Real>(inner * inner);
    pruneDisplacementSquared = 0.25 * pruneSkin * pruneSkin;
    for (Domain& domain : domains) domain.pruned = false;
  }

  void resize(std::size_t numberOfDomains) { domains.resize(numberOfDomains); }

  std::size_t paddedLocalAtoms(std::size_t d) const { return domains[d].padded; }
  bool pruning() const { return pruneSkin > 0.0; }

  /// Builds the clusters and the block list of a sub-domain from its CellList lists (after CellList::buildLists).
  void buildLists(std::size_t d, const CellList::DomainLists& lists)
  {
    Domain& domain = domains[d];
    domain.owned = lists.ownedAtoms.size();
    domain.local = lists.numberOfLocalAtoms();
    domain.iClusters = (domain.owned + clusterI - 1) / clusterI;
    domain.jClusters = (domain.local + clusterJ - 1) / clusterJ;
    domain.padded = std::max(domain.iClusters * clusterI, domain.jClusters * clusterJ);
    domain.pruned = false;

    // spatial order of the owned atoms and of the images, separately (the owned atoms stay in front)
    domain.order.clear();
    domain.order.reserve(domain.local);
    ClusterKernelDetail::appendSpatialOrder(lists.positions.data(), 0, domain.owned, domain.sortScratch, domain.order);
    ClusterKernelDetail::appendSpatialOrder(lists.positions.data(), domain.owned, domain.local, domain.sortScratch,
                                            domain.order);
    domain.rank.resize(domain.local);
    for (std::size_t k = 0; k < domain.local; ++k) domain.rank[domain.order[k]] = static_cast<std::uint32_t>(k);

    domain.x.assign(domain.padded, Real(0));
    domain.y.assign(domain.padded, Real(0));
    domain.z.assign(domain.padded, Real(0));
    domain.charge.assign(domain.padded, Real(0));
    domain.typeOffset.assign(domain.padded, 0u);
    domain.typeIndex.assign(domain.padded, 0u);
    domain.force.resize(domain.padded);
    for (std::size_t k = 0; k < domain.local; ++k)
    {
      const std::uint32_t l = domain.order[k];
      domain.charge[k] = static_cast<Real>(lists.positions[l].charge);
      domain.typeIndex[k] = lists.localType[l];
      domain.typeOffset[k] = static_cast<std::uint32_t>(lists.localType[l] * numberOfTypes);
    }

    // the origin: the mean of the owned positions at the build
    double3 sum(0.0, 0.0, 0.0);
    for (std::size_t l = 0; l < domain.owned; ++l)
    {
      sum += double3(lists.positions[l].x, lists.positions[l].y, lists.positions[l].z);
    }
    domain.origin = domain.owned > 0 ? sum / static_cast<double>(domain.owned) : double3(0.0, 0.0, 0.0);

    // fold the per-atom half lists into blocks: one entry per (i-cluster, j-cluster) with a bit per listed pair
    domain.pairStart.assign(domain.iClusters + 1, 0u);
    domain.pairCluster.clear();
    domain.pairMask.clear();
    constexpr std::uint32_t noEntry = std::numeric_limits<std::uint32_t>::max();
    domain.lastEntryOfCluster.assign(domain.jClusters, noEntry);
    for (std::size_t ic = 0; ic < domain.iClusters; ++ic)
    {
      const std::uint32_t start = static_cast<std::uint32_t>(domain.pairCluster.size());
      domain.pairStart[ic] = start;
      for (std::size_t a = 0; a < clusterI; ++a)
      {
        const std::size_t i = ic * clusterI + a;
        if (i >= domain.owned) break;
        const std::uint32_t li = domain.order[i];
        for (std::uint32_t n = lists.neighbourStart[li]; n < lists.neighbourStart[li + 1]; ++n)
        {
          const std::uint32_t j = domain.rank[lists.neighbourList[n]];
          const std::uint32_t jc = j / static_cast<std::uint32_t>(clusterJ);
          const std::uint32_t lane = j % static_cast<std::uint32_t>(clusterJ);
          // the last block created for this j-cluster is reusable only if it belongs to the current i-cluster
          std::uint32_t entry = domain.lastEntryOfCluster[jc];
          if (entry == noEntry || entry < start)
          {
            entry = static_cast<std::uint32_t>(domain.pairCluster.size());
            domain.pairCluster.push_back(jc);
            domain.pairMask.push_back(0u);
            domain.lastEntryOfCluster[jc] = entry;
          }
          domain.pairMask[entry] |= 1u << (a * clusterJ + lane);
        }
      }
    }
    domain.pairStart[domain.iClusters] = static_cast<std::uint32_t>(domain.pairCluster.size());

    refreshPositions(d, lists);
  }

  /// Copies the current local positions (after CellList::gatherPositions) into the cluster arrays.
  void refreshPositions(std::size_t d, const CellList::DomainLists& lists)
  {
    Domain& domain = domains[d];
    const LocalAtom* positions = lists.positions.data();
    const std::uint32_t* order = domain.order.data();
    const double ox = domain.origin.x;
    const double oy = domain.origin.y;
    const double oz = domain.origin.z;
    Real* x = domain.x.data();
    Real* y = domain.y.data();
    Real* z = domain.z.data();
    for (std::size_t k = 0; k < domain.local; ++k)
    {
      const LocalAtom& p = positions[order[k]];
      x[k] = static_cast<Real>(p.x - ox);
      y[k] = static_cast<Real>(p.y - oy);
      z[k] = static_cast<Real>(p.z - oz);
    }
  }

  /// Re-prunes the list of a sub-domain when it has none or a local atom has moved more than pruneSkin / 2 since
  /// the last prune (call after refreshPositions). Returns whether a prune was done.
  bool pruneIfNeeded(std::size_t d, const CellList::DomainLists& lists)
  {
    if (!pruning()) return false;
    Domain& domain = domains[d];
    const LocalAtom* positions = lists.positions.data();
    if (domain.pruned)
    {
      const double* px = domain.pruneX.data();
      const double* py = domain.pruneY.data();
      const double* pz = domain.pruneZ.data();
      double maximum = 0.0;
      for (std::size_t l = 0; l < domain.local; ++l)
      {
        const double dx = positions[l].x - px[l];
        const double dy = positions[l].y - py[l];
        const double dz = positions[l].z - pz[l];
        maximum = std::max(maximum, dx * dx + dy * dy + dz * dz);
      }
      if (maximum <= pruneDisplacementSquared) return false;
    }

    domain.pruneX.resize(domain.local);
    domain.pruneY.resize(domain.local);
    domain.pruneZ.resize(domain.local);
    for (std::size_t l = 0; l < domain.local; ++l)
    {
      domain.pruneX[l] = positions[l].x;
      domain.pruneY[l] = positions[l].y;
      domain.pruneZ[l] = positions[l].z;
    }
    prune(domain);
    domain.pruned = true;
    ++domain.prunes;
    return true;
  }

  std::size_t numberOfBlocks() const
  {
    std::size_t total = 0;
    for (const Domain& domain : domains) total += domain.pairCluster.size();
    return total;
  }
  std::size_t numberOfPrunedBlocks() const
  {
    std::size_t total = 0;
    for (const Domain& domain : domains)
      total += domain.pruned ? domain.innerCluster.size() : domain.pairCluster.size();
    return total;
  }
  std::size_t numberOfPrunes() const
  {
    std::size_t total = 0;
    for (const Domain& domain : domains) total += domain.prunes;
    return total;
  }

  /// Evaluates all blocks of a sub-domain: adds the forces to `force` (indexed like the local atoms, at least
  /// paddedLocalAtoms long), the energies and, with `withVirial`, the strain derivative sum f (x) dr.
  void compute(std::size_t d, LocalForce* force, bool withVirial, double& energyVDW, double& energyCharge,
               double3x3& strain)
  {
    Domain& domain = domains[d];
    const Real* x = domain.x.data();
    const Real* y = domain.y.data();
    const Real* z = domain.z.data();
    const Real* q = domain.charge.data();
    const std::uint32_t* typeOffset = domain.typeOffset.data();
    const std::uint32_t* typeIndex = domain.typeIndex.data();
    const Real* eps4 = epsilon4.data();
    const Real* sig6 = sigma6.data();
    const Real* shf = shift.data();
    const Real* tableV = tableValue.data();
    const Real* tableS = tableSlope.data();
    const bool usePruned = domain.pruned;
    const std::uint32_t* blockStart = usePruned ? domain.innerStart.data() : domain.pairStart.data();
    const std::uint32_t* blockCluster = usePruned ? domain.innerCluster.data() : domain.pairCluster.data();
    const std::uint32_t* blockMask = usePruned ? domain.innerMask.data() : domain.pairMask.data();
    constexpr std::uint32_t laneBits = (clusterJ >= 32) ? 0xFFFFFFFFu : ((1u << clusterJ) - 1u);
    // pairs are never this close between molecules; the clamp keeps the masked (padding, beyond-cutoff) lanes
    // finite
    constexpr Real minimumRR = Real(1e-4);
    // Abramowitz-Stegun 7.1.26: erfc(x) = t (a1 + t (a2 + t (a3 + t (a4 + t a5)))) exp(-x^2), t = 1 / (1 + p x)
    constexpr Real asP = Real(0.3275911);
    constexpr Real asA1 = Real(0.254829592);
    constexpr Real asA2 = Real(-0.284496736);
    constexpr Real asA3 = Real(1.421413741);
    constexpr Real asA4 = Real(-1.453152027);
    constexpr Real asA5 = Real(1.061405429);
    const Real alpha = ewaldAlpha;
    const Real alphaSquared = alpha * alpha;
    const Real alphaOverSqrtPi = static_cast<Real>(static_cast<double>(alpha) / std::sqrt(std::numbers::pi));

    LocalForce* kernelForce = domain.force.data();
    std::fill(kernelForce, kernelForce + domain.padded, LocalForce{});

    double totalVDW = 0.0;
    double totalCharge = 0.0;

    for (std::size_t ic = 0; ic < domain.iClusters; ++ic)
    {
      const std::size_t i0 = ic * clusterI;
      Real xi[clusterI], yi[clusterI], zi[clusterI], qi[clusterI];
      std::uint32_t rowOffset[clusterI];
      Vec fix[clusterI], fiy[clusterI], fiz[clusterI];
      for (std::size_t a = 0; a < clusterI; ++a)
      {
        xi[a] = x[i0 + a];
        yi[a] = y[i0 + a];
        zi[a] = z[i0 + a];
        qi[a] = q[i0 + a];
        rowOffset[a] = typeOffset[i0 + a];
        fix[a] = Vec::zero();
        fiy[a] = Vec::zero();
        fiz[a] = Vec::zero();
      }
      Vec eVDW = Vec::zero();
      Vec eCharge = Vec::zero();

      for (std::uint32_t e = blockStart[ic]; e < blockStart[ic + 1]; ++e)
      {
        const std::size_t j0 = static_cast<std::size_t>(blockCluster[e]) * clusterJ;
        const std::uint32_t mask = blockMask[e];
        const Vec xj = Vec::load(x + j0);
        const Vec yj = Vec::load(y + j0);
        const Vec zj = Vec::load(z + j0);
        const Vec qj = Vec::load(q + j0);
        const std::uint32_t* tj = typeIndex + j0;
        Vec fjx = Vec::zero(), fjy = Vec::zero(), fjz = Vec::zero();

        for (std::size_t a = 0; a < clusterI; ++a)
        {
          const std::uint32_t bits = (mask >> (a * clusterJ)) & laneBits;
          if (bits == 0) continue;
          const Vec inList = Vec::fromBits(bits);

          const Vec dx = xi[a] - xj;
          const Vec dy = yi[a] - yj;
          const Vec dz = zi[a] - zj;
          const Vec rr = dx * dx + dy * dy + dz * dz;
          const Vec rrSafe = max(rr, minimumRR);
          const Vec invRR = reciprocal(rrSafe);

          // Lennard-Jones: parameters of the pair types gathered per lane
          std::uint32_t index[clusterJ];
          for (std::size_t n = 0; n < clusterJ; ++n) index[n] = rowOffset[a] + tj[n];
          const Vec e4 = Vec::gather(eps4, index);
          const Vec s6 = Vec::gather(sig6, index);
          const Vec sh = Vec::gather(shf, index);
          const Vec maskVDW = inList * less(rr, cutOffVDWSquared);
          const Vec invRR3 = invRR * invRR * invRR;
          const Vec rri3 = s6 * invRR3;
          const Vec rri6 = rri3 * rri3;
          eVDW += maskVDW * (e4 * (rri6 - rri3) - sh);
          Vec factor = maskVDW * (Real(12) * e4 * rri3 * (Real(0.5) - rri3) * invRR);

          if (useCharge)
          {
            const Vec qq = qi[a] * qj;
            const Vec maskCharge = inList * less(rr, cutOffChargeSquared) * nonZero(qq);
            Vec u, dudrr;
            if constexpr (analyticEwald)
            {
              // erfc(alpha r) / r and its derivative to r^2 in closed form
              const Vec r = sqrt(rrSafe);
              const Vec invR = r * invRR;
              const Vec ar = alpha * r;
              const Vec t = reciprocal(Real(1) + asP * ar);
              const Vec gauss = ClusterKernelDetail::expNonPositive(-(alphaSquared * rrSafe));
              const Vec erfcValue = t * (asA1 + t * (asA2 + t * (asA3 + t * (asA4 + t * asA5)))) * gauss;
              u = erfcValue * invR;
              dudrr = -((alphaOverSqrtPi * gauss + Real(0.5) * erfcValue * invR) * invRR);
            }
            else
            {
              // tabulated erfc(alpha r) / r in r^2: cubic Hermite interpolation with gathered nodes
              Vec position = (rrSafe - tableRRMin) * tableInverseSpacing;
              position = min(max(position, Real(0)), tableMaximumPosition);
              std::uint32_t node[clusterJ], next[clusterJ];
              position.truncate(node);
              for (std::size_t n = 0; n < clusterJ; ++n) next[n] = node[n] + 1;
              Vec t = position;
              for (std::size_t n = 0; n < clusterJ; ++n) t.v[n] -= static_cast<Real>(node[n]);
              const Vec p0 = Vec::gather(tableV, node);
              const Vec p1 = Vec::gather(tableV, next);
              const Vec m0 = Vec::gather(tableS, node);
              const Vec m1 = Vec::gather(tableS, next);
              const Vec t2 = t * t;
              const Vec t3 = t2 * t;
              u = (Real(2) * t3 - Real(3) * t2 + Vec::broadcast(Real(1))) * p0 + (t3 - Real(2) * t2 + t) * m0 +
                  (Real(3) * t2 - Real(2) * t3) * p1 + (t3 - t2) * m1;
              dudrr = tableInverseSpacing *
                      ((Real(6) * t2 - Real(6) * t) * p0 + (Real(3) * t2 - Real(4) * t + Vec::broadcast(Real(1))) * m0 +
                       (Real(6) * t - Real(6) * t2) * p1 + (Real(3) * t2 - Real(2) * t) * m1);
              // interacting lanes below the table range (closer than 1 Angstrom): the exact term
              const Vec below = maskCharge * less(rr, tableRRMin);
              if (below.sum() != Real(0)) [[unlikely]]
              {
                for (std::size_t n = 0; n < clusterJ; ++n)
                {
                  if (below.v[n] == Real(0)) continue;
                  double exactU, exactD;
                  exactEwald(static_cast<double>(rr.v[n]), exactU, exactD);
                  u.v[n] = static_cast<Real>(exactU);
                  dudrr.v[n] = static_cast<Real>(exactD);
                }
              }
            }
            const Vec prefactor = coulombFactor * (maskCharge * qq);
            eCharge += prefactor * u;
            factor += Real(2) * (prefactor * dudrr);
          }

          const Vec gx = factor * dx;
          const Vec gy = factor * dy;
          const Vec gz = factor * dz;
          fix[a] += gx;
          fiy[a] += gy;
          fiz[a] += gz;
          fjx -= gx;
          fjy -= gy;
          fjz -= gz;
        }

        LocalForce* fj = kernelForce + j0;
        for (std::size_t n = 0; n < clusterJ; ++n)
        {
          fj[n].x += static_cast<double>(fjx.v[n]);
          fj[n].y += static_cast<double>(fjy.v[n]);
          fj[n].z += static_cast<double>(fjz.v[n]);
        }
      }

      for (std::size_t a = 0; a < clusterI; ++a)
      {
        kernelForce[i0 + a].x += static_cast<double>(fix[a].sum());
        kernelForce[i0 + a].y += static_cast<double>(fiy[a].sum());
        kernelForce[i0 + a].z += static_cast<double>(fiz[a].sum());
      }
      totalVDW += static_cast<double>(eVDW.sum());
      totalCharge += static_cast<double>(eCharge.sum());
    }

    energyVDW += totalVDW;
    energyCharge += totalCharge;

    // back to the local order of the caller
    const std::uint32_t* order = domain.order.data();
    for (std::size_t k = 0; k < domain.local; ++k)
    {
      LocalForce& f = force[order[k]];
      f.x += kernelForce[k].x;
      f.y += kernelForce[k].y;
      f.z += kernelForce[k].z;
    }

    if (withVirial)
    {
      // Both ends of every listed pair are local atoms (owned or ghost image) whose forces are accumulated in
      // kernelForce, so sum_pairs g_ij (x) (x_i - x_j) = sum_k f_k (x) x_k with the local (image) positions; the
      // position offset drops out because the pair forces sum to zero.
      double sxx = 0.0, syx = 0.0, szx = 0.0, sxy = 0.0, syy = 0.0, szy = 0.0, sxz = 0.0, syz = 0.0, szz = 0.0;
      for (std::size_t k = 0; k < domain.local; ++k)
      {
        const double px = static_cast<double>(x[k]);
        const double py = static_cast<double>(y[k]);
        const double pz = static_cast<double>(z[k]);
        sxx += kernelForce[k].x * px;
        syx += kernelForce[k].y * px;
        szx += kernelForce[k].z * px;
        sxy += kernelForce[k].x * py;
        syy += kernelForce[k].y * py;
        szy += kernelForce[k].z * py;
        sxz += kernelForce[k].x * pz;
        syz += kernelForce[k].y * pz;
        szz += kernelForce[k].z * pz;
      }
      strain.ax += sxx;
      strain.bx += syx;
      strain.cx += szx;
      strain.ay += sxy;
      strain.by += syy;
      strain.cy += szy;
      strain.az += sxz;
      strain.bz += syz;
      strain.cz += szz;
    }
  }

 private:
  /// Derives the pruned list from the outer list at the current positions: the blocks with at least one pair
  /// within cutoff + pruneSkin, with the bits of the farther pairs cleared.
  void prune(Domain& domain) const
  {
    const Real* x = domain.x.data();
    const Real* y = domain.y.data();
    const Real* z = domain.z.data();
    constexpr std::uint32_t laneBits = (clusterJ >= 32) ? 0xFFFFFFFFu : ((1u << clusterJ) - 1u);
    const Real threshold = innerCutoffSquared;

    domain.innerStart.assign(domain.iClusters + 1, 0u);
    domain.innerCluster.clear();
    domain.innerMask.clear();
    domain.innerCluster.reserve(domain.pairCluster.size());
    domain.innerMask.reserve(domain.pairMask.size());
    for (std::size_t ic = 0; ic < domain.iClusters; ++ic)
    {
      const std::size_t i0 = ic * clusterI;
      domain.innerStart[ic] = static_cast<std::uint32_t>(domain.innerCluster.size());
      for (std::uint32_t e = domain.pairStart[ic]; e < domain.pairStart[ic + 1]; ++e)
      {
        const std::size_t j0 = static_cast<std::size_t>(domain.pairCluster[e]) * clusterJ;
        const std::uint32_t mask = domain.pairMask[e];
        const Vec xj = Vec::load(x + j0);
        const Vec yj = Vec::load(y + j0);
        const Vec zj = Vec::load(z + j0);
        std::uint32_t kept = 0;
        for (std::size_t a = 0; a < clusterI; ++a)
        {
          const std::uint32_t bits = (mask >> (a * clusterJ)) & laneBits;
          if (bits == 0) continue;
          const Vec dx = Vec::broadcast(x[i0 + a]) - xj;
          const Vec dy = Vec::broadcast(y[i0 + a]) - yj;
          const Vec dz = Vec::broadcast(z[i0 + a]) - zj;
          const Vec rr = dx * dx + dy * dy + dz * dz;
          kept |= (bits & less(rr, threshold).toBits()) << (a * clusterJ);
        }
        if (kept != 0)
        {
          domain.innerCluster.push_back(domain.pairCluster[e]);
          domain.innerMask.push_back(kept);
        }
      }
    }
    domain.innerStart[domain.iClusters] = static_cast<std::uint32_t>(domain.innerCluster.size());
  }

  /// erfc(alpha r) / r and its derivative to r^2 from the library functions (EwaldRealSpaceTable::exact).
  void exactEwald(double rr, double& u, double& dudrr) const
  {
    const double r = std::sqrt(rr);
    const double inverseR = 1.0 / r;
    const double erfcTerm = std::erfc(ewaldAlphaExact * r);
    const double gaussian = std::exp(-ewaldAlphaExact * ewaldAlphaExact * rr) * std::numbers::inv_sqrtpi_v<double>;
    u = erfcTerm * inverseR;
    dudrr = -0.5 * (erfcTerm * inverseR * inverseR + 2.0 * ewaldAlphaExact * gaussian * inverseR) * inverseR;
  }

  std::size_t numberOfTypes{0};
  std::vector<Real> epsilon4{}, sigma6{}, shift{};  ///< Per pair of types (row = type_i * numberOfTypes).
  Real cutOffVDWSquared{0};
  Real cutOffChargeSquared{0};
  Real coulombFactor{1};
  bool useCharge{false};
  Real ewaldAlpha{0};
  double ewaldAlphaExact{0.0};
  std::vector<Real> tableValue{}, tableSlope{};  ///< u and du/d(r^2) * spacing at the nodes.
  Real tableRRMin{1};
  Real tableInverseSpacing{1};
  Real tableMaximumPosition{0};

  double pruneSkin{0.0};  ///< 0: no pruning, the outer list is evaluated.
  Real innerCutoffSquared{0};
  double pruneDisplacementSquared{0.0};

  std::vector<Domain> domains{};
};
