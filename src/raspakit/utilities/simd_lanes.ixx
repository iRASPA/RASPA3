module;

export module simd_lanes;

import std;

/**
 * \brief The portable SIMD vocabulary shared by the pair kernels: the cluster kernel of the spatial-decomposition
 * MD engine (ClusterPairKernel) and the neighbourhood kernel of the Monte Carlo trial energies
 * (Interactions::computeNeighbourhoodEnergy).
 */
export namespace SIMD
{
/// Width of the SIMD registers of the target in bytes (AVX-512 is treated as two 256-bit halves: the kernels
/// cap their lane counts so that a block fits a 32-bit interaction mask).
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
  /// Four consecutive table entries per lane: a = table[4 i], b = table[4 i + 1], c = table[4 i + 2],
  /// d = table[4 i + 3] (the coefficient records of a piecewise cubic).
  static void gather4(const Real* table, const std::uint32_t* index, Lanes& a, Lanes& b, Lanes& c, Lanes& d)
  {
    for (std::size_t n = 0; n < N; ++n)
    {
      const Real* record = table + 4 * static_cast<std::size_t>(index[n]);
      a.v[n] = record[0];
      b.v[n] = record[1];
      c.v[n] = record[2];
      d.v[n] = record[3];
    }
  }

  Real sum() const
  {
    Real s = Real(0);
    for (std::size_t n = 0; n < N; ++n) s += v[n];
    return s;
  }
  Real maximum() const
  {
    Real m = v[0];
    for (std::size_t n = 1; n < N; ++n) m = v[n] > m ? v[n] : m;
    return m;
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
  /// Rounds every lane to the nearest integer, halves away from zero (the convention of
  /// SimulationBox::applyPeriodicBoundaryConditions).
  friend Lanes round(Lanes a)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = std::round(a.v[n]);
    return a;
  }
  friend Lanes floor(Lanes a)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = std::floor(a.v[n]);
    return a;
  }
  friend Lanes max(Lanes a, Real b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = a.v[n] > b ? a.v[n] : b;
    return a;
  }
  friend Lanes max(Lanes a, const Lanes& b)
  {
    for (std::size_t n = 0; n < N; ++n) a.v[n] = a.v[n] > b.v[n] ? a.v[n] : b.v[n];
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
  /// 1 where a != b, 0 elsewhere.
  friend Lanes notEqual(const Lanes& a, Real b)
  {
    Lanes r;
    for (std::size_t n = 0; n < N; ++n) r.v[n] = a.v[n] != b ? Real(1) : Real(0);
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
}  // namespace SIMD
