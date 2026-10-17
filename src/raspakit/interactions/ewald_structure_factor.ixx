module;

export module ewald_structure_factor;

import std;

import archive;
import running_energy;

/**
 * \brief The Ewald structure factor S(k) = sum_i q_i exp(i k.r_i) over the wave vectors in use, with the
 * per-group derivative sums of the dU/dlambda machinery, stored as a structure of arrays.
 *
 * The Monte Carlo energy differences stream through every wave vector per move: they read the stored
 * sum, write the trial sum, and -- only when dU/dlambda groups are present (fractional molecules) -- carry the
 * group sums along. Keeping the sums and the group sums in separate arrays, and leaving the group sums
 * untouched while no group-tagged atom exists ('hasGroupSums'), reduces the traffic per wave vector from 160
 * bytes (a pair of the sum and four group sums, read and written) to 32.
 *
 * Invariants: 'sum' and 'groupSums' have the same size. When 'hasGroupSums' is false the contents of
 * 'groupSums' are stale and read as zero through 'groups()'; a write through 'set()' with a non-zero group
 * sum raises the flag. Entries past the number of wave vectors in use are unused.
 */
export struct EwaldStructureFactor
{
  using GroupSums = std::array<std::complex<double>, maximumNumberOfDUDlambdaGroups>;

  std::vector<std::complex<double>> sum{};
  std::vector<GroupSums> groupSums{};
  bool hasGroupSums{false};

  EwaldStructureFactor() = default;
  EwaldStructureFactor(const EwaldStructureFactor&) = default;
  EwaldStructureFactor(EwaldStructureFactor&&) noexcept = default;
  EwaldStructureFactor& operator=(EwaldStructureFactor&&) noexcept = default;

  /// Copies the sums; the group sums only when they carry information.
  EwaldStructureFactor& operator=(const EwaldStructureFactor& other)
  {
    if (this == &other) return *this;
    sum = other.sum;
    hasGroupSums = other.hasGroupSums;
    if (hasGroupSums)
    {
      groupSums = other.groupSums;
    }
    else
    {
      groupSums.resize(sum.size());
    }
    return *this;
  }

  [[nodiscard]] std::size_t size() const { return sum.size(); }
  [[nodiscard]] bool empty() const { return sum.empty(); }

  friend bool operator==(const EwaldStructureFactor& a, const EwaldStructureFactor& b)
  {
    if (a.sum != b.sum) return false;
    if (a.hasGroupSums == b.hasGroupSums)
    {
      return !a.hasGroupSums || a.groupSums == b.groupSums;
    }
    for (std::size_t k = 0; k < a.sum.size(); ++k)
    {
      if (a.groups(k) != b.groups(k)) return false;
    }
    return true;
  }

  void clear()
  {
    sum.clear();
    groupSums.clear();
    hasGroupSums = false;
  }

  /// Room for at least 'numberOfWaveVectors' entries (the contents are kept).
  void reserveWaveVectors(std::size_t numberOfWaveVectors)
  {
    if (sum.size() < numberOfWaveVectors)
    {
      sum.resize(numberOfWaveVectors);
      groupSums.resize(numberOfWaveVectors);
    }
  }

  /// Room for 'numberOfWaveVectors' entries and no group sums yet: before a computation from scratch that
  /// writes every entry through 'set'.
  void beginFill(std::size_t numberOfWaveVectors)
  {
    reserveWaveVectors(numberOfWaveVectors);
    hasGroupSums = false;
  }

  /// Before a difference computation that reads 'stored' and writes every entry of this through 'set' (or the
  /// plain arrays). The chained moves pass one object as both: every entry is then read before it is written,
  /// so the flag is kept and the stale group sums are zeroed, which keeps 'groups()' exact once 'set' raises
  /// the flag part-way through. With two distinct objects this is 'beginFill'.
  void beginFill(std::size_t numberOfWaveVectors, const EwaldStructureFactor& stored)
  {
    if (this != &stored)
    {
      beginFill(numberOfWaveVectors);
      return;
    }
    reserveWaveVectors(numberOfWaveVectors);
    if (!hasGroupSums) std::fill_n(groupSums.begin(), numberOfWaveVectors, GroupSums{});
  }

  /// The group sums of wave vector 'k' (zeros while no group-tagged atom exists).
  [[nodiscard]] const GroupSums& groups(std::size_t k) const
  {
    static constexpr GroupSums zero{};
    return hasGroupSums ? groupSums[k] : zero;
  }

  /// Writes wave vector 'k'.
  void set(std::size_t k, const std::complex<double>& value, const GroupSums& groupValues)
  {
    sum[k] = value;
    if (!hasGroupSums)
    {
      bool nonZero = false;
      for (const std::complex<double>& g : groupValues) nonZero |= (g.real() != 0.0 || g.imag() != 0.0);
      if (!nonZero) return;
      // the first group sum: the stale entries before 'k' are zero from here on
      std::fill_n(groupSums.begin(), k, GroupSums{});
      hasGroupSums = true;
    }
    groupSums[k] = groupValues;
  }

  /// The sums as one bulk array; the group sums only when they carry information.
  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EwaldStructureFactor& s)
  {
    archive << s.sum;
    archive << s.hasGroupSums;
    if (s.hasGroupSums) archive << s.groupSums;
    return archive;
  }
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EwaldStructureFactor& s)
  {
    archive >> s.sum;
    archive >> s.hasGroupSums;
    if (s.hasGroupSums)
    {
      archive >> s.groupSums;
      if (s.groupSums.size() != s.sum.size()) throw std::runtime_error("malformed Ewald structure factor");
    }
    else
    {
      s.groupSums.assign(s.sum.size(), GroupSums{});
    }
    return archive;
  }
};
