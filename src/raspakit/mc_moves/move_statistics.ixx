module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <functional>
#include <map>
#include <print>
#include <source_location>
#include <string>
#endif

export module move_statistics;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

export template <typename T>
/**
 * \brief Collects and manages statistical data for move operations.
 *
 * The `MoveStatistics` struct accumulates statistics related to move operations in simulations.
 * It keeps track of counts, acceptance rates, and allows for optimization of move parameters
 * based on these statistics.
 *
 * \tparam T The numeric type used for counting and statistical calculations (e.g., `double`).
 */
struct MoveStatistics
{
  std::uint64_t versionNumber{2};  ///< Version number for serialization purposes.

  bool operator==(MoveStatistics<T> const &) const = default;

  T counts{};               ///< Number of move attempts.
  T constructed{};          ///< Number of moves constructed.
  T accepted{};             ///< Number of moves accepted.
  std::size_t allCounts{};  ///< Total number of counts across all types.
  T totalCounts{};          ///< Total move attempts across all simulations.
  T totalConstructed{};     ///< Total moves constructed across all simulations.
  T totalAccepted{};        ///< Total moves accepted across all simulations.
  T maxChange{};            ///< Maximum allowed change in move parameters.
  T targetAcceptance{0.5};  ///< Target acceptance rate for moves.
  T lowerLimit{};
  T upperLimit{};
  bool optimize{true};

  /**
   * \brief Resets the statistical counters.
   *
   * Clears all the counts and statistics to start fresh.
   */
  void clear()
  {
    counts = T{};
    constructed = T{};
    accepted = T{};
    allCounts = std::size_t{};
    totalCounts = T{};
    totalConstructed = T{};
    totalAccepted = T{};
  }

  /**
   * \brief Optimizes the acceptance rate of moves.
   *
   * Adjusts the maximum change allowed in move parameters based on the acceptance rate to approach the target
   * acceptance rate.
   */
  void optimizeAcceptance()
  {
    if (!optimize) return;
    T ratio = accepted / (counts + T(1.0));
    if constexpr (std::is_same_v<double, T>)
    {
      T scaling = std::clamp(ratio / targetAcceptance, T(0.5), T(1.5));
      maxChange = std::clamp(maxChange * scaling, lowerLimit, upperLimit);
    }
    else
    {
      T scaling = clamp(ratio / targetAcceptance, T(0.5), T(1.5));
      maxChange = clamp(maxChange * scaling, lowerLimit, upperLimit);
    }
    counts = T{};
    constructed = T{};
    accepted = T{};
  }

  template <class U>
  inline MoveStatistics &operator+=(const MoveStatistics<U> &b)
  {
    counts += b.counts;
    constructed += b.constructed;
    accepted += b.accepted;
    allCounts += b.allCounts;
    totalCounts += b.totalCounts;
    totalConstructed += b.totalConstructed;
    totalAccepted += b.totalAccepted;
    maxChange = 0.5 * (maxChange + b.maxChange);
    targetAcceptance = 0.5 * (targetAcceptance + b.targetAcceptance);
    lowerLimit = 0.5 * (lowerLimit + b.lowerLimit);
    upperLimit = 0.5 * (upperLimit + b.upperLimit);
    optimize = (optimize && b.optimize);
    return *this;
  }

  template <class U>
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MoveStatistics<U> &m);

  template <class U>
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MoveStatistics<U> &m);
};

export template <class T>
inline MoveStatistics<T> operator+(const MoveStatistics<T> &a, const MoveStatistics<T> &b)
{
  MoveStatistics<T> c;
  c.counts = a.counts + b.counts;
  c.constructed = a.constructed + b.constructed;
  c.accepted = a.accepted + b.accepted;
  c.allCounts = a.allCounts + b.allCounts;
  c.totalCounts = a.totalCounts + b.totalCounts;
  c.totalConstructed = a.totalConstructed + b.totalConstructed;
  c.totalAccepted = a.totalAccepted + b.totalAccepted;
  c.maxChange = 0.5 * (a.maxChange + b.maxChange);
  c.targetAcceptance = 0.5 * (a.targetAcceptance + b.targetAcceptance);
  c.lowerLimit = 0.5 * (a.lowerLimit + b.lowerLimit);
  c.upperLimit = 0.5 * (a.upperLimit + b.upperLimit);
  c.optimize = (a.optimize && b.optimize);
  return c;
}

export template <class T>
Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MoveStatistics<T> &m)
{
  archive << m.versionNumber;

  archive << m.counts;
  archive << m.constructed;
  archive << m.accepted;
  archive << m.allCounts;
  archive << m.totalCounts;
  archive << m.totalConstructed;
  archive << m.totalAccepted;
  archive << m.maxChange;
  archive << m.targetAcceptance;
  archive << m.lowerLimit;
  archive << m.upperLimit;
  archive << m.optimize;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

export template <class T>
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MoveStatistics<T> &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MoveStatistics' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> m.counts;
  archive >> m.constructed;
  archive >> m.accepted;
  archive >> m.allCounts;
  archive >> m.totalCounts;
  archive >> m.totalConstructed;
  archive >> m.totalAccepted;
  archive >> m.maxChange;
  archive >> m.targetAcceptance;
  archive >> m.lowerLimit;
  archive >> m.upperLimit;
  archive >> m.optimize;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("SimulationBox: Error in binary restart\n"));
  }
#endif

  return archive;
}
