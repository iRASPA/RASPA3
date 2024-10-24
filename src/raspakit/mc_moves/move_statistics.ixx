module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
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
import <string>;
import <algorithm>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <map>;
import <print>;
import <functional>;
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
  uint64_t versionNumber{1};  ///< Version number for serialization purposes.

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
   *
   * \param lowerLimit The lower limit for `maxChange` adjustment.
   * \param upperLimit The upper limit for `maxChange` adjustment.
   */
  void optimizeAcceptance(T lowerLimit = T(0.0), T upperLimit = T(1.0))
  {
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
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MoveStatistics<U> &m);

  template <class U>
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MoveStatistics<U> &m);
};

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

  return archive;
}

export template <class T>
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MoveStatistics<T> &m)
{
  uint64_t versionNumber;
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

  return archive;
}
