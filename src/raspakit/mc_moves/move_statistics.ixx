module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <algorithm>
#include <fstream>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <map>
#include <source_location>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#include <functional>
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
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
import <functional>;
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import archive;

export template <typename T>
struct MoveStatistics
{
  uint64_t versionNumber{1};

  bool operator==(MoveStatistics<T> const&) const = default;

  T counts{};
  T constructed{};
  T accepted{};
  std::size_t allCounts{};
  T totalCounts{};
  T totalConstructed{};
  T totalAccepted{};
  T maxChange{};
  T targetAcceptance{0.5};

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
