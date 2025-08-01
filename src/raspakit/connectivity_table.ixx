module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <optional>
#endif

export module connectivity_table;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;

export struct ConnectivityTable
{
  std::uint64_t versionNumber{1};
  std::size_t numberOfBeads;
  std::vector<bool> table;

  ConnectivityTable() {};
  ConnectivityTable(std::size_t numberOfBeads)
      : numberOfBeads(numberOfBeads), table(numberOfBeads * numberOfBeads) {};

  bool operator[](std::size_t i, std::size_t j) const
  {
    return table[i * numberOfBeads + j];
  }
  std::vector<bool>::reference operator[](std::size_t i, std::size_t j)
  {
    return table[i * numberOfBeads + j];
  }

  std::string print(const std::string &prestring) const;

  std::tuple<std::optional<std::size_t>, std::size_t, std::vector<std::size_t>> nextBeads(
      const std::vector<std::size_t> &placedBeads) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ConnectivityTable &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ConnectivityTable &b);
};
