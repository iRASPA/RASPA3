module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <type_traits>
#include <vector>
#include <tuple>
#include <utility>
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
  ConnectivityTable(std::size_t numberOfBeads): numberOfBeads(numberOfBeads), table(numberOfBeads * (numberOfBeads + 1) / 2) {};

  bool operator[](std::size_t i, std::size_t j) const
  {
    if (i < j) std::swap(i, j);
    return table[i * (i + 1) / 2 + j];
  }
  std::vector<bool>::reference operator[](std::size_t i, std::size_t j) 
  {
    if (i < j) std::swap(i, j);
    return table[i * (i + 1) / 2 + j];
  }

  std::string print(const std::string &prestring) const;

  std::tuple<std::optional<std::size_t>, std::size_t, std::vector<std::size_t>> nextBeads(const std::vector<std::size_t> &placedBeads);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ConnectivityTable &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ConnectivityTable &b);
};
