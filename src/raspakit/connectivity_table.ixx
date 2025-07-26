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
import <string>;
import <map>;
import <vector>;
import <array>;
import <fstream>;
import <type_traits>;
import <print>;
#endif

import stringutils;
import archive;

export struct ConnectivityTable
{
  uint64_t versionNumber{1};
  std::size_t numberOfBeads;
  std::vector<bool> table;

  ConnectivityTable() {};
  ConnectivityTable(size_t numberOfBeads): numberOfBeads(numberOfBeads), table(numberOfBeads * (numberOfBeads + 1) / 2) {};

  bool operator[](size_t i, size_t j) const
  {
    if (i < j) std::swap(i, j);
    return table[i * (i + 1) / 2 + j];
  }
  std::vector<bool>::reference operator[](size_t i, size_t j) 
  {
    if (i < j) std::swap(i, j);
    return table[i * (i + 1) / 2 + j];
  }

  std::string print(const std::string &prestring) const;

  std::tuple<std::optional<size_t>, size_t, std::vector<size_t>> nextBeads(const std::vector<size_t> &placedBeads);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ConnectivityTable &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ConnectivityTable &b);
};
