module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <type_traits>
#include <vector>
#endif

export module chiral_center;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

export struct ChiralCenter
{
  enum class Chirality : std::size_t
  {
    S_Chiral = 0,
    R_Chiral = 1
  };

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  Chirality type;
  std::array<std::size_t, 4> ids;

  ChiralCenter() : type(Chirality::S_Chiral), ids({0, 0, 0, 0}) {}

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ChiralCenter &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ChiralCenter &b);
};
