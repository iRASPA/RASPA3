module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cmath>
#include <array>
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
import randomnumbers;
import double3;
import units;

export struct ChiralCenter
{
  enum class Chirality : size_t
  {
    S_Chiral = 0,
    R_Chiral = 1
  };

  uint64_t versionNumber{1};  ///< Version number for serialization.

  Chirality type;
  std::array<size_t, 4> ids;

  ChiralCenter() : type(Chirality::S_Chiral), ids({0, 0, 0, 0}) {}


  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ChiralCenter &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ChiralCenter &b);
};
