module;

export module chiral_center;

import std;

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
