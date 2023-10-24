export module reactions;

import <vector>;
import <numbers>;
import <string>;
import <sstream>;
import <fstream>;

import archive;
import reaction;

export struct Reactions
{
  uint64_t versionNumber{ 1 };

  bool operator==(Reactions const&) const = default;

  std::vector<Reaction> list;

  std::string printStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reactions &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reactions &r);
};
