module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>
#endif

export module reactions;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <numbers>;
import <string>;
import <sstream>;
import <fstream>;
#endif

import archive;
import reaction;
import json;

export struct Reactions
{
  uint64_t versionNumber{1};

  bool operator==(Reactions const &) const = default;

  std::vector<Reaction> list;

  std::string printStatus() const;
  nlohmann::json jsonStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reactions &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reactions &r);
};
