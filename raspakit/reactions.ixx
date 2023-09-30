export module reactions;

import <vector>;
import <numbers>;
import <string>;
import <sstream>;

import reaction;

export struct Reactions
{
  std::vector<Reaction> list;


  std::string printStatus() const;
};
