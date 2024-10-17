module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>
#endif

export module reaction;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <numbers>;
import <string>;
import <sstream>;
import <fstream>;
#endif

import archive;
import property_lambda_probability_histogram;
import json;

export struct Reaction
{
  Reaction() {};
  Reaction(size_t id, std::vector<size_t> reactantStoichiometry, std::vector<size_t> productStoichiometry)
      : id(id), reactantStoichiometry(reactantStoichiometry), productStoichiometry(productStoichiometry), lambda(5, 21)
  {
  }

  bool operator==(Reaction const &) const = default;

  uint64_t versionNumber{1};

  size_t id;
  std::vector<size_t> reactantStoichiometry;
  std::vector<size_t> productStoichiometry;

  PropertyLambdaProbabilityHistogram lambda;

  std::string printStatus() const;
  nlohmann::json jsonStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reaction &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reaction &r);
};
