export module reaction;

import <vector>;
import <numbers>;
import <string>;
import <sstream>;

import print;

export struct Reaction
{
  Reaction(size_t id, std::vector<size_t> reactantStoichiometry, std::vector<size_t> productStoichiometry): 
            id(id), reactantStoichiometry(reactantStoichiometry), productStoichiometry(productStoichiometry) {}

  size_t id;
  std::vector<size_t> reactantStoichiometry;
  std::vector<size_t> productStoichiometry;

  std::string printStatus() const;
};
