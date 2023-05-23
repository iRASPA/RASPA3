export module reaction;

import <vector>;
import <numbers>;
import <string>;
import <sstream>;

import print;
import property_lambda_probability_histogram;

export struct Reaction
{
  Reaction(size_t id, std::vector<size_t> reactantStoichiometry, std::vector<size_t> productStoichiometry): 
            id(id), reactantStoichiometry(reactantStoichiometry), productStoichiometry(productStoichiometry),
            lambda(5,21)
 {}

  size_t id;
  std::vector<size_t> reactantStoichiometry;
  std::vector<size_t> productStoichiometry;

  PropertyLambdaProbabilityHistogram lambda;

  std::string printStatus() const;
};
