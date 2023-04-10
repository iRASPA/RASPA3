export module reaction;

import <numbers>;
import <string>;
import <sstream>;

import print;

export struct Reaction
{
  Reaction(size_t numberOfComponents): reactant_components(numberOfComponents), product_components(numberOfComponents)
  std::vector<size_t> reactant_components;
  std::vector<size_t> product_components;
}
