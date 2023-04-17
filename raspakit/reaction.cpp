module;

module reaction;

import <string>;
import <sstream>;

import print;

std::string Reaction::printStatus() const
{
  std::ostringstream stream;

  std::string r,p;
  for (const size_t & i : reactantStoichiometry)
    r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const size_t& i : productStoichiometry)
    p += (p.empty() ? "" : ",") + std::to_string(i);
  std::print(stream, "    reaction [{}]: Stoichiometry reactants {} -> products {}\n", id, r, p);

  return stream.str();
}