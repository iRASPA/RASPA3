module;

module pair_interactions;

import std;

std::optional<std::size_t> PairInteractions::findType(const std::string& name) const
{
  std::vector<std::string>::const_iterator match = std::find(this->names.begin(), this->names.end(), name);
  if (match == this->names.end()) return std::nullopt;

  return static_cast<std::size_t>(std::distance(this->names.begin(), match));
}
