module;

module pair_interactions;

import std;

std::optional<std::size_t> PairInteractions::findType(const std::string& name) const
{
  std::vector<std::string>::const_iterator match = std::find(this->names.begin(), this->names.end(), name);
  if (match == this->names.end()) return std::nullopt;

  return static_cast<std::size_t>(std::distance(this->names.begin(), match));
}

namespace
{
double lennardJonesShift(double strength, double size, double cutOff)
{
  if (!(cutOff > 0.0) || !(size > 0.0) || !(strength > 0.0)) return 0.0;
  const double ratio = size / cutOff;
  const double ratio3 = ratio * ratio * ratio;
  const double ratio6 = ratio3 * ratio3;
  return 4.0 * strength * (ratio6 * ratio6 - ratio6);
}
}  // namespace

std::size_t PairInteractions::addSphericalProbe(std::string name, double strength, double size)
{
  if (name.empty()) name = "custom";

  const bool shifted =
      std::ranges::any_of(this->parameters, [](const PairParameters &pair) { return pair.shift != 0.0; });

  auto mix = [&](double strengthOther, double sizeOther) -> PairParameters
  {
    PairParameters pair;
    pair.sizeParameter = 0.5 * (size + sizeOther);
    pair.strengthParameter = std::sqrt(std::max(0.0, strength * strengthOther));
    if (shifted) pair.shift = lennardJonesShift(pair.strengthParameter, pair.sizeParameter, this->cutOffVDW);
    return pair;
  };

  if (std::optional<std::size_t> existing = this->findType(name); existing.has_value())
  {
    const std::size_t probe = existing.value();
    this->charges[probe] = 0.0;
    for (std::size_t other = 0; other < this->numberOfTypes; ++other)
    {
      PairParameters pair =
          (other == probe) ? PairParameters{size, strength, shifted ? lennardJonesShift(strength, size, this->cutOffVDW) : 0.0}
                           : mix(this->parameters[other * this->numberOfTypes + other].strengthParameter,
                                 this->parameters[other * this->numberOfTypes + other].sizeParameter);
      this->parameters[probe * this->numberOfTypes + other] = pair;
      this->parameters[other * this->numberOfTypes + probe] = pair;
    }
    return probe;
  }

  const std::size_t oldTypes = this->numberOfTypes;
  const std::size_t newTypes = oldTypes + 1;
  std::vector<PairParameters> expanded(newTypes * newTypes);

  for (std::size_t row = 0; row < oldTypes; ++row)
  {
    for (std::size_t column = 0; column < oldTypes; ++column)
    {
      expanded[row * newTypes + column] = this->parameters[row * oldTypes + column];
    }
  }

  this->names.push_back(name);
  this->charges.push_back(0.0);
  this->numberOfTypes = newTypes;
  this->parameters = std::move(expanded);

  const std::size_t probe = oldTypes;
  for (std::size_t other = 0; other < newTypes; ++other)
  {
    PairParameters pair =
        (other == probe)
            ? PairParameters{size, strength, shifted ? lennardJonesShift(strength, size, this->cutOffVDW) : 0.0}
            : mix(this->parameters[other * newTypes + other].strengthParameter,
                  this->parameters[other * newTypes + other].sizeParameter);
    this->parameters[probe * newTypes + other] = pair;
    this->parameters[other * newTypes + probe] = pair;
  }
  return probe;
}
