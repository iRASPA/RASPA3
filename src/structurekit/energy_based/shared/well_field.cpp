module;

module energy_shared_well_field;

import std;

WellField::WellField() {}

WellField::~WellField() {}

double WellField::deepestEnergy() const
{
  if (this->energy.empty()) return 0.0;
  return static_cast<double>(*std::ranges::min_element(this->energy));
}

double ScreenedCoulomb::exactly(double alpha, double rr)
{
  double r = std::sqrt(std::max(rr, 1.0e-8));
  return std::erfc(alpha * r) / r;
}

void ScreenedCoulomb::build(double alphaValue, double largestSquared)
{
  this->alpha = alphaValue;
  double range = std::max(largestSquared, smallestSquared + 1.0) - smallestSquared;
  this->scale = static_cast<double>(bins) / range;
  this->table.resize(bins + 2);
  for (std::size_t i = 0; i < this->table.size(); ++i)
  {
    this->table[i] = exactly(alphaValue, smallestSquared + static_cast<double>(i) / this->scale);
  }
}

double ScreenedCoulomb::at(double rr) const
{
  if (rr < smallestSquared) return exactly(this->alpha, rr);
  double x = (rr - smallestSquared) * this->scale;
  std::size_t i = static_cast<std::size_t>(x);
  if (i + 1 >= this->table.size()) return this->table.back();
  double f = x - static_cast<double>(i);
  return this->table[i] + f * (this->table[i + 1] - this->table[i]);
}
