export module energy_factor;

import <cmath>;
import <fstream>;

import archive;

export struct EnergyFactor
{
  double energy;
  double dUdlambda;

  EnergyFactor(double energy, double dUdlambda):
      energy(energy),
      dUdlambda(dUdlambda) {}

  bool operator==(EnergyFactor const&) const = default;

  inline EnergyFactor& operator+=(const EnergyFactor& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    return *this;
  }

  inline EnergyFactor& operator-=(const EnergyFactor& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    return *this;
  }

  inline EnergyFactor operator-() const
  {
    EnergyFactor v(0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    return v;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyFactor &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyFactor &e);
};

export inline EnergyFactor operator+(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;

  return m;
}

export inline EnergyFactor operator-(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;

  return m;
}

export inline EnergyFactor operator*(const EnergyFactor& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;

  return m;
}

export inline EnergyFactor operator*(const double& a, const EnergyFactor& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;

  return m;
}

export inline EnergyFactor operator*(const EnergyFactor& a, const double& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;

  return m;
}

export inline EnergyFactor operator/(const EnergyFactor& a, const double& b)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;

  return m;
}

export inline EnergyFactor sqrt(const EnergyFactor & a)
{
  EnergyFactor m(0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);

  return m;
}
