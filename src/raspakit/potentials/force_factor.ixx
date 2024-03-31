export module force_factor;

import <cmath>;
import <fstream>;

import archive;

export struct ForceFactor
{
  double energy;
  double forceFactor;
  double dUdlambda;

  ForceFactor(double energy, double forceFactor, double dUdlambda):
      energy(energy),
      forceFactor(forceFactor),
      dUdlambda(dUdlambda) {}

  bool operator==(ForceFactor const&) const = default;

  inline ForceFactor& operator+=(const ForceFactor& b)
  {
    energy += b.energy;
    forceFactor += b.forceFactor;
    dUdlambda += b.dUdlambda;
    return *this;
  }

  inline ForceFactor& operator-=(const ForceFactor& b)
  {
    energy -= b.energy;
    forceFactor -= b.forceFactor;
    dUdlambda -= b.dUdlambda;
    return *this;
  }

  inline ForceFactor operator-() const
  {
    ForceFactor v(0.0, 0.0, 0.0);
    v.energy = -energy;
    v.forceFactor = -forceFactor;
    v.dUdlambda = -dUdlambda;
    return v;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceFactor &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceFactor &e);
};

export inline ForceFactor operator+(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.forceFactor = a.forceFactor +b.forceFactor;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;

  return m;
}

export inline ForceFactor operator-(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.forceFactor = a.forceFactor - b.forceFactor;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const ForceFactor& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.forceFactor = a.forceFactor * b.forceFactor;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const double& a, const ForceFactor& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a * b.energy;
  m.forceFactor = a * b.forceFactor;
  m.dUdlambda = a * b.dUdlambda;

  return m;
}

export inline ForceFactor operator*(const ForceFactor& a, const double& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy * b;
  m.forceFactor = a.forceFactor * b;
  m.dUdlambda = a.dUdlambda * b;

  return m;
}

export inline ForceFactor operator/(const ForceFactor& a, const double& b)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = a.energy / b;
  m.forceFactor = a.forceFactor / b;
  m.dUdlambda = a.dUdlambda / b;

  return m;
}

export inline ForceFactor sqrt(const ForceFactor & a)
{
  ForceFactor m(0.0, 0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.forceFactor = std::sqrt(a.forceFactor);
  m.dUdlambda = std::sqrt(a.dUdlambda);

  return m;
}
