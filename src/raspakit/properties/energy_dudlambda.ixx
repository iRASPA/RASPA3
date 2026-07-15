module;

export module energy_dudlambda;

import std;

import archive;

/**
 * \brief Stores an energy term together with its thermodynamic-integration derivative dU/dlambda.
 *
 * EnergyDuDlambda is the accumulator type of the energy bookkeeping (EnergyStatus and its
 * per-component sub-structures). The element-wise arithmetic and sqrt support the
 * block-averaging statistics; the Archive operators store it in binary restart files.
 */
export struct EnergyDuDlambda
{
  double energy;     ///< The energy value.
  double dUdlambda;  ///< The derivative of energy with respect to lambda.

  /**
   * \brief Constructs an EnergyDuDlambda with specified energy and derivative.
   *
   * \param energy The energy value.
   * \param dUdlambda The derivative of energy with respect to lambda.
   */
  EnergyDuDlambda(double energy, double dUdlambda) : energy(energy), dUdlambda(dUdlambda) {}

  bool operator==(EnergyDuDlambda const&) const = default;

  inline EnergyDuDlambda& operator+=(const EnergyDuDlambda& b)
  {
    energy += b.energy;
    dUdlambda += b.dUdlambda;
    return *this;
  }

  inline EnergyDuDlambda& operator-=(const EnergyDuDlambda& b)
  {
    energy -= b.energy;
    dUdlambda -= b.dUdlambda;
    return *this;
  }

  inline EnergyDuDlambda operator-() const
  {
    EnergyDuDlambda v(0.0, 0.0);
    v.energy = -energy;
    v.dUdlambda = -dUdlambda;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnergyDuDlambda& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnergyDuDlambda& e);
};

export Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnergyDuDlambda& e);
export Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnergyDuDlambda& e);

export inline EnergyDuDlambda operator+(const EnergyDuDlambda& a, const EnergyDuDlambda& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a.energy + b.energy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;

  return m;
}

export inline EnergyDuDlambda operator-(const EnergyDuDlambda& a, const EnergyDuDlambda& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a.energy - b.energy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;

  return m;
}

export inline EnergyDuDlambda operator*(const EnergyDuDlambda& a, const EnergyDuDlambda& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a.energy * b.energy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;

  return m;
}

export inline EnergyDuDlambda operator*(const double& a, const EnergyDuDlambda& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a * b.energy;
  m.dUdlambda = a * b.dUdlambda;

  return m;
}

export inline EnergyDuDlambda operator*(const EnergyDuDlambda& a, const double& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a.energy * b;
  m.dUdlambda = a.dUdlambda * b;

  return m;
}

export inline EnergyDuDlambda operator/(const EnergyDuDlambda& a, const double& b)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = a.energy / b;
  m.dUdlambda = a.dUdlambda / b;

  return m;
}

export inline EnergyDuDlambda sqrt(const EnergyDuDlambda& a)
{
  EnergyDuDlambda m(0.0, 0.0);
  m.energy = std::sqrt(a.energy);
  m.dUdlambda = std::sqrt(a.dUdlambda);

  return m;
}
