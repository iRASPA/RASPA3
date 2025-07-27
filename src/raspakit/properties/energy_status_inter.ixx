module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <print>
#include <sstream>
#include <string>
#endif

export module energy_status_inter;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import energy_factor;

export struct EnergyInter
{
  std::uint64_t versionNumber{1};

  Potentials::EnergyFactor VanDerWaals;
  Potentials::EnergyFactor VanDerWaalsTailCorrection;
  Potentials::EnergyFactor CoulombicReal;
  Potentials::EnergyFactor CoulombicFourier;
  Potentials::EnergyFactor totalInter;

  EnergyInter()
      : VanDerWaals(0.0, 0.0),
        VanDerWaalsTailCorrection(0.0, 0.0),
        CoulombicReal(0.0, 0.0),
        CoulombicFourier(0.0, 0.0),
        totalInter(0.0, 0.0)
  {
  }

  bool operator==(EnergyInter const&) const = default;

  void zero()
  {
    VanDerWaals = Potentials::EnergyFactor(0.0, 0.0);
    VanDerWaalsTailCorrection = Potentials::EnergyFactor(0.0, 0.0);
    CoulombicReal = Potentials::EnergyFactor(0.0, 0.0);
    CoulombicFourier = Potentials::EnergyFactor(0.0, 0.0);
    totalInter = Potentials::EnergyFactor(0.0, 0.0);
  }

  void sumTotal() { totalInter = VanDerWaals + VanDerWaalsTailCorrection + CoulombicReal + CoulombicFourier; }

  inline Potentials::EnergyFactor totalEnergyFactor() const
  {
    return VanDerWaals + VanDerWaalsTailCorrection + CoulombicReal + CoulombicFourier;
  }

  inline Potentials::EnergyFactor total() const
  {
    return VanDerWaals + VanDerWaalsTailCorrection + CoulombicReal + CoulombicFourier;
  }

  inline EnergyInter& operator+=(const EnergyInter& b)
  {
    this->VanDerWaals += b.VanDerWaals;
    this->VanDerWaalsTailCorrection += b.VanDerWaalsTailCorrection;
    this->CoulombicReal += b.CoulombicReal;
    this->CoulombicFourier += b.CoulombicFourier;
    this->totalInter += b.totalInter;
    return *this;
  }

  inline EnergyInter& operator-=(const EnergyInter& b)
  {
    this->VanDerWaals -= b.VanDerWaals;
    this->VanDerWaalsTailCorrection -= b.VanDerWaalsTailCorrection;
    this->CoulombicReal -= b.CoulombicReal;
    this->CoulombicFourier -= b.CoulombicFourier;
    this->totalInter -= b.totalInter;
    return *this;
  }

  inline EnergyInter operator-() const
  {
    EnergyInter v;
    v.VanDerWaals = -VanDerWaals;
    v.VanDerWaalsTailCorrection = -VanDerWaalsTailCorrection;
    v.CoulombicReal = -CoulombicReal;
    v.CoulombicFourier = -CoulombicFourier;
    v.totalInter = -totalInter;
    return v;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnergyInter& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnergyInter& e);
};

export inline EnergyInter operator+(const EnergyInter& a, const EnergyInter& b)
{
  EnergyInter m{};
  m.VanDerWaals = a.VanDerWaals + b.VanDerWaals;
  m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection + b.VanDerWaalsTailCorrection;
  m.CoulombicReal = a.CoulombicReal + b.CoulombicReal;
  m.CoulombicFourier = a.CoulombicFourier + b.CoulombicFourier;
  m.totalInter = a.totalInter + b.totalInter;
  return m;
}

export inline EnergyInter operator-(const EnergyInter& a, const EnergyInter& b)
{
  EnergyInter m{};
  m.VanDerWaals = a.VanDerWaals - b.VanDerWaals;
  m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection - b.VanDerWaalsTailCorrection;
  m.CoulombicReal = a.CoulombicReal - b.CoulombicReal;
  m.CoulombicFourier = a.CoulombicFourier - b.CoulombicFourier;
  m.totalInter = a.totalInter - b.totalInter;
  return m;
}

export inline EnergyInter operator*(const EnergyInter& a, const EnergyInter& b)
{
  EnergyInter m{};
  m.VanDerWaals = a.VanDerWaals * b.VanDerWaals;
  m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection * b.VanDerWaalsTailCorrection;
  m.CoulombicReal = a.CoulombicReal * b.CoulombicReal;
  m.CoulombicFourier = a.CoulombicFourier * b.CoulombicFourier;
  m.totalInter = a.totalInter * b.totalInter;
  return m;
}

export inline EnergyInter operator*(const double& a, const EnergyInter& b)
{
  EnergyInter m{};
  m.VanDerWaals = a * b.VanDerWaals;
  m.VanDerWaalsTailCorrection = a * b.VanDerWaalsTailCorrection;
  m.CoulombicReal = a * b.CoulombicReal;
  m.CoulombicFourier = a * b.CoulombicFourier;
  m.totalInter = a * b.totalInter;
  return m;
}

export inline EnergyInter operator/(const EnergyInter& a, const double& b)
{
  EnergyInter m{};
  m.VanDerWaals = a.VanDerWaals / b;
  m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection / b;
  m.CoulombicReal = a.CoulombicReal / b;
  m.CoulombicFourier = a.CoulombicFourier / b;
  m.totalInter = a.totalInter / b;
  return m;
}

export inline EnergyInter sqrt(const EnergyInter& a)
{
  EnergyInter m{};
  m.VanDerWaals = sqrt(a.VanDerWaals);
  m.VanDerWaalsTailCorrection = sqrt(a.VanDerWaalsTailCorrection);
  m.CoulombicReal = sqrt(a.CoulombicReal);
  m.CoulombicFourier = sqrt(a.CoulombicFourier);
  m.totalInter = sqrt(a.totalInter);
  return m;
}
