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

export module energy_status_intra;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;
import energy_factor;

export struct EnergyIntra
{
  std::uint64_t versionNumber{1};

  double bond;
  double bend;
  double inversionBend;
  double ureyBradley;
  double torsion;
  double improperTorsion;
  double bondBond;
  double bondBend;
  double bondTorsion;
  double bendBend;
  double bendTorsion;
  double intraVDW;
  double intraChargeCharge;

  EnergyIntra()
      : bond(0.0),
        bend(0.0),
        inversionBend(0.0),
        ureyBradley(0.0),
        torsion(0.0),
        improperTorsion(0.0),
        bondBond(0.0),
        bondBend(0.0),
        bondTorsion(0.0),
        bendBend(0.0),
        bendTorsion(0.0),
        intraVDW(0.0),
        intraChargeCharge(0.0)
  {
  }

  bool operator==(EnergyIntra const&) const = default;

  inline Potentials::EnergyFactor total() const
  {
    return Potentials::EnergyFactor(bond + bend + inversionBend + ureyBradley + torsion + improperTorsion + bondBond +
                                        bondBend + bondTorsion + bendBend + bendTorsion + intraVDW + intraChargeCharge,
                                    0.0);
  }

  void zero()
  {
    this->bond = 0.0;
    this->bend = 0.0;
    this->inversionBend = 0.0;
    this->ureyBradley = 0.0;
    this->torsion = 0.0;
    this->improperTorsion = 0.0;
    this->bondBond = 0.0;
    this->bondBend = 0.0;
    this->bondTorsion = 0.0;
    this->bendBend = 0.0;
    this->bendTorsion = 0.0;
    this->intraVDW = 0.0;
    this->intraChargeCharge = 0.0;
  }

  inline EnergyIntra& operator+=(const EnergyIntra& b)
  {
    this->bond += b.bond;
    this->bend += b.bend;
    this->inversionBend += b.inversionBend;
    this->ureyBradley += b.ureyBradley;
    this->torsion += b.torsion;
    this->improperTorsion += b.improperTorsion;
    this->bondBond += b.bondBond;
    this->bondBend += b.bondBend;
    this->bondTorsion += b.bondTorsion;
    this->bendBend += b.bendBend;
    this->bendTorsion += b.bendTorsion;
    this->intraVDW += b.intraVDW;
    this->intraChargeCharge += b.intraChargeCharge;
    return *this;
  }

  inline EnergyIntra& operator-=(const EnergyIntra& b)
  {
    this->bond -= b.bond;
    this->bend -= b.bend;
    this->inversionBend -= b.inversionBend;
    this->ureyBradley -= b.ureyBradley;
    this->torsion -= b.torsion;
    this->improperTorsion -= b.improperTorsion;
    this->bondBond -= b.bondBond;
    this->bondBend -= b.bondBend;
    this->bondTorsion -= b.bondTorsion;
    this->bendBend -= b.bendBend;
    this->bendTorsion -= b.bendTorsion;
    this->intraVDW -= b.intraVDW;
    this->intraChargeCharge -= b.intraChargeCharge;
    return *this;
  }

  inline EnergyIntra operator-() const
  {
    EnergyIntra v;
    v.bond = -bond;
    v.bend = -bend;
    v.inversionBend = -inversionBend;
    v.ureyBradley = -ureyBradley;
    v.torsion = -torsion;
    v.improperTorsion = -improperTorsion;
    v.bondBond = -bondBond;
    v.bondBend = -bondBend;
    v.bondTorsion = -bondTorsion;
    v.bendBend = -bendBend;
    v.bendTorsion = -bendTorsion;
    v.intraVDW = -intraVDW;
    v.intraChargeCharge = -intraChargeCharge;
    return v;
  }

  std::string printEnergyStatus([[maybe_unused]] int i)
  {
    std::ostringstream stream;

    std::print(stream, "    bond: {}", bond);

    return stream.str();
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnergyIntra& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnergyIntra& e);
};

export inline EnergyIntra operator+(const EnergyIntra& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a.bond + b.bond;
  m.bend = a.bend + b.bend;
  m.inversionBend = a.inversionBend + b.inversionBend;
  m.ureyBradley = a.ureyBradley + b.ureyBradley;
  m.torsion = a.torsion + b.torsion;
  m.improperTorsion = a.improperTorsion + b.improperTorsion;
  m.bondBond = a.bondBond + b.bondBond;
  m.bondBend = a.bondBend + b.bondBend;
  m.bondTorsion = a.bondTorsion + b.bondTorsion;
  m.bendBend = a.bendBend + b.bendBend;
  m.bendTorsion = a.bendTorsion + b.bendTorsion;
  m.intraVDW = a.intraVDW + b.intraVDW;
  m.intraChargeCharge = a.intraChargeCharge + b.intraChargeCharge;
  return m;
}

export inline EnergyIntra operator-(const EnergyIntra& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a.bond - b.bond;
  m.bend = a.bend - b.bend;
  m.inversionBend = a.inversionBend - b.inversionBend;
  m.ureyBradley = a.ureyBradley - b.ureyBradley;
  m.torsion = a.torsion - b.torsion;
  m.improperTorsion = a.improperTorsion - b.improperTorsion;
  m.bondBond = a.bondBond - b.bondBond;
  m.bondBend = a.bondBend - b.bondBend;
  m.bondTorsion = a.bondTorsion - b.bondTorsion;
  m.bendBend = a.bendBend - b.bendBend;
  m.bendTorsion = a.bendTorsion - b.bendTorsion;
  m.intraVDW = a.intraVDW - b.intraVDW;
  m.intraChargeCharge = a.intraChargeCharge - b.intraChargeCharge;
  return m;
}

export inline EnergyIntra operator*(const EnergyIntra& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a.bond * b.bond;
  m.bend = a.bend * b.bend;
  m.inversionBend = a.inversionBend * b.inversionBend;
  m.ureyBradley = a.ureyBradley * b.ureyBradley;
  m.torsion = a.torsion * b.torsion;
  m.improperTorsion = a.improperTorsion * b.improperTorsion;
  m.bondBond = a.bondBond * b.bondBond;
  m.bondBend = a.bondBend * b.bondBend;
  m.bondTorsion = a.bondTorsion * b.bondTorsion;
  m.bendBend = a.bendBend * b.bendBend;
  m.bendTorsion = a.bendTorsion * b.bendTorsion;
  m.intraVDW = a.intraVDW * b.intraVDW;
  m.intraChargeCharge = a.intraChargeCharge * b.intraChargeCharge;
  return m;
}

export inline EnergyIntra operator*(const double& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a * b.bond;
  m.bend = a * b.bend;
  m.inversionBend = a * b.inversionBend;
  m.ureyBradley = a * b.ureyBradley;
  m.torsion = a * b.torsion;
  m.improperTorsion = a * b.improperTorsion;
  m.bondBond = a * b.bondBond;
  m.bondBend = a * b.bondBend;
  m.bondTorsion = a * b.bondTorsion;
  m.bendBend = a * b.bendBend;
  m.bendTorsion = a * b.bendTorsion;
  m.intraVDW = a * b.intraVDW;
  m.intraChargeCharge = a * b.intraChargeCharge;
  return m;
}

export inline EnergyIntra operator/(const EnergyIntra& a, const double& b)
{
  EnergyIntra m{};
  m.bond = a.bond / b;
  m.bend = a.bend / b;
  m.inversionBend = a.inversionBend / b;
  m.ureyBradley = a.ureyBradley / b;
  m.torsion = a.torsion / b;
  m.improperTorsion = a.improperTorsion / b;
  m.bondBond = a.bondBond / b;
  m.bondBend = a.bondBend / b;
  m.bondTorsion = a.bondTorsion / b;
  m.bendBend = a.bendBend / b;
  m.bendTorsion = a.bendTorsion / b;
  m.intraVDW = a.intraVDW / b;
  m.intraChargeCharge = a.intraChargeCharge / b;
  return m;
}

export inline EnergyIntra sqrt(const EnergyIntra& a)
{
  EnergyIntra m{};
  m.bond = std::sqrt(a.bond);
  m.bend = std::sqrt(a.bend);
  m.inversionBend = std::sqrt(a.inversionBend);
  m.ureyBradley = std::sqrt(a.ureyBradley);
  m.torsion = std::sqrt(a.torsion);
  m.improperTorsion = std::sqrt(a.improperTorsion);
  m.bondBond = std::sqrt(a.bondBond);
  m.bondBend = std::sqrt(a.bondBend);
  m.bondTorsion = std::sqrt(a.bondTorsion);
  m.bendBend = std::sqrt(a.bendBend);
  m.bendTorsion = std::sqrt(a.bendTorsion);
  m.intraVDW = std::sqrt(a.intraVDW);
  m.intraChargeCharge = std::sqrt(a.intraChargeCharge);
  return m;
}
