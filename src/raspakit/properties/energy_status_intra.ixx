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
  double ureyBradley;
  double bend;
  double inversionBend;
  double outOfPlaneBend;
  double torsion;
  double improperTorsion;
  double bondBond;
  double bondBend;
  double bondTorsion;
  double bendBend;
  double bendTorsion;
  double vanDerWaals;
  double coulomb;

  EnergyIntra()
      : bond(0.0),
        ureyBradley(0.0),
        bend(0.0),
        inversionBend(0.0),
        outOfPlaneBend(0.0),
        torsion(0.0),
        improperTorsion(0.0),
        bondBond(0.0),
        bondBend(0.0),
        bondTorsion(0.0),
        bendBend(0.0),
        bendTorsion(0.0),
        vanDerWaals(0.0),
        coulomb(0.0)
  {
  }

  bool operator==(EnergyIntra const&) const = default;

  inline Potentials::EnergyFactor total() const
  {
    return Potentials::EnergyFactor(bond + ureyBradley + bend + inversionBend + outOfPlaneBend + torsion +
                                        improperTorsion + bondBond + bondBend + bondTorsion + bendBend + bendTorsion +
                                        vanDerWaals + coulomb,
                                    0.0);
  }

  void zero()
  {
    this->bond = 0.0;
    this->ureyBradley = 0.0;
    this->bend = 0.0;
    this->inversionBend = 0.0;
    this->outOfPlaneBend = 0.0;
    this->torsion = 0.0;
    this->improperTorsion = 0.0;
    this->bondBond = 0.0;
    this->bondBend = 0.0;
    this->bondTorsion = 0.0;
    this->bendBend = 0.0;
    this->bendTorsion = 0.0;
    this->vanDerWaals = 0.0;
    this->coulomb = 0.0;
  }

  inline EnergyIntra& operator+=(const EnergyIntra& b)
  {
    this->bond += b.bond;
    this->ureyBradley += b.ureyBradley;
    this->bend += b.bend;
    this->inversionBend += b.inversionBend;
    this->outOfPlaneBend += b.outOfPlaneBend;
    this->torsion += b.torsion;
    this->improperTorsion += b.improperTorsion;
    this->bondBond += b.bondBond;
    this->bondBend += b.bondBend;
    this->bondTorsion += b.bondTorsion;
    this->bendBend += b.bendBend;
    this->bendTorsion += b.bendTorsion;
    this->vanDerWaals += b.vanDerWaals;
    this->coulomb += b.coulomb;
    return *this;
  }

  inline EnergyIntra& operator-=(const EnergyIntra& b)
  {
    this->bond -= b.bond;
    this->ureyBradley -= b.ureyBradley;
    this->bend -= b.bend;
    this->inversionBend -= b.inversionBend;
    this->outOfPlaneBend -= b.outOfPlaneBend;
    this->torsion -= b.torsion;
    this->improperTorsion -= b.improperTorsion;
    this->bondBond -= b.bondBond;
    this->bondBend -= b.bondBend;
    this->bondTorsion -= b.bondTorsion;
    this->bendBend -= b.bendBend;
    this->bendTorsion -= b.bendTorsion;
    this->vanDerWaals -= b.vanDerWaals;
    this->coulomb -= b.coulomb;
    return *this;
  }

  inline EnergyIntra operator-() const
  {
    EnergyIntra v;
    v.bond = -bond;
    v.ureyBradley = -ureyBradley;
    v.bend = -bend;
    v.inversionBend = -inversionBend;
    v.outOfPlaneBend = -outOfPlaneBend;
    v.torsion = -torsion;
    v.improperTorsion = -improperTorsion;
    v.bondBond = -bondBond;
    v.bondBend = -bondBend;
    v.bondTorsion = -bondTorsion;
    v.bendBend = -bendBend;
    v.bendTorsion = -bendTorsion;
    v.vanDerWaals = -vanDerWaals;
    v.coulomb = -coulomb;
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
  m.ureyBradley = a.ureyBradley + b.ureyBradley;
  m.bend = a.bend + b.bend;
  m.inversionBend = a.inversionBend + b.inversionBend;
  m.outOfPlaneBend = a.outOfPlaneBend + b.outOfPlaneBend;
  m.torsion = a.torsion + b.torsion;
  m.improperTorsion = a.improperTorsion + b.improperTorsion;
  m.bondBond = a.bondBond + b.bondBond;
  m.bondBend = a.bondBend + b.bondBend;
  m.bondTorsion = a.bondTorsion + b.bondTorsion;
  m.bendBend = a.bendBend + b.bendBend;
  m.bendTorsion = a.bendTorsion + b.bendTorsion;
  m.vanDerWaals = a.vanDerWaals + b.vanDerWaals;
  m.coulomb = a.coulomb + b.coulomb;
  return m;
}

export inline EnergyIntra operator-(const EnergyIntra& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a.bond - b.bond;
  m.ureyBradley = a.ureyBradley - b.ureyBradley;
  m.bend = a.bend - b.bend;
  m.inversionBend = a.inversionBend - b.inversionBend;
  m.outOfPlaneBend = a.outOfPlaneBend - b.outOfPlaneBend;
  m.torsion = a.torsion - b.torsion;
  m.improperTorsion = a.improperTorsion - b.improperTorsion;
  m.bondBond = a.bondBond - b.bondBond;
  m.bondBend = a.bondBend - b.bondBend;
  m.bondTorsion = a.bondTorsion - b.bondTorsion;
  m.bendBend = a.bendBend - b.bendBend;
  m.bendTorsion = a.bendTorsion - b.bendTorsion;
  m.vanDerWaals = a.vanDerWaals - b.vanDerWaals;
  m.coulomb = a.coulomb - b.coulomb;
  return m;
}

export inline EnergyIntra operator*(const EnergyIntra& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a.bond * b.bond;
  m.ureyBradley = a.ureyBradley * b.ureyBradley;
  m.bend = a.bend * b.bend;
  m.inversionBend = a.inversionBend * b.inversionBend;
  m.outOfPlaneBend = a.outOfPlaneBend * b.outOfPlaneBend;
  m.torsion = a.torsion * b.torsion;
  m.improperTorsion = a.improperTorsion * b.improperTorsion;
  m.bondBond = a.bondBond * b.bondBond;
  m.bondBend = a.bondBend * b.bondBend;
  m.bondTorsion = a.bondTorsion * b.bondTorsion;
  m.bendBend = a.bendBend * b.bendBend;
  m.bendTorsion = a.bendTorsion * b.bendTorsion;
  m.vanDerWaals = a.vanDerWaals * b.vanDerWaals;
  m.coulomb = a.coulomb * b.coulomb;
  return m;
}

export inline EnergyIntra operator*(const double& a, const EnergyIntra& b)
{
  EnergyIntra m{};
  m.bond = a * b.bond;
  m.ureyBradley = a * b.ureyBradley;
  m.bend = a * b.bend;
  m.inversionBend = a * b.inversionBend;
  m.outOfPlaneBend = a * b.outOfPlaneBend;
  m.torsion = a * b.torsion;
  m.improperTorsion = a * b.improperTorsion;
  m.bondBond = a * b.bondBond;
  m.bondBend = a * b.bondBend;
  m.bondTorsion = a * b.bondTorsion;
  m.bendBend = a * b.bendBend;
  m.bendTorsion = a * b.bendTorsion;
  m.vanDerWaals = a * b.vanDerWaals;
  m.coulomb = a * b.coulomb;
  return m;
}

export inline EnergyIntra operator/(const EnergyIntra& a, const double& b)
{
  EnergyIntra m{};
  m.bond = a.bond / b;
  m.ureyBradley = a.ureyBradley / b;
  m.bend = a.bend / b;
  m.inversionBend = a.inversionBend / b;
  m.outOfPlaneBend = a.outOfPlaneBend / b;
  m.torsion = a.torsion / b;
  m.improperTorsion = a.improperTorsion / b;
  m.bondBond = a.bondBond / b;
  m.bondBend = a.bondBend / b;
  m.bondTorsion = a.bondTorsion / b;
  m.bendBend = a.bendBend / b;
  m.bendTorsion = a.bendTorsion / b;
  m.vanDerWaals = a.vanDerWaals / b;
  m.coulomb = a.coulomb / b;
  return m;
}

export inline EnergyIntra sqrt(const EnergyIntra& a)
{
  EnergyIntra m{};
  m.bond = std::sqrt(a.bond);
  m.ureyBradley = std::sqrt(a.ureyBradley);
  m.bend = std::sqrt(a.bend);
  m.inversionBend = std::sqrt(a.inversionBend);
  m.outOfPlaneBend = std::sqrt(a.outOfPlaneBend);
  m.torsion = std::sqrt(a.torsion);
  m.improperTorsion = std::sqrt(a.improperTorsion);
  m.bondBond = std::sqrt(a.bondBond);
  m.bondBend = std::sqrt(a.bondBend);
  m.bondTorsion = std::sqrt(a.bondTorsion);
  m.bendBend = std::sqrt(a.bendBend);
  m.bendTorsion = std::sqrt(a.bendTorsion);
  m.vanDerWaals = std::sqrt(a.vanDerWaals);
  m.coulomb = std::sqrt(a.coulomb);
  return m;
}
