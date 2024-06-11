module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#endif

export module running_energy;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <map>;
import <cmath>;
import <algorithm>;
import <iostream>;
import <numeric>;
import <string>;
import <sstream>;
import <ostream>;
import <fstream>;
import <functional>;
#endif

import archive;
import scaling;
import hdf5;

export struct RunningEnergy
{
  RunningEnergy()
      : externalFieldVDW(0.0),
        frameworkMoleculeVDW(0.0),
        moleculeMoleculeVDW(0.0),
        externalFieldCharge(0.0),
        frameworkMoleculeCharge(0.0),
        moleculeMoleculeCharge(0.0),
        ewald(0.0),
        intraVDW(0.0),
        intraCoul(0.0),
        tail(0.0),
        polarization(0.0),
        dudlambdaVDW(0.0),
        dudlambdaCharge(0.0),
        dudlambdaEwald(0.0),
        translationalKineticEnergy(0.0),
        rotationalKineticEnergy(0.0)
  {
  }

  std::string printMC() const;
  std::string printMD() const;
  std::string printMC(const std::string& label) const;
  std::string printMD(const std::string& label, double referenceEnergy) const;

  void logMC(HDF5Writer& hdf5, const std::string& label) const;
  void logMD(HDF5Writer& hdf5, const std::string& label, double referenceEnergy) const;

  inline double potentialEnergy() const
  {
    return externalFieldVDW + frameworkMoleculeVDW + moleculeMoleculeVDW + externalFieldCharge +
           frameworkMoleculeCharge + moleculeMoleculeCharge + ewald + intraVDW + intraCoul + tail + polarization;
  }

  inline double kineticEnergy() const { return translationalKineticEnergy + rotationalKineticEnergy; }

  inline double conservedEnergy() const
  {
    return externalFieldVDW + frameworkMoleculeVDW + moleculeMoleculeVDW + externalFieldCharge +
           frameworkMoleculeCharge + moleculeMoleculeCharge + ewald + intraVDW + intraCoul + tail + polarization +
           translationalKineticEnergy + rotationalKineticEnergy;
  }

  inline double dudlambda(double lambda) const
  {
    return Scaling::scalingVDWDerivative(lambda) * dudlambdaVDW +
           Scaling::scalingCoulombDerivative(lambda) * (dudlambdaCharge + dudlambdaEwald);
  }

  inline void zero()
  {
    externalFieldVDW = 0.0;
    frameworkMoleculeVDW = 0.0;
    moleculeMoleculeVDW = 0.0;
    externalFieldCharge = 0.0;
    frameworkMoleculeCharge = 0.0;
    moleculeMoleculeCharge = 0.0;
    ewald = 0.0;
    intraVDW = 0.0;
    intraCoul = 0.0;
    tail = 0.0;
    polarization = 0.0;
    dudlambdaVDW = 0.0;
    dudlambdaCharge = 0.0;
    dudlambdaEwald = 0.0;
    translationalKineticEnergy = 0.0;
    rotationalKineticEnergy = 0.0;
  }

  inline RunningEnergy& operator+=(const RunningEnergy& b)
  {
    externalFieldVDW += b.externalFieldVDW;
    frameworkMoleculeVDW += b.frameworkMoleculeVDW;
    moleculeMoleculeVDW += b.moleculeMoleculeVDW;
    externalFieldCharge += b.externalFieldCharge;
    frameworkMoleculeCharge += b.frameworkMoleculeCharge;
    moleculeMoleculeCharge += b.moleculeMoleculeCharge;
    ewald += b.ewald;
    intraVDW += b.intraVDW;
    intraCoul += b.intraCoul;
    tail += b.tail;
    polarization += b.polarization;
    dudlambdaVDW += b.dudlambdaVDW;
    dudlambdaCharge += b.dudlambdaCharge;
    dudlambdaEwald += b.dudlambdaEwald;
    translationalKineticEnergy += b.translationalKineticEnergy;
    rotationalKineticEnergy += b.rotationalKineticEnergy;

    return *this;
  }

  inline RunningEnergy& operator-=(const RunningEnergy& b)
  {
    externalFieldVDW -= b.externalFieldVDW;
    frameworkMoleculeVDW -= b.frameworkMoleculeVDW;
    moleculeMoleculeVDW -= b.moleculeMoleculeVDW;
    externalFieldCharge -= b.externalFieldCharge;
    frameworkMoleculeCharge -= b.frameworkMoleculeCharge;
    moleculeMoleculeCharge -= b.moleculeMoleculeCharge;
    ewald -= b.ewald;
    intraVDW -= b.intraVDW;
    intraCoul -= b.intraCoul;
    tail -= b.tail;
    polarization -= b.polarization;
    dudlambdaVDW -= b.dudlambdaVDW;
    dudlambdaCharge -= b.dudlambdaCharge;
    dudlambdaEwald -= b.dudlambdaEwald;
    translationalKineticEnergy -= b.translationalKineticEnergy;
    rotationalKineticEnergy -= b.rotationalKineticEnergy;

    return *this;
  }

  inline RunningEnergy operator-() const
  {
    RunningEnergy v{};
    v.externalFieldVDW = -externalFieldVDW;
    v.frameworkMoleculeVDW = -frameworkMoleculeVDW;
    v.moleculeMoleculeVDW = -moleculeMoleculeVDW;
    v.externalFieldCharge = -externalFieldCharge;
    v.frameworkMoleculeCharge = -frameworkMoleculeCharge;
    v.moleculeMoleculeCharge = -moleculeMoleculeCharge;
    v.ewald = -ewald;
    v.intraVDW = -intraVDW;
    v.intraCoul = -intraCoul;
    v.tail = -tail;
    v.polarization = -polarization;
    v.dudlambdaVDW = -dudlambdaVDW;
    v.dudlambdaCharge = -dudlambdaCharge;
    v.dudlambdaEwald = -dudlambdaEwald;
    v.translationalKineticEnergy = -translationalKineticEnergy;
    v.rotationalKineticEnergy = -rotationalKineticEnergy;

    return v;
  }

  uint64_t versionNumber{1};

  double externalFieldVDW;
  double frameworkMoleculeVDW;
  double moleculeMoleculeVDW;
  double externalFieldCharge;
  double frameworkMoleculeCharge;
  double moleculeMoleculeCharge;
  double ewald;
  double intraVDW;
  double intraCoul;
  double tail;
  double polarization;
  double dudlambdaVDW;
  double dudlambdaCharge;
  double dudlambdaEwald;
  double translationalKineticEnergy;
  double rotationalKineticEnergy;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const RunningEnergy& c);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, RunningEnergy& c);
};

export inline RunningEnergy operator+(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a.externalFieldVDW + b.externalFieldVDW;
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW + b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW + b.moleculeMoleculeVDW;
  m.externalFieldCharge = a.externalFieldCharge + b.externalFieldCharge;
  m.frameworkMoleculeCharge = a.frameworkMoleculeCharge + b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a.moleculeMoleculeCharge + b.moleculeMoleculeCharge;
  m.ewald = a.ewald + b.ewald;
  m.intraVDW = a.intraVDW + b.intraVDW;
  m.intraCoul = a.intraCoul + b.intraCoul;
  m.tail = a.tail + b.tail;
  m.polarization = a.polarization + b.polarization;
  m.dudlambdaVDW = a.dudlambdaVDW + b.dudlambdaVDW;
  m.dudlambdaCharge = a.dudlambdaCharge + b.dudlambdaCharge;
  m.dudlambdaEwald = a.dudlambdaEwald + b.dudlambdaEwald;
  m.translationalKineticEnergy = a.translationalKineticEnergy + b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a.rotationalKineticEnergy + b.rotationalKineticEnergy;

  return m;
}

export inline RunningEnergy operator-(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a.externalFieldVDW - b.externalFieldVDW;
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW - b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW - b.moleculeMoleculeVDW;
  m.externalFieldCharge = a.externalFieldCharge - b.externalFieldCharge;
  m.frameworkMoleculeCharge = a.frameworkMoleculeCharge - b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a.moleculeMoleculeCharge - b.moleculeMoleculeCharge;
  m.ewald = a.ewald - b.ewald;
  m.intraVDW = a.intraVDW - b.intraVDW;
  m.intraCoul = a.intraCoul - b.intraCoul;
  m.tail = a.tail - b.tail;
  m.polarization = a.polarization - b.polarization;
  m.dudlambdaVDW = a.dudlambdaVDW - b.dudlambdaVDW;
  m.dudlambdaCharge = a.dudlambdaCharge - b.dudlambdaCharge;
  m.dudlambdaEwald = a.dudlambdaEwald - b.dudlambdaEwald;
  m.translationalKineticEnergy = a.translationalKineticEnergy - b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a.rotationalKineticEnergy - b.rotationalKineticEnergy;
  return m;
}

export inline RunningEnergy operator*(double a, const RunningEnergy b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a * b.externalFieldVDW;
  m.frameworkMoleculeVDW = a * b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a * b.moleculeMoleculeVDW;
  m.externalFieldCharge = a * b.externalFieldCharge;
  m.frameworkMoleculeCharge = a * b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a * b.moleculeMoleculeCharge;
  m.ewald = a * b.ewald;
  m.intraVDW = a * b.intraVDW;
  m.intraCoul = a * b.intraCoul;
  m.tail = a * b.tail;
  m.polarization = a * b.polarization;
  m.dudlambdaVDW = a * b.dudlambdaVDW;
  m.dudlambdaCharge = a * b.dudlambdaCharge;
  m.dudlambdaEwald = a * b.dudlambdaEwald;
  m.translationalKineticEnergy = a * b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a * b.rotationalKineticEnergy;

  return m;
}

export inline RunningEnergy operator*(const RunningEnergy a, double b)
{
  RunningEnergy m{};
  m.externalFieldVDW = b * a.externalFieldVDW;
  m.frameworkMoleculeVDW = b * a.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = b * a.moleculeMoleculeVDW;
  m.externalFieldCharge = b * a.externalFieldCharge;
  m.frameworkMoleculeCharge = b * a.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = b * a.moleculeMoleculeCharge;
  m.ewald = b * a.ewald;
  m.intraVDW = b * a.intraVDW;
  m.intraCoul = b * a.intraCoul;
  m.tail = b * a.tail;
  m.polarization = b * a.polarization;
  m.dudlambdaVDW = b * a.dudlambdaVDW;
  m.dudlambdaCharge = b * a.dudlambdaCharge;
  m.dudlambdaEwald = b * a.dudlambdaEwald;
  m.translationalKineticEnergy = b * a.translationalKineticEnergy;
  m.rotationalKineticEnergy = b * a.rotationalKineticEnergy;

  return m;
}
