export module running_energy;

import <string>;
import <vector>;
import <cmath>;
import <algorithm>;
import <iostream>;
import <numeric>;
import <string>;
import <sstream>;
import <ostream>;
import <fstream>;

import archive;
import scaling;


export struct RunningEnergy
{
  RunningEnergy() : frameworkMoleculeVDW(0.0), moleculeMoleculeVDW(0.0), 
                    frameworkMoleculeCharge(0.0), moleculeMoleculeCharge(0.0),
                    ewald(0.0),
                    intraVDW(0.0), intraCoul(0.0),
                    tail(0.0),
                    polarization(0.0),
                    dudlambdaVDW(0.0),
                    dudlambdaCharge(0.0),
                    dudlambdaEwald(0.0)
  {
  }

  void print(std::ostream &stream, const std::string& label);
  
  inline double total() const
  {
      return frameworkMoleculeVDW + moleculeMoleculeVDW + 
             frameworkMoleculeCharge + moleculeMoleculeCharge +
             ewald + 
             intraVDW + intraCoul + 
             tail + 
             polarization;
  }

  inline double dudlambda(double lambda) const
  {
    return Scaling::scalingVDWDerivative(lambda) * dudlambdaVDW + 
           Scaling::scalingCoulombDerivative(lambda) * (dudlambdaCharge + dudlambdaEwald);
  }

  inline void zero()
  {
    frameworkMoleculeVDW = 0.0;
    moleculeMoleculeVDW = 0.0;
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
  }

  inline RunningEnergy& operator+=(const RunningEnergy& b)
  {
    frameworkMoleculeVDW += b.frameworkMoleculeVDW;
    moleculeMoleculeVDW += b.moleculeMoleculeVDW;
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
    
    return *this;
  }

  inline RunningEnergy& operator-=(const RunningEnergy& b)
  {
    frameworkMoleculeVDW -= b.frameworkMoleculeVDW;
    moleculeMoleculeVDW -= b.moleculeMoleculeVDW;
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

    return *this;
  }

  inline RunningEnergy operator-() const
  {
    RunningEnergy v{};
    v.frameworkMoleculeVDW = -frameworkMoleculeVDW;
    v.moleculeMoleculeVDW = -moleculeMoleculeVDW;
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

    return v;
  }

  uint64_t versionNumber{ 1 };

  double frameworkMoleculeVDW;
  double moleculeMoleculeVDW;
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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RunningEnergy &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RunningEnergy &c);
};

export inline bool operator==(const RunningEnergy& a, const RunningEnergy& b) 
{
  return (std::abs(b.frameworkMoleculeVDW - a.frameworkMoleculeVDW) < 1e-3) &&
         (std::abs(b.moleculeMoleculeVDW - a.moleculeMoleculeVDW) < 1e-3) &&
         (std::abs(b.moleculeMoleculeCharge - a.moleculeMoleculeCharge) < 1e-3) &&
         (std::abs(b.ewald - a.ewald) < 1e-3) &&
         (std::abs(b.intraVDW - a.intraVDW) < 1e-3) &&
         (std::abs(b.intraCoul - a.intraCoul) < 1e-3) &&
         (std::abs(b.tail - a.tail) < 1e-3) &&
         (std::abs(b.polarization - a.polarization) < 1e-3) &&
         (std::abs(b.dudlambdaVDW - a.dudlambdaVDW) < 1e-3) &&
         (std::abs(b.dudlambdaCharge - a.dudlambdaCharge) < 1e-3) &&
         (std::abs(b.dudlambdaEwald - a.dudlambdaEwald) < 1e-3);
}


export inline RunningEnergy operator+(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW + b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW + b.moleculeMoleculeVDW;
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

  return m;
}

export inline RunningEnergy operator-(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW - b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW - b.moleculeMoleculeVDW;
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
  return m;
}


export inline RunningEnergy operator*(double a, const RunningEnergy b)
{
  RunningEnergy m{};
  m.frameworkMoleculeVDW = a * b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a * b.moleculeMoleculeVDW;
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

  return m;
}

export inline RunningEnergy operator*(const RunningEnergy a, double b)
{
  RunningEnergy m{};
  m.frameworkMoleculeVDW = b * a.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = b * a.moleculeMoleculeVDW;
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

  return m;
}
