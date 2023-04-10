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
};


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
