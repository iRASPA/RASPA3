export module energy_status;

import <string>;
import <vector>;
import <cmath>;
import <algorithm>;
import <iostream>;
import <numeric>;
import <string>;
import <sstream>;
import <fstream>;

import archive;
import energy_factor;
import energy_status_intra;
import energy_status_inter;
import component;


export struct EnergyStatus
{
  EnergyStatus(): totalEnergy(0.0, 0.0) {};

  EnergyStatus(size_t size) : size(size), totalEnergy(0.0, 0.0), intraEnergy({}), interEnergy({}),
      intraEnergies(std::vector<EnergyIntra>(size)),
      interEnergies(std::vector<EnergyInter>(size * size)),
      dUdlambda(0.0)
  {
  }

  bool operator==(EnergyStatus const&) const = default;

  void resize(size_t numberOfComponents)
  {
    size = numberOfComponents;
    intraEnergies.resize(size);
    interEnergies.resize(size * size);
  }

  inline EnergyIntra& operator() (size_t compA) { return intraEnergies[compA]; }
  inline EnergyInter& operator() (size_t compA, size_t compB) { return interEnergies[compA * size + compB]; }


  EnergyFactor interEnergyComponent(size_t compA) 
  {
    EnergyFactor sum(0.0, 0.0);
    for(size_t i=0; i < size; i++)
    {
        sum += interEnergies[compA * size + i].totalInter;
    }
    return sum;
  }

  EnergyFactor interEnergyComponentZero() {return interEnergies[0 * size + 0].VanDerWaals + 
                                                  interEnergies[0 * size + 1].VanDerWaals +
                                                  interEnergies[0 * size + 2].VanDerWaals;}
  EnergyFactor interEnergyComponentOne() {return interEnergies[1 * size + 0].VanDerWaals + 
                                                 interEnergies[1 * size + 1].VanDerWaals +
                                                 interEnergies[1 * size + 2].VanDerWaals;}
  EnergyFactor interEnergyComponentTwo() {return interEnergies[2 * size + 0].VanDerWaals + 
                                                 interEnergies[2 * size + 1].VanDerWaals +
                                                 interEnergies[2 * size + 2].VanDerWaals;}

  void zero()
  {
    totalEnergy = EnergyFactor(0.0, 0.0);
    dUdlambda = 0.0;
    intraEnergy.zero();
    interEnergy.zero();
    std::fill(intraEnergies.begin(), intraEnergies.end(), EnergyIntra());
    std::fill(interEnergies.begin(), interEnergies.end(), EnergyInter());
  }



  void sumTotal()
  {
    for (size_t i = 0; i < this->intraEnergies.size(); ++i)
    {
        intraEnergy += intraEnergies[i];
    }
    for (size_t i = 0; i < this->interEnergies.size(); ++i)
    {
        interEnergies[i].sumTotal();
        interEnergy += interEnergies[i];
    }
    totalEnergy = intraEnergy.total() + interEnergy.total();
  }

  std::string printEnergyStatus(const std::vector<Component> &components, const std::string& label);
  
  inline EnergyStatus& operator+=(const EnergyStatus& b)
  {
    totalEnergy += b.totalEnergy;
    dUdlambda += b.dUdlambda;
    intraEnergy += b.intraEnergy;
    interEnergy += b.interEnergy;
    for (size_t i = 0; i < this->intraEnergies.size(); ++i)
    {
        intraEnergies[i] += b.intraEnergies[i];
    }
    for (size_t i = 0; i < this->interEnergies.size(); ++i)
    {
        interEnergies[i] += b.interEnergies[i];
    }
    
    return *this;
  }

  inline EnergyStatus& operator-=(const EnergyStatus& b)
  {
    totalEnergy -= b.totalEnergy;
    dUdlambda -= b.dUdlambda;
    intraEnergy -= b.intraEnergy;
    interEnergy -= b.interEnergy;
    for (size_t i = 0; i < this->intraEnergies.size(); ++i)
    {
        intraEnergies[i] -= b.intraEnergies[i];
    }
    for (size_t i = 0; i < this->interEnergies.size(); ++i)
    {
        interEnergies[i] -= b.interEnergies[i];
    }

    return *this;
  }

  inline EnergyStatus operator-() const
  {
    EnergyStatus v(size);
    v.totalEnergy = -totalEnergy;
    v.dUdlambda = -dUdlambda;
    v.intraEnergy = -intraEnergy;
    v.interEnergy = -interEnergy;
    for (size_t i = 0; i < this->intraEnergies.size(); ++i)
    {
        v.intraEnergies[i] = -intraEnergies[i];
    }
    for (size_t i = 0; i < this->interEnergies.size(); ++i)
    {
        v.interEnergies[i] = -interEnergies[i];
    }

    return v;
  }

  uint64_t versionNumber{ 1 };
  size_t size;
  EnergyFactor totalEnergy;
  EnergyIntra intraEnergy;
  EnergyInter interEnergy;
  std::vector < EnergyIntra > intraEnergies;
  std::vector < EnergyInter > interEnergies;
  double dUdlambda;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyStatus &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyStatus &e);
};


export inline EnergyStatus operator+(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.size);
  m.totalEnergy = a.totalEnergy + b.totalEnergy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.intraEnergy = a.intraEnergy + b.intraEnergy;
  m.interEnergy = a.interEnergy + b.interEnergy;
  for (size_t i = 0; i < a.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = a.intraEnergies[i] + b.intraEnergies[i];
  }
  for (size_t i = 0; i < a.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = a.interEnergies[i] + b.interEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator-(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.size);
  m.totalEnergy = a.totalEnergy - b.totalEnergy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.intraEnergy = a.intraEnergy - b.intraEnergy;
  m.interEnergy = a.interEnergy - b.interEnergy;
  for (size_t i = 0; i < a.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = a.intraEnergies[i] - b.intraEnergies[i];
  }
  for (size_t i = 0; i < a.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = a.interEnergies[i] - b.interEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator*(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.size);
  m.totalEnergy = a.totalEnergy * b.totalEnergy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.intraEnergy = a.intraEnergy * b.intraEnergy;
  m.interEnergy = a.interEnergy * b.interEnergy;
  for (size_t i = 0; i < a.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = a.intraEnergies[i] * b.intraEnergies[i];
  }
  for (size_t i = 0; i < a.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = a.interEnergies[i] * b.interEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator*(const double& a, const EnergyStatus& b)
{
  EnergyStatus m(b.size);
  m.totalEnergy = a * b.totalEnergy;
  m.dUdlambda = a * b.dUdlambda;
  m.intraEnergy = a * b.intraEnergy;
  m.interEnergy = a * b.interEnergy;
  for (size_t i = 0; i < b.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = a * b.intraEnergies[i];
  }
  for (size_t i = 0; i < b.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = a * b.interEnergies[i];
  }

  return m;
}


export inline EnergyStatus operator/(const EnergyStatus& a, const double& b)
{
  EnergyStatus m(a.size);
  m.totalEnergy = a.totalEnergy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.intraEnergy = a.intraEnergy / b;
  m.interEnergy = a.interEnergy / b;
  for (size_t i = 0; i < a.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = a.intraEnergies[i] / b;
  }
  for (size_t i = 0; i < a.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = a.interEnergies[i] / b;
  }

  return m;
}


export inline EnergyStatus sqrt(const EnergyStatus& a)
{
  EnergyStatus m(a.size);
  m.totalEnergy = sqrt(a.totalEnergy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.intraEnergy = sqrt(a.intraEnergy);
  m.interEnergy = sqrt(a.interEnergy);
  for (size_t i = 0; i < a.intraEnergies.size(); ++i)
  {
      m.intraEnergies[i] = sqrt(a.intraEnergies[i]);
  }
  for (size_t i = 0; i < a.interEnergies.size(); ++i)
  {
      m.interEnergies[i] = sqrt(a.interEnergies[i]);
  }

  return m;
}

