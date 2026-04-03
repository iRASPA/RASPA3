module;

export module average_energy_type;

import std;

import archive;

export struct AverageEnergyType
{
  double totalEnergy;
  double VanDerWaalsEnergy;
  double CoulombEnergy;
  double polarizationEnergy;

  AverageEnergyType& operator+=(const AverageEnergyType& b)
  {
    this->totalEnergy += b.totalEnergy;
    this->VanDerWaalsEnergy += b.VanDerWaalsEnergy,
    this->CoulombEnergy += b.CoulombEnergy;
    this->polarizationEnergy += b.polarizationEnergy;
    return *this;
  };

  AverageEnergyType& operator-=(const AverageEnergyType& b)
  {
    this->totalEnergy -= b.totalEnergy;
    this->VanDerWaalsEnergy -= b.VanDerWaalsEnergy;
    this->CoulombEnergy -= b.CoulombEnergy;
    this->polarizationEnergy -= b.polarizationEnergy;
    return *this;
  };

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const AverageEnergyType& vec);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, AverageEnergyType& vec);
};

export inline AverageEnergyType operator+(const AverageEnergyType& a, const AverageEnergyType& b)
{
  return AverageEnergyType(a.totalEnergy + b.totalEnergy, a.VanDerWaalsEnergy + b.VanDerWaalsEnergy, 
                           a.CoulombEnergy + b.CoulombEnergy, a.polarizationEnergy + b.polarizationEnergy);
}

export inline AverageEnergyType operator-(const AverageEnergyType& a, const AverageEnergyType& b)
{
  return AverageEnergyType(a.totalEnergy - b.totalEnergy, a.VanDerWaalsEnergy - b.VanDerWaalsEnergy, 
                           a.CoulombEnergy - b.CoulombEnergy, a.polarizationEnergy - b.polarizationEnergy);
}

export inline AverageEnergyType operator*(const AverageEnergyType& a, const AverageEnergyType& b)
{
  return AverageEnergyType(a.totalEnergy * b.totalEnergy, a.VanDerWaalsEnergy * b.VanDerWaalsEnergy, 
                           a.CoulombEnergy * b.CoulombEnergy, a.polarizationEnergy * b.polarizationEnergy);
}

export inline AverageEnergyType operator*(const AverageEnergyType& a, const double& b)
{
  return AverageEnergyType(a.totalEnergy * b, a.VanDerWaalsEnergy * b, a.CoulombEnergy * b, a.polarizationEnergy * b);
}

export inline AverageEnergyType operator*(const double& a, const AverageEnergyType& b)
{
  return AverageEnergyType(a * b.totalEnergy, a * b.VanDerWaalsEnergy, a * b.CoulombEnergy, a * b.polarizationEnergy);
}

export inline AverageEnergyType operator/(const AverageEnergyType& a, double b) 
{ 
  return AverageEnergyType(a.totalEnergy / b, a.VanDerWaalsEnergy / b, a.CoulombEnergy / b, a.polarizationEnergy / b); 
}

export inline AverageEnergyType sqrt(const AverageEnergyType& a)
{
  return AverageEnergyType(std::sqrt(a.totalEnergy), std::sqrt(a.VanDerWaalsEnergy), std::sqrt(a.CoulombEnergy), std::sqrt(a.polarizationEnergy));
}


