module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#endif

export module energy_status;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <cmath>;
import <algorithm>;
import <iostream>;
import <numeric>;
import <string>;
import <sstream>;
import <fstream>;
#endif

import archive;
import energy_factor;
import energy_status_intra;
import energy_status_inter;
import component;
import json;

export struct EnergyStatus
{
  EnergyStatus() : totalEnergy(0.0, 0.0) {};

  EnergyStatus(size_t numberOfExternalFields, size_t numberOfFrameworks, size_t numberOfComponents)
      : numberOfExternalFields(numberOfExternalFields),
        numberOfFrameworks(numberOfFrameworks),
        numberOfComponents(numberOfComponents),
        totalEnergy(0.0, 0.0),
        intraEnergy({}),
        externalFieldMoleculeEnergy({}),
        frameworkMoleculeEnergy({}),
        interEnergy({}),
        intraComponentEnergies(std::vector<EnergyIntra>(numberOfComponents)),
        externalFieldComponentEnergies(
            std::vector<EnergyInter>(std::max(1uz, numberOfExternalFields) * numberOfComponents)),
        frameworkComponentEnergies(std::vector<EnergyInter>(std::max(1uz, numberOfFrameworks) * numberOfComponents)),
        interComponentEnergies(std::vector<EnergyInter>(numberOfComponents * numberOfComponents)),
        dUdlambda(0.0)
  {
  }

  inline EnergyInter& externalFieldComponentEnergy(size_t compA, size_t compB)
  {
    return externalFieldComponentEnergies[compA * numberOfComponents + compB];
  }

  inline EnergyInter& frameworkComponentEnergy(size_t compA, size_t compB)
  {
    return frameworkComponentEnergies[compA * numberOfComponents + compB];
  }

  inline EnergyInter& componentEnergy(size_t compA, size_t compB)
  {
    return interComponentEnergies[compA * numberOfComponents + compB];
  }

  Potentials::EnergyFactor interEnergyComponent(size_t compA)
  {
    Potentials::EnergyFactor sum(0.0, 0.0);
    for (size_t i = 0; i < numberOfComponents; i++)
    {
      sum += interComponentEnergies[compA * numberOfComponents + i].totalInter;
    }
    return sum;
  }

  void zero()
  {
    totalEnergy = Potentials::EnergyFactor(0.0, 0.0);
    dUdlambda = 0.0;
    intraEnergy.zero();
    externalFieldMoleculeEnergy.zero();
    frameworkMoleculeEnergy.zero();
    interEnergy.zero();
    std::fill(intraComponentEnergies.begin(), intraComponentEnergies.end(), EnergyIntra());
    std::fill(externalFieldComponentEnergies.begin(), externalFieldComponentEnergies.end(), EnergyInter());
    std::fill(frameworkComponentEnergies.begin(), frameworkComponentEnergies.end(), EnergyInter());
    std::fill(interComponentEnergies.begin(), interComponentEnergies.end(), EnergyInter());
  }

  void sumTotal()
  {
    intraEnergy.zero();
    externalFieldMoleculeEnergy.zero();
    frameworkMoleculeEnergy.zero();
    interEnergy.zero();
    for (size_t i = 0; i < this->intraComponentEnergies.size(); ++i)
    {
      intraEnergy += intraComponentEnergies[i];
    }
    for (size_t i = 0; i < this->externalFieldComponentEnergies.size(); ++i)
    {
      externalFieldComponentEnergies[i].sumTotal();
      externalFieldMoleculeEnergy += externalFieldComponentEnergies[i];
    }
    for (size_t i = 0; i < this->frameworkComponentEnergies.size(); ++i)
    {
      frameworkComponentEnergies[i].sumTotal();
      frameworkMoleculeEnergy += frameworkComponentEnergies[i];
    }
    for (size_t i = 0; i < this->interComponentEnergies.size(); ++i)
    {
      interComponentEnergies[i].sumTotal();
      interEnergy += interComponentEnergies[i];
    }
    totalEnergy = intraEnergy.total() + externalFieldMoleculeEnergy.total() + frameworkMoleculeEnergy.total() +
                  interEnergy.total();
  }

  std::string printEnergyStatus(const std::vector<Component>& components, const std::string& label);

  inline EnergyStatus& operator+=(const EnergyStatus& b)
  {
    totalEnergy += b.totalEnergy;
    dUdlambda += b.dUdlambda;
    intraEnergy += b.intraEnergy;
    externalFieldMoleculeEnergy += b.externalFieldMoleculeEnergy;
    frameworkMoleculeEnergy += b.frameworkMoleculeEnergy;
    interEnergy += b.interEnergy;
    for (size_t i = 0; i < this->intraComponentEnergies.size(); ++i)
    {
      intraComponentEnergies[i] += b.intraComponentEnergies[i];
    }
    for (size_t i = 0; i < this->externalFieldComponentEnergies.size(); ++i)
    {
      externalFieldComponentEnergies[i] += b.externalFieldComponentEnergies[i];
    }
    for (size_t i = 0; i < this->frameworkComponentEnergies.size(); ++i)
    {
      frameworkComponentEnergies[i] += b.frameworkComponentEnergies[i];
    }
    for (size_t i = 0; i < this->interComponentEnergies.size(); ++i)
    {
      interComponentEnergies[i] += b.interComponentEnergies[i];
    }

    return *this;
  }

  inline EnergyStatus& operator-=(const EnergyStatus& b)
  {
    totalEnergy -= b.totalEnergy;
    dUdlambda -= b.dUdlambda;
    intraEnergy -= b.intraEnergy;
    externalFieldMoleculeEnergy -= b.externalFieldMoleculeEnergy;
    frameworkMoleculeEnergy -= b.frameworkMoleculeEnergy;
    interEnergy -= b.interEnergy;
    for (size_t i = 0; i < this->intraComponentEnergies.size(); ++i)
    {
      intraComponentEnergies[i] -= b.intraComponentEnergies[i];
    }
    for (size_t i = 0; i < this->externalFieldComponentEnergies.size(); ++i)
    {
      externalFieldComponentEnergies[i] -= b.externalFieldComponentEnergies[i];
    }
    for (size_t i = 0; i < this->frameworkComponentEnergies.size(); ++i)
    {
      frameworkComponentEnergies[i] -= b.frameworkComponentEnergies[i];
    }
    for (size_t i = 0; i < this->interComponentEnergies.size(); ++i)
    {
      interComponentEnergies[i] -= b.interComponentEnergies[i];
    }

    return *this;
  }

  inline EnergyStatus operator-() const
  {
    EnergyStatus v(numberOfExternalFields, numberOfFrameworks, numberOfComponents);
    v.totalEnergy = -totalEnergy;
    v.dUdlambda = -dUdlambda;
    v.intraEnergy = -intraEnergy;
    v.externalFieldMoleculeEnergy = -externalFieldMoleculeEnergy;
    v.frameworkMoleculeEnergy = -frameworkMoleculeEnergy;
    v.interEnergy = -interEnergy;
    for (size_t i = 0; i < this->intraComponentEnergies.size(); ++i)
    {
      v.intraComponentEnergies[i] = -intraComponentEnergies[i];
    }
    for (size_t i = 0; i < this->externalFieldComponentEnergies.size(); ++i)
    {
      v.externalFieldComponentEnergies[i] = -externalFieldComponentEnergies[i];
    }
    for (size_t i = 0; i < this->frameworkComponentEnergies.size(); ++i)
    {
      v.frameworkComponentEnergies[i] = -frameworkComponentEnergies[i];
    }
    for (size_t i = 0; i < this->interComponentEnergies.size(); ++i)
    {
      v.interComponentEnergies[i] = -interComponentEnergies[i];
    }

    return v;
  }

  uint64_t versionNumber{1};
  size_t numberOfExternalFields;
  size_t numberOfFrameworks;
  size_t numberOfComponents;
  Potentials::EnergyFactor totalEnergy;
  EnergyIntra intraEnergy;
  EnergyInter externalFieldMoleculeEnergy;
  EnergyInter frameworkMoleculeEnergy;
  EnergyInter interEnergy;
  std::vector<EnergyIntra> intraComponentEnergies;
  std::vector<EnergyInter> externalFieldComponentEnergies;
  std::vector<EnergyInter> frameworkComponentEnergies;
  std::vector<EnergyInter> interComponentEnergies;
  double dUdlambda;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnergyStatus& e);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnergyStatus& e);
};

export inline EnergyStatus operator+(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.numberOfExternalFields, a.numberOfFrameworks, a.numberOfComponents);
  m.totalEnergy = a.totalEnergy + b.totalEnergy;
  m.dUdlambda = a.dUdlambda + b.dUdlambda;
  m.intraEnergy = a.intraEnergy + b.intraEnergy;
  m.externalFieldMoleculeEnergy = a.externalFieldMoleculeEnergy + b.externalFieldMoleculeEnergy;
  m.frameworkMoleculeEnergy = a.frameworkMoleculeEnergy + b.frameworkMoleculeEnergy;
  m.interEnergy = a.interEnergy + b.interEnergy;
  for (size_t i = 0; i < a.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = a.intraComponentEnergies[i] + b.intraComponentEnergies[i];
  }
  for (size_t i = 0; i < a.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = a.externalFieldComponentEnergies[i] + b.externalFieldComponentEnergies[i];
  }
  for (size_t i = 0; i < a.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = a.frameworkComponentEnergies[i] + b.frameworkComponentEnergies[i];
  }
  for (size_t i = 0; i < a.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = a.interComponentEnergies[i] + b.interComponentEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator-(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.numberOfExternalFields, a.numberOfFrameworks, a.numberOfComponents);
  m.totalEnergy = a.totalEnergy - b.totalEnergy;
  m.dUdlambda = a.dUdlambda - b.dUdlambda;
  m.intraEnergy = a.intraEnergy - b.intraEnergy;
  m.externalFieldMoleculeEnergy = a.externalFieldMoleculeEnergy - b.externalFieldMoleculeEnergy;
  m.frameworkMoleculeEnergy = a.frameworkMoleculeEnergy - b.frameworkMoleculeEnergy;
  m.interEnergy = a.interEnergy - b.interEnergy;
  for (size_t i = 0; i < a.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = a.intraComponentEnergies[i] - b.intraComponentEnergies[i];
  }
  for (size_t i = 0; i < a.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = a.externalFieldComponentEnergies[i] - b.externalFieldComponentEnergies[i];
  }
  for (size_t i = 0; i < a.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = a.frameworkComponentEnergies[i] - b.frameworkComponentEnergies[i];
  }
  for (size_t i = 0; i < a.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = a.interComponentEnergies[i] - b.interComponentEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator*(const EnergyStatus& a, const EnergyStatus& b)
{
  EnergyStatus m(a.numberOfExternalFields, a.numberOfFrameworks, a.numberOfComponents);
  m.totalEnergy = a.totalEnergy * b.totalEnergy;
  m.dUdlambda = a.dUdlambda * b.dUdlambda;
  m.intraEnergy = a.intraEnergy * b.intraEnergy;
  m.externalFieldMoleculeEnergy = a.externalFieldMoleculeEnergy * b.externalFieldMoleculeEnergy;
  m.frameworkMoleculeEnergy = a.frameworkMoleculeEnergy * b.frameworkMoleculeEnergy;
  m.interEnergy = a.interEnergy * b.interEnergy;
  for (size_t i = 0; i < a.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = a.intraComponentEnergies[i] * b.intraComponentEnergies[i];
  }
  for (size_t i = 0; i < a.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = a.externalFieldComponentEnergies[i] * b.externalFieldComponentEnergies[i];
  }
  for (size_t i = 0; i < a.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = a.frameworkComponentEnergies[i] * b.frameworkComponentEnergies[i];
  }
  for (size_t i = 0; i < a.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = a.interComponentEnergies[i] * b.interComponentEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator*(const double& a, const EnergyStatus& b)
{
  EnergyStatus m(b.numberOfExternalFields, b.numberOfFrameworks, b.numberOfComponents);
  m.totalEnergy = a * b.totalEnergy;
  m.dUdlambda = a * b.dUdlambda;
  m.intraEnergy = a * b.intraEnergy;
  m.externalFieldMoleculeEnergy = a * b.externalFieldMoleculeEnergy;
  m.frameworkMoleculeEnergy = a * b.frameworkMoleculeEnergy;
  m.interEnergy = a * b.interEnergy;
  for (size_t i = 0; i < b.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = a * b.intraComponentEnergies[i];
  }
  for (size_t i = 0; i < b.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = a * b.externalFieldComponentEnergies[i];
  }
  for (size_t i = 0; i < b.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = a * b.frameworkComponentEnergies[i];
  }
  for (size_t i = 0; i < b.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = a * b.interComponentEnergies[i];
  }

  return m;
}

export inline EnergyStatus operator/(const EnergyStatus& a, const double& b)
{
  EnergyStatus m(a.numberOfExternalFields, a.numberOfFrameworks, a.numberOfComponents);
  m.totalEnergy = a.totalEnergy / b;
  m.dUdlambda = a.dUdlambda / b;
  m.intraEnergy = a.intraEnergy / b;
  m.externalFieldMoleculeEnergy = a.externalFieldMoleculeEnergy / b;
  m.frameworkMoleculeEnergy = a.frameworkMoleculeEnergy / b;
  m.interEnergy = a.interEnergy / b;
  for (size_t i = 0; i < a.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = a.intraComponentEnergies[i] / b;
  }
  for (size_t i = 0; i < a.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = a.externalFieldComponentEnergies[i] / b;
  }
  for (size_t i = 0; i < a.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = a.frameworkComponentEnergies[i] / b;
  }
  for (size_t i = 0; i < a.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = a.interComponentEnergies[i] / b;
  }

  return m;
}

export inline EnergyStatus sqrt(const EnergyStatus& a)
{
  EnergyStatus m(a.numberOfExternalFields, a.numberOfFrameworks, a.numberOfComponents);
  m.totalEnergy = sqrt(a.totalEnergy);
  m.dUdlambda = std::sqrt(a.dUdlambda);
  m.intraEnergy = sqrt(a.intraEnergy);
  m.externalFieldMoleculeEnergy = sqrt(a.externalFieldMoleculeEnergy);
  m.frameworkMoleculeEnergy = sqrt(a.frameworkMoleculeEnergy);
  m.interEnergy = sqrt(a.interEnergy);
  for (size_t i = 0; i < a.intraComponentEnergies.size(); ++i)
  {
    m.intraComponentEnergies[i] = sqrt(a.intraComponentEnergies[i]);
  }
  for (size_t i = 0; i < a.externalFieldComponentEnergies.size(); ++i)
  {
    m.externalFieldComponentEnergies[i] = sqrt(a.externalFieldComponentEnergies[i]);
  }
  for (size_t i = 0; i < a.frameworkComponentEnergies.size(); ++i)
  {
    m.frameworkComponentEnergies[i] = sqrt(a.frameworkComponentEnergies[i]);
  }
  for (size_t i = 0; i < a.interComponentEnergies.size(); ++i)
  {
    m.interComponentEnergies[i] = sqrt(a.interComponentEnergies[i]);
  }

  return m;
}
