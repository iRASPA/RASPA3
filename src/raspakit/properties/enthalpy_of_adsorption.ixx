module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#if !defined(_WIN32)
#include <assert.h>
#endif
#endif

export module enthalpy_of_adsorption;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <numeric>;
import <fstream>;
import <utility>;
import <string>;
import <cmath>;
import <iostream>;
#if defined(_WIN32)
import <cassert>;
#endif
#endif

import matrix;
import archive;
import energy_status;
import averages;
import units;

export struct EnthalpyOfAdsorption
{
  EnthalpyOfAdsorption(size_t size) : size(size), values(size) {}

  EnthalpyOfAdsorption(std::vector<double> values) : size(values.size()), values(values) {}

  bool operator==(EnthalpyOfAdsorption const&) const = default;

  size_t size;
  std::vector<double> values;

  inline EnthalpyOfAdsorption& operator+=(const EnthalpyOfAdsorption& b)
  {
    for (size_t i = 0; i < size; ++i)
    {
      values[i] += b.values[i];
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnthalpyOfAdsorption& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnthalpyOfAdsorption& p);
};

export inline EnthalpyOfAdsorption operator+(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
  EnthalpyOfAdsorption m(a.size);

  for (size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] + b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorption operator-(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
  EnthalpyOfAdsorption m(a.size);

  for (size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] - b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorption operator*(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
  EnthalpyOfAdsorption m(a.size);

  for (size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] * b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorption operator*(const double& a, const EnthalpyOfAdsorption& b)
{
  EnthalpyOfAdsorption m(b.size);

  for (size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a * b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorption sqrt(const EnthalpyOfAdsorption& a)
{
  EnthalpyOfAdsorption m(a.size);

  for (size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = std::sqrt(a.values[i]);
  }

  return m;
}

export struct EnthalpyOfAdsorptionTerms
{
  EnthalpyOfAdsorptionTerms(size_t size)
      : size(size),
        swappableComponents(size),
        totalEnergyTimesNumberOfMolecules(size),
        numberOfMoleculesSquared(size, std::vector<double>(size)),
        numberOfMolecules(size),
        temperature(0.0),
        totalEnergy(0.0)
  {
  }

  EnthalpyOfAdsorptionTerms(const EnthalpyOfAdsorptionTerms& a) noexcept = default;
  EnthalpyOfAdsorptionTerms& operator=(const EnthalpyOfAdsorptionTerms& a) noexcept = default;

  EnthalpyOfAdsorptionTerms(const std::vector<size_t> swappableComponents,
                            const std::vector<size_t> numberOfIntegerMolecules, double totalEnergy, double temperature)
      : size(swappableComponents.size()),
        swappableComponents(swappableComponents),
        totalEnergyTimesNumberOfMolecules(swappableComponents.size()),
        numberOfMoleculesSquared(swappableComponents.size(), std::vector<double>(swappableComponents.size())),
        numberOfMolecules(swappableComponents.size()),
        temperature(temperature * Units::KelvinToEnergy),
        totalEnergy(totalEnergy)
  {
    for (size_t i = 0; i < swappableComponents.size(); ++i)
    {
      size_t index_i = swappableComponents[i];
      numberOfMolecules[i] = static_cast<double>(numberOfIntegerMolecules[index_i]);

      totalEnergyTimesNumberOfMolecules[i] = totalEnergy * static_cast<double>(numberOfIntegerMolecules[index_i]);
      for (size_t j = 0; j < swappableComponents.size(); ++j)
      {
        size_t index_j = swappableComponents[j];
        numberOfMoleculesSquared[i][j] = static_cast<double>(numberOfIntegerMolecules[index_i]) *
                                         static_cast<double>(numberOfIntegerMolecules[index_j]);
      }
    }
  }

  EnthalpyOfAdsorptionTerms() = default;

  bool operator==(EnthalpyOfAdsorptionTerms const&) const = default;

  size_t size;
  std::vector<size_t> swappableComponents;
  std::vector<double> totalEnergyTimesNumberOfMolecules;
  std::vector<std::vector<double>> numberOfMoleculesSquared;
  std::vector<double> numberOfMolecules;
  double temperature;
  double totalEnergy;

  inline EnthalpyOfAdsorptionTerms& operator+=(const EnthalpyOfAdsorptionTerms& b)
  {
    totalEnergy += b.totalEnergy;
    temperature += b.temperature;
    for (size_t i = 0; i < size; ++i)
    {
      numberOfMolecules[i] += b.numberOfMolecules[i];
      totalEnergyTimesNumberOfMolecules[i] += b.totalEnergyTimesNumberOfMolecules[i];
      for (size_t j = 0; j < size; ++j)
      {
        numberOfMoleculesSquared[i][j] += b.numberOfMoleculesSquared[i][j];
      }
    }
    return *this;
  }

  inline EnthalpyOfAdsorption compositeProperty [[nodiscard]] () const
  {
    EnthalpyOfAdsorption v(size);
    if (size > 0)
    {
      // Symmetric matrix
      Matrix m(size, size, 0.0);
      for (size_t i = 0; i < size; ++i)
      {
        for (size_t j = 0; j < size; ++j)
        {
          m(i, j) = numberOfMoleculesSquared[i][j] - numberOfMolecules[i] * numberOfMolecules[j];
        }
      }

      m.inverse();

      for (size_t i = 0; i < size; ++i)
      {
        v.values[i] = 0.0;
        for (size_t j = 0; j < size; ++j)
        {
          v.values[i] += (totalEnergyTimesNumberOfMolecules[j] - totalEnergy * numberOfMolecules[j]) * m(j, i);
        }
        v.values[i] -= temperature;
      }
    }
    return v;
  }
  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnthalpyOfAdsorptionTerms& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnthalpyOfAdsorptionTerms& p);
};

export inline EnthalpyOfAdsorptionTerms operator+(const EnthalpyOfAdsorptionTerms& a,
                                                  const EnthalpyOfAdsorptionTerms& b)
{
  EnthalpyOfAdsorptionTerms m(a.size);

  m.totalEnergy = a.totalEnergy + b.totalEnergy;
  m.temperature = a.temperature + b.temperature;
  for (size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] + b.numberOfMolecules[i];
    m.totalEnergyTimesNumberOfMolecules[i] =
        a.totalEnergyTimesNumberOfMolecules[i] + b.totalEnergyTimesNumberOfMolecules[i];
    for (size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] + b.numberOfMoleculesSquared[i][j];
    }
  }

  return m;
}

export inline EnthalpyOfAdsorptionTerms operator*(const double& a, const EnthalpyOfAdsorptionTerms& b)
{
  EnthalpyOfAdsorptionTerms m(b.size);

  m.totalEnergy = a * b.totalEnergy;
  m.temperature = a * b.temperature;
  for (size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a * b.numberOfMolecules[i];
    m.totalEnergyTimesNumberOfMolecules[i] = a * b.totalEnergyTimesNumberOfMolecules[i];
    for (size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a * b.numberOfMoleculesSquared[i][j];
    }
  }

  return m;
}

export inline EnthalpyOfAdsorptionTerms operator/(const EnthalpyOfAdsorptionTerms& a, const double& b)
{
  EnthalpyOfAdsorptionTerms m(a.size);

  m.totalEnergy = a.totalEnergy / b;
  m.temperature = a.temperature / b;
  for (size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] / b;
    m.totalEnergyTimesNumberOfMolecules[i] = a.totalEnergyTimesNumberOfMolecules[i] / b;
    for (size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] / b;
    }
  }

  return m;
}
