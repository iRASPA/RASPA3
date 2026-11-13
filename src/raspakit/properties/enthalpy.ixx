module;

export module property_enthalpy;

import std;

import archive;
import averages;
import component;
import json;
import units;
export import property_block_average;
import matrix;

export struct EnthalpyOfAdsorptionData
{
  EnthalpyOfAdsorptionData(std::size_t size) : size(size), values(size) {}

  EnthalpyOfAdsorptionData(std::vector<double> values) : size(values.size()), values(values) {}

  bool operator==(EnthalpyOfAdsorptionData const&) const = default;

  double& operator[](std::size_t index) {return values[index];}

  std::size_t size;
  std::vector<double> values;
  /// False for a block whose particle-number covariance is singular: it has no enthalpy (the code base is
  /// built with fast math, so a not-a-number cannot be tested for; this flag is what marks the block).
  bool defined{true};

  inline EnthalpyOfAdsorptionData& operator+=(const EnthalpyOfAdsorptionData& b)
  {
    for (std::size_t i = 0; i < size; ++i)
    {
      values[i] += b.values[i];
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EnthalpyOfAdsorptionData& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EnthalpyOfAdsorptionData& p);
};

/// Whether a block's enthalpies are defined (a block with a singular particle-number covariance has none).
export inline bool isFinite(const EnthalpyOfAdsorptionData& data) { return data.defined; }

export inline EnthalpyOfAdsorptionData operator+(const EnthalpyOfAdsorptionData& a, const EnthalpyOfAdsorptionData& b)
{
  EnthalpyOfAdsorptionData m(a.size);

  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] + b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorptionData operator-(const EnthalpyOfAdsorptionData& a, const EnthalpyOfAdsorptionData& b)
{
  EnthalpyOfAdsorptionData m(a.size);

  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] - b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorptionData operator*(const EnthalpyOfAdsorptionData& a, const EnthalpyOfAdsorptionData& b)
{
  EnthalpyOfAdsorptionData m(a.size);

  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a.values[i] * b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorptionData operator*(const double& a, const EnthalpyOfAdsorptionData& b)
{
  EnthalpyOfAdsorptionData m(b.size);

  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = a * b.values[i];
  }

  return m;
}

export inline EnthalpyOfAdsorptionData sqrt(const EnthalpyOfAdsorptionData& a)
{
  EnthalpyOfAdsorptionData m(a.size);

  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.values[i] = std::sqrt(a.values[i]);
  }

  return m;
}

export struct EnthalpyOfAdsorptionTerms
{
  EnthalpyOfAdsorptionTerms(std::size_t size)
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

  EnthalpyOfAdsorptionTerms(const std::vector<std::size_t> swappableComponents,
                            const std::vector<std::size_t> numberOfIntegerMolecules, double totalEnergy,
                            double temperature)
      : size(swappableComponents.size()),
        swappableComponents(swappableComponents),
        totalEnergyTimesNumberOfMolecules(swappableComponents.size()),
        numberOfMoleculesSquared(swappableComponents.size(), std::vector<double>(swappableComponents.size())),
        numberOfMolecules(swappableComponents.size()),
        temperature(temperature * Units::KelvinToEnergy),
        totalEnergy(totalEnergy)
  {
    for (std::size_t i = 0; i < swappableComponents.size(); ++i)
    {
      std::size_t index_i = swappableComponents[i];
      numberOfMolecules[i] = static_cast<double>(numberOfIntegerMolecules[index_i]);

      totalEnergyTimesNumberOfMolecules[i] = totalEnergy * static_cast<double>(numberOfIntegerMolecules[index_i]);
      for (std::size_t j = 0; j < swappableComponents.size(); ++j)
      {
        std::size_t index_j = swappableComponents[j];
        numberOfMoleculesSquared[i][j] = static_cast<double>(numberOfIntegerMolecules[index_i]) *
                                         static_cast<double>(numberOfIntegerMolecules[index_j]);
      }
    }
  }

  EnthalpyOfAdsorptionTerms() = default;

  bool operator==(EnthalpyOfAdsorptionTerms const&) const = default;

  std::size_t size;
  std::vector<std::size_t> swappableComponents;
  std::vector<double> totalEnergyTimesNumberOfMolecules;
  std::vector<std::vector<double>> numberOfMoleculesSquared;
  std::vector<double> numberOfMolecules;
  double temperature;
  double totalEnergy;

  inline EnthalpyOfAdsorptionTerms& operator+=(const EnthalpyOfAdsorptionTerms& b)
  {
    totalEnergy += b.totalEnergy;
    temperature += b.temperature;
    for (std::size_t i = 0; i < size; ++i)
    {
      numberOfMolecules[i] += b.numberOfMolecules[i];
      totalEnergyTimesNumberOfMolecules[i] += b.totalEnergyTimesNumberOfMolecules[i];
      for (std::size_t j = 0; j < size; ++j)
      {
        numberOfMoleculesSquared[i][j] += b.numberOfMoleculesSquared[i][j];
      }
    }
    return *this;
  }

  /// Correctly-shaped zero of the derived type; the composite of zero terms is not a valid zero.
  inline EnthalpyOfAdsorptionData zeroCompositeProperty [[nodiscard]] () const
  {
    return EnthalpyOfAdsorptionData(size);
  }

  inline EnthalpyOfAdsorptionData compositeProperty [[nodiscard]] () const
  {
    EnthalpyOfAdsorptionData v(size);
    if (size > 0)
    {
      // Symmetric matrix
      Matrix m(size, size, 0.0);
      for (std::size_t i = 0; i < size; ++i)
      {
        for (std::size_t j = 0; j < size; ++j)
        {
          m(i, j) = numberOfMoleculesSquared[i][j] - numberOfMolecules[i] * numberOfMolecules[j];
        }
      }

      // a singular covariance (a component whose count did not fluctuate in this block, or counts that
      // are collinear, as with pair or group swaps) has no enthalpy: the block is marked undefined and
      // skipped by the block statistics (see isFinite); using the failed inverse gave -RT for one
      // component and the LU factors for several
      if (!m.inverse())
      {
        v.defined = false;
        return v;
      }

      for (std::size_t i = 0; i < size; ++i)
      {
        v.values[i] = 0.0;
        for (std::size_t j = 0; j < size; ++j)
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
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] + b.numberOfMolecules[i];
    m.totalEnergyTimesNumberOfMolecules[i] =
        a.totalEnergyTimesNumberOfMolecules[i] + b.totalEnergyTimesNumberOfMolecules[i];
    for (std::size_t j = 0; j < m.size; ++j)
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
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a * b.numberOfMolecules[i];
    m.totalEnergyTimesNumberOfMolecules[i] = a * b.totalEnergyTimesNumberOfMolecules[i];
    for (std::size_t j = 0; j < m.size; ++j)
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
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] / b;
    m.totalEnergyTimesNumberOfMolecules[i] = a.totalEnergyTimesNumberOfMolecules[i] / b;
    for (std::size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] / b;
    }
  }

  return m;
}

/**
 * \brief Block-averaged enthalpy of adsorption.
 *
 * A composite (fluctuation-derived) property: EnthalpyOfAdsorptionTerms are accumulated per block
 * and reduced to the per-component enthalpy through the matrix fluctuation formula. All accumulation,
 * averaging and error estimation is inherited from BlockAverage; only the enthalpy-specific
 * output formatting lives here.
 */
export struct PropertyEnthalpy : BlockAverage<EnthalpyOfAdsorptionTerms>
{
  PropertyEnthalpy() = default;
  PropertyEnthalpy(std::size_t numberOfBlocks, std::size_t numberOfComponents)
      : BlockAverage<EnthalpyOfAdsorptionTerms>(numberOfBlocks, EnthalpyOfAdsorptionTerms(numberOfComponents))
  {
  }

  bool operator==(PropertyEnthalpy const &) const = default;

  void resize(std::size_t newNumberOfComponents)
  {
    BlockAverage<EnthalpyOfAdsorptionTerms>::resize(EnthalpyOfAdsorptionTerms(newNumberOfComponents));
  }

  std::vector<std::pair<double, double>> result() const
  {
    std::vector<std::pair<double, double>> data(zeroSample.size);
    std::pair<EnthalpyOfAdsorptionData, EnthalpyOfAdsorptionData> enthalpy = average();
    for (std::size_t k = 0; k < zeroSample.size; k++)
    {
      data[k] = {enthalpy.first.values[k], enthalpy.second[k]};
    }
    return data;
  }

  /// The adsorbed-phase mole fractions <N_k> / sum_j <N_j> of the swappable components in a block (all zero
  /// when nothing was adsorbed): the weights of the total enthalpy of adsorption, which is per adsorbed
  /// molecule; the bulk gas fractions, used before, weigh the reservoir instead.
  std::vector<double> blockAdsorbedMoleFractions(std::size_t blockIndex) const
  {
    const EnthalpyOfAdsorptionTerms normalized =
        bookKeeping[blockIndex].first / std::max(1.0, bookKeeping[blockIndex].second);
    std::vector<double> fractions(normalized.size, 0.0);
    const double total = std::accumulate(normalized.numberOfMolecules.begin(), normalized.numberOfMolecules.end(), 0.0);
    if (total > 0.0)
    {
      for (std::size_t k = 0; k < normalized.size; ++k) fractions[k] = normalized.numberOfMolecules[k] / total;
    }
    return fractions;
  }

  /// The total enthalpy of adsorption per block, in internal energy units; nothing for a block whose
  /// enthalpies are undefined (see enthalpy.cpp).
  std::vector<std::optional<double>> blockTotalEnthalpies(const std::vector<std::size_t>& swappableComponents,
                                                          const std::vector<Component>& components) const;
  std::string writeAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                      std::vector<Component> &components) const;
  nlohmann::json jsonAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                        std::vector<Component> &components) const;
};
