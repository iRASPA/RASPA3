module;

export module partial_molar_properties_data;

import std;

import matrix;
import archive;
import units;

/**
 * \brief Result of a partial-molar-property calculation.
 *
 * Holds, per swappable component, the partial molar internal energy
 * \f$\bar{U}_i=(\partial U/\partial N_i)_{T,V,N_{j\neq i}}\f$ and the partial molar
 * volume \f$\bar{V}_i=(\partial V/\partial N_i)_{T,P,N_{j\neq i}}\f$, both obtained
 * from grand-canonical/osmotic energy(volume)-particle fluctuations.
 */
export struct PartialMolarPropertiesData
{
  PartialMolarPropertiesData(std::size_t size) : size(size), partialMolarEnergy(size), partialMolarVolume(size) {}

  bool operator==(PartialMolarPropertiesData const&) const = default;

  std::size_t size;
  std::vector<double> partialMolarEnergy;
  std::vector<double> partialMolarVolume;

  inline PartialMolarPropertiesData& operator+=(const PartialMolarPropertiesData& b)
  {
    for (std::size_t i = 0; i < size; ++i)
    {
      partialMolarEnergy[i] += b.partialMolarEnergy[i];
      partialMolarVolume[i] += b.partialMolarVolume[i];
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const PartialMolarPropertiesData& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, PartialMolarPropertiesData& p);
};

export inline PartialMolarPropertiesData operator-(const PartialMolarPropertiesData& a,
                                                   const PartialMolarPropertiesData& b)
{
  PartialMolarPropertiesData m(a.size);
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.partialMolarEnergy[i] = a.partialMolarEnergy[i] - b.partialMolarEnergy[i];
    m.partialMolarVolume[i] = a.partialMolarVolume[i] - b.partialMolarVolume[i];
  }
  return m;
}

export inline PartialMolarPropertiesData operator*(const PartialMolarPropertiesData& a,
                                                   const PartialMolarPropertiesData& b)
{
  PartialMolarPropertiesData m(a.size);
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.partialMolarEnergy[i] = a.partialMolarEnergy[i] * b.partialMolarEnergy[i];
    m.partialMolarVolume[i] = a.partialMolarVolume[i] * b.partialMolarVolume[i];
  }
  return m;
}

export inline PartialMolarPropertiesData operator*(const double& a, const PartialMolarPropertiesData& b)
{
  PartialMolarPropertiesData m(b.size);
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.partialMolarEnergy[i] = a * b.partialMolarEnergy[i];
    m.partialMolarVolume[i] = a * b.partialMolarVolume[i];
  }
  return m;
}

export inline PartialMolarPropertiesData sqrt(const PartialMolarPropertiesData& a)
{
  PartialMolarPropertiesData m(a.size);
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.partialMolarEnergy[i] = std::sqrt(a.partialMolarEnergy[i]);
    m.partialMolarVolume[i] = std::sqrt(a.partialMolarVolume[i]);
  }
  return m;
}

/**
 * \brief Raw fluctuation terms accumulated for partial-molar-property calculations.
 *
 * Each production sample stores the instantaneous \f$N_i\f$, \f$N_iN_j\f$, total internal
 * energy \f$U\f$, \f$N_iU\f$, volume \f$V\f$, and \f$N_iV\f$. After block averaging, these
 * become the ensemble averages needed by the matrix fluctuation formula in
 * compositeProperty().
 */
export struct PartialMolarPropertiesTerms
{
  PartialMolarPropertiesTerms(std::size_t size)
      : size(size),
        swappableComponents(size),
        totalEnergyTimesNumberOfMolecules(size),
        volumeTimesNumberOfMolecules(size),
        numberOfMoleculesSquared(size, std::vector<double>(size)),
        numberOfMolecules(size),
        totalEnergy(0.0),
        volume(0.0)
  {
  }

  PartialMolarPropertiesTerms(const PartialMolarPropertiesTerms& a) noexcept = default;
  PartialMolarPropertiesTerms& operator=(const PartialMolarPropertiesTerms& a) noexcept = default;

  PartialMolarPropertiesTerms(const std::vector<std::size_t> swappableComponents,
                              const std::vector<std::size_t> numberOfIntegerMolecules, double totalEnergy,
                              double volume)
      : size(swappableComponents.size()),
        swappableComponents(swappableComponents),
        totalEnergyTimesNumberOfMolecules(swappableComponents.size()),
        volumeTimesNumberOfMolecules(swappableComponents.size()),
        numberOfMoleculesSquared(swappableComponents.size(), std::vector<double>(swappableComponents.size())),
        numberOfMolecules(swappableComponents.size()),
        totalEnergy(totalEnergy),
        volume(volume)
  {
    for (std::size_t i = 0; i < swappableComponents.size(); ++i)
    {
      std::size_t index_i = swappableComponents[i];
      numberOfMolecules[i] = static_cast<double>(numberOfIntegerMolecules[index_i]);

      totalEnergyTimesNumberOfMolecules[i] = totalEnergy * static_cast<double>(numberOfIntegerMolecules[index_i]);
      volumeTimesNumberOfMolecules[i] = volume * static_cast<double>(numberOfIntegerMolecules[index_i]);
      for (std::size_t j = 0; j < swappableComponents.size(); ++j)
      {
        std::size_t index_j = swappableComponents[j];
        numberOfMoleculesSquared[i][j] = static_cast<double>(numberOfIntegerMolecules[index_i]) *
                                         static_cast<double>(numberOfIntegerMolecules[index_j]);
      }
    }
  }

  PartialMolarPropertiesTerms() = default;

  bool operator==(PartialMolarPropertiesTerms const&) const = default;

  std::size_t size;
  std::vector<std::size_t> swappableComponents;
  std::vector<double> totalEnergyTimesNumberOfMolecules;
  std::vector<double> volumeTimesNumberOfMolecules;
  std::vector<std::vector<double>> numberOfMoleculesSquared;
  std::vector<double> numberOfMolecules;
  double totalEnergy;
  double volume;

  inline PartialMolarPropertiesTerms& operator+=(const PartialMolarPropertiesTerms& b)
  {
    totalEnergy += b.totalEnergy;
    volume += b.volume;
    for (std::size_t i = 0; i < size; ++i)
    {
      numberOfMolecules[i] += b.numberOfMolecules[i];
      totalEnergyTimesNumberOfMolecules[i] += b.totalEnergyTimesNumberOfMolecules[i];
      volumeTimesNumberOfMolecules[i] += b.volumeTimesNumberOfMolecules[i];
      for (std::size_t j = 0; j < size; ++j)
      {
        numberOfMoleculesSquared[i][j] += b.numberOfMoleculesSquared[i][j];
      }
    }
    return *this;
  }

  /**
   * \brief Compute partial molar energy and volume from the (block-averaged) fluctuation terms.
   *
   * Builds the symmetric particle-number covariance matrix
   * \f$\mathrm{Cov}(N)_{ij}=\langle N_iN_j\rangle-\langle N_i\rangle\langle N_j\rangle\f$, inverts it,
   * and contracts it with the energy-number and volume-number cross-correlation vectors:
   * \f[\bar{U}_i=\sum_j(\langle N_jU\rangle-\langle N_j\rangle\langle U\rangle)[\mathrm{Cov}(N)^{-1}]_{ji},\qquad
   *    \bar{V}_i=\sum_j(\langle N_jV\rangle-\langle N_j\rangle\langle V\rangle)[\mathrm{Cov}(N)^{-1}]_{ji}.\f]
   */
  inline PartialMolarPropertiesData compositeProperty [[nodiscard]] () const
  {
    PartialMolarPropertiesData v(size);
    if (size > 0)
    {
      // Symmetric particle-number covariance matrix
      Matrix m(size, size, 0.0);
      for (std::size_t i = 0; i < size; ++i)
      {
        for (std::size_t j = 0; j < size; ++j)
        {
          m(i, j) = numberOfMoleculesSquared[i][j] - numberOfMolecules[i] * numberOfMolecules[j];
        }
      }

      m.inverse();

      for (std::size_t i = 0; i < size; ++i)
      {
        v.partialMolarEnergy[i] = 0.0;
        v.partialMolarVolume[i] = 0.0;
        for (std::size_t j = 0; j < size; ++j)
        {
          v.partialMolarEnergy[i] += (totalEnergyTimesNumberOfMolecules[j] - totalEnergy * numberOfMolecules[j]) * m(j, i);
          v.partialMolarVolume[i] += (volumeTimesNumberOfMolecules[j] - volume * numberOfMolecules[j]) * m(j, i);
        }
      }
    }
    return v;
  }
  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const PartialMolarPropertiesTerms& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, PartialMolarPropertiesTerms& p);
};

export inline PartialMolarPropertiesTerms operator*(const double& a, const PartialMolarPropertiesTerms& b)
{
  PartialMolarPropertiesTerms m(b.size);

  m.totalEnergy = a * b.totalEnergy;
  m.volume = a * b.volume;
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a * b.numberOfMolecules[i];
    m.totalEnergyTimesNumberOfMolecules[i] = a * b.totalEnergyTimesNumberOfMolecules[i];
    m.volumeTimesNumberOfMolecules[i] = a * b.volumeTimesNumberOfMolecules[i];
    for (std::size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a * b.numberOfMoleculesSquared[i][j];
    }
  }

  return m;
}

export inline PartialMolarPropertiesTerms operator/(const PartialMolarPropertiesTerms& a, const double& b)
{
  PartialMolarPropertiesTerms m(a.size);

  m.totalEnergy = a.totalEnergy / b;
  m.volume = a.volume / b;
  for (std::size_t i = 0; i < m.size; ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] / b;
    m.totalEnergyTimesNumberOfMolecules[i] = a.totalEnergyTimesNumberOfMolecules[i] / b;
    m.volumeTimesNumberOfMolecules[i] = a.volumeTimesNumberOfMolecules[i] / b;
    for (std::size_t j = 0; j < m.size; ++j)
    {
      m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] / b;
    }
  }

  return m;
}
