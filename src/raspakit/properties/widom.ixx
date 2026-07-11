module;

export module property_widom;

import std;

import archive;
import int3;
import averages;
export import property_block_average;

export struct WidomData
{
  WidomData():
    total(0.0),
    excess(0.0),
    idealGas(0.0)
  {
  };

  WidomData(double total, double excess, double idealGas):
    total(total),
    excess(excess),
    idealGas(idealGas)
  {
  };

  inline WidomData& operator+=(const WidomData& b)
  {
    total += b.total;
    excess += b.excess;
    idealGas += b.idealGas;

    return *this;
  }

  std::uint64_t versionNumber{1};

  double total{};
  double excess{};
  double idealGas{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const WidomData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, WidomData& l);
};

export inline WidomData operator+(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total + b.total;
  m.excess = a.excess + b.excess;
  m.idealGas = a.idealGas + b.idealGas;

  return m;
}

export inline WidomData operator-(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total - b.total;
  m.excess = a.excess - b.excess;
  m.idealGas = a.idealGas - b.idealGas;

  return m;
}

export inline WidomData operator*(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total * b.total;
  m.excess = a.excess * b.excess;
  m.idealGas = a.idealGas * b.idealGas;

  return m;
}

export inline WidomData operator*(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a * b.total;
  m.excess = a * b.excess;
  m.idealGas = a * b.idealGas;

  return m;
}

export inline WidomData operator/(const WidomData& a, const double& b)
{
  WidomData m{}; 

  double inv_b = 1.0 / b;
  m.total = inv_b * a.total;
  m.excess = inv_b * a.excess;
  m.idealGas = inv_b * a.idealGas;

  return m;
}

export inline WidomData operator/(const WidomData& a, const std::array<double,3>& b)
{
  WidomData m{}; 

  m.total = a.total / b[0];
  m.excess = a.excess / b[1];
  m.idealGas = a.idealGas / b[2];

  return m;
}

export inline WidomData operator/(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a / b.total;
  m.excess = a / b.excess;
  m.idealGas = a / b.idealGas;

  return m;
}

export inline WidomData sqrt(const WidomData& a)
{
  WidomData m{};

  m.total = std::sqrt(a.total);
  m.excess = std::sqrt(a.excess);
  m.idealGas = std::sqrt(a.idealGas);

  return m;
}

export inline WidomData log(const WidomData& a)
{
  WidomData m{};

  m.total = std::log(a.total);
  m.excess = std::log(a.excess);
  m.idealGas = std::log(a.idealGas);

  return m;
}

/**
 * \brief Widom-insertion statistics: Rosenbluth weight, chemical potential and fugacity.
 *
 * Two block-averaged channels are sampled: the bare Rosenbluth weight and the raw chemical
 * potential terms (excess Rosenbluth weight and ideal-gas density). The chemical potential and
 * fugacity are non-linear functions of those averages, propagated through
 * BlockAverage::statistics().
 */
export struct PropertyWidom
{
  PropertyWidom() = default;

  PropertyWidom(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), rosenbluthWeight(numberOfBlocks), chemicalPotentialTerms(numberOfBlocks)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  BlockAverage<double> rosenbluthWeight;
  BlockAverage<WidomData> chemicalPotentialTerms;

  std::string writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                      std::optional<double> frameworkMass,
                                                      std::optional<int3> number_of_unit_cells) const;
  std::string writeAveragesChemicalPotentialStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                                       std::optional<double> imposedFugacity) const;

  inline void addWidomSample(std::size_t blockIndex, double RosenbluthValue, std::size_t N, double V, double weight)
  {
    rosenbluthWeight.addSample(blockIndex, RosenbluthValue, weight);
    chemicalPotentialTerms.addSample(blockIndex, WidomData(0.0, RosenbluthValue, static_cast<double>(N) / V), weight);
  }

  //====================================================================================================================

  double averagedRosenbluthWeight() const { return rosenbluthWeight.averaged(); }
  double averagedRosenbluthWeight(std::size_t blockIndex) const { return rosenbluthWeight.averaged(blockIndex); }

  std::pair<double, double> result() const { return rosenbluthWeight.average(); }

  //====================================================================================================================

  /// Chemical potential as a non-linear function of the averaged raw terms.
  static WidomData chemicalPotentialTransform(const WidomData &terms, double beta)
  {
    return WidomData(-(1.0 / beta) * std::log(terms.excess) + (1.0 / beta) * std::log(terms.idealGas),
                     -(1.0 / beta) * std::log(terms.excess),
                      (1.0 / beta) * std::log(terms.idealGas));
  }

  WidomData averagedChemicalPotential(double beta) const
  {
    return chemicalPotentialTransform(chemicalPotentialTerms.averaged(), beta);
  }

  WidomData averagedChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return chemicalPotentialTransform(chemicalPotentialTerms.averaged(blockIndex), beta);
  }

  std::pair<WidomData, WidomData> chemicalPotentialResult(double beta) const
  {
    return chemicalPotentialTerms.statistics(
        [beta](const WidomData &terms) { return chemicalPotentialTransform(terms, beta); });
  }

  //====================================================================================================================

  static double fugacityTransform(const WidomData &terms, double beta)
  {
    return std::exp(beta * chemicalPotentialTransform(terms, beta).total) / beta;
  }

  double averagedFugacity(double beta) const { return fugacityTransform(chemicalPotentialTerms.averaged(), beta); }

  double averagedFugacity(std::size_t blockIndex, double beta) const
  {
    return fugacityTransform(chemicalPotentialTerms.averaged(blockIndex), beta);
  }

  std::pair<double, double> fugacityResult(double beta) const
  {
    return chemicalPotentialTerms.statistics([beta](const WidomData &terms) { return fugacityTransform(terms, beta); });
  }

  //====================================================================================================================

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w);
};
