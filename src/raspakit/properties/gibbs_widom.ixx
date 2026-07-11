module;

export module property_gibbs_widom;

import std;

import archive;
import int3;
import averages;
export import property_widom;
export import property_block_average;

/**
 * \brief Gibbs-ensemble Widom-insertion statistics.
 *
 * Mirrors PropertyWidom, but the raw chemical potential terms include the fluctuating volume per
 * molecule of the Gibbs ensemble. The chemical potential and fugacity are non-linear functions of
 * the averaged terms, propagated through BlockAverage::statistics().
 */
export struct PropertyGibbsWidom
{
  PropertyGibbsWidom() = default;

  PropertyGibbsWidom(std::size_t numberOfBlocks)
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
    double volumePerMolecule = V / static_cast<double>(std::max(N, 1uz));
    rosenbluthWeight.addSample(blockIndex, RosenbluthValue, weight);
    chemicalPotentialTerms.addSample(
        blockIndex, WidomData(volumePerMolecule * RosenbluthValue, RosenbluthValue, volumePerMolecule), weight);
  }

  //====================================================================================================================

  double averagedRosenbluthWeight() const { return rosenbluthWeight.averaged(); }
  double averagedRosenbluthWeight(std::size_t blockIndex) const { return rosenbluthWeight.averaged(blockIndex); }

  std::pair<double, double> result() const { return rosenbluthWeight.average(); }

  //====================================================================================================================

  static WidomData chemicalPotentialTransform(const WidomData &terms, double beta)
  {
    return -(1.0 / beta) * log(terms);
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

  static double fugacityTransform(const WidomData &terms, double beta) { return 1.0 / (beta * terms.total); }

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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyGibbsWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyGibbsWidom &w);
};
