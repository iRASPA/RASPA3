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
    // Deliberately V / N, NOT the textbook V / (N + 1) of Sindzingre et al. (1987) / Smit & Frenkel (1989):
    //   1. V / (N + 1) only cancels the volume-fluctuation bias in a pure single-component closed NpT fluid where
    //      <V> = (N + 1) / (beta P). Whenever V is fixed or set by a host framework or solvent (NVT, flexible-host
    //      NpT, or a dilute solute in NpT), -ln(<V / (N + 1) W>) / beta is the forward step F(N + 1) - F(N) at
    //      density (N + 1) / V rather than the box density N / V: at N = 1 it reports twice the fugacity (+ln 2 in
    //      beta mu) and halves the Henry coefficient K_H = rho / f, whereas V / N gives F(1) - F(0) and the exact
    //      Henry limit.
    //   2. In an open or Gibbs box that visits N = 0, summing Q(N + 1, V) = (V / (N + 1)) Q(N, V) <W> from N = 0
    //      omits the M = 0 state and gives <V / (N + 1) W> = (1 - P(N = 0)) / (beta f), which plateaus at a finite
    //      fugacity 1 / (beta <V W>_0) instead of vanishing as f -> 0 (for boxes that empty frequently, see
    //      PropertyWidom, whose <N / V> / <W> = beta f has no N = 0 boundary term).
    // An empty box (N = 0) has no volume per molecule: it gives no chemical-potential sample (counted as one
    // molecule, it gave V W, a finite fugacity for a box without molecules).
    rosenbluthWeight.addSample(blockIndex, RosenbluthValue, weight);
    if (N == 0uz) return;
    const double volumePerMolecule = V / static_cast<double>(N);
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
    // a term without samples (zero) has no estimate: left at zero rather than the infinity of log(0)
    return WidomData(terms.total > 0.0 ? -(1.0 / beta) * std::log(terms.total) : 0.0,
                     terms.excess > 0.0 ? -(1.0 / beta) * std::log(terms.excess) : 0.0,
                     terms.idealGas > 0.0 ? -(1.0 / beta) * std::log(terms.idealGas) : 0.0);
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

  static double fugacityTransform(const WidomData& terms, double beta)
  {
    return terms.total > 0.0 ? 1.0 / (beta * terms.total) : 0.0;
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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyGibbsWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyGibbsWidom &w);
};
