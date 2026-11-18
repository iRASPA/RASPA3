module;

export module property_lambda_probability_histogram;

import std;

import randomnumbers;
import archive;
import averages;
import json;
import property_block_average;

/**
 * \brief CFCMC lambda statistics: probability histogram, dU/dlambda and fluid density.
 *
 * Three statistical channels are sampled per block:
 * - the lambda visit counts (per-bin histogram), from which the Landau free-energy profile and
 *   the excess chemical potential (-ln(P(1)/P(0))/beta) follow;
 * - the thermodynamic-integration terms dU/dlambda (per-bin weighted averages), integrated with
 *   Simpson's rule to an alternative excess chemical potential;
 * - the weighted fluid density, a scalar BlockAverage whose non-linear transform
 *   ln(rho)/beta yields the ideal-gas chemical potential.
 *
 * The Wang-Landau biasing state (bias factors, scaling factor, visit histogram) is kept alongside.
 */
export struct PropertyLambdaProbabilityHistogram
{
  enum class WangLandauPhase : std::size_t
  {
    Initialize = 0,
    Sample = 1,
    AdjustBiasingFactors = 2,
    Finalize = 3
  };

  PropertyLambdaProbabilityHistogram() {};

  PropertyLambdaProbabilityHistogram(std::size_t numberOfBlocks, std::size_t numberOfSamplePoints)
      : numberOfBlocks(numberOfBlocks),
        numberOfSamplePoints(numberOfSamplePoints),
        currentBin(0),
        // a single point has no spacing: its lambda is zero (0 times infinity would be a NaN)
        delta(numberOfSamplePoints > 1 ? 1.0 / static_cast<double>(numberOfSamplePoints - 1) : 0.0),
        histogram(numberOfSamplePoints),
        biasFactor(numberOfSamplePoints),
        bookKeepingLambda(std::vector<std::vector<double>>(numberOfBlocks, std::vector<double>(numberOfSamplePoints))),
        density(numberOfBlocks),
        computeDUdlambda(false),
        bookKeepingDUdlambda(std::vector<std::vector<std::pair<double, double>>>(
            numberOfBlocks, std::vector<std::pair<double, double>>(numberOfSamplePoints)))
  {
  }

  bool operator==(PropertyLambdaProbabilityHistogram const &) const = default;

  std::uint64_t versionNumber{1};
  // zero for a default-constructed histogram (Python constructs one): read uninitialized they sized the results
  std::size_t numberOfBlocks{0};

  std::size_t numberOfSamplePoints{0};
  std::size_t currentBin{0};
  double delta{0.0};

  double WangLandauScalingFactor{1.0};

  std::vector<double> histogram;
  std::vector<double> biasFactor;

  // lambda-histogram: per-block visit counts per bin
  std::vector<std::vector<double>> bookKeepingLambda;

  // weighted fluid density
  BlockAverage<double> density;

  // dU/dlambda-histogram
  bool computeDUdlambda{false};
  // 1-based thermodynamic-integration group id assigned by System::assignDUdlambdaGroups();
  // 0 when this lambda coordinate is not tracked. Atoms of the corresponding fractional molecule
  // carry this id in Atom::groupId, and RunningEnergy accumulates dU/dlambda per group.
  std::uint8_t dUdlambdaGroupId{0};
  // per-block, per-bin (sum of dU/dlambda, number of samples)
  std::vector<std::vector<std::pair<double, double>>> bookKeepingDUdlambda;

  // fractional molecule occupancy
  std::size_t occupancyCount{0};
  std::size_t occupancyTotal{0};

  void clear()
  {
    std::fill(histogram.begin(), histogram.end(), 0.0);

    density = BlockAverage<double>(numberOfBlocks);
    for (std::size_t i = 0; i < numberOfBlocks; ++i)
    {
      std::fill(bookKeepingLambda[i].begin(), bookKeepingLambda[i].end(), 0.0);
      std::fill(bookKeepingDUdlambda[i].begin(), bookKeepingDUdlambda[i].end(),
                std::make_pair<double, double>(0.0, 0.0));
    }

    occupancyCount = 0;
    occupancyTotal = 0;
  }

  inline double lambdaValue() const { return static_cast<double>(currentBin) * delta; }

  /// The grid point nearest to 'lambda': the inverse of 'lambdaValue' (exact for a lambda on the grid). The CFCMC
  /// moves step along the grid, so the lambda of a fractional molecule is always that of a bin.
  inline std::size_t binOfLambda(double lambda) const
  {
    if (numberOfSamplePoints < 2) return 0;
    const double points = static_cast<double>(numberOfSamplePoints - 1);
    return static_cast<std::size_t>(std::clamp(std::round(lambda * points), 0.0, points));
  }

  inline int selectNewBin(RandomNumber &random, double scale) const
  {
    // a single point: every step is a direct crossing of one end of the lambda path (a fraction of a bin would
    // truncate to no step at all)
    if (numberOfSamplePoints <= 1) return static_cast<int>(currentBin) + (random.uniform() < 0.5 ? -1 : 1);
    // A step of 1 .. K bins either way, K = max(1, floor(scale * points)), each with the same probability (a
    // symmetric proposal). The truncated real step it replaces was zero for every draw once scale * points fell
    // to one or below (two points from the start, ten points at the smallest scale): lambda then never moved, and
    // above that a share of the moves evaluated the full energy of a step to the bin they were in and were
    // counted as accepted, which the step-size optimization read as room to shrink further.
    const int maximumStep = std::max(1, static_cast<int>(scale * static_cast<double>(numberOfSamplePoints)));
    const int draw = std::min(2 * maximumStep - 1,
                              static_cast<int>(random.uniform() * static_cast<double>(2 * maximumStep)));
    const int step = draw < maximumStep ? draw - maximumStep : draw - maximumStep + 1;
    return static_cast<int>(currentBin) + step;
  }

  inline void setCurrentBin(std::size_t index) { currentBin = index; }

  inline void updateHistogram() { histogram[currentBin] += 1.0; }

  void sampleOccupancy(bool containsTheFractionalMolecule)
  {
    ++occupancyTotal;
    if (containsTheFractionalMolecule)
    {
      ++occupancyCount;
    }
  }

  double occupancy() const
  {
    return static_cast<double>(occupancyCount) / static_cast<double>(std::max(std::size_t{1}, occupancyTotal));
  }

  void normalize(double normalizationFactor)
  {
    for (double &bias : biasFactor)
    {
      bias -= normalizationFactor;
    }
  }

  /// Shifts the bias so that its smallest value is zero. The weight exp(-bias) of every bin is then at most
  /// one: Wang-Landau only ever lowers the bias, and exp(-bias) of the raw values overflows after a few
  /// thousand cycles, turning every weighted average into inf / inf.
  void normalizeToMinimum()
  {
    if (biasFactor.empty()) return;
    normalize(*std::min_element(biasFactor.cbegin(), biasFactor.cend()));
  }

  /// 'w' is the weight of the configuration, exp(-sum of the biases of every fractional molecule)
  /// ('System::weight'), and 'ownBiasInWeight' says whether the bias of this histogram is one of them.
  ///
  /// The histogram and dU/dlambda are kept under the bias of this histogram itself (the results take it out: the
  /// free-energy profile and the excess chemical potential add the bias back, and dU/dlambda at one lambda does not
  /// feel it), but the biases of the other fractional molecules in the box have to go: a sample counts with
  /// exp(-sum of their biases), which is 'w' over the factor of this histogram. Counted with one, the
  /// distribution along this lambda and the average dU/dlambda at each lambda were those of a system in which the
  /// other fractional molecules are biased, as soon as two fractional molecules shared a box. With one fractional
  /// molecule the factor is exactly one.
  void sampleHistogram(std::size_t blockIndex, double fluidDensity, double dUdlambda,
                       bool containsTheFractionalMolecule, double w, bool ownBiasInWeight = true)
  {
    if (containsTheFractionalMolecule)
    {
      double weightOfTheOthers = w;
      if (ownBiasInWeight)
      {
        const double own = weight();
        weightOfTheOthers = (own > 0.0 && std::isfinite(own) && std::isfinite(w)) ? w / own : 1.0;
      }
      bookKeepingLambda[blockIndex][currentBin] += weightOfTheOthers;

      bookKeepingDUdlambda[blockIndex][currentBin].first += weightOfTheOthers * dUdlambda;
      bookKeepingDUdlambda[blockIndex][currentBin].second += weightOfTheOthers;
    }

    density.addSample(blockIndex, fluidDensity, w);
  }

  inline double weight() const { return std::exp(-biasFactor[currentBin]); }

  void WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase,
                           bool containsTheFractionalMolecule = true, double value = 1.0);

  std::pair<std::vector<double>, std::vector<double>> result();

  std::string writeAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                      std::optional<double> imposedFugacity) const;

  std::string writeDUdLambdaStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                       std::optional<double> imposedFugacity) const;

  nlohmann::json jsonAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                        std::optional<double> imposedFugacity) const;

  nlohmann::json jsonDUdLambdaStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                         std::optional<double> imposedFugacity) const;

  //====================================================================================================================
  // density and ideal-gas chemical potential (scalar BlockAverage channel)

  double averagedDensity(std::size_t blockIndex) const { return density.averaged(blockIndex); }
  double averagedDensity() const { return density.averaged(); }
  std::pair<double, double> averageDensity() const { return density.average(); }

  static double idealGasChemicalPotentialTransform(double fluidDensity, double beta)
  {
    // without molecules (no density) there is no ideal-gas term: zero rather than the infinity of log(0)
    return fluidDensity > 0.0 ? std::log(fluidDensity) / beta : 0.0;
  }

  double averagedIdealGasChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return idealGasChemicalPotentialTransform(density.averaged(blockIndex), beta);
  }

  double averagedIdealGasChemicalPotential(double beta) const
  {
    return idealGasChemicalPotentialTransform(density.averaged(), beta);
  }

  std::pair<double, double> averageIdealGasChemicalPotential(double beta) const
  {
    return density.statistics([beta](double d) { return idealGasChemicalPotentialTransform(d, beta); });
  }

  //====================================================================================================================
  // lambda probability histogram (per-bin visit counts)

  /// Visit counts pooled over all blocks.
  std::vector<double> summedLambdaBins() const
  {
    std::vector<double> summedBlocks(numberOfSamplePoints);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingLambda[blockIndex].begin(),
                     summedBlocks.begin(), std::plus<double>());
    }
    return summedBlocks;
  }

  /// The visit counts normalised to a probability (unit sum); zero when nothing was sampled.
  static std::vector<double> normalizedHistogram(std::vector<double> counts)
  {
    const double total = std::accumulate(counts.begin(), counts.end(), 0.0);
    if (total > 0.0)
    {
      for (double& value : counts) value /= total;
    }
    return counts;
  }

  /// The probability histogram of one block: its own counts normalised, so blocks of different size compare bin by
  /// bin (unnormalised counts shifted every block's free-energy profile by the logarithm of its size).
  std::vector<double> averagedProbabilityHistogram(std::size_t blockIndex) const
  {
    return normalizedHistogram(bookKeepingLambda[blockIndex]);
  }

  /// The probability histogram pooled over all blocks.
  std::vector<double> averagedProbabilityHistogram() const { return normalizedHistogram(summedLambdaBins()); }

  /// The blocks that hold samples (and, for the chemical potentials, samples at both ends of the lambda path): a
  /// block that was never reached has no estimate, and taking it as zero pulled the averages and errors towards it.
  bool blockSampled(std::size_t blockIndex) const
  {
    return std::ranges::any_of(bookKeepingLambda[blockIndex], [](double count) { return count > 0.0; });
  }
  bool blockSampledAtBothEnds(std::size_t blockIndex) const
  {
    return bookKeepingLambda[blockIndex].front() > 0.0 && bookKeepingLambda[blockIndex].back() > 0.0;
  }

  /// Whether the run as a whole has samples at both ends of the lambda path.
  bool sampledAtBothEnds() const
  {
    const std::vector<double> bins = summedLambdaBins();
    return !bins.empty() && bins.front() > 0.0 && bins.back() > 0.0;
  }

  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram() const
  {
    std::vector<double> average = averagedProbabilityHistogram();

    std::vector<std::vector<double>> blockAverages;
    blockAverages.reserve(numberOfBlocks);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampled(blockIndex)) blockAverages.push_back(averagedProbabilityHistogram(blockIndex));

    std::vector<double> confidenceIntervalError = blockErrorEstimate(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================
  // Landau free-energy profile: -ln(P(lambda))/beta per bin. The histogram is sampled under the bias eta_i
  // (P_biased(i) ~ P(i) exp(eta_i), see weight()), so the unbiased profile is -ln(count_i)/beta + eta_i/beta;
  // the bias was left out, which made the profile flat (the Wang-Landau bias flattens the sampled
  // histogram) instead of showing the barrier. The excess chemical potential below adds the same
  // eta_last - eta_0 difference.

  static std::vector<double> landauFreeEnergyTransform(const std::vector<double>& bins, const std::vector<double>& bias,
                                                       double beta)
  {
    std::vector<double> freeEnergy(bins.size());
    for (std::size_t i = 0; i < bins.size(); ++i)
    {
      // an unsampled bin has no estimate: left at zero rather than the infinity of log(0)
      freeEnergy[i] = bins[i] > 0.0 ? (-std::log(bins[i]) + bias[i]) / beta : 0.0;
    }
    return freeEnergy;
  }

  std::vector<double> averagedLandauFreeEnergyHistogram(std::size_t blockIndex, double beta) const
  {
    return landauFreeEnergyTransform(averagedProbabilityHistogram(blockIndex), biasFactor, beta);
  }

  std::vector<double> averagedLandauFreeEnergyHistogram(double beta) const
  {
    return landauFreeEnergyTransform(averagedProbabilityHistogram(), biasFactor, beta);
  }

  std::pair<std::vector<double>, std::vector<double>> averageLandauFreeEnergyHistogram(double beta) const
  {
    std::vector<double> average = averagedLandauFreeEnergyHistogram(beta);

    std::vector<std::vector<double>> blockAverages;
    blockAverages.reserve(numberOfBlocks);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampled(blockIndex)) blockAverages.push_back(averagedLandauFreeEnergyHistogram(blockIndex, beta));

    std::vector<double> confidenceIntervalError = blockErrorEstimate(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================
  // excess chemical potential from the histogram: -ln(P(lambda=1)/P(lambda=0))/beta of the sampled (biased)
  // histogram; the callers add the bias difference (eta_last - eta_0)/beta (see landauFreeEnergyTransform)

  static double excessChemicalPotentialTransform(const std::vector<double> &bins, double beta)
  {
    // no estimate without samples at both ends of the lambda path
    if (!(bins.front() > 0.0) || !(bins.back() > 0.0)) return 0.0;
    return -std::log(bins.back() / bins.front()) / beta;
  }

  double averagedExcessChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return excessChemicalPotentialTransform(bookKeepingLambda[blockIndex], beta);
  }

  double averagedExcessChemicalPotential(double beta) const
  {
    return excessChemicalPotentialTransform(summedLambdaBins(), beta);
  }

  std::pair<double, double> averageExcessChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotential(beta);

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtBothEnds(blockIndex)) blockAverages.push_back(averagedExcessChemicalPotential(blockIndex, beta));

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta, double bias) const
  {
    double average = averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta) + bias;

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtBothEnds(blockIndex))
        blockAverages.push_back(averagedExcessChemicalPotential(blockIndex, beta) +
                                averagedIdealGasChemicalPotential(blockIndex, beta) + bias);

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  /// The fugacity of a total chemical potential. Without molecules (no density) there is none: the ideal-gas
  /// chemical potential of an empty box is minus infinity, and with that term left at zero the excess term alone
  /// made an empty box report a fugacity of exp(beta mu_ex) / beta (as PropertyWidom::fugacityTransform did).
  static double fugacityOfChemicalPotential(double fluidDensity, double totalChemicalPotential, double beta)
  {
    return fluidDensity > 0.0 ? std::exp(beta * totalChemicalPotential) / beta : 0.0;
  }

  std::pair<double, double> averageFugacity(double beta, double bias) const
  {
    double average = fugacityOfChemicalPotential(
        density.averaged(),
        averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta) + bias, beta);

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtBothEnds(blockIndex))
        blockAverages.push_back(fugacityOfChemicalPotential(
            density.averaged(blockIndex),
            averagedExcessChemicalPotential(blockIndex, beta) + averagedIdealGasChemicalPotential(blockIndex, beta) +
                bias,
            beta));

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================
  // dU/dlambda (thermodynamic integration; per-bin weighted averages)

  /// Per-bin average of (sum, count) accumulators.
  static std::vector<double> binAverages(const std::vector<std::pair<double, double>> &bins)
  {
    std::vector<double> averagedData(bins.size());
    std::transform(bins.cbegin(), bins.cend(), averagedData.begin(), [](const std::pair<double, double> &sample)
                   { return sample.first / (sample.second > 0.0 ? sample.second : 1.0); });
    return averagedData;
  }

  /// dU/dlambda accumulators pooled over all blocks.
  std::vector<std::pair<double, double>> summedDUdlambdaBins() const
  {
    std::vector<std::pair<double, double>> summedBlocks(numberOfSamplePoints);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingDUdlambda[blockIndex].begin(),
                     summedBlocks.begin(), pairSum<double, double>);
    }
    return summedBlocks;
  }

  std::vector<double> averagedDUdlambda(std::size_t blockIndex) const
  {
    return binAverages(bookKeepingDUdlambda[blockIndex]);
  }

  std::vector<double> averagedDUdlambda() const { return binAverages(summedDUdlambdaBins()); }

  std::pair<std::vector<double>, std::vector<double>> averageDuDlambda() const
  {
    std::vector<double> average = averagedDUdlambda();

    std::vector<std::vector<double>> blockAverages;
    blockAverages.reserve(numberOfBlocks);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampled(blockIndex)) blockAverages.push_back(averagedDUdlambda(blockIndex));

    std::vector<double> confidenceIntervalError = blockErrorEstimate(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================
  // excess chemical potential from thermodynamic integration: Simpson's rule over <dU/dlambda>

  /// Composite Simpson's rule over the lambda grid (spacing delta). The 1/3 rule requires an even
  /// number of intervals; for an odd number of intervals the last three are integrated with the
  /// 3/8 rule (same fourth-order accuracy), so any grid size is handled exactly.
  double simpsonIntegral(const std::vector<double> &data) const
  {
    const std::size_t n = data.size();
    if (n < 2uz) return 0.0;
    if (n == 2uz) return 0.5 * delta * (data[0] + data[1]);

    const std::size_t intervals = n - 1uz;
    const std::size_t oneThirdIntervals = (intervals % 2uz == 0uz) ? intervals : intervals - 3uz;

    double sum = 0.0;
    if (oneThirdIntervals > 0uz)
    {
      double weighted = data[0] + data[oneThirdIntervals];
      for (std::size_t i = 1; i < oneThirdIntervals; ++i)
      {
        weighted += (i % 2uz != 0uz) ? 4.0 * data[i] : 2.0 * data[i];
      }
      sum += weighted * (delta / 3.0);
    }

    if (intervals % 2uz != 0uz)
    {
      const std::size_t s = oneThirdIntervals;
      sum += (3.0 * delta / 8.0) * (data[s] + 3.0 * data[s + 1] + 3.0 * data[s + 2] + data[s + 3]);
    }

    return sum;
  }

  /// The blocks whose integral over lambda is an estimate: those with a dU/dlambda sample in every bin. A bin
  /// without one enters the integral as zero, so a block that visited part of the lambda path gave part of the
  /// integral (the blocks were taken when any bin was visited, or, for the total chemical potential and the
  /// fugacity, when the density was sampled).
  bool blockSampledAtEveryLambda(std::size_t blockIndex) const
  {
    return std::ranges::all_of(bookKeepingDUdlambda[blockIndex],
                               [](const std::pair<double, double> &sample) { return sample.second > 0.0; });
  }

  double averagedExcessChemicalPotentialDUdlambda(std::size_t blockIndex) const
  {
    return simpsonIntegral(averagedDUdlambda(blockIndex));
  }

  double averagedExcessChemicalPotentialDUdlambda() const { return simpsonIntegral(averagedDUdlambda()); }

  std::pair<double, double> averageExcessChemicalPotentialDUdlambda() const
  {
    double average = averagedExcessChemicalPotentialDUdlambda();

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtEveryLambda(blockIndex))
        blockAverages.push_back(averagedExcessChemicalPotentialDUdlambda(blockIndex));

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotentialDUdlambda() + averagedIdealGasChemicalPotential(beta);

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtEveryLambda(blockIndex))
        blockAverages.push_back(averagedExcessChemicalPotentialDUdlambda(blockIndex) +
                                averagedIdealGasChemicalPotential(blockIndex, beta));

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacityDUdlambda(double beta) const
  {
    double average = fugacityOfChemicalPotential(
        density.averaged(), averagedExcessChemicalPotentialDUdlambda() + averagedIdealGasChemicalPotential(beta),
        beta);

    std::vector<double> blockAverages;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
      if (blockSampledAtEveryLambda(blockIndex))
        blockAverages.push_back(fugacityOfChemicalPotential(
            density.averaged(blockIndex),
            averagedExcessChemicalPotentialDUdlambda(blockIndex) + averagedIdealGasChemicalPotential(blockIndex, beta),
            beta));

    double confidenceIntervalError = blockErrorEstimate<double>(blockAverages, average);

    return std::make_pair(average, confidenceIntervalError);
  }

  void readBiasingFile(std::filesystem::path path);
  void writeBiasingFile(std::filesystem::path path);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyLambdaProbabilityHistogram &p);

  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyLambdaProbabilityHistogram &p);
};
