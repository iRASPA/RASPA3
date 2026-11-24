module;

export module property_conventional_rdf;

import std;

import archive;
import atom;
import simulationbox;
import forcefield;
export import property_block_average;

// Computes Radial Distribution Function
// Also works correctly for a small number of molecules (RDF still goes to unity)

export struct PropertyConventionalRadialDistributionFunction
{
  PropertyConventionalRadialDistributionFunction() {};

  PropertyConventionalRadialDistributionFunction(std::size_t numberOfBins, double range, std::size_t sampleEvery,
                                                 std::optional<std::size_t> writeEvery):
        numberOfBins(numberOfBins),
        range(range),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
    // an interval of zero is a division by zero at the first sample
    if (sampleEvery == 0uz)
      throw std::runtime_error("PropertyConventionalRadialDistributionFunction: the sampling interval must be at least 1\n");
    if (writeEvery.has_value() && writeEvery.value() == 0uz)
      throw std::runtime_error("PropertyConventionalRadialDistributionFunction: the writing interval must be at least 1\n");
  }

  PropertyConventionalRadialDistributionFunction(std::size_t numberOfBlocks, std::size_t numberOfPseudoAtoms,
                                                 std::size_t numberOfBins, double range, std::size_t sampleEvery,
                                                 std::optional<std::size_t> writeEvery)
      : numberOfPseudoAtoms(numberOfPseudoAtoms),
        numberOfBins(numberOfBins),
        range(range),
        deltaR(range / static_cast<double>(numberOfBins)),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        histogram(numberOfBlocks, numberOfPseudoAtoms * numberOfPseudoAtoms, numberOfBins),
        pairCount(numberOfPseudoAtoms * numberOfPseudoAtoms),
        pairWeight(numberOfPseudoAtoms * numberOfPseudoAtoms),
        pairWeightPerBlock(numberOfBlocks, std::vector<double>(numberOfPseudoAtoms * numberOfPseudoAtoms))
  {
    // an interval of zero is a division by zero at the first sample
    if (sampleEvery == 0uz)
      throw std::runtime_error("PropertyConventionalRadialDistributionFunction: the sampling interval must be at least 1\n");
    if (writeEvery.has_value() && writeEvery.value() == 0uz)
      throw std::runtime_error("PropertyConventionalRadialDistributionFunction: the writing interval must be at least 1\n");
  }

  std::uint64_t versionNumber{3};

  /// Channel index of an ordered pseudo-atom pair.
  std::size_t channel(std::size_t atomTypeA, std::size_t atomTypeB) const
  {
    return atomTypeB + atomTypeA * numberOfPseudoAtoms;
  }

  std::vector<double> averagedProbabilityHistogram(std::size_t blockIndex, std::size_t atomTypeA,
                                                   std::size_t atomTypeB) const
  {
    return histogram.averaged(blockIndex, channel(atomTypeA, atomTypeB));
  }
  std::vector<double> averagedProbabilityHistogram(std::size_t atomTypeA, std::size_t atomTypeB) const
  {
    return histogram.averaged(channel(atomTypeA, atomTypeB));
  }
  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram(std::size_t atomTypeA,
                                                                                  std::size_t atomTypeB) const
  {
    return histogram.average(channel(atomTypeA, atomTypeB));
  }

  /// g(r) of a pair of types and its error, per bin, for the given volume.
  std::pair<std::vector<double>, std::vector<double>> normalized(std::size_t atomTypeA, std::size_t atomTypeB,
                                                                 double averageVolume) const;

  std::vector<std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>> result() const;

  std::size_t numberOfPseudoAtoms{};
  std::size_t numberOfBins{};
  double range{};
  double deltaR{};
  std::size_t sampleEvery{};
  std::optional<std::size_t> writeEvery;
  BlockHistogram<double> histogram;
  std::vector<std::size_t> pairCount;
  /// The ordered pairs sampled per channel, each with the weight of its configuration: the denominator of g(r).
  std::vector<double> pairWeight;
  /// The same per block, [block][channel]: a block is normalized by its own pairs (see 'normalized').
  std::vector<std::vector<double>> pairWeightPerBlock;
  double volumeCumulative{};
  double volumeSampleCount{};

  void sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms, std::span<Atom> moleculeAtoms,
              std::size_t currentCycle, std::size_t block, double weight = 1.0);

  void writeOutput(const ForceField &forceField, std::size_t systemId, double volume,
                   std::vector<std::size_t> &numberOfPseudoAtomsType, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyConventionalRadialDistributionFunction &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive,
                                            PropertyConventionalRadialDistributionFunction &temp);
};
