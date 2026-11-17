module;

export module property_rdf;

import std;

import archive;
import atom;
import double3;
import atom_dynamics;
import molecule;
import simulationbox;
import forcefield;
export import property_block_average;

// Force-based (Borgis) radial distribution function from known site gradients.
// Does not evaluate forces: callers must provide current full-U gradients in atomDynamics
// (MD integrator, or System::sampleForceBasedRDFWithFullGradients on the MC path).

export struct PropertyRadialDistributionFunction
{
  PropertyRadialDistributionFunction() {};

  PropertyRadialDistributionFunction(std::size_t numberOfBlocks, std::size_t numberOfPseudoAtoms,
                                     std::size_t numberOfBins, double range, std::size_t sampleEvery,
                                     std::size_t writeEvery)
      : numberOfPseudoAtoms(numberOfPseudoAtoms),
        numberOfBins(numberOfBins),
        range(range),
        deltaR(range / static_cast<double>(numberOfBins)),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        histogram(numberOfBlocks, numberOfPseudoAtoms * numberOfPseudoAtoms, numberOfBins),
        pairCount(numberOfPseudoAtoms * numberOfPseudoAtoms),
        pairCountPerBlock(numberOfBlocks, std::vector<double>(numberOfPseudoAtoms * numberOfPseudoAtoms))
  {
  }

  std::uint64_t versionNumber{3};

  /// Channel index of an ordered pseudo-atom pair.
  std::size_t channel(std::size_t atomTypeA, std::size_t atomTypeB) const
  {
    return atomTypeB + atomTypeA * numberOfPseudoAtoms;
  }

  /// g(r) of one block and of all blocks: the accumulated force sums over the accumulated number of pairs.
  std::vector<double> averagedProbabilityHistogram(std::size_t blockIndex, std::size_t atomTypeA,
                                                   std::size_t atomTypeB) const;
  std::vector<double> averagedProbabilityHistogram(std::size_t atomTypeA, std::size_t atomTypeB) const;
  /// g(r) and its error: g = <beta V sum / 8 pi> / <P>, the force integral of the pairs closer than r (Borgis et
  /// al.) over the number of pairs, each averaged over the samples. Integrated from g(0) = 0: no assumption on g at
  /// the end of the range.
  std::pair<std::vector<double>, std::vector<double>> result(std::size_t atomTypeA, std::size_t atomTypeB) const;

  std::size_t numberOfPseudoAtoms{};
  std::size_t numberOfBins{};
  double range{};
  double deltaR{};
  std::size_t sampleEvery{};
  std::size_t writeEvery{};
  BlockHistogram<double> histogram;
  std::vector<std::size_t> pairCount;
  /// The number of ordered pairs sampled per block and channel: the denominator of g(r).
  std::vector<std::vector<double>> pairCountPerBlock;

  /// 'frameworkGradients' and 'moleculeGradients' are the gradients of the bodies the sites move with (the atom
  /// itself, or its rigid molecule, fragment or group); a framework site whose 'frameworkMobile' entry is zero does
  /// not move and contributes the force on its partner only.
  void sample(const SimulationBox &simulationBox, double beta, std::span<const Atom> frameworkAtoms,
              std::span<const double3> frameworkGradients, std::span<const std::uint8_t> frameworkMobile,
              std::span<const Atom> moleculeAtoms, std::span<const double3> moleculeGradients,
              std::size_t currentCycle, std::size_t block);
  void writeOutput(const ForceField &forceField, std::size_t systemId, double volume,
                   std::vector<std::size_t> &numberOfPseudoAtomsType, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyRadialDistributionFunction &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadialDistributionFunction &temp);
};
