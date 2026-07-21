module;

export module property_rdf;

import std;

import archive;
import atom;
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
        pairCount(numberOfPseudoAtoms * numberOfPseudoAtoms)
  {
  }

  std::uint64_t versionNumber{1};

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
  std::pair<std::vector<double>, std::vector<double>> result(std::size_t atomTypeA, std::size_t atomTypeB) const
  {
    return histogram.average(channel(atomTypeA, atomTypeB));
  }

  std::size_t numberOfPseudoAtoms{};
  std::size_t numberOfBins{};
  double range{};
  double deltaR{};
  std::size_t sampleEvery{};
  std::size_t writeEvery{};
  BlockHistogram<double> histogram;
  std::vector<std::size_t> pairCount;

  void sample(const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
              std::span<const AtomDynamics> frameworkDynamics, const std::vector<Molecule> &molecules,
              std::span<const Atom> moleculeAtoms, std::span<const AtomDynamics> moleculeDynamics,
              std::size_t currentCycle, std::size_t block);
  void writeOutput(const ForceField &forceField, std::size_t systemId, double volume,
                   std::vector<std::size_t> &numberOfPseudoAtomsType, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyRadialDistributionFunction &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadialDistributionFunction &temp);
};
