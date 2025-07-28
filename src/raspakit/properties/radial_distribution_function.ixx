module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_rdf;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import atom;
import molecule;
import molecule;
import simulationbox;
import forcefield;

// Computes Radial Distribution Function
// Also works correctly for a small number of molecules (RDF still goes to unity)

export struct PropertyRadialDistributionFunction
{
  PropertyRadialDistributionFunction() {};

  PropertyRadialDistributionFunction(std::size_t numberOfBlocks, std::size_t numberOfPseudoAtoms,
                                     std::size_t numberOfBins, double range, std::size_t sampleEvery,
                                     std::size_t writeEvery)
      : numberOfBlocks(numberOfBlocks),
        numberOfPseudoAtoms(numberOfPseudoAtoms),
        numberOfBins(numberOfBins),
        range(range),
        deltaR(range / static_cast<double>(numberOfBins)),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        sumProperty(
            std::vector(numberOfBlocks * numberOfPseudoAtoms * numberOfPseudoAtoms, std::vector<double>(numberOfBins))),
        totalNumberOfCounts(0uz),
        numberOfCounts(numberOfBlocks),
        pairCount(numberOfPseudoAtoms * numberOfPseudoAtoms)
  {
  }

  std::uint64_t versionNumber{1};

  std::vector<double> averagedProbabilityHistogram(std::size_t blockIndex, std::size_t atomTypeA,
                                                   std::size_t atomTypeB) const;
  std::vector<double> averagedProbabilityHistogram(std::size_t atomTypeA, std::size_t atomTypeB) const;
  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram(std::size_t atomTypeA,
                                                                                  std::size_t atomTypeB) const;

  std::size_t numberOfBlocks;
  std::size_t numberOfPseudoAtoms;
  std::size_t numberOfPseudoAtomsSymmetricMatrix;
  std::size_t numberOfBins;
  double range;
  double deltaR;
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::vector<std::vector<double>> sumProperty;
  std::size_t totalNumberOfCounts;
  std::vector<std::size_t> numberOfCounts;
  std::vector<std::size_t> pairCount;

  void sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
              const std::vector<Molecule> &molecules, std::span<Atom> moleculeAtoms, std::size_t currentCycle,
              std::size_t block);
  void writeOutput(const ForceField &forceField, std::size_t systemId, double volume,
                   std::vector<std::size_t> &numberOfPseudoAtomsType, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyRadialDistributionFunction &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadialDistributionFunction &temp);
};
