module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_conventional_rdf;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <span>;
import <tuple>;
#endif

import archive;
import atom;
import simulationbox;
import forcefield;

// Computes Radial Distribution Function
// Also works correctly for a small number of molecules (RDF still goes to unity)

export struct PropertyConventionalRadialDistributionFunction
{
  PropertyConventionalRadialDistributionFunction() {};

  PropertyConventionalRadialDistributionFunction(size_t numberOfBlocks, size_t numberOfPseudoAtoms, size_t numberOfBins,
                                                 double range, size_t sampleEvery, size_t writeEvery)
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

  uint64_t versionNumber{1};

  std::vector<double> averagedProbabilityHistogram(size_t blockIndex, size_t atomTypeA, size_t atomTypeB) const;
  std::vector<double> averagedProbabilityHistogram(size_t atomTypeA, size_t atomTypeB) const;
  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram(size_t atomTypeA,
                                                                                  size_t atomTypeB) const;

  size_t numberOfBlocks;
  size_t numberOfPseudoAtoms;
  size_t numberOfPseudoAtomsSymmetricMatrix;
  size_t numberOfBins;
  double range;
  double deltaR;
  size_t sampleEvery;
  size_t writeEvery;
  std::vector<std::vector<double>> sumProperty;
  size_t totalNumberOfCounts;
  std::vector<size_t> numberOfCounts;
  std::vector<size_t> pairCount;

  void sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms, std::span<Atom> moleculeAtoms,
              size_t currentCycle, size_t block);
  void writeOutput(const ForceField &forceField, size_t systemId, double volume,
                   std::vector<size_t> &numberOfPseudoAtomsType, size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyConventionalRadialDistributionFunction &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive,
                                            PropertyConventionalRadialDistributionFunction &temp);
};
