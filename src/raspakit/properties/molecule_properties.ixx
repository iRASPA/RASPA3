module;

export module property_molecule_properties;

import std;

import archive;
import atom;
import component;

// Samples intra-molecular geometry histograms (bond lengths, bend angles, and
// torsion/dihedral angles) for the flexible connectivity of every component.
//
// This mirrors the "molecule properties" analysis from RASPA2
// (src/molecule_properties.c): for each component the bond, bend and torsion
// distributions are accumulated into probability histograms and written to disk.

export struct PropertyMoleculeProperties
{
  PropertyMoleculeProperties() {};

  PropertyMoleculeProperties(std::size_t numberOfBlocks, const std::vector<Component> &components,
                             std::size_t numberOfBins, double bondRange, std::size_t sampleEvery,
                             std::optional<std::size_t> writeEvery);

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks{0};
  std::size_t numberOfBins{0};
  std::size_t numberOfComponents{0};

  double bondRange{4.0};    ///< Upper limit of the bond-length histogram [Angstrom].
  double bendRange{180.0};  ///< Upper limit of the bend-angle histogram [degrees].
  double torsionRange{360.0};  ///< Full width of the torsion-angle histogram [degrees] (range [-180, 180]).

  double deltaBond{0.0};     ///< Bin width for bond-length histograms [Angstrom].
  double deltaBend{0.0};     ///< Bin width for bend-angle histograms [degrees].
  double deltaTorsion{0.0};  ///< Bin width for torsion-angle histograms [degrees].

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};

  std::vector<std::size_t> numberOfBondsPerComponent{};
  std::vector<std::size_t> numberOfBendsPerComponent{};
  std::vector<std::size_t> numberOfTorsionsPerComponent{};

  // Histograms indexed as [block][component][potentialIndex][bin].
  std::vector<std::vector<std::vector<std::vector<double>>>> bondHistogram{};
  std::vector<std::vector<std::vector<std::vector<double>>>> bendHistogram{};
  std::vector<std::vector<std::vector<std::vector<double>>>> torsionHistogram{};

  // Number of molecule-samples per [block][component]; identical for every
  // potential of a given component, used as the normalization factor.
  std::vector<std::vector<double>> numberOfCounts{};
  double totalNumberOfCounts{0.0};

  void sample(const std::vector<Component> &components,
              const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
              std::size_t currentCycle, std::size_t block);

  // Returns {binCenters, averageProbabilityDensity, confidenceIntervalError} for
  // one potential instance of one component.
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> result(
      const std::vector<std::vector<std::vector<std::vector<double>>>> &histogram, std::size_t component,
      std::size_t index, double delta, double rangeStart) const;

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyMoleculeProperties &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMoleculeProperties &p);
};
