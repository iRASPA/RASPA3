export module property_rdf;

import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <span>;

import atom;
import simulationbox;


export struct PropertyRadialDistributionFunction
{
  PropertyRadialDistributionFunction() {};

  PropertyRadialDistributionFunction(size_t numberOfBlocks, size_t numberOfPseudoAtoms, size_t numberOfBins, double range) :
    numberOfBlocks(numberOfBlocks),
    numberOfPseudoAtoms(numberOfPseudoAtoms),
    numberOfPseudoAtomsSymmetricMatrix(numberOfPseudoAtoms * (numberOfPseudoAtoms + 1) / 2),
    numberOfBins(numberOfBins),
    range(range),
    deltaR(range / static_cast<double>(numberOfBins)),
    sumProperty(numberOfBlocks * numberOfPseudoAtomsSymmetricMatrix * numberOfBins),
    numberOfCounts(0uz)
  {
  }

  double& operator[](size_t block, size_t i, size_t j, size_t bin)
  {
    if (i < j) std::swap(i, j);
    size_t index_pseudo_atoms = i * (i + 1) / 2 + j;
    size_t index = bin + numberOfBins * (index_pseudo_atoms + block * numberOfPseudoAtomsSymmetricMatrix);
    return sumProperty[index];
  }

  const double& operator[](size_t block, size_t i, size_t j, size_t bin) const
  {
    if (i < j) std::swap(i, j);
    size_t index_pseudo_atoms = i * (i + 1) / 2 + j;
    size_t index = bin + numberOfBins * (index_pseudo_atoms + block * numberOfPseudoAtomsSymmetricMatrix);
    return sumProperty[index];
  }

  size_t numberOfBlocks;
  size_t numberOfPseudoAtoms;
  size_t numberOfPseudoAtomsSymmetricMatrix;
  size_t numberOfBins;
  double range;
  double deltaR;
  std::vector<double> sumProperty;
  size_t numberOfCounts;

  void sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms, std::span<Atom> moleculeAtoms, size_t block);
  void writeOutput(size_t systemId);
};
