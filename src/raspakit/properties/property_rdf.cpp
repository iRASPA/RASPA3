module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#endif

module property_rdf;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <string>;
import <iostream>;
import <fstream>;
import <sstream>;
import <tuple>;
import <vector>;
import <algorithm>;
import <complex>;
import <format>;
import <numbers>;
import <span>;
import <array>;
import <cmath>;
import <exception>;
import <source_location>;
import <print>;
#endif

import archive;
import double3;
import atom;
import simulationbox;
import forcefield;
import averages;

void PropertyRadialDistributionFunction::sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
                                                [[maybe_unused]] const std::vector<Molecule> &molecules,
                                                std::span<Atom> moleculeAtoms, size_t currentCycle, size_t block)
{
  double3 dr, posA, posB, f;
  double3 gradientA, gradientB;
  double rr, r;

  if (currentCycle % sampleEvery != 0uz) return;

  if (moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    gradientA = it1->gradient;
    size_t typeA = static_cast<size_t>(it1->type);
    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      gradientB = it2->gradient;
      size_t typeB = static_cast<size_t>(it2->type);

      pairCount[typeB + typeA * numberOfPseudoAtoms]++;
      pairCount[typeA + typeB * numberOfPseudoAtoms]++;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      r = std::sqrt(rr);

      double value = double3::dot(gradientA - gradientB, dr) / (rr * r);

      size_t bin = static_cast<size_t>(r / deltaR);
      for (size_t i = 0; i < bin; ++i)
      {
        size_t index = typeB + typeA * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
        sumProperty[index][bin] += value;
        index = typeA + typeB * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
        sumProperty[index][bin] += value;
      }
    }
  }

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    gradientA = it1->gradient;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);
      gradientB = it2->gradient;

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);

        pairCount[typeB + typeA * numberOfPseudoAtoms]++;
        pairCount[typeA + typeB * numberOfPseudoAtoms]++;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);
        r = std::sqrt(rr);

        if (r < range)
        {
          double value = double3::dot(gradientA - gradientB, dr) / (rr * r);

          for (size_t i = 0; i < numberOfBins; ++i)
          {
            if ((static_cast<double>(i) + 0.5) * deltaR < r)
            {
              size_t index = typeB + typeA * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
              sumProperty[index][i] += value;
              index = typeA + typeB * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
              sumProperty[index][i] += value;
            }
          }
        }
      }
    }
  }

  /*
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    gradientA = it1->gradient;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);

    (atom.position.x - com.x) * atom.gradient.x
  }
  */

  totalNumberOfCounts += 2;
  numberOfCounts[block] += 2;
}

std::vector<double> PropertyRadialDistributionFunction::averagedProbabilityHistogram(size_t blockIndex,
                                                                                     size_t atomTypeA,
                                                                                     size_t atomTypeB) const
{
  size_t index_pseudo_atoms = atomTypeB + atomTypeA * numberOfPseudoAtoms;

  std::vector<double> averagedData(numberOfBins);
  std::transform(sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].begin(),
                 sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].end(),
                 averagedData.begin(),
                 [&](const double &sample) { return sample / static_cast<double>(numberOfCounts[blockIndex]); });
  return averagedData;
}

std::vector<double> PropertyRadialDistributionFunction::averagedProbabilityHistogram(size_t atomTypeA,
                                                                                     size_t atomTypeB) const
{
  size_t index_pseudo_atoms = atomTypeB + atomTypeA * numberOfPseudoAtoms;

  std::vector<double> summedBlocks(numberOfBins);
  for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(),
                   sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].begin(),
                   summedBlocks.begin(), [](const double &a, const double &b) { return a + b; });
  }
  std::vector<double> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double &sample) { return sample / static_cast<double>(totalNumberOfCounts); });

  return average;
}

std::pair<std::vector<double>, std::vector<double>> PropertyRadialDistributionFunction::averageProbabilityHistogram(
    size_t atomTypeA, size_t atomTypeB) const
{
  size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  std::vector<double> average = averagedProbabilityHistogram(atomTypeA, atomTypeB);

  std::vector<double> sumOfSquares(numberOfBins);
  for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex, atomTypeA, atomTypeB);
    for (size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
    {
      double value = blockAverage[binIndex] - average[binIndex];
      sumOfSquares[binIndex] += value * value;
    }
  }
  std::vector<double> standardDeviation(numberOfBins);
  std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(), [&](const double &sumofsquares)
                 { return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

  std::vector<double> standardError(numberOfBins);
  std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                 [&](const double &sigma) { return sigma / std::sqrt(static_cast<double>(numberOfBlocks)); });

  std::vector<double> confidenceIntervalError(numberOfBins);
  std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                 [&](const double &error) { return intermediateStandardNormalDeviate * error; });

  return std::make_pair(average, confidenceIntervalError);
}

void PropertyRadialDistributionFunction::writeOutput(const ForceField &forceField, size_t systemId,
                                                     [[maybe_unused]] double volume,
                                                     [[maybe_unused]] std::vector<size_t> &numberOfPseudoAtomsType,
                                                     size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("rdf");

  for (size_t atomTypeA = 0; atomTypeA < numberOfPseudoAtoms; ++atomTypeA)
  {
    for (size_t atomTypeB = atomTypeA; atomTypeB < numberOfPseudoAtoms; ++atomTypeB)
    {
      if (pairCount[atomTypeA + atomTypeB * numberOfPseudoAtoms] > 0)
      {
        std::ofstream stream_rdf_output(std::format("rdf/rdf_{}_{}.s{}.data", forceField.pseudoAtoms[atomTypeA].name,
                                                    forceField.pseudoAtoms[atomTypeB].name, systemId));

        stream_rdf_output << std::format("# rdf, number of counts: {}\n", totalNumberOfCounts);
        stream_rdf_output << "# column 1: distance []\n";
        stream_rdf_output << "# column 2: normalize rdf []\n";
        stream_rdf_output << "# column 3: error normalize rdf []\n";

        auto [average, error] = averageProbabilityHistogram(atomTypeA, atomTypeB);

        // n_pairs is the number of unique pairs of atoms where one atom is from each of two sets
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
        [[maybe_unused]] double avg_n_pairs =
            static_cast<double>(pairCount[atomTypeB + atomTypeA * numberOfPseudoAtoms] +
                                pairCount[atomTypeA + atomTypeB * numberOfPseudoAtoms]) /
            static_cast<double>(totalNumberOfCounts);
        double normalization = 1.0 / std::abs(average[0]);

        for (size_t bin = 0; bin != numberOfBins; ++bin)
        {
          stream_rdf_output << std::format("{} {} {}\n", (static_cast<double>(bin) + 0.5) * deltaR,
                                           1.0 + average[bin] * normalization, error[bin] * normalization);
        }
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyRadialDistributionFunction &rdf)
{
  archive << rdf.versionNumber;

  archive << rdf.numberOfBlocks;
  archive << rdf.numberOfPseudoAtoms;
  archive << rdf.numberOfPseudoAtomsSymmetricMatrix;
  archive << rdf.numberOfBins;
  archive << rdf.range;
  archive << rdf.deltaR;
  archive << rdf.sampleEvery;
  archive << rdf.writeEvery;
  archive << rdf.sumProperty;
  archive << rdf.totalNumberOfCounts;
  archive << rdf.numberOfCounts;
  archive << rdf.pairCount;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadialDistributionFunction &rdf)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > rdf.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyRadialDistributionFunction' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> rdf.numberOfBlocks;
  archive >> rdf.numberOfPseudoAtoms;
  archive >> rdf.numberOfPseudoAtomsSymmetricMatrix;
  archive >> rdf.numberOfBins;
  archive >> rdf.range;
  archive >> rdf.deltaR;
  archive >> rdf.sampleEvery;
  archive >> rdf.writeEvery;
  archive >> rdf.sumProperty;
  archive >> rdf.totalNumberOfCounts;
  archive >> rdf.numberOfCounts;
  archive >> rdf.pairCount;

  return archive;
}
