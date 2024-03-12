module;

module property_conventional_rdf;

import <cstddef>;
import <string>;
import <iostream>;
import <fstream>;
import <sstream>;
import <tuple>;
import <vector>;
import <algorithm>;
import <format>;
import <numbers>;

import double3;
import atom;
import simulationbox;
import forcefield;
import averages;

void PropertyConventionalRadialDistributionFunction::sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms, std::span<Atom> moleculeAtoms, size_t block)
{
  double3 dr, posA, posB, f;
  double rr, r;

  if(moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);

      pairCount[typeB + typeA * numberOfPseudoAtoms]++;
      pairCount[typeA + typeB * numberOfPseudoAtoms]++;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      r = std::sqrt(rr);

      size_t bin = static_cast<size_t>(r / deltaR);
      if(bin < numberOfBins)
      {
        size_t index = typeB + typeA * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
        sumProperty[index][bin] += 1.0;
        index = typeA + typeB * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
        sumProperty[index][bin] += 1.0;
      }
    }
  }

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);

    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);

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

        size_t bin = static_cast<size_t>(r / deltaR);
        if(bin < numberOfBins)
        {
          size_t index = typeB + typeA * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
          sumProperty[index][bin] += 1.0;
          index = typeA + typeB * numberOfPseudoAtoms + block * numberOfPseudoAtoms * numberOfPseudoAtoms;
          sumProperty[index][bin] += 1.0;
        }
      }
    }
  }

  totalNumberOfCounts += 2;
  numberOfCounts[block] += 2;
}


std::vector<double> PropertyConventionalRadialDistributionFunction::averagedProbabilityHistogram(size_t blockIndex, size_t atomTypeA, size_t atomTypeB) const
{
  size_t index_pseudo_atoms = atomTypeB + atomTypeA * numberOfPseudoAtoms;

  std::vector<double> averagedData(numberOfBins);
  std::transform(sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].begin(), 
                 sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].end(), averagedData.begin(),
                 [&](const double &sample){return sample / static_cast<double>(numberOfCounts[blockIndex]);});
  return averagedData;
}

std::vector<double> PropertyConventionalRadialDistributionFunction::averagedProbabilityHistogram(size_t atomTypeA, size_t atomTypeB) const
{
  size_t index_pseudo_atoms = atomTypeB + atomTypeA * numberOfPseudoAtoms;

  std::vector<double> summedBlocks(numberOfBins);
  for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::transform(summedBlocks.begin(), summedBlocks.end(), 
                   sumProperty[index_pseudo_atoms + blockIndex * numberOfPseudoAtoms * numberOfPseudoAtoms].begin(), summedBlocks.begin(),
                   [](const double & a, const double & b){ return a + b; });
  }
  std::vector<double> average(numberOfBins);
  std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double &sample){return sample / static_cast<double>(totalNumberOfCounts);});

  return average;
}

std::pair<std::vector<double>, std::vector<double>> 
PropertyConventionalRadialDistributionFunction::averageProbabilityHistogram(size_t atomTypeA, size_t atomTypeB) const
{
  size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  std::vector<double> average = averagedProbabilityHistogram(atomTypeA, atomTypeB);

  std::vector<double> sumOfSquares(numberOfBins);
  for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex, atomTypeA, atomTypeB);
    for(size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
    {
      double value = blockAverage[binIndex] - average[binIndex];
      sumOfSquares[binIndex] += value * value;
    }
  }
  std::vector<double> standardDeviation(numberOfBins);
  std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                 [&](const double &sumofsquares)
                 {return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom));});

  std::vector<double> standardError(numberOfBins);
  std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                 [&](const double &sigma){return sigma / sqrt(static_cast<double>(numberOfBlocks));});

  std::vector<double> confidenceIntervalError(numberOfBins);
  std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                 [&](const double &error){return intermediateStandardNormalDeviate * error;});

  return std::make_pair(average, confidenceIntervalError);
}


void PropertyConventionalRadialDistributionFunction::writeOutput(const ForceField &forceField, size_t systemId, double volume,
                                                                 [[maybe_unused]]std::vector<size_t> &numberOfPseudoAtomsType)
{
  std::filesystem::create_directory("conventional_rdf");

  for(size_t atomTypeA = 0; atomTypeA < numberOfPseudoAtoms; ++atomTypeA)
  {
    for(size_t atomTypeB = atomTypeA; atomTypeB < numberOfPseudoAtoms; ++atomTypeB)
    {
      if(pairCount[atomTypeA + atomTypeB * numberOfPseudoAtoms] > 0)
      {
        std::ofstream stream_rdf_output(std::format("conventional_rdf/rdf_{}_{}_{}.data", 
              forceField.pseudoAtoms[atomTypeA].name,
              forceField.pseudoAtoms[atomTypeB].name,
              systemId));

        stream_rdf_output << std::format("# rdf, number of counts: {}\n", totalNumberOfCounts);
        stream_rdf_output << "# column 1: distance []\n";
        stream_rdf_output << "# column 2: normalize rdf []\n";
        stream_rdf_output << "# column 3: error normalize rdf []\n";

        auto [average, error] = averageProbabilityHistogram(atomTypeA, atomTypeB);

        // n_pairs is the number of unique pairs of atoms where one atom is from each of two sets
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
        double avg_n_pairs = static_cast<double>(pairCount[atomTypeB + atomTypeA * numberOfPseudoAtoms] +
                                                 pairCount[atomTypeA + atomTypeB * numberOfPseudoAtoms]) / static_cast<double>(totalNumberOfCounts);
        double normalization = volume / (2.0 * std::numbers::pi * deltaR * deltaR * deltaR * avg_n_pairs);

        for(size_t bin = 0; bin != numberOfBins; ++bin)
        {
          stream_rdf_output << std::format("{} {} {}\n", (static_cast<double>(bin) + 0.5) * deltaR, 
              average[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5)),
              error[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5)));
        }
      }
    }
  }
}

