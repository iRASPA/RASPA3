module;

module property_conventional_rdf;

import std;

import archive;
import double3;
import atom;
import simulationbox;
import forcefield;
import averages;

void PropertyConventionalRadialDistributionFunction::sample(const SimulationBox &simulationBox,
                                                            std::span<Atom> frameworkAtoms,
                                                            std::span<Atom> moleculeAtoms, std::size_t currentCycle,
                                                            std::size_t block)
{
  double3 dr, posA, posB, f;
  double rr, r;

  if (currentCycle % sampleEvery != 0uz) return;

  if (moleculeAtoms.empty()) return;

  volumeCumulative += simulationBox.volume;
  ++volumeSampleCount;


  for (std::span<Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    for (std::span<Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);

      pairCount[channel(typeA, typeB)]++;
      pairCount[channel(typeB, typeA)]++;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      r = std::sqrt(rr);

      std::size_t bin = static_cast<std::size_t>(r / deltaR);
      if (bin < numberOfBins)
      {
        histogram(block, channel(typeA, typeB), bin) += 1.0;
        histogram(block, channel(typeB, typeA), bin) += 1.0;
      }
    }
  }

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    std::size_t compA = static_cast<std::size_t>(it1->componentId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);

    for (std::span<Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);
      std::size_t compB = static_cast<std::size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);

        pairCount[channel(typeA, typeB)]++;
        pairCount[channel(typeB, typeA)]++;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);
        r = std::sqrt(rr);

        std::size_t bin = static_cast<std::size_t>(r / deltaR);
        if (bin < numberOfBins)
        {
          histogram(block, channel(typeA, typeB), bin) += 1.0;
          histogram(block, channel(typeB, typeA), bin) += 1.0;
        }
      }
    }
  }

  histogram.addCount(block, 2.0);
}

std::vector<std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>>
PropertyConventionalRadialDistributionFunction::result() const
{
   std::vector<std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>> 
     results(numberOfPseudoAtoms, std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>(numberOfPseudoAtoms));

  for (std::size_t atomTypeA = 0; atomTypeA < numberOfPseudoAtoms; ++atomTypeA)
  {
    for (std::size_t atomTypeB = atomTypeA; atomTypeB < numberOfPseudoAtoms; ++atomTypeB)
    {
      if (pairCount[channel(atomTypeB, atomTypeA)] > 0)
      {
        auto [average, error] = averageProbabilityHistogram(atomTypeA, atomTypeB);
        
        // n_pairs is the number of unique pairs of atoms where one atom is from each of two sets
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
        double avg_n_pairs = static_cast<double>(pairCount[channel(atomTypeA, atomTypeB)] +
                                                 pairCount[channel(atomTypeB, atomTypeA)]) /
                             histogram.totalNumberOfCounts;
        double normalization = (volumeCumulative / volumeSampleCount) / (2.0 * std::numbers::pi * deltaR * deltaR * deltaR * avg_n_pairs);
        
        std::vector<double> x(numberOfBins);
        std::vector<double> y(numberOfBins);
        std::vector<double> error_y(numberOfBins);
        for (std::size_t bin = 0; bin != numberOfBins; ++bin)
        {
          x[bin] = (static_cast<double>(bin) + 0.5) * deltaR,
          y[bin] = average[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5)),
          error_y[bin] = error[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5));
        }

        results[atomTypeA][atomTypeB] = {x, y, error_y};
        results[atomTypeB][atomTypeA] = {x, y, error_y};
      }
    }
  }
   return results;
}

void PropertyConventionalRadialDistributionFunction::writeOutput(
    const ForceField &forceField, std::size_t systemId, double volume,
    [[maybe_unused]] std::vector<std::size_t> &numberOfPseudoAtomsType, std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  std::filesystem::create_directory("conventional_rdf");

  for (std::size_t atomTypeA = 0; atomTypeA < numberOfPseudoAtoms; ++atomTypeA)
  {
    for (std::size_t atomTypeB = atomTypeA; atomTypeB < numberOfPseudoAtoms; ++atomTypeB)
    {
      if (pairCount[channel(atomTypeB, atomTypeA)] > 0)
      {
        std::ofstream stream_rdf_output(std::format("conventional_rdf/rdf_{}_{}.s{}.txt",
                                                    forceField.pseudoAtoms[atomTypeA].name,
                                                    forceField.pseudoAtoms[atomTypeB].name, systemId));

        stream_rdf_output << std::format("# rdf, number of counts: {}\n", histogram.totalNumberOfCounts);
        stream_rdf_output << "# column 1: distance []\n";
        stream_rdf_output << "# column 2: normalize rdf []\n";
        stream_rdf_output << "# column 3: error normalize rdf []\n";

        auto [average, error] = averageProbabilityHistogram(atomTypeA, atomTypeB);

        // n_pairs is the number of unique pairs of atoms where one atom is from each of two sets
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
        double avg_n_pairs = static_cast<double>(pairCount[channel(atomTypeA, atomTypeB)] +
                                                 pairCount[channel(atomTypeB, atomTypeA)]) /
                             histogram.totalNumberOfCounts;
        double normalization = volume / (2.0 * std::numbers::pi * deltaR * deltaR * deltaR * avg_n_pairs);

        for (std::size_t bin = 0; bin != numberOfBins; ++bin)
        {
          stream_rdf_output << std::format(
              "{} {} {}\n", (static_cast<double>(bin) + 0.5) * deltaR,
              average[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5)),
              error[bin] * normalization / ((static_cast<double>(bin) + 0.5) * (static_cast<double>(bin) + 0.5)));
        }
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                   const PropertyConventionalRadialDistributionFunction &rdf)
{
  archive << rdf.versionNumber;

  archive << rdf.numberOfPseudoAtoms;
  archive << rdf.numberOfBins;
  archive << rdf.range;
  archive << rdf.deltaR;
  archive << rdf.sampleEvery;
  archive << rdf.writeEvery;
  archive << rdf.histogram;
  archive << rdf.pairCount;
  archive << rdf.volumeCumulative;
  archive << rdf.volumeSampleCount;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyConventionalRadialDistributionFunction &rdf)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > rdf.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyConventionalRadialDistributionFunction' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> rdf.numberOfPseudoAtoms;
  archive >> rdf.numberOfBins;
  archive >> rdf.range;
  archive >> rdf.deltaR;
  archive >> rdf.sampleEvery;
  archive >> rdf.writeEvery;
  archive >> rdf.histogram;
  archive >> rdf.pairCount;
  archive >> rdf.volumeCumulative;
  archive >> rdf.volumeSampleCount;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyConventionalRadialDistributionFunction: Error in binary restart\n"));
  }
#endif

  return archive;
}
