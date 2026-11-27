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
                                                            std::size_t block, double weight)
{
  double3 dr, posA, posB, f;
  double rr, r;

  if (currentCycle % sampleEvery != 0uz) return;

  // a configuration without integer molecules is no sample (a fractional molecule is no member of the ensemble)
  if (std::ranges::none_of(moleculeAtoms, [](const Atom &atom) { return !atom.isFractional; })) return;

  // the counts, the pairs and the volume carry the weight of the configuration (one, unless a biased walk is
  // reweighted)
  volumeCumulative += weight * simulationBox.volume;
  volumeSampleCount += weight;


  for (std::span<Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    for (std::span<Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      // the distribution of the integer molecules: a fractional molecule is no member of the ensemble
      if (it2->isFractional) continue;
      posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);

      pairCount[channel(typeA, typeB)]++;
      pairCount[channel(typeB, typeA)]++;
      pairWeight[channel(typeA, typeB)] += weight;
      pairWeight[channel(typeB, typeA)] += weight;
      pairWeightPerBlock[block][channel(typeA, typeB)] += weight;
      pairWeightPerBlock[block][channel(typeB, typeA)] += weight;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      r = std::sqrt(rr);

      std::size_t bin = static_cast<std::size_t>(r / deltaR);
      if (bin < numberOfBins)
      {
        histogram(block, channel(typeA, typeB), bin) += weight;
        histogram(block, channel(typeB, typeA), bin) += weight;
      }
    }
  }

  for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    if (it1->isFractional) continue;
    posA = it1->position;
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    std::size_t compA = static_cast<std::size_t>(it1->componentId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);

    for (std::span<Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      if (it2->isFractional) continue;
      std::size_t molB = static_cast<std::size_t>(it2->moleculeId);
      std::size_t compB = static_cast<std::size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        std::size_t typeB = static_cast<std::size_t>(it2->type);

        pairCount[channel(typeA, typeB)]++;
        pairCount[channel(typeB, typeA)]++;
        pairWeight[channel(typeA, typeB)] += weight;
        pairWeight[channel(typeB, typeA)] += weight;
        pairWeightPerBlock[block][channel(typeA, typeB)] += weight;
        pairWeightPerBlock[block][channel(typeB, typeA)] += weight;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);
        r = std::sqrt(rr);

        std::size_t bin = static_cast<std::size_t>(r / deltaR);
        if (bin < numberOfBins)
        {
          histogram(block, channel(typeA, typeB), bin) += weight;
          histogram(block, channel(typeB, typeA), bin) += weight;
        }
      }
    }
  }

  histogram.addCount(block, 2.0 * weight);
}

std::pair<std::vector<double>, std::vector<double>> PropertyConventionalRadialDistributionFunction::normalized(
    std::size_t atomTypeA, std::size_t atomTypeB, double averageVolume) const
{
  const std::size_t ab = channel(atomTypeA, atomTypeB);
  const std::size_t ba = channel(atomTypeB, atomTypeA);

  // g(r) of a bin is the pairs counted in it over the pairs an ideal gas of the same density has there: the pairs
  // sampled, n_pairs, times the volume of the shell over the volume of the box.
  // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
  // The volume of the shell is (4 pi / 3) ((b + 1)^3 - b^3) dr^3 = 4 pi (b^2 + b + 1/3) dr^3. The mid-point value
  // 4 pi (b + 1/2)^2 dr^3 is smaller by dr^3 pi / 3: g(r) was too large by a third in the first bin and by 3.7% in
  // the second.
  const double prefactor = averageVolume / (2.0 * std::numbers::pi * deltaR * deltaR * deltaR);
  auto shell = [](std::size_t bin)
  {
    const double b = static_cast<double>(bin);
    return b * b + b + 1.0 / 3.0;
  };

  std::vector<double> mean(numberOfBins);
  const double pairs = pairWeight[ab] + pairWeight[ba];
  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    double counts = 0.0;
    for (std::size_t block = 0; block != histogram.numberOfBlocks; ++block) counts += histogram.bookKeeping[block][ab][bin];
    mean[bin] = pairs > 0.0 ? prefactor * counts / (pairs * shell(bin)) : 0.0;
  }

  // Each block is normalized by its own pairs. The blocks were all divided by the pairs of the whole run: where the
  // number of molecules changes (an open ensemble), a block with more of them than the average then lay above the
  // mean in every bin, and the spread of the number of pairs was reported as the error of g(r), for an ideal gas as
  // well.
  std::vector<std::vector<double>> blocks;
  for (std::size_t block = 0; block != histogram.numberOfBlocks; ++block)
  {
    if (!(histogram.numberOfCounts[block] > 0.0)) continue;
    const double blockPairs = pairWeightPerBlock[block][ab] + pairWeightPerBlock[block][ba];
    if (!(blockPairs > 0.0)) continue;
    std::vector<double> value(numberOfBins);
    for (std::size_t bin = 0; bin != numberOfBins; ++bin)
    {
      value[bin] = prefactor * histogram.bookKeeping[block][ab][bin] / (blockPairs * shell(bin));
    }
    blocks.push_back(std::move(value));
  }

  return {mean, blockErrorEstimate(blocks, mean)};
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
        // (samples of no weight: no volume, rather than 0 / 0)
        const double averageVolume = volumeSampleCount > 0.0 ? volumeCumulative / volumeSampleCount : 0.0;
        const auto [y, error_y] = normalized(atomTypeA, atomTypeB, averageVolume);
        std::vector<double> x(numberOfBins);
        for (std::size_t bin = 0; bin != numberOfBins; ++bin) x[bin] = (static_cast<double>(bin) + 0.5) * deltaR;

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

        // the average volume of the samples: the current volume differs from it under NPT
        const double averageVolume = volumeSampleCount > 0.0 ? volumeCumulative / volumeSampleCount : volume;
        const auto [y, error_y] = normalized(atomTypeA, atomTypeB, averageVolume);

        for (std::size_t bin = 0; bin != numberOfBins; ++bin)
        {
          stream_rdf_output << std::format("{} {} {}\n", (static_cast<double>(bin) + 0.5) * deltaR, y[bin],
                                           error_y[bin]);
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
  archive << rdf.pairWeight;
  archive << rdf.pairWeightPerBlock;

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
  archive >> rdf.pairWeight;
  archive >> rdf.pairWeightPerBlock;

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
