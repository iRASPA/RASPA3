module;

module property_rdf;

import std;

import archive;
import double3;
import atom;
import atom_dynamics;
import simulationbox;
import forcefield;
import averages;

void PropertyRadialDistributionFunction::sample(const SimulationBox& simulationBox, double beta,
                                                std::span<const Atom> frameworkAtoms,
                                                std::span<const double3> frameworkGradients,
                                                std::span<const std::uint8_t> frameworkMobile,
                                                std::span<const Atom> moleculeAtoms,
                                                std::span<const double3> moleculeGradients, std::size_t currentCycle,
                                                std::size_t block)
{
  if (currentCycle % sampleEvery != 0uz) return;

  if (moleculeAtoms.empty()) return;

  // Site forces must already be present in the gradients (MD integrator, or
  // System::sampleForceBasedRDFWithFullGradients for MC). This sampler does not evaluate forces.
  //
  // Borgis et al., Mol. Phys. 111, 3486 (2013): integrating the Yvon-Born-Green relation from g(0) = 0,
  //
  //   g_AB(r) = beta V / (8 pi P_AB) < sum_{i in A, j in B} H(r - r_ij) (F_i - F_j) . r_ij / r_ij^3 >
  //
  // over the ordered pairs, with P_AB their number (the normalisation of the conventional RDF). A fixed site i does
  // not move: its force is replaced by -F_j, the one-body relation of the site j around it.
  const std::size_t numberOfChannels = numberOfPseudoAtoms * numberOfPseudoAtoms;
  std::vector<double> sums(numberOfChannels * numberOfBins);
  std::vector<std::size_t> pairs(numberOfChannels);

  auto accumulatePair =
      [&](std::size_t typeA, std::size_t typeB, const double3& gradA, const double3& gradB, const double3& separation)
  {
    pairs[channel(typeA, typeB)]++;
    pairs[channel(typeB, typeA)]++;

    const double3 dr = simulationBox.applyPeriodicBoundaryConditions(separation);
    const double rr = double3::dot(dr, dr);
    const double r = std::sqrt(rr);
    if (!(r < range) || r < 1.0e-12) return;

    // (F_A - F_B) . r_AB / r^3 with F = -gradient, added to every bin beyond the pair distance
    const double value = double3::dot(gradB - gradA, dr) / (rr * r);
    const std::size_t first = static_cast<std::size_t>(std::max(0.0, std::ceil(r / deltaR - 0.5)));
    for (std::size_t i = first; i < numberOfBins; ++i)
    {
      sums[channel(typeA, typeB) * numberOfBins + i] += value;
      sums[channel(typeB, typeA) * numberOfBins + i] += value;
    }
  };

  for (std::size_t ia = 0; ia < frameworkAtoms.size(); ++ia)
  {
    const double3 posA = frameworkAtoms[ia].position;
    const bool mobile = ia < frameworkMobile.size() && frameworkMobile[ia] != 0;
    const std::size_t typeA = static_cast<std::size_t>(frameworkAtoms[ia].type);
    for (std::size_t ib = 0; ib < moleculeAtoms.size(); ++ib)
    {
      // the distribution of the integer molecules: a fractional molecule is no member of the ensemble
      if (moleculeAtoms[ib].isFractional) continue;
      const double3 gradientB = moleculeGradients[ib];
      const double3 gradientA = mobile ? frameworkGradients[ia] : -gradientB;
      const std::size_t typeB = static_cast<std::size_t>(moleculeAtoms[ib].type);
      accumulatePair(typeA, typeB, gradientA, gradientB, posA - moleculeAtoms[ib].position);
    }
  }

  for (std::size_t ia = 0; ia + 1 < moleculeAtoms.size(); ++ia)
  {
    if (moleculeAtoms[ia].isFractional) continue;
    const double3 posA = moleculeAtoms[ia].position;
    const std::size_t molA = static_cast<std::size_t>(moleculeAtoms[ia].moleculeId);
    const std::size_t compA = static_cast<std::size_t>(moleculeAtoms[ia].componentId);
    const std::size_t typeA = static_cast<std::size_t>(moleculeAtoms[ia].type);

    for (std::size_t ib = ia + 1; ib < moleculeAtoms.size(); ++ib)
    {
      if (moleculeAtoms[ib].isFractional) continue;
      const std::size_t molB = static_cast<std::size_t>(moleculeAtoms[ib].moleculeId);
      const std::size_t compB = static_cast<std::size_t>(moleculeAtoms[ib].componentId);

      // skip interactions within the same molecule
      if ((compA == compB) && (molA == molB)) continue;

      const std::size_t typeB = static_cast<std::size_t>(moleculeAtoms[ib].type);
      accumulatePair(typeA, typeB, moleculeGradients[ia], moleculeGradients[ib], posA - moleculeAtoms[ib].position);
    }
  }

  for (std::size_t c = 0; c != numberOfChannels; ++c)
  {
    pairCount[c] += pairs[c];
    if (pairs[c] == 0uz) continue;
    const double normalization = beta * simulationBox.volume / (8.0 * std::numbers::pi * static_cast<double>(pairs[c]));
    for (std::size_t i = 0; i != numberOfBins; ++i)
      histogram(block, c, i) += normalization * sums[c * numberOfBins + i];
  }
  histogram.addCount(block, 1.0);
}

void PropertyRadialDistributionFunction::writeOutput(const ForceField &forceField, std::size_t systemId,
                                                     [[maybe_unused]] double volume,
                                                     [[maybe_unused]] std::vector<std::size_t> &numberOfPseudoAtomsType,
                                                     std::size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("rdf");

  for (std::size_t atomTypeA = 0; atomTypeA < numberOfPseudoAtoms; ++atomTypeA)
  {
    for (std::size_t atomTypeB = atomTypeA; atomTypeB < numberOfPseudoAtoms; ++atomTypeB)
    {
      if (pairCount[channel(atomTypeB, atomTypeA)] > 0)
      {
        std::ofstream stream_rdf_output(std::format("rdf/rdf_{}_{}.s{}.txt", forceField.pseudoAtoms[atomTypeA].name,
                                                    forceField.pseudoAtoms[atomTypeB].name, systemId));

        stream_rdf_output << std::format("# rdf, number of counts: {}\n", histogram.totalNumberOfCounts);
        stream_rdf_output << "# column 1: distance []\n";
        stream_rdf_output << "# column 2: normalize rdf []\n";
        stream_rdf_output << "# column 3: error normalize rdf []\n";

        auto [average, error] = result(atomTypeA, atomTypeB);

        // n_pairs is the number of unique pairs of atoms where one atom is from each of two sets
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3085256/
        [[maybe_unused]] double avg_n_pairs =
            static_cast<double>(pairCount[channel(atomTypeA, atomTypeB)] + pairCount[channel(atomTypeB, atomTypeA)]) /
            histogram.totalNumberOfCounts;
        for (std::size_t bin = 0; bin != numberOfBins; ++bin)
        {
          stream_rdf_output << std::format("{} {} {}\n", (static_cast<double>(bin) + 0.5) * deltaR, average[bin],
                                           error[bin]);
        }
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyRadialDistributionFunction &rdf)
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

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadialDistributionFunction &rdf)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > rdf.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyRadialDistributionFunction' at line {} in file {}\n",
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

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyRadialDistributionFunction: Error in binary restart\n"));
  }
#endif

  return archive;
}
