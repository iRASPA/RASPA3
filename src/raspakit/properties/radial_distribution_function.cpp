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

void PropertyRadialDistributionFunction::sample(const SimulationBox &simulationBox,
                                                std::span<const Atom> frameworkAtoms,
                                                std::span<const AtomDynamics> frameworkDynamics,
                                                [[maybe_unused]] const std::vector<Molecule> &molecules,
                                                std::span<const Atom> moleculeAtoms,
                                                std::span<const AtomDynamics> moleculeDynamics,
                                                std::size_t currentCycle, std::size_t block)
{
  double3 dr, posA, posB;
  double3 gradientA, gradientB;
  double rr, r;

  if (currentCycle % sampleEvery != 0uz) return;

  if (moleculeAtoms.empty()) return;

  // Site forces must already be present in frameworkDynamics / moleculeDynamics (MD integrator, or
  // System::sampleForceBasedRDFWithFullGradients for MC). This sampler does not evaluate forces.

  auto accumulatePair = [&](std::size_t typeA, std::size_t typeB, const double3 &gradA, const double3 &gradB,
                            const double3 &separation)
  {
    pairCount[channel(typeA, typeB)]++;
    pairCount[channel(typeB, typeA)]++;

    dr = simulationBox.applyPeriodicBoundaryConditions(separation);
    rr = double3::dot(dr, dr);
    r = std::sqrt(rr);
    if (!(r < range) || r < 1.0e-12) return;

    // Borgis estimator: (∇U_A - ∇U_B)·(r_A - r_B) / r^3 accumulated for all bins with r_bin < r.
    const double value = double3::dot(gradA - gradB, dr) / (rr * r);
    for (std::size_t i = 0; i < numberOfBins; ++i)
    {
      if ((static_cast<double>(i) + 0.5) * deltaR < r)
      {
        histogram(block, channel(typeA, typeB), i) += value;
        histogram(block, channel(typeB, typeA), i) += value;
      }
    }
  };

  for (std::size_t ia = 0; ia < frameworkAtoms.size(); ++ia)
  {
    posA = frameworkAtoms[ia].position;
      gradientA = ia < frameworkDynamics.size() ? frameworkDynamics[ia].gradient : double3();
    std::size_t typeA = static_cast<std::size_t>(frameworkAtoms[ia].type);
    for (std::size_t ib = 0; ib < moleculeAtoms.size(); ++ib)
    {
      posB = moleculeAtoms[ib].position;
      gradientB = moleculeDynamics[ib].gradient;
      std::size_t typeB = static_cast<std::size_t>(moleculeAtoms[ib].type);
      accumulatePair(typeA, typeB, gradientA, gradientB, posA - posB);
    }
  }

  for (std::size_t ia = 0; ia + 1 < moleculeAtoms.size(); ++ia)
  {
    posA = moleculeAtoms[ia].position;
    gradientA = moleculeDynamics[ia].gradient;
    std::size_t molA = static_cast<std::size_t>(moleculeAtoms[ia].moleculeId);
    std::size_t compA = static_cast<std::size_t>(moleculeAtoms[ia].componentId);
    std::size_t typeA = static_cast<std::size_t>(moleculeAtoms[ia].type);

    for (std::size_t ib = ia + 1; ib < moleculeAtoms.size(); ++ib)
    {
      std::size_t molB = static_cast<std::size_t>(moleculeAtoms[ib].moleculeId);
      std::size_t compB = static_cast<std::size_t>(moleculeAtoms[ib].componentId);

      // skip interactions within the same molecule
      if ((compA == compB) && (molA == molB)) continue;

      posB = moleculeAtoms[ib].position;
      gradientB = moleculeDynamics[ib].gradient;
      std::size_t typeB = static_cast<std::size_t>(moleculeAtoms[ib].type);
      accumulatePair(typeA, typeB, gradientA, gradientB, posA - posB);
    }
  }

  histogram.addCount(block, 2.0);
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
        double normalization = 1.0 / std::abs(average[0]);

        for (std::size_t bin = 0; bin != numberOfBins; ++bin)
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
