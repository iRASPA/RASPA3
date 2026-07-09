module;

module property_molecule_properties;

import std;

import archive;
import double3;
import atom;
import component;
import averages;

PropertyMoleculeProperties::PropertyMoleculeProperties(std::size_t numberOfBlocks,
                                                       const std::vector<Component> &components,
                                                       std::size_t numberOfBins, double bondRange,
                                                       std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
    : numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      numberOfComponents(components.size()),
      bondRange(bondRange),
      deltaBond(bondRange / static_cast<double>(numberOfBins)),
      deltaBend(bendRange / static_cast<double>(numberOfBins)),
      deltaTorsion(torsionRange / static_cast<double>(numberOfBins)),
      sampleEvery(sampleEvery),
      writeEvery(writeEvery),
      numberOfBondsPerComponent(components.size()),
      numberOfBendsPerComponent(components.size()),
      numberOfTorsionsPerComponent(components.size()),
      bondHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      bendHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      torsionHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      numberOfCounts(numberOfBlocks, std::vector<double>(components.size()))
{
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    numberOfBondsPerComponent[c] = components[c].intraMolecularPotentials.bonds.size();
    numberOfBendsPerComponent[c] = components[c].intraMolecularPotentials.bends.size();
    numberOfTorsionsPerComponent[c] = components[c].intraMolecularPotentials.torsions.size();

    for (std::size_t b = 0; b < numberOfBlocks; ++b)
    {
      bondHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfBondsPerComponent[c], std::vector<double>(numberOfBins));
      bendHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfBendsPerComponent[c], std::vector<double>(numberOfBins));
      torsionHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfTorsionsPerComponent[c], std::vector<double>(numberOfBins));
    }
  }
}

void PropertyMoleculeProperties::sample(const std::vector<Component> &components,
                                        const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                        std::span<const Atom> moleculeAtoms, std::size_t currentCycle,
                                        std::size_t block)
{
  if (currentCycle % sampleEvery != 0uz) return;
  if (moleculeAtoms.empty()) return;

  std::size_t offset{0};
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    std::size_t numberOfAtoms = components[c].atoms.size();
    std::size_t numberOfMolecules = numberOfMoleculesPerComponent[c];

    const auto &bonds = components[c].intraMolecularPotentials.bonds;
    const auto &bends = components[c].intraMolecularPotentials.bends;
    const auto &torsions = components[c].intraMolecularPotentials.torsions;

    for (std::size_t m = 0; m < numberOfMolecules; ++m)
    {
      std::span<const Atom> molecule = moleculeAtoms.subspan(offset, numberOfAtoms);
      offset += numberOfAtoms;

      // Bond lengths [Angstrom].
      for (std::size_t i = 0; i < bonds.size(); ++i)
      {
        double3 dr = molecule[bonds[i].identifiers[0]].position - molecule[bonds[i].identifiers[1]].position;
        double r = dr.length();
        std::size_t bin = static_cast<std::size_t>(r / deltaBond);
        if (bin < numberOfBins)
        {
          bondHistogram[block][c][i][bin] += 1.0;
        }
      }

      // Bend angles [degrees].
      for (std::size_t i = 0; i < bends.size(); ++i)
      {
        double3 dr_ab = molecule[bends[i].identifiers[0]].position - molecule[bends[i].identifiers[1]].position;
        double3 dr_cb = molecule[bends[i].identifiers[2]].position - molecule[bends[i].identifiers[1]].position;
        double cos_theta = double3::dot(dr_ab, dr_cb) / (dr_ab.length() * dr_cb.length());
        cos_theta = std::clamp(cos_theta, -1.0, 1.0);
        double theta = std::acos(cos_theta) * (180.0 / std::numbers::pi);
        std::size_t bin = static_cast<std::size_t>(theta / deltaBend);
        if (bin < numberOfBins)
        {
          bendHistogram[block][c][i][bin] += 1.0;
        }
      }

      // Torsion (dihedral) angles [degrees], protein convention with sign in [-180, 180].
      for (std::size_t i = 0; i < torsions.size(); ++i)
      {
        double3 posA = molecule[torsions[i].identifiers[0]].position;
        double3 posB = molecule[torsions[i].identifiers[1]].position;
        double3 posC = molecule[torsions[i].identifiers[2]].position;
        double3 posD = molecule[torsions[i].identifiers[3]].position;

        double3 Dab = posA - posB;
        double3 Dcb = (posC - posB).normalized();
        double3 Ddc = posD - posC;

        double dot_ab = double3::dot(Dab, Dcb);
        double dot_dc = double3::dot(Ddc, Dcb);

        double3 dr = (Dab - dot_ab * Dcb).normalized();
        double3 ds = (Ddc - dot_dc * Dcb).normalized();

        double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
        double sign = double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc)));
        double phi = std::copysign(std::acos(cos_phi), sign) * (180.0 / std::numbers::pi);

        std::size_t bin = static_cast<std::size_t>((phi + 0.5 * torsionRange) / deltaTorsion);
        if (bin < numberOfBins)
        {
          torsionHistogram[block][c][i][bin] += 1.0;
        }
      }

      numberOfCounts[block][c] += 1.0;
    }
  }

  totalNumberOfCounts += 1.0;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> PropertyMoleculeProperties::result(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &histogram, std::size_t component,
    std::size_t index, double delta, double rangeStart) const
{
  std::vector<double> bins(numberOfBins);
  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    bins[bin] = rangeStart + (static_cast<double>(bin) + 0.5) * delta;
  }

  // Overall (block-combined) probability density.
  double totalSamples{0.0};
  std::vector<double> summedBlocks(numberOfBins);
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    totalSamples += numberOfCounts[blockIndex][component];
    for (std::size_t bin = 0; bin != numberOfBins; ++bin)
    {
      summedBlocks[bin] += histogram[blockIndex][component][index][bin];
    }
  }

  std::vector<double> average(numberOfBins);
  if (totalSamples > 0.0)
  {
    for (std::size_t bin = 0; bin != numberOfBins; ++bin)
    {
      average[bin] = summedBlocks[bin] / (totalSamples * delta);
    }
  }

  // Block averages for the confidence interval.
  std::size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];

  std::vector<double> sumOfSquares(numberOfBins);
  std::size_t numberOfSamples{0};
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    if (numberOfCounts[blockIndex][component] > 0.0)
    {
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        double blockAverage =
            histogram[blockIndex][component][index][bin] / (numberOfCounts[blockIndex][component] * delta);
        double value = blockAverage - average[bin];
        sumOfSquares[bin] += value * value;
      }
      ++numberOfSamples;
    }
  }

  std::vector<double> confidenceIntervalError(numberOfBins);
  if (numberOfSamples >= 3)
  {
    for (std::size_t bin = 0; bin != numberOfBins; ++bin)
    {
      double standardDeviation = std::sqrt(sumOfSquares[bin] / static_cast<double>(degreesOfFreedom));
      double standardError = standardDeviation / std::sqrt(static_cast<double>(numberOfBlocks));
      confidenceIntervalError[bin] = intermediateStandardNormalDeviate * standardError;
    }
  }

  return {bins, average, confidenceIntervalError};
}

void PropertyMoleculeProperties::writeOutput(std::size_t systemId, const std::vector<Component> &components,
                                             std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  bool anything = std::any_of(numberOfBondsPerComponent.begin(), numberOfBondsPerComponent.end(),
                              [](std::size_t n) { return n > 0; }) ||
                  std::any_of(numberOfBendsPerComponent.begin(), numberOfBendsPerComponent.end(),
                              [](std::size_t n) { return n > 0; }) ||
                  std::any_of(numberOfTorsionsPerComponent.begin(), numberOfTorsionsPerComponent.end(),
                              [](std::size_t n) { return n > 0; });
  if (!anything) return;

  std::filesystem::create_directory("molecule_properties");

  for (std::size_t c = 0; c < components.size() && c < numberOfComponents; ++c)
  {
    const auto &bonds = components[c].intraMolecularPotentials.bonds;
    for (std::size_t i = 0; i < bonds.size(); ++i)
    {
      std::ofstream stream(std::format("molecule_properties/bond_{}_{}_{}.s{}.txt", components[c].name,
                                       bonds[i].identifiers[0], bonds[i].identifiers[1], systemId));
      stream << std::format("# bond-length histogram, component: {}, number of counts: {}\n", components[c].name,
                            totalNumberOfCounts);
      stream << "# column 1: bond length [Angstrom]\n";
      stream << "# column 2: probability density [1/Angstrom]\n";
      stream << "# column 3: probability density error [1/Angstrom]\n";

      auto [values, average, error] = result(bondHistogram, c, i, deltaBond, 0.0);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }

    const auto &bends = components[c].intraMolecularPotentials.bends;
    for (std::size_t i = 0; i < bends.size(); ++i)
    {
      std::ofstream stream(std::format("molecule_properties/bend_{}_{}_{}_{}.s{}.txt", components[c].name,
                                       bends[i].identifiers[0], bends[i].identifiers[1], bends[i].identifiers[2],
                                       systemId));
      stream << std::format("# bend-angle histogram, component: {}, number of counts: {}\n", components[c].name,
                            totalNumberOfCounts);
      stream << "# column 1: bend angle [degrees]\n";
      stream << "# column 2: probability density [1/degrees]\n";
      stream << "# column 3: probability density error [1/degrees]\n";

      auto [values, average, error] = result(bendHistogram, c, i, deltaBend, 0.0);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }

    const auto &torsions = components[c].intraMolecularPotentials.torsions;
    for (std::size_t i = 0; i < torsions.size(); ++i)
    {
      std::ofstream stream(std::format("molecule_properties/torsion_{}_{}_{}_{}_{}.s{}.txt", components[c].name,
                                       torsions[i].identifiers[0], torsions[i].identifiers[1],
                                       torsions[i].identifiers[2], torsions[i].identifiers[3], systemId));
      stream << std::format("# torsion-angle histogram, component: {}, number of counts: {}\n", components[c].name,
                            totalNumberOfCounts);
      stream << "# column 1: torsion angle [degrees]\n";
      stream << "# column 2: probability density [1/degrees]\n";
      stream << "# column 3: probability density error [1/degrees]\n";

      auto [values, average, error] = result(torsionHistogram, c, i, deltaTorsion, -0.5 * torsionRange);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }
  }
}

std::string PropertyMoleculeProperties::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "Molecule-properties histograms:\n");
  std::print(stream, "    sample every: {}\n", sampleEvery);
  if (writeEvery.has_value())
  {
    std::print(stream, "    write every: {}\n", writeEvery.value());
  }
  std::print(stream, "    number of bins: {}\n", numberOfBins);
  std::print(stream, "    bond range: 0 - {} [Angstrom]\n", bondRange);
  std::print(stream, "    bend range: 0 - {} [degrees]\n", bendRange);
  std::print(stream, "    torsion range: {} - {} [degrees]\n", -0.5 * torsionRange, 0.5 * torsionRange);
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyMoleculeProperties &p)
{
  archive << p.versionNumber;

  archive << p.numberOfBlocks;
  archive << p.numberOfBins;
  archive << p.numberOfComponents;
  archive << p.bondRange;
  archive << p.bendRange;
  archive << p.torsionRange;
  archive << p.deltaBond;
  archive << p.deltaBend;
  archive << p.deltaTorsion;
  archive << p.sampleEvery;
  archive << p.writeEvery;
  archive << p.numberOfBondsPerComponent;
  archive << p.numberOfBendsPerComponent;
  archive << p.numberOfTorsionsPerComponent;
  archive << p.bondHistogram;
  archive << p.bendHistogram;
  archive << p.torsionHistogram;
  archive << p.numberOfCounts;
  archive << p.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMoleculeProperties &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyMoleculeProperties' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.numberOfBlocks;
  archive >> p.numberOfBins;
  archive >> p.numberOfComponents;
  archive >> p.bondRange;
  archive >> p.bendRange;
  archive >> p.torsionRange;
  archive >> p.deltaBond;
  archive >> p.deltaBend;
  archive >> p.deltaTorsion;
  archive >> p.sampleEvery;
  archive >> p.writeEvery;
  archive >> p.numberOfBondsPerComponent;
  archive >> p.numberOfBendsPerComponent;
  archive >> p.numberOfTorsionsPerComponent;
  archive >> p.bondHistogram;
  archive >> p.bendHistogram;
  archive >> p.torsionHistogram;
  archive >> p.numberOfCounts;
  archive >> p.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyMoleculeProperties: Error in binary restart\n"));
  }
#endif

  return archive;
}
