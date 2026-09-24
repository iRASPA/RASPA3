module;

module property_molecule_properties;

import std;

import archive;
import double3;
import atom;
import component;
import bond_potential;
import averages;

namespace
{

// Equilibrium length of a bond potential: the declared length for a Fixed bond, otherwise the
// grid minimum of the potential (analytic equilibria are type-specific; a fine grid is exact
// enough for a histogram range).
double bondEquilibriumLength(const BondPotential &bond)
{
  if (bond.type == BondType::Fixed) return bond.parameters[0];

  constexpr std::size_t numberOfGridPoints = 1024;
  constexpr double maximumLength = 5.0;
  double bestLength = 1.54;
  double bestEnergy = std::numeric_limits<double>::max();
  for (std::size_t i = 1; i <= numberOfGridPoints; ++i)
  {
    double r = maximumLength * static_cast<double>(i) / static_cast<double>(numberOfGridPoints);
    double energy = bond.calculateEnergy(double3(0.0, 0.0, 0.0), double3(r, 0.0, 0.0));
    if (energy < bestEnergy)
    {
      bestEnergy = energy;
      bestLength = r;
    }
  }
  return bestLength;
}

// Contour length between the end-to-end atoms: the sum of equilibrium bond lengths along the
// shortest topological path. Bonds without a potential (e.g. inside a rigid fragment) use the
// reference geometry, falling back to a generic 1.54 Angstrom when that is degenerate.
double contourLength(const Component &component, const std::array<std::size_t, 2> &ends)
{
  std::vector<std::size_t> path = component.connectivityTable.shortestPath(ends[0], ends[1]);
  double length = 0.0;
  for (std::size_t i = 1; i < path.size(); ++i)
  {
    std::size_t atomA = path[i - 1];
    std::size_t atomB = path[i];
    auto match = std::find_if(component.intraMolecularPotentials.bonds.begin(),
                              component.intraMolecularPotentials.bonds.end(),
                              [&](const BondPotential &bond)
                              {
                                return (bond.identifiers[0] == atomA && bond.identifiers[1] == atomB) ||
                                       (bond.identifiers[0] == atomB && bond.identifiers[1] == atomA);
                              });
    if (match != component.intraMolecularPotentials.bonds.end())
    {
      length += bondEquilibriumLength(*match);
    }
    else
    {
      double referenceLength = (component.atoms[atomA].position - component.atoms[atomB].position).length();
      length += referenceLength > 1e-6 ? referenceLength : 1.54;
    }
  }
  return length;
}

}  // namespace

PropertyMoleculeProperties::PropertyMoleculeProperties(std::size_t numberOfBlocks,
                                                       const std::vector<Component> &components,
                                                       std::size_t numberOfBins, double bondRange,
                                                       std::size_t sampleEvery, std::optional<std::size_t> writeEvery,
                                                       std::optional<double> endToEndRangeOverride)
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
      endToEndAtomsPerComponent(components.size()),
      endToEndRangePerComponent(components.size()),
      deltaEndToEndPerComponent(components.size()),
      bondHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      bendHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      torsionHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      endToEndHistogram(numberOfBlocks, std::vector<std::vector<std::vector<double>>>(components.size())),
      endToEndSum(numberOfBlocks, std::vector<double>(components.size())),
      endToEndSquaredSum(numberOfBlocks, std::vector<double>(components.size())),
      numberOfCounts(numberOfBlocks, std::vector<double>(components.size()))
{
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    numberOfBondsPerComponent[c] = components[c].intraMolecularPotentials.bonds.size();
    numberOfBendsPerComponent[c] = components[c].intraMolecularPotentials.bends.size();
    numberOfTorsionsPerComponent[c] = components[c].intraMolecularPotentials.torsions.size();

    endToEndAtomsPerComponent[c] = components[c].endToEndAtoms;
    if (endToEndAtomsPerComponent[c].has_value())
    {
      // Default range: the contour length between the ends (every reachable distance fits; bonds
      // stretching beyond their summed equilibria are astronomically rare and land in the last bin's
      // out-of-range discard, exactly like the other histograms).
      endToEndRangePerComponent[c] = endToEndRangeOverride.has_value()
                                         ? endToEndRangeOverride.value()
                                         : 1.05 * contourLength(components[c], endToEndAtomsPerComponent[c].value());
      deltaEndToEndPerComponent[c] = endToEndRangePerComponent[c] / static_cast<double>(numberOfBins);
    }

    for (std::size_t b = 0; b < numberOfBlocks; ++b)
    {
      bondHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfBondsPerComponent[c], std::vector<double>(numberOfBins));
      bendHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfBendsPerComponent[c], std::vector<double>(numberOfBins));
      torsionHistogram[b][c] =
          std::vector<std::vector<double>>(numberOfTorsionsPerComponent[c], std::vector<double>(numberOfBins));
      endToEndHistogram[b][c] = std::vector<std::vector<double>>(endToEndAtomsPerComponent[c].has_value() ? 1 : 0,
                                                                 std::vector<double>(numberOfBins));
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

      // End-to-end distance [Angstrom]. Positions are stored unwrapped (as the bond lengths above
      // rely on), so the plain difference is the physical distance.
      if (endToEndAtomsPerComponent[c].has_value())
      {
        const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();
        double r = (molecule[ends[0]].position - molecule[ends[1]].position).length();
        std::size_t bin = static_cast<std::size_t>(r / deltaEndToEndPerComponent[c]);
        if (bin < numberOfBins)
        {
          endToEndHistogram[block][c][0][bin] += 1.0;
        }
        endToEndSum[block][c] += r;
        endToEndSquaredSum[block][c] += r * r;
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
                              [](std::size_t n) { return n > 0; }) ||
                  std::any_of(endToEndAtomsPerComponent.begin(), endToEndAtomsPerComponent.end(),
                              [](const auto &pair) { return pair.has_value(); });
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

    if (endToEndAtomsPerComponent[c].has_value())
    {
      const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();

      // Block average and confidence interval of a per-molecule moment accumulator.
      auto momentStatistics = [&](const std::vector<std::vector<double>> &accumulator) -> std::pair<double, double>
      {
        double totalSamples{0.0}, totalSum{0.0};
        for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
        {
          totalSamples += numberOfCounts[blockIndex][c];
          totalSum += accumulator[blockIndex][c];
        }
        double mean = totalSamples > 0.0 ? totalSum / totalSamples : 0.0;

        std::size_t degreesOfFreedom = numberOfBlocks - 1;
        double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
        double sumOfSquares{0.0};
        std::size_t numberOfSamples{0};
        for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
        {
          if (numberOfCounts[blockIndex][c] > 0.0)
          {
            double value = accumulator[blockIndex][c] / numberOfCounts[blockIndex][c] - mean;
            sumOfSquares += value * value;
            ++numberOfSamples;
          }
        }
        double confidenceIntervalError{0.0};
        if (numberOfSamples >= 3)
        {
          double standardDeviation = std::sqrt(sumOfSquares / static_cast<double>(degreesOfFreedom));
          double standardError = standardDeviation / std::sqrt(static_cast<double>(numberOfBlocks));
          confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
        }
        return {mean, confidenceIntervalError};
      };
      auto [meanR, errorR] = momentStatistics(endToEndSum);
      auto [meanR2, errorR2] = momentStatistics(endToEndSquaredSum);

      std::ofstream stream(std::format("molecule_properties/end_to_end_{}_{}_{}.s{}.txt", components[c].name, ends[0],
                                       ends[1], systemId));
      stream << std::format("# end-to-end distance histogram, component: {}, atoms: ({}, {}), number of counts: {}\n",
                            components[c].name, ends[0], ends[1], totalNumberOfCounts);
      stream << std::format("# <R>    = {:g} +/- {:g} [Angstrom]\n", meanR, errorR);
      stream << std::format("# <R^2>  = {:g} +/- {:g} [Angstrom^2]\n", meanR2, errorR2);
      stream << std::format("# sqrt(<R^2>) = {:g} [Angstrom]\n", meanR2 > 0.0 ? std::sqrt(meanR2) : 0.0);
      stream << "# column 1: end-to-end distance [Angstrom]\n";
      stream << "# column 2: probability density [1/Angstrom]\n";
      stream << "# column 3: probability density error (95% confidence) [1/Angstrom]\n";

      auto [values, average, error] = result(endToEndHistogram, c, 0, deltaEndToEndPerComponent[c], 0.0);
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
  for (std::size_t c = 0; c < endToEndAtomsPerComponent.size(); ++c)
  {
    if (endToEndAtomsPerComponent[c].has_value())
    {
      std::print(stream, "    end-to-end component {}: atoms ({}, {}), range: 0 - {:g} [Angstrom]\n", c,
                 endToEndAtomsPerComponent[c]->at(0), endToEndAtomsPerComponent[c]->at(1),
                 endToEndRangePerComponent[c]);
    }
  }
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
  archive << p.endToEndAtomsPerComponent;
  archive << p.endToEndRangePerComponent;
  archive << p.deltaEndToEndPerComponent;
  archive << p.bondHistogram;
  archive << p.bendHistogram;
  archive << p.torsionHistogram;
  archive << p.endToEndHistogram;
  archive << p.endToEndSum;
  archive << p.endToEndSquaredSum;
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
  archive >> p.endToEndAtomsPerComponent;
  archive >> p.endToEndRangePerComponent;
  archive >> p.deltaEndToEndPerComponent;
  archive >> p.bondHistogram;
  archive >> p.bendHistogram;
  archive >> p.torsionHistogram;
  archive >> p.endToEndHistogram;
  archive >> p.endToEndSum;
  archive >> p.endToEndSquaredSum;
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
