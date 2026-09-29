module;

module property_polymer_backbone;

import std;

import archive;
import double3;
import atom;
import component;
import averages;
import property_molecule_properties;

namespace
{

// Weighted least-squares slope of y against x (unit weights when 'w' is empty).
std::optional<double> slope(std::span<const double> x, std::span<const double> y, std::span<const double> w = {})
{
  if (x.size() < 2) return std::nullopt;
  double sumW{0.0}, meanX{0.0}, meanY{0.0};
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    double wi = w.empty() ? 1.0 : w[i];
    sumW += wi;
    meanX += wi * x[i];
    meanY += wi * y[i];
  }
  if (sumW <= 0.0) return std::nullopt;
  meanX /= sumW;
  meanY /= sumW;
  double sxx{0.0}, sxy{0.0};
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    double wi = w.empty() ? 1.0 : w[i];
    sxx += wi * (x[i] - meanX) * (x[i] - meanX);
    sxy += wi * (x[i] - meanX) * (y[i] - meanY);
  }
  if (sxx <= 0.0) return std::nullopt;
  return sxy / sxx;
}

}  // namespace

PropertyPolymerBackbone::PropertyPolymerBackbone(std::size_t numberOfBlocks, const std::vector<Component> &components,
                                                 std::size_t numberOfWaveVectors, double waveVectorLowerLimit,
                                                 double waveVectorUpperLimit, std::size_t sampleEvery,
                                                 std::optional<std::size_t> writeEvery)
    : numberOfBlocks(numberOfBlocks),
      numberOfComponents(components.size()),
      sampleEvery(sampleEvery),
      writeEvery(writeEvery),
      numberOfWaveVectors(numberOfWaveVectors),
      waveVectorLowerLimit(waveVectorLowerLimit),
      waveVectorUpperLimit(waveVectorUpperLimit),
      waveVectors(numberOfWaveVectors),
      backbonePerComponent(components.size()),
      contourLengthPerComponent(components.size()),
      internalDistanceSquaredSum(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      bondCorrelationSum(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      formFactorSum(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      sums(numberOfBlocks, std::vector<Moments>(components.size(), Moments{})),
      numberOfCounts(numberOfBlocks, std::vector<double>(components.size()))
{
  // Logarithmic grid; a single point sits at the lower limit.
  for (std::size_t i = 0; i < numberOfWaveVectors; ++i)
  {
    double fraction = numberOfWaveVectors > 1 ? static_cast<double>(i) / static_cast<double>(numberOfWaveVectors - 1)
                                              : 0.0;
    waveVectors[i] = waveVectorLowerLimit * std::pow(waveVectorUpperLimit / waveVectorLowerLimit, fraction);
  }

  for (std::size_t c = 0; c < components.size(); ++c)
  {
    std::vector<std::size_t> backbone = components[c].backboneAtoms();
    // Two backbone beads give a single bond: no internal structure to speak of.
    if (backbone.size() < 3) continue;

    backbonePerComponent[c] = backbone;
    contourLengthPerComponent[c] = contourLength(components[c], {backbone.front(), backbone.back()});

    std::size_t numberOfBeads = backbone.size();
    for (std::size_t b = 0; b < numberOfBlocks; ++b)
    {
      internalDistanceSquaredSum[b][c] = std::vector<double>(numberOfBeads - 1);
      bondCorrelationSum[b][c] = std::vector<double>(numberOfBeads - 1);
      formFactorSum[b][c] = std::vector<double>(numberOfWaveVectors);
    }
  }
}

void PropertyPolymerBackbone::accumulateInternalDistances(std::span<const Atom> molecule,
                                                          std::span<const std::size_t> backbone,
                                                          std::span<double> internalDistanceSquared)
{
  std::size_t numberOfBeads = backbone.size();
  for (std::size_t k = 1; k < numberOfBeads; ++k)
  {
    double sum = 0.0;
    for (std::size_t i = 0; i + k < numberOfBeads; ++i)
    {
      double3 dr = molecule[backbone[i + k]].position - molecule[backbone[i]].position;
      sum += double3::dot(dr, dr);
    }
    internalDistanceSquared[k - 1] += sum / static_cast<double>(numberOfBeads - k);
  }
}

void PropertyPolymerBackbone::accumulateBondCorrelation(std::span<const Atom> molecule,
                                                        std::span<const std::size_t> backbone,
                                                        std::span<double> bondCorrelation)
{
  std::size_t numberOfBonds = backbone.size() - 1;
  std::vector<double3> unitBonds(numberOfBonds);
  for (std::size_t i = 0; i < numberOfBonds; ++i)
  {
    unitBonds[i] = (molecule[backbone[i + 1]].position - molecule[backbone[i]].position).normalized();
  }
  for (std::size_t k = 0; k < numberOfBonds; ++k)
  {
    double sum = 0.0;
    for (std::size_t i = 0; i + k < numberOfBonds; ++i)
    {
      sum += double3::dot(unitBonds[i], unitBonds[i + k]);
    }
    bondCorrelation[k] += sum / static_cast<double>(numberOfBonds - k);
  }
}

void PropertyPolymerBackbone::accumulateFormFactor(std::span<const Atom> molecule, std::span<const double> waveVectors,
                                                   std::span<double> formFactor)
{
  std::size_t n = molecule.size();
  double inverseN2 = 1.0 / (static_cast<double>(n) * static_cast<double>(n));

  // Pair distances once; the self terms contribute N.
  std::vector<double> distances;
  distances.reserve(n * (n - 1) / 2);
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    for (std::size_t j = i + 1; j < n; ++j)
    {
      distances.push_back((molecule[i].position - molecule[j].position).length());
    }
  }

  for (std::size_t iq = 0; iq < waveVectors.size(); ++iq)
  {
    double q = waveVectors[iq];
    double sum = 0.0;
    for (double r : distances)
    {
      double qr = q * r;
      sum += qr > 1e-12 ? std::sin(qr) / qr : 1.0;
    }
    formFactor[iq] += (static_cast<double>(n) + 2.0 * sum) * inverseN2;
  }
}

PropertyPolymerBackbone::Moments PropertyPolymerBackbone::computeMoments(std::span<const Atom> molecule,
                                                                         std::span<const std::size_t> backbone)
{
  Moments moments{};
  std::size_t numberOfBonds = backbone.size() - 1;

  std::vector<double3> bonds(numberOfBonds);
  double sumLength{0.0}, sumLengthSquared{0.0};
  for (std::size_t i = 0; i < numberOfBonds; ++i)
  {
    bonds[i] = molecule[backbone[i + 1]].position - molecule[backbone[i]].position;
    double l2 = double3::dot(bonds[i], bonds[i]);
    sumLength += std::sqrt(l2);
    sumLengthSquared += l2;
  }
  moments[BondLength] = sumLength / static_cast<double>(numberOfBonds);
  moments[BondLengthSquared] = sumLengthSquared / static_cast<double>(numberOfBonds);

  double3 endToEnd = molecule[backbone.back()].position - molecule[backbone.front()].position;
  moments[EndToEndSquared] = double3::dot(endToEnd, endToEnd);

  // Projection of the chain onto its first bond and (by symmetry of the reversed chain) onto its
  // last bond: both equal <sum_j b^_end . b_j> for the respective end.
  double3 firstUnit = bonds.front().normalized();
  double3 lastUnit = bonds.back().normalized();
  moments[Projection] = 0.5 * (double3::dot(firstUnit, endToEnd) + double3::dot(lastUnit, endToEnd));

  double3 center{};
  for (const Atom &atom : molecule) center += atom.position;
  center = center / static_cast<double>(molecule.size());
  double rg2 = 0.0;
  for (const Atom &atom : molecule)
  {
    double3 d = atom.position - center;
    rg2 += double3::dot(d, d);
  }
  moments[RadiusOfGyrationSquared] = rg2 / static_cast<double>(molecule.size());

  return moments;
}

void PropertyPolymerBackbone::sample(const std::vector<Component> &components,
                                     const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                     std::span<const Atom> moleculeAtoms, std::size_t currentCycle, std::size_t block)
{
  if (currentCycle % sampleEvery != 0uz) return;
  if (moleculeAtoms.empty()) return;

  std::size_t offset{0};
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    std::size_t numberOfAtoms = components[c].atoms.size();
    std::size_t numberOfMolecules = numberOfMoleculesPerComponent[c];

    if (!isSampled(c))
    {
      offset += numberOfAtoms * numberOfMolecules;
      continue;
    }

    const std::vector<std::size_t> &backbone = backbonePerComponent[c];

    for (std::size_t m = 0; m < numberOfMolecules; ++m)
    {
      // Positions are stored unwrapped, so plain differences are the physical intramolecular vectors.
      std::span<const Atom> molecule = moleculeAtoms.subspan(offset, numberOfAtoms);
      offset += numberOfAtoms;

      accumulateInternalDistances(molecule, backbone, internalDistanceSquaredSum[block][c]);
      accumulateBondCorrelation(molecule, backbone, bondCorrelationSum[block][c]);
      accumulateFormFactor(molecule, waveVectors, formFactorSum[block][c]);

      Moments moments = computeMoments(molecule, backbone);
      for (std::size_t i = 0; i < NumberOfMoments; ++i) sums[block][c][i] += moments[i];

      numberOfCounts[block][c] += 1.0;
    }
  }

  totalNumberOfCounts += 1.0;
}

PropertyPolymerBackbone::Averages PropertyPolymerBackbone::blockAverages(std::size_t block,
                                                                         std::size_t component) const
{
  Averages averages{};
  averages.internalDistanceSquared = internalDistanceSquaredSum[block][component];
  averages.bondCorrelation = bondCorrelationSum[block][component];
  averages.formFactor = formFactorSum[block][component];
  averages.moments = sums[block][component];

  double count = numberOfCounts[block][component];
  if (count > 0.0)
  {
    for (double &v : averages.internalDistanceSquared) v /= count;
    for (double &v : averages.bondCorrelation) v /= count;
    for (double &v : averages.formFactor) v /= count;
    for (double &v : averages.moments) v /= count;
  }
  return averages;
}

PropertyPolymerBackbone::Averages PropertyPolymerBackbone::overallAverages(std::size_t component) const
{
  Averages averages{};
  averages.internalDistanceSquared = std::vector<double>(internalDistanceSquaredSum[0][component].size());
  averages.bondCorrelation = std::vector<double>(bondCorrelationSum[0][component].size());
  averages.formFactor = std::vector<double>(formFactorSum[0][component].size());

  double total{0.0};
  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    total += numberOfCounts[block][component];
    for (std::size_t i = 0; i < averages.internalDistanceSquared.size(); ++i)
      averages.internalDistanceSquared[i] += internalDistanceSquaredSum[block][component][i];
    for (std::size_t i = 0; i < averages.bondCorrelation.size(); ++i)
      averages.bondCorrelation[i] += bondCorrelationSum[block][component][i];
    for (std::size_t i = 0; i < averages.formFactor.size(); ++i)
      averages.formFactor[i] += formFactorSum[block][component][i];
    for (std::size_t i = 0; i < NumberOfMoments; ++i) averages.moments[i] += sums[block][component][i];
  }
  if (total > 0.0)
  {
    for (double &v : averages.internalDistanceSquared) v /= total;
    for (double &v : averages.bondCorrelation) v /= total;
    for (double &v : averages.formFactor) v /= total;
    for (double &v : averages.moments) v /= total;
  }
  return averages;
}

std::pair<double, double> PropertyPolymerBackbone::statistics(
    std::size_t component, const std::function<double(const Averages &)> &function) const
{
  double totalSamples{0.0};
  for (std::size_t block = 0; block < numberOfBlocks; ++block) totalSamples += numberOfCounts[block][component];
  if (totalSamples <= 0.0) return {0.0, 0.0};

  double mean = function(overallAverages(component));

  std::size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  double sumOfSquares{0.0};
  std::size_t numberOfSamples{0};
  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    if (numberOfCounts[block][component] > 0.0)
    {
      double value = function(blockAverages(block, component)) - mean;
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
}

double PropertyPolymerBackbone::persistenceLengthFromProjection(const Averages &averages)
{
  return averages.moments[Projection];
}

// Least-squares fit of ln C(k) = a - k <l> / l_p over the initial decay: consecutive points from
// k = 0 up to the first C(k) <= 0.05 (the noise floor), weighted by C(k)^2 so that the fit is that
// of C(k) itself and the log of a small anticorrelation dip (a chain whose backbone alternates bond
// angles has C(k) oscillating with the repeat period) does not dominate. Requires at least three
// points and a decaying slope; otherwise zero is returned.
double PropertyPolymerBackbone::persistenceLengthFromFit(const Averages &averages)
{
  std::vector<double> x{}, y{}, w{};
  for (std::size_t k = 0; k < averages.bondCorrelation.size(); ++k)
  {
    double c = averages.bondCorrelation[k];
    if (c <= 0.05) break;
    x.push_back(static_cast<double>(k));
    y.push_back(std::log(c));
    w.push_back(c * c);
  }
  if (x.size() < 3) return 0.0;
  std::optional<double> s = slope(x, y, w);
  if (!s.has_value() || s.value() >= 0.0) return 0.0;
  return -averages.moments[BondLength] / s.value();
}

// Log-log slope of <r^2(k)> against k over Nb/8 <= k <= Nb/2 divided by two; the window skips the
// stiff short-range regime and the finite-size roll-off near the full chain. Requires at least four
// points (Nb >= 16); otherwise zero is returned.
double PropertyPolymerBackbone::floryExponent(const Averages &averages)
{
  std::size_t numberOfBeads = averages.internalDistanceSquared.size() + 1;
  std::size_t kMin = std::max<std::size_t>(2, numberOfBeads / 8);
  std::size_t kMax = numberOfBeads / 2;
  if (kMax < kMin + 3) return 0.0;

  std::vector<double> x{}, y{};
  for (std::size_t k = kMin; k <= kMax; ++k)
  {
    double r2 = averages.internalDistanceSquared[k - 1];
    if (r2 <= 0.0) continue;
    x.push_back(std::log(static_cast<double>(k)));
    y.push_back(std::log(r2));
  }
  if (x.size() < 4) return 0.0;
  std::optional<double> s = slope(x, y);
  return s.has_value() ? 0.5 * s.value() : 0.0;
}

double PropertyPolymerBackbone::debyeFunction(double x)
{
  if (x < 1e-6) return 1.0 - x / 3.0;
  return 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
}

void PropertyPolymerBackbone::writeOutput(std::size_t systemId, const std::vector<Component> &components,
                                          std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  bool anything = false;
  for (std::size_t c = 0; c < numberOfComponents; ++c) anything = anything || isSampled(c);
  if (!anything) return;

  std::filesystem::create_directory("polymer_backbone");

  for (std::size_t c = 0; c < components.size() && c < numberOfComponents; ++c)
  {
    if (!isSampled(c)) continue;
    const std::string &name = components[c].name;
    const std::vector<std::size_t> &backbone = backbonePerComponent[c];
    std::size_t numberOfBeads = backbone.size();
    std::size_t numberOfBonds = numberOfBeads - 1;
    double rMax = contourLengthPerComponent[c];

    auto [meanL, errorL] = statistics(c, [](const Averages &a) { return a.moments[BondLength]; });
    auto [meanL2, errorL2] = statistics(c, [](const Averages &a) { return a.moments[BondLengthSquared]; });
    auto [meanR2, errorR2] = statistics(c, [](const Averages &a) { return a.moments[EndToEndSquared]; });
    auto [meanRg2, errorRg2] = statistics(c, [](const Averages &a) { return a.moments[RadiusOfGyrationSquared]; });
    auto [ratioR2Rg2, errorRatioR2Rg2] = statistics(
        c, [](const Averages &a)
        {
          return a.moments[RadiusOfGyrationSquared] > 0.0 ? a.moments[EndToEndSquared] / a.moments[RadiusOfGyrationSquared]
                                                          : 0.0;
        });
    auto [characteristicRatio, errorCharacteristicRatio] = statistics(
        c, [numberOfBonds](const Averages &a)
        {
          double l = a.moments[BondLength];
          return l > 0.0 ? a.moments[EndToEndSquared] / (static_cast<double>(numberOfBonds) * l * l) : 0.0;
        });
    auto [kuhnLength, errorKuhnLength] =
        statistics(c, [rMax](const Averages &a) { return rMax > 0.0 ? a.moments[EndToEndSquared] / rMax : 0.0; });
    auto [kuhnSegments, errorKuhnSegments] = statistics(
        c, [rMax](const Averages &a)
        { return a.moments[EndToEndSquared] > 0.0 ? rMax * rMax / a.moments[EndToEndSquared] : 0.0; });
    auto [persistenceProjection, errorPersistenceProjection] = statistics(c, persistenceLengthFromProjection);
    auto [persistenceFit, errorPersistenceFit] = statistics(c, persistenceLengthFromFit);
    auto [nu, errorNu] = statistics(c, floryExponent);

    std::size_t kMin = std::max<std::size_t>(2, numberOfBeads / 8);
    std::size_t kMax = numberOfBeads / 2;

    {
      std::ofstream summary(std::format("polymer_backbone/polymer_backbone_{}.s{}.txt", name, systemId));
      summary << std::format("# backbone chain statistics, component: {}, number of counts: {}\n", name,
                             totalNumberOfCounts);
      summary << std::format("# backbone: {} beads, {} bonds, atoms {} .. {}; errors are 95% confidence intervals\n",
                             numberOfBeads, numberOfBonds, backbone.front(), backbone.back());
      summary << "#\n";
      summary << std::format("R_max contour length            {:<14.6g}                   [Angstrom]\n", rMax);
      summary << std::format("<l>   bond length                {:<14.6g} +/- {:<12.6g} [Angstrom]\n", meanL, errorL);
      summary << std::format("<l^2>                           {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanL2, errorL2);
      summary << std::format("<R^2> end-to-end                {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanR2, errorR2);
      summary << std::format("<Rg^2> (all atoms)              {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanRg2,
                             errorRg2);
      summary << std::format("<R^2>/<Rg^2>                    {:<14.6g} +/- {:<12.6g} [-]   (ideal chain: 6)\n",
                             ratioR2Rg2, errorRatioR2Rg2);
      summary << "#\n";
      summary << std::format("C_N = <R^2>/((Nb-1)<l>^2)       {:<14.6g} +/- {:<12.6g} [-]   (freely jointed: 1)\n",
                             characteristicRatio, errorCharacteristicRatio);
      summary << std::format("b_K = <R^2>/R_max Kuhn length   {:<14.6g} +/- {:<12.6g} [Angstrom]\n", kuhnLength,
                             errorKuhnLength);
      summary << std::format("N_K = R_max/b_K Kuhn segments   {:<14.6g} +/- {:<12.6g} [-]\n", kuhnSegments,
                             errorKuhnSegments);
      summary << std::format("l_p  persistence (projection)   {:<14.6g} +/- {:<12.6g} [Angstrom]   (<sum_j b^_end . b_j>)\n",
                             persistenceProjection, errorPersistenceProjection);
      summary << std::format("l_p  persistence (fit)          {:<14.6g} +/- {:<12.6g} [Angstrom]   (ln C(k) = -k<l>/l_p, C(k) > 0.05, weights C(k)^2; 0 if undetermined)\n",
                             persistenceFit, errorPersistenceFit);
      summary << std::format("nu   Flory exponent             {:<14.6g} +/- {:<12.6g} [-]   (<r^2(k)> ~ k^2nu, {} <= k <= {}; ideal 0.5, SAW 0.588; 0 if undetermined)\n",
                             nu, errorNu, kMin, kMax);
    }

    {
      std::ofstream stream(std::format("polymer_backbone/internal_distances_{}.s{}.txt", name, systemId));
      stream << std::format("# mean squared internal distances along the backbone, component: {}, number of counts: {}\n",
                            name, totalNumberOfCounts);
      stream << std::format("# <l> = {:g} [Angstrom]\n", meanL);
      stream << "# column 1: separation k [bonds]\n";
      stream << "# column 2: contour separation k <l> [Angstrom]\n";
      stream << "# column 3: <r^2(k)> [Angstrom^2]\n";
      stream << "# column 4: <r^2(k)> error (95% confidence) [Angstrom^2]\n";
      stream << "# column 5: <r^2(k)> / (k <l>^2) [-]  (flat for an ideal chain)\n";
      for (std::size_t k = 1; k < numberOfBeads; ++k)
      {
        auto [value, error] = statistics(c, [k](const Averages &a) { return a.internalDistanceSquared[k - 1]; });
        double reduced = meanL > 0.0 ? value / (static_cast<double>(k) * meanL * meanL) : 0.0;
        stream << std::format("{} {} {} {} {}\n", k, static_cast<double>(k) * meanL, value, error, reduced);
      }
    }

    {
      std::ofstream stream(std::format("polymer_backbone/bond_correlation_{}.s{}.txt", name, systemId));
      stream << std::format("# bond-vector correlation along the backbone, component: {}, number of counts: {}\n",
                            name, totalNumberOfCounts);
      stream << std::format("# l_p (projection) = {:g} [Angstrom], l_p (fit) = {:g} [Angstrom]\n",
                            persistenceProjection, persistenceFit);
      stream << "# column 1: separation k [bonds]\n";
      stream << "# column 2: C(k) = <b^_i . b^_(i+k)> [-]\n";
      stream << "# column 3: C(k) error (95% confidence) [-]\n";
      for (std::size_t k = 0; k < numberOfBonds; ++k)
      {
        auto [value, error] = statistics(c, [k](const Averages &a) { return a.bondCorrelation[k]; });
        stream << std::format("{} {} {}\n", k, value, error);
      }
    }

    {
      std::ofstream stream(std::format("polymer_backbone/form_factor_{}.s{}.txt", name, systemId));
      stream << std::format("# single-chain form factor (all atoms, uniform weights), component: {}, number of counts: {}\n",
                            name, totalNumberOfCounts);
      stream << std::format("# <Rg^2> = {:g} [Angstrom^2]; Debye function evaluated at x = q^2 <Rg^2>\n", meanRg2);
      stream << "# column 1: q [1/Angstrom]\n";
      stream << "# column 2: q sqrt(<Rg^2>) [-]\n";
      stream << "# column 3: P(q) [-]\n";
      stream << "# column 4: P(q) error (95% confidence) [-]\n";
      stream << "# column 5: Debye P(q) for a Gaussian chain with the same <Rg^2> [-]\n";
      double rg = std::sqrt(std::max(meanRg2, 0.0));
      for (std::size_t iq = 0; iq < numberOfWaveVectors; ++iq)
      {
        auto [value, error] = statistics(c, [iq](const Averages &a) { return a.formFactor[iq]; });
        double q = waveVectors[iq];
        stream << std::format("{} {} {} {} {}\n", q, q * rg, value, error, debyeFunction(q * q * meanRg2));
      }
    }
  }
}

std::string PropertyPolymerBackbone::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "Polymer-backbone chain statistics:\n");
  std::print(stream, "    sample every: {}\n", sampleEvery);
  if (writeEvery.has_value())
  {
    std::print(stream, "    write every: {}\n", writeEvery.value());
  }
  std::print(stream, "    wave vectors: {} log-spaced in {:g} - {:g} [1/Angstrom]\n", numberOfWaveVectors,
             waveVectorLowerLimit, waveVectorUpperLimit);
  for (std::size_t c = 0; c < backbonePerComponent.size(); ++c)
  {
    if (isSampled(c))
    {
      std::print(stream, "    component {}: backbone of {} beads (atoms {} .. {}), contour length {:g} [Angstrom]\n", c,
                 backbonePerComponent[c].size(), backbonePerComponent[c].front(), backbonePerComponent[c].back(),
                 contourLengthPerComponent[c]);
    }
  }
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerBackbone &p)
{
  archive << p.versionNumber;

  archive << p.numberOfBlocks;
  archive << p.numberOfComponents;
  archive << p.sampleEvery;
  archive << p.writeEvery;
  archive << p.numberOfWaveVectors;
  archive << p.waveVectorLowerLimit;
  archive << p.waveVectorUpperLimit;
  archive << p.waveVectors;
  archive << p.backbonePerComponent;
  archive << p.contourLengthPerComponent;
  archive << p.internalDistanceSquaredSum;
  archive << p.bondCorrelationSum;
  archive << p.formFactorSum;
  archive << p.sums;
  archive << p.numberOfCounts;
  archive << p.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerBackbone &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPolymerBackbone' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.numberOfBlocks;
  archive >> p.numberOfComponents;
  archive >> p.sampleEvery;
  archive >> p.writeEvery;
  archive >> p.numberOfWaveVectors;
  archive >> p.waveVectorLowerLimit;
  archive >> p.waveVectorUpperLimit;
  archive >> p.waveVectors;
  archive >> p.backbonePerComponent;
  archive >> p.contourLengthPerComponent;
  archive >> p.internalDistanceSquaredSum;
  archive >> p.bondCorrelationSum;
  archive >> p.formFactorSum;
  archive >> p.sums;
  archive >> p.numberOfCounts;
  archive >> p.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyPolymerBackbone: Error in binary restart\n"));
  }
#endif

  return archive;
}
