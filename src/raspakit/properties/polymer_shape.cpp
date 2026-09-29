module;

module property_polymer_shape;

import std;

import archive;
import double3;
import double3x3;
import atom;
import forcefield;
import component;
import averages;
import property_molecule_properties;

namespace
{

// Eigenvalues of a real symmetric 3x3 matrix in descending order, by the trigonometric (Cardano)
// solution of the characteristic polynomial. Exact for diagonal input and stable for the
// positive-semidefinite gyration tensors this is used on.
double3 symmetricEigenvaluesDescending(double a11, double a22, double a33, double a12, double a13, double a23)
{
  double offDiagonal = a12 * a12 + a13 * a13 + a23 * a23;
  if (offDiagonal < 1e-300)
  {
    std::array<double, 3> diagonal{a11, a22, a33};
    std::sort(diagonal.begin(), diagonal.end(), std::greater<double>());
    return double3(diagonal[0], diagonal[1], diagonal[2]);
  }

  double q = (a11 + a22 + a33) / 3.0;
  double d11 = a11 - q, d22 = a22 - q, d33 = a33 - q;
  double p = std::sqrt((d11 * d11 + d22 * d22 + d33 * d33 + 2.0 * offDiagonal) / 6.0);

  // B = (A - qI) / p has det(B)/2 = cos(3 phi).
  double b11 = d11 / p, b22 = d22 / p, b33 = d33 / p;
  double b12 = a12 / p, b13 = a13 / p, b23 = a23 / p;
  double determinant = b11 * (b22 * b33 - b23 * b23) - b12 * (b12 * b33 - b23 * b13) + b13 * (b12 * b23 - b22 * b13);
  double r = std::clamp(0.5 * determinant, -1.0, 1.0);
  double phi = std::acos(r) / 3.0;

  double lambda1 = q + 2.0 * p * std::cos(phi);
  double lambda3 = q + 2.0 * p * std::cos(phi + 2.0 * std::numbers::pi / 3.0);
  double lambda2 = 3.0 * q - lambda1 - lambda3;
  return double3(lambda1, lambda2, lambda3);
}

// Normalized weights for a set of atoms: uniform, or proportional to the pseudo-atom masses.
std::vector<double> normalizedWeights(const ForceField &forceField, const Component &component,
                                      std::span<const std::size_t> atoms, bool massWeighted)
{
  std::vector<double> weights(atoms.size(), 1.0);
  if (massWeighted)
  {
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      weights[i] = forceField.pseudoAtoms[component.atoms[atoms[i]].type].mass;
    }
  }
  double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  if (total <= 0.0)
  {
    std::fill(weights.begin(), weights.end(), 1.0);
    total = static_cast<double>(atoms.size());
  }
  for (double &w : weights) w /= total;
  return weights;
}

// Gyration-tensor descriptors of n weighted points; 'positionAt(i)' yields the i-th position.
template <typename PositionAt>
PropertyPolymerShape::Descriptors descriptorsOf(std::size_t n, PositionAt positionAt, std::span<const double> weights)
{
  PropertyPolymerShape::Descriptors descriptors{};

  double3 center{};
  for (std::size_t i = 0; i < n; ++i) center += weights[i] * positionAt(i);

  double sxx{0.0}, syy{0.0}, szz{0.0}, sxy{0.0}, sxz{0.0}, syz{0.0};
  for (std::size_t i = 0; i < n; ++i)
  {
    double3 d = positionAt(i) - center;
    double w = weights[i];
    sxx += w * d.x * d.x;
    syy += w * d.y * d.y;
    szz += w * d.z * d.z;
    sxy += w * d.x * d.y;
    sxz += w * d.x * d.z;
    syz += w * d.y * d.z;
  }

  descriptors.tensor = double3x3(sxx, sxy, sxz, sxy, syy, syz, sxz, syz, szz);

  double3 lambda = symmetricEigenvaluesDescending(sxx, syy, szz, sxy, sxz, syz);
  // Round-off can push the smallest eigenvalue of a planar/linear conformation slightly negative.
  lambda.x = std::max(lambda.x, 0.0);
  lambda.y = std::max(lambda.y, 0.0);
  lambda.z = std::max(lambda.z, 0.0);
  descriptors.eigenvalues = lambda;

  double rg2 = lambda.x + lambda.y + lambda.z;
  descriptors.radiusOfGyrationSquared = rg2;
  descriptors.asphericity = lambda.x - 0.5 * (lambda.y + lambda.z);
  descriptors.acylindricity = lambda.y - lambda.z;

  if (rg2 > 0.0)
  {
    double b = descriptors.asphericity;
    double c = descriptors.acylindricity;
    descriptors.shapeAnisotropy = std::clamp((b * b + 0.75 * c * c) / (rg2 * rg2), 0.0, 1.0);

    double mean = rg2 / 3.0;
    descriptors.prolateness =
        std::clamp(27.0 * (lambda.x - mean) * (lambda.y - mean) * (lambda.z - mean) / (rg2 * rg2 * rg2), -0.25, 2.0);
  }

  return descriptors;
}

// Unit eigenvector of a symmetric 3x3 matrix for a known eigenvalue: the cross product of the two
// rows of (S - lambda I) with the largest product norm. Returns nullopt when the eigenvalue is
// (numerically) degenerate and the rows are parallel.
std::optional<double3> principalAxis(const double3x3 &s, double lambda)
{
  double3 r0(s.ax - lambda, s.ay, s.az);
  double3 r1(s.ay, s.by - lambda, s.bz);
  double3 r2(s.az, s.bz, s.cz - lambda);
  double3 c01 = double3::cross(r0, r1), c02 = double3::cross(r0, r2), c12 = double3::cross(r1, r2);
  double n01 = c01.length(), n02 = c02.length(), n12 = c12.length();
  double best = std::max({n01, n02, n12});
  double scale = std::max({r0.length(), r1.length(), r2.length(), 1e-300});
  if (best < 1e-10 * scale * scale) return std::nullopt;
  if (best == n01) return c01 / n01;
  if (best == n02) return c02 / n02;
  return c12 / n12;
}

// Radius of gyration of the reference geometry with the given normalized weights.
double referenceRadiusOfGyration(const Component &component, std::span<const double> weights)
{
  double3 center{};
  for (std::size_t i = 0; i < component.atoms.size(); ++i) center += weights[i] * component.atoms[i].position;
  double sum = 0.0;
  for (std::size_t i = 0; i < component.atoms.size(); ++i)
  {
    double3 d = component.atoms[i].position - center;
    sum += weights[i] * double3::dot(d, d);
  }
  return std::sqrt(std::max(sum, 0.0));
}

}  // namespace

PropertyPolymerShape::PropertyPolymerShape(std::size_t numberOfBlocks, const ForceField &forceField,
                                           const std::vector<Component> &components, std::size_t numberOfBins,
                                           bool massWeighted, std::size_t sampleEvery,
                                           std::optional<std::size_t> writeEvery,
                                           std::optional<double> radiusOfGyrationRangeOverride)
    : numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      numberOfComponents(components.size()),
      massWeighted(massWeighted),
      sampleEvery(sampleEvery),
      writeEvery(writeEvery),
      weightsPerComponent(components.size()),
      endToEndAtomsPerComponent(components.size()),
      radiusOfGyrationRangePerComponent(components.size()),
      deltaRadiusOfGyrationPerComponent(components.size()),
      deltaShapeAnisotropy(shapeAnisotropyRange / static_cast<double>(numberOfBins)),
      deltaProlateness(prolatenessRange / static_cast<double>(numberOfBins)),
      radiusOfGyrationHistogram(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      shapeAnisotropyHistogram(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      prolatenessHistogram(numberOfBlocks, std::vector<std::vector<double>>(components.size())),
      sums(numberOfBlocks, std::vector<Moments>(components.size(), Moments{})),
      numberOfCounts(numberOfBlocks, std::vector<double>(components.size())),
      unitAtomsPerComponent(components.size()),
      unitWeightsPerComponent(components.size()),
      unitSums(numberOfBlocks, std::vector<std::vector<Moments>>(components.size())),
      unitCounts(numberOfBlocks, std::vector<std::vector<double>>(components.size()))
{
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    const Component &component = components[c];
    std::size_t numberOfAtoms = component.atoms.size();
    if (numberOfAtoms < 2) continue;

    std::vector<std::size_t> allAtoms(numberOfAtoms);
    std::iota(allAtoms.begin(), allAtoms.end(), 0uz);
    std::vector<double> weights = normalizedWeights(forceField, component, allAtoms, massWeighted);
    weightsPerComponent[c] = weights;

    endToEndAtomsPerComponent[c] = component.endToEndAtoms;

    // Per-monomer sampling for chains declared with repeat units.
    unitAtomsPerComponent[c] = component.repeatUnits;
    unitWeightsPerComponent[c].resize(component.repeatUnits.size());
    for (std::size_t u = 0; u < component.repeatUnits.size(); ++u)
    {
      unitWeightsPerComponent[c][u] = normalizedWeights(forceField, component, component.repeatUnits[u], massWeighted);
    }
    for (std::size_t b = 0; b < numberOfBlocks; ++b)
    {
      unitSums[b][c] = std::vector<Moments>(component.repeatUnits.size(), Moments{});
      unitCounts[b][c] = std::vector<double>(component.repeatUnits.size());
    }

    // Default histogram range: 0.6 times the contour length of the bond-graph diameter. Every atom is
    // within that path length of every other, so the Euclidean extent of the molecule is bounded by
    // it and Rg by 0.61 of it (Jung's theorem); the rod value Rg = L/sqrt(12) sits well inside. When
    // there is no bond graph (a rigid molecule) Rg is constant and twice the reference value is used.
    if (radiusOfGyrationRangeOverride.has_value())
    {
      radiusOfGyrationRangePerComponent[c] = radiusOfGyrationRangeOverride.value();
    }
    else
    {
      std::optional<std::array<std::size_t, 2>> diameter = component.connectivityTable.graphDiameterEndpoints();
      double range = 0.0;
      if (diameter.has_value() && diameter->at(0) != diameter->at(1))
      {
        range = 0.6 * contourLength(component, diameter.value());
      }
      if (range <= 0.0)
      {
        range = 2.0 * referenceRadiusOfGyration(component, weights);
      }
      radiusOfGyrationRangePerComponent[c] = std::max(range, 1.0);
    }
    deltaRadiusOfGyrationPerComponent[c] = radiusOfGyrationRangePerComponent[c] / static_cast<double>(numberOfBins);

    for (std::size_t b = 0; b < numberOfBlocks; ++b)
    {
      radiusOfGyrationHistogram[b][c] = std::vector<double>(numberOfBins);
      shapeAnisotropyHistogram[b][c] = std::vector<double>(numberOfBins);
      prolatenessHistogram[b][c] = std::vector<double>(numberOfBins);
    }
  }
}

PropertyPolymerShape::Descriptors PropertyPolymerShape::computeDescriptors(std::span<const Atom> molecule,
                                                                           std::span<const double> weights)
{
  return descriptorsOf(molecule.size(), [&](std::size_t i) { return molecule[i].position; }, weights);
}

PropertyPolymerShape::Descriptors PropertyPolymerShape::computeDescriptors(std::span<const Atom> molecule,
                                                                           std::span<const std::size_t> atoms,
                                                                           std::span<const double> weights)
{
  return descriptorsOf(atoms.size(), [&](std::size_t i) { return molecule[atoms[i]].position; }, weights);
}

void PropertyPolymerShape::accumulateShapeMoments(Moments &moments, const Descriptors &d)
{
  double rg2 = d.radiusOfGyrationSquared;
  moments[RadiusOfGyration] += std::sqrt(rg2);
  moments[RadiusOfGyrationSquared] += rg2;
  moments[RadiusOfGyrationFourth] += rg2 * rg2;
  moments[Lambda1] += d.eigenvalues.x;
  moments[Lambda2] += d.eigenvalues.y;
  moments[Lambda3] += d.eigenvalues.z;
  moments[Asphericity] += d.asphericity;
  moments[Acylindricity] += d.acylindricity;
  moments[ShapeAnisotropy] += d.shapeAnisotropy;
  moments[SecondInvariant] +=
      d.eigenvalues.x * d.eigenvalues.y + d.eigenvalues.y * d.eigenvalues.z + d.eigenvalues.z * d.eigenvalues.x;
  moments[Prolateness] += d.prolateness;
  moments[TensorXX] += d.tensor.ax;
  moments[TensorYY] += d.tensor.by;
  moments[TensorZZ] += d.tensor.cz;
  moments[TensorXY] += d.tensor.ay;
  moments[TensorXZ] += d.tensor.az;
  moments[TensorYZ] += d.tensor.bz;
}

double PropertyPolymerShape::computeInverseHydrodynamicRadius(std::span<const Atom> molecule)
{
  std::size_t n = molecule.size();
  if (n < 2) return 0.0;

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    for (std::size_t j = i + 1; j < n; ++j)
    {
      double r = (molecule[i].position - molecule[j].position).length();
      if (r > 1e-10) sum += 1.0 / r;
    }
  }
  // The ordered double sum over i != j is twice the sum over unordered pairs.
  return 2.0 * sum / (static_cast<double>(n) * static_cast<double>(n));
}

void PropertyPolymerShape::sample(const std::vector<Component> &components,
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

    const std::vector<double> &weights = weightsPerComponent[c];
    Moments &moments = sums[block][c];

    for (std::size_t m = 0; m < numberOfMolecules; ++m)
    {
      // Positions are stored unwrapped, so plain differences are the physical intramolecular vectors.
      std::span<const Atom> molecule = moleculeAtoms.subspan(offset, numberOfAtoms);
      offset += numberOfAtoms;

      Descriptors d = computeDescriptors(molecule, weights);
      double rg = std::sqrt(d.radiusOfGyrationSquared);

      accumulateShapeMoments(moments, d);
      moments[InverseHydrodynamicRadius] += computeInverseHydrodynamicRadius(molecule);

      // Per-monomer descriptors.
      const std::vector<std::vector<std::size_t>> &units = unitAtomsPerComponent[c];
      for (std::size_t u = 0; u < units.size(); ++u)
      {
        if (units[u].size() < 2) continue;
        accumulateShapeMoments(unitSums[block][c][u], computeDescriptors(molecule, units[u], unitWeightsPerComponent[c][u]));
        unitCounts[block][c][u] += 1.0;
      }

      if (endToEndAtomsPerComponent[c].has_value())
      {
        const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();
        double3 dr = molecule[ends[0]].position - molecule[ends[1]].position;
        moments[EndToEndSquared] += double3::dot(dr, dr);
      }

      std::size_t binRg = static_cast<std::size_t>(rg / deltaRadiusOfGyrationPerComponent[c]);
      if (binRg < numberOfBins) radiusOfGyrationHistogram[block][c][binRg] += 1.0;

      std::size_t binKappa = static_cast<std::size_t>(d.shapeAnisotropy / deltaShapeAnisotropy);
      if (binKappa < numberOfBins) shapeAnisotropyHistogram[block][c][binKappa] += 1.0;

      std::size_t binS = static_cast<std::size_t>((d.prolateness - prolatenessLowerLimit) / deltaProlateness);
      if (binS < numberOfBins) prolatenessHistogram[block][c][binS] += 1.0;

      numberOfCounts[block][c] += 1.0;
    }
  }

  totalNumberOfCounts += 1.0;
}

std::pair<double, double> PropertyPolymerShape::statistics(
    std::size_t component, const std::function<double(const Moments &)> &function) const
{
  return blockStatistics([&](std::size_t block) { return sums[block][component]; },
                         [&](std::size_t block) { return numberOfCounts[block][component]; }, function);
}

std::pair<double, double> PropertyPolymerShape::unitStatistics(
    std::size_t component, std::optional<std::size_t> unit,
    const std::function<double(const Moments &)> &function) const
{
  if (unit.has_value())
  {
    return blockStatistics([&](std::size_t block) { return unitSums[block][component][unit.value()]; },
                           [&](std::size_t block) { return unitCounts[block][component][unit.value()]; }, function);
  }
  // Pooled over all units: every sampled unit of every molecule counts as one sample.
  return blockStatistics(
      [&](std::size_t block)
      {
        Moments pooled{};
        for (const Moments &m : unitSums[block][component])
          for (std::size_t i = 0; i != NumberOfMoments; ++i) pooled[i] += m[i];
        return pooled;
      },
      [&](std::size_t block)
      {
        const std::vector<double> &counts = unitCounts[block][component];
        return std::accumulate(counts.begin(), counts.end(), 0.0);
      },
      function);
}

std::pair<double, double> PropertyPolymerShape::blockStatistics(
    const std::function<Moments(std::size_t)> &sumOf, const std::function<double(std::size_t)> &countOf,
    const std::function<double(const Moments &)> &function) const
{
  Moments total{};
  double totalSamples{0.0};
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    totalSamples += countOf(blockIndex);
    Moments blockSum = sumOf(blockIndex);
    for (std::size_t i = 0; i != NumberOfMoments; ++i) total[i] += blockSum[i];
  }
  if (totalSamples <= 0.0) return {0.0, 0.0};

  Moments overall{};
  for (std::size_t i = 0; i != NumberOfMoments; ++i) overall[i] = total[i] / totalSamples;
  double mean = function(overall);

  std::size_t degreesOfFreedom = numberOfBlocks - 1;
  double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
  double sumOfSquares{0.0};
  std::size_t numberOfSamples{0};
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    double count = countOf(blockIndex);
    if (count > 0.0)
    {
      Moments blockSum = sumOf(blockIndex);
      Moments blockAverage{};
      for (std::size_t i = 0; i != NumberOfMoments; ++i) blockAverage[i] = blockSum[i] / count;
      double value = function(blockAverage) - mean;
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

std::pair<double, double> PropertyPolymerShape::momentStatistics(std::size_t component, Moment moment) const
{
  return statistics(component, [moment](const Moments &m) { return m[moment]; });
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> PropertyPolymerShape::result(
    const std::vector<std::vector<std::vector<double>>> &histogram, std::size_t component, double delta,
    double rangeStart) const
{
  std::vector<double> bins(numberOfBins);
  for (std::size_t bin = 0; bin != numberOfBins; ++bin)
  {
    bins[bin] = rangeStart + (static_cast<double>(bin) + 0.5) * delta;
  }

  double totalSamples{0.0};
  std::vector<double> summedBlocks(numberOfBins);
  for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
  {
    totalSamples += numberOfCounts[blockIndex][component];
    for (std::size_t bin = 0; bin != numberOfBins; ++bin)
    {
      summedBlocks[bin] += histogram[blockIndex][component][bin];
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
        double blockAverage = histogram[blockIndex][component][bin] / (numberOfCounts[blockIndex][component] * delta);
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

void PropertyPolymerShape::writeOutput(std::size_t systemId, const std::vector<Component> &components,
                                       std::size_t currentCycle)
{
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;

  bool anything = false;
  for (std::size_t c = 0; c < numberOfComponents; ++c) anything = anything || isSampled(c);
  if (!anything) return;

  std::filesystem::create_directory("polymer_shape");

  for (std::size_t c = 0; c < components.size() && c < numberOfComponents; ++c)
  {
    if (!isSampled(c)) continue;
    const std::string &name = components[c].name;

    auto [meanRg, errorRg] = momentStatistics(c, RadiusOfGyration);
    auto [meanRg2, errorRg2] = momentStatistics(c, RadiusOfGyrationSquared);
    auto [rmsRg, errorRmsRg] =
        statistics(c, [](const Moments &m) { return std::sqrt(std::max(m[RadiusOfGyrationSquared], 0.0)); });
    auto [meanL1, errorL1] = momentStatistics(c, Lambda1);
    auto [meanL2, errorL2] = momentStatistics(c, Lambda2);
    auto [meanL3, errorL3] = momentStatistics(c, Lambda3);
    auto [ratioL1L3, errorRatioL1L3] =
        statistics(c, [](const Moments &m) { return m[Lambda3] > 0.0 ? m[Lambda1] / m[Lambda3] : 0.0; });
    auto [ratioL2L3, errorRatioL2L3] =
        statistics(c, [](const Moments &m) { return m[Lambda3] > 0.0 ? m[Lambda2] / m[Lambda3] : 0.0; });
    auto [meanB, errorB] = momentStatistics(c, Asphericity);
    auto [meanC, errorC] = momentStatistics(c, Acylindricity);
    auto [relativeB, errorRelativeB] = statistics(
        c, [](const Moments &m) { return m[RadiusOfGyrationSquared] > 0.0 ? m[Asphericity] / m[RadiusOfGyrationSquared] : 0.0; });
    auto [relativeC, errorRelativeC] = statistics(
        c, [](const Moments &m) { return m[RadiusOfGyrationSquared] > 0.0 ? m[Acylindricity] / m[RadiusOfGyrationSquared] : 0.0; });
    auto [meanKappa2, errorKappa2] = momentStatistics(c, ShapeAnisotropy);
    auto [ensembleKappa2, errorEnsembleKappa2] = statistics(
        c, [](const Moments &m)
        { return m[RadiusOfGyrationFourth] > 0.0 ? 1.0 - 3.0 * m[SecondInvariant] / m[RadiusOfGyrationFourth] : 0.0; });
    auto [meanS, errorS] = momentStatistics(c, Prolateness);
    auto [hydrodynamicRadius, errorHydrodynamicRadius] = statistics(
        c, [](const Moments &m) { return m[InverseHydrodynamicRadius] > 0.0 ? 1.0 / m[InverseHydrodynamicRadius] : 0.0; });
    auto [ratioRgRh, errorRatioRgRh] = statistics(
        c, [](const Moments &m)
        { return std::sqrt(std::max(m[RadiusOfGyrationSquared], 0.0)) * m[InverseHydrodynamicRadius]; });

    std::ofstream summary(std::format("polymer_shape/polymer_shape_{}.s{}.txt", name, systemId));
    summary << std::format("# gyration-tensor shape descriptors, component: {}, number of counts: {}\n", name,
                           totalNumberOfCounts);
    summary << std::format("# weights: {}\n", massWeighted ? "pseudo-atom masses" : "uniform per bead");
    summary << "# eigenvalues of the gyration tensor ordered l1 >= l2 >= l3; errors are 95% confidence intervals\n";
    summary << "#\n";
    summary << std::format("<Rg>                       {:<14.6g} +/- {:<12.6g} [Angstrom]\n", meanRg, errorRg);
    summary << std::format("<Rg^2>                     {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanRg2, errorRg2);
    summary << std::format("sqrt(<Rg^2>)               {:<14.6g} +/- {:<12.6g} [Angstrom]\n", rmsRg, errorRmsRg);
    summary << "#\n";
    summary << std::format("<l1>                       {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanL1, errorL1);
    summary << std::format("<l2>                       {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanL2, errorL2);
    summary << std::format("<l3>                       {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanL3, errorL3);
    summary << std::format("<l1>/<l3>                  {:<14.6g} +/- {:<12.6g} [-]   (ideal chain: 11.8)\n",
                           ratioL1L3, errorRatioL1L3);
    summary << std::format("<l2>/<l3>                  {:<14.6g} +/- {:<12.6g} [-]   (ideal chain: 2.7)\n",
                           ratioL2L3, errorRatioL2L3);
    summary << "#\n";
    summary << std::format("<b>   asphericity          {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanB, errorB);
    summary << std::format("<c>   acylindricity        {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", meanC, errorC);
    summary << std::format("<b>/<Rg^2>                 {:<14.6g} +/- {:<12.6g} [-]   (0 sphere, 1 rod)\n",
                           relativeB, errorRelativeB);
    summary << std::format("<c>/<Rg^2>                 {:<14.6g} +/- {:<12.6g} [-]\n", relativeC, errorRelativeC);
    summary << std::format("<k^2> shape anisotropy     {:<14.6g} +/- {:<12.6g} [-]   (0 sphere, 1 rod; per molecule)\n",
                           meanKappa2, errorKappa2);
    summary << std::format("k^2 = 1 - 3<I2>/<Rg^4>     {:<14.6g} +/- {:<12.6g} [-]   (ensemble form, ideal chain: 0.39)\n",
                           ensembleKappa2, errorEnsembleKappa2);
    summary << std::format("<S>   prolateness          {:<14.6g} +/- {:<12.6g} [-]   (-1/4 disk, 0 sphere, 2 rod)\n",
                           meanS, errorS);
    summary << "#\n";
    summary << std::format("Rh    hydrodynamic radius  {:<14.6g} +/- {:<12.6g} [Angstrom]   (Kirkwood, 1/Rh = <(1/N^2) sum 1/r_ij>)\n",
                           hydrodynamicRadius, errorHydrodynamicRadius);
    summary << std::format("sqrt(<Rg^2>)/Rh            {:<14.6g} +/- {:<12.6g} [-]   (Gaussian chain: 1.5, hard sphere: 0.78)\n",
                           ratioRgRh, errorRatioRgRh);

    // The total (ensemble-averaged, lab-frame) gyration tensor and its principal axes.
    {
      std::array<Moment, 6> tensorMoments{TensorXX, TensorYY, TensorZZ, TensorXY, TensorXZ, TensorYZ};
      std::array<double, 6> t{}, e{};
      for (std::size_t i = 0; i < 6; ++i) std::tie(t[i], e[i]) = momentStatistics(c, tensorMoments[i]);
      double3x3 s(t[0], t[3], t[4], t[3], t[1], t[5], t[4], t[5], t[2]);
      double3 lambda = symmetricEigenvaluesDescending(t[0], t[1], t[2], t[3], t[4], t[5]);
      double trace = lambda.x + lambda.y + lambda.z;

      summary << "#\n";
      summary << "# total gyration tensor <S_ab> [Angstrom^2], lab frame (isotropic bulk: <Rg^2>/3 on the diagonal)\n";
      summary << std::format("<S_xx> <S_xy> <S_xz>       {:<14.6g} {:<14.6g} {:<14.6g}   +/- {:<12.6g} {:<12.6g} {:<12.6g}\n",
                             t[0], t[3], t[4], e[0], e[3], e[4]);
      summary << std::format("<S_yx> <S_yy> <S_yz>       {:<14.6g} {:<14.6g} {:<14.6g}   +/- {:<12.6g} {:<12.6g} {:<12.6g}\n",
                             t[3], t[1], t[5], e[3], e[1], e[5]);
      summary << std::format("<S_zx> <S_zy> <S_zz>       {:<14.6g} {:<14.6g} {:<14.6g}   +/- {:<12.6g} {:<12.6g} {:<12.6g}\n",
                             t[4], t[5], t[2], e[4], e[5], e[2]);
      summary << std::format("eigenvalues of <S>         {:<14.6g} {:<14.6g} {:<14.6g}   [Angstrom^2]\n", lambda.x,
                             lambda.y, lambda.z);
      if (trace > 0.0)
      {
        summary << std::format("eigenvalues of <S> / tr    {:<14.6g} {:<14.6g} {:<14.6g}   [-]   (isotropic: 1/3 each)\n",
                               lambda.x / trace, lambda.y / trace, lambda.z / trace);
        summary << std::format("lab-frame anisotropy       {:<14.6g}                                 [-]   ((l1 - l3)/tr of <S>; 0 isotropic)\n",
                               (lambda.x - lambda.z) / trace);
      }
      std::optional<double3> axis1 = principalAxis(s, lambda.x);
      std::optional<double3> axis3 = principalAxis(s, lambda.z);
      if (axis1.has_value() && axis3.has_value())
      {
        double3 axis2 = double3::cross(axis3.value(), axis1.value());
        summary << std::format("principal axis 1 of <S>    ({:.6g}, {:.6g}, {:.6g})\n", axis1->x, axis1->y, axis1->z);
        summary << std::format("principal axis 2 of <S>    ({:.6g}, {:.6g}, {:.6g})\n", axis2.x, axis2.y, axis2.z);
        summary << std::format("principal axis 3 of <S>    ({:.6g}, {:.6g}, {:.6g})\n", axis3->x, axis3->y, axis3->z);
      }
      else
      {
        summary << "# principal axes of <S> undetermined (degenerate eigenvalues: isotropic to within noise)\n";
      }
    }

    if (endToEndAtomsPerComponent[c].has_value())
    {
      const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();
      auto [meanR2, errorR2] = momentStatistics(c, EndToEndSquared);
      auto [ratioR2Rg2, errorRatioR2Rg2] = statistics(
          c, [](const Moments &m)
          { return m[RadiusOfGyrationSquared] > 0.0 ? m[EndToEndSquared] / m[RadiusOfGyrationSquared] : 0.0; });
      summary << "#\n";
      summary << std::format("<R^2> end-to-end ({}, {})   {:<14.6g} +/- {:<12.6g} [Angstrom^2]\n", ends[0], ends[1],
                             meanR2, errorR2);
      summary << std::format("<R^2>/<Rg^2>               {:<14.6g} +/- {:<12.6g} [-]   (ideal chain: 6, good solvent: ~6.3)\n",
                             ratioR2Rg2, errorRatioR2Rg2);
    }

    if (numberOfUnits(c) > 0)
    {
      std::ofstream stream(std::format("polymer_shape/monomer_shape_{}.s{}.txt", name, systemId));
      stream << std::format(
          "# per-repeat-unit gyration-tensor shape descriptors, component: {}, units: {}, number of counts: {}\n", name,
          numberOfUnits(c), totalNumberOfCounts);
      stream << std::format("# weights: {} (renormalized within each unit); errors are 95% confidence intervals\n",
                            massWeighted ? "pseudo-atom masses" : "uniform per bead");
      stream << "# column 1: unit index along the chain ('all' = pooled over units)\n";
      stream << "# column 2: number of atoms in the unit\n";
      stream << "# columns 3, 4: <Rg>, error [Angstrom]\n";
      stream << "# columns 5, 6: <Rg^2>, error [Angstrom^2]\n";
      stream << "# columns 7, 8: <l1>, error [Angstrom^2]\n";
      stream << "# columns 9, 10: <l2>, error [Angstrom^2]\n";
      stream << "# columns 11, 12: <l3>, error [Angstrom^2]\n";
      stream << "# columns 13, 14: <b>/<Rg^2> relative asphericity, error [-]\n";
      stream << "# columns 15, 16: <c>/<Rg^2> relative acylindricity, error [-]\n";
      stream << "# columns 17, 18: <k^2> shape anisotropy, error [-]\n";
      stream << "# columns 19, 20: <S> prolateness, error [-]\n";

      auto writeRow = [&](const std::string &label, std::size_t numberOfAtoms, std::optional<std::size_t> unit)
      {
        auto [rgU, eRgU] = unitStatistics(c, unit, [](const Moments &m) { return m[RadiusOfGyration]; });
        auto [rg2U, eRg2U] = unitStatistics(c, unit, [](const Moments &m) { return m[RadiusOfGyrationSquared]; });
        auto [l1U, eL1U] = unitStatistics(c, unit, [](const Moments &m) { return m[Lambda1]; });
        auto [l2U, eL2U] = unitStatistics(c, unit, [](const Moments &m) { return m[Lambda2]; });
        auto [l3U, eL3U] = unitStatistics(c, unit, [](const Moments &m) { return m[Lambda3]; });
        auto [bU, eBU] = unitStatistics(
            c, unit, [](const Moments &m)
            { return m[RadiusOfGyrationSquared] > 0.0 ? m[Asphericity] / m[RadiusOfGyrationSquared] : 0.0; });
        auto [cU, eCU] = unitStatistics(
            c, unit, [](const Moments &m)
            { return m[RadiusOfGyrationSquared] > 0.0 ? m[Acylindricity] / m[RadiusOfGyrationSquared] : 0.0; });
        auto [kU, eKU] = unitStatistics(c, unit, [](const Moments &m) { return m[ShapeAnisotropy]; });
        auto [sU, eSU] = unitStatistics(c, unit, [](const Moments &m) { return m[Prolateness]; });
        stream << std::format("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n", label, numberOfAtoms,
                              rgU, eRgU, rg2U, eRg2U, l1U, eL1U, l2U, eL2U, l3U, eL3U, bU, eBU, cU, eCU, kU, eKU, sU,
                              eSU);
      };

      for (std::size_t u = 0; u < numberOfUnits(c); ++u)
      {
        writeRow(std::to_string(u), unitAtomsPerComponent[c][u].size(), u);
      }
      std::size_t pooledAtoms = 0;
      for (const auto &unit : unitAtomsPerComponent[c]) pooledAtoms += unit.size();
      writeRow("all", pooledAtoms / numberOfUnits(c), std::nullopt);
    }

    {
      std::ofstream stream(std::format("polymer_shape/radius_of_gyration_{}.s{}.txt", name, systemId));
      stream << std::format("# radius-of-gyration histogram, component: {}, number of counts: {}\n", name,
                            totalNumberOfCounts);
      stream << std::format("# <Rg> = {:g} +/- {:g} [Angstrom], sqrt(<Rg^2>) = {:g} [Angstrom]\n", meanRg, errorRg,
                            rmsRg);
      stream << "# column 1: radius of gyration [Angstrom]\n";
      stream << "# column 2: probability density [1/Angstrom]\n";
      stream << "# column 3: probability density error (95% confidence) [1/Angstrom]\n";
      auto [values, average, error] =
          result(radiusOfGyrationHistogram, c, deltaRadiusOfGyrationPerComponent[c], 0.0);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }

    {
      std::ofstream stream(std::format("polymer_shape/shape_anisotropy_{}.s{}.txt", name, systemId));
      stream << std::format("# relative-shape-anisotropy histogram, component: {}, number of counts: {}\n", name,
                            totalNumberOfCounts);
      stream << std::format("# <k^2> = {:g} +/- {:g} [-]\n", meanKappa2, errorKappa2);
      stream << "# column 1: relative shape anisotropy k^2 [-]\n";
      stream << "# column 2: probability density [-]\n";
      stream << "# column 3: probability density error (95% confidence) [-]\n";
      auto [values, average, error] = result(shapeAnisotropyHistogram, c, deltaShapeAnisotropy, 0.0);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }

    {
      std::ofstream stream(std::format("polymer_shape/prolateness_{}.s{}.txt", name, systemId));
      stream << std::format("# prolateness histogram, component: {}, number of counts: {}\n", name,
                            totalNumberOfCounts);
      stream << std::format("# <S> = {:g} +/- {:g} [-]\n", meanS, errorS);
      stream << "# column 1: prolateness S [-]\n";
      stream << "# column 2: probability density [-]\n";
      stream << "# column 3: probability density error (95% confidence) [-]\n";
      auto [values, average, error] = result(prolatenessHistogram, c, deltaProlateness, prolatenessLowerLimit);
      for (std::size_t bin = 0; bin != numberOfBins; ++bin)
      {
        stream << std::format("{} {} {}\n", values[bin], average[bin], error[bin]);
      }
    }
  }
}

std::string PropertyPolymerShape::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "Polymer-shape (gyration tensor) sampling:\n");
  std::print(stream, "    sample every: {}\n", sampleEvery);
  if (writeEvery.has_value())
  {
    std::print(stream, "    write every: {}\n", writeEvery.value());
  }
  std::print(stream, "    number of bins: {}\n", numberOfBins);
  std::print(stream, "    weights: {}\n", massWeighted ? "pseudo-atom masses" : "uniform per bead");
  for (std::size_t c = 0; c < weightsPerComponent.size(); ++c)
  {
    if (isSampled(c))
    {
      std::print(stream, "    component {}: Rg range: 0 - {:g} [Angstrom]\n", c, radiusOfGyrationRangePerComponent[c]);
      if (numberOfUnits(c) > 0)
      {
        std::print(stream, "    component {}: per-monomer shape over {} repeat units\n", c, numberOfUnits(c));
      }
    }
  }
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerShape &p)
{
  archive << p.versionNumber;

  archive << p.numberOfBlocks;
  archive << p.numberOfBins;
  archive << p.numberOfComponents;
  archive << p.massWeighted;
  archive << p.sampleEvery;
  archive << p.writeEvery;
  archive << p.weightsPerComponent;
  archive << p.endToEndAtomsPerComponent;
  archive << p.radiusOfGyrationRangePerComponent;
  archive << p.deltaRadiusOfGyrationPerComponent;
  archive << p.shapeAnisotropyRange;
  archive << p.prolatenessLowerLimit;
  archive << p.prolatenessRange;
  archive << p.deltaShapeAnisotropy;
  archive << p.deltaProlateness;
  archive << p.radiusOfGyrationHistogram;
  archive << p.shapeAnisotropyHistogram;
  archive << p.prolatenessHistogram;
  archive << p.sums;
  archive << p.numberOfCounts;
  archive << p.totalNumberOfCounts;
  archive << p.unitAtomsPerComponent;
  archive << p.unitWeightsPerComponent;
  archive << p.unitSums;
  archive << p.unitCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerShape &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPolymerShape' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.numberOfBlocks;
  archive >> p.numberOfBins;
  archive >> p.numberOfComponents;
  archive >> p.massWeighted;
  archive >> p.sampleEvery;
  archive >> p.writeEvery;
  archive >> p.weightsPerComponent;
  archive >> p.endToEndAtomsPerComponent;
  archive >> p.radiusOfGyrationRangePerComponent;
  archive >> p.deltaRadiusOfGyrationPerComponent;
  archive >> p.shapeAnisotropyRange;
  archive >> p.prolatenessLowerLimit;
  archive >> p.prolatenessRange;
  archive >> p.deltaShapeAnisotropy;
  archive >> p.deltaProlateness;
  archive >> p.radiusOfGyrationHistogram;
  archive >> p.shapeAnisotropyHistogram;
  archive >> p.prolatenessHistogram;
  archive >> p.sums;
  archive >> p.numberOfCounts;
  archive >> p.totalNumberOfCounts;
  archive >> p.unitAtomsPerComponent;
  archive >> p.unitWeightsPerComponent;
  archive >> p.unitSums;
  archive >> p.unitCounts;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyPolymerShape: Error in binary restart\n"));
  }
#endif

  return archive;
}
