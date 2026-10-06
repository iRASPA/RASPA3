module;

module property_end_to_end_acf;

import std;

import archive;
import double3;
import atom;
import molecule;
import component;
import molecule_property_settings;

namespace
{
// exact integer power, avoids floating-point round-off of std::pow
std::size_t integerPower(std::size_t base, std::size_t exponent)
{
  std::size_t result{1uz};
  for (std::size_t i = 0; i < exponent; ++i) result *= base;
  return result;
}
}  // namespace

PropertyEndToEndAutoCorrelationFunction::PropertyEndToEndAutoCorrelationFunction(
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::size_t numberOfParticles, double timeStep)
    : PropertyEndToEndAutoCorrelationFunction(
          numberOfMoleculesPerComponent,
          [&]
          {
            std::vector<std::optional<std::array<std::size_t, 2>>> ends;
            ends.reserve(components.size());
            for (const Component &component : components) ends.push_back(component.endToEndAtoms);
            return ends;
          }(),
          [&]
          {
            std::vector<std::optional<EndToEndACFSettings>> settings;
            settings.reserve(components.size());
            for (const Component &component : components) settings.push_back(component.endToEndACFSettings);
            return settings;
          }(),
          numberOfParticles, timeStep)
{
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    if (components[c].endToEndACFSettings.has_value() && !components[c].endToEndAtoms.has_value())
    {
      throw std::runtime_error(std::format(
          "[Input reader]: 'ComputeEndToEndACF' is set for component '{}', which has no end-to-end atoms (set "
          "'EndToEndAtoms' in the molecule definition, or use a chain molecule with two ends)\n",
          components[c].name));
    }
  }
}

PropertyEndToEndAutoCorrelationFunction::PropertyEndToEndAutoCorrelationFunction(
    const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    const std::vector<std::optional<std::array<std::size_t, 2>>> &endToEndAtomsPerComponent,
    const std::vector<std::optional<EndToEndACFSettings>> &settingsPerComponent, std::size_t numberOfParticles,
    double timeStep)
    : numberOfMoleculesPerComponent(numberOfMoleculesPerComponent),
      moleculeOffsetPerComponent(numberOfMoleculesPerComponent.size()),
      endToEndAtomsPerComponent(endToEndAtomsPerComponent),
      settingsPerComponent(settingsPerComponent),
      numberOfComponents(numberOfMoleculesPerComponent.size()),
      numberOfParticles(numberOfParticles),
      timeStep(timeStep),
      dataPerComponent(numberOfMoleculesPerComponent.size())
{
  if (endToEndAtomsPerComponent.size() != numberOfComponents || settingsPerComponent.size() != numberOfComponents)
  {
    throw std::runtime_error(
        "PropertyEndToEndAutoCorrelationFunction: one end-to-end atom pair and one settings entry per component are "
        "required\n");
  }

  std::size_t offset{0};
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    moleculeOffsetPerComponent[c] = offset;
    offset += numberOfMoleculesPerComponent[c];

    if (!this->settingsPerComponent[c].has_value()) continue;
    // a component without end-to-end atoms has nothing to correlate
    if (!endToEndAtomsPerComponent[c].has_value())
    {
      this->settingsPerComponent[c].reset();
      continue;
    }

    EndToEndACFSettings &settings = this->settingsPerComponent[c].value();
    settings.sampleEvery = std::max<std::size_t>(1, settings.sampleEvery);
    settings.numberOfBlockElements = std::max<std::size_t>(2, settings.numberOfBlockElements);

    const std::size_t n = settings.numberOfBlockElements;
    ComponentData &data = dataPerComponent[c];
    data.maxNumberOfBlocks = 1;
    data.blockLength.assign(1, 0uz);
    data.acfCount.assign(1, std::vector<std::size_t>(n, 0uz));
    data.blockData.assign(1, std::vector<std::vector<double3>>(numberOfMoleculesPerComponent[c],
                                                               std::vector<double3>(n, double3())));
    data.acf.assign(1, std::vector<double>(n, 0.0));
    data.blockScalars.assign(1, std::vector<std::vector<Scalars>>(numberOfMoleculesPerComponent[c],
                                                                  std::vector<Scalars>(n, Scalars{})));
    data.acfScalars.assign(1, std::vector<Scalars>(n, Scalars{}));
    data.scalarSampled.fill(true);  // cleared for a scalar as soon as a sample omits it
  }
  if (offset != numberOfParticles)
  {
    throw std::runtime_error(std::format(
        "PropertyEndToEndAutoCorrelationFunction: the numbers of molecules per component sum to {}, not {}\n", offset,
        numberOfParticles));
  }
}

bool PropertyEndToEndAutoCorrelationFunction::hasData(std::size_t component) const
{
  return isSampled(component) && numberOfMoleculesPerComponent[component] > 0;
}

double PropertyEndToEndAutoCorrelationFunction::lagOf(std::size_t component, std::size_t block, std::size_t k) const
{
  const double unit = usesCycles() ? 1.0 : timeStep;
  return static_cast<double>(k) * static_cast<double>(sampleEvery(component)) * unit *
         static_cast<double>(integerPower(numberOfBlockElements(component), block));
}

void PropertyEndToEndAutoCorrelationFunction::addSample(std::size_t currentCycle,
                                                        const std::vector<Molecule> &molecules,
                                                        std::span<const Atom> atoms)
{
  std::vector<std::size_t> due;
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    if (isSampled(c) && currentCycle % sampleEvery(c) == 0uz) due.push_back(c);
  }
  if (due.empty()) return;

  // the molecules are tracked by index and the number of molecules must stay fixed
  if (molecules.size() != numberOfParticles)
  {
    throw std::runtime_error(std::format(
        "PropertyEndToEndAutoCorrelationFunction: the number of molecules changed from {} to {}; computing the "
        "end-to-end autocorrelation function requires a fixed number of molecules (do not combine "
        "'ComputeEndToEndACF' with insertion/deletion moves)\n",
        numberOfParticles, molecules.size()));
  }

  std::vector<double3> vectors(numberOfParticles, double3(0.0, 0.0, 0.0));
  std::vector<double> radiiOfGyrationSquared(numberOfParticles, 0.0);
  for (std::size_t c : due)
  {
    const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();
    for (std::size_t m = 0; m < numberOfMoleculesPerComponent[c]; ++m)
    {
      const std::size_t moleculeIndex = moleculeOffsetPerComponent[c] + m;
      const Molecule &molecule = molecules[moleculeIndex];
      if (ends[0] >= molecule.numberOfAtoms || ends[1] >= molecule.numberOfAtoms)
      {
        throw std::runtime_error(std::format(
            "PropertyEndToEndAutoCorrelationFunction: end-to-end atoms ({}, {}) out of range for a molecule of {} "
            "atoms\n",
            ends[0], ends[1], molecule.numberOfAtoms));
      }
      // positions are stored unwrapped, so the plain difference is the physical end-to-end vector
      vectors[moleculeIndex] =
          atoms[molecule.atomIndex + ends[1]].position - atoms[molecule.atomIndex + ends[0]].position;
      radiiOfGyrationSquared[moleculeIndex] = radiusOfGyrationSquared(molecule, atoms);
    }
  }
  for (std::size_t c : due) addSampleVectors(c, vectors, radiiOfGyrationSquared);
}

double PropertyEndToEndAutoCorrelationFunction::radiusOfGyrationSquared(const Molecule &molecule,
                                                                        std::span<const Atom> atoms)
{
  if (molecule.numberOfAtoms == 0) return 0.0;
  const std::span<const Atom> moleculeAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);

  double3 center(0.0, 0.0, 0.0);
  for (const Atom &atom : moleculeAtoms) center += atom.position;
  center /= static_cast<double>(moleculeAtoms.size());

  double sum = 0.0;
  for (const Atom &atom : moleculeAtoms)
  {
    const double3 d = atom.position - center;
    sum += double3::dot(d, d);
  }
  return sum / static_cast<double>(moleculeAtoms.size());
}

void PropertyEndToEndAutoCorrelationFunction::addSampleVectors(std::span<const double3> endToEndVectors,
                                                               std::span<const double> radiusOfGyrationSquared)
{
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    if (isSampled(c)) addSampleVectors(c, endToEndVectors, radiusOfGyrationSquared);
  }
}

void PropertyEndToEndAutoCorrelationFunction::addSampleVectors(std::size_t component,
                                                               std::span<const double3> endToEndVectors,
                                                               std::span<const double> radiusOfGyrationSquared)
{
  if (endToEndVectors.size() != numberOfParticles)
  {
    throw std::runtime_error(std::format(
        "PropertyEndToEndAutoCorrelationFunction: {} end-to-end vectors given for {} molecules\n",
        endToEndVectors.size(), numberOfParticles));
  }
  if (!radiusOfGyrationSquared.empty() && radiusOfGyrationSquared.size() != numberOfParticles)
  {
    throw std::runtime_error(std::format(
        "PropertyEndToEndAutoCorrelationFunction: {} squared radii of gyration given for {} molecules\n",
        radiusOfGyrationSquared.size(), numberOfParticles));
  }
  if (!isSampled(component)) return;

  const std::size_t n = numberOfBlockElements(component);
  const std::size_t numberOfMolecules = numberOfMoleculesPerComponent[component];
  const std::size_t offset = moleculeOffsetPerComponent[component];
  ComponentData &data = dataPerComponent[component];

  // number of blocks needed for the current count
  data.numberOfBlocks = 1;
  std::size_t p = data.count / n;
  while (p != 0)
  {
    ++data.numberOfBlocks;
    p /= n;
  }

  if (data.numberOfBlocks > data.maxNumberOfBlocks)
  {
    data.blockLength.resize(data.numberOfBlocks, 0uz);
    data.acfCount.resize(data.numberOfBlocks, std::vector<std::size_t>(n, 0uz));
    data.blockData.resize(data.numberOfBlocks,
                          std::vector<std::vector<double3>>(numberOfMolecules, std::vector<double3>(n, double3())));
    data.acf.resize(data.numberOfBlocks, std::vector<double>(n, 0.0));
    data.blockScalars.resize(data.numberOfBlocks,
                             std::vector<std::vector<Scalars>>(numberOfMolecules, std::vector<Scalars>(n, Scalars{})));
    data.acfScalars.resize(data.numberOfBlocks, std::vector<Scalars>(n, Scalars{}));
    data.maxNumberOfBlocks = data.numberOfBlocks;
  }

  // the rotation-invariant scalars of this sample; an Rg^2 that is not supplied enters as zero and marks the
  // channel as unavailable (NaN is not used as a marker because the code may be compiled with fast-math)
  if (radiusOfGyrationSquared.empty()) data.scalarSampled[RadiusOfGyrationSquared] = false;
  auto scalarsOf = [&](std::size_t m) -> Scalars
  {
    const double3 R = endToEndVectors[offset + m];
    Scalars scalars{};
    scalars[EndToEndSquared] = double3::dot(R, R);
    scalars[RadiusOfGyrationSquared] = radiusOfGyrationSquared.empty() ? 0.0 : radiusOfGyrationSquared[offset + m];
    return scalars;
  };

  for (std::size_t block = 0; block < data.numberOfBlocks; ++block)
  {
    // block 'block' takes a sample when the count is a multiple of n^block
    if (data.count % integerPower(n, block) != 0) continue;

    ++data.blockLength[block];
    const std::size_t currentBlockLength = std::min(data.blockLength[block], n);

    for (std::size_t m = 0; m < numberOfMolecules; ++m)
    {
      const double3 value = endToEndVectors[offset + m];
      std::vector<double3> &history = data.blockData[block][m];
      std::shift_right(history.begin(), history.end(), 1);
      history[0] = value;

      const Scalars scalars = scalarsOf(m);
      std::vector<Scalars> &scalarHistory = data.blockScalars[block][m];
      std::shift_right(scalarHistory.begin(), scalarHistory.end(), 1);
      scalarHistory[0] = scalars;

      // k = 0 is the zero lag, <R^2>; it is kept (needed for the normalization) but only reported from block 0
      for (std::size_t k = 0; k < currentBlockLength; ++k)
      {
        ++data.acfCount[block][k];
        data.acf[block][k] += double3::dot(history[k], value);
        for (std::size_t s = 0; s < NumberOfScalars; ++s)
        {
          data.acfScalars[block][k][s] += scalarHistory[k][s] * scalars[s];
        }
      }
    }
  }

  // running mean of the scalars (every sample enters once, through block 0)
  for (std::size_t m = 0; m < numberOfMolecules; ++m)
  {
    const Scalars scalars = scalarsOf(m);
    for (std::size_t s = 0; s < NumberOfScalars; ++s) data.sumScalars[s] += scalars[s];
  }
  data.numberOfScalarSamples += numberOfMolecules;

  ++data.count;
}

namespace
{
// mean and variance of scalar s of a component; the variance is not available (nullopt) when the scalar does not
// fluctuate or was not sampled
struct ScalarMoments
{
  double mean{0.0};
  std::optional<double> variance{};
};

ScalarMoments scalarMoments(const PropertyEndToEndAutoCorrelationFunction::ComponentData &d, std::size_t s)
{
  ScalarMoments moments{};
  if (!d.scalarSampled[s] || d.numberOfScalarSamples == 0 || d.acfCount.empty() || d.acfCount[0][0] == 0)
  {
    return moments;
  }
  moments.mean = d.sumScalars[s] / static_cast<double>(d.numberOfScalarSamples);
  const double meanSquare = d.acfScalars[0][0][s] / static_cast<double>(d.acfCount[0][0]);
  const double variance = meanSquare - moments.mean * moments.mean;
  // a relative tolerance guards against round-off for a frozen conformation (variance ~ eps * <x^2>)
  if (variance > 1e-12 * std::max(meanSquare, 1.0)) moments.variance = variance;
  return moments;
}
}  // namespace

std::vector<EndToEndAutoCorrelationFunctionData> PropertyEndToEndAutoCorrelationFunction::result(
    std::size_t component) const
{
  std::vector<EndToEndAutoCorrelationFunctionData> data;
  if (!hasData(component)) return data;

  const ComponentData &d = dataPerComponent[component];
  if (d.count == 0 || d.acfCount[0][0] == 0) return data;
  const std::size_t n = numberOfBlockElements(component);
  const double meanSquared = d.acf[0][0] / static_cast<double>(d.acfCount[0][0]);
  const double inverseMeanSquared = meanSquared > 0.0 ? 1.0 / meanSquared : 0.0;

  // fluctuation functions (<x(0)x(t)> - <x>^2) / (<x^2> - <x>^2) of the scalars
  std::array<ScalarMoments, NumberOfScalars> moments{};
  for (std::size_t s = 0; s < NumberOfScalars; ++s) moments[s] = scalarMoments(d, s);

  for (std::size_t block = 0; block < d.numberOfBlocks; ++block)
  {
    const std::size_t currentBlockLength = std::min(d.blockLength[block], n);
    for (std::size_t k = (block == 0 ? 0 : 1); k < currentBlockLength; ++k)
    {
      const std::size_t count = d.acfCount[block][k];
      if (count == 0) continue;
      const double value = d.acf[block][k] / static_cast<double>(count);
      EndToEndAutoCorrelationFunctionData point{lagOf(component, block, k), value, value * inverseMeanSquared,
                                                static_cast<double>(count)};
      std::array<std::optional<double>, NumberOfScalars> scalarACF{};
      for (std::size_t s = 0; s < NumberOfScalars; ++s)
      {
        if (!moments[s].variance.has_value()) continue;
        const double correlation = d.acfScalars[block][k][s] / static_cast<double>(count);
        scalarACF[s] = (correlation - moments[s].mean * moments[s].mean) / moments[s].variance.value();
      }
      point.endToEndSquaredACF = scalarACF[EndToEndSquared];
      point.radiusOfGyrationSquaredACF = scalarACF[RadiusOfGyrationSquared];
      data.push_back(point);
    }
  }
  // the blocks are already in increasing order of lag (block b+1 starts at n^(b+1) > (n-1) n^b)
  std::ranges::sort(data, {}, &EndToEndAutoCorrelationFunctionData::time);
  return data;
}

EndToEndRelaxationTimes PropertyEndToEndAutoCorrelationFunction::relaxationTimes(std::size_t component) const
{
  EndToEndRelaxationTimes times = relaxationTimes(result(component));
  if (!hasData(component)) return times;

  // the static version only sees the normalized functions; add the moments of the scalars
  const ComponentData &d = dataPerComponent[component];
  for (std::size_t s = 0; s < NumberOfScalars; ++s)
  {
    ScalarRelaxation &relaxation = (s == EndToEndSquared) ? times.endToEndSquared : times.radiusOfGyrationSquared;
    const ScalarMoments moments = scalarMoments(d, s);
    relaxation.mean = moments.mean;
    relaxation.variance = moments.variance.value_or(0.0);
  }
  return times;
}

EndToEndRelaxationTimes PropertyEndToEndAutoCorrelationFunction::relaxationTimes(
    const std::vector<EndToEndAutoCorrelationFunctionData> &data)
{
  EndToEndRelaxationTimes times{};
  if (data.empty()) return times;

  times.meanSquaredEndToEnd = data.front().acf;
  times.longestLag = data.back().time;

  std::vector<double> lags;
  lags.reserve(data.size());
  for (const EndToEndAutoCorrelationFunctionData &point : data) lags.push_back(point.time);

  if (data.front().acf > 0.0)
  {
    std::vector<double> normalized;
    normalized.reserve(data.size());
    for (const EndToEndAutoCorrelationFunctionData &point : data) normalized.push_back(point.normalized);
    const RelaxationTimeEstimates estimates = estimateRelaxationTimes(lags, normalized);
    times.integrated = estimates.integrated;
    times.integratedIsLowerBound = estimates.integratedIsLowerBound;
    times.oneOverE = estimates.oneOverE;
    times.exponentialFit = estimates.exponentialFit;
  }

  auto scalarRelaxation = [&](std::optional<double> EndToEndAutoCorrelationFunctionData::*column) -> ScalarRelaxation
  {
    ScalarRelaxation relaxation{};
    if (!(data.front().*column).has_value()) return relaxation;
    std::vector<double> normalized;
    normalized.reserve(data.size());
    for (const EndToEndAutoCorrelationFunctionData &point : data) normalized.push_back((point.*column).value_or(0.0));
    relaxation.available = true;
    relaxation.times = estimateRelaxationTimes(lags, normalized);
    return relaxation;
  };
  times.endToEndSquared = scalarRelaxation(&EndToEndAutoCorrelationFunctionData::endToEndSquaredACF);
  times.radiusOfGyrationSquared = scalarRelaxation(&EndToEndAutoCorrelationFunctionData::radiusOfGyrationSquaredACF);

  return times;
}

RelaxationTimeEstimates PropertyEndToEndAutoCorrelationFunction::estimateRelaxationTimes(
    std::span<const double> lags, std::span<const double> normalized)
{
  RelaxationTimeEstimates times{};
  const std::size_t size = std::min(lags.size(), normalized.size());
  if (size < 2) return times;

  // integrated correlation time up to the first zero crossing (trapezoid rule on the normalized function)
  {
    double integral = 0.0;
    bool crossed = false;
    for (std::size_t i = 1; i < size; ++i)
    {
      const double t0 = lags[i - 1], t1 = lags[i];
      const double c0 = normalized[i - 1], c1 = normalized[i];
      if (c1 <= 0.0)
      {
        // linear interpolation of the zero crossing; the area of the remaining triangle
        const double tc = (c0 > 0.0 && c0 != c1) ? t0 + (t1 - t0) * c0 / (c0 - c1) : t0;
        integral += 0.5 * c0 * (tc - t0);
        crossed = true;
        break;
      }
      integral += 0.5 * (c0 + c1) * (t1 - t0);
    }
    times.integrated = integral;
    times.integratedIsLowerBound = !crossed;
  }

  // 1/e time
  {
    const double threshold = std::exp(-1.0);
    for (std::size_t i = 1; i < size; ++i)
    {
      const double c0 = normalized[i - 1], c1 = normalized[i];
      if (c1 <= threshold)
      {
        const double t0 = lags[i - 1], t1 = lags[i];
        times.oneOverE = (c0 != c1) ? t0 + (t1 - t0) * (c0 - threshold) / (c0 - c1) : t1;
        break;
      }
    }
  }

  // single-exponential fit of the tail: ln C(t)/C(0) = a - t / tau over the leading, contiguous stretch of lags
  // with 0.05 < C/C(0) <= 0.5; the fit stops at the first lag below 0.05 so that the statistical noise of the
  // long lags (which fluctuates around zero and can re-enter the window) does not enter
  {
    double sumT = 0.0, sumY = 0.0, sumTT = 0.0, sumTY = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
      const double t = lags[i], c = normalized[i];
      if (t <= 0.0 || c > 0.5) continue;
      if (c <= 0.05) break;
      const double y = std::log(c);
      sumT += t;
      sumY += y;
      sumTT += t * t;
      sumTY += t * y;
      ++n;
    }
    if (n >= 3)
    {
      const double denominator = static_cast<double>(n) * sumTT - sumT * sumT;
      if (denominator > 0.0)
      {
        const double slope = (static_cast<double>(n) * sumTY - sumT * sumY) / denominator;
        if (slope < 0.0) times.exponentialFit = -1.0 / slope;
      }
    }
  }

  return times;
}

void PropertyEndToEndAutoCorrelationFunction::writeOutput(std::size_t systemId, const std::vector<Component> &components,
                                                          std::size_t currentCycle) const
{
  bool anything = false;
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    anything = anything || (hasData(c) && writeEvery(c).has_value() && currentCycle % writeEvery(c).value() == 0uz &&
                            dataPerComponent[c].count > 0uz);
  }
  if (!anything) return;

  std::filesystem::create_directory("end_to_end_acf");

  const std::string unit = timeUnit();
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    if (!hasData(c)) continue;
    if (!writeEvery(c).has_value() || currentCycle % writeEvery(c).value() != 0uz) continue;
    const std::vector<EndToEndAutoCorrelationFunctionData> data = result(c);
    if (data.empty()) continue;
    const EndToEndRelaxationTimes times = relaxationTimes(c);
    const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();

    std::ofstream stream(std::format("end_to_end_acf/end_to_end_acf_{}.s{}.txt", components[c].name, systemId));

    stream << std::format(
        "# end-to-end vector autocorrelation function <R(0).R(t)> and conformational (R^2, Rg^2) autocorrelation "
        "functions of component '{}'\n",
        components[c].name);
    stream << std::format("# end-to-end atoms: {} and {}; {} molecules; {} time origins sampled every {} cycles\n",
                          ends[0], ends[1], numberOfMoleculesPerComponent[c], dataPerComponent[c].count,
                          sampleEvery(c));
    stream << std::format("# order-N blocking with {} elements per block\n", numberOfBlockElements(c));
    stream << std::format("# <R^2> = C(0) = {:.6f} [A^2], <R^2>^(1/2) = {:.6f} [A]\n", times.meanSquaredEndToEnd,
                          std::sqrt(std::max(0.0, times.meanSquaredEndToEnd)));
    stream << std::format("# longest lag: {} [{}]\n", times.longestLag, unit);
    stream << "# end-to-end relaxation time tau_R (spacing of independent samples of R):\n";
    if (times.integrated.has_value())
    {
      stream << std::format("#   integrated  tau_int = int C(t)/C(0) dt = {:.6g} [{}]{}\n", times.integrated.value(),
                            unit,
                            times.integratedIsLowerBound
                                ? " (lower bound: C(t) has not crossed zero within the longest lag)"
                                : " (up to the first zero crossing)");
    }
    else
    {
      stream << "#   integrated  tau_int = n/a\n";
    }
    if (times.oneOverE.has_value())
    {
      stream << std::format("#   1/e time    tau_e   = {:.6g} [{}]\n", times.oneOverE.value(), unit);
    }
    else
    {
      stream << std::format("#   1/e time    tau_e   > {} [{}] (C(t)/C(0) has not decayed to 1/e)\n",
                            times.longestLag, unit);
    }
    if (times.exponentialFit.has_value())
    {
      stream << std::format("#   exp. fit    tau_fit = {:.6g} [{}] (single exponential over 0.05 < C/C(0) <= 0.5)\n",
                            times.exponentialFit.value(), unit);
    }
    else
    {
      stream << "#   exp. fit    tau_fit = n/a (fewer than three lags with 0.05 < C/C(0) <= 0.5)\n";
    }

    // conformational relaxation: the vector function above also decays when an unchanged conformation tumbles;
    // the fluctuation functions of the rotation-invariant R^2 and Rg^2 only decay when the chain changes shape
    stream << "# conformational relaxation (rotation-invariant; stays at 1 while the conformation is unchanged, "
              "e.g. a closed hairpin that only tumbles):\n";
    auto writeScalar = [&](const std::string &symbol, const ScalarRelaxation &relaxation)
    {
      if (!relaxation.available)
      {
        stream << std::format("#   {}: n/a (no fluctuations sampled)\n", symbol);
        return;
      }
      stream << std::format("#   <{}> = {:.6f} [A^2], <{}>^(1/2) = {:.6f} [A], std.dev. = {:.6f} [A^2]\n", symbol,
                            relaxation.mean, symbol, std::sqrt(std::max(0.0, relaxation.mean)),
                            std::sqrt(std::max(0.0, relaxation.variance)));
      const RelaxationTimeEstimates &estimates = relaxation.times;
      if (estimates.integrated.has_value())
      {
        stream << std::format("#   {} integrated  tau_int = {:.6g} [{}]{}\n", symbol, estimates.integrated.value(),
                              unit, estimates.integratedIsLowerBound ? " (lower bound: no zero crossing yet)" : "");
      }
      if (estimates.oneOverE.has_value())
      {
        stream << std::format("#   {} 1/e time    tau_e   = {:.6g} [{}]\n", symbol, estimates.oneOverE.value(), unit);
      }
      else
      {
        stream << std::format("#   {} 1/e time    tau_e   > {} [{}] (not decayed to 1/e: the conformation has not "
                              "relaxed within the longest lag)\n",
                              symbol, times.longestLag, unit);
      }
      if (estimates.exponentialFit.has_value())
      {
        stream << std::format("#   {} exp. fit    tau_fit = {:.6g} [{}]\n", symbol, estimates.exponentialFit.value(),
                              unit);
      }
    };
    writeScalar("R^2", times.endToEndSquared);
    writeScalar("Rg^2", times.radiusOfGyrationSquared);

    stream << std::format("# column 1: time [{}]\n", unit);
    stream << "# column 2: <R(0).R(t)> [A^2]\n";
    stream << "# column 3: <R(0).R(t)>/<R^2> [-]\n";
    stream << "# column 4: <dR^2(0) dR^2(t)>/<(dR^2)^2> [-], dR^2 = R^2 - <R^2> (end-to-end distance "
              "fluctuations)\n";
    stream << "# column 5: <dRg^2(0) dRg^2(t)>/<(dRg^2)^2> [-], dRg^2 = Rg^2 - <Rg^2> (radius-of-gyration "
              "fluctuations)\n";
    stream << "# column 6: number of samples [-]\n";

    auto column = [](const std::optional<double> &value)
    { return value.has_value() ? std::format("{}", value.value()) : std::string("n/a"); };
    for (const EndToEndAutoCorrelationFunctionData &point : data)
    {
      stream << std::format("{} {} {} {} {} (count: {})\n", point.time, point.acf, point.normalized,
                            column(point.endToEndSquaredACF), column(point.radiusOfGyrationSquaredACF),
                            static_cast<std::size_t>(point.numberOfSamples));
    }
  }
}

std::string PropertyEndToEndAutoCorrelationFunction::printSettings() const
{
  std::ostringstream stream;

  std::print(stream, "End-to-end vector and conformational (R^2, Rg^2) autocorrelation functions (order-N):\n");
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    if (!isSampled(c)) continue;
    const EndToEndACFSettings &settings = settingsPerComponent[c].value();
    std::print(stream, "    component {}: sample every {}", c, settings.sampleEvery);
    if (settings.writeEvery.has_value())
    {
      std::print(stream, ", write every {}", settings.writeEvery.value());
    }
    std::print(stream, ", {} elements per block, end-to-end atoms ({}, {})\n", settings.numberOfBlockElements,
               endToEndAtomsPerComponent[c]->at(0), endToEndAtomsPerComponent[c]->at(1));
  }
  std::print(stream, "\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                   const PropertyEndToEndAutoCorrelationFunction::ComponentData &d)
{
  archive << d.count;
  archive << d.numberOfBlocks;
  archive << d.maxNumberOfBlocks;
  archive << d.blockLength;
  archive << d.acfCount;
  archive << d.blockData;
  archive << d.acf;
  archive << d.blockScalars;
  archive << d.acfScalars;
  archive << d.sumScalars;
  archive << d.numberOfScalarSamples;
  archive << d.scalarSampled;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEndToEndAutoCorrelationFunction::ComponentData &d)
{
  archive >> d.count;
  archive >> d.numberOfBlocks;
  archive >> d.maxNumberOfBlocks;
  archive >> d.blockLength;
  archive >> d.acfCount;
  archive >> d.blockData;
  archive >> d.acf;
  archive >> d.blockScalars;
  archive >> d.acfScalars;
  archive >> d.sumScalars;
  archive >> d.numberOfScalarSamples;
  archive >> d.scalarSampled;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEndToEndAutoCorrelationFunction &p)
{
  archive << p.versionNumber;

  archive << p.numberOfMoleculesPerComponent;
  archive << p.moleculeOffsetPerComponent;
  archive << p.endToEndAtomsPerComponent;
  archive << p.settingsPerComponent;
  archive << p.numberOfComponents;
  archive << p.numberOfParticles;
  archive << p.timeStep;
  archive << p.dataPerComponent;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEndToEndAutoCorrelationFunction &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber != p.versionNumber)
  {
    // version 1 held one system-wide schedule and block structure; it cannot be split per component;
    // version 2 lacks the conformational (R^2, Rg^2) accumulators
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format(
        "Invalid version {} reading 'PropertyEndToEndAutoCorrelationFunction' (expected {}) at line {} in file {}\n",
        versionNumber, p.versionNumber, location.line(), location.file_name()));
  }

  archive >> p.numberOfMoleculesPerComponent;
  archive >> p.moleculeOffsetPerComponent;
  archive >> p.endToEndAtomsPerComponent;
  archive >> p.settingsPerComponent;
  archive >> p.numberOfComponents;
  archive >> p.numberOfParticles;
  archive >> p.timeStep;
  archive >> p.dataPerComponent;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyEndToEndAutoCorrelationFunction: Error in binary restart\n"));
  }
#endif

  return archive;
}
