module;

module property_end_to_end_acf;

import std;

import archive;
import double3;
import atom;
import molecule;
import component;

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
    const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    const std::vector<std::optional<std::array<std::size_t, 2>>> &endToEndAtomsPerComponent,
    std::size_t numberOfParticles, double timeStep, std::size_t numberOfBlockElements, std::size_t sampleEvery,
    std::optional<std::size_t> writeEvery)
    : numberOfMoleculesPerComponent(numberOfMoleculesPerComponent),
      endToEndAtomsPerComponent(endToEndAtomsPerComponent),
      numberOfComponents(numberOfMoleculesPerComponent.size()),
      numberOfParticles(numberOfParticles),
      timeStep(timeStep),
      numberOfBlockElements(std::max<std::size_t>(2, numberOfBlockElements)),
      sampleEvery(std::max<std::size_t>(1, sampleEvery)),
      writeEvery(writeEvery),
      maxNumberOfBlocks(1),
      blockLength(1, 0uz),
      acfCount(1, std::vector<std::vector<std::size_t>>(
                      numberOfComponents, std::vector<std::size_t>(this->numberOfBlockElements, 0uz))),
      blockData(1, std::vector<std::vector<double3>>(numberOfParticles,
                                                     std::vector<double3>(this->numberOfBlockElements, double3()))),
      acf(1, std::vector<std::vector<double>>(numberOfComponents, std::vector<double>(this->numberOfBlockElements, 0.0)))
{
  if (endToEndAtomsPerComponent.size() != numberOfComponents)
  {
    throw std::runtime_error(
        "PropertyEndToEndAutoCorrelationFunction: one end-to-end atom pair per component is required\n");
  }
}

bool PropertyEndToEndAutoCorrelationFunction::hasData(std::size_t component) const
{
  return component < numberOfComponents && endToEndAtomsPerComponent[component].has_value() &&
         numberOfMoleculesPerComponent[component] > 0;
}

double PropertyEndToEndAutoCorrelationFunction::lagOf(std::size_t block, std::size_t k) const
{
  const double unit = usesCycles() ? 1.0 : timeStep;
  return static_cast<double>(k) * static_cast<double>(sampleEvery) * unit *
         static_cast<double>(integerPower(numberOfBlockElements, block));
}

void PropertyEndToEndAutoCorrelationFunction::addSample(std::size_t currentCycle,
                                                        const std::vector<Molecule> &molecules,
                                                        std::span<const Atom> atoms)
{
  if (currentCycle % sampleEvery != 0uz) return;

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
  std::size_t moleculeIndex{0};
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    const std::optional<std::array<std::size_t, 2>> &ends = endToEndAtomsPerComponent[c];
    for (std::size_t m = 0; m < numberOfMoleculesPerComponent[c]; ++m, ++moleculeIndex)
    {
      if (!ends.has_value()) continue;
      const Molecule &molecule = molecules[moleculeIndex];
      if (ends->at(0) >= molecule.numberOfAtoms || ends->at(1) >= molecule.numberOfAtoms)
      {
        throw std::runtime_error(std::format(
            "PropertyEndToEndAutoCorrelationFunction: end-to-end atoms ({}, {}) out of range for a molecule of {} "
            "atoms\n",
            ends->at(0), ends->at(1), molecule.numberOfAtoms));
      }
      // positions are stored unwrapped, so the plain difference is the physical end-to-end vector
      vectors[moleculeIndex] =
          atoms[molecule.atomIndex + ends->at(1)].position - atoms[molecule.atomIndex + ends->at(0)].position;
    }
  }
  addSampleVectors(vectors);
}

void PropertyEndToEndAutoCorrelationFunction::addSampleVectors(std::span<const double3> endToEndVectors)
{
  if (endToEndVectors.size() != numberOfParticles)
  {
    throw std::runtime_error(std::format(
        "PropertyEndToEndAutoCorrelationFunction: {} end-to-end vectors given for {} molecules\n",
        endToEndVectors.size(), numberOfParticles));
  }

  // number of blocks needed for the current count
  numberOfBlocks = 1;
  std::size_t p = count / numberOfBlockElements;
  while (p != 0)
  {
    ++numberOfBlocks;
    p /= numberOfBlockElements;
  }

  if (numberOfBlocks > maxNumberOfBlocks)
  {
    blockLength.resize(numberOfBlocks, 0uz);
    acfCount.resize(numberOfBlocks, std::vector<std::vector<std::size_t>>(
                                        numberOfComponents, std::vector<std::size_t>(numberOfBlockElements, 0uz)));
    blockData.resize(numberOfBlocks, std::vector<std::vector<double3>>(
                                         numberOfParticles, std::vector<double3>(numberOfBlockElements, double3())));
    acf.resize(numberOfBlocks,
               std::vector<std::vector<double>>(numberOfComponents, std::vector<double>(numberOfBlockElements, 0.0)));
    maxNumberOfBlocks = numberOfBlocks;
  }

  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    // block 'block' takes a sample when the count is a multiple of n^block
    if (count % integerPower(numberOfBlockElements, block) != 0) continue;

    ++blockLength[block];
    const std::size_t currentBlockLength = std::min(blockLength[block], numberOfBlockElements);

    std::size_t moleculeIndex{0};
    for (std::size_t c = 0; c < numberOfComponents; ++c)
    {
      const bool active = endToEndAtomsPerComponent[c].has_value();
      for (std::size_t m = 0; m < numberOfMoleculesPerComponent[c]; ++m, ++moleculeIndex)
      {
        if (!active) continue;
        const double3 value = endToEndVectors[moleculeIndex];
        std::vector<double3> &history = blockData[block][moleculeIndex];
        std::shift_right(history.begin(), history.end(), 1);
        history[0] = value;

        // k = 0 is the zero lag, <R^2>; it is kept (needed for the normalization) but only reported from block 0
        for (std::size_t k = 0; k < currentBlockLength; ++k)
        {
          ++acfCount[block][c][k];
          acf[block][c][k] += double3::dot(history[k], value);
        }
      }
    }
  }

  ++count;
}

std::vector<EndToEndAutoCorrelationFunctionData> PropertyEndToEndAutoCorrelationFunction::result(
    std::size_t component) const
{
  std::vector<EndToEndAutoCorrelationFunctionData> data;
  if (!hasData(component) || count == 0) return data;

  const std::size_t c = component;
  if (acfCount[0][c][0] == 0) return data;
  const double meanSquared = acf[0][c][0] / static_cast<double>(acfCount[0][c][0]);
  const double inverseMeanSquared = meanSquared > 0.0 ? 1.0 / meanSquared : 0.0;

  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    const std::size_t currentBlockLength = std::min(blockLength[block], numberOfBlockElements);
    for (std::size_t k = (block == 0 ? 0 : 1); k < currentBlockLength; ++k)
    {
      const std::size_t n = acfCount[block][c][k];
      if (n == 0) continue;
      const double value = acf[block][c][k] / static_cast<double>(n);
      data.push_back(EndToEndAutoCorrelationFunctionData{lagOf(block, k), value, value * inverseMeanSquared,
                                                         static_cast<double>(n)});
    }
  }
  // the blocks are already in increasing order of lag (block b+1 starts at n^(b+1) > (n-1) n^b)
  std::ranges::sort(data, {}, &EndToEndAutoCorrelationFunctionData::time);
  return data;
}

EndToEndRelaxationTimes PropertyEndToEndAutoCorrelationFunction::relaxationTimes(std::size_t component) const
{
  return relaxationTimes(result(component));
}

EndToEndRelaxationTimes PropertyEndToEndAutoCorrelationFunction::relaxationTimes(
    const std::vector<EndToEndAutoCorrelationFunctionData> &data)
{
  EndToEndRelaxationTimes times{};
  if (data.empty()) return times;

  times.meanSquaredEndToEnd = data.front().acf;
  times.longestLag = data.back().time;
  if (!(data.front().acf > 0.0) || data.size() < 2) return times;

  // integrated correlation time up to the first zero crossing (trapezoid rule on the normalized function)
  {
    double integral = 0.0;
    bool crossed = false;
    for (std::size_t i = 1; i < data.size(); ++i)
    {
      const double t0 = data[i - 1].time, t1 = data[i].time;
      const double c0 = data[i - 1].normalized, c1 = data[i].normalized;
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
    for (std::size_t i = 1; i < data.size(); ++i)
    {
      const double c0 = data[i - 1].normalized, c1 = data[i].normalized;
      if (c1 <= threshold)
      {
        const double t0 = data[i - 1].time, t1 = data[i].time;
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
    for (const EndToEndAutoCorrelationFunctionData &point : data)
    {
      if (point.time <= 0.0 || point.normalized > 0.5) continue;
      if (point.normalized <= 0.05) break;
      const double y = std::log(point.normalized);
      sumT += point.time;
      sumY += y;
      sumTT += point.time * point.time;
      sumTY += point.time * y;
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
  if (!writeEvery.has_value()) return;
  if (currentCycle % writeEvery.value() != 0uz) return;
  if (count == 0uz) return;

  std::filesystem::create_directory("end_to_end_acf");

  const std::string unit = timeUnit();
  for (std::size_t c = 0; c < numberOfComponents; ++c)
  {
    if (!hasData(c)) continue;
    const std::vector<EndToEndAutoCorrelationFunctionData> data = result(c);
    if (data.empty()) continue;
    const EndToEndRelaxationTimes times = relaxationTimes(data);
    const std::array<std::size_t, 2> &ends = endToEndAtomsPerComponent[c].value();

    std::ofstream stream(std::format("end_to_end_acf/end_to_end_acf_{}.s{}.txt", components[c].name, systemId));

    stream << std::format("# end-to-end vector autocorrelation function <R(0).R(t)> of component '{}'\n",
                          components[c].name);
    stream << std::format("# end-to-end atoms: {} and {}; {} molecules; {} time origins sampled every {} cycles\n",
                          ends[0], ends[1], numberOfMoleculesPerComponent[c], count, sampleEvery);
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
    stream << std::format("# column 1: time [{}]\n", unit);
    stream << "# column 2: <R(0).R(t)> [A^2]\n";
    stream << "# column 3: <R(0).R(t)>/<R^2> [-]\n";
    stream << "# column 4: number of samples [-]\n";

    for (const EndToEndAutoCorrelationFunctionData &point : data)
    {
      stream << std::format("{} {} {} (count: {})\n", point.time, point.acf, point.normalized,
                            static_cast<std::size_t>(point.numberOfSamples));
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEndToEndAutoCorrelationFunction &p)
{
  archive << p.versionNumber;

  archive << p.numberOfMoleculesPerComponent;
  archive << p.endToEndAtomsPerComponent;
  archive << p.numberOfComponents;
  archive << p.numberOfParticles;
  archive << p.timeStep;
  archive << p.numberOfBlockElements;
  archive << p.sampleEvery;
  archive << p.writeEvery;

  archive << p.count;
  archive << p.numberOfBlocks;
  archive << p.maxNumberOfBlocks;
  archive << p.blockLength;
  archive << p.acfCount;
  archive << p.blockData;
  archive << p.acf;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEndToEndAutoCorrelationFunction &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyEndToEndAutoCorrelationFunction' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> p.numberOfMoleculesPerComponent;
  archive >> p.endToEndAtomsPerComponent;
  archive >> p.numberOfComponents;
  archive >> p.numberOfParticles;
  archive >> p.timeStep;
  archive >> p.numberOfBlockElements;
  archive >> p.sampleEvery;
  archive >> p.writeEvery;

  archive >> p.count;
  archive >> p.numberOfBlocks;
  archive >> p.maxNumberOfBlocks;
  archive >> p.blockLength;
  archive >> p.acfCount;
  archive >> p.blockData;
  archive >> p.acf;

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
