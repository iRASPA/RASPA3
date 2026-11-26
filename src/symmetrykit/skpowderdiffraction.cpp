module;

module skpowderdiffraction;

import std;

import int3;
import double3;
import double3x3;

import skelement;
import skscatteringfactors;

namespace
{
constexpr double degreesToRadians = std::numbers::pi / 180.0;
constexpr double radiansToDegrees = 180.0 / std::numbers::pi;

// The atoms that scatter alike: one element at one temperature factor. The scattering factor is worked out once
// per reflection for each of these rather than once per atom.
struct Species
{
  std::size_t atomicNumber;
  double isotropicB;
  std::vector<double3> positions;
  std::vector<double> occupancies;
};

struct Reflection
{
  int3 hkl;
  double reciprocalLength;  // 1/d
  double structureFactorSquared;
  std::size_t count;  // 2 where the reflection stands for its Friedel mate as well
};

// The reflections of one d-spacing.
struct Shell
{
  int3 hkl;
  double reciprocalLength;
  double summedStructureFactorSquared;
  std::size_t multiplicity;
};

// A peak as the pattern is drawn from it, before anything is scaled.
struct Bragg
{
  double twoTheta;  // degrees
  double fullWidth;  // degrees
  double intensity;
};

std::string elementSymbol(std::size_t atomicNumber)
{
  if (atomicNumber < PredefinedElements::predefinedElements.size())
  {
    return PredefinedElements::predefinedElements[atomicNumber]._chemicalSymbol;
  }
  return std::format("of atomic number {}", atomicNumber);
}

// Whether a number is finite, read off its bits: the library is compiled with finite mathematics only, under which
// a comparison with infinity is whatever the optimiser makes of it.
// (Taken by reference and read as bytes: passed by value, the number is one the compiler holds to be finite, and
// the test was folded away.)
bool isFiniteNumber(const double& value)
{
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return (bits & 0x7ff0000000000000ull) != 0x7ff0000000000000ull;
}

void validate(const SKPowderDiffraction::Settings& settings)
{
  if (!(settings.twoThetaMin >= 0.0) || !(settings.twoThetaMax > settings.twoThetaMin) ||
      !(settings.twoThetaMax < 180.0))
  {
    throw std::invalid_argument(
        std::format("powder diffraction: the two-theta range [{}, {}] is not within 0 and 180 degrees",
                    settings.twoThetaMin, settings.twoThetaMax));
  }
  if (!(settings.twoThetaStep > 0.0))
  {
    throw std::invalid_argument(
        std::format("powder diffraction: the two-theta step {} is not positive", settings.twoThetaStep));
  }
  if (settings.lines.empty())
  {
    throw std::invalid_argument("powder diffraction: no wavelength given");
  }
  for (const SKPowderDiffraction::Line& line : settings.lines)
  {
    // (an infinite weight made the strongest peak infinite, and the scaling to 100 every intensity not a number)
    if (!(line.wavelength > 0.0) || !(line.weight >= 0.0) || !isFiniteNumber(line.wavelength) ||
        !isFiniteNumber(line.weight))
    {
      throw std::invalid_argument(std::format(
          "powder diffraction: a wavelength of {} Angstrom with weight {} is not a line of the radiation",
          line.wavelength, line.weight));
    }
  }
  if (!(settings.mixing >= 0.0) || !(settings.mixing <= 1.0))
  {
    throw std::invalid_argument(
        std::format("powder diffraction: the pseudo-Voigt mixing {} is not between 0 and 1", settings.mixing));
  }
}

std::vector<Species> gatherSpecies(std::span<const SKPowderDiffraction::Atom> atoms)
{
  std::vector<Species> species;
  for (const SKPowderDiffraction::Atom& atom : atoms)
  {
    auto found = std::ranges::find_if(species, [&](const Species& s)
                                      { return s.atomicNumber == atom.atomicNumber && s.isotropicB == atom.isotropicB; });
    if (found == species.end())
    {
      species.push_back(Species{atom.atomicNumber, atom.isotropicB, {}, {}});
      found = std::prev(species.end());
    }
    found->positions.push_back(atom.fractionalPosition);
    found->occupancies.push_back(atom.occupancy);
  }
  return species;
}

// Every reflection out to `maximumReciprocalLength`, with its |F|^2. Without anomalous scattering the scattering
// factors are real and a reflection and its Friedel mate are equally strong, so half of them are worked out and
// counted twice.
std::vector<Reflection> reflections(const double3x3& unitCell, const std::vector<Species>& species,
                                    const SKPowderDiffraction::Settings& settings, double maximumReciprocalLength)
{
  const double3x3 inverseCell = unitCell.inverse();

  const bool anomalous =
      settings.radiation == SKPowderDiffraction::Radiation::XRay && settings.anomalousScattering.has_value();

  std::vector<std::complex<double>> dispersion(species.size(), std::complex<double>(0.0, 0.0));
  for (std::size_t i = 0; i < species.size(); ++i)
  {
    // asked once here, so that an element the table lacks is reported whether or not a reflection is in range
    if (!SKScatteringFactors::scatteringFactor(settings.radiation, species[i].atomicNumber, 0.0).has_value())
    {
      throw std::invalid_argument(std::format("powder diffraction: no scattering factor for element {}",
                                              elementSymbol(species[i].atomicNumber)));
    }
    if (anomalous)
    {
      dispersion[i] = SKScatteringFactors::anomalous(species[i].atomicNumber, *settings.anomalousScattering);
    }
  }

  // h is the component of the reciprocal vector along a, so |h| is at most |G| |a|, in a cell of any shape
  const double3 a(unitCell.ax, unitCell.ay, unitCell.az);
  const double3 b(unitCell.bx, unitCell.by, unitCell.bz);
  const double3 c(unitCell.cx, unitCell.cy, unitCell.cz);
  const double reach = maximumReciprocalLength * (1.0 + 1.0e-9);
  // more indices than can be counted, let alone summed over: a wavelength of 1e-12 Angstrom overflowed the int
  constexpr double mostIndices = 1.0e4;
  if (!(reach * a.length() < mostIndices) || !(reach * b.length() < mostIndices) || !(reach * c.length() < mostIndices))
  {
    throw std::invalid_argument(
        "powder diffraction: the wavelength is too short for the cell (more than 10^4 Miller indices along an axis)");
  }
  const int hmax = static_cast<int>(std::floor(reach * a.length()));
  const int kmax = static_cast<int>(std::floor(reach * b.length()));
  const int lmax = static_cast<int>(std::floor(reach * c.length()));

  std::vector<Reflection> result;
  for (int h = anomalous ? -hmax : 0; h <= hmax; ++h)
  {
    for (int k = -kmax; k <= kmax; ++k)
    {
      for (int l = -lmax; l <= lmax; ++l)
      {
        if (h == 0 && k == 0 && l == 0) continue;
        // the half of the reflections whose first index that is not zero is positive
        if (!anomalous && h == 0 && (k < 0 || (k == 0 && l < 0))) continue;

        const double3 index(static_cast<double>(h), static_cast<double>(k), static_cast<double>(l));
        const double reciprocalLength = transposedMultiply(inverseCell, index).length();
        if (reciprocalLength > reach) continue;

        const double sSquared = 0.25 * reciprocalLength * reciprocalLength;

        std::complex<double> structureFactor(0.0, 0.0);
        for (std::size_t i = 0; i < species.size(); ++i)
        {
          double real = 0.0;
          double imaginary = 0.0;
          for (std::size_t j = 0; j < species[i].positions.size(); ++j)
          {
            const double phase = 2.0 * std::numbers::pi * double3::dot(index, species[i].positions[j]);
            real += species[i].occupancies[j] * std::cos(phase);
            imaginary += species[i].occupancies[j] * std::sin(phase);
          }
          const double factor =
              SKScatteringFactors::scatteringFactor(settings.radiation, species[i].atomicNumber, sSquared).value();
          structureFactor += (factor + dispersion[i]) * std::exp(-species[i].isotropicB * sSquared) *
                             std::complex<double>(real, imaginary);
        }

        result.push_back(Reflection{int3(h, k, l), reciprocalLength, std::norm(structureFactor),
                                    anomalous ? std::size_t{1} : std::size_t{2}});
      }
    }
  }

  // A reflection the symmetry (or the repetition of a supercell) extinguishes is no part of a peak: its |F|^2 is
  // zero, or the rounding left of sums that cancel. Left in, it was gathered with the allowed reflections of the same
  // d-spacing: in P2_13 the absent (3,0,0) with (2,2,1), in Fd-3m (6,0,0) with (4,4,2), in a 3x3x3 supercell (4,1,1)
  // with (3,3,0). The peak then carried the index of the absent reflection, its multiplicity as well, and a mean
  // |F|^2 diluted by it.
  double strongest = 0.0;
  for (const Reflection& reflection : result) strongest = std::max(strongest, reflection.structureFactorSquared);
  std::erase_if(result, [&](const Reflection& reflection)
                { return !(reflection.structureFactorSquared > 1.0e-12 * strongest); });
  return result;
}

bool laterIndex(const int3& left, const int3& right)
{
  return std::tie(left.x, left.y, left.z) > std::tie(right.x, right.y, right.z);
}

// The reflections gathered by d-spacing, by rising reciprocal length.
std::vector<Shell> shells(std::vector<Reflection> list)
{
  std::ranges::sort(list, {}, &Reflection::reciprocalLength);

  std::vector<Shell> result;
  for (const Reflection& reflection : list)
  {
    if (result.empty() ||
        reflection.reciprocalLength - result.back().reciprocalLength > 1.0e-8 * reflection.reciprocalLength)
    {
      result.push_back(Shell{reflection.hkl, reflection.reciprocalLength, 0.0, 0});
    }
    Shell& shell = result.back();
    if (laterIndex(reflection.hkl, shell.hkl)) shell.hkl = reflection.hkl;
    shell.summedStructureFactorSquared += static_cast<double>(reflection.count) * reflection.structureFactorSquared;
    shell.multiplicity += reflection.count;
  }
  return result;
}

// The profile of a peak of unit area and full width at half maximum `fullWidth`, `distance` from its centre.
double profile(SKPowderDiffraction::PeakShape shape, double mixing, double fullWidth, double distance)
{
  const double ratioSquared = (distance * distance) / (fullWidth * fullWidth);
  const double fourLogTwo = 4.0 * std::numbers::ln2;

  auto gaussian = [&]()
  {
    // beyond seven widths it is below 1e-58 of its height
    if (ratioSquared > 49.0) return 0.0;
    return std::sqrt(fourLogTwo / std::numbers::pi) * std::exp(-fourLogTwo * ratioSquared) / fullWidth;
  };
  auto lorentzian = [&]() { return 2.0 / (std::numbers::pi * fullWidth * (1.0 + 4.0 * ratioSquared)); };

  switch (shape)
  {
    case SKPowderDiffraction::PeakShape::Gaussian:
      return gaussian();
    case SKPowderDiffraction::PeakShape::Lorentzian:
      return lorentzian();
    case SKPowderDiffraction::PeakShape::PseudoVoigt:
      return mixing * lorentzian() + (1.0 - mixing) * gaussian();
  }
  return 0.0;
}
}  // namespace

SKPowderDiffraction::Settings SKPowderDiffraction::Settings::forAnode(Anode anode, bool doublet)
{
  const SKScatteringFactors::Wavelengths wavelengths = SKScatteringFactors::wavelengths(anode);

  Settings settings;
  settings.radiation = Radiation::XRay;
  settings.anomalousScattering = anode;
  if (doublet)
  {
    settings.lines = {Line{wavelengths.alpha1, 2.0 / 3.0}, Line{wavelengths.alpha2, 1.0 / 3.0}};
  }
  else
  {
    settings.lines = {Line{wavelengths.mean, 1.0}};
  }
  return settings;
}

SKPowderDiffraction::Settings SKPowderDiffraction::Settings::forWavelength(Radiation radiation, double wavelength)
{
  Settings settings;
  settings.radiation = radiation;
  settings.anomalousScattering = std::nullopt;
  settings.lines = {Line{wavelength, 1.0}};
  return settings;
}

SKPowderDiffraction SKPowderDiffraction::compute(const double3x3& unitCell, std::span<const Atom> atoms)
{
  return compute(unitCell, atoms, Settings{});
}

SKPowderDiffraction SKPowderDiffraction::compute(const double3x3& unitCell, std::span<const Atom> atoms,
                                                 const Settings& settings)
{
  validate(settings);
  if (!(std::abs(unitCell.determinant()) > 0.0))
  {
    throw std::invalid_argument("powder diffraction: the unit cell has no volume");
  }

  SKPowderDiffraction result;

  if (!((settings.twoThetaMax - settings.twoThetaMin) / settings.twoThetaStep < 1.0e8))
  {
    throw std::invalid_argument(std::format("powder diffraction: a two-theta step of {} gives more than 10^8 points",
                                            settings.twoThetaStep));
  }
  const std::size_t numberOfPoints =
      static_cast<std::size_t>(std::floor((settings.twoThetaMax - settings.twoThetaMin) / settings.twoThetaStep + 1.0e-9)) +
      1;
  result.twoTheta.resize(numberOfPoints);
  result.intensity.assign(numberOfPoints, 0.0);
  for (std::size_t i = 0; i < numberOfPoints; ++i)
  {
    result.twoTheta[i] = settings.twoThetaMin + static_cast<double>(i) * settings.twoThetaStep;
  }

  // A peak just past the end of the range still has a flank inside it, so the reflections are taken a degree
  // further than the range and drawn; the table lists the range alone.
  const double margin = 1.0;
  const double twoThetaReach = std::min(settings.twoThetaMax + margin, 179.0);
  const double twoThetaLimit = std::max(twoThetaReach, settings.twoThetaMax);
  const double shortestWavelength = std::ranges::min(settings.lines, {}, &Line::wavelength).wavelength;
  const double maximumReciprocalLength = 2.0 * std::sin(0.5 * twoThetaLimit * degreesToRadians) / shortestWavelength;

  const std::vector<Species> species = gatherSpecies(atoms);
  const std::vector<Shell> gathered = shells(reflections(unitCell, species, settings, maximumReciprocalLength));

  std::vector<Bragg> drawn;
  for (const Line& line : settings.lines)
  {
    for (const Shell& shell : gathered)
    {
      const double sinTheta = 0.5 * line.wavelength * shell.reciprocalLength;
      if (sinTheta >= 1.0) continue;
      const double theta = std::asin(sinTheta);
      const double twoTheta = 2.0 * theta * radiansToDegrees;
      if (twoTheta > twoThetaLimit) continue;
      // the Lorentz factor has no value in the forward direction, and no detector stands there
      if (!(sinTheta > 0.0)) continue;

      const double cosTheta = std::cos(theta);
      const double lorentzPolarization =
          settings.radiation == Radiation::XRay
              ? (1.0 + std::cos(2.0 * theta) * std::cos(2.0 * theta)) / (sinTheta * sinTheta * cosTheta)
              : 1.0 / (std::sin(2.0 * theta) * sinTheta);
      const double peakIntensity = line.weight * lorentzPolarization * shell.summedStructureFactorSquared;

      const double tanTheta = sinTheta / cosTheta;
      const double fullWidth =
          std::sqrt(std::max(0.0, settings.u * tanTheta * tanTheta + settings.v * tanTheta + settings.w));
      drawn.push_back(Bragg{twoTheta, fullWidth, peakIntensity});

      if (twoTheta >= settings.twoThetaMin - 1.0e-9 && twoTheta <= settings.twoThetaMax + 1.0e-9)
      {
        result.peaks.push_back(Peak{
            .hkl = shell.hkl,
            .multiplicity = shell.multiplicity,
            .dSpacing = 1.0 / shell.reciprocalLength,
            .twoTheta = twoTheta,
            .wavelength = line.wavelength,
            .lorentzPolarization = lorentzPolarization,
            .structureFactorSquared = shell.summedStructureFactorSquared / static_cast<double>(shell.multiplicity),
            .intensity = peakIntensity});
      }
    }
  }

  // the table: the strongest peak of the range is 100, and what is left without intensity is not a peak
  double strongest = 0.0;
  for (const Peak& peak : result.peaks) strongest = std::max(strongest, peak.intensity);
  // (every reflection that is drawn, also those past the range whose tails reach into it: one that is not finite made
  // the whole pattern not a number)
  if (!isFiniteNumber(strongest) ||
      !std::ranges::all_of(drawn, [](const Bragg& bragg) { return isFiniteNumber(bragg.intensity); }))
  {
    // exp(-B s^2) of a large negative B overflows, and the scaling below made every intensity not a number
    throw std::invalid_argument("powder diffraction: the intensities are not finite (check the displacement factors)");
  }
  if (strongest > 0.0)
  {
    for (Peak& peak : result.peaks) peak.intensity *= 100.0 / strongest;
  }
  std::erase_if(result.peaks, [](const Peak& peak) { return !(peak.intensity > 1.0e-10); });
  std::ranges::stable_sort(result.peaks, {}, &Peak::twoTheta);

  // the pattern
  for (const Bragg& bragg : drawn)
  {
    if (!(bragg.intensity > 0.0)) continue;

    if (bragg.fullWidth < 1.0e-5)
    {
      // a line only inside the range, as in the table of peaks: one within half a step outside it was rounded onto
      // the first or the last point, a spike in the pattern with no peak in the table
      const double position = std::round((bragg.twoTheta - settings.twoThetaMin) / settings.twoThetaStep);
      if (bragg.twoTheta >= settings.twoThetaMin - 1.0e-9 && bragg.twoTheta <= settings.twoThetaMax + 1.0e-9 &&
          position >= 0.0 && position < static_cast<double>(numberOfPoints))
      {
        result.intensity[static_cast<std::size_t>(position)] += bragg.intensity / settings.twoThetaStep;
      }
      continue;
    }

    for (std::size_t i = 0; i < numberOfPoints; ++i)
    {
      result.intensity[i] += bragg.intensity * profile(settings.peakShape, settings.mixing, bragg.fullWidth,
                                                       result.twoTheta[i] - bragg.twoTheta);
    }
  }

  const double highest = result.intensity.empty() ? 0.0 : std::ranges::max(result.intensity);
  if (highest > 0.0)
  {
    for (double& value : result.intensity) value *= 100.0 / highest;
  }

  return result;
}
