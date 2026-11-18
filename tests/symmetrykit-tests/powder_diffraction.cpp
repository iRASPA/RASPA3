#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;

import skscatteringfactors;
import skpowderdiffraction;

// The powder pattern of cells small enough that every peak can be worked out by hand: where the peaks of a cubic
// cell fall and how many reflections each gathers, which of them a centring removes, and what the profile drawn
// of a single one looks like.

namespace
{
using Atom = SKPowderDiffraction::Atom;
using Settings = SKPowderDiffraction::Settings;
using Radiation = SKPowderDiffraction::Radiation;
using Anode = SKPowderDiffraction::Anode;
using PeakShape = SKPowderDiffraction::PeakShape;

constexpr double degrees = std::numbers::pi / 180.0;

// Build a cell (columns are lattice vectors) from cell parameters (angles in degrees).
double3x3 cell(double a, double b, double c, double alpha, double beta, double gamma)
{
  double ca = std::cos(alpha * degrees);
  double cb = std::cos(beta * degrees);
  double cg = std::cos(gamma * degrees);
  double sg = std::sin(gamma * degrees);
  double vf = std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
  return double3x3(double3(a, 0.0, 0.0), double3(b * cg, b * sg, 0.0),
                   double3(c * cb, c * (ca - cb * cg) / sg, c * vf / sg));
}

double3x3 cubic(double a) { return double3x3(double3(a, 0.0, 0.0), double3(0.0, a, 0.0), double3(0.0, 0.0, a)); }

// h^2 + k^2 + l^2 of each peak, read back off its d-spacing.
std::vector<int> squaredIndices(const SKPowderDiffraction& pattern, double a)
{
  std::vector<int> result;
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks)
  {
    result.push_back(static_cast<int>(std::lround((a * a) / (peak.dSpacing * peak.dSpacing))));
  }
  return result;
}

Settings copperTo90()
{
  Settings settings = Settings::forWavelength(Radiation::XRay, 1.5418);
  settings.twoThetaMin = 5.0;
  settings.twoThetaMax = 90.0;
  settings.twoThetaStep = 0.02;
  return settings;
}
}  // namespace

TEST(scattering_factors, xray_factor_at_zero_angle_is_the_number_of_electrons)
{
  for (std::size_t atomicNumber = 1; atomicNumber <= SKScatteringFactors::maximumAtomicNumber; ++atomicNumber)
  {
    const std::optional<double> factor = SKScatteringFactors::xray(atomicNumber, 0.0);
    ASSERT_TRUE(factor.has_value());
    EXPECT_NEAR(*factor, static_cast<double>(atomicNumber), 0.05) << "atomic number " << atomicNumber;
  }
}

TEST(scattering_factors, elements_outside_a_table_are_nothing)
{
  EXPECT_FALSE(SKScatteringFactors::xray(0, 0.0).has_value());
  EXPECT_FALSE(SKScatteringFactors::xray(SKScatteringFactors::maximumAtomicNumber + 1, 0.0).has_value());
  // Doyle and Turner have no hydrogen, and there is no scattering length for polonium
  EXPECT_FALSE(SKScatteringFactors::electron(1, 0.0).has_value());
  EXPECT_FALSE(SKScatteringFactors::neutron(84).has_value());

  EXPECT_NEAR(SKScatteringFactors::neutron(1).value(), -3.7390, 1e-12);
  EXPECT_NEAR(SKScatteringFactors::neutron(14).value(), 4.1491, 1e-12);
  EXPECT_NEAR(SKScatteringFactors::electron(6, 0.0).value(), 0.7307 + 1.1951 + 0.4563 + 0.1247, 1e-12);
}

TEST(scattering_factors, anomalous_scattering_follows_the_anode)
{
  const std::complex<double> copper = SKScatteringFactors::anomalous(14, Anode::Copper);
  EXPECT_NEAR(copper.real(), 0.244, 1e-12);
  EXPECT_NEAR(copper.imag(), 0.330, 1e-12);
  const std::complex<double> molybdenum = SKScatteringFactors::anomalous(14, Anode::Molybdenum);
  EXPECT_NEAR(molybdenum.real(), 0.072, 1e-12);
  EXPECT_NEAR(molybdenum.imag(), 0.071, 1e-12);

  EXPECT_NEAR(SKScatteringFactors::wavelengths(Anode::Copper).mean, 1.5418, 1e-12);
  EXPECT_NEAR(SKScatteringFactors::wavelengths(Anode::Copper).alpha1, 1.54056, 1e-12);
  EXPECT_NEAR(SKScatteringFactors::wavelengths(Anode::Copper).alpha2, 1.54439, 1e-12);
}

TEST(powder_diffraction, simple_cubic_peaks_and_multiplicities)
{
  const double a = 4.0;
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 29}};
  const Settings settings = copperTo90();
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), atoms, settings);

  // every sum of three squares up to 13: 7 is not one
  EXPECT_EQ(squaredIndices(pattern, a), (std::vector<int>{1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13}));

  // 300 and 221 share a d-spacing, as do 310 alone, 311 alone, 222 alone and 320 alone
  const std::vector<std::size_t> multiplicities{6, 12, 8, 6, 24, 24, 12, 30, 24, 24, 8, 24};
  const std::vector<int3> indices{int3(1, 0, 0), int3(1, 1, 0), int3(1, 1, 1), int3(2, 0, 0),
                                  int3(2, 1, 0), int3(2, 1, 1), int3(2, 2, 0), int3(3, 0, 0),
                                  int3(3, 1, 0), int3(3, 1, 1), int3(2, 2, 2), int3(3, 2, 0)};
  ASSERT_EQ(pattern.peaks.size(), multiplicities.size());
  for (std::size_t i = 0; i < pattern.peaks.size(); ++i)
  {
    const SKPowderDiffraction::Peak& peak = pattern.peaks[i];
    EXPECT_EQ(peak.multiplicity, multiplicities[i]) << "peak " << i;
    EXPECT_EQ(peak.hkl, indices[i]) << "peak " << i;
    // Bragg's law
    EXPECT_NEAR(peak.twoTheta, 2.0 * std::asin(0.5 * 1.5418 / peak.dSpacing) / degrees, 1e-10);
    // one atom at the origin: |F|^2 is its scattering factor squared
    const double factor = SKScatteringFactors::xray(29, 0.25 / (peak.dSpacing * peak.dSpacing)).value();
    EXPECT_NEAR(peak.structureFactorSquared, factor * factor, 1e-9 * factor * factor);
  }

  double strongest = 0.0;
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks) strongest = std::max(strongest, peak.intensity);
  EXPECT_NEAR(strongest, 100.0, 1e-10);
  EXPECT_NEAR(std::ranges::max(pattern.intensity), 100.0, 1e-10);
  EXPECT_EQ(pattern.twoTheta.size(), pattern.intensity.size());
  EXPECT_EQ(pattern.twoTheta.size(), 4251uz);
  EXPECT_NEAR(pattern.twoTheta.front(), 5.0, 1e-12);
  EXPECT_NEAR(pattern.twoTheta.back(), 90.0, 1e-9);
}

TEST(powder_diffraction, centring_removes_reflections)
{
  const double a = 4.0;
  const Settings settings = copperTo90();

  // body-centred: h + k + l even
  const std::vector<Atom> body{Atom{double3(0.0, 0.0, 0.0), 26}, Atom{double3(0.5, 0.5, 0.5), 26}};
  EXPECT_EQ(squaredIndices(SKPowderDiffraction::compute(cubic(a), body, settings), a),
            (std::vector<int>{2, 4, 6, 8, 10, 12}));

  // face-centred: h, k and l all even or all odd
  const std::vector<Atom> face{Atom{double3(0.0, 0.0, 0.0), 29}, Atom{double3(0.5, 0.5, 0.0), 29},
                               Atom{double3(0.5, 0.0, 0.5), 29}, Atom{double3(0.0, 0.5, 0.5), 29}};
  EXPECT_EQ(squaredIndices(SKPowderDiffraction::compute(cubic(a), face, settings), a),
            (std::vector<int>{3, 4, 8, 11, 12}));
}

TEST(powder_diffraction, two_elements_scatter_with_their_difference)
{
  // caesium chloride: the reflections a body centring would remove are left with the difference of the two
  const double a = 4.123;
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 55}, Atom{double3(0.5, 0.5, 0.5), 17}};
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), atoms, copperTo90());

  ASSERT_GE(pattern.peaks.size(), 2uz);
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks)
  {
    const double sSquared = 0.25 / (peak.dSpacing * peak.dSpacing);
    const double caesium = SKScatteringFactors::xray(55, sSquared).value();
    const double chlorine = SKScatteringFactors::xray(17, sSquared).value();
    const int sum = peak.hkl.x + peak.hkl.y + peak.hkl.z;
    // 300 and 221 are both odd, so the peak they share is one or the other throughout
    const double expected = (sum % 2 == 0) ? (caesium + chlorine) : (caesium - chlorine);
    EXPECT_NEAR(peak.structureFactorSquared, expected * expected, 1e-9 * expected * expected);
  }
}

TEST(powder_diffraction, temperature_factor_and_occupancy_scale_the_structure_factor)
{
  const double a = 4.0;
  const double B = 1.7;
  const std::vector<Atom> plain{Atom{double3(0.0, 0.0, 0.0), 29}};
  const std::vector<Atom> damped{Atom{double3(0.0, 0.0, 0.0), 29, 0.5, B}};
  const SKPowderDiffraction reference = SKPowderDiffraction::compute(cubic(a), plain, copperTo90());
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), damped, copperTo90());

  ASSERT_EQ(pattern.peaks.size(), reference.peaks.size());
  for (std::size_t i = 0; i < pattern.peaks.size(); ++i)
  {
    const double sSquared = 0.25 / (pattern.peaks[i].dSpacing * pattern.peaks[i].dSpacing);
    const double expected = 0.25 * std::exp(-2.0 * B * sSquared) * reference.peaks[i].structureFactorSquared;
    EXPECT_NEAR(pattern.peaks[i].structureFactorSquared, expected, 1e-9 * expected);
  }
}

TEST(powder_diffraction, oblique_cell_has_every_reflection_in_range)
{
  // the index range of an oblique cell is wider than one over its reciprocal lengths
  const double3x3 unitCell = cell(5.0, 6.0, 7.0, 70.0, 80.0, 65.0);
  const double3x3 inverseCell = unitCell.inverse();
  const std::vector<Atom> atoms{Atom{double3(0.1, 0.2, 0.3), 14}};
  Settings settings = Settings::forWavelength(Radiation::XRay, 1.5418);
  settings.twoThetaMin = 1.0;
  settings.twoThetaMax = 70.0;
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(unitCell, atoms, settings);

  std::size_t expected = 0;
  for (int h = -20; h <= 20; ++h)
    for (int k = -20; k <= 20; ++k)
      for (int l = -20; l <= 20; ++l)
      {
        if (h == 0 && k == 0 && l == 0) continue;
        const double reciprocalLength = transposedMultiply(inverseCell, double3(h, k, l)).length();
        const double sinTheta = 0.5 * 1.5418 * reciprocalLength;
        if (sinTheta >= 1.0) continue;
        const double twoTheta = 2.0 * std::asin(sinTheta) / degrees;
        if (twoTheta >= 1.0 && twoTheta <= 70.0) ++expected;
      }

  std::size_t found = 0;
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks) found += peak.multiplicity;
  EXPECT_EQ(found, expected);
  // a triclinic cell has no reflections that coincide but a reflection and its Friedel mate
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks) EXPECT_EQ(peak.multiplicity, 2uz);
}

TEST(powder_diffraction, pattern_does_not_depend_on_the_origin_or_on_the_choice_of_cell)
{
  const double3x3 unitCell = cell(5.0, 6.0, 7.0, 90.0, 105.0, 90.0);
  const std::vector<Atom> atoms{Atom{double3(0.10, 0.20, 0.30), 14}, Atom{double3(0.55, 0.15, 0.80), 8},
                                Atom{double3(0.35, 0.70, 0.45), 30}};
  const Settings settings = Settings::forAnode(Anode::Copper);
  const SKPowderDiffraction reference = SKPowderDiffraction::compute(unitCell, atoms, settings);

  std::vector<Atom> shifted = atoms;
  for (Atom& atom : shifted) atom.fractionalPosition = atom.fractionalPosition + double3(0.13, 0.27, 0.41);
  const SKPowderDiffraction moved = SKPowderDiffraction::compute(unitCell, shifted, settings);

  // the same crystal described by a cell twice as long along a
  const double3x3 doubled(double3(2.0 * unitCell.ax, 2.0 * unitCell.ay, 2.0 * unitCell.az),
                          double3(unitCell.bx, unitCell.by, unitCell.bz),
                          double3(unitCell.cx, unitCell.cy, unitCell.cz));
  std::vector<Atom> repeated;
  for (const Atom& atom : atoms)
  {
    for (double offset : {0.0, 0.5})
    {
      Atom copy = atom;
      copy.fractionalPosition.x = 0.5 * atom.fractionalPosition.x + offset;
      repeated.push_back(copy);
    }
  }
  const SKPowderDiffraction supercell = SKPowderDiffraction::compute(doubled, repeated, settings);

  ASSERT_EQ(moved.peaks.size(), reference.peaks.size());
  ASSERT_EQ(supercell.peaks.size(), reference.peaks.size());
  for (std::size_t i = 0; i < reference.peaks.size(); ++i)
  {
    EXPECT_NEAR(moved.peaks[i].twoTheta, reference.peaks[i].twoTheta, 1e-9);
    EXPECT_NEAR(moved.peaks[i].intensity, reference.peaks[i].intensity, 1e-8);
    EXPECT_NEAR(supercell.peaks[i].twoTheta, reference.peaks[i].twoTheta, 1e-9);
    EXPECT_NEAR(supercell.peaks[i].intensity, reference.peaks[i].intensity, 1e-8);
  }
  ASSERT_EQ(supercell.intensity.size(), reference.intensity.size());
  for (std::size_t i = 0; i < reference.intensity.size(); ++i)
  {
    EXPECT_NEAR(moved.intensity[i], reference.intensity[i], 1e-8);
    EXPECT_NEAR(supercell.intensity[i], reference.intensity[i], 1e-8);
  }
}

TEST(powder_diffraction, profile_of_a_single_peak)
{
  // a cubic cell whose 100 reflection falls on 22 degrees exactly, and nothing else within a degree of the range
  const double wavelength = 1.5418;
  const double a = wavelength / (2.0 * std::sin(11.0 * degrees));
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 29}};
  const double fullWidth = 0.2;

  Settings settings = Settings::forWavelength(Radiation::XRay, wavelength);
  settings.twoThetaMin = 20.0;
  settings.twoThetaMax = 24.0;
  settings.twoThetaStep = 0.001;
  settings.u = 0.0;
  settings.v = 0.0;
  settings.w = fullWidth * fullWidth;
  settings.mixing = 0.3;

  const std::size_t centre = 2000;
  const std::size_t halfWidth = 100;  // half the full width, in steps
  const double gaussianHeight = std::sqrt(4.0 * std::numbers::ln2 / std::numbers::pi);
  const double lorentzianHeight = 2.0 / std::numbers::pi;

  // the height one full width from the centre, over the height at the centre
  const std::map<PeakShape, double> expected{
      {PeakShape::Gaussian, 1.0 / 16.0},
      {PeakShape::Lorentzian, 1.0 / 5.0},
      {PeakShape::PseudoVoigt, (0.3 * lorentzianHeight / 5.0 + 0.7 * gaussianHeight / 16.0) /
                                   (0.3 * lorentzianHeight + 0.7 * gaussianHeight)}};

  for (const auto& [shape, ratio] : expected)
  {
    settings.peakShape = shape;
    const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), atoms, settings);

    ASSERT_EQ(pattern.peaks.size(), 1uz);
    EXPECT_NEAR(pattern.peaks[0].twoTheta, 22.0, 1e-9);
    ASSERT_EQ(pattern.intensity.size(), 4001uz);
    EXPECT_NEAR(pattern.intensity[centre], 100.0, 1e-9);
    EXPECT_NEAR(pattern.intensity[centre - halfWidth], 50.0, 1e-6);
    EXPECT_NEAR(pattern.intensity[centre + halfWidth], 50.0, 1e-6);
    EXPECT_NEAR(pattern.intensity[centre + 2 * halfWidth], 100.0 * ratio, 1e-6);
  }
}

TEST(powder_diffraction, peak_is_drawn_at_its_own_angle)
{
  // the peak falls between two points of a coarse grid; the centre of what is drawn is the peak and not a point
  const double wavelength = 1.5418;
  const double a = wavelength / (2.0 * std::sin(0.5 * 22.013 * degrees));
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 29}};

  Settings settings = Settings::forWavelength(Radiation::XRay, wavelength);
  settings.twoThetaMin = 20.0;
  settings.twoThetaMax = 24.0;
  settings.twoThetaStep = 0.05;
  settings.peakShape = PeakShape::Gaussian;
  settings.u = 0.0;
  settings.v = 0.0;
  settings.w = 0.04;
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), atoms, settings);

  double sum = 0.0;
  double moment = 0.0;
  for (std::size_t i = 0; i < pattern.intensity.size(); ++i)
  {
    sum += pattern.intensity[i];
    moment += pattern.intensity[i] * pattern.twoTheta[i];
  }
  EXPECT_NEAR(moment / sum, 22.013, 1e-6);
}

TEST(powder_diffraction, doublet_draws_both_lines_two_to_one)
{
  const double a = 4.0;
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 29}};
  const SKPowderDiffraction single = SKPowderDiffraction::compute(cubic(a), atoms, Settings::forAnode(Anode::Copper));
  const SKPowderDiffraction doublet =
      SKPowderDiffraction::compute(cubic(a), atoms, Settings::forAnode(Anode::Copper, true));

  ASSERT_EQ(doublet.peaks.size(), 2 * single.peaks.size());
  for (std::size_t i = 0; i < single.peaks.size(); ++i)
  {
    const SKPowderDiffraction::Peak& first = doublet.peaks[2 * i];
    const SKPowderDiffraction::Peak& second = doublet.peaks[2 * i + 1];
    EXPECT_NEAR(first.wavelength, 1.54056, 1e-12);
    EXPECT_NEAR(second.wavelength, 1.54439, 1e-12);
    EXPECT_NEAR(first.dSpacing, second.dSpacing, 1e-12);
    EXPECT_NEAR(first.twoTheta, 2.0 * std::asin(0.5 * 1.54056 / first.dSpacing) / degrees, 1e-10);
    EXPECT_NEAR(second.twoTheta, 2.0 * std::asin(0.5 * 1.54439 / second.dSpacing) / degrees, 1e-10);
    EXPECT_NEAR((first.intensity / first.lorentzPolarization) / (second.intensity / second.lorentzPolarization), 2.0,
                1e-9);
  }
}

TEST(powder_diffraction, neutrons_see_the_scattering_lengths)
{
  // a hydride: hydrogen scatters neutrons in antiphase, and hardly scatters X-rays at all
  const double a = 4.0;
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 3}, Atom{double3(0.5, 0.5, 0.5), 1}};
  Settings settings = Settings::forWavelength(Radiation::Neutron, 1.5418);
  settings.twoThetaMin = 5.0;
  settings.twoThetaMax = 60.0;
  const SKPowderDiffraction pattern = SKPowderDiffraction::compute(cubic(a), atoms, settings);

  const double lithium = SKScatteringFactors::neutron(3).value();
  const double hydrogen = SKScatteringFactors::neutron(1).value();
  ASSERT_FALSE(pattern.peaks.empty());
  for (const SKPowderDiffraction::Peak& peak : pattern.peaks)
  {
    const int sum = peak.hkl.x + peak.hkl.y + peak.hkl.z;
    const double expected = (sum % 2 == 0) ? (lithium + hydrogen) : (lithium - hydrogen);
    EXPECT_NEAR(peak.structureFactorSquared, expected * expected, 1e-9 * expected * expected);
    const double theta = 0.5 * peak.twoTheta * degrees;
    EXPECT_NEAR(peak.lorentzPolarization, 1.0 / (std::sin(2.0 * theta) * std::sin(theta)), 1e-9);
  }
}

TEST(powder_diffraction, what_cannot_be_measured_is_refused)
{
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 1}};

  // no electron scattering factor for hydrogen
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), atoms, Settings::forWavelength(Radiation::Electron, 0.025)),
               std::invalid_argument);
  // not an element
  const std::vector<Atom> unknown{Atom{double3(0.0, 0.0, 0.0), 0}};
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), unknown), std::invalid_argument);

  Settings settings;
  settings.twoThetaMax = 180.0;
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), atoms, settings), std::invalid_argument);
  settings = Settings{};
  settings.twoThetaStep = 0.0;
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), atoms, settings), std::invalid_argument);
  settings = Settings{};
  settings.lines.clear();
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), atoms, settings), std::invalid_argument);
  settings = Settings{};
  settings.lines = {SKPowderDiffraction::Line{0.0, 1.0}};
  EXPECT_THROW(SKPowderDiffraction::compute(cubic(4.0), atoms, settings), std::invalid_argument);

  EXPECT_THROW(SKPowderDiffraction::compute(double3x3(double3(1.0, 0.0, 0.0), double3(2.0, 0.0, 0.0),
                                                      double3(0.0, 0.0, 1.0)),
                                            atoms),
               std::invalid_argument);

  // an empty cell scatters nothing, and that is a pattern
  const SKPowderDiffraction empty = SKPowderDiffraction::compute(cubic(4.0), std::span<const Atom>{});
  EXPECT_TRUE(empty.peaks.empty());
  EXPECT_EQ(empty.intensity.size(), 4701uz);
  EXPECT_EQ(std::ranges::max(empty.intensity), 0.0);
}
