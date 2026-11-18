module;

export module skpowderdiffraction;

import std;

import int3;
import double3;
import double3x3;

import skscatteringfactors;

// The powder diffraction pattern of a crystal: where its reflections fall in the scattering angle, how strong
// each is, and the curve a diffractometer would draw of them.
//
// What it is given is one unit cell: the cell matrix, with the lattice vectors as its columns, and the atoms in
// it by fractional position and atomic number. Symmetry has been applied by then; nothing here expands an
// asymmetric unit, and nothing here needs a force field, a charge or a type.
//
// The structure factor of a reflection hkl is
//
//   F(hkl) = sum_j o_j (f_j(s) + f'_j + i f''_j) exp(-B_j s^2) exp(2 pi i (h x_j + k y_j + l z_j)),
//
// with s = sin(theta)/lambda = 1/(2d), o the occupancy and B the isotropic temperature factor. The reflections
// at one d-spacing fall on one angle and are one peak of the powder, whether symmetry makes them equivalent or
// they only happen to coincide (333 and 511 of a cubic cell), so they are gathered by d-spacing alone. The peak
// has the summed |F|^2 of its reflections times the Lorentz and polarization factor of its angle.
//
// This follows the powder diffraction of RASPA2 (spectra.c), and departs from it where that went wrong:
//
//   - the range of h, k and l is the length of the cell vector over the smallest d-spacing, which is the bound
//     for a cell of any shape; one over the reciprocal length is that only for a rectangular cell, and misses
//     reflections of an oblique one;
//   - a peak is drawn at its own angle and not at the nearest point of the grid;
//   - the two lines of a K-alpha doublet are both drawn, two to one;
//   - the Gaussian part of the pseudo-Voigt profile is normalized as the Gaussian profile is;
//   - the anomalous scattering is that of the anode in use rather than always that of copper;
//   - the peak table lists d-spacings, in Angstrom, where the other listed their reciprocals.
export struct SKPowderDiffraction
{
  using Radiation = SKScatteringFactors::Radiation;
  using Anode = SKScatteringFactors::Anode;

  enum class PeakShape : std::size_t
  {
    Gaussian = 0,
    Lorentzian = 1,
    PseudoVoigt = 2
  };

  struct Atom
  {
    double3 fractionalPosition{};
    std::size_t atomicNumber{0};
    double occupancy{1.0};
    double isotropicB{0.0};  // B = 8 pi^2 <u^2>, Angstrom^2
  };

  // One wavelength of the radiation, in Angstrom, and its share of the intensity.
  struct Line
  {
    double wavelength{1.5418};
    double weight{1.0};
  };

  struct Settings
  {
    Radiation radiation{Radiation::XRay};

    // The anode whose anomalous scattering is added to the X-ray scattering factors, or none: the corrections
    // are tabulated at the five K-alpha wavelengths only, so radiation of another wavelength goes without.
    std::optional<Anode> anomalousScattering{Anode::Copper};

    std::vector<Line> lines{Line{}};

    // The range and the step of the pattern, in degrees of two theta.
    double twoThetaMin{3.0};
    double twoThetaMax{50.0};
    double twoThetaStep{0.01};

    PeakShape peakShape{PeakShape::PseudoVoigt};

    // The full width at half maximum of a peak at Bragg angle theta, in degrees (Caglioti):
    // FWHM^2 = U tan^2(theta) + V tan(theta) + W. A peak of no width is drawn into the one point nearest to it.
    double u{0.01};
    double v{-0.001};
    double w{0.002};

    // The Lorentzian share of the pseudo-Voigt profile; the rest is Gaussian.
    double mixing{0.5};

    // X-rays of a tube with the given anode: the single K-alpha wavelength, or the two lines of the doublet.
    static Settings forAnode(Anode anode, bool doublet = false);

    // Radiation of a wavelength of the caller's choosing, as from a synchrotron, a reactor or a microscope.
    static Settings forWavelength(Radiation radiation, double wavelength);
  };

  struct Peak
  {
    int3 hkl{};                          // one of the reflections of the peak: the largest, h before k before l
    std::size_t multiplicity{0};         // how many reflections fall on it
    double dSpacing{0.0};                // Angstrom
    double twoTheta{0.0};                // degrees
    double wavelength{0.0};              // Angstrom, the line this peak is of
    double lorentzPolarization{0.0};
    double structureFactorSquared{0.0};  // |F|^2, the mean over the reflections of the peak
    double intensity{0.0};               // the strongest peak of the range is 100
  };

  // The peaks within the range, by rising angle. Those that symmetry or chance leaves without intensity are
  // not among them.
  std::vector<Peak> peaks;

  // The pattern: the angles, in degrees, and the intensity at each, the highest point being 100.
  std::vector<double> twoTheta;
  std::vector<double> intensity;

  // Throws std::invalid_argument for settings that describe no measurement, a cell without volume, or an atom
  // of an element the table of the radiation does not have.
  static SKPowderDiffraction compute(const double3x3& unitCell, std::span<const Atom> atoms,
                                     const Settings& settings);

  // The same with the settings as they stand above: copper K-alpha, 3 to 50 degrees, pseudo-Voigt peaks.
  static SKPowderDiffraction compute(const double3x3& unitCell, std::span<const Atom> atoms);
};
