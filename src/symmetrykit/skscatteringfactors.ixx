module;

export module skscatteringfactors;

import std;

// What one atom of an element scatters, for the three kinds of radiation a powder pattern is measured with.
//
// The tables are those of RASPA2 (scattering_factors.c), indexed here by atomic number so that they sit beside
// the element table of `skelement` and are read with the same index:
//
//   X-ray:     D. Waasmaier and A. Kirfel, Acta Cryst. A51, 416-431 (1995), neutral atoms, five Gaussians and a
//              constant in s = sin(theta)/lambda, valid up to s = 6 1/Angstrom;
//   anomalous: D. T. Cromer and J. A. Ibers, International Tables for X-ray Crystallography vol. 4 (1974), f' and
//              f'' at the K-alpha radiation of the five common anodes;
//   electron:  P. A. Doyle and P. S. Turner, Acta Cryst. A24, 390-397 (1968), four Gaussians;
//   neutron:   bound coherent scattering lengths, in fm, which do not depend on the angle.
//
// The electron and the neutron table do not cover every element. An element they do not have is reported as
// nothing rather than as an atom that scatters nothing, which is what a row of zeros would have been read as.
export struct SKScatteringFactors
{
  enum class Radiation : std::size_t
  {
    XRay = 0,
    Neutron = 1,
    Electron = 2
  };

  // The anode of an X-ray tube, in the order of the columns of the anomalous-scattering table.
  enum class Anode : std::size_t
  {
    Chromium = 0,
    Iron = 1,
    Copper = 2,
    Molybdenum = 3,
    Silver = 4
  };

  // The K-alpha radiation of an anode, in Angstrom: the two lines of the doublet, and the single wavelength
  // that stands for both where they are not resolved (their mean, weighted two to one).
  struct Wavelengths
  {
    double mean;
    double alpha1;
    double alpha2;
  };

  // The highest atomic number the tables reach (Cf).
  static constexpr std::size_t maximumAtomicNumber = 98;

  static Wavelengths wavelengths(Anode anode);

  // The scattering factor of a neutral atom at s^2 = (sin(theta)/lambda)^2, in 1/Angstrom^2: in electrons for
  // X-rays, in Angstrom for electrons, and in fm, whatever the angle, for neutrons. Nothing for an atomic
  // number that is not an element of the table asked.
  static std::optional<double> xray(std::size_t atomicNumber, double sSquared);
  static std::optional<double> electron(std::size_t atomicNumber, double sSquared);
  static std::optional<double> neutron(std::size_t atomicNumber);
  static std::optional<double> scatteringFactor(Radiation radiation, std::size_t atomicNumber, double sSquared);

  // The dispersion correction f' + i f'' to the X-ray scattering factor at the K-alpha radiation of an anode.
  // Zero for an atomic number outside the table.
  static std::complex<double> anomalous(std::size_t atomicNumber, Anode anode);
};
