module;

export module lammps_styles;

import std;

import double3;
import forcefield;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import torsion_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bend_torsion_potential;

/**
 * \brief RASPA <-> LAMMPS functional-form mapping.
 *
 * Pure functions: one RASPA intramolecular term in, one or more LAMMPS coefficient lines out (in LAMMPS
 * 'real' units: kcal/mol, Angstrom, degrees), and the inverse for the styles a LAMMPS data file can carry.
 * Nothing here knows about molecules or type numbering; that is the topology builder's job.
 *
 * Conventions that the mapping relies on:
 *   - RASPA and LAMMPS use the same dihedral angle (phi = 180 degrees for trans); both compute the angle
 *     between the ABC and BCD planes, so every function of cos(n phi) maps exactly.
 *   - RASPA's harmonic forms carry a 1/2 that LAMMPS's harmonic styles do not: E_LAMMPS = K (x - x0)^2.
 *   - 'dihedral_style nharmonic' is E = sum_{i=1}^{n} A_i cos^(i-1)(phi); every RASPA cosine series is
 *     rewritten as a polynomial in cos(phi) through the Chebyshev identities cos(n phi) = T_n(cos phi).
 *   - Forms without an analytic LAMMPS equivalent are exported as 'table' styles; the topology builder
 *     tabulates them from RASPA's own energy functions (see the *TableBody functions below).
 */
export namespace LAMMPS
{
/// One LAMMPS coefficient line: the style it belongs to and its arguments, already formatted in real units.
struct Term
{
  std::string style{};         ///< LAMMPS sub-style name ("harmonic", "nharmonic", "table", "zero", ...).
  std::string coefficients{};  ///< Formatted coefficient string, without the type index.
  std::string tableKeyword{};  ///< For style "table": the keyword of the table section (its body is built later).
  bool exact{true};            ///< False when the LAMMPS form only approximates the RASPA form (or drops a constant).
  std::string note{};          ///< Human-readable remark written as a trailing comment.

  /// Optional class2 cross-term coefficient blocks keyed by LAMMPS keyword ("bb", "ba", "aat", "aa", ...).
  std::map<std::string, std::string> extra{};

  /// Permutation applied to the RASPA identifiers to obtain LAMMPS's I,J,K,L ordering for this style.
  std::array<std::size_t, 4> atomOrder{0, 1, 2, 3};

  /// Key used to merge identical terms into one LAMMPS type.
  std::string key() const;
};

std::string formatValue(double value);
std::string formatValues(std::span<const double> values);

// ---------------------------------------------------------------------------------------------------
// Cosine-series torsions
// ---------------------------------------------------------------------------------------------------

/// Fourier coefficients a_n (internal energy units) with E = sum_n a_n cos(n phi), for every RASPA torsion
/// type that is a pure cosine series; nullopt otherwise (Harmonic, Fixed, phase-shifted ModifiedTraPPE, ...).
std::optional<std::vector<double>> torsionFourierCoefficients(const TorsionPotential &torsion);

/// Rewrites sum_n a_n cos(n phi) as sum_k c_k cos^k(phi) (Chebyshev T_n expansion).
std::vector<double> chebyshevToPolynomial(std::span<const double> fourier);

// ---------------------------------------------------------------------------------------------------
// RASPA -> LAMMPS, one class at a time
// ---------------------------------------------------------------------------------------------------
Term bondTerm(const BondPotential &bond);
Term ureyBradleyTerm(const UreyBradleyPotential &ureyBradley);
Term bendTerm(const BendPotential &bend);
Term torsionTerm(const TorsionPotential &torsion);
std::vector<Term> improperTorsionTerms(const TorsionPotential &torsion);
Term inversionBendTerm(const InversionBendPotential &inversion);
Term outOfPlaneBendTerm(const OutOfPlaneBendPotential &outOfPlane);

/// class2 cross-term coefficient strings; nullopt when the RASPA form is not the class2 form.
std::optional<std::string> bondBondClass2(const BondBondPotential &bondBond);
std::optional<std::string> bondBendClass2(const BondBendPotential &bondBend, double bendTheta0);
std::optional<std::string> bendTorsionClass2(const BendTorsionPotential &bendTorsion);
std::optional<std::string> bendBendClass2(const BendBendPotential &bendBend);

/// The dihedral in 'dihedral_style class2' form (K_n, phi_n for n = 1..3); nullopt if not representable.
std::optional<std::string> torsionClass2(const TorsionPotential &torsion);
/// The bend in 'angle_style class2' form (theta0 K2 K3 K4); nullopt if not representable.
std::optional<std::string> bendClass2(const BendPotential &bend);

// ---------------------------------------------------------------------------------------------------
// Pair styles
// ---------------------------------------------------------------------------------------------------
struct PairTerm
{
  std::string style{};         ///< "lj/cut", "buck", "morse", "mie/cut", "born", "lj/smooth/linear", "table", "zero"
  std::string coefficients{};  ///< without i j and without the style-specific cut-off
  bool exact{true};
  std::string note{};
};
PairTerm pairTerm(const ForceField &forceField, std::size_t typeA, std::size_t typeB);

// ---------------------------------------------------------------------------------------------------
// Table bodies (LAMMPS table-file sections) tabulated from RASPA's own energy functions
// ---------------------------------------------------------------------------------------------------
struct TableGrid
{
  std::size_t bondPoints{901};      ///< bond tables: r in [bondLow, bondHigh]
  double bondLow{0.5};
  double bondHigh{5.0};
  std::size_t anglePoints{1801};    ///< angle tables: theta in [0, 180] degrees
  std::size_t dihedralPoints{720};  ///< dihedral tables: phi in [-180, 180) degrees
  std::size_t pairPoints{4000};     ///< pair tables: r in [pairLow, cut-off]
  double pairLow{0.5};
};

std::string bondTableBody(const BondPotential &bond, std::string_view keyword, const TableGrid &grid);
std::string angleTableBody(const BendPotential &bend, std::string_view keyword, const TableGrid &grid);
std::string dihedralTableBody(const TorsionPotential &torsion, std::string_view keyword, const TableGrid &grid);
std::string pairTableBody(const ForceField &forceField, std::size_t typeA, std::size_t typeB, double cutOff,
                          std::string_view keyword, const TableGrid &grid);

// ---------------------------------------------------------------------------------------------------
// LAMMPS -> RASPA (used by the data-file reader). Coefficients are the LAMMPS real-unit numbers.
// Returned parameters are in RASPA's JSON input units (K, Angstrom, degrees).
// ---------------------------------------------------------------------------------------------------
struct RaspaTerm
{
  std::string type{};               ///< RASPA type name as used in the component JSON ("HARMONIC", ...)
  std::vector<double> parameters{};
  std::string note{};
};
std::optional<RaspaTerm> bondFromLammps(std::string_view style, std::span<const double> coefficients);
std::optional<RaspaTerm> bendFromLammps(std::string_view style, std::span<const double> coefficients);
std::optional<RaspaTerm> torsionFromLammps(std::string_view style, std::span<const double> coefficients);
/// Returns {epsilon [K], sigma [A]} for lj/cut-like pair styles; nullopt otherwise.
std::optional<std::array<double, 2>> lennardJonesFromLammps(std::string_view style,
                                                            std::span<const double> coefficients);
}  // namespace LAMMPS
