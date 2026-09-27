module;

module lammps_styles;

import std;

import double3;
import units;
import forcefield;
import vdwparameters;
import potential_pair_vdw;
import potential_pair_derivatives;
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

namespace
{
// The unit system is a runtime setting of RASPA, so the energy conversions are read on use.
const double &toKCal = Units::EnergyToKCalPerMol;
double kelvinPerKCal() { return Units::KCalPerMolToEnergy * Units::EnergyToKelvin; }
constexpr double toDeg = Units::RadiansToDegrees;

LAMMPS::Term zero(std::string note, std::string coefficients = {})
{
  LAMMPS::Term term{};
  term.style = "zero";
  term.coefficients = std::move(coefficients);
  term.exact = false;
  term.note = std::move(note);
  return term;
}

LAMMPS::Term table(std::string_view prefix, std::string note = {})
{
  // The keyword is completed by the topology builder (it appends a running index); here it only carries
  // the class prefix so the builder knows which table generator to call.
  LAMMPS::Term term{};
  term.style = "table";
  term.tableKeyword = std::string(prefix);
  term.exact = true;
  term.note = std::move(note);
  return term;
}

// Signed dihedral angle with RASPA's convention (phi = pi for trans), for the table generator.
double raspaDihedral(const double3 &A, const double3 &B, const double3 &C, const double3 &D)
{
  double3 Dab = A - B;
  double3 Dcb = (C - B).normalized();
  double3 Ddc = D - C;
  double dot_ab = double3::dot(Dab, Dcb);
  double dot_dc = double3::dot(Ddc, Dcb);
  double3 dr = (Dab - dot_ab * Dcb).normalized();
  double3 ds = (Ddc - dot_dc * Dcb).normalized();
  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  double sign = double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc)));
  return std::copysign(std::acos(cos_phi), sign);
}
}  // namespace

namespace LAMMPS
{
std::string formatValue(double value)
{
  if (std::abs(value) < 1e-14) value = 0.0;
  return std::format("{:.10g}", value);
}

std::string formatValues(std::span<const double> values)
{
  std::string out{};
  for (std::size_t i = 0; i != values.size(); ++i)
  {
    if (i) out += " ";
    out += formatValue(values[i]);
  }
  return out;
}

std::string Term::key() const
{
  std::string k = style + "|" + coefficients + "|" + tableKeyword;
  for (const auto &[name, value] : extra) k += "|" + name + ":" + value;
  return k;
}

// ---------------------------------------------------------------------------------------------------
// Cosine series
// ---------------------------------------------------------------------------------------------------
std::optional<std::vector<double>> torsionFourierCoefficients(const TorsionPotential &torsion)
{
  const std::array<double, maximumNumberOfTorsionParameters> &p = torsion.parameters;
  switch (torsion.type)
  {
    case TorsionType::ThreeCosine:
    case TorsionType::MM3:
      // (1/2) p0 (1 + cos phi) + (1/2) p1 (1 - cos 2phi) + (1/2) p2 (1 + cos 3phi)
      return std::vector<double>{0.5 * (p[0] + p[1] + p[2]), 0.5 * p[0], -0.5 * p[1], 0.5 * p[2]};
    case TorsionType::TraPPE:
      // p0 + p1 (1 + cos phi) + p2 (1 - cos 2phi) + p3 (1 + cos 3phi)
      return std::vector<double>{p[0] + p[1] + p[2] + p[3], p[1], -p[2], p[3]};
    case TorsionType::TraPPE_Extended:
      // sum_{n=0}^{4} p_n cos(n phi)
      return std::vector<double>{p[0], p[1], p[2], p[3], p[4]};
    case TorsionType::ModifiedTraPPE:
      if (std::abs(p[4]) < 1e-12)
      {
        return std::vector<double>{p[0] + p[1] + p[2] + p[3], p[1], -p[2], p[3]};
      }
      return std::nullopt;
    case TorsionType::CFF:
      // p0 (1 - cos phi) + p1 (1 - cos 2phi) + p2 (1 - cos 3phi)
      return std::vector<double>{p[0] + p[1] + p[2], -p[0], -p[1], -p[2]};
    case TorsionType::CFF2:
      // p0 (1 + cos phi) + p1 (1 + cos 2phi) + p2 (1 + cos 3phi)
      return std::vector<double>{p[0] + p[1] + p[2], p[0], p[1], p[2]};
    case TorsionType::OPLS:
      // (1/2) p0 + (1/2) p1 (1 + cos phi) + (1/2) p2 (1 - cos 2phi) + (1/2) p3 (1 + cos 3phi)
      return std::vector<double>{0.5 * (p[0] + p[1] + p[2] + p[3]), 0.5 * p[1], -0.5 * p[2], 0.5 * p[3]};
    case TorsionType::FourierSeries:
      // (1/2)[p0 (1 + cos phi) + p1 (1 - cos 2phi) + p2 (1 + cos 3phi) + p3 (1 - cos 4phi) + p4 (1 + cos 5phi)
      //       + p5 (1 - cos 6phi)]
      // RASPA's implementation (as RASPA2's) alternates the sign for every even harmonic, including n = 6;
      // the header comment in torsion_potential.cpp that says (1 + cos 6phi) does not match the code.
      return std::vector<double>{0.5 * (p[0] + p[1] + p[2] + p[3] + p[4] + p[5]), 0.5 * p[0], -0.5 * p[1],
                                 0.5 * p[2], -0.5 * p[3], 0.5 * p[4], -0.5 * p[5]};
    case TorsionType::FourierSeries2:
      return std::vector<double>{0.5 * (p[0] + p[1] + p[2] + p[3] + p[4] + p[5]), 0.5 * p[0], -0.5 * p[1],
                                 0.5 * p[2], 0.5 * p[3], 0.5 * p[4], 0.5 * p[5]};
    case TorsionType::CVFF:
    {
      // p0 (1 + cos(p1 phi - p2)) is a cosine series only for integer p1 and p2 in {0, pi}
      const long n = std::lround(p[1]);
      if (n < 0 || n > 6 || std::abs(p[1] - static_cast<double>(n)) > 1e-9) return std::nullopt;
      const double c = std::cos(p[2]);
      if (std::abs(std::abs(c) - 1.0) > 1e-9) return std::nullopt;
      std::vector<double> a(static_cast<std::size_t>(n) + 1, 0.0);
      a[0] += p[0];
      a[static_cast<std::size_t>(n)] += p[0] * c;
      return a;
    }
    default:
      return std::nullopt;
  }
}

std::vector<double> chebyshevToPolynomial(std::span<const double> fourier)
{
  // T_0 = 1, T_1 = c, T_{n+1} = 2 c T_n - T_{n-1}
  std::vector<double> polynomial(fourier.size(), 0.0);
  std::vector<double> previous{1.0};
  std::vector<double> current{0.0, 1.0};
  for (std::size_t n = 0; n != fourier.size(); ++n)
  {
    const std::vector<double> &Tn = (n == 0) ? previous : current;
    for (std::size_t k = 0; k != Tn.size(); ++k) polynomial[k] += fourier[n] * Tn[k];
    if (n >= 1)
    {
      std::vector<double> next(n + 2, 0.0);
      for (std::size_t k = 0; k != current.size(); ++k) next[k + 1] += 2.0 * current[k];
      for (std::size_t k = 0; k != previous.size(); ++k) next[k] -= previous[k];
      previous = std::move(current);
      current = std::move(next);
    }
  }
  return polynomial;
}

// ---------------------------------------------------------------------------------------------------
// Bonds
// ---------------------------------------------------------------------------------------------------
Term bondTerm(const BondPotential &bond)
{
  const auto &p = bond.parameters;
  Term term{};
  switch (bond.type)
  {
    case BondType::Harmonic:
      // (1/2) p0 (r - p1)^2  ->  harmonic: K r0
      term.style = "harmonic";
      term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, p[1]});
      return term;
    case BondType::CoreShellSpring:
      term.style = "harmonic";
      term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, 0.0});
      return term;
    case BondType::Morse:
      // p0 [(1 - exp(-p1 (r - p2)))^2 - 1]  ->  morse: D0 alpha r0 (LAMMPS omits the constant -D0)
      term.style = "morse";
      term.coefficients = formatValues(std::array{p[0] * toKCal, p[1], p[2]});
      term.exact = false;
      term.note = "constant -D0 dropped";
      return term;
    case BondType::Quartic:
      // (1/2) p0 d^2 + (1/3) p2 d^3 + (1/4) p3 d^4, d = r - p1  ->  class2: r0 K2 K3 K4
      term.style = "class2";
      term.coefficients =
          formatValues(std::array{p[1], 0.5 * p[0] * toKCal, p[2] * toKCal / 3.0, 0.25 * p[3] * toKCal});
      return term;
    case BondType::CFF_Quartic:
      term.style = "class2";
      term.coefficients = formatValues(std::array{p[1], p[0] * toKCal, p[2] * toKCal, p[3] * toKCal});
      return term;
    case BondType::MM3:
    case BondType::LJ_12_6:
    case BondType::LennardJones:
    case BondType::Buckingham:
    case BondType::RestrainedHarmonic:
      return table("B");
    case BondType::Fixed:
      return zero("fixed bond length: constrained with fix shake (see input script)", formatValue(p[0]));
    case BondType::None:
    default:
      return zero("no LAMMPS equivalent");
  }
}

// ---------------------------------------------------------------------------------------------------
// Urey-Bradley: 'angle_style charmm' with K = 0 carries a pure 1-3 harmonic term.
// ---------------------------------------------------------------------------------------------------
Term ureyBradleyTerm(const UreyBradleyPotential &ureyBradley)
{
  const auto &p = ureyBradley.parameters;
  Term term{};
  switch (ureyBradley.type)
  {
    case UreyBradleyType::Harmonic:
      // (1/2) p0 (r13 - p1)^2  ->  charmm: K theta0 K_ub r_ub with K = 0
      term.style = "charmm";
      term.coefficients = formatValues(std::array{0.0, 0.0, 0.5 * p[0] * toKCal, p[1]});
      return term;
    case UreyBradleyType::CoreShellSpring:
      term.style = "charmm";
      term.coefficients = formatValues(std::array{0.0, 0.0, 0.5 * p[0] * toKCal, 0.0});
      return term;
    case UreyBradleyType::Fixed:
      return zero("fixed 1-3 distance: no LAMMPS angle constraint; use fix shake on the angle", "");
    default:
      return zero(std::format("Urey-Bradley type {} has no LAMMPS equivalent",
                              static_cast<std::size_t>(ureyBradley.type)));
  }
}

// ---------------------------------------------------------------------------------------------------
// Bends
// ---------------------------------------------------------------------------------------------------
Term bendTerm(const BendPotential &bend)
{
  const auto &p = bend.parameters;
  Term term{};
  switch (bend.type)
  {
    case BendType::Harmonic:
    case BendType::CoreShell:
      // (1/2) p0 (theta - p1)^2  ->  harmonic: K theta0
      term.style = "harmonic";
      term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg});
      return term;
    case BendType::Quartic:
      term.style = "quartic";
      term.coefficients =
          formatValues(std::array{p[1] * toDeg, 0.5 * p[0] * toKCal, p[2] * toKCal / 3.0, 0.25 * p[3] * toKCal});
      return term;
    case BendType::CFF_Quartic:
      term.style = "quartic";
      term.coefficients = formatValues(std::array{p[1] * toDeg, p[0] * toKCal, p[2] * toKCal, p[3] * toKCal});
      return term;
    case BendType::HarmonicCosine:
      // (1/2) p0 (cos theta - p1)^2, p1 = cos theta0  ->  cosine/squared: K theta0
      term.style = "cosine/squared";
      term.coefficients =
          formatValues(std::array{0.5 * p[0] * toKCal, std::acos(std::clamp(p[1], -1.0, 1.0)) * toDeg});
      return term;
    case BendType::Cosine:
    case BendType::Tafipolsky:
    case BendType::MM3:
    case BendType::MM3_inplane:
      return table("A");
    case BendType::Fixed:
      return zero("fixed bend angle: constrained with fix shake (see input script)", formatValue(p[0] * toDeg));
    case BendType::Rigid:
      return zero("bend inside a rigid fragment: handled by fix rigid", "");
    default:
      return zero("no LAMMPS equivalent");
  }
}

std::optional<std::string> bendClass2(const BendPotential &bend)
{
  const auto &p = bend.parameters;
  switch (bend.type)
  {
    case BendType::Harmonic:
      return formatValues(std::array{p[1] * toDeg, 0.5 * p[0] * toKCal, 0.0, 0.0});
    case BendType::Quartic:
      return formatValues(std::array{p[1] * toDeg, 0.5 * p[0] * toKCal, p[2] * toKCal / 3.0, 0.25 * p[3] * toKCal});
    case BendType::CFF_Quartic:
      return formatValues(std::array{p[1] * toDeg, p[0] * toKCal, p[2] * toKCal, p[3] * toKCal});
    default:
      return std::nullopt;
  }
}

// ---------------------------------------------------------------------------------------------------
// Torsions
// ---------------------------------------------------------------------------------------------------
namespace
{
/// 'nharmonic n A1 .. An' with trailing zero coefficients removed (n >= 1), so identical series merge.
Term nharmonicTerm(std::vector<double> polynomial)
{
  while (polynomial.size() > 1 && std::abs(polynomial.back()) < 1e-14) polynomial.pop_back();
  Term term{};
  term.style = "nharmonic";
  term.coefficients = std::format("{} {}", polynomial.size(), formatValues(polynomial));
  return term;
}
}  // namespace

Term torsionTerm(const TorsionPotential &torsion)
{
  const auto &p = torsion.parameters;
  Term term{};
  switch (torsion.type)
  {
    case TorsionType::Harmonic:
      // (1/2) p0 (phi - p1)^2  ->  quadratic: K phi0
      term.style = "quadratic";
      term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg});
      return term;
    case TorsionType::HarmonicCosine:
    {
      // (1/2) p0 (cos phi - c0)^2  ->  nharmonic 3: (1/2) p0 c0^2, -p0 c0, (1/2) p0
      const double k = p[0] * toKCal, c0 = p[1];
      return nharmonicTerm({0.5 * k * c0 * c0, -k * c0, 0.5 * k});
    }
    case TorsionType::RyckaertBellemans:
    {
      // sum_i (-1)^i p_i cos^i phi
      std::vector<double> polynomial(6);
      for (std::size_t i = 0; i != 6; ++i) polynomial[i] = ((i % 2) ? -1.0 : 1.0) * p[i] * toKCal;
      return nharmonicTerm(polynomial);
    }
    case TorsionType::Polynomial:
    {
      std::vector<double> polynomial(6);
      for (std::size_t i = 0; i != 6; ++i) polynomial[i] = p[i] * toKCal;
      return nharmonicTerm(polynomial);
    }
    case TorsionType::CVFF:
    {
      // p0 (1 + cos(p1 phi - p2))  ->  charmm: K n d w  (integer n, integer d in degrees)
      const long n = std::lround(p[1]);
      const long d = std::lround(p[2] * toDeg);
      term.style = "charmm";
      term.coefficients = std::format("{} {} {} 0.0", formatValue(p[0] * toKCal), n, d);
      if (std::abs(p[1] - static_cast<double>(n)) > 1e-9 || std::abs(p[2] * toDeg - static_cast<double>(d)) > 1e-6)
      {
        term.exact = false;
        term.note = "charmm needs integer n and d (degrees); values were rounded";
      }
      return term;
    }
    case TorsionType::Fixed:
      return zero("fixed torsion: no LAMMPS constraint; keep the fragment rigid", "");
    case TorsionType::CVFFBlocked:
      return zero("CVFF_BLOCKED carries no energy", "");
    default:
      break;
  }

  std::optional<std::vector<double>> fourier = torsionFourierCoefficients(torsion);
  if (fourier.has_value())
  {
    std::vector<double> a = *fourier;
    for (double &value : a) value *= toKCal;
    return nharmonicTerm(chebyshevToPolynomial(a));
  }

  // phase-shifted or otherwise non-polynomial forms: tabulate
  return table("D", "tabulated; assumes LAMMPS and RASPA share the sign convention of phi");
}

std::optional<std::string> torsionClass2(const TorsionPotential &torsion)
{
  // dihedral_style class2: E = sum_{n=1}^{3} K_n [1 - cos(n phi - phi_n)]
  const auto &p = torsion.parameters;
  switch (torsion.type)
  {
    case TorsionType::CFF:
      // p_i (1 - cos(i phi))
      return formatValues(std::array{p[0] * toKCal, 0.0, p[1] * toKCal, 0.0, p[2] * toKCal, 0.0});
    case TorsionType::CFF2:
      // p_i (1 + cos(i phi)) = p_i (1 - cos(i phi - 180))
      return formatValues(std::array{p[0] * toKCal, 180.0, p[1] * toKCal, 180.0, p[2] * toKCal, 180.0});
    default:
      return std::nullopt;
  }
}

// ---------------------------------------------------------------------------------------------------
// Improper torsions (dihedral-angle based): every cosine series becomes one 'improper_style cvff' term
// per harmonic, E = K [1 + d cos(n phi)], with the constant collected and reported.
// ---------------------------------------------------------------------------------------------------
std::vector<Term> improperTorsionTerms(const TorsionPotential &torsion)
{
  const auto &p = torsion.parameters;
  if (torsion.type == TorsionType::Harmonic)
  {
    // improper_style harmonic: K (chi - chi0)^2 with chi the (unsigned) angle between the ijk and jkl planes
    Term term{};
    term.style = "harmonic";
    term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg});
    const double phi0 = std::fmod(std::abs(p[1] * toDeg), 360.0);
    if (std::abs(phi0) > 1e-6 && std::abs(phi0 - 180.0) > 1e-6)
    {
      term.exact = false;
      term.note = "LAMMPS improper harmonic uses an unsigned angle; exact only for chi0 = 0 or 180";
    }
    return {term};
  }
  if (torsion.type == TorsionType::Fixed || torsion.type == TorsionType::CVFFBlocked)
  {
    return {zero("carries no energy", "")};
  }

  std::optional<std::vector<double>> fourier = torsionFourierCoefficients(torsion);
  if (!fourier.has_value())
  {
    return {zero("improper torsion form has no LAMMPS equivalent", "")};
  }

  std::vector<Term> terms{};
  double constant = (*fourier)[0];
  for (std::size_t n = 1; n < fourier->size(); ++n)
  {
    const double a = (*fourier)[n];
    if (std::abs(a) < 1e-14) continue;
    Term term{};
    term.style = "cvff";
    term.coefficients = std::format("{} {} {}", formatValue(std::abs(a) * toKCal), (a < 0.0) ? -1 : 1, n);
    constant -= std::abs(a);
    terms.push_back(term);
  }
  if (terms.empty())
  {
    return {zero(std::format("constant improper energy {} kcal/mol dropped", formatValue(constant * toKCal)), "")};
  }
  if (std::abs(constant) > 1e-12)
  {
    terms.front().exact = false;
    terms.front().note = std::format("constant {} kcal/mol dropped", formatValue(constant * toKCal));
  }
  return terms;
}

// ---------------------------------------------------------------------------------------------------
// Inversion bends (Wilson angle chi between bond B-A and the plane B-C-D, B central)
// -> improper_style umbrella with LAMMPS ordering I = B (central), J = C, K = D, L = A
// ---------------------------------------------------------------------------------------------------
Term inversionBendTerm(const InversionBendPotential &inversion)
{
  const auto &p = inversion.parameters;
  Term term{};
  term.atomOrder = {1, 2, 3, 0};
  switch (inversion.type)
  {
    case InversionBendType::Planar:
      // p0 (1 - cos chi)  ->  umbrella: K omega0 with omega0 = 0: E = K (1 - cos omega)
      term.style = "umbrella";
      term.coefficients = formatValues(std::array{p[0] * toKCal, 0.0});
      return term;
    case InversionBendType::HarmonicCosine:
    {
      // (1/2) p0 (cos chi - c0)^2  ->  umbrella: (1/2) K / sin^2(omega0) (cos omega - cos omega0)^2
      const double c0 = std::clamp(p[1], -1.0, 1.0);
      const double sin2 = 1.0 - c0 * c0;
      if (sin2 < 1e-12)
      {
        return zero("HARMONIC_COSINE inversion with chi0 = 0 has no umbrella equivalent", "");
      }
      term.style = "umbrella";
      term.coefficients = formatValues(std::array{p[0] * sin2 * toKCal, std::acos(c0) * toDeg});
      return term;
    }
    case InversionBendType::Harmonic:
      // (1/2) p0 (chi - chi0)^2  ->  improper_style inversion/harmonic (YAFF package): K (omega - omega0)^2
      term.style = "inversion/harmonic";
      term.coefficients = formatValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg});
      term.exact = false;
      term.note = "verify the central-atom ordering of improper_style inversion/harmonic for your LAMMPS version";
      return term;
    default:
      return zero("inversion-bend type has no LAMMPS equivalent (plane through A-C-D or MM3 form)", "");
  }
}

Term outOfPlaneBendTerm([[maybe_unused]] const OutOfPlaneBendPotential &outOfPlane)
{
  return zero("out-of-plane bend carries no energy in RASPA", "");
}

// ---------------------------------------------------------------------------------------------------
// class2 cross terms
// ---------------------------------------------------------------------------------------------------
std::optional<std::string> bondBondClass2(const BondBondPotential &bondBond)
{
  const auto &p = bondBond.parameters;
  switch (bondBond.type)
  {
    case BondBondType::CVFF:
    case BondBondType::CFF:
      // p0 (r_ab - p1)(r_cb - p2)  ->  bb: M r1 r2
      return formatValues(std::array{p[0] * toKCal, p[1], p[2]});
    default:
      return std::nullopt;
  }
}

std::optional<std::string> bondBendClass2(const BondBendPotential &bondBend, [[maybe_unused]] double bendTheta0)
{
  const auto &p = bondBend.parameters;
  switch (bondBend.type)
  {
    case BondBendType::CVFF:
    case BondBendType::CFF:
      // (theta - p0)(p1 (r_ab - p2) + p3 (r_cb - p4))  ->  ba: N1 N2 r1 r2  (theta0 taken from the angle)
      return formatValues(std::array{p[1] * toKCal, p[3] * toKCal, p[2], p[4]});
    default:
      return std::nullopt;
  }
}

std::optional<std::string> bendTorsionClass2(const BendTorsionPotential &bendTorsion)
{
  const auto &p = bendTorsion.parameters;
  switch (bendTorsion.type)
  {
    case BendTorsionType::CVFF:
    case BendTorsionType::CFF:
      // p0 (theta1 - p1)(theta2 - p2) cos phi  ->  aat: M theta1 theta2
      return formatValues(std::array{p[0] * toKCal, p[1] * toDeg, p[2] * toDeg});
    default:
      return std::nullopt;
  }
}

std::optional<std::string> bendBendClass2(const BendBendPotential &bendBend)
{
  const auto &p = bendBend.parameters;
  switch (bendBend.type)
  {
    case BendBendType::CVFF:
    case BendBendType::CFF:
      // p0 (theta_ABC - p1)(theta_ABD - p2), B central  ->  improper class2 aa: M1 M2 M3 theta1 theta2 theta3
      // with LAMMPS I=A, J=B (central), K=C, L=D: theta_ijk = ABC, theta_ijl = ABD, so only M2 is non-zero.
      return formatValues(std::array{0.0, p[0] * toKCal, 0.0, p[1] * toDeg, p[2] * toDeg, 0.0});
    default:
      return std::nullopt;
  }
}

// ---------------------------------------------------------------------------------------------------
// Pair styles
// ---------------------------------------------------------------------------------------------------
PairTerm pairTerm(const ForceField &forceField, std::size_t typeA, std::size_t typeB)
{
  const VDWParameters &vdw = forceField(typeA, typeB);
  const double x = vdw.parameters.x, y = vdw.parameters.y, z = vdw.parameters.z, w = vdw.parameters.w;
  PairTerm term{};
  switch (vdw.type)
  {
    case VDWParameters::Type::None:
      term.style = "zero";
      return term;
    case VDWParameters::Type::LennardJones:
      term.style = "lj/cut";
      term.coefficients = formatValues(std::array{x * toKCal, y});
      return term;
    case VDWParameters::Type::LennardJonesShiftedForce:
      term.style = "lj/smooth/linear";
      term.coefficients = formatValues(std::array{x * toKCal, y});
      return term;
    case VDWParameters::Type::BuckingHam:
      // x exp(-y r) - z / r^6  ->  buck: A rho C
      term.style = "buck";
      term.coefficients = formatValues(std::array{x * toKCal, 1.0 / y, z * toKCal});
      return term;
    case VDWParameters::Type::Morse:
      // x [(1 - exp(-y (r - z)))^2 - 1] = x e^{-2y(r-z)} - 2x e^{-y(r-z)}  ->  morse: D0 alpha r0
      term.style = "morse";
      term.coefficients = formatValues(std::array{x * toKCal, y, z});
      return term;
    case VDWParameters::Type::BornHugginsMeyer:
      // x exp(y (z - r)) - w/r^6 - p2.x/r^8 - p2.y/r^10  ->  born: A rho sigma C D  (E = ... - C/r^6 + D/r^8)
      if (std::abs(vdw.parameters2.y) > 0.0) break;
      term.style = "born";
      term.coefficients = formatValues(std::array{x * toKCal, 1.0 / y, z, w * toKCal, -vdw.parameters2.x * toKCal});
      return term;
    case VDWParameters::Type::Potential12_6:
      // x / r^12 - y / r^6  ->  lj/cut with sigma = (x/y)^(1/6), epsilon = y^2 / (4x)
      if (x > 0.0 && y > 0.0)
      {
        const double sigma = std::pow(x / y, 1.0 / 6.0);
        term.style = "lj/cut";
        term.coefficients = formatValues(std::array{(y * y / (4.0 * x)) * toKCal, sigma});
        return term;
      }
      break;
    case VDWParameters::Type::CFFEpsilonSigma:
      // x [2 (y/r)^9 - 3 (y/r)^6]  ->  lj/class2: epsilon sigma
      term.style = "lj/class2";
      term.coefficients = formatValues(std::array{x * toKCal, y});
      return term;
    case VDWParameters::Type::CFF9_6:
      // x / r^9 - y / r^6 = eps [2 (s/r)^9 - 3 (s/r)^6]: s^3 = 3x / (2y), eps = y / (3 s^6)
      if (x > 0.0 && y > 0.0)
      {
        const double sigma = std::cbrt(1.5 * x / y);
        term.style = "lj/class2";
        term.coefficients = formatValues(std::array{(y / (3.0 * std::pow(sigma, 6.0))) * toKCal, sigma});
        return term;
      }
      break;
    case VDWParameters::Type::Mie:
    {
      // x / r^y - z / r^w  ->  mie/cut: epsilon sigma gammaR gammaA with
      // E = C eps [(s/r)^gR - (s/r)^gA], C = gR/(gR-gA) (gR/gA)^(gA/(gR-gA))
      const double gR = y, gA = w;
      if (x > 0.0 && z > 0.0 && gR > gA && gA > 0.0)
      {
        const double sigma = std::pow(x / z, 1.0 / (gR - gA));
        const double C = gR / (gR - gA) * std::pow(gR / gA, gA / (gR - gA));
        const double epsilon = z / (C * std::pow(sigma, gA));
        term.style = "mie/cut";
        term.coefficients = formatValues(std::array{epsilon * toKCal, sigma, gR, gA});
        return term;
      }
      break;
    }
    default:
      break;
  }
  term.style = "table";
  term.note = std::format("{} pair potential tabulated", VDWParameters::nameOfType(vdw.type));
  return term;
}

// ---------------------------------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------------------------------
std::string bondTableBody(const BondPotential &bond, std::string_view keyword, const TableGrid &grid)
{
  std::ostringstream out;
  const std::size_t n = grid.bondPoints;
  const double dr = (grid.bondHigh - grid.bondLow) / static_cast<double>(n - 1);
  auto energy = [&](double r) { return bond.calculateEnergy(double3{0.0, 0.0, 0.0}, double3{r, 0.0, 0.0}) * toKCal; };
  std::print(out, "# RASPA bond type {} tabulated: r [A], E [kcal/mol], F = -dE/dr [kcal/mol/A]\n",
             static_cast<std::size_t>(bond.type));
  std::print(out, "{}\nN {}\n\n", keyword, n);
  const double h = 1e-5;
  for (std::size_t i = 0; i < n; ++i)
  {
    const double r = grid.bondLow + static_cast<double>(i) * dr;
    const double force = -(energy(r + h) - energy(r - h)) / (2.0 * h);
    std::print(out, "{} {} {} {}\n", i + 1, formatValue(r), formatValue(energy(r)), formatValue(force));
  }
  return out.str();
}

std::string angleTableBody(const BendPotential &bend, std::string_view keyword, const TableGrid &grid)
{
  std::ostringstream out;
  const std::size_t n = grid.anglePoints;
  const double dtheta = 180.0 / static_cast<double>(n - 1);
  auto energy = [&](double thetaDeg)
  {
    const double theta = thetaDeg * Units::DegreesToRadians;
    return bend.calculateEnergy(double3{std::cos(theta), std::sin(theta), 0.0}, double3{0.0, 0.0, 0.0},
                                double3{1.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}) *
           toKCal;
  };
  std::print(out,
             "# RASPA bend type {} tabulated: theta [degrees], E [kcal/mol], F = -dE/dtheta [kcal/mol/degree]\n",
             static_cast<std::size_t>(bend.type));
  std::print(out, "{}\nN {}\n\n", keyword, n);
  const double h = 1e-3;  // degrees
  for (std::size_t i = 0; i < n; ++i)
  {
    const double theta = static_cast<double>(i) * dtheta;
    const double lo = std::max(0.0, theta - h), hi = std::min(180.0, theta + h);
    const double force = -(energy(hi) - energy(lo)) / (hi - lo);
    std::print(out, "{} {} {} {}\n", i + 1, formatValue(theta), formatValue(energy(theta)), formatValue(force));
  }
  return out.str();
}

std::string dihedralTableBody(const TorsionPotential &torsion, std::string_view keyword, const TableGrid &grid)
{
  std::ostringstream out;
  const std::size_t n = grid.dihedralPoints;
  const double dphi = 360.0 / static_cast<double>(n);
  std::print(out, "# RASPA torsion type {} tabulated: phi [degrees], E [kcal/mol]\n",
             static_cast<std::size_t>(torsion.type));
  std::print(out, "{}\nN {} DEGREES NOF\n\n", keyword, n);
  for (std::size_t i = 0; i < n; ++i)
  {
    const double phiDeg = -180.0 + static_cast<double>(i) * dphi;
    const double phi = phiDeg * Units::DegreesToRadians;
    // B-C along x; A in the xy-plane; D rotated by phi about the B-C axis relative to A's side
    const double3 B{0.0, 0.0, 0.0}, C{1.5, 0.0, 0.0};
    const double3 A{-0.5, 1.0, 0.0};
    double3 D{2.0, std::cos(phi), std::sin(phi)};
    // Make sure the tabulated abscissa is RASPA's signed angle for this geometry (flip D if needed).
    if (std::abs(raspaDihedral(A, B, C, D) - phi) > 1e-6 && std::abs(std::abs(phi) - std::numbers::pi) > 1e-6)
    {
      D = double3{2.0, std::cos(phi), -std::sin(phi)};
    }
    std::print(out, "{} {} {}\n", i + 1, formatValue(phiDeg), formatValue(torsion.calculateEnergy(A, B, C, D) * toKCal));
  }
  return out.str();
}

std::string pairTableBody(const ForceField &forceField, std::size_t typeA, std::size_t typeB, double cutOff,
                          std::string_view keyword, const TableGrid &grid)
{
  std::ostringstream out;
  const std::size_t n = grid.pairPoints;
  const double dr = (cutOff - grid.pairLow) / static_cast<double>(n - 1);
  std::print(out, "# RASPA {} pair potential between types {} and {}: r [A], E [kcal/mol], F [kcal/mol/A]\n",
             VDWParameters::nameOfType(forceField(typeA, typeB).type), typeA + 1, typeB + 1);
  std::print(out, "{}\nN {} R {} {}\n\n", keyword, n, formatValue(grid.pairLow), formatValue(cutOff));
  for (std::size_t i = 0; i < n; ++i)
  {
    const double r = grid.pairLow + static_cast<double>(i) * dr;
    const Potentials::PairDerivatives<1> d = Potentials::potentialVDW<1>(forceField, 1.0, 1.0, r * r, typeA, typeB);
    // firstDerivativeFactor = (1/r) dU/dr  ->  F = -dU/dr = -factor * r
    std::print(out, "{} {} {} {}\n", i + 1, formatValue(r), formatValue(d.energy * toKCal),
               formatValue(-d.firstDerivativeFactor * r * toKCal));
  }
  return out.str();
}

// ---------------------------------------------------------------------------------------------------
// LAMMPS -> RASPA
// ---------------------------------------------------------------------------------------------------
std::optional<RaspaTerm> bondFromLammps(std::string_view style, std::span<const double> c)
{
  const double K = kelvinPerKCal();
  if (style == "harmonic" && c.size() >= 2)
  {
    // K (r - r0)^2 = (1/2)(2K)(r - r0)^2
    return RaspaTerm{"HARMONIC", {2.0 * c[0] * K, c[1]}};
  }
  if (style == "morse" && c.size() >= 3)
  {
    return RaspaTerm{"MORSE", {c[0] * K, c[1], c[2]}, "RASPA's Morse bond includes the constant -D0"};
  }
  if (style == "class2" && c.size() >= 4)
  {
    // r0 K2 K3 K4  ->  CFF_QUARTIC: K2 r0 K3 K4
    return RaspaTerm{"CFF_QUARTIC", {c[1] * K, c[0], c[2] * K, c[3] * K}};
  }
  if (style == "zero")
  {
    if (!c.empty()) return RaspaTerm{"FIXED", {c[0]}};
    return RaspaTerm{"FIXED", {}, "zero bond without r0: length taken from the coordinates"};
  }
  return std::nullopt;
}

std::optional<RaspaTerm> bendFromLammps(std::string_view style, std::span<const double> c)
{
  const double K = kelvinPerKCal();
  if (style == "harmonic" && c.size() >= 2)
  {
    return RaspaTerm{"HARMONIC", {2.0 * c[0] * K, c[1]}};
  }
  if (style == "cosine/squared" && c.size() >= 2)
  {
    return RaspaTerm{"HARMONIC_COSINE", {2.0 * c[0] * K, c[1]}};
  }
  if (style == "quartic" && c.size() >= 4)
  {
    // theta0 K2 K3 K4  ->  CFF_QUARTIC: K2 theta0 K3 K4
    return RaspaTerm{"CFF_QUARTIC", {c[1] * K, c[0], c[2] * K, c[3] * K}};
  }
  if (style == "charmm" && c.size() >= 4)
  {
    // K theta0 K_ub r_ub: the Urey-Bradley part is reported separately by the reader
    return RaspaTerm{"HARMONIC", {2.0 * c[0] * K, c[1]}, (std::abs(c[2]) > 0.0) ? "Urey-Bradley term present" : ""};
  }
  if (style == "zero")
  {
    if (!c.empty()) return RaspaTerm{"FIXED", {c[0]}};
    return RaspaTerm{"FIXED", {}, "zero angle without theta0: angle taken from the coordinates"};
  }
  return std::nullopt;
}

std::optional<RaspaTerm> torsionFromLammps(std::string_view style, std::span<const double> c)
{
  const double K = kelvinPerKCal();
  if (style == "nharmonic" && c.size() >= 2)
  {
    const std::size_t n = static_cast<std::size_t>(std::lround(c[0]));
    if (n > 6 || c.size() < n + 1) return std::nullopt;
    std::vector<double> parameters(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) parameters[i] = c[i + 1] * K;
    return RaspaTerm{"POLYNOMIAL", parameters};
  }
  if (style == "multi/harmonic" && c.size() >= 5)
  {
    std::vector<double> parameters(5);
    for (std::size_t i = 0; i < 5; ++i) parameters[i] = c[i] * K;
    return RaspaTerm{"POLYNOMIAL", parameters};
  }
  if (style == "trappe" && c.size() >= 5)
  {
    // Custom style (dihedral_trappe.cpp, Torres-Knoop et al.): U = sum_n (1/2) K_n cos(n phi), n = 0..4
    return RaspaTerm{"TRAPPE_EXTENDED", {0.5 * c[0] * K, 0.5 * c[1] * K, 0.5 * c[2] * K, 0.5 * c[3] * K, 0.5 * c[4] * K}};
  }
  if (style == "opls" && c.size() >= 4)
  {
    // (1/2) K1 (1 + cos phi) + (1/2) K2 (1 - cos 2phi) + (1/2) K3 (1 + cos 3phi) + (1/2) K4 (1 - cos 4phi)
    // = FOURIER_SERIES with p0..p3 = K1..K4
    return RaspaTerm{"FOURIER_SERIES", {c[0] * K, c[1] * K, c[2] * K, c[3] * K, 0.0, 0.0}};
  }
  if (style == "harmonic" && c.size() >= 3)
  {
    // K [1 + d cos(n phi)]  ->  TRAPPE_EXTENDED with p0 = K, p_n = d K
    const long n = std::lround(c[2]);
    if (n < 1 || n > 4) return std::nullopt;
    std::vector<double> parameters(5, 0.0);
    parameters[0] = c[0] * K;
    parameters[static_cast<std::size_t>(n)] = c[1] * c[0] * K;
    return RaspaTerm{"TRAPPE_EXTENDED", parameters};
  }
  if (style == "charmm" && c.size() >= 3)
  {
    // K [1 + cos(n phi - d)]  ->  CVFF: p0 = K, p1 = n, p2 = d
    return RaspaTerm{"CVFF", {c[0] * K, c[1], c[2]}};
  }
  if (style == "quadratic" && c.size() >= 2)
  {
    return RaspaTerm{"HARMONIC", {2.0 * c[0] * K, c[1]}};
  }
  if (style == "zero")
  {
    return RaspaTerm{"FIXED", {}};
  }
  return std::nullopt;
}

std::optional<std::array<double, 2>> lennardJonesFromLammps(std::string_view style, std::span<const double> c)
{
  if (c.size() < 2) return std::nullopt;
  if (style.starts_with("lj/cut") || style.starts_with("lj/charmm") || style == "lj/smooth/linear" ||
      style.starts_with("lj/long"))
  {
    return std::array{c[0] * kelvinPerKCal(), c[1]};
  }
  return std::nullopt;
}
}  // namespace LAMMPS
