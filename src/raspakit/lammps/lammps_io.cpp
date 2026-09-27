module;

module lammps_io;

import std;

import double3;
import double4;
import component;
import atom;
import atom_dynamics;
import molecule;
import simulationbox;
import forcefield;
import vdwparameters;
import pseudo_atom;
import units;
import framework;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import connectivity_table;

namespace
{
// ---------------------------------------------------------------------------------------------------
// One coefficient line of a Bond/Angle/Dihedral Coeffs section: the LAMMPS style it belongs to, the
// coefficients in LAMMPS 'real' units, and an optional remark (a dropped constant, an unsupported form).
// A class whose lines span several styles is written for 'hybrid' (style name on every line).
// ---------------------------------------------------------------------------------------------------
struct CoefficientLine
{
  std::string style{};
  std::string coefficients{};
  std::string note{};
};

std::string fmt(double value) { return std::format("{:.10g}", value); }

std::string joinValues(std::span<const double> values)
{
  std::string out{};
  for (std::size_t i = 0; i != values.size(); ++i)
  {
    if (i) out += " ";
    out += fmt(values[i]);
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------
// Bonds. RASPA's harmonic forms carry a 1/2 that LAMMPS's do not (E_lammps = K (r - r0)^2).
// ---------------------------------------------------------------------------------------------------
CoefficientLine bondCoefficients(const BondPotential &bond)
{
  const double toKCal = Units::EnergyToKCalPerMol;
  const std::array<double, maximumNumberOfBondParameters> &p = bond.parameters;
  switch (bond.type)
  {
    case BondType::Harmonic:
      // (1/2) p0 (r - p1)^2  ->  harmonic: K r0
      return {"harmonic", joinValues(std::array{0.5 * p[0] * toKCal, p[1]}), ""};
    case BondType::CoreShellSpring:
      // (1/2) p0 r^2  ->  harmonic with r0 = 0
      return {"harmonic", joinValues(std::array{0.5 * p[0] * toKCal, 0.0}), ""};
    case BondType::Morse:
      // p0 [(1 - exp(-p1 (r - p2)))^2 - 1]  ->  morse: D0 alpha r0 (the constant -p0 is dropped)
      return {"morse", joinValues(std::array{p[0] * toKCal, p[1], p[2]}), "constant -D0 dropped"};
    case BondType::Quartic:
      // (1/2) p0 d^2 + (1/3) p2 d^3 + (1/4) p3 d^4, d = r - p1  ->  class2: r0 K2 K3 K4
      return {"class2", joinValues(std::array{p[1], 0.5 * p[0] * toKCal, p[2] * toKCal / 3.0, 0.25 * p[3] * toKCal}),
              ""};
    case BondType::CFF_Quartic:
      // p0 d^2 + p2 d^3 + p3 d^4  ->  class2: r0 K2 K3 K4
      return {"class2", joinValues(std::array{p[1], p[0] * toKCal, p[2] * toKCal, p[3] * toKCal}), ""};
    case BondType::Fixed:
      return {"zero", fmt(p[0]), "fixed bond length in RASPA: constrain it in LAMMPS (fix shake/rattle)"};
    case BondType::None:
      return {"zero", "", ""};
    default:
      return {"zero", "", std::format("RASPA bond type {} has no LAMMPS equivalent; written as zero",
                                      static_cast<std::size_t>(bond.type))};
  }
}

// ---------------------------------------------------------------------------------------------------
// Bends. Angles in degrees; RASPA's (1/2) k conventions become LAMMPS's K.
// ---------------------------------------------------------------------------------------------------
CoefficientLine bendCoefficients(const BendPotential &bend)
{
  const double toKCal = Units::EnergyToKCalPerMol;
  const double toDeg = Units::RadiansToDegrees;
  const std::array<double, maximumNumberOfBendParameters> &p = bend.parameters;
  switch (bend.type)
  {
    case BendType::Harmonic:
    case BendType::CoreShell:
      // (1/2) p0 (theta - p1)^2  ->  harmonic: K theta0
      return {"harmonic", joinValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg}), ""};
    case BendType::Quartic:
      // (1/2) p0 d^2 + (1/3) p2 d^3 + (1/4) p3 d^4  ->  quartic: theta0 K2 K3 K4
      return {"quartic",
              joinValues(std::array{p[1] * toDeg, 0.5 * p[0] * toKCal, p[2] * toKCal / 3.0, 0.25 * p[3] * toKCal}),
              ""};
    case BendType::CFF_Quartic:
      // p0 d^2 + p2 d^3 + p3 d^4  ->  quartic: theta0 K2 K3 K4
      return {"quartic", joinValues(std::array{p[1] * toDeg, p[0] * toKCal, p[2] * toKCal, p[3] * toKCal}), ""};
    case BendType::HarmonicCosine:
      // (1/2) p0 (cos theta - cos theta0)^2, p1 stored as cos theta0  ->  cosine/squared: K theta0
      return {"cosine/squared",
              joinValues(std::array{0.5 * p[0] * toKCal, std::acos(std::clamp(p[1], -1.0, 1.0)) * toDeg}), ""};
    case BendType::Fixed:
      return {"zero", fmt(p[0] * toDeg), "fixed bend angle in RASPA: constrain it in LAMMPS (fix shake/rattle)"};
    case BendType::Rigid:
      return {"zero", "", "bend inside a rigid fragment: keep the fragment rigid in LAMMPS (fix rigid)"};
    default:
      return {"zero", "", std::format("RASPA bend type {} has no LAMMPS equivalent; written as zero",
                                      static_cast<std::size_t>(bend.type))};
  }
}

// ---------------------------------------------------------------------------------------------------
// Torsions. RASPA and LAMMPS share the dihedral convention (phi = 180 degrees for trans), so every
// cosine series maps exactly. A Fourier series sum_n a_n cos(n phi) is written as a polynomial in
// cos(phi) for 'dihedral_style nharmonic' (E = sum_i A_i cos^(i-1) phi) via the Chebyshev identities
// cos(n phi) = T_n(cos phi); this keeps constant terms, so the energy matches RASPA's exactly.
// ---------------------------------------------------------------------------------------------------
std::vector<double> chebyshevToPolynomial(std::span<const double> fourier)
{
  // T_n coefficient tables in powers of c: T_0 = 1, T_1 = c, T_{n+1} = 2 c T_n - T_{n-1}.
  std::vector<double> polynomial(fourier.size(), 0.0);
  std::vector<double> previous{1.0};      // T_0
  std::vector<double> current{0.0, 1.0};  // T_1
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

CoefficientLine nharmonicLine(std::span<const double> fourierKelvinToEnergy)
{
  std::vector<double> fourier(fourierKelvinToEnergy.begin(), fourierKelvinToEnergy.end());
  for (double &a : fourier) a *= Units::EnergyToKCalPerMol;
  const std::vector<double> polynomial = chebyshevToPolynomial(fourier);
  return {"nharmonic", std::format("{} {}", polynomial.size(), joinValues(polynomial)), ""};
}

CoefficientLine torsionCoefficients(const TorsionPotential &torsion)
{
  const double toKCal = Units::EnergyToKCalPerMol;
  const double toDeg = Units::RadiansToDegrees;
  const std::array<double, maximumNumberOfTorsionParameters> &p = torsion.parameters;
  switch (torsion.type)
  {
    case TorsionType::Harmonic:
      // (1/2) p0 (phi - p1)^2  ->  quadratic: K phi0
      return {"quadratic", joinValues(std::array{0.5 * p[0] * toKCal, p[1] * toDeg}), ""};
    case TorsionType::HarmonicCosine:
    {
      // (1/2) p0 (cos phi - c0)^2, p1 stored as c0  ->  polynomial (1/2) p0 c0^2 - p0 c0 c + (1/2) p0 c^2
      const double k = p[0] * toKCal, c0 = p[1];
      return {"nharmonic", std::format("3 {}", joinValues(std::array{0.5 * k * c0 * c0, -k * c0, 0.5 * k})), ""};
    }
    case TorsionType::ThreeCosine:
    case TorsionType::MM3:
      // (1/2) p0 (1 + cos phi) + (1/2) p1 (1 - cos 2phi) + (1/2) p2 (1 + cos 3phi)
      return nharmonicLine(std::array{0.5 * (p[0] + p[1] + p[2]), 0.5 * p[0], -0.5 * p[1], 0.5 * p[2]});
    case TorsionType::RyckaertBellemans:
    {
      // sum_i (-1)^i p_i cos^i phi (RASPA evaluates the alternating-sign form)
      std::vector<double> polynomial(6);
      for (std::size_t i = 0; i != 6; ++i) polynomial[i] = ((i % 2) ? -1.0 : 1.0) * p[i] * toKCal;
      return {"nharmonic", std::format("6 {}", joinValues(polynomial)), ""};
    }
    case TorsionType::TraPPE:
      // p0 + p1 (1 + cos phi) + p2 (1 - cos 2phi) + p3 (1 + cos 3phi)
      return nharmonicLine(std::array{p[0] + p[1] + p[2] + p[3], p[1], -p[2], p[3]});
    case TorsionType::TraPPE_Extended:
      // sum_{n=0}^{4} p_n cos(n phi)
      return nharmonicLine(std::array{p[0], p[1], p[2], p[3], p[4]});
    case TorsionType::ModifiedTraPPE:
      // TraPPE with a phase shift p4; only the unshifted case is a polynomial in cos(phi)
      if (std::abs(p[4]) < 1e-12)
      {
        return nharmonicLine(std::array{p[0] + p[1] + p[2] + p[3], p[1], -p[2], p[3]});
      }
      return {"zero", "", "MODIFIED_TRAPPE with a non-zero phase shift has no LAMMPS equivalent; written as zero"};
    case TorsionType::CVFF:
    {
      // p0 (1 + cos(p1 phi - p2))  ->  charmm: K n d(degrees) w, with integer n and d
      const long n = std::lround(p[1]);
      const long d = std::lround(p[2] * toDeg);
      std::string note{};
      if (std::abs(p[1] - static_cast<double>(n)) > 1e-9 || std::abs(p[2] * toDeg - static_cast<double>(d)) > 1e-6)
      {
        note = "charmm needs integer n and d (degrees); values were rounded";
      }
      return {"charmm", std::format("{} {} {} 0.0", fmt(p[0] * toKCal), n, d), note};
    }
    case TorsionType::CFF:
      // p0 (1 - cos phi) + p1 (1 - cos 2phi) + p2 (1 - cos 3phi)
      return nharmonicLine(std::array{p[0] + p[1] + p[2], -p[0], -p[1], -p[2]});
    case TorsionType::CFF2:
      // p0 (1 + cos phi) + p1 (1 + cos 2phi) + p2 (1 + cos 3phi)
      return nharmonicLine(std::array{p[0] + p[1] + p[2], p[0], p[1], p[2]});
    case TorsionType::OPLS:
      // (1/2) p0 + (1/2) p1 (1 + cos phi) + (1/2) p2 (1 - cos 2phi) + (1/2) p3 (1 + cos 3phi)
      return nharmonicLine(std::array{0.5 * (p[0] + p[1] + p[2] + p[3]), 0.5 * p[1], -0.5 * p[2], 0.5 * p[3]});
    case TorsionType::FourierSeries:
      // (1/2) [p0 (1 + cos phi) + p1 (1 - cos 2phi) + p2 (1 + cos 3phi) + p3 (1 - cos 4phi) + p4 (1 + cos 5phi)
      //        + p5 (1 + cos 6phi)]
      return nharmonicLine(std::array{0.5 * (p[0] + p[1] + p[2] + p[3] + p[4] + p[5]), 0.5 * p[0], -0.5 * p[1],
                                      0.5 * p[2], -0.5 * p[3], 0.5 * p[4], 0.5 * p[5]});
    case TorsionType::FourierSeries2:
      // as FourierSeries but with (1 + cos 4phi)
      return nharmonicLine(std::array{0.5 * (p[0] + p[1] + p[2] + p[3] + p[4] + p[5]), 0.5 * p[0], -0.5 * p[1],
                                      0.5 * p[2], 0.5 * p[3], 0.5 * p[4], 0.5 * p[5]});
    case TorsionType::Fixed:
      return {"zero", "", "fixed torsion in RASPA: constrain it in LAMMPS"};
    case TorsionType::CVFFBlocked:
      return {"zero", "", "CVFF_BLOCKED carries no energy"};
    default:
      return {"zero", "", std::format("RASPA torsion type {} has no LAMMPS equivalent; written as zero",
                                      static_cast<std::size_t>(torsion.type))};
  }
}

// ---------------------------------------------------------------------------------------------------
// A coefficient section: the collected lines of one class over all components, in the same order as
// the per-molecule type numbering of the corresponding topology section.
// ---------------------------------------------------------------------------------------------------
struct CoefficientSection
{
  std::vector<CoefficientLine> lines{};

  std::vector<std::string> styles() const
  {
    std::vector<std::string> result{};
    for (const CoefficientLine &line : lines)
    {
      if (std::find(result.begin(), result.end(), line.style) == result.end()) result.push_back(line.style);
    }
    return result;
  }

  bool hybrid() const { return styles().size() > 1; }

  std::string styleCommand(std::string_view command) const
  {
    const std::vector<std::string> list = styles();
    if (list.empty()) return std::format("{} none", command);
    if (list.size() == 1) return std::format("{} {}", command, list.front());
    std::string out = std::format("{} hybrid", command);
    for (const std::string &style : list) out += " " + style;
    return out;
  }

  void write(std::ostringstream &out, std::string_view title) const
  {
    if (lines.empty()) return;
    std::print(out, "\n{}\n\n", title);
    const bool useHybrid = hybrid();
    for (std::size_t i = 0; i != lines.size(); ++i)
    {
      const CoefficientLine &line = lines[i];
      std::print(out, "{}", i + 1);
      if (useHybrid) std::print(out, " {}", line.style);
      if (!line.coefficients.empty()) std::print(out, " {}", line.coefficients);
      if (!line.note.empty()) std::print(out, "  # {}", line.note);
      std::print(out, "\n");
    }
  }
};

// The 1-4 scaling a component uses for its intramolecular van der Waals and Coulomb pairs, read off
// the pair lists (a pair whose two atoms are the ends of a torsion is a 1-4 pair). Nullopt when the
// component has no 1-4 pair at all.
struct Scaling14
{
  double vanDerWaals{0.0};
  double coulomb{0.0};
};

Scaling14 scaling14(const Component &component)
{
  Scaling14 result{};
  const std::vector<std::array<std::size_t, 4>> torsions = component.connectivityTable.findAllTorsions();
  auto is14 = [&](std::size_t A, std::size_t B)
  {
    for (const std::array<std::size_t, 4> &t : torsions)
    {
      if ((t[0] == A && t[3] == B) || (t[0] == B && t[3] == A)) return true;
    }
    return false;
  };
  for (const VanDerWaalsPotential &pair : component.intraMolecularPotentials.vanDerWaals)
  {
    if (is14(pair.identifiers[0], pair.identifiers[1]))
    {
      result.vanDerWaals = pair.scaling;
      break;
    }
  }
  for (const CoulombPotential &pair : component.intraMolecularPotentials.coulombs)
  {
    if (is14(pair.identifiers[0], pair.identifiers[1]))
    {
      result.coulomb = pair.scaling;
      break;
    }
  }
  return result;
}
}  // namespace

std::string IO::WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                    std::span<const AtomDynamics> atomDynamics,
                                    std::span<const Molecule> moleculeData, const SimulationBox simulationBox,
                                    const ForceField forceField,
                                    std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                    std::optional<Framework> framework)
{
  std::ostringstream out;

  // LAMMPS 'real' units: Angstrom, kcal/mol, fs, e, g/mol.
  const double toAngstrom = Units::LengthConversionFactor * 1e10;
  const double toAngstromPerFemtosecond = Units::VelocityConversionFactor * 1e-5;
  const double toKCal = Units::EnergyToKCalPerMol;

  // ------------------------------------------------------------------------------------------------
  // Coefficient sections (one type per bonded term of each component, matching the topology sections).
  // ------------------------------------------------------------------------------------------------
  CoefficientSection bondSection{}, bendSection{}, torsionSection{};
  std::size_t numberOfBonds = 0uz, numberOfBends = 0uz, numberOfTorsions = 0uz;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    const Potentials::IntraMolecularPotentials &intra = components[i].intraMolecularPotentials;
    numberOfBonds += intra.bonds.size() * numberOfIntegerMoleculesPerComponent[i];
    numberOfBends += intra.bends.size() * numberOfIntegerMoleculesPerComponent[i];
    numberOfTorsions += intra.torsions.size() * numberOfIntegerMoleculesPerComponent[i];
    for (const BondPotential &bond : intra.bonds) bondSection.lines.push_back(bondCoefficients(bond));
    for (const BendPotential &bend : intra.bends) bendSection.lines.push_back(bendCoefficients(bend));
    for (const TorsionPotential &torsion : intra.torsions) torsionSection.lines.push_back(torsionCoefficients(torsion));
  }

  // ------------------------------------------------------------------------------------------------
  // Header comment: the input-script settings a data file cannot carry, derived from the force field.
  // (read_data skips the first line and treats '#' lines as blank.)
  // ------------------------------------------------------------------------------------------------
  std::print(out, "LAMMPS data file written by RASPA3\n\n");
  std::print(out, "# Companion input-script settings derived from the RASPA force field (a data file cannot carry\n");
  std::print(out, "# them); paste before read_data:\n");
  std::print(out, "#   units real\n");
  std::print(out, "#   atom_style full\n");
  std::print(out, "#   {}\n", bondSection.styleCommand("bond_style"));
  std::print(out, "#   {}\n", bendSection.styleCommand("angle_style"));
  std::print(out, "#   {}\n", torsionSection.styleCommand("dihedral_style"));

  const double cutOffVDW = forceField.cutOffMoleculeVDW;
  if (forceField.useCharge)
  {
    if (forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      std::print(out, "#   pair_style lj/cut/coul/long {} {}\n", fmt(cutOffVDW), fmt(forceField.cutOffCoulomb));
      std::print(out, "#   kspace_style ewald {:.1e}\n", forceField.EwaldPrecision);
    }
    else
    {
      std::print(out, "#   pair_style lj/cut/coul/cut {} {}   # RASPA charge method {} is not Ewald; choose the\n",
                 fmt(cutOffVDW), fmt(forceField.cutOffCoulomb), static_cast<int>(forceField.chargeMethod));
      std::print(out, "#                                        # matching LAMMPS coul/* variant\n");
    }
  }
  else
  {
    std::print(out, "#   pair_style lj/cut {}\n", fmt(cutOffVDW));
  }
  {
    const bool anyTail = std::any_of(forceField.tailCorrections.begin(), forceField.tailCorrections.end(),
                                     [](bool b) { return b; });
    const bool anyShift =
        std::any_of(forceField.shiftPotentials.begin(), forceField.shiftPotentials.end(), [](bool b) { return b; });
    if (anyTail || anyShift)
    {
      std::print(out, "#   pair_modify{}{}\n", anyTail ? " tail yes" : "", anyShift ? " shift yes" : "");
    }
  }
  if (framework.has_value() && std::abs(forceField.cutOffFrameworkVDW - cutOffVDW) > 1e-10)
  {
    std::print(out, "#   # RASPA uses a framework-molecule VDW cut-off of {} A; LAMMPS pair_style has one cut-off\n",
               fmt(forceField.cutOffFrameworkVDW));
  }

  // special_bonds from the components' 1-4 scaling (1-2 and 1-3 are always excluded in RASPA; from 1-5 on
  // the full interaction acts). LAMMPS applies one setting to all molecules.
  {
    std::optional<Scaling14> common{};
    bool consistent = true;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      if (components[i].intraMolecularPotentials.bonds.empty()) continue;
      const Scaling14 s = scaling14(components[i]);
      if (!common.has_value())
      {
        common = s;
      }
      else if (std::abs(common->vanDerWaals - s.vanDerWaals) > 1e-12 || std::abs(common->coulomb - s.coulomb) > 1e-12)
      {
        consistent = false;
      }
    }
    if (common.has_value())
    {
      std::print(out, "#   special_bonds lj 0 0 {} coul 0 0 {}\n", fmt(common->vanDerWaals),
                 fmt(common->coulomb));
      if (!consistent)
      {
        std::print(out, "#   # WARNING: the components use different 1-4 scalings; LAMMPS special_bonds is global\n");
      }
    }
  }
  std::print(out, "# Cross pair terms are listed in 'PairIJ Coeffs' (no pair_modify mix needed).\n");
  {
    bool nonLJ = false;
    for (std::size_t i = 0; i < forceField.numberOfPseudoAtoms; ++i)
    {
      for (std::size_t j = i; j < forceField.numberOfPseudoAtoms; ++j)
      {
        const VDWParameters::Type type = forceField(i, j).type;
        if (type != VDWParameters::Type::LennardJones && type != VDWParameters::Type::None) nonLJ = true;
      }
    }
    if (nonLJ)
    {
      std::print(out, "# WARNING: some pair potentials are not Lennard-Jones; their lines list the raw p0 (kcal/mol) "
                      "and p1 and\n#          need a matching pair_style.\n");
    }
  }
  std::print(out, "\n");

  // ------------------------------------------------------------------------------------------------
  // Header.
  // ------------------------------------------------------------------------------------------------
  std::print(out, "{} atoms\n", atomData.size());
  std::print(out, "{} bonds\n", numberOfBonds);
  std::print(out, "{} angles\n", numberOfBends);
  std::print(out, "{} dihedrals\n\n", numberOfTorsions);

  std::print(out, "{} atom types\n", forceField.numberOfPseudoAtoms);
  std::print(out, "{} bond types\n", bondSection.lines.size());
  std::print(out, "{} angle types\n", bendSection.lines.size());
  std::print(out, "{} dihedral types\n", torsionSection.lines.size());
  std::print(out, "0 improper types\n\n");

  double3 lengths = simulationBox.lengths() * toAngstrom;
  double3 angles = simulationBox.angles();

  double xy = lengths.y * std::cos(angles.z);
  double xz = lengths.z * std::cos(angles.y);
  double yz = (lengths.y * lengths.z * std::cos(angles.x) - xy * xz) / (lengths.y * std::sin(angles.z));

  // LAMMPS box edge lengths are the projections lx, ly, lz for a triclinic cell.
  const bool triclinic = (std::abs(xy) > 1e-10) || (std::abs(xz) > 1e-10) || (std::abs(yz) > 1e-10);
  const double ly = triclinic ? std::sqrt(lengths.y * lengths.y - xy * xy) : lengths.y;
  const double lz = triclinic ? std::sqrt(lengths.z * lengths.z - xz * xz - yz * yz) : lengths.z;
  std::print(out, "{} {} xlo xhi\n", 0.0, fmt(lengths.x));
  std::print(out, "{} {} ylo yhi\n", 0.0, fmt(ly));
  std::print(out, "{} {} zlo zhi\n", 0.0, fmt(lz));
  if (triclinic)
  {
    std::print(out, "{} {} {} xy xz yz\n", fmt(xy), fmt(xz), fmt(yz));
  }

  // ------------------------------------------------------------------------------------------------
  // Masses and pair coefficients (all i <= j, so the file is exact for any RASPA mixing rule or
  // explicitly overridden cross term).
  // ------------------------------------------------------------------------------------------------
  std::print(out, "\nMasses\n\n");
  std::size_t idx = 1;  // lammps works with 1-indexing
  for (const PseudoAtom &pseudoAtom : forceField.pseudoAtoms)
  {
    std::print(out, "  {} {}  # {}\n", idx++, fmt(pseudoAtom.mass), pseudoAtom.name);
  }

  std::print(out, "\nPairIJ Coeffs\n\n");
  for (std::size_t i = 0; i < forceField.numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = i; j < forceField.numberOfPseudoAtoms; ++j)
    {
      const VDWParameters &vdw = forceField(i, j);
      // Lennard-Jones: p0 = epsilon (energy), p1 = sigma (Angstrom). Other types: raw p0 (in kcal/mol), p1.
      std::print(out, "  {} {} {} {}\n", i + 1, j + 1, fmt(vdw.parameters.x * toKCal), fmt(vdw.parameters.y));
    }
  }

  bondSection.write(out, "Bond Coeffs");
  bendSection.write(out, "Angle Coeffs");
  torsionSection.write(out, "Dihedral Coeffs");

  // ------------------------------------------------------------------------------------------------
  // Atoms (atom_style full: id molecule type q x y z). Framework atoms, when present, form molecule 1
  // and the component molecules follow; without a framework the molecules start at 1.
  // ------------------------------------------------------------------------------------------------
  std::vector<std::size_t> molAtomOffset(numberOfIntegerMoleculesPerComponent.size());
  molAtomOffset[0] = 0;
  for (std::size_t i = 1; i < components.size(); ++i)
  {
    molAtomOffset[i] = molAtomOffset[i - 1] + numberOfIntegerMoleculesPerComponent[i - 1];
  }

  std::print(out, "\nAtoms\n\n");
  idx = 1;  // lammps works with 1-indexing
  const std::size_t numberOfFrameworkAtoms = (framework.has_value()) ? framework->atoms.size() : 0uz;
  for (const Atom &atom : atomData)
  {
    const bool isFrameworkAtom = framework.has_value() && idx <= numberOfFrameworkAtoms;
    const std::size_t moleculeIndex =
        isFrameworkAtom ? 1uz
                        : atom.moleculeId + molAtomOffset[atom.componentId] + (framework.has_value() ? 2uz : 1uz);
    const double3 position = atom.position * toAngstrom;
    std::print(out, "  {} {} {} {} {} {} {}\n", idx++, moleculeIndex, atom.type + 1, fmt(atom.charge),
               fmt(position.x), fmt(position.y), fmt(position.z));
  }

  std::print(out, "\nVelocities\n\n");
  idx = 1;  // lammps works with 1-indexing
  for (std::size_t i = 0; i < atomData.size(); ++i)
  {
    const Atom &atom = atomData[i];
    double3 velocity{};
    if (!components[atom.componentId].rigid)
    {
      velocity = atomDynamics[i].velocity;
    }
    else
    {
      velocity = moleculeData[atom.moleculeId].velocity;
    }
    velocity = velocity * toAngstromPerFemtosecond;
    std::print(out, "  {} {} {} {}\n", idx++, fmt(velocity.x), fmt(velocity.y), fmt(velocity.z));
  }

  // ------------------------------------------------------------------------------------------------
  // Topology. Type numbering: term b of component i has type (number of terms of the preceding
  // components) + b + 1, the same order in which the coefficient sections were filled.
  // ------------------------------------------------------------------------------------------------
  if (numberOfBonds)
  {
    std::print(out, "\nBonds\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t bondIndex = 1;
    std::size_t bondCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t b = 0; b < components[i].intraMolecularPotentials.bonds.size(); b++)
        {
          auto &bond = components[i].intraMolecularPotentials.bonds[b];
          std::print(out, "  {} {} {} {}\n", bondIndex++, bondCoeffCount + b, bond.identifiers[0] + idx,
                     bond.identifiers[1] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      bondCoeffCount += components[i].intraMolecularPotentials.bonds.size();
    }
  }

  if (numberOfBends)
  {
    std::print(out, "\nAngles\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t bendIndex = 1;
    std::size_t bendCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t b = 0; b < components[i].intraMolecularPotentials.bends.size(); b++)
        {
          auto &bend = components[i].intraMolecularPotentials.bends[b];
          std::print(out, "  {} {} {} {} {}\n", bendIndex++, bendCoeffCount + b, bend.identifiers[0] + idx,
                     bend.identifiers[1] + idx, bend.identifiers[2] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      bendCoeffCount += components[i].intraMolecularPotentials.bends.size();
    }
  }

  if (numberOfTorsions)
  {
    std::print(out, "\nDihedrals\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t torsionIndex = 1;
    std::size_t torsionCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t t = 0; t < components[i].intraMolecularPotentials.torsions.size(); t++)
        {
          auto &torsion = components[i].intraMolecularPotentials.torsions[t];
          std::print(out, "  {} {} {} {} {} {}\n", torsionIndex++, torsionCoeffCount + t, torsion.identifiers[0] + idx,
                     torsion.identifiers[1] + idx, torsion.identifiers[2] + idx, torsion.identifiers[3] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      torsionCoeffCount += components[i].intraMolecularPotentials.torsions.size();
    }
  }

  return out.str();
}
