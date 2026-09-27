module;

module lammps_io;

import std;

import double3;
import int3;
import component;
import atom;
import atom_dynamics;
import simulationbox;
import forcefield;
import framework;
import molecule;
import lammps_styles;
import lammps_topology;

namespace
{
using namespace LAMMPS;

std::string styleArguments(std::string_view style, const Topology &topology, std::string_view bondedClass)
{
  if (style != "table") return std::string(style);
  const TableGrid &grid = topology.options.grid;
  if (bondedClass == "bond") return std::format("table linear {}", grid.bondPoints);
  if (bondedClass == "angle") return std::format("table linear {}", grid.anglePoints);
  if (bondedClass == "dihedral") return std::format("table linear {}", grid.dihedralPoints);
  return "table";
}

/// 'bond_style harmonic', 'angle_style hybrid harmonic table linear 1801', ...
std::string styleCommand(std::string_view command, std::string_view bondedClass, const TypeTable &table,
                         const Topology &topology)
{
  const std::vector<std::string> styles = table.styles();
  if (styles.empty()) return std::format("{} none", command);
  if (styles.size() == 1) return std::format("{} {}", command, styleArguments(styles.front(), topology, bondedClass));
  std::string out = std::format("{} hybrid", command);
  for (const std::string &style : styles) out += " " + styleArguments(style, topology, bondedClass);
  return out;
}

/// The coefficient line of one type: "[style] coefficients" (style only under hybrid).
std::string coefficientLine(const Term &term, bool hybrid, const Topology &topology)
{
  std::string out{};
  if (hybrid) out += term.style;
  if (term.style == "table")
  {
    out += std::format("{}{} {}", hybrid ? " " : "", topology.options.tableFile, term.tableKeyword);
  }
  else if (!term.coefficients.empty())
  {
    out += std::format("{}{}", hybrid ? " " : "", term.coefficients);
  }
  return out;
}

void writeCoefficientSection(std::ostringstream &out, std::string_view title, const TypeTable &table,
                             const Topology &topology)
{
  if (table.empty()) return;
  std::print(out, "\n{}\n\n", title);
  const bool hybrid = table.hybrid();
  for (std::size_t i = 0; i < table.types.size(); ++i)
  {
    const Term &term = table.types[i];
    std::string line = coefficientLine(term, hybrid, topology);
    std::print(out, "{}{}{}", i + 1, line.empty() ? "" : " ", line);
    if (!term.note.empty()) std::print(out, "  # {}", term.note);
    std::print(out, "\n");
  }
}

void writeTopologySection(std::ostringstream &out, std::string_view title, std::span<const BondedEntry> entries,
                          std::size_t arity)
{
  if (entries.empty()) return;
  std::print(out, "\n{}\n\n", title);
  for (std::size_t i = 0; i < entries.size(); ++i)
  {
    std::print(out, "{} {}", i + 1, entries[i].type);
    for (std::size_t k = 0; k < arity; ++k) std::print(out, " {}", entries[i].atoms[k]);
    std::print(out, "\n");
  }
}

std::string coulombStyle(const Topology &topology, std::vector<std::string> &notes)
{
  const std::string rc = formatValue(topology.cutOffCoulomb);
  switch (topology.chargeMethod)
  {
    case ForceField::ChargeMethod::Ewald:
      return "coul/long " + rc;
    case ForceField::ChargeMethod::Coulomb:
      return "coul/cut " + rc;
    case ForceField::ChargeMethod::Wolf:
      return std::format("coul/wolf {} {}", formatValue(topology.alpha), rc);
    case ForceField::ChargeMethod::DampedShiftedForce:
      return std::format("coul/dsf {} {}", formatValue(topology.alpha), rc);
    default:
      notes.push_back("RASPA charge method has no LAMMPS pair style; coul/cut written instead");
      return "coul/cut " + rc;
  }
}

/// lj/cut + coul/X and buck + coul/X have fused LAMMPS styles; anything else is hybrid/overlay.
std::optional<std::string> fusedPairStyle(std::string_view vdw, const Topology &topology)
{
  const std::string rvdw = formatValue(topology.cutOffVDW);
  const std::string rc = formatValue(topology.cutOffCoulomb);
  const std::string alpha = formatValue(topology.alpha);
  if (vdw == "lj/cut")
  {
    switch (topology.chargeMethod)
    {
      case ForceField::ChargeMethod::Ewald:
        return std::format("lj/cut/coul/long {} {}", rvdw, rc);
      case ForceField::ChargeMethod::Coulomb:
        return std::format("lj/cut/coul/cut {} {}", rvdw, rc);
      case ForceField::ChargeMethod::Wolf:
        return std::format("lj/cut/coul/wolf {} {} {}", alpha, rvdw, rc);
      case ForceField::ChargeMethod::DampedShiftedForce:
        return std::format("lj/cut/coul/dsf {} {} {}", alpha, rvdw, rc);
      default:
        return std::format("lj/cut/coul/cut {} {}", rvdw, rc);
    }
  }
  if (vdw == "buck")
  {
    switch (topology.chargeMethod)
    {
      case ForceField::ChargeMethod::Ewald:
        return std::format("buck/coul/long {} {}", rvdw, rc);
      case ForceField::ChargeMethod::Coulomb:
        return std::format("buck/coul/cut {} {}", rvdw, rc);
      default:
        return std::nullopt;
    }
  }
  return std::nullopt;
}

std::string vdwStyleArguments(std::string_view style, const Topology &topology)
{
  if (style == "table") return std::format("table linear {}", topology.options.grid.pairPoints);
  return std::format("{} {}", style, formatValue(topology.cutOffVDW));
}

/// Compresses sorted ids into LAMMPS range syntax ("1:40 42 50:60").
std::string idRanges(std::vector<std::size_t> ids)
{
  std::sort(ids.begin(), ids.end());
  std::string out{};
  std::size_t i = 0;
  while (i < ids.size())
  {
    std::size_t j = i;
    while (j + 1 < ids.size() && ids[j + 1] == ids[j] + 1) ++j;
    if (!out.empty()) out += " ";
    out += (j > i) ? std::format("{}:{}", ids[i], ids[j]) : std::format("{}", ids[i]);
    i = j + 1;
  }
  return out;
}
}  // namespace

namespace LAMMPS
{
std::string writeDataFile(const Topology &topology)
{
  std::ostringstream out;
  std::print(out, "LAMMPS data file written by RASPA3 (units real, atom_style full); input script: see the companion "
                  ".in file\n\n");

  std::print(out, "{} atoms\n", topology.atoms.size());
  std::print(out, "{} bonds\n", topology.bondList.size());
  std::print(out, "{} angles\n", topology.angleList.size());
  std::print(out, "{} dihedrals\n", topology.dihedralList.size());
  std::print(out, "{} impropers\n\n", topology.improperList.size());

  std::print(out, "{} atom types\n", topology.masses.size());
  std::print(out, "{} bond types\n", topology.bonds.types.size());
  std::print(out, "{} angle types\n", topology.angles.types.size());
  std::print(out, "{} dihedral types\n", topology.dihedrals.types.size());
  std::print(out, "{} improper types\n\n", topology.impropers.types.size());

  std::print(out, "0.0 {} xlo xhi\n", formatValue(topology.lengths.x));
  std::print(out, "0.0 {} ylo yhi\n", formatValue(topology.lengths.y));
  std::print(out, "0.0 {} zlo zhi\n", formatValue(topology.lengths.z));
  if (topology.triclinic)
  {
    std::print(out, "{} {} {} xy xz yz\n", formatValue(topology.xy), formatValue(topology.xz),
               formatValue(topology.yz));
  }

  std::print(out, "\nMasses\n\n");
  for (std::size_t i = 0; i < topology.masses.size(); ++i)
  {
    std::print(out, "{} {}  # {}\n", i + 1, formatValue(topology.masses[i].second), topology.masses[i].first);
  }

  if (topology.pairCoefficientsInDataFile())
  {
    std::print(out, "\nPairIJ Coeffs\n\n");
    for (const PairEntry &pair : topology.pairs)
    {
      std::print(out, "{} {} {}\n", pair.i, pair.j, pair.term.coefficients);
    }
  }

  writeCoefficientSection(out, "Bond Coeffs", topology.bonds, topology);
  writeCoefficientSection(out, "Angle Coeffs", topology.angles, topology);
  writeCoefficientSection(out, "Dihedral Coeffs", topology.dihedrals, topology);
  writeCoefficientSection(out, "Improper Coeffs", topology.impropers, topology);

  std::print(out, "\nAtoms # full\n\n");
  for (const AtomEntry &atom : topology.atoms)
  {
    std::print(out, "{} {} {} {} {} {} {} {} {} {}\n", atom.id, atom.molecule, atom.type, formatValue(atom.charge),
               formatValue(atom.position.x), formatValue(atom.position.y), formatValue(atom.position.z),
               atom.image.x, atom.image.y, atom.image.z);
  }

  std::print(out, "\nVelocities\n\n");
  for (const AtomEntry &atom : topology.atoms)
  {
    std::print(out, "{} {} {} {}\n", atom.id, formatValue(atom.velocity.x), formatValue(atom.velocity.y),
               formatValue(atom.velocity.z));
  }

  writeTopologySection(out, "Bonds", topology.bondList, 2);
  writeTopologySection(out, "Angles", topology.angleList, 3);
  writeTopologySection(out, "Dihedrals", topology.dihedralList, 4);
  writeTopologySection(out, "Impropers", topology.improperList, 4);

  if (topology.numberOfFragments > 0)
  {
    std::print(out, "\nFragments\n\n");
    for (const AtomEntry &atom : topology.atoms) std::print(out, "{} {}\n", atom.id, atom.fragment);
  }

  return out.str();
}

std::string writeTableFile(const Topology &topology)
{
  if (topology.tables.empty()) return {};
  std::ostringstream out;
  std::print(out, "# Tabulated potentials written by RASPA3 for {}\n", topology.options.dataFile);
  for (const auto &[keyword, body] : topology.tables) std::print(out, "\n{}", body);
  return out.str();
}

std::string writePairListFile(const Topology &topology)
{
  if (topology.pairList.empty()) return {};
  std::ostringstream out;
  std::print(out, "# 1-4 pairs with non-uniform scaling for pair_style list (RASPA3)\n");
  for (const std::string &line : topology.pairList) std::print(out, "{}\n", line);
  return out.str();
}

std::string writeInputScript(const Topology &topology)
{
  std::ostringstream out;
  std::vector<std::string> notes = topology.warnings;

  std::print(out, "# LAMMPS input script written by RASPA3 for {}\n", topology.options.dataFile);
  std::print(out, "# Reads the exported configuration and evaluates the energy once (run 0); append your own\n");
  std::print(out, "# integrator after the '# integration' marker.\n\n");

  std::print(out, "units real\natom_style full\nboundary p p p\n\n");

  // ---- bonded styles -------------------------------------------------------------------------------
  std::print(out, "{}\n", styleCommand("bond_style", "bond", topology.bonds, topology));
  std::print(out, "{}\n", styleCommand("angle_style", "angle", topology.angles, topology));
  std::print(out, "{}\n", styleCommand("dihedral_style", "dihedral", topology.dihedrals, topology));
  std::print(out, "{}\n", styleCommand("improper_style", "improper", topology.impropers, topology));

  // ---- pair style ----------------------------------------------------------------------------------
  const std::vector<std::string> vdwStyles = topology.pairStyles();
  const bool anyZeroPair =
      std::any_of(topology.pairs.begin(), topology.pairs.end(), [](const PairEntry &p) { return p.term.style == "zero"; });
  bool overlay = false;
  bool hybrid = false;
  std::string pairStyle{};
  std::string coulomb{};
  if (topology.useCharge)
  {
    coulomb = coulombStyle(topology, notes);
    std::optional<std::string> fused = (vdwStyles.size() == 1 && !anyZeroPair)
                                           ? fusedPairStyle(vdwStyles.front(), topology)
                                           : std::nullopt;
    if (fused.has_value())
    {
      pairStyle = *fused;
    }
    else
    {
      overlay = hybrid = true;
      pairStyle = "hybrid/overlay " + coulomb;
      for (const std::string &style : vdwStyles) pairStyle += " " + vdwStyleArguments(style, topology);
    }
  }
  else
  {
    if (vdwStyles.size() == 1 && !anyZeroPair)
    {
      pairStyle = vdwStyleArguments(vdwStyles.front(), topology);
    }
    else if (vdwStyles.empty())
    {
      pairStyle = std::format("zero {}", formatValue(topology.cutOffVDW));
    }
    else
    {
      hybrid = true;
      pairStyle = "hybrid";
      for (const std::string &style : vdwStyles) pairStyle += " " + vdwStyleArguments(style, topology);
      if (anyZeroPair) pairStyle += std::format(" zero {}", formatValue(topology.cutOffVDW));
    }
  }
  std::print(out, "pair_style {}\n", pairStyle);
  if (topology.special.present)
  {
    std::print(out, "special_bonds lj 0.0 0.0 {} coul 0.0 0.0 {}\n", formatValue(topology.special.vdw14),
               formatValue(topology.special.coul14));
  }
  if (!topology.pairList.empty())
  {
    notes.push_back("scaled 1-4 VDW pairs are supplied through pair_style list (see the .pairs file); add "
                    "'list " + topology.options.pairListFile + " " + formatValue(topology.cutOffVDW) +
                    "' to the hybrid/overlay pair_style and 'pair_coeff * * list'");
  }
  if (topology.hasFramework && std::abs(topology.cutOffFrameworkVDW - topology.cutOffVDW) > 1e-10)
  {
    notes.push_back(std::format("RASPA uses a framework-molecule VDW cut-off of {} A; LAMMPS uses one cut-off ({})",
                                formatValue(topology.cutOffFrameworkVDW), formatValue(topology.cutOffVDW)));
  }

  // ---- data file -----------------------------------------------------------------------------------
  std::print(out, "\n");
  if (topology.numberOfFragments > 0)
  {
    std::print(out, "fix fragments all property/atom i_fragment\n");
    std::print(out, "read_data {} fix fragments NULL Fragments\n", topology.options.dataFile);
  }
  else
  {
    std::print(out, "read_data {}\n", topology.options.dataFile);
  }

  // ---- pair coefficients (when not in the data file) -----------------------------------------------
  if (!topology.pairCoefficientsInDataFile() || hybrid)
  {
    std::print(out, "\n");
    if (overlay)
    {
      std::print(out, "pair_coeff * * {}\n", coulomb.substr(0, coulomb.find(' ')));
    }
    for (const PairEntry &pair : topology.pairs)
    {
      if (pair.term.style == "zero")
      {
        if (hybrid && !overlay) std::print(out, "pair_coeff {} {} zero\n", pair.i, pair.j);
        continue;
      }
      std::string line = std::format("pair_coeff {} {}", pair.i, pair.j);
      if (hybrid) line += " " + pair.term.style;
      if (pair.term.style == "table")
      {
        line += std::format(" {} {} {}", topology.options.tableFile, pair.tableKeyword,
                            formatValue(topology.cutOffVDW));
      }
      else
      {
        line += " " + pair.term.coefficients;
      }
      if (!pair.term.note.empty()) line += "  # " + pair.term.note;
      std::print(out, "{}\n", line);
    }
  }
  {
    // The mixing rule is informative only (every pair is written explicitly), but it lets a reader recover it.
    const char *mix = topology.mixingRule == ForceField::MixingRule::Lorentz_Berthelot ? "arithmetic"
                      : topology.mixingRule == ForceField::MixingRule::SixthPower     ? "sixthpower"
                                                                                       : "geometric";
    std::print(out, "pair_modify mix {}{}{}\n", mix, topology.anyTail ? " tail yes" : "",
               topology.anyShift ? " shift yes" : "");
  }
  if (topology.useCharge && topology.chargeMethod == ForceField::ChargeMethod::Ewald)
  {
    std::print(out, "kspace_style ewald {:.1e}\n", topology.ewaldPrecision);
  }

  // ---- class2 cross coefficients ---------------------------------------------------------------------
  auto writeExtras = [&](std::string_view command, const TypeTable &table, std::span<const std::string_view> keys)
  {
    for (std::size_t i = 0; i < table.types.size(); ++i)
    {
      const Term &term = table.types[i];
      if (term.extra.empty()) continue;
      for (std::string_view key : keys)
      {
        auto it = term.extra.find(std::string(key));
        if (it == term.extra.end()) continue;
        std::print(out, "{} {}{} {} {}\n", command, i + 1, table.hybrid() ? " class2" : "", key, it->second);
      }
    }
  };
  {
    bool any = false;
    for (const TypeTable *t : {&topology.angles, &topology.dihedrals, &topology.impropers})
    {
      any = any || std::any_of(t->types.begin(), t->types.end(), [](const Term &x) { return !x.extra.empty(); });
    }
    if (any) std::print(out, "\n# class2 cross terms\n");
    constexpr std::array<std::string_view, 2> angleKeys{"bb", "ba"};
    constexpr std::array<std::string_view, 5> dihedralKeys{"mbt", "ebt", "at", "aat", "bb13"};
    constexpr std::array<std::string_view, 1> improperKeys{"aa"};
    writeExtras("angle_coeff", topology.angles, angleKeys);
    writeExtras("dihedral_coeff", topology.dihedrals, dihedralKeys);
    writeExtras("improper_coeff", topology.impropers, improperKeys);
  }

  // ---- groups, exclusions, constraints, rigid bodies --------------------------------------------------
  std::print(out, "\n# groups\n");
  for (const MoleculeGroup &group : topology.groups)
  {
    if (group.firstMolecule == group.lastMolecule)
      std::print(out, "group {} molecule {}\n", group.name, group.firstMolecule);
    else
      std::print(out, "group {} molecule {}:{}\n", group.name, group.firstMolecule, group.lastMolecule);
  }
  for (const MoleculeGroup &group : topology.groups)
  {
    if (group.rigid)
    {
      // RASPA has no intramolecular interactions inside a rigid molecule
      std::print(out, "neigh_modify exclude molecule/intra {}\n", group.name);
    }
  }

  std::string mobile = "all";
  if (topology.hasFixedAtoms)
  {
    std::vector<std::size_t> frozenIds{};
    for (const AtomEntry &atom : topology.atoms)
    {
      if (atom.fixed) frozenIds.push_back(atom.id);
    }
    std::print(out, "\n# frozen framework atoms\n");
    std::print(out, "group frozen id {}\n", idRanges(frozenIds));
    std::print(out, "group mobile subtract all frozen\n");
    std::print(out, "neigh_modify exclude group frozen frozen\n");
    std::print(out, "velocity frozen set 0.0 0.0 0.0\n");
    std::print(out, "fix freeze frozen setforce 0.0 0.0 0.0\n");
    mobile = "mobile";
  }

  if (!topology.shakeBondTypes.empty() || !topology.shakeAngleTypes.empty())
  {
    std::string line = std::format("fix constraints {} shake 1.0e-6 100 0", mobile);
    if (!topology.shakeBondTypes.empty())
    {
      line += " b";
      for (std::size_t t : topology.shakeBondTypes) line += std::format(" {}", t);
    }
    if (!topology.shakeAngleTypes.empty())
    {
      line += " a";
      for (std::size_t t : topology.shakeAngleTypes) line += std::format(" {}", t);
    }
    std::print(out, "\n# fixed bond lengths / bend angles of RASPA\n{}\n", line);
  }

  std::string integrate = mobile;
  if (topology.numberOfFragments > 0)
  {
    std::print(out, "\n# rigid bodies (RASPA rigid molecules and rigid fragments)\n");
    std::print(out, "variable inbody atom \"i_fragment > 0\"\n");
    std::print(out, "group bodies variable inbody\n");
    std::print(out, "group flexible subtract {} bodies\n", mobile);
    std::print(out, "fix bodies bodies rigid/small custom i_fragment\n");
    integrate = "flexible";
  }

  // ---- diagnostics and run --------------------------------------------------------------------------
  if (!notes.empty())
  {
    std::print(out, "\n# notes from the RASPA export:\n");
    for (const std::string &note : notes) std::print(out, "#   {}\n", note);
  }
  std::print(out, "\nthermo_style custom step temp pe ebond eangle edihed eimp evdwl ecoul elong etail press\n");
  std::print(out, "thermo_modify format float %.10g\n");
  std::print(out, "thermo 1\n");
  std::print(out, "run 0\n");
  std::print(out, "\n# integration (example)\n");
  std::print(out, "# timestep 1.0\n");
  std::print(out, "# fix integrate {} nvt temp 300.0 300.0 100.0\n", integrate);
  std::print(out, "# run 100000\n");

  return out.str();
}

ExportFiles exportSystem(std::span<const Component> components, std::span<const Atom> atomData,
                         std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                         const SimulationBox &simulationBox, const ForceField &forceField,
                         std::span<const std::size_t> numberOfIntegerMoleculesPerComponent,
                         const std::optional<Framework> &framework, ExportOptions options)
{
  Topology topology = buildTopology(components, atomData, atomDynamics, moleculeData, simulationBox, forceField,
                                    numberOfIntegerMoleculesPerComponent, framework, std::move(options));
  ExportFiles files{};
  files.data = writeDataFile(topology);
  files.input = writeInputScript(topology);
  files.table = writeTableFile(topology);
  files.pairList = writePairListFile(topology);
  files.warnings = topology.warnings;
  return files;
}
}  // namespace LAMMPS

std::string IO::WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                    std::span<const AtomDynamics> atomDynamics,
                                    std::span<const Molecule> moleculeData, const SimulationBox simulationBox,
                                    const ForceField forceField,
                                    std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                    std::optional<Framework> framework)
{
  return LAMMPS::exportSystem(components, atomData, atomDynamics, moleculeData, simulationBox, forceField,
                              numberOfIntegerMoleculesPerComponent, framework)
      .data;
}
