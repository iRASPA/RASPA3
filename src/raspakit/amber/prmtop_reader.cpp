module;

module amber_prmtop_reader;

import std;

import json;
import units;
import skelement;

namespace
{
constexpr double amberChargeUnit = 18.2223;  ///< prmtop CHARGE = q [e] * 18.2223 (sqrt of 332.0522 kcal Å/mol)
const double chamberChargeUnit = std::sqrt(332.0716);  ///< CHAMBER prmtop CHARGE = q [e] * sqrt(CCELEC)

double kelvinPerKCal() { return Units::KCalPerMolToEnergy * Units::EnergyToKelvin; }

std::string trim(std::string_view text)
{
  std::size_t begin = text.find_first_not_of(" \t\r\n");
  if (begin == std::string_view::npos) return {};
  std::size_t end = text.find_last_not_of(" \t\r\n");
  return std::string(text.substr(begin, end - begin + 1));
}

std::string toUpper(std::string text)
{
  for (char& c : text) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return text;
}

std::string sanitizeName(std::string name)
{
  std::string result{};
  for (char c : name)
  {
    if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')
      result += c;
    else if (c == '+')
      result += "_plus";
    else if (c == '*')
      result += "_star";
    else if (c == '\'')
      result += "_prime";
    else
      result += '_';
  }
  if (result.empty()) result = "X";
  return result;
}

// ---------------------------------------------------------------------------------------------------
// prmtop sections
// ---------------------------------------------------------------------------------------------------

struct Section
{
  char kind{'a'};         // a (string), I (integer), E/F/G (floating point)
  std::size_t width{};    // field width
  std::size_t perLine{};  // fields per line
  std::vector<std::string> fields{};
};

std::map<std::string, Section> readSections(const std::filesystem::path& file)
{
  std::ifstream stream(file);
  if (!stream) throw std::runtime_error(std::format("[prmtop reader] cannot open '{}'", file.string()));

  std::map<std::string, Section> sections{};
  std::string line{};
  std::string current{};
  bool firstLine = true;
  while (std::getline(stream, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (firstLine)
    {
      firstLine = false;
      if (!line.starts_with("%VERSION"))
        throw std::runtime_error(
            std::format("[prmtop reader] '{}' does not start with %VERSION: only AMBER7-format prmtop files are "
                        "supported",
                        file.string()));
      continue;
    }
    if (line.starts_with("%FLAG"))
    {
      current = trim(line.substr(5));
      sections[current] = Section{};
      continue;
    }
    if (line.starts_with("%FORMAT"))
    {
      // %FORMAT(20a4), %FORMAT(10I8), %FORMAT(5E16.8), %FORMAT(3I8), %FORMAT(1a80)
      std::size_t open = line.find('(');
      std::size_t close = line.find(')', open == std::string::npos ? 0 : open);
      if (open == std::string::npos || close == std::string::npos || current.empty())
        throw std::runtime_error(std::format("[prmtop reader] malformed format line '{}'", line));
      std::string spec = line.substr(open + 1, close - open - 1);
      std::erase_if(spec, [](char c) { return c == '(' || c == ')' || c == ' '; });  // CHAMBER: 5(E16.8)
      Section& section = sections[current];
      std::size_t position = 0;
      while (position < spec.size() && std::isdigit(static_cast<unsigned char>(spec[position]))) ++position;
      section.perLine = position == 0 ? 1 : static_cast<std::size_t>(std::stoul(spec.substr(0, position)));
      if (position >= spec.size()) throw std::runtime_error(std::format("[prmtop reader] malformed format '{}'", spec));
      section.kind = static_cast<char>(std::toupper(static_cast<unsigned char>(spec[position])));
      if (section.kind == 'A') section.kind = 'a';
      ++position;
      std::size_t widthEnd = position;
      while (widthEnd < spec.size() && std::isdigit(static_cast<unsigned char>(spec[widthEnd]))) ++widthEnd;
      section.width = static_cast<std::size_t>(std::stoul(spec.substr(position, widthEnd - position)));
      continue;
    }
    if (line.starts_with("%COMMENT") || line.starts_with("%")) continue;
    if (current.empty()) continue;

    Section& section = sections[current];
    if (section.width == 0) continue;
    if (section.kind == 'a')
    {
      // fixed-width strings, kept untrimmed (a title keeps its blanks; names are trimmed on use)
      for (std::size_t offset = 0; offset < line.size(); offset += section.width)
      {
        std::string field = line.substr(offset, section.width);
        if (trim(field).empty() && offset + section.width >= line.size()) break;  // trailing padding
        section.fields.push_back(field);
      }
    }
    else
    {
      for (std::size_t offset = 0; offset < line.size(); offset += section.width)
      {
        std::string field = trim(line.substr(offset, section.width));
        if (field.empty()) continue;
        section.fields.push_back(field);
      }
    }
  }
  return sections;
}

std::vector<std::int64_t> integers(const std::map<std::string, Section>& sections, const std::string& flag)
{
  std::vector<std::int64_t> result{};
  auto it = sections.find(flag);
  if (it == sections.end()) return result;
  result.reserve(it->second.fields.size());
  for (const std::string& field : it->second.fields)
  {
    try
    {
      result.push_back(std::stoll(field));
    }
    catch (std::exception&)
    {
      throw std::runtime_error(std::format("[prmtop reader] '{}' is not an integer in section {}", field, flag));
    }
  }
  return result;
}

std::vector<double> reals(const std::map<std::string, Section>& sections, const std::string& flag)
{
  std::vector<double> result{};
  auto it = sections.find(flag);
  if (it == sections.end()) return result;
  result.reserve(it->second.fields.size());
  for (std::string field : it->second.fields)
  {
    // Fortran prints 'D' exponents occasionally
    for (char& c : field)
      if (c == 'D' || c == 'd') c = 'E';
    try
    {
      result.push_back(std::stod(field));
    }
    catch (std::exception&)
    {
      throw std::runtime_error(std::format("[prmtop reader] '{}' is not a number in section {}", field, flag));
    }
  }
  return result;
}

std::vector<std::string> strings(const std::map<std::string, Section>& sections, const std::string& flag)
{
  auto it = sections.find(flag);
  if (it == sections.end()) return {};
  std::vector<std::string> result{};
  result.reserve(it->second.fields.size());
  for (const std::string& field : it->second.fields) result.push_back(trim(field));
  return result;
}

std::string joinedText(const std::map<std::string, Section>& sections, const std::string& flag)
{
  auto it = sections.find(flag);
  if (it == sections.end()) return {};
  std::string text{};
  for (const std::string& field : it->second.fields) text += field;
  return trim(text);
}

// ---------------------------------------------------------------------------------------------------
// topology helpers
// ---------------------------------------------------------------------------------------------------

std::string elementSymbol(std::int64_t atomicNumber, double mass)
{
  if (atomicNumber > 0 && static_cast<std::size_t>(atomicNumber) < PredefinedElements::predefinedElements.size())
    return PredefinedElements::predefinedElements[static_cast<std::size_t>(atomicNumber)]._chemicalSymbol;
  // fall back on the nearest mass (hydrogen mass repartitioning moves masses by up to 3 amu, so prefer the
  // explicit atomic number whenever the file carries it)
  std::size_t best = 0;
  double bestDifference = std::numeric_limits<double>::max();
  for (std::size_t z = 1; z < PredefinedElements::predefinedElements.size(); ++z)
  {
    double difference = std::abs(PredefinedElements::predefinedElements[z]._mass - mass);
    if (difference < bestDifference)
    {
      bestDifference = difference;
      best = z;
    }
  }
  return best == 0 ? "C" : PredefinedElements::predefinedElements[best]._chemicalSymbol;
}

struct DihedralTerm
{
  double forceConstant{};  // kcal/mol (V_n / 2 in the AMBER functional form)
  double periodicity{};
  double phase{};  // radians
  bool operator==(const DihedralTerm&) const = default;
};

struct PseudoAtomKey
{
  std::int64_t typeIndex{};
  std::string typeName{};
  double mass{};
  auto operator<=>(const PseudoAtomKey&) const = default;
};

struct PseudoAtomDefinition
{
  PseudoAtomKey key{};
  std::string name{};
  std::string element{};
  std::optional<double> commonCharge{};  // nullopt when the atoms of the pseudo-atom carry different charges
  bool chargeSeen{false};
};

/// A molecule as compared between instances: everything but the coordinates.
struct MoleculeTemplate
{
  std::vector<std::size_t> pseudoAtomIndex{};
  std::vector<double> charges{};
  std::vector<std::string> residues{};               // residue label per atom
  std::vector<std::array<std::int64_t, 3>> bonds{};  // a, b, type (local 0-based atoms, 1-based type)
  std::vector<std::array<std::int64_t, 4>> angles{};
  std::vector<std::array<std::int64_t, 6>> dihedrals{};  // a, b, c, d, type, flags (bit0: no 1-4, bit1: improper)
  std::vector<std::array<std::int64_t, 6>> cmaps{};      // a, b, c, d, e, map (0-based)
  bool operator==(const MoleculeTemplate&) const = default;
};

std::vector<std::vector<std::size_t>> findMolecules(std::size_t numberOfAtoms,
                                                    const std::vector<std::array<std::size_t, 2>>& bonds)
{
  std::vector<std::size_t> parent(numberOfAtoms);
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  auto find = [&](std::size_t a)
  {
    while (parent[a] != a)
    {
      parent[a] = parent[parent[a]];
      a = parent[a];
    }
    return a;
  };
  for (const std::array<std::size_t, 2>& bond : bonds)
  {
    std::size_t ra = find(bond[0]);
    std::size_t rb = find(bond[1]);
    if (ra != rb) parent[std::max(ra, rb)] = std::min(ra, rb);
  }
  std::map<std::size_t, std::vector<std::size_t>> groups{};
  for (std::size_t atom = 0; atom < numberOfAtoms; ++atom) groups[find(atom)].push_back(atom);
  std::vector<std::vector<std::size_t>> molecules{};
  for (auto& [_, atoms] : groups) molecules.push_back(std::move(atoms));
  std::ranges::sort(molecules, {}, [](const std::vector<std::size_t>& m) { return m.front(); });
  return molecules;
}

// Chebyshev polynomials T_n(x) = cos(n phi) for x = cos(phi), n <= 5, as coefficients of x^k.
constexpr std::array<std::array<double, 6>, 6> chebyshev{{{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                          {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                                          {-1.0, 0.0, 2.0, 0.0, 0.0, 0.0},
                                                          {0.0, -3.0, 0.0, 4.0, 0.0, 0.0},
                                                          {1.0, 0.0, -8.0, 0.0, 8.0, 0.0},
                                                          {0.0, 5.0, 0.0, -20.0, 0.0, 16.0}}};

/// Whether an AMBER dihedral term can be folded into a polynomial in cos(phi): integer periodicity <= 5 and
/// a phase of 0 or 180 degrees.
bool isPolynomialTerm(const DihedralTerm& term)
{
  double n = std::round(term.periodicity);
  if (std::abs(term.periodicity - n) > 1e-6 || n < 0.0 || n > 5.0) return false;
  // prmtop files carry 8 significant digits (pi is written as 3.14159400), so allow a few 1e-6 in the phase
  return std::abs(std::sin(term.phase)) < 1e-5;
}

std::array<double, 6> polynomialCoefficients(std::span<const DihedralTerm> terms, double energyFactor)
{
  // sum_i (V_i/2) (1 + cos(gamma_i) T_{n_i}(cos phi)) with cos(gamma_i) = +-1
  std::array<double, 6> coefficients{};
  for (const DihedralTerm& term : terms)
  {
    std::size_t n = static_cast<std::size_t>(std::round(term.periodicity));
    double sign = std::cos(term.phase) < 0.0 ? -1.0 : 1.0;
    double half = term.forceConstant * energyFactor;
    coefficients[0] += half;
    for (std::size_t k = 0; k < 6; ++k) coefficients[k] += half * sign * chebyshev[n][k];
  }
  return coefficients;
}

std::array<double, 3> centreOfMass(std::span<const std::array<double, 3>> positions, std::span<const double> masses)
{
  std::array<double, 3> com{};
  double total = 0.0;
  for (std::size_t i = 0; i < positions.size(); ++i)
  {
    double mass = masses[i] > 0.0 ? masses[i] : 1.0;
    for (std::size_t k = 0; k < 3; ++k) com[k] += mass * positions[i][k];
    total += mass;
  }
  for (double& c : com) c /= total;
  return com;
}

const std::set<std::string> aminoAcids{"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN", "GLU", "GLH",
                                       "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
                                       "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ACE", "NME", "NHE"};
const std::set<std::string> nucleicAcids{"DA", "DC", "DG", "DT", "DA5", "DC5", "DG5", "DT5", "DA3", "DC3", "DG3", "DT3",
                                         "A",  "C",  "G",  "U",  "A5",  "C5",  "G5",  "U5",  "A3",  "C3",  "G3",  "U3"};

std::string baseResidue(std::string label)
{
  // N-/C-terminal variants: NALA, CALA
  if (label.size() == 4 && (label[0] == 'N' || label[0] == 'C') && aminoAcids.contains(label.substr(1)))
    return label.substr(1);
  return label;
}
}  // namespace

namespace AMBER
{
Prmtop parsePrmtop(const std::filesystem::path& file)
{
  std::map<std::string, Section> sections = readSections(file);

  Prmtop prmtop{};
  for (const auto& [flag, _] : sections) prmtop.flags.insert(flag);

  prmtop.title = joinedText(sections, "TITLE");
  if (prmtop.title.empty()) prmtop.title = joinedText(sections, "CTITLE");
  prmtop.pointers = integers(sections, "POINTERS");
  if (prmtop.pointers.size() < 31)
    throw std::runtime_error(std::format("[prmtop reader] '{}': POINTERS section missing or too short", file.string()));

  prmtop.atomNames = strings(sections, "ATOM_NAME");
  prmtop.charges = reals(sections, "CHARGE");
  // CHAMBER prmtops (CTITLE) store q * sqrt(332.0716) (the CHARMM Coulomb constant), AMBER ones q * 18.2223
  const double chargeUnit = sections.contains("CTITLE") ? chamberChargeUnit : amberChargeUnit;
  for (double& charge : prmtop.charges) charge /= chargeUnit;
  prmtop.atomicNumbers = integers(sections, "ATOMIC_NUMBER");
  prmtop.masses = reals(sections, "MASS");
  prmtop.atomTypeIndex = integers(sections, "ATOM_TYPE_INDEX");
  prmtop.numberExcludedAtoms = integers(sections, "NUMBER_EXCLUDED_ATOMS");
  prmtop.nonbondedParmIndex = integers(sections, "NONBONDED_PARM_INDEX");
  prmtop.residueLabels = strings(sections, "RESIDUE_LABEL");
  prmtop.residuePointer = integers(sections, "RESIDUE_POINTER");
  prmtop.bondForceConstant = reals(sections, "BOND_FORCE_CONSTANT");
  prmtop.bondEquilValue = reals(sections, "BOND_EQUIL_VALUE");
  prmtop.angleForceConstant = reals(sections, "ANGLE_FORCE_CONSTANT");
  prmtop.angleEquilValue = reals(sections, "ANGLE_EQUIL_VALUE");
  prmtop.dihedralForceConstant = reals(sections, "DIHEDRAL_FORCE_CONSTANT");
  prmtop.dihedralPeriodicity = reals(sections, "DIHEDRAL_PERIODICITY");
  prmtop.dihedralPhase = reals(sections, "DIHEDRAL_PHASE");
  prmtop.sceeScaleFactor = reals(sections, "SCEE_SCALE_FACTOR");
  prmtop.scnbScaleFactor = reals(sections, "SCNB_SCALE_FACTOR");
  prmtop.lennardJonesACoef = reals(sections, "LENNARD_JONES_ACOEF");
  prmtop.lennardJonesBCoef = reals(sections, "LENNARD_JONES_BCOEF");
  prmtop.lennardJones14ACoef = reals(sections, "LENNARD_JONES_14_ACOEF");
  prmtop.lennardJones14BCoef = reals(sections, "LENNARD_JONES_14_BCOEF");
  prmtop.bondsIncHydrogen = integers(sections, "BONDS_INC_HYDROGEN");
  prmtop.bondsWithoutHydrogen = integers(sections, "BONDS_WITHOUT_HYDROGEN");
  prmtop.anglesIncHydrogen = integers(sections, "ANGLES_INC_HYDROGEN");
  prmtop.anglesWithoutHydrogen = integers(sections, "ANGLES_WITHOUT_HYDROGEN");
  prmtop.dihedralsIncHydrogen = integers(sections, "DIHEDRALS_INC_HYDROGEN");
  prmtop.dihedralsWithoutHydrogen = integers(sections, "DIHEDRALS_WITHOUT_HYDROGEN");
  prmtop.excludedAtomsList = integers(sections, "EXCLUDED_ATOMS_LIST");
  prmtop.amberAtomTypes = strings(sections, "AMBER_ATOM_TYPE");
  prmtop.solventPointers = integers(sections, "SOLVENT_POINTERS");
  prmtop.atomsPerMolecule = integers(sections, "ATOMS_PER_MOLECULE");
  prmtop.boxDimensions = reals(sections, "BOX_DIMENSIONS");

  // CMAP (ff19SB: CMAP_*; CHAMBER: CHARMM_CMAP_*)
  {
    const std::string prefix = sections.contains("CHARMM_CMAP_COUNT") ? "CHARMM_CMAP_" : "CMAP_";
    const std::vector<std::int64_t> count = integers(sections, prefix + "COUNT");
    if (count.size() >= 2)
    {
      const std::size_t numberOfTerms = static_cast<std::size_t>(std::max<std::int64_t>(count[0], 0));
      const std::size_t numberOfMaps = static_cast<std::size_t>(std::max<std::int64_t>(count[1], 0));
      prmtop.cmapResolution = integers(sections, prefix + "RESOLUTION");
      if (prmtop.cmapResolution.size() != numberOfMaps)
        throw std::runtime_error(std::format("[prmtop reader] '{}': {}RESOLUTION has {} entries, expected {} maps",
                                             file.string(), prefix, prmtop.cmapResolution.size(), numberOfMaps));
      for (std::size_t map = 0; map < numberOfMaps; ++map)
      {
        const std::string flag = std::format("{}PARAMETER_{:02d}", prefix, map + 1);
        std::vector<double> values = reals(sections, flag);
        const std::size_t resolution = static_cast<std::size_t>(prmtop.cmapResolution[map]);
        if (values.size() != resolution * resolution)
          throw std::runtime_error(std::format("[prmtop reader] '{}': section {} has {} entries, expected {} ({}^2)",
                                               file.string(), flag, values.size(), resolution * resolution,
                                               resolution));
        prmtop.cmapParameters.push_back(std::move(values));
      }
      prmtop.cmapIndex = integers(sections, prefix + "INDEX");
      if (prmtop.cmapIndex.size() != 6 * numberOfTerms)
        throw std::runtime_error(std::format("[prmtop reader] '{}': {}INDEX has {} entries, expected {} (6 x {})",
                                             file.string(), prefix, prmtop.cmapIndex.size(), 6 * numberOfTerms,
                                             numberOfTerms));
    }
  }

  std::size_t numberOfAtoms = prmtop.numberOfAtoms();
  auto require = [&](std::size_t size, std::string_view flag)
  {
    if (size != numberOfAtoms)
      throw std::runtime_error(std::format("[prmtop reader] '{}': section {} has {} entries, expected {} atoms",
                                           file.string(), flag, size, numberOfAtoms));
  };
  require(prmtop.atomNames.size(), "ATOM_NAME");
  require(prmtop.charges.size(), "CHARGE");
  require(prmtop.masses.size(), "MASS");
  require(prmtop.atomTypeIndex.size(), "ATOM_TYPE_INDEX");
  require(prmtop.amberAtomTypes.size(), "AMBER_ATOM_TYPE");
  if (!prmtop.atomicNumbers.empty()) require(prmtop.atomicNumbers.size(), "ATOMIC_NUMBER");

  std::size_t numberOfTypes = prmtop.numberOfTypes();
  if (prmtop.nonbondedParmIndex.size() != numberOfTypes * numberOfTypes)
    throw std::runtime_error(std::format("[prmtop reader] '{}': NONBONDED_PARM_INDEX has {} entries, expected {}",
                                         file.string(), prmtop.nonbondedParmIndex.size(),
                                         numberOfTypes * numberOfTypes));
  if (prmtop.residueLabels.size() != prmtop.residuePointer.size())
    throw std::runtime_error(
        std::format("[prmtop reader] '{}': RESIDUE_LABEL / RESIDUE_POINTER sizes differ", file.string()));
  if (prmtop.bondsIncHydrogen.size() % 3 != 0 || prmtop.bondsWithoutHydrogen.size() % 3 != 0 ||
      prmtop.anglesIncHydrogen.size() % 4 != 0 || prmtop.anglesWithoutHydrogen.size() % 4 != 0 ||
      prmtop.dihedralsIncHydrogen.size() % 5 != 0 || prmtop.dihedralsWithoutHydrogen.size() % 5 != 0)
    throw std::runtime_error(
        std::format("[prmtop reader] '{}': a bond/angle/dihedral list has an unexpected length", file.string()));

  return prmtop;
}

Coordinates readCoordinates(const std::filesystem::path& file, std::size_t numberOfAtoms)
{
  std::ifstream stream(file, std::ios::binary);
  if (!stream) throw std::runtime_error(std::format("[inpcrd reader] cannot open '{}'", file.string()));

  Coordinates coordinates{};
  std::string line{};
  if (!std::getline(stream, line))
    throw std::runtime_error(std::format("[inpcrd reader] '{}' is empty", file.string()));
  if (line.starts_with("CDF") || line.starts_with("\x89HDF"))
    throw std::runtime_error(std::format(
        "[inpcrd reader] '{}' is a NetCDF restart; convert it to ASCII first (cpptraj: trajout file.rst7 restart)",
        file.string()));
  coordinates.title = trim(line);

  if (!std::getline(stream, line))
    throw std::runtime_error(std::format("[inpcrd reader] '{}' has no atom-count line", file.string()));
  std::istringstream header(line);
  std::size_t fileAtoms{};
  if (!(header >> fileAtoms))
    throw std::runtime_error(std::format("[inpcrd reader] '{}': cannot read the number of atoms", file.string()));
  if (numberOfAtoms != 0 && fileAtoms != numberOfAtoms)
    throw std::runtime_error(
        std::format("[inpcrd reader] '{}' has {} atoms, the topology has {}", file.string(), fileAtoms, numberOfAtoms));

  // 6F12.7 fixed-width fields; adjacent negative numbers may touch, so slice rather than split
  std::vector<double> values{};
  while (std::getline(stream, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    for (std::size_t offset = 0; offset < line.size(); offset += 12)
    {
      std::string field = trim(line.substr(offset, 12));
      if (field.empty()) continue;
      try
      {
        values.push_back(std::stod(field));
      }
      catch (std::exception&)
      {
        throw std::runtime_error(std::format("[inpcrd reader] '{}': '{}' is not a number", file.string(), field));
      }
    }
  }

  std::size_t coordinateCount = 3 * fileAtoms;
  if (values.size() < coordinateCount)
    throw std::runtime_error(std::format("[inpcrd reader] '{}': expected {} coordinates, found {}", file.string(),
                                         coordinateCount, values.size()));
  coordinates.positions.resize(fileAtoms);
  for (std::size_t i = 0; i < fileAtoms; ++i)
    coordinates.positions[i] = {values[3 * i], values[3 * i + 1], values[3 * i + 2]};

  std::size_t remaining = values.size() - coordinateCount;
  std::size_t offset = coordinateCount;
  if (remaining == coordinateCount || remaining == coordinateCount + 6)
  {
    coordinates.velocities.resize(fileAtoms);
    for (std::size_t i = 0; i < fileAtoms; ++i)
      coordinates.velocities[i] = {values[offset + 3 * i], values[offset + 3 * i + 1], values[offset + 3 * i + 2]};
    offset += coordinateCount;
    remaining -= coordinateCount;
  }
  if (remaining >= 6)
  {
    coordinates.box = std::array<double, 6>{values[offset],     values[offset + 1], values[offset + 2],
                                            values[offset + 3], values[offset + 4], values[offset + 5]};
    remaining -= 6;
  }
  if (remaining != 0)
    throw std::runtime_error(
        std::format("[inpcrd reader] '{}': {} trailing numbers could not be interpreted", file.string(), remaining));
  return coordinates;
}

ReadResult readPrmtop(const std::filesystem::path& prmtopFile,
                      const std::optional<std::filesystem::path>& coordinateFile, const ReadOptions& options)
{
  ReadResult result{};
  Prmtop prmtop = parsePrmtop(prmtopFile);
  const std::size_t numberOfAtoms = prmtop.numberOfAtoms();
  const std::size_t numberOfTypes = prmtop.numberOfTypes();
  const double energyFactor = kelvinPerKCal();

  std::optional<Coordinates> coordinates{};
  if (coordinateFile) coordinates = readCoordinates(*coordinateFile, numberOfAtoms);

  if (std::ranges::any_of(prmtop.nonbondedParmIndex, [](std::int64_t i) { return i < 0; }))
    result.warnings.push_back("[prmtop reader] 10-12 hydrogen-bond terms are not converted (set to zero)");

  // residue of every atom
  std::vector<std::string> residueOfAtom(numberOfAtoms, "UNK");
  for (std::size_t r = 0; r < prmtop.residuePointer.size(); ++r)
  {
    std::size_t begin = static_cast<std::size_t>(prmtop.residuePointer[r] - 1);
    std::size_t end = r + 1 < prmtop.residuePointer.size() ? static_cast<std::size_t>(prmtop.residuePointer[r + 1] - 1)
                                                           : numberOfAtoms;
    for (std::size_t atom = begin; atom < std::min(end, numberOfAtoms); ++atom)
      residueOfAtom[atom] = prmtop.residueLabels[r];
  }

  // pseudo-atoms
  std::vector<PseudoAtomDefinition> pseudoAtoms{};
  std::map<PseudoAtomKey, std::size_t> pseudoAtomOf{};
  std::vector<std::size_t> pseudoAtomOfAtom(numberOfAtoms);
  for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
  {
    PseudoAtomKey key{prmtop.atomTypeIndex[atom], prmtop.amberAtomTypes[atom], prmtop.masses[atom]};
    auto it = pseudoAtomOf.find(key);
    if (it == pseudoAtomOf.end())
    {
      PseudoAtomDefinition definition{};
      definition.key = key;
      definition.name = sanitizeName(key.typeName);
      std::int64_t atomicNumber = atom < prmtop.atomicNumbers.size() ? prmtop.atomicNumbers[atom] : 0;
      definition.element = elementSymbol(atomicNumber, key.mass);
      pseudoAtoms.push_back(definition);
      it = pseudoAtomOf.emplace(key, pseudoAtoms.size() - 1).first;
    }
    PseudoAtomDefinition& definition = pseudoAtoms[it->second];
    if (!definition.chargeSeen)
    {
      definition.commonCharge = prmtop.charges[atom];
      definition.chargeSeen = true;
    }
    else if (definition.commonCharge && std::abs(*definition.commonCharge - prmtop.charges[atom]) > 1e-8)
    {
      definition.commonCharge.reset();
    }
    pseudoAtomOfAtom[atom] = it->second;
  }
  {
    std::map<std::string, std::size_t> seen{};
    for (PseudoAtomDefinition& definition : pseudoAtoms)
    {
      std::size_t n = seen[definition.name]++;
      if (n > 0) definition.name = std::format("{}_{}", definition.name, n);
    }
  }

  // Lennard-Jones parameters per AMBER type pair
  auto lennardJonesOf = [&](const std::vector<double>& ACoef, const std::vector<double>& BCoef, std::int64_t typeA,
                            std::int64_t typeB) -> std::array<double, 2>
  {
    std::size_t slot = static_cast<std::size_t>(typeA - 1) * numberOfTypes + static_cast<std::size_t>(typeB - 1);
    std::int64_t index = prmtop.nonbondedParmIndex[slot];
    if (index <= 0) return {0.0, 1.0};
    double A = ACoef[static_cast<std::size_t>(index - 1)];
    double B = BCoef[static_cast<std::size_t>(index - 1)];
    if (A <= 0.0 || B <= 0.0) return {0.0, 1.0};
    return {B * B / (4.0 * A), std::pow(A / B, 1.0 / 6.0)};  // epsilon [kcal/mol], sigma [Å]
  };
  auto lennardJones = [&](std::int64_t typeA, std::int64_t typeB)
  { return lennardJonesOf(prmtop.lennardJonesACoef, prmtop.lennardJonesBCoef, typeA, typeB); };
  // CHAMBER prmtops carry the Lennard-Jones coefficients of the 1-4 pairs (CHARMM gives a type separate 1-4
  // parameters); the pairs then use the 1-4 table of the force field with the scaling SCNB = 1
  const bool has14Table = prmtop.lennardJones14ACoef.size() == prmtop.lennardJonesACoef.size() &&
                          prmtop.lennardJones14BCoef.size() == prmtop.lennardJonesBCoef.size() &&
                          !prmtop.lennardJones14ACoef.empty();
  auto lennardJones14 = [&](std::int64_t typeA, std::int64_t typeB)
  { return lennardJonesOf(prmtop.lennardJones14ACoef, prmtop.lennardJones14BCoef, typeA, typeB); };

  // box
  bool periodic = false;
  std::array<double, 3> boxLengths{};
  std::array<double, 3> boxAngles{90.0, 90.0, 90.0};
  if (coordinates && coordinates->box)
  {
    periodic = true;
    boxLengths = {(*coordinates->box)[0], (*coordinates->box)[1], (*coordinates->box)[2]};
    boxAngles = {(*coordinates->box)[3], (*coordinates->box)[4], (*coordinates->box)[5]};
  }
  else if (prmtop.hasBox() && prmtop.boxDimensions.size() >= 4)
  {
    periodic = true;
    boxLengths = {prmtop.boxDimensions[1], prmtop.boxDimensions[2], prmtop.boxDimensions[3]};
    boxAngles = {prmtop.boxDimensions[0], prmtop.boxDimensions[0], prmtop.boxDimensions[0]};
  }

  // bonded lists as (atoms..., type)
  std::vector<std::array<std::size_t, 2>> allBonds{};
  std::vector<std::array<std::int64_t, 3>> bondList{};
  for (const std::vector<std::int64_t>* list : {&prmtop.bondsIncHydrogen, &prmtop.bondsWithoutHydrogen})
    for (std::size_t i = 0; i + 2 < list->size(); i += 3)
    {
      std::int64_t a = (*list)[i] / 3, b = (*list)[i + 1] / 3;
      bondList.push_back({a, b, (*list)[i + 2]});
      allBonds.push_back({static_cast<std::size_t>(a), static_cast<std::size_t>(b)});
    }
  std::vector<std::array<std::int64_t, 4>> angleList{};
  for (const std::vector<std::int64_t>* list : {&prmtop.anglesIncHydrogen, &prmtop.anglesWithoutHydrogen})
    for (std::size_t i = 0; i + 3 < list->size(); i += 4)
      angleList.push_back({(*list)[i] / 3, (*list)[i + 1] / 3, (*list)[i + 2] / 3, (*list)[i + 3]});
  std::vector<std::array<std::int64_t, 6>> dihedralList{};
  for (const std::vector<std::int64_t>* list : {&prmtop.dihedralsIncHydrogen, &prmtop.dihedralsWithoutHydrogen})
    for (std::size_t i = 0; i + 4 < list->size(); i += 5)
    {
      std::int64_t flags = ((*list)[i + 2] < 0 ? 1 : 0) | ((*list)[i + 3] < 0 ? 2 : 0);
      dihedralList.push_back({std::abs((*list)[i]) / 3, std::abs((*list)[i + 1]) / 3, std::abs((*list)[i + 2]) / 3,
                              std::abs((*list)[i + 3]) / 3, (*list)[i + 4], flags});
    }

  // 1-4 scaling (uniform in RASPA; AMBER allows it per dihedral type)
  double scee = 1.2, scnb = 2.0;
  {
    std::set<double> sceeValues{}, scnbValues{};
    for (const std::array<std::int64_t, 6>& dihedral : dihedralList)
    {
      if (dihedral[5] != 0) continue;  // only terms that carry the 1-4 interaction
      std::size_t type = static_cast<std::size_t>(dihedral[4] - 1);
      if (type < prmtop.sceeScaleFactor.size() && prmtop.sceeScaleFactor[type] > 0.0)
        sceeValues.insert(prmtop.sceeScaleFactor[type]);
      if (type < prmtop.scnbScaleFactor.size() && prmtop.scnbScaleFactor[type] > 0.0)
        scnbValues.insert(prmtop.scnbScaleFactor[type]);
    }
    if (!sceeValues.empty()) scee = *sceeValues.begin();
    if (!scnbValues.empty()) scnb = *scnbValues.begin();
    if (sceeValues.size() > 1 || scnbValues.size() > 1)
      result.warnings.push_back(std::format(
          "[prmtop reader] non-uniform 1-4 scale factors (SCEE {} values, SCNB {} values); RASPA applies one pair of "
          "factors per component, using 1/{} and 1/{}",
          sceeValues.size(), scnbValues.size(), scee, scnb));
  }

  // molecules
  std::vector<std::vector<std::size_t>> molecules = findMolecules(numberOfAtoms, allBonds);
  for (const std::vector<std::size_t>& molecule : molecules)
    if (molecule.back() - molecule.front() + 1 != molecule.size())
      throw std::runtime_error(std::format(
          "[prmtop reader] molecule starting at atom {} is not a contiguous atom range; RASPA components require "
          "contiguous molecules",
          molecule.front() + 1));
  if (!prmtop.atomsPerMolecule.empty())
  {
    std::size_t total = 0;
    for (std::int64_t n : prmtop.atomsPerMolecule) total += static_cast<std::size_t>(n);
    if (prmtop.atomsPerMolecule.size() != molecules.size() && total == numberOfAtoms)
      result.warnings.push_back(std::format(
          "[prmtop reader] ATOMS_PER_MOLECULE lists {} molecules, the bond graph has {}; using the bond graph",
          prmtop.atomsPerMolecule.size(), molecules.size()));
  }

  // per-molecule topology, keyed on the first atom of the molecule
  std::vector<std::size_t> moleculeOfAtom(numberOfAtoms);
  for (std::size_t m = 0; m < molecules.size(); ++m)
    for (std::size_t atom : molecules[m]) moleculeOfAtom[atom] = m;

  std::vector<MoleculeTemplate> moleculeTemplates(molecules.size());
  for (std::size_t m = 0; m < molecules.size(); ++m)
  {
    MoleculeTemplate& tpl = moleculeTemplates[m];
    for (std::size_t atom : molecules[m])
    {
      tpl.pseudoAtomIndex.push_back(pseudoAtomOfAtom[atom]);
      tpl.charges.push_back(prmtop.charges[atom]);
      tpl.residues.push_back(residueOfAtom[atom]);
    }
  }
  auto local = [&](std::int64_t atom, std::size_t molecule) -> std::int64_t
  { return atom - static_cast<std::int64_t>(molecules[molecule].front()); };
  for (const std::array<std::int64_t, 3>& bond : bondList)
  {
    std::size_t m = moleculeOfAtom[static_cast<std::size_t>(bond[0])];
    moleculeTemplates[m].bonds.push_back({local(bond[0], m), local(bond[1], m), bond[2]});
  }
  for (const std::array<std::int64_t, 4>& angle : angleList)
  {
    std::size_t m = moleculeOfAtom[static_cast<std::size_t>(angle[0])];
    moleculeTemplates[m].angles.push_back({local(angle[0], m), local(angle[1], m), local(angle[2], m), angle[3]});
  }
  for (const std::array<std::int64_t, 6>& dihedral : dihedralList)
  {
    std::size_t m = moleculeOfAtom[static_cast<std::size_t>(dihedral[0])];
    moleculeTemplates[m].dihedrals.push_back({local(dihedral[0], m), local(dihedral[1], m), local(dihedral[2], m),
                                              local(dihedral[3], m), dihedral[4], dihedral[5]});
  }
  for (std::size_t i = 0; i + 5 < prmtop.cmapIndex.size(); i += 6)
  {
    std::array<std::int64_t, 6> cmap{};
    for (std::size_t k = 0; k < 6; ++k) cmap[k] = prmtop.cmapIndex[i + k] - 1;  // 0-based atoms and map
    for (std::size_t k = 0; k < 5; ++k)
      if (cmap[k] < 0 || static_cast<std::size_t>(cmap[k]) >= numberOfAtoms)
        throw std::runtime_error(std::format("[prmtop reader] CMAP_INDEX refers to atom {} (out of range)", cmap[k] + 1));
    if (cmap[5] < 0 || static_cast<std::size_t>(cmap[5]) >= prmtop.cmapParameters.size())
      throw std::runtime_error(std::format("[prmtop reader] CMAP_INDEX refers to map {} (out of range)", cmap[5] + 1));
    std::size_t m = moleculeOfAtom[static_cast<std::size_t>(cmap[0])];
    for (std::size_t k = 1; k < 5; ++k)
      if (moleculeOfAtom[static_cast<std::size_t>(cmap[k])] != m)
        throw std::runtime_error("[prmtop reader] a CMAP term spans two molecules");
    moleculeTemplates[m].cmaps.push_back(
        {local(cmap[0], m), local(cmap[1], m), local(cmap[2], m), local(cmap[3], m), local(cmap[4], m), cmap[5]});
  }
  for (MoleculeTemplate& tpl : moleculeTemplates)
  {
    std::ranges::sort(tpl.bonds);
    std::ranges::sort(tpl.angles);
    std::ranges::sort(tpl.dihedrals);
    std::ranges::sort(tpl.cmaps);
  }

  // distinct templates -> components (in order of first appearance)
  std::vector<std::size_t> templateIndexOfMolecule(molecules.size());
  std::vector<std::size_t> representative{};  // molecule index of the first instance of every component
  std::vector<std::vector<std::size_t>> instances{};
  for (std::size_t m = 0; m < molecules.size(); ++m)
  {
    std::optional<std::size_t> found{};
    for (std::size_t c = 0; c < representative.size(); ++c)
    {
      if (moleculeTemplates[representative[c]] == moleculeTemplates[m])
      {
        found = c;
        break;
      }
    }
    if (!found)
    {
      representative.push_back(m);
      instances.push_back({});
      found = representative.size() - 1;
    }
    templateIndexOfMolecule[m] = *found;
    instances[*found].push_back(m);
  }

  // positions (shifted to the centre of a vacuum box when the system is not periodic)
  std::vector<std::array<double, 3>> positions{};
  if (coordinates) positions = coordinates->positions;
  double cutOff = options.cutOff;
  if (!periodic)
  {
    double extent = 0.0;
    std::array<double, 3> centre{};
    if (!positions.empty())
    {
      for (const std::array<double, 3>& p : positions)
        for (std::size_t k = 0; k < 3; ++k) centre[k] += p[k] / static_cast<double>(positions.size());
      for (const std::array<double, 3>& p : positions)
      {
        double rr = 0.0;
        for (std::size_t k = 0; k < 3; ++k) rr += (p[k] - centre[k]) * (p[k] - centre[k]);
        extent = std::max(extent, std::sqrt(rr));
      }
    }
    cutOff = std::max(options.cutOff, 2.0 * extent + 2.0);
    double length = 2.0 * cutOff + 2.0;
    boxLengths = {length, length, length};
    for (std::array<double, 3>& p : positions)
      for (std::size_t k = 0; k < 3; ++k) p[k] = p[k] - centre[k] + 0.5 * length;
    result.warnings.push_back(std::format(
        "[prmtop reader] no box: built a {:.1f} Å vacuum box with direct Coulomb and a {:.1f} Å cut-off spanning "
        "every atom pair",
        length, cutOff));
  }
  result.boxLengths = boxLengths;
  result.boxAngles = boxAngles;
  result.positions = positions;

  // force field
  nlohmann::json forceField{};
  forceField["MixingRule"] = "Lorentz-Berthelot";
  forceField["TruncationMethod"] = options.truncationMethod;
  if (options.switchingDistance > 0.0) forceField["SwitchingDistance"] = options.switchingDistance;
  forceField["TailCorrections"] = false;
  forceField["CutOffVDW"] = cutOff;
  forceField["CutOffCoulomb"] = cutOff;
  if (periodic) forceField["EwaldPrecision"] = 1e-6;
  forceField["PseudoAtoms"] = nlohmann::json::array();
  for (const PseudoAtomDefinition& definition : pseudoAtoms)
  {
    forceField["PseudoAtoms"].push_back({{"name", definition.name},
                                         {"framework", false},
                                         {"print_to_output", true},
                                         {"element", definition.element},
                                         {"print_as", definition.element},
                                         {"mass", definition.key.mass},
                                         {"charge", definition.commonCharge.value_or(0.0)}});
  }
  forceField["SelfInteractions"] = nlohmann::json::array();
  std::vector<std::array<double, 2>> selfParameters(pseudoAtoms.size());
  std::vector<std::array<double, 2>> selfParameters14(pseudoAtoms.size());
  for (std::size_t i = 0; i < pseudoAtoms.size(); ++i)
  {
    selfParameters[i] = lennardJones(pseudoAtoms[i].key.typeIndex, pseudoAtoms[i].key.typeIndex);
    nlohmann::json self = {{"name", pseudoAtoms[i].name},
                           {"type", "lennard-jones"},
                           {"parameters", {selfParameters[i][0] * energyFactor, selfParameters[i][1]}}};
    if (has14Table)
    {
      selfParameters14[i] = lennardJones14(pseudoAtoms[i].key.typeIndex, pseudoAtoms[i].key.typeIndex);
      self["parameters14"] = {selfParameters14[i][0] * energyFactor, selfParameters14[i][1]};
    }
    forceField["SelfInteractions"].push_back(self);
  }
  auto mixes = [](const std::array<double, 2>& pair, const std::array<double, 2>& selfA,
                  const std::array<double, 2>& selfB)
  {
    double mixedEpsilon = std::sqrt(selfA[0] * selfB[0]);
    double mixedSigma = 0.5 * (selfA[1] + selfB[1]);
    bool epsilonMatches = std::abs(pair[0] - mixedEpsilon) <= 1e-6 * std::max(1.0, std::abs(pair[0]));
    bool sigmaMatches = pair[0] == 0.0 || std::abs(pair[1] - mixedSigma) <= 1e-6 * std::max(1.0, pair[1]);
    return epsilonMatches && sigmaMatches;
  };
  nlohmann::json binary = nlohmann::json::array();
  for (std::size_t i = 0; i < pseudoAtoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j < pseudoAtoms.size(); ++j)
    {
      std::array<double, 2> pair = lennardJones(pseudoAtoms[i].key.typeIndex, pseudoAtoms[j].key.typeIndex);
      const bool pairMixes = mixes(pair, selfParameters[i], selfParameters[j]);
      std::array<double, 2> pair14{};
      bool pair14Mixes = true;
      if (has14Table)
      {
        pair14 = lennardJones14(pseudoAtoms[i].key.typeIndex, pseudoAtoms[j].key.typeIndex);
        pair14Mixes = mixes(pair14, selfParameters14[i], selfParameters14[j]);
      }
      if (pairMixes && pair14Mixes) continue;
      nlohmann::json entry = {{"names", {pseudoAtoms[i].name, pseudoAtoms[j].name}}, {"type", "lennard-jones"}};
      if (!pairMixes) entry["parameters"] = {pair[0] * energyFactor, pair[1]};
      if (!pair14Mixes) entry["parameters14"] = {pair14[0] * energyFactor, pair14[1]};
      binary.push_back(entry);
    }
  }
  if (!binary.empty()) forceField["BinaryInteractions"] = binary;

  // CMAP maps: the prmtop grid order (first angle slow, angles from -180 degrees) is RASPA's; kcal/mol -> K
  const auto cmapName = [](std::size_t map) { return std::format("CMAP_{}", map + 1); };
  if (!prmtop.cmapParameters.empty())
  {
    forceField["CMAPs"] = nlohmann::json::array();
    for (std::size_t map = 0; map < prmtop.cmapParameters.size(); ++map)
    {
      std::vector<double> energies = prmtop.cmapParameters[map];
      for (double& energy : energies) energy *= energyFactor;
      forceField["CMAPs"].push_back(
          {{"Name", cmapName(map)}, {"Resolution", prmtop.cmapResolution[map]}, {"Energies", energies}});
    }
  }
  result.forceField = forceField;

  // components
  std::map<std::string, std::size_t> namesSeen{};
  std::size_t unnamedCounter = 0;
  for (std::size_t c = 0; c < representative.size(); ++c)
  {
    const std::size_t m = representative[c];
    const MoleculeTemplate& tpl = moleculeTemplates[m];
    const std::vector<std::size_t>& atoms = molecules[m];
    const std::size_t n = atoms.size();

    ReadComponent component{};
    component.count = instances[c].size();
    component.atomsPerMolecule = n;
    for (std::size_t instance : instances[c]) component.firstAtoms.push_back(molecules[instance].front());

    // name
    std::set<std::string> residues(tpl.residues.begin(), tpl.residues.end());
    if (residues.size() == 1)
      component.name = sanitizeName(*residues.begin());
    else if (std::ranges::all_of(residues,
                                 [](const std::string& r) { return aminoAcids.contains(baseResidue(toUpper(r))); }))
      component.name = "protein";
    else if (std::ranges::all_of(residues, [](const std::string& r) { return nucleicAcids.contains(toUpper(r)); }))
      component.name = "nucleic_acid";
    else
      component.name = std::format("molecule{}", ++unnamedCounter);
    if (std::size_t k = namesSeen[component.name]++; k > 0) component.name = std::format("{}_{}", component.name, k);

    // rigid: single atoms and the residues asked for (water)
    bool rigid =
        n == 1 || (residues.size() == 1 && std::ranges::any_of(options.rigidResidues, [&](const std::string& r)
                                                               { return toUpper(r) == toUpper(*residues.begin()); }));
    component.rigid = rigid;

    // reference geometry: the first instance, centred at its centre of mass
    std::vector<std::array<double, 3>> geometry(n, std::array<double, 3>{});
    std::vector<double> masses(n);
    for (std::size_t k = 0; k < n; ++k) masses[k] = prmtop.masses[atoms[k]];
    if (!positions.empty())
      for (std::size_t k = 0; k < n; ++k) geometry[k] = positions[atoms[k]];
    else if (n > 1 && !rigid)
      result.warnings.push_back(std::format(
          "[prmtop reader] component '{}': no coordinates given, the reference geometry is all zeros", component.name));

    // rigid three-site water: ideal geometry from the equilibrium bond lengths and angle
    bool idealWater = false;
    if (rigid && n == 3 && tpl.bonds.size() >= 2)
    {
      std::array<std::size_t, 3> bondCount{};
      for (const std::array<std::int64_t, 3>& bond : tpl.bonds)
      {
        ++bondCount[static_cast<std::size_t>(bond[0])];
        ++bondCount[static_cast<std::size_t>(bond[1])];
      }
      std::optional<std::size_t> centralAtom{};
      for (std::size_t k = 0; k < 3; ++k)
        if (bondCount[k] == 2 && masses[k] > masses[(k + 1) % 3] && masses[k] > masses[(k + 2) % 3]) centralAtom = k;
      if (!centralAtom)
        for (std::size_t k = 0; k < 3; ++k)
          if (masses[k] >= masses[(k + 1) % 3] && masses[k] >= masses[(k + 2) % 3]) centralAtom = k;
      if (centralAtom)
      {
        std::size_t o = *centralAtom;
        std::size_t h1 = (o + 1) % 3, h2 = (o + 2) % 3;
        if (h1 > h2) std::swap(h1, h2);
        std::optional<double> rOH{}, rHH{}, theta{};
        for (const std::array<std::int64_t, 3>& bond : tpl.bonds)
        {
          std::size_t a = static_cast<std::size_t>(bond[0]), b = static_cast<std::size_t>(bond[1]);
          double r0 = prmtop.bondEquilValue[static_cast<std::size_t>(bond[2] - 1)];
          if (a == o || b == o)
            rOH = r0;
          else
            rHH = r0;
        }
        for (const std::array<std::int64_t, 4>& angle : tpl.angles)
          if (static_cast<std::size_t>(angle[1]) == o)
            theta = prmtop.angleEquilValue[static_cast<std::size_t>(angle[3] - 1)];
        if (rOH && !theta && rHH) theta = 2.0 * std::asin(std::clamp(*rHH / (2.0 * *rOH), 0.0, 1.0));
        if (rOH && theta)
        {
          geometry[o] = {0.0, 0.0, 0.0};
          geometry[h1] = {*rOH * std::sin(0.5 * *theta), *rOH * std::cos(0.5 * *theta), 0.0};
          geometry[h2] = {-*rOH * std::sin(0.5 * *theta), *rOH * std::cos(0.5 * *theta), 0.0};
          idealWater = true;
        }
      }
    }
    std::array<double, 3> com = centreOfMass(geometry, masses);
    for (std::array<double, 3>& p : geometry)
      for (std::size_t k = 0; k < 3; ++k) p[k] -= com[k];

    nlohmann::json definition{};
    definition["CriticalTemperature"] = 500.0;
    definition["CriticalPressure"] = 5.0e6;
    definition["AcentricFactor"] = 0.3;
    definition["_comment"] =
        std::format("converted from AMBER prmtop '{}'{}; critical constants are placeholders{}",
                    prmtopFile.filename().string(), prmtop.title.empty() ? "" : std::format(" ({})", prmtop.title),
                    rigid ? (idealWater ? "; rigid body with the equilibrium water geometry" : "; rigid body")
                          : "; geometry is the first instance");
    definition["PseudoAtoms"] = nlohmann::json::array();
    for (std::size_t k = 0; k < n; ++k)
    {
      const PseudoAtomDefinition& pseudoAtom = pseudoAtoms[tpl.pseudoAtomIndex[k]];
      nlohmann::json entry = {pseudoAtom.name, {geometry[k][0], geometry[k][1], geometry[k][2]}};
      if (!pseudoAtom.commonCharge || std::abs(*pseudoAtom.commonCharge - tpl.charges[k]) > 1e-8)
        entry.push_back(tpl.charges[k]);
      definition["PseudoAtoms"].push_back(entry);
    }

    if (!rigid && !tpl.bonds.empty())
    {
      // adjacency for the bend/torsion enumeration RASPA performs from 'Connectivity'
      std::vector<std::set<std::size_t>> neighbours(n);
      definition["Connectivity"] = nlohmann::json::array();
      definition["Bonds"] = nlohmann::json::array();
      for (const std::array<std::int64_t, 3>& bond : tpl.bonds)
      {
        std::size_t a = static_cast<std::size_t>(bond[0]), b = static_cast<std::size_t>(bond[1]);
        neighbours[a].insert(b);
        neighbours[b].insert(a);
        definition["Connectivity"].push_back({a, b});
        std::size_t type = static_cast<std::size_t>(bond[2] - 1);
        // AMBER: k (r - r0)^2 ; RASPA HARMONIC: (1/2) p0 (r - p1)^2
        definition["Bonds"].push_back(
            {{a, b}, "HARMONIC", {2.0 * prmtop.bondForceConstant[type] * energyFactor, prmtop.bondEquilValue[type]}});
      }

      // bends: every i-j-k of the bond graph needs a definition
      std::map<std::array<std::size_t, 3>, std::size_t> angleType{};
      for (const std::array<std::int64_t, 4>& angle : tpl.angles)
      {
        std::array<std::size_t, 3> ids{static_cast<std::size_t>(angle[0]), static_cast<std::size_t>(angle[1]),
                                       static_cast<std::size_t>(angle[2])};
        if (ids[0] > ids[2]) std::swap(ids[0], ids[2]);
        angleType[ids] = static_cast<std::size_t>(angle[3] - 1);
      }
      definition["Bends"] = nlohmann::json::array();
      std::size_t missingBends = 0;
      for (std::size_t j = 0; j < n; ++j)
      {
        for (std::size_t i : neighbours[j])
          for (std::size_t k : neighbours[j])
          {
            if (i >= k) continue;
            auto it = angleType.find({i, j, k});
            if (it == angleType.end())
            {
              ++missingBends;
              definition["Bends"].push_back({{i, j, k}, "HARMONIC", {0.0, 109.5}});
              continue;
            }
            // AMBER: k (theta - theta0)^2 [theta0 in rad]; RASPA HARMONIC: (1/2) p0 (theta - p1)^2 [p1 in degrees]
            definition["Bends"].push_back({{i, j, k},
                                           "HARMONIC",
                                           {2.0 * prmtop.angleForceConstant[it->second] * energyFactor,
                                            prmtop.angleEquilValue[it->second] * Units::RadiansToDegrees}});
          }
      }
      if (missingBends > 0)
        result.warnings.push_back(std::format(
            "[prmtop reader] component '{}': {} bend(s) of the bond graph have no AMBER angle term (set to zero)",
            component.name, missingBends));

      // torsions: AMBER terms grouped per quadruple
      auto canonical = [](std::array<std::size_t, 4> ids)
      {
        if (ids[1] > ids[2] || (ids[1] == ids[2] && ids[0] > ids[3])) std::ranges::reverse(ids);
        return ids;
      };
      std::map<std::array<std::size_t, 4>, std::vector<DihedralTerm>> properTerms{};
      std::vector<std::pair<std::array<std::size_t, 4>, DihedralTerm>> improperTerms{};
      for (const std::array<std::int64_t, 6>& dihedral : tpl.dihedrals)
      {
        std::size_t type = static_cast<std::size_t>(dihedral[4] - 1);
        DihedralTerm term{prmtop.dihedralForceConstant[type], prmtop.dihedralPeriodicity[type],
                          prmtop.dihedralPhase[type]};
        std::array<std::size_t, 4> ids{static_cast<std::size_t>(dihedral[0]), static_cast<std::size_t>(dihedral[1]),
                                       static_cast<std::size_t>(dihedral[2]), static_cast<std::size_t>(dihedral[3])};
        if (dihedral[5] & 2)
        {
          if (term.forceConstant != 0.0) improperTerms.emplace_back(ids, term);
          continue;
        }
        properTerms[canonical(ids)].push_back(term);
      }

      definition["Torsions"] = nlohmann::json::array();
      nlohmann::json impropers = nlohmann::json::array();
      std::set<std::array<std::size_t, 4>> graphTorsions{};
      for (std::size_t b = 0; b < n; ++b)
        for (std::size_t c : neighbours[b])
        {
          if (b >= c) continue;
          for (std::size_t a : neighbours[b])
          {
            if (a == c) continue;
            for (std::size_t d : neighbours[c])
            {
              if (d == b) continue;
              graphTorsions.insert(canonical({a, b, c, d}));
            }
          }
        }
      std::size_t missingTorsions = 0, offPhaseTerms = 0;
      for (const std::array<std::size_t, 4>& ids : graphTorsions)
      {
        std::vector<DihedralTerm> polynomialTerms{};
        auto it = properTerms.find(ids);
        if (it != properTerms.end())
        {
          for (const DihedralTerm& term : it->second)
          {
            if (term.forceConstant == 0.0) continue;
            if (isPolynomialTerm(term))
              polynomialTerms.push_back(term);
            else
            {
              ++offPhaseTerms;
              impropers.push_back(
                  {{ids[0], ids[1], ids[2], ids[3]},
                   "CVFF",
                   {term.forceConstant * energyFactor, term.periodicity, term.phase * Units::RadiansToDegrees}});
            }
          }
        }
        else
        {
          ++missingTorsions;
        }
        std::array<double, 6> coefficients = polynomialCoefficients(polynomialTerms, energyFactor);
        definition["Torsions"].push_back({{ids[0], ids[1], ids[2], ids[3]}, "POLYNOMIAL", coefficients});
      }
      for (const auto& [ids, term] : improperTerms)
      {
        // AMBER impropers are ordinary dihedrals over the four listed atoms (third atom central)
        impropers.push_back(
            {{ids[0], ids[1], ids[2], ids[3]},
             "CVFF",
             {term.forceConstant * energyFactor, term.periodicity, term.phase * Units::RadiansToDegrees}});
      }
      if (!impropers.empty()) definition["ImproperTorsions"] = impropers;
      if (missingTorsions > 0)
        result.warnings.push_back(std::format(
            "[prmtop reader] component '{}': {} torsion(s) of the bond graph have no AMBER dihedral term (set to zero)",
            component.name, missingTorsions));
      if (offPhaseTerms > 0)
        result.warnings.push_back(std::format(
            "[prmtop reader] component '{}': {} dihedral term(s) with a phase other than 0/180 degrees or a "
            "periodicity above 5 are written as CVFF entries under 'ImproperTorsions'",
            component.name, offPhaseTerms));
      std::size_t properQuadruplesOutsideGraph = 0;
      for (const auto& [ids, _] : properTerms)
        if (!graphTorsions.contains(ids)) ++properQuadruplesOutsideGraph;
      if (properQuadruplesOutsideGraph > 0)
        result.warnings.push_back(std::format(
            "[prmtop reader] component '{}': {} AMBER dihedral(s) do not follow a bonded path and were dropped",
            component.name, properQuadruplesOutsideGraph));

      if (!tpl.cmaps.empty())
      {
        definition["CMAPTorsions"] = nlohmann::json::array();
        for (const std::array<std::int64_t, 6>& cmap : tpl.cmaps)
        {
          definition["CMAPTorsions"].push_back(
              {{cmap[0], cmap[1], cmap[2], cmap[3], cmap[4]}, cmapName(static_cast<std::size_t>(cmap[5]))});
        }
      }

      definition["Intra14VanDerWaalsScalingValue"] = 1.0 / scnb;
      definition["Intra14ChargeChargeScalingValue"] = 1.0 / scee;
    }
    else if (!rigid && tpl.bonds.empty() && n > 1)
    {
      result.warnings.push_back(std::format(
          "[prmtop reader] component '{}' has {} atoms but no bonds; modelled as a rigid body", component.name, n));
      component.rigid = true;
    }

    component.definition = definition;
    result.components.push_back(component);
  }

  // restart (coordinates)
  if (!positions.empty())
  {
    nlohmann::json restart{};
    restart["SimulationBox"] = {{"length-a", boxLengths[0]},  {"length-b", boxLengths[1]},
                                {"length-c", boxLengths[2]},  {"angle-alpha", boxAngles[0]},
                                {"angle-beta", boxAngles[1]}, {"angle-gamma", boxAngles[2]}};
    for (const ReadComponent& component : result.components)
    {
      nlohmann::json list = nlohmann::json::array();
      for (std::size_t first : component.firstAtoms)
        for (std::size_t k = 0; k < component.atomsPerMolecule; ++k)
        {
          const std::array<double, 3>& p = positions[first + k];
          list.push_back({p[0], p[1], p[2]});
        }
      restart[component.name] = list;
    }
    result.restart = restart;
  }

  // simulation skeleton
  nlohmann::json simulation{};
  simulation["SimulationType"] = "MolecularDynamics";
  simulation["NumberOfInitializationCycles"] = 0;
  simulation["NumberOfEquilibrationCycles"] = 0;
  simulation["NumberOfProductionCycles"] = 1000;
  simulation["PrintEvery"] = 100;
  simulation["ForceField"] = ".";
  nlohmann::json system{};
  system["Type"] = "Box";
  system["BoxLengths"] = boxLengths;
  if (boxAngles != std::array<double, 3>{90.0, 90.0, 90.0}) system["BoxAngles"] = boxAngles;
  system["ExternalTemperature"] = options.temperature;
  system["Ensemble"] = "NVT";
  system["TimeStep"] = 0.001;
  system["ChargeMethod"] = periodic ? "Ewald" : "Coulomb";
  if (!result.restart.empty()) system["RestartFileName"] = "restart.json";
  simulation["Systems"] = nlohmann::json::array({system});
  simulation["Components"] = nlohmann::json::array();
  for (const ReadComponent& component : result.components)
  {
    nlohmann::json entry{};
    entry["Name"] = component.name;
    entry["CreateNumberOfMolecules"] = result.restart.empty() ? component.count : 0;
    simulation["Components"].push_back(entry);
  }
  result.simulation = simulation;

  return result;
}

void writeRaspaInput(const ReadResult& result, const std::filesystem::path& directory)
{
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory / "force_field.json");
    file << result.forceField.dump(2) << '\n';
  }
  for (const ReadComponent& component : result.components)
  {
    std::ofstream file(directory / (component.name + ".json"));
    file << component.definition.dump(2) << '\n';
  }
  {
    std::ofstream file(directory / "simulation.json");
    file << result.simulation.dump(2) << '\n';
  }
  if (!result.restart.empty())
  {
    std::ofstream file(directory / "restart.json");
    file << result.restart.dump(1) << '\n';
  }
  if (!result.warnings.empty())
  {
    std::ofstream file(directory / "conversion_warnings.txt");
    for (const std::string& warning : result.warnings) file << warning << '\n';
  }
}
}  // namespace AMBER
