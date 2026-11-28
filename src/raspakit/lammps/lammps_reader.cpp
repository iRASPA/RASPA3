module;

module lammps_reader;

import std;

import json;
import lammps_styles;

namespace LAMMPS
{
namespace
{
// ---------------------------------------------------------------------------------------------------
// Small text helpers
// ---------------------------------------------------------------------------------------------------
std::string stripComment(std::string_view line)
{
  std::size_t hash = line.find('#');
  std::string s(hash == std::string_view::npos ? line : line.substr(0, hash));
  std::size_t first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  std::size_t last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

std::string trailingComment(std::string_view line)
{
  std::size_t hash = line.find('#');
  if (hash == std::string_view::npos) return {};
  std::string s(line.substr(hash + 1));
  std::size_t first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  std::size_t last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

std::vector<std::string> tokens(std::string_view line)
{
  std::vector<std::string> result;
  std::istringstream stream{std::string(line)};
  std::string word;
  while (stream >> word) result.push_back(word);
  return result;
}

double toDouble(const std::string &s) { return std::stod(s); }

bool isNumber(const std::string &s)
{
  if (s.empty()) return false;
  char *end{};
  std::strtod(s.c_str(), &end);
  return end != s.c_str() && *end == '\0';
}

std::string elementForMass(double mass)
{
  struct Entry
  {
    double mass;
    const char *symbol;
  };
  static constexpr std::array<Entry, 30> table{{{1.008, "H"},   {4.0026, "He"}, {6.94, "Li"},   {9.0122, "Be"},
                                                {10.81, "B"},   {12.011, "C"},  {14.007, "N"},  {15.999, "O"},
                                                {18.998, "F"},  {20.180, "Ne"}, {22.990, "Na"}, {24.305, "Mg"},
                                                {26.982, "Al"}, {28.085, "Si"}, {30.974, "P"},  {32.06, "S"},
                                                {35.45, "Cl"},  {39.948, "Ar"}, {39.098, "K"},  {40.078, "Ca"},
                                                {47.867, "Ti"}, {55.845, "Fe"}, {58.693, "Ni"}, {63.546, "Cu"},
                                                {65.38, "Zn"},  {79.904, "Br"}, {83.798, "Kr"}, {126.90, "I"},
                                                {131.29, "Xe"}, {91.224, "Zr"}}};
  // united-atom carbons (CH, CH2, CH3, CH4): 13 - 16.05 -> carbon
  if (mass > 12.5 && mass < 16.1 && std::abs(mass - 14.007) > 0.005 && std::abs(mass - 15.999) > 0.005) return "C";
  const Entry *best = nullptr;
  double bestDiff = std::numeric_limits<double>::max();
  for (const Entry &entry : table)
  {
    double diff = std::abs(entry.mass - mass);
    if (diff < bestDiff)
    {
      bestDiff = diff;
      best = &entry;
    }
  }
  if (best && bestDiff < 0.5) return best->symbol;
  return "C";
}

// ---------------------------------------------------------------------------------------------------
// Parsed representation of the LAMMPS files
// ---------------------------------------------------------------------------------------------------
struct StyleSpec
{
  std::string name{"none"};
  std::vector<std::string> arguments{};
  std::vector<std::string> hybridStyles{};  ///< sub-styles for hybrid / hybrid/overlay

  bool hybrid() const { return name == "hybrid" || name == "hybrid/overlay"; }
};

struct CoefficientLine
{
  std::string style{};  ///< empty when not given (non-hybrid)
  std::vector<double> values{};
};

struct DataAtom
{
  std::size_t id{};
  std::size_t molecule{};
  std::size_t type{};
  double charge{};
  std::array<double, 3> position{};
  std::array<int, 3> image{0, 0, 0};
};

struct DataBonded
{
  std::size_t type{};
  std::vector<std::size_t> atoms{};
};

struct DataFile
{
  std::size_t numberOfAtomTypes{}, numberOfBondTypes{}, numberOfAngleTypes{}, numberOfDihedralTypes{},
      numberOfImproperTypes{};
  std::array<double, 3> low{}, high{};
  double xy{}, xz{}, yz{};
  bool triclinic{false};
  std::string atomStyle{"full"};
  std::map<std::size_t, std::array<double, 3>> velocities{};  ///< 'Velocities' section: by atom id, units of the file
  std::map<std::size_t, double> masses{};
  std::map<std::size_t, std::string> typeNames{};
  std::map<std::size_t, CoefficientLine> pairCoefficients{};
  std::map<std::pair<std::size_t, std::size_t>, CoefficientLine> pairIJCoefficients{};
  std::map<std::size_t, CoefficientLine> bondCoefficients{}, angleCoefficients{}, dihedralCoefficients{},
      improperCoefficients{};
  std::vector<DataAtom> atoms{};
  std::vector<DataBonded> bonds{}, angles{}, dihedrals{}, impropers{};
};

struct InputScript
{
  std::string units{"real"};
  std::optional<std::string> atomStyle{};  ///< 'atom_style': the columns of the Atoms section without a style comment
  StyleSpec bond{}, angle{}, dihedral{}, improper{}, pair{};
  bool haveBond{false}, haveAngle{false}, haveDihedral{false}, haveImproper{false}, havePair{false};
  std::string mixing{};  ///< geometric (LAMMPS default), arithmetic, sixthpower
  bool shift{false};
  bool tail{false};
  std::optional<std::array<double, 3>> specialLJ{};
  std::optional<std::array<double, 3>> specialCoulomb{};
  std::optional<std::string> kspace{};
  std::optional<double> kspaceAccuracy{};
  std::map<std::pair<std::size_t, std::size_t>, CoefficientLine> pairCoefficients{};
  std::map<std::size_t, CoefficientLine> bondCoefficients{}, angleCoefficients{}, dihedralCoefficients{},
      improperCoefficients{};

  /// 'fix shake' (or 'fix rattle'): the bonds and angles LAMMPS holds as constraints. They become FIXED terms,
  /// which the molecular dynamics of RASPA holds with SHAKE/RATTLE and its Monte Carlo through its moves.
  struct Shake
  {
    double tolerance{};
    std::size_t iterations{};
    std::set<std::size_t> bondTypes{}, angleTypes{}, atomTypes{};
    std::vector<double> masses{};
  };
  std::optional<Shake> shake{};

  /// The run the script asks for: 'timestep', 'run' and the fixes that integrate, thermostat and barostat. Times,
  /// temperatures and pressures are in the units of the script, converted when the simulation is written.
  struct Dynamics
  {
    std::optional<double> timeStep{};
    std::vector<std::size_t> runs{};  ///< the steps of every 'run' with at least one step, in order
    std::string integrator{};         ///< "nve", "nvt" or "npt": the fix that moves the atoms; empty without one
    std::optional<double> temperature{}, temperatureDamping{};
    std::optional<double> pressure{}, pressureDamping{};
    bool flexibleCell{false};  ///< 'aniso', 'tri' or separate 'x', 'y', 'z' of fix npt
    std::optional<std::size_t> thermostatChain{}, barostatChain{};
    std::optional<double> langevinTemperature{}, langevinDamping{};
  };
  Dynamics dynamics{};
  /// 'velocity ... create': the script draws new velocities, so those of the data file are not what it starts from
  bool createsVelocities{false};
};

CoefficientLine parseCoefficientLine(const std::vector<std::string> &words, std::size_t firstValue)
{
  CoefficientLine line;
  std::size_t i = firstValue;
  if (i < words.size() && !isNumber(words[i]))
  {
    line.style = words[i];
    ++i;
  }
  for (; i < words.size(); ++i)
  {
    if (isNumber(words[i])) line.values.push_back(toDouble(words[i]));
  }
  return line;
}

StyleSpec parseStyle(const std::vector<std::string> &words)
{
  StyleSpec spec;
  if (words.size() < 2) return spec;
  spec.name = words[1];
  spec.arguments.assign(words.begin() + 2, words.end());
  if (spec.hybrid())
  {
    for (const std::string &word : spec.arguments)
    {
      if (!isNumber(word)) spec.hybridStyles.push_back(word);
    }
  }
  return spec;
}

// The fixes that integrate the equations of motion and hold the temperature and the pressure: 'fix nve', 'fix nvt',
// 'fix npt' and 'fix langevin' (with an accelerator suffix or without). Returns whether 'words' is such a fix (or one
// of their kind that is not converted, which is reported).
bool readDynamicsFix(const std::vector<std::string> &words, InputScript::Dynamics &dynamics,
                     std::vector<std::string> &warnings)
{
  std::string style = words[3];
  for (const std::string_view suffix : {"/omp", "/kk", "/gpu", "/intel", "/opt"})
  {
    if (style.ends_with(suffix)) style.erase(style.size() - suffix.size());
  }
  const auto number = [&](std::size_t i) -> std::optional<double>
  {
    if (i < words.size() && isNumber(words[i])) return toDouble(words[i]);
    return std::nullopt;
  };
  // three numbers after a keyword: the start and stop values (a ramp is not converted) and the damping time
  const auto triple = [&](std::size_t i, std::string_view what, std::optional<double> &value,
                          std::optional<double> &damping)
  {
    const std::optional<double> start = number(i), stop = number(i + 1), damp = number(i + 2);
    if (!start || !stop || !damp)
    {
      warnings.push_back(std::format("[LAMMPS reader] fix {}: the {} is not given by numbers and is not converted",
                                     words[3], what));
      return;
    }
    if (*start != *stop)
    {
      warnings.push_back(std::format("[LAMMPS reader] fix {}: the {} ramp from {} to {} is not converted, {} is used",
                                     words[3], what, *start, *stop, *stop));
    }
    value = *stop;
    damping = *damp;
  };

  if (style != "nve" && style != "nvt" && style != "npt" && style != "langevin")
  {
    const bool ofTheirKind = style.starts_with("nve/") || style.starts_with("nvt/") || style.starts_with("npt/") ||
                             style.starts_with("nph") || style.starts_with("temp/") ||
                             style.starts_with("press/") || style == "gld" || style == "nvk";
    if (ofTheirKind)
    {
      warnings.push_back(std::format(
          "[LAMMPS reader] fix {} is not converted (the integrators read are nve, nvt and npt, the thermostat next "
          "to them langevin)",
          words[3]));
    }
    return ofTheirKind;
  }
  if (words[2] != "all")
  {
    warnings.push_back(std::format(
        "[LAMMPS reader] fix {} acts on group '{}': it is converted as acting on all atoms", words[3], words[2]));
  }

  if (style == "langevin")
  {
    triple(4, "temperature", dynamics.langevinTemperature, dynamics.langevinDamping);
    return true;
  }
  if (!dynamics.integrator.empty() && dynamics.integrator != style)
  {
    warnings.push_back(std::format("[LAMMPS reader] fix {} after fix {}: the last one is converted", words[3],
                                   dynamics.integrator));
  }
  dynamics.integrator = style;
  for (std::size_t i = 4; i < words.size(); ++i)
  {
    if (words[i] == "temp")
    {
      triple(i + 1, "temperature", dynamics.temperature, dynamics.temperatureDamping);
    }
    else if (words[i] == "iso" || words[i] == "aniso" || words[i] == "tri" || words[i] == "x" || words[i] == "y" ||
             words[i] == "z")
    {
      if (words[i] != "iso") dynamics.flexibleCell = true;
      triple(i + 1, "pressure", dynamics.pressure, dynamics.pressureDamping);
    }
    else if (words[i] == "tchain" && number(i + 1))
    {
      dynamics.thermostatChain = static_cast<std::size_t>(*number(i + 1));
    }
    else if (words[i] == "pchain" && number(i + 1))
    {
      dynamics.barostatChain = static_cast<std::size_t>(*number(i + 1));
    }
  }
  return true;
}

InputScript readInputScript(const std::filesystem::path &path, std::vector<std::string> &warnings)
{
  InputScript script;
  std::ifstream file(path);
  if (!file) throw std::runtime_error(std::format("[LAMMPS reader] cannot open input script '{}'", path.string()));
  std::map<std::string, std::string> variables;

  std::string raw;
  std::string logical;
  while (std::getline(file, raw))
  {
    std::string line = stripComment(raw);
    if (line.ends_with('&'))
    {
      logical += line.substr(0, line.size() - 1) + " ";
      continue;
    }
    logical += line;
    std::vector<std::string> words = tokens(logical);
    logical.clear();
    if (words.empty()) continue;
    // 'variable name equal|index|string value' with one value, used as ${name} or v_name further on
    for (std::size_t i = 1; i < words.size(); ++i)
    {
      for (const auto &[name, value] : variables)
      {
        if (words[i] == "v_" + name) words[i] = value;
        // the search goes on after the value put in: a value that holds the reference itself ('${dir}/out' for
        // 'dir') was expanded forever. A variable of one letter is also referenced as '$x'.
        for (const std::string& reference :
             name.size() == 1 ? std::vector<std::string>{"${" + name + "}", "$" + name}
                              : std::vector<std::string>{"${" + name + "}"})
        {
          for (std::size_t at = words[i].find(reference); at != std::string::npos;
               at = words[i].find(reference, at + value.size()))
          {
            words[i].replace(at, reference.size(), value);
          }
        }
      }
    }
    const std::string &command = words[0];
    if (command == "variable" && words.size() == 4 &&
        (words[2] == "equal" || words[2] == "index" || words[2] == "string"))
    {
      variables[words[1]] = words[3];
      continue;
    }
    if (command == "atom_style" && words.size() > 1)
    {
      script.atomStyle = words[1];
      continue;
    }
    if (command == "timestep" && words.size() > 1 && isNumber(words[1]))
    {
      script.dynamics.timeStep = toDouble(words[1]);
      continue;
    }
    if (command == "run" && words.size() > 1 && isNumber(words[1]))
    {
      // 'run 0' evaluates the energy of the configuration: no step of a run
      if (toDouble(words[1]) >= 1.0) script.dynamics.runs.push_back(static_cast<std::size_t>(toDouble(words[1])));
      continue;
    }
    if (command == "fix" && words.size() > 3 && readDynamicsFix(words, script.dynamics, warnings)) continue;
    if (command == "velocity" && words.size() > 2)
    {
      if (words[2] == "create")
      {
        script.createsVelocities = true;
      }
      else if (words[2] != "zero")  // RASPA takes the drift of the centre of mass out itself
      {
        warnings.push_back(std::format(
            "[LAMMPS reader] 'velocity {} {}' is not applied: the velocities of the data file are kept as they are",
            words[1], words[2]));
      }
      continue;
    }

    if (command == "units" && words.size() > 1)
    {
      script.units = words[1];
    }
    else if (command == "bond_style")
    {
      script.bond = parseStyle(words);
      script.haveBond = true;
    }
    else if (command == "angle_style")
    {
      script.angle = parseStyle(words);
      script.haveAngle = true;
    }
    else if (command == "dihedral_style")
    {
      script.dihedral = parseStyle(words);
      script.haveDihedral = true;
    }
    else if (command == "improper_style")
    {
      script.improper = parseStyle(words);
      script.haveImproper = true;
    }
    else if (command == "pair_style")
    {
      script.pair = parseStyle(words);
      script.havePair = true;
    }
    else if (command == "pair_modify")
    {
      for (std::size_t i = 1; i + 1 < words.size(); ++i)
      {
        if (words[i] == "mix") script.mixing = words[i + 1];
        if (words[i] == "shift") script.shift = (words[i + 1] == "yes");
        if (words[i] == "tail") script.tail = (words[i + 1] == "yes");
      }
    }
    else if (command == "special_bonds")
    {
      for (std::size_t i = 1; i < words.size(); ++i)
      {
        if (words[i] == "lj/coul" && i + 3 < words.size())
        {
          std::array<double, 3> f{toDouble(words[i + 1]), toDouble(words[i + 2]), toDouble(words[i + 3])};
          script.specialLJ = f;
          script.specialCoulomb = f;
          i += 3;
        }
        else if (words[i] == "lj" && i + 3 < words.size())
        {
          script.specialLJ = std::array{toDouble(words[i + 1]), toDouble(words[i + 2]), toDouble(words[i + 3])};
          i += 3;
        }
        else if (words[i] == "coul" && i + 3 < words.size())
        {
          script.specialCoulomb = std::array{toDouble(words[i + 1]), toDouble(words[i + 2]), toDouble(words[i + 3])};
          i += 3;
        }
        else if (words[i] == "amber")
        {
          script.specialLJ = std::array{0.0, 0.0, 0.5};
          script.specialCoulomb = std::array{0.0, 0.0, 5.0 / 6.0};
        }
        else if (words[i] == "charmm")
        {
          script.specialLJ = std::array{0.0, 0.0, 0.0};
          script.specialCoulomb = std::array{0.0, 0.0, 0.0};
        }
        else if (words[i] == "dreiding")
        {
          script.specialLJ = std::array{0.0, 0.0, 1.0};
          script.specialCoulomb = std::array{0.0, 0.0, 1.0};
        }
        else if (words[i] == "fene")
        {
          script.specialLJ = std::array{0.0, 1.0, 1.0};
          script.specialCoulomb = std::array{0.0, 1.0, 1.0};
        }
      }
    }
    else if (command == "fix" && words.size() > 6 && (words[3] == "shake" || words[3] == "rattle"))
    {
      // fix ID group shake|rattle tolerance iterations N keyword values ... (b: bond types, a: angle types,
      // t: atom types, m: masses)
      if (script.shake) warnings.push_back("[LAMMPS reader] more than one 'fix shake': their constraints are merged");
      InputScript::Shake shake = script.shake.value_or(InputScript::Shake{});
      shake.tolerance = toDouble(words[4]);
      shake.iterations = static_cast<std::size_t>(std::stoul(words[5]));
      char keyword{};
      for (std::size_t i = 7; i < words.size(); ++i)
      {
        if (!isNumber(words[i]))
        {
          keyword = words[i].size() == 1 ? words[i][0] : '\0';
          // 'mol template' and other keywords carry no constraint list
          if (keyword != 'b' && keyword != 'a' && keyword != 't' && keyword != 'm')
          {
            keyword = '\0';
            if (words[i] != "mol")
              warnings.push_back(std::format("[LAMMPS reader] fix {}: keyword '{}' is not read", words[3], words[i]));
          }
          continue;
        }
        switch (keyword)
        {
          case 'b':
            shake.bondTypes.insert(static_cast<std::size_t>(std::stoul(words[i])));
            break;
          case 'a':
            shake.angleTypes.insert(static_cast<std::size_t>(std::stoul(words[i])));
            break;
          case 't':
            shake.atomTypes.insert(static_cast<std::size_t>(std::stoul(words[i])));
            break;
          case 'm':
            shake.masses.push_back(toDouble(words[i]));
            break;
          default:
            break;
        }
      }
      script.shake = shake;
    }
    else if (command == "kspace_style" && words.size() > 1)
    {
      script.kspace = words[1];
      if (words.size() > 2 && isNumber(words[2])) script.kspaceAccuracy = toDouble(words[2]);
    }
    else if (command == "pair_coeff" && words.size() >= 3)
    {
      if (words[1] == "*" || words[2] == "*" || words[1].contains('*') || words[2].contains('*'))
      {
        warnings.push_back(std::format("[LAMMPS reader] pair_coeff with wildcard '{}' is not expanded", raw));
        continue;
      }
      std::size_t i = static_cast<std::size_t>(std::stoul(words[1]));
      std::size_t j = static_cast<std::size_t>(std::stoul(words[2]));
      script.pairCoefficients[{std::min(i, j), std::max(i, j)}] = parseCoefficientLine(words, 3);
    }
    else if ((command == "bond_coeff" || command == "angle_coeff" || command == "dihedral_coeff" ||
              command == "improper_coeff") &&
             words.size() >= 2)  // 'bond_coeff 1' of the zero style carries no coefficients
    {
      if (!isNumber(words[1]))
      {
        warnings.push_back(std::format("[LAMMPS reader] {} with wildcard '{}' is not expanded", command, raw));
        continue;
      }
      std::size_t type = static_cast<std::size_t>(std::stoul(words[1]));
      CoefficientLine line = parseCoefficientLine(words, 2);
      if (command == "bond_coeff") script.bondCoefficients[type] = line;
      if (command == "angle_coeff") script.angleCoefficients[type] = line;
      if (command == "dihedral_coeff") script.dihedralCoefficients[type] = line;
      if (command == "improper_coeff") script.improperCoefficients[type] = line;
    }
  }
  return script;
}

DataFile readData(const std::filesystem::path &path, std::vector<std::string> &warnings,
                  const std::optional<std::string> &atomStyle)
{
  DataFile data;
  // the style of the script, unless the Atoms section names its own
  if (atomStyle.has_value()) data.atomStyle = *atomStyle;
  std::ifstream file(path);
  if (!file) throw std::runtime_error(std::format("[LAMMPS reader] cannot open data file '{}'", path.string()));

  static const std::set<std::string> knownSections{
      "Masses",         "Pair Coeffs",      "PairIJ Coeffs",  "Bond Coeffs",      "Angle Coeffs",
      "Dihedral Coeffs", "Improper Coeffs",  "Atoms",          "Velocities",       "Bonds",
      "Angles",         "Dihedrals",        "Impropers",      "BondBond Coeffs",  "BondAngle Coeffs",
      "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs", "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs",
      "BondBond13 Coeffs", "AngleAngle Coeffs", "Ellipsoids", "Lines", "Triangles", "Bodies", "Fragments"};

  std::string raw;
  std::getline(file, raw);  // title line

  // header
  std::string section;
  std::string sectionComment;
  std::streampos sectionStart;
  while (std::getline(file, raw))
  {
    std::string line = stripComment(raw);
    if (line.empty()) continue;
    std::vector<std::string> words = tokens(line);

    // section keyword?
    std::string candidate = line;
    if (knownSections.contains(candidate))
    {
      section = candidate;
      sectionComment = trailingComment(raw);
      break;
    }
    if (words.size() >= 2 && words[1] == "atoms") continue;
    if (words.size() >= 2 && (words[1] == "bonds" || words[1] == "angles" || words[1] == "dihedrals" ||
                              words[1] == "impropers"))
      continue;
    if (words.size() >= 3 && words[2] == "types")
    {
      std::size_t n = static_cast<std::size_t>(std::stoul(words[0]));
      if (words[1] == "atom") data.numberOfAtomTypes = n;
      if (words[1] == "bond") data.numberOfBondTypes = n;
      if (words[1] == "angle") data.numberOfAngleTypes = n;
      if (words[1] == "dihedral") data.numberOfDihedralTypes = n;
      if (words[1] == "improper") data.numberOfImproperTypes = n;
      continue;
    }
    if (words.size() >= 4 && words[2] == "xlo")
    {
      data.low[0] = toDouble(words[0]);
      data.high[0] = toDouble(words[1]);
      continue;
    }
    if (words.size() >= 4 && words[2] == "ylo")
    {
      data.low[1] = toDouble(words[0]);
      data.high[1] = toDouble(words[1]);
      continue;
    }
    if (words.size() >= 4 && words[2] == "zlo")
    {
      data.low[2] = toDouble(words[0]);
      data.high[2] = toDouble(words[1]);
      continue;
    }
    if (words.size() >= 6 && words[3] == "xy")
    {
      data.xy = toDouble(words[0]);
      data.xz = toDouble(words[1]);
      data.yz = toDouble(words[2]);
      data.triclinic = true;
      continue;
    }
    // "extra ... per atom" and similar header lines are ignored
  }

  auto readSection = [&](const std::string &name, const std::string &comment, auto &&onLine)
  {
    (void)name;
    (void)comment;
    std::string next;
    while (std::getline(file, next))
    {
      std::string line = stripComment(next);
      if (line.empty()) continue;
      if (knownSections.contains(line))
      {
        return std::make_pair(line, trailingComment(next));
      }
      onLine(tokens(line), next);
    }
    return std::make_pair(std::string{}, std::string{});
  };

  while (!section.empty())
  {
    std::pair<std::string, std::string> nextSection;
    if (section == "Masses")
    {
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &original)
                                {
                                  std::size_t type = static_cast<std::size_t>(std::stoul(words[0]));
                                  data.masses[type] = toDouble(words[1]);
                                  std::string name = trailingComment(original);
                                  if (!name.empty()) data.typeNames[type] = tokens(name)[0];
                                });
    }
    else if (section == "Pair Coeffs")
    {
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &)
                                {
                                  std::size_t type = static_cast<std::size_t>(std::stoul(words[0]));
                                  data.pairCoefficients[type] = parseCoefficientLine(words, 1);
                                });
    }
    else if (section == "PairIJ Coeffs")
    {
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &)
                                {
                                  std::size_t i = static_cast<std::size_t>(std::stoul(words[0]));
                                  std::size_t j = static_cast<std::size_t>(std::stoul(words[1]));
                                  data.pairIJCoefficients[{std::min(i, j), std::max(i, j)}] =
                                      parseCoefficientLine(words, 2);
                                });
    }
    else if (section == "Bond Coeffs" || section == "Angle Coeffs" || section == "Dihedral Coeffs" ||
             section == "Improper Coeffs")
    {
      std::map<std::size_t, CoefficientLine> &target =
          section == "Bond Coeffs"    ? data.bondCoefficients
          : section == "Angle Coeffs" ? data.angleCoefficients
          : section == "Dihedral Coeffs" ? data.dihedralCoefficients
                                         : data.improperCoefficients;
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &)
                                {
                                  std::size_t type = static_cast<std::size_t>(std::stoul(words[0]));
                                  target[type] = parseCoefficientLine(words, 1);
                                });
    }
    else if (section == "Atoms")
    {
      if (!sectionComment.empty()) data.atomStyle = tokens(sectionComment)[0];
      nextSection = readSection(
          section, sectionComment,
          [&](const std::vector<std::string> &words, const std::string &)
          {
            DataAtom atom;
            atom.id = static_cast<std::size_t>(std::stoul(words[0]));
            std::size_t k = 1;
            if (data.atomStyle == "full" || data.atomStyle == "molecular" || data.atomStyle == "bond" ||
                data.atomStyle == "angle")
            {
              atom.molecule = static_cast<std::size_t>(std::stoul(words[k++]));
            }
            else
            {
              atom.molecule = atom.id;  // atomic / charge styles: every atom its own molecule
            }
            atom.type = static_cast<std::size_t>(std::stoul(words[k++]));
            if (data.atomStyle == "full" || data.atomStyle == "charge") atom.charge = toDouble(words[k++]);
            if (k + 2 >= words.size())
              throw std::runtime_error(
                  std::format("[LAMMPS reader] cannot parse atom line '{}' (atom_style {})", words[0], data.atomStyle));
            atom.position = {toDouble(words[k]), toDouble(words[k + 1]), toDouble(words[k + 2])};
            k += 3;
            if (k + 2 < words.size() && isNumber(words[k]) && isNumber(words[k + 1]) && isNumber(words[k + 2]))
            {
              atom.image = {std::stoi(words[k]), std::stoi(words[k + 1]), std::stoi(words[k + 2])};
            }
            data.atoms.push_back(atom);
          });
    }
    else if (section == "Bonds" || section == "Angles" || section == "Dihedrals" || section == "Impropers")
    {
      std::vector<DataBonded> &target = section == "Bonds"       ? data.bonds
                                        : section == "Angles"    ? data.angles
                                        : section == "Dihedrals" ? data.dihedrals
                                                                 : data.impropers;
      std::size_t count = section == "Bonds" ? 2 : section == "Angles" ? 3 : 4;
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &)
                                {
                                  if (words.size() < 2 + count) return;
                                  DataBonded entry;
                                  entry.type = static_cast<std::size_t>(std::stoul(words[1]));
                                  for (std::size_t i = 0; i < count; ++i)
                                    entry.atoms.push_back(static_cast<std::size_t>(std::stoul(words[2 + i])));
                                  target.push_back(entry);
                                });
    }
    else if (section == "Velocities")
    {
      // 'id vx vy vz', followed by the angular velocities of the extended particles of some atom styles
      nextSection = readSection(section, sectionComment,
                                [&](const std::vector<std::string> &words, const std::string &)
                                {
                                  if (words.size() < 4) return;
                                  data.velocities[static_cast<std::size_t>(std::stoul(words[0]))] = {
                                      toDouble(words[1]), toDouble(words[2]), toDouble(words[3])};
                                });
    }
    else
    {
      warnings.push_back(std::format("[LAMMPS reader] section '{}' is ignored", section));
      nextSection = readSection(section, sectionComment, [](const std::vector<std::string> &, const std::string &) {});
    }
    section = nextSection.first;
    sectionComment = nextSection.second;
  }

  std::sort(data.atoms.begin(), data.atoms.end(), [](const DataAtom &a, const DataAtom &b) { return a.id < b.id; });
  return data;
}

// ---------------------------------------------------------------------------------------------------
// Conversion
// ---------------------------------------------------------------------------------------------------
double energyToKelvinFactor(const std::string& units, [[maybe_unused]] std::vector<std::string>& warnings)
{
  // The lammps_styles converters assume kcal/mol ('real'); 'metal' differs from it only in the energy (eV), with
  // lengths in Angstrom, masses in g/mol and charges in e. The other unit systems (si, cgs, lj, electron, micro,
  // nano) have lengths, masses and charges in other units as well, which the reader does not convert: scaling the
  // energies alone gave positions in metres and masses in kilograms read as Angstrom and g/mol.
  if (units == "real") return 1.0;
  if (units == "metal") return 23.060548;  // eV -> kcal/mol
  throw std::runtime_error(std::format(
      "[LAMMPS reader] units '{}' are not supported (only 'real' and 'metal'): its lengths, masses and charges would "
      "not be converted\n",
      units));
}

void scaleEnergies(RaspaTerm &term, std::span<const std::size_t> energyIndices, double factor)
{
  if (factor == 1.0) return;
  for (std::size_t index : energyIndices)
  {
    if (index < term.parameters.size()) term.parameters[index] *= factor;
  }
}

std::string styleFor(const StyleSpec &spec, const CoefficientLine &line)
{
  if (spec.hybrid()) return line.style.empty() ? (spec.hybridStyles.empty() ? "none" : spec.hybridStyles[0]) : line.style;
  return spec.name;
}

struct Template
{
  std::vector<std::size_t> pseudoAtomIndex{};  // per local atom
  std::vector<double> charges{};               // per local atom
  std::vector<long> chargeKeys{};              // per local atom, the charge rounded to 1e-6 (for the key)
  std::vector<std::array<std::size_t, 3>> bonds{};        // local a, local b, type
  std::vector<std::array<std::size_t, 4>> angles{};
  std::vector<std::array<std::size_t, 5>> dihedrals{};
  std::vector<std::array<std::size_t, 5>> impropers{};
  std::vector<std::array<double, 3>> positions{};  // first instance, unwrapped
  std::size_t count{};
  std::vector<std::size_t> firstAtomIds{};  // global id of the first atom of every instance
  std::vector<std::vector<std::size_t>> instances{};  // the atoms (indices in the data file) of every instance

  auto key() const { return std::tie(pseudoAtomIndex, chargeKeys, bonds, angles, dihedrals, impropers); }
};

std::string sanitizeName(std::string name)
{
  for (char &c : name)
  {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
  }
  return name;
}
}  // namespace

ReadResult readDataFile(const std::filesystem::path &dataFile, const std::optional<std::filesystem::path> &inputScript)
{
  ReadResult result;
  InputScript script;
  if (inputScript) script = readInputScript(*inputScript, result.warnings);
  DataFile data = readData(dataFile, result.warnings, script.atomStyle);

  const double unitFactor = energyToKelvinFactor(script.units, result.warnings);

  // box
  for (std::size_t k = 0; k < 3; ++k) result.boxLengths[k] = data.high[k] - data.low[k];
  if (data.triclinic)
  {
    double lx = result.boxLengths[0], ly = result.boxLengths[1], lz = result.boxLengths[2];
    double b = std::sqrt(ly * ly + data.xy * data.xy);
    double c = std::sqrt(lz * lz + data.xz * data.xz + data.yz * data.yz);
    result.boxLengths = {lx, b, c};
    const double toDeg = 180.0 / std::numbers::pi;
    result.boxAngles = {std::acos(std::clamp((data.xy * data.xz + ly * data.yz) / (b * c), -1.0, 1.0)) * toDeg,
                        std::acos(std::clamp(data.xz / c, -1.0, 1.0)) * toDeg,
                        std::acos(std::clamp(data.xy / b, -1.0, 1.0)) * toDeg};
  }

  // default styles
  if (!script.haveBond) script.bond.name = data.bondCoefficients.empty() ? "none" : "harmonic";
  if (!script.haveAngle) script.angle.name = data.angleCoefficients.empty() ? "none" : "harmonic";
  if (!script.haveDihedral) script.dihedral.name = data.dihedralCoefficients.empty() ? "none" : "nharmonic";
  if (!script.havePair) script.pair.name = "lj/cut";
  if (!inputScript)
  {
    result.warnings.push_back(
        "[LAMMPS reader] no input script given: assuming bond_style harmonic, angle_style harmonic, "
        "dihedral_style nharmonic, pair_style lj/cut, units real, special_bonds lj/coul 0 0 0");
  }
  for (const auto &[type, line] : script.bondCoefficients) data.bondCoefficients[type] = line;
  for (const auto &[type, line] : script.angleCoefficients) data.angleCoefficients[type] = line;
  for (const auto &[type, line] : script.dihedralCoefficients) data.dihedralCoefficients[type] = line;
  for (const auto &[type, line] : script.improperCoefficients) data.improperCoefficients[type] = line;
  for (const auto &[ij, line] : script.pairCoefficients)
  {
    if (ij.first == ij.second)
      data.pairCoefficients[ij.first] = line;
    else
      data.pairIJCoefficients[ij] = line;
  }

  // unwrap positions
  result.positions.resize(data.atoms.size());
  std::unordered_map<std::size_t, std::size_t> indexOfId;
  for (std::size_t i = 0; i < data.atoms.size(); ++i)
  {
    const DataAtom &atom = data.atoms[i];
    indexOfId[atom.id] = i;
    for (std::size_t k = 0; k < 3; ++k)
    {
      result.positions[i][k] = atom.position[k] + static_cast<double>(atom.image[k]) * (data.high[k] - data.low[k]);
    }
    if (data.triclinic)
    {
      result.positions[i][0] += static_cast<double>(atom.image[1]) * data.xy + static_cast<double>(atom.image[2]) * data.xz;
      result.positions[i][1] += static_cast<double>(atom.image[2]) * data.yz;
    }
  }

  // A molecule without image flags (all zero) that straddles a periodic boundary is broken in two by the wrap of
  // the data file, and its centre of mass and internal geometry come out wrong. Its atoms are joined along the
  // bonds, each bonded atom put at the nearest image of the atom it is bonded to (a bond is far shorter than half
  // the box), which holds for a molecule of any size; an atom without bonds is put nearest the first atom.
  {
    // the cell vectors a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz) of the LAMMPS box
    const double lx = data.high[0] - data.low[0];
    const double ly = data.high[1] - data.low[1];
    const double lz = data.high[2] - data.low[2];
    auto nearestImage = [&](std::array<double, 3> d)
    {
      const double sz = std::round(d[2] / lz);
      d[0] -= sz * data.xz;
      d[1] -= sz * data.yz;
      d[2] -= sz * lz;
      const double sy = std::round(d[1] / ly);
      d[0] -= sy * data.xy;
      d[1] -= sy * ly;
      d[0] -= std::round(d[0] / lx) * lx;
      return d;
    };

    std::map<std::size_t, std::vector<std::size_t>> moleculeAtoms;
    for (std::size_t i = 0; i < data.atoms.size(); ++i) moleculeAtoms[data.atoms[i].molecule].push_back(i);
    std::vector<std::vector<std::size_t>> neighbours(data.atoms.size());
    for (const DataBonded& bond : data.bonds)
    {
      if (bond.atoms.size() != 2 || !indexOfId.contains(bond.atoms[0]) || !indexOfId.contains(bond.atoms[1])) continue;
      const std::size_t a = indexOfId.at(bond.atoms[0]);
      const std::size_t b = indexOfId.at(bond.atoms[1]);
      neighbours[a].push_back(b);
      neighbours[b].push_back(a);
    }
    std::vector<bool> placed(data.atoms.size(), false);
    for (const auto& [molecule, members] : moleculeAtoms)
    {
      if (members.size() < 2) continue;
      const bool unflagged = std::ranges::all_of(
          members, [&](std::size_t i) { return data.atoms[i].image == std::array<int, 3>{0, 0, 0}; });
      if (!unflagged) continue;
      const std::size_t first = members.front();
      for (std::size_t start : members)
      {
        if (placed[start]) continue;
        if (start != first)
        {
          // a new bonded piece (or a lone atom): nearest the first atom of the molecule
          std::array<double, 3> d{};
          for (std::size_t k = 0; k < 3; ++k) d[k] = result.positions[start][k] - result.positions[first][k];
          d = nearestImage(d);
          for (std::size_t k = 0; k < 3; ++k) result.positions[start][k] = result.positions[first][k] + d[k];
        }
        placed[start] = true;
        std::vector<std::size_t> stack{start};
        while (!stack.empty())
        {
          const std::size_t i = stack.back();
          stack.pop_back();
          for (std::size_t j : neighbours[i])
          {
            if (placed[j] || data.atoms[j].molecule != molecule) continue;
            std::array<double, 3> d{};
            for (std::size_t k = 0; k < 3; ++k) d[k] = result.positions[j][k] - result.positions[i][k];
            d = nearestImage(d);
            for (std::size_t k = 0; k < 3; ++k) result.positions[j][k] = result.positions[i][k] + d[k];
            placed[j] = true;
            stack.push_back(j);
          }
        }
      }
    }
  }

  // pseudo-atoms: one per LAMMPS atom type. The charge of the pseudo-atom is the charge shared by all atoms
  // of the type; when the atoms of a type carry different charges the pseudo-atom gets charge zero and the
  // charges are written per atom in the component definitions (RASPA's per-atom charge overrides the
  // pseudo-atom default).
  struct PseudoAtomDefinition
  {
    std::size_t type;
    double charge;
    std::string name;
  };
  std::vector<PseudoAtomDefinition> pseudoAtoms;
  std::map<std::size_t, std::size_t> pseudoAtomIndexOf;
  auto chargeKey = [](double q) { return std::lround(q * 1.0e6); };
  std::vector<std::size_t> pseudoAtomOfAtom(data.atoms.size());
  std::map<std::size_t, std::set<long>> chargeVariantsOfType;
  for (std::size_t i = 0; i < data.atoms.size(); ++i)
  {
    const DataAtom &atom = data.atoms[i];
    auto it = pseudoAtomIndexOf.find(atom.type);
    if (it == pseudoAtomIndexOf.end())
    {
      std::string base = data.typeNames.contains(atom.type) ? data.typeNames.at(atom.type)
                                                             : std::format("T{}", atom.type);
      pseudoAtoms.push_back({atom.type, atom.charge, sanitizeName(base)});
      it = pseudoAtomIndexOf.emplace(atom.type, pseudoAtoms.size() - 1).first;
    }
    chargeVariantsOfType[atom.type].insert(chargeKey(atom.charge));
    pseudoAtomOfAtom[i] = it->second;
  }
  for (const auto &[type, variants] : chargeVariantsOfType)
  {
    if (variants.size() > 1)
    {
      pseudoAtoms[pseudoAtomIndexOf.at(type)].charge = 0.0;
      result.warnings.push_back(std::format(
          "[LAMMPS reader] atom type {} carries {} different charges; the pseudo-atom charge is set to zero and "
          "the charges are listed per atom in the component definitions",
          type, variants.size()));
    }
  }
  // pair coefficients for types that appear only in the coefficient tables (no atoms)
  for (const auto &[type, mass] : data.masses)
  {
    if (!pseudoAtomIndexOf.contains(type))
    {
      std::string base = data.typeNames.contains(type) ? data.typeNames.at(type) : std::format("T{}", type);
      pseudoAtoms.push_back({type, 0.0, sanitizeName(base)});
      pseudoAtomIndexOf.emplace(type, pseudoAtoms.size() - 1);
    }
  }
  // ensure unique names
  {
    std::map<std::string, std::size_t> seen;
    for (PseudoAtomDefinition &definition : pseudoAtoms)
    {
      std::size_t n = seen[definition.name]++;
      if (n > 0) definition.name = std::format("{}_{}", definition.name, n);
    }
  }

  // force field JSON
  nlohmann::json forceField;
  std::string mixing = script.mixing.empty() ? "geometric" : script.mixing;
  if (mixing == "arithmetic")
    forceField["MixingRule"] = "Lorentz-Berthelot";
  else if (mixing == "sixthpower")
    forceField["MixingRule"] = "SixthPower";
  else
    forceField["MixingRule"] = "Jorgensen";
  if (script.mixing.empty() && data.pairIJCoefficients.empty())
    result.warnings.push_back(
        "[LAMMPS reader] no 'pair_modify mix' given: LAMMPS defaults to geometric mixing (RASPA 'Jorgensen')");
  forceField["TruncationMethod"] = script.shift ? "shifted" : "truncated";
  forceField["TailCorrections"] = script.tail;
  double cutOffVDW = 12.0, cutOffCoulomb = 12.0;
  // The numbers of a pair style belong to its (sub-)styles: under hybrid/overlay every sub-style takes its own
  // ('hybrid/overlay coul/wolf 0.25 12.0 morse 12.0'). The damped Coulomb styles (coul/wolf, coul/dsf and the
  // fused lj/cut/coul/wolf, lj/cut/coul/dsf) take the damping alpha first and the cut-offs after it, a fused style
  // the VDW and then the Coulomb cut-off, a table style its interpolation keyword and table length.
  std::vector<std::pair<std::string, std::vector<double>>> pairStyleNumbers{};
  if (script.pair.hybrid())
  {
    for (const std::string &argument : script.pair.arguments)
    {
      if (isNumber(argument))
      {
        if (!pairStyleNumbers.empty()) pairStyleNumbers.back().second.push_back(toDouble(argument));
      }
      else if (pairStyleNumbers.empty() || pairStyleNumbers.back().first != "table" ||
               !(argument == "linear" || argument == "spline" || argument == "lookup" || argument == "bitmap"))
      {
        pairStyleNumbers.push_back({argument, {}});
      }
    }
  }
  else
  {
    pairStyleNumbers.push_back({script.pair.name, {}});
    for (const std::string &argument : script.pair.arguments)
      if (isNumber(argument)) pairStyleNumbers.back().second.push_back(toDouble(argument));
  }
  bool wolf = false;
  bool dampedShiftedForce = false;
  std::optional<double> dampingAlpha{};
  std::optional<double> styleCutOffVDW{};
  std::optional<double> styleCutOffCoulomb{};
  for (auto [style, numbers] : pairStyleNumbers)
  {
    if (style == "table" || numbers.empty()) continue;
    if (style.contains("coul/wolf") || style.contains("coul/dsf"))
    {
      (style.contains("coul/wolf") ? wolf : dampedShiftedForce) = true;
      dampingAlpha = numbers.front();
      numbers.erase(numbers.begin());
      if (numbers.empty()) continue;
    }
    if (style.starts_with("coul"))
    {
      if (!styleCutOffCoulomb) styleCutOffCoulomb = numbers[0];
    }
    else if (style.contains("coul"))
    {
      if (!styleCutOffVDW) styleCutOffVDW = numbers[0];
      if (!styleCutOffCoulomb) styleCutOffCoulomb = numbers.size() > 1 ? numbers[1] : numbers[0];
    }
    else if (!styleCutOffVDW)
    {
      styleCutOffVDW = numbers[0];
    }
  }
  if (styleCutOffVDW) cutOffVDW = *styleCutOffVDW;
  cutOffCoulomb = styleCutOffCoulomb.value_or(styleCutOffVDW.value_or(cutOffCoulomb));
  if (!styleCutOffVDW && styleCutOffCoulomb) cutOffVDW = *styleCutOffCoulomb;
  forceField["CutOffVDW"] = cutOffVDW;
  forceField["CutOffCoulomb"] = cutOffCoulomb;

  std::string pairStyleName = script.pair.name;
  bool anyCharge = std::any_of(data.atoms.begin(), data.atoms.end(), [](const DataAtom &a) { return a.charge != 0.0; });
  bool coulomb = pairStyleName.contains("coul") || script.pair.hybrid() &&
                                                       std::any_of(script.pair.hybridStyles.begin(),
                                                                   script.pair.hybridStyles.end(),
                                                                   [](const std::string &s) { return s.contains("coul"); });
  if (anyCharge && !coulomb)
    result.warnings.push_back("[LAMMPS reader] atoms carry charges but the pair style has no Coulomb part");

  forceField["PseudoAtoms"] = nlohmann::json::array();
  for (const PseudoAtomDefinition &definition : pseudoAtoms)
  {
    double mass = data.masses.contains(definition.type) ? data.masses.at(definition.type) : 0.0;
    std::string element = elementForMass(mass);
    forceField["PseudoAtoms"].push_back({{"name", definition.name},
                                         {"framework", false},
                                         {"print_to_output", true},
                                         {"element", element},
                                         {"print_as", element},
                                         {"mass", mass},
                                         {"charge", definition.charge}});
  }

  auto pairStyleOf = [&](const CoefficientLine &line) -> std::string
  {
    if (script.pair.hybrid()) return line.style.empty() ? (script.pair.hybridStyles.empty() ? "" : script.pair.hybridStyles[0]) : line.style;
    return pairStyleName;
  };

  forceField["SelfInteractions"] = nlohmann::json::array();
  std::set<std::size_t> missingSelf;
  for (const PseudoAtomDefinition &definition : pseudoAtoms)
  {
    auto it = data.pairCoefficients.find(definition.type);
    if (it == data.pairCoefficients.end())
    {
      auto ij = data.pairIJCoefficients.find({definition.type, definition.type});
      if (ij != data.pairIJCoefficients.end()) it = data.pairCoefficients.emplace(definition.type, ij->second).first;
    }
    if (it == data.pairCoefficients.end())
    {
      missingSelf.insert(definition.type);
      forceField["SelfInteractions"].push_back(
          {{"name", definition.name}, {"type", "lennard-jones"}, {"parameters", {0.0, 1.0}}});
      continue;
    }
    std::optional<std::array<double, 2>> lj = lennardJonesFromLammps(pairStyleOf(it->second), it->second.values);
    if (!lj)
    {
      result.warnings.push_back(std::format("[LAMMPS reader] pair style '{}' for type {} is not converted (set to zero)",
                                            pairStyleOf(it->second), definition.type));
      forceField["SelfInteractions"].push_back(
          {{"name", definition.name}, {"type", "lennard-jones"}, {"parameters", {0.0, 1.0}}});
      continue;
    }
    forceField["SelfInteractions"].push_back({{"name", definition.name},
                                              {"type", "lennard-jones"},
                                              {"parameters", {(*lj)[0] * unitFactor, (*lj)[1]}}});
  }
  for (std::size_t type : missingSelf)
    result.warnings.push_back(std::format("[LAMMPS reader] no pair coefficients for atom type {}", type));

  if (!data.pairIJCoefficients.empty())
  {
    forceField["BinaryInteractions"] = nlohmann::json::array();
    for (const auto &[ij, line] : data.pairIJCoefficients)
    {
      if (ij.first == ij.second) continue;
      std::optional<std::array<double, 2>> lj = lennardJonesFromLammps(pairStyleOf(line), line.values);
      if (!lj) continue;
      for (const PseudoAtomDefinition &a : pseudoAtoms)
      {
        if (a.type != ij.first) continue;
        for (const PseudoAtomDefinition &b : pseudoAtoms)
        {
          if (b.type != ij.second) continue;
          forceField["BinaryInteractions"].push_back({{"names", {a.name, b.name}},
                                                      {"type", "lennard-jones"},
                                                      {"parameters", {(*lj)[0] * unitFactor, (*lj)[1]}}});
        }
      }
    }
  }
  result.forceField = forceField;

  // bonded terms per type
  auto convertAll = [&](const std::map<std::size_t, CoefficientLine> &lines, const StyleSpec &spec, auto converter,
                        std::span<const std::size_t> energyIndices, const char *label)
  {
    std::map<std::size_t, RaspaTerm> converted;
    for (const auto &[type, line] : lines)
    {
      std::string style = styleFor(spec, line);
      std::optional<RaspaTerm> term = converter(style, std::span<const double>(line.values));
      if (!term)
      {
        result.warnings.push_back(
            std::format("[LAMMPS reader] {} type {} with style '{}' cannot be converted; term dropped", label, type, style));
        continue;
      }
      // a FIXED term holds a length or an angle and no energy: the energy units (eV under 'units metal') do not apply,
      // and scaling turned a 1.54 A bond into 35.5 A
      if (term->type != "FIXED") scaleEnergies(*term, energyIndices, unitFactor);
      if (!term->note.empty()) result.warnings.push_back(std::format("[LAMMPS reader] {} type {}: {}", label, type, term->note));
      converted[type] = *term;
    }
    return converted;
  };
  static constexpr std::array<std::size_t, 1> bondEnergy{0};
  static constexpr std::array<std::size_t, 3> bendEnergy{0, 2, 3};
  static constexpr std::array<std::size_t, 6> torsionEnergyAll{0, 1, 2, 3, 4, 5};
  std::map<std::size_t, RaspaTerm> bondTerms = convertAll(data.bondCoefficients, script.bond, bondFromLammps, bondEnergy, "bond");
  std::map<std::size_t, RaspaTerm> bendTerms = convertAll(data.angleCoefficients, script.angle, bendFromLammps, bendEnergy, "angle");
  // angle_style charmm (K theta0 K_ub r_ub) carries a Urey-Bradley 1-3 term next to the bend: it becomes a RASPA
  // Urey-Bradley potential on the outer atoms. A charmm angle with K = 0 (how RASPA exports a Urey-Bradley term next
  // to the bend of the same triple) has no bend.
  std::map<std::size_t, RaspaTerm> ureyBradleyTerms;
  for (const auto &[type, line] : data.angleCoefficients)
  {
    if (styleFor(script.angle, line) != "charmm" || line.values.size() < 4) continue;
    if (line.values[2] != 0.0)
    {
      std::optional<RaspaTerm> term = bondFromLammps("harmonic", std::span<const double>(line.values).subspan(2, 2));
      scaleEnergies(*term, bondEnergy, unitFactor);
      ureyBradleyTerms[type] = *term;
    }
    if (line.values[0] == 0.0) bendTerms.erase(type);
  }
  std::map<std::size_t, RaspaTerm> torsionTerms;
  for (const auto &[type, line] : data.dihedralCoefficients)
  {
    std::string style = styleFor(script.dihedral, line);
    std::optional<RaspaTerm> term = torsionFromLammps(style, line.values);
    if (!term)
    {
      result.warnings.push_back(std::format(
          "[LAMMPS reader] dihedral type {} with style '{}' cannot be converted; term dropped", type, style));
      continue;
    }
    if (unitFactor != 1.0)
    {
      // CVFF carries (K, n, d), HARMONIC (K, phi_0 in degrees): only their first parameter is an energy
      if (term->type == "CVFF" || term->type == "HARMONIC")
        term->parameters[0] *= unitFactor;
      else
        scaleEnergies(*term, torsionEnergyAll, unitFactor);
    }
    torsionTerms[type] = *term;
  }
  if (!data.improperCoefficients.empty())
    result.warnings.push_back("[LAMMPS reader] improper terms are not converted");

  // molecules -> templates
  std::map<std::size_t, std::vector<std::size_t>> atomsOfMolecule;  // molecule id -> atom indices (sorted by id)
  for (std::size_t i = 0; i < data.atoms.size(); ++i) atomsOfMolecule[data.atoms[i].molecule].push_back(i);

  std::map<std::size_t, std::vector<const DataBonded *>> bondsOfMolecule, anglesOfMolecule, dihedralsOfMolecule,
      impropersOfMolecule;
  auto assign = [&](const std::vector<DataBonded> &list, std::map<std::size_t, std::vector<const DataBonded *>> &target,
                    const char *label)
  {
    for (const DataBonded &entry : list)
    {
      std::size_t molecule = data.atoms[indexOfId.at(entry.atoms[0])].molecule;
      for (std::size_t id : entry.atoms)
      {
        if (data.atoms[indexOfId.at(id)].molecule != molecule)
        {
          result.warnings.push_back(std::format("[LAMMPS reader] {} between different molecule ids ignored", label));
          molecule = std::numeric_limits<std::size_t>::max();
          break;
        }
      }
      if (molecule != std::numeric_limits<std::size_t>::max()) target[molecule].push_back(&entry);
    }
  };
  assign(data.bonds, bondsOfMolecule, "bond");
  assign(data.angles, anglesOfMolecule, "angle");
  assign(data.dihedrals, dihedralsOfMolecule, "dihedral");
  assign(data.impropers, impropersOfMolecule, "improper");

  std::vector<Template> templates;
  for (const auto &[moleculeId, atomIndices] : atomsOfMolecule)
  {
    Template candidate;
    std::unordered_map<std::size_t, std::size_t> local;  // global id -> local index
    for (std::size_t n = 0; n < atomIndices.size(); ++n)
    {
      local[data.atoms[atomIndices[n]].id] = n;
      candidate.pseudoAtomIndex.push_back(pseudoAtomOfAtom[atomIndices[n]]);
      candidate.charges.push_back(data.atoms[atomIndices[n]].charge);
      candidate.chargeKeys.push_back(chargeKey(data.atoms[atomIndices[n]].charge));
      candidate.positions.push_back(result.positions[atomIndices[n]]);
    }
    auto localOf = [&](std::size_t id) { return local.at(id); };
    for (const DataBonded *bond : bondsOfMolecule[moleculeId])
    {
      std::size_t a = localOf(bond->atoms[0]), b = localOf(bond->atoms[1]);
      candidate.bonds.push_back({std::min(a, b), std::max(a, b), bond->type});
    }
    for (const DataBonded *angle : anglesOfMolecule[moleculeId])
    {
      std::size_t a = localOf(angle->atoms[0]), b = localOf(angle->atoms[1]), c = localOf(angle->atoms[2]);
      if (a > c) std::swap(a, c);
      candidate.angles.push_back({a, b, c, angle->type});
    }
    for (const DataBonded *dihedral : dihedralsOfMolecule[moleculeId])
    {
      std::size_t a = localOf(dihedral->atoms[0]), b = localOf(dihedral->atoms[1]), c = localOf(dihedral->atoms[2]),
                  d = localOf(dihedral->atoms[3]);
      if (a > d)
      {
        std::swap(a, d);
        std::swap(b, c);
      }
      candidate.dihedrals.push_back({a, b, c, d, dihedral->type});
    }
    for (const DataBonded *improper : impropersOfMolecule[moleculeId])
      candidate.impropers.push_back({localOf(improper->atoms[0]), localOf(improper->atoms[1]), localOf(improper->atoms[2]),
                                     localOf(improper->atoms[3]), improper->type});
    std::sort(candidate.bonds.begin(), candidate.bonds.end());
    std::sort(candidate.angles.begin(), candidate.angles.end());
    std::sort(candidate.dihedrals.begin(), candidate.dihedrals.end());
    std::sort(candidate.impropers.begin(), candidate.impropers.end());

    bool found = false;
    for (Template &existing : templates)
    {
      if (existing.key() == candidate.key())
      {
        existing.count++;
        existing.firstAtomIds.push_back(data.atoms[atomIndices[0]].id);
        existing.instances.push_back(atomIndices);
        found = true;
        break;
      }
    }
    if (!found)
    {
      candidate.count = 1;
      candidate.firstAtomIds.push_back(data.atoms[atomIndices[0]].id);
      candidate.instances.push_back(atomIndices);
      templates.push_back(std::move(candidate));
    }
  }

  // components
  double lj14 = script.specialLJ ? (*script.specialLJ)[2] : 0.0;
  double coul14 = script.specialCoulomb ? (*script.specialCoulomb)[2] : 0.0;
  if ((script.specialLJ && ((*script.specialLJ)[0] != 0.0 || (*script.specialLJ)[1] != 0.0)) ||
      (script.specialCoulomb && ((*script.specialCoulomb)[0] != 0.0 || (*script.specialCoulomb)[1] != 0.0)))
    result.warnings.push_back(
        "[LAMMPS reader] non-zero 1-2 / 1-3 special_bonds factors cannot be represented in RASPA and are ignored");

  std::size_t componentCounter = 0;
  std::map<std::string, std::size_t> componentNamesSeen;
  for (const Template &tpl : templates)
  {
    ReadComponent component;
    component.count = tpl.count;
    component.atomsPerMolecule = tpl.pseudoAtomIndex.size();
    component.name = templates.size() == 1 ? "molecule" : std::format("molecule{}", componentCounter + 1);
    if (tpl.pseudoAtomIndex.size() == 1) component.name = pseudoAtoms[tpl.pseudoAtomIndex[0]].name;
    // single-atom templates of one type but different charges would share a name
    if (std::size_t n = componentNamesSeen[component.name]++; n > 0)
      component.name = std::format("{}_{}", component.name, n);
    ++componentCounter;

    nlohmann::json definition;
    definition["CriticalTemperature"] = 500.0;
    definition["CriticalPressure"] = 5.0e6;
    definition["AcentricFactor"] = 0.3;
    definition["_comment"] = std::format(
        "converted from LAMMPS data file '{}'; critical constants are placeholders, geometry is the first instance",
        dataFile.filename().string());

    // centre the template geometry
    std::array<double, 3> centre{};
    for (const std::array<double, 3> &p : tpl.positions)
      for (std::size_t k = 0; k < 3; ++k) centre[k] += p[k] / static_cast<double>(tpl.positions.size());
    definition["PseudoAtoms"] = nlohmann::json::array();
    for (std::size_t n = 0; n < tpl.pseudoAtomIndex.size(); ++n)
    {
      const PseudoAtomDefinition &pseudoAtom = pseudoAtoms[tpl.pseudoAtomIndex[n]];
      nlohmann::json entry = {
          pseudoAtom.name,
          {tpl.positions[n][0] - centre[0], tpl.positions[n][1] - centre[1], tpl.positions[n][2] - centre[2]}};
      // the per-atom charge overrides the pseudo-atom default; written only when it differs
      if (tpl.chargeKeys[n] != chargeKey(pseudoAtom.charge)) entry.push_back(tpl.charges[n]);
      definition["PseudoAtoms"].push_back(entry);
    }

    // geometry of the first instance: the reference of the FIXED terms of a 'zero' style without coefficients
    auto separation = [&](std::size_t from, std::size_t to)
    {
      return std::array<double, 3>{tpl.positions[to][0] - tpl.positions[from][0],
                                   tpl.positions[to][1] - tpl.positions[from][1],
                                   tpl.positions[to][2] - tpl.positions[from][2]};
    };
    auto dot = [](const std::array<double, 3> &u, const std::array<double, 3> &v)
    { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
    auto cross = [](const std::array<double, 3> &u, const std::array<double, 3> &v)
    {
      return std::array<double, 3>{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                   u[0] * v[1] - u[1] * v[0]};
    };
    constexpr double degrees = 180.0 / std::numbers::pi;
    auto fixedReference = [](RaspaTerm term, double value)
    {
      if (term.type == "FIXED" && term.parameters.empty()) term.parameters = {value};
      return term;
    };

    // 'fix shake' constrains a bond by its type, or by the type or the mass of either of its atoms (LAMMPS matches
    // a mass to within 0.1), and an angle by its type when its two bonds are constrained. The constrained length
    // and angle are the equilibrium values of the harmonic term, as in LAMMPS, or else those of the geometry.
    std::set<std::pair<std::size_t, std::size_t>> shakenBonds;
    const auto shakenAtom = [&](std::size_t local)
    {
      const std::size_t type = pseudoAtoms[tpl.pseudoAtomIndex[local]].type;
      if (script.shake->atomTypes.contains(type)) return true;
      if (!data.masses.contains(type)) return false;
      return std::ranges::any_of(script.shake->masses,
                                 [&](double mass) { return std::abs(mass - data.masses.at(type)) <= 0.1; });
    };
    if (!tpl.bonds.empty())
    {
      definition["Connectivity"] = nlohmann::json::array();
      definition["Bonds"] = nlohmann::json::array();
      for (const std::array<std::size_t, 3> &bond : tpl.bonds)
      {
        definition["Connectivity"].push_back({bond[0], bond[1]});
        auto it = bondTerms.find(bond[2]);
        const std::array<double, 3> r = separation(bond[0], bond[1]);
        if (script.shake &&
            (script.shake->bondTypes.contains(bond[2]) || shakenAtom(bond[0]) || shakenAtom(bond[1])))
        {
          const bool harmonic = it != bondTerms.end() && it->second.type == "HARMONIC" && it->second.parameters.size() > 1;
          const bool fixed = it != bondTerms.end() && it->second.type == "FIXED" && !it->second.parameters.empty();
          const double length =
              harmonic ? it->second.parameters[1] : fixed ? it->second.parameters[0] : std::sqrt(dot(r, r));
          definition["Bonds"].push_back({{bond[0], bond[1]}, "FIXED", {length}});
          shakenBonds.insert(std::minmax(bond[0], bond[1]));
          continue;
        }
        if (it == bondTerms.end()) continue;
        const RaspaTerm term = fixedReference(it->second, std::sqrt(dot(r, r)));
        definition["Bonds"].push_back({{bond[0], bond[1]}, term.type, term.parameters});
      }
    }
    if (!tpl.angles.empty())
    {
      // RASPA takes one bend per triple: several LAMMPS angle types on one triple (a bend and a separate Urey-Bradley
      // angle) keep the first bend
      definition["Bends"] = nlohmann::json::array();
      std::set<std::array<std::size_t, 3>> bentTriples;
      for (const std::array<std::size_t, 4> &angle : tpl.angles)
      {
        if (auto it = ureyBradleyTerms.find(angle[3]); it != ureyBradleyTerms.end())
        {
          if (!definition.contains("UreyBradleys")) definition["UreyBradleys"] = nlohmann::json::array();
          definition["UreyBradleys"].push_back({{angle[0], angle[2]}, it->second.type, it->second.parameters});
        }
        auto it = bendTerms.find(angle[3]);
        if (it == bendTerms.end()) continue;
        if (!bentTriples.insert({angle[0], angle[1], angle[2]}).second)
        {
          result.warnings.push_back(std::format("[LAMMPS reader] angle {}-{}-{}: more than one bend; the first is kept",
                                                angle[0], angle[1], angle[2]));
          continue;
        }
        const std::array<double, 3> u = separation(angle[1], angle[0]), v = separation(angle[1], angle[2]);
        if (script.shake && script.shake->angleTypes.contains(angle[3]))
        {
          if (shakenBonds.contains(std::minmax(angle[0], angle[1])) &&
              shakenBonds.contains(std::minmax(angle[1], angle[2])))
          {
            const bool harmonic = it->second.type == "HARMONIC" && it->second.parameters.size() > 1;
            const bool fixed = it->second.type == "FIXED" && !it->second.parameters.empty();
            const double theta =
                harmonic ? it->second.parameters[1]
                         : fixed ? it->second.parameters[0]
                                 : degrees * std::acos(std::clamp(dot(u, v) / std::sqrt(dot(u, u) * dot(v, v)), -1.0, 1.0));
            definition["Bends"].push_back({{angle[0], angle[1], angle[2]}, "FIXED", {theta}});
            continue;
          }
          result.warnings.push_back(std::format(
              "[LAMMPS reader] fix shake: angle {}-{}-{} of type {} is not constrained, its two bonds are not both",
              angle[0], angle[1], angle[2], angle[3]));
        }
        const RaspaTerm term =
            fixedReference(it->second, degrees * std::acos(std::clamp(dot(u, v) / std::sqrt(dot(u, u) * dot(v, v)),
                                                                      -1.0, 1.0)));
        definition["Bends"].push_back({{angle[0], angle[1], angle[2]}, term.type, term.parameters});
      }
    }
    if (!tpl.dihedrals.empty())
    {
      // RASPA takes one torsion per quadruple: the terms of a multi-term LAMMPS dihedral (several dihedral types on
      // one quadruple, as CHARMM/AMBER files write them) are summed into one cosine series c_0 + sum_n c_n cos(n phi)
      // (TRAPPE_EXTENDED, n <= 4) when every term is one (TRAPPE_EXTENDED, or CVFF with a phase of 0 or 180 degrees)
      auto cosineSeries = [](const RaspaTerm &term) -> std::optional<std::array<double, 5>>
      {
        std::array<double, 5> c{};
        if (term.type == "TRAPPE_EXTENDED")
        {
          for (std::size_t n = 0; n < 5 && n < term.parameters.size(); ++n) c[n] = term.parameters[n];
          return c;
        }
        if (term.type == "CVFF" && term.parameters.size() >= 3)
        {
          const long n = std::lround(term.parameters[1]);
          const double phase = std::remainder(term.parameters[2], 360.0);
          if (n < 1 || n > 4 || std::abs(term.parameters[1] - static_cast<double>(n)) > 1e-9) return std::nullopt;
          if (std::abs(phase) > 1e-9 && std::abs(std::abs(phase) - 180.0) > 1e-9) return std::nullopt;
          c[0] = term.parameters[0];
          c[static_cast<std::size_t>(n)] = std::abs(phase) < 1e-9 ? term.parameters[0] : -term.parameters[0];
          return c;
        }
        return std::nullopt;
      };
      std::map<std::array<std::size_t, 4>, std::vector<RaspaTerm>> termsOfQuadruple;
      std::vector<std::array<std::size_t, 4>> quadruples;
      for (const std::array<std::size_t, 5> &dihedral : tpl.dihedrals)
      {
        auto it = torsionTerms.find(dihedral[4]);
        if (it == torsionTerms.end()) continue;
        const std::array<std::size_t, 4> quadruple{dihedral[0], dihedral[1], dihedral[2], dihedral[3]};
        if (!termsOfQuadruple.contains(quadruple)) quadruples.push_back(quadruple);
        termsOfQuadruple[quadruple].push_back(it->second);
      }
      definition["Torsions"] = nlohmann::json::array();
      for (const std::array<std::size_t, 4> &quadruple : quadruples)
      {
        const std::vector<RaspaTerm> &terms = termsOfQuadruple[quadruple];
        RaspaTerm term = terms.front();
        if (terms.size() > 1)
        {
          std::array<double, 5> sum{};
          bool summable = true;
          for (const RaspaTerm &part : terms)
          {
            std::optional<std::array<double, 5>> c = cosineSeries(part);
            if (!c)
            {
              summable = false;
              break;
            }
            for (std::size_t n = 0; n < 5; ++n) sum[n] += (*c)[n];
          }
          if (summable)
          {
            term = RaspaTerm{"TRAPPE_EXTENDED", std::vector<double>(sum.begin(), sum.end())};
          }
          else
          {
            result.warnings.push_back(std::format(
                "[LAMMPS reader] dihedral {}-{}-{}-{}: {} terms that do not sum to one cosine series; the first is kept",
                quadruple[0], quadruple[1], quadruple[2], quadruple[3], terms.size()));
          }
        }
        // IUPAC: the angle between the planes (a, b, c) and (b, c, d), positive clockwise looking along b -> c
        const std::array<double, 3> b1 = separation(quadruple[0], quadruple[1]),
                                    b2 = separation(quadruple[1], quadruple[2]),
                                    b3 = separation(quadruple[2], quadruple[3]);
        const std::array<double, 3> n1 = cross(b1, b2), n2 = cross(b2, b3);
        const double phi = degrees * std::atan2(std::sqrt(dot(b2, b2)) * dot(b1, n2), dot(n1, n2));
        term = fixedReference(term, phi);
        definition["Torsions"].push_back(
            {{quadruple[0], quadruple[1], quadruple[2], quadruple[3]}, term.type, term.parameters});
      }
    }
    // The component reader asks a potential of every bend and torsion of the bond graph; LAMMPS has none for an
    // angle or dihedral the file leaves out, which is no force at all. Those get a potential without energy (left
    // out, the component of the file was refused).
    if (!tpl.bonds.empty())
    {
      std::vector<std::set<std::size_t>> neighbours(tpl.pseudoAtomIndex.size());
      for (const std::array<std::size_t, 3> &bond : tpl.bonds)
      {
        neighbours[bond[0]].insert(bond[1]);
        neighbours[bond[1]].insert(bond[0]);
      }
      const auto defined = [&](const char *key)
      {
        std::set<std::vector<std::size_t>> listed;
        if (definition.contains(key))
        {
          for (const auto &item : definition[key])
          {
            std::vector<std::size_t> ids = item[0].get<std::vector<std::size_t>>();
            listed.insert(ids);
            std::ranges::reverse(ids);
            listed.insert(ids);
          }
        }
        return listed;
      };
      const std::set<std::vector<std::size_t>> bends = defined("Bends");
      for (std::size_t b = 0; b < neighbours.size(); ++b)
      {
        for (std::size_t a : neighbours[b])
        {
          for (std::size_t c : neighbours[b])
          {
            if (a >= c || bends.contains({a, b, c})) continue;
            if (!definition.contains("Bends")) definition["Bends"] = nlohmann::json::array();
            definition["Bends"].push_back({{a, b, c}, "HARMONIC", {0.0, 0.0}});
          }
        }
      }
      const std::set<std::vector<std::size_t>> torsions = defined("Torsions");
      for (const std::array<std::size_t, 3> &bond : tpl.bonds)
      {
        const std::size_t b = bond[0], c = bond[1];
        for (std::size_t a : neighbours[b])
        {
          for (std::size_t d : neighbours[c])
          {
            if (a == c || d == b || a == d || torsions.contains({a, b, c, d})) continue;
            if (!definition.contains("Torsions")) definition["Torsions"] = nlohmann::json::array();
            definition["Torsions"].push_back({{a, b, c, d}, "THREE_COSINE", {0.0, 0.0, 0.0}});
          }
        }
      }
    }
    if (!tpl.bonds.empty())
    {
      if (lj14 != 0.0) definition["Intra14VanDerWaalsScalingValue"] = lj14;
      if (coul14 != 0.0) definition["Intra14ChargeChargeScalingValue"] = coul14;
    }
    component.definition = definition;
    // The configuration of the data file: the positions of the atoms of every molecule of the component, molecule
    // after molecule in the atom order of the template, each molecule in one piece, relative to the origin of the
    // box. It is what the simulation starts from ('RestartFileName').
    nlohmann::json positions = nlohmann::json::array();
    for (const std::vector<std::size_t> &instance : tpl.instances)
    {
      for (const std::size_t atom : instance)
      {
        positions.push_back({result.positions[atom][0] - data.low[0], result.positions[atom][1] - data.low[1],
                             result.positions[atom][2] - data.low[2]});
      }
    }
    result.configuration[component.name] = std::move(positions);

    result.components.push_back(component);
  }

  // The velocities of the data file, in the same order, next to the positions: molecular dynamics starts from
  // them. LAMMPS 'real' units have Angstrom per femtosecond, 'metal' units Angstrom per picosecond, as RASPA.
  if (!data.velocities.empty())
  {
    const bool complete = std::ranges::all_of(
        data.atoms, [&](const DataAtom &atom) { return data.velocities.contains(atom.id); });
    if (script.createsVelocities)
    {
      result.warnings.push_back(
          "[LAMMPS reader] the script draws new velocities ('velocity ... create'): those of the data file are not "
          "carried over, molecular dynamics draws them at the temperature of the run");
    }
    else if (!complete)
    {
      result.warnings.push_back(
          "[LAMMPS reader] the 'Velocities' section does not list every atom: the velocities are not carried over");
    }
    else
    {
      const double perPicosecond = script.units == "metal" ? 1.0 : 1.0e3;
      nlohmann::json velocities = nlohmann::json::object();
      for (std::size_t c = 0; c < templates.size(); ++c)
      {
        nlohmann::json list = nlohmann::json::array();
        for (const std::vector<std::size_t> &instance : templates[c].instances)
        {
          for (const std::size_t atom : instance)
          {
            const std::array<double, 3> &v = data.velocities.at(data.atoms[atom].id);
            list.push_back({v[0] * perPicosecond, v[1] * perPicosecond, v[2] * perPicosecond});
          }
        }
        velocities[result.components[c].name] = std::move(list);
      }
      result.configuration["Velocities"] = std::move(velocities);
    }
  }

  // simulation skeleton
  nlohmann::json simulation;
  simulation["SimulationType"] = "MonteCarlo";
  simulation["NumberOfProductionCycles"] = 1000;
  simulation["NumberOfInitializationCycles"] = 0;
  simulation["PrintEvery"] = 100;
  simulation["ForceField"] = ".";
  nlohmann::json system;
  system["Type"] = "Box";
  system["BoxLengths"] = result.boxLengths;
  if (data.triclinic) system["BoxAngles"] = result.boxAngles;
  system["ExternalTemperature"] = 300.0;
  system["RestartFileName"] = "initial_configuration";
  system["ChargeMethod"] = !(coulomb && anyCharge) ? "None"
                           : wolf                  ? "Wolf"
                           : dampedShiftedForce    ? "DampedShiftedForce"
                                                   : "Ewald";
  if (dampingAlpha && coulomb && anyCharge)
  {
    // the damping parameter of Wolf and DSF is the Ewald alpha of the force field
    result.forceField["EwaldParameters"] = nlohmann::json::array({*dampingAlpha, 1, 1, 1});
  }
  if (script.kspaceAccuracy && coulomb && anyCharge) result.forceField["EwaldPrecision"] = *script.kspaceAccuracy;
  // The run of the script: molecular dynamics in the ensemble of its integrating fix. LAMMPS 'real' units have
  // femtoseconds and atmospheres, 'metal' units picoseconds and bar; RASPA takes picoseconds and Pascal.
  const InputScript::Dynamics &dynamics = script.dynamics;
  if (!dynamics.integrator.empty())
  {
    const bool metal = script.units == "metal";
    const double picoseconds = metal ? 1.0 : 1.0e-3;
    const double pascal = metal ? 1.0e5 : 101325.0;
    const bool langevin = dynamics.langevinTemperature.has_value();

    simulation["SimulationType"] = "MolecularDynamics";
    simulation["PrintEvery"] = 1000;
    if (!dynamics.runs.empty())
    {
      // the last run is the production, those before it the equilibration
      simulation["NumberOfProductionCycles"] = dynamics.runs.back();
      simulation["NumberOfEquilibrationCycles"] =
          std::accumulate(dynamics.runs.begin(), dynamics.runs.end() - 1, std::size_t{0});
    }
    system["TimeStep"] = dynamics.timeStep.value_or(metal ? 0.001 : 1.0) * picoseconds;

    if (dynamics.integrator == "nve")
    {
      system["Ensemble"] = langevin ? "NVT" : "NVE";
    }
    else
    {
      system["Ensemble"] = dynamics.integrator == "nvt" ? "NVT" : dynamics.flexibleCell ? "NPTPR" : "NPT";
      if (dynamics.temperature) system["ExternalTemperature"] = *dynamics.temperature;
      if (dynamics.temperatureDamping)
        system["TimeScaleParameterThermostat"] = *dynamics.temperatureDamping * picoseconds;
      if (dynamics.thermostatChain) system["ThermostatChainLength"] = *dynamics.thermostatChain;
      if (dynamics.integrator == "npt")
      {
        if (dynamics.pressure) system["ExternalPressure"] = *dynamics.pressure * pascal;
        if (dynamics.pressureDamping) system["TimeScaleParameterBarostat"] = *dynamics.pressureDamping * picoseconds;
        if (dynamics.barostatChain) system["BarostatChainLength"] = *dynamics.barostatChain;
        if (dynamics.flexibleCell)
        {
          result.warnings.push_back(
              "[LAMMPS reader] fix npt with aniso, tri or separate x, y, z is converted to the flexible cell of "
              "'NPTPR': check the cell degrees of freedom it should have");
        }
      }
      if (langevin)
      {
        result.warnings.push_back(std::format(
            "[LAMMPS reader] fix langevin next to fix {}: the Langevin thermostat replaces the Nose-Hoover chain "
            "of the particles",
            dynamics.integrator));
      }
    }
    if (langevin)
    {
      // the 'damp' of fix langevin is the time over which a velocity is forgotten: the inverse of the friction
      system["Thermostat"] = "Langevin";
      system["ExternalTemperature"] = *dynamics.langevinTemperature;
      if (*dynamics.langevinDamping > 0.0)
        system["LangevinFriction"] = 1.0 / (*dynamics.langevinDamping * picoseconds);
    }
  }
  else if (dynamics.langevinTemperature.has_value())
  {
    result.warnings.push_back(
        "[LAMMPS reader] fix langevin without a fix that integrates (nve, nvt, npt): no dynamics is converted");
  }
  simulation["Systems"] = nlohmann::json::array({system});
  simulation["Components"] = nlohmann::json::array();
  for (const ReadComponent &component : result.components)
  {
    nlohmann::json entry;
    entry["Name"] = component.name;
    entry["TranslationProbability"] = 1.0;
    entry["RotationProbability"] = 1.0;
    if (component.definition.contains("Bonds")) entry["ReinsertionProbability"] = 0.5;
    // the molecules are those of the data file, read from the configuration: none are created next to them
    entry["CreateNumberOfMolecules"] = 0;
    if (script.shake)
    {
      const bool constrained = component.definition.contains("Bonds") &&
                               std::ranges::any_of(component.definition["Bonds"],
                                                   [](const nlohmann::json &bond) { return bond[1] == "FIXED"; });
      if (constrained && script.shake->tolerance > 0.0) entry["ConstraintTolerance"] = script.shake->tolerance;
      if (constrained && script.shake->iterations > 0)
        entry["ConstraintMaximumIterations"] = script.shake->iterations;
    }
    simulation["Components"].push_back(entry);
  }
  result.simulation = simulation;

  return result;
}

void writeRaspaInput(const ReadResult &result, const std::filesystem::path &directory)
{
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory / "force_field.json");
    file << result.forceField.dump(2) << '\n';
  }
  for (const ReadComponent &component : result.components)
  {
    std::ofstream file(directory / (component.name + ".json"));
    file << component.definition.dump(2) << '\n';
  }
  {
    std::ofstream file(directory / "simulation.json");
    file << result.simulation.dump(2) << '\n';
  }
  {
    std::ofstream file(directory / "initial_configuration.json");
    file << result.configuration.dump() << '\n';
  }
  if (!result.warnings.empty())
  {
    std::ofstream file(directory / "conversion_warnings.txt");
    for (const std::string &warning : result.warnings) file << warning << '\n';
  }
}
}  // namespace LAMMPS
