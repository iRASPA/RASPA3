module;

module janaf;

import std;

import units;

// Embedded data from the NIST-JANAF Thermochemical Tables, Fourth Edition
// (M.W. Chase, J. Phys. Chem. Ref. Data, Monograph 9, 1998; https://janaf.nist.gov).
// See janaf.ixx for the equations used (appendix A.2 of the PhD thesis of A. Rahbari).

namespace JANAF
{
namespace
{
struct GibbsEnergyFunctionEntry
{
  double temperature;          // [K]
  double gibbsEnergyFunction;  // -[G°(T) - H°(298.15 K)]/T  [J mol⁻¹ K⁻¹]
};

struct AtomCount
{
  std::string_view element;
  std::size_t count;
};

struct MoleculeData
{
  std::string_view name;                                  // canonical name (chemical formula)
  std::vector<std::string_view> aliases;                  // common names
  std::vector<AtomCount> composition;                     // atom counts for Eq. A78
  double enthalpyOfFormation0K;                           // ΔfH°(0 K) [kJ/mol]
  double enthalpyShift;                                   // H°(0 K) - H°(298.15 K) [kJ/mol]
  std::vector<GibbsEnergyFunctionEntry> gibbsEnergyFunction;
};

// ΔfH°(0 K) of the gaseous atoms [kJ/mol], used for the atomization energy (Eq. A78).

using namespace std::literals;

const std::map<std::string_view, double> atomEnthalpyOfFormation0K = {{"H"sv, 216.035}, {"C"sv, 711.185}, {"N"sv, 470.82}, {"O"sv, 246.79}, {"S"sv, 274.735}, {"Cl"sv, 119.621}};


// JANAF table codes:
//   H2    H-050
//   N2    N-023
//   O2    O-029
//   Cl2   Cl-073
//   H2O   H-064
//   NH3   H-083
//   H2S   H-080
//   CO    C-093
//   CO2   C-095
//   CH4   C-067
//   C2H2  C-127
//   C2H4  C-128
//   NO    N-005
//   NO2   N-007
//   N2O   N-026
//   SO2   O-034
//   HCl   Cl-026
const std::vector<MoleculeData> moleculeDatabase = {
    MoleculeData{
        "H2"sv,
        {"hydrogen"sv},
        {{"H"sv, 2}},
        0,
        -8.467,
        {{100, 155.408}, {200, 133.284}, {250, 131.152}, {298.15, 130.68}, {300, 130.68},
         {350, 131.032}, {400, 131.817}, {450, 132.834}, {500, 133.973}, {600, 136.392},
         {700, 138.822}, {800, 141.171}, {900, 143.411}, {1000, 145.536}, {1100, 147.549},
         {1200, 149.459}, {1300, 151.274}, {1400, 153.003}, {1500, 154.652}, {1600, 156.231},
         {1700, 157.743}, {1800, 159.197}, {1900, 160.595}, {2000, 161.943}, {2100, 163.244},
         {2200, 164.501}, {2300, 165.719}, {2400, 166.899}, {2500, 168.044}, {2600, 169.155},
         {2700, 170.236}, {2800, 171.288}, {2900, 172.313}, {3000, 173.311}}},
    MoleculeData{
        "N2"sv,
        {"nitrogen"sv},
        {{"N"sv, 2}},
        0,
        -8.67,
        {{100, 217.49}, {200, 194.272}, {250, 192.088}, {298.15, 191.609}, {300, 191.61},
         {350, 191.964}, {400, 192.753}, {450, 193.774}, {500, 194.917}, {600, 197.353},
         {700, 199.813}, {800, 202.209}, {900, 204.51}, {1000, 206.708}, {1100, 208.804},
         {1200, 210.802}, {1300, 212.71}, {1400, 214.533}, {1500, 216.277}, {1600, 217.948},
         {1700, 219.552}, {1800, 221.094}, {1900, 222.577}, {2000, 224.006}, {2100, 225.385},
         {2200, 226.717}, {2300, 228.004}, {2400, 229.25}, {2500, 230.458}, {2600, 231.629},
         {2700, 232.765}, {2800, 233.869}, {2900, 234.942}, {3000, 235.986}}},
    MoleculeData{
        "O2"sv,
        {"oxygen"sv},
        {{"O"sv, 2}},
        0,
        -8.683,
        {{100, 231.094}, {200, 207.823}, {250, 205.63}, {298.15, 205.147}, {300, 205.148},
         {350, 205.506}, {400, 206.308}, {450, 207.35}, {500, 208.524}, {600, 211.044},
         {700, 213.611}, {800, 216.126}, {900, 218.552}, {1000, 220.875}, {1100, 223.093},
         {1200, 225.209}, {1300, 227.229}, {1400, 229.158}, {1500, 231.002}, {1600, 232.768},
         {1700, 234.462}, {1800, 236.089}, {1900, 237.653}, {2000, 239.16}, {2100, 240.613},
         {2200, 242.017}, {2300, 243.374}, {2400, 244.687}, {2500, 245.959}, {2600, 247.194},
         {2700, 248.393}, {2800, 249.558}, {2900, 250.691}, {3000, 251.795}}},
    MoleculeData{
        "Cl2"sv,
        {"chlorine"sv},
        {{"Cl"sv, 2}},
        0,
        -9.181,
        {{100, 251.696}, {200, 226.12}, {250, 223.633}, {298.15, 223.079}, {300, 223.08},
         {350, 223.496}, {400, 224.431}, {450, 225.648}, {500, 227.02}, {600, 229.956},
         {700, 232.926}, {800, 235.814}, {900, 238.577}, {1000, 241.203}, {1100, 243.692},
         {1200, 246.051}, {1300, 248.289}, {1400, 250.415}, {1500, 252.438}, {1600, 254.366},
         {1700, 256.207}, {1800, 257.969}, {1900, 259.657}, {2000, 261.277}, {2100, 262.834},
         {2200, 264.333}, {2300, 265.778}, {2400, 267.173}, {2500, 268.522}, {2600, 269.827},
         {2700, 271.091}, {2800, 272.317}, {2900, 273.508}, {3000, 274.665}}},
    MoleculeData{
        "H2O"sv,
        {"water"sv},
        {{"H"sv, 2}, {"O"sv, 1}},
        -238.921,
        -9.904,
        {{100, 218.534}, {200, 191.896}, {298.15, 188.834}, {300, 188.835}, {400, 190.159},
         {500, 192.685}, {600, 195.55}, {700, 198.465}, {800, 201.322}, {900, 204.084},
         {1000, 206.738}, {1100, 209.285}, {1200, 211.73}, {1300, 214.08}, {1400, 216.341},
         {1500, 218.52}, {1600, 220.623}, {1700, 222.655}, {1800, 224.621}, {1900, 226.526},
         {2000, 228.374}, {2100, 230.167}, {2200, 231.909}, {2300, 233.604}, {2400, 235.253},
         {2500, 236.86}, {2600, 238.425}, {2700, 239.952}, {2800, 241.443}, {2900, 242.899},
         {3000, 244.321}}},
    MoleculeData{
        "NH3"sv,
        {"ammonia"sv},
        {{"N"sv, 1}, {"H"sv, 3}},
        -38.907,
        -10.045,
        {{100, 223.211}, {200, 195.962}, {298.15, 192.774}, {300, 192.775}, {400, 194.209},
         {500, 197.021}, {600, 200.302}, {700, 203.727}, {800, 207.16}, {900, 210.543},
         {1000, 213.849}, {1100, 217.069}, {1200, 220.197}, {1300, 223.236}, {1400, 226.187},
         {1500, 229.054}, {1600, 231.84}, {1700, 234.549}, {1800, 237.184}, {1900, 239.748},
         {2000, 242.244}, {2100, 244.677}, {2200, 247.048}, {2300, 249.36}, {2400, 251.616},
         {2500, 253.818}, {2600, 255.969}, {2700, 258.07}, {2800, 260.124}, {2900, 262.132},
         {3000, 264.097}}},
    MoleculeData{
        "H2S"sv,
        {"hydrogen sulfide"sv},
        {{"H"sv, 2}, {"S"sv, 1}},
        -17.584,
        -9.962,
        {{100, 235.33}, {200, 208.586}, {298.15, 205.757}, {300, 205.758}, {400, 207.115},
         {500, 209.726}, {600, 212.713}, {700, 215.777}, {800, 218.804}, {900, 221.75},
         {1000, 224.599}, {1100, 227.346}, {1200, 229.993}, {1300, 232.544}, {1400, 235.003},
         {1500, 237.375}, {1600, 239.666}, {1700, 241.879}, {1800, 244.019}, {1900, 246.091},
         {2000, 248.098}, {2100, 250.044}, {2200, 251.932}, {2300, 253.766}, {2400, 255.549},
         {2500, 257.282}, {2600, 258.97}, {2700, 260.613}, {2800, 262.214}, {2900, 263.776},
         {3000, 265.3}}},
    MoleculeData{
        "CO"sv,
        {"carbon monoxide"sv},
        {{"C"sv, 1}, {"O"sv, 1}},
        -113.805,
        -8.671,
        {{100, 223.539}, {200, 200.317}, {298.15, 197.653}, {300, 197.653}, {400, 198.798},
         {500, 200.968}, {600, 203.415}, {700, 205.89}, {800, 208.305}, {900, 210.628},
         {1000, 212.848}, {1100, 214.967}, {1200, 216.988}, {1300, 218.917}, {1400, 220.761},
         {1500, 222.526}, {1600, 224.216}, {1700, 225.839}, {1800, 227.398}, {1900, 228.897},
         {2000, 230.342}, {2100, 231.736}, {2200, 233.081}, {2300, 234.382}, {2400, 235.641},
         {2500, 236.86}, {2600, 238.043}, {2700, 239.19}, {2800, 240.304}, {2900, 241.387},
         {3000, 242.441}}},
    MoleculeData{
        "CO2"sv,
        {"carbon dioxide"sv},
        {{"C"sv, 1}, {"O"sv, 2}},
        -393.151,
        -9.364,
        {{100, 243.568}, {200, 217.046}, {298.15, 213.795}, {300, 213.795}, {400, 215.307},
         {500, 218.29}, {600, 221.772}, {700, 225.388}, {800, 228.986}, {900, 232.5},
         {1000, 235.901}, {1100, 239.178}, {1200, 242.329}, {1300, 245.356}, {1400, 248.265},
         {1500, 251.062}, {1600, 253.753}, {1700, 256.343}, {1800, 258.84}, {1900, 261.248},
         {2000, 263.574}, {2100, 265.822}, {2200, 267.996}, {2300, 270.102}, {2400, 272.144},
         {2500, 274.124}, {2600, 276.046}, {2700, 277.914}, {2800, 279.73}, {2900, 281.497},
         {3000, 283.218}}},
    MoleculeData{
        "CH4"sv,
        {"methane"sv},
        {{"C"sv, 1}, {"H"sv, 4}},
        -66.911,
        -10.024,
        {{100, 216.485}, {200, 189.418}, {250, 186.829}, {298.15, 186.251}, {300, 186.252},
         {350, 186.694}, {400, 187.704}, {450, 189.053}, {500, 190.614}, {600, 194.103},
         {700, 197.84}, {800, 201.675}, {900, 205.532}, {1000, 209.37}, {1100, 213.162},
         {1200, 216.895}, {1300, 220.558}, {1400, 224.148}, {1500, 227.66}, {1600, 231.095},
         {1700, 234.45}, {1800, 237.728}, {1900, 240.93}, {2000, 244.057}, {2100, 247.11},
         {2200, 250.093}, {2300, 253.007}, {2400, 255.854}, {2500, 258.638}, {2600, 261.359},
         {2700, 264.02}, {2800, 266.623}, {2900, 269.171}, {3000, 271.664}}},
    MoleculeData{
        "C2H2"sv,
        {"ethyne"sv, "acetylene"sv},
        {{"C"sv, 2}, {"H"sv, 2}},
        227.288,
        -10.012,
        {{100, 234.338}, {200, 204.72}, {298.15, 200.958}, {300, 200.959}, {400, 202.774},
         {500, 206.393}, {600, 210.64}, {700, 215.064}, {800, 219.476}, {900, 223.794},
         {1000, 227.984}, {1100, 232.034}, {1200, 235.941}, {1300, 239.709}, {1400, 243.344},
         {1500, 246.853}, {1600, 250.242}, {1700, 253.518}, {1800, 256.688}, {1900, 259.758},
         {2000, 262.733}, {2100, 265.62}, {2200, 268.422}, {2300, 271.145}, {2400, 273.793},
         {2500, 276.369}, {2600, 278.878}, {2700, 281.322}, {2800, 283.706}, {2900, 286.031},
         {3000, 288.301}}},
    MoleculeData{
        "C2H4"sv,
        {"ethene"sv, "ethylene"sv},
        {{"C"sv, 2}, {"H"sv, 4}},
        60.986,
        -10.518,
        {{100, 252.466}, {200, 222.975}, {250, 220.011}, {298.15, 219.33}, {300, 219.331},
         {350, 219.873}, {400, 221.138}, {450, 222.858}, {500, 224.879}, {600, 229.456},
         {700, 234.408}, {800, 239.511}, {900, 244.644}, {1000, 249.742}, {1100, 254.768},
         {1200, 259.698}, {1300, 264.522}, {1400, 269.233}, {1500, 273.827}, {1600, 278.306},
         {1700, 282.67}, {1800, 286.922}, {1900, 291.064}, {2000, 295.101}, {2100, 299.035},
         {2200, 302.87}, {2300, 306.61}, {2400, 310.258}, {2500, 313.818}, {2600, 317.294},
         {2700, 320.689}, {2800, 324.006}, {2900, 327.248}, {3000, 330.418}}},
    MoleculeData{
        "NO"sv,
        {"nitric oxide"sv, "nitrogen oxide"sv},
        {{"N"sv, 1}, {"O"sv, 1}},
        89.775,
        -9.192,
        {{100, 237.757}, {200, 213.501}, {250, 211.251}, {298.15, 210.758}, {300, 210.759},
         {350, 211.122}, {400, 211.929}, {450, 212.974}, {500, 214.145}, {600, 216.646},
         {700, 219.179}, {800, 221.652}, {900, 224.031}, {1000, 226.307}, {1100, 228.478},
         {1200, 230.549}, {1300, 232.525}, {1400, 234.412}, {1500, 236.217}, {1600, 237.945},
         {1700, 239.603}, {1800, 241.195}, {1900, 242.725}, {2000, 244.199}, {2100, 245.619},
         {2200, 246.99}, {2300, 248.315}, {2400, 249.596}, {2500, 250.837}, {2600, 252.039},
         {2700, 253.205}, {2800, 254.338}, {2900, 255.438}, {3000, 256.508}}},
    MoleculeData{
        "NO2"sv,
        {"nitrogen dioxide"sv},
        {{"N"sv, 1}, {"O"sv, 2}},
        35.927,
        -10.186,
        {{100, 271.168}, {200, 243.325}, {250, 240.634}, {298.15, 240.034}, {300, 240.034},
         {350, 240.491}, {400, 241.524}, {450, 242.886}, {500, 244.44}, {600, 247.83},
         {700, 251.345}, {800, 254.84}, {900, 258.25}, {1000, 261.545}, {1100, 264.717},
         {1200, 267.761}, {1300, 270.683}, {1400, 273.485}, {1500, 276.175}, {1600, 278.759},
         {1700, 281.244}, {1800, 283.634}, {1900, 285.937}, {2000, 288.158}, {2100, 290.302},
         {2200, 292.373}, {2300, 294.377}, {2400, 296.316}, {2500, 298.196}, {2600, 300.018},
         {2700, 301.787}, {2800, 303.505}, {2900, 305.176}, {3000, 306.8}}},
    MoleculeData{
        "N2O"sv,
        {"nitrous oxide"sv},
        {{"N"sv, 2}, {"O"sv, 1}},
        85.481,
        -9.579,
        {{100, 250.829}, {200, 223.335}, {250, 220.58}, {298.15, 219.957}, {300, 219.958},
         {350, 220.437}, {400, 221.526}, {450, 222.967}, {500, 224.613}, {600, 228.204},
         {700, 231.924}, {800, 235.618}, {900, 239.217}, {1000, 242.694}, {1100, 246.038},
         {1200, 249.248}, {1300, 252.328}, {1400, 255.283}, {1500, 258.12}, {1600, 260.846},
         {1700, 263.467}, {1800, 265.991}, {1900, 268.423}, {2000, 270.769}, {2100, 273.035},
         {2200, 275.225}, {2300, 277.344}, {2400, 279.396}, {2500, 281.385}, {2600, 283.314},
         {2700, 285.188}, {2800, 287.008}, {2900, 288.778}, {3000, 290.5}}},
    MoleculeData{
        "SO2"sv,
        {"sulfur dioxide"sv},
        {{"S"sv, 1}, {"O"sv, 2}},
        -294.299,
        -10.552,
        {{100, 281.199}, {200, 251.714}, {298.15, 248.212}, {300, 248.213}, {400, 249.824},
         {500, 252.979}, {600, 256.641}, {700, 260.427}, {800, 264.178}, {900, 267.825},
         {1000, 271.339}, {1100, 274.71}, {1200, 277.937}, {1300, 281.026}, {1400, 283.983},
         {1500, 286.816}, {1600, 289.533}, {1700, 292.141}, {1800, 294.647}, {1900, 297.059},
         {2000, 299.383}, {2100, 301.624}, {2200, 303.787}, {2300, 305.878}, {2400, 307.902},
         {2500, 309.861}, {2600, 311.761}, {2700, 313.604}, {2800, 315.394}, {2900, 317.133},
         {3000, 318.825}}},
    MoleculeData{
        "HCl"sv,
        {"hydrogen chloride"sv, "hydrochloric acid"sv},
        {{"H"sv, 1}, {"Cl"sv, 1}},
        -92.127,
        -8.64,
        {{100, 212.797}, {200, 189.566}, {250, 187.381}, {298.15, 186.901}, {300, 186.902},
         {350, 187.256}, {400, 188.045}, {450, 189.064}, {500, 190.205}, {600, 192.629},
         {700, 195.068}, {800, 197.434}, {900, 199.699}, {1000, 201.857}, {1100, 203.91},
         {1200, 205.866}, {1300, 207.73}, {1400, 209.51}, {1500, 211.214}, {1600, 212.846},
         {1700, 214.413}, {1800, 215.92}, {1900, 217.371}, {2000, 218.769}, {2100, 220.119},
         {2200, 221.424}, {2300, 222.687}, {2400, 223.91}, {2500, 225.097}, {2600, 226.248},
         {2700, 227.366}, {2800, 228.453}, {2900, 229.511}, {3000, 230.54}}},
};

// case-insensitive matching, ignoring spaces, hyphens, and underscores
std::string normalizeName(std::string_view name)
{
  std::string normalized;
  normalized.reserve(name.size());
  for (char c : name)
  {
    if (c == ' ' || c == '-' || c == '_')
    {
      continue;
    }
    normalized.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  return normalized;
}

const MoleculeData *findMolecule(std::string_view moleculeName)
{
  const std::string normalized = normalizeName(moleculeName);
  for (const MoleculeData &molecule : moleculeDatabase)
  {
    if (normalizeName(molecule.name) == normalized)
    {
      return &molecule;
    }
    for (std::string_view alias : molecule.aliases)
    {
      if (normalizeName(alias) == normalized)
      {
        return &molecule;
      }
    }
  }
  return nullptr;
}

const MoleculeData &getMolecule(std::string_view moleculeName)
{
  const MoleculeData *molecule = findMolecule(moleculeName);
  if (!molecule)
  {
    throw std::runtime_error(
        std::format("JANAF: no thermochemical data available for molecule '{}'\n", moleculeName));
  }
  return *molecule;
}

// -[G°(T) - H°(0 K)]/T in [J mol⁻¹ K⁻¹], obtained by linear interpolation of the tabulated
// Gibbs energy function and shifting the reference enthalpy from 298.15 K to 0 K (Eq. A82)
double gibbsEnergyFunctionReferencedTo0K(const MoleculeData &molecule, double temperature)
{
  const std::vector<GibbsEnergyFunctionEntry> &table = molecule.gibbsEnergyFunction;
  if (temperature < table.front().temperature || temperature > table.back().temperature)
  {
    throw std::runtime_error(
        std::format("JANAF: temperature {} K for molecule '{}' outside the tabulated range [{}, {}] K\n", temperature,
                    molecule.name, table.front().temperature, table.back().temperature));
  }

  auto upper = std::ranges::lower_bound(table, temperature, {}, &GibbsEnergyFunctionEntry::temperature);
  if (upper == table.begin())
  {
    ++upper;
  }
  auto lower = std::prev(upper);
  const double fraction = (temperature - lower->temperature) / (upper->temperature - lower->temperature);
  const double gibbsEnergyFunction =
      lower->gibbsEnergyFunction + fraction * (upper->gibbsEnergyFunction - lower->gibbsEnergyFunction);

  return gibbsEnergyFunction + molecule.enthalpyShift * 1000.0 / temperature;
}
}  // namespace

bool contains(std::string_view moleculeName) { return findMolecule(moleculeName) != nullptr; }

std::vector<std::string_view> availableMolecules()
{
  std::vector<std::string_view> names;
  names.reserve(moleculeDatabase.size());
  for (const MoleculeData &molecule : moleculeDatabase)
  {
    names.push_back(molecule.name);
  }
  return names;
}

double logPartitionFunctionGroundStateReference(std::string_view moleculeName, double temperature)
{
  const MoleculeData &molecule = getMolecule(moleculeName);

  // Eq. A82: -[G°(T) - H°(0 K)]/T = R ln[(q₀/V)(k_B T / P°)] with P° = 1 bar
  constexpr double standardPressure = 1.0e5;  // [Pa]
  constexpr double cubicMetersToCubicAngstrom = 1.0e30;
  const double thermalVolume =
      Units::BoltzmannConstant * temperature / standardPressure * cubicMetersToCubicAngstrom;  // k_B T / P° in [Å³]

  return gibbsEnergyFunctionReferencedTo0K(molecule, temperature) / Units::MolarGasConstant - std::log(thermalVolume);
}

double atomizationEnergy(std::string_view moleculeName)
{
  const MoleculeData &molecule = getMolecule(moleculeName);

  // Eq. A78: D₀ = Σᵢ yᵢ ΔfH°ᵢ(0 K) - ΔfH°(0 K)
  double atomizationEnergy = -molecule.enthalpyOfFormation0K;
  for (const AtomCount &atom : molecule.composition)
  {
    atomizationEnergy += static_cast<double>(atom.count) * atomEnthalpyOfFormation0K.at(atom.element);
  }
  return atomizationEnergy;
}

double logPartitionFunction(std::string_view moleculeName, double temperature)
{
  // Eq. A74: q/V = (q₀/V) exp[D₀/(k_B T)], referencing the energy to the dissociated atoms
  return logPartitionFunctionGroundStateReference(moleculeName, temperature) +
         atomizationEnergy(moleculeName) * 1000.0 / (Units::MolarGasConstant * temperature);
}
}  // namespace JANAF
