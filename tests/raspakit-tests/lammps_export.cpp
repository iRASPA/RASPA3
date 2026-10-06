#include <gtest/gtest.h>

import std;

import atom;
import bend_potential;
import bond_potential;
import component;
import connectivity_table;
import double3;
import forcefield;
import framework;
import inversion_bend_potential;
import json;
import lammps_io;
import lammps_reader;
import lammps_styles;
import lammps_topology;
import molecule;
import simulationbox;
import system;
import torsion_potential;
import units;
import vdwparameters;

// Every LAMMPS functional form is re-implemented here, straight from the LAMMPS documentation, and evaluated
// on random geometries against RASPA's own calculateEnergy of the term that was mapped onto it.

namespace
{
constexpr double kDegrees = std::numbers::pi / 180.0;

std::vector<double> numbers(std::string_view text)
{
  std::vector<double> values;
  std::istringstream stream{std::string(text)};
  std::string word;
  while (stream >> word)
  {
    try
    {
      values.push_back(std::stod(word));
    }
    catch (...)
    {
    }
  }
  return values;
}

double kcal(double internalEnergy) { return internalEnergy * Units::EnergyToKCalPerMol; }

struct Geometry
{
  double3 a, b, c, d;
};

Geometry randomGeometry(std::mt19937 &rng)
{
  std::uniform_real_distribution<double> unit(-1.0, 1.0);
  auto point = [&] { return double3(unit(rng), unit(rng), unit(rng)); };
  // chain-like: keep consecutive atoms 1.2 - 1.8 apart so all angles are well defined
  double3 a = point();
  auto next = [&](const double3 &from)
  {
    double3 direction = point();
    direction = direction / std::sqrt(double3::dot(direction, direction));
    return from + direction * (1.2 + 0.6 * (unit(rng) + 1.0) * 0.5);
  };
  double3 b = next(a), c = next(b), d = next(c);
  return {a, b, c, d};
}

double distance(const double3 &a, const double3 &b)
{
  double3 r = a - b;
  return std::sqrt(double3::dot(r, r));
}

double angle(const double3 &a, const double3 &b, const double3 &c)
{
  double3 u = a - b, v = c - b;
  return std::acos(std::clamp(double3::dot(u, v) / std::sqrt(double3::dot(u, u) * double3::dot(v, v)), -1.0, 1.0));
}

// LAMMPS dihedral: angle between planes ABC and BCD, trans = 180 degrees
double dihedral(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  double3 b1 = b - a, b2 = c - b, b3 = d - c;
  double3 n1 = double3::cross(b1, b2), n2 = double3::cross(b2, b3);
  double3 m = double3::cross(n1, b2 / std::sqrt(double3::dot(b2, b2)));
  double x = double3::dot(n1, n2), y = double3::dot(m, n2);
  return std::atan2(y, x);
}

// LAMMPS improper umbrella: omega is the angle between the I-L axis and the I-J-K plane
double umbrellaAngle(const double3 &i, const double3 &j, const double3 &k, const double3 &l)
{
  double3 axis = l - i;
  double3 normal = double3::cross(j - i, k - i);
  double s = double3::dot(axis, normal) / std::sqrt(double3::dot(axis, axis) * double3::dot(normal, normal));
  return std::asin(std::clamp(s, -1.0, 1.0));
}

double lammpsBond(std::string_view style, std::span<const double> c, double r)
{
  if (style == "harmonic") return c[0] * (r - c[1]) * (r - c[1]);
  if (style == "class2")
  {
    double dr = r - c[0];
    return c[1] * dr * dr + c[2] * dr * dr * dr + c[3] * dr * dr * dr * dr;
  }
  if (style == "morse") return c[0] * std::pow(1.0 - std::exp(-c[1] * (r - c[2])), 2.0);
  if (style == "zero") return 0.0;
  throw std::runtime_error(std::format("no evaluator for bond style {}", style));
}

double lammpsAngle(std::string_view style, std::span<const double> c, double theta)
{
  if (style == "harmonic") return c[0] * std::pow(theta - c[1] * kDegrees, 2.0);
  if (style == "cosine/squared") return c[0] * std::pow(std::cos(theta) - std::cos(c[1] * kDegrees), 2.0);
  if (style == "quartic")
  {
    double dt = theta - c[0] * kDegrees;
    return c[1] * dt * dt + c[2] * dt * dt * dt + c[3] * dt * dt * dt * dt;
  }
  if (style == "charmm") return c[0] * std::pow(theta - c[1] * kDegrees, 2.0);  // plus K_ub (r13 - r_ub)^2, tested separately
  if (style == "zero") return 0.0;
  throw std::runtime_error(std::format("no evaluator for angle style {}", style));
}

double lammpsDihedral(std::string_view style, std::span<const double> c, double phi)
{
  if (style == "nharmonic")
  {
    std::size_t n = static_cast<std::size_t>(std::lround(c[0]));
    double energy = 0.0, power = 1.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      energy += c[1 + i] * power;
      power *= std::cos(phi);
    }
    return energy;
  }
  if (style == "quadratic")
  {
    double dphi = phi - c[1] * kDegrees;
    // periodic image of the difference, as RASPA's harmonic torsion does
    dphi = std::remainder(dphi, 2.0 * std::numbers::pi);
    return c[0] * dphi * dphi;
  }
  if (style == "charmm") return c[0] * (1.0 + std::cos(c[1] * phi - c[2] * kDegrees));
  if (style == "zero") return 0.0;
  throw std::runtime_error(std::format("no evaluator for dihedral style {}", style));
}

double lammpsUmbrella(std::span<const double> c, double omega)
{
  double omega0 = c[1] * kDegrees;
  if (std::abs(omega0) < 1e-12) return c[0] * (1.0 - std::cos(omega));
  return 0.5 * c[0] / std::pow(std::sin(omega0), 2.0) * std::pow(std::cos(omega) - std::cos(omega0), 2.0);
}

ForceField chainForceField()
{
  return ForceField({{"C", false, 12.0, 0.0, 0.0, 6, true}, {"O", false, 16.0, -0.4, 0.0, 8, true}},
                    {{120.0, 3.4}, {80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false,
                    false, true);
}

class TemporaryDirectory
{
 public:
  TemporaryDirectory()
  {
    static std::size_t counter{};
    path = std::filesystem::temp_directory_path() /
           std::format("raspa3-lammps-test-{}-{}", std::random_device{}(), counter++);
    std::filesystem::create_directories(path);
  }
  ~TemporaryDirectory()
  {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
  std::filesystem::path path;
};

Component readComponent(const ForceField &forceField, const nlohmann::json &definition, const std::string &name)
{
  TemporaryDirectory directory;
  std::filesystem::path file = directory.path / (name + ".json");
  {
    std::ofstream stream(file);
    stream << definition.dump(2);
  }
  std::filesystem::path stem = file;
  stem.replace_extension();
  return Component(Component::Type::Adsorbate, 0, forceField, name, stem.string(), 5, 21);
}

nlohmann::json chainDefinition()
{
  return nlohmann::json::parse(R"({
    "CriticalTemperature": 500.0, "CriticalPressure": 4.0e6, "AcentricFactor": 0.3,
    "PseudoAtoms": [
      ["C", [0.0, 0.0, 0.0]], ["C", [1.5, 0.0, 0.0]], ["O", [2.2, 1.2, 0.0]], ["C", [3.7, 1.2, 0.0]], ["C", [4.4, 2.4, 0.3]]
    ],
    "Connectivity": [[0, 1], [1, 2], [2, 3], [3, 4]],
    "Intra14ChargeChargeScalingValue": 0.5,
    "Bonds": [
      [[0, 1], "HARMONIC", [301931.64, 1.54]],
      [[1, 2], "HARMONIC", [250000.0, 1.41]],
      [[2, 3], "HARMONIC", [250000.0, 1.41]],
      [[3, 4], "HARMONIC", [301931.64, 1.54]]
    ],
    "Bends": [
      [[0, 1, 2], "HARMONIC", [62500.0, 112.0]],
      [[1, 2, 3], "HARMONIC", [60400.0, 112.0]],
      [[2, 3, 4], "HARMONIC", [62500.0, 112.0]]
    ],
    "Torsions": [
      [[0, 1, 2, 3], "TRAPPE_EXTENDED", [1078.16, 355.03, 68.19, 791.32, 0.0]],
      [[1, 2, 3, 4], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
    ]
  })");
}
}  // namespace

// ---------------------------------------------------------------------------------------------------
// Functional forms
// ---------------------------------------------------------------------------------------------------

TEST(lammps_export, bond_styles_match_raspa_energies)
{
  std::mt19937 rng(7);
  struct Case
  {
    BondType type;
    std::vector<double> parameters;
    std::string expectedStyle;
  };
  std::vector<Case> cases{{BondType::Harmonic, {300000.0, 1.52}, "harmonic"},
                          {BondType::Quartic, {300000.0, 1.52, -8000.0, 2000.0}, "class2"},
                          {BondType::CFF_Quartic, {300000.0, 1.52, -8000.0, 2000.0}, "class2"}};
  for (const Case &testCase : cases)
  {
    BondPotential bond({0, 1}, testCase.type, testCase.parameters);
    LAMMPS::Term term = LAMMPS::bondTerm(bond);
    ASSERT_EQ(term.style, testCase.expectedStyle);
    EXPECT_TRUE(term.exact);
    std::vector<double> c = numbers(term.coefficients);
    for (int sample = 0; sample < 20; ++sample)
    {
      Geometry g = randomGeometry(rng);
      double expected = kcal(bond.calculateEnergy(g.a, g.b));
      double actual = lammpsBond(term.style, c, distance(g.a, g.b));
      EXPECT_NEAR(actual, expected, 1e-6 * std::max(1.0, std::abs(expected))) << term.style << " " << term.coefficients;
    }
  }
}

TEST(lammps_export, bend_styles_match_raspa_energies)
{
  std::mt19937 rng(11);
  struct Case
  {
    BendType type;
    std::vector<double> parameters;
    std::string expectedStyle;
  };
  std::vector<Case> cases{{BendType::Harmonic, {62500.0, 114.0}, "harmonic"},
                          {BendType::HarmonicCosine, {40000.0, 109.5}, "cosine/squared"},
                          {BendType::Quartic, {62500.0, 114.0, -3000.0, 500.0}, "quartic"},
                          {BendType::CFF_Quartic, {62500.0, 114.0, -3000.0, 500.0}, "quartic"}};
  for (const Case &testCase : cases)
  {
    BendPotential bend({0, 1, 2}, testCase.type, testCase.parameters);
    LAMMPS::Term term = LAMMPS::bendTerm(bend);
    ASSERT_EQ(term.style, testCase.expectedStyle);
    EXPECT_TRUE(term.exact);
    std::vector<double> c = numbers(term.coefficients);
    for (int sample = 0; sample < 20; ++sample)
    {
      Geometry g = randomGeometry(rng);
      double expected = kcal(bend.calculateEnergy(g.a, g.b, g.c, std::nullopt));
      double actual = lammpsAngle(term.style, c, angle(g.a, g.b, g.c));
      EXPECT_NEAR(actual, expected, 1e-6 * std::max(1.0, std::abs(expected))) << term.style << " " << term.coefficients;
    }
  }
}

TEST(lammps_export, torsion_styles_match_raspa_energies)
{
  std::mt19937 rng(13);
  struct Case
  {
    TorsionType type;
    std::vector<double> parameters;
    std::string expectedStyle;
  };
  std::vector<Case> cases{
      {TorsionType::TraPPE_Extended, {1078.16, 355.03, 68.19, 791.32, -51.27}, "nharmonic"},
      {TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32}, "nharmonic"},
      {TorsionType::RyckaertBellemans, {1116.0, 1462.0, -1578.0, -368.0, 3156.0, -3788.0}, "nharmonic"},
      {TorsionType::OPLS, {0.0, 355.03, -68.19, 791.32}, "nharmonic"},
      {TorsionType::FourierSeries, {100.0, 200.0, 300.0, 400.0, 500.0, 600.0}, "nharmonic"},
      {TorsionType::HarmonicCosine, {5000.0, 60.0}, "nharmonic"},
      {TorsionType::ThreeCosine, {100.0, 200.0, 300.0}, "nharmonic"},
      {TorsionType::CVFF, {800.0, 3.0, 0.0}, "charmm"},
      {TorsionType::CVFF, {800.0, 2.0, 180.0}, "charmm"},
      {TorsionType::Harmonic, {5000.0, 180.0}, "quadratic"},
      {TorsionType::Polynomial, {10.0, 20.0, 30.0, 40.0, 50.0, 60.0}, "nharmonic"},
  };
  for (const Case &testCase : cases)
  {
    TorsionPotential torsion({0, 1, 2, 3}, testCase.type, testCase.parameters);
    LAMMPS::Term term = LAMMPS::torsionTerm(torsion);
    ASSERT_EQ(term.style, testCase.expectedStyle) << std::to_underlying(testCase.type);
    EXPECT_TRUE(term.exact) << term.note;
    std::vector<double> c = numbers(term.coefficients);
    for (int sample = 0; sample < 30; ++sample)
    {
      Geometry g = randomGeometry(rng);
      double expected = kcal(torsion.calculateEnergy(g.a, g.b, g.c, g.d));
      double actual = lammpsDihedral(term.style, c, dihedral(g.a, g.b, g.c, g.d));
      EXPECT_NEAR(actual, expected, 1e-6 * std::max(1.0, std::abs(expected)))
          << std::to_underlying(testCase.type) << ": " << term.style << " " << term.coefficients;
    }
  }
}

TEST(lammps_export, inversion_bends_map_onto_improper_umbrella)
{
  std::mt19937 rng(17);
  struct Case
  {
    InversionBendType type;
    std::vector<double> parameters;
  };
  std::vector<Case> cases{{InversionBendType::Planar, {2000.0}}, {InversionBendType::HarmonicCosine, {3000.0, 30.0}}};
  for (const Case &testCase : cases)
  {
    InversionBendPotential inversion({0, 1, 2, 3}, testCase.type, testCase.parameters);
    LAMMPS::Term term = LAMMPS::inversionBendTerm(inversion);
    ASSERT_EQ(term.style, "umbrella");
    EXPECT_TRUE(term.exact);
    std::vector<double> c = numbers(term.coefficients);
    for (int sample = 0; sample < 30; ++sample)
    {
      Geometry g = randomGeometry(rng);
      std::array<double3, 4> raspaOrder{g.a, g.b, g.c, g.d};
      double expected = kcal(inversion.calculateEnergy(g.a, g.b, g.c, g.d));
      // LAMMPS I J K L = RASPA identifiers permuted by term.atomOrder
      double omega = umbrellaAngle(raspaOrder[term.atomOrder[0]], raspaOrder[term.atomOrder[1]],
                                   raspaOrder[term.atomOrder[2]], raspaOrder[term.atomOrder[3]]);
      double actual = lammpsUmbrella(c, omega);
      EXPECT_NEAR(actual, expected, 1e-6 * std::max(1.0, std::abs(expected))) << term.coefficients;
    }
  }
}

TEST(lammps_export, polynomial_torsion_equals_cosine_series_and_gradient)
{
  std::mt19937 rng(19);
  TorsionPotential trappe({0, 1, 2, 3}, TorsionType::TraPPE_Extended, {1078.16, 355.03, 68.19, 791.32, -51.27});
  std::optional<std::vector<double>> fourier = LAMMPS::torsionFourierCoefficients(trappe);
  ASSERT_TRUE(fourier.has_value());
  std::vector<double> polynomial = LAMMPS::chebyshevToPolynomial(*fourier);
  ASSERT_EQ(polynomial.size(), 5u);
  for (double &value : polynomial) value *= Units::EnergyToKelvin;
  TorsionPotential poly({0, 1, 2, 3}, TorsionType::Polynomial, polynomial);

  for (int sample = 0; sample < 30; ++sample)
  {
    Geometry g = randomGeometry(rng);
    EXPECT_NEAR(poly.calculateEnergy(g.a, g.b, g.c, g.d), trappe.calculateEnergy(g.a, g.b, g.c, g.d), 1e-9);

    auto [energy, gradient, strain] = poly.potentialEnergyGradientStrain(g.a, g.b, g.c, g.d);
    EXPECT_NEAR(energy, trappe.calculateEnergy(g.a, g.b, g.c, g.d), 1e-9);
    // finite-difference check of the gradient of atom A and D
    const double h = 1e-5;
    for (std::size_t k = 0; k < 3; ++k)
    {
      double3 plus = g.a, minus = g.a;
      plus[k] += h;
      minus[k] -= h;
      double numerical = (poly.calculateEnergy(plus, g.b, g.c, g.d) - poly.calculateEnergy(minus, g.b, g.c, g.d)) / (2.0 * h);
      EXPECT_NEAR(gradient[0][k], numerical, 1e-5 * std::max(1.0, std::abs(numerical)));
      double3 dplus = g.d, dminus = g.d;
      dplus[k] += h;
      dminus[k] -= h;
      numerical = (poly.calculateEnergy(g.a, g.b, g.c, dplus) - poly.calculateEnergy(g.a, g.b, g.c, dminus)) / (2.0 * h);
      EXPECT_NEAR(gradient[3][k], numerical, 1e-5 * std::max(1.0, std::abs(numerical)));
    }
  }
}

TEST(lammps_export, pair_lennard_jones_round_trip)
{
  ForceField forceField = chainForceField();
  LAMMPS::PairTerm term = LAMMPS::pairTerm(forceField, 0, 1);
  EXPECT_EQ(term.style, "lj/cut");
  EXPECT_TRUE(term.exact);
  std::vector<double> c = numbers(term.coefficients);
  ASSERT_EQ(c.size(), 2u);
  // Lorentz-Berthelot: epsilon = sqrt(120 * 80) K, sigma = 3.2
  EXPECT_NEAR(c[0], std::sqrt(120.0 * 80.0) * Units::KelvinToEnergy * Units::EnergyToKCalPerMol, 1e-9);
  EXPECT_NEAR(c[1], 3.2, 1e-12);

  std::optional<std::array<double, 2>> back = LAMMPS::lennardJonesFromLammps("lj/cut", c);
  ASSERT_TRUE(back.has_value());
  EXPECT_NEAR((*back)[0], std::sqrt(120.0 * 80.0), 1e-6);
  EXPECT_NEAR((*back)[1], 3.2, 1e-12);
}

TEST(lammps_export, reader_inverses_recover_raspa_parameters)
{
  BondPotential bond({0, 1}, BondType::Harmonic, {301931.64, 1.54});
  LAMMPS::Term bondTerm = LAMMPS::bondTerm(bond);
  std::optional<LAMMPS::RaspaTerm> bondBack = LAMMPS::bondFromLammps(bondTerm.style, numbers(bondTerm.coefficients));
  ASSERT_TRUE(bondBack.has_value());
  EXPECT_EQ(bondBack->type, "HARMONIC");
  EXPECT_NEAR(bondBack->parameters[0], 301931.64, 1e-3);
  EXPECT_NEAR(bondBack->parameters[1], 1.54, 1e-9);

  BendPotential bend({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0});
  LAMMPS::Term bendTerm = LAMMPS::bendTerm(bend);
  std::optional<LAMMPS::RaspaTerm> bendBack = LAMMPS::bendFromLammps(bendTerm.style, numbers(bendTerm.coefficients));
  ASSERT_TRUE(bendBack.has_value());
  EXPECT_EQ(bendBack->type, "HARMONIC");
  EXPECT_NEAR(bendBack->parameters[0], 62500.0, 1e-3);
  EXPECT_NEAR(bendBack->parameters[1], 114.0, 1e-9);

  std::mt19937 rng(23);
  TorsionPotential torsion({0, 1, 2, 3}, TorsionType::TraPPE_Extended, {1078.16, 355.03, 68.19, 791.32, -51.27});
  LAMMPS::Term torsionTerm = LAMMPS::torsionTerm(torsion);
  std::optional<LAMMPS::RaspaTerm> torsionBack =
      LAMMPS::torsionFromLammps(torsionTerm.style, numbers(torsionTerm.coefficients));
  ASSERT_TRUE(torsionBack.has_value());
  EXPECT_EQ(torsionBack->type, "POLYNOMIAL");
  TorsionPotential recovered({0, 1, 2, 3}, TorsionType::Polynomial, torsionBack->parameters);
  for (int sample = 0; sample < 20; ++sample)
  {
    Geometry g = randomGeometry(rng);
    EXPECT_NEAR(recovered.calculateEnergy(g.a, g.b, g.c, g.d), torsion.calculateEnergy(g.a, g.b, g.c, g.d), 1e-6);
  }

  // the authors' custom 'trappe' style stores (1/2) K_n
  std::optional<LAMMPS::RaspaTerm> trappe = LAMMPS::torsionFromLammps("trappe", std::array{4.0, 2.0, 0.0, 0.0, 0.0});
  ASSERT_TRUE(trappe.has_value());
  EXPECT_EQ(trappe->type, "TRAPPE_EXTENDED");
  EXPECT_NEAR(trappe->parameters[0], 2.0 * Units::KCalPerMolToEnergy * Units::EnergyToKelvin, 1e-6);
  EXPECT_NEAR(trappe->parameters[1], 1.0 * Units::KCalPerMolToEnergy * Units::EnergyToKelvin, 1e-6);
}

// ---------------------------------------------------------------------------------------------------
// Topology, writer and reader on a small chain system
// ---------------------------------------------------------------------------------------------------

TEST(lammps_export, chain_system_export_and_read_back)
{
  ForceField forceField = chainForceField();
  Component chain = readComponent(forceField, chainDefinition(), "chain");
  ASSERT_FALSE(chain.rigid);

  const SimulationBox box(20.0, 20.0, 20.0);
  System system(forceField, box, false, 300.0, 1e5, 1.0, {}, {chain}, {}, {2}, 5);
  ASSERT_EQ(system.atomData.size(), 10u);

  // place the two molecules (component positions are centred); the second one straddles the boundary in x
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t m = 0; m < 2; ++m)
  {
    double3 offset = m == 0 ? double3(5.0, 5.0, 5.0) : double3(20.0, 12.0, 8.0);
    for (std::size_t i = 0; i < 5; ++i)
    {
      atoms[m * 5 + i].position = offset + chain.atoms[i].position;
    }
  }

  LAMMPS::ExportOptions options{};
  options.dataFile = "chain.data";
  LAMMPS::ExportFiles files =
      LAMMPS::exportSystem(system.components, system.atomData, system.atomDynamics, system.moleculeData,
                           system.simulationBox, system.forceField, system.numberOfIntegerMoleculesPerComponent,
                           system.framework, options);

  // header and deduplicated types
  EXPECT_NE(files.data.find("10 atoms\n"), std::string::npos);
  EXPECT_NE(files.data.find("8 bonds\n"), std::string::npos);
  EXPECT_NE(files.data.find("6 angles\n"), std::string::npos);
  EXPECT_NE(files.data.find("4 dihedrals\n"), std::string::npos);
  EXPECT_NE(files.data.find("2 bond types\n"), std::string::npos);
  EXPECT_NE(files.data.find("2 angle types\n"), std::string::npos);
  // the two torsions are the same cosine series -> one nharmonic type
  EXPECT_NE(files.data.find("1 dihedral types\n"), std::string::npos);
  EXPECT_NE(files.data.find("Bond Coeffs"), std::string::npos);
  EXPECT_NE(files.data.find("PairIJ Coeffs"), std::string::npos);
  EXPECT_NE(files.data.find("Atoms # full"), std::string::npos);

  // molecule ids start at 1 and the second molecule carries an image flag in x
  EXPECT_NE(files.data.find("\n1 1 1 "), std::string::npos);
  EXPECT_NE(files.data.find("\n6 2 1 "), std::string::npos);
  {
    std::istringstream stream(files.data);
    std::string line;
    bool inAtoms = false;
    int imageX = 0;
    while (std::getline(stream, line))
    {
      if (line.starts_with("Atoms")) { inAtoms = true; continue; }
      if (inAtoms && line.starts_with("Velocities")) break;
      if (!inAtoms || line.empty()) continue;
      std::vector<double> v = numbers(line);
      ASSERT_EQ(v.size(), 10u) << line;
      EXPECT_GE(v[4], 0.0);
      EXPECT_LT(v[4], 20.0);
      if (v[0] >= 6) imageX += static_cast<int>(v[7]) != 0;
    }
    EXPECT_GT(imageX, 0);
  }

  // input script
  EXPECT_NE(files.input.find("bond_style harmonic"), std::string::npos);
  EXPECT_NE(files.input.find("angle_style harmonic"), std::string::npos);
  EXPECT_NE(files.input.find("dihedral_style nharmonic"), std::string::npos);
  EXPECT_NE(files.input.find("special_bonds lj 0.0 0.0 0 coul 0.0 0.0 0.5"), std::string::npos);
  EXPECT_NE(files.input.find("read_data chain.data"), std::string::npos);
  EXPECT_NE(files.input.find("run 0"), std::string::npos);
  EXPECT_TRUE(files.table.empty());
  EXPECT_TRUE(files.pairList.empty());

  // read back
  TemporaryDirectory directory;
  std::filesystem::path dataPath = directory.path / "chain.data";
  std::filesystem::path inputPath = directory.path / "chain.in";
  std::ofstream(dataPath) << files.data;
  std::ofstream(inputPath) << files.input;

  LAMMPS::ReadResult read = LAMMPS::readDataFile(dataPath, inputPath);
  ASSERT_EQ(read.components.size(), 1u);
  EXPECT_EQ(read.components[0].count, 2u);
  EXPECT_EQ(read.components[0].atomsPerMolecule, 5u);
  EXPECT_EQ(read.positions.size(), 10u);
  EXPECT_NEAR(read.boxLengths[0], 20.0, 1e-9);
  EXPECT_EQ(read.forceField["MixingRule"], "Lorentz-Berthelot");
  EXPECT_NEAR(read.forceField["CutOffVDW"].get<double>(), 12.0, 1e-9);

  const nlohmann::json &definition = read.components[0].definition;
  ASSERT_EQ(definition["Bonds"].size(), 4u);
  EXPECT_EQ(definition["Bonds"][0][1], "HARMONIC");
  EXPECT_NEAR(definition["Bonds"][0][2][0].get<double>(), 301931.64, 0.05);
  EXPECT_NEAR(definition["Bonds"][0][2][1].get<double>(), 1.54, 1e-6);
  ASSERT_EQ(definition["Bends"].size(), 3u);
  EXPECT_NEAR(definition["Bends"][0][2][0].get<double>(), 62500.0, 0.05);
  ASSERT_EQ(definition["Torsions"].size(), 2u);
  EXPECT_EQ(definition["Torsions"][0][1], "POLYNOMIAL");
  EXPECT_NEAR(definition["Intra14ChargeChargeScalingValue"].get<double>(), 0.5, 1e-12);
  EXPECT_FALSE(definition.contains("Intra14VanDerWaalsScalingValue"));

  // pseudo-atoms keep their names, masses and charges
  bool foundOxygen = false;
  for (const auto &pseudoAtom : read.forceField["PseudoAtoms"])
  {
    if (pseudoAtom["name"] == "O")
    {
      foundOxygen = true;
      EXPECT_NEAR(pseudoAtom["mass"].get<double>(), 16.0, 1e-6);
      EXPECT_NEAR(pseudoAtom["charge"].get<double>(), -0.4, 1e-6);
    }
  }
  EXPECT_TRUE(foundOxygen);

  // the converted torsion reproduces RASPA's TRAPPE_EXTENDED energy
  std::mt19937 rng(29);
  std::vector<double> polynomial = definition["Torsions"][0][2].get<std::vector<double>>();
  TorsionPotential recovered({0, 1, 2, 3}, TorsionType::Polynomial, polynomial);
  TorsionPotential original({0, 1, 2, 3}, TorsionType::TraPPE_Extended, {1078.16, 355.03, 68.19, 791.32, 0.0});
  for (int sample = 0; sample < 20; ++sample)
  {
    Geometry g = randomGeometry(rng);
    EXPECT_NEAR(recovered.calculateEnergy(g.a, g.b, g.c, g.d), original.calculateEnergy(g.a, g.b, g.c, g.d), 1e-5);
  }

  // and the converted component can be read by RASPA again
  Component reread = readComponent(forceField, definition, "chain-reread");
  EXPECT_EQ(reread.atoms.size(), 5u);
  EXPECT_EQ(reread.intraMolecularPotentials.bonds.size(), 4u);
  EXPECT_EQ(reread.intraMolecularPotentials.torsions.size(), 2u);
}

TEST(lammps_export, reader_keeps_one_pseudo_atom_per_type_with_per_atom_charges)
{
  // one LAMMPS type (1) with two different charges in a 3-atom chain, and a type (2) with one charge
  const std::string data = R"(LAMMPS data file

3 atoms
2 bonds
1 angles
2 atom types
1 bond types
1 angle types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 12.0
2 16.0

Pair Coeffs # lj/cut/coul/long

1 0.1 3.4
2 0.2 3.0

Bond Coeffs # harmonic

1 300.0 1.5

Angle Coeffs # harmonic

1 60.0 110.0

Atoms # full

1 1 1 -0.3 5.0 5.0 5.0
2 1 2  0.6 6.5 5.0 5.0
3 1 1 -0.3 8.0 5.0 5.0

Bonds

1 1 1 2
2 1 2 3

Angles

1 1 1 2 3
)";
  const std::string input = R"(units real
atom_style full
pair_style lj/cut/coul/long 12.0
bond_style harmonic
angle_style harmonic
kspace_style pppm 1e-5
read_data mixed.data
)";
  TemporaryDirectory directory;
  std::filesystem::path dataPath = directory.path / "mixed.data";
  std::filesystem::path inputPath = directory.path / "mixed.in";
  std::ofstream(dataPath) << data;
  std::ofstream(inputPath) << input;

  // first: every atom of a type carries the same charge -> pseudo-atom charge, no per-atom entries
  {
    LAMMPS::ReadResult read = LAMMPS::readDataFile(dataPath, inputPath);
    ASSERT_EQ(read.forceField["PseudoAtoms"].size(), 2u);
    EXPECT_NEAR(read.forceField["PseudoAtoms"][0]["charge"].get<double>(), -0.3, 1e-9);
    EXPECT_NEAR(read.forceField["PseudoAtoms"][1]["charge"].get<double>(), 0.6, 1e-9);
    ASSERT_EQ(read.components.size(), 1u);
    for (const auto &atom : read.components[0].definition["PseudoAtoms"]) EXPECT_EQ(atom.size(), 2u);
  }

  // second: the two type-1 atoms carry different charges -> still one pseudo-atom (charge 0), charges per atom
  std::string mixed = data;
  mixed.replace(mixed.find("3 1 1 -0.3"), std::string("3 1 1 -0.3").size(), "3 1 1 -0.1");
  std::ofstream(dataPath, std::ios::trunc) << mixed;
  LAMMPS::ReadResult read = LAMMPS::readDataFile(dataPath, inputPath);

  ASSERT_EQ(read.forceField["PseudoAtoms"].size(), 2u);
  EXPECT_NEAR(read.forceField["PseudoAtoms"][0]["charge"].get<double>(), 0.0, 1e-12);
  EXPECT_NEAR(read.forceField["PseudoAtoms"][1]["charge"].get<double>(), 0.6, 1e-9);
  EXPECT_TRUE(std::any_of(read.warnings.begin(), read.warnings.end(),
                          [](const std::string &w) { return w.contains("different charges"); }));

  ASSERT_EQ(read.components.size(), 1u);
  const nlohmann::json &definition = read.components[0].definition;
  ASSERT_EQ(definition["PseudoAtoms"].size(), 3u);
  ASSERT_EQ(definition["PseudoAtoms"][0].size(), 3u);
  EXPECT_NEAR(definition["PseudoAtoms"][0][2].get<double>(), -0.3, 1e-9);
  EXPECT_EQ(definition["PseudoAtoms"][1].size(), 2u);  // equals the pseudo-atom charge
  ASSERT_EQ(definition["PseudoAtoms"][2].size(), 3u);
  EXPECT_NEAR(definition["PseudoAtoms"][2][2].get<double>(), -0.1, 1e-9);

  // the converted input reads back into RASPA with the per-atom charges
  EXPECT_EQ(read.forceField["PseudoAtoms"][0]["name"], "T1");
  EXPECT_EQ(read.forceField["PseudoAtoms"][1]["name"], "T2");
  ForceField forceField({{"T1", false, 12.0, 0.0, 0.0, 6, true}, {"T2", false, 16.0, 0.6, 0.0, 8, true}},
                        {{50.0, 3.4}, {100.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                        false, false, true);
  Component reread = readComponent(forceField, definition, "mixed-reread");
  ASSERT_EQ(reread.atoms.size(), 3u);
  EXPECT_NEAR(reread.atoms[0].charge, -0.3, 1e-9);
  EXPECT_NEAR(reread.atoms[1].charge, 0.6, 1e-9);
  EXPECT_NEAR(reread.atoms[2].charge, -0.1, 1e-9);
  EXPECT_EQ(reread.atoms[0].type, reread.atoms[2].type);
  EXPECT_NEAR(reread.netCharge, 0.2, 1e-9);
}

TEST(lammps_export, rigid_molecule_gets_fragment_ids_and_rigid_small)
{
  ForceField forceField = chainForceField();
  nlohmann::json definition = nlohmann::json::parse(R"({
    "CriticalTemperature": 300.0, "CriticalPressure": 7.0e6, "AcentricFactor": 0.2,
    "PseudoAtoms": [["O", [-1.16, 0.0, 0.0]], ["C", [0.0, 0.0, 0.0]], ["O", [1.16, 0.0, 0.0]]]
  })");
  Component rigid = readComponent(forceField, definition, "rigid3");
  ASSERT_TRUE(rigid.rigid);

  const SimulationBox box(15.0, 15.0, 15.0);
  System system(forceField, box, false, 300.0, 1e5, 1.0, {}, {rigid}, {}, {3}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t m = 0; m < 3; ++m)
    for (std::size_t i = 0; i < 3; ++i)
      atoms[m * 3 + i].position = double3(3.0 + 4.0 * static_cast<double>(m), 7.0, 7.0) + rigid.atoms[i].position;

  LAMMPS::Topology topology =
      LAMMPS::buildTopology(system.components, system.atomData, system.atomDynamics, system.moleculeData,
                            system.simulationBox, system.forceField, system.numberOfIntegerMoleculesPerComponent,
                            system.framework);
  EXPECT_EQ(topology.numberOfFragments, 3u);
  EXPECT_EQ(topology.atoms.size(), 9u);
  EXPECT_EQ(topology.atoms[0].fragment, 1u);
  EXPECT_EQ(topology.atoms[8].fragment, 3u);
  EXPECT_TRUE(topology.bondList.empty());

  std::string data = LAMMPS::writeDataFile(topology);
  std::string input = LAMMPS::writeInputScript(topology);
  EXPECT_NE(data.find("Fragments"), std::string::npos);
  EXPECT_NE(input.find("fix fragments all property/atom i_fragment"), std::string::npos);
  EXPECT_NE(input.find("rigid/small custom i_fragment"), std::string::npos);
  EXPECT_NE(input.find("neigh_modify exclude molecule/intra"), std::string::npos);
  EXPECT_NE(input.find("read_data raspa.data fix fragments NULL Fragments"), std::string::npos);
}
