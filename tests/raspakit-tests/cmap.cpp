#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import units;
import atom;
import atom_dynamics;
import forcefield;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import cmap_potential;
import running_energy;
import randomnumbers;
import system;
import simulationbox;
import integrators_update;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;

// The CMAP correction of the biomolecular force fields (CHARMM, ff19SB): a periodic bicubic-spline map E(phi, psi)
// over two consecutive dihedral angles; the map, its evaluation, the five-atom term, and its integration in the
// component reader, the spatial-decomposition engine and the analytic Hessian.
namespace
{
constexpr double pi = std::numbers::pi;

// a smooth periodic test surface and its derivatives
double surface(double phi, double psi) { return 100.0 * std::cos(phi) + 50.0 * std::sin(2.0 * psi) + 30.0 * std::cos(phi) * std::sin(psi); }
double surfacePhi(double phi, double psi) { return -100.0 * std::sin(phi) - 30.0 * std::sin(phi) * std::sin(psi); }
double surfacePsi(double phi, double psi) { return 100.0 * std::cos(2.0 * psi) + 30.0 * std::cos(phi) * std::cos(psi); }

CMAPMap makeSurfaceMap(std::size_t n, std::string name = "surface")
{
  std::vector<double> energies(n * n);
  const double h = 2.0 * pi / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      energies[i * n + j] = surface(-pi + static_cast<double>(i) * h, -pi + static_cast<double>(j) * h);
  return CMAPMap(std::move(name), n, std::move(energies));
}

std::array<double3, 5> randomPositions(RandomNumber& random)
{
  // a chain with bonds of about 1.5 Angstrom and well-defined dihedral angles
  std::array<double3, 5> positions{};
  positions[0] = double3(0.0, 0.0, 0.0);
  for (std::size_t k = 1; k < 5; ++k)
  {
    const double3 step(1.0 + 0.5 * random.uniform(), 0.8 * (random.uniform() - 0.5), 0.8 * (random.uniform() - 0.5));
    positions[k] = positions[k - 1] + step;
  }
  return positions;
}

class TemporaryDirectory
{
 public:
  explicit TemporaryDirectory(std::string name) : path(std::filesystem::temp_directory_path() / std::move(name))
  {
    std::filesystem::create_directories(path);
  }
  ~TemporaryDirectory()
  {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
  void write(const std::string& file, std::string_view contents) const
  {
    std::ofstream stream(path / file);
    stream << contents;
  }
  std::filesystem::path path;
};

std::string energiesJson(std::size_t n, double scale)
{
  std::string result = "[";
  const double h = 2.0 * pi / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
    {
      if (i + j > 0) result += ", ";
      result += std::format("{:.10f}", scale * surface(-pi + static_cast<double>(i) * h, -pi + static_cast<double>(j) * h));
    }
  return result + "]";
}

std::string forceFieldJson()
{
  return R"({
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "truncated",
  "TailCorrections": false,
  "CutOff": 9.0,
  "ChargeMethod": "None",
  "PseudoAtoms": [
    {"name": "C", "framework": false, "element": "C", "mass": 12.0, "charge": 0.0}
  ],
  "SelfInteractions": [
    {"name": "C", "type": "lennard-jones", "parameters": [40.0, 3.5]}
  ]
})";
}
}  // namespace

TEST(cmap, spline_reproduces_the_grid_and_converges_to_a_smooth_surface)
{
  // the spline passes through the grid values and is periodic in both angles
  const std::size_t n = 12;
  const CMAPMap map = makeSurfaceMap(n);
  const double h = 2.0 * pi / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
    {
      const double phi = -pi + static_cast<double>(i) * h;
      const double psi = -pi + static_cast<double>(j) * h;
      EXPECT_NEAR(map.evaluate(phi, psi).energy, map.energies[i * n + j], 1e-10);
      EXPECT_NEAR(map.evaluate(phi + 2.0 * pi, psi - 4.0 * pi).energy, map.energies[i * n + j], 1e-10);
    }

  // against an independent implementation of the same scheme (periodic C^2 splines along both angles, the cross
  // derivative from a spline through dE/dphi along psi, bicubic patch alpha = A F A^T)
  struct Reference
  {
    double phi, psi, energy, dPhi, dPsi;
  };
  for (const Reference& r : {Reference{0.3, -1.1, 29.552288929161, -21.623005907789, -46.152817882865},
                             Reference{2.9, 2.0, -161.378627942722, -30.485503527010, -53.812352940578},
                             Reference{-3.0, 0.1, -92.086755473415, 14.649034908596, 67.891105803625},
                             Reference{1.0, 1.0, 113.098147891247, -105.294095726641, -31.908632902830}})
  {
    const CMAPMap::Evaluation value = map.evaluate(r.phi, r.psi);
    EXPECT_NEAR(value.energy, r.energy, 1e-9);
    EXPECT_NEAR(value.dPhi, r.dPhi, 1e-9);
    EXPECT_NEAR(value.dPsi, r.dPsi, 1e-9);
  }

  // the derivatives of the patch are those of the surface it interpolates
  for (const auto [phi, psi] : {std::pair{0.3, -1.1}, std::pair{2.9, 2.0}, std::pair{-3.0, 0.1}})
  {
    const CMAPMap::Evaluation value = map.evaluate(phi, psi);
    const double d = 1e-5;
    EXPECT_NEAR(value.dPhi, (map.evaluate(phi + d, psi).energy - map.evaluate(phi - d, psi).energy) / (2 * d), 1e-6);
    EXPECT_NEAR(value.dPsi, (map.evaluate(phi, psi + d).energy - map.evaluate(phi, psi - d).energy) / (2 * d), 1e-6);
    EXPECT_NEAR(value.dPhiPhi, (map.evaluate(phi + d, psi).dPhi - map.evaluate(phi - d, psi).dPhi) / (2 * d), 1e-5);
    EXPECT_NEAR(value.dPsiPsi, (map.evaluate(phi, psi + d).dPsi - map.evaluate(phi, psi - d).dPsi) / (2 * d), 1e-5);
    EXPECT_NEAR(value.dPhiPsi, (map.evaluate(phi + d, psi).dPsi - map.evaluate(phi - d, psi).dPsi) / (2 * d), 1e-5);
    EXPECT_NEAR(value.dPhiPsi, (map.evaluate(phi, psi + d).dPhi - map.evaluate(phi, psi - d).dPhi) / (2 * d), 1e-5);
  }

  // C^1 across a cell boundary (the gradient is continuous at the nodes)
  {
    const double phi = -pi + 3.0 * h;
    const double psi = -pi + 7.0 * h;
    const CMAPMap::Evaluation left = map.evaluate(phi - 1e-9, psi + 0.37 * h);
    const CMAPMap::Evaluation right = map.evaluate(phi + 1e-9, psi + 0.37 * h);
    EXPECT_NEAR(left.dPhi, right.dPhi, 1e-6);
    EXPECT_NEAR(left.dPsi, right.dPsi, 1e-6);
    const CMAPMap::Evaluation below = map.evaluate(phi + 0.61 * h, psi - 1e-9);
    const CMAPMap::Evaluation above = map.evaluate(phi + 0.61 * h, psi + 1e-9);
    EXPECT_NEAR(below.dPhi, above.dPhi, 1e-6);
    EXPECT_NEAR(below.dPsi, above.dPsi, 1e-6);
  }

  // fourth-order convergence to the smooth surface
  double previous = 0.0;
  for (const std::size_t resolution : {12uz, 24uz, 48uz})
  {
    const CMAPMap fine = makeSurfaceMap(resolution);
    double maxError = 0.0, maxGradientError = 0.0;
    for (double phi = -3.1; phi < 3.1; phi += 0.37)
      for (double psi = -3.0; psi < 3.1; psi += 0.41)
      {
        const CMAPMap::Evaluation value = fine.evaluate(phi, psi);
        maxError = std::max(maxError, std::abs(value.energy - surface(phi, psi)));
        maxGradientError = std::max({maxGradientError, std::abs(value.dPhi - surfacePhi(phi, psi)),
                                     std::abs(value.dPsi - surfacePsi(phi, psi))});
      }
    if (previous > 0.0) EXPECT_LT(maxError, previous / 12.0) << "resolution " << resolution;
    if (resolution == 48) EXPECT_LT(maxGradientError, 2e-2);  // relative 1.5e-4 of the ~130 amplitude
    previous = maxError;
  }
  EXPECT_LT(previous, 1e-3);
}

TEST(cmap, map_validation_scaling_and_equality)
{
  EXPECT_THROW(CMAPMap("bad", 2, std::vector<double>(4, 0.0)), std::runtime_error);
  EXPECT_THROW(CMAPMap("bad", 4, std::vector<double>(15, 0.0)), std::runtime_error);

  CMAPMap map = makeSurfaceMap(8);
  const CMAPMap copy = map;
  EXPECT_TRUE(map == copy);
  const CMAPMap::Evaluation before = map.evaluate(0.7, -2.1);
  map.scaleEnergy(0.25);
  const CMAPMap::Evaluation after = map.evaluate(0.7, -2.1);
  EXPECT_NEAR(after.energy, 0.25 * before.energy, 1e-12);
  EXPECT_NEAR(after.dPhi, 0.25 * before.dPhi, 1e-12);
  EXPECT_NEAR(after.dPsiPsi, 0.25 * before.dPsiPsi, 1e-12);
  EXPECT_FALSE(map == copy);
}

TEST(cmap, dihedral_angle_follows_the_iupac_convention)
{
  // looking along B->C, D is rotated counter-clockwise from A by +90 degrees
  const double3 A(1.0, 0.0, 0.0), B(0.0, 0.0, 0.0), C(0.0, 0.0, 1.0), D(0.0, 1.0, 1.0);
  EXPECT_NEAR(CMAPPotential::dihedralAngle(A, B, C, D), 0.5 * pi, 1e-12);
  EXPECT_NEAR(CMAPPotential::dihedralAngle(D, C, B, A), 0.5 * pi, 1e-12);  // the reversed chain
  EXPECT_NEAR(CMAPPotential::dihedralAngle(A, B, C, double3(0.0, -1.0, 1.0)), -0.5 * pi, 1e-12);
  EXPECT_NEAR(CMAPPotential::dihedralAngle(A, B, C, double3(-1.0, 0.0, 1.0)), pi, 1e-12);  // trans
  EXPECT_NEAR(CMAPPotential::dihedralAngle(A, B, C, double3(1.0, 0.0, 1.0)), 0.0, 1e-12);  // cis

  // the gradient of the angle by central differences
  RandomNumber random(17);
  for (std::size_t trial = 0; trial < 5; ++trial)
  {
    std::array<double3, 5> positions = randomPositions(random);
    const auto [angle, gradient] =
        CMAPPotential::dihedralAngleAndGradient(positions[0], positions[1], positions[2], positions[3]);
    EXPECT_NEAR(angle, CMAPPotential::dihedralAngle(positions[0], positions[1], positions[2], positions[3]), 1e-14);
    const double d = 1e-6;
    for (std::size_t k = 0; k < 4; ++k)
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        std::array<double3, 5> plus = positions, minus = positions;
        plus[k][axis] += d;
        minus[k][axis] -= d;
        const double numerical = (CMAPPotential::dihedralAngle(plus[0], plus[1], plus[2], plus[3]) -
                                  CMAPPotential::dihedralAngle(minus[0], minus[1], minus[2], minus[3])) /
                                 (2.0 * d);
        EXPECT_NEAR(gradient[k][axis], numerical, 1e-7) << "atom " << k << " axis " << axis;
      }
  }
}

TEST(cmap, gradient_and_strain_match_finite_differences)
{
  const CMAPMap map = makeSurfaceMap(16);
  const CMAPPotential term({0, 1, 2, 3, 4}, 0);
  RandomNumber random(23);
  for (std::size_t trial = 0; trial < 6; ++trial)
  {
    const std::array<double3, 5> positions = randomPositions(random);
    const auto [energy, gradient, strain] =
        term.potentialEnergyGradientStrain(map, positions[0], positions[1], positions[2], positions[3], positions[4]);
    EXPECT_NEAR(energy, term.calculateEnergy(map, positions[0], positions[1], positions[2], positions[3], positions[4]),
                1e-12);

    const double d = 1e-6;
    double3 total{};
    for (std::size_t k = 0; k < 5; ++k)
    {
      total += gradient[k];
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        std::array<double3, 5> plus = positions, minus = positions;
        plus[k][axis] += d;
        minus[k][axis] -= d;
        const double numerical = (term.calculateEnergy(map, plus[0], plus[1], plus[2], plus[3], plus[4]) -
                                  term.calculateEnergy(map, minus[0], minus[1], minus[2], minus[3], minus[4])) /
                                 (2.0 * d);
        EXPECT_NEAR(gradient[k][axis], numerical, 1e-5 * std::max(1.0, std::abs(numerical)))
            << "trial " << trial << " atom " << k << " axis " << axis;
      }
    }
    // no net force
    EXPECT_NEAR(total.x, 0.0, 1e-10);
    EXPECT_NEAR(total.y, 0.0, 1e-10);
    EXPECT_NEAR(total.z, 0.0, 1e-10);

    // the strain derivative is dU/d(epsilon) of a homogeneous deformation r -> (1 + epsilon) r
    const auto energyUnderStrain = [&](std::size_t row, std::size_t column, double epsilon)
    {
      std::array<double3, 5> deformed = positions;
      for (double3& p : deformed) p[row] += epsilon * p[column];
      return term.calculateEnergy(map, deformed[0], deformed[1], deformed[2], deformed[3], deformed[4]);
    };
    const std::array<std::array<double, 3>, 3> analytic{{{strain.ax, strain.bx, strain.cx},
                                                          {strain.ay, strain.by, strain.cy},
                                                          {strain.az, strain.bz, strain.cz}}};
    for (std::size_t row = 0; row < 3; ++row)
      for (std::size_t column = 0; column < 3; ++column)
      {
        const double numerical = (energyUnderStrain(row, column, d) - energyUnderStrain(row, column, -d)) / (2.0 * d);
        // strain(ij) = sum_k (r_k)_j (g_k)_i: the (position component, gradient component) convention
        EXPECT_NEAR(analytic[row][column], numerical, 1e-5 * std::max(1.0, std::abs(numerical)))
            << "trial " << trial << " strain " << row << column;
      }
  }
}

TEST(cmap, component_reads_maps_and_terms_from_json)
{
  TemporaryDirectory directory("raspa_cmap_component_test");
  directory.write("force_field.json", forceFieldJson());
  const ForceField forceField((directory.path / "force_field.json").string());

  // the force-field map 'shared', the component map 'local' (in Kelvin)
  const std::string shared = energiesJson(8, 1.0);
  const std::string local = energiesJson(6, 2.0);
  directory.write("force_field.json", std::format(R"({{
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "truncated",
  "TailCorrections": false,
  "CutOff": 9.0,
  "ChargeMethod": "None",
  "PseudoAtoms": [ {{"name": "C", "framework": false, "element": "C", "mass": 12.0, "charge": 0.0}} ],
  "SelfInteractions": [ {{"name": "C", "type": "lennard-jones", "parameters": [40.0, 3.5]}} ],
  "CMAPs": [ {{"Name": "shared", "Resolution": 8, "Energies": {}}} ]
}})",
                                                  shared));
  const ForceField withMaps((directory.path / "force_field.json").string());
  ASSERT_EQ(withMaps.cmapMaps.size(), 1uz);
  EXPECT_EQ(withMaps.cmapMaps[0].name, "shared");
  EXPECT_EQ(withMaps.cmapMaps[0].resolution, 8uz);
  // the grid is converted from Kelvin to internal energy units
  EXPECT_NEAR(withMaps.cmapMaps[0].energies[0] * Units::EnergyToKelvin, surface(-pi, -pi), 1e-8);

  const std::string chain = std::format(R"({{
  "CriticalTemperature": 400.0, "CriticalPressure": 4.0e6, "AcentricFactor": 0.1,
  "PseudoAtoms": [["C", [0.0, 0.0, 0.0]], ["C", [1.5, 0.0, 0.0]], ["C", [2.3, 1.2, 0.0]], ["C", [3.8, 1.3, 0.4]],
                  ["C", [4.5, 2.5, 1.0]], ["C", [6.0, 2.6, 1.1]]],
  "Connectivity": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]],
  "Bonds": [[[0, 1], "HARMONIC", [1000.0, 1.5]], [[1, 2], "HARMONIC", [1000.0, 1.5]], [[2, 3], "HARMONIC", [1000.0, 1.5]],
            [[3, 4], "HARMONIC", [1000.0, 1.5]], [[4, 5], "HARMONIC", [1000.0, 1.5]]],
  "Bends": [[["C", "C", "C"], "HARMONIC", [500.0, 110.0]]],
  "Torsions": [[["C", "C", "C", "C"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]],
  "CMAPs": [ {{"Name": "local", "Resolution": 6, "Energies": {}}} ],
  "CMAPTorsions": [ [[0, 1, 2, 3, 4], "shared"], [[1, 2, 3, 4, 5], "local"], [[1, 2, 3, 4, 5], 0] ]
}})",
                                        local);
  directory.write("chain.json", chain);
  const Component component(Component::Type::Adsorbate, 0, withMaps, "chain", (directory.path / "chain").string(), 5,
                            21);
  const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
  ASSERT_EQ(potentials.cmapMaps.size(), 2uz);
  ASSERT_EQ(potentials.cmaps.size(), 3uz);
  // the component map is listed first, so index 0 is 'local'
  EXPECT_EQ(potentials.cmapMaps[potentials.cmaps[0].mapIndex].name, "shared");
  EXPECT_EQ(potentials.cmapMaps[potentials.cmaps[1].mapIndex].name, "local");
  EXPECT_EQ(potentials.cmapMaps[potentials.cmaps[2].mapIndex].name, "local");
  EXPECT_EQ(potentials.cmaps[0].identifiers, (std::array<std::size_t, 5>{0, 1, 2, 3, 4}));
  EXPECT_EQ(potentials.cmaps[1].identifiers, (std::array<std::size_t, 5>{1, 2, 3, 4, 5}));
  EXPECT_NEAR(potentials.cmapMaps[potentials.cmaps[1].mapIndex].energies[7] * Units::EnergyToKelvin,
              2.0 * surface(-pi + 2.0 * pi / 6.0, -pi + 2.0 * pi / 6.0), 1e-8);

  // the CMAP energy of the component geometry
  const RunningEnergy energy = potentials.computeInternalCMAPEnergies(component.atoms);
  double expected = 0.0;
  for (const CMAPPotential& term : potentials.cmaps)
  {
    const auto& a = component.atoms;
    expected += term.calculateEnergy(potentials.cmapMaps[term.mapIndex], a[term.identifiers[0]].position,
                                     a[term.identifiers[1]].position, a[term.identifiers[2]].position,
                                     a[term.identifiers[3]].position, a[term.identifiers[4]].position);
  }
  EXPECT_NE(expected, 0.0);
  EXPECT_NEAR(energy.cmap, expected, 1e-10 * std::abs(expected));
  EXPECT_NEAR(potentials.computeInternalEnergies(withMaps, SimulationBox(30.0, 30.0, 30.0), component.atoms).cmap,
              expected, 1e-10 * std::abs(expected));

  // errors: an unknown map, an atom out of range, a term with four atoms
  const auto rejects = [&](std::string_view terms)
  {
    directory.write("bad.json", std::format(R"({{
  "CriticalTemperature": 400.0, "CriticalPressure": 4.0e6, "AcentricFactor": 0.1,
  "PseudoAtoms": [["C", [0.0, 0.0, 0.0]], ["C", [1.5, 0.0, 0.0]], ["C", [2.3, 1.2, 0.0]], ["C", [3.8, 1.3, 0.4]],
                  ["C", [4.5, 2.5, 1.0]]],
  "Connectivity": [[0, 1], [1, 2], [2, 3], [3, 4]],
  "CMAPTorsions": [ {} ]
}})",
                                            terms));
    EXPECT_THROW(Component(Component::Type::Adsorbate, 0, withMaps, "bad", (directory.path / "bad").string(), 5, 21),
                 std::runtime_error)
        << terms;
  };
  rejects(R"([[0, 1, 2, 3, 4], "missing"])");
  rejects(R"([[0, 1, 2, 3, 5], "shared"])");
  rejects(R"([[0, 1, 2, 3], "shared"])");
  rejects(R"([[0, 1, 2, 3, 4], 3])");
}

namespace
{
/// A chain of n^3 Lennard-Jones beads on a boustrophedon path through an n x n x n lattice with harmonic bonds and
/// bends and a CMAP term on every five consecutive beads (two maps, alternating).
System makeCMAPChainSystem(RandomNumber& random, std::size_t n, std::size_t numberOfMolecules)
{
  ForceField forceField = ForceField({{"CH2", false, 14.02658, 0.0, 0.0, 6, false}}, {{46.0, 3.95}},
                                     ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, false);

  const std::size_t numberOfBeads = n * n * n;
  ConnectivityTable connectivityTable(numberOfBeads);
  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i + 1 < numberOfBeads; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
    potentials.bonds.push_back(BondPotential({i, i + 1}, BondType::Harmonic, {96500.0, 3.8}));
  }
  for (std::size_t i = 0; i + 2 < numberOfBeads; ++i)
  {
    potentials.bends.push_back(BendPotential({i, i + 1, i + 2}, BendType::Harmonic, {62500.0, 100.0}));
  }
  potentials.cmapMaps = {makeSurfaceMap(12, "alpha"), makeSurfaceMap(24, "beta")};
  potentials.cmapMaps[1].scaleEnergy(-0.5);
  for (std::size_t i = 0; i + 4 < numberOfBeads; ++i)
  {
    potentials.addCMAP({i, i + 1, i + 2, i + 3, i + 4}, (i % 2 == 0) ? "alpha" : "beta");
  }

  std::vector<Atom> beads{};
  const double spacing = 3.8;
  const double offset = 0.5 * spacing * static_cast<double>(n - 1);
  for (std::size_t z = 0; z < n; ++z)
    for (std::size_t row = 0; row < n; ++row)
    {
      const std::size_t y = (z % 2 == 0) ? row : n - 1 - row;
      for (std::size_t column = 0; column < n; ++column)
      {
        const std::size_t x = ((z * n + row) % 2 == 0) ? column : n - 1 - column;
        beads.push_back(Atom({spacing * static_cast<double>(x) - offset, spacing * static_cast<double>(y) - offset,
                              spacing * static_cast<double>(z) - offset},
                             0.0, 1.0, 0, 0, 0, false, false));
      }
    }
  Component chain =
      Component(forceField, "cmap-chain", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  std::vector<double3> positions{};
  for (std::size_t m = 0; m < numberOfMolecules; ++m)
    for (const Atom& bead : beads) positions.push_back(bead.position);
  System system = System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain},
                         {positions}, {0}, 5);

  // random placement and jitter (some molecules straddle the box boundary)
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t m = 0; m < numberOfMolecules; ++m)
  {
    const double3 center = system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
    for (std::size_t k = 0; k < numberOfBeads; ++k)
    {
      const double3 local =
          beads[k].position + 0.3 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
      atoms[m * numberOfBeads + k].position = center + rotation * local;
    }
  }
  return system;
}

RunningEnergy referenceGradients(System& system)
{
  return Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
}

std::vector<double3> gradientsOf(const System& system)
{
  std::vector<double3> result;
  for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) result.push_back(dynamics.gradient);
  return result;
}

double rmsDifference(const std::vector<double3>& a, const std::vector<double3>& b)
{
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += double3::dot(a[i] - b[i], a[i] - b[i]);
  return std::sqrt(sum / static_cast<double>(a.size()));
}

double rmsNorm(const std::vector<double3>& a)
{
  double sum = 0.0;
  for (const double3& g : a) sum += double3::dot(g, g);
  return std::sqrt(sum / static_cast<double>(a.size()));
}

double maxAbsDifference(const double3x3& a, const double3x3& b)
{
  return std::max({std::abs(a.ax - b.ax), std::abs(a.ay - b.ay), std::abs(a.az - b.az), std::abs(a.bx - b.bx),
                   std::abs(a.by - b.by), std::abs(a.bz - b.bz), std::abs(a.cx - b.cx), std::abs(a.cy - b.cy),
                   std::abs(a.cz - b.cz)});
}

double maxAbs(const double3x3& a)
{
  return std::max({std::abs(a.ax), std::abs(a.ay), std::abs(a.az), std::abs(a.bx), std::abs(a.by), std::abs(a.bz),
                   std::abs(a.cx), std::abs(a.cy), std::abs(a.cz)});
}
}  // namespace

TEST(cmap, molecule_energy_gradient_and_strain_are_consistent)
{
  // the per-molecule bonded routines: the CMAP energy is booked, the gradient is the derivative of the energy, and
  // the strain derivative that of the molecule under a homogeneous deformation
  RandomNumber random(5);
  System system = makeCMAPChainSystem(random, 3, 1);
  const Component& component = system.components[0];
  const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
  ASSERT_EQ(potentials.cmaps.size(), 23uz);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  for (AtomDynamics& d : dynamics) d.gradient = double3{};
  const auto [energy, strain] = potentials.computeInternalBondedStrainDerivative(system.simulationBox, atoms, dynamics);
  EXPECT_NE(energy.cmap, 0.0);
  EXPECT_NEAR(energy.cmap, potentials.computeInternalCMAPEnergies(atoms).cmap, 1e-10 * std::abs(energy.cmap));
  EXPECT_NEAR(energy.cmap, potentials.computeInternalEnergies(system.forceField, system.simulationBox, atoms).cmap,
              1e-10 * std::abs(energy.cmap));

  // the gradient and the strain by central differences of the total bonded energy
  const auto bondedEnergy = [&]()
  {
    const RunningEnergy e = potentials.computeInternalEnergies(system.forceField, system.simulationBox, atoms);
    return e.bond + e.bend + e.cmap;
  };
  const double d = 1e-5;
  for (std::size_t k = 0; k < atoms.size(); k += 5)
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const double original = atoms[k].position[axis];
      atoms[k].position[axis] = original + d;
      const double plus = bondedEnergy();
      atoms[k].position[axis] = original - d;
      const double minus = bondedEnergy();
      atoms[k].position[axis] = original;
      EXPECT_NEAR(dynamics[k].gradient[axis], (plus - minus) / (2.0 * d), 1e-4 * std::max(1.0, std::abs(plus - minus) / (2.0 * d)))
          << "atom " << k << " axis " << axis;
    }
  const std::vector<double3> saved = [&]
  {
    std::vector<double3> p;
    for (const Atom& atom : atoms) p.push_back(atom.position);
    return p;
  }();
  const std::array<std::array<double, 3>, 3> analytic{
      {{strain.ax, strain.bx, strain.cx}, {strain.ay, strain.by, strain.cy}, {strain.az, strain.bz, strain.cz}}};
  for (std::size_t row = 0; row < 3; ++row)
    for (std::size_t column = 0; column < 3; ++column)
    {
      double values[2];
      for (std::size_t s = 0; s < 2; ++s)
      {
        const double epsilon = (s == 0) ? d : -d;
        for (std::size_t k = 0; k < atoms.size(); ++k)
        {
          atoms[k].position = saved[k];
          atoms[k].position[row] += epsilon * saved[k][column];
        }
        values[s] = bondedEnergy();
      }
      for (std::size_t k = 0; k < atoms.size(); ++k) atoms[k].position = saved[k];
      const double numerical = (values[0] - values[1]) / (2.0 * d);
      EXPECT_NEAR(analytic[row][column], numerical, 1e-4 * std::max(1.0, std::abs(numerical))) << row << column;
    }
}

TEST(cmap, spatial_decomposition_engine_matches_the_molecule_reference)
{
  // 125-bead chains: the host bonded path of the engine partitions the CMAP terms over the threads
  RandomNumber random(31);
  System system = makeCMAPChainSystem(random, 5, 3);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;
  EXPECT_NE(reference.cmap, 0.0);

  for (std::size_t threads : {1uz, 3uz, 8uz})
  {
    SpatialDecompositionSettings settings;
    settings.numberOfThreads = threads;
    settings.verletSkin = 1.5;
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    const std::string label = std::format("threads {}", threads);

    EXPECT_NEAR(energy.cmap, reference.cmap, 1e-10 * std::abs(reference.cmap)) << label;
    EXPECT_NEAR(energy.bond, reference.bond, 1e-10 * std::abs(reference.bond)) << label;
    EXPECT_NEAR(energy.bend, reference.bend, 1e-10 * std::abs(reference.bend)) << label;
    EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 1e-8 * std::abs(reference.potentialEnergy()))
        << label;
    EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-8 * rmsNorm(referenceGradient)) << label;
    EXPECT_LT(maxAbsDifference(engine.molecularPressureTensor(), referencePressure), 1e-8 * maxAbs(referencePressure))
        << label;
  }
}

TEST(cmap, analytic_hessian_matches_finite_differences_of_the_gradient)
{
  RandomNumber random(41);
  System system = makeCMAPChainSystem(random, 2, 1);  // 8 beads, 4 CMAP terms
  ASSERT_EQ(system.components[0].intraMolecularPotentials.cmaps.size(), 4uz);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 3 * atoms.size());

  const auto evaluate = [&](bool withHessian, GeneralizedHessian& hessian, std::vector<double>& gradient)
  {
    std::ranges::fill(gradient, 0.0);
    DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = withHessian};
    DerivativeResults results{.gradient = gradient, .hessian = hessian};
    evaluateDerivatives(system, layout, capabilities, results);
    return results.energy;
  };

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  const double energy = evaluate(true, hessian, gradient);
  const RunningEnergy components = system.components[0].intraMolecularPotentials.computeInternalEnergies(
      system.forceField, system.simulationBox, atoms);
  EXPECT_NE(components.cmap, 0.0);
  EXPECT_NEAR(energy, components.potentialEnergy(), 1e-8 * std::abs(components.potentialEnergy()));

  GeneralizedHessian scratch(layout.numDofs(), 0);
  std::vector<double> plus(layout.numDofs()), minus(layout.numDofs());
  const double d = 1e-5;
  double maxAbsHessian = 0.0;
  for (const double value : hessian.positionPosition()) maxAbsHessian = std::max(maxAbsHessian, std::abs(value));
  ASSERT_GT(maxAbsHessian, 0.0);
  for (std::size_t k = 0; k < atoms.size(); ++k)
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const std::size_t column = layout.flexibleAtomDofBase(0, k) + axis;
      const double original = atoms[k].position[axis];
      atoms[k].position[axis] = original + d;
      evaluate(false, scratch, plus);
      atoms[k].position[axis] = original - d;
      evaluate(false, scratch, minus);
      atoms[k].position[axis] = original;
      for (std::size_t row = 0; row < layout.numDofs(); ++row)
      {
        const double numerical = (plus[row] - minus[row]) / (2.0 * d);
        EXPECT_NEAR(hessian(row, column), numerical, 1e-5 * maxAbsHessian) << "row " << row << " column " << column;
      }
    }
}
