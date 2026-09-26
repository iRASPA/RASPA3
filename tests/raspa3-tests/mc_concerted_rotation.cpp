#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import double3x3;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import units;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_pivot;
import mc_moves_concerted_rotation;
import mc_moves_concerted_rotation_geometry;

// Tests for the concerted-rotation (ConRot) move of Dodd, Boone and Theodorou, Mol. Phys. 78, 961
// (1993). The rebridging kernel is checked on its own (the solver must reproduce the original
// trimer when the driver is not rotated, every solution must satisfy the bond-length and bend-angle
// constraints, and the closure Jacobian must agree with finite differences), and the complete move
// is checked as a sampler: on a chain with fixed bonds and bends the torsion angles are independent
// with density proportional to exp(-beta U_torsion), so both the U = 0 case (uniform torsions) and
// the TraPPE case have exact reference distributions. The Jacobian ratio in the acceptance rule is
// essential for the sampled distribution, so these tests are sensitive to it.

namespace
{

constexpr double kBondLength = 1.54;
constexpr double kBendAngleDegrees = 114.0;

// Linear 10-bead chain with FIXED bonds and bends: the moves must preserve them exactly, and CBMC
// growth at creation then puts every bond and bend exactly at its fixed value.
std::string decaneJson(std::string_view torsionParameters)
{
  return std::format(R"({{
  "CriticalTemperature" : 617.7,
  "CriticalPressure" : 2110000.0,
  "AcentricFactor" : 0.492,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH2", [3.08, 0.0, 0.0]], ["CH2", [4.62, 0.0, 0.0]],
      ["CH2", [6.16, 0.0, 0.0]], ["CH2", [7.70, 0.0, 0.0]], ["CH2", [9.24, 0.0, 0.0]], ["CH2", [10.78, 0.0, 0.0]],
      ["CH2", [12.32, 0.0, 0.0]], ["CH2", [13.86, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "FIXED", [114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [{}]]
  ],
  "VanDerWaals" : "auto"
}}
)",
                     torsionParameters);
}

// Branched chain: a CH3 side group on backbone atoms 2 and 5 of an 8-bead backbone. Side groups
// inside a window are carried rigidly, so all bond lengths must be preserved exactly.
constexpr std::string_view kBranchedJson =
R"({
  "CriticalTemperature" : 617.7,
  "CriticalPressure" : 2110000.0,
  "AcentricFactor" : 0.492,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH2", [3.08, 0.0, 0.0]], ["CH2", [4.62, 0.0, 0.0]],
      ["CH2", [6.16, 0.0, 0.0]], ["CH2", [7.70, 0.0, 0.0]], ["CH2", [9.24, 0.0, 0.0]], ["CH2", [10.78, 0.0, 0.0]],
      ["CH3", [3.08, 1.54, 0.0]], ["CH3", [7.70, 1.54, 0.0]]
    ],
  "Connectivity" : [
    [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [2, 8], [5, 9]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "FIXED", [1.54]],
    [["CH2", "CH3"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [62500.0, 114.0]],
    [["CH2", "CH2", "CH3"], "HARMONIC", [62500.0, 112.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH3", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}, {"CH3", false, 15.04, 0.0, 0.0, 6, false}},
                    {{0.0, 3.95}, {0.0, 3.75}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

Component makeComponent(const ForceField& forceField, std::string name, std::string_view json)
{
  TemporaryFile file(name + ".json", json);
  return Component(Component::Type::Adsorbate, 0, forceField, name, file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false);
}

// Signed dihedral angle of the atom sequence a-b-c-d, in [-pi, pi].
double dihedralAngle(const double3& a, const double3& b, const double3& c, const double3& d)
{
  double3 b1 = b - a;
  double3 b2 = c - b;
  double3 b3 = d - c;
  double3 n1 = double3::cross(b1, b2);
  double3 n2 = double3::cross(b2, b3);
  double3 m1 = double3::cross(n1, b2.normalized());
  return std::atan2(double3::dot(m1, n2), double3::dot(n1, n2));
}

double cosBendAngle(const double3& a, const double3& b, const double3& c)
{
  return double3::dot((a - b).normalized(), (c - b).normalized());
}

// Builds a chain with constant bond length and bend angle and the given torsion angles (one per
// interior bond), starting from a fixed reference placement of the first three atoms.
std::vector<double3> buildChain(double bondLength, double bendAngle, const std::vector<double>& torsions)
{
  std::vector<double3> chain;
  chain.emplace_back(0.0, 0.0, 0.0);
  chain.emplace_back(bondLength, 0.0, 0.0);
  chain.emplace_back(bondLength - bondLength * std::cos(bendAngle), bondLength * std::sin(bendAngle), 0.0);
  for (double phi : torsions)
  {
    const double3& a = chain[chain.size() - 3];
    const double3& b = chain[chain.size() - 2];
    const double3& c = chain[chain.size() - 1];
    double3 bc = (c - b).normalized();
    double3 n = double3::cross(b - a, bc).normalized();
    double3 m = double3::cross(n, bc);
    double3 d2(-bondLength * std::cos(bendAngle), bondLength * std::sin(bendAngle) * std::cos(phi),
               bondLength * std::sin(bendAngle) * std::sin(phi));
    chain.push_back(c + bc * d2.x + m * d2.y + n * d2.z);
  }
  return chain;
}

// Pearson chi-squared statistic of 'counts' against the expected bin probabilities.
double chiSquared(const std::vector<std::size_t>& counts, const std::vector<double>& probabilities)
{
  std::size_t total = std::accumulate(counts.begin(), counts.end(), 0uz);
  double statistic = 0.0;
  for (std::size_t bin = 0; bin != counts.size(); ++bin)
  {
    double expected = static_cast<double>(total) * probabilities[bin];
    double deviation = static_cast<double>(counts[bin]) - expected;
    statistic += deviation * deviation / expected;
  }
  return statistic;
}

struct TorsionStatistics
{
  std::vector<std::vector<std::size_t>> counts;
  std::vector<double> sumCos;
  std::vector<double> sumSin;
  std::size_t samples{};

  TorsionStatistics(std::size_t torsions, std::size_t bins)
      : counts(torsions, std::vector<std::size_t>(bins, 0uz)), sumCos(torsions, 0.0), sumSin(torsions, 0.0)
  {
  }

  void sample(std::span<const Atom> atoms)
  {
    ++samples;
    const std::size_t bins = counts.front().size();
    for (std::size_t torsion = 0; torsion != counts.size(); ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      sumCos[torsion] += std::cos(phi);
      sumSin[torsion] += std::sin(phi);
      std::size_t bin = std::min(
          bins - 1, static_cast<std::size_t>((phi + std::numbers::pi) / (2.0 * std::numbers::pi) *
                                             static_cast<double>(bins)));
      ++counts[torsion][bin];
    }
  }
};

// Largest deviation of the backbone bond lengths and bend cosines of a linear chain from the fixed
// values.
std::pair<double, double> backboneDeviation(std::span<const Atom> atoms)
{
  double bondDeviation = 0.0, bendDeviation = 0.0;
  const double cosFixed = std::cos(kBendAngleDegrees * std::numbers::pi / 180.0);
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    bondDeviation = std::max(bondDeviation, std::abs((atoms[i + 1].position - atoms[i].position).length() - kBondLength));
  }
  for (std::size_t i = 0; i + 2 < atoms.size(); ++i)
  {
    bendDeviation = std::max(
        bendDeviation, std::abs(cosBendAngle(atoms[i].position, atoms[i + 1].position, atoms[i + 2].position) - cosFixed));
  }
  return {bondDeviation, bendDeviation};
}

}  // namespace

// The rebridging solver must (i) return the original trimer among its solutions when called with the
// original driver position, and (ii) return only solutions that satisfy all four bond lengths and
// five bend angles of the window.
TEST(MC_CONCERTED_ROTATION, rebridge_reproduces_original_trimer_and_satisfies_constraints)
{
  RandomNumber random(7);
  const double theta = kBendAngleDegrees * std::numbers::pi / 180.0;

  std::size_t reproduced = 0;
  std::size_t totalSolutions = 0;
  constexpr std::size_t trials = 300;
  double maxConstraintViolation = 0.0;

  for (std::size_t trial = 0; trial != trials; ++trial)
  {
    std::vector<double> torsions(6);
    for (double& phi : torsions) phi = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
    std::vector<double3> chain = buildChain(kBondLength, theta, torsions);  // 9 atoms: a0 .. a8
    ASSERT_EQ(chain.size(), 9uz);

    std::span<const double3, 7> a1_to_a7(chain.data() + 1, 7);
    ConcertedRotation::BackboneGeometry geometry = ConcertedRotation::BackboneGeometry::fromPositions(a1_to_a7);
    std::vector<ConcertedRotation::Trimer> solutions =
        ConcertedRotation::rebridge(chain[1], chain[2], chain[6], chain[7], geometry);

    bool found = false;
    for (const ConcertedRotation::Trimer& trimer : solutions)
    {
      ++totalSolutions;
      if ((trimer.a3 - chain[3]).length() < 1e-6 && (trimer.a4 - chain[4]).length() < 1e-6 &&
          (trimer.a5 - chain[5]).length() < 1e-6)
      {
        found = true;
      }

      const std::array<double3, 7> window{chain[1], chain[2], trimer.a3, trimer.a4, trimer.a5, chain[6], chain[7]};
      for (std::size_t i = 1; i + 1 < 6; ++i)
      {
        maxConstraintViolation =
            std::max(maxConstraintViolation, std::abs((window[i + 1] - window[i]).length() - kBondLength));
      }
      for (std::size_t i = 0; i + 2 < 7; ++i)
      {
        maxConstraintViolation = std::max(
            maxConstraintViolation, std::abs(cosBendAngle(window[i], window[i + 1], window[i + 2]) - std::cos(theta)));
      }
    }
    if (found) ++reproduced;
  }

  // Root finding on a grid can in principle miss a double root between two grid points; that must
  // be extremely rare (the move rejects such proposals).
  EXPECT_GE(reproduced, trials - 1);
  EXPECT_GE(totalSolutions, trials);  // at least the identity solution per window
  EXPECT_LT(maxConstraintViolation, 1e-8);
}

// The analytic closure Jacobian must agree with a central finite difference of the chart
// (a6, the two perpendicular components of a7, the out-of-plane component of a8) with respect to
// the six torsion-like rotations about the bonds a1a2 ... a6a7.
TEST(MC_CONCERTED_ROTATION, closure_jacobian_matches_finite_differences)
{
  RandomNumber random(11);
  const double theta = kBendAngleDegrees * std::numbers::pi / 180.0;

  for (std::size_t trial = 0; trial != 20; ++trial)
  {
    std::vector<double> torsions(6);
    for (double& phi : torsions) phi = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
    std::vector<double3> chain = buildChain(kBondLength, theta, torsions);  // a0 .. a8

    for (bool withA8 : {true, false})
    {
      std::span<const double3, 7> a1_to_a7(chain.data() + 1, 7);
      std::optional<double3> a8 = withA8 ? std::optional<double3>(chain[8]) : std::nullopt;
      const double analytic = ConcertedRotation::closureJacobian(a1_to_a7, a8);

      const std::size_t n = withA8 ? 6 : 5;
      const double3 u67 = (chain[7] - chain[6]).normalized();
      double3 helper = std::abs(u67.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
      const double3 e1 = double3::cross(u67, helper).normalized();
      const double3 e2 = double3::cross(u67, e1);
      const double3 normal = double3::cross(u67, chain[8] - chain[7]).normalized();

      auto chart = [&](const std::vector<double3>& c) -> std::vector<double>
      {
        std::vector<double> values{c[6].x, c[6].y, c[6].z, double3::dot(e1, c[7]), double3::dot(e2, c[7])};
        if (withA8) values.push_back(double3::dot(normal, c[8]));
        return values;
      };

      constexpr double epsilon = 1e-6;
      std::vector<double> matrix(n * n, 0.0);
      for (std::size_t k = 1; k <= n; ++k)
      {
        auto rotated = [&](double angle)
        {
          std::vector<double3> c = chain;
          const double3 axis = (chain[k + 1] - chain[k]).normalized();
          for (std::size_t m = k + 2; m < c.size(); ++m)
          {
            c[m] = ConcertedRotation::rotateAboutAxis(chain[k], axis, angle, chain[m]);
          }
          return chart(c);
        };
        std::vector<double> plus = rotated(epsilon), minus = rotated(-epsilon);
        for (std::size_t row = 0; row != n; ++row)
        {
          matrix[row * n + (k - 1)] = (plus[row] - minus[row]) / (2.0 * epsilon);
        }
      }
      const double numeric = std::abs(ConcertedRotation::determinant(matrix, n));
      EXPECT_NEAR(analytic, numeric, 1e-6 * std::max(1.0, numeric)) << "trial " << trial << " a8 " << withA8;
    }
  }
}

// A linear 10-bead chain has six windows: a1 = 1, 2, 3 in the forward direction (a0 = 0, 1, 2) and
// their mirror images in the backward direction; a8 exists for all but the two windows that reach a
// chain end.
TEST(MC_CONCERTED_ROTATION, window_enumeration_linear_and_branched)
{
  const ForceField forceField = makeZeroForceField();
  Component decane = makeComponent(forceField, "conrot-decane", decaneJson("0.0, 0.0, 0.0, 0.0"));

  const std::vector<Component::ConcertedRotationWindow>& windows = decane.concertedRotationWindows();
  EXPECT_EQ(windows.size(), 6uz);
  std::size_t withA8 = 0;
  for (const Component::ConcertedRotationWindow& window : windows)
  {
    if (window.a8.has_value()) ++withA8;
    for (const std::vector<std::size_t>& group : window.substituents) EXPECT_TRUE(group.empty());
  }
  EXPECT_EQ(withA8, 4uz);

  // Branched chain: side groups on atoms 2 and 5 must be carried by their backbone atom whenever it
  // is one of a2 ... a5, and windows in which the side-group atom itself is a backbone member are
  // valid as well (the backbone may run through the branch as a chain end).
  Component branched = makeComponent(forceField, "conrot-branched", kBranchedJson);
  const std::vector<Component::ConcertedRotationWindow>& branchedWindows = branched.concertedRotationWindows();
  EXPECT_FALSE(branchedWindows.empty());
  for (const Component::ConcertedRotationWindow& window : branchedWindows)
  {
    for (std::size_t i = 2; i <= 5; ++i)
    {
      const std::size_t atom = window.backbone[i];
      const std::vector<std::size_t>& group = window.substituents[i - 2];
      if (atom == 2)
      {
        EXPECT_TRUE(std::find(group.begin(), group.end(), 8uz) != group.end() ||
                    std::find(window.backbone.begin(), window.backbone.end(), 8uz) != window.backbone.end());
      }
      if (atom == 5)
      {
        EXPECT_TRUE(std::find(group.begin(), group.end(), 9uz) != group.end() ||
                    std::find(window.backbone.begin(), window.backbone.end(), 9uz) != window.backbone.end());
      }
    }
  }
}

// Number of closure solutions of the window a0 = 1 ... a7 = 8 at the current driver angle (the
// current trimer is one of them). Its equilibrium distribution is a sensitive probe of the
// solution-count ratio in the acceptance rule.
std::size_t solutionCountOfWindow(std::span<const Atom> atoms)
{
  std::array<double3, 7> window{};
  for (std::size_t k = 0; k != 7; ++k) window[k] = atoms[k + 2].position;
  const std::span<const double3, 7> span(window);
  return ConcertedRotation::rebridge(window[0], window[1], window[5], window[6],
                                     ConcertedRotation::BackboneGeometry::fromPositions(span))
      .size();
}

// U = 0 sampler test: with zero torsion parameters, fixed bonds and bends and no non-bonded
// interactions, all seven torsion angles must be uniform on [-pi, pi]. The concerted rotation is
// mixed with a few pivot moves (exact at U = 0) only to guarantee ergodicity of the chain ends; the
// concerted rotation dominates, so a missing or inverted Jacobian ratio biases the torsion
// histograms strongly (chi-squared in the hundreds to thousands with these settings). The
// solution-count ratio hardly affects the torsion marginals, but it does shift the distribution
// of the number of closure solutions, which is therefore compared against a pivot-only reference.
TEST(MC_CONCERTED_ROTATION, u0_torsions_uniform_and_constraints_preserved)
{
  const ForceField forceField = makeZeroForceField();
  Component decane = makeComponent(forceField, "conrot-decane-u0", decaneJson("0.0, 0.0, 0.0, 0.0"));

  constexpr std::size_t numberOfMoves = 250000;
  constexpr std::size_t burnIn = 5000;
  constexpr std::size_t sampleEvery = 6;
  constexpr std::size_t bins = 12;
  constexpr std::size_t countBins = 10;  // last bin: >= 9 solutions

  auto run = [&](bool withConcertedRotation, std::uint64_t seed, TorsionStatistics& statistics,
                 std::array<double, countBins>& countFractions, double& acceptance, double& maxBond, double& maxBend)
  {
    System system =
        System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {decane}, {}, {1}, 5);
    system.runningEnergies = system.computeTotalEnergies();
    system.components[0].concertedRotationRandomizationFraction = 0.5;
    system.components[0].pivotRandomizationFraction = 1.0;

    RandomNumber random(seed);
    std::span<Atom> atoms = system.spanOfMolecule(0, 0);
    ASSERT_EQ(atoms.size(), 10uz);

    std::array<std::size_t, countBins> counts{};
    std::size_t attempts = 0, accepted = 0;
    for (std::size_t i = 0; i != numberOfMoves; ++i)
    {
      std::optional<RunningEnergy> energyDifference;
      if (!withConcertedRotation || i % 10 == 0)
      {
        energyDifference = MC_Moves::pivotMove(random, system, 0, 0);
      }
      else
      {
        ++attempts;
        energyDifference = MC_Moves::concertedRotationMove(random, system, 0, 0);
        if (energyDifference.has_value()) ++accepted;
      }
      if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

      if (i < burnIn || i % sampleEvery != 0) continue;
      statistics.sample(atoms);
      ++counts[std::min(countBins - 1, solutionCountOfWindow(atoms))];
      auto [bond, bend] = backboneDeviation(atoms);
      maxBond = std::max(maxBond, bond);
      maxBend = std::max(maxBend, bend);
    }
    for (std::size_t n = 0; n != countBins; ++n)
    {
      countFractions[n] = static_cast<double>(counts[n]) / static_cast<double>(statistics.samples);
    }
    acceptance = attempts == 0 ? 0.0 : static_cast<double>(accepted) / static_cast<double>(attempts);

    RunningEnergy recomputed = system.computeTotalEnergies();
    EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
  };

  TorsionStatistics statistics(7, bins), referenceStatistics(7, bins);
  std::array<double, countBins> countFractions{}, referenceCountFractions{};
  double acceptance = 0.0, referenceAcceptance = 0.0;
  double maxBond = 0.0, maxBend = 0.0, referenceMaxBond = 0.0, referenceMaxBend = 0.0;
  run(true, 2024, statistics, countFractions, acceptance, maxBond, maxBend);
  run(false, 4048, referenceStatistics, referenceCountFractions, referenceAcceptance, referenceMaxBond,
      referenceMaxBend);

  // The move is not rejection-free even at U = 0 (Jacobian ratio, solution counts), but it must
  // be accepted often enough to be useful.
  EXPECT_GT(acceptance, 0.2);
  EXPECT_LT(acceptance, 1.0);

  EXPECT_LT(maxBond, 1e-9);
  EXPECT_LT(maxBend, 1e-9);

  // Uniform torsions; the chi-squared threshold (11 degrees of freedom) allows for the residual
  // correlation between successive samples.
  const std::vector<double> uniform(bins, 1.0 / static_cast<double>(bins));
  for (std::size_t torsion = 0; torsion != 7; ++torsion)
  {
    EXPECT_LT(chiSquared(statistics.counts[torsion], uniform), 60.0) << "torsion " << torsion;
    EXPECT_NEAR(statistics.sumCos[torsion] / static_cast<double>(statistics.samples), 0.0, 0.03)
        << "torsion " << torsion;
    EXPECT_NEAR(statistics.sumSin[torsion] / static_cast<double>(statistics.samples), 0.0, 0.03)
        << "torsion " << torsion;
  }

  // Distribution of the number of closure solutions: the fractions of states with few (<= 2) and
  // many (>= 8) solutions must agree with the exact pivot-only sampler. Dropping the solution-count
  // ratio shifts the <= 2 fraction by about +0.03, dropping the Jacobian ratio by about -0.04.
  auto fewSolutions = [](const std::array<double, countBins>& f) { return f[0] + f[1] + f[2]; };
  auto manySolutions = [](const std::array<double, countBins>& f) { return f[8] + f[9]; };
  EXPECT_NEAR(fewSolutions(countFractions), fewSolutions(referenceCountFractions), 0.015);
  EXPECT_NEAR(manySolutions(countFractions), manySolutions(referenceCountFractions), 0.015);
}

// Boltzmann sampler test with the TraPPE alkane torsion potential at 600 K: each torsion angle is
// independently distributed as exp(-U(phi)/T), with <cos phi> = -0.2633 and a trans fraction
// (|phi| > 120 degrees) of 0.4946. This checks the Metropolis factor together with the Jacobian and
// solution-count ratios (without the Jacobian ratio <cos phi> averaged over the torsions shifts by
// about +0.035, with the interior torsions off by up to +0.06), and the internal-energy bookkeeping
// of the move. Barrier crossings make successive samples strongly correlated, so the checks are on
// means with tolerances calibrated on repeated runs rather than on histogram chi-squared values.
TEST(MC_CONCERTED_ROTATION, trappe_torsions_follow_boltzmann_distribution)
{
  constexpr double c1 = 355.03, c2 = -68.19, c3 = 791.32;
  constexpr double temperature = 600.0;
  constexpr double exactCos = -0.263285;
  constexpr double exactTrans = 0.494568;

  const ForceField forceField = makeZeroForceField();
  Component decane =
      makeComponent(forceField, "conrot-decane-trappe", decaneJson(std::format("0.0, {}, {}, {}", c1, c2, c3)));
  System system = System(forceField, SimulationBox(200.0, 200.0, 200.0), false, temperature, 1e4, 1.0, {}, {decane},
                         {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].concertedRotationRandomizationFraction = 0.2;
  system.components[0].pivotRandomizationFraction = 0.3;

  RandomNumber random(31337);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);

  constexpr std::size_t numberOfMoves = 400000;
  constexpr std::size_t burnIn = 20000;
  constexpr std::size_t sampleEvery = 10;
  constexpr std::size_t bins = 12;
  TorsionStatistics statistics(7, bins);
  std::array<std::size_t, 7> transCounts{};

  std::size_t conrotAttempts = 0, conrotAccepted = 0;
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    // Adapt the small-step channels during burn-in only (the step sizes must be fixed while
    // sampling to keep the chain reversible).
    if (i < burnIn && i % 2000 == 1999) system.components[0].mc_moves_statistics.optimizeMCMoves();

    // The concerted rotation never moves the first and last two beads of a linear chain, so it is
    // combined with pivot moves for ergodicity (as in any production use).
    std::optional<RunningEnergy> energyDifference;
    if (i % 10 == 0)
    {
      energyDifference = MC_Moves::pivotMove(random, system, 0, 0);
    }
    else
    {
      energyDifference = MC_Moves::concertedRotationMove(random, system, 0, 0);
      if (i >= burnIn)
      {
        ++conrotAttempts;
        if (energyDifference.has_value()) ++conrotAccepted;
      }
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnIn || i % sampleEvery != 0) continue;
    statistics.sample(atoms);
    for (std::size_t torsion = 0; torsion != 7; ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      if (std::abs(phi) > 2.0 * std::numbers::pi / 3.0) ++transCounts[torsion];
    }
  }

  EXPECT_GT(conrotAccepted, conrotAttempts / 10);

  double meanCos = 0.0, meanTrans = 0.0;
  for (std::size_t torsion = 0; torsion != 7; ++torsion)
  {
    const double cosine = statistics.sumCos[torsion] / static_cast<double>(statistics.samples);
    const double trans = static_cast<double>(transCounts[torsion]) / static_cast<double>(statistics.samples);
    meanCos += cosine / 7.0;
    meanTrans += trans / 7.0;
    EXPECT_NEAR(cosine, exactCos, 0.05) << "torsion " << torsion;
    EXPECT_NEAR(trans, exactTrans, 0.05) << "torsion " << torsion;
  }
  EXPECT_NEAR(meanCos, exactCos, 0.02);
  EXPECT_NEAR(meanTrans, exactTrans, 0.02);

  // Energy bookkeeping: the accumulated running energy must match a full recomputation.
  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy() * Units::EnergyToKelvin,
              recomputed.potentialEnergy() * Units::EnergyToKelvin, 1e-6);
}

// Branched chain: side groups inside the window are carried rigidly, so every bond length of the
// molecule is preserved exactly, and the running energy stays consistent.
TEST(MC_CONCERTED_ROTATION, branched_chain_preserves_all_bond_lengths)
{
  const ForceField forceField = makeZeroForceField();
  Component branched = makeComponent(forceField, "conrot-branched-run", kBranchedJson);
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 400.0, 1e4, 1.0, {}, {branched}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(99);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  ASSERT_EQ(atoms.size(), 10uz);

  const std::vector<std::pair<std::size_t, std::size_t>> bonds{
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {2, 8}, {5, 9}};
  std::vector<double> initial;
  for (auto [i, j] : bonds) initial.push_back((atoms[i].position - atoms[j].position).length());

  std::size_t accepted = 0;
  double maxBondDeviation = 0.0;
  for (std::size_t move = 0; move != 20000; ++move)
  {
    std::optional<RunningEnergy> energyDifference = MC_Moves::concertedRotationMove(random, system, 0, 0);
    if (energyDifference.has_value())
    {
      ++accepted;
      system.runningEnergies += energyDifference.value();
    }
    for (std::size_t b = 0; b != bonds.size(); ++b)
    {
      auto [i, j] = bonds[b];
      maxBondDeviation =
          std::max(maxBondDeviation, std::abs((atoms[i].position - atoms[j].position).length() - initial[b]));
    }
  }

  EXPECT_GT(accepted, 100uz);
  EXPECT_LT(maxBondDeviation, 1e-9);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy() * Units::EnergyToKelvin,
              recomputed.potentialEnergy() * Units::EnergyToKelvin, 1e-6);
}


