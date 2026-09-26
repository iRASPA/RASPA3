#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
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
import mc_moves_translation;
import mc_moves_double_bridging;
import mc_moves_double_rebridging;

// Tests for the connectivity-altering bridging moves of Karayiannis, Mavrantzas and Theodorou:
// double bridging (DB, two chains exchange tails) and intramolecular double rebridging (IDR, a chain
// segment is reversed). On chains with fixed bonds and bends and no non-bonded interactions the
// torsion angles are independent and uniform, and the mean square end-to-end distance is that of
// the freely rotating chain, which gives exact references that are sensitive to the closure
// Jacobians and the partner-selection factor in the acceptance rules. The energy bookkeeping of the
// tail exchange (pairs changing between intra- and intermolecular, Ewald exclusions) is checked
// against full recomputation in a dense system with Lennard-Jones and charged beads.

namespace
{

constexpr double kBondLength = 1.54;
constexpr double kBendAngleDegrees = 114.0;

// Linear chain of 'beads' united atoms with FIXED bonds and bends. With 'alternating' the beads
// alternate between the pseudo-atoms CHA and CHB (used to give the chain alternating charges).
std::string chainJson(std::size_t beads, std::string_view torsionParameters, bool alternating = false)
{
  std::string pseudoAtoms{};
  std::string connectivity{};
  for (std::size_t i = 0; i != beads; ++i)
  {
    const std::string type = alternating ? (i % 2 == 0 ? "CHA" : "CHB") : "CH2";
    pseudoAtoms +=
        std::format("{}[\"{}\", [{}, 0.0, 0.0]]", i == 0 ? "" : ", ", type, kBondLength * static_cast<double>(i));
    if (i + 1 != beads) connectivity += std::format("{}[{}, {}]", i == 0 ? "" : ", ", i, i + 1);
  }
  std::string bonds, bends, torsions;
  if (alternating)
  {
    bonds = R"([["CHA", "CHB"], "FIXED", [1.54]])";
    bends = R"([["CHA", "CHB", "CHA"], "FIXED", [114.0]], [["CHB", "CHA", "CHB"], "FIXED", [114.0]])";
    torsions = std::format(R"([["CHA", "CHB", "CHA", "CHB"], "TRAPPE", [{}]])", torsionParameters);
  }
  else
  {
    bonds = R"([["CH2", "CH2"], "FIXED", [1.54]])";
    bends = R"([["CH2", "CH2", "CH2"], "FIXED", [114.0]])";
    torsions = std::format(R"([["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [{}]])", torsionParameters);
  }
  return std::format(R"({{
  "CriticalTemperature" : 617.7,
  "CriticalPressure" : 2110000.0,
  "AcentricFactor" : 0.492,
  "pseudoAtoms" : [{}],
  "Connectivity" : [{}],
  "Bonds" : [{}],
  "Bends" : [{}],
  "Torsions" : [{}],
  "VanDerWaals" : "auto"
}}
)",
                     pseudoAtoms, connectivity, bonds, bends, torsions);
}

// Branched chain: a CH3 side group on backbone atoms 2 and 5 of an 8-bead backbone.
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
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CHA", false, 14.03, 0.0, 0.0, 6, false},
                     {"CHB", false, 14.03, 0.0, 0.0, 6, false}},
                    {{0.0, 3.95}, {0.0, 3.75}, {0.0, 3.95}, {0.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot,
                    12.0, 12.0, 12.0, true, true, false);
}

// Small Lennard-Jones beads with alternating charges (the chain is neutral) and Ewald summation.
ForceField makeChargedForceField()
{
  return ForceField({{"CHA", false, 14.03, 0.05, 0.0, 6, false}, {"CHB", false, 14.03, -0.05, 0.0, 6, false}},
                    {{40.0, 2.5}, {40.0, 2.5}}, ForceField::MixingRule::Lorentz_Berthelot, 6.0, 6.0, 6.0, true, false,
                    true);
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

// Largest deviation of the backbone bond lengths and bend cosines of a linear chain from the fixed
// values.
std::pair<double, double> backboneDeviation(std::span<const Atom> atoms)
{
  double bondDeviation = 0.0, bendDeviation = 0.0;
  const double cosFixed = std::cos(kBendAngleDegrees * std::numbers::pi / 180.0);
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    bondDeviation =
        std::max(bondDeviation, std::abs((atoms[i + 1].position - atoms[i].position).length() - kBondLength));
  }
  for (std::size_t i = 0; i + 2 < atoms.size(); ++i)
  {
    bendDeviation =
        std::max(bendDeviation,
                 std::abs(cosBendAngle(atoms[i].position, atoms[i + 1].position, atoms[i + 2].position) - cosFixed));
  }
  return {bondDeviation, bendDeviation};
}

// Torsion-angle statistics collected in blocks. Successive samples of a Metropolis chain are
// strongly correlated (the bridging moves have acceptance ratios of a few percent), so a plain
// Pearson chi-squared against the expected bin probabilities is inflated by an unknown factor. The
// blocked statistic uses the block-to-block variance of every bin fraction instead and is
// approximately chi-squared with (bins - 1) degrees of freedom for an unbiased sampler when the
// blocks are long compared with the correlation time.
struct BlockedTorsionStatistics
{
  std::size_t torsions, bins, blocks;
  std::vector<std::vector<std::vector<std::size_t>>> counts;  // [block][torsion][bin]
  std::vector<std::vector<double>> sumCos, sumSin;            // [block][torsion]
  std::vector<std::size_t> samples;                           // [block]

  BlockedTorsionStatistics(std::size_t torsions_, std::size_t bins_, std::size_t blocks_)
      : torsions(torsions_),
        bins(bins_),
        blocks(blocks_),
        counts(blocks_, std::vector<std::vector<std::size_t>>(torsions_, std::vector<std::size_t>(bins_, 0uz))),
        sumCos(blocks_, std::vector<double>(torsions_, 0.0)),
        sumSin(blocks_, std::vector<double>(torsions_, 0.0)),
        samples(blocks_, 0uz)
  {
  }

  void sample(std::size_t block, std::span<const Atom> atoms)
  {
    ++samples[block];
    for (std::size_t torsion = 0; torsion != torsions; ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      sumCos[block][torsion] += std::cos(phi);
      sumSin[block][torsion] += std::sin(phi);
      std::size_t bin = std::min(
          bins - 1,
          static_cast<std::size_t>((phi + std::numbers::pi) / (2.0 * std::numbers::pi) * static_cast<double>(bins)));
      ++counts[block][torsion][bin];
    }
  }

  std::size_t totalSamples() const { return std::accumulate(samples.begin(), samples.end(), 0uz); }

  // Chi-squared of the mean bin fractions against 'probabilities', with the variance of every bin
  // fraction estimated from the block values.
  double blockedChiSquared(std::size_t torsion, const std::vector<double>& probabilities) const
  {
    double statistic = 0.0;
    for (std::size_t bin = 0; bin != bins; ++bin)
    {
      double mean = 0.0, meanSquare = 0.0;
      for (std::size_t block = 0; block != blocks; ++block)
      {
        double fraction = static_cast<double>(counts[block][torsion][bin]) / static_cast<double>(samples[block]);
        mean += fraction;
        meanSquare += fraction * fraction;
      }
      mean /= static_cast<double>(blocks);
      meanSquare /= static_cast<double>(blocks);
      double variance = (meanSquare - mean * mean) * static_cast<double>(blocks) / static_cast<double>(blocks - 1);
      double deviation = mean - probabilities[bin];
      statistic += deviation * deviation / (variance / static_cast<double>(blocks));
    }
    return statistic;
  }

  double meanCos(std::size_t torsion) const
  {
    double sum = 0.0;
    for (std::size_t block = 0; block != blocks; ++block) sum += sumCos[block][torsion];
    return sum / static_cast<double>(totalSamples());
  }

  double meanSin(std::size_t torsion) const
  {
    double sum = 0.0;
    for (std::size_t block = 0; block != blocks; ++block) sum += sumSin[block][torsion];
    return sum / static_cast<double>(totalSamples());
  }

  // Fraction of samples in bins whose centres satisfy 'predicate' (angle in [-pi, pi]).
  double fraction(std::size_t torsion, std::predicate<double> auto predicate) const
  {
    std::size_t selected = 0;
    for (std::size_t block = 0; block != blocks; ++block)
    {
      for (std::size_t bin = 0; bin != bins; ++bin)
      {
        double centre =
            -std::numbers::pi + (static_cast<double>(bin) + 0.5) * 2.0 * std::numbers::pi / static_cast<double>(bins);
        if (predicate(centre)) selected += counts[block][torsion][bin];
      }
    }
    return static_cast<double>(selected) / static_cast<double>(totalSamples());
  }
};

// Bin probabilities and <cos phi> of a single torsion with the TraPPE potential
// U(phi) = p1 (1 + cos phi) + p2 (1 - cos 2 phi) + p3 (1 + cos 3 phi) at temperature T (midpoint
// quadrature; the integrand is smooth and periodic, so the rule converges exponentially).
struct TorsionReference
{
  std::vector<double> binProbabilities;
  double meanCos;
};

TorsionReference torsionReference(double p1, double p2, double p3, double temperature, std::size_t bins)
{
  constexpr std::size_t points = 36000;
  std::vector<double> weights(bins, 0.0);
  double partition = 0.0, cosine = 0.0;
  for (std::size_t i = 0; i != points; ++i)
  {
    double phi =
        -std::numbers::pi + (static_cast<double>(i) + 0.5) * 2.0 * std::numbers::pi / static_cast<double>(points);
    double energy = p1 * (1.0 + std::cos(phi)) + p2 * (1.0 - std::cos(2.0 * phi)) + p3 * (1.0 + std::cos(3.0 * phi));
    double weight = std::exp(-energy / temperature);
    partition += weight;
    cosine += weight * std::cos(phi);
    weights[std::min(bins - 1, i * bins / points)] += weight;
  }
  for (double& weight : weights) weight /= partition;
  return {weights, cosine / partition};
}

// Mean square end-to-end distance of the freely rotating chain with 'bonds' bonds of length l and
// bend angle theta (uniform torsions): <R^2> = l^2 [n (1+a)/(1-a) - 2a (1-a^n)/(1-a)^2], a = -cos theta.
double freelyRotatingMeanSquareEndToEnd(std::size_t bonds)
{
  const double a = -std::cos(kBendAngleDegrees * std::numbers::pi / 180.0);
  const double n = static_cast<double>(bonds);
  return kBondLength * kBondLength *
         (n * (1.0 + a) / (1.0 - a) - 2.0 * a * (1.0 - std::pow(a, n)) / ((1.0 - a) * (1.0 - a)));
}

double squareEndToEnd(std::span<const Atom> atoms)
{
  const double3 r = atoms.back().position - atoms.front().position;
  return double3::dot(r, r);
}

}  // namespace

// A linear chain of M beads has the sites s = 1 ... M-6 (the trimer s+1 .. s+3 needs the anchors
// s-1, s and s+4, s+5) and the pairs (a, b) with b >= a + 5. Side groups belong to the unit of the
// backbone atom they hang off.
TEST(MC_DOUBLE_BRIDGING, topology_sites_pairs_and_units)
{
  const ForceField forceField = makeZeroForceField();

  Component decane = makeComponent(forceField, "db-decane", chainJson(10, "0.0, 0.0, 0.0, 0.0"));
  const Component::BridgingTopology& decaneTopology = decane.bridgingTopology();
  ASSERT_EQ(decaneTopology.units.size(), 10uz);
  EXPECT_EQ(decaneTopology.sites, (std::vector<std::size_t>{1, 2, 3, 4}));
  EXPECT_TRUE(decaneTopology.sitePairs.empty());
  EXPECT_NEAR(decaneTopology.maximumBridgeDistance, 4.0 * kBondLength, 1e-12);
  for (const Component::BridgingTopology::Unit& unit : decaneTopology.units) EXPECT_TRUE(unit.sideAtoms.empty());

  Component c16 = makeComponent(forceField, "db-c16", chainJson(16, "0.0, 0.0, 0.0, 0.0"));
  const Component::BridgingTopology& c16Topology = c16.bridgingTopology();
  EXPECT_EQ(c16Topology.sites.size(), 10uz);
  std::size_t expectedPairs = 0;
  for (std::size_t a : c16Topology.sites)
    for (std::size_t b : c16Topology.sites)
      if (b >= a + 5) ++expectedPairs;
  EXPECT_EQ(c16Topology.sitePairs.size(), expectedPairs);
  EXPECT_GT(expectedPairs, 0uz);

  // Alternating bead types: a segment of units a+4 .. b is congruent under reversal only when the
  // mirror unit has the same type, i.e. when a + b is even.
  Component alternating = makeComponent(forceField, "db-alt", chainJson(16, "0.0, 0.0, 0.0, 0.0", true));
  for (const std::array<std::size_t, 2>& pair : alternating.bridgingTopology().sitePairs)
  {
    EXPECT_EQ((pair[0] + pair[1]) % 2, 0uz);
  }
  EXPECT_LT(alternating.bridgingTopology().sitePairs.size(), expectedPairs);
  EXPECT_GT(alternating.bridgingTopology().sitePairs.size(), 0uz);

  Component branched = makeComponent(forceField, "db-branched", kBranchedJson);
  const Component::BridgingTopology& branchedTopology = branched.bridgingTopology();
  ASSERT_EQ(branchedTopology.units.size(), 8uz);
  EXPECT_EQ(branchedTopology.units[2].sideAtoms, (std::vector<std::size_t>{8}));
  EXPECT_EQ(branchedTopology.units[5].sideAtoms, (std::vector<std::size_t>{9}));
  EXPECT_EQ(branchedTopology.sites, (std::vector<std::size_t>{1, 2}));
}

// U = 0 sampler test of the intramolecular double rebridging on a 16-bead chain: torsions uniform,
// bonds and bends preserved exactly, <R^2> of the freely rotating chain. Pivot moves (exact at
// U = 0) are mixed in for ergodicity of the chain ends; the rebridging refreshes the torsions about
// three times as often as the pivots. The acceptance ratio of the rebridging is only ~2% (the two
// new closures must both have solutions), so the samples are strongly correlated and the blocked
// chi-squared is used. Dropping the closure Jacobians from the acceptance rule raises the mean
// blocked chi-squared from ~14 to ~80 with these settings; the solution-count ratio hardly affects
// the torsion marginals (as for the concerted rotation).
TEST(MC_DOUBLE_BRIDGING, idr_u0_torsions_uniform_and_constraints_preserved)
{
  const ForceField forceField = makeZeroForceField();
  Component chain = makeComponent(forceField, "idr-c16-u0", chainJson(16, "0.0, 0.0, 0.0, 0.0"));

  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {chain}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotRandomizationFraction = 1.0;

  RandomNumber random(2024);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  ASSERT_EQ(atoms.size(), 16uz);

  constexpr std::size_t numberOfMoves = 1000000;
  constexpr std::size_t burnIn = 10000;
  constexpr std::size_t sampleEvery = 5;
  constexpr std::size_t pivotEvery = 20;
  constexpr std::size_t bins = 12;
  constexpr std::size_t blocks = 20;
  constexpr std::size_t torsions = 13;
  constexpr std::size_t blockLength = (numberOfMoves - burnIn) / blocks;
  BlockedTorsionStatistics statistics(torsions, bins, blocks);
  std::size_t attempts = 0, accepted = 0;
  double maxBond = 0.0, maxBend = 0.0, sumR2 = 0.0;
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    std::optional<RunningEnergy> energyDifference;
    if (i % pivotEvery == 0)
    {
      energyDifference = MC_Moves::pivotMove(random, system, 0, 0);
    }
    else
    {
      ++attempts;
      energyDifference = MC_Moves::intramolecularDoubleRebridgingMove(random, system, 0, 0);
      if (energyDifference.has_value()) ++accepted;
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnIn || i % sampleEvery != 0) continue;
    statistics.sample(std::min(blocks - 1, (i - burnIn) / blockLength), atoms);
    sumR2 += squareEndToEnd(atoms);
    auto [bond, bend] = backboneDeviation(atoms);
    maxBond = std::max(maxBond, bond);
    maxBend = std::max(maxBend, bend);
  }

  EXPECT_GT(accepted, attempts / 100);
  EXPECT_LT(maxBond, 1e-9);
  EXPECT_LT(maxBend, 1e-9);

  const std::vector<double> uniform(bins, 1.0 / static_cast<double>(bins));
  double meanChiSquared = 0.0;
  for (std::size_t torsion = 0; torsion != torsions; ++torsion)
  {
    const double chiSquared = statistics.blockedChiSquared(torsion, uniform);
    meanChiSquared += chiSquared / static_cast<double>(torsions);
    EXPECT_LT(chiSquared, 60.0) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanCos(torsion), 0.0, 0.04) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanSin(torsion), 0.0, 0.04) << "torsion " << torsion;
  }
  EXPECT_LT(meanChiSquared, 30.0);

  const double meanR2 = sumR2 / static_cast<double>(statistics.totalSamples());
  EXPECT_NEAR(meanR2, freelyRotatingMeanSquareEndToEnd(15), 0.04 * freelyRotatingMeanSquareEndToEnd(15));

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
}

// Boltzmann sampler test of the rebridging with a torsion potential U = p1 (1 + cos phi), p1 = T/2
// (a single well at trans, <cos phi> = -I1(1/2)/I0(1/2) = -0.2425): every torsion is independently
// Boltzmann distributed. The full TraPPE potential at 600 K gives an acceptance ratio below 0.1%
// (six to twelve torsions are re-drawn at once), which leaves too few accepted moves for a test; the
// weaker potential keeps ~1%. The pivots alone sample the same distribution, so the test is
// sensitive to the energy difference entering the rebridging acceptance only in proportion to the
// fraction of torsion refreshes it performs (roughly one half here): an IDR ignoring the torsion
// energy shifts <cos phi> by about +0.1.
TEST(MC_DOUBLE_BRIDGING, idr_torsion_potential_follows_boltzmann_distribution)
{
  constexpr double temperature = 600.0;
  constexpr double p1 = 300.0;
  constexpr std::size_t bins = 12;
  const TorsionReference reference = torsionReference(p1, 0.0, 0.0, temperature, bins);
  ASSERT_NEAR(reference.meanCos, -0.242504, 1e-5);

  const ForceField forceField = makeZeroForceField();
  Component chain = makeComponent(forceField, "idr-c16-torsion", chainJson(16, std::format("0.0, {}, 0.0, 0.0", p1)));
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, temperature, 1e4, 1.0, {}, {chain}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotRandomizationFraction = 1.0;

  RandomNumber random(31337);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);

  constexpr std::size_t numberOfMoves = 1000000;
  constexpr std::size_t burnIn = 10000;
  constexpr std::size_t sampleEvery = 5;
  constexpr std::size_t pivotEvery = 20;
  constexpr std::size_t blocks = 20;
  constexpr std::size_t torsions = 13;
  constexpr std::size_t blockLength = (numberOfMoves - burnIn) / blocks;
  BlockedTorsionStatistics statistics(torsions, bins, blocks);
  std::size_t attempts = 0, accepted = 0;
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    std::optional<RunningEnergy> energyDifference;
    if (i % pivotEvery == 0)
    {
      energyDifference = MC_Moves::pivotMove(random, system, 0, 0);
    }
    else
    {
      ++attempts;
      energyDifference = MC_Moves::intramolecularDoubleRebridgingMove(random, system, 0, 0);
      if (energyDifference.has_value()) ++accepted;
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnIn || i % sampleEvery != 0) continue;
    statistics.sample(std::min(blocks - 1, (i - burnIn) / blockLength), atoms);
  }

  EXPECT_GT(accepted, attempts / 400);

  double meanCos = 0.0, meanChiSquared = 0.0;
  for (std::size_t torsion = 0; torsion != torsions; ++torsion)
  {
    const double chiSquared = statistics.blockedChiSquared(torsion, reference.binProbabilities);
    meanChiSquared += chiSquared / static_cast<double>(torsions);
    meanCos += statistics.meanCos(torsion) / static_cast<double>(torsions);
    EXPECT_LT(chiSquared, 60.0) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanCos(torsion), reference.meanCos, 0.04) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanSin(torsion), 0.0, 0.04) << "torsion " << torsion;
  }
  EXPECT_LT(meanChiSquared, 30.0);
  EXPECT_NEAR(meanCos, reference.meanCos, 0.015);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy() * Units::EnergyToKelvin,
              recomputed.potentialEnergy() * Units::EnergyToKelvin, 1e-6);
}

// U = 0 sampler test of the double bridging: six decane chains in a small periodic box (so that
// bridging partners exist) without any non-bonded interaction. Every accepted move must change
// exactly two chains, all bonds and bends stay fixed, the torsions of every chain are uniform and
// <R^2> is that of the freely rotating chain. Pivot and translation moves (exact at U = 0) provide
// ergodicity of the chain ends and of the positions; the double bridging dominates.
TEST(MC_DOUBLE_BRIDGING, db_u0_torsions_uniform_and_tails_exchanged)
{
  const ForceField forceField = makeZeroForceField();
  Component decane = makeComponent(forceField, "db-decane-u0", chainJson(10, "0.0, 0.0, 0.0, 0.0"));

  constexpr std::size_t numberOfChains = 6;
  System system = System(forceField, SimulationBox(14.0, 14.0, 14.0), false, 300.0, 1e4, 1.0, {}, {decane}, {},
                         {numberOfChains}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotRandomizationFraction = 1.0;

  RandomNumber random(4711);

  constexpr std::size_t numberOfMoves = 600000;
  constexpr std::size_t burnIn = 10000;
  constexpr std::size_t sampleEvery = 10;
  constexpr std::size_t bins = 12;
  constexpr std::size_t blocks = 20;
  constexpr std::size_t torsions = 7;
  constexpr std::size_t blockLength = (numberOfMoves - burnIn) / blocks;
  BlockedTorsionStatistics statistics(torsions, bins, blocks);
  std::size_t attempts = 0, accepted = 0;
  double maxBond = 0.0, maxBend = 0.0, sumR2 = 0.0;
  std::size_t r2Samples = 0;
  std::vector<Atom> snapshot{};
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    const std::size_t molecule = static_cast<std::size_t>(random.uniform() * static_cast<double>(numberOfChains));
    std::optional<RunningEnergy> energyDifference;
    if (i % 20 == 0)
    {
      energyDifference = MC_Moves::pivotMove(random, system, 0, molecule);
    }
    else if (i % 20 == 10)
    {
      energyDifference = MC_Moves::translationMove(random, system, 0, molecule);
    }
    else
    {
      ++attempts;
      const std::span<const Atom> all = system.spanOfMoleculeAtoms();
      snapshot.assign(all.begin(), all.end());
      energyDifference = MC_Moves::doubleBridgingMove(random, system, 0, molecule);
      if (energyDifference.has_value())
      {
        ++accepted;
        std::size_t changedChains = 0;
        for (std::size_t m = 0; m != numberOfChains; ++m)
        {
          std::span<const Atom> atoms = system.spanOfMolecule(0, m);
          bool changed = false;
          for (std::size_t k = 0; k != atoms.size(); ++k)
          {
            if ((atoms[k].position - snapshot[m * 10 + k].position).length() > 1e-12) changed = true;
            EXPECT_EQ(atoms[k].moleculeId, snapshot[m * 10 + k].moleculeId);
            EXPECT_EQ(atoms[k].type, snapshot[m * 10 + k].type);
          }
          if (changed) ++changedChains;
        }
        EXPECT_EQ(changedChains, 2uz);
      }
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnIn || i % sampleEvery != 0) continue;
    const std::size_t block = std::min(blocks - 1, (i - burnIn) / blockLength);
    for (std::size_t m = 0; m != numberOfChains; ++m)
    {
      std::span<const Atom> atoms = system.spanOfMolecule(0, m);
      statistics.sample(block, atoms);
      sumR2 += squareEndToEnd(atoms);
      ++r2Samples;
      auto [bond, bend] = backboneDeviation(atoms);
      maxBond = std::max(maxBond, bond);
      maxBend = std::max(maxBend, bend);
    }
  }

  EXPECT_GT(accepted, attempts / 100);
  EXPECT_LT(maxBond, 1e-9);
  EXPECT_LT(maxBend, 1e-9);

  const std::vector<double> uniform(bins, 1.0 / static_cast<double>(bins));
  double meanChiSquared = 0.0;
  for (std::size_t torsion = 0; torsion != torsions; ++torsion)
  {
    const double chiSquared = statistics.blockedChiSquared(torsion, uniform);
    meanChiSquared += chiSquared / static_cast<double>(torsions);
    EXPECT_LT(chiSquared, 60.0) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanCos(torsion), 0.0, 0.03) << "torsion " << torsion;
    EXPECT_NEAR(statistics.meanSin(torsion), 0.0, 0.03) << "torsion " << torsion;
  }
  EXPECT_LT(meanChiSquared, 30.0);
  const double meanR2 = sumR2 / static_cast<double>(r2Samples);
  EXPECT_NEAR(meanR2, freelyRotatingMeanSquareEndToEnd(9), 0.03 * freelyRotatingMeanSquareEndToEnd(9));

  system.checkMoleculeIds();
  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
}

// Energy bookkeeping of both moves in a system of twelve interacting chains with Lennard-Jones
// beads carrying alternating charges (Ewald summation) and a weak torsion potential: the tail
// exchange turns intermolecular pairs into intramolecular ones and vice versa (real-space Coulomb
// exclusions, Ewald exclusion terms, intramolecular van der Waals), so the accumulated running
// energy is compared term by term with a full recomputation. The beads are small (sigma = 2.5 A) so
// that the unbiased trimer placement of the bridging moves is accepted often enough at this
// density; with alkane-sized beads the acceptance ratio in a melt is well below 0.1%.
TEST(MC_DOUBLE_BRIDGING, energy_bookkeeping_with_lennard_jones_and_ewald)
{
  const ForceField forceField = makeChargedForceField();
  Component chain = makeComponent(forceField, "db-c14-charged", chainJson(14, "0.0, 100.0, 0.0, 0.0", true));

  constexpr std::size_t numberOfChains = 12;
  System system =
      System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 600.0, 1e4, 1.0, {}, {chain}, {}, {numberOfChains}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotRandomizationFraction = 0.5;

  RandomNumber random(8080);
  std::size_t acceptedBridging = 0, acceptedRebridging = 0;
  for (std::size_t i = 0; i != 60000; ++i)
  {
    const std::size_t molecule = static_cast<std::size_t>(random.uniform() * static_cast<double>(numberOfChains));
    std::optional<RunningEnergy> energyDifference;
    switch (i % 4)
    {
      case 0:
        energyDifference = MC_Moves::pivotMove(random, system, 0, molecule);
        break;
      case 1:
        energyDifference = MC_Moves::translationMove(random, system, 0, molecule);
        break;
      case 2:
        energyDifference = MC_Moves::doubleBridgingMove(random, system, 0, molecule);
        if (energyDifference.has_value()) ++acceptedBridging;
        break;
      default:
        energyDifference = MC_Moves::intramolecularDoubleRebridgingMove(random, system, 0, molecule);
        if (energyDifference.has_value()) ++acceptedRebridging;
        break;
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();
  }

  EXPECT_GT(acceptedBridging, 20uz);
  EXPECT_GT(acceptedRebridging, 20uz);

  system.checkMoleculeIds();
  const RunningEnergy running = system.runningEnergies;
  const RunningEnergy recomputed = system.computeTotalEnergies();
  const double toKelvin = Units::EnergyToKelvin;
  EXPECT_NEAR(running.moleculeMoleculeVDW * toKelvin, recomputed.moleculeMoleculeVDW * toKelvin, 1e-4);
  EXPECT_NEAR(running.moleculeMoleculeCharge * toKelvin, recomputed.moleculeMoleculeCharge * toKelvin, 1e-4);
  EXPECT_NEAR(running.ewald_fourier * toKelvin, recomputed.ewald_fourier * toKelvin, 1e-4);
  EXPECT_NEAR(running.ewald_exclusion * toKelvin, recomputed.ewald_exclusion * toKelvin, 1e-4);
  EXPECT_NEAR(running.intraVDW * toKelvin, recomputed.intraVDW * toKelvin, 1e-4);
  EXPECT_NEAR(running.torsion * toKelvin, recomputed.torsion * toKelvin, 1e-4);
  EXPECT_NEAR(running.potentialEnergy() * toKelvin, recomputed.potentialEnergy() * toKelvin, 1e-4);
}
