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
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_pivot;
import mc_moves_crankshaft;

// 'U = 0' sanity test for the conformational Monte-Carlo moves (Vitalis & Pappu, Methods 46
// (2009)): with every conformation-dependent potential switched off, the exact answers are known.
// Every pivot/crankshaft proposal must be accepted (the moves preserve all bond lengths, so with
// zero-strength bends/torsions, zero-epsilon Lennard-Jones and no charges the energy difference is
// exactly zero), and the sampled distribution must be the ideal one: with only bond-length
// constraints the configurational measure carries no conformation-dependent Jacobian, so the bond
// direction vectors are i.i.d. uniform on the sphere. Hence every torsion angle is uniform on
// [-pi, pi] and the cosine of every bend angle is uniform on [-1, 1]. Any hidden asymmetry in the
// proposals (biased axis selection, a missing Jacobian, an asymmetric angle window) would show up
// as a deviation from these distributions, independent of any force field.

namespace
{

// Fully flexible 6-bead chain: harmonic bonds (constant under pivot/crankshaft, which preserve all
// bond lengths), zero-strength bends and torsions, and 'auto' intramolecular van-der-Waals pairs
// that vanish because the force field carries epsilon = 0.
constexpr std::string_view kFlexibleHexaneJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.16637443033673, 1.40686000476961, 0.0]],
      ["CH2", [3.70637443033673, 1.40686000476961, 0.0]],
      ["CH2", [4.33274886067346, 2.81372000953922, 0.0]],
      ["CH2", [5.87274886067346, 2.81372000953922, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [0.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto"
}
)";

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}}, {{0.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeFlexibleHexane(const ForceField& forceField)
{
  TemporaryFile file("flexible-hexane-u0.json", kFlexibleHexaneJson);
  return Component(Component::Type::Adsorbate, 0, forceField, "flexible-hexane-u0", file.stemPath().string(), 5, 21,
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

// Cosine of the bend angle a-b-c.
double cosBendAngle(const double3& a, const double3& b, const double3& c)
{
  return double3::dot((a - b).normalized(), (c - b).normalized());
}

// Pearson chi-squared statistic of 'counts' against the uniform distribution.
double chiSquaredUniform(const std::vector<std::size_t>& counts)
{
  std::size_t total = std::accumulate(counts.begin(), counts.end(), 0uz);
  double expected = static_cast<double>(total) / static_cast<double>(counts.size());
  double chiSquared = 0.0;
  for (std::size_t count : counts)
  {
    double deviation = static_cast<double>(count) - expected;
    chiSquared += deviation * deviation / expected;
  }
  return chiSquared;
}

}  // namespace

// Structural check of the precomputed move geometry on the 6-bead chain: three interior pivot
// bonds (terminal bonds have nothing to rotate) and ten crankshaft units (four single-bead, three
// two-bead, two three-bead and one four-bead segment at the default cap of four).
TEST(MC_U0_SANITY, chain_pivot_bonds_and_crankshaft_units)
{
  const ForceField forceField = makeZeroForceField();
  Component hexane = makeFlexibleHexane(forceField);

  EXPECT_EQ(hexane.pivotBonds().size(), 3uz);
  for (const Component::PivotBond& pivotBond : hexane.pivotBonds())
  {
    EXPECT_FALSE(pivotBond.rotatedAtoms.empty());
  }

  EXPECT_EQ(hexane.crankshaftUnits().size(), 10uz);
  for (const Component::CrankshaftUnit& unit : hexane.crankshaftUnits())
  {
    EXPECT_GE(unit.rotatedAtoms.size(), 1uz);
    EXPECT_LE(unit.rotatedAtoms.size(), 4uz);
  }
}

TEST(MC_U0_SANITY, pivot_crankshaft_accept_all_and_sample_ideal_distributions)
{
  const ForceField forceField = makeZeroForceField();
  Component hexane = makeFlexibleHexane(forceField);
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {hexane}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  // Full randomization on both moves gives the fastest decorrelation; at U = 0 the acceptance is
  // independent of the angle window anyway.
  system.components[0].pivotRandomizationFraction = 1.0;
  system.components[0].crankshaftRandomizationFraction = 1.0;

  RandomNumber random(1984);

  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  constexpr std::size_t numberOfAtoms = 6;
  ASSERT_EQ(atoms.size(), numberOfAtoms);

  // The moves must preserve the bond lengths frozen in at the initial CBMC growth.
  std::array<double, numberOfAtoms - 1> initialBondLengths{};
  for (std::size_t i = 0; i + 1 < numberOfAtoms; ++i)
  {
    initialBondLengths[i] = (atoms[i + 1].position - atoms[i].position).length();
  }

  // A bend angle is altered by only three of the ten crankshaft units (rotations that move all
  // three bend atoms rigidly leave it invariant), i.e. roughly every seventh move, so samples are
  // taken sparsely enough to be nearly independent and keep the chi-squared statistic calibrated.
  constexpr std::size_t numberOfMoves = 600000;
  constexpr std::size_t sampleEvery = 12;
  constexpr std::size_t burnInMoves = 3000;
  constexpr std::size_t torsionBins = 12;
  constexpr std::size_t bendBins = 10;

  std::array<std::vector<std::size_t>, 3> torsionCounts;
  torsionCounts.fill(std::vector<std::size_t>(torsionBins, 0uz));
  std::array<std::vector<std::size_t>, 4> cosBendCounts;
  cosBendCounts.fill(std::vector<std::size_t>(bendBins, 0uz));
  std::array<double, 3> sumCosTorsion{};
  std::array<double, 3> sumSinTorsion{};
  std::array<double, 4> sumCosBend{};
  std::size_t samples = 0;
  std::size_t rejected = 0;
  double maxBondDeviation = 0.0;

  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    std::optional<RunningEnergy> energyDifference = (i % 2 == 0)
                                                        ? MC_Moves::pivotMove(random, system, 0, 0)
                                                        : MC_Moves::crankshaftMove(random, system, 0, 0);
    if (!energyDifference.has_value()) ++rejected;

    if (i < burnInMoves || i % sampleEvery != 0) continue;
    ++samples;

    for (std::size_t bond = 0; bond + 1 < numberOfAtoms; ++bond)
    {
      double length = (atoms[bond + 1].position - atoms[bond].position).length();
      maxBondDeviation = std::max(maxBondDeviation, std::abs(length - initialBondLengths[bond]));
    }

    for (std::size_t torsion = 0; torsion != 3; ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      sumCosTorsion[torsion] += std::cos(phi);
      sumSinTorsion[torsion] += std::sin(phi);
      std::size_t bin = std::min(
          torsionBins - 1, static_cast<std::size_t>((phi + std::numbers::pi) / (2.0 * std::numbers::pi) *
                                                    static_cast<double>(torsionBins)));
      ++torsionCounts[torsion][bin];
    }

    for (std::size_t bend = 0; bend != 4; ++bend)
    {
      double cosTheta =
          cosBendAngle(atoms[bend].position, atoms[bend + 1].position, atoms[bend + 2].position);
      sumCosBend[bend] += cosTheta;
      std::size_t bin =
          std::min(bendBins - 1, static_cast<std::size_t>((cosTheta + 1.0) / 2.0 * static_cast<double>(bendBins)));
      ++cosBendCounts[bend][bin];
    }
  }

  // At U = 0 every proposal must be accepted.
  EXPECT_EQ(rejected, 0uz);

  // The moves are pure rotations: bond lengths are preserved to machine precision.
  EXPECT_LT(maxBondDeviation, 1e-9);

  // Every torsion angle must be uniform on [-pi, pi]. The chi-squared statistic against the
  // uniform histogram has 11 degrees of freedom (expectation 11); the generous threshold accounts
  // for the residual correlation between successive samples while still rejecting a systematic
  // bin bias of a few percent.
  for (std::size_t torsion = 0; torsion != 3; ++torsion)
  {
    EXPECT_LT(chiSquaredUniform(torsionCounts[torsion]), 60.0) << "torsion " << torsion;
    EXPECT_NEAR(sumCosTorsion[torsion] / static_cast<double>(samples), 0.0, 0.03) << "torsion " << torsion;
    EXPECT_NEAR(sumSinTorsion[torsion] / static_cast<double>(samples), 0.0, 0.03) << "torsion " << torsion;
  }

  // The cosine of every bend angle must be uniform on [-1, 1] (i.i.d. bond directions on the
  // sphere); this part is sampled by the crankshaft move alone, since the pivot preserves bends.
  for (std::size_t bend = 0; bend != 4; ++bend)
  {
    EXPECT_LT(chiSquaredUniform(cosBendCounts[bend]), 60.0) << "bend " << bend;
    EXPECT_NEAR(sumCosBend[bend] / static_cast<double>(samples), 0.0, 0.03) << "bend " << bend;
  }

  // Standard bookkeeping check: the running energies must match a full recomputation.
  RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
}
