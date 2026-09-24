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
import mc_moves_reptation;

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

// The same chain declared as two three-bead repeat units, making it eligible for reptation.
constexpr std::string_view kPeriodicHexaneJson =
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
  "VanDerWaals" : "auto",
  "RepeatUnits" : [
    [0, 1, 2],
    [3, 4, 5]
  ]
}
)";

// A CH3-capped chain is chemically NOT shift-periodic (the end slots differ in type from the
// interior), so declaring repeat units on it must be rejected at parse time.
constexpr std::string_view kCappedChainJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.16637443033673, 1.40686000476961, 0.0]],
      ["CH2", [3.70637443033673, 1.40686000476961, 0.0]],
      ["CH2", [4.33274886067346, 2.81372000953922, 0.0]],
      ["CH3", [5.87274886067346, 2.81372000953922, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]],
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [0.0, 114.0]],
    [["CH2", "CH2", "CH2"], "HARMONIC", [0.0, 114.0]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]],
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto",
  "RepeatUnits" : [
    [0, 1, 2],
    [3, 4, 5]
  ]
}
)";

// A periodic chain with a direction-asymmetric repeat unit (a CH3 side group on the second
// backbone bead): types, connectivity and potentials are shift-periodic, but the unit grows with a
// branch step from one chain end and single-bead steps from the other. With the exact flexible-bead
// base sampler this is valid for reptation (no rigid/ring steps involved).
constexpr std::string_view kBranchedPolymerJson =
R"({
  "CriticalTemperature" : 600.0,
  "CriticalPressure" : 3000000.0,
  "AcentricFactor" : 0.3,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],  ["CH2", [1.54, 0.0, 0.0]],  ["CH3", [2.17, 1.41, 0.0]],
      ["CH2", [3.08, 0.0, 0.0]], ["CH2", [4.62, 0.0, 0.0]],  ["CH3", [5.25, 1.41, 0.0]],
      ["CH2", [6.16, 0.0, 0.0]], ["CH2", [7.70, 0.0, 0.0]],  ["CH3", [8.33, 1.41, 0.0]]
    ],
  "Connectivity" : [
    [0, 1], [1, 2], [1, 3],
    [3, 4], [4, 5], [4, 6],
    [6, 7], [7, 8]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]],
    [["CH2", "CH3"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH3"], "HARMONIC", [0.0, 112.0]],
    [["CH2", "CH2", "CH2"], "HARMONIC", [0.0, 112.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]],
    [["CH3", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto",
  "RepeatUnits" : [
    [0, 1, 2],
    [3, 4, 5],
    [6, 7, 8]
  ]
}
)";

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}, {"CH3", false, 15.04, 0.0, 0.0, 6, false}},
                    {{0.0, 3.95}, {0.0, 3.75}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
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
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

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

// Repeat-unit validation: a homopolymer chain declared as equal monomer blocks passes the
// shift-periodicity checks; the same chain with CH3 end caps is chemically not shift-periodic
// (the end slots differ in type from the interior) and must be rejected at parse time.
TEST(MC_U0_SANITY, repeat_units_validation)
{
  const ForceField forceField = makeZeroForceField();

  {
    TemporaryFile file("periodic-hexane-u0.json", kPeriodicHexaneJson);
    Component periodic = Component(Component::Type::Adsorbate, 0, forceField, "periodic-hexane-u0",
                                   file.stemPath().string(), 5, 21, MCMoveProbabilities(), std::nullopt, false);
    EXPECT_EQ(periodic.repeatUnits.size(), 2uz);
  }

  {
    TemporaryFile file("capped-chain-u0.json", kCappedChainJson);
    EXPECT_THROW(Component(Component::Type::Adsorbate, 0, forceField, "capped-chain-u0", file.stemPath().string(), 5,
                           21, MCMoveProbabilities(), std::nullopt, false),
                 std::runtime_error);
  }

  {
    // Direction-asymmetric repeat unit (branch step from one chain end, single-bead steps from the
    // other): valid for reptation, because purely flexible acyclic steps draw their trial
    // conformations from the exact bonded-Boltzmann base sampler, which makes the cross-plan
    // grow/retrace pairing of reptation exact for arbitrary plans.
    TemporaryFile file("branched-polymer-u0.json", kBranchedPolymerJson);
    Component component(Component::Type::Adsorbate, 0, forceField, "branched-polymer-u0", file.stemPath().string(), 5,
                        21, MCMoveProbabilities(), std::nullopt, false);
    EXPECT_EQ(component.repeatUnits.size(), 3);
  }
}

// U = 0 sanity test for the reptation move alone. Every reptation regrows one repeat unit at a
// chain end with CBMC; at U = 0 every trial direction carries unit weight, so the Rosenbluth ratio
// is exactly one and every proposal must be accepted. Unlike the pivot and crankshaft moves,
// reptation resamples the bond lengths of the regrown unit (from the bond potential), the bend
// angles (uniform in the cosine at zero bend strength) and the torsions (uniform), so the ideal
// distributions must emerge from reptation alone.
TEST(MC_U0_SANITY, reptation_accepts_all_and_samples_ideal_distributions)
{
  const ForceField forceField = makeZeroForceField();
  TemporaryFile file("periodic-hexane-u0.json", kPeriodicHexaneJson);
  Component hexane = Component(Component::Type::Adsorbate, 0, forceField, "periodic-hexane-u0",
                               file.stemPath().string(), 5, 21, MCMoveProbabilities(), std::nullopt, false);
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {hexane}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(1868);

  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  constexpr std::size_t numberOfAtoms = 6;
  ASSERT_EQ(atoms.size(), numberOfAtoms);

  std::array<double, numberOfAtoms - 1> initialBondLengths{};
  for (std::size_t i = 0; i + 1 < numberOfAtoms; ++i)
  {
    initialBondLengths[i] = (atoms[i + 1].position - atoms[i].position).length();
  }

  constexpr std::size_t numberOfMoves = 40000;
  constexpr std::size_t sampleEvery = 2;
  constexpr std::size_t burnInMoves = 1000;
  constexpr std::size_t torsionBins = 12;
  constexpr std::size_t bendBins = 10;

  std::array<std::vector<std::size_t>, 3> torsionCounts;
  torsionCounts.fill(std::vector<std::size_t>(torsionBins, 0uz));
  std::array<std::vector<std::size_t>, 4> cosBendCounts;
  cosBendCounts.fill(std::vector<std::size_t>(bendBins, 0uz));
  std::size_t samples = 0;
  std::size_t rejected = 0;
  double maxBondDeviation = 0.0;
  double sumBondLength = 0.0;
  std::size_t bondSamples = 0;

  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    std::optional<RunningEnergy> energyDifference = MC_Moves::reptationMove(random, system, 0, 0);
    if (!energyDifference.has_value()) ++rejected;
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnInMoves || i % sampleEvery != 0) continue;
    ++samples;

    for (std::size_t bond = 0; bond + 1 < numberOfAtoms; ++bond)
    {
      double length = (atoms[bond + 1].position - atoms[bond].position).length();
      maxBondDeviation = std::max(maxBondDeviation, std::abs(length - initialBondLengths[bond]));
      sumBondLength += length;
      ++bondSamples;
    }

    for (std::size_t torsion = 0; torsion != 3; ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      std::size_t bin = std::min(
          torsionBins - 1, static_cast<std::size_t>((phi + std::numbers::pi) / (2.0 * std::numbers::pi) *
                                                    static_cast<double>(torsionBins)));
      ++torsionCounts[torsion][bin];
    }

    for (std::size_t bend = 0; bend != 4; ++bend)
    {
      double cosTheta = cosBendAngle(atoms[bend].position, atoms[bend + 1].position, atoms[bend + 2].position);
      std::size_t bin =
          std::min(bendBins - 1, static_cast<std::size_t>((cosTheta + 1.0) / 2.0 * static_cast<double>(bendBins)));
      ++cosBendCounts[bend][bin];
    }
  }

  // At U = 0 the Rosenbluth ratio is exactly one: every proposal must be accepted.
  EXPECT_EQ(rejected, 0uz);

  // Reptation must resample the bond lengths (thermal width of the harmonic bond at 300 K is about
  // 0.06 Angstrom), in contrast to the pivot and crankshaft moves which preserve them exactly.
  EXPECT_GT(maxBondDeviation, 0.01);
  EXPECT_NEAR(sumBondLength / static_cast<double>(bondSamples), 1.54, 0.02);

  for (std::size_t torsion = 0; torsion != 3; ++torsion)
  {
    EXPECT_LT(chiSquaredUniform(torsionCounts[torsion]), 60.0) << "torsion " << torsion;
  }
  for (std::size_t bend = 0; bend != 4; ++bend)
  {
    EXPECT_LT(chiSquaredUniform(cosBendCounts[bend]), 60.0) << "bend " << bend;
  }

  // Standard bookkeeping check: the running energies must match a full recomputation.
  RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
}
