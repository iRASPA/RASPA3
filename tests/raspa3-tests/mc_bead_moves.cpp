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
import move_statistics;
import mc_moves_probabilities;
import mc_moves_bead_displacement;
import mc_moves_bead_flip;

// Tests for the single-bead moves: the bead displacement (random Cartesian displacement of one
// flexible bead) and the bead flip (bond-length-preserving rotation of one bead: kink jump for an
// interior bead, end rotation for a terminal bead). Both are symmetric proposals accepted with the
// plain Metropolis criterion, so the sampled distribution of a chain with simple intramolecular
// potentials and no non-bonded interactions is known exactly:
//  - harmonic bonds, free bends, no torsions: the bond vectors are independent with density
//    proportional to exp(-beta k (l - l0)^2 / 2) d^3l, so the bond length has the marginal
//    l^2 exp(-beta k (l - l0)^2 / 2) and bends and torsions are uniform;
//  - fixed bonds, free bends, no torsions: uniform bond directions, so cos(bend) is uniform on
//    [-1, 1] and the torsions are uniform on [-pi, pi].

namespace
{

constexpr double kBondLength = 1.54;

// Linear 6-bead chain with harmonic bonds of the given force constant (K/A^2), free bends and a
// zero torsional potential. With zero epsilon there is no other energy.
std::string harmonicHexaneJson(double bondForceConstant)
{
  return std::format(R"({{
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH2", [3.08, 0.0, 0.0]],
      ["CH2", [4.62, 0.0, 0.0]], ["CH2", [6.16, 0.0, 0.0]], ["CH2", [7.70, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [{}, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [0.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto"
}}
)",
                     bondForceConstant);
}

// Linear 10-bead chain with FIXED bonds and the given bend type (FIXED or free HARMONIC), zero
// torsional potential.
std::string fixedBondDecaneJson(std::string_view bendDefinition)
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
    [["CH2", "CH2", "CH2"], {}]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
  ],
  "VanDerWaals" : "auto"
}}
)",
                     bendDefinition);
}

// Branched chain with FIXED bonds: CH3 side groups on backbone atoms 2 and 5 of an 8-bead backbone.
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

// Fully flexible hexane with realistic intramolecular potentials, used for the energy bookkeeping.
constexpr std::string_view kTrappeHexaneJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH2", [3.08, 0.0, 0.0]],
      ["CH2", [4.62, 0.0, 0.0]], ["CH2", [6.16, 0.0, 0.0]], ["CH3", [7.70, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]],
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114.0]],
    [["CH2", "CH2", "CH2"], "HARMONIC", [62500.0, 114.0]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

ForceField makeForceField(double epsilonCH2, double epsilonCH3)
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}, {"CH3", false, 15.04, 0.0, 0.0, 6, false}},
                    {{epsilonCH2, 3.95}, {epsilonCH3, 3.75}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                    12.0, true, true, false);
}

Component makeComponent(const ForceField& forceField, std::string name, std::string_view json)
{
  TemporaryFile file(name + ".json", json);
  return Component(Component::Type::Adsorbate, 0, forceField, name, file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false);
}

System makeSingleMoleculeSystem(const ForceField& forceField, const Component& component, double temperature)
{
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, temperature, 1e4, 1.0, {}, {component}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
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

// Running means of the bond lengths, bend cosines and torsion angles of a linear chain.
struct ChainStatistics
{
  std::size_t samples{};
  double sumBond{}, sumBondSquared{};
  double sumCosBend{}, sumCosBendSquared{};
  double sumCosTorsion{}, sumSinTorsion{};
  double maxBondDeviation{};
  std::size_t bonds{}, bends{}, torsions{};

  void sample(std::span<const Atom> atoms)
  {
    ++samples;
    bonds = atoms.size() - 1;
    bends = atoms.size() - 2;
    torsions = atoms.size() - 3;
    for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
    {
      double l = (atoms[i + 1].position - atoms[i].position).length();
      sumBond += l;
      sumBondSquared += l * l;
      maxBondDeviation = std::max(maxBondDeviation, std::abs(l - kBondLength));
    }
    for (std::size_t i = 0; i + 2 < atoms.size(); ++i)
    {
      double c = cosBendAngle(atoms[i].position, atoms[i + 1].position, atoms[i + 2].position);
      sumCosBend += c;
      sumCosBendSquared += c * c;
    }
    for (std::size_t i = 0; i + 3 < atoms.size(); ++i)
    {
      double phi = dihedralAngle(atoms[i].position, atoms[i + 1].position, atoms[i + 2].position, atoms[i + 3].position);
      sumCosTorsion += std::cos(phi);
      sumSinTorsion += std::sin(phi);
    }
  }

  double meanBond() const { return sumBond / static_cast<double>(samples * bonds); }
  double meanBondSquared() const { return sumBondSquared / static_cast<double>(samples * bonds); }
  double meanCosBend() const { return sumCosBend / static_cast<double>(samples * bends); }
  double meanCosBendSquared() const { return sumCosBendSquared / static_cast<double>(samples * bends); }
  double meanCosTorsion() const { return sumCosTorsion / static_cast<double>(samples * torsions); }
  double meanSinTorsion() const { return sumSinTorsion / static_cast<double>(samples * torsions); }
};

// Moments <l^n> of the bond-length marginal l^2 exp(-beta k (l - l0)^2 / 2) on [0, infinity) by
// composite Simpson integration (the integrand is negligible beyond l0 + 10 sigma).
double bondLengthMoment(double betaK, double l0, int n)
{
  const double sigma = 1.0 / std::sqrt(betaK);
  const double upper = l0 + 10.0 * sigma;
  constexpr std::size_t intervals = 200000;
  const double h = upper / static_cast<double>(intervals);
  auto weight = [&](double l) { return l * l * std::exp(-0.5 * betaK * (l - l0) * (l - l0)); };
  double numerator = 0.0, denominator = 0.0;
  for (std::size_t i = 0; i <= intervals; ++i)
  {
    double l = static_cast<double>(i) * h;
    double simpson = (i == 0 || i == intervals) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    double w = simpson * weight(l);
    denominator += w;
    numerator += w * std::pow(l, n);
  }
  return numerator / denominator;
}

}  // namespace

// The lists of movable beads follow from the connectivity and the holonomic constraints:
//  - fully flexible chain: every bead can be displaced, and every bead with one or two neighbours
//    can be flipped;
//  - FIXED bonds: no bead can be displaced (that would change a bond), but all beads can still be
//    flipped when the bends are free;
//  - FIXED bends as well: no bead can be flipped either (a terminal bead sits at the end of a bend,
//    an interior bead is the end atom of the bend centred on its neighbour);
//  - a branch point (three neighbours) cannot be flipped, its side group can.
TEST(MC_BEAD_MOVES, movable_bead_enumeration_respects_constraints)
{
  const ForceField forceField = makeForceField(0.0, 0.0);

  Component flexible = makeComponent(forceField, "bead-hexane-flexible", harmonicHexaneJson(96500.0));
  EXPECT_EQ(flexible.displaceableBeads().size(), 6uz);
  ASSERT_EQ(flexible.flipBeads().size(), 6uz);
  EXPECT_EQ(flexible.flipBeads()[0].neighbours, (std::vector<std::size_t>{1}));
  EXPECT_EQ(flexible.flipBeads()[2].neighbours, (std::vector<std::size_t>{1, 3}));
  EXPECT_EQ(flexible.flipBeads()[5].neighbours, (std::vector<std::size_t>{4}));

  Component fixedBonds =
      makeComponent(forceField, "bead-decane-fixed-bonds", fixedBondDecaneJson("\"HARMONIC\", [0.0, 114.0]"));
  EXPECT_TRUE(fixedBonds.displaceableBeads().empty());
  EXPECT_EQ(fixedBonds.flipBeads().size(), 10uz);

  Component fixedBondsAndBends =
      makeComponent(forceField, "bead-decane-fixed-all", fixedBondDecaneJson("\"FIXED\", [114.0]"));
  EXPECT_TRUE(fixedBondsAndBends.displaceableBeads().empty());
  EXPECT_TRUE(fixedBondsAndBends.flipBeads().empty());

  Component branched = makeComponent(forceField, "bead-branched", kBranchedJson);
  EXPECT_TRUE(branched.displaceableBeads().empty());
  std::vector<std::size_t> flippable;
  for (const Component::FlipBead& bead : branched.flipBeads()) flippable.push_back(bead.atom);
  EXPECT_EQ(flippable, (std::vector<std::size_t>{0, 1, 3, 4, 6, 7, 8, 9}));

  // A component without movable beads rejects the move without counting a trial.
  System system = makeSingleMoleculeSystem(forceField, fixedBondsAndBends, 300.0);
  RandomNumber random(1);
  for (std::size_t i = 0; i != 50; ++i)
  {
    EXPECT_FALSE(MC_Moves::beadDisplacementMove(random, system, 0, 0).has_value());
    EXPECT_FALSE(MC_Moves::beadFlipMove(random, system, 0, 0).has_value());
  }
  for (Move::Types move : {Move::Types::BeadDisplacement, Move::Types::BeadFlip})
  {
    const MoveStatistics<double3>& moveStatistics =
        std::get<MoveStatistics<double3>>(system.components[0].mc_moves_statistics[move]);
    EXPECT_EQ(moveStatistics.counts.x + moveStatistics.counts.y + moveStatistics.counts.z, 0.0);
  }
}

// Bead displacement alone is ergodic for a fully flexible chain. With harmonic bonds, free bends and
// no other interactions the bond length has the exact marginal l^2 exp(-beta k (l - l0)^2 / 2); a
// soft bond (sigma = 0.24 A) makes the l^2 factor shift <l> by 0.08 A above l0, so the test
// distinguishes the correct distribution from a Gaussian in l and from any proposal asymmetry. The
// bends and torsions must be uniform.
TEST(MC_BEAD_MOVES, displacement_samples_harmonic_bond_distribution)
{
  constexpr double bondForceConstant = 5000.0;  // K/A^2
  constexpr double temperature = 300.0;
  const double betaK = bondForceConstant / temperature;

  const ForceField forceField = makeForceField(0.0, 0.0);
  Component hexane = makeComponent(forceField, "bead-hexane-soft", harmonicHexaneJson(bondForceConstant));
  System system = makeSingleMoleculeSystem(forceField, hexane, temperature);

  RandomNumber random(12345);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  ASSERT_EQ(atoms.size(), 6uz);

  constexpr std::size_t numberOfMoves = 600000;
  constexpr std::size_t burnIn = 20000;
  constexpr std::size_t sampleEvery = 10;
  ChainStatistics statistics;
  std::size_t accepted = 0;

  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    if (i < burnIn && i % 2000 == 1999) system.components[0].mc_moves_statistics.optimizeMCMoves();

    std::optional<RunningEnergy> energyDifference = MC_Moves::beadDisplacementMove(random, system, 0, 0);
    if (energyDifference.has_value())
    {
      system.runningEnergies += energyDifference.value();
      if (i >= burnIn) ++accepted;
    }
    if (i >= burnIn && i % sampleEvery == 0) statistics.sample(atoms);
  }

  const double acceptance = static_cast<double>(accepted) / static_cast<double>(numberOfMoves - burnIn);
  EXPECT_GT(acceptance, 0.2);
  EXPECT_LT(acceptance, 0.9);

  const double exactMeanBond = bondLengthMoment(betaK, kBondLength, 1);
  const double exactMeanBondSquared = bondLengthMoment(betaK, kBondLength, 2);
  EXPECT_GT(exactMeanBond, kBondLength + 0.05);  // the l^2 factor matters at this force constant
  EXPECT_NEAR(statistics.meanBond(), exactMeanBond, 0.006);
  EXPECT_NEAR(statistics.meanBondSquared(), exactMeanBondSquared, 0.02);

  EXPECT_NEAR(statistics.meanCosBend(), 0.0, 0.03);
  EXPECT_NEAR(statistics.meanCosBendSquared(), 1.0 / 3.0, 0.02);
  EXPECT_NEAR(statistics.meanCosTorsion(), 0.0, 0.03);
  EXPECT_NEAR(statistics.meanSinTorsion(), 0.0, 0.03);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
}

// Bead flips (kink jumps for the interior beads, end rotations for the terminal beads) preserve
// every bond length exactly. With FIXED bonds, free bends and no other interactions the bond
// directions are independent and uniform, so cos(bend) is uniform on [-1, 1] (mean 0, mean square
// 1/3) and the torsions are uniform on [-pi, pi].
TEST(MC_BEAD_MOVES, flip_preserves_bonds_and_samples_uniform_bends_and_torsions)
{
  const ForceField forceField = makeForceField(0.0, 0.0);
  Component decane =
      makeComponent(forceField, "bead-decane-flip", fixedBondDecaneJson("\"HARMONIC\", [0.0, 114.0]"));
  System system = makeSingleMoleculeSystem(forceField, decane, 300.0);
  system.components[0].beadFlipRandomizationFraction = 0.3;

  RandomNumber random(777);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  ASSERT_EQ(atoms.size(), 10uz);

  constexpr std::size_t numberOfMoves = 400000;
  constexpr std::size_t burnIn = 20000;
  constexpr std::size_t sampleEvery = 10;
  ChainStatistics statistics;
  std::size_t accepted = 0;

  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    if (i < burnIn && i % 2000 == 1999) system.components[0].mc_moves_statistics.optimizeMCMoves();

    std::optional<RunningEnergy> energyDifference = MC_Moves::beadFlipMove(random, system, 0, 0);
    if (energyDifference.has_value())
    {
      system.runningEnergies += energyDifference.value();
      if (i >= burnIn) ++accepted;
    }
    if (i >= burnIn && i % sampleEvery == 0) statistics.sample(atoms);
  }

  // At U = 0 every proposal is accepted.
  EXPECT_EQ(accepted, numberOfMoves - burnIn);

  EXPECT_LT(statistics.maxBondDeviation, 1e-9);
  EXPECT_NEAR(statistics.meanCosBend(), 0.0, 0.03);
  EXPECT_NEAR(statistics.meanCosBendSquared(), 1.0 / 3.0, 0.02);
  EXPECT_NEAR(statistics.meanCosTorsion(), 0.0, 0.03);
  EXPECT_NEAR(statistics.meanSinTorsion(), 0.0, 0.03);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
}

// With realistic bonded potentials and intramolecular Lennard-Jones interactions both moves must
// keep the running energies consistent with a full recomputation, and neither may be trivially
// accepted or rejected.
TEST(MC_BEAD_MOVES, energy_bookkeeping_with_interactions)
{
  const ForceField forceField = makeForceField(46.0, 98.0);
  Component hexane = makeComponent(forceField, "bead-hexane-trappe", kTrappeHexaneJson);
  System system = makeSingleMoleculeSystem(forceField, hexane, 300.0);

  RandomNumber random(99);
  std::size_t displacementAccepted = 0, flipAccepted = 0;
  constexpr std::size_t numberOfMoves = 40000;
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    if (i % 2000 == 1999) system.components[0].mc_moves_statistics.optimizeMCMoves();

    std::optional<RunningEnergy> energyDifference;
    if (i % 2 == 0)
    {
      energyDifference = MC_Moves::beadDisplacementMove(random, system, 0, 0);
      if (energyDifference.has_value()) ++displacementAccepted;
    }
    else
    {
      energyDifference = MC_Moves::beadFlipMove(random, system, 0, 0);
      if (energyDifference.has_value()) ++flipAccepted;
    }
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();
  }

  EXPECT_GT(displacementAccepted, numberOfMoves / 20);
  EXPECT_LT(displacementAccepted, numberOfMoves / 2);
  EXPECT_GT(flipAccepted, numberOfMoves / 20);
  EXPECT_LT(flipAccepted, numberOfMoves / 2);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6);
  EXPECT_NEAR(system.runningEnergies.bond, recomputed.bond, 1e-6);
  EXPECT_NEAR(system.runningEnergies.bend, recomputed.bend, 1e-6);
  EXPECT_NEAR(system.runningEnergies.torsion, recomputed.torsion, 1e-6);
  EXPECT_NEAR(system.runningEnergies.intraVDW, recomputed.intraVDW, 1e-6);
}
