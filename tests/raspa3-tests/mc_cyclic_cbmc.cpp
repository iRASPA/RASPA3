#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import double3;
import atom;
import molecule;
import forcefield;
import component;
import system;
import monte_carlo;
import molecular_dynamics;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import move_statistics;
import fragment;
import fragment_graph;
import cbmc_growth_plan;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import mc_moves_reinsertion;
import intra_molecular_potentials;

// Tests for flexible rings grown with ring-closure CBMC. Rings are inferred from cycles in the
// molecule's connectivity graph. The reference molecules are united-atom cyclohexane (six cyclic
// CH2 beads, one simple ring) and trans-decalin (two fused six-membered rings).

namespace
{

ForceField makeCyclohexaneForceField()
{
  return ForceField({{"CH2_c", false, 14.02658, 0.0, 0.0, 6, false}}, {{52.5, 3.91}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeCyclohexane(const ForceField& forceField, std::size_t componentId,
                          const MCMoveProbabilities& probabilities)
{
  TemporaryFile file("cyclohexane.json", molecule_fixtures::kCyclohexaneJson);
  return Component(Component::Type::Adsorbate, componentId, forceField, "cyclohexane", file.stemPath().string(), 5, 21,
                   probabilities, std::nullopt, false);
}

constexpr double kRingBondLength = 1.54;

// Every consecutive pair around the ring (including the closure bond 5-0) must stay close to the
// stiff harmonic equilibrium length; the ring must remain closed throughout the sampling.
void expectRingClosedAndReasonable(const System& system)
{
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  constexpr std::size_t moleculeSize = 6;
  ASSERT_EQ(atoms.size() % moleculeSize, 0uz);

  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    for (std::size_t i = 0; i != 6; ++i)
    {
      double3 p0 = atoms[offset + i].position;
      double3 p1 = atoms[offset + (i + 1) % 6].position;
      EXPECT_NEAR((p1 - p0).length(), kRingBondLength, 0.4);
    }

    // The ring bends stay in a physically reasonable window around their harmonic equilibrium.
    for (std::size_t i = 0; i != 6; ++i)
    {
      double angle = double3::angle(atoms[offset + (i + 5) % 6].position, atoms[offset + i].position,
                                    atoms[offset + (i + 1) % 6].position);
      EXPECT_GT(angle, 90.0 * (std::numbers::pi / 180.0));
      EXPECT_LT(angle, 135.0 * (std::numbers::pi / 180.0));
    }
  }
}

void expectNoEnergyDrift(System& system)
{
  RunningEnergy recomputedEnergies = system.computeTotalEnergies();
  RunningEnergy drift = system.runningEnergies - recomputedEnergies;

  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
  EXPECT_NEAR(drift.intraVDW, 0.0, 1e-6);
  EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
  EXPECT_NEAR(drift.tail, 0.0, 1e-6);
}

// Mean ring bond length and bend angle over all molecules; used to cross-check the ring geometry
// distribution sampled by CBMC against the one sampled by molecular dynamics.
std::pair<double, double> meanRingGeometry(const System& system)
{
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  constexpr std::size_t moleculeSize = 6;

  double bondSum{};
  double bendSum{};
  std::size_t count{};
  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    for (std::size_t i = 0; i != 6; ++i)
    {
      double3 pm = atoms[offset + (i + 5) % 6].position;
      double3 p0 = atoms[offset + i].position;
      double3 pp = atoms[offset + (i + 1) % 6].position;
      bondSum += (pp - p0).length();
      bendSum += double3::angle(pm, p0, pp);
      ++count;
    }
  }
  return {bondSum / static_cast<double>(count), bendSum / static_cast<double>(count)};
}

constexpr std::array<std::array<std::size_t, 2>, 11> kDecalinBonds{
    {{0, 1}, {0, 5}, {0, 9}, {1, 2}, {1, 6}, {2, 3}, {3, 4}, {4, 5}, {6, 7}, {7, 8}, {8, 9}}};

constexpr std::array<std::array<std::size_t, 2>, 8> kNorbornaneBonds{
    {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 6}, {3, 6}}};

constexpr std::array<std::array<std::size_t, 2>, 7> kMethylcyclohexaneBonds{
    {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 6}}};

// Every listed bond of every molecule must stay close to the stiff harmonic equilibrium length:
// all rings of the molecule (simple, fused, or bridged) remain closed throughout the sampling.
template <std::size_t N>
void expectBondsIntact(const System& system, std::size_t moleculeSize,
                       const std::array<std::array<std::size_t, 2>, N>& bonds)
{
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  ASSERT_EQ(atoms.size() % moleculeSize, 0uz);

  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    for (const std::array<std::size_t, 2>& bond : bonds)
    {
      double3 p0 = atoms[offset + bond[0]].position;
      double3 p1 = atoms[offset + bond[1]].position;
      EXPECT_NEAR((p1 - p0).length(), kRingBondLength, 0.4);
    }
  }
}

void expectDecalinRingsClosed(const System& system) { expectBondsIntact(system, 10, kDecalinBonds); }

// Mean bond length and bend angle over all listed bonds and all bonded triples derived from them;
// used to cross-check the geometry distribution sampled by CBMC against molecular dynamics.
template <std::size_t N>
std::pair<double, double> meanMoleculeGeometry(const System& system, std::size_t moleculeSize,
                                               const std::array<std::array<std::size_t, 2>, N>& bonds)
{
  // All bonded triples i-j-k (i < k) with j bonded to both.
  std::vector<std::array<std::size_t, 3>> bends{};
  for (std::size_t j = 0; j != moleculeSize; ++j)
  {
    std::vector<std::size_t> neighbors{};
    for (const std::array<std::size_t, 2>& bond : bonds)
    {
      if (bond[0] == j) neighbors.push_back(bond[1]);
      if (bond[1] == j) neighbors.push_back(bond[0]);
    }
    for (std::size_t a = 0; a != neighbors.size(); ++a)
    {
      for (std::size_t b = a + 1; b != neighbors.size(); ++b)
      {
        bends.push_back({neighbors[a], j, neighbors[b]});
      }
    }
  }

  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  double bondSum{};
  double bendSum{};
  std::size_t bondCount{};
  std::size_t bendCount{};
  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    for (const std::array<std::size_t, 2>& bond : bonds)
    {
      bondSum += (atoms[offset + bond[1]].position - atoms[offset + bond[0]].position).length();
      ++bondCount;
    }
    for (const std::array<std::size_t, 3>& bend : bends)
    {
      bendSum += double3::angle(atoms[offset + bend[0]].position, atoms[offset + bend[1]].position,
                                atoms[offset + bend[2]].position);
      ++bendCount;
    }
  }
  return {bondSum / static_cast<double>(bondCount), bendSum / static_cast<double>(bendCount)};
}

ForceField makeMethylcyclohexaneForceField()
{
  return ForceField({{"CH2_c", false, 14.02658, 0.0, 0.0, 6, false},
                     {"CH_c", false, 13.01864, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.03452, 0.0, 0.0, 6, false}},
                    {{52.5, 3.91}, {10.0, 4.65}, {98.0, 3.75}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                    12.0, true, true, false);
}

MonteCarlo makeShortMonteCarlo(std::vector<System> systems)
{
  return MonteCarlo({20, 0, 5, 5, 1000, 10000, 5000, 5000}, std::move(systems), 42uz, 5, false);
}

MonteCarlo makeMonteCarlo(std::vector<System> systems, std::size_t production)
{
  return MonteCarlo({production, 0, 200, 200, 100000, 1000000, 5000, 5000}, std::move(systems), 42uz, 5, false);
}

MolecularDynamics makeMolecularDynamics(std::vector<System> systems, std::size_t production)
{
  return MolecularDynamics({production, 0, 200, 200, 100000, 1000000, 5000, 5000}, std::move(systems), 42uz, 5, false);
}

}  // namespace

TEST(MC_CYCLIC_CBMC, fragment_graph_parsing)
{
  const ForceField forceField = makeCyclohexaneForceField();
  Component cyclohexane = makeCyclohexane(forceField, 0, MCMoveProbabilities());

  const FragmentGraph& graph = cyclohexane.fragmentGraph;

  // Six single-atom flexible fragments; the ring is inferred from the connectivity as one cyclic
  // cluster covering all beads, with one closure bond (a simple cycle has |E| - |V| + 1 = 1).
  ASSERT_EQ(graph.fragments.size(), 6uz);
  for (std::size_t bead = 0; bead != 6; ++bead)
  {
    EXPECT_FALSE(graph.fragments[graph.fragmentContaining(bead)].isRigidBody());
  }

  ASSERT_EQ(graph.cyclicClusters.size(), 1uz);
  EXPECT_EQ(graph.cyclicClusters[0], (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));
  EXPECT_EQ(graph.closureBonds.size(), 1uz);

  // A flexible ring keeps all intra-ring interactions: six bonds, six bends and six torsions (one
  // per edge / apex / central bond of the ring).
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.bonds.size(), 6uz);
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.bends.size(), 6uz);
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.torsions.size(), 6uz);

  // A cyclic component with no rigid bodies is fully flexible (integrated per atom).
  EXPECT_FALSE(cyclohexane.isSemiFlexible());
  EXPECT_EQ(graph.numberOfRigidFragments(), 0uz);
}

TEST(MC_CYCLIC_CBMC, growth_plan)
{
  const ForceField forceField = makeCyclohexaneForceField();
  Component cyclohexane = makeCyclohexane(forceField, 0, MCMoveProbabilities());

  std::vector<CBMC::GrowStep> plan = CBMC::buildGrowthPlan(cyclohexane, {cyclohexane.startingBead});

  // The whole ring is grown as one ring-closure step hinged on the seed bead: the remaining five
  // beads are grown together, with the closure bond keeping the ring closed.
  ASSERT_EQ(plan.size(), 1uz);

  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::CloseRing);
  EXPECT_FALSE(plan[0].rigidBody);
  EXPECT_EQ(plan[0].currentBead, cyclohexane.startingBead);
  EXPECT_FALSE(plan[0].previousBead.has_value());
  EXPECT_EQ(plan[0].nextBeads, (std::vector<std::size_t>{1, 2, 3, 4, 5}));

  // The ring-closure step carries all intra-ring bond/bend/torsion potentials (including the
  // closure bond 5-0 and the closure bends/torsions), so the ring is closed during growth.
  EXPECT_EQ(plan[0].intra.bonds.size(), 6uz);
  EXPECT_EQ(plan[0].intra.bends.size(), 6uz);
  EXPECT_EQ(plan[0].intra.torsions.size(), 6uz);
}

// Fused rings are valid input: the whole fused system is one cyclic cluster with two closure bonds,
// grown as a single ring-closure step.
TEST(MC_CYCLIC_CBMC, fused_ring_fragment_graph_and_growth_plan)
{
  const ForceField forceField = makeCyclohexaneForceField();
  TemporaryFile file("decalin.json", molecule_fixtures::kDecalinJson);
  Component decalin = Component(Component::Type::Adsorbate, 0, forceField, "decalin", file.stemPath().string(), 5, 21,
                                MCMoveProbabilities(), std::nullopt, false);

  const FragmentGraph& graph = decalin.fragmentGraph;
  ASSERT_EQ(graph.fragments.size(), 10uz);
  ASSERT_EQ(graph.cyclicClusters.size(), 1uz);
  EXPECT_EQ(graph.cyclicClusters[0], (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
  // Two independent cycles: |E| - |V| + 1 = 11 - 10 + 1 = 2 closure bonds.
  EXPECT_EQ(graph.closureBonds.size(), 2uz);

  std::vector<CBMC::GrowStep> plan = CBMC::buildGrowthPlan(decalin, {decalin.startingBead});
  ASSERT_EQ(plan.size(), 1uz);
  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::CloseRing);
  EXPECT_EQ(plan[0].nextBeads.size(), 9uz);

  // All eleven bonds (nine tree bonds plus two closure bonds) are sampled during ring closure.
  EXPECT_EQ(plan[0].intra.bonds.size(), 11uz);

  // The minimal fused system (two three-membered rings sharing an edge) was a parse error under the
  // old 'Cycle' group model; under the fragment graph it is a valid molecule.
  TemporaryFile fusedFile("fused_ring.json", molecule_fixtures::kFusedRingJson);
  Component fused = Component(Component::Type::Adsorbate, 0, forceField, "fused_ring", fusedFile.stemPath().string(),
                              5, 21, MCMoveProbabilities(), std::nullopt, false);
  ASSERT_EQ(fused.fragmentGraph.cyclicClusters.size(), 1uz);
  EXPECT_EQ(fused.fragmentGraph.cyclicClusters[0], (std::vector<std::size_t>{0, 1, 2, 3}));
  EXPECT_EQ(fused.fragmentGraph.closureBonds.size(), 2uz);
}

// Fused-ring growth: NVT sampling of trans-decalin with CBMC reinsertion. The running energies must
// stay consistent with a full recomputation and both fused rings must remain closed.
TEST(MC_CYCLIC_CBMC, fused_ring_nvt_drift_and_geometry)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  TemporaryFile file("decalin.json", molecule_fixtures::kDecalinJson);
  Component decalin = Component(Component::Type::Adsorbate, 0, forceField, "decalin", file.stemPath().string(), 5, 21,
                                probabilities, std::nullopt, false);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {decalin}, {}, {10}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectDecalinRingsClosed(s);
  }
}

TEST(MC_CYCLIC_CBMC, nvt_drift_and_ring_geometry_cbmc)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component cyclohexane = makeCyclohexane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {cyclohexane}, {}, {20}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingClosedAndReasonable(s);
  }
}

TEST(MC_CYCLIC_CBMC, nvt_drift_and_ring_geometry_recoil_growth)
{
  ForceField forceField = makeCyclohexaneForceField();
  forceField.useRecoilGrowth = true;

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component cyclohexane = makeCyclohexane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {cyclohexane}, {}, {20}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingClosedAndReasonable(s);
  }
}

// Grand-canonical swap (CBMC insertion/deletion) of whole flexible rings: the running energies must
// stay consistent with a full recomputation (grow and retrace weights are mirrored) and both an
// insertion and a deletion must be accepted so that both directions of the ring growth are exercised.
TEST(MC_CYCLIC_CBMC, muvt_swap_cbmc_drift)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);

  Component cyclohexane = makeCyclohexane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {cyclohexane}, {}, {20}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingClosedAndReasonable(s);
  }
}

// Widom test-particle insertion of a flexible ring: the ghost ring is grown with ring-closure CBMC
// but never inserted, so the state is untouched and the sampled Rosenbluth weight is finite/positive.
TEST(MC_CYCLIC_CBMC, widom_rosenbluth_weight)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::Widom, 1.0);

  Component cyclohexane = makeCyclohexane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {cyclohexane}, {}, {20}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    const double rosenbluthWeight = s.components[0].averageRosenbluthWeights.averagedRosenbluthWeight();
    EXPECT_TRUE(std::isfinite(rosenbluthWeight));
    EXPECT_GT(rosenbluthWeight, 0.0);

    expectNoEnergyDrift(s);
    expectRingClosedAndReasonable(s);
  }
}

TEST(MC_CYCLIC_CBMC, nve_ring_geometry_molecular_dynamics)
{
  const ForceField forceField = makeCyclohexaneForceField();

  Component cyclohexane = makeCyclohexane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {cyclohexane}, {}, {10}, 5);
  system.timeStep = 1e-4;

  MolecularDynamics md = MolecularDynamics({20, 0, 5, 5, 1000, 10000, 5000, 5000}, {std::move(system)}, 42uz, 5, false);
  md.run();

  for (System& s : md.systems)
  {
    expectRingClosedAndReasonable(s);
    EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
  }
}

// Bridged rings are valid input: norbornane (bicyclo[2.2.1]heptane) is one cyclic cluster covering
// all seven beads with two closure bonds, grown as a single ring-closure step.
TEST(MC_CYCLIC_CBMC, norbornane_fragment_graph_and_growth_plan)
{
  const ForceField forceField = makeCyclohexaneForceField();
  TemporaryFile file("norbornane.json", molecule_fixtures::kNorbornaneJson);
  Component norbornane = Component(Component::Type::Adsorbate, 0, forceField, "norbornane", file.stemPath().string(),
                                   5, 21, MCMoveProbabilities(), std::nullopt, false);

  const FragmentGraph& graph = norbornane.fragmentGraph;
  ASSERT_EQ(graph.fragments.size(), 7uz);
  ASSERT_EQ(graph.cyclicClusters.size(), 1uz);
  EXPECT_EQ(graph.cyclicClusters[0], (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6}));
  // Two independent cycles: |E| - |V| + 1 = 8 - 7 + 1 = 2 closure bonds.
  EXPECT_EQ(graph.closureBonds.size(), 2uz);

  std::vector<CBMC::GrowStep> plan = CBMC::buildGrowthPlan(norbornane, {norbornane.startingBead});
  ASSERT_EQ(plan.size(), 1uz);
  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::CloseRing);
  EXPECT_EQ(plan[0].nextBeads.size(), 6uz);

  // All eight bonds (six tree bonds plus two closure bonds) are sampled during ring closure.
  EXPECT_EQ(plan[0].intra.bonds.size(), 8uz);
}

// Bridged-ring growth: NVT sampling of norbornane with CBMC reinsertion. The running energies must
// stay consistent with a full recomputation and both bridged rings must remain closed.
TEST(MC_CYCLIC_CBMC, bridged_ring_nvt_drift_and_geometry)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  TemporaryFile file("norbornane.json", molecule_fixtures::kNorbornaneJson);
  Component norbornane = Component(Component::Type::Adsorbate, 0, forceField, "norbornane", file.stemPath().string(),
                                   5, 21, probabilities, std::nullopt, false);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {norbornane}, {}, {10}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectBondsIntact(s, 7, kNorbornaneBonds);
  }
}

// A ring with a flexible tail: the ring is grown first as a ring-closure step from the seed bead;
// the tail is then grown off the closed ring as an ordinary flexible attachment whose anchor (the
// ring atom) has two placed ring neighbors.
TEST(MC_CYCLIC_CBMC, ring_with_tail_growth_plan)
{
  const ForceField forceField = makeMethylcyclohexaneForceField();
  TemporaryFile file("methylcyclohexane.json", molecule_fixtures::kMethylcyclohexaneJson);
  Component methylcyclohexane = Component(Component::Type::Adsorbate, 0, forceField, "methylcyclohexane",
                                          file.stemPath().string(), 5, 21, MCMoveProbabilities(), std::nullopt, false);

  const FragmentGraph& graph = methylcyclohexane.fragmentGraph;
  ASSERT_EQ(graph.fragments.size(), 7uz);
  ASSERT_EQ(graph.cyclicClusters.size(), 1uz);
  EXPECT_EQ(graph.cyclicClusters[0], (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));
  EXPECT_EQ(graph.closureBonds.size(), 1uz);

  std::vector<CBMC::GrowStep> plan = CBMC::buildGrowthPlan(methylcyclohexane, {methylcyclohexane.startingBead});
  ASSERT_EQ(plan.size(), 2uz);

  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::CloseRing);
  EXPECT_EQ(plan[0].currentBead, 0uz);
  EXPECT_EQ(plan[0].nextBeads, (std::vector<std::size_t>{1, 2, 3, 4, 5}));
  EXPECT_EQ(plan[0].intra.bonds.size(), 6uz);

  EXPECT_EQ(plan[1].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[1].currentBead, 0uz);
  EXPECT_EQ(plan[1].previousBead, std::optional<std::size_t>{1});
  EXPECT_EQ(plan[1].nextBeads, (std::vector<std::size_t>{6}));
}

// The central consistency check: a flexible ring is a genuine flexible molecule, so the internal
// ring geometry (bonds and bends) sampled by ring-closure CBMC reinsertion must agree with the
// geometry sampled by molecular dynamics of the same molecule.
TEST(MC_CYCLIC_CBMC, ring_geometry_matches_molecular_dynamics)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component cyclohexaneMC = makeCyclohexane(forceField, 0, probabilities);
  System mcSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {cyclohexaneMC}, {}, {8}, 5);
  MonteCarlo mc = makeMonteCarlo({mcSystem}, 400);
  mc.run();

  Component cyclohexaneMD = makeCyclohexane(forceField, 0, MCMoveProbabilities());
  System mdSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {cyclohexaneMD}, {}, {8}, 5);
  mdSystem.timeStep = 5e-4;
  MolecularDynamics md = makeMolecularDynamics({mdSystem}, 2000);
  md.run();

  auto [mcBond, mcBend] = meanRingGeometry(mc.systems[0]);
  auto [mdBond, mdBend] = meanRingGeometry(md.systems[0]);

  // The instantaneous end-of-run geometry is only a single sample of a stiff distribution, so the
  // tolerances are generous; the point is that CBMC does not systematically distort the ring.
  EXPECT_NEAR(mcBond, mdBond, 0.05);
  EXPECT_NEAR(mcBond, kRingBondLength, 0.1);
  EXPECT_NEAR(mcBend, mdBend, 8.0 * (std::numbers::pi / 180.0));
}

// Grow/retrace weight symmetry (ring-closure Jacobian sanity): the average Rosenbluth weight of
// freshly grown ghost rings must match the average retrace weight of Boltzmann-sampled ring
// configurations — both estimate the same ideal-gas Rosenbluth weight of the molecule. A
// bookkeeping asymmetry between the grow and retrace branches of the ring-closure operator (e.g. a
// bias factor entering only one direction) would systematically shift the two estimators apart.
TEST(MC_CYCLIC_CBMC, grow_retrace_weight_symmetry)
{
  const ForceField forceField = makeCyclohexaneForceField();

  // One molecule in a huge box: external (inter-molecular) energies are negligible, so both
  // averages reduce to the internal-sampling normalization of the ring-closure growth.
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  Component cyclohexane = makeCyclohexane(forceField, 0, probabilities);
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {cyclohexane}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(12);
  double growSum{};
  double retraceSum{};
  std::size_t growCount{};
  std::size_t retraceCount{};

  constexpr std::size_t samples = 400;
  for (std::size_t i = 0; i != samples; ++i)
  {
    // Advance the Boltzmann sampler: regrow the molecule with CBMC reinsertion.
    (void)MC_Moves::reinsertionMove(random, system, 0, 0);

    const CBMC::GrowContext context{system.hasExternalField,
                                    system.forceField,
                                    system.simulationBox,
                                    system.interpolationGrids,
                                    system.externalFieldInterpolationGrid,
                                    system.framework,
                                    system.spanOfFrameworkAtoms(),
                                    system.spanOfMoleculeAtoms(),
                                    system.beta,
                                    system.forceField.cutOffFrameworkVDW,
                                    system.forceField.cutOffMoleculeVDW,
                                    system.forceField.cutOffCoulomb};

    // Grow a ghost ring (never inserted).
    std::optional<ChainGrowData> growData =
        CBMC::growMoleculeSwapInsertion(random, context, system.components[0], 0, 1, 1.0, std::uint8_t{0}, false);
    if (growData.has_value())
    {
      growSum += growData->RosenbluthWeight;
      ++growCount;
    }

    // Retrace the current (Boltzmann-sampled) configuration.
    ChainRetraceData retraceData =
        CBMC::retraceMoleculeSwapDeletion(random, context, system.components[0], system.spanOfMolecule(0, 0));
    retraceSum += retraceData.RosenbluthWeight;
    ++retraceCount;
  }

  ASSERT_GT(growCount, samples / 2);
  const double growMean = growSum / static_cast<double>(growCount);
  const double retraceMean = retraceSum / static_cast<double>(retraceCount);

  EXPECT_TRUE(std::isfinite(growMean));
  EXPECT_TRUE(std::isfinite(retraceMean));
  EXPECT_GT(growMean, 0.0);
  EXPECT_GT(retraceMean, 0.0);
  EXPECT_NEAR(growMean / retraceMean, 1.0, 0.15);
}

// Same CBMC-vs-MD ensemble comparison for a fused-ring molecule: trans-decalin regrown with
// ring-closure CBMC must sample the same internal geometry distribution as molecular dynamics.
TEST(MC_CYCLIC_CBMC, fused_ring_geometry_matches_molecular_dynamics)
{
  const ForceField forceField = makeCyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  TemporaryFile mcFile("decalin.json", molecule_fixtures::kDecalinJson);
  Component decalinMC = Component(Component::Type::Adsorbate, 0, forceField, "decalin", mcFile.stemPath().string(), 5,
                                  21, probabilities, std::nullopt, false);
  System mcSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {decalinMC}, {}, {8}, 5);
  MonteCarlo mc = makeMonteCarlo({mcSystem}, 400);
  mc.run();

  TemporaryFile mdFile("decalin.json", molecule_fixtures::kDecalinJson);
  Component decalinMD = Component(Component::Type::Adsorbate, 0, forceField, "decalin", mdFile.stemPath().string(), 5,
                                  21, MCMoveProbabilities(), std::nullopt, false);
  System mdSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {decalinMD}, {}, {8}, 5);
  mdSystem.timeStep = 5e-4;
  MolecularDynamics md = makeMolecularDynamics({mdSystem}, 2000);
  md.run();

  expectDecalinRingsClosed(mc.systems[0]);
  expectDecalinRingsClosed(md.systems[0]);

  auto [mcBond, mcBend] = meanMoleculeGeometry(mc.systems[0], 10, kDecalinBonds);
  auto [mdBond, mdBend] = meanMoleculeGeometry(md.systems[0], 10, kDecalinBonds);

  EXPECT_NEAR(mcBond, mdBond, 0.05);
  EXPECT_NEAR(mcBond, kRingBondLength, 0.1);
  EXPECT_NEAR(mcBend, mdBend, 8.0 * (std::numbers::pi / 180.0));
}

// Same CBMC-vs-MD ensemble comparison for a ring with a flexible tail: methylcyclohexane regrown
// with ring-closure CBMC plus a flexible tail attachment must sample the same internal geometry
// distribution as molecular dynamics.
TEST(MC_CYCLIC_CBMC, ring_with_tail_geometry_matches_molecular_dynamics)
{
  const ForceField forceField = makeMethylcyclohexaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  TemporaryFile mcFile("methylcyclohexane.json", molecule_fixtures::kMethylcyclohexaneJson);
  Component moleculeMC = Component(Component::Type::Adsorbate, 0, forceField, "methylcyclohexane",
                                   mcFile.stemPath().string(), 5, 21, probabilities, std::nullopt, false);
  System mcSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {moleculeMC}, {}, {8}, 5);
  MonteCarlo mc = makeMonteCarlo({mcSystem}, 400);
  mc.run();

  TemporaryFile mdFile("methylcyclohexane.json", molecule_fixtures::kMethylcyclohexaneJson);
  Component moleculeMD = Component(Component::Type::Adsorbate, 0, forceField, "methylcyclohexane",
                                   mdFile.stemPath().string(), 5, 21, MCMoveProbabilities(), std::nullopt, false);
  System mdSystem =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e4, 1.0, {}, {moleculeMD}, {}, {8}, 5);
  mdSystem.timeStep = 5e-4;
  MolecularDynamics md = makeMolecularDynamics({mdSystem}, 2000);
  md.run();

  expectBondsIntact(mc.systems[0], 7, kMethylcyclohexaneBonds);
  expectBondsIntact(md.systems[0], 7, kMethylcyclohexaneBonds);

  auto [mcBond, mcBend] = meanMoleculeGeometry(mc.systems[0], 7, kMethylcyclohexaneBonds);
  auto [mdBond, mdBend] = meanMoleculeGeometry(md.systems[0], 7, kMethylcyclohexaneBonds);

  EXPECT_NEAR(mcBond, mdBond, 0.05);
  EXPECT_NEAR(mcBond, kRingBondLength, 0.1);
  EXPECT_NEAR(mcBend, mdBend, 8.0 * (std::numbers::pi / 180.0));
}
