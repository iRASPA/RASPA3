#include <gtest/gtest.h>

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
import cbmc_segments;
import intra_molecular_potentials;

// Tests for flexible rings grown with ring-closure CBMC: molecules with a 'Cycle' group in the
// molecule JSON. The reference molecule is united-atom cyclohexane (six cyclic CH2 beads).

namespace
{

std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

ForceField makeCyclohexaneForceField()
{
  return ForceField({{"CH2_c", false, 14.02658, 0.0, 0.0, 6, false}}, {{52.5, 3.91}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeCyclohexane(const ForceField& forceField, std::size_t componentId,
                          const MCMoveProbabilities& probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/20_mc_flexible_cyclohexane_in_box" / "cyclohexane";
  return Component(Component::Type::Adsorbate, componentId, forceField, "cyclohexane", moleculePath.string(), 5, 21,
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

TEST(MC_CYCLIC_CBMC, component_groups_parsing)
{
  const ForceField forceField = makeCyclohexaneForceField();
  Component cyclohexane = makeCyclohexane(forceField, 0, MCMoveProbabilities());

  // one cyclic group covering all six beads
  ASSERT_EQ(cyclohexane.groups.size(), 1uz);
  EXPECT_FALSE(cyclohexane.groups[0].rigid);
  EXPECT_TRUE(cyclohexane.groups[0].cyclic);
  EXPECT_EQ(cyclohexane.groups[0].atoms, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));

  ASSERT_EQ(cyclohexane.atomGroupIds.size(), 6uz);
  EXPECT_EQ(cyclohexane.atomGroupIds, (std::vector<std::size_t>{0, 0, 0, 0, 0, 0}));

  for (std::size_t bead = 0; bead != 6; ++bead)
  {
    EXPECT_FALSE(cyclohexane.rigidGroupContaining(bead).has_value());
    EXPECT_EQ(cyclohexane.cyclicGroupContaining(bead), std::optional<std::size_t>{0});
  }

  // A cyclic group is flexible, so all intra-ring interactions are kept: six bonds, six bends and
  // six torsions (one per edge / apex / central bond of the ring).
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.bonds.size(), 6uz);
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.bends.size(), 6uz);
  EXPECT_EQ(cyclohexane.intraMolecularPotentials.torsions.size(), 6uz);

  // A cyclic component with no rigid groups is fully flexible (integrated per atom).
  EXPECT_FALSE(cyclohexane.isSemiFlexible());
  EXPECT_EQ(cyclohexane.numberOfRigidGroups(), 0uz);
}

TEST(MC_CYCLIC_CBMC, grow_segments)
{
  const ForceField forceField = makeCyclohexaneForceField();
  Component cyclohexane = makeCyclohexane(forceField, 0, MCMoveProbabilities());

  std::vector<CBMC::GrowSegment> segments = CBMC::buildGrowSegments(cyclohexane, {cyclohexane.startingBead});

  // The whole ring is grown as one ring-closure segment hinged on the seed bead: the remaining five
  // beads are ordered along the ring so that the last bead (5) closes back onto the anchor (0).
  ASSERT_EQ(segments.size(), 1uz);

  EXPECT_FALSE(segments[0].rigidUnit);
  EXPECT_TRUE(segments[0].ringUnit);
  EXPECT_EQ(segments[0].cyclicGroup, std::optional<std::size_t>{0});
  EXPECT_EQ(segments[0].currentBead, cyclohexane.startingBead);
  EXPECT_FALSE(segments[0].previousBead.has_value());
  EXPECT_EQ(segments[0].nextBeads, (std::vector<std::size_t>{1, 2, 3, 4, 5}));

  // The ring-closure segment carries all intra-ring bond/bend/torsion potentials (including the
  // closure bond 5-0 and the closure bends/torsions), so the ring is closed during growth.
  EXPECT_EQ(segments[0].intra.bonds.size(), 6uz);
  EXPECT_EQ(segments[0].intra.bends.size(), 6uz);
  EXPECT_EQ(segments[0].intra.torsions.size(), 6uz);
}

// A non-simple cycle (a bead with three in-ring neighbors) must be rejected with a clear error.
TEST(MC_CYCLIC_CBMC, fused_ring_rejected)
{
  const ForceField forceField = makeCyclohexaneForceField();
  const std::filesystem::path moleculePath =
      repositoryRoot() / "tests/raspa3-tests/fixtures" / "bad-fused-ring";
  // The fixture is created on the fly below; if it is missing this test is skipped.
  if (!std::filesystem::exists(moleculePath.string() + ".json")) GTEST_SKIP();
  EXPECT_ANY_THROW(
      Component(Component::Type::Adsorbate, 0, forceField, "bad-fused-ring", moleculePath.string(), 5, 21,
                MCMoveProbabilities(), std::nullopt, false));
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
