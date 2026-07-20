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
import simulation_schedule;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import move_statistics;
import property_widom;
import fragment;
import fragment_graph;
import cbmc_growth_plan;
import generalized_hessian;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_generalized_coordinates;
import minimization_options;
import minimization;
import thermostat;
import thermobarostat;
import intra_molecular_potentials;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bond_torsion_potential;
import bend_torsion_potential;
import inversion_bend_potential;

// Tests for semi-flexible molecules: molecules mixing rigid bodies ('RigidBodies' in the molecule
// JSON) with flexible beads, grown with CBMC or recoil growth.

namespace
{

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

// Semi-flexible pentane: flexible CH3 ends (beads 0 and 4), rigid CH2-CH2-CH2 core (beads 1,2,3)
// with bond length 1.54 and a fixed 114-degree bend.
Component makeSemiFlexiblePentane(const ForceField& forceField, std::size_t componentId,
                                  const MCMoveProbabilities& probabilities)
{
  TemporaryFile file("semi-flexible-pentane.json", molecule_fixtures::kSemiFlexiblePentaneJson);
  return Component(Component::Type::Adsorbate, componentId, forceField, "semi-flexible-pentane",
                   file.stemPath().string(), 5, 21, probabilities, std::nullopt, false);
}

// United-atom 4,4'-diethyl-biphenyl: flexible ethyl tails (beads 0-1 and 14-15), two rigid aromatic
// 6-rings (beads 2-7 and 8-13) connected by a flexible biphenyl bond (5-8).
ForceField makeBiphenylForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false},
                     {"C_ar", false, 12.011, 0.0, 0.0, 6, false},
                     {"CH_ar", false, 13.019, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}, {21.0, 3.88}, {50.5, 3.695}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeDiethylBiphenyl(const ForceField& forceField, std::size_t componentId,
                              const MCMoveProbabilities& probabilities)
{
  TemporaryFile file("diethyl-biphenyl.json", molecule_fixtures::kDiethylBiphenylJson);
  return Component(Component::Type::Adsorbate, componentId, forceField, "diethyl-biphenyl", file.stemPath().string(), 5,
                   21, probabilities, std::nullopt, false);
}

// Same diethyl-biphenyl fixture, including the four 'Partial-reinsertion' fixed-atom sets that
// exercise group-aware partial reinsertion (regrow one tail, the other, or half the molecule).
Component makeDiethylBiphenylFromBasicExample(const ForceField& forceField, std::size_t componentId,
                                              const MCMoveProbabilities& probabilities)
{
  return makeDiethylBiphenyl(forceField, componentId, probabilities);
}

// Fixed distances inside a rigid aromatic ring (regular hexagon, side 1.40).
constexpr double kAromaticBondLength = 1.40;
const double kAromaticMetaDistance = 1.40 * std::sqrt(3.0);
constexpr double kAromaticParaDistance = 2.80;

void expectRingGeometryPreserved(const System& system)
{
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  constexpr std::size_t moleculeSize = 16;
  ASSERT_EQ(atoms.size() % moleculeSize, 0uz);

  constexpr std::array<std::size_t, 6> ringA{2, 3, 4, 5, 6, 7};
  constexpr std::array<std::size_t, 6> ringB{8, 9, 10, 11, 12, 13};

  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    for (const std::array<std::size_t, 6>& ring : {ringA, ringB})
    {
      for (std::size_t i = 0; i != 6; ++i)
      {
        double3 p0 = atoms[offset + ring[i]].position;
        double3 p1 = atoms[offset + ring[(i + 1) % 6]].position;
        double3 p2 = atoms[offset + ring[(i + 2) % 6]].position;
        double3 p3 = atoms[offset + ring[(i + 3) % 6]].position;

        EXPECT_NEAR((p1 - p0).length(), kAromaticBondLength, 1e-6);
        EXPECT_NEAR((p2 - p0).length(), kAromaticMetaDistance, 1e-6);
        EXPECT_NEAR((p3 - p0).length(), kAromaticParaDistance, 1e-6);
      }
    }

    // The flexible junction bonds are sampled from stiff harmonic potentials; they must stay close
    // to their equilibrium lengths.
    EXPECT_NEAR((atoms[offset + 2].position - atoms[offset + 1].position).length(), 1.51, 0.5);
    EXPECT_NEAR((atoms[offset + 8].position - atoms[offset + 5].position).length(), 1.48, 0.5);
    EXPECT_NEAR((atoms[offset + 14].position - atoms[offset + 11].position).length(), 1.51, 0.5);

    // Both ipso bends on each side must be sampled (stiff harmonic around 120 degrees). This guards
    // against under-relaxed rigid-body orientations: the pair of ipso bends forces the exocyclic
    // bond towards the ring plane, identically for the CH2 grown before ring A and the CH2 grown
    // after ring B.
    auto bendAngle = [&](std::size_t a, std::size_t apex, std::size_t c)
    { return double3::angle(atoms[offset + a].position, atoms[offset + apex].position, atoms[offset + c].position); };

    constexpr double kEquilibriumBend = 120.0 * (std::numbers::pi / 180.0);
    constexpr double kBendTolerance = 30.0 * (std::numbers::pi / 180.0);
    EXPECT_NEAR(bendAngle(1, 2, 3), kEquilibriumBend, kBendTolerance);
    EXPECT_NEAR(bendAngle(1, 2, 7), kEquilibriumBend, kBendTolerance);
    EXPECT_NEAR(bendAngle(10, 11, 14), kEquilibriumBend, kBendTolerance);
    EXPECT_NEAR(bendAngle(12, 11, 14), kEquilibriumBend, kBendTolerance);
  }
}

// Fixed distances inside the rigid core (from the reference geometry in the molecule JSON).
constexpr double kRigidBondLength = 1.54;
const double kRigid13Distance =
    std::sqrt(2.16637443033673 * 2.16637443033673 + 1.40686000476961 * 1.40686000476961);

void expectRigidCoreGeometryPreserved(const System& system)
{
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  constexpr std::size_t moleculeSize = 5;
  ASSERT_EQ(atoms.size() % moleculeSize, 0uz);

  for (std::size_t offset = 0; offset < atoms.size(); offset += moleculeSize)
  {
    double3 p1 = atoms[offset + 1].position;
    double3 p2 = atoms[offset + 2].position;
    double3 p3 = atoms[offset + 3].position;

    EXPECT_NEAR((p2 - p1).length(), kRigidBondLength, 1e-6);
    EXPECT_NEAR((p3 - p2).length(), kRigidBondLength, 1e-6);
    EXPECT_NEAR((p3 - p1).length(), kRigid13Distance, 1e-5);

    // The flexible CH3-CH2 junction bonds are sampled from a stiff harmonic potential
    // (96500 K/A^2 at 573 K, sigma ~ 0.08 A); they must stay close to 1.54 A.
    double3 p0 = atoms[offset + 0].position;
    double3 p4 = atoms[offset + 4].position;
    EXPECT_NEAR((p1 - p0).length(), kRigidBondLength, 0.5);
    EXPECT_NEAR((p4 - p3).length(), kRigidBondLength, 0.5);
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

MonteCarlo makeShortMonteCarlo(std::vector<System> systems)
{
  std::size_t numberOfProductionCycles{20};
  std::size_t numberOfInitializationCycles{5};
  std::size_t numberOfEquilibrationCycles{5};
  std::size_t printEvery{1000};
  std::size_t writeBinaryRestartEvery{10000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t numberOfBlocks{5};
  bool outputToFiles{false};

  return MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles,
                     printEvery, writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery},
                    std::move(systems), 42uz, numberOfBlocks, outputToFiles);
}

MolecularDynamics makeShortMolecularDynamics(std::vector<System> systems)
{
  std::size_t numberOfProductionCycles{20};
  std::size_t numberOfInitializationCycles{5};
  std::size_t numberOfEquilibrationCycles{5};
  std::size_t printEvery{1000};
  std::size_t writeBinaryRestartEvery{10000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t numberOfBlocks{5};
  bool outputToFiles{false};

  return MolecularDynamics({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles,
                            printEvery, writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery},
                           std::move(systems), 42uz, numberOfBlocks, outputToFiles);
}

}  // namespace

TEST(MC_SEMI_FLEXIBLE_CBMC, component_fragment_graph_parsing)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  const FragmentGraph& graph = pentane.fragmentGraph;

  // three fragments (ordered by lowest atom index): flexible {0}, rigid {1,2,3}, flexible {4}
  ASSERT_EQ(graph.fragments.size(), 3uz);
  EXPECT_FALSE(graph.fragments[0].isRigidBody());
  EXPECT_TRUE(graph.fragments[1].isRigidBody());
  EXPECT_FALSE(graph.fragments[2].isRigidBody());
  EXPECT_EQ(graph.fragments[1].atoms, (std::vector<std::size_t>{1, 2, 3}));

  ASSERT_EQ(graph.atomFragmentIds.size(), 5uz);
  EXPECT_EQ(graph.atomFragmentIds, (std::vector<std::size_t>{0, 1, 1, 1, 2}));

  EXPECT_FALSE(pentane.rigidFragmentContaining(0).has_value());
  EXPECT_EQ(pentane.rigidFragmentContaining(1), std::optional<std::size_t>{1});
  EXPECT_EQ(pentane.rigidFragmentContaining(2), std::optional<std::size_t>{1});
  EXPECT_EQ(pentane.rigidFragmentContaining(3), std::optional<std::size_t>{1});
  EXPECT_FALSE(pentane.rigidFragmentContaining(4).has_value());

  // interactions entirely inside the rigid body are dropped: only the two junction bonds
  // (0-1, 3-4), the two junction bends (0-1-2, 2-3-4), and the two junction torsions remain.
  EXPECT_EQ(pentane.intraMolecularPotentials.bonds.size(), 2uz);
  EXPECT_EQ(pentane.intraMolecularPotentials.bends.size(), 2uz);
  EXPECT_EQ(pentane.intraMolecularPotentials.torsions.size(), 2uz);
}

TEST(MC_SEMI_FLEXIBLE_CBMC, growth_plan)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  std::vector<CBMC::GrowStep> plan = pentane.growthPlan({pentane.startingBead});

  // Starting from flexible bead 0, the rigid core is entered from outside, so its connecting atom
  // (bead 1) is peeled off as an ordinary flexible bond first; then the rest of the rigid core
  // {2,3} is grown as one rigid body hinged on bead 1; finally the other CH3 end (bead 4) is grown.
  ASSERT_EQ(plan.size(), 3uz);

  EXPECT_FALSE(plan[0].rigidBody);
  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::PlaceSeedFragment);
  EXPECT_EQ(plan[0].currentBead, 0uz);
  EXPECT_FALSE(plan[0].previousBead.has_value());
  EXPECT_EQ(plan[0].nextBeads, (std::vector<std::size_t>{1}));

  EXPECT_TRUE(plan[1].rigidBody);
  EXPECT_EQ(plan[1].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[1].currentBead, 1uz);
  EXPECT_EQ(plan[1].previousBead, std::optional<std::size_t>{0});
  EXPECT_EQ(plan[1].nextBeads, (std::vector<std::size_t>{2, 3}));

  EXPECT_FALSE(plan[2].rigidBody);
  EXPECT_EQ(plan[2].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[2].currentBead, 3uz);
  EXPECT_EQ(plan[2].previousBead, std::optional<std::size_t>{2});
  EXPECT_EQ(plan[2].nextBeads, (std::vector<std::size_t>{4}));
}

TEST(MC_SEMI_FLEXIBLE_CBMC, nvt_drift_and_rigid_geometry_cbmc)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {pentane}, {}, {30}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, biphenyl_component_parsing)
{
  const ForceField forceField = makeBiphenylForceField();
  Component biphenyl = makeDiethylBiphenyl(forceField, 0, MCMoveProbabilities());

  // six fragments (ordered by lowest atom index): flexible {0}, {1}, rigid ring A {2..7},
  // rigid ring B {8..13}, flexible {14}, {15}
  const FragmentGraph& graph = biphenyl.fragmentGraph;
  ASSERT_EQ(graph.fragments.size(), 6uz);
  EXPECT_FALSE(graph.fragments[0].isRigidBody());
  EXPECT_FALSE(graph.fragments[1].isRigidBody());
  EXPECT_TRUE(graph.fragments[2].isRigidBody());
  EXPECT_TRUE(graph.fragments[3].isRigidBody());
  EXPECT_FALSE(graph.fragments[4].isRigidBody());
  EXPECT_FALSE(graph.fragments[5].isRigidBody());
  EXPECT_EQ(graph.fragments[2].atoms, (std::vector<std::size_t>{2, 3, 4, 5, 6, 7}));
  EXPECT_EQ(graph.fragments[3].atoms, (std::vector<std::size_t>{8, 9, 10, 11, 12, 13}));

  // interactions entirely inside a ring are dropped; only the junction terms remain:
  // bonds 0-1, 1-2, 5-8 (biphenyl), 11-14, 14-15
  EXPECT_EQ(biphenyl.intraMolecularPotentials.bonds.size(), 5uz);
  // bends 0-1-2, {1-2-3, 1-2-7}, {4-5-8, 6-5-8}, {5-8-9, 5-8-13}, {10-11-14, 12-11-14}, 11-14-15
  EXPECT_EQ(biphenyl.intraMolecularPotentials.bends.size(), 10uz);
  // 16 crossing torsions, including the four biphenyl-twist torsions around the 5-8 bond
  EXPECT_EQ(biphenyl.intraMolecularPotentials.torsions.size(), 16uz);
}

TEST(MC_SEMI_FLEXIBLE_CBMC, biphenyl_growth_plan)
{
  const ForceField forceField = makeBiphenylForceField();
  Component biphenyl = makeDiethylBiphenyl(forceField, 0, MCMoveProbabilities());

  std::vector<CBMC::GrowStep> plan = biphenyl.growthPlan({biphenyl.startingBead});

  // CH3(0) -> CH2(1) -> peel ipso(2) -> ring A {3,4,5,6,7} hinged on 2 -> peel ipso(8) across the
  // biphenyl bond -> ring B {9,10,11,12,13} hinged on 8 -> CH2(14) -> CH3(15)
  ASSERT_EQ(plan.size(), 7uz);

  EXPECT_FALSE(plan[0].rigidBody);
  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::PlaceSeedFragment);
  EXPECT_EQ(plan[0].currentBead, 0uz);
  EXPECT_EQ(plan[0].nextBeads, (std::vector<std::size_t>{1}));

  EXPECT_FALSE(plan[1].rigidBody);
  EXPECT_EQ(plan[1].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[1].currentBead, 1uz);
  EXPECT_EQ(plan[1].previousBead, std::optional<std::size_t>{0});
  EXPECT_EQ(plan[1].nextBeads, (std::vector<std::size_t>{2}));

  EXPECT_TRUE(plan[2].rigidBody);
  EXPECT_EQ(plan[2].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[2].currentBead, 2uz);
  EXPECT_EQ(plan[2].previousBead, std::optional<std::size_t>{1});
  EXPECT_EQ(plan[2].nextBeads, (std::vector<std::size_t>{3, 4, 5, 6, 7}));

  EXPECT_FALSE(plan[3].rigidBody);
  EXPECT_EQ(plan[3].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[3].currentBead, 5uz);
  EXPECT_EQ(plan[3].previousBead, std::optional<std::size_t>{4});
  EXPECT_EQ(plan[3].nextBeads, (std::vector<std::size_t>{8}));

  EXPECT_TRUE(plan[4].rigidBody);
  EXPECT_EQ(plan[4].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[4].currentBead, 8uz);
  EXPECT_EQ(plan[4].previousBead, std::optional<std::size_t>{5});
  EXPECT_EQ(plan[4].nextBeads, (std::vector<std::size_t>{9, 10, 11, 12, 13}));

  EXPECT_FALSE(plan[5].rigidBody);
  EXPECT_EQ(plan[5].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[5].currentBead, 11uz);
  EXPECT_EQ(plan[5].previousBead, std::optional<std::size_t>{10});
  EXPECT_EQ(plan[5].nextBeads, (std::vector<std::size_t>{14}));

  EXPECT_FALSE(plan[6].rigidBody);
  EXPECT_EQ(plan[6].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[6].currentBead, 14uz);
  EXPECT_EQ(plan[6].previousBead, std::optional<std::size_t>{11});
  EXPECT_EQ(plan[6].nextBeads, (std::vector<std::size_t>{15}));
}

TEST(MC_SEMI_FLEXIBLE_CBMC, biphenyl_nvt_drift_and_ring_geometry_cbmc)
{
  const ForceField forceField = makeBiphenylForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component biphenyl = makeDiethylBiphenyl(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {biphenyl}, {}, {10}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingGeometryPreserved(s);
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, biphenyl_nvt_drift_and_ring_geometry_recoil_growth)
{
  ForceField forceField = makeBiphenylForceField();
  forceField.useRecoilGrowth = true;

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component biphenyl = makeDiethylBiphenyl(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {biphenyl}, {}, {10}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingGeometryPreserved(s);
  }
}

// The basic example defines four 'Partial-reinsertion' fixed-atom sets; parsing must accept them
// (the group-aware validation allows regrowing a whole tail or one half of the molecule while the
// rest, including the rigid rings, stays fixed).
TEST(MC_SEMI_FLEXIBLE_CBMC, basic_example_biphenyl_partial_reinsertion_parsing)
{
  const ForceField forceField = makeBiphenylForceField();
  Component biphenyl = makeDiethylBiphenylFromBasicExample(forceField, 0, MCMoveProbabilities());

  ASSERT_EQ(biphenyl.partialReinsertionFixedAtoms.size(), 4uz);
  EXPECT_EQ(biphenyl.partialReinsertionFixedAtoms[0],
            (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}));
  EXPECT_EQ(biphenyl.partialReinsertionFixedAtoms[1],
            (std::vector<std::size_t>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}));
  EXPECT_EQ(biphenyl.partialReinsertionFixedAtoms[2], (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7}));
  EXPECT_EQ(biphenyl.partialReinsertionFixedAtoms[3], (std::vector<std::size_t>{8, 9, 10, 11, 12, 13, 14, 15}));
}

// Full basic-example move set (translation, rotation, reinsertion, partial reinsertion) with the
// group-aware CBMC growth (box only, 298 K, no charges). The rigid aromatic rings must stay
// exactly rigid and the running energy must not drift.
TEST(MC_SEMI_FLEXIBLE_CBMC, basic_example_biphenyl_nvt_partial_reinsertion_cbmc)
{
  const ForceField forceField = makeBiphenylForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);

  Component biphenyl = makeDiethylBiphenylFromBasicExample(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 298.0, 1e4, 1.0, {}, {biphenyl}, {}, {4}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingGeometryPreserved(s);
  }
}

// Same basic-example move set with recoil growth as the underlying chain-growth scheme.
TEST(MC_SEMI_FLEXIBLE_CBMC, basic_example_biphenyl_nvt_partial_reinsertion_recoil_growth)
{
  ForceField forceField = makeBiphenylForceField();
  forceField.useRecoilGrowth = true;

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);

  Component biphenyl = makeDiethylBiphenylFromBasicExample(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 298.0, 1e4, 1.0, {}, {biphenyl}, {}, {4}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRingGeometryPreserved(s);
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, nvt_drift_and_rigid_geometry_recoil_growth)
{
  ForceField forceField = makeAlkaneForceField();
  forceField.useRecoilGrowth = true;

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {pentane}, {}, {30}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, fragment_rigid_properties)
{
  const ForceField forceField = makeBiphenylForceField();
  Component biphenyl = makeDiethylBiphenyl(forceField, 0, MCMoveProbabilities());

  ASSERT_EQ(biphenyl.fragmentGraph.fragments.size(), 6uz);
  EXPECT_EQ(biphenyl.numberOfRigidFragments(), 2uz);
  EXPECT_TRUE(biphenyl.isSemiFlexible());

  for (std::size_t g : {2uz, 3uz})
  {
    const Fragment& ring = biphenyl.fragmentGraph.fragments[g];
    EXPECT_TRUE(ring.isRigidBody());
    EXPECT_EQ(ring.bodyFixedPositions.size(), 6uz);
    EXPECT_EQ(ring.rotationalDegreesOfFreedom, 3uz);
    EXPECT_EQ(ring.shapeType, 0uz);
    EXPECT_GT(ring.mass, 0.0);
    EXPECT_GT(ring.inertiaVector.x, 0.0);
    EXPECT_GT(ring.inertiaVector.y, 0.0);
    EXPECT_GT(ring.inertiaVector.z, 0.0);
  }

  // Semi-flexible degree-of-freedom count: four flexible atoms (0,1,14,15) contribute 3 each, plus
  // three center-of-mass translations per rigid fragment; rotational DOF is 3 per ring.
  EXPECT_EQ(biphenyl.translationalDegreesOfFreedom, 3uz * 4uz + 3uz * 2uz);
  EXPECT_EQ(biphenyl.rotationalDegreesOfFreedom, 6uz);
}

TEST(MC_SEMI_FLEXIBLE_CBMC, derive_and_regenerate_fragment_state_round_trip)
{
  const ForceField forceField = makeBiphenylForceField();
  Component biphenyl = makeDiethylBiphenyl(forceField, 0, MCMoveProbabilities());

  // Use the component reference geometry as a laboratory configuration.
  std::vector<Atom> atoms;
  atoms.reserve(biphenyl.definedAtoms.size());
  for (const auto& [atom, mass] : biphenyl.definedAtoms)
  {
    atoms.push_back(atom);
  }

  for (std::size_t g = 0; g != biphenyl.fragmentGraph.fragments.size(); ++g)
  {
    if (!biphenyl.fragmentGraph.fragments[g].isRigidBody()) continue;

    GroupState state = biphenyl.deriveFragmentState(g, atoms);

    std::vector<Atom> rebuilt = atoms;
    biphenyl.regenerateFragmentAtoms(state, g, rebuilt);

    for (std::size_t atomIndex : biphenyl.fragmentGraph.fragments[g].atoms)
    {
      double3 delta = rebuilt[atomIndex].position - atoms[atomIndex].position;
      EXPECT_NEAR(delta.length(), 0.0, 1e-9);
    }
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, biphenyl_nve_ring_geometry_molecular_dynamics)
{
  const ForceField forceField = makeBiphenylForceField();

  Component biphenyl = makeDiethylBiphenyl(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {biphenyl}, {}, {8}, 5);
  system.timeStep = 1e-4;

  MolecularDynamics md = makeShortMolecularDynamics({system});
  md.run();

  for (System& s : md.systems)
  {
    // The rigid rings are integrated as rigid bodies, so their internal geometry must be preserved
    // to machine precision throughout the dynamics.
    expectRingGeometryPreserved(s);
  }
}

TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_nve_geometry_molecular_dynamics)
{
  const ForceField forceField = makeAlkaneForceField();

  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {16}, 5);
  system.timeStep = 5e-5;

  MolecularDynamics md = makeShortMolecularDynamics({system});
  md.run();

  for (System& s : md.systems)
  {
    // The rigid core geometry is preserved exactly by the rigid-body integration.
    expectRigidCoreGeometryPreserved(s);
    EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
  }
}

namespace
{
// Evaluate the total minimization energy of the current configuration for a fixed-cell layout.
// The Hessian capability is requested because evaluateDerivatives only includes the intermolecular
// contributions on the Hessian path; this keeps the finite-difference energy consistent with the
// analytic derivatives under test.
double evaluateGeneralizedEnergy(System& system, const MinimizationDofLayout& layout)
{
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);
  return results.energy;
}
}  // namespace

// The generalized gradient for a semi-flexible molecule mixes Cartesian degrees of freedom (flexible
// atoms) with rigid-group center-of-mass and orientation degrees of freedom. Validate the analytic
// gradient (assembled from the per-atom Cartesian gradients projected onto the group degrees of
// freedom) against a central finite difference of the energy taken in generalized coordinates.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_minimization_gradient_matches_finite_difference)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {6}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  ASSERT_GT(layout.numDofs(), 0u);

  // Snap the rigid groups onto an exactly rigid configuration so the finite-difference reference is
  // consistent with the rigid-body regeneration used by applyGeneralizedDisplacement.
  applyGeneralizedDisplacement(system, layout, std::vector<double>(layout.numDofs(), 0.0));

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const double h = 1e-5;
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    std::vector<double> step(layout.numDofs(), 0.0);

    System plus = system;
    step[dof] = h;
    applyGeneralizedDisplacement(plus, layout, step);
    const double energyPlus = evaluateGeneralizedEnergy(plus, layout);

    System minus = system;
    step[dof] = -h;
    applyGeneralizedDisplacement(minus, layout, step);
    const double energyMinus = evaluateGeneralizedEnergy(minus, layout);

    const double numerical = (energyPlus - energyMinus) / (2.0 * h);
    const double analytic = gradient[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, 1e-3 * scale) << "dof=" << dof;
  }
}

// The position-position generalized Hessian projects the per-interaction Cartesian second derivatives
// onto the group degrees of freedom (J^T H J plus the gradient-curvature term on rigid orientations).
// Validate it against a central second difference of the energy taken with joint displacements in a
// single chart anchored at the base configuration. (Finite-differencing the analytic gradient would
// re-center the rotation chart at each displaced configuration and pick up a spurious antisymmetric
// torque term on the orientation blocks.)
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_minimization_hessian_matches_finite_difference)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {2}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  const std::size_t numberOfDofs = layout.numDofs();
  ASSERT_GT(numberOfDofs, 0u);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(numberOfDofs, 0.0));

  GeneralizedHessian hessian(numberOfDofs, 0);
  std::vector<double> gradient(numberOfDofs, 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const auto energyAtDisplacement = [&](std::span<const double> step)
  {
    System displaced = system;
    applyGeneralizedDisplacement(displaced, layout, step);
    return evaluateGeneralizedEnergy(displaced, layout);
  };

  std::vector<double> step(numberOfDofs, 0.0);
  const auto secondDifference = [&](std::size_t row, std::size_t column, double h)
  {
    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] += h;
    const double ePP = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] -= h;
    const double ePM = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] += h;
    const double eMP = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] -= h;
    const double eMM = energyAtDisplacement(step);

    return (ePP - ePM - eMP + eMM) / (4.0 * h * h);
  };

  const double h = 1e-3;
  for (std::size_t row = 0; row < numberOfDofs; ++row)
  {
    for (std::size_t column = row; column < numberOfDofs; ++column)
    {
      // Richardson extrapolation removes the leading O(h^2) truncation error.
      const double numerical = (4.0 * secondDifference(row, column, h) - secondDifference(row, column, 2.0 * h)) / 3.0;
      const double analytic = hessian(row, column);
      const double scale = std::max({10.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, 5e-3 * scale) << "row=" << row << " column=" << column;
    }
  }
}

// Constant-pressure (NPT) molecular dynamics couples the barostat to the rigid-group centers of
// mass and the flexible atoms: the cell fluctuates while the rigid core geometry is preserved
// exactly by the rigid-body integration.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_npt_geometry_molecular_dynamics)
{
  const ForceField forceField = makeAlkaneForceField();

  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {16}, 5);
  system.timeStep = 5e-5;
  system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::NPT;
  system.thermostat = Thermostat(300.0, system.timeStep, system.translationalDegreesOfFreedom,
                                 system.rotationalDegreesOfFreedom, 3, 3, 0.15);
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta, 300.0,
                     system.pressure, system.timeStep, system.translationalDegreesOfFreedom);

  MolecularDynamics md = makeShortMolecularDynamics({system});
  md.run();

  for (System& s : md.systems)
  {
    expectRigidCoreGeometryPreserved(s);
    EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
    EXPECT_GT(s.simulationBox.volume, 0.0);
  }
}

namespace
{
// The stored group states must stay aligned with the molecule layout: one GroupState per rigid
// fragment, in molecule order, and regenerating the rigid-fragment atoms from each stored state
// must reproduce the current atom positions of that molecule.
void expectGroupStatesAligned(const System& system)
{
  std::size_t expectedGroups{};
  for (const Molecule& molecule : system.moleculeData)
  {
    expectedGroups += system.components[molecule.componentId].numberOfRigidFragments();
  }
  ASSERT_EQ(system.groupData.size(), expectedGroups);

  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::size_t atomIndex{};
  std::size_t groupIndex{};
  for (const Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    if (component.isSemiFlexible())
    {
      std::vector<Atom> regenerated(atoms.begin() + static_cast<std::ptrdiff_t>(atomIndex),
                                    atoms.begin() + static_cast<std::ptrdiff_t>(atomIndex + molecule.numberOfAtoms));
      std::size_t rigidRank{};
      for (std::size_t g = 0; g != component.fragmentGraph.fragments.size(); ++g)
      {
        if (component.fragmentGraph.fragments[g].isRigidBody())
        {
          component.regenerateFragmentAtoms(system.groupData[groupIndex + rigidRank], g, regenerated);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        EXPECT_NEAR((regenerated[i].position - atoms[atomIndex + i].position).length(), 0.0, 1e-8);
      }
      groupIndex += component.numberOfRigidFragments();
    }
    atomIndex += molecule.numberOfAtoms;
  }
}
}  // namespace

// Grand-canonical MD (MuVT): the particle-exchange step inserts and deletes whole molecules; for
// semi-flexible components the per-group rigid-body states (GroupState) are spliced in and out of
// System::groupData so that the rigid-body integration stays aligned with the molecule layout. The
// alignment is verified every cycle through the stage callbacks while exchanges are happening.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_muvt_geometry_molecular_dynamics)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);
  system.timeStep = 5e-5;
  system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::MuVT;
  system.thermostat = Thermostat(300.0, system.timeStep, system.translationalDegreesOfFreedom,
                                 system.rotationalDegreesOfFreedom, 3, 3, 0.15);

  // More production cycles than the other MD tests: particle exchange is only attempted every
  // third cycle and the acceptance assertions below need accepted exchanges within the MD stages.
  MolecularDynamics md = MolecularDynamics({60, 0, 5, 5, 1000, 10000, 5000, 5000}, {std::move(system)}, 42uz, 5, false);

  const auto verifyGroupStates = [&md]()
  {
    for (System& s : md.systems)
    {
      expectGroupStatesAligned(s);
      expectRigidCoreGeometryPreserved(s);
      EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
    }
  };

  md.setup();
  md.preInitialize();
  md.initialize();
  md.equilibrate(verifyGroupStates, 1);
  md.production(verifyGroupStates, 1);

  for (System& s : md.systems)
  {
    // both an insertion and a deletion must have been accepted, otherwise the group-state splicing
    // was not exercised (the run is deterministic through the fixed seed)
    const auto& swapStatistics =
        std::get<MoveStatistics<double3>>(s.components[0].mc_moves_statistics[Move::Types::SwapCBMC]);
    const double3 acceptedCounts = swapStatistics.accepted + swapStatistics.totalAccepted;
    EXPECT_GT(acceptedCounts.x, 0.0);
    EXPECT_GT(acceptedCounts.y, 0.0);

    expectGroupStatesAligned(s);
    expectRigidCoreGeometryPreserved(s);
    EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
  }
}

// Grand-canonical MD at constant pressure (MuPT): particle exchange (group-state splicing) combined
// with the barostat propagation of the rigid-group centers of mass and the flexible atoms.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_mupt_geometry_molecular_dynamics)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);
  system.timeStep = 5e-5;
  system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::MuPT;
  system.thermostat = Thermostat(300.0, system.timeStep, system.translationalDegreesOfFreedom,
                                 system.rotationalDegreesOfFreedom, 3, 3, 0.15);
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::MuPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta, 300.0,
                     system.pressure, system.timeStep, system.translationalDegreesOfFreedom);

  MolecularDynamics md = MolecularDynamics({60, 0, 5, 5, 1000, 10000, 5000, 5000}, {std::move(system)}, 42uz, 5, false);

  const auto verifyGroupStates = [&md]()
  {
    for (System& s : md.systems)
    {
      expectGroupStatesAligned(s);
      expectRigidCoreGeometryPreserved(s);
      EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
    }
  };

  md.setup();
  md.preInitialize();
  md.initialize();
  md.equilibrate(verifyGroupStates, 1);
  md.production(verifyGroupStates, 1);

  for (System& s : md.systems)
  {
    const auto& swapStatistics =
        std::get<MoveStatistics<double3>>(s.components[0].mc_moves_statistics[Move::Types::SwapCBMC]);
    const double3 acceptedCounts = swapStatistics.accepted + swapStatistics.totalAccepted;
    EXPECT_GT(acceptedCounts.x + acceptedCounts.y, 0.0);

    expectGroupStatesAligned(s);
    expectRigidCoreGeometryPreserved(s);
    EXPECT_TRUE(std::isfinite(s.runningEnergies.conservedEnergy()));
    EXPECT_GT(s.simulationBox.volume, 0.0);
  }
}

// Under a symmetric logarithmic cell strain the rigid-group centers of mass follow the cell while
// the internal group geometry does not scale; the flexible atoms deform affinely.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_variable_cell_strain_preserves_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {4}, 5);
  system.cellMinimizationType = CellMinimizationType::Regular;

  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  // Snap the rigid groups onto an exactly rigid configuration first.
  applyGeneralizedDisplacement(system, layout, std::vector<double>(layout.numDofs(), 0.0));

  std::vector<double> displacement(layout.numDofs(), 0.0);
  displacement[*layout.cellDof(0)] = -0.15;
  displacement[*layout.cellDof(1)] = 0.06;
  displacement[*layout.cellDof(5)] = 0.09;
  applyGeneralizedDisplacement(system, layout, displacement);

  EXPECT_GT(system.simulationBox.volume, 0.0);
  expectRigidCoreGeometryPreserved(system);
}

// Variable-cell derivatives for semi-flexible molecules: the cell-strain gradient (molecular virial
// with group center-of-mass arms) and the mixed position-cell and cell-cell Hessian blocks (affine
// contraction over group and flexible degrees of freedom) must match finite differences of the
// energy taken in generalized coordinates.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_variable_cell_derivatives_match_finite_difference)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {2}, 5);
  system.cellMinimizationType = CellMinimizationType::Regular;

  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  ASSERT_EQ(layout.numberOfCellDofs(), 6uz);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(layout.numDofs(), 0.0));

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const auto energyAtDisplacement = [&](const std::vector<double>& step)
  {
    System displaced = system;
    applyGeneralizedDisplacement(displaced, layout, step);
    return evaluateGeneralizedEnergy(displaced, layout);
  };
  const double referenceEnergy = energyAtDisplacement(std::vector<double>(layout.numDofs(), 0.0));

  constexpr double delta = 2.0e-5;
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const std::size_t column = *layout.cellDof(a);
    std::vector<double> plus(layout.numDofs(), 0.0);
    std::vector<double> minus(layout.numDofs(), 0.0);
    plus[column] = delta;
    minus[column] = -delta;
    const double energyPlus = energyAtDisplacement(plus);
    const double energyMinus = energyAtDisplacement(minus);

    const double numericalGradient = (energyPlus - energyMinus) / (2.0 * delta);
    const double gradientScale = std::max({1.0, std::abs(numericalGradient), std::abs(gradient[column])});
    EXPECT_NEAR(gradient[column], numericalGradient, 1.0e-4 * gradientScale) << "cell coordinate=" << a;

    for (std::size_t row = 0; row < layout.numDofs(); ++row)
    {
      double numericalHessian{};
      if (row == column)
      {
        numericalHessian = (energyPlus - 2.0 * referenceEnergy + energyMinus) / (delta * delta);
      }
      else
      {
        std::vector<double> pp(layout.numDofs(), 0.0);
        std::vector<double> pm(layout.numDofs(), 0.0);
        std::vector<double> mp(layout.numDofs(), 0.0);
        std::vector<double> mm(layout.numDofs(), 0.0);
        pp[row] = pp[column] = delta;
        pm[row] = delta;
        pm[column] = -delta;
        mp[row] = -delta;
        mp[column] = delta;
        mm[row] = mm[column] = -delta;
        numericalHessian = (energyAtDisplacement(pp) - energyAtDisplacement(pm) - energyAtDisplacement(mp) +
                            energyAtDisplacement(mm)) /
                           (4.0 * delta * delta);
      }
      const double analyticHessian = hessian(row, column);
      const double hessianScale = std::max({1.0, std::abs(numericalHessian), std::abs(analyticHessian)});
      EXPECT_NEAR(analyticHessian, numericalHessian, 2.0e-3 * hessianScale) << "row=" << row << " cell coordinate=" << a;
    }
  }
}

// End-to-end Baker eigenvector-following minimization of semi-flexible pentane molecules: the energy
// must decrease to a local minimum with a vanishing generalized gradient, and the rigid core must be
// preserved exactly since it only moves through its six rigid-body degrees of freedom.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_minimization_reaches_local_minimum)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {4}, 5);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 500;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-6;
  options.maxGradientTolerance = 1e-5;
  options.printEvery = 500;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());

  ASSERT_EQ(minimization.results.size(), 1uz);
  EXPECT_TRUE(minimization.results[0].converged);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);

  expectRigidCoreGeometryPreserved(minimization.systems[0]);
}

// Variable-cell Baker minimization of semi-flexible pentane: the driver takes combined position and
// cell steps; the rigid core geometry must stay exact because the strain only drives the group
// centers of mass.
TEST(MC_SEMI_FLEXIBLE_CBMC, pentane_variable_cell_minimization_preserves_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {4}, 5);
  system.cellMinimizationType = CellMinimizationType::Isotropic;

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 2000;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-6;
  options.maxGradientTolerance = 1e-5;
  options.printEvery = 2000;
  Minimization minimization(options, {system}, false);
  // A dilute gas with a volume degree of freedom converges slowly; the driver throws when the step
  // budget is exhausted, which is acceptable for this smoke test as long as the state stays sane.
  try
  {
    minimization.run();
  }
  catch (const std::runtime_error&)
  {
  }

  ASSERT_EQ(minimization.results.size(), 1uz);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);
  EXPECT_GT(minimization.systems[0].simulationBox.volume, 0.0);
  expectRigidCoreGeometryPreserved(minimization.systems[0]);
}

namespace
{
// Charged variant of the alkane force field (net-neutral pentane: 2 x CH3(-0.3) + 3 x CH2(+0.2)).
// Exercises the Ewald reciprocal-space and exclusion Hessians for semi-flexible molecules: the
// rigid-group atoms couple through the group center-of-mass and orientation degrees of freedom.
ForceField makeChargedAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, -0.3, 0.0, 6, false}, {"CH2", false, 14.03, 0.2, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, true);
}
}  // namespace

// With electrostatics enabled, the generalized gradient additionally contains the Ewald
// reciprocal-space forces (projected per site onto group center-of-mass and orientation degrees of
// freedom) and the intramolecular exclusion terms (energy-only inside a rigid group, full coupling
// across group boundaries). Validate against a central finite difference of the energy.
TEST(MC_SEMI_FLEXIBLE_CBMC, charged_pentane_minimization_gradient_matches_finite_difference)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {4}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  ASSERT_GT(layout.numDofs(), 0u);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(layout.numDofs(), 0.0));

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const double h = 1e-5;
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    std::vector<double> step(layout.numDofs(), 0.0);

    System plus = system;
    step[dof] = h;
    applyGeneralizedDisplacement(plus, layout, step);
    const double energyPlus = evaluateGeneralizedEnergy(plus, layout);

    System minus = system;
    step[dof] = -h;
    applyGeneralizedDisplacement(minus, layout, step);
    const double energyMinus = evaluateGeneralizedEnergy(minus, layout);

    const double numerical = (energyPlus - energyMinus) / (2.0 * h);
    const double analytic = gradient[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, 1e-3 * scale) << "dof=" << dof;
  }
}

// Full generalized Hessian with electrostatics: the Ewald reciprocal-space part contributes
// site-projected phase curvature (including orientation-orientation curvature terms) and the
// exclusion part contributes mixed rigid-body / Cartesian pair blocks. Validated against a central
// second difference of the energy in a single rotation chart anchored at the base configuration.
TEST(MC_SEMI_FLEXIBLE_CBMC, charged_pentane_minimization_hessian_matches_finite_difference)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {2}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  const std::size_t numberOfDofs = layout.numDofs();
  ASSERT_GT(numberOfDofs, 0u);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(numberOfDofs, 0.0));

  GeneralizedHessian hessian(numberOfDofs, 0);
  std::vector<double> gradient(numberOfDofs, 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const auto energyAtDisplacement = [&](std::span<const double> step)
  {
    System displaced = system;
    applyGeneralizedDisplacement(displaced, layout, step);
    return evaluateGeneralizedEnergy(displaced, layout);
  };

  std::vector<double> step(numberOfDofs, 0.0);
  const auto secondDifference = [&](std::size_t row, std::size_t column, double h)
  {
    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] += h;
    const double ePP = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] -= h;
    const double ePM = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] += h;
    const double eMP = energyAtDisplacement(step);

    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] -= h;
    const double eMM = energyAtDisplacement(step);

    return (ePP - ePM - eMP + eMM) / (4.0 * h * h);
  };

  const double h = 1e-3;
  for (std::size_t row = 0; row < numberOfDofs; ++row)
  {
    for (std::size_t column = row; column < numberOfDofs; ++column)
    {
      // Richardson extrapolation removes the leading O(h^2) truncation error.
      const double numerical = (4.0 * secondDifference(row, column, h) - secondDifference(row, column, 2.0 * h)) / 3.0;
      const double analytic = hessian(row, column);
      const double scale = std::max({10.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, 5e-3 * scale) << "row=" << row << " column=" << column;
    }
  }
}

// Variable-cell derivatives with electrostatics: the Ewald reciprocal-space cell blocks use the
// group center-of-mass strain phases (internal offsets do not scale with the cell) and the
// exclusion pairs crossing a rigid-group boundary use the center-of-mass virial arm. Both the
// cell gradient and the mixed position-cell / cell-cell Hessian blocks must match finite
// differences of the energy in generalized coordinates.
TEST(MC_SEMI_FLEXIBLE_CBMC, charged_pentane_variable_cell_derivatives_match_finite_difference)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {2}, 5);
  system.cellMinimizationType = CellMinimizationType::Regular;

  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  ASSERT_EQ(layout.numberOfCellDofs(), 6uz);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(layout.numDofs(), 0.0));

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const auto energyAtDisplacement = [&](const std::vector<double>& step)
  {
    System displaced = system;
    applyGeneralizedDisplacement(displaced, layout, step);
    return evaluateGeneralizedEnergy(displaced, layout);
  };
  const double referenceEnergy = energyAtDisplacement(std::vector<double>(layout.numDofs(), 0.0));

  constexpr double delta = 2.0e-5;
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const std::size_t column = *layout.cellDof(a);
    std::vector<double> plus(layout.numDofs(), 0.0);
    std::vector<double> minus(layout.numDofs(), 0.0);
    plus[column] = delta;
    minus[column] = -delta;
    const double energyPlus = energyAtDisplacement(plus);
    const double energyMinus = energyAtDisplacement(minus);

    const double numericalGradient = (energyPlus - energyMinus) / (2.0 * delta);
    const double gradientScale = std::max({1.0, std::abs(numericalGradient), std::abs(gradient[column])});
    EXPECT_NEAR(gradient[column], numericalGradient, 1.0e-4 * gradientScale) << "cell coordinate=" << a;

    for (std::size_t row = 0; row < layout.numDofs(); ++row)
    {
      double numericalHessian{};
      if (row == column)
      {
        numericalHessian = (energyPlus - 2.0 * referenceEnergy + energyMinus) / (delta * delta);
      }
      else
      {
        std::vector<double> pp(layout.numDofs(), 0.0);
        std::vector<double> pm(layout.numDofs(), 0.0);
        std::vector<double> mp(layout.numDofs(), 0.0);
        std::vector<double> mm(layout.numDofs(), 0.0);
        pp[row] = pp[column] = delta;
        pm[row] = delta;
        pm[column] = -delta;
        mp[row] = -delta;
        mp[column] = delta;
        mm[row] = mm[column] = -delta;
        numericalHessian = (energyAtDisplacement(pp) - energyAtDisplacement(pm) - energyAtDisplacement(mp) +
                            energyAtDisplacement(mm)) /
                           (4.0 * delta * delta);
      }
      const double analyticHessian = hessian(row, column);
      const double hessianScale = std::max({1.0, std::abs(numericalHessian), std::abs(analyticHessian)});
      EXPECT_NEAR(analyticHessian, numericalHessian, 2.0e-3 * hessianScale) << "row=" << row << " cell coordinate=" << a;
    }
  }
}

// End-to-end Baker minimization with electrostatics: the energy decreases to a local minimum while
// the rigid core geometry stays exact.
TEST(MC_SEMI_FLEXIBLE_CBMC, charged_pentane_minimization_reaches_local_minimum)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {4}, 5);

  MinimizationOptions options{};
  options.maximumNumberOfSteps = 1000;
  options.maximumStepLength = 0.2;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1e-6;
  options.maxGradientTolerance = 1e-5;
  options.printEvery = 1000;
  Minimization minimization(options, {system}, false);
  ASSERT_NO_THROW(minimization.run());

  ASSERT_EQ(minimization.results.size(), 1uz);
  EXPECT_TRUE(minimization.results[0].converged);
  EXPECT_LT(minimization.results[0].finalEnergy, minimization.results[0].initialEnergy);

  expectRigidCoreGeometryPreserved(minimization.systems[0]);
}

namespace
{
// Validate the analytic generalized gradient and position-position Hessian of a semi-flexible
// pentane carrying a higher-order cross term against finite differences of the energy taken in
// generalized coordinates. The injected term spans the rigid core / flexible boundary, so it
// couples the flexible-atom Cartesian degrees of freedom to the rigid-group center-of-mass and
// orientation degrees of freedom -- exactly the projection under test.
void runSemiFlexibleCrossTermTest(const Potentials::IntraMolecularPotentials& extra)
{
  const ForceField forceField = makeAlkaneForceField();
  Component pentane = makeSemiFlexiblePentane(forceField, 0, MCMoveProbabilities());

  Potentials::IntraMolecularPotentials& p = pentane.intraMolecularPotentials;
  p.bondBonds.insert(p.bondBonds.end(), extra.bondBonds.begin(), extra.bondBonds.end());
  p.bondBends.insert(p.bondBends.end(), extra.bondBends.begin(), extra.bondBends.end());
  p.bendBends.insert(p.bendBends.end(), extra.bendBends.begin(), extra.bendBends.end());
  p.bondTorsions.insert(p.bondTorsions.end(), extra.bondTorsions.begin(), extra.bondTorsions.end());
  p.bendTorsions.insert(p.bendTorsions.end(), extra.bendTorsions.begin(), extra.bendTorsions.end());
  p.inversionBends.insert(p.inversionBends.end(), extra.inversionBends.begin(), extra.inversionBends.end());

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {pentane}, {}, {1}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  const std::size_t numberOfDofs = layout.numDofs();
  ASSERT_GT(numberOfDofs, 0u);

  applyGeneralizedDisplacement(system, layout, std::vector<double>(numberOfDofs, 0.0));

  GeneralizedHessian hessian(numberOfDofs, 0);
  std::vector<double> gradient(numberOfDofs, 0.0);
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout,
                      DerivativeCapabilities{.energy = true, .gradient = true, .hessianPositionPosition = true},
                      results);

  const auto energyAtDisplacement = [&](std::span<const double> step)
  {
    System displaced = system;
    applyGeneralizedDisplacement(displaced, layout, step);
    return evaluateGeneralizedEnergy(displaced, layout);
  };

  // Gradient: central difference.
  const double hGradient = 1e-5;
  for (std::size_t dof = 0; dof < numberOfDofs; ++dof)
  {
    std::vector<double> step(numberOfDofs, 0.0);
    step[dof] = hGradient;
    const double energyPlus = energyAtDisplacement(step);
    step[dof] = -hGradient;
    const double energyMinus = energyAtDisplacement(step);
    const double numerical = (energyPlus - energyMinus) / (2.0 * hGradient);
    const double scale = std::max({1.0, std::abs(numerical), std::abs(gradient[dof])});
    EXPECT_NEAR(numerical, gradient[dof], 1e-3 * scale) << "gradient dof=" << dof;
  }

  // Hessian: central second difference in a single chart anchored at the base configuration, with
  // Richardson extrapolation to remove the leading O(h^2) truncation error.
  std::vector<double> step(numberOfDofs, 0.0);
  const auto secondDifference = [&](std::size_t row, std::size_t column, double h)
  {
    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] += h;
    const double ePP = energyAtDisplacement(step);
    std::ranges::fill(step, 0.0);
    step[row] += h;
    step[column] -= h;
    const double ePM = energyAtDisplacement(step);
    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] += h;
    const double eMP = energyAtDisplacement(step);
    std::ranges::fill(step, 0.0);
    step[row] -= h;
    step[column] -= h;
    const double eMM = energyAtDisplacement(step);
    return (ePP - ePM - eMP + eMM) / (4.0 * h * h);
  };

  const double h = 1e-3;
  for (std::size_t row = 0; row < numberOfDofs; ++row)
  {
    for (std::size_t column = row; column < numberOfDofs; ++column)
    {
      const double numerical = (4.0 * secondDifference(row, column, h) - secondDifference(row, column, 2.0 * h)) / 3.0;
      const double analytic = hessian(row, column);
      const double scale = std::max({10.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, 5e-3 * scale) << "row=" << row << " column=" << column;
    }
  }
}
}  // namespace

// Bond-bend cross term across the rigid core / flexible boundary (atoms 0 flexible, 1-2 rigid).
TEST(MC_SEMI_FLEXIBLE_CBMC, semi_flexible_bond_bend_cross_hessian_matches_finite_difference)
{
  Potentials::IntraMolecularPotentials extra{};
  extra.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::CVFF, {109.5, 40000.0, 1.5, 40000.0, 1.5})};
  runSemiFlexibleCrossTermTest(extra);
}

// Bend-bend cross term (angles 0-1-2 and 0-1-3): atom 0 flexible, 1-2-3 rigid core.
TEST(MC_SEMI_FLEXIBLE_CBMC, semi_flexible_bend_bend_cross_hessian_matches_finite_difference)
{
  Potentials::IntraMolecularPotentials extra{};
  extra.bendBends = {BendBendPotential({0, 1, 2, 3}, BendBendType::CVFF, {50000.0, 109.5, 109.5})};
  runSemiFlexibleCrossTermTest(extra);
}

// Bond-torsion cross term (bond 1-2 inside the rigid core, torsion 0-1-2-3 crossing the boundary).
TEST(MC_SEMI_FLEXIBLE_CBMC, semi_flexible_bond_torsion_cross_hessian_matches_finite_difference)
{
  Potentials::IntraMolecularPotentials extra{};
  extra.bondTorsions = {BondTorsionPotential({0, 1, 2, 3}, BondTorsionType::MM3, {2.0, 1.0, 0.5, 1.54})};
  runSemiFlexibleCrossTermTest(extra);
}

// Bend-torsion cross term (angles 0-1-2 and 1-2-3, torsion 0-1-2-3) crossing the boundary.
TEST(MC_SEMI_FLEXIBLE_CBMC, semi_flexible_bend_torsion_cross_hessian_matches_finite_difference)
{
  Potentials::IntraMolecularPotentials extra{};
  extra.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::CVFF, {5000.0, 109.5, 109.5})};
  runSemiFlexibleCrossTermTest(extra);
}

// Inversion-bend cross term (dense four-body Cartesian Hessian projected onto the group DOFs).
TEST(MC_SEMI_FLEXIBLE_CBMC, semi_flexible_inversion_bend_cross_hessian_matches_finite_difference)
{
  Potentials::IntraMolecularPotentials extra{};
  extra.inversionBends = {InversionBendPotential({0, 1, 2, 3}, InversionBendType::Harmonic, {60000.0, 25.0})};
  runSemiFlexibleCrossTermTest(extra);
}

// Grand-canonical CFCMC with a semi-flexible component: the fractional molecule is created and
// regrown through the group-aware CBMC machinery and the lambda moves rescale its interactions.
// The running energies must stay consistent with a full recomputation and the rigid cores (of the
// integer and the fractional molecules alike) must stay exactly rigid.
TEST(MC_SEMI_FLEXIBLE_CBMC, muvt_swap_cfcmc_drift_and_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  // 10 bar at 300 K keeps an equilibrium loading of a handful of molecules in the box, so both
  // insertions and deletions are sampled during the run.
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfFractionalMoleculesPerComponent[0], 1uz);
    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

// Same as above with the configurational-bias flavor (SwapCBCFCMC), whose insertion/deletion and
// lambda-shuffle sub-moves regrow the fractional molecule with the group-aware CBMC growth.
TEST(MC_SEMI_FLEXIBLE_CBMC, muvt_swap_cbcfcmc_drift_and_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfFractionalMoleculesPerComponent[0], 1uz);
    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

// Gibbs ensemble with a semi-flexible component: molecules are exchanged between the two boxes with
// the group-aware CBMC growth while the volumes fluctuate. Both boxes must keep consistent running
// energies and exactly rigid cores.
TEST(MC_SEMI_FLEXIBLE_CBMC, gibbs_swap_cbmc_drift_and_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCBMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {pentane}, {},
                              {20}, 5, systemProbabilities);
  System systemLiquid = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 573.0, 1e4, 1.0, {}, {pentane}, {},
                               {20}, 5, systemProbabilities);

  MonteCarlo mc = makeShortMonteCarlo({systemVapor, systemLiquid});
  mc.run();

  for (System& s : mc.systems)
  {
    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

// Widom test-particle insertions with a semi-flexible component: the ghost molecule is grown with
// the group-aware CBMC machinery (rigid cores hinged on flexible junctions) but never inserted, so
// the system state must be untouched while the sampled Rosenbluth weight is finite and positive.
TEST(MC_SEMI_FLEXIBLE_CBMC, widom_rosenbluth_weight_and_rigid_geometry)
{
  const ForceField forceField = makeAlkaneForceField();

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::Widom, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    const double rosenbluthWeight = s.components[0].averageRosenbluthWeights.averagedRosenbluthWeight();
    EXPECT_TRUE(std::isfinite(rosenbluthWeight));
    EXPECT_GT(rosenbluthWeight, 0.0);

    expectNoEnergyDrift(s);
    expectRigidCoreGeometryPreserved(s);
  }
}

namespace
{
// Charged, polarizable united-atom alkane force field (net-neutral pentane: 2 x CH3(-0.3) +
// 3 x CH2(+0.2)) with Ewald summation and molecule-molecule polarization enabled.
ForceField makePolarizableAlkaneForceField()
{
  ForceField forceField({{"CH3", false, 15.04, -0.3, 1.0, 6, false}, {"CH2", false, 14.03, 0.2, 0.8, 6, false}},
                        {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                        true, true, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  return forceField;
}
}  // namespace

// Polarization combined with semi-flexible molecules: the electric-field bookkeeping of the
// translation/rotation moves and the CBMC regrowth/exchange must stay consistent for molecules with
// rigid groups. The running polarization energy must agree with a full recomputation and the rigid
// cores must remain exactly rigid.
TEST(MC_SEMI_FLEXIBLE_CBMC, nvt_polarization_drift_and_rigid_geometry)
{
  const ForceField forceField = makePolarizableAlkaneForceField();

  // Fixed loading (NVT): particle-exchange moves would deplete the box (the ideal-gas Rosenbluth
  // reference of pentane is far below one), leaving no molecules to carry a polarization energy.
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component pentane = makeSemiFlexiblePentane(forceField, 0, probabilities);

  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e6, 1.0, {}, {pentane}, {}, {16}, 5);

  MonteCarlo mc = makeShortMonteCarlo({system});
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);

    // the fixture must actually exercise a non-trivial polarization energy
    EXPECT_LT(recomputedEnergies.polarization, -1e-8);

    expectRigidCoreGeometryPreserved(s);
  }
}
