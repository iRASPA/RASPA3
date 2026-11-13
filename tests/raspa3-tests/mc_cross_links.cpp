#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import units;
import forcefield;
import vdwparameters;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import bond_potential;
import bend_potential;
import cross_links;
import interactions_cross_link;
import mc_moves_probabilities;
import mc_moves_move_types;
import move_statistics;
import mc_moves_cross_link_swap;
import mc_moves_cross_link_exchange;
import mc_moves_cross_link_formation;
import mc_moves_reinsertion;
import mc_moves_translation;
import mc_moves_partial_reinsertion;
import mc_moves_reaction_common;
import monte_carlo;

// Tests for the cross-links: inter-molecular bonds between reactive sites owned by the system
// (CrossLinkTable), their energy (bond, junction bends, and the exclusion corrections that turn the
// linked pairs into bonded pairs), and the two topology moves (bond swap, formation/scission).
//
// The energy is checked against the one thing it has to reproduce: two molecules joined by a
// cross-link must have exactly the energy of the single molecule with the same bond in its topology.
// The topology moves are checked at frozen positions, where the Boltzmann distribution over the
// link topologies is known exactly.

namespace
{

// Single-bead monomer of pseudo-atom 'atomName' with one reactive site of type 'siteType'.
std::string monomerJson(std::string_view atomName, std::string_view siteType)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["{}", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [[0, "{}"]]
}}
)",
                     atomName, siteType);
}

// Two-bead molecule B-A with a harmonic bond; atom 1 (A) is the reactive site 'X'.
std::string dimerBAJson(double bondK, double bondLength)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["B", [0.0, 0.0, 0.0]], ["A", [{1}, 0.0, 0.0]]],
  "Connectivity" : [[0, 1]],
  "Bonds" : [[["B", "A"], "HARMONIC", [{0}, {1}]]],
  "ReactiveSites" : [{{"Atom": 1, "Type": "X", "Valence": 1}}],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}}
)",
                     bondK, bondLength);
}

// The molecule A-B with a harmonic bond: the reference for two linked monomers.
std::string bondedABJson(double bondK, double bondLength)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["A", [0.0, 0.0, 0.0]], ["B", [{1}, 0.0, 0.0]]],
  "Connectivity" : [[0, 1]],
  "Bonds" : [[["A", "B"], "HARMONIC", [{0}, {1}]]],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}}
)",
                     bondK, bondLength);
}

// The chain B-A-A-B: the reference for two linked B-A dimers with junction bends at the link.
std::string chainBAABJson(double bondK, double bondLength, double bendK, double bendAngle)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["B", [0.0, 0.0, 0.0]], ["A", [{1}, 0.0, 0.0]], ["A", [{4}, 0.0, 0.0]], ["B", [{5}, 0.0, 0.0]]],
  "Connectivity" : [[0, 1], [1, 2], [2, 3]],
  "Bonds" : [[["B", "A"], "HARMONIC", [{0}, {1}]], [["A", "A"], "HARMONIC", [{0}, {1}]]],
  "Bends" : [[["B", "A", "A"], "HARMONIC", [{2}, {3}]]],
  "Torsions" : [[["B", "A", "A", "B"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}}
)",
                     bondK, bondLength, bendK, bendAngle, 2.0 * bondLength, 3.0 * bondLength);
}

ForceField makeForceField(double epsilon, double chargeA, double chargeB, bool useCharge)
{
  return ForceField({{"A", false, 14.03, chargeA, 0.0, 6, false}, {"B", false, 15.04, chargeB, 0.0, 6, false}},
                    {{epsilon, 3.4}, {epsilon, 3.6}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                    true, false, useCharge);
}

Component makeComponent(const ForceField& forceField, std::size_t id, std::string name, std::string_view json)
{
  TemporaryFile file(name + ".json", json);
  return Component(Component::Type::Adsorbate, id, forceField, name, file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false);
}

System makeSystem(const ForceField& forceField, std::vector<Component> components,
                  std::vector<std::vector<double3>> positions, double boxLength = 30.0)
{
  std::vector<std::size_t> create(components.size(), 0);
  System system = System(forceField, SimulationBox(boxLength, boxLength, boxLength), false, 300.0, 1e5, 1.0, {},
                         std::move(components), std::move(positions), create, 5);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}

CrossLinkBondType harmonicBondType(std::string a, std::string b, double k, double r0, double captureRadius,
                                   double formationEnergyKelvin = 0.0)
{
  CrossLinkBondType type{};
  type.siteTypeA = std::move(a);
  type.siteTypeB = std::move(b);
  type.bond = BondPotential({0, 1}, BondType::Harmonic, {k, r0});
  type.captureRadius = captureRadius;
  type.formationEnergy = formationEnergyKelvin * Units::KelvinToEnergy;
  return type;
}

CrossLinkSite site(std::size_t c, std::size_t m, std::size_t a)
{
  return CrossLinkSite{static_cast<std::uint32_t>(c), static_cast<std::uint32_t>(m), static_cast<std::uint32_t>(a)};
}

constexpr double kBondK = 96500.0;
constexpr double kBondLength = 1.54;
constexpr double kBendK = 62500.0;
constexpr double kBendAngle = 114.0;

}  // namespace

// Two charged LJ monomers joined by a cross-link have the energy of the bonded molecule A-B: the
// inter-molecular LJ and real-space Coulomb of the pair are removed, the Ewald reciprocal part of the
// pair is compensated by the same exclusion term the bonded molecule receives, and the bond is added.
TEST(MC_CROSS_LINKS, two_linked_monomers_equal_bonded_molecule_with_ewald)
{
  const ForceField forceField = makeForceField(100.0, 0.4, -0.4, true);

  const double3 posA(10.0, 10.0, 10.0);
  const double3 posB(10.9, 11.1, 10.6);  // stretched, off-axis bond

  Component monomerA = makeComponent(forceField, 0, "xl-monomer-a", monomerJson("A", "X"));
  Component monomerB = makeComponent(forceField, 1, "xl-monomer-b", monomerJson("B", "Y"));
  ASSERT_EQ(monomerA.reactiveSites.size(), 1uz);
  EXPECT_EQ(monomerA.reactiveSites[0].siteType, "X");
  EXPECT_EQ(monomerA.reactiveSites[0].valence, 1uz);

  System linked = makeSystem(forceField, {monomerA, monomerB}, {{posA}, {posB}});
  linked.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 3.0));
  linked.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  const RunningEnergy linkedEnergy = linked.computeTotalEnergies();

  Component bonded = makeComponent(forceField, 0, "xl-bonded-ab", bondedABJson(kBondK, kBondLength));
  System reference = makeSystem(forceField, {bonded}, {{posA, posB}});
  const RunningEnergy referenceEnergy = reference.computeTotalEnergies();

  EXPECT_NEAR(linkedEnergy.potentialEnergy(), referenceEnergy.potentialEnergy(),
              1e-8 * std::abs(referenceEnergy.potentialEnergy()) + 1e-10);
  // The bond shows up in the cross-link slot, and the pair interactions cancel out exactly.
  EXPECT_NEAR(linkedEnergy.crossLink, referenceEnergy.bond, 1e-9 * std::abs(referenceEnergy.bond));
  EXPECT_NEAR(linkedEnergy.moleculeMoleculeVDW, 0.0, 1e-9);
  EXPECT_NEAR(linkedEnergy.moleculeMoleculeCharge, 0.0, 1e-9);
  EXPECT_NEAR(linkedEnergy.ewald_exclusion, referenceEnergy.ewald_exclusion,
              1e-9 * std::abs(referenceEnergy.ewald_exclusion));
  EXPECT_GT(std::abs(referenceEnergy.ewald_exclusion), 1.0);  // the charges actually contribute
}

// Two B-A dimers linked at their A atoms with junction bends have the energy of the chain B-A-A-B.
TEST(MC_CROSS_LINKS, two_linked_dimers_with_junction_bends_equal_chain)
{
  const ForceField forceField = makeForceField(0.0, 0.0, 0.0, false);

  // A bent, stretched geometry so that every term is non-trivial.
  const double3 b0(10.0, 10.0, 10.0);
  const double3 a1(11.5, 10.2, 10.1);
  const double3 a2(12.4, 11.4, 10.5);
  const double3 b3(12.9, 12.6, 11.5);

  Component dimer = makeComponent(forceField, 0, "xl-dimer-ba", dimerBAJson(kBondK, kBondLength));
  ASSERT_EQ(dimer.reactiveSites.size(), 1uz);
  EXPECT_EQ(dimer.reactiveSites[0].atom, 1uz);

  System linked = makeSystem(forceField, {dimer}, {{b0, a1, b3, a2}});
  CrossLinkBondType type = harmonicBondType("X", "X", kBondK, kBondLength, 3.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {kBendK, kBendAngle});
  linked.crossLinks.bondTypes.push_back(type);
  linked.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 1), 0});
  const RunningEnergy linkedEnergy = linked.computeTotalEnergies();

  Component chain = makeComponent(forceField, 0, "xl-chain-baab", chainBAABJson(kBondK, kBondLength, kBendK, kBendAngle));
  System reference = makeSystem(forceField, {chain}, {{b0, a1, a2, b3}});
  const RunningEnergy referenceEnergy = reference.computeTotalEnergies();

  EXPECT_NEAR(linkedEnergy.potentialEnergy(), referenceEnergy.potentialEnergy(), 1e-8);
  EXPECT_NEAR(linkedEnergy.bond + linkedEnergy.crossLink, referenceEnergy.bond + referenceEnergy.bend, 1e-8);
  EXPECT_GT(referenceEnergy.bend, 1.0);
}

// The incremental energy difference of a moved molecule equals the difference of full recomputes,
// including the exclusion corrections (LJ, real-space Coulomb and the Ewald exclusion) of the 1-2
// and 1-3 pairs across the link.
TEST(MC_CROSS_LINKS, energy_difference_of_moved_molecule_matches_recompute)
{
  const ForceField forceField = makeForceField(80.0, 0.3, -0.3, true);

  Component dimer = makeComponent(forceField, 0, "xl-dimer-ba-diff", dimerBAJson(kBondK, kBondLength));
  System system = makeSystem(
      forceField, {dimer},
      {{double3(10.0, 10.0, 10.0), double3(11.5, 10.2, 10.1), double3(12.9, 12.6, 11.5), double3(12.4, 11.4, 10.5),
        double3(4.0, 4.0, 4.0), double3(5.5, 4.0, 4.0)}});
  CrossLinkBondType type = harmonicBondType("X", "X", kBondK, kBondLength, 3.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {kBendK, kBendAngle});
  system.crossLinks.bondTypes.push_back(type);
  system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 1), 0});

  const RunningEnergy before = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());

  std::span<Atom> molecule = system.spanOfMolecule(0, 1);
  std::vector<Atom> oldAtoms(molecule.begin(), molecule.end());
  std::vector<Atom> newAtoms = oldAtoms;
  newAtoms[0].position += double3(0.3, -0.2, 0.4);
  newAtoms[1].position += double3(-0.1, 0.25, 0.15);

  const RunningEnergy difference = system.crossLinkEnergyDifference(0, 1, newAtoms, oldAtoms);
  std::copy(newAtoms.begin(), newAtoms.end(), molecule.begin());
  const RunningEnergy after = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(difference.potentialEnergy(), (after - before).potentialEnergy(), 1e-9);
  EXPECT_NEAR(difference.crossLink, after.crossLink - before.crossLink, 1e-9);
  EXPECT_NEAR(difference.moleculeMoleculeVDW, after.moleculeMoleculeVDW - before.moleculeMoleculeVDW, 1e-9);
  EXPECT_NEAR(difference.ewald_exclusion, after.ewald_exclusion - before.ewald_exclusion, 1e-9);
  EXPECT_GT(std::abs(difference.crossLink), 1e-3);
  EXPECT_GT(std::abs(difference.moleculeMoleculeVDW), 1e-6);
  EXPECT_GT(std::abs(difference.ewald_exclusion), 1e-6);

  // An unlinked molecule contributes nothing.
  std::span<Atom> free = system.spanOfMolecule(0, 2);
  std::vector<Atom> freeOld(free.begin(), free.end());
  std::vector<Atom> freeNew = freeOld;
  freeNew[0].position += double3(0.5, 0.0, 0.0);
  EXPECT_EQ(system.crossLinkEnergyDifference(0, 2, freeNew, freeOld).potentialEnergy(), 0.0);
}

// The analytic gradients of the cross-link energy agree with central finite differences.
TEST(MC_CROSS_LINKS, gradients_match_finite_differences)
{
  const ForceField forceField = makeForceField(80.0, 0.3, -0.3, true);

  Component dimer = makeComponent(forceField, 0, "xl-dimer-ba-grad", dimerBAJson(kBondK, kBondLength));
  System system = makeSystem(
      forceField, {dimer},
      {{double3(10.0, 10.0, 10.0), double3(11.5, 10.2, 10.1), double3(12.9, 12.6, 11.5), double3(12.4, 11.4, 10.5)}});
  CrossLinkBondType type = harmonicBondType("X", "X", kBondK, kBondLength, 3.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {kBendK, kBendAngle});
  system.crossLinks.bondTypes.push_back(type);
  system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 1), 0});

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::vector<AtomDynamics> dynamics(atoms.size());
  const RunningEnergy energy = Interactions::computeCrossLinkGradient(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms, dynamics,
      system.crossLinks);
  std::vector<AtomDynamics> strainDynamics(atoms.size());
  const auto [strainEnergy, strain] = Interactions::computeCrossLinkEnergyStrainDerivative(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms,
      strainDynamics, system.crossLinks);
  EXPECT_NEAR(energy.potentialEnergy(), strainEnergy.potentialEnergy(), 1e-10);

  const double h = 1e-5;
  double3 totalGradient{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    totalGradient += dynamics[i].gradient;
    for (std::size_t k = 0; k < 3; ++k)
    {
      auto shifted = [&](double delta)
      {
        std::vector<Atom> copy(atoms.begin(), atoms.end());
        double3 d{};
        (k == 0 ? d.x : (k == 1 ? d.y : d.z)) = delta;
        copy[i].position += d;
        return Interactions::computeCrossLinkEnergy(forceField, system.simulationBox, system.components,
                                                    system.numberOfMoleculesPerComponent, copy, system.crossLinks)
            .potentialEnergy();
      };
      const double numerical = (shifted(h) - shifted(-h)) / (2.0 * h);
      const double analytic = (k == 0 ? dynamics[i].gradient.x : (k == 1 ? dynamics[i].gradient.y : dynamics[i].gradient.z));
      EXPECT_NEAR(analytic, numerical, 1e-5 * std::max(1.0, std::abs(numerical))) << "atom " << i << " axis " << k;
    }
  }
  // Internal forces sum to zero.
  EXPECT_NEAR(totalGradient.x, 0.0, 1e-8);
  EXPECT_NEAR(totalGradient.y, 0.0, 1e-8);
  EXPECT_NEAR(totalGradient.z, 0.0, 1e-8);
}

namespace
{

// Enumeration of the link topologies of 'sites' with valence 1: all matchings of the sites.
struct Topology
{
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
};

// Every matching exactly once: the smallest unused site is either left unmatched or paired with a
// larger unused site.
void enumerateMatchings(std::size_t numberOfSites, std::size_t from, std::vector<bool>& used, Topology& current,
                        std::vector<Topology>& out)
{
  std::size_t i = from;
  while (i < numberOfSites && used[i]) ++i;
  if (i == numberOfSites)
  {
    out.push_back(current);
    return;
  }
  used[i] = true;
  enumerateMatchings(numberOfSites, i + 1, used, current, out);  // i unmatched
  for (std::size_t j = i + 1; j < numberOfSites; ++j)
  {
    if (used[j]) continue;
    used[j] = true;
    current.pairs.emplace_back(i, j);
    enumerateMatchings(numberOfSites, i + 1, used, current, out);
    current.pairs.pop_back();
    used[j] = false;
  }
  used[i] = false;
}

std::vector<Topology> allMatchings(std::size_t numberOfSites)
{
  std::vector<Topology> out;
  std::vector<bool> used(numberOfSites, false);
  Topology current;
  enumerateMatchings(numberOfSites, 0, used, current, out);
  for (Topology& t : out) std::ranges::sort(t.pairs);
  return out;
}

std::vector<std::pair<std::size_t, std::size_t>> currentTopology(const CrossLinkTable& table)
{
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  for (const CrossLink& link : table.links)
  {
    pairs.emplace_back(std::min(link.a.moleculeIndex, link.b.moleculeIndex),
                       std::max(link.a.moleculeIndex, link.b.moleculeIndex));
  }
  std::ranges::sort(pairs);
  return pairs;
}

}  // namespace

// At frozen positions the formation/scission move (with and without the swap move) samples the exact
// Boltzmann distribution over the link topologies of four valence-1 monomers: the empty topology,
// the six single links and the three perfect matchings, weighted by exp(-beta sum of link energies)
// where every link energy includes the formation energy.
TEST(MC_CROSS_LINKS, topology_moves_sample_exact_distribution_at_frozen_positions)
{
  const ForceField forceField = makeForceField(0.0, 0.0, 0.0, false);
  Component monomer = makeComponent(forceField, 0, "xl-monomer-frozen", monomerJson("A", "X"));

  // Irregular tetrahedron-like arrangement: all pair distances differ, all within the capture radius.
  const std::vector<double3> positions{double3(10.0, 10.0, 10.0), double3(11.9, 10.3, 10.1), double3(10.4, 11.7, 10.6),
                                       double3(11.2, 11.0, 12.0)};
  const std::size_t numberOfSites = positions.size();

  auto run = [&](double swapFraction, unsigned seed)
  {
    System system = makeSystem(forceField, {monomer}, {positions});
    // Soft bond so that several topologies have comparable weight; a negative formation energy favours links.
    system.crossLinks.bondTypes.push_back(harmonicBondType("X", "X", 800.0, 1.8, 6.0, -250.0));

    // Exact weights.
    const std::vector<Topology> topologies = allMatchings(numberOfSites);
    std::vector<double> logWeights;
    for (const Topology& topology : topologies)
    {
      std::vector<CrossLink> links;
      for (const auto& [i, j] : topology.pairs) links.push_back(CrossLink{site(0, i, 0), site(0, j, 0), 0});
      const RunningEnergy energy = Interactions::computeCrossLinkEnergyOfLinks(
          forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
          system.spanOfMoleculeAtoms(), system.crossLinks, links);
      logWeights.push_back(-system.beta * energy.potentialEnergy());
    }
    const double maxLog = *std::ranges::max_element(logWeights);
    double normalization = 0.0;
    for (double lw : logWeights) normalization += std::exp(lw - maxLog);
    std::vector<double> exact;
    for (double lw : logWeights) exact.push_back(std::exp(lw - maxLog) / normalization);

    // Sample.
    RandomNumber random(seed);
    std::map<std::vector<std::pair<std::size_t, std::size_t>>, std::size_t> counts;
    const std::size_t numberOfMoves = 400000;
    RunningEnergy running = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
    for (std::size_t step = 0; step < numberOfMoves; ++step)
    {
      std::optional<RunningEnergy> difference;
      if (random.uniform() < swapFraction)
      {
        difference = MC_Moves::crossLinkSwapMove(random, system);
      }
      else
      {
        difference = MC_Moves::crossLinkFormationScissionMove(random, system);
      }
      if (difference) running += difference.value();
      counts[currentTopology(system.crossLinks)] += 1;

      // Invariants: valence, no duplicates, no intramolecular links.
      for (const CrossLink& link : system.crossLinks.links)
      {
        ASSERT_FALSE(link.a.sameMolecule(link.b));
        ASSERT_LE(system.crossLinks.linkCount(link.a), 1u);
        ASSERT_LE(system.crossLinks.linkCount(link.b), 1u);
      }
    }
    // The running energy tracks the topology changes.
    const RunningEnergy recomputed = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
    EXPECT_NEAR(running.potentialEnergy(), recomputed.potentialEnergy(), 1e-8 * std::max(1.0, std::abs(recomputed.potentialEnergy())));

    ASSERT_EQ(topologies.size(), 10uz);
    for (std::size_t t = 0; t < topologies.size(); ++t)
    {
      const double sampled = static_cast<double>(counts[topologies[t].pairs]) / static_cast<double>(numberOfMoves);
      EXPECT_NEAR(sampled, exact[t], 0.012) << "topology " << t << " (swap fraction " << swapFraction << ")";
    }
    // The distribution is genuinely spread over the topologies.
    EXPECT_LT(*std::ranges::max_element(exact), 0.7);
    if (swapFraction > 0.0)
    {
      EXPECT_GT(std::get<MoveStatistics<double>>(system.mc_moves_statistics[Move::Types::CrossLinkSwap]).totalAccepted,
                0.0);
    }
  };

  run(0.0, 7);
  run(0.5, 11);
}

// Sites with a higher valence carry several links; the moves never exceed the valence, never link a
// site to itself or its own molecule, and never create the same link twice.
TEST(MC_CROSS_LINKS, valence_limits_are_respected)
{
  const ForceField forceField = makeForceField(0.0, 0.0, 0.0, false);

  // Monomer with valence 2: the shorthand form with an explicit valence.
  const std::string json = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["A", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [[0, "X", 2]]
}
)";
  Component monomer = makeComponent(forceField, 0, "xl-monomer-valence2", json);
  ASSERT_EQ(monomer.reactiveSites[0].valence, 2uz);

  std::vector<double3> positions;
  for (std::size_t i = 0; i < 6; ++i)
  {
    positions.emplace_back(10.0 + 1.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / 6.0),
                           10.0 + 1.5 * std::sin(2.0 * std::numbers::pi * static_cast<double>(i) / 6.0), 10.0);
  }
  System system = makeSystem(forceField, {monomer}, {positions});
  system.crossLinks.bondTypes.push_back(harmonicBondType("X", "X", 500.0, 1.6, 6.0, -400.0));

  RandomNumber random(3);
  std::size_t maximumLinks = 0;
  for (std::size_t step = 0; step < 50000; ++step)
  {
    if (random.uniform() < 0.3)
    {
      MC_Moves::crossLinkSwapMove(random, system);
    }
    else
    {
      MC_Moves::crossLinkFormationScissionMove(random, system);
    }
    maximumLinks = std::max(maximumLinks, system.crossLinks.links.size());
    std::set<std::pair<std::uint64_t, std::uint64_t>> seen;
    for (const CrossLink& link : system.crossLinks.links)
    {
      ASSERT_FALSE(link.a.sameMolecule(link.b));
      ASSERT_LE(system.crossLinks.linkCount(link.a), 2u);
      ASSERT_LE(system.crossLinks.linkCount(link.b), 2u);
      const std::uint64_t keyA = link.a.siteKey();
      const std::uint64_t keyB = link.b.siteKey();
      ASSERT_TRUE(seen.insert(std::pair{std::min(keyA, keyB), std::max(keyA, keyB)}).second) << "duplicate link";
    }
  }
  // Six valence-2 sites admit up to six links; the strongly favourable formation energy reaches beyond one per site.
  EXPECT_GT(maximumLinks, 3uz);
  EXPECT_LE(maximumLinks, 6uz);
}

// The bond-exchange move (a,b)+(c,d) -> (a,c)+(b,d) conserves the number of links of every site. At
// frozen positions it must therefore sample the exact Boltzmann distribution over all link topologies
// with the degree sequence of the initial state. Two cases: four valence-1 sites (the three perfect
// matchings) and one valence-2 site with four valence-1 sites (six graphs of degree sequence
// 2,1,1,1,1), where the n_links(c)/n_links(b) proposal factor of the move is exercised.
TEST(MC_CROSS_LINKS, exchange_move_samples_exact_distribution_over_fixed_degree_topologies)
{
  const ForceField forceField = makeForceField(0.0, 0.0, 0.0, false);
  const std::string valence2Json = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["A", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [[0, "X", 2]]
}
)";

  // Global site index g -> (component, molecule).
  struct SiteRef
  {
    std::size_t component;
    std::size_t molecule;
  };
  using Edge = std::pair<std::size_t, std::size_t>;
  using Graph = std::vector<Edge>;

  auto run = [&](std::vector<Component> components, std::vector<std::vector<double3>> positions,
                 std::vector<SiteRef> sites, Graph initial, unsigned seed, std::size_t expectedGraphs,
                 double bondK, double bondLength, double captureRadius)
  {
    System system = makeSystem(forceField, std::move(components), std::move(positions));
    system.crossLinks.bondTypes.push_back(harmonicBondType("X", "X", bondK, bondLength, captureRadius, -250.0));
    auto siteOf = [&](std::size_t g) { return site(sites[g].component, sites[g].molecule, 0); };
    auto indexOf = [&](const CrossLinkSite& cs)
    {
      for (std::size_t g = 0; g < sites.size(); ++g)
      {
        if (sites[g].component == cs.componentId && sites[g].molecule == cs.moleculeIndex) return g;
      }
      throw std::runtime_error("unknown site");
    };
    for (const Edge& e : initial) system.crossLinks.addLink(CrossLink{siteOf(e.first), siteOf(e.second), 0});

    // Degree sequence of the initial state, and all simple graphs on the sites with that sequence.
    const std::size_t n = sites.size();
    std::vector<std::size_t> degree(n, 0);
    for (const Edge& e : initial)
    {
      ++degree[e.first];
      ++degree[e.second];
    }
    std::vector<Edge> allEdges;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j) allEdges.emplace_back(i, j);
    std::vector<Graph> graphs;
    for (std::size_t mask = 0; mask < (1uz << allEdges.size()); ++mask)
    {
      if (static_cast<std::size_t>(std::popcount(mask)) != initial.size()) continue;
      std::vector<std::size_t> d(n, 0);
      Graph g;
      for (std::size_t k = 0; k < allEdges.size(); ++k)
      {
        if (!(mask & (1uz << k))) continue;
        g.push_back(allEdges[k]);
        ++d[allEdges[k].first];
        ++d[allEdges[k].second];
      }
      if (d == degree) graphs.push_back(g);
    }
    ASSERT_EQ(graphs.size(), expectedGraphs);

    std::vector<double> logWeights;
    for (const Graph& g : graphs)
    {
      std::vector<CrossLink> links;
      for (const Edge& e : g) links.push_back(CrossLink{siteOf(e.first), siteOf(e.second), 0});
      const RunningEnergy energy = Interactions::computeCrossLinkEnergyOfLinks(
          forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
          system.spanOfMoleculeAtoms(), system.crossLinks, links);
      logWeights.push_back(-system.beta * energy.potentialEnergy());
    }
    const double maxLog = *std::ranges::max_element(logWeights);
    double normalization = 0.0;
    for (double lw : logWeights) normalization += std::exp(lw - maxLog);
    std::vector<double> exact;
    for (double lw : logWeights) exact.push_back(std::exp(lw - maxLog) / normalization);
    EXPECT_LT(*std::ranges::max_element(exact), 0.7);

    // Sample with the exchange move only.
    RandomNumber random(seed);
    std::map<Graph, std::size_t> counts;
    const std::size_t numberOfMoves = 400000;
    RunningEnergy running = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
    for (std::size_t step = 0; step < numberOfMoves; ++step)
    {
      std::optional<RunningEnergy> difference = MC_Moves::crossLinkExchangeMove(random, system);
      if (difference) running += difference.value();

      Graph current;
      for (const CrossLink& link : system.crossLinks.links)
      {
        const std::size_t i = indexOf(link.a), j = indexOf(link.b);
        current.emplace_back(std::min(i, j), std::max(i, j));
      }
      std::ranges::sort(current);
      counts[current] += 1;

      // Invariants: the degree of every site is conserved, no duplicate or intramolecular links.
      std::vector<std::size_t> d(n, 0);
      for (const Edge& e : current)
      {
        ASSERT_NE(e.first, e.second);
        ++d[e.first];
        ++d[e.second];
      }
      ASSERT_EQ(d, degree);
      ASSERT_EQ(std::ranges::adjacent_find(current), current.end()) << "duplicate link";
    }
    const RunningEnergy recomputed = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
    EXPECT_NEAR(running.potentialEnergy(), recomputed.potentialEnergy(),
                1e-8 * std::max(1.0, std::abs(recomputed.potentialEnergy())));

    for (std::size_t t = 0; t < graphs.size(); ++t)
    {
      const double sampled = static_cast<double>(counts[graphs[t]]) / static_cast<double>(numberOfMoves);
      EXPECT_NEAR(sampled, exact[t], 0.012) << "graph " << t << " (seed " << seed << ")";
    }
    EXPECT_GT(std::get<MoveStatistics<double>>(system.mc_moves_statistics[Move::Types::CrossLinkExchange]).totalAccepted,
              1000.0);
  };

  // Four valence-1 monomers: perfect matchings.
  {
    Component monomer = makeComponent(forceField, 0, "xl-monomer-exchange", monomerJson("A", "X"));
    run({monomer},
        {{double3(10.0, 10.0, 10.0), double3(11.9, 10.3, 10.1), double3(10.4, 11.7, 10.6), double3(11.2, 11.0, 12.0)}},
        {{0, 0}, {0, 1}, {0, 2}, {0, 3}}, {{0, 1}, {2, 3}}, 5, 3uz, 800.0, 1.8, 6.0);
  }
  // One valence-2 site (component 0) and four valence-1 sites (component 1): degree sequence 2,1,1,1,1.
  // The capture radius leaves some pairs out of reach of each other (while every graph stays reachable),
  // so that the candidate sets of the four routes to the same exchange differ in size and the total
  // proposal probability is NOT symmetric: only the per-route n_links(c)/n_links(b) factor makes the
  // move exact here (a plain Metropolis rule is off by up to 0.028 in these probabilities).
  {
    Component hub = makeComponent(forceField, 0, "xl-hub-exchange", valence2Json);
    Component monomer = makeComponent(forceField, 1, "xl-monomer-exchange-b", monomerJson("A", "X"));
    run({hub, monomer},
        {{double3(10.0, 10.0, 10.0)},
         {double3(11.9, 10.3, 10.1), double3(10.4, 11.7, 10.6), double3(11.2, 11.0, 12.0), double3(8.7, 9.2, 11.1)}},
        {{0, 0}, {1, 0}, {1, 1}, {1, 2}, {1, 3}}, {{0, 1}, {0, 2}, {3, 4}}, 9, 6uz, 200.0, 2.5, 3.0);
  }
}

// Deleting a molecule renumbers the sites of the molecules above it; a linked molecule can not be deleted.
TEST(MC_CROSS_LINKS, table_renumbering_on_deletion_and_insertion)
{
  CrossLinkTable table;
  table.bondTypes.push_back(harmonicBondType("X", "X", 1.0, 1.0, 1.0));
  table.addLink(CrossLink{site(0, 3, 0), site(0, 5, 0), 0});
  table.addLink(CrossLink{site(1, 2, 1), site(0, 1, 0), 0});

  EXPECT_TRUE(table.moleculeIsLinked(0, 3));
  EXPECT_TRUE(table.moleculeIsLinked(0, 5));
  EXPECT_TRUE(table.moleculeIsLinked(1, 2));
  EXPECT_FALSE(table.moleculeIsLinked(0, 4));

  // Deleting molecule 4 of component 0 shifts molecule 5 down to 4.
  table.moleculeDeleted(0, 4);
  EXPECT_TRUE(table.isLinked(site(0, 3, 0), site(0, 4, 0)));
  EXPECT_FALSE(table.moleculeIsLinked(0, 5));
  EXPECT_TRUE(table.isLinked(site(1, 2, 1), site(0, 1, 0)));

  // Deleting molecule 0 of component 1 shifts molecule 2 down to 1.
  table.moleculeDeleted(1, 0);
  EXPECT_TRUE(table.isLinked(site(1, 1, 1), site(0, 1, 0)));

  // Inserting at index 1 of component 0 shifts 1, 3 and 4 up.
  table.moleculeInserted(0, 1);
  EXPECT_TRUE(table.isLinked(site(0, 4, 0), site(0, 5, 0)));
  EXPECT_TRUE(table.isLinked(site(1, 1, 1), site(0, 2, 0)));
  EXPECT_EQ(table.linkCount(site(0, 4, 0)), 1u);
  EXPECT_EQ(table.linkCount(site(0, 1, 0)), 0u);

  // A linked molecule can not be deleted.
  EXPECT_THROW(table.moleculeDeleted(0, 4), std::runtime_error);

  // Removing a link frees its sites.
  table.removeLink(0);
  EXPECT_EQ(table.links.size(), 1uz);
  EXPECT_EQ(table.linkCount(site(0, 4, 0)) + table.linkCount(site(1, 1, 1)), 1u);
}

// The reaction moves draw the whole molecules they demote uniformly from all of them and refuse a linked one,
// so that every tuple of unlinked molecules keeps the proposal probability the acceptance rule divides out.
// Leaving the linked molecules out of the draw would favour the others by N! (N - L)! / ((N - nu)! N - L)!.
TEST(MC_CROSS_LINKS, reaction_moves_draw_from_every_whole_molecule_and_refuse_a_linked_one)
{
  const ForceField forceField = makeForceField(80.0, 0.0, 0.0, false);
  Component dimer = makeComponent(forceField, 0, "xl-dimer-draw", dimerBAJson(kBondK, kBondLength));
  std::vector<std::vector<double3>> positions(1);
  for (std::size_t m = 0; m != 6; ++m)
  {
    positions[0].push_back(double3(4.0 * static_cast<double>(m) + 3.0, 5.0, 5.0));
    positions[0].push_back(double3(4.0 * static_cast<double>(m) + 4.5, 5.0, 5.0));
  }
  System system = makeSystem(forceField, {dimer}, positions);
  system.crossLinks.bondTypes.push_back(harmonicBondType("X", "X", kBondK, kBondLength, 3.0));
  // molecules 0-1 and 2-3 linked: two of the six are available for a reaction that demotes one
  system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 0), 0});
  system.crossLinks.addLink(CrossLink{site(0, 2, 1), site(0, 3, 0), 0});

  RandomNumber random(99);
  const std::vector<std::size_t> stoichiometry{1};
  std::vector<std::pair<std::size_t, std::size_t>> selected;
  std::array<std::size_t, 6> picked{};
  std::size_t refused = 0;
  constexpr std::size_t draws = 60000;
  for (std::size_t i = 0; i != draws; ++i)
  {
    if (!MC_Moves::ReactionCommon::selectRandomIntegerMolecules(random, system, stoichiometry, selected))
    {
      ++refused;
      continue;
    }
    ASSERT_EQ(selected.size(), 1uz);
    ASSERT_FALSE(system.crossLinks.moleculeIsLinked(0, selected.front().second));
    ++picked[selected.front().second];
  }
  // four of six draws meet a linked molecule and refuse; the two free molecules are drawn one sixth each
  EXPECT_NEAR(static_cast<double>(refused) / draws, 4.0 / 6.0, 0.02);
  EXPECT_NEAR(static_cast<double>(picked[4]) / draws, 1.0 / 6.0, 0.02);
  EXPECT_NEAR(static_cast<double>(picked[5]) / draws, 1.0 / 6.0, 0.02);
  EXPECT_EQ(picked[0] + picked[1] + picked[2] + picked[3], 0uz);
}

TEST(MC_CROSS_LINKS, table_renumbering_on_exchange)
{
  CrossLinkTable table;
  table.bondTypes.push_back(harmonicBondType("X", "X", 1.0, 1.0, 1.0));
  table.addLink(CrossLink{site(0, 3, 0), site(0, 5, 0), 0});
  table.addLink(CrossLink{site(1, 2, 1), site(0, 1, 0), 0});

  // Exchanging molecules 1 and 5 of component 0 takes each molecule's links along.
  table.moleculesExchanged(0, 1, 5);
  EXPECT_TRUE(table.isLinked(site(0, 3, 0), site(0, 1, 0)));
  EXPECT_TRUE(table.isLinked(site(1, 2, 1), site(0, 5, 0)));
  EXPECT_FALSE(table.moleculeIsLinked(0, 4));
  EXPECT_EQ(table.linkCount(site(0, 1, 0)), 1u);
  EXPECT_EQ(table.linkCount(site(0, 5, 0)), 1u);

  // An exchange with an unlinked molecule moves the links to its index; the other component is untouched.
  table.moleculesExchanged(0, 3, 4);
  EXPECT_TRUE(table.isLinked(site(0, 4, 0), site(0, 1, 0)));
  EXPECT_FALSE(table.moleculeIsLinked(0, 3));
  EXPECT_TRUE(table.moleculeIsLinked(1, 2));
  table.moleculesExchanged(1, 0, 1);
  EXPECT_TRUE(table.moleculeIsLinked(1, 2));
}

// The exchange of two molecules of a system takes the cross-links along with the molecules, so the linked pair
// stays linked at its new indices and the cross-link energy is the same; a linked molecule can not leave for
// another system.
TEST(MC_CROSS_LINKS, swapping_molecules_takes_their_links_along)
{
  const ForceField forceField = makeForceField(80.0, 0.3, -0.3, true);
  Component dimer = makeComponent(forceField, 0, "xl-dimer-ba-swap", dimerBAJson(kBondK, kBondLength));
  const std::vector<std::vector<double3>> positions{{double3(10.0, 10.0, 10.0), double3(11.5, 10.2, 10.1),
                                                     double3(12.9, 12.6, 11.5), double3(12.4, 11.4, 10.5),
                                                     double3(4.0, 4.0, 4.0), double3(5.5, 4.0, 4.0)}};
  System system = makeSystem(forceField, {dimer}, positions);
  CrossLinkBondType type = harmonicBondType("X", "X", kBondK, kBondLength, 3.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {kBendK, kBendAngle});
  system.crossLinks.bondTypes.push_back(type);
  system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 1), 0});
  const RunningEnergy before = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  const double3 linkedPosition = system.spanOfMolecule(0, 1)[1].position;

  // the linked molecule 1 changes places with the free molecule 2
  system.swapMolecules(0, 1, 2);
  EXPECT_EQ(system.spanOfMolecule(0, 2)[1].position, linkedPosition);
  EXPECT_TRUE(system.crossLinks.isLinked(site(0, 0, 1), site(0, 2, 1)));
  EXPECT_FALSE(system.crossLinks.moleculeIsLinked(0, 1));
  const RunningEnergy after = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  EXPECT_NEAR(after.potentialEnergy(), before.potentialEnergy(), 1e-9);
  EXPECT_NEAR(after.crossLink, before.crossLink, 1e-9);
  EXPECT_GT(std::abs(after.crossLink), 1e-6);

  // both linked molecules change places with each other: the link is between the same two molecules
  system.swapMolecules(0, 0, 2);
  EXPECT_TRUE(system.crossLinks.isLinked(site(0, 2, 1), site(0, 0, 1)));
  EXPECT_NEAR(system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms()).potentialEnergy(),
              before.potentialEnergy(), 1e-9);

  // a linked molecule can not be exchanged with a molecule of another system
  System other = makeSystem(forceField, {dimer}, positions);
  other.crossLinks.bondTypes.push_back(type);
  EXPECT_THROW(system.exchangeMolecule(other, 0, 0, 1), std::runtime_error);
  EXPECT_THROW(other.exchangeMolecule(system, 0, 1, 2), std::runtime_error);
  EXPECT_NO_THROW(system.exchangeMolecule(other, 0, 1, 1));
}

// A full Monte Carlo run of charged, flexible dimers with cross-links: the coordinate moves
// (translation, rotation, bead displacement, CBMC reinsertion, grand-canonical swaps), the volume move
// and the two topology moves together keep the running energy consistent with a full recompute,
// including the cross-link slot and the exclusion corrections booked in the pair and Ewald slots.
TEST(MC_CROSS_LINKS, monte_carlo_drift_with_cross_links)
{
  const ForceField forceField = makeForceField(60.0, 0.25, -0.25, true);

  MCMoveProbabilities componentMoves;
  componentMoves.setProbability(Move::Types::Translation, 1.0);
  componentMoves.setProbability(Move::Types::Rotation, 0.5);
  componentMoves.setProbability(Move::Types::BeadDisplacement, 1.0);
  componentMoves.setProbability(Move::Types::ReinsertionCBMC, 0.5);
  componentMoves.setProbability(Move::Types::SwapCBMC, 0.5);

  const std::string json = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["B", [0.0, 0.0, 0.0]], ["A", [1.54, 0.0, 0.0]]],
  "Connectivity" : [[0, 1]],
  "Bonds" : [[["B", "A"], "HARMONIC", [96500.0, 1.54]]],
  "ReactiveSites" : [[1, "X"]],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}
)";
  TemporaryFile file("xl-dimer-mc.json", json);
  Component dimer = Component(Component::Type::Adsorbate, 0, forceField, "xl-dimer-mc", file.stemPath().string(), 5,
                              21, componentMoves, std::nullopt, false);

  MCMoveProbabilities systemMoveProbabilities;
  systemMoveProbabilities.setProbability(Move::Types::CrossLinkSwap, 0.5);
  systemMoveProbabilities.setProbability(Move::Types::CrossLinkFormationScission, 1.0);
  systemMoveProbabilities.setProbability(Move::Types::CrossLinkExchange, 0.5);
  systemMoveProbabilities.setProbability(Move::Types::VolumeChange, 0.05);

  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 5e5, 1.0, {}, {dimer}, {}, {30}, 5,
                         systemMoveProbabilities);
  // Soft link at contact distance so that links actually form and break during the run.
  CrossLinkBondType type = harmonicBondType("X", "X", 100.0, 3.8, 6.5, -2500.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {100.0, 120.0});
  system.crossLinks.bondTypes.push_back(type);
  system.runningEnergies = system.computeTotalEnergies();

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({400, 0, 20, 20, 100000, 1000000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    const RunningEnergy recomputed = s.computeTotalEnergies();
    const RunningEnergy drift = s.runningEnergies - recomputed;
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.crossLink, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.bond, 0.0, 1e-6);

    // The run actually exercised the moves: links were formed, broken and swapped, and every link is
    // consistent with the topology rules.
    const MoveStatistics<double3>& formation =
        std::get<MoveStatistics<double3>>(s.mc_moves_statistics[Move::Types::CrossLinkFormationScission]);
    const MoveStatistics<double>& swap =
        std::get<MoveStatistics<double>>(s.mc_moves_statistics[Move::Types::CrossLinkSwap]);
    const MoveStatistics<double>& exchange =
        std::get<MoveStatistics<double>>(s.mc_moves_statistics[Move::Types::CrossLinkExchange]);
    const MoveStatistics<double3>& translation =
        std::get<MoveStatistics<double3>>(s.components[0].mc_moves_statistics[Move::Types::Translation]);
    EXPECT_GT(translation.totalCounts.x + translation.totalCounts.y + translation.totalCounts.z, 100.0);
    EXPECT_GT(formation.totalCounts.x + formation.totalCounts.y, 100.0);
    EXPECT_GT(formation.totalAccepted.x, 5.0);
    EXPECT_GT(formation.totalAccepted.y, 5.0);
    EXPECT_GT(swap.totalAccepted, 0.0);
    EXPECT_GT(exchange.totalAccepted, 0.0);
    for (const CrossLink& link : s.crossLinks.links)
    {
      EXPECT_FALSE(link.a.sameMolecule(link.b));
      EXPECT_LT(link.a.moleculeIndex, s.numberOfMoleculesPerComponent[0]);
      EXPECT_LT(link.b.moleculeIndex, s.numberOfMoleculesPerComponent[0]);
      EXPECT_LE(s.crossLinks.linkCount(link.a), 1u);
      EXPECT_LE(s.crossLinks.linkCount(link.b), 1u);
    }
  }
}

// CBMC regrowth of a cross-linked molecule: the linked site stays in place and the rest of the
// molecule is regrown with the link's terms (here the junction bend B-A-A') entering the Rosenbluth
// weights through the tethers of the grow context. Two B-A dimers linked at their A atoms, no
// non-bonded interactions: the regrown B bead of each dimer must follow
//   p(r, theta) ~ r^2 sin(theta) exp(-beta [u_bond(r) + u_bend(theta)])
// exactly, with theta the junction angle. The sampled <cos theta> and <u_bend> are compared with the
// quadrature of this density (a soft bend, so that the sin(theta) Jacobian is visible).
namespace
{
// Two B-A dimers linked at their A atoms; 'numberOfMoves' reinsertion moves regrow the B atoms, the junction angle
// must follow its Boltzmann distribution (with the sin(theta) Jacobian) and the running energy must stay consistent.
void expectTetheredRegrowthSamplesJunctionAngle(ForceField forceField, std::size_t numberOfMoves,
                                                double toleranceCos, double toleranceBend, std::string_view label)
{
  constexpr double bendK = 3000.0;      // K/rad^2
  constexpr double bendAngle = 100.0;   // degrees
  constexpr double temperature = 300.0;
  const double theta0 = bendAngle * std::numbers::pi / 180.0;

  // Quadrature of sin(theta) exp(-u_bend / T) over [0, pi].
  double norm = 0.0;
  double expectedCos = 0.0;
  double expectedBend = 0.0;
  constexpr std::size_t n = 20000;
  for (std::size_t k = 0; k <= n; ++k)
  {
    const double theta = std::numbers::pi * static_cast<double>(k) / static_cast<double>(n);
    const double u = 0.5 * bendK * (theta - theta0) * (theta - theta0);
    const double weight =
        ((k == 0 || k == n) ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0)) * std::sin(theta) * std::exp(-u / temperature);
    norm += weight;
    expectedCos += weight * std::cos(theta);
    expectedBend += weight * u;
  }
  expectedCos /= norm;
  expectedBend /= norm;
  // The Jacobian shifts the angle away from the minimum: a sampler that ignored it (or the bend) would
  // sit at cos(theta0) or at <u_bend> = T/2 instead.
  EXPECT_GT(std::abs(expectedCos - std::cos(theta0)), 0.02);

  // Both chain schemes carry the tethers.
  for (const bool useRecoilGrowth : {false, true})
  {
    forceField.settings.useRecoilGrowth = useRecoilGrowth;

    MCMoveProbabilities componentMoves;
    componentMoves.setProbability(Move::Types::ReinsertionCBMC, 1.0);
    TemporaryFile file(std::format("xl-dimer-regrow-{}.json", label), dimerBAJson(kBondK, kBondLength));
    Component dimer = Component(Component::Type::Adsorbate, 0, forceField, std::format("xl-dimer-regrow-{}", label),
                                file.stemPath().string(), 5, 21, componentMoves, std::nullopt, false);

    // Molecule 0: B0 A0, molecule 1: B1 A1; the A atoms (the sites) sit a bond length apart and both
    // junction angles start at theta0. (Recoil growth tests openness against the absolute energy of
    // the step, which includes the partner's junction-bend strain through the tether; a start far up
    // the bend potential would hold the pair in place for thousands of moves and bias the averages.)
    System system = makeSystem(forceField, {dimer},
                               {{double3(9.733, 11.516, 10.0), double3(10.0, 10.0, 10.0),
                                 double3(11.807, 11.516, 10.0), double3(11.54, 10.0, 10.0)}});
    CrossLinkBondType type = harmonicBondType("X", "X", kBondK, kBondLength, 3.0);
    type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {bendK, bendAngle});
    system.crossLinks.bondTypes.push_back(type);
    system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(0, 1, 1), 0});
    system.runningEnergies = system.computeTotalEnergies();
    const double3 siteA0 = system.spanOfMolecule(0, 0)[1].position;
    const double3 siteA1 = system.spanOfMolecule(0, 1)[1].position;

    RandomNumber random(7);
    double sumCos = 0.0;
    double sumBend = 0.0;
    double samples = 0.0;
    const std::size_t burnIn = 5000;
    for (std::size_t i = 0; i < numberOfMoves; ++i)
    {
      const std::size_t molecule = random.uniform_integer(0, 1);
      std::optional<RunningEnergy> energy = MC_Moves::reinsertionMove(random, system, 0, molecule);
      if (energy.has_value()) system.runningEnergies += energy.value();
      if (i < burnIn) continue;  // burn-in

      for (std::size_t m = 0; m < 2; ++m)
      {
        std::span<const Atom> atoms = system.spanOfMolecule(0, m);
        const double3 partner = (m == 0) ? siteA1 : siteA0;
        const double3 u = atoms[0].position - atoms[1].position;
        const double3 v = partner - atoms[1].position;
        const double cosTheta = double3::dot(u, v) / (u.length() * v.length());
        const double theta = std::acos(std::clamp(cosTheta, -1.0, 1.0));
        sumCos += cosTheta;
        sumBend += 0.5 * bendK * (theta - theta0) * (theta - theta0);
        samples += 1.0;
      }
    }

    // The sites never moved, the moves were accepted, and the running energy stayed consistent.
    EXPECT_EQ(system.spanOfMolecule(0, 0)[1].position.x, siteA0.x);
    EXPECT_EQ(system.spanOfMolecule(0, 1)[1].position.x, siteA1.x);
    const MoveStatistics<double>& reinsertion =
        std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::ReinsertionCBMC]);
    EXPECT_GT(reinsertion.totalAccepted, 0.1 * static_cast<double>(numberOfMoves))
        << label << ", recoil growth: " << useRecoilGrowth;
    const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6) << label << ", recoil growth: " << useRecoilGrowth;
    EXPECT_NEAR(drift.crossLink, 0.0, 1e-6);

    EXPECT_NEAR(sumCos / samples, expectedCos, toleranceCos) << label << ", recoil growth: " << useRecoilGrowth;
    EXPECT_NEAR(sumBend / samples, expectedBend, toleranceBend * expectedBend)
        << label << ", recoil growth: " << useRecoilGrowth;
  }
}
}  // namespace

TEST(MC_CROSS_LINKS, tethered_regrowth_samples_junction_angle_distribution)
{
  expectTetheredRegrowthSamplesJunctionAngle(makeForceField(0.0, 0.0, 0.0, false), 105000, 0.01, 0.03, "plain");
}

namespace
{
// Lennard-Jones only between the pairs the link excludes: A-A (the linked sites, at bond distance: 5e6 K) and A-B
// (the 1-3 pairs across the link, about 2.5 A apart: 1.6e4 K). B-B, the 1-4 pair across the link that is not
// excluded, is zero, so the junction angle keeps the distribution of the bend alone.
ForceField makeExcludedPairsOnlyForceField()
{
  ForceField forceField = makeForceField(100.0, 0.0, 0.0, false);
  forceField(1, 1) = VDWParameters(0.0, 3.6);
  forceField(0, 1) = VDWParameters(100.0, 3.5);
  forceField(1, 0) = VDWParameters(100.0, 3.5);
  forceField.preComputeCutOffDependentParameters();
  return forceField;
}
}  // namespace

// The dual cut-off scheme: the linked site atom of the regrown molecule sits at bond distance from its partner,
// far above the overlap criterion, and the kept site atoms are not part of the regrowth; the correction to the
// full cut-offs must neither reject every grow nor throw on the retrace, and must leave the sampling exact.
TEST(MC_CROSS_LINKS, tethered_regrowth_with_dual_cut_off)
{
  ForceField forceField = makeExcludedPairsOnlyForceField();
  forceField.settings.useDualCutOff = true;
  forceField.settings.dualCutOff = 6.0;
  expectTetheredRegrowthSamplesJunctionAngle(forceField, 45000, 0.015, 0.05, "dual-cut-off");
}

// An overlap criterion below the non-bonded energy of the 1-3 pair across the link (the regrown B with the
// partner's site A, about 2.5 A apart: 1.6e4 K here): the pair is excluded by the link, so it is no overlap.
TEST(MC_CROSS_LINKS, tethered_regrowth_with_excluded_pair_above_overlap_criterion)
{
  ForceField forceField = makeExcludedPairsOnlyForceField();
  forceField.energyOverlapCriteria = 5000.0;
  expectTetheredRegrowthSamplesJunctionAngle(forceField, 45000, 0.015, 0.05, "overlap-criterion");
}

// Fixed-endpoint regrowth of cross-linked chains: two B-A-A-B chains linked end to end at both ends
// (a ring of two chains) have both terminal atoms fixed by links, so the reinsertion regrows the two
// interior atoms with the bridge-closure step, and the partial reinsertion adds its own fixed atoms.
// The running energy stays consistent and the moves are accepted.
TEST(MC_CROSS_LINKS, fixed_endpoint_regrowth_of_doubly_linked_chains)
{
  const ForceField forceField = makeForceField(40.0, 0.0, 0.0, false);

  MCMoveProbabilities componentMoves;
  componentMoves.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentMoves.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);
  const std::string json = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["B", [0.0, 0.0, 0.0]], ["A", [1.54, 0.0, 0.0]], ["A", [3.08, 0.0, 0.0]], ["B", [4.62, 0.0, 0.0]]],
  "Connectivity" : [[0, 1], [1, 2], [2, 3]],
  "Bonds" : [[["B", "A"], "HARMONIC", [96500.0, 1.54]], [["A", "A"], "HARMONIC", [96500.0, 1.54]]],
  "Bends" : [[["B", "A", "A"], "HARMONIC", [62500.0, 114.0]]],
  "Torsions" : [[["B", "A", "A", "B"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]],
  "ReactiveSites" : [[0, "X"], [3, "X"]],
  "Partial-reinsertion" : [[0, 1], [2, 3]],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}
)";
  TemporaryFile file("xl-chain-ring.json", json);
  Component chain = Component(Component::Type::Adsorbate, 0, forceField, "xl-chain-ring", file.stemPath().string(), 5,
                              21, componentMoves, std::nullopt, false);

  // Two trans zigzag chains (1.54 A bonds, 114 degree bends) stacked 3.8 A apart in z, linked B0-B0'
  // and B3-B3' (the junction angles A-B-B' are 90 degrees).
  const std::vector<double3> zigzag{double3(10.0, 10.0, 10.0), double3(11.54, 10.0, 10.0), double3(12.166, 11.407, 10.0),
                                    double3(13.706, 11.407, 10.0)};
  std::vector<double3> positions = zigzag;
  for (const double3& p : zigzag) positions.push_back(p + double3(0.0, 0.0, 3.8));
  System system = makeSystem(forceField, {chain}, {positions});
  CrossLinkBondType type = harmonicBondType("X", "X", 1000.0, 3.8, 6.0);
  type.junctionBend = BendPotential({0, 1, 2}, BendType::Harmonic, {2000.0, 90.0});
  system.crossLinks.bondTypes.push_back(type);
  system.crossLinks.addLink(CrossLink{site(0, 0, 0), site(0, 1, 0), 0});
  system.crossLinks.addLink(CrossLink{site(0, 0, 3), site(0, 1, 3), 0});
  system.runningEnergies = system.computeTotalEnergies();

  // Both ends are linked: the whole-molecule regrowth keeps atoms 0 and 3, the partial one 0,1 or 2,3
  // plus the linked ends.
  EXPECT_EQ(system.crossLinkedSiteAtoms(0, 0), (std::vector<std::size_t>{0, 3}));
  EXPECT_EQ(system.crossLinkRegrowthPlacedSet(0, 0).value(), (std::vector<std::size_t>{0, 3}));
  EXPECT_EQ(system.crossLinkRegrowthPlacedSet(0, 0, std::vector<std::size_t>{0, 1}).value(),
            (std::vector<std::size_t>{0, 1, 3}));
  EXPECT_FALSE(system.crossLinkRegrowthPlacedSet(0, 0, std::vector<std::size_t>{1, 2}).has_value());
  EXPECT_EQ(system.crossLinkTethers(0, 1).size(), 2uz);

  const std::vector<double3> endsBefore{system.spanOfMolecule(0, 0)[0].position, system.spanOfMolecule(0, 0)[3].position,
                                        system.spanOfMolecule(0, 1)[0].position, system.spanOfMolecule(0, 1)[3].position};

  RandomNumber random(11);
  for (std::size_t i = 0; i < 20000; ++i)
  {
    const std::size_t molecule = random.uniform_integer(0, 1);
    std::optional<RunningEnergy> energy = (i % 2 == 0) ? MC_Moves::reinsertionMove(random, system, 0, molecule)
                                                       : MC_Moves::partialReinsertionMove(random, system, 0, molecule);
    if (energy.has_value()) system.runningEnergies += energy.value();
  }

  const std::vector<double3> endsAfter{system.spanOfMolecule(0, 0)[0].position, system.spanOfMolecule(0, 0)[3].position,
                                       system.spanOfMolecule(0, 1)[0].position, system.spanOfMolecule(0, 1)[3].position};
  for (std::size_t k = 0; k < 4; ++k)
  {
    EXPECT_EQ(endsAfter[k].x, endsBefore[k].x);
    EXPECT_EQ(endsAfter[k].y, endsBefore[k].y);
    EXPECT_EQ(endsAfter[k].z, endsBefore[k].z);
  }
  const MoveStatistics<double>& reinsertion =
      std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::ReinsertionCBMC]);
  const MoveStatistics<double>& partial =
      std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::PartialReinsertionCBMC]);
  EXPECT_GT(reinsertion.totalAccepted, 100.0);
  EXPECT_GT(partial.totalAccepted, 100.0);

  const RunningEnergy recomputed = system.computeTotalEnergies();
  const RunningEnergy drift = system.runningEnergies - recomputed;
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  EXPECT_NEAR(drift.crossLink, 0.0, 1e-6);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
  EXPECT_NEAR(drift.bond, 0.0, 1e-6);
  EXPECT_NEAR(drift.bend, 0.0, 1e-6);
  EXPECT_NEAR(drift.torsion, 0.0, 1e-6);
  // The interior atoms actually moved.
  const std::span<const Atom> atoms = system.spanOfMolecule(0, 0);
  EXPECT_GT((atoms[1].position - zigzag[1]).length() + (atoms[2].position - zigzag[2]).length(), 0.05);
}

// Several links can exclude the same atom pair (links between adjacent sites of the same two molecules), while the
// inter-molecular sum counted that pair once: every pair is excluded once, however many links reach it.
namespace
{
// Linear molecule of 'numberOfAtoms' beads of 'atomName' with the given reactive atoms, all of site type 'siteType'.
std::string multiSiteChainJson(std::string_view atomName, std::size_t numberOfAtoms,
                               std::vector<std::size_t> reactiveAtoms, std::string_view siteType, std::size_t valence)
{
  std::string atoms{};
  std::string connectivity{};
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    atoms += std::format("{}[\"{}\", [{}, 0.0, 0.0]]", i == 0 ? "" : ", ", atomName, static_cast<double>(i) * kBondLength);
    if (i > 0) connectivity += std::format("{}[{}, {}]", i == 1 ? "" : ", ", i - 1, i);
  }
  std::string sites{};
  for (std::size_t i = 0; i < reactiveAtoms.size(); ++i)
  {
    sites += std::format("{}{{\"Atom\": {}, \"Type\": \"{}\", \"Valence\": {}}}", i == 0 ? "" : ", ", reactiveAtoms[i],
                         siteType, valence);
  }
  std::string topology{};
  if (numberOfAtoms > 1)
  {
    topology = std::format(R"("Connectivity" : [{0}],
  "Bonds" : [[["{1}", "{1}"], "HARMONIC", [{2}, {3}]]],
  {4}"VanDerWaals" : "auto",
  "Coulomb" : "auto",
)",
                           connectivity, atomName, kBondK, kBondLength,
                           numberOfAtoms > 2 ? std::format("\"Bends\" : [[[\"{0}\", \"{0}\", \"{0}\"], \"HARMONIC\", [0.0, 114.0]]],\n  ", atomName)
                                             : std::string{});
  }
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [{}],
  {}"ReactiveSites" : [{}]
}}
)",
                     atoms, topology, sites);
}

// A closed ring of the given atoms in which every bond is the same harmonic bond: the reference of a linked network.
std::string ringJson(std::vector<std::string> atomNames)
{
  std::string atoms{};
  std::string connectivity{};
  std::set<std::pair<std::string, std::string>> bondTypes{};
  for (std::size_t i = 0; i < atomNames.size(); ++i)
  {
    const std::size_t next = (i + 1) % atomNames.size();
    atoms += std::format("{}[\"{}\", [{}, 0.0, 0.0]]", i == 0 ? "" : ", ", atomNames[i], static_cast<double>(i) * kBondLength);
    connectivity += std::format("{}[{}, {}]", i == 0 ? "" : ", ", i, next);
    bondTypes.insert(std::minmax(atomNames[i], atomNames[next]));
  }
  std::string bonds{};
  for (const auto& [a, b] : bondTypes)
  {
    bonds += std::format("{}[[\"{}\", \"{}\"], \"HARMONIC\", [{}, {}]]", bonds.empty() ? "" : ", ", a, b, kBondK, kBondLength);
  }
  // the bends and torsions of a ring have to be defined; they are zero, as a link has none here
  const std::size_t n = atomNames.size();
  std::set<std::vector<std::string>> bendTypes{};
  std::set<std::vector<std::string>> torsionTypes{};
  for (std::size_t i = 0; i < n; ++i)
  {
    // a definition matches in both directions
    auto canonical = [](std::vector<std::string> names)
    {
      std::vector<std::string> reversed(names.rbegin(), names.rend());
      return std::min(names, reversed);
    };
    bendTypes.insert(canonical({atomNames[i], atomNames[(i + 1) % n], atomNames[(i + 2) % n]}));
    torsionTypes.insert(
        canonical({atomNames[i], atomNames[(i + 1) % n], atomNames[(i + 2) % n], atomNames[(i + 3) % n]}));
  }
  std::string bends{};
  for (const std::vector<std::string>& t : bendTypes)
  {
    bends += std::format("{}[[\"{}\", \"{}\", \"{}\"], \"HARMONIC\", [0.0, 114.0]]", bends.empty() ? "" : ", ", t[0], t[1], t[2]);
  }
  std::string torsions{};
  for (const std::vector<std::string>& t : torsionTypes)
  {
    torsions += std::format("{}[[\"{}\", \"{}\", \"{}\", \"{}\"], \"TRAPPE\", [0.0, 0.0, 0.0, 0.0]]",
                            torsions.empty() ? "" : ", ", t[0], t[1], t[2], t[3]);
  }
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [{}],
  "Connectivity" : [{}],
  "Bonds" : [{}],
  "Bends" : [{}],
  "Torsions" : [{}],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}}
)",
                     atoms, connectivity, bonds, bends, torsions);
}

void expectSameEnergies(const RunningEnergy& linked, const RunningEnergy& reference)
{
  const double tolerance = 1e-8 * std::max(1.0, std::abs(reference.potentialEnergy()));
  EXPECT_NEAR(linked.potentialEnergy(), reference.potentialEnergy(), tolerance);
  EXPECT_NEAR(linked.moleculeMoleculeVDW, 0.0, 1e-8);
  EXPECT_NEAR(linked.moleculeMoleculeCharge, 0.0, 1e-8);
  EXPECT_NEAR(linked.ewald_exclusion, reference.ewald_exclusion, tolerance);
  EXPECT_NEAR(linked.crossLink + linked.bond, reference.bond + reference.bend, tolerance);
}

// Two dimers A1-A2 and B1-B2, linked A1-B1 and A2-B2: the 1-3 pairs (A2,B1) and (A1,B2) are reached by both links.
System makeLadder(const ForceField& forceField, const std::string& name)
{
  Component dimerA = makeComponent(forceField, 0, name + "-a", multiSiteChainJson("A", 2, {0, 1}, "X", 1));
  Component dimerB = makeComponent(forceField, 1, name + "-b", multiSiteChainJson("B", 2, {0, 1}, "Y", 1));
  System system = makeSystem(forceField, {dimerA, dimerB},
                             {{double3(10.0, 10.0, 10.0), double3(11.5, 10.2, 10.1)},
                              {double3(10.3, 11.6, 10.4), double3(11.8, 11.9, 10.9)}});
  system.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 6.0));
  system.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  system.crossLinks.addLink(CrossLink{site(0, 0, 1), site(1, 0, 1), 0});
  return system;
}
}  // namespace

TEST(MC_CROSS_LINKS, ladder_of_two_links_equals_four_membered_ring)
{
  const ForceField forceField = makeForceField(100.0, 0.4, -0.4, true);
  const System linked = makeLadder(forceField, "xl-ladder");
  const RunningEnergy linkedEnergy = const_cast<System&>(linked).computeTotalEnergies();

  const std::span<const Atom> atoms = linked.spanOfMoleculeAtoms();
  Component ring = makeComponent(forceField, 0, "xl-ring-aabb", ringJson({"A", "A", "B", "B"}));
  System reference = makeSystem(forceField, {ring},
                                {{atoms[0].position, atoms[1].position, atoms[3].position, atoms[2].position}});
  expectSameEnergies(linkedEnergy, reference.computeTotalEnergies());
}

TEST(MC_CROSS_LINKS, two_links_to_one_site_from_adjacent_sites_equal_three_membered_ring)
{
  const ForceField forceField = makeForceField(100.0, 0.4, -0.4, true);
  Component dimer = makeComponent(forceField, 0, "xl-shared-a", multiSiteChainJson("A", 2, {0, 1}, "X", 1));
  Component monomer = makeComponent(forceField, 1, "xl-shared-b", multiSiteChainJson("B", 1, {0}, "Y", 2));
  const double3 a1(10.0, 10.0, 10.0);
  const double3 a2(11.5, 10.2, 10.1);
  const double3 b(10.6, 11.4, 10.5);
  System linked = makeSystem(forceField, {dimer, monomer}, {{a1, a2}, {b}});
  linked.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 6.0));
  linked.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  linked.crossLinks.addLink(CrossLink{site(0, 0, 1), site(1, 0, 0), 0});

  Component ring = makeComponent(forceField, 0, "xl-ring-aab", ringJson({"A", "A", "B"}));
  System reference = makeSystem(forceField, {ring}, {{a1, a2, b}});
  expectSameEnergies(linked.computeTotalEnergies(), reference.computeTotalEnergies());
}

TEST(MC_CROSS_LINKS, two_links_sharing_a_neighbour_of_their_sites_equal_four_membered_ring)
{
  // the tethers are compared with the total energy below: no excluded pair may be above the overlap criterion
  // (a tether leaves such a pair out of the exclusion, as the regrowth leaves it out of the external energy)
  ForceField forceField = makeForceField(100.0, 0.4, -0.4, true);
  forceField.energyOverlapCriteria = 1e300;
  Component trimer = makeComponent(forceField, 0, "xl-bridge-a", multiSiteChainJson("A", 3, {0, 2}, "X", 1));
  Component monomer = makeComponent(forceField, 1, "xl-bridge-b", multiSiteChainJson("B", 1, {0}, "Y", 2));
  const double3 a0(10.0, 10.0, 10.0);
  const double3 a1(11.4, 10.3, 10.2);
  const double3 a2(11.9, 11.7, 10.6);
  const double3 b(10.4, 11.5, 10.9);
  System linked = makeSystem(forceField, {trimer, monomer}, {{a0, a1, a2}, {b}});
  linked.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 6.0));
  linked.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  linked.crossLinks.addLink(CrossLink{site(0, 0, 2), site(1, 0, 0), 0});

  Component ring = makeComponent(forceField, 0, "xl-ring-aaab", ringJson({"A", "A", "A", "B"}));
  System reference = makeSystem(forceField, {ring}, {{a0, a1, a2, b}});
  expectSameEnergies(linked.computeTotalEnergies(), reference.computeTotalEnergies());

  // the regrown middle atom is reached by both links: the tethers of the trimer hold its pair with b once
  const std::vector<CrossLinkTether> tethers = linked.crossLinkTethers(0, 0);
  ASSERT_EQ(tethers.size(), 2uz);
  EXPECT_EQ(tethers[0].exclusions.size() + tethers[1].exclusions.size(), 3uz);
  const std::array<std::size_t, 2> all{0, 1};
  const RunningEnergy tetherEnergy = Interactions::computeCrossLinkTetherEnergy(
      forceField, linked.simulationBox, linked.components[0], linked.spanOfMolecule(0, 0), tethers, all);
  const RunningEnergy total = linked.computeCrossLinkEnergy(linked.simulationBox, linked.spanOfMoleculeAtoms());
  EXPECT_NEAR(tetherEnergy.potentialEnergy(), total.potentialEnergy(), 1e-8);
}

TEST(MC_CROSS_LINKS, shared_exclusions_gradient_strain_difference_and_tethers_are_consistent)
{
  // the tethers are compared with the total energy below: no excluded pair may be above the overlap criterion
  ForceField forceField = makeForceField(100.0, 0.4, -0.4, true);
  forceField.energyOverlapCriteria = 1e300;
  System system = makeLadder(forceField, "xl-ladder-derivatives");
  const std::span<Atom> atoms = system.spanOfMoleculeAtoms();

  auto total = [&](std::span<const Atom> positions)
  {
    return Interactions::computeCrossLinkEnergy(forceField, system.simulationBox, system.components,
                                                system.numberOfMoleculesPerComponent, positions, system.crossLinks);
  };

  std::vector<AtomDynamics> dynamics(atoms.size());
  const RunningEnergy energy =
      Interactions::computeCrossLinkGradient(forceField, system.simulationBox, system.components,
                                             system.numberOfMoleculesPerComponent, atoms, dynamics, system.crossLinks);
  EXPECT_NEAR(energy.potentialEnergy(), total(atoms).potentialEnergy(), 1e-10);
  std::vector<AtomDynamics> strainDynamics(atoms.size());
  const auto [strainEnergy, strain] = Interactions::computeCrossLinkEnergyStrainDerivative(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms, strainDynamics,
      system.crossLinks);
  EXPECT_NEAR(energy.potentialEnergy(), strainEnergy.potentialEnergy(), 1e-10);

  const double h = 1e-5;
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    EXPECT_NEAR((dynamics[i].gradient - strainDynamics[i].gradient).length(), 0.0, 1e-10);
    for (std::size_t k = 0; k < 3; ++k)
    {
      auto shifted = [&](double delta)
      {
        std::vector<Atom> copy(atoms.begin(), atoms.end());
        double3 d{};
        (k == 0 ? d.x : (k == 1 ? d.y : d.z)) = delta;
        copy[i].position += d;
        return total(copy).potentialEnergy();
      };
      const double numerical = (shifted(h) - shifted(-h)) / (2.0 * h);
      const double analytic =
          (k == 0 ? dynamics[i].gradient.x : (k == 1 ? dynamics[i].gradient.y : dynamics[i].gradient.z));
      EXPECT_NEAR(analytic, numerical, 1e-5 * std::max(1.0, std::abs(numerical))) << "atom " << i << " axis " << k;
    }
  }

  // The energy difference of a moved molecule (all links that can share a pair are links of that molecule).
  const RunningEnergy before = total(atoms);
  std::vector<Atom> oldAtoms(atoms.begin(), atoms.begin() + 2);
  std::vector<Atom> newAtoms = oldAtoms;
  newAtoms[0].position += double3(0.11, -0.07, 0.05);
  newAtoms[1].position += double3(-0.06, 0.09, -0.04);
  const RunningEnergy difference = Interactions::computeCrossLinkEnergyDifference(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms,
      system.crossLinks, 0, 0, newAtoms, oldAtoms);
  std::vector<Atom> movedAtoms(atoms.begin(), atoms.end());
  movedAtoms[0] = newAtoms[0];
  movedAtoms[1] = newAtoms[1];
  EXPECT_NEAR(difference.potentialEnergy(), total(movedAtoms).potentialEnergy() - before.potentialEnergy(), 1e-8);

  // The tethers of a molecule hold each excluded pair once, also when only some of them are evaluated at a time.
  const std::vector<CrossLinkTether> tethers = system.crossLinkTethers(0, 0);
  ASSERT_EQ(tethers.size(), 2uz);
  const std::span<const Atom> molecule = system.spanOfMolecule(0, 0);
  const std::array<std::size_t, 1> first{0};
  const std::array<std::size_t, 1> second{1};
  const RunningEnergy tetherEnergy = Interactions::computeCrossLinkTetherEnergy(
      forceField, system.simulationBox, system.components[0], molecule, tethers, first) +
                                     Interactions::computeCrossLinkTetherEnergy(forceField, system.simulationBox,
                                                                                system.components[0], molecule, tethers,
                                                                                second);
  EXPECT_NEAR(tetherEnergy.potentialEnergy(), before.potentialEnergy(), 1e-8);
}

// The topology moves score the change of the set of excluded pairs: breaking one link of a ladder keeps the pairs
// the other link still excludes. Frozen positions: the distribution over the 7 topologies of two dimers is exact.
TEST(MC_CROSS_LINKS, topology_moves_with_shared_exclusions_sample_exact_distribution)
{
  // a small epsilon: the pairs at these short distances are then of the order of the thermal energy
  const ForceField forceField = makeForceField(0.01, 0.0, 0.0, false);
  System system = makeLadder(forceField, "xl-ladder-moves");
  system.crossLinks.links.clear();
  system.crossLinks.rebuildIndex();
  system.crossLinks.bondTypes[0] = harmonicBondType("X", "Y", 800.0, 1.8, 6.0, -250.0);

  using Key = std::vector<std::pair<std::uint32_t, std::uint32_t>>;
  auto keyOf = [](const std::vector<CrossLink>& links)
  {
    Key key{};
    for (const CrossLink& link : links)
    {
      key.push_back(link.a.componentId == 0 ? std::make_pair(link.a.atomIndex, link.b.atomIndex)
                                            : std::make_pair(link.b.atomIndex, link.a.atomIndex));
    }
    std::ranges::sort(key);
    return key;
  };
  const std::vector<Key> topologies{{},
                                    {{0, 0}},
                                    {{0, 1}},
                                    {{1, 0}},
                                    {{1, 1}},
                                    {{0, 0}, {1, 1}},
                                    {{0, 1}, {1, 0}}};
  std::vector<double> logWeights{};
  std::vector<double> energies{};
  for (const Key& topology : topologies)
  {
    std::vector<CrossLink> links{};
    for (const auto& [i, j] : topology) links.push_back(CrossLink{site(0, 0, i), site(1, 0, j), 0});
    const RunningEnergy energy = Interactions::computeCrossLinkEnergyOfLinks(
        forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms(), system.crossLinks, links);
    energies.push_back(energy.potentialEnergy());
    logWeights.push_back(-system.beta * energy.potentialEnergy());
  }
  // the two links of a ladder share pairs: the energy of the pair of links is not the sum of the single links
  EXPECT_GT(system.beta * std::abs(energies[5] - energies[1] - energies[4]), 0.2);
  const double maxLog = *std::ranges::max_element(logWeights);
  std::vector<double> weights{};
  for (double logWeight : logWeights) weights.push_back(std::exp(logWeight - maxLog));
  const double normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  RandomNumber random(2024);
  std::map<Key, std::size_t> counts{};
  const std::size_t numberOfMoves = 600000;
  RunningEnergy running = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  for (std::size_t step = 0; step < numberOfMoves; ++step)
  {
    const double choice = random.uniform();
    std::optional<RunningEnergy> difference;
    if (choice < 0.5)
      difference = MC_Moves::crossLinkFormationScissionMove(random, system);
    else if (choice < 0.8)
      difference = MC_Moves::crossLinkSwapMove(random, system);
    else
      difference = MC_Moves::crossLinkExchangeMove(random, system);
    if (difference) running += difference.value();
    counts[keyOf(system.crossLinks.links)] += 1;
  }
  const RunningEnergy recomputed = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  EXPECT_NEAR(running.potentialEnergy(), recomputed.potentialEnergy(),
              1e-8 * std::max(1.0, std::abs(recomputed.potentialEnergy())));

  double largest = 0.0;
  for (std::size_t t = 0; t < topologies.size(); ++t)
  {
    const double exact = weights[t] / normalization;
    largest = std::max(largest, exact);
    const double sampled = static_cast<double>(counts[topologies[t]]) / static_cast<double>(numberOfMoves);
    EXPECT_NEAR(sampled, exact, 0.015) << "topology " << t;
  }
  // the distribution is spread over the topologies, also the ones whose links share pairs
  EXPECT_LT(largest, 0.8);
  EXPECT_GT(weights[5] / normalization + weights[6] / normalization, 0.02);
}

// A site that carries two links: the two partners form a 1-3 pair through the site and the bend of the site, as the
// two bonds of an atom in a molecule do.
namespace
{
// Single-bead monomer with one reactive site of valence 2 and the bend between its two links.
std::string valenceTwoMonomerJson(std::string_view atomName, std::string_view siteType, double bendK, double bendAngle)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["{}", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [{{"Atom": 0, "Type": "{}", "Valence": 2, "LinkBend": ["HARMONIC", [{}, {}]]}}]
}}
)",
                     atomName, siteType, bendK, bendAngle);
}

// The molecule A-B-A with harmonic bonds and the bend at B: the reference for two monomers A linked to one B.
std::string chainABAJson(double bondK, double bondLength, double bendK, double bendAngle)
{
  return std::format(R"({{
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["A", [0.0, 0.0, 0.0]], ["B", [{1}, 0.0, 0.0]], ["A", [{4}, 0.0, 0.0]]],
  "Connectivity" : [[0, 1], [1, 2]],
  "Bonds" : [[["A", "B"], "HARMONIC", [{0}, {1}]]],
  "Bends" : [[["A", "B", "A"], "HARMONIC", [{2}, {3}]]],
  "VanDerWaals" : "auto",
  "Coulomb" : "auto"
}}
)",
                     bondK, bondLength, bendK, bendAngle, 2.0 * bondLength);
}

// Two monomers A (valence 1) linked to one monomer B (valence 2, with a link bend), in a bent geometry.
System makeValenceTwoChain(const ForceField& forceField, const std::string& name, double bendK = kBendK,
                           std::vector<double3> positionsA = {double3(10.0, 10.0, 10.0), double3(12.4, 11.4, 10.5)},
                           double3 positionB = double3(11.5, 10.2, 10.1))
{
  Component monomerA = makeComponent(forceField, 0, name + "-a", monomerJson("A", "X"));
  Component monomerB = makeComponent(forceField, 1, name + "-b", valenceTwoMonomerJson("B", "Y", bendK, kBendAngle));
  System system = makeSystem(forceField, {monomerA, monomerB}, {std::move(positionsA), {positionB}});
  system.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 6.0));
  return system;
}
}  // namespace

TEST(MC_CROSS_LINKS, two_links_at_one_site_equal_bonded_chain_with_bend)
{
  const ForceField forceField = makeForceField(100.0, 0.2, -0.4, true);
  System linked = makeValenceTwoChain(forceField, "xl-valence-two");
  linked.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  linked.crossLinks.addLink(CrossLink{site(0, 1, 0), site(1, 0, 0), 0});
  const RunningEnergy linkedEnergy = linked.computeTotalEnergies();

  const std::span<const Atom> atoms = linked.spanOfMoleculeAtoms();
  Component chain = makeComponent(forceField, 0, "xl-chain-aba", chainABAJson(kBondK, kBondLength, kBendK, kBendAngle));
  System reference = makeSystem(forceField, {chain}, {{atoms[0].position, atoms[2].position, atoms[1].position}});
  const RunningEnergy referenceEnergy = reference.computeTotalEnergies();

  expectSameEnergies(linkedEnergy, referenceEnergy);
  EXPECT_NEAR(linkedEnergy.crossLink, referenceEnergy.bond + referenceEnergy.bend, 1e-8);
  EXPECT_GT(referenceEnergy.bend, 1.0);
}

TEST(MC_CROSS_LINKS, two_links_at_one_site_derivatives_and_differences_are_consistent)
{
  const ForceField forceField = makeForceField(100.0, 0.2, -0.4, true);
  System system = makeValenceTwoChain(forceField, "xl-valence-two-derivatives");
  system.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  system.crossLinks.addLink(CrossLink{site(0, 1, 0), site(1, 0, 0), 0});
  const std::span<Atom> atoms = system.spanOfMoleculeAtoms();

  auto total = [&](std::span<const Atom> positions)
  {
    return Interactions::computeCrossLinkEnergy(forceField, system.simulationBox, system.components,
                                                system.numberOfMoleculesPerComponent, positions, system.crossLinks);
  };

  std::vector<AtomDynamics> dynamics(atoms.size());
  const RunningEnergy energy =
      Interactions::computeCrossLinkGradient(forceField, system.simulationBox, system.components,
                                             system.numberOfMoleculesPerComponent, atoms, dynamics, system.crossLinks);
  EXPECT_NEAR(energy.potentialEnergy(), total(atoms).potentialEnergy(), 1e-10);
  std::vector<AtomDynamics> strainDynamics(atoms.size());
  const auto [strainEnergy, strain] = Interactions::computeCrossLinkEnergyStrainDerivative(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms, strainDynamics,
      system.crossLinks);
  EXPECT_NEAR(energy.potentialEnergy(), strainEnergy.potentialEnergy(), 1e-10);

  const double h = 1e-5;
  double3 totalGradient{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    totalGradient += dynamics[i].gradient;
    EXPECT_NEAR((dynamics[i].gradient - strainDynamics[i].gradient).length(), 0.0, 1e-10);
    for (std::size_t k = 0; k < 3; ++k)
    {
      auto shifted = [&](double delta)
      {
        std::vector<Atom> copy(atoms.begin(), atoms.end());
        double3 d{};
        (k == 0 ? d.x : (k == 1 ? d.y : d.z)) = delta;
        copy[i].position += d;
        return total(copy).potentialEnergy();
      };
      const double numerical = (shifted(h) - shifted(-h)) / (2.0 * h);
      const double analytic =
          (k == 0 ? dynamics[i].gradient.x : (k == 1 ? dynamics[i].gradient.y : dynamics[i].gradient.z));
      EXPECT_NEAR(analytic, numerical, 1e-5 * std::max(1.0, std::abs(numerical))) << "atom " << i << " axis " << k;
    }
  }
  EXPECT_NEAR(totalGradient.length(), 0.0, 1e-6);

  // the strain derivative against a finite difference of the energy under a homogeneous deformation
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      auto deformed = [&](double delta)
      {
        double3x3 deformation = double3x3::identity();
        deformation[row][column] += delta;
        const SimulationBox box(deformation * system.simulationBox.cell);
        std::vector<Atom> copy(atoms.begin(), atoms.end());
        for (Atom& atom : copy) atom.position = deformation * atom.position;
        return Interactions::computeCrossLinkEnergy(forceField, box, system.components,
                                                    system.numberOfMoleculesPerComponent, copy, system.crossLinks)
            .potentialEnergy();
      };
      const double numerical = (deformed(h) - deformed(-h)) / (2.0 * h);
      EXPECT_NEAR(strain[row][column], numerical, 1e-5 * std::max(1.0, std::abs(numerical)))
          << "strain " << row << " " << column;
    }
  }

  // The energy difference of a moved molecule: the pair of the two partners and the bend at the site involve the
  // other partner, whose link is not a link of the moved molecule.
  const RunningEnergy before = total(atoms);
  for (const auto [component, molecule, first] : {std::tuple{0uz, 1uz, 1uz}, std::tuple{1uz, 0uz, 2uz}})
  {
    std::vector<Atom> oldAtoms{atoms[first]};
    std::vector<Atom> newAtoms = oldAtoms;
    newAtoms[0].position += double3(0.13, -0.09, 0.07);
    const RunningEnergy difference = Interactions::computeCrossLinkEnergyDifference(
        forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms,
        system.crossLinks, component, molecule, newAtoms, oldAtoms);
    std::vector<Atom> movedAtoms(atoms.begin(), atoms.end());
    movedAtoms[first] = newAtoms[0];
    EXPECT_NEAR(difference.potentialEnergy(), total(movedAtoms).potentialEnergy() - before.potentialEnergy(), 1e-8)
        << "component " << component << " molecule " << molecule;
  }

  // The topology change: adding the second link adds the terms it forms with the first one.
  const CrossLink second = system.crossLinks.links[1];
  system.crossLinks.removeLink(1);
  const RunningEnergy change = Interactions::computeCrossLinkTopologyChangeEnergy(
      forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atoms,
      system.crossLinks, {}, std::span<const CrossLink>(&second, 1));
  EXPECT_NEAR(change.potentialEnergy(), before.potentialEnergy() - total(atoms).potentialEnergy(), 1e-8);
}

// At frozen positions the topology moves sample the exact distribution over the four topologies of two valence-1
// monomers A and one valence-2 monomer B; the doubly linked topology carries the bend and the A-A exclusion.
TEST(MC_CROSS_LINKS, topology_moves_with_two_links_at_one_site_sample_exact_distribution)
{
  // the monomers at the Lennard-Jones contact distance of the bond, where the pair energies are thermal
  const ForceField forceField = makeForceField(100.0, 0.0, 0.0, false);
  const double3 b(10.0, 10.0, 10.0);
  const double angle = 100.0 * Units::DegreesToRadians;
  System system = makeValenceTwoChain(forceField, "xl-valence-two-moves", 2000.0,
                                      {b + 3.4 * double3(1.0, 0.0, 0.05), b + 3.4 * double3(std::cos(angle), std::sin(angle), -0.05)}, b);
  system.crossLinks.bondTypes[0] = harmonicBondType("X", "Y", 800.0, 3.4, 6.0, -250.0);

  using Key = std::vector<std::uint32_t>;  // the molecules of component A linked to B
  auto keyOf = [](const std::vector<CrossLink>& links)
  {
    Key key{};
    for (const CrossLink& link : links) key.push_back(link.a.componentId == 0 ? link.a.moleculeIndex : link.b.moleculeIndex);
    std::ranges::sort(key);
    return key;
  };
  const std::vector<Key> topologies{{}, {0}, {1}, {0, 1}};
  std::vector<double> energies{};
  for (const Key& topology : topologies)
  {
    std::vector<CrossLink> links{};
    for (std::uint32_t m : topology) links.push_back(CrossLink{site(0, m, 0), site(1, 0, 0), 0});
    energies.push_back(Interactions::computeCrossLinkEnergyOfLinks(
                           forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
                           system.spanOfMoleculeAtoms(), system.crossLinks, links)
                           .potentialEnergy());
  }
  // the two links at B form a bend and exclude the A-A pair: the energy is not the sum of the single links
  EXPECT_GT(system.beta * std::abs(energies[3] - energies[1] - energies[2]), 0.2);
  std::vector<double> weights{};
  for (double energy : energies) weights.push_back(std::exp(-system.beta * (energy - energies[3])));
  const double normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  RandomNumber random(4711);
  std::map<Key, std::size_t> counts{};
  const std::size_t numberOfMoves = 400000;
  RunningEnergy running = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  for (std::size_t step = 0; step < numberOfMoves; ++step)
  {
    const double choice = random.uniform();
    std::optional<RunningEnergy> difference;
    if (choice < 0.6)
      difference = MC_Moves::crossLinkFormationScissionMove(random, system);
    else if (choice < 0.8)
      difference = MC_Moves::crossLinkSwapMove(random, system);
    else
      difference = MC_Moves::crossLinkExchangeMove(random, system);
    if (difference) running += difference.value();
    counts[keyOf(system.crossLinks.links)] += 1;
    for (const CrossLink& link : system.crossLinks.links)
    {
      ASSERT_LE(system.crossLinks.linkCount(link.a), 2u);
      ASSERT_LE(system.crossLinks.linkCount(link.b), 2u);
    }
  }
  const RunningEnergy recomputed = system.computeCrossLinkEnergy(system.simulationBox, system.spanOfMoleculeAtoms());
  EXPECT_NEAR(running.potentialEnergy(), recomputed.potentialEnergy(),
              1e-8 * std::max(1.0, std::abs(recomputed.potentialEnergy())));

  double largest = 0.0;
  for (std::size_t t = 0; t < topologies.size(); ++t)
  {
    const double exact = weights[t] / normalization;
    largest = std::max(largest, exact);
    const double sampled = static_cast<double>(counts[topologies[t]]) / static_cast<double>(numberOfMoves);
    EXPECT_NEAR(sampled, exact, 0.015) << "topology " << t;
  }
  EXPECT_LT(largest, 0.8);
  EXPECT_GT(weights[3] / normalization, 0.05);
}

// The pair of atoms joined by a covalent cross-link sits at bond distance, where its unexcluded Lennard-Jones energy
// is far above the overlap criterion. The inter-molecular difference of a displacement move must not reject the move
// as a hard-core overlap on that pair: the cross-link term removes the pair again, and the move goes to the
// Metropolis test with the bond energy.
TEST(MC_CROSS_LINKS, displacing_a_linked_monomer_is_not_rejected_as_an_overlap)
{
  const ForceField forceField = makeForceField(100.0, 0.0, 0.0, false);
  const double3 posA(10.0, 10.0, 10.0);
  const double3 posB(10.0 + kBondLength, 10.0, 10.0);
  Component monomerA = makeComponent(forceField, 0, "xl-monomer-a-move", monomerJson("A", "X"));
  Component monomerB = makeComponent(forceField, 1, "xl-monomer-b-move", monomerJson("B", "Y"));
  System system = makeSystem(forceField, {monomerA, monomerB}, {{posA}, {posB}});
  system.crossLinks.bondTypes.push_back(harmonicBondType("X", "Y", kBondK, kBondLength, 3.0));
  system.crossLinks.addLink(CrossLink{site(0, 0, 0), site(1, 0, 0), 0});
  system.runningEnergies = system.computeTotalEnergies();

  // the bare pair energy at the bond distance is far above the overlap criterion
  const double sigma = 0.5 * (3.4 + 3.6);
  const double bareLennardJones = 4.0 * 100.0 * (std::pow(sigma / kBondLength, 12) - std::pow(sigma / kBondLength, 6));
  ASSERT_GT(bareLennardJones, forceField.energyOverlapCriteria);

  // tiny displacements: the bond energy changes by less than a kelvin, so nearly every attempt is accepted
  system.components[0].mc_moves_statistics.setMaxChange(Move::Types::Translation, double3(0.001, 0.001, 0.001));
  RandomNumber random(11);
  std::size_t accepted = 0;
  for (std::size_t attempt = 0; attempt < 200; ++attempt)
  {
    const std::optional<RunningEnergy> difference = MC_Moves::translationMove(random, system, 0, 0);
    if (difference.has_value())
    {
      ++accepted;
      system.runningEnergies += difference.value();
    }
  }
  EXPECT_GT(accepted, 150uz);
  // the running energy follows the accepted moves exactly
  const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
  EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6 * std::max(1.0, std::abs(system.runningEnergies.potentialEnergy())));
}
