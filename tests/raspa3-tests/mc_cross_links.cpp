#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import atom_dynamics;
import units;
import forcefield;
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
import mc_moves_partial_reinsertion;
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
TEST(MC_CROSS_LINKS, tethered_regrowth_samples_junction_angle_distribution)
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
    ForceField forceField = makeForceField(0.0, 0.0, 0.0, false);
    forceField.useRecoilGrowth = useRecoilGrowth;

    MCMoveProbabilities componentMoves;
    componentMoves.setProbability(Move::Types::ReinsertionCBMC, 1.0);
    TemporaryFile file("xl-dimer-regrow.json", dimerBAJson(kBondK, kBondLength));
    Component dimer = Component(Component::Type::Adsorbate, 0, forceField, "xl-dimer-regrow", file.stemPath().string(),
                                5, 21, componentMoves, std::nullopt, false);

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
    for (std::size_t i = 0; i < 105000; ++i)
    {
      const std::size_t molecule = random.uniform_integer(0, 1);
      std::optional<RunningEnergy> energy = MC_Moves::reinsertionMove(random, system, 0, molecule);
      if (energy.has_value()) system.runningEnergies += energy.value();
      if (i < 5000) continue;  // burn-in

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
    EXPECT_GT(reinsertion.totalAccepted, 10000.0) << "recoil growth: " << useRecoilGrowth;
    const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.crossLink, 0.0, 1e-6);

    EXPECT_NEAR(sumCos / samples, expectedCos, 0.01) << "recoil growth: " << useRecoilGrowth;
    EXPECT_NEAR(sumBend / samples, expectedBend, 0.03 * expectedBend) << "recoil growth: " << useRecoilGrowth;
  }
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
