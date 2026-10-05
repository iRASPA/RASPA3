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
import mc_moves_pivot_cbmc;
import mc_moves_translation;

// Tests for the configurational-bias pivot move. On a linear chain with fixed bonds and bends every
// pivot bond is a torsion, so the move alone is ergodic in the torsions and their exact Boltzmann
// distribution is a sharp test of the Rosenbluth acceptance rule (W_new / W_old): the torsion energy
// enters the move only through the trial weights. The bookkeeping with intermolecular Lennard-Jones,
// real-space and Fourier-space Coulomb terms is checked against a full recomputation in a dense
// system of chains, where the biased selection must also raise the acceptance of fully randomized
// pivots above that of the plain pivot while sampling the same distribution.

namespace
{

constexpr double kBondLength = 1.54;

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

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}}, {{0.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
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

// Largest deviation of the backbone bond lengths and bend cosines of a linear chain from the fixed
// values.
std::pair<double, double> backboneDeviation(std::span<const Atom> atoms)
{
  double bondDeviation = 0.0, bendDeviation = 0.0;
  const double cosFixed = std::cos(114.0 * std::numbers::pi / 180.0);
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    bondDeviation = std::max(bondDeviation, std::abs((atoms[i + 1].position - atoms[i].position).length() - kBondLength));
  }
  for (std::size_t i = 0; i + 2 < atoms.size(); ++i)
  {
    const double cosine = double3::dot((atoms[i].position - atoms[i + 1].position).normalized(),
                                       (atoms[i + 2].position - atoms[i + 1].position).normalized());
    bendDeviation = std::max(bendDeviation, std::abs(cosine - cosFixed));
  }
  return {bondDeviation, bendDeviation};
}

// Mean cosine of all backbone torsions of all chains of component 0.
double meanTorsionCosine(const System& system)
{
  double sum = 0.0;
  std::size_t count = 0;
  for (std::size_t m = 0; m != system.numberOfMoleculesPerComponent[0]; ++m)
  {
    std::span<const Atom> atoms = system.spanOfMolecule(0, m);
    for (std::size_t i = 0; i + 3 < atoms.size(); ++i)
    {
      sum += std::cos(
          dihedralAngle(atoms[i].position, atoms[i + 1].position, atoms[i + 2].position, atoms[i + 3].position));
      ++count;
    }
  }
  return sum / static_cast<double>(count);
}

}  // namespace

// With U = 0 every trial has unit weight, W_new = W_old = k, and every constructed proposal must be
// accepted; the torsions are then uniform and the fixed bonds and bends are preserved exactly.
TEST(MC_PIVOT_CBMC, u0_accepts_every_proposal_and_preserves_constraints)
{
  const ForceField forceField = makeZeroForceField();
  Component decane = makeComponent(forceField, "cbpivot-decane-u0", chainJson(10, "0.0, 0.0, 0.0, 0.0"));
  System system =
      System(forceField, SimulationBox(200.0, 200.0, 200.0), false, 300.0, 1e4, 1.0, {}, {decane}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotCBMCNumberOfTrialAngles = 6;
  system.components[0].pivotCBMCRandomizationFraction = 0.5;

  RandomNumber random(12345);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);

  constexpr std::size_t numberOfMoves = 20000;
  std::size_t accepted = 0;
  double sumCos = 0.0;
  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    std::optional<RunningEnergy> energyDifference = MC_Moves::pivotCBMCMove(random, system, 0, 0);
    if (energyDifference.has_value())
    {
      ++accepted;
      system.runningEnergies += energyDifference.value();
      EXPECT_NEAR(energyDifference->potentialEnergy(), 0.0, 1e-10);
    }
    if (i % 10 == 0)
    {
      sumCos += std::cos(dihedralAngle(atoms[3].position, atoms[4].position, atoms[5].position, atoms[6].position));
    }
  }
  EXPECT_EQ(accepted, numberOfMoves);
  EXPECT_NEAR(sumCos / static_cast<double>(numberOfMoves / 10), 0.0, 0.05);

  const auto [bondDeviation, bendDeviation] = backboneDeviation(atoms);
  EXPECT_LT(bondDeviation, 1e-8);
  EXPECT_LT(bendDeviation, 1e-8);
}

// Boltzmann sampler test with the TraPPE alkane torsion potential at 600 K on a 10-bead chain in
// vacuum: each torsion angle is independently distributed as exp(-U(phi)/T), with
// <cos phi> = -0.2633 and a trans fraction (|phi| > 120 degrees) of 0.4946. The torsion energy
// enters the move only through the Rosenbluth weights, so this checks the W_new / W_old rule
// (including the reverse trial set drawn about the new configuration) on both channels. Without
// the W_old factor the sampled distribution would be exp(-2 beta U) (too few gauche states).
TEST(MC_PIVOT_CBMC, trappe_torsions_follow_boltzmann_distribution)
{
  constexpr double temperature = 600.0;
  constexpr double exactCos = -0.263285;
  constexpr double exactTrans = 0.494568;

  const ForceField forceField = makeZeroForceField();
  Component decane =
      makeComponent(forceField, "cbpivot-decane-trappe", chainJson(10, "0.0, 355.03, -68.19, 791.32"));
  System system = System(forceField, SimulationBox(200.0, 200.0, 200.0), false, temperature, 1e4, 1.0, {}, {decane},
                         {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotCBMCNumberOfTrialAngles = 8;
  system.components[0].pivotCBMCRandomizationFraction = 0.5;

  RandomNumber random(31337);
  std::span<Atom> atoms = system.spanOfMolecule(0, 0);

  constexpr std::size_t numberOfMoves = 200000;
  constexpr std::size_t burnIn = 20000;
  constexpr std::size_t sampleEvery = 5;
  std::array<double, 7> sumCos{};
  std::array<std::size_t, 7> transCounts{};
  std::size_t samples = 0;

  for (std::size_t i = 0; i != numberOfMoves; ++i)
  {
    // Adapt the small-step channel during burn-in only (the step sizes must be fixed while
    // sampling to keep the chain reversible).
    if (i < burnIn && i % 2000 == 1999) system.components[0].mc_moves_statistics.optimizeMCMoves();

    std::optional<RunningEnergy> energyDifference = MC_Moves::pivotCBMCMove(random, system, 0, 0);
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();

    if (i < burnIn || i % sampleEvery != 0) continue;
    ++samples;
    for (std::size_t torsion = 0; torsion != 7; ++torsion)
    {
      double phi = dihedralAngle(atoms[torsion].position, atoms[torsion + 1].position, atoms[torsion + 2].position,
                                 atoms[torsion + 3].position);
      sumCos[torsion] += std::cos(phi);
      if (std::abs(phi) > 2.0 * std::numbers::pi / 3.0) ++transCounts[torsion];
    }
  }

  double meanCos = 0.0, meanTrans = 0.0;
  for (std::size_t torsion = 0; torsion != 7; ++torsion)
  {
    const double cosine = sumCos[torsion] / static_cast<double>(samples);
    const double trans = static_cast<double>(transCounts[torsion]) / static_cast<double>(samples);
    meanCos += cosine / 7.0;
    meanTrans += trans / 7.0;
    EXPECT_NEAR(cosine, exactCos, 0.04) << "torsion " << torsion;
    EXPECT_NEAR(trans, exactTrans, 0.04) << "torsion " << torsion;
  }
  EXPECT_NEAR(meanCos, exactCos, 0.015);
  EXPECT_NEAR(meanTrans, exactTrans, 0.015);

  const auto [bondDeviation, bendDeviation] = backboneDeviation(atoms);
  EXPECT_LT(bondDeviation, 1e-8);
  EXPECT_LT(bendDeviation, 1e-8);

  // Energy bookkeeping: the accumulated running energy must match a full recomputation.
  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy() * Units::EnergyToKelvin,
              recomputed.potentialEnergy() * Units::EnergyToKelvin, 1e-6);
}

// Twelve interacting chains with Lennard-Jones beads carrying alternating charges (Ewald summation)
// and a TraPPE torsion potential at 600 K. (i) The accumulated running energy must match a full
// recomputation term by term (the bias contains the real-space terms, the Fourier part enters as a
// correction). (ii) Fully randomized CBMC pivots must be accepted more often than fully randomized
// plain pivots in this dense state. (iii) Both samplers must agree on the mean torsion cosine.
TEST(MC_PIVOT_CBMC, dense_chains_bookkeeping_acceptance_and_agreement_with_plain_pivot)
{
  const ForceField forceField = makeChargedForceField();
  Component chain =
      makeComponent(forceField, "cbpivot-c14-charged", chainJson(14, "0.0, 355.03, -68.19, 791.32", true));

  constexpr std::size_t numberOfChains = 12;
  System system =
      System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 600.0, 1e4, 1.0, {}, {chain}, {}, {numberOfChains}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  system.components[0].pivotRandomizationFraction = 1.0;
  system.components[0].pivotCBMCRandomizationFraction = 1.0;
  system.components[0].pivotCBMCNumberOfTrialAngles = 8;

  RandomNumber random(8080);
  auto randomChain = [&] { return static_cast<std::size_t>(random.uniform() * static_cast<double>(numberOfChains)); };

  // Sampling with one pivot flavour (alternated with translations); returns the mean torsion cosine
  // and the acceptance ratio of the pivot flavour.
  auto run = [&](bool biased, std::size_t numberOfMoves, std::size_t burnIn) -> std::pair<double, double>
  {
    std::size_t attempts = 0, accepted = 0, samples = 0;
    double sumCos = 0.0;
    for (std::size_t i = 0; i != numberOfMoves; ++i)
    {
      std::optional<RunningEnergy> energyDifference;
      if (i % 2 == 0)
      {
        energyDifference = biased ? MC_Moves::pivotCBMCMove(random, system, 0, randomChain())
                                  : MC_Moves::pivotMove(random, system, 0, randomChain());
        if (i >= burnIn)
        {
          ++attempts;
          if (energyDifference.has_value()) ++accepted;
        }
      }
      else
      {
        energyDifference = MC_Moves::translationMove(random, system, 0, randomChain());
      }
      if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();
      if (i >= burnIn && i % 20 == 0)
      {
        ++samples;
        sumCos += meanTorsionCosine(system);
      }
    }
    return {sumCos / static_cast<double>(samples), static_cast<double>(accepted) / static_cast<double>(attempts)};
  };

  const auto [plainCos, plainAcceptance] = run(false, 60000, 10000);
  const auto [biasedCos, biasedAcceptance] = run(true, 60000, 10000);

  // Measured with this seed: plain 0.34, biased 0.71; mean cosines agree to 0.005.
  EXPECT_GT(biasedAcceptance, 1.5 * plainAcceptance);
  EXPECT_NEAR(biasedCos, plainCos, 0.03);

  system.checkMoleculeIds();
  const RunningEnergy running = system.runningEnergies;
  const RunningEnergy recomputed = system.computeTotalEnergies();
  const double toKelvin = Units::EnergyToKelvin;
  EXPECT_NEAR(running.moleculeMoleculeVDW * toKelvin, recomputed.moleculeMoleculeVDW * toKelvin, 1e-4);
  EXPECT_NEAR(running.moleculeMoleculeCharge * toKelvin, recomputed.moleculeMoleculeCharge * toKelvin, 1e-4);
  EXPECT_NEAR(running.ewald_fourier * toKelvin, recomputed.ewald_fourier * toKelvin, 1e-4);
  EXPECT_NEAR(running.ewald_exclusion * toKelvin, recomputed.ewald_exclusion * toKelvin, 1e-4);
  EXPECT_NEAR(running.intraVDW * toKelvin, recomputed.intraVDW * toKelvin, 1e-4);
  EXPECT_NEAR(running.intraCoul * toKelvin, recomputed.intraCoul * toKelvin, 1e-4);
  EXPECT_NEAR(running.torsion * toKelvin, recomputed.torsion * toKelvin, 1e-4);
  EXPECT_NEAR(running.potentialEnergy() * toKelvin, recomputed.potentialEnergy() * toKelvin, 1e-4);
}
