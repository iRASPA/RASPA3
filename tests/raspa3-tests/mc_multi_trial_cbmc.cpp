#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import int3;
import double3;
import double3x3;
import randomnumbers;
import units;
import atom;
import molecule;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_statistics;
import move_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import simd_quatd;
import spatial_decomposition_settings;
import cbmc;
import cbmc_external_energy;
import cbmc_multi_trial;
import cbmc_multi_trial_cpu;
import cbmc_device_environment;
import cbmc_multi_trial_device;

// Multiple-try CBMC moves (swap, reinsertion, Gibbs swap): the M trial molecules are grown with the long-range
// terms folded into every weight before one is selected. These tests check that
//   * the CPU grower with M = 1 is the plain CBMC grow (same selected molecule, same weight),
//   * the running energies do not drift for rigid (Ewald) and flexible (tail-corrected) molecules,
//   * the moves are exercised and accepted,
//   * the sampled loading agrees with the single-trial swap move.

namespace
{

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

Component makeAlkaneFromExample(const ForceField &forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities &probabilities)
{
  TemporaryFile file(std::string(name) + ".json", molecule_fixtures::alkaneJson(name));
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), file.stemPath().string(), 5,
                   21, probabilities, std::nullopt, false);
}

void expectNoEnergyDrift(System &system, double tolerance = 1e-6)
{
  RunningEnergy recomputedEnergies = system.computeTotalEnergies();
  RunningEnergy drift = system.runningEnergies - recomputedEnergies;

  EXPECT_NEAR(drift.potentialEnergy(), 0.0, tolerance);
  EXPECT_NEAR(drift.externalFieldVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.externalFieldCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_fourier, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_self, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_exclusion, 0.0, tolerance);
  EXPECT_NEAR(drift.intraVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.intraCoul, 0.0, tolerance);
  EXPECT_NEAR(drift.tail, 0.0, tolerance);
  EXPECT_NEAR(drift.polarization, 0.0, tolerance);
}

MonteCarlo makeMonteCarlo(std::vector<System> systems, std::size_t production, std::size_t initialization,
                          std::size_t equilibration)
{
  return MonteCarlo({production, 0, initialization, equilibration, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
}

/// The first device kind the multiple-try grower can run on here (Metal on macOS, else CUDA or OpenCL).
std::optional<PairDevice> availableCBMCDevice()
{
  for (PairDevice device : {PairDevice::Metal, PairDevice::CUDA, PairDevice::OpenCL})
  {
    if (CBMC::DeviceGrowEnvironment::available(device)) return device;
  }
  return std::nullopt;
}

Component makeCO2(const ForceField &forceField, const MCMoveProbabilities &probabilities, std::uint8_t componentId = 0)
{
  return Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                   {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, componentId, false, false),
                    Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, componentId, false, false),
                    Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, componentId, false, false)},
                   {}, {}, 5, 21, probabilities, std::nullopt, false);
}

Component makeWater(const ForceField &forceField, const MCMoveProbabilities &probabilities, std::uint8_t componentId)
{
  return Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, componentId, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, componentId, false,
            false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, componentId, false,
            false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, componentId, false,
            false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, componentId, false,
            false)},
      {}, {}, 5, 21, probabilities, std::nullopt, false);
}

}  // namespace

// With M = 1 the multi-trial grower is the plain CBMC grow: the same random stream gives the same molecule and
// log W-bar = log W (the long-range callback returns zero here).
TEST(MC_MULTI_TRIAL_CBMC, cpu_grower_single_trial_is_plain_cbmc)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, MCMoveProbabilities(), std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2}, {}, {8}, 5);

  const CBMC::GrowContext context = system.makeGrowContext();
  const CBMC::MultiTrialGrower &grower = CBMC::cpuMultiTrialGrower();
  const CBMC::NewMoleculeIdentity identity{.componentId = 0, .moleculeId = system.numberOfMolecules()};
  const CBMC::LongRangeCorrection zero = [](std::span<const Atom>) { return 0.0; };

  RandomNumber randomA(7);
  RandomNumber randomB(7);

  std::optional<CBMC::GrowResult> plain = CBMC::growNewMolecule(randomA, context, system.components[0], identity);
  std::optional<CBMC::MultiTrialGrowResult> multi =
      grower.growNew(randomB, context, system.components[0], identity, 1, zero);

  ASSERT_TRUE(plain.has_value());
  ASSERT_TRUE(multi.has_value());
  EXPECT_EQ(multi->numberOfConstructed, 1uz);
  EXPECT_NEAR(multi->logMeanWeight, plain->logRosenbluthWeight, 1e-12);
  EXPECT_NEAR(multi->selected.logRosenbluthWeight, plain->logRosenbluthWeight, 1e-12);
  EXPECT_NEAR(multi->selected.energies.potentialEnergy(), plain->energies.potentialEnergy(), 1e-10);
  ASSERT_EQ(multi->selected.atoms.size(), plain->atoms.size());
  for (std::size_t i = 0; i < plain->atoms.size(); ++i)
  {
    EXPECT_NEAR(multi->selected.atoms[i].position.x, plain->atoms[i].position.x, 1e-12);
    EXPECT_NEAR(multi->selected.atoms[i].position.y, plain->atoms[i].position.y, 1e-12);
    EXPECT_NEAR(multi->selected.atoms[i].position.z, plain->atoms[i].position.z, 1e-12);
  }
}

// log W-bar = log((1/M) sum_i W_i e^{-beta dU_i}): the mean lies between the smallest and the largest trial
// weight, and a constant long-range shift moves it by exactly -beta * shift.
TEST(MC_MULTI_TRIAL_CBMC, cpu_grower_mean_weight_shift_invariance)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 0, false, false)}, {}, {}, 5, 21,
                                MCMoveProbabilities(), std::nullopt, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {methane}, {}, {8}, 5);

  const CBMC::GrowContext context = system.makeGrowContext();
  const CBMC::MultiTrialGrower &grower = CBMC::cpuMultiTrialGrower();
  const CBMC::NewMoleculeIdentity identity{.componentId = 0, .moleculeId = system.numberOfMolecules()};

  const double shift = 250.0 * Units::KelvinToEnergy;
  const CBMC::LongRangeCorrection zero = [](std::span<const Atom>) { return 0.0; };
  const CBMC::LongRangeCorrection shifted = [shift](std::span<const Atom>) { return shift; };

  RandomNumber randomA(11);
  RandomNumber randomB(11);
  std::optional<CBMC::MultiTrialGrowResult> a = grower.growNew(randomA, context, system.components[0], identity, 16, zero);
  std::optional<CBMC::MultiTrialGrowResult> b =
      grower.growNew(randomB, context, system.components[0], identity, 16, shifted);

  ASSERT_TRUE(a.has_value());
  ASSERT_TRUE(b.has_value());
  EXPECT_EQ(a->numberOfConstructed, 16uz);
  EXPECT_NEAR(b->logMeanWeight, a->logMeanWeight - system.beta * shift, 1e-9);
  EXPECT_NEAR(b->selectedLongRangeEnergy, shift, 1e-12);
  EXPECT_NEAR(a->selectedLongRangeEnergy, 0.0, 1e-12);

  // the selection is unchanged by a constant shift (same random stream, same relative weights)
  EXPECT_NEAR(a->selected.logRosenbluthWeight, b->selected.logRosenbluthWeight, 1e-12);
}

TEST(MC_MULTI_TRIAL_CBMC, log_sum_exp_and_selection)
{
  const std::array<double, 4> logWeights{std::log(1.0), std::log(2.0), std::log(3.0), std::log(4.0)};
  EXPECT_NEAR(CBMC::logSumExp(logWeights), std::log(10.0), 1e-12);

  // selection probabilities proportional to the weights
  RandomNumber random(3);
  std::array<std::size_t, 4> counts{};
  const std::size_t samples = 200000;
  for (std::size_t i = 0; i < samples; ++i) ++counts[CBMC::selectByLogWeight(random, logWeights)];
  for (std::size_t i = 0; i < 4; ++i)
  {
    const double expected = static_cast<double>(i + 1) / 10.0;
    EXPECT_NEAR(static_cast<double>(counts[i]) / static_cast<double>(samples), expected, 0.01);
  }
}

TEST(MC_MULTI_TRIAL_CBMC, swap_and_reinsertion_drift_rigid_ewald_fau)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionMultiTrialCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapMultiTrialCBMC, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities, std::nullopt, false);
  co2.numberOfTrialMolecules = 8;

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21, probabilities,
                                std::nullopt, false);
  methane.numberOfTrialMolecules = 5;

  MCMoveProbabilities probabilities_water = MCMoveProbabilities();
  probabilities_water.setProbability(Move::Types::Translation, 1.0);
  probabilities_water.setProbability(Move::Types::Rotation, 1.0);
  probabilities_water.setProbability(Move::Types::ReinsertionMultiTrialCBMC, 1.0);

  Component water = Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, 2, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), 0.241, 1.0, 0, 8, 2, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), -0.241, 1.0, 0, 9, 2, false, false)},
      {}, {}, 5, 21, probabilities_water, std::nullopt, false);
  water.numberOfTrialMolecules = 4;

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);

  MonteCarlo mc = makeMonteCarlo({system}, 20, 5, 5);
  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);

    const MoveStatistics<double3> &swapStatistics =
        std::get<MoveStatistics<double3>>(s.components[0].mc_moves_statistics[Move::Types::SwapMultiTrialCBMC]);
    EXPECT_GT(swapStatistics.totalCounts.x, 0.0);
    EXPECT_GT(swapStatistics.totalCounts.y, 0.0);
    EXPECT_GT(swapStatistics.totalAccepted.x + swapStatistics.totalAccepted.y, 0.0);

    const MoveStatistics<double> &reinsertionStatistics = std::get<MoveStatistics<double>>(
        s.components[2].mc_moves_statistics[Move::Types::ReinsertionMultiTrialCBMC]);
    EXPECT_GT(reinsertionStatistics.totalCounts, 0.0);
    EXPECT_GT(reinsertionStatistics.totalAccepted, 0.0);
  }
}

TEST(MC_MULTI_TRIAL_CBMC, swap_and_reinsertion_drift_flexible_tail_corrections)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionMultiTrialCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapMultiTrialCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
  propane.numberOfTrialMolecules = 6;
  butane.numberOfTrialMolecules = 3;
  propane.idealGasRosenbluthWeight = 1.0;
  butane.idealGasRosenbluthWeight = 1.0;

  System system = System(forceField, box, false, 400.0, 1e6, 1.0, {}, {propane, butane}, {}, {30, 30}, 5);

  MonteCarlo mc = makeMonteCarlo({system}, 20, 5, 5);
  mc.run();

  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
    for (std::size_t i = 0; i < 2; ++i)
    {
      const MoveStatistics<double3> &swapStatistics =
          std::get<MoveStatistics<double3>>(s.components[i].mc_moves_statistics[Move::Types::SwapMultiTrialCBMC]);
      EXPECT_GT(swapStatistics.totalAccepted.x + swapStatistics.totalAccepted.y, 0.0);
      const MoveStatistics<double> &reinsertionStatistics = std::get<MoveStatistics<double>>(
          s.components[i].mc_moves_statistics[Move::Types::ReinsertionMultiTrialCBMC]);
      EXPECT_GT(reinsertionStatistics.totalAccepted, 0.0);
    }
  }
}

TEST(MC_MULTI_TRIAL_CBMC, gibbs_swap_drift_propane_butane_tail_corrections)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapMultiTrialCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", probabilities);
  propane.numberOfTrialMolecules = 4;
  butane.numberOfTrialMolecules = 4;

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::GibbsVolume, 0.01);

  System systemVapor =
      System(forceField, box, false, 400.0, 1e4, 1.0, {}, {propane, butane}, {}, {40, 40}, 5, systemProbabilities);
  System systemLiquid =
      System(forceField, box, false, 400.0, 1e4, 1.0, {}, {propane, butane}, {}, {40, 40}, 5, systemProbabilities);

  MonteCarlo mc = makeMonteCarlo({systemVapor, systemLiquid}, 50, 10, 10);
  mc.run();

  std::size_t total{0};
  double accepted{0.0};
  for (System &s : mc.systems)
  {
    expectNoEnergyDrift(s);
    total += s.numberOfIntegerMoleculesPerComponent[0] + s.numberOfIntegerMoleculesPerComponent[1];
    for (std::size_t i = 0; i < 2; ++i)
    {
      accepted += std::get<MoveStatistics<double>>(
                      s.components[i].mc_moves_statistics[Move::Types::GibbsSwapMultiTrialCBMC])
                      .totalAccepted;
    }
  }
  EXPECT_EQ(total, 160uz);
  EXPECT_GT(accepted, 0.0);
}

// The multiple-try swap samples the same grand-canonical distribution as the single-trial swap: the average
// loading of methane in FAU must agree within the statistical error.
TEST(MC_MULTI_TRIAL_CBMC, loading_matches_single_trial_swap_methane_fau)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  auto run = [&](Move::Types swapMove, std::size_t numberOfTrialMolecules) -> std::pair<double, double>
  {
    MCMoveProbabilities probabilities = MCMoveProbabilities();
    probabilities.setProbability(Move::Types::Translation, 1.0);
    probabilities.setProbability(swapMove, 2.0);

    Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                  {Atom({0, 0, 0}, 0.0, 1.0, 0, 6, 0, false, false)}, {}, {}, 5, 21, probabilities,
                                  std::nullopt, false);
    methane.numberOfTrialMolecules = numberOfTrialMolecules;

    System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {methane}, {}, {10}, 5);
    MonteCarlo mc = makeMonteCarlo({system}, 400, 100, 0);
    mc.run();
    return mc.systems[0].averageLoadings.averageLoadingNumberOfMolecules(0);
  };

  const auto [meanSingle, errorSingle] = run(Move::Types::SwapCBMC, 1);
  const auto [meanMulti, errorMulti] = run(Move::Types::SwapMultiTrialCBMC, 16);

  EXPECT_GT(meanSingle, 1.0);
  EXPECT_GT(meanMulti, 1.0);
  EXPECT_NEAR(meanMulti, meanSingle, 3.0 * (errorSingle + errorMulti) + 0.05 * meanSingle)
      << "single-trial: " << meanSingle << " +/- " << errorSingle << ", multi-trial: " << meanMulti << " +/- "
      << errorMulti;
}

// ---------------------------------------------------------------------------------------------------------------
// Device grower (skipped when no OpenCL/Metal/CUDA device is available)

// The device kernel's single-precision external energies (Lennard-Jones + real-space Ewald, framework and
// molecules, same-molecule exclusion) agree with the host's double-precision CBMC external energy for random
// rigid CO2 and water trial molecules, and for the molecules already present (the retrace evaluation).
TEST(MC_MULTI_TRIAL_CBMC, device_energies_match_host_double)
{
  const std::optional<PairDevice> device = availableCBMCDevice();
  if (!device.has_value()) GTEST_SKIP() << "no CBMC device available";

  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component co2 = makeCO2(forceField, MCMoveProbabilities(), 0);
  Component water = makeWater(forceField, MCMoveProbabilities(), 1);
  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, water}, {}, {12, 8}, 5);

  const CBMC::GrowContext context = system.makeGrowContext(CBMC::CutOffMode::Full);
  CBMC::DeviceGrowEnvironment environment(device.value());
  environment.setForceField(context.forceField, context.simulationBox);
  environment.setBackground(context.frameworkAtoms, context.moleculeAtoms);

  RandomNumber random(5);
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double tolerance = 2e-4;  // relative, single precision over ~10^4 pair terms

  auto compare = [&](const Component &component, std::span<const Atom> batch, std::uint32_t trialMoleculeId,
                     std::size_t &compared, std::size_t &flagMismatches)
  {
    const std::size_t n = component.atoms.size();
    std::vector<CBMC::DeviceTrialEnergy> energies;
    environment.evaluate(batch, n, trialMoleculeId, std::nullopt, energies);
    ASSERT_EQ(energies.size(), batch.size() / n);
    for (std::size_t set = 0; set < energies.size(); ++set)
    {
      const std::span<const Atom> atoms = batch.subspan(set * n, n);
      const std::optional<RunningEnergy> host = computeExternalNonOverlappingEnergy(context, component, atoms);
      if (!host.has_value() || energies[set].overlap)
      {
        // the flags may only differ when a pair energy lies at the overlap threshold
        if (host.has_value() != !energies[set].overlap)
        {
          if (host.has_value() && host->potentialEnergy() < 0.5 * overlapCriteria) ++flagMismatches;
          if (!host.has_value() && energies[set].total() < 0.5 * overlapCriteria) ++flagMismatches;
        }
        continue;
      }
      const double hostVDW = host->frameworkMoleculeVDW + host->moleculeMoleculeVDW;
      const double hostCharge = host->frameworkMoleculeCharge + host->moleculeMoleculeCharge;
      const double scale = std::max({std::abs(hostVDW), std::abs(hostCharge), 1000.0 * Units::KelvinToEnergy});
      EXPECT_NEAR(energies[set].vdw, hostVDW, tolerance * scale) << "set " << set;
      EXPECT_NEAR(energies[set].charge, hostCharge, tolerance * scale) << "set " << set;
      ++compared;
    }
  };

  // random rigid trial molecules of a new molecule
  for (const Component &component : {std::cref(co2), std::cref(water)})
  {
    const std::size_t n = component.atoms.size();
    const std::size_t sets = 512;
    std::vector<Atom> batch(sets * n);
    for (std::size_t set = 0; set < sets; ++set)
    {
      const double3 center = context.simulationBox.randomPosition(random);
      const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
      for (std::size_t i = 0; i < n; ++i)
      {
        Atom atom = component.atoms[i];
        atom.position = center + rotation * component.atoms[i].position;
        atom.moleculeId = static_cast<std::uint32_t>(system.numberOfMolecules());
        batch[set * n + i] = atom;
      }
    }
    std::size_t compared = 0, flagMismatches = 0;
    compare(component, batch, static_cast<std::uint32_t>(system.numberOfMolecules()), compared, flagMismatches);
    EXPECT_GT(compared, 50uz) << component.name;
    EXPECT_EQ(flagMismatches, 0uz) << component.name;
  }

  // the molecules present, evaluated with their own id left out of the background (the retrace side)
  for (std::size_t componentId = 0; componentId < 2; ++componentId)
  {
    const Component &component = system.components[componentId];
    for (std::size_t m = 0; m < system.numberOfMoleculesPerComponent[componentId]; ++m)
    {
      const std::span<const Atom> atoms = system.spanOfMolecule(componentId, m);
      std::size_t compared = 0, flagMismatches = 0;
      compare(component, atoms, atoms.front().moleculeId, compared, flagMismatches);
      EXPECT_EQ(compared, 1uz) << component.name << " molecule " << m;
      EXPECT_EQ(flagMismatches, 0uz);
    }
  }
}

// The device grower reproduces the CPU grower's bookkeeping: with the system's CBMC device set, the multiple-try
// swap and reinsertion of rigid Ewald molecules in FAU leave no energy drift and are accepted.
TEST(MC_MULTI_TRIAL_CBMC, device_swap_and_reinsertion_drift_rigid_ewald_fau)
{
  const std::optional<PairDevice> device = availableCBMCDevice();
  if (!device.has_value()) GTEST_SKIP() << "no CBMC device available";

  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionMultiTrialCBMC, 1.0);
  probabilities.setProbability(Move::Types::SwapMultiTrialCBMC, 1.0);

  Component co2 = makeCO2(forceField, probabilities, 0);
  co2.numberOfTrialMolecules = 8;
  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)}, {}, {}, 5, 21, probabilities,
                                std::nullopt, false);
  methane.numberOfTrialMolecules = 5;
  Component water = makeWater(forceField, probabilities, 2);
  water.numberOfTrialMolecules = 4;

  System system =
      System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {co2, methane, water}, {}, {10, 15, 8}, 5);
  system.cbmcDevice = device.value();

  // all three components are grown on the device
  const CBMC::GrowContext context = system.makeGrowContext();
  for (const Component &component : system.components)
  {
    std::string reason;
    EXPECT_TRUE(CBMC::DeviceMultiTrialGrower::supports(context, component, reason)) << component.name << ": " << reason;
  }

  MonteCarlo mc = makeMonteCarlo({system}, 20, 5, 5);
  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.cbmcDevice, device.value());
    expectNoEnergyDrift(s);
    for (std::size_t componentId = 0; componentId < 3; ++componentId)
    {
      const MoveStatistics<double3> &swapStatistics = std::get<MoveStatistics<double3>>(
          s.components[componentId].mc_moves_statistics[Move::Types::SwapMultiTrialCBMC]);
      EXPECT_GT(swapStatistics.totalCounts.x, 0.0);
      EXPECT_GT(swapStatistics.totalCounts.y, 0.0);
      EXPECT_GT(swapStatistics.totalAccepted.x + swapStatistics.totalAccepted.y, 0.0) << componentId;

      const MoveStatistics<double> &reinsertionStatistics = std::get<MoveStatistics<double>>(
          s.components[componentId].mc_moves_statistics[Move::Types::ReinsertionMultiTrialCBMC]);
      EXPECT_GT(reinsertionStatistics.totalCounts, 0.0);
      EXPECT_GT(reinsertionStatistics.totalAccepted, 0.0) << componentId;
    }
  }
}

// The device kernel evaluates unscaled pairs: a fractional (CFCMC) molecule among the background atoms, not only
// in the component being grown, sends the grow to the CPU grower. The 'supports' check is static: no device needed.
TEST(MC_MULTI_TRIAL_CBMC, device_grower_declines_a_background_with_a_fractional_molecule)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities multiTrial = MCMoveProbabilities();
  multiTrial.setProbability(Move::Types::Translation, 1.0);
  multiTrial.setProbability(Move::Types::SwapMultiTrialCBMC, 1.0);
  MCMoveProbabilities cfcmc = MCMoveProbabilities();
  cfcmc.setProbability(Move::Types::Translation, 1.0);
  cfcmc.setProbability(Move::Types::SwapCFCMC, 1.0);

  auto makeMethane = [&](const MCMoveProbabilities &probabilities)
  {
    return Component(forceField, "methane", 190.564, 45599200, 0.01142,
                     {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                     false);
  };

  // CO2 with a CFCMC fractional molecule: the methane grows see it in their background
  {
    System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f},
                           {makeMethane(multiTrial), makeCO2(forceField, cfcmc, 1)}, {}, {10, 5}, 5);
    ASSERT_TRUE(system.components[1].hasFractionalMolecule);
    const CBMC::GrowContext context = system.makeGrowContext();
    const bool backgroundIsFractional = std::ranges::any_of(
        context.moleculeAtoms, [](const Atom &atom) { return atom.isFractional || atom.scalingVDW != 1.0; });
    ASSERT_TRUE(backgroundIsFractional);

    std::string reason;
    EXPECT_FALSE(CBMC::DeviceMultiTrialGrower::supports(context, system.components[0], reason));
    EXPECT_NE(reason.find("fractional"), std::string::npos) << reason;
  }

  // the same mixture without a fractional molecule is grown on the device
  {
    System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f},
                           {makeMethane(multiTrial), makeCO2(forceField, multiTrial, 1)}, {}, {10, 5}, 5);
    ASSERT_FALSE(system.components[1].hasFractionalMolecule);
    const CBMC::GrowContext context = system.makeGrowContext();
    std::string reason;
    EXPECT_TRUE(CBMC::DeviceMultiTrialGrower::supports(context, system.components[0], reason)) << reason;
  }
}

// With the CBMC device set and a fractional molecule present, the multiple-try moves run through the CPU grower:
// the simulation proceeds with the moves accepted and the running energies without drift.
TEST(MC_MULTI_TRIAL_CBMC, device_grower_falls_back_to_the_cpu_grower_next_to_a_fractional_molecule)
{
  const std::optional<PairDevice> device = availableCBMCDevice();
  if (!device.has_value()) GTEST_SKIP() << "no CBMC device available";

  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities multiTrial = MCMoveProbabilities();
  multiTrial.setProbability(Move::Types::Translation, 1.0);
  multiTrial.setProbability(Move::Types::ReinsertionMultiTrialCBMC, 1.0);
  multiTrial.setProbability(Move::Types::SwapMultiTrialCBMC, 1.0);
  MCMoveProbabilities cfcmc = MCMoveProbabilities();
  cfcmc.setProbability(Move::Types::Translation, 1.0);
  cfcmc.setProbability(Move::Types::Rotation, 1.0);
  cfcmc.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 0, false, false)}, {}, {}, 5, 21, multiTrial,
                                std::nullopt, false);
  methane.numberOfTrialMolecules = 5;
  Component co2 = makeCO2(forceField, cfcmc, 1);

  System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {methane, co2}, {}, {10, 5}, 5);
  system.cbmcDevice = device.value();

  MonteCarlo mc = makeMonteCarlo({system}, 20, 5, 5);
  mc.run();

  for (System &s : mc.systems)
  {
    EXPECT_EQ(s.cbmcDevice, device.value());
    expectNoEnergyDrift(s);
    const MoveStatistics<double3> &swapStatistics =
        std::get<MoveStatistics<double3>>(s.components[0].mc_moves_statistics[Move::Types::SwapMultiTrialCBMC]);
    EXPECT_GT(swapStatistics.totalCounts.x, 0.0);
    EXPECT_GT(swapStatistics.totalCounts.y, 0.0);
    EXPECT_GT(swapStatistics.totalAccepted.x + swapStatistics.totalAccepted.y, 0.0);
    const MoveStatistics<double> &reinsertionStatistics =
        std::get<MoveStatistics<double>>(s.components[0].mc_moves_statistics[Move::Types::ReinsertionMultiTrialCBMC]);
    EXPECT_GT(reinsertionStatistics.totalCounts, 0.0);
    EXPECT_GT(reinsertionStatistics.totalAccepted, 0.0);
  }
}

// The device proposal with its double-precision correction samples the same grand-canonical distribution as the
// CPU grower: the average loading of methane in FAU agrees within the statistical error.
TEST(MC_MULTI_TRIAL_CBMC, device_loading_matches_cpu_grower_methane_fau)
{
  const std::optional<PairDevice> device = availableCBMCDevice();
  if (!device.has_value()) GTEST_SKIP() << "no CBMC device available";

  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  auto run = [&](PairDevice cbmcDevice) -> std::pair<double, double>
  {
    MCMoveProbabilities probabilities = MCMoveProbabilities();
    probabilities.setProbability(Move::Types::Translation, 1.0);
    probabilities.setProbability(Move::Types::SwapMultiTrialCBMC, 2.0);

    Component methane = Component(forceField, "methane", 190.564, 45599200, 0.01142,
                                  {Atom({0, 0, 0}, 0.0, 1.0, 0, 6, 0, false, false)}, {}, {}, 5, 21, probabilities,
                                  std::nullopt, false);
    methane.numberOfTrialMolecules = 16;

    System system = System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {f}, {methane}, {}, {10}, 5);
    system.cbmcDevice = cbmcDevice;
    MonteCarlo mc = makeMonteCarlo({system}, 400, 100, 0);
    mc.run();
    return mc.systems[0].averageLoadings.averageLoadingNumberOfMolecules(0);
  };

  const auto [meanCPU, errorCPU] = run(PairDevice::CPU);
  const auto [meanDevice, errorDevice] = run(device.value());

  EXPECT_GT(meanCPU, 1.0);
  EXPECT_GT(meanDevice, 1.0);
  EXPECT_NEAR(meanDevice, meanCPU, 3.0 * (errorCPU + errorDevice) + 0.05 * meanCPU)
      << "CPU: " << meanCPU << " +/- " << errorCPU << ", device: " << meanDevice << " +/- " << errorDevice;
}
