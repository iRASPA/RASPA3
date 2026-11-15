#include <gtest/gtest.h>

#include "molecule_fixtures.hpp"

import std;

import int3;
import double3;
import simd_quatd;
import units;
import json;
import simulationbox;
import molecule;
import fragment;
import mc_moves_probabilities;
import forcefield;
import framework;
import component;
import system;
import running_energy;
import dctst_free_energy;
import dctst_free_energy_settings;

// The biased moves of the tagged molecule of the dcTST free-energy driver with a polarizable force field. The
// window cycle hands every other molecule to the standard moves, which track the electric field and the
// polarization energy incrementally; the tagged molecule's own translation and reinsertion have to do the same,
// or the window samples the wrong distribution, the stored field goes stale, and the running energy drifts
// from a full recomputation.

namespace
{

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath()
  {
    std::error_code ignored;
    const std::filesystem::path used = std::filesystem::current_path();
    std::filesystem::current_path(original_);
    std::filesystem::remove_all(used, ignored);
  }

 private:
  std::filesystem::path original_;
};

// CO2 in ITQ-29 with polarization on, molecule-molecule polarization included (the same fixture as the
// polarization drift tests of the standard moves)
System makePolarizableSystem(std::size_t numberOfMolecules)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  const Component co2 = Component::makeCO2(forceField, 0, true);
  return System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {numberOfMolecules}, 5);
}

DcTSTFreeEnergySettings makeSettings()
{
  DcTSTFreeEnergySettings settings;
  settings.taggedComponentName = "CO2";
  settings.taggedComponent = 0;
  settings.barrierPosition = double3(0.25, 0.0, 0.0);  // an 8-ring window of the 2x2x2 supercell
  settings.barrierNormal = double3(1.0, 0.0, 0.0);
  settings.latticeDistance = 11.8671;  // cage-to-cage distance (the unit cell)
  settings.numberOfWindows = 2;        // one window on the barrier, one in the cage
  settings.flatteningBias = false;
  settings.numberOfPilotCycles = 0;
  settings.maximumTranslationTagged = 0.5;
  settings.numberOfReplicas = 1;
  return settings;
}

// runs 'steps' tagged moves in the given window, 'reinsertionFraction' of them reinsertions, and returns the
// drift of the running energy from a full recomputation
RunningEnergy driftOfTaggedMoves(DynamicallyCorrectedTSTFreeEnergy& driver, std::size_t windowIndex,
                                 double reinsertionFraction, std::size_t steps)
{
  System& system = driver.systems[0];
  driver.startWindow(0, windowIndex);
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();

  for (std::size_t step = 0; step < steps; ++step)
  {
    const std::optional<RunningEnergy> energyDifference = (driver.randoms[0].uniform() < reinsertionFraction)
                                                              ? driver.biasedReinsertionMove(0, windowIndex)
                                                              : driver.biasedTranslationMove(0, windowIndex);
    if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();
  }
  return system.runningEnergies - system.computeTotalEnergies();
}

}  // namespace

TEST(DCTST_POLARIZATION_DRIFT, tagged_translation_and_reinsertion_keep_the_polarization_energy_consistent)
{
  ScopedCurrentPath workspace(
      std::filesystem::temp_directory_path() /
      std::format("raspa3_dctst_polarization_{}", std::chrono::steady_clock::now().time_since_epoch().count()));

  DcTSTFreeEnergyParameters parameters;
  parameters.numberOfBlocks = 5;
  parameters.randomSeed = 17;
  DynamicallyCorrectedTSTFreeEnergy driver(makePolarizableSystem(16), makeSettings(), parameters);
  driver.setup();
  driver.buildWindows();
  ASSERT_EQ(driver.windows.size(), 2uz);

  // the fixture must exercise a non-trivial polarization energy
  driver.systems[0].precomputeTotalRigidEnergy();
  EXPECT_LT(driver.systems[0].computeTotalEnergies().polarization, -1e-8);

  for (std::size_t windowIndex = 0; windowIndex < driver.windows.size(); ++windowIndex)
  {
    // translations only, then reinsertions only, then the mix of the window cycle
    for (const double reinsertionFraction : {0.0, 1.0, 0.3})
    {
      const RunningEnergy drift = driftOfTaggedMoves(driver, windowIndex, reinsertionFraction, 400);
      EXPECT_NEAR(drift.polarization, 0.0, 1e-6)
          << "window " << windowIndex << ", reinsertion fraction " << reinsertionFraction;
      EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6) << "window " << windowIndex;
      EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6)
          << "window " << windowIndex << ", reinsertion fraction " << reinsertionFraction;
    }
  }

  // the moves were actually accepted (a test of rejected moves only would pass trivially)
  std::size_t translationsAccepted = 0;
  std::size_t reinsertionsAccepted = 0;
  for (const DcTSTWindowAccumulator& accumulator : driver.windowData)
  {
    translationsAccepted += accumulator.translationAccepted;
    reinsertionsAccepted += accumulator.reinsertionAccepted;
  }
  EXPECT_GT(translationsAccepted, 0uz);
  EXPECT_GT(reinsertionsAccepted, 0uz);
}

// The reaction bead of a semi-flexible molecule may sit inside a rigid fragment. Its thermal speed along the barrier
// normal is then that of the fragment (translation plus the lever of the bead about the fragment's axes), as it is for
// a bead of a rigid molecule, and not the speed of a free atom of the bead's mass.
TEST(DCTST_FREE_ENERGY, the_pilot_bead_speed_of_a_rigid_fragment_bead_uses_the_fragment_effective_mass)
{
  ScopedCurrentPath workspace(
      std::filesystem::temp_directory_path() /
      std::format("raspa3_dctst_pilot_fragment_{}", std::chrono::steady_clock::now().time_since_epoch().count()));

  // semi-flexible pentane with the reaction bead on the middle atom of its rigid CH2-CH2-CH2 core
  nlohmann::json definition = nlohmann::json::parse(molecule_fixtures::kSemiFlexiblePentaneJson);
  definition["StartingBead"] = 2;
  {
    std::ofstream stream("pentane.json");
    stream << definition.dump(2);
  }
  const ForceField forceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                              {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                              true, false, false);
  const Component pentane(Component::Type::Adsorbate, 0, forceField, "pentane", "pentane", 5, 21, MCMoveProbabilities(),
                          std::nullopt, false);
  ASSERT_EQ(pentane.startingBead, 2uz);
  ASSERT_TRUE(pentane.isSemiFlexible());
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {pentane}, {}, {1}, 5);

  DcTSTFreeEnergySettings settings;
  settings.taggedComponentName = "pentane";
  settings.taggedComponent = 0;
  settings.barrierPosition = double3(0.5, 0.0, 0.0);
  settings.barrierNormal = double3(1.0, 0.0, 0.0);
  settings.latticeDistance = 30.0;
  settings.numberOfWindows = 2;
  settings.flatteningBias = false;
  settings.numberOfPilotCycles = 0;
  settings.numberOfReplicas = 1;
  DcTSTFreeEnergyParameters parameters;
  parameters.numberOfBlocks = 5;
  parameters.randomSeed = 17;
  DynamicallyCorrectedTSTFreeEnergy driver(system, settings, parameters);
  driver.setup();

  System& replica = driver.systems[0];
  driver.samplePilot(0, 0);
  ASSERT_EQ(driver.pilot[0].beadSpeedCount, 1.0);

  // the effective mass of the bead along the normal, from the fragment's rigid-body state
  const Component& component = replica.components[0];
  const std::optional<std::size_t> fragmentIndex = component.rigidFragmentContaining(component.startingBead);
  ASSERT_TRUE(fragmentIndex.has_value());
  const Fragment& fragment = component.fragmentGraph.fragments[*fragmentIndex];
  const GroupState state = component.deriveFragmentState(*fragmentIndex, replica.spanOfMolecule(0, 0));
  const std::size_t k =
      static_cast<std::size_t>(std::ranges::find(fragment.atoms, component.startingBead) - fragment.atoms.begin());
  const double3 normalBody = state.orientation.inverse() * driver.barrier.normal;
  const double3 lever = double3::cross(fragment.bodyFixedPositions[k], normalBody);
  const double3 invI = fragment.inverseInertiaVector;
  const double inverseMass =
      1.0 / fragment.mass + lever.x * lever.x * invI.x + lever.y * lever.y * invI.y + lever.z * lever.z * invI.z;
  const double expected = std::sqrt(2.0 * Units::KB * replica.temperature * inverseMass / std::numbers::pi);
  EXPECT_NEAR(driver.pilot[0].beadSpeedSum, expected, 1e-9 * expected);
  // (a free CH2 atom would be sqrt(42.08 / 14.03) faster than the fragment's translation alone)
  EXPECT_LT(driver.pilot[0].beadSpeedSum, std::sqrt(2.0 * Units::KB * replica.temperature / 14.03 / std::numbers::pi));
}

namespace
{
System makeMethaneSystem(std::size_t numberOfCations, std::size_t numberOfMolecules)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component cation = Component::makeMethane(forceField, 0);
  cation.name = "cation";
  cation.type = Component::Type::Cation;
  Component methane = Component::makeMethane(forceField, 1);
  methane.name = "methane";
  return System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {cation, methane}, {},
                {numberOfCations, numberOfMolecules}, 5);
}

DcTSTFreeEnergySettings makeMethaneSettings()
{
  DcTSTFreeEnergySettings settings = makeSettings();
  settings.taggedComponentName = "methane";
  settings.taggedComponent = 1;
  return settings;
}

DcTSTFreeEnergyParameters makeParameters()
{
  DcTSTFreeEnergyParameters parameters;
  parameters.numberOfBlocks = 5;
  parameters.randomSeed = 17;
  return parameters;
}
}  // namespace

// The frame of cage A runs over [-lambda/2, 3 lambda/4) with the dividing surface at lambda/2: four fifths of the
// way, a bin edge for one bin count in five only. In the others the pilot filled part of a bin the reweighting
// took as forbidden to it (or as wholly its own).
TEST(DCTST_FREE_ENERGY, the_dividing_surface_of_the_occupancy_frame_is_a_bin_edge)
{
  for (const double fineBinWidth : {0.01, 0.013, 0.02, 0.031, 0.05, 0.3})
  {
    DcTSTFreeEnergySettings settings = makeMethaneSettings();
    settings.cageOccupancies = std::make_pair(1uz, 0uz);
    settings.fineBinWidth = fineBinWidth;
    const DynamicallyCorrectedTSTFreeEnergy driver(makeMethaneSystem(0, 8), settings, makeParameters());
    const double lambda = driver.barrier.lambda;
    const std::size_t numberOfBins = driver.grid.numberOfBins;
    ASSERT_EQ(numberOfBins % 5uz, 0uz) << fineBinWidth;
    const std::size_t first = 4uz * numberOfBins / 5uz;  // the first bin of cage B
    EXPECT_NEAR(driver.grid.begin + static_cast<double>(first) * driver.grid.width, 0.5 * lambda, 1e-9);
    EXPECT_GT(driver.grid.center(first), 0.5 * lambda);
    EXPECT_LT(driver.grid.center(first - 1uz), 0.5 * lambda);
    EXPECT_NEAR(driver.grid.width, fineBinWidth, 0.5 * fineBinWidth);
  }
}

TEST(DCTST_FREE_ENERGY, settings_that_cannot_be_run_are_refused)
{
  auto build = [](auto&& change)
  {
    DcTSTFreeEnergySettings settings = makeMethaneSettings();
    change(settings);
    DynamicallyCorrectedTSTFreeEnergy driver(makeMethaneSystem(0, 8), settings, makeParameters());
  };
  EXPECT_NO_THROW(build([](DcTSTFreeEnergySettings&) {}));
  // n_A counts the tagged molecule: zero wrapped n_A - 1 around
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.cageOccupancies = std::make_pair(0uz, 1uz); }),
               std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.fineBinWidth = 0.0; }), std::runtime_error);
  // a window off the coordinate it biases: every sample of it would be rejected
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s)
                     {
                       s.mirrorSymmetry = true;
                       s.windowCenters = {-1.0, 1.0};
                     }),
               std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s)
                     {
                       s.cageOccupancies = std::make_pair(1uz, 0uz);
                       s.windowCenters = {0.0, 12.0};
                     }),
               std::runtime_error);
  // the plain coordinate is periodic: a center outside the period is folded
  EXPECT_NO_THROW(build([](DcTSTFreeEnergySettings& s) { s.windowCenters = {0.0, 7.0}; }));
}

// A center given twice made the smallest spacing zero: an infinite stiffness with a support of zero width, and
// 0.5 * inf * 0 = NaN in the bias of the bin on the center.
TEST(DCTST_FREE_ENERGY, repeated_window_centers_keep_a_finite_bias)
{
  ScopedCurrentPath workspace(
      std::filesystem::temp_directory_path() /
      std::format("raspa3_dctst_windows_{}", std::chrono::steady_clock::now().time_since_epoch().count()));
  DcTSTFreeEnergySettings settings = makeMethaneSettings();
  settings.windowCenters = {0.0, 0.0, 1.0};
  DynamicallyCorrectedTSTFreeEnergy driver(makeMethaneSystem(0, 8), settings, makeParameters());
  driver.setup();
  driver.buildWindows();

  ASSERT_GE(driver.windows.size(), 3uz);
  for (const auto& window : driver.windows)
  {
    EXPECT_NEAR(window.stiffness, 4.0, 1e-9);  // 4 / spacing^2 with the spacing of the distinct centers
    for (const double bias : window.bias)
    {
      EXPECT_GE(bias, 0.0);
      EXPECT_LT(bias, 1.0e300);
    }
  }

  // the thermal speed is that of the samples on the dividing surface once there are enough of them
  const double cageSpeed = driver.beadThermalSpeed();
  EXPECT_GT(cageSpeed, 0.0);
  driver.windowData[0].barrierSpeedSum = 250.0;
  driver.windowData[0].barrierSpeedCount = 50.0;
  EXPECT_DOUBLE_EQ(driver.beadThermalSpeed(), cageSpeed);
  driver.windowData[0].barrierSpeedCount = 125.0;
  EXPECT_DOUBLE_EQ(driver.beadThermalSpeed(), 2.0);
}

// The cage occupancies are those of the tagged species: a cation in the cage was counted, and taken off its site
// to fill or empty a cage.
TEST(DCTST_FREE_ENERGY, cage_occupancies_count_the_tagged_species_only)
{
  ScopedCurrentPath workspace(
      std::filesystem::temp_directory_path() /
      std::format("raspa3_dctst_occupancy_{}", std::chrono::steady_clock::now().time_since_epoch().count()));
  DcTSTFreeEnergySettings settings = makeMethaneSettings();
  settings.cageOccupancies = std::make_pair(2uz, 1uz);
  DynamicallyCorrectedTSTFreeEnergy driver(makeMethaneSystem(6, 12), settings, makeParameters());
  driver.setup();
  driver.buildWindows();

  System& system = driver.systems[0];
  const std::vector<std::size_t> occupancy = driver.countOccupancies(system);
  EXPECT_EQ(std::accumulate(occupancy.begin(), occupancy.end(), 0uz), 12uz);

  std::vector<double3> before;
  for (std::size_t m = 0; m < 6uz; ++m) before.push_back(system.spanOfMolecule(0, m).front().position);
  driver.startWindow(0, 0);
  for (std::size_t m = 0; m < 6uz; ++m)
  {
    const double3 after = system.spanOfMolecule(0, m).front().position;
    EXPECT_DOUBLE_EQ(after.x, before[m].x);
    EXPECT_DOUBLE_EQ(after.y, before[m].y);
    EXPECT_DOUBLE_EQ(after.z, before[m].z);
  }
  // the requested occupancies hold for the tagged species (the tagged molecule itself is the +1 of cage A)
  const DcTSTWindowState& state = driver.windowStates[0];
  const std::vector<std::size_t> counted = driver.countOccupancies(system);
  EXPECT_EQ(counted[state.cageA], 2uz);
  EXPECT_EQ(counted[state.cageB], 1uz);
}
