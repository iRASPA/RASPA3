#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "input_reader_fixtures.hpp"
#include "molecule_fixtures.hpp"

import std;

import int3;
import double3;
import json;
import simulationbox;
import randomnumbers;
import forcefield;
import framework;
import component;
import molecule;
import atom;
import system;
import cbmc;
import input_reader;
import dctst;
import dctst_settings;
import dctst_barrier;
import dctst_free_energy;
import dctst_free_energy_settings;
import cross_links;

// The seam of the box on the barrier top, the committor of an asymmetric profile, the components the trajectories
// hold in place, and the settings the drivers check for themselves.

namespace
{
class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::current_path(path);
  }
  ~ScopedCurrentPath() { std::filesystem::current_path(original_); }

 private:
  std::filesystem::path original_;
};

nlohmann::json transmissionInput()
{
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "DynamicallyCorrectedTST";
  input["Components"][0].erase("SwapProbability");
  input["RandomSeed"] = 17;
  input["Systems"][0]["dcTST"] = {{"TaggedComponent", "methane"},
                                  {"BarrierPosition", {0.25, 0.0, 0.0}},
                                  {"BarrierNormal", {1.0, 0.0, 0.0}},
                                  {"LatticeDistance", 30.69253},
                                  {"TrajectoryTime", 0.1}};
  return input;
}

void writeInput(const nlohmann::json& input)
{
  std::ofstream stream("simulation.json");
  stream << input.dump(2);
}

System makeCO2System(std::size_t numberOfMolecules)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  const Component co2 = Component::makeCO2(forceField, 0, true);
  return System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {numberOfMolecules}, 5);
}

DcTSTFreeEnergySettings makeFreeEnergySettings()
{
  DcTSTFreeEnergySettings settings;
  settings.taggedComponentName = "CO2";
  settings.taggedComponent = 0;
  settings.barrierPosition = double3(0.25, 0.0, 0.0);
  settings.barrierNormal = double3(1.0, 0.0, 0.0);
  settings.latticeDistance = 11.8671;
  settings.flatteningBias = false;
  settings.numberOfPilotCycles = 0;
  settings.numberOfReplicas = 1;
  return settings;
}
}  // namespace

// The lattice distance is the whole fraction of the repeat distance of the box, not the number as typed: with
// 12.28 for a period of 49.128 the seam of the box at +-period/2 folded 0.004 A beside the barrier plane, and the
// slab between them had five preimages where every other point of the reduced coordinate has four.
TEST(DCTST_SEAM_AND_COMMITTOR, the_lattice_distance_is_a_whole_fraction_of_the_period)
{
  const double period = 49.128;
  const SimulationBox box(period, period, period);
  for (const double typed : {12.28, 12.284})
  {
    const DcTSTBarrier barrier(box, double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), typed);
    EXPECT_DOUBLE_EQ(barrier.lambda, period / 4.0);

    // the reduced coordinate is continuous across the seam of the box, which is a barrier plane
    const double below = barrier.reducedCoordinate(barrier.projection(box, double3(0.5 * period - 1e-6, 3.0, 4.0)));
    const double above = barrier.reducedCoordinate(barrier.projection(box, double3(0.5 * period + 1e-6, 3.0, 4.0)));
    EXPECT_NEAR(below, -1e-6, 1e-9);
    EXPECT_NEAR(above, 1e-6, 1e-9);

    // every reduced coordinate has the same number of preimages: a uniform scan of the box fills the bins evenly
    constexpr std::size_t numberOfBins = 3071;  // 0.004 A per bin
    std::vector<std::size_t> histogram(numberOfBins, 0uz);
    constexpr std::size_t numberOfPoints = 3071uz * 4uz * 16uz;
    for (std::size_t i = 0; i < numberOfPoints; ++i)
    {
      const double x = (static_cast<double>(i) + 0.5) * period / static_cast<double>(numberOfPoints);
      ++histogram[barrier.bin(barrier.projection(box, double3(x, 0.0, 0.0)), numberOfBins)];
    }
    const auto [lowest, highest] = std::minmax_element(histogram.begin(), histogram.end());
    EXPECT_LE(*highest - *lowest, 2uz) << "typed lattice distance " << typed;
  }
}

// A profile of 100 points fills every fourth bin of the committor table. The bins between took the mirror image
// across the barrier before the interpolation, so that an asymmetric barrier was symmetrized under the exponential.
TEST(DCTST_SEAM_AND_COMMITTOR, the_committor_of_an_asymmetric_profile_keeps_its_asymmetry)
{
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  workspace.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);
  ScopedCurrentPath currentPath(workspace.path());
  writeInput(transmissionInput());
  InputReader reader("simulation.json");
  DynamicallyCorrectedTST driver(reader);

  // a barrier of 12 kT that is narrow towards cage A and wide towards cage B
  const double lambda = driver.barrier.lambda;
  const auto betaF = [&](double q)
  {
    const double width = (q < 0.0) ? 0.04 * lambda : 0.16 * lambda;
    return 12.0 * std::exp(-0.5 * q * q / (width * width));
  };
  for (const std::size_t numberOfProfilePoints : {100uz, 200uz, 400uz})
  {
    std::vector<double> q, f;
    for (std::size_t m = 0; m < numberOfProfilePoints; ++m)
    {
      q.push_back(-0.5 * lambda + static_cast<double>(m) * lambda / static_cast<double>(numberOfProfilePoints));
      f.push_back(betaF(q.back()));
    }
    driver.buildCommittorTable(q, f, 0.0);
    ASSERT_EQ(driver.committorTable.size(), 400uz);

    // the committor from the function itself
    constexpr std::size_t fine = 40000;
    std::vector<double> cumulative(fine + 1uz, 0.0);
    for (std::size_t i = 0; i < fine; ++i)
    {
      const double x = -0.5 * lambda + (static_cast<double>(i) + 0.5) * lambda / static_cast<double>(fine);
      cumulative[i + 1uz] = cumulative[i] + std::exp(betaF(x) - 12.0);
    }
    for (std::size_t i = 0; i < 400uz; ++i)
    {
      const double reference = cumulative[(2uz * i + 1uz) * fine / 800uz] / cumulative.back();
      EXPECT_NEAR(driver.committorTable[i], reference, 0.02) << "bin " << i << " of " << numberOfProfilePoints;
    }
    // at the top: the smaller part of the weight lies on the narrow side (1/2 for the symmetrized profile)
    EXPECT_LT(driver.committorTable[200], 0.3);
  }

  // a profile that ends beyond the barrier still takes the missing bins from its mirror image
  {
    std::vector<double> q, f;
    for (std::size_t m = 0; m < 70uz; ++m)
    {
      q.push_back(-0.5 * lambda + static_cast<double>(m) * lambda / 100.0);
      f.push_back(12.0 * std::exp(-0.5 * q.back() * q.back() / (0.01 * lambda * lambda)));
    }
    driver.buildCommittorTable(q, f, 0.0);
    ASSERT_EQ(driver.committorTable.size(), 400uz);
    EXPECT_NEAR(driver.committorTable[200], 0.5, 0.02);
    EXPECT_NEAR(driver.committorTable[399], 1.0, 1e-3);
  }
}

// A component without a Monte Carlo move is part of the fixed environment of the constrained ensemble, and of the
// trajectories: it was given thermal velocities and integrated. The trajectories conserve their energy with such a
// component held in place.
TEST(DCTST_SEAM_AND_COMMITTOR, trajectories_hold_the_components_without_moves_in_place)
{
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  workspace.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);
  ScopedCurrentPath currentPath(workspace.path());
  writeInput(transmissionInput());
  InputReader reader("simulation.json");
  DynamicallyCorrectedTST driver(reader);

  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component cation = Component::makeMethane(forceField, 0);
  cation.name = "cation";
  cation.type = Component::Type::Cation;
  Component methane = Component::makeMethane(forceField, 1);
  driver.systems[0] =
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {cation, methane}, {}, {8, 8}, 5);
  driver.settings.taggedComponent = 1;
  driver.barrier =
      DcTSTBarrier(driver.systems[0].simulationBox, double3(0.25, 0.0, 0.0), double3(1.0, 0.0, 0.0), 11.8671);
  driver.cageLattice =
      DcTSTCageLattice(driver.systems[0].framework->simulationBox, driver.systems[0].framework->numberOfUnitCells,
                       driver.barrier.point + 0.5 * driver.barrier.lambda * driver.barrier.normal);
  driver.setup();
  driver.createConstrainedConfiguration(0);

  System& system = driver.systems[0];
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
  std::vector<double3> before;
  for (std::size_t m = 0; m < 8uz; ++m) before.push_back(system.spanOfMolecule(0, m).front().position);

  driver.launchTrajectories(0, 0, 0);

  const DcTSTTransmissionAccumulator& accumulator = driver.transmission[0];
  ASSERT_GT(accumulator.numberOfTrajectories, 0uz);
  EXPECT_LT(accumulator.maximumRelativeEnergyDrift, 1e-3);
  // (the trajectories run on copies: the Monte Carlo state is untouched either way)
  for (std::size_t m = 0; m < 8uz; ++m)
  {
    EXPECT_DOUBLE_EQ(system.spanOfMolecule(0, m).front().position.x, before[m].x);
  }
}

// The handover times belong to the event-terminated trajectories: settings filled in by hand kept the default list
// with "FixedTime", which switched the hybrid estimator on with accumulators that are never filled.
TEST(DCTST_SEAM_AND_COMMITTOR, handover_times_are_dropped_for_fixed_time_trajectories)
{
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  workspace.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);
  ScopedCurrentPath currentPath(workspace.path());
  writeInput(transmissionInput());
  InputReader reader("simulation.json");
  reader.dctstSettings.trajectoryTermination = DcTSTSettings::TrajectoryTermination::FixedTime;
  reader.dctstSettings.handoverTimes = {0.02, 0.05};
  DynamicallyCorrectedTST driver(reader);
  EXPECT_TRUE(driver.settings.handoverTimes.empty());
  EXPECT_TRUE(driver.transmission[0].hybridOccupancySigned.empty());
}

// A grown molecule carries the component it was grown for, and the molecule records the component of their place in
// the system, whatever id the component was constructed with.
TEST(DCTST_SEAM_AND_COMMITTOR, grown_molecules_and_molecule_records_carry_their_component)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  // both components constructed with id 0
  const Component first = Component::makeMethane(forceField, 0);
  Component second = Component::makeMethane(forceField, 0);
  second.name = "second";
  System system(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {first, second}, {}, {2, 3}, 5);

  for (const Molecule& molecule : system.moleculeData)
  {
    const std::size_t expected = (&molecule - system.moleculeData.data() < 2) ? 0uz : 1uz;
    EXPECT_EQ(molecule.componentId, expected);
  }

  RandomNumber random(3);
  std::optional<CBMC::GrowResult> grown;
  for (std::size_t attempt = 0; attempt < 1000uz && !grown.has_value(); ++attempt)
  {
    grown = CBMC::growNewMolecule(
        random, system.makeGrowContext(CBMC::CutOffMode::Full), system.components[1],
        CBMC::NewMoleculeIdentity{.componentId = 1uz, .moleculeId = system.numberOfMolecules()});
  }
  ASSERT_TRUE(grown.has_value());
  for (const Atom& atom : grown->atoms) EXPECT_EQ(static_cast<std::size_t>(atom.componentId), 1uz);
  EXPECT_EQ(grown->molecule.componentId, 1uz);
}

// The settings the input reader checks, checked by the free-energy driver for a construction by hand: zero blocks
// indexed a histogram that does not exist, zero bins took a modulo by zero.
TEST(DCTST_FREE_ENERGY_CHECKS, programmatic_settings_are_held_to_the_rules_of_the_input_reader)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  auto build = [](auto&& change, std::size_t numberOfBlocks = 5uz)
  {
    DcTSTFreeEnergySettings settings = makeFreeEnergySettings();
    change(settings);
    DcTSTFreeEnergyParameters parameters;
    parameters.numberOfBlocks = numberOfBlocks;
    parameters.randomSeed = 17;
    DynamicallyCorrectedTSTFreeEnergy driver(makeCO2System(8), settings, parameters);
  };
  EXPECT_NO_THROW(build([](DcTSTFreeEnergySettings&) {}));
  EXPECT_THROW(build([](DcTSTFreeEnergySettings&) {}, 0uz), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.numberOfFreeEnergyBins = 0; }), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.whamTolerance = 0.0; }), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.whamMaximumIterations = 0; }), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.reinsertionTaggedProbability = 1.5; }), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s) { s.taggedMoveFraction = -0.1; }), std::runtime_error);
  EXPECT_THROW(build([](DcTSTFreeEnergySettings& s)
                     {
                       s.mirrorSymmetry = true;
                       s.cageOccupancies = std::make_pair(1uz, 0uz);
                     }),
               std::runtime_error);
}

// The plain reaction coordinate is periodic: explicit window centers are points of a circle (the gap across
// +-lambda/2 counts for the stiffness), and a 'WindowRange' of the whole period is that period, not an interval
// with a window on both of its ends.
TEST(DCTST_FREE_ENERGY_CHECKS, windows_of_the_plain_coordinate_are_laid_out_on_the_period)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());
  DcTSTFreeEnergyParameters parameters;
  parameters.numberOfBlocks = 5;
  parameters.randomSeed = 17;

  {
    DcTSTFreeEnergySettings settings = makeFreeEnergySettings();
    settings.windowCenters = {-5.5, 5.5};
    DynamicallyCorrectedTSTFreeEnergy driver(makeCO2System(8), settings, parameters);
    driver.setup();
    driver.buildWindows();
    const double gap = driver.barrier.lambda - 11.0;
    ASSERT_EQ(driver.windows.size(), 2uz);
    EXPECT_NEAR(driver.windows[0].stiffness, 4.0 / (gap * gap), 1e-9);
    EXPECT_TRUE(driver.windowsClosePeriod);
  }
  {
    // centers outside the period are folded into it
    DcTSTFreeEnergySettings settings = makeFreeEnergySettings();
    settings.windowCenters = {0.0, 7.0, 10.0};
    DynamicallyCorrectedTSTFreeEnergy driver(makeCO2System(8), settings, parameters);
    driver.setup();
    driver.buildWindows();
    const double lambda = driver.barrier.lambda;
    ASSERT_EQ(driver.windows.size(), 3uz);
    EXPECT_NEAR(driver.windows[1].center, 7.0 - lambda, 1e-9);
    EXPECT_NEAR(driver.windows[2].center, 10.0 - lambda, 1e-9);
    EXPECT_NEAR(driver.windows[0].stiffness, 4.0 / ((lambda - 10.0) * (lambda - 10.0)), 1e-9);
  }
  {
    // the whole period as a range: no window twice on the same plane
    DcTSTFreeEnergySettings settings = makeFreeEnergySettings();
    settings.numberOfWindows = 4;
    DynamicallyCorrectedTSTFreeEnergy probe(makeCO2System(8), settings, parameters);
    const double lambda = probe.barrier.lambda;
    settings.windowRange = std::make_pair(0.0, lambda);
    DynamicallyCorrectedTSTFreeEnergy driver(makeCO2System(8), settings, parameters);
    driver.setup();
    driver.buildWindows();
    ASSERT_EQ(driver.windows.size(), 4uz);
    for (std::size_t i = 0; i < 4uz; ++i)
    {
      for (std::size_t j = i + 1uz; j < 4uz; ++j)
      {
        EXPECT_GT(std::abs(driver.barrier.reducedCoordinate(driver.windows[i].center - driver.windows[j].center)),
                  0.2 * lambda);
      }
    }
    EXPECT_TRUE(driver.windowsClosePeriod);
  }
  {
    // windows over a part of the period do not go around it
    DcTSTFreeEnergySettings settings = makeFreeEnergySettings();
    settings.numberOfWindows = 5;
    settings.windowRange = std::make_pair(-1.0, 1.0);
    DynamicallyCorrectedTSTFreeEnergy driver(makeCO2System(8), settings, parameters);
    driver.setup();
    driver.buildWindows();
    ASSERT_EQ(driver.windows.size(), 5uz);
    EXPECT_FALSE(driver.windowsClosePeriod);
  }
}

// The transmission coefficient takes the speed of the tagged bead from its effective mass along the normal: a FIXED
// bond of the molecule (RATTLE) moves its partner along, which that mass of the free bead left out. And bonds between
// molecules tie the tagged molecule to others. Both are refused by the free-energy driver.
TEST(DCTST_FREE_ENERGY_CHECKS, constrained_or_linked_tagged_molecules_are_refused)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());
  DcTSTFreeEnergyParameters parameters;
  parameters.numberOfBlocks = 5uz;
  parameters.randomSeed = 17;

  System constrained = makeCO2System(8);
  constrained.components[0].distanceConstraints.push_back(Component::DistanceConstraint{0, 1, 1.149});
  EXPECT_THROW(DynamicallyCorrectedTSTFreeEnergy(constrained, makeFreeEnergySettings(), parameters),
               std::runtime_error);

  System linked = makeCO2System(8);
  linked.crossLinks.bondTypes.push_back(CrossLinkBondType{});
  linked.crossLinks.links.push_back(CrossLink{{0, 0, 0}, {0, 1, 0}, 0});
  linked.crossLinks.rebuildIndex();
  EXPECT_THROW(DynamicallyCorrectedTSTFreeEnergy(linked, makeFreeEnergySettings(), parameters), std::runtime_error);
}
