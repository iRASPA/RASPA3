#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "input_reader_fixtures.hpp"
#include "molecule_fixtures.hpp"

import std;

import int3;
import double3;
import json;
import forcefield;
import framework;
import component;
import system;
import running_energy;
import input_reader;
import dctst;
import dctst_barrier;

// The in-plane moves of the tagged molecule of the dcTST transmission driver with a polarizable force field. Every
// other molecule moves with the standard moves, which track the electric field and the polarization energy
// incrementally; the tagged molecule's in-plane translation and reinsertion have to do the same, or the constrained
// ensemble is sampled with the wrong energy, the stored field goes stale, and the running energy drifts from a full
// recomputation.

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
}  // namespace

TEST(DCTST_TRANSMISSION_POLARIZATION_DRIFT,
     in_plane_translation_and_reinsertion_keep_the_polarization_energy_consistent)
{
  // the driver is built from an input file (methane in tobacco-667) and then handed the polarizable system: CO2 in
  // ITQ-29 with molecule-molecule polarization, the fixture of the polarization drift tests of the standard moves
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  workspace.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "DynamicallyCorrectedTST";
  input["Components"][0].erase("SwapProbability");
  input["RandomSeed"] = 17;
  input["Systems"][0]["dcTST"] = {{"TaggedComponent", "methane"},
                                  {"BarrierPosition", {0.25, 0.0, 0.0}},
                                  {"BarrierNormal", {1.0, 0.0, 0.0}},
                                  // the (100) plane spacing of the (nearly orthogonal, triclinic) tobacco-667 cell:
                                  // the box must repeat a whole number of lattice distances along the normal
                                  {"LatticeDistance", 30.69253}};
  {
    std::ofstream stream("simulation.json");
    stream << input.dump(2);
  }
  InputReader reader("simulation.json");
  DynamicallyCorrectedTST driver(reader);

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;
  const Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  const Component co2 = Component::makeCO2(forceField, 0, true);
  driver.systems[0] = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {16}, 5);
  driver.settings.taggedComponent = 0;
  driver.barrier =
      DcTSTBarrier(driver.systems[0].simulationBox, double3(0.25, 0.0, 0.0), double3(1.0, 0.0, 0.0), 11.8671);
  driver.cageLattice =
      DcTSTCageLattice(driver.systems[0].framework->simulationBox, driver.systems[0].framework->numberOfUnitCells,
                       driver.barrier.point + 0.5 * driver.barrier.lambda * driver.barrier.normal);
  driver.constrainedMoves[0].maximumTranslation = 0.5;
  driver.setup();
  driver.createConstrainedConfiguration(0);

  System& system = driver.systems[0];
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
  // the fixture must exercise a non-trivial polarization energy
  ASSERT_LT(system.runningEnergies.polarization, -1e-8);

  // translations only, then reinsertions only, then a mix
  for (const double reinsertionFraction : {0.0, 1.0, 0.3})
  {
    for (std::size_t step = 0; step < 400uz; ++step)
    {
      const std::optional<RunningEnergy> energyDifference = (driver.randoms[0].uniform() < reinsertionFraction)
                                                                ? driver.reinsertionInPlaneMove(0)
                                                                : driver.translationInPlaneMove(0);
      if (energyDifference.has_value()) system.runningEnergies += energyDifference.value();
    }
    const RunningEnergy drift = system.runningEnergies - system.computeTotalEnergies();
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6) << "reinsertion fraction " << reinsertionFraction;
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6) << "reinsertion fraction " << reinsertionFraction;
  }

  // the moves were actually accepted (a test of rejected moves only would pass trivially)
  EXPECT_GT(driver.constrainedMoves[0].translationAccepted, 0uz);
  EXPECT_GT(driver.constrainedMoves[0].reinsertionAccepted, 0uz);
}

namespace
{
nlohmann::json transmissionInput()
{
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "DynamicallyCorrectedTST";
  input["Components"][0].erase("SwapProbability");
  input["RandomSeed"] = 17;
  input["Systems"][0]["dcTST"] = {{"TaggedComponent", "methane"},
                                  {"BarrierPosition", {0.25, 0.0, 0.0}},
                                  {"BarrierNormal", {1.0, 0.0, 0.0}},
                                  {"LatticeDistance", 30.69253}};
  return input;
}

void writeInput(const nlohmann::json& input)
{
  std::ofstream stream("simulation.json");
  stream << input.dump(2);
}
}  // namespace

// The extra-framework cations, and any component without a move, stay where the system has them when the
// constrained configuration is built: every component was deleted and regrown at random positions.
TEST(DCTST_TRANSMISSION_POLARIZATION_DRIFT, constrained_configuration_keeps_the_cations_in_place)
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
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {cation, methane}, {}, {4, 8}, 5);
  driver.settings.taggedComponent = 1;
  driver.barrier =
      DcTSTBarrier(driver.systems[0].simulationBox, double3(0.25, 0.0, 0.0), double3(1.0, 0.0, 0.0), 11.8671);
  driver.cageLattice =
      DcTSTCageLattice(driver.systems[0].framework->simulationBox, driver.systems[0].framework->numberOfUnitCells,
                       driver.barrier.point + 0.5 * driver.barrier.lambda * driver.barrier.normal);
  driver.setup();

  System& system = driver.systems[0];
  std::vector<double3> before;
  for (std::size_t m = 0; m < 4uz; ++m) before.push_back(system.spanOfMolecule(0, m).front().position);

  driver.createConstrainedConfiguration(0);

  ASSERT_EQ(system.numberOfMoleculesPerComponent[0], 4uz);
  ASSERT_EQ(system.numberOfMoleculesPerComponent[1], 8uz);
  for (std::size_t m = 0; m < 4uz; ++m)
  {
    const double3 after = system.spanOfMolecule(0, m).front().position;
    EXPECT_DOUBLE_EQ(after.x, before[m].x);
    EXPECT_DOUBLE_EQ(after.y, before[m].y);
    EXPECT_DOUBLE_EQ(after.z, before[m].z);
  }
  // the tagged molecule is on the plane
  const double3 bead = system.spanOfMolecule(1, 0).front().position;
  EXPECT_NEAR(driver.barrier.reducedCoordinate(driver.barrier.projection(system.simulationBox, bead)), 0.0, 1e-8);

  // the cage occupancies are those of the tagged species
  const std::vector<std::size_t> occupancy = driver.countOccupancies(system);
  EXPECT_EQ(std::accumulate(occupancy.begin(), occupancy.end(), 0uz), 8uz);
}

// A commit distance beyond the center of the cage, and a plateau window outside the trajectories, were accepted:
// the first committed the tagged molecule to the wrong cage, the second was replaced by the last time point.
TEST(DCTST_TRANSMISSION_POLARIZATION_DRIFT, commit_distance_and_plateau_window_are_held_to_their_ranges)
{
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  workspace.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);
  ScopedCurrentPath currentPath(workspace.path());

  auto build = [](auto&& change)
  {
    nlohmann::json input = transmissionInput();
    change(input["Systems"][0]["dcTST"]);
    writeInput(input);
    InputReader reader("simulation.json");
    DynamicallyCorrectedTST driver(reader);
  };
  EXPECT_NO_THROW(build([](nlohmann::json&) {}));
  EXPECT_NO_THROW(build([](nlohmann::json& block) { block["CommitDistance"] = 10.0; }));
  EXPECT_THROW(build([](nlohmann::json& block) { block["CommitDistance"] = 20.0; }), std::runtime_error);
  EXPECT_THROW(build([](nlohmann::json& block) { block["PlateauWindow"] = {5.0, 2.0}; }), std::runtime_error);
  EXPECT_THROW(build([](nlohmann::json& block) { block["PlateauWindow"] = {2.0, 1.0e6}; }), std::runtime_error);
}
