#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "input_reader_fixtures.hpp"
#include "molecule_fixtures.hpp"

import std;

import input_reader;
import json;
import thermobarostat;

namespace
{

class TemporaryInput
{
 public:
  TemporaryInput(nlohmann::json input, std::string_view suffix, const std::filesystem::path& directory)
      : path_(directory / std::format("simulation_tmmc_{}_{}.json", suffix,
                                      std::chrono::steady_clock::now().time_since_epoch().count()))
  {
    std::ofstream stream(path_);
    stream << input.dump(2);
  }

  ~TemporaryInput()
  {
    std::error_code error;
    std::filesystem::remove(path_, error);
  }

  const std::filesystem::path& path() const { return path_; }

 private:
  std::filesystem::path path_;
};

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

TemporaryDirectory makeTmmcWorkspace()
{
  TemporaryDirectory dir;
  dir.write("force_field.json", input_reader_fixtures::kTmmcForceFieldJson);
  dir.write("tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  dir.write("methane.json", molecule_fixtures::kMethaneJson);
  return dir;
}

TemporaryDirectory makeBoxMdWorkspace()
{
  TemporaryDirectory dir;
  dir.write("force_field.json", input_reader_fixtures::kBoxForceFieldJson);
  dir.write("methane.json", molecule_fixtures::kMethaneJson);
  return dir;
}

nlohmann::json readTMMCExample() { return nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson); }

nlohmann::json readNVTExample() { return nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson); }

}  // namespace

TEST(INPUT_READER_TMMC, rejects_multidimensional_component_moves)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["Components"][0]["IdentityChangeProbability"] = 1.0;
  TemporaryInput temporary(std::move(input), "identity", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}

TEST(INPUT_READER_TMMC, rejects_non_nearest_neighbor_system_moves)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["Systems"][0]["ParallelTemperingSwapProbability"] = 1.0;
  TemporaryInput temporary(std::move(input), "parallel_tempering", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}

TEST(INPUT_READER_TMMC, accepts_elastic_constant_minimization_options)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["ComputeElasticConstants"] = true;
  input["ElasticEigenvalueTolerance"] = 2.5e-9;
  TemporaryInput temporary(std::move(input), "elastic_constants", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.minimizationOptions.computeElasticConstants);
  EXPECT_DOUBLE_EQ(reader.minimizationOptions.elasticEigenvalueTolerance, 2.5e-9);
}

TEST(INPUT_READER_MD, accepts_nvt_stress_fluctuation_elastic_constants)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = readNVTExample();
  input["NumberOfBlocks"] = 7;
  input["Systems"][0]["ComputeElasticConstantsFromFluctuations"] = true;
  input["Systems"][0]["ElasticConstantsSampleEvery"] = 25;
  TemporaryInput temporary(std::move(input), "elastic_fluctuation", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_TRUE(reader.systems[0].propertyElasticConstantsFluctuation.has_value());
  EXPECT_EQ(reader.numberOfBlocks, 7u);
  EXPECT_EQ(reader.systems[0].propertyElasticConstantsFluctuation->numberOfBlocks, 7u);
  EXPECT_EQ(reader.systems[0].elasticConstantsSampleEvery, 25u);
}

TEST(INPUT_READER_MD, rejects_non_nvt_stress_fluctuation_elastic_constants)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = readNVTExample();
  input["Systems"][0]["Ensemble"] = "NVE";
  input["Systems"][0]["ComputeElasticConstantsFromFluctuations"] = true;
  TemporaryInput temporary(std::move(input), "elastic_fluctuation_nve", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}

TEST(INPUT_READER_MD, accepts_grand_canonical_molecular_dynamics_ensembles)
{
  struct Case
  {
    std::string name;
    MolecularDynamicsEnsemble ensemble;
    bool hasBarostat;
  };
  const std::array cases{
      Case{"MuVT", MolecularDynamicsEnsemble::MuVT, false},
      Case{"MuPT", MolecularDynamicsEnsemble::MuPT, true},
      Case{"MuPTPR", MolecularDynamicsEnsemble::MuPTPR, true},
  };

  for (const Case& test : cases)
  {
    TemporaryDirectory workspace = makeBoxMdWorkspace();
    nlohmann::json input = readNVTExample();
    input["Systems"][0]["Ensemble"] = test.name;
    input["Systems"][0]["ExternalPressure"] = 1.0e5;
    if (test.ensemble == MolecularDynamicsEnsemble::MuPTPR)
      input["Systems"][0]["CellType"] = "RegularUpperTriangle";
    input["Components"][0]["SwapProbability"] = 1.0;
    TemporaryInput temporary(std::move(input), test.name, workspace.path());
    ScopedCurrentPath currentPath(workspace.path());

    InputReader reader(temporary.path().filename().string());
    EXPECT_EQ(reader.systems[0].molecularDynamicsEnsemble, test.ensemble);
    EXPECT_TRUE(reader.systems[0].thermostat.has_value());
    EXPECT_EQ(reader.systems[0].thermobarostat.has_value(), test.hasBarostat);
  }
}

TEST(INPUT_READER_MD, rejects_grand_canonical_md_without_reservoir_or_swap)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();

  nlohmann::json missingPressure = readNVTExample();
  missingPressure["Systems"][0]["Ensemble"] = "MuVT";
  missingPressure["Components"][0]["SwapProbability"] = 1.0;
  TemporaryInput pressureInput(std::move(missingPressure), "muvt_pressure", workspace.path());

  nlohmann::json missingSwap = readNVTExample();
  missingSwap["Systems"][0]["Ensemble"] = "MuVT";
  missingSwap["Systems"][0]["ExternalPressure"] = 1.0e5;
  TemporaryInput swapInput(std::move(missingSwap), "muvt_swap", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader pressureReader(pressureInput.path().filename().string()), std::runtime_error);
  EXPECT_THROW(InputReader swapReader(swapInput.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER_MD, reads_grand_canonical_md_examples)
{
  const std::array simulations{
      input_reader_fixtures::kMuvtSimulationJson,
      input_reader_fixtures::kMuptSimulationJson,
      input_reader_fixtures::kMuptprSimulationJson,
  };
  for (std::string_view simulation : simulations)
  {
    TemporaryDirectory workspace = makeBoxMdWorkspace();
    workspace.write("simulation.json", simulation);
    ScopedCurrentPath currentPath(workspace.path());
    InputReader reader("simulation.json");
    EXPECT_TRUE(molecularDynamicsHasParticleExchange(reader.systems[0].molecularDynamicsEnsemble));
  }
}
