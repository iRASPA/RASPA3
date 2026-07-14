#include <gtest/gtest.h>

import std;

import input_reader;
import json;

namespace
{

const std::filesystem::path tmmcExample = "examples/advanced/5_tmmc_methane_in_tobacco_667/0/simulation.json";
const std::filesystem::path nvtExample = "examples/basic/5_md_methane_in_box_msd/simulation.json";

class TemporaryInput
{
 public:
  TemporaryInput(nlohmann::json input, std::string_view suffix,
                 const std::filesystem::path& directory = tmmcExample.parent_path())
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

nlohmann::json readTMMCExample()
{
  std::ifstream stream(tmmcExample);
  return nlohmann::json::parse(stream);
}

nlohmann::json readNVTExample()
{
  std::ifstream stream(nvtExample);
  return nlohmann::json::parse(stream);
}

}  // namespace

TEST(INPUT_READER_TMMC, rejects_multidimensional_component_moves)
{
  nlohmann::json input = readTMMCExample();
  input["Components"][0]["IdentityChangeProbability"] = 1.0;
  TemporaryInput temporary(std::move(input), "identity");
  ScopedCurrentPath currentPath(tmmcExample.parent_path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}

TEST(INPUT_READER_TMMC, rejects_non_nearest_neighbor_system_moves)
{
  nlohmann::json input = readTMMCExample();
  input["Systems"][0]["ParallelTemperingSwapProbability"] = 1.0;
  TemporaryInput temporary(std::move(input), "parallel_tempering");
  ScopedCurrentPath currentPath(tmmcExample.parent_path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}

TEST(INPUT_READER_TMMC, accepts_elastic_constant_minimization_options)
{
  nlohmann::json input = readTMMCExample();
  input["ComputeElasticConstants"] = true;
  input["ElasticEigenvalueTolerance"] = 2.5e-9;
  TemporaryInput temporary(std::move(input), "elastic_constants");
  ScopedCurrentPath currentPath(tmmcExample.parent_path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.minimizationOptions.computeElasticConstants);
  EXPECT_DOUBLE_EQ(reader.minimizationOptions.elasticEigenvalueTolerance, 2.5e-9);
}

TEST(INPUT_READER_MD, accepts_nvt_stress_fluctuation_elastic_constants)
{
  nlohmann::json input = readNVTExample();
  input["NumberOfBlocks"] = 7;
  input["Systems"][0]["ComputeElasticConstantsFromFluctuations"] = true;
  input["Systems"][0]["ElasticConstantsSampleEvery"] = 25;
  TemporaryInput temporary(std::move(input), "elastic_fluctuation", nvtExample.parent_path());
  ScopedCurrentPath currentPath(nvtExample.parent_path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_TRUE(reader.systems[0].propertyElasticConstantsFluctuation.has_value());
  EXPECT_EQ(reader.numberOfBlocks, 7u);
  EXPECT_EQ(reader.systems[0].propertyElasticConstantsFluctuation->numberOfBlocks, 7u);
  EXPECT_EQ(reader.systems[0].elasticConstantsSampleEvery, 25u);
}

TEST(INPUT_READER_MD, rejects_non_nvt_stress_fluctuation_elastic_constants)
{
  nlohmann::json input = readNVTExample();
  input["Systems"][0]["Ensemble"] = "NVE";
  input["Systems"][0]["ComputeElasticConstantsFromFluctuations"] = true;
  TemporaryInput temporary(std::move(input), "elastic_fluctuation_nve", nvtExample.parent_path());
  ScopedCurrentPath currentPath(nvtExample.parent_path());

  EXPECT_THROW(
      {
        InputReader reader(temporary.path().filename().string());
        static_cast<void>(reader);
      },
      std::runtime_error);
}
