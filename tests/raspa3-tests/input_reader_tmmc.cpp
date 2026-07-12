#include <gtest/gtest.h>

import std;

import input_reader;
import json;

namespace
{

const std::filesystem::path tmmcExample = "examples/advanced/5_tmmc_methane_in_tobacco_667/0/simulation.json";

class TemporaryInput
{
 public:
  TemporaryInput(nlohmann::json input, std::string_view suffix)
      : path_(tmmcExample.parent_path() / std::format("simulation_tmmc_{}_{}.json", suffix,
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
