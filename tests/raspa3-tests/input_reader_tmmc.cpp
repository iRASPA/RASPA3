#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "input_reader_fixtures.hpp"
#include "molecule_fixtures.hpp"

import std;

import input_reader;
import int3;
import json;
import thermobarostat;

namespace
{

void setBETProbeProperties(nlohmann::json& component)
{
  component["CrossSection"] = 16.2;
  component["LiquidVolume"] = 57.7;
  component["SaturationPressure"] = 101325.0;
}

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

// The transition-matrix example run as 'ParallelTMMC' instead: the same framework, component and macrostate
// range, two windows, and nothing created (the walkers are grown into their windows).
nlohmann::json readParallelTMMCExample()
{
  nlohmann::json input = readTMMCExample();
  input["SimulationType"] = "ParallelTMMC";
  input["NumberOfWindows"] = 2;
  input["Components"][0].erase("CreateNuMberofmolecules");
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  return input;
}

// What reading the input says when it refuses it, or nothing when it does not.
std::string refusalOf(nlohmann::json input, std::string_view suffix, const std::filesystem::path& directory)
{
  TemporaryInput temporary(std::move(input), suffix, directory);
  try
  {
    InputReader reader(temporary.path().filename().string());
    static_cast<void>(reader);
  }
  catch (const std::runtime_error& error)
  {
    return error.what();
  }
  return {};
}

}  // namespace

TEST(INPUT_READER, rejects_component_name_with_directory_separator)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  workspace.write("top/methane.json", molecule_fixtures::kMethaneJson);
  nlohmann::json input = readTMMCExample();
  input["Components"][0]["Name"] = "top/methane";
  TemporaryInput temporary(std::move(input), "component_path", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  try
  {
    InputReader reader(temporary.path().filename().string());
    FAIL() << "expected InputReader to reject a component Name that is a path";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    EXPECT_NE(message.find("directory separator"), std::string::npos);
  }
}

TEST(INPUT_READER, rejects_framework_name_with_directory_separator)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["Systems"][0]["Name"] = "top/tobacco-667";
  TemporaryInput temporary(std::move(input), "framework_path", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  try
  {
    InputReader reader(temporary.path().filename().string());
    FAIL() << "expected InputReader to reject a framework Name that is a path";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    EXPECT_NE(message.find("directory separator"), std::string::npos);
  }
}

TEST(INPUT_READER, component_file_name_loads_path_and_derives_name_from_stem)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  workspace.write("top/methane.json", molecule_fixtures::kMethaneJson);
  nlohmann::json input = readTMMCExample();
  input["Components"][0].erase("Name");
  input["Components"][0]["FileName"] = "top/methane";
  TemporaryInput temporary(std::move(input), "component_filename", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_EQ(reader.systems[0].components.size(), 1uz);
  EXPECT_EQ(reader.systems[0].components[0].name, "methane");
}

TEST(INPUT_READER, component_name_overrides_stem_of_file_name)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  workspace.write("top/methane.json", molecule_fixtures::kMethaneJson);
  nlohmann::json input = readTMMCExample();
  input["Components"][0]["Name"] = "CH4-alias";
  input["Components"][0]["FileName"] = "top/methane.json";
  TemporaryInput temporary(std::move(input), "component_name_override", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_EQ(reader.systems[0].components.size(), 1uz);
  EXPECT_EQ(reader.systems[0].components[0].name, "CH4-alias");
}

TEST(INPUT_READER, framework_file_name_loads_path_and_derives_name_from_stem)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  workspace.write("top/tobacco-667.cif", input_reader_fixtures::kTobacco667Cif);
  nlohmann::json input = readTMMCExample();
  input["Systems"][0].erase("Name");
  input["Systems"][0]["FileName"] = "top/tobacco-667";
  TemporaryInput temporary(std::move(input), "framework_filename", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_TRUE(reader.systems[0].framework.has_value());
  EXPECT_EQ(reader.systems[0].framework->name, "tobacco-667");
}

TEST(INPUT_READER, rejects_duplicate_component_names)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["Components"].push_back(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "duplicate_names", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  try
  {
    InputReader reader(temporary.path().filename().string());
    FAIL() << "expected InputReader to reject duplicate component names";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    EXPECT_NE(message.find("duplicate component name"), std::string::npos);
  }
}

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

// 'ParallelTMMC' samples the same nearest-neighbour transition matrix, and the moves that do not fit it are
// refused for it as they are for 'MonteCarloTransitionMatrix'. It turns the transition matrix on after the
// moves are looked at, and they used to go through unlooked at.
TEST(INPUT_READER_TMMC, parallel_tmmc_rejects_the_moves_the_transition_matrix_cannot_hold)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_EQ(refusalOf(readParallelTMMCExample(), "ptmmc_plain", workspace.path()), "");

  nlohmann::json identity = readParallelTMMCExample();
  identity["Components"][0]["IdentityChangeProbability"] = 1.0;
  const std::string identityRefusal = refusalOf(std::move(identity), "ptmmc_identity", workspace.path());
  EXPECT_NE(identityRefusal.find("IdentityChangeProbability"), std::string::npos) << identityRefusal;
  EXPECT_NE(identityRefusal.find("must be zero"), std::string::npos) << identityRefusal;

  nlohmann::json tempering = readParallelTMMCExample();
  tempering["Systems"][0]["ParallelTemperingSwapProbability"] = 1.0;
  const std::string temperingRefusal = refusalOf(std::move(tempering), "ptmmc_tempering", workspace.path());
  EXPECT_NE(temperingRefusal.find("ParallelTemperingSwapProbability"), std::string::npos) << temperingRefusal;
  EXPECT_NE(temperingRefusal.find("must be zero"), std::string::npos) << temperingRefusal;
}

// Its walkers are grown into their windows, so what is created may lie below the macrostate range, which
// for a single transition-matrix run it may not.
TEST(INPUT_READER_TMMC, parallel_tmmc_may_create_fewer_molecules_than_the_range_starts_at)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json parallel = readParallelTMMCExample();
  parallel["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 5;
  EXPECT_EQ(refusalOf(std::move(parallel), "ptmmc_below_range", workspace.path()), "");

  nlohmann::json single = readTMMCExample();
  single["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 5;
  single["Components"][0].erase("CreateNuMberofmolecules");
  single["Components"][0]["CreateNumberOfMolecules"] = 0;
  EXPECT_NE(refusalOf(std::move(single), "tmmc_below_range", workspace.path()).find("macrostate"),
            std::string::npos);
}

TEST(INPUT_READER_TMMC, enables_lambda_chain_for_cfcmc_swap)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["NumberOfLambdaBins"] = 7;
  input["Components"][0]["SwapProbability"] = 0.0;
  input["Components"][0]["CFCMC_CBMC_SwapProbability"] = 1.0;
  TemporaryInput temporary(std::move(input), "lambda_chain", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_EQ(reader.systems.size(), 1uz);
  EXPECT_EQ(reader.systems[0].components[0].lambdaGC.numberOfSamplePoints, 7uz);
  EXPECT_EQ(reader.systems[0].tmmc.numberOfLambdaBins, 7uz);
}

// An integer swap beside the CFCMC swap on a lambda chain jumps (N, k) -> (N +/- 1, k), which the (N, lambda)
// matrix cannot hold; the mix is refused. With a single lambda bin the chain is plain N and the mix is fine.
TEST(INPUT_READER_TMMC, rejects_integer_swaps_beside_a_cfcmc_swap_on_a_lambda_chain)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  {
    nlohmann::json input = readTMMCExample();
    input["NumberOfLambdaBins"] = 7;
    input["Components"][0]["SwapProbability"] = 1.0;
    input["Components"][0]["CFCMC_CBMC_SwapProbability"] = 1.0;
    TemporaryInput temporary(std::move(input), "lambda_chain_mix", workspace.path());
    ScopedCurrentPath currentPath(workspace.path());

    EXPECT_THROW(
        {
          InputReader reader(temporary.path().filename().string());
          static_cast<void>(reader);
        },
        std::runtime_error);
  }
  {
    nlohmann::json input = readTMMCExample();
    input["NumberOfLambdaBins"] = 1;
    input["Components"][0]["SwapProbability"] = 1.0;
    input["Components"][0]["CFCMC_CBMC_SwapProbability"] = 1.0;
    TemporaryInput temporary(std::move(input), "plain_chain_mix", workspace.path());
    ScopedCurrentPath currentPath(workspace.path());

    InputReader reader(temporary.path().filename().string());
    ASSERT_EQ(reader.systems.size(), 1uz);
    EXPECT_FALSE(reader.systems[0].tmmc.lambdaChain());
  }
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

TEST(INPUT_READER, compute_bet_is_accepted_for_parallel_tmmc)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ComputeBET"] = true;
  input["NumberOfWindows"] = 2;
  input["NumberOfInitializationCycles"] = 10;
  input["NumberOfEquilibrationCycles"] = 10;
  input["NumberOfProductionCycles"] = 10;
  input["Systems"][0].erase("Ensemble");
  input["Systems"][0].erase("ComputeMSD");
  input["Systems"][0].erase("SampleMSDEvery");
  input["Systems"][0].erase("WriteMSDEvery");
  input["Systems"][0]["ExternalPressure"] = 1.0e5;
  input["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 0;
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = 10;
  input["Components"][0]["SwapProbability"] = 1.0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "compute_bet", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.computeBET);
  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::ParallelTMMC);
}

TEST(INPUT_READER, compute_bte_alias_sets_compute_bet)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ComputeBTE"] = true;
  input["NumberOfWindows"] = 2;
  input["Systems"][0].erase("Ensemble");
  input["Systems"][0].erase("ComputeMSD");
  input["Systems"][0].erase("SampleMSDEvery");
  input["Systems"][0].erase("WriteMSDEvery");
  input["Systems"][0]["ExternalPressure"] = 1.0e5;
  input["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 0;
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = 10;
  input["Components"][0]["SwapProbability"] = 1.0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "compute_bte", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.computeBET);
}

TEST(INPUT_READER, compute_bet_rejected_for_molecular_dynamics)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson);
  input["ComputeBET"] = true;
  TemporaryInput temporary(std::move(input), "compute_bet_md", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, compute_bet_accepts_auto_pressures_for_wham)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["ReweightingPressureRange"] = "auto";
  input["NumberOfThreads"] = 16;
  input["Systems"][0]["ExternalPressures"] = "auto";
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "wham_auto", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.computeBET);
  EXPECT_TRUE(reader.autoExternalPressures);
  EXPECT_TRUE(reader.autoReweightingPressureRange);
  EXPECT_TRUE(reader.parallelTemperingPressures.empty());
  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::ReweightedHistogram);
}

TEST(INPUT_READER, compute_bet_omitted_pressures_imply_auto_for_wham)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "wham_omit", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.autoExternalPressures);
  EXPECT_TRUE(reader.autoReweightingPressureRange);
  EXPECT_EQ(reader.parallelTemperingTemperatures.size(), 1uz);
}

TEST(INPUT_READER, auto_pressures_rejected_without_compute_bet)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ReweightingPressureRange"] = "auto";
  input["Systems"][0]["ExternalTemperatures"] = nlohmann::json::array({111.58, 120.0});
  input["Systems"][0]["ExternalPressures"] = "auto";
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  TemporaryInput temporary(std::move(input), "auto_no_bet", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, compute_bet_omitted_range_implies_auto_for_tmmc)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ComputeBET"] = true;
  input["NumberOfWindows"] = 2;
  input["Components"][0].erase("CreateNuMberofmolecules");
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "tmmc_omit_range", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.computeBET);
  EXPECT_TRUE(reader.autoReweightingPressureRange);
  EXPECT_FALSE(reader.reweightingPressureRange.has_value());
}

// With 'ComputeBET' on a framework a macrostate maximum that is left out is automatic, placed by a scout
// once the system exists. The checks of the range have to know that when they run: they used to run first,
// against the transition matrix's default maximum of a hundred, and refused more than a hundred windows, or
// a minimum of a hundred or more, in an input that had given no maximum to hold them to.
TEST(INPUT_READER, compute_bet_omitted_maximum_is_automatic_before_the_range_is_checked)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  auto withoutMaximum = []()
  {
    nlohmann::json input = readParallelTMMCExample();
    input["ComputeBET"] = true;
    input["Systems"][0].erase("MacroStateMaximumNumberOfMolecules");
    setBETProbeProperties(input["Components"][0]);
    return input;
  };

  {
    nlohmann::json input = withoutMaximum();
    input["NumberOfWindows"] = 150;
    TemporaryInput temporary(std::move(input), "ptmmc_many_windows", workspace.path());
    InputReader reader(temporary.path().filename().string());
    EXPECT_TRUE(reader.autoMacroStateMaximum);
    EXPECT_EQ(reader.tmmcNumberOfWindows, 150uz);
  }
  {
    nlohmann::json input = withoutMaximum();
    input["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 120;
    TemporaryInput temporary(std::move(input), "ptmmc_high_minimum", workspace.path());
    InputReader reader(temporary.path().filename().string());
    EXPECT_TRUE(reader.autoMacroStateMaximum);
  }

  // A maximum that is given is still held to: more windows than it has macrostates for are refused.
  nlohmann::json given = readParallelTMMCExample();
  given["ComputeBET"] = true;
  setBETProbeProperties(given["Components"][0]);
  given["NumberOfWindows"] = 150;
  EXPECT_NE(refusalOf(std::move(given), "ptmmc_given_maximum", workspace.path()).find("NumberOfWindows"),
            std::string::npos);
}

TEST(INPUT_READER, auto_reweighting_range_rejected_without_compute_bet)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ReweightingPressureRange"] = "auto";
  input["NumberOfWindows"] = 2;
  input["Systems"][0].erase("Ensemble");
  input["Systems"][0].erase("ComputeMSD");
  input["Systems"][0].erase("SampleMSDEvery");
  input["Systems"][0].erase("WriteMSDEvery");
  input["Systems"][0]["ExternalPressure"] = 1.0e5;
  input["Systems"][0]["MacroStateMinimumNumberOfMolecules"] = 0;
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = 10;
  input["Components"][0]["SwapProbability"] = 1.0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  TemporaryInput temporary(std::move(input), "tmmc_auto_no_bet", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, compute_bet_accepts_auto_macrostate_for_wham)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = "auto";
  TemporaryInput temporary(std::move(input), "wham_auto_nmax", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.autoMacroStateMaximum);
}

TEST(INPUT_READER, wham_tolerance_defaults_and_accepts_override)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporaryDefault(input, "wham_tol_default", workspace.path());
  {
    ScopedCurrentPath currentPath(workspace.path());
    InputReader reader(temporaryDefault.path().filename().string());
    EXPECT_DOUBLE_EQ(reader.whamTolerance, 1.0e-6);
    EXPECT_EQ(reader.whamMaximumIterations, 1000000uz);
    EXPECT_EQ(reader.betScoutMaximumCycles, 15000uz);
  }

  input["WHAMTolerance"] = 1.0e-5;
  input["WHAMIterations"] = 500000;
  input["BETScoutMaximumCycles"] = 30000;
  TemporaryInput temporaryOverride(std::move(input), "wham_tol_override", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());
  InputReader reader(temporaryOverride.path().filename().string());
  EXPECT_DOUBLE_EQ(reader.whamTolerance, 1.0e-5);
  EXPECT_EQ(reader.whamMaximumIterations, 500000uz);
  EXPECT_EQ(reader.betScoutMaximumCycles, 30000uz);
}

TEST(INPUT_READER, wham_tolerance_rejects_non_positive)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["WHAMTolerance"] = 0.0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "wham_tol_bad", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, wham_iterations_rejects_zero)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["WHAMIterations"] = 0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "wham_iters_bad", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, bet_scout_maximum_cycles_rejects_zero)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ReweightedHistogram";
  input["ComputeBET"] = true;
  input["BETScoutMaximumCycles"] = 0;
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  TemporaryInput temporary(std::move(input), "bet_scout_bad", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, compute_bet_accepts_auto_macrostate_for_tmmc)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ComputeBET"] = true;
  input["NumberOfWindows"] = 2;
  input["Components"][0].erase("CreateNuMberofmolecules");
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = "auto";
  TemporaryInput temporary(std::move(input), "tmmc_auto_nmax", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.autoMacroStateMaximum);
}

TEST(INPUT_READER, compute_bet_omitted_macrostate_implies_auto)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["ComputeBET"] = true;
  input["NumberOfWindows"] = 2;
  input["Components"][0].erase("CreateNuMberofmolecules");
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  setBETProbeProperties(input["Components"][0]);
  input["Systems"][0].erase("MacroStateMaximumNumberOfMolecules");
  TemporaryInput temporary(std::move(input), "tmmc_omit_nmax", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  EXPECT_TRUE(reader.autoMacroStateMaximum);
}

TEST(INPUT_READER, auto_macrostate_rejected_without_compute_bet)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kTmmcSimulationJson);
  input["SimulationType"] = "ParallelTMMC";
  input["NumberOfWindows"] = 2;
  input["Components"][0].erase("CreateNuMberofmolecules");
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  input["Systems"][0]["MacroStateMaximumNumberOfMolecules"] = "auto";
  TemporaryInput temporary(std::move(input), "auto_nmax_no_bet", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

TEST(INPUT_READER, number_of_unit_cells_auto_satisfies_minimum_image_convention)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  // Unit cell ~30.6 Å; CutOff 20 Å → need ceil(40 / 30.6) = 2 replicas per axis.
  // The cut-off is a force-field option, so it is set in force_field.json.
  nlohmann::json forceField = nlohmann::json::parse(input_reader_fixtures::kTmmcForceFieldJson);
  forceField["CutOffVDW"] = 20.0;
  workspace.write("force_field.json", forceField.dump(2));
  nlohmann::json input = readTMMCExample();
  input["Systems"][0]["NumberOfUnitCells"] = "auto";
  TemporaryInput temporary(std::move(input), "unit_cells_auto", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_TRUE(reader.systems[0].framework.has_value());
  const int3 cells = reader.systems[0].framework->numberOfUnitCells;
  EXPECT_EQ(cells.x, 2);
  EXPECT_EQ(cells.y, 2);
  EXPECT_EQ(cells.z, 2);
}

TEST(INPUT_READER, number_of_unit_cells_rejects_invalid_string)
{
  TemporaryDirectory workspace = makeTmmcWorkspace();
  nlohmann::json input = readTMMCExample();
  input["Systems"][0]["NumberOfUnitCells"] = "automatic";
  TemporaryInput temporary(std::move(input), "unit_cells_bad", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  EXPECT_THROW(InputReader reader(temporary.path().filename().string()), std::runtime_error);
}

namespace
{
std::string dcTSTInputError(bool crossLinked)
{
  TemporaryDirectory workspace;
  workspace.write("force_field.json", input_reader_fixtures::kBoxForceFieldJson);
  nlohmann::json molecule = nlohmann::json::parse(molecule_fixtures::kMethaneJson);
  molecule["ReactiveSites"] = nlohmann::json::array({nlohmann::json::array({0, "X"})});
  workspace.write("methane.json", molecule.dump());

  nlohmann::json input = readNVTExample();
  input["SimulationType"] = "DynamicallyCorrectedTST";
  nlohmann::json& system = input["Systems"][0];
  system.erase("Ensemble");
  system.erase("ComputeMSD");
  system.erase("SampleMSDEvery");
  system.erase("WriteMSDEvery");
  system["dcTST"] = {{"TaggedComponent", "methane"},
                     {"BarrierPosition", {0.0, 0.5, 0.5}},
                     {"BarrierNormal", {1.0, 0.0, 0.0}},
                     {"LatticeDistance", 30.0}};
  input["Components"][0]["CreateNumberOfMolecules"] = 4;
  if (crossLinked)
  {
    system["CrossLinkBonds"] = nlohmann::json::array(
        {{{"Sites", {"X", "X"}}, {"Bond", nlohmann::json::array({"HARMONIC", {100.0, 3.8}})}}});
    system["InitialCrossLinks"] = nlohmann::json::array({{{0, 0, 0}, {0, 1, 0}}});
  }
  TemporaryInput temporary(std::move(input), "dctst_cross_links", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());
  try
  {
    InputReader reader(temporary.path().filename().string());
  }
  catch (const std::exception& error)
  {
    return error.what();
  }
  return {};
}
}  // namespace

// The tagged-molecule moves of dcTST leave the link energy out, while the trajectories include the link forces.
TEST(INPUT_READER, dctst_rejects_cross_links)
{
  EXPECT_EQ(dcTSTInputError(false).find("cross-link"), std::string::npos) << dcTSTInputError(false);
  EXPECT_NE(dcTSTInputError(true).find("cross-link"), std::string::npos) << dcTSTInputError(true);
}

namespace
{
constexpr std::string_view kCrossLinkMonomerJson = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["CH4", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [{"Atom": 0, "Type": "X", "Valence": 2}]
}
)";

constexpr std::string_view kCrossLinkMonomerWithBendJson = R"({
  "CriticalTemperature" : 190.6,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.011,
  "pseudoAtoms" : [["CH4", [0.0, 0.0, 0.0]]],
  "ReactiveSites" : [{"Atom": 0, "Type": "X", "Valence": 2, "LinkBend": ["HARMONIC", [0.0, 120.0]]}]
}
)";

nlohmann::json crossLinkInput()
{
  nlohmann::json input = readNVTExample();
  input["SimulationType"] = "MonteCarlo";
  input["Systems"][0]["CrossLinkBonds"] = {{{"Sites", {"X", "X"}}, {"Bond", {"HARMONIC", {1000.0, 1.5}}}}};
  input["Components"][0]["Name"] = "monomer";
  input["Components"][0]["CreateNumberOfMolecules"] = 4;
  input["Systems"][0]["CrossLinkFormationProbability"] = 0.5;
  return input;
}
}  // namespace

// A reactive site that can carry two links forms the bend between its partners: the bend has to be given.
TEST(INPUT_READER, cross_link_site_with_valence_above_one_needs_a_link_bend)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  workspace.write("monomer.json", kCrossLinkMonomerJson);
  TemporaryInput temporary(crossLinkInput(), "cross_link_valence", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  try
  {
    InputReader reader(temporary.path().filename().string());
    FAIL() << "expected InputReader to reject a valence-2 site without a 'LinkBend'";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    EXPECT_NE(message.find("LinkBend"), std::string::npos) << message;
    EXPECT_NE(message.find("valence 2"), std::string::npos) << message;
  }
}

TEST(INPUT_READER, cross_link_site_with_valence_above_one_and_a_link_bend_is_accepted)
{
  TemporaryDirectory workspace = makeBoxMdWorkspace();
  workspace.write("monomer.json", kCrossLinkMonomerWithBendJson);
  TemporaryInput temporary(crossLinkInput(), "cross_link_valence_bend", workspace.path());
  ScopedCurrentPath currentPath(workspace.path());

  InputReader reader(temporary.path().filename().string());
  ASSERT_EQ(reader.systems.size(), 1uz);
  ASSERT_EQ(reader.systems[0].components.size(), 1uz);
  ASSERT_EQ(reader.systems[0].components[0].reactiveSites.size(), 1uz);
  EXPECT_EQ(reader.systems[0].components[0].reactiveSites[0].valence, 2uz);
  EXPECT_TRUE(reader.systems[0].components[0].reactiveSites[0].linkBend.has_value());
}

// The electric field of the induced-dipole model does not exclude the pairs a bond between molecules excludes.
TEST(INPUT_READER, cross_links_and_polymerization_reject_polarization)
{
  for (const bool polymerization : {false, true})
  {
    TemporaryDirectory workspace = makeBoxMdWorkspace();
    workspace.write("monomer.json", kCrossLinkMonomerWithBendJson);
    nlohmann::json input = crossLinkInput();
    input["Systems"][0]["ComputePolarization"] = true;
    if (polymerization)
    {
      input["Systems"][0].erase("CrossLinkBonds");
      input["Systems"][0].erase("CrossLinkFormationProbability");
      input["Systems"][0]["Polymerization"] = {{"StartCycle", 0}};
    }
    TemporaryInput temporary(std::move(input), polymerization ? "polymerization_polarization" : "cross_link_polarization",
                             workspace.path());
    ScopedCurrentPath currentPath(workspace.path());

    try
    {
      InputReader reader(temporary.path().filename().string());
      FAIL() << "expected InputReader to reject polarization with bonds between molecules (polymerization: "
             << polymerization << ")";
    }
    catch (const std::runtime_error& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("polarization is not supported for bonds between molecules"), std::string::npos)
          << message;
      EXPECT_NE(message.find(polymerization ? "'Polymerization'" : "'CrossLinkBonds'"), std::string::npos) << message;
    }
  }
}

TEST(INPUT_READER, rejects_an_unknown_or_non_string_use_charges_from)
{
  // a 'UseChargesFrom' that matched none of the three values, or was not a string, fell back to the pseudo-atom
  // charges without a word, which in a charge-equilibration set-up are all zero
  TemporaryDirectory workspace = makeTmmcWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json misspelt = readTMMCExample();
  misspelt["Systems"][0]["UseChargesFrom"] = "Charge_Equilibration";
  const std::string misspeltMessage = refusalOf(std::move(misspelt), "charges_from_misspelt", workspace.path());
  EXPECT_NE(misspeltMessage.find("unknown 'UseChargesFrom' 'Charge_Equilibration'"), std::string::npos)
      << misspeltMessage;

  nlohmann::json number = readTMMCExample();
  number["Systems"][0]["UseChargesFrom"] = 2;
  const std::string numberMessage = refusalOf(std::move(number), "charges_from_number", workspace.path());
  EXPECT_NE(numberMessage.find("'UseChargesFrom' must be a string"), std::string::npos) << numberMessage;

  // the three values, in any case, and the misspelling that was accepted for a while
  for (const std::string_view accepted : {"CIF_File", "cif_file", "PseudoAtoms", "PsuedoAtoms"})
  {
    nlohmann::json input = readTMMCExample();
    input["Systems"][0]["UseChargesFrom"] = accepted;
    EXPECT_EQ(refusalOf(std::move(input), "charges_from_accepted", workspace.path()), "") << accepted;
  }
}
