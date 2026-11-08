#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import archive;
import characterset;
import getopt;
import indexpath;
import scanner;
import stringutils;
import threadpool;
import units;
import json;

TEST(FOUNDATIONKIT_FIXES, scanner_reads_numbers_and_the_last_token_and_line)
{
  Scanner numbers("1.5 -7 2.25", CharacterSet::whitespaceAndNewlineCharacterSet());
  double x{};
  int n{};
  EXPECT_TRUE(numbers.scanDouble(x));
  EXPECT_DOUBLE_EQ(x, 1.5);
  EXPECT_TRUE(numbers.scanInt(n));
  EXPECT_EQ(n, -7);
  // the last token has no delimiter after it
  EXPECT_TRUE(numbers.scanDouble(x));
  EXPECT_DOUBLE_EQ(x, 2.25);
  EXPECT_FALSE(numbers.scanDouble(x));
  EXPECT_FALSE(Scanner("abc", CharacterSet::whitespaceAndNewlineCharacterSet()).scanDouble(x));

  Scanner tokens("alpha\tbeta", CharacterSet::whitespaceAndNewlineCharacterSet());
  std::string token;
  EXPECT_TRUE(tokens.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), token));
  EXPECT_EQ(token, "alpha");  // the tab is whitespace
  EXPECT_TRUE(tokens.scanUpToCharacters(CharacterSet::whitespaceAndNewlineCharacterSet(), token));
  EXPECT_EQ(token, "beta");
  EXPECT_TRUE(tokens.isAtEnd());

  Scanner lines("first\nlast", CharacterSet::whitespaceCharacterSet());
  std::string line;
  EXPECT_TRUE(lines.scanLine(line));
  EXPECT_EQ(line, "first");
  EXPECT_TRUE(lines.scanLine(line));
  EXPECT_EQ(line, "last");
  EXPECT_FALSE(lines.scanLine(line));
}

TEST(FOUNDATIONKIT_FIXES, whitespace_and_newline_set_holds_the_tab)
{
  const std::string set = CharacterSet::whitespaceAndNewlineCharacterSet().string();
  EXPECT_TRUE(set.contains('\t'));
  EXPECT_TRUE(set.contains(' '));
  EXPECT_TRUE(set.contains('\n'));
}

TEST(FOUNDATIONKIT_FIXES, equal_index_paths_compare_equal)
{
  EXPECT_TRUE(IndexPath(3) == IndexPath(3));
  EXPECT_FALSE(IndexPath(3) == IndexPath(4));
}

TEST(FOUNDATIONKIT_FIXES, an_optional_argument_takes_the_next_value_but_not_the_next_option)
{
  std::array<std::string, 5> storage{"raspa3", "-o", "value", "-o", "-v"};
  std::array<char*, 5> argv{};
  for (std::size_t i = 0; i < storage.size(); ++i) argv[i] = storage[i].data();

  std::vector<std::string> optionalValues;
  bool verbose = false;
  argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
  parser.reg({"-o"}, "value", argparser::argparser::optional_argument, "",
             [&](const std::string& value) { optionalValues.push_back(value); });
  parser.reg({"-v"}, argparser::argparser::no_argument, "", [&](const std::string&) { verbose = true; });
  parser();

  ASSERT_EQ(optionalValues.size(), 2uz);
  EXPECT_EQ(optionalValues[0], "value");
  EXPECT_EQ(optionalValues[1], "");
  EXPECT_TRUE(verbose);
}

TEST(FOUNDATIONKIT_FIXES, trim_handles_empty_and_blank_strings)
{
  EXPECT_EQ(trim(""), "");
  EXPECT_EQ(trim("   "), "");
  EXPECT_EQ(trim("  a b \t"), "a b");
}

TEST(FOUNDATIONKIT_FIXES, reading_a_map_replaces_its_contents)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "map.bin";
  {
    std::ofstream stream(file, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << std::map<int, double>{{1, 1.5}};
    archive << std::unordered_map<int, double>{{2, 2.5}};
  }
  std::map<int, double> map{{7, 7.0}};
  std::unordered_map<int, double> unorderedMap{{8, 8.0}};
  {
    std::ifstream stream(file, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> map;
    archive >> unorderedMap;
  }
  EXPECT_EQ(map, (std::map<int, double>{{1, 1.5}}));
  EXPECT_EQ(unorderedMap, (std::unordered_map<int, double>{{2, 2.5}}));
}

TEST(FOUNDATIONKIT_FIXES, acceleration_unit_is_length_per_time_squared)
{
  // the acceleration factor and strings carried the diffusion unit m^2/s
  Units::setUnits(Units::System::RASPA);
  EXPECT_DOUBLE_EQ(Units::AccelerationConversionFactor, Units::LengthUnit / (Units::TimeUnit * Units::TimeUnit));
  EXPECT_DOUBLE_EQ(Units::DiffusionConversionFactor, Units::LengthUnit * Units::LengthUnit / Units::TimeUnit);
  EXPECT_EQ(Units::unitOfAccelerationString, "m.s⁻²");
  EXPECT_EQ(Units::unitOfDiffusionString, "m².s⁻¹");

  const nlohmann::json status = Units::jsonStatus();
  EXPECT_TRUE(status.contains("acceleration [m/s²]"));
  EXPECT_FALSE(status.contains("acceleration [m²/s]"));
  EXPECT_DOUBLE_EQ(status["acceleration [m/s²]"].get<double>(), Units::AccelerationConversionFactor);
}

TEST(FOUNDATIONKIT_FIXES, reduced_units_set_the_diffusion_and_energy_per_mol_strings)
{
  Units::setUnits(Units::System::ReducedUnits);
  EXPECT_EQ(Units::unitOfDiffusionString, "σ².𝜏⁻¹");
  EXPECT_EQ(Units::unitOfAccelerationString, "ε.σ⁻¹.m⁻¹");
  // used to keep the SI string of the previous unit system
  EXPECT_NE(Units::unitOfEnergyPerMolString, "J.mol⁻¹");
  EXPECT_EQ(Units::unitOfEnergyPerMolString, Units::unitOfEnthalpyPerMolString);

  Units::setUnits(Units::System::RASPA);
  EXPECT_EQ(Units::unitOfEnergyPerMolString, "J.mol⁻¹");
  EXPECT_EQ(Units::unitOfAccelerationString, "m.s⁻²");
}

TEST(FOUNDATIONKIT_FIXES, serial_thread_pool_runs_enqueued_tasks_inline)
{
  // one thread requested: no worker; the task used to be dropped and the future threw broken_promise
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  pool.init(1, ThreadPool::ThreadingType::ThreadPool);
  EXPECT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::Serial);
  EXPECT_EQ(pool.getThreadCount(), 0uz);

  std::future<int> future = pool.enqueue([](int a, int b) { return a * b; }, 6, 7);
  EXPECT_EQ(future.get(), 42);

  std::atomic<int> counter{0};
  pool.enqueue_detach([&counter] { counter.fetch_add(1); });
  pool.wait_for_tasks();
  EXPECT_EQ(counter.load(), 1);

  std::future<void> failing = pool.enqueue([] { throw std::runtime_error("inside the task"); });
  EXPECT_THROW(failing.get(), std::runtime_error);
}

TEST(FOUNDATIONKIT_FIXES, thread_pool_completes_tasks_queued_faster_than_its_workers_wake)
{
  // many tasks land on a worker before it wakes: the signal used to be a binary semaphore released past its
  // single permit
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  pool.init(3, ThreadPool::ThreadingType::ThreadPool);
  ASSERT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::ThreadPool);
  ASSERT_EQ(pool.getThreadCount(), 2uz);

  std::vector<std::future<std::size_t>> futures;
  for (std::size_t i = 0; i < 2000; ++i) futures.push_back(pool.enqueue([i] { return i * i; }));
  pool.wait_for_tasks();
  for (std::size_t i = 0; i < futures.size(); ++i) EXPECT_EQ(futures[i].get(), i * i);

  // asking for one thread afterwards switches to serial work, even though the workers stay alive
  pool.init(1, ThreadPool::ThreadingType::ThreadPool);
  EXPECT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::Serial);
  EXPECT_EQ(pool.getThreadCount(), 0uz);
  EXPECT_EQ(pool.enqueue([] { return 5; }).get(), 5);

  pool.init(3, ThreadPool::ThreadingType::ThreadPool);
  EXPECT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::ThreadPool);
  EXPECT_EQ(pool.getThreadCount(), 2uz);
}

TEST(FOUNDATIONKIT_FIXES, archive_chains_optional_and_variant_and_rejects_a_bad_variant_index)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "chain.bin";
  {
    // the optional and variant insertions returned a const archive, which ended the chain
    std::ofstream stream(file, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << std::optional<int>{3} << std::variant<int, double>{2.5} << std::optional<int>{} << 7;
    archive << std::size_t{2};  // not a valid variant index
  }
  std::ifstream stream(file, std::ios::binary);
  Archive<std::ifstream> archive(stream);
  std::optional<int> first, third;
  std::variant<int, double> second;
  int fourth = 0;
  archive >> first >> second >> third >> fourth;
  EXPECT_EQ(first, std::optional<int>{3});
  EXPECT_EQ(second, (std::variant<int, double>{2.5}));
  EXPECT_FALSE(third.has_value());
  EXPECT_EQ(fourth, 7);
  std::variant<int, double> bad{1};
  EXPECT_THROW(archive >> bad, std::runtime_error);
}

TEST(FOUNDATIONKIT_FIXES, a_restart_file_shorter_than_its_header_is_reported_as_truncated)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "restart_data.bin";
  {
    std::ofstream stream(file, std::ios::binary);
    const std::uint32_t marker = 0x01020304;
    stream.write(reinterpret_cast<const char*>(&marker), sizeof(marker));
    stream.write("abcdef", 6);  // 10 bytes: past the marker but short of the 20-byte header
  }
  int object = 0;
  try
  {
    readBinaryRestartFile(object, file.string());
    FAIL() << "a 10-byte restart file was accepted";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    EXPECT_NE(message.find(file.string()), std::string::npos) << message;
    EXPECT_NE(message.find("truncated"), std::string::npos) << message;
  }
}
