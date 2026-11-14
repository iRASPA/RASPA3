#include <H5Cpp.h>
#include <gtest/gtest.h>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

#include <signal.h>
#include <sys/resource.h>

#include "../test_support.hpp"

import std;

import archive;
import characterset;
import getopt;
import indexpath;
import scanner;
import stringutils;
import threadpool;
import hdf5;
import opencl;
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

TEST(FOUNDATIONKIT_FIXES, hdf5_writer_refuses_short_buffers_and_writes_bools_as_hbool)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "data.h5";
  {
    HDF5Writer writer(file.string());
    writer.createGroup("/g");
    writer.createDataset<double>("/g", "d", {3}, {});
    writer.createDataset<bool>("/g", "b", {3}, {});
    writer.createDataset<int>("/g", "m", {2, 2}, {});

    // a whole-dataset write reads as many elements as the dataset holds: a shorter vector was read past its end
    EXPECT_THROW(writer.writeVector<double>("/g", "d", std::vector<double>{1.0, 2.0}), std::invalid_argument);
    EXPECT_THROW(writer.writeVector<bool>("/g", "b", std::vector<bool>{true}), std::invalid_argument);
    writer.writeVector<double>("/g", "d", std::vector<double>{1.0, 2.0, 3.0});
    writer.writeVector<bool>("/g", "b", std::vector<bool>{true, false, true});

    // one coordinate per dimension: fewer indices were read past the end of the index vector
    EXPECT_THROW(writer.writeSingleValue<int>("/g", "m", std::vector<hsize_t>{1}, 5), std::invalid_argument);
    EXPECT_THROW(writer.writeSingleValue<int>("/g", "m", std::vector<hsize_t>{}, 5), std::invalid_argument);
    writer.writeSingleValue<int>("/g", "m", std::vector<hsize_t>{1, 0}, 5);
  }

  H5::H5File reader(file.string(), H5F_ACC_RDONLY);
  std::vector<double> d(3);
  reader.openDataSet("/g/d").read(d.data(), H5::PredType::NATIVE_DOUBLE);
  EXPECT_EQ(d, (std::vector<double>{1.0, 2.0, 3.0}));
  hbool_t b[3] = {false, false, false};
  reader.openDataSet("/g/b").read(b, H5::PredType::NATIVE_HBOOL);
  EXPECT_TRUE(b[0]);
  EXPECT_FALSE(b[1]);
  EXPECT_TRUE(b[2]);
  std::vector<int> m(4);
  reader.openDataSet("/g/m").read(m.data(), H5::PredType::NATIVE_INT);
  EXPECT_EQ(m, (std::vector<int>{0, 0, 5, 0}));
}

TEST(FOUNDATIONKIT_FIXES, opencl_device_queries_tolerate_absent_devices)
{
  // the counts start at zero and the errors are checked, so a device type that no platform offers (or no
  // platform at all) yields no device instead of sizing a vector from an uninitialized count
  EXPECT_FALSE(OpenCL::bestOpenCLDevice(CL_DEVICE_TYPE_ACCELERATOR).has_value());

  // the status report prints zero-initialized work-item and image sizes when a query is not answered
  OpenCL::initialize();
  const std::string status = OpenCL::printBestOpenCLDevice();
  EXPECT_FALSE(status.empty());
}

TEST(FOUNDATIONKIT_FIXES, simplified_trims_and_collapses_whitespace)
{
  // it used to return its argument unchanged
  EXPECT_EQ(simplified("  P 21/c   "), "P 21/c");
  EXPECT_EQ(simplified("\tF\n m  -3   m\t"), "F m -3 m");
  EXPECT_EQ(simplified("   "), "");
  EXPECT_EQ(simplified(""), "");
  EXPECT_EQ(simplified("single"), "single");
}

TEST(FOUNDATIONKIT_FIXES, case_insensitive_compare_accepts_bytes_above_ascii)
{
  // bytes above 0x7f are negative as char and were passed to std::tolower as such
  const std::string angstrom = "\xC3\x85ngstr\xC3\xB6m";
  EXPECT_TRUE(caseInSensStringCompare(angstrom, angstrom));
  EXPECT_TRUE(caseInSensStringCompare("\xC3\x85NGSTR\xC3\xB6M", angstrom));
  EXPECT_FALSE(caseInSensStringCompare("\xC3\x85ngstr\xC3\xB6m", "\xC3\x85ngstr\xC3\xB6n"));
  EXPECT_TRUE(caseInSensStringCompare("Methane", "mEtHaNe"));
  EXPECT_FALSE(caseInSensStringCompare("Methane", "Methan"));
}

TEST(FOUNDATIONKIT_FIXES, read_file_content_rejects_a_directory)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "content.txt";
  {
    std::ofstream stream(file);
    stream << "hello\n";
  }
  EXPECT_EQ(readFileContent(file.string(), ""), "hello\n");
  EXPECT_EQ(readFileContent((workspace.path() / "content").string(), ".txt"), "hello\n");
  // a directory exists, and reading it as a stream used to return an empty string silently
  EXPECT_THROW(readFileContent(workspace.path().string(), ""), std::runtime_error);
  EXPECT_THROW(readFileContent((workspace.path() / "missing").string(), ".txt"), std::runtime_error);
}

TEST(FOUNDATIONKIT_FIXES, index_path_decrement_stops_at_zero_and_const_paths_compose)
{
  IndexPath path(0);
  path.decreaseValueAtLastIndex();  // used to wrap to the largest std::size_t
  EXPECT_EQ(path.lastIndex(), 0uz);
  IndexPath two(2);
  two.decreaseValueAtLastIndex();
  EXPECT_EQ(two.lastIndex(), 1uz);

  // size, appending and + are usable on a const path
  const IndexPath base = IndexPath(1).appending(2);
  EXPECT_EQ(base.size(), 2uz);
  EXPECT_EQ(base.appending(3).size(), 3uz);
  EXPECT_EQ((base + IndexPath(7)).lastIndex(), 7uz);
}

TEST(FOUNDATIONKIT_FIXES, argparser_tolerates_a_negative_count_and_shows_help_from_the_parsed_instance)
{
  // a negative argc wrapped into a huge reservation
  argparser::argparser empty(-1, nullptr);
  EXPECT_TRUE(empty());

  std::array<std::string, 2> storage{"raspa3", "--help"};
  std::array<char*, 2> argv{};
  for (std::size_t i = 0; i < storage.size(); ++i) argv[i] = storage[i].data();
  argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
  parser.info("test", "raspa3").help({"-h", "--help"}, "show this help");
  // the help option no longer binds the address of the registering instance: a copy works as well
  argparser::argparser copy = parser;
  testing::internal::CaptureStdout();
  EXPECT_THROW(copy(), argparser::help_requested_exception);
  const std::string shown = testing::internal::GetCapturedStdout();
  EXPECT_NE(shown.find("Usage:"), std::string::npos);
}

// A restart checkpoint that runs out of space mid-payload (simulated with a file-size limit, which fails the write
// with EFBIG like ENOSPC or a quota) is dropped: the stream buffer swap after the payload cleared the failure, so
// the truncated payload was sealed with a valid checksum and renamed over the latest checkpoint, and the next
// checkpoint rotated it over the previous one as well.
TEST(FOUNDATIONKIT_FIXES, restart_checkpoint_that_runs_out_of_space_keeps_both_previous_checkpoints)
{
  TemporaryDirectory directory;
  const std::string fileName = (directory.path() / "restart_data.bin").string();
  auto payload = [](double value) { return std::vector<double>(1uz << 18uz, value); };  // 2 MiB
  auto contents = [](const std::string& name)
  {
    std::ifstream stream(name, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(stream), {});
  };

  writeBinaryRestartFile(payload(1.0), fileName);
  writeBinaryRestartFile(payload(2.0), fileName);
  const std::string latest = contents(fileName);
  const std::string previous = contents(fileName + ".prev");
  ASSERT_GT(latest.size(), 1uz << 21uz);
  ASSERT_GT(previous.size(), 1uz << 21uz);

  {
    rlimit original{};
    ASSERT_EQ(::getrlimit(RLIMIT_FSIZE, &original), 0);
    rlimit limited = original;
    limited.rlim_cur = 1u << 20u;
    // the limit signals SIGXFSZ on the write past it; ignored, the write fails with EFBIG instead
    auto* const originalHandler = ::signal(SIGXFSZ, SIG_IGN);
    ASSERT_EQ(::setrlimit(RLIMIT_FSIZE, &limited), 0);
    writeBinaryRestartFile(payload(3.0), fileName);
    ::setrlimit(RLIMIT_FSIZE, &original);
    ::signal(SIGXFSZ, originalHandler);
  }

  EXPECT_TRUE(contents(fileName) == latest) << "the latest checkpoint was replaced";
  EXPECT_TRUE(contents(fileName + ".prev") == previous) << "the previous checkpoint was replaced";
  EXPECT_FALSE(std::filesystem::exists(fileName + "_temp"));

  std::vector<double> restored;
  readBinaryRestartFile(restored, fileName);
  EXPECT_TRUE(restored == payload(2.0));
}

// A pool switched to serial after its workers started keeps them alive, but runs every task on the calling thread:
// the dispatch used to hand tasks to the workers whatever the threading type.
TEST(FOUNDATIONKIT_FIXES, a_pool_switched_to_serial_runs_tasks_on_the_calling_thread)
{
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  pool.init(3, ThreadPool::ThreadingType::ThreadPool);
  ASSERT_EQ(pool.getThreadCount(), 2uz);

  pool.init(3, ThreadPool::ThreadingType::Serial);
  ASSERT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::Serial);
  const std::thread::id caller = std::this_thread::get_id();
  for (std::size_t i = 0; i < 64; ++i)
  {
    EXPECT_EQ(pool.enqueue([] { return std::this_thread::get_id(); }).get(), caller);
  }

  pool.init(3, ThreadPool::ThreadingType::ThreadPool);
  EXPECT_EQ(pool.getThreadingType(), ThreadPool::ThreadingType::ThreadPool);
}

// No character of the set at the scan location is no match, as for the other scans (it was true with an empty
// string, so a loop over tokens never ended).
TEST(FOUNDATIONKIT_FIXES, scanning_characters_without_a_match_is_false)
{
  Scanner scanner("abc def", CharacterSet::newlineCharacterSet());
  std::string into = "unchanged";
  EXPECT_FALSE(scanner.scanCharacters(CharacterSet::whitespaceCharacterSet(), into));
  EXPECT_TRUE(into.empty());
  EXPECT_TRUE(scanner.scanUpToCharacters(CharacterSet::whitespaceCharacterSet(), into));
  EXPECT_EQ(into, "abc");
  EXPECT_TRUE(scanner.scanCharacters(CharacterSet::whitespaceCharacterSet(), into));
  EXPECT_EQ(into, " ");
}

// A single bool is written through an hbool_t, the size of the memory type, and reads back as written.
TEST(FOUNDATIONKIT_FIXES, hdf5_writer_writes_a_single_bool_as_hbool)
{
  TemporaryDirectory workspace;
  const std::filesystem::path file = workspace.path() / "single.h5";
  {
    HDF5Writer writer(file.string());
    writer.createGroup("/g");
    writer.createDataset<bool>("/g", "b", {3}, {});
    writer.writeVector<bool>("/g", "b", std::vector<bool>{false, false, false});
    writer.writeSingleValue<bool>("/g", "b", std::vector<hsize_t>{1}, true);
  }
  H5::H5File reader(file.string(), H5F_ACC_RDONLY);
  hbool_t b[3] = {true, false, true};
  reader.openDataSet("/g/b").read(b, H5::PredType::NATIVE_HBOOL);
  EXPECT_FALSE(b[0]);
  EXPECT_TRUE(b[1]);
  EXPECT_FALSE(b[2]);
}

// A misspelled switch is not a positional argument: '--orientation' (for '--orientations') was taken for the input
// file name when the positional slot was still free, and the switch meant was silently missing. A negative number or
// a plain file name is still a positional argument.
TEST(FOUNDATIONKIT_FIXES, argparser_refuses_an_unknown_switch_instead_of_taking_it_for_a_file)
{
  for (const std::string unknown : {"--orientation", "-x"})
  {
    std::array<std::string, 2> storage{"raspa3", unknown};
    std::array<char*, 2> argv{};
    for (std::size_t i = 0; i < storage.size(); ++i) argv[i] = storage[i].data();
    std::string file;
    argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
    parser.pos("FILE", "the input", [&](const std::string& value) { file = value; });
    EXPECT_THROW(parser(), argparser::unknown_option_exception) << unknown;
    EXPECT_TRUE(file.empty());
  }
  for (const std::string positional : {"my structure.cif", "-5"})
  {
    std::array<std::string, 2> storage{"raspa3", positional};
    std::array<char*, 2> argv{};
    for (std::size_t i = 0; i < storage.size(); ++i) argv[i] = storage[i].data();
    std::string file;
    argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
    parser.pos("FILE", "the input", [&](const std::string& value) { file = value; });
    EXPECT_NO_THROW(parser());
    EXPECT_EQ(file, positional);
  }
}
