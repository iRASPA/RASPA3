#include <gtest/gtest.h>

import std;

import archive;
import randomnumbers;

// The restart file stores the state of the std::mt19937_64 engine in its standard textual form
// (312 decimal words, [rand.eng.mers]); restoring is O(1) instead of replaying 'count' draws with
// discard(). Files in the legacy 'seed, count' format must still read.

namespace
{
std::filesystem::path temporaryFile(const char* name) { return std::filesystem::temp_directory_path() / name; }

template <typename T>
void writeArchive(const std::filesystem::path& path, const T& object)
{
  std::ofstream stream(path, std::ios::binary);
  Archive<std::ofstream> archive(stream);
  archive << object;
}

template <typename T>
void readArchive(const std::filesystem::path& path, T& object)
{
  std::ifstream stream(path, std::ios::binary);
  Archive<std::ifstream> archive(stream);
  archive >> object;
}
}  // namespace

TEST(RANDOM_ARCHIVE, engine_state_round_trip_continues_the_sequence)
{
  RandomNumber original(std::size_t{20240928});
  for (std::size_t i = 0; i < 12345; ++i) original.uniform();
  for (std::size_t i = 0; i < 77; ++i) original.Gaussian();
  const std::size_t countAtWrite = original.count;

  const std::filesystem::path path = temporaryFile("raspa3_random_archive.bin");
  writeArchive(path, original);  // also drops the cached Box-Muller value of 'original'
  RandomNumber restored;
  readArchive(path, restored);
  std::filesystem::remove(path);

  EXPECT_EQ(restored, original);
  EXPECT_EQ(restored.seed, std::size_t{20240928});
  EXPECT_EQ(restored.count, countAtWrite);

  // uninterrupted run and restarted run stay on the same sequence, including the Gaussian draws
  for (std::size_t i = 0; i < 1000; ++i)
  {
    EXPECT_EQ(restored.uniform(), original.uniform());
    EXPECT_EQ(restored.Gaussian(), original.Gaussian());
    EXPECT_EQ(restored.uniform_integer(0, 1000), original.uniform_integer(0, 1000));
  }
}

TEST(RANDOM_ARCHIVE, legacy_seed_count_format_still_reads)
{
  RandomNumber reference(std::size_t{4242});
  for (std::size_t i = 0; i < 5000; ++i) reference.uniform();

  const std::filesystem::path path = temporaryFile("raspa3_random_archive_legacy.bin");
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << std::size_t{4242};  // seed (the first word is not the marker, so this is the legacy layout)
    archive << std::size_t{5000};  // count: one engine call per uniform()
  }
  RandomNumber restored;
  readArchive(path, restored);
  std::filesystem::remove(path);

  EXPECT_EQ(restored.seed, std::size_t{4242});
  EXPECT_EQ(restored.count, std::size_t{5000});
  for (std::size_t i = 0; i < 100; ++i) EXPECT_EQ(restored.uniform(), reference.uniform());
}

TEST(RANDOM_ARCHIVE, textual_engine_state_is_standard_form)
{
  RandomNumber random(std::size_t{7});
  for (std::size_t i = 0; i < 1000; ++i) random.uniform();

  const std::string state = random.engineState();

  // 312 white-space separated unsigned decimal integers, nothing else
  std::istringstream stream(state);
  std::size_t words = 0;
  std::string word;
  while (stream >> word)
  {
    ++words;
    EXPECT_TRUE(std::ranges::all_of(word, [](char c) { return c >= '0' && c <= '9'; })) << word;
  }
  EXPECT_EQ(words, std::size_t{312});

  RandomNumber other;
  other.setEngineState(state);
  EXPECT_EQ(other.mt, random.mt);
  EXPECT_EQ(other.mt(), random.mt());

  EXPECT_THROW(other.setEngineState("not a state"), std::runtime_error);
}
