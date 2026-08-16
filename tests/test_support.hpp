#pragma once

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <random>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace test_support
{
inline std::string uniqueSuffix()
{
  static thread_local std::mt19937 rng{std::random_device{}()};
  return std::format("{}_{}", std::chrono::steady_clock::now().time_since_epoch().count(), rng());
}
}  // namespace test_support

// RAII helpers for materializing embedded fixture strings to the filesystem so that
// production loaders (Component / ForceField / InputReader) can consume them.

class TemporaryFile
{
 public:
  TemporaryFile(std::string name, std::string_view contents)
      : path_(std::filesystem::temp_directory_path() /
              std::format("raspa3_test_{}_{}", test_support::uniqueSuffix(), std::move(name)))
  {
    std::ofstream stream(path_);
    stream << contents;
  }

  TemporaryFile(const TemporaryFile&) = delete;
  TemporaryFile& operator=(const TemporaryFile&) = delete;

  TemporaryFile(TemporaryFile&& other) noexcept : path_(std::move(other.path_)) { other.path_.clear(); }

  TemporaryFile& operator=(TemporaryFile&& other) noexcept
  {
    if (this != &other)
    {
      remove();
      path_ = std::move(other.path_);
      other.path_.clear();
    }
    return *this;
  }

  ~TemporaryFile() { remove(); }

  const std::filesystem::path& path() const { return path_; }

  // Path without the .json extension, matching the Component(fileName) convention.
  std::filesystem::path stemPath() const
  {
    std::filesystem::path stem = path_;
    stem.replace_extension();
    return stem;
  }

 private:
  void remove()
  {
    if (!path_.empty())
    {
      std::error_code ignored;
      std::filesystem::remove(path_, ignored);
      path_.clear();
    }
  }

  std::filesystem::path path_;
};

class TemporaryDirectory
{
 public:
  TemporaryDirectory()
      : path_(std::filesystem::temp_directory_path() /
              std::format("raspa3_test_dir_{}", test_support::uniqueSuffix()))
  {
    std::filesystem::create_directories(path_);
  }

  TemporaryDirectory(const TemporaryDirectory&) = delete;
  TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

  TemporaryDirectory(TemporaryDirectory&& other) noexcept : path_(std::move(other.path_)) { other.path_.clear(); }

  TemporaryDirectory& operator=(TemporaryDirectory&& other) noexcept
  {
    if (this != &other)
    {
      remove();
      path_ = std::move(other.path_);
      other.path_.clear();
    }
    return *this;
  }

  ~TemporaryDirectory() { remove(); }

  const std::filesystem::path& path() const { return path_; }

  void write(std::string_view relativeName, std::string_view contents) const
  {
    const std::filesystem::path filePath = path_ / relativeName;
    if (filePath.has_parent_path())
    {
      std::filesystem::create_directories(filePath.parent_path());
    }
    std::ofstream stream(filePath);
    stream << contents;
  }

 private:
  void remove()
  {
    if (!path_.empty())
    {
      std::error_code ignored;
      std::filesystem::remove_all(path_, ignored);
      path_.clear();
    }
  }

  std::filesystem::path path_;
};
