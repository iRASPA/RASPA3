module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module reactions;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;
import reaction;
import json;

std::string Reactions::printStatus() const
{
  std::ostringstream stream;

  if (list.empty()) return stream.str();

  std::print(stream, "Reactions:\n");
  std::print(stream, "===============================================================================\n");

  std::print(stream, "{} reactions\n", list.size());
  for (const Reaction &reaction : list)
  {
    std::print(stream, "{}", reaction.printStatus());
  }
  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json Reactions::jsonStatus() const
{
  nlohmann::json status;

  if (list.empty()) return status;
  status["n_reactions"] = list.size();

  nlohmann::json reactions(list.size());
  for (std::size_t i = 0; i < list.size(); i++) reactions[i] = list[i].jsonStatus();
  status["reactions"] = reactions;
  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reactions &r)
{
  archive << r.versionNumber;
  archive << r.list;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reactions &r)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > r.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Reactions' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> r.list;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Reactions: Error in binary restart\n"));
  }
#endif

  return archive;
}
