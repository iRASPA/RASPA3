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

module reaction;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;

std::string Reaction::printStatus() const
{
  std::ostringstream stream;

  std::string r, p;
  for (const std::size_t &i : reactantStoichiometry) r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const std::size_t &i : productStoichiometry) p += (p.empty() ? "" : ",") + std::to_string(i);
  std::print(stream, "    reaction [{}]: Stoichiometry reactants {} -> products {}\n", id, r, p);

  return stream.str();
}

nlohmann::json Reaction::jsonStatus() const
{
  nlohmann::json status;
  status["id"] = id;

  std::string r, p;
  for (const std::size_t &i : reactantStoichiometry) r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const std::size_t &i : productStoichiometry) p += (p.empty() ? "" : ",") + std::to_string(i);
  status["reactants"] = r;
  status["products"] = p;

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reaction &r)
{
  archive << r.versionNumber;

  archive << r.id;
  archive << r.reactantStoichiometry;
  archive << r.productStoichiometry;
  archive << r.lambda;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reaction &r)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > r.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Reaction' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> r.id;
  archive >> r.reactantStoichiometry;
  archive >> r.productStoichiometry;
  archive >> r.lambda;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Reaction: Error in binary restart\n"));
  }
#endif

  return archive;
}
