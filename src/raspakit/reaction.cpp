module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <format>
#include <exception>
#include <source_location>
#include <fstream>
#include <complex>
#include <print>
#endif

module reaction;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <vector>;
import <map>;
import <format>;
import <exception>;
import <source_location>;
import <fstream>;
import <complex>;
import <print>;
#endif


import archive;
import stringutils;


std::string Reaction::printStatus() const
{
  std::ostringstream stream;

  std::string r,p;
  for (const size_t & i : reactantStoichiometry)
    r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const size_t& i : productStoichiometry)
    p += (p.empty() ? "" : ",") + std::to_string(i);
  std::print(stream, "    reaction [{}]: Stoichiometry reactants {} -> products {}\n", id, r, p);

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reaction &r)
{
  archive << r.versionNumber;

  archive << r.id;
  archive << r.reactantStoichiometry;
  archive << r.productStoichiometry;
  archive << r.lambda;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reaction &r)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > r.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Reaction' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> r.id;
  archive >> r.reactantStoichiometry;
  archive >> r.productStoichiometry;
  archive >> r.lambda;

  return archive;
}
