module;

module reactions;

import <string>;
import <sstream>;
import <ostream>;
import <vector>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <fstream>;
import <complex>;

import archive;
import stringutils;
import reaction;


std::string Reactions::printStatus() const
{
  std::ostringstream stream;

  if (list.empty()) return stream.str();

  std::print(stream, "Reactions:\n");
  std::print(stream, "===============================================================================\n");

  std::print(stream,"{} reactions\n", list.size());
  for (const Reaction& reaction : list)
  {
    std::print(stream, "{}", reaction.printStatus());
  }
  std::print(stream, "\n\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reactions &r)
{
  archive << r.versionNumber;
  archive << r.list;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reactions &r)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > r.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Reactions' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> r.list;

  return archive;
}
