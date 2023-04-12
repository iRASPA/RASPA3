module;

module reactions;

import <string>;
import <sstream>;

import print;
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
    std::print(stream, reaction.printStatus());
  }
  std::print(stream, "\n\n");

  return stream.str();
}
