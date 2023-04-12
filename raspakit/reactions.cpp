module;

module reactions;

import <string>;
import <sstream>;

std::string Reactions::printStatus() const
{
  std::ostringstream stream;

  if (list.empty()) return stream.str();

  std::print(stream, "Reactions:\n");
  std::print(stream, "===============================================================================\n");

  std::print(stream,"{} reactions\n", list.size());
  for (size_t reactionId = 0; const Reaction& reaction : list)
  {
    std::print(stream, reaction.printStatus());
    ++reactionId;
  }
  std::print(stream, "\n\n");

  return stream.str();
}