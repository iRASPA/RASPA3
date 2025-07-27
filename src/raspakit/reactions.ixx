module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>
#endif

export module reactions;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import reaction;
import json;

/**
 * \brief Manages a collection of reactions within the simulation system.
 *
 * The Reactions struct encapsulates a list of Reaction objects and provides methods
 * to retrieve the current status in both string and JSON formats. It also includes
 * serialization operators for reading from and writing to archives.
 */
export struct Reactions
{
  ///< The version number of the Reactions data structure.
  std::uint64_t versionNumber{1};

  bool operator==(Reactions const &) const = default;

  ///< A list of Reaction objects managed by this struct.
  std::vector<Reaction> list;

  /**
   * \brief Retrieves the current status of all reactions as a string.
   *
   * Generates a formatted string that includes the number of reactions and the status
   * of each individual reaction in the list.
   *
   * \return A string representation of the current reactions status.
   */
  std::string printStatus() const;

  /**
   * \brief Retrieves the current status of all reactions in JSON format.
   *
   * Constructs a JSON object containing the number of reactions and an array
   * representing each reaction's status.
   *
   * \return A JSON object representing the current reactions status.
   */
  nlohmann::json jsonStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reactions &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reactions &r);
};
