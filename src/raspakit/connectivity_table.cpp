module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <sstream>
#include <map>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>
#endif

module connectivity_table;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;

std::string ConnectivityTable::print(const std::string &prestring) const
{
  std::ostringstream stream;

  for(std::size_t i = 0; i != numberOfBeads - 1; ++i)
  {
    for(std::size_t j = i + 1; j != numberOfBeads; ++j)
    {
      if(table[j * (j + 1) / 2 + i])
      {
        std::print(stream, "{}{} - {}\n", prestring, i, j);
      }
    }
  }

  return stream.str();
}

std::tuple<std::optional<std::size_t>, std::size_t, std::vector<std::size_t>> ConnectivityTable::nextBeads(const std::vector<std::size_t> &placedBeads)
{
  // copy the connectvity from the molecule
  ConnectivityTable grown_connectivity = *this;

  // remove the already grown beads
  for (std::size_t i = 0; i != placedBeads.size(); ++i)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      grown_connectivity[i, j] = false;
    }
  }

  // search for next-bonds, i.e. everything connected to 'placedBeads'
  std::vector<std::pair<std::size_t, std::size_t>> nextBonds{};
  nextBonds.reserve(16);
  for (std::size_t i = 0; i != placedBeads.size(); ++i)
  {
    std::size_t k = placedBeads[i];
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (grown_connectivity[k, j])
      {
        nextBonds.push_back(std::make_pair(k, j));
      }
    }
  }

  if(nextBonds.empty())
  {
    throw std::runtime_error(std::format("Error in CBMMC: No bead can be grown\n"));
  }

  // always select the first
  std::optional<std::size_t> previous_bead{};
  std::size_t current_bead = nextBonds[0].first;
  std::size_t next_bead = nextBonds[0].second;

  // there can be more bonds connected to 'current_bead', so search for these
  std::vector<std::size_t> nextBeads{};
  nextBeads.reserve(nextBonds.size());
  std::size_t number_of_previous_beads{};

  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    if((*this)[i, current_bead])
    {
      if(grown_connectivity[i, current_bead])
      {
        nextBeads.push_back(i);
      }
      else
      {
        ++number_of_previous_beads;
        previous_bead = i;
      }
    }
  }

  // if there are no previous beads, than grow the chosen, single bond
  if(number_of_previous_beads == 0)
  {
    return {std::nullopt, current_bead, {next_bead}};
  }
  else
  {
    if(!previous_bead)
    {
      throw std::runtime_error(std::format("Error in CBMC: No previous bead\n"));
    }
    else if (number_of_previous_beads > 1)
    {
      throw std::runtime_error(std::format("Error in CBMC: Multiple previous beads\n"));
    }
  }

  return {previous_bead, current_bead, nextBeads};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ConnectivityTable &b)
{
  archive << b.versionNumber;

  archive << b.numberOfBeads;
  archive << b.table;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ConnectivityTable &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ConnectivityTable' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.numberOfBeads;
  archive >> b.table;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("ConnectivityTable: Error in binary restart\n"));
  }
#endif

  return archive;
}
