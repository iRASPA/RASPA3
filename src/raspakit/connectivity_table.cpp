module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module connectivity_table;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import stringutils;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

#ifdef BLAS_ILP64
typedef long long blas_int;
#else
typedef int blas_int;
#endif

/* DSYEV prototype */
extern "C"
{
  void dsyev_(char *jobz, char *uplo, blas_int *n, double *a, blas_int *lda, double *w, double *work, blas_int *lwork,
              blas_int *info);
}

std::string ConnectivityTable::print(const std::string &prestring) const
{
  std::ostringstream stream;

  for (std::size_t i = 0; i != numberOfBeads - 1; ++i)
  {
    for (std::size_t j = i + 1; j != numberOfBeads; ++j)
    {
      if (table[i * numberOfBeads + j])
      {
        std::print(stream, "{}{} - {}\n", prestring, i, j);
      }
    }
  }

  return stream.str();
}

std::tuple<std::optional<std::size_t>, std::size_t, std::vector<std::size_t>> ConnectivityTable::nextBeads(
    const std::vector<std::size_t> &placedBeads) const
{
  // copy the connectvity from the molecule
  ConnectivityTable to_do_connectivity = *this;

  // remove the already grown beads
  for (const std::size_t placed_bead : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      // note: update is asymmetric (only [placed_bead, j], not [j, placed_bead])
      to_do_connectivity[placed_bead, j] = false;
    }
  }

  // search for next-bonds, i.e. everything connected to 'placedBeads'
  std::vector<std::pair<std::size_t, std::size_t>> nextBonds{};
  nextBonds.reserve(16);
  for (const std::size_t k : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (to_do_connectivity[j, k])
      {
        nextBonds.push_back(std::make_pair(k, j));
      }
    }
  }

  if (nextBonds.empty())
  {
    throw std::runtime_error(std::format("Error in CBMC: No bead can be grown\n"));
  }

  // always select the first for reversibility in coupled/decoupled
  std::optional<std::size_t> previous_bead{};
  std::size_t current_bead = nextBonds[0].first;
  std::size_t next_bead = nextBonds[0].second;

  // there can be more bonds connected to 'current_bead', so search for these
  std::vector<std::size_t> nextBeads{};
  nextBeads.reserve(nextBonds.size());
  std::size_t number_of_previous_beads{};

  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    if ((*this)[i, current_bead])
    {
      if (to_do_connectivity[i, current_bead])
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
  if (number_of_previous_beads == 0)
  {
    return {std::nullopt, current_bead, {next_bead}};
  }
  else
  {
    if (!previous_bead)
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

std::vector<std::size_t> ConnectivityTable::findAllNeighbors(std::size_t currentBead) const
{
  std::vector<std::size_t> neighbors{};
  neighbors.reserve(numberOfBeads);

  for (std::size_t i = 0; i < numberOfBeads; ++i)
  {
    if (((*this)[currentBead, i]) && (i != currentBead))
    {
      neighbors.push_back(i);
    }
  }
  return neighbors;
}

// https://en.wikipedia.org/wiki/Algebraic_connectivity
// https://en.wikipedia.org/wiki/Laplacian_matrix
// The algebraic connectivity (also known as Fiedler value or Fiedler eigenvalue after Miroslav Fiedler)
// of a graph G is the second-smallest eigenvalue (counting multiple eigenvalues separately) of the
// Laplacian matrix of G.[1] This eigenvalue is greater than 0 if and only if G is a connected graph.
// the number of times 0 appears as an eigenvalue in the Laplacian is the number of connected components in the graph
bool ConnectivityTable::checkIsConnectedSubgraph(const std::vector<std::size_t> &set)
{
  if (set.empty() || (set.size() == 1))
  {
    return true;
  }

  std::vector<double> laplacian_matrix(set.size() * set.size());
  std::mdspan<double, std::dextents<std::size_t, 2>, std::layout_left> data_laplacian_matrix(laplacian_matrix.data(),
                                                                                             set.size(), set.size());

  for (std::size_t i = 0; i < set.size(); ++i)
  {
    for (std::size_t j = 0; j < set.size(); ++j)
    {
      if (i == j)
      {
        for (std::size_t k = 0; k < set.size(); ++k)
        {
          if (table[set[i] * numberOfBeads + set[k]])
          {
            data_laplacian_matrix[i, i] += 1.0;
          }
        }
      }
      else
      {
        data_laplacian_matrix[i, j] = table[set[i] * numberOfBeads + set[j]] ? -1.0 : 0.0;
      }
    }
  }

  char decompositionJobV = 'N';
  char upload = 'U';
  std::vector<double> work(9 * set.size());
  blas_int lwork = 9 * static_cast<blas_int>(set.size());
  std::vector<double> e = std::vector<double>(set.size());
  blas_int error = 0;
  blas_int N = static_cast<blas_int>(set.size());
  blas_int M = static_cast<blas_int>(set.size());

  dsyev_(&decompositionJobV, &upload, &M, laplacian_matrix.data(), &N, e.data(), work.data(), &lwork, &error);

  return std::fabs(e[1]) > 1e-8;
}

std::optional<std::vector<std::size_t>> ConnectivityTable::checkValidityNextBeads(
    const std::vector<std::size_t> &placedBeads) const
{
  // copy the connectvity from the molecule
  ConnectivityTable to_do_connectivity = *this;

  // remove the already grown beads
  for (const std::size_t placed_bead : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      // note: update is asymmetric (only [placed_bead, j], not [j, placed_bead])
      to_do_connectivity[placed_bead, j] = false;
    }
  }

  // search for next-bonds, i.e. everything connected to 'placedBeads'
  std::vector<std::pair<std::size_t, std::size_t>> nextBonds{};
  nextBonds.reserve(16);
  for (const std::size_t k : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (to_do_connectivity[j, k])
      {
        nextBonds.push_back(std::make_pair(k, j));
      }
    }
  }

  if (nextBonds.empty())
  {
    return std::nullopt;
  }

  // always select the first for reversibility in coupled/decoupled
  std::optional<std::size_t> previous_bead{};
  std::size_t current_bead = nextBonds[0].first;
  std::size_t next_bead = nextBonds[0].second;

  // there can be more bonds connected to 'current_bead', so search for these
  std::vector<std::size_t> nextBeads{};
  nextBeads.reserve(nextBonds.size());
  std::size_t number_of_previous_beads{};

  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    if ((*this)[i, current_bead])
    {
      if (to_do_connectivity[i, current_bead])
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
  if (number_of_previous_beads == 0)
  {
    return {{next_bead}};
  }
  else
  {
    if (!previous_bead)
    {
      return std::nullopt;
    }
    else if (number_of_previous_beads > 1)
    {
      return std::nullopt;
    }
  }

  return {nextBeads};
}

std::vector<std::array<std::size_t, 2>> ConnectivityTable::findAllBonds() const
{
  std::vector<std::array<std::size_t, 2>> result{};

  for (std::size_t i = 0; i < numberOfBeads - 1; ++i)
  {
    for (std::size_t j = i + 1; j < numberOfBeads; ++j)
    {
      if ((*this)[i, j])
      {
        result.push_back({i, j});
      }
    }
  }

  return result;
}

std::vector<std::array<std::size_t, 3>> ConnectivityTable::findAllBends() const
{
  std::vector<std::array<std::size_t, 3>> result{};

  std::vector<std::array<std::size_t, 2>> found_bonds = findAllBonds();

  for (const std::array<std::size_t, 2> found_bond : found_bonds)
  {
    // extend to the left
    for (std::size_t left_bead : findAllNeighbors(found_bond[0]))
    {
      if (left_bead != found_bond[1])
      {
        std::array<std::size_t, 3> bend = {std::min(left_bead, found_bond[1]), found_bond[0],
                                           std::max(left_bead, found_bond[1])};
        if (std::find(result.begin(), result.end(), bend) == result.end())
        {
          result.push_back(bend);
        }
      }
    }

    // extend to the right
    for (std::size_t right_bead : findAllNeighbors(found_bond[1]))
    {
      if (right_bead != found_bond[0])
      {
        std::array<std::size_t, 3> bend = {std::min(found_bond[0], right_bead), found_bond[1],
                                           std::max(found_bond[0], right_bead)};
        if (std::find(result.begin(), result.end(), bend) == result.end())
        {
          result.push_back(bend);
        }
      }
    }
  }

  return result;
}

std::vector<std::array<std::size_t, 4>> ConnectivityTable::findAllTorsions() const
{
  std::vector<std::array<std::size_t, 4>> result{};

  std::vector<std::array<std::size_t, 3>> found_bends = findAllBends();

  for (const std::array<std::size_t, 3> found_bend : found_bends)
  {
    // extend to the left
    for (std::size_t left_bead : findAllNeighbors(found_bend[0]))
    {
      if (left_bead != found_bend[1])
      {
        std::array<std::size_t, 4> torsion = {left_bead, found_bend[0], found_bend[1], found_bend[2]};
        std::array<std::size_t, 4> torsion_mirrored = {found_bend[2], found_bend[1], found_bend[0], left_bead};
        if (std::find(result.begin(), result.end(), torsion) == result.end() &&
            std::find(result.begin(), result.end(), torsion_mirrored) == result.end())
        {
          result.push_back(torsion);
        }
      }
    }

    // extend to the right
    for (std::size_t right_bead : findAllNeighbors(found_bend[2]))
    {
      if (right_bead != found_bend[1])
      {
        std::array<std::size_t, 4> torsion = {found_bend[0], found_bend[1], found_bend[2], right_bead};
        std::array<std::size_t, 4> torsion_mirrored = {right_bead, found_bend[2], found_bend[1], found_bend[0]};

        if (std::find(result.begin(), result.end(), torsion) == result.end() &&
            std::find(result.begin(), result.end(), torsion_mirrored) == result.end())
        {
          result.push_back(torsion);
        }
      }
    }
  }

  return result;
}

std::vector<std::array<std::size_t, 2>> ConnectivityTable::findAllVanDerWaals() const
{
  std::vector<std::array<std::size_t, 2>> result{};
  result.reserve(numberOfBeads * (numberOfBeads - 1) / 2);

  for (std::size_t i = 0; i < numberOfBeads - 1; ++i)
  {
    for (std::size_t j = i + 1; j < numberOfBeads; ++j)
    {
      result.push_back({i, j});
    }
  }

  std::vector<std::array<std::size_t, 2>> found_bonds = findAllBonds();
  std::vector<std::array<std::size_t, 3>> found_bends = findAllBends();
  std::vector<std::array<std::size_t, 4>> found_torsions = findAllTorsions();

  for (std::array<std::size_t, 2> &found_bond : found_bonds)
  {
    result.erase(std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_bond[0], found_bond[1]}),
                 result.end());
    result.erase(std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_bond[1], found_bond[0]}),
                 result.end());
  }

  for (std::array<std::size_t, 3> &found_bend : found_bends)
  {
    result.erase(std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_bend[0], found_bend[2]}),
                 result.end());
    result.erase(std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_bend[2], found_bend[0]}),
                 result.end());
  }

  for (std::array<std::size_t, 4> &found_torsion : found_torsions)
  {
    result.erase(
        std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_torsion[0], found_torsion[3]}),
        result.end());
    result.erase(
        std::remove(result.begin(), result.end(), std::array<std::size_t, 2>{found_torsion[3], found_torsion[0]}),
        result.end());
  }

  return result;
}

const std::string ConnectivityTable::printConnectivityTable() const
{
  std::ostringstream stream;
  std::print(stream, "Connectivity Table\n");
  std::print(stream, "------------------------------------------------\n");
  for (std::size_t i = 0; i < numberOfBeads; i++)
  {
    for (std::size_t j = 0; j < numberOfBeads; j++)
    {
      std::print(stream, "{} ", table[i * numberOfBeads + j] ? "1" : "0");
    }
    std::print(stream, "\n");
  }
  std::print(stream, "\n");
  return stream.str();
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
