module;

#ifdef USE_LEGACY_HEADERS
#include <tuple>
#include <vector>
#endif

export module cbmc_growing_status;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
#endif

import component;
import atom;
import bond_potential;

export struct GrowingStatus
{
  GrowingStatus(const Component& component, const std::vector<bool>& placedBeads)
      : size(component.atoms.size()), placement(size, false), connectivity(size * size, false)
  {
    // for (size_t index = 0; const bool& bead : placedBeads)
    //{
    //     placement[index] = true;
    //     ++index;
    // }

    for (size_t i = 0; i != component.bonds.size(); ++i)
    {
      // connectivity[component.bonds[i].first + size * component.bonds[i].second] = true;
      // connectivity[component.bonds[i].second + size * component.bonds[i].first] = true;
    }

    for (size_t i = 0; i != placedBeads.size(); ++i)
    {
      for (size_t j = 0; j != component.atoms.size(); ++j)
      {
        connectivity[i + size * j] = false;
        connectivity[j + size * i] = false;
      }
    }

    for (size_t i = 0; i != placedBeads.size(); ++i)
    {
      size_t k = placedBeads[i];
      for (size_t j = 0; j != component.atoms.size(); ++j)
      {
        if (connectivity[k + size * j])
        {
          currentNextBeadPairs.push_back(std::make_pair(k, j));
        }
      }
    }
  }

  size_t size;
  std::vector<bool> placement;
  std::vector<bool> connectivity;
  std::vector<std::pair<size_t, size_t>> currentNextBeadPairs{};
};
