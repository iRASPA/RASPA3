module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <iostream>
#include <vector>
#endif

export module multi_site_isotherm;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import hashcombine;
import randomnumbers;
import isotherm;

export struct MultiSiteIsotherm
{
  enum class PredictionMethod
  {
    IAST = 0,
    SIAST = 1,
    EI = 2,
    SEI = 3
  };

  bool operator==(MultiSiteIsotherm const &) const = default;

  std::size_t numberOfSites{0};
  std::vector<Isotherm> sites{};

  std::size_t numberOfParameters{0};
  std::vector<std::pair<std::size_t, std::size_t>> parameterIndices{};
  std::vector<std::size_t> siteParameterIndex{};

  double &parameters(std::size_t i)
  {
    std::pair<std::size_t, std::size_t> index = parameterIndices[i];
    return sites[index.first].parameters[index.second];
  }
  const double &parameters(std::size_t i) const
  {
    std::pair<std::size_t, std::size_t> index = parameterIndices[i];
    return sites[index.first].parameters[index.second];
  }

  void add(const Isotherm &isotherm);

  std::string print() const;
  std::string printAsInputFormat() const;

  MultiSiteIsotherm randomized(RandomNumber &random, double maximumLoading)
  {
    MultiSiteIsotherm copy(*this);
    for (std::size_t i = 0; i < numberOfSites; ++i)
    {
      copy.sites[i].randomize(random, maximumLoading);
    }
    return copy;
  }

  inline double value(double pressure) const
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].value(pressure);
    }
    return sum;
  }

  inline double value(std::size_t site, double pressure) const
  {
    if (site < numberOfSites)
    {
      return sites[site].value(pressure);
    }
    return 0.0;
  }

  // computed reduced grand potential for pressure
  inline double psiForPressure(double pressure) const
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].psiForPressure(pressure);
    }
    return sum;
  }

  // computed reduced grand potential for pressure
  inline double psiForPressure(std::size_t site, double pressure) const
  {
    if (site < numberOfSites)
    {
      return sites[site].psiForPressure(pressure);
    }
    return 0.0;
  }

  double inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const;

  double inversePressureForPsi(std::size_t site, double reduced_grand_potential, double &cachedP0) const
  {
    if (site < numberOfSites)
    {
      return sites[site].inversePressureForPsi(reduced_grand_potential, cachedP0);
    }
    return 0.0;
  }

  double fitness() const;
  std::string gnuplotFunctionString(char s) const;
};

namespace std
{
export template <>
struct hash<MultiSiteIsotherm>
{
  std::size_t operator()(const MultiSiteIsotherm &k) const
  {
    std::size_t h = 0;
    for (const Isotherm &isotherm : k.sites)
    {
      for (std::size_t i = 0; i < isotherm.numberOfParameters; ++i)
      {
        hash_combine(h, isotherm.parameters[i]);
      }
    }
    return h;
  }
};
}  // namespace std
