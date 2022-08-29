export module multi_site_isotherm;

import <cstddef>;
import <array>;
import <vector>;
import <iostream>;

import isotherm;
import hashcombine;

export struct MultiSiteIsotherm
{
  MultiSiteIsotherm() noexcept = default;
  MultiSiteIsotherm(const MultiSiteIsotherm &a) noexcept = default;
  MultiSiteIsotherm& operator=(const MultiSiteIsotherm& a) noexcept = default;
  MultiSiteIsotherm(MultiSiteIsotherm&& a) noexcept = default;
  MultiSiteIsotherm& operator=(MultiSiteIsotherm&& a) noexcept = default;
  ~MultiSiteIsotherm() noexcept = default;

  size_t numberOfSites{ 0 };
  std::vector<Isotherm> sites{};

  std::vector<double> parameters{};
  std::vector<size_t> siteParameterIndex{};

  std::string print() const;
  std::string printAsInputFormat() const;
  void add(const Isotherm &isotherm, const std::vector<double> &params);

  MultiSiteIsotherm randomized(double maximumLoading)
  {
    MultiSiteIsotherm copy(*this);
    for(size_t i = 0; i < numberOfSites; ++i)
    {
      copy.sites[i].randomize(&parameters[siteParameterIndex[i]], maximumLoading);
    }
    return copy;
  }

  inline double value(double pressure) const
  {
    double sum = 0.0;
    for(size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].value(pressure, &parameters[siteParameterIndex[i]]);
    }
    return sum;
  }

  // computed reduced grand potential for pressure
  inline double psiForPressure(double pressure) const
  {
    double sum = 0.0;
    for(size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].psiForPressure(pressure, &parameters[siteParameterIndex[i]]);
    }
    return sum;
  }

  double inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const;

  double fitness() const;
  std::string gnuplotFunctionString(char s) const;
};

export namespace std
{
  template <> struct hash<MultiSiteIsotherm>
  {
    size_t operator()(const MultiSiteIsotherm& k) const
    {
      std::size_t h=0;
      for(const double &parameter: k.parameters)
      {
        hash_combine(h, parameter);
      }
      return h;
    }
  };
}
