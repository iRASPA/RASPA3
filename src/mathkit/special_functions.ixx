module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>
#endif

export module special_functions;

#ifdef USE_STD_IMPORT
import std;
#endif

export extern double li2(double x);
export extern double hypergeometric2F1(double a, double b, double c, double z);
export extern double hypergeometric(double a, double b, double c, double x);

export template <typename T>
std::vector<std::size_t> sort_indexes(const std::vector<T> &v)
{
  // initialize original index locations
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });

  return idx;
}
