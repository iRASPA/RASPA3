module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>
#endif

export module special_functions;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <numeric>;
import <algorithm>;
import <string>;
#endif

export extern double li2(double x);
export extern double hypergeometric2F1(double a, double b, double c, double z);
export extern double hypergeometric(double a, double b, double c, double x);

export template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}
