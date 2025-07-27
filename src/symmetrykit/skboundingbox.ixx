module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>
#endif

export module skboundingbox;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4x4;
import double3;

export struct SKBoundingBox
{
  SKBoundingBox() : _minimum(double3(0, 0, 0)), _maximum(double3(20, 20, 20)) {}

  SKBoundingBox(double3 minimum, double3 maximum);
  SKBoundingBox(const double3 center, const double3 width, const double scale);
  double3 widths() const;
  std::array<double3, 8> const corners() const;
  std::array<std::pair<double3, double3>, 12> const sides() const;
  double3 center();
  double volume();
  double shortestEdge();
  double longestEdge();
  double3 aspectRatio();
  double boundingSphereRadius();
  double3 maximum() const { return _maximum; }
  double3 minimum() const { return _minimum; }
  SKBoundingBox adjustForTransformation(double4x4 transformation);

  std::int64_t _versionNumber{1};
  double3 _minimum = double3(0.0, 0.0, 0.0);
  double3 _maximum = double3(0.0, 0.0, 0.0);
};

inline SKBoundingBox operator+(const SKBoundingBox left, const SKBoundingBox right)
{
  return SKBoundingBox(left._minimum + right.minimum(), left._maximum + right._maximum);
}

inline SKBoundingBox operator+(const SKBoundingBox left, double3 right)
{
  return SKBoundingBox(left._minimum + right, left._maximum + right);
}

inline SKBoundingBox operator-(const SKBoundingBox left, double3 right)
{
  return SKBoundingBox(left._minimum - right, left._maximum - right);
}
