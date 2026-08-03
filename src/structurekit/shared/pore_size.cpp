module;

#ifdef _OPENMP
#include <omp.h>
#endif

module grid_pore_size;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;


double coveringSlack(const UnitCell &unitCell, uint3 gridSize)
{
  double3 step(1.0 / static_cast<double>(gridSize.x), 1.0 / static_cast<double>(gridSize.y),
               1.0 / static_cast<double>(gridSize.z));
  double3 halfStep = 0.5 * (unitCell.cell * step);
  return 2.0 * halfStep.length();
}


std::vector<float> distanceToIsosurface(uint3 gridSize, const UnitCell &unitCell,
                                        std::span<const float> openness, double lowestOpenness)
{
  const std::int64_t nx = static_cast<std::int64_t>(gridSize.x);
  const std::int64_t ny = static_cast<std::int64_t>(gridSize.y);
  const std::int64_t nz = static_cast<std::int64_t>(gridSize.z);
  const std::size_t numberOfVoxels = openness.size();

  std::vector<float> distance(numberOfVoxels, -1.0f);
  if (numberOfVoxels == 0) return distance;

  const double3x3 cell = unitCell.cell;
  const float level = static_cast<float>(lowestOpenness);

  // A step of one grid point along an axis carries the point at least this far towards the opposite face,
  // whatever the shape of the cell. Reaching a Chebyshev radius of r therefore costs at least r of these, so
  // it bounds from below what the next shell out could possibly find and is what lets the search stop.
  double3 widths = unitCell.perpendicularWidths();
  const double shortestStep = std::min({widths.x / static_cast<double>(nx), widths.y / static_cast<double>(ny),
                                        widths.z / static_cast<double>(nz)});

  // The three grid steps as displacements, and the longest of them. A crossing sits on an edge, so it lies
  // within one step of the grid point that found it, and the stopping bound has to allow for that.
  const std::array<double3, 3> axisStep{cell * double3(1.0 / static_cast<double>(nx), 0.0, 0.0),
                                        cell * double3(0.0, 1.0 / static_cast<double>(ny), 0.0),
                                        cell * double3(0.0, 0.0, 1.0 / static_cast<double>(nz))};
  const double longestStep =
      std::max({axisStep[0].length(), axisStep[1].length(), axisStep[2].length()});

  // A void filling the whole cell has no boundary to find, and the search would otherwise walk the grid.
  const std::int64_t largestReach = std::max({nx, ny, nz}) / 2 + 1;

#pragma omp parallel for schedule(dynamic, 8)
  for (std::int64_t voxel = 0; voxel < static_cast<std::int64_t>(numberOfVoxels); ++voxel)
  {
    const float here = openness[static_cast<std::size_t>(voxel)];
    if (here < level) continue;

    const std::int64_t k = voxel / (nx * ny);
    const std::int64_t j = (voxel / nx) % ny;
    const std::int64_t i = voxel % nx;

    double best = std::numeric_limits<double>::max();

    for (std::int64_t reach = 1; reach <= largestReach; ++reach)
    {
      // Nothing at this radius or beyond can beat what is already in hand. A crossing found from a point at
      // this radius can lie a step nearer than the point itself, which the bound has to allow for.
      if (static_cast<double>(reach - 1) * shortestStep - longestStep >= best) break;

      for (std::int64_t dk = -reach; dk <= reach; ++dk)
      {
        for (std::int64_t dj = -reach; dj <= reach; ++dj)
        {
          // Only the shell, since everything inside it was walked on an earlier pass.
          const bool onFace = (std::abs(dk) == reach) || (std::abs(dj) == reach);
          const std::int64_t stride = onFace ? 1 : 2 * reach;
          for (std::int64_t di = -reach; di <= reach; di += stride)
          {
            const std::int64_t ii = ((i + di) % nx + nx) % nx;
            const std::int64_t jj = ((j + dj) % ny + ny) % ny;
            const std::int64_t kk = ((k + dk) % nz + nz) % nz;
            const float outside = openness[static_cast<std::size_t>((kk * ny + jj) * nx + ii)];
            if (outside >= level) continue;

            double3 fractional(static_cast<double>(di) / static_cast<double>(nx),
                               static_cast<double>(dj) / static_cast<double>(ny),
                               static_cast<double>(dk) / static_cast<double>(nz));
            double3 offset = cell * fractional;
            if (offset.length() - longestStep >= best) continue;

            // The surface is placed where it crosses the edges of the grid, and only ever between two points
            // that straddle it, which is where marching cubes puts its vertices too. Interpolating over any
            // longer a span would be worse than useless on a field of this kind: an energy climbs as the
            // twelfth power of separation, so a point in the wall is not larger than a point in the pore by
            // some modest factor but by ten orders of magnitude, and a straight line drawn to it puts the
            // crossing almost on top of the pore point. Distances measured that way come out far too small,
            // and, worse, too small by an amount that depends on how deep the wall behind them happens to
            // be.
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
              for (const std::int64_t step : {std::int64_t{-1}, std::int64_t{1}})
              {
                std::int64_t ni = ii, nj = jj, nk = kk;
                if (axis == 0) ni = ((ii + step) % nx + nx) % nx;
                if (axis == 1) nj = ((jj + step) % ny + ny) % ny;
                if (axis == 2) nk = ((kk + step) % nz + nz) % nz;

                const float inside = openness[static_cast<std::size_t>((nk * ny + nj) * nx + ni)];
                if (inside < level) continue;

                const double climb = static_cast<double>(inside) - static_cast<double>(outside);
                const double fraction =
                    (climb > 0.0) ? (static_cast<double>(level) - static_cast<double>(outside)) / climb : 0.0;
                double3 crossing = offset + std::clamp(fraction, 0.0, 1.0) * static_cast<double>(step) *
                                                axisStep[axis];
                best = std::min(best, crossing.length());
              }
            }
          }
        }
      }
    }

    distance[static_cast<std::size_t>(voxel)] =
        (best == std::numeric_limits<double>::max()) ? std::numeric_limits<float>::max() : static_cast<float>(best);
  }

  return distance;
}


std::vector<float> poreRadiusField(uint3 gridSize, const UnitCell &unitCell,
                                   std::span<const float> distance, double slack)
{
  const std::int64_t nx = static_cast<std::int64_t>(gridSize.x);
  const std::int64_t ny = static_cast<std::int64_t>(gridSize.y);
  const std::int64_t nz = static_cast<std::int64_t>(gridSize.z);
  const std::size_t numberOfVoxels = distance.size();

  std::vector<float> poreRadius(numberOfVoxels, 0.0f);
  if (numberOfVoxels == 0) return poreRadius;

  const double3x3 cell = unitCell.cell;

  // How many grid steps a ball of a given radius can possibly reach along each axis. The perpendicular width
  // is the right measure for an oblique cell: a step along one axis moves at least that far towards the
  // opposite face, so the box built from it holds the ball.
  double3 widths = unitCell.perpendicularWidths();
  const double stepsPerLengthX = static_cast<double>(nx) / widths.x;
  const double stepsPerLengthY = static_cast<double>(ny) / widths.y;
  const double stepsPerLengthZ = static_cast<double>(nz) / widths.z;

#pragma omp parallel for schedule(dynamic, 8)
  for (std::int64_t voxel = 0; voxel < static_cast<std::int64_t>(numberOfVoxels); ++voxel)
  {
    const float radius = distance[static_cast<std::size_t>(voxel)];
    if (!(radius > 0.0f) || radius == std::numeric_limits<float>::max()) continue;

    const std::int64_t k = voxel / (nx * ny);
    const std::int64_t j = (voxel / nx) % ny;
    const std::int64_t i = voxel % nx;

    const double reach = static_cast<double>(radius) + slack;
    const std::int64_t reachX = static_cast<std::int64_t>(std::ceil(reach * stepsPerLengthX));
    const std::int64_t reachY = static_cast<std::int64_t>(std::ceil(reach * stepsPerLengthY));
    const std::int64_t reachZ = static_cast<std::int64_t>(std::ceil(reach * stepsPerLengthZ));

    for (std::int64_t dk = -reachZ; dk <= reachZ; ++dk)
    {
      for (std::int64_t dj = -reachY; dj <= reachY; ++dj)
      {
        for (std::int64_t di = -reachX; di <= reachX; ++di)
        {
          double3 fractional(static_cast<double>(di) / static_cast<double>(nx),
                             static_cast<double>(dj) / static_cast<double>(ny),
                             static_cast<double>(dk) / static_cast<double>(nz));
          double3 offset = cell * fractional;
          if (offset.length_squared() > reach * reach) continue;

          const std::int64_t ii = ((i + di) % nx + nx) % nx;
          const std::int64_t jj = ((j + dj) % ny + ny) % ny;
          const std::int64_t kk = ((k + dk) % nz + nz) % nz;

          // Spheres from different centres reach the same point, and which of them gets there first is not
          // fixed, so the largest has to be kept without assuming an order.
          std::atomic_ref<float> covered(poreRadius[static_cast<std::size_t>((kk * ny + jj) * nx + ii)]);
          float current = covered.load(std::memory_order_relaxed);
          while (current < radius &&
                 !covered.compare_exchange_weak(current, radius, std::memory_order_relaxed,
                                                std::memory_order_relaxed))
          {
          }
        }
      }
    }
  }

  return poreRadius;
}


PoreSizeCurve binPoreSizes(std::span<const float> poreRadius, std::span<const float> distance,
                           std::span<const std::uint8_t> accessible, std::span<const double> weight,
                           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  PoreSizeCurve curve;
  const std::size_t numberOfVoxels = distance.size();

  auto weightOf = [&](std::size_t v) { return weight.empty() ? 1.0 : weight[v]; };

  double weightedDiameter = 0.0;
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (distance[v] < 0.0f) continue;
    const bool reached = !accessible.empty() && accessible[v] != 0;
    const double w = weightOf(v);

    ++curve.numberOfVoidVoxels;
    curve.totalWeight += w;
    if (reached)
    {
      ++curve.numberOfAccessibleVoxels;
      curve.totalAccessibleWeight += w;
    }

    const double diameter = 2.0 * static_cast<double>(poreRadius[v]);
    curve.largestDiameter = std::max(curve.largestDiameter, diameter);
    weightedDiameter += w * diameter;
  }
  if (curve.totalWeight > 0.0) curve.meanDiameter = weightedDiameter / curve.totalWeight;

  const std::size_t bins = numberOfBins.value_or(200);
  const double largest = maximumDiameter.value_or(std::max(curve.largestDiameter, 1.0));
  const double step = largest / static_cast<double>(bins);
  curve.binWidth = step;

  std::vector<double> voidAtLeast(bins + 1, 0.0);
  std::vector<double> reachedAtLeast(bins + 1, 0.0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (distance[v] < 0.0f) continue;
    const double diameter = 2.0 * static_cast<double>(poreRadius[v]);
    std::size_t bin = static_cast<std::size_t>(std::floor(diameter / step));
    if (bin > bins) bin = bins;

    const double w = weightOf(v);
    voidAtLeast[bin] += w;
    if (!accessible.empty() && accessible[v] != 0) reachedAtLeast[bin] += w;
  }

  // Turn the weight per size into the weight at or above a size, from the top down.
  for (std::size_t bin = bins; bin-- > 0;)
  {
    voidAtLeast[bin] += voidAtLeast[bin + 1];
    reachedAtLeast[bin] += reachedAtLeast[bin + 1];
  }

  const double voidNormalisation = (curve.totalWeight > 0.0) ? curve.totalWeight : 1.0;
  const double reachedNormalisation = (curve.totalAccessibleWeight > 0.0) ? curve.totalAccessibleWeight : 1.0;

  curve.points.resize(bins);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    PoreSizeCurvePoint point;
    point.diameter = (static_cast<double>(bin) + 0.5) * step;
    point.cumulativeVoidFraction = voidAtLeast[bin] / voidNormalisation;
    point.cumulativeAccessibleFraction = reachedAtLeast[bin] / reachedNormalisation;
    point.distribution = (voidAtLeast[bin] - voidAtLeast[bin + 1]) / (voidNormalisation * step);
    point.accessibleDistribution = (reachedAtLeast[bin] - reachedAtLeast[bin + 1]) / (reachedNormalisation * step);
    curve.points[bin] = point;
  }

  return curve;
}
