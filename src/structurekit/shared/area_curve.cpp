module;

#ifdef _OPENMP
#include <omp.h>
#endif

module grid_area_curve;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;


double AreaCurve::areaAt(double level) const
{
  if (this->points.empty()) return 0.0;
  if (level < this->points.front().level || level > this->points.back().level) return 0.0;

  const double position = (level - this->points.front().level) / this->binWidth;
  const std::size_t below = static_cast<std::size_t>(std::floor(position));
  if (below + 1 >= this->points.size()) return this->points.back().area;

  const double fraction = position - static_cast<double>(below);
  return (1.0 - fraction) * this->points[below].area + fraction * this->points[below + 1].area;
}


AreaCurve areaAgainstLevel(uint3 gridSize, const UnitCell &unitCell, std::span<const float> field,
                           double lowestLevel, double highestLevel, std::size_t numberOfBins)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  AreaCurve curve;
  curve.lowestLevel = lowestLevel;
  curve.highestLevel = highestLevel;

  const std::int64_t nx = static_cast<std::int64_t>(gridSize.x);
  const std::int64_t ny = static_cast<std::int64_t>(gridSize.y);
  const std::int64_t nz = static_cast<std::int64_t>(gridSize.z);
  const std::size_t numberOfVoxels = field.size();
  if (numberOfVoxels == 0 || numberOfBins == 0 || highestLevel <= lowestLevel) return curve;

  const std::size_t bins = numberOfBins;
  const double step = (highestLevel - lowestLevel) / static_cast<double>(bins);
  curve.binWidth = step;

  // A gradient taken on the grid is a gradient with respect to the fractional coordinate, and the field is
  // sampled at intervals of 1/n along each of them. Turning it into a gradient with respect to position needs
  // the inverse cell transposed, since s = H^-1 x gives d/dx = (H^-1)^T d/ds, and on an oblique cell that is
  // not the same as dividing each component by its own spacing.
  const double3x3 inverseCell = unitCell.inverseCell;
  const double3 fractionalStep(1.0 / static_cast<double>(nx), 1.0 / static_cast<double>(ny),
                               1.0 / static_cast<double>(nz));
  const double voxelVolume = unitCell.volume / static_cast<double>(numberOfVoxels);

  // A voxel whose value lies outside the range may still have a band reaching into it, and dropping those
  // starves the bins at either end of the curve, by as much as a third in the last one. They are let through
  // and their bands clipped. The guard is wide enough that no voxel with anything to contribute is turned
  // away, and narrow enough to go on rejecting the wall of an energy field, whose values are ten orders of
  // magnitude out and whose gradients would otherwise be smeared across the whole curve.
  const double guard = highestLevel - lowestLevel;

  std::vector<double> gathered(bins, 0.0);
  std::vector<double> occupied(bins, 0.0);
  double underRange = 0.0;
  double overRange = 0.0;
  double totalSpan = 0.0;
  std::size_t spanned = 0;

#pragma omp parallel
  {
    std::vector<double> localGathered(bins, 0.0);
    std::vector<double> localOccupied(bins, 0.0);
    double localUnder = 0.0;
    double localOver = 0.0;
    double localSpan = 0.0;
    std::size_t localSpanned = 0;

#pragma omp for schedule(static)
    for (std::int64_t voxel = 0; voxel < static_cast<std::int64_t>(numberOfVoxels); ++voxel)
    {
      const double here = static_cast<double>(field[static_cast<std::size_t>(voxel)]);

      if (here < lowestLevel - guard)
      {
        localUnder += 1.0;
        continue;
      }
      if (here >= highestLevel + guard)
      {
        localOver += 1.0;
        continue;
      }

      const std::int64_t k = voxel / (nx * ny);
      const std::int64_t j = (voxel / nx) % ny;
      const std::int64_t i = voxel % nx;

      auto at = [&](std::int64_t a, std::int64_t b, std::int64_t c)
      {
        const std::int64_t aa = ((a % nx) + nx) % nx;
        const std::int64_t bb = ((b % ny) + ny) % ny;
        const std::int64_t cc = ((c % nz) + nz) % nz;
        return static_cast<double>(field[static_cast<std::size_t>((cc * ny + bb) * nx + aa)]);
      };

      const double alongX = 0.5 * (at(i + 1, j, k) - at(i - 1, j, k));
      const double alongY = 0.5 * (at(i, j + 1, k) - at(i, j - 1, k));
      const double alongZ = 0.5 * (at(i, j, k + 1) - at(i, j, k - 1));

      double3 fractionalGradient(alongX / fractionalStep.x, alongY / fractionalStep.y, alongZ / fractionalStep.z);
      const double magnitude = transposedMultiply(inverseCell, fractionalGradient).length();

      // How far the field moves across this one voxel, which is the projection of the voxel onto the gradient
      // and works out to the sum of the three differences already in hand.
      const double span = std::abs(alongX) + std::abs(alongY) + std::abs(alongZ);
      if (here >= lowestLevel && here < highestLevel)
      {
        localSpan += span;
        ++localSpanned;
      }

      // The voxel is not a point. Its value is a sample of a field that varies across it, and dropping the
      // whole of its weight into the single bin its center happens to fall in is what makes this estimator
      // alias: bins narrower than `span` are then filled by whether a grid point landed in them rather than
      // by the field, and neighbouring bins come out at wildly different heights. Taking the field as linear
      // across the voxel, which is the same approximation the gradient already is, spreads that weight evenly
      // over the values the voxel actually covers. Nothing is lost, since the integral over the band is what
      // it was, and the curve becomes smooth in the only sense available to it: the field cannot resolve
      // structure finer than it moves in one voxel, and it now says so rather than inventing some.
      const double lowest = here - 0.5 * span;
      const double highest = here + 0.5 * span;

      if (span <= 0.0)
      {
        if (here < lowestLevel)
        {
          localUnder += 1.0;
        }
        else if (here >= highestLevel)
        {
          localOver += 1.0;
        }
        else
        {
          std::size_t bin = static_cast<std::size_t>((here - lowestLevel) / step);
          if (bin >= bins) bin = bins - 1;
          localGathered[bin] += magnitude;
          localOccupied[bin] += 1.0;
        }
        continue;
      }

      // Whatever of the band falls outside the range is not lost, it is counted as lying below or above.
      localUnder += std::clamp(lowestLevel - lowest, 0.0, span) / span;
      localOver += std::clamp(highest - highestLevel, 0.0, span) / span;

      const std::int64_t first = std::max<std::int64_t>(
          0, static_cast<std::int64_t>(std::floor((lowest - lowestLevel) / step)));
      const std::int64_t last = std::min<std::int64_t>(
          static_cast<std::int64_t>(bins) - 1, static_cast<std::int64_t>(std::floor((highest - lowestLevel) / step)));

      for (std::int64_t bin = first; bin <= last; ++bin)
      {
        const double binLow = lowestLevel + static_cast<double>(bin) * step;
        const double overlap = std::min(highest, binLow + step) - std::max(lowest, binLow);
        if (overlap <= 0.0) continue;

        const double share = overlap / span;
        localGathered[bin] += magnitude * share;
        localOccupied[bin] += share;
      }
    }

#pragma omp critical
    {
      for (std::size_t bin = 0; bin < bins; ++bin)
      {
        gathered[bin] += localGathered[bin];
        occupied[bin] += localOccupied[bin];
      }
      underRange += localUnder;
      overRange += localOver;
      totalSpan += localSpan;
      spanned += localSpanned;
    }
  }

  const double total = static_cast<double>(numberOfVoxels);
  curve.fractionBelowRange = underRange / total;
  curve.fractionAboveRange = overRange / total;
  curve.meanSpanPerVoxel = (spanned > 0) ? totalSpan / static_cast<double>(spanned) : 0.0;

  // Everything under the range is under every level of it, so the running share starts there rather than at
  // nothing.
  double running = underRange;

  curve.points.resize(bins);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    running += occupied[bin];
    curve.points[bin].level = lowestLevel + (static_cast<double>(bin) + 0.5) * step;
    curve.points[bin].area = gathered[bin] * voxelVolume / step;
    curve.points[bin].fractionBelow = running / total;
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  curve.seconds = elapsed.count();

  return curve;
}
