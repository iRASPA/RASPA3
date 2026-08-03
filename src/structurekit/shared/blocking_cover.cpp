module;

module grid_blocking_cover;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import grid_connected_components;


GridBlockingCover coverPockets(uint3 gridSize, const UnitCell &unitCell,
                           std::span<const std::int32_t> voxelPore, const std::vector<GridPore> &pores,
                           std::span<const std::uint8_t> needsCover, std::span<const float> width,
                           std::span<const std::uint8_t> reachable)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridBlockingCover cover;

  const std::int32_t nx = static_cast<std::int32_t>(gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(gridSize.z);
  const std::size_t numberOfVoxels = voxelPore.size();
  const double3x3 cell = unitCell.cell;

  // The distance a step of the grid covers, so that a set of grid steps can be turned into a length.
  auto lengthOfSteps = [&cell, nx, ny, nz](int3 steps)
  {
    return (cell * double3(static_cast<double>(steps.x) / static_cast<double>(nx),
                           static_cast<double>(steps.y) / static_cast<double>(ny),
                           static_cast<double>(steps.z) / static_cast<double>(nz)))
        .length();
  };

  // Half a voxel diagonal: how far a point may be from the nearest grid point, and so how much further than
  // the last grid point of the pocket a sphere has to reach to have covered the pocket itself.
  const double halfDiagonal = 0.5 * lengthOfSteps(int3(1, 1, 1));

  double3 widths = unitCell.perpendicularWidths();

  auto voxelAt = [nx, ny, nz](std::int32_t i, std::int32_t j, std::int32_t k)
  {
    std::int32_t wi = ((i % nx) + nx) % nx;
    std::int32_t wj = ((j % ny) + ny) % ny;
    std::int32_t wk = ((k % nz) + nz) % nz;
    return static_cast<std::size_t>((wk * ny + wj) * nx + wi);
  };

  bool anyReachable = false;
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (!reachable.empty() && reachable[voxel] != 0)
    {
      anyReachable = true;
      break;
    }
  }

  // No infinities: the build turns them off.
  const double unreachable = std::numeric_limits<double>::max();

  // The least a single step of the grid can carry a point away from the plane of the other two axes, the
  // perpendicular width being the right measure of that in an oblique cell. A set of steps of which the
  // largest is n spans at least n of these, which is what lets the search below stop before it has looked
  // everywhere.
  const double shortestStep = std::min({widths.x / static_cast<double>(nx), widths.y / static_cast<double>(ny),
                                        widths.z / static_cast<double>(nz)});

  // How far it is from a grid point to the nearest reachable point, looked for no further than `furthest`
  // because a sphere is never grown further than that anyway.
  //
  // The looking is done in shells about the point and stops as soon as no shell still to come could hold
  // anything nearer than what has been found. A pocket that all but touches a channel is answered at once,
  // and the whole of the reach is only ever swept for a pocket with nothing reachable near it, which is the
  // case where the answer does not matter.
  auto distanceToReachable = [&](int3 centre, double furthest)
  {
    if (!anyReachable) return unreachable;

    auto look = [&](std::int32_t di, std::int32_t dj, std::int32_t dk, double &nearest)
    {
      if (reachable[voxelAt(centre.x + di, centre.y + dj, centre.z + dk)] == 0) return;
      nearest = std::min(nearest, lengthOfSteps(int3(di, dj, dk)));
    };

    const std::int32_t shells = static_cast<std::int32_t>(std::ceil(furthest / shortestStep));
    double nearest = unreachable;
    for (std::int32_t shell = 0; shell <= shells; ++shell)
    {
      if (nearest <= static_cast<double>(shell) * shortestStep) break;

      for (std::int32_t di = -shell; di <= shell; ++di)
      {
        for (std::int32_t dj = -shell; dj <= shell; ++dj)
        {
          if (std::abs(di) == shell || std::abs(dj) == shell)
          {
            for (std::int32_t dk = -shell; dk <= shell; ++dk) look(di, dj, dk, nearest);
          }
          else
          {
            look(di, dj, -shell, nearest);
            look(di, dj, shell, nearest);
          }
        }
      }
    }
    return nearest;
  };

  // Kept between pockets and cleared only where it was written, there being far more grid points than there
  // are points in any pocket.
  std::vector<std::uint8_t> seen(numberOfVoxels, 0);

  for (std::size_t poreIndex = 0; poreIndex < pores.size(); ++poreIndex)
  {
    if (poreIndex >= needsCover.size() || needsCover[poreIndex] == 0) continue;
    const GridPore &pore = pores[poreIndex];

    // The pocket is followed out of the cell, so that a pocket lying across a cell boundary is one connected
    // body of points rather than two halves at opposite ends of the cell. It closes on nothing, that being
    // what makes it a pocket, so following it can never contradict itself.
    std::vector<int3> points;
    points.reserve(pore.numberOfVoxels);

    std::vector<int3> stack;
    stack.push_back(pore.exampleVoxel);
    seen[voxelAt(pore.exampleVoxel.x, pore.exampleVoxel.y, pore.exampleVoxel.z)] = 1;
    while (!stack.empty())
    {
      int3 here = stack.back();
      stack.pop_back();
      points.push_back(here);

      for (std::int32_t dk = -1; dk <= 1; ++dk)
      {
        for (std::int32_t dj = -1; dj <= 1; ++dj)
        {
          for (std::int32_t di = -1; di <= 1; ++di)
          {
            if (di == 0 && dj == 0 && dk == 0) continue;
            int3 there(here.x + di, here.y + dj, here.z + dk);
            std::size_t voxel = voxelAt(there.x, there.y, there.z);
            if (seen[voxel] != 0) continue;
            if (voxelPore[voxel] != static_cast<std::int32_t>(poreIndex)) continue;
            seen[voxel] = 1;
            stack.push_back(there);
          }
        }
      }
    }
    cover.numberOfPocketVoxels += points.size();
    for (int3 point : points) seen[voxelAt(point.x, point.y, point.z)] = 0;

    // Widest point first, so that when a pocket does need more than one sphere the roomiest part of it is
    // dealt with while there is still the most left to take.
    std::ranges::sort(points, [&](int3 a, int3 b)
                      { return width[voxelAt(a.x, a.y, a.z)] > width[voxelAt(b.x, b.y, b.z)]; });

    std::vector<std::uint8_t> covered(points.size(), 0);
    std::size_t left = points.size();

    // What a sphere about a given point would be. It reaches the furthest point of the pocket still
    // uncovered, and half a voxel past it because a grid point stands for its voxel and not only for itself.
    auto sphereAbout = [&](int3 centre)
    {
      double wanted = halfDiagonal;
      for (std::size_t p = 0; p < points.size(); ++p)
      {
        if (covered[p] != 0) continue;
        wanted = std::max(wanted, halfDiagonal + lengthOfSteps(int3(points[p].x - centre.x, points[p].y - centre.y,
                                                                    points[p].z - centre.z)));
      }

      GridBlockingSphere sphere;
      double admissible = distanceToReachable(centre, wanted);
      sphere.radius = std::min(wanted, admissible);
      sphere.clipped = admissible < wanted;
      sphere.pocket = poreIndex;
      sphere.centreVoxel = centre;
      sphere.centreFractional =
          double3::fract(double3(static_cast<double>(centre.x) / static_cast<double>(nx),
                                 static_cast<double>(centre.y) / static_cast<double>(ny),
                                 static_cast<double>(centre.z) / static_cast<double>(nz)));
      return sphere;
    };

    while (left > 0)
    {
      // The widest point still uncovered. A sphere about it is the fallback, and covering it is what makes
      // each round of this get somewhere.
      std::size_t anchor = 0;
      while (covered[anchor] != 0) ++anchor;

      // The middle of what is left is the better place to put the sphere. From the widest point instead, a
      // pocket drawn out to one side is covered by a sphere reaching all the way across it, which blocks
      // nothing it should not but blocks a great deal more than it needs to and runs into channels that were
      // never in the way.
      double3 middle{0.0, 0.0, 0.0};
      for (std::size_t p = 0; p < points.size(); ++p)
      {
        if (covered[p] != 0) continue;
        middle += double3(static_cast<double>(points[p].x), static_cast<double>(points[p].y),
                          static_cast<double>(points[p].z));
      }
      middle /= static_cast<double>(left);
      int3 nearestToMiddle(static_cast<std::int32_t>(std::lround(middle.x)),
                           static_cast<std::int32_t>(std::lround(middle.y)),
                           static_cast<std::int32_t>(std::lround(middle.z)));

      // The middle is only taken if its sphere reaches the widest point, since otherwise there is no telling
      // that this round covers anything at all and the pocket would never be finished.
      GridBlockingSphere sphere = sphereAbout(nearestToMiddle);
      double toAnchor = lengthOfSteps(int3(points[anchor].x - nearestToMiddle.x, points[anchor].y - nearestToMiddle.y,
                                           points[anchor].z - nearestToMiddle.z));
      if (toAnchor > sphere.radius) sphere = sphereAbout(points[anchor]);

      for (std::size_t p = 0; p < points.size(); ++p)
      {
        if (covered[p] != 0) continue;
        double distance = lengthOfSteps(int3(points[p].x - sphere.centreVoxel.x, points[p].y - sphere.centreVoxel.y,
                                             points[p].z - sphere.centreVoxel.z));
        if (distance > sphere.radius && p != anchor) continue;
        covered[p] = 1;
        --left;
        ++sphere.voxelsCovered;
      }

      // A pocket point nearer to somewhere reachable than the grid can resolve leaves nothing to place a
      // sphere in. It is a neck the grid has closed rather than a pocket, and holds no volume, so it is
      // passed over and counted.
      if (sphere.radius < halfDiagonal)
      {
        cover.numberOfRefusedPoints += sphere.voxelsCovered;
        continue;
      }
      if (sphere.clipped) ++cover.numberOfClippedSpheres;
      cover.spheres.push_back(sphere);
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  cover.seconds = elapsed.count();

  return cover;
}
