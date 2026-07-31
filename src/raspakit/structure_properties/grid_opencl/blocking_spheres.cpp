module;

module opencl_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;


GridBlockingSpheres::GridBlockingSpheres() {}


GridBlockingSpheres::~GridBlockingSpheres() {}


void GridBlockingSpheres::run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
                              uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(forceField, framework, gridSize);
  this->run(forceField, framework, probePseudoAtom, grid);
}


void GridBlockingSpheres::run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
                              const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridBlockingSpheres: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();
  this->gridSeconds = grid.seconds;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridComponents components = GridComponents::compute(grid, this->probeRadius);

  const std::int32_t nx = static_cast<std::int32_t>(grid.gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(grid.gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(grid.gridSize.z);
  const std::size_t numberOfVoxels = grid.numberOfVoxels();
  const double3x3 cell = framework.simulationBox.cell;

  this->numberOfPockets = components.numberOfPockets;
  this->numberOfSinglePointPockets = components.numberOfSinglePointPockets;

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

  double3 widths = framework.simulationBox.perpendicularWidths();

  auto voxelAt = [nx, ny, nz](std::int32_t i, std::int32_t j, std::int32_t k)
  {
    std::int32_t wi = ((i % nx) + nx) % nx;
    std::int32_t wj = ((j % ny) + ny) % ny;
    std::int32_t wk = ((k % nz) + nz) % nz;
    return static_cast<std::size_t>((wk * ny + wj) * nx + wi);
  };

  // Where the probe both fits and can get out, which is what a sphere must not reach into.
  std::vector<std::uint8_t> isChannel(numberOfVoxels, 0);
  bool anyChannel = false;
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    std::int32_t pore = components.voxelPore[voxel];
    if (pore < 0) continue;
    if (components.pores[static_cast<std::size_t>(pore)].isChannel)
    {
      isChannel[voxel] = 1;
      anyChannel = true;
    }
  }

  // The nearest point of a channel to a given grid point, looked for no further than a given distance because
  // a sphere is never grown further than that anyway. Doing it here, once per sphere that is actually placed,
  // rather than once per grid point of every pocket, is what keeps it cheap: there are a handful of spheres and
  // there are a great many grid points.
  // No infinities: the build turns them off.
  const double unreachable = std::numeric_limits<double>::max();

  // The least a single step of the grid can carry a point away from the plane of the other two axes, the
  // perpendicular width being the right measure of that in an oblique cell. A set of steps of which the largest
  // is n spans at least n of these, which is what lets the search below stop before it has looked everywhere.
  const double shortestStep =
      std::min({widths.x / static_cast<double>(nx), widths.y / static_cast<double>(ny),
                widths.z / static_cast<double>(nz)});

  // How far it is from a grid point to the nearest point of a channel, looked for no further than `furthest`
  // because a sphere is never grown further than that anyway.
  //
  // The looking is done in shells about the point and stops as soon as no shell still to come could hold
  // anything nearer than what has been found. A pocket that all but touches a channel is answered at once, and
  // the whole of the reach is only ever swept for a pocket with no channel anywhere near it, which is the case
  // where the answer does not matter.
  auto distanceToChannel = [&](int3 centre, double furthest)
  {
    if (!anyChannel) return unreachable;

    auto look = [&](std::int32_t di, std::int32_t dj, std::int32_t dk, double &nearest)
    {
      if (isChannel[voxelAt(centre.x + di, centre.y + dj, centre.z + dk)] == 0) return;
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

  const double voxelVolume = grid.voxelVolume();
  std::size_t pocketVoxels = 0;

  // Kept between pockets and cleared only where it was written, there being far more grid points than there are
  // points in any pocket.
  std::vector<std::uint8_t> seen(numberOfVoxels, 0);

  for (std::size_t poreIndex = 0; poreIndex < components.pores.size(); ++poreIndex)
  {
    const GridPore &pore = components.pores[poreIndex];
    if (pore.isChannel) continue;

    // The pocket is followed out of the cell, so that a pocket lying across a cell boundary is one connected
    // body of points rather than two halves at opposite ends of the cell. It closes on nothing, that being what
    // makes it a pocket, so following it can never contradict itself.
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
            if (components.voxelPore[voxel] != static_cast<std::int32_t>(poreIndex)) continue;
            seen[voxel] = 1;
            stack.push_back(there);
          }
        }
      }
    }
    pocketVoxels += points.size();
    for (int3 point : points) seen[voxelAt(point.x, point.y, point.z)] = 0;

    // Widest point first, so that when a pocket does need more than one sphere the roomiest part of it is dealt
    // with while there is still the most left to take.
    std::ranges::sort(points, [&](int3 a, int3 b)
                      { return grid.clearance[voxelAt(a.x, a.y, a.z)] > grid.clearance[voxelAt(b.x, b.y, b.z)]; });

    std::vector<std::uint8_t> covered(points.size(), 0);
    std::size_t left = points.size();

    // What a sphere about a given point would be. It reaches the furthest point of the pocket still uncovered,
    // and half a voxel past it because a grid point stands for its voxel and not only for itself.
    //
    // Reaching that far takes in framework, and other pockets, and void too narrow for the probe. None of that
    // costs anything: a simulation compares this sphere against the centre of an atom it is trying to place, and
    // nowhere in any of it was a centre going to be allowed. What does cost something is reaching a channel, so
    // the sphere is stopped short of the nearest point the probe can both occupy and leave.
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
      double admissible = distanceToChannel(centre, wanted);
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
      // The widest point still uncovered. A sphere about it is the fallback, and covering it is what makes each
      // round of this get somewhere.
      std::size_t anchor = 0;
      while (covered[anchor] != 0) ++anchor;

      // The middle of what is left is the better place to put the sphere, and is where the Voronoi and
      // Apollonius routes put theirs. From the widest point instead, a pocket drawn out to one side is covered
      // by a sphere reaching all the way across it, which blocks nothing it should not but blocks a great deal
      // more than it needs to and runs into channels that were never in the way.
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

      // A pocket point nearer to a channel than the grid can resolve leaves nothing to place a sphere in. It is
      // a neck the grid has closed rather than a pocket, and holds no volume, so it is passed over and counted.
      if (sphere.radius < halfDiagonal)
      {
        this->numberOfRefusedPoints += sphere.voxelsCovered;
        continue;
      }
      if (sphere.clipped) ++this->numberOfClippedSpheres;
      this->spheres.push_back(sphere);
    }
  }

  this->numberOfPocketVoxels = pocketVoxels;
  this->pocketFraction = static_cast<double>(pocketVoxels) / static_cast<double>(numberOfVoxels);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  // What a simulation reads: nothing but the numbers, a comment being more than the format allows.
  std::ofstream blockFile;
  blockFile.open(framework.name + ".block");
  std::print(blockFile, "{}\n", this->spheres.size());
  for (const GridBlockingSphere &sphere : this->spheres)
  {
    std::print(blockFile, "{} {} {} {}\n", sphere.centreFractional.x, sphere.centreFractional.y,
               sphere.centreFractional.z, sphere.radius);
  }
  blockFile.close();

  double3 spacing = grid.spacing();

  std::ofstream report;
  report.open(framework.name + ".grid.block.gpu.txt");
  std::print(report, "# Blocking spheres (clearance grid)\n");
  std::print(report, "# Framework: {}\n", framework.name);
  std::print(report, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, this->probeRadius);
  std::print(report, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(report, "# Channels: {}, pockets: {}, of those {} holding a single grid point\n",
             components.numberOfChannels, this->numberOfPockets, this->numberOfSinglePointPockets);
  std::print(report, "# Grid points in the pockets: {} of {}, {} [Å³], a fraction {} of the cell\n", pocketVoxels,
             numberOfVoxels, static_cast<double>(pocketVoxels) * voxelVolume, this->pocketFraction);
  std::print(report, "# Spheres written to {}.block: {}\n", framework.name, this->spheres.size());
  std::print(report, "# Of those, {} were stopped by a channel rather than by the pocket's own extent\n",
             this->numberOfClippedSpheres);
  std::print(report, "# {} points of a pocket had no sphere placed on them, lying nearer a channel than the grid\n",
             this->numberOfRefusedPoints);
  std::print(report, "# can resolve: a neck the grid has closed rather than a pocket, holding no volume\n");
  std::print(report, "# GPU Timing: {} [s] for the clearance field\n", this->gridSeconds);
  std::print(report, "# CPU Timing: {} [s] for the pockets and the spheres\n", this->seconds);
  std::print(report, "#\n");
  std::print(report, "# A pocket is a piece of the region where the probe fits that does not run anywhere, so what\n");
  std::print(report, "# is blocked here is decided by the shape of that region and not by sampling it. Each pocket\n");
  std::print(report, "# is covered from its widest point outwards, and a sphere stops at the nearest point the\n");
  std::print(report, "# probe can actually reach: a sphere past that would block pore volume the probe is owed,\n");
  std::print(report, "# and a simulation given such a file would report an uptake too low for a reason nothing in\n");
  std::print(report, "# its output would explain.\n");
  std::print(report, "#\n");
  std::print(report, "# The grid decides where the necks are. A passage narrower than the spacing reads as closed,\n");
  std::print(report, "# which turns a piece of channel into a pocket and blocks it; a finer grid opens it again.\n");
  std::print(report, "# So the count of pockets is the part of this to check against a finer grid, the spheres\n");
  std::print(report, "# themselves being no more than a covering of whatever the pockets turned out to be.\n");
  std::print(report, "#\n");
  std::print(report, "# column 1-3: centre, fractional\n");
  std::print(report, "# column 4: radius [Å]\n");
  std::print(report, "# column 5: which pocket\n");
  std::print(report, "# column 6: grid points of that pocket this sphere covered\n");
  std::print(report, "# column 7: what stopped it growing\n");
  std::print(report, "#       s_a         s_b         s_c   radius [Å]  pocket    points  stopped by\n");
  for (const GridBlockingSphere &sphere : this->spheres)
  {
    std::print(report, "Sphere: {:11.7f} {:11.7f} {:11.7f} {:11.5f} {:7} {:9}  {}\n", sphere.centreFractional.x,
               sphere.centreFractional.y, sphere.centreFractional.z, sphere.radius, sphere.pocket,
               sphere.voxelsCovered,
               sphere.clipped ? "a channel" : "the pocket");
  }
  report.close();
}
