module;

module mc_blocking_pockets;

import std;

import double3;
import unit_cell;
import voronoi_network;
import blocking_spheres;
import sampled_structure;
import sampled_roadmap;
import mc_backend;

namespace
{
double periodicDistance(const UnitCell &unitCell, const double3 &a, const double3 &b)
{
  double3 dr = unitCell.applyPeriodicBoundaryConditions(a - b);
  return dr.length();
}
}  // namespace

void MC_BlockingPockets::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                             std::optional<std::size_t> numberOfInnerSteps)
{
  run(structure, SampledRoadmap::build(structure, samplingBackendCPU(), numberOfIterations, numberOfInnerSteps));
}

void MC_BlockingPockets::run(const SampledStructure &structure, const SampledRoadmap &roadmap)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->roadmap = roadmap;
  this->spheres.clear();

  const UnitCell &unitCell = structure.unitCell;
  const VoronoiNetwork &network = roadmap.network;

  this->numberOfPockets = roadmap.numberOfPockets;
  this->numberOfUnresolvedPieces = roadmap.numberOfUnresolvedPieces;

  // The points of each piece of the void, sorted into the ones a probe can reach from outside and the ones
  // it cannot. Every node counts here, the pocket centres as well as the samples: a cavity only just wide
  // enough to hold the probe's centre has almost no volume for a point to land in, and where one did land
  // the walk uphill from it has already found the deepest point of the cavity, which is the best centre a
  // sphere over it could have.
  //
  // A piece the sample could not resolve is left alone: nothing is written over it, and nothing about it is
  // taken as reachable either, since either would be an assertion the sample cannot support and one of the
  // two would put a sphere in a working pore.
  struct PocketPoint
  {
    double3 position;
    double clearance;
  };

  std::vector<double3> reachable;
  std::map<std::int32_t, std::vector<PocketPoint>> shutIn;

  for (std::size_t node = 0; node < network.nodes.size(); ++node)
  {
    if (roadmap.isReachable(node))
      reachable.push_back(network.nodes[node].position);
    else if (roadmap.isShutIn(node))
      shutIn[roadmap.components.nodePoreId[node]].push_back(
          {network.nodes[node].position, network.nodes[node].radius});
  }

  // How much of the void the pockets hold, counted over the sample points alone: the pocket centres stand
  // at places the sample already covers and counting them would count that volume twice.
  std::size_t volumeNodes = roadmap.numberOfVolumeNodes();
  std::size_t blockedSamples = 0;
  for (std::size_t node = 0; node < volumeNodes; ++node)
  {
    if (roadmap.isShutIn(node)) ++blockedSamples;
  }
  this->blockedFractionOfVoid =
      volumeNodes > 0 ? static_cast<double>(blockedSamples) / static_cast<double>(volumeNodes) : 0.0;

  // The finest detail the sample can resolve. A sphere reaches this far past the outermost point it covers,
  // so that the point is inside it rather than on it, and a sphere thinner than this is not written at all,
  // it would block less than the volume of a single sample.
  this->resolution =
      roadmap.numberOfSamples > 0
          ? std::cbrt(unitCell.volume / static_cast<double>(roadmap.numberOfSamples))
          : 0.0;

  // How far the distance to the nearest reachable point can overstate the distance to the reachable region
  // itself, which is the spacing of the reachable points alone. It is held off the cap below.
  double channelMargin =
      reachable.empty() ? 0.0 : std::cbrt(unitCell.volume / static_cast<double>(reachable.size()));

  for (const auto &[poreId, points] : shutIn)
  {
    std::vector<PocketPoint> pocket = points;

    // Deepest point first, so each pocket is covered from the inside out and the widest sphere is placed
    // while there is still the most left to cover.
    std::sort(pocket.begin(), pocket.end(),
              [](const PocketPoint &a, const PocketPoint &b) { return a.clearance > b.clearance; });

    std::vector<char> covered(pocket.size(), 0);
    bool first = true;
    bool wrote = false;

    for (std::size_t i = 0; i < pocket.size(); ++i)
    {
      if (covered[i] != 0) continue;
      const double3 center = pocket[i].position;

      double furthestPocket = 0.0;
      for (std::size_t j = i; j < pocket.size(); ++j)
      {
        if (covered[j] == 0)
          furthestPocket = std::max(furthestPocket, periodicDistance(unitCell, center, pocket[j].position));
      }

      double closestChannel = std::numeric_limits<double>::max();
      for (const double3 &point : reachable)
        closestChannel = std::min(closestChannel, periodicDistance(unitCell, center, point));

      // Reaching the outermost point left to cover is all that is wanted, and there is nothing to gain by
      // going further.
      double wanted = furthestPocket + this->resolution;

      // A sphere may not reach a position the probe can occupy from a channel, or the simulation loses part
      // of its pore rather than a pocket. Nothing else constrains it: swallowing the atoms around the
      // pocket is harmless, and so is swallowing a neighbouring pocket, since no molecule belongs in
      // either. That is why a pocket in a structure with no channel at all gets one sphere over the whole
      // of it.
      double admissible =
          reachable.empty() ? std::numeric_limits<double>::max() : closestChannel - channelMargin;

      // Whatever those two say, a sphere of the clearance at the centre is always safe, by the argument in
      // the interface. That floor is what keeps a pocket from being abandoned when the sample puts a
      // reachable point improbably close to its centre.
      double radius = std::max(pocket[i].clearance, std::min(wanted, admissible));

      // Near the rim of a pocket there is no room for a sphere worth writing, and the classification there
      // is the least certain anyway. Leaving those points uncovered is what bounds the count: the clearance
      // falls to zero at the boundary, so insisting on covering the last of it would take unboundedly many
      // ever smaller spheres. What is given up are positions where the probe is already touching an atom.
      // The first sphere of a pocket is exempt, a pocket narrow throughout being still a pocket.
      if (radius < this->resolution && !first)
      {
        covered[i] = 1;
        continue;
      }
      first = false;

      this->spheres.push_back(
          BlockingSphere{double3::fract(unitCell.inverseCell * center), radius});
      wrote = true;

      for (std::size_t j = i; j < pocket.size(); ++j)
      {
        if (covered[j] == 0 && periodicDistance(unitCell, center, pocket[j].position) < radius) covered[j] = 1;
      }
    }

    if (wrote) ++this->numberOfCoveredPockets;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  this->seconds = timing.count();

  writeBlockFile(structure.name, this->spheres);

  std::ofstream report;
  report.open(std::format("{}.mc.block.{}.txt", structure.name, roadmap.backend));
  std::print(report, "# Blocking spheres by sampling the void\n");
  structure.writeHeader(report);
  roadmap.writeHeader(report);
  std::print(report, "# Timing (blocking, on the processor either way): {} [s]\n", this->seconds);
  std::print(report, "# Spheres written to {}.block: {}\n", structure.name, this->spheres.size());
  std::print(report, "# Sample resolution: {} [Å]\n", this->resolution);
  std::print(report, "# Every point below is a position for the probe's centre, so the radii need no further\n");
  std::print(report, "# allowance for the probe's size. A sphere reaches as far as the furthest point of its\n");
  std::print(report, "# pocket left uncovered and no nearer a channel than the sample can tell, and never\n");
  std::print(report, "# less than the clearance at its own centre, which is inside the pocket by\n");
  std::print(report, "# construction.\n");
  std::print(report, "#\n");
  std::print(report, "# A pocket no point landed in is a pocket left open, which is what this route can do\n");
  std::print(report, "# that the surface routes cannot. The counts below are what says whether the sample was\n");
  std::print(report, "# equal to the structure.\n");
  std::print(report, "#     s_a         s_b         s_c      radius [Å]\n");
  for (const BlockingSphere &sphere : this->spheres)
  {
    std::print(report, "Sphere: {:11.7f} {:11.7f} {:11.7f} {:11.5f}\n", sphere.centerFractional.x,
               sphere.centerFractional.y, sphere.centerFractional.z, sphere.radius);
  }
  std::print(report, "Number of pockets: {}\n", this->numberOfPockets);
  std::print(report, "Number of pockets covered: {}\n", this->numberOfCoveredPockets);
  std::print(report, "Number of pieces left open as too thinly sampled: {}\n", this->numberOfUnresolvedPieces);
  std::print(report, "Share of the sampled void shut in a pocket: {}\n", this->blockedFractionOfVoid);
  report.close();
}
