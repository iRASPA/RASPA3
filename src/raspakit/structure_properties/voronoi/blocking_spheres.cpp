module;

module voronoi_blocking_spheres;

import std;

import double3;
import simulationbox;
import randomnumbers;
import atom;
import framework;
import forcefield;
import pore_accessibility;
import exact_void_split;

std::vector<BlockingSphere> exactBlockingSpheres(const ExactVoidSplit& split)
{
  std::vector<BlockingSphere> spheres;
  spheres.reserve(split.pockets.size());
  for (const PocketGeometry& pocket : split.pockets)
  {
    if (pocket.blockingRadius() <= 0.0) continue;
    spheres.push_back(BlockingSphere{pocket.centreFractional, pocket.blockingRadius()});
  }
  return spheres;
}

std::string measuredSpheresRefused(const ExactVoidSplit& split)
{
  if (!split.reliable) return split.rejection;

  // Both radii are distances from the centroid, so they say what they say about the region the centroid is in.
  // Where that is not the pocket but the void beyond it, the sphere would be centred in a channel, and there is
  // no radius small enough to make that safe.
  std::size_t adrift = static_cast<std::size_t>(
      std::ranges::count_if(split.pockets, [](const PocketGeometry& pocket) { return pocket.centreInChannel; }));
  if (adrift > 0)
  {
    return std::format("{} of {} pockets have their centroid in a channel", adrift, split.pockets.size());
  }
  return {};
}

double periodicDistance(const SimulationBox& simulationBox, const double3& a, const double3& b)
{
  double3 delta = simulationBox.applyPeriodicBoundaryConditions(a - b);
  return delta.length();
}

std::vector<BlockingSphere> computeBlockingSpheres(const PoreAccessibility& accessibility,
                                                   std::size_t numberOfSamples)
{
  RandomNumber random{samplingSeed};
  const SimulationBox& simulationBox = accessibility.simulationBox;
  std::vector<BlockingSphere> spheres;

  std::size_t samples = std::max<std::size_t>(1, numberOfSamples);

  // Sample the void: points a probe can reach on one side, points it cannot grouped by pocket on the
  // other. Each pocket point carries the room left for the probe's centre there, which sets both the
  // order they are covered in and a radius that is safe whatever the sampling.
  struct PocketPoint
  {
    double3 position;
    double clearance;
  };
  std::vector<double3> accessiblePoints;
  std::map<std::int32_t, std::vector<PocketPoint>> pocketPoints;

  // Every pocket starts from the diagram's own account of it, its nodes, so that no pocket depends on
  // the sampling to be discovered at all. A cavity only just wide enough to hold the probe's centre has
  // almost no volume for a sample to land in -- MON's eight cages leave 0.04 Å of room and would be
  // missed every time -- but the diagram has a node sitting in each one.
  for (std::size_t pore = 0; pore < accessibility.channels.pores.size(); ++pore)
  {
    if (accessibility.channels.pores[pore].isChannel) continue;
    for (std::size_t node : accessibility.channels.pores[pore].nodeIndices)
    {
      pocketPoints[static_cast<std::int32_t>(pore)].push_back(
          {accessibility.network.nodes[node].position, accessibility.nodeClearance[node]});
    }
  }

  for (std::size_t s = 0; s < samples; ++s)
  {
    double3 point;
    PointClassification classification;
    for (std::size_t attempt = 0; attempt < resampleLimit; ++attempt)
    {
      point = simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      classification = accessibility.classify(point);
      if (!classification.resample) break;
    }
    if (classification.inside || classification.resample) continue;
    if (classification.accessible)
      accessiblePoints.push_back(point);
    else if (classification.poreId >= 0)
      pocketPoints[classification.poreId].push_back({point, accessibility.clearance(point)});
  }

  // Every sample here is a position for the probe's *centre*, so distances between them are already in
  // those terms and a radius read off them needs no further allowance for the probe's size. What does
  // need allowing for is the resolution of the sampling, and these are the two lengths that set it.
  //
  // The first is the mean spacing of the samples, the finest detail the sampling can resolve. A sphere
  // reaches this far past the outermost point it covers, so that the point is inside it rather than on
  // it, and a sphere thinner than this is not written at all: it would block less than the volume of a
  // single sample.
  double resolution = std::cbrt(simulationBox.volume / static_cast<double>(samples));

  // The second is the mean spacing of the accessible samples alone, which is how far the distance to
  // the nearest of them can overstate the distance to the accessible region itself. It is held off the
  // cap below. Both lengths fall away as the sampling densifies, which is what keeps the radii a
  // property of the structure rather than of the number of samples drawn.
  double channelMargin = accessiblePoints.empty()
                             ? 0.0
                             : std::cbrt(simulationBox.volume / static_cast<double>(accessiblePoints.size()));

  for (auto& [poreId, points] : pocketPoints)
  {
    // Deepest point first, so each pocket is covered from the inside out and the widest sphere is placed
    // while there is still the most left to cover.
    std::sort(points.begin(), points.end(),
              [](const PocketPoint& a, const PocketPoint& b) { return a.clearance > b.clearance; });

    std::vector<char> covered(points.size(), 0);
    bool first = true;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      if (covered[i] != 0) continue;
      const double3 center = points[i].position;

      // Refuse the sphere if this point turns out to be in a channel after all. The pocket it was
      // assigned to came from the diagram's connectivity, and where that is wrong -- a strongly
      // degenerate diagram can leave a stretch of channel joined to nothing and so counted as a pocket
      // -- a sphere here would block part of the pore. The free ball of an accessible node holding the
      // point settles it against the diagram, since the ball is free space and connected. Only points
      // that would otherwise become sphere centres are checked, which is a few hundred at most.
      if (accessibility.provablyAccessible(center))
      {
        covered[i] = 1;
        continue;
      }

      double furthestPocket = 0.0;
      for (std::size_t j = i; j < points.size(); ++j)
      {
        if (covered[j] == 0)
          furthestPocket = std::max(furthestPocket, periodicDistance(simulationBox, center, points[j].position));
      }

      double closestChannel = std::numeric_limits<double>::max();
      for (const double3& p : accessiblePoints)
        closestChannel = std::min(closestChannel, periodicDistance(simulationBox, center, p));

      // Reaching the outermost point left to cover is all that is wanted, and there is nothing to gain
      // by going further.
      double wanted = furthestPocket + resolution;

      // A sphere may not reach a position the probe can occupy from a channel, or the simulation loses
      // part of its pore rather than a pocket. Nothing else constrains it: swallowing the atoms around
      // the pocket is harmless, and so is swallowing a neighbouring pocket, since no molecule belongs in
      // either. That is why a pocket with no channel anywhere in the structure gets a single sphere over
      // the whole of it.
      double admissible = accessiblePoints.empty() ? std::numeric_limits<double>::max()
                                                   : closestChannel - channelMargin;

      // Whatever the two above say, a sphere of the clearance at the centre is always safe: it holds no
      // atom, so it lies in the free space, and being connected it cannot cross the neck that made this
      // a pocket, so it is inside the pocket by construction. That floor is what keeps a pocket from
      // being abandoned when the sampling puts an accessible point improbably close to its centre.
      double radius = std::max(points[i].clearance, std::min(wanted, admissible));

      // Near the rim of a pocket there is no room for a sphere worth writing, and the classification
      // there is the least certain anyway. Leaving those points uncovered is what bounds the count: the
      // clearance falls to zero at the boundary, so insisting on covering the last of it would take
      // unboundedly many ever smaller spheres. What is given up are positions where the probe is already
      // touching an atom. The first sphere of a pocket is exempt, since a pocket narrow throughout is
      // still a pocket and closing it is the whole point.
      if (radius < resolution && !first)
      {
        covered[i] = 1;
        continue;
      }
      first = false;

      spheres.push_back(BlockingSphere{double3::fract(simulationBox.inverseCell * center), radius});

      for (std::size_t j = i; j < points.size(); ++j)
      {
        if (covered[j] == 0 && periodicDistance(simulationBox, center, points[j].position) < radius) covered[j] = 1;
      }
    }
  }

  return spheres;
}

void writeBlockingSpheres(const std::string& frameworkName, const std::string& diagramName,
                          const std::string& probeName, double probeRadius,
                          const std::vector<BlockingSphere>& spheres, const std::vector<PocketGeometry>& pockets,
                          const std::string& fallbackReason)
{
  // What a simulation reads: nothing but the numbers, a comment being more than the format allows.
  std::ofstream blockFile;
  blockFile.open(frameworkName + ".block");
  std::print(blockFile, "{}\n", spheres.size());
  for (const BlockingSphere& sphere : spheres)
  {
    std::print(blockFile, "{} {} {} {}\n", sphere.centerFractional.x, sphere.centerFractional.y,
               sphere.centerFractional.z, sphere.radius);
  }
  blockFile.close();

  std::ofstream report;
  report.open(std::format("{}.{}.block.txt", frameworkName, diagramName));
  std::print(report, "# Blocking spheres ({}, {})\n", diagramName,
             fallbackReason.empty() ? "measured" : "sampled");
  std::print(report, "# Framework: {}\n", frameworkName);
  std::print(report, "# Probe atom: {} radius: {} [Å]\n", probeName, probeRadius);
  std::print(report, "# Spheres written to {}.block: {}\n", frameworkName, spheres.size());
  if (!fallbackReason.empty())
  {
    std::print(report, "# The measured pockets were refused ({}); the spheres are sampled\n", fallbackReason);
  }
  else
  {
    std::print(report, "# One sphere per pocket, at the centroid of the pocket, of the lesser of two radii:\n");
    std::print(report, "#   `covering`, how far the pocket's own wall gets from that centre, which holds the\n");
    std::print(report, "#     whole pocket and nothing outside its walls; and\n");
    std::print(report, "#   `channel`, how near the boundary of the accessible void comes, past which a sphere\n");
    std::print(report, "#     would block part of a pore. `none` where the structure has no channel at all.\n");
    std::print(report, "# Both are extrema over the same patches of the same surfaces, in closed form, so the\n");
    std::print(report, "# radius is neither a sampling estimate nor a bound: it is attained.\n");
    std::print(report, "#     s_a         s_b         s_c      volume [Å³]  covering [Å]  channel [Å]  radius [Å]\n");
    for (const PocketGeometry& pocket : pockets)
    {
      std::print(report, "Sphere: {:11.7f} {:11.7f} {:11.7f} {:11.5f} {:9.4f} {:>12} {:9.4f}{}\n",
                 pocket.centreFractional.x, pocket.centreFractional.y, pocket.centreFractional.z, pocket.volume,
                 pocket.coveringRadius, pocket.hasChannel() ? std::format("{:9.4f}", pocket.channelRadius) : "none",
                 pocket.blockingRadius(), pocket.coversPocket() ? "" : "  (capped short of the pocket)");
    }
  }
  report.close();
}


void VoronoiBlockingSpheres::run(const ForceField& forceField, const Framework& framework,
                                 std::string probePseudoAtom, std::optional<std::size_t> numberOfSamples)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiBlockingSpheres: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  PoreAccessibility accessibility =
      PoreAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);

  double volume = framework.simulationBox.volume;

  // The pockets come from the surface around them, which also says how far each may be blocked. Sampling is
  // what is left where that has to be refused: for the reasons the void fraction's division is refused, which
  // are surface the decomposition could not place and boundaries that did not close, or for the one reason
  // particular to a sphere, a centre that turns out to lie in a channel.
  ExactVoidSplit split = exactVoidSplitByComponents(accessibility, volume);
  fallbackReason = measuredSpheresRefused(split);
  if (fallbackReason.empty())
  {
    measured = true;
    pockets = split.pockets;
    spheres = exactBlockingSpheres(split);
  }
  else
  {
    std::size_t samples = numberOfSamples.value_or(static_cast<std::size_t>(200.0 * volume));
    spheres = computeBlockingSpheres(accessibility, samples);
  }

  writeBlockingSpheres(framework.name, "voronoi", probePseudoAtom, probeRadius, spheres, pockets, fallbackReason);
}
