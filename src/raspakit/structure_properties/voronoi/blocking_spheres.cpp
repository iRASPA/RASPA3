module;

module voronoi_blocking_spheres;

import std;

import double3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import voronoi_accessibility;

constexpr double overshoot = 0.1;  // Å; zeo++ sphere_radius_overshoot

double periodicDistance(const Framework& framework, const double3& a, const double3& b)
{
  double3 delta = framework.simulationBox.applyPeriodicBoundaryConditions(a - b);
  return delta.length();
}

// Index of the point with the highest local density (Gaussian-weighted neighbour count),
// evaluated on a capped subsample to keep the cost bounded.
std::size_t mostDenseIndex(const Framework& framework, const std::vector<double3>& points)
{
  std::size_t count = std::min<std::size_t>(points.size(), 1000);
  if (count <= 1) return 0;

  double meanDistance = 0.0;
  std::size_t pairs = 0;
  for (std::size_t i = 0; i < count; ++i)
  {
    for (std::size_t j = i + 1; j < count; ++j)
    {
      meanDistance += periodicDistance(framework, points[i], points[j]);
      ++pairs;
    }
  }
  meanDistance /= static_cast<double>(std::max<std::size_t>(1, pairs));
  double sigmaSquared = std::max(1.0e-6, meanDistance * meanDistance);

  std::size_t best = 0;
  double bestDensity = -1.0;
  for (std::size_t i = 0; i < count; ++i)
  {
    double density = 0.0;
    for (std::size_t j = 0; j < count; ++j)
    {
      if (i == j) continue;
      double d = periodicDistance(framework, points[i], points[j]);
      density += std::exp(-d * d / sigmaSquared);
    }
    if (density > bestDensity)
    {
      bestDensity = density;
      best = i;
    }
  }
  return best;
}

void VoronoiBlockingSpheres::run(const ForceField& forceField, const Framework& framework,
                                 std::string probePseudoAtom, std::optional<std::size_t> numberOfSamples)
{
  RandomNumber random{std::nullopt};

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

  VoronoiAccessibility accessibility =
      VoronoiAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);

  double volume = framework.simulationBox.volume;
  std::size_t samples = numberOfSamples.value_or(static_cast<std::size_t>(200.0 * volume));
  samples = std::max<std::size_t>(1, samples);

  // Sample the void; separate accessible points from inaccessible points grouped by pocket.
  std::vector<double3> accessiblePoints;
  std::map<std::int32_t, std::vector<double3>> pocketPoints;
  for (std::size_t s = 0; s < samples; ++s)
  {
    double3 point = framework.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
    PointClassification classification = accessibility.classify(point);
    if (classification.inside || classification.resample) continue;
    if (classification.accessible)
      accessiblePoints.push_back(point);
    else if (classification.poreId >= 0)
      pocketPoints[classification.poreId].push_back(point);
  }

  // Greedily cover each pocket with spheres.
  for (auto& [poreId, points] : pocketPoints)
  {
    std::vector<double3> remaining = points;
    while (!remaining.empty())
    {
      std::size_t centerIndex = mostDenseIndex(framework, remaining);
      double3 center = remaining[centerIndex];

      double furthestPocket = 0.0;
      for (const double3& p : remaining) furthestPocket = std::max(furthestPocket, periodicDistance(framework, center, p));

      double closestChannel = std::numeric_limits<double>::max();
      for (const double3& p : accessiblePoints)
        closestChannel = std::min(closestChannel, periodicDistance(framework, center, p));

      double radius;
      if (accessiblePoints.empty())
      {
        radius = furthestPocket + probeRadius + overshoot;
      }
      else if (furthestPocket < closestChannel)
      {
        radius = std::min(furthestPocket + probeRadius + overshoot, closestChannel - (probeRadius + overshoot));
      }
      else
      {
        radius = std::max(overshoot, closestChannel - (probeRadius + overshoot));
      }

      double3 centerFractional = double3::fract(framework.simulationBox.inverseCell * center);
      spheres.push_back(BlockingSphere{centerFractional, radius});

      std::vector<double3> survivors;
      survivors.reserve(remaining.size());
      for (const double3& p : remaining)
      {
        if (periodicDistance(framework, center, p) >= radius) survivors.push_back(p);
      }
      // Guard against non-termination if the sphere covers nothing.
      if (survivors.size() == remaining.size()) break;
      remaining = std::move(survivors);
    }
  }

  std::ofstream myfile;
  myfile.open(framework.name + ".block");
  std::print(myfile, "{}\n", spheres.size());
  for (const BlockingSphere& sphere : spheres)
  {
    std::print(myfile, "{} {} {} {}\n", sphere.centerFractional.x, sphere.centerFractional.y,
               sphere.centerFractional.z, sphere.radius);
  }
  myfile.close();
}
