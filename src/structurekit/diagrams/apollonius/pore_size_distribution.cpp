module;

module apollonius_pore_size_distribution;

import std;

import double3;
import crystal;
import pair_interactions;
import apollonius_accessibility;
import exact_pore_size_distribution;
import exact_void_split;
import voronoi_blocking_spheres;
import voronoi_pore_size_distribution;

void ApolloniusPoreSizeDistribution::run(const PairInteractions& interactions, const Crystal& framework,
                                         std::string probePseudoAtom, std::optional<double> maximumDiameter,
                                         std::optional<std::size_t> numberOfBins, std::size_t subdivisions,
                                         double floorRadius, bool peaksOnly)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusPoreSizeDistribution: Unknown probe-atom type\n");
  }
  const double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.atoms.size());
  radii.reserve(framework.atoms.size());
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  auto build = [&](double inflation)
  {
    return ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, inflation)
        .accessibility;
  };

  const double cellVolume = framework.unitCell.volume;
  const double maxDiameter = maximumDiameter.value_or(20.0);
  const std::size_t bins = numberOfBins.value_or(100);

  this->sweptSpheres.clear();
  this->pocketsFilled = false;

  // The same division of the routes as the radical diagram's: see there.
  if (floorRadius >= probeRadius)
  {
    curve = peaksOnly ? exactPoreSizePeaks(build, cellVolume, subdivisions, probeRadius, floorRadius)
                      : exactPoreSizeDistribution(build, cellVolume, maxDiameter, bins, subdivisions, probeRadius,
                                                  floorRadius);
  }
  else
  {
    ApolloniusAccessibility classifier =
        ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius);
    ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, cellVolume, subdivisions);
    std::vector<BlockingSphere> blockedBalls =
        hybridBlockedBalls(classifier.accessibility, split, fractionalPositions, radii, probeRadius,
                           static_cast<std::size_t>(200.0 * cellVolume));

    this->sweptSpheres = blockedBalls;

    // The atoms come back in another order, widest first, so the two lists are compared sorted alike.
    auto sorted = [](std::vector<std::array<double, 4>> spheres)
    {
      std::ranges::sort(spheres);
      return spheres;
    };
    std::vector<std::array<double, 4>> ofAtoms;
    std::vector<std::array<double, 4>> ofBalls;
    ofAtoms.reserve(radii.size());
    ofBalls.reserve(blockedBalls.size());
    for (std::size_t i = 0; i < radii.size(); ++i)
    {
      ofAtoms.push_back({radii[i], fractionalPositions[i].x, fractionalPositions[i].y, fractionalPositions[i].z});
    }
    for (const BlockingSphere& ball : blockedBalls)
    {
      ofBalls.push_back({ball.radius, ball.centerFractional.x, ball.centerFractional.y, ball.centerFractional.z});
    }
    this->pocketsFilled = sorted(std::move(ofAtoms)) != sorted(std::move(ofBalls));

    std::vector<double3> blockedPositions;
    std::vector<double> blockedRadii;
    blockedPositions.reserve(blockedBalls.size());
    blockedRadii.reserve(blockedBalls.size());
    for (const BlockingSphere& ball : blockedBalls)
    {
      blockedPositions.push_back(ball.centerFractional);
      blockedRadii.push_back(ball.radius);
    }

    auto buildBlocked = [&](double inflation)
    {
      return ApolloniusAccessibility::create(framework.unitCell, blockedPositions, blockedRadii, inflation)
          .accessibility;
    };

    PoreSizeDistributionCurve blocked =
        peaksOnly ? exactPoreSizePeaks(buildBlocked, cellVolume, subdivisions, 0.0, floorRadius)
                  : exactPoreSizeDistribution(buildBlocked, cellVolume, maxDiameter, bins, subdivisions, 0.0,
                                              floorRadius);
    curve = hybridFromBlockedCurve(std::move(blocked), probeRadius, floorRadius);
  }

  if (peaksOnly)
    writePoreSizePeaks(framework, "apollonius", probePseudoAtom, curve);
  else
    writePoreSizeDistribution(framework, "apollonius", probePseudoAtom, curve);
}
