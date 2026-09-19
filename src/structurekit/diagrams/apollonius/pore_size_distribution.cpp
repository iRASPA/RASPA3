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

  if (floorRadius > 0.0)
  {
    curve = peaksOnly ? exactPoreSizePeaks(build, cellVolume, subdivisions, probeRadius, floorRadius)
                      : exactPoreSizeDistribution(build, cellVolume, maxDiameter, bins, subdivisions, probeRadius,
                                                  floorRadius);
  }
  else
  {
    // Hybrid: fill helium-sealed pockets with their blocking spheres, then take the bare PSD of what remains.
    ApolloniusAccessibility classifier =
        ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius);
    ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, cellVolume, subdivisions);
    std::vector<BlockingSphere> spheres;
    if (measuredSpheresRefused(split).empty())
    {
      spheres = exactBlockingSpheres(split);
    }
    else
    {
      spheres = computeBlockingSpheres(classifier.accessibility, static_cast<std::size_t>(200.0 * cellVolume));
    }

    std::vector<double3> blockedPositions = fractionalPositions;
    std::vector<double> blockedRadii = radii;
    blockedPositions.reserve(fractionalPositions.size() + spheres.size());
    blockedRadii.reserve(radii.size() + spheres.size());
    for (const BlockingSphere& sphere : spheres)
    {
      blockedPositions.push_back(sphere.centerFractional);
      blockedRadii.push_back(sphere.radius);
    }

    auto buildBlocked = [&](double inflation)
    {
      return ApolloniusAccessibility::create(framework.unitCell, blockedPositions, blockedRadii, inflation)
          .accessibility;
    };

    PoreSizeDistributionCurve blocked =
        peaksOnly ? exactPoreSizePeaks(buildBlocked, cellVolume, subdivisions, 0.0, 0.0)
                  : exactPoreSizeDistribution(buildBlocked, cellVolume, maxDiameter, bins, subdivisions, 0.0, 0.0);
    curve = hybridFromBlockedCurve(std::move(blocked), probeRadius);
  }

  if (peaksOnly)
    writePoreSizePeaks(framework, "apollonius", probePseudoAtom, curve);
  else
    writePoreSizeDistribution(framework, "apollonius", probePseudoAtom, curve);
}
