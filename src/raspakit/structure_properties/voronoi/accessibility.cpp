module;

module voronoi_accessibility;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
import voronoi_network;
import voronoi_channels;

// Splits a possibly out-of-range bin coordinate into a wrapped bin index and the periodic
// image it came from.
std::pair<int, int> binAndImage(int coordinate, int gridExtent)
{
  int image = (coordinate >= 0) ? coordinate / gridExtent : -((-coordinate + gridExtent - 1) / gridExtent);
  return {coordinate - image * gridExtent, image};
}

int3 binOfFractional(const double3& fractional, const int3& gridSize)
{
  return int3(std::min(gridSize.x - 1, static_cast<int>(fractional.x * static_cast<double>(gridSize.x))),
              std::min(gridSize.y - 1, static_cast<int>(fractional.y * static_cast<double>(gridSize.y))),
              std::min(gridSize.z - 1, static_cast<int>(fractional.z * static_cast<double>(gridSize.z))));
}

VoronoiAccessibility VoronoiAccessibility::create(const SimulationBox& simulationBox,
                                                  const std::vector<double3>& fractionalPositions,
                                                  const std::vector<double>& radii, double probeRadius)
{
  VoronoiAccessibility accessibility;
  accessibility.simulationBox = simulationBox;

  // Inflate atoms by the probe radius; the network's node radii then measure the room
  // available to the probe's centre.
  std::vector<double> inflatedRadii(radii.size());
  for (std::size_t i = 0; i < radii.size(); ++i) inflatedRadii[i] = radii[i] + probeRadius;
  accessibility.atomRadii = inflatedRadii;
  accessibility.maximumAtomRadius = *std::max_element(inflatedRadii.begin(), inflatedRadii.end());

  accessibility.network = VoronoiNetwork::create(simulationBox, fractionalPositions, inflatedRadii);

  // With inflated atoms, any node with positive radius has room for the probe centre, so
  // channels are detected at probe radius zero.
  accessibility.channels = ChannelAnalysis::compute(accessibility.network, 0.0);

  accessibility.nodeAccessible.assign(accessibility.network.nodes.size(), 0);
  for (std::size_t i = 0; i < accessibility.network.nodes.size(); ++i)
  {
    std::int32_t poreId = accessibility.channels.nodePoreId[i];
    if (poreId >= 0 && accessibility.channels.pores[static_cast<std::size_t>(poreId)].isChannel)
    {
      accessibility.nodeAccessible[i] = 1;
    }
  }

  accessibility.atomPositions.reserve(accessibility.network.atomPositionsFractional.size());
  for (const double3& fractional : accessibility.network.atomPositionsFractional)
  {
    accessibility.atomPositions.push_back(simulationBox.cell * fractional);
  }

  // Cell list over the fractional unit cube (same construction as SKVoronoi): roughly
  // four atoms per bin, bin counts proportional to the perpendicular widths so that the
  // bins are approximately metrically cubic.
  const double3x3 cell = simulationBox.cell;
  double volume = simulationBox.volume;
  double3 a = double3(cell[0].x, cell[0].y, cell[0].z);
  double3 b = double3(cell[1].x, cell[1].y, cell[1].z);
  double3 c = double3(cell[2].x, cell[2].y, cell[2].z);
  double3 perpendicularWidths = double3(volume / double3::cross(b, c).length(),
                                        volume / double3::cross(c, a).length(),
                                        volume / double3::cross(a, b).length());
  std::size_t numberOfAtoms = accessibility.atomPositions.size();
  double targetBinSize = std::cbrt(volume / std::max(1.0, static_cast<double>(numberOfAtoms) / 4.0));
  accessibility.gridSize = int3(std::max(1, static_cast<int>(perpendicularWidths.x / targetBinSize)),
                                std::max(1, static_cast<int>(perpendicularWidths.y / targetBinSize)),
                                std::max(1, static_cast<int>(perpendicularWidths.z / targetBinSize)));
  accessibility.minimumBinWidth =
      std::min({perpendicularWidths.x / static_cast<double>(accessibility.gridSize.x),
                perpendicularWidths.y / static_cast<double>(accessibility.gridSize.y),
                perpendicularWidths.z / static_cast<double>(accessibility.gridSize.z)});

  const int3 gridSize = accessibility.gridSize;
  accessibility.bins.assign(static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) *
                                static_cast<std::size_t>(gridSize.z),
                            {});
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    int3 bin = binOfFractional(accessibility.network.atomPositionsFractional[i], gridSize);
    accessibility.bins[static_cast<std::size_t>((bin.z * gridSize.y + bin.y) * gridSize.x + bin.x)].push_back(i);
  }

  return accessibility;
}

PointClassification VoronoiAccessibility::classify(const double3& point) const
{
  PointClassification classification;

  // Work with the wrapped point; the configuration is periodic so the classification of
  // the wrapped point equals that of the original.
  double3 fractional = double3::fract(simulationBox.inverseCell * point);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  // Nearest atom via the cell list: walk bin shells outward (wrapping periodically); all
  // atoms in shell k are at least (k-1)·(minimum bin width) away, so the walk stops once
  // that bound exceeds the best distance found.
  std::size_t nearestAtom = 0;
  double nearestDistanceSquared = std::numeric_limits<double>::max();
  double3 nearestDelta(0.0, 0.0, 0.0);
  bool found = false;
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (k > 0 && found && lowerBound * lowerBound > nearestDistanceSquared) break;

    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(pointBin.x + ox, gridSize.x);
          auto [by, ly] = binAndImage(pointBin.y + oy, gridSize.y);
          auto [bz, lz] = binAndImage(pointBin.z + oz, gridSize.z);

          double3 imageShift =
              simulationBox.cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz)) -
              wrappedPoint;

          for (std::size_t j : bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
          {
            double3 delta = atomPositions[j] + imageShift;
            double distanceSquared = double3::dot(delta, delta);
            if (distanceSquared < nearestDistanceSquared)
            {
              nearestDistanceSquared = distanceSquared;
              nearestAtom = j;
              nearestDelta = delta;
              found = true;
            }
          }
        }
      }
    }
  }
  double3 nearestAtomImage = wrappedPoint + nearestDelta;  // nearest periodic image of the atom

  double nearestDistance = std::sqrt(nearestDistanceSquared);
  if (nearestDistance < atomRadii[nearestAtom] - 1.0e-8)
  {
    classification.inside = true;
    return classification;
  }

  // Line-of-sight test against the Voronoi nodes of the nearest atom's cell.
  double3 sampleRay = wrappedPoint - nearestAtomImage;
  double bestDistanceSquared = std::numeric_limits<double>::max();
  bool decided = false;
  for (const auto& [nodeIndex, vertexRelative] : network.atomNodeVectors[nearestAtom])
  {
    double3 nodePosition = nearestAtomImage + vertexRelative;
    double3 otherRay = wrappedPoint - nodePosition;
    if (double3::dot(sampleRay, otherRay) <= 0.0)
    {
      double distanceSquared = double3::dot(otherRay, otherRay);
      if (distanceSquared < bestDistanceSquared)
      {
        bestDistanceSquared = distanceSquared;
        classification.accessible = nodeAccessible[nodeIndex] != 0;
        classification.poreId = channels.nodePoreId[nodeIndex];
        decided = true;
      }
    }
  }

  if (!decided) classification.resample = true;
  return classification;
}

bool VoronoiAccessibility::overlapsAtom(const double3& point, std::size_t excludedAtom) const
{
  double3 fractional = double3::fract(simulationBox.inverseCell * point);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  // Only atoms within the largest inflated radius can contain the point.
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (lowerBound > maximumAtomRadius) break;

    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(pointBin.x + ox, gridSize.x);
          auto [by, ly] = binAndImage(pointBin.y + oy, gridSize.y);
          auto [bz, lz] = binAndImage(pointBin.z + oz, gridSize.z);

          double3 imageShift =
              simulationBox.cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz)) -
              wrappedPoint;

          for (std::size_t j : bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
          {
            if (j == excludedAtom) continue;
            double3 delta = atomPositions[j] + imageShift;
            if (double3::dot(delta, delta) < atomRadii[j] * atomRadii[j]) return true;
          }
        }
      }
    }
  }
  return false;
}
