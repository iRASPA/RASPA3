module;

module pore_accessibility;

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

PoreAccessibility PoreAccessibility::create(const SimulationBox& simulationBox,
                                                  const std::vector<double3>& fractionalPositions,
                                                  const std::vector<double>& radii, double probeRadius)
{
  // Inflate atoms by the probe radius; the network's node radii then measure the room
  // available to the probe's centre.
  std::vector<double> inflatedRadii(radii.size());
  for (std::size_t i = 0; i < radii.size(); ++i) inflatedRadii[i] = radii[i] + probeRadius;

  return createFromNetwork(VoronoiNetwork::create(simulationBox, fractionalPositions, inflatedRadii), Metric::Power,
                           probeRadius);
}

PoreAccessibility PoreAccessibility::createFromNetwork(VoronoiNetwork network, Metric metric, double probeRadius)
{
  PoreAccessibility accessibility;
  const SimulationBox simulationBox = network.simulationBox;
  accessibility.simulationBox = simulationBox;
  accessibility.metric = metric;
  accessibility.probeRadius = probeRadius;
  accessibility.atomRadii = network.atomRadii;
  accessibility.maximumAtomRadius = *std::max_element(network.atomRadii.begin(), network.atomRadii.end());
  accessibility.network = std::move(network);

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

  // Needs the cell list, so it comes last.
  accessibility.nodeClearance.reserve(accessibility.network.nodes.size());
  for (const VoronoiNode& node : accessibility.network.nodes)
  {
    accessibility.nodeClearance.push_back(accessibility.clearance(node.position));
  }

  // The nodes into the same bins, for the containment test. Nodes too narrow for the probe are left out:
  // they carry no pore and so can say nothing about a point that lands in one of them.
  accessibility.nodeBins.assign(accessibility.bins.size(), {});
  accessibility.nodeBinMaximumClearance.assign(accessibility.bins.size(), 0.0);
  for (std::size_t i = 0; i < accessibility.network.nodes.size(); ++i)
  {
    if (accessibility.channels.nodePoreId[i] < 0) continue;
    double3 fractional = double3::fract(simulationBox.inverseCell * accessibility.network.nodes[i].position);
    int3 bin = binOfFractional(fractional, gridSize);
    std::size_t index = static_cast<std::size_t>((bin.z * gridSize.y + bin.y) * gridSize.x + bin.x);
    accessibility.nodeBins[index].push_back(i);
    accessibility.nodeBinMaximumClearance[index] =
        std::max(accessibility.nodeBinMaximumClearance[index], accessibility.nodeClearance[i]);
    accessibility.maximumNodeClearance = std::max(accessibility.maximumNodeClearance, accessibility.nodeClearance[i]);
  }

  return accessibility;
}

PointClassification PoreAccessibility::classify(const double3& point) const
{
  PointClassification classification;

  // Work with the wrapped point; the configuration is periodic so the classification of
  // the wrapped point equals that of the original. The translation that wrapping removed is kept, since
  // the frame the answer refers to is the caller's point and not the wrapped one.
  double3 unwrappedFractional = simulationBox.inverseCell * point;
  const int3 wrap(static_cast<int>(std::floor(unwrappedFractional.x)),
                  static_cast<int>(std::floor(unwrappedFractional.y)),
                  static_cast<int>(std::floor(unwrappedFractional.z)));
  double3 fractional = double3::fract(unwrappedFractional);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  // The pore of the node that decides the point, and the lift of that pore the point sits next to. The
  // lift is taken from the nearest node of the pore rather than from the deciding node's own image: the
  // deciding node is chosen to settle which pore, and the tests below that do that by proximity or by
  // reach may well pick a node in another copy of it, which is no error in an answer about pores but is
  // the wrong frame to bring the point into.
  auto decideFrom = [&](std::size_t nodeIndex)
  {
    classification.accessible = nodeAccessible[nodeIndex] != 0;
    classification.poreId = channels.nodePoreId[nodeIndex];
    classification.latticeOffset = -nearestLift(wrappedPoint, classification.poreId) - wrap;
  };

  // The cell containing the point is the one nearest in the metric the network's cells are cut by:
  // the power distance |x-x_i|^2 - r_i^2 for a radical diagram, the clearance |x-x_i| - r_i for an
  // Apollonius diagram. Picking the Euclidean-nearest atom instead can land on a site whose cell is
  // empty, which either diagram allows, leaving no vertices for the line-of-sight test below.
  //
  // Search via the cell list: walk bin shells outward (wrapping periodically); all atoms in
  // shell k are at least (k-1)·(minimum bin width) away, so their distance in either metric is at
  // least what that bound gives, and the walk stops once that lower bound exceeds the best found.
  std::size_t nearestAtom = 0;
  double nearestDistance = std::numeric_limits<double>::max();
  double3 nearestDelta(0.0, 0.0, 0.0);
  bool found = false;
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    double reachable = (metric == Metric::Power) ? lowerBound * lowerBound - maximumAtomRadius * maximumAtomRadius
                                                 : lowerBound - maximumAtomRadius;
    if (k > 0 && found && lowerBound > 0.0 && reachable > nearestDistance) break;

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
            double distance = (metric == Metric::Power) ? double3::dot(delta, delta) - atomRadii[j] * atomRadii[j]
                                                        : delta.length() - atomRadii[j];
            if (distance < nearestDistance)
            {
              nearestDistance = distance;
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

  // A negative distance in either metric means the point lies inside that inflated atom.
  if (nearestDistance < -1.0e-8)
  {
    classification.inside = true;
    return classification;
  }

  // A node whose free ball holds the point settles the question outright, so that is asked first. It is
  // the only test here that proves rather than infers, and asking it second was how the Apollonius route
  // came to write blocking spheres in the middle of a channel in YFI: a point with over an \si{\angstrom}
  // of room to spare, held by any number of channel nodes, was handed to a pocket because the diagram's
  // account of which nodes surround the nearest atom is incomplete on a strongly degenerate structure.
  if (std::optional<std::size_t> holder = containingNode(wrappedPoint); holder.has_value())
  {
    decideFrom(holder.value());
    return classification;
  }

  // Line-of-sight test against the Voronoi nodes of the nearest atom's cell, for the points no ball
  // holds: those nearer to an atom than any node has room, which is most of a surface-area sweep.
  double3 sampleRay = wrappedPoint - nearestAtomImage;
  double bestDistanceSquared = std::numeric_limits<double>::max();
  bool decided = false;
  std::size_t bestNode = 0;
  for (const auto& [nodeIndex, vertexRelative] : network.atomNodeVectors[nearestAtom])
  {
    double3 nodePosition = nearestAtomImage + vertexRelative;
    double3 otherRay = wrappedPoint - nodePosition;
    double distanceSquared = double3::dot(otherRay, otherRay);

    if (double3::dot(sampleRay, otherRay) <= 0.0 && distanceSquared < bestDistanceSquared)
    {
      bestDistanceSquared = distanceSquared;
      bestNode = nodeIndex;
      decided = true;
    }
  }
  if (decided)
  {
    decideFrom(bestNode);
    return classification;
  }

  // In a compact cage every node of the cell can sit forward of the sample point, and then nothing
  // passes the test above. What settles it then is the node whose own free ball comes nearest to holding
  // the point: distance to the node, less the room for the probe's centre there. A point inside such a
  // ball is provably in the same channel or pocket as its node, the ball being free and connected, so
  // the measure is negative exactly when the answer is certain and otherwise says by how much it falls
  // short. Every node is a candidate, not just those of the nearest atom's cell, because on a strongly
  // degenerate diagram the cell of an atom lining a cage need not carry any of that cage's nodes -- and
  // then a search confined to the cell can only return a node in the channel next door.
  //
  // What this replaces is the nearest node of the cell, taken with no regard for how much room it has.
  // That guess went wrong in one direction: the nearest node of an atom in the wall between a cage and a
  // channel is as likely to be on the channel side as not, so a sealed cage was reported as reachable
  // and its pocket never blocked. It is the whole reason the Apollonius route lost LTA's sodalite cage.
  // Drawing another point instead is no answer either: the points that cannot be decided are the ones
  // inside cages, so replacing them costs a cage the volume it should have been counted for.
  //
  // This runs over every node and only points that get this far pay for it, which is few enough that it
  // does not show in the cost of a sweep.
  double bestReach = std::numeric_limits<double>::max();
  std::size_t reachNode = 0;
  for (std::size_t nodeIndex = 0; nodeIndex < network.nodes.size(); ++nodeIndex)
  {
    if (channels.nodePoreId[nodeIndex] < 0) continue;  // too narrow for the probe, so it says nothing
    double3 delta = simulationBox.applyPeriodicBoundaryConditions(wrappedPoint - network.nodes[nodeIndex].position);
    double reach = delta.length() - nodeClearance[nodeIndex];
    if (reach < bestReach)
    {
      bestReach = reach;
      reachNode = nodeIndex;
      decided = true;
    }
  }

  if (decided)
  {
    decideFrom(reachNode);
    return classification;
  }

  // No node anywhere has room for the probe, so there is nothing here to decide it with.
  classification.resample = true;
  return classification;
}

std::optional<std::size_t> PoreAccessibility::containingNode(const double3& point, bool accessibleOnly) const
{
  double3 fractional = double3::fract(simulationBox.inverseCell * point);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  std::optional<std::size_t> best;
  double deepest = 0.0;
  for (int k = 0;; ++k)
  {
    // Nothing in shell k is nearer than this, so once even the largest free ball anywhere falls short of
    // it there is nothing left to find.
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (k > 0 && lowerBound >= maximumNodeClearance) break;

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
          std::size_t index = static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx);
          if (nodeBinMaximumClearance[index] <= lowerBound) continue;

          double3 imageShift =
              simulationBox.cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz)) -
              wrappedPoint;

          for (std::size_t nodeIndex : nodeBins[index])
          {
            if (accessibleOnly && nodeAccessible[nodeIndex] == 0) continue;
            double3 delta = network.nodes[nodeIndex].position + imageShift;
            double reach = delta.length() - nodeClearance[nodeIndex];
            if (reach < deepest)
            {
              deepest = reach;
              best = nodeIndex;
            }
          }
        }
      }
    }
  }
  return best;
}

int3 PoreAccessibility::nearestLift(const double3& point, std::int32_t poreId) const
{
  if (poreId < 0) return int3(0, 0, 0);
  const VoronoiPore& pore = channels.pores[static_cast<std::size_t>(poreId)];

  int3 best(0, 0, 0);
  double nearest = std::numeric_limits<double>::max();
  for (std::size_t nodeIndex : pore.nodeIndices)
  {
    // Where this node sits in the lift the pore is assembled in, and the translation carrying that to the
    // copy of it nearest the point, which is the minimum image and needs no search.
    const int3& offset = channels.nodeLatticeOffset[nodeIndex];
    double3 assembled = network.nodes[nodeIndex].position +
                        simulationBox.cell * double3(static_cast<double>(offset.x), static_cast<double>(offset.y),
                                                     static_cast<double>(offset.z));
    double3 fractional = simulationBox.inverseCell * (point - assembled);
    int3 translation(static_cast<int>(std::lround(fractional.x)), static_cast<int>(std::lround(fractional.y)),
                     static_cast<int>(std::lround(fractional.z)));
    double3 shift = simulationBox.cell * double3(static_cast<double>(translation.x),
                                                static_cast<double>(translation.y),
                                                static_cast<double>(translation.z));
    double distance = (point - (assembled + shift)).length();
    if (distance < nearest)
    {
      nearest = distance;
      best = translation;
    }
  }
  return best;
}

bool PoreAccessibility::provablyAccessible(const double3& point) const
{
  return containingNode(point, true).has_value();
}

double PoreAccessibility::clearance(const double3& point) const
{
  double3 fractional = double3::fract(simulationBox.inverseCell * point);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  double best = std::numeric_limits<double>::max();
  bool found = false;
  for (int k = 0;; ++k)
  {
    // Atoms in shell k are at least (k-1)·(minimum bin width) away, so the smallest clearance they can
    // contribute is that less the largest radius; stop once even that cannot beat what is in hand.
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (k > 0 && found && lowerBound - maximumAtomRadius > best) break;

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
            best = std::min(best, delta.length() - atomRadii[j]);
            found = true;
          }
        }
      }
    }
  }
  return found ? best : std::numeric_limits<double>::max();
}

std::vector<std::pair<double3, double>> PoreAccessibility::neighbourAtoms(const double3& point, double reach) const
{
  double3 fractional = double3::fract(simulationBox.inverseCell * point);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  std::vector<std::pair<double3, double>> neighbours;
  for (int k = 0;; ++k)
  {
    // Nothing in shell k is nearer to the point than this, so once that passes the reach there is
    // nothing left to collect. The walk is over the shells of the same cell list the other queries
    // use, and wraps, so a reach larger than the cell returns the images as separate entries.
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (k > 0 && lowerBound > reach) break;

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
            if (delta.length() <= reach) neighbours.emplace_back(delta, atomRadii[j]);
          }
        }
      }
    }
  }
  return neighbours;
}


std::vector<NeighbourImage> PoreAccessibility::neighbourAtomImages(const double3& point, double reach) const
{
  // The same walk as `neighbourAtoms`, keeping what that one drops: which atom each image is of, and
  // which image. The bin walk knows both -- the atom by its index in the bin, the image by how far the
  // shell reached outside the cell -- and the only correction is the wrap the query point itself went
  // through on the way in, since the translations are wanted relative to the point as given.
  double3 unwrappedFractional = simulationBox.inverseCell * point;
  const int3 wrap(static_cast<int>(std::floor(unwrappedFractional.x)),
                  static_cast<int>(std::floor(unwrappedFractional.y)),
                  static_cast<int>(std::floor(unwrappedFractional.z)));
  double3 fractional = double3::fract(unwrappedFractional);
  double3 wrappedPoint = simulationBox.cell * fractional;
  const int3 pointBin = binOfFractional(fractional, gridSize);

  std::vector<NeighbourImage> neighbours;
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
    if (k > 0 && lowerBound > reach) break;

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
            if (delta.length() > reach) continue;

            NeighbourImage neighbour;
            neighbour.delta = delta;
            neighbour.radius = atomRadii[j];
            neighbour.index = j;
            neighbour.image = int3(lx + wrap.x, ly + wrap.y, lz + wrap.z);
            neighbours.push_back(neighbour);
          }
        }
      }
    }
  }
  return neighbours;
}

bool PoreAccessibility::overlapsAtom(const double3& point, std::size_t excludedAtom) const
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
