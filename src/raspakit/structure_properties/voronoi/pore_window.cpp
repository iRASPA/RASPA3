module;

module pore_window;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
import voronoi_network;
import voronoi_channels;
import voronoi_pore_diameters;

// A disc: what one atom leaves of the plane of the window. Its centre is in the plane's own two
// coordinates, relative to the bottleneck, and `gap` is the squared distance from the bottleneck to the
// disc's rim along the line through its centre, which is positive exactly when the bottleneck is outside
// the disc and is what the ray cast below needs.
struct WindowDisc
{
  double first;
  double second;
  double radius;
  double gap;
};

// The perpendicular widths of the cell: the distance between the pair of faces normal to each axis.
// A displacement of d along any direction crosses at most ceil(d / width) cells of the corresponding
// index, whichever way the cell is sheared.
static double3 perpendicularWidthsOf(const SimulationBox& simulationBox)
{
  const double3x3 cell = simulationBox.cell;
  double3 a(cell[0].x, cell[0].y, cell[0].z);
  double3 b(cell[1].x, cell[1].y, cell[1].z);
  double3 c(cell[2].x, cell[2].y, cell[2].z);
  double volume = simulationBox.volume;
  return double3(volume / double3::cross(b, c).length(), volume / double3::cross(c, a).length(),
                 volume / double3::cross(a, b).length());
}

PoreWindow PoreWindow::measure(const VoronoiNetwork& network, const double3& position, const double3& normal)
{
  PoreWindow window;
  window.position = position;

  double normalLength = normal.length();
  if (normalLength <= 0.0 || network.atomPositionsFractional.empty()) return window;
  double3 axisNormal = normal / normalLength;
  window.normal = axisNormal;

  // Two directions spanning the plane. The first is the smallest component of the normal turned into a
  // vector that cannot be parallel to it.
  double3 helper = (std::abs(axisNormal.x) <= std::abs(axisNormal.y) && std::abs(axisNormal.x) <= std::abs(axisNormal.z))
                       ? double3(1.0, 0.0, 0.0)
                       : ((std::abs(axisNormal.y) <= std::abs(axisNormal.z)) ? double3(0.0, 1.0, 0.0)
                                                                            : double3(0.0, 0.0, 1.0));
  double3 firstAxis = double3::cross(axisNormal, helper);
  firstAxis = firstAxis / firstAxis.length();
  double3 secondAxis = double3::cross(axisNormal, firstAxis);

  const double3x3 cell = network.simulationBox.cell;
  const double3x3 inverseCell = network.simulationBox.inverseCell;
  const std::vector<double>& radii = network.atomRadii;
  double3 perpendicularWidths = perpendicularWidthsOf(network.simulationBox);
  double maximumRadius = *std::max_element(radii.begin(), radii.end());
  double3 positionFractional = inverseCell * position;

  // Visits every atom image whose surface comes within `searchRadius` of the bottleneck. Each atom is
  // first brought to its image nearest the bottleneck, so the loop over lattice offsets is symmetric
  // around it however the cell is sheared.
  auto forEachAtomImage = [&](double searchRadius, auto&& visit)
  {
    double extent = searchRadius + maximumRadius;
    int3 span(static_cast<std::int32_t>(std::ceil(extent / perpendicularWidths.x)) + 1,
              static_cast<std::int32_t>(std::ceil(extent / perpendicularWidths.y)) + 1,
              static_cast<std::int32_t>(std::ceil(extent / perpendicularWidths.z)) + 1);

    for (std::size_t j = 0; j < network.atomPositionsFractional.size(); ++j)
    {
      double3 relative = network.atomPositionsFractional[j] - positionFractional;
      double3 nearest(relative.x - std::round(relative.x), relative.y - std::round(relative.y),
                      relative.z - std::round(relative.z));

      for (std::int32_t ox = -span.x; ox <= span.x; ++ox)
      {
        for (std::int32_t oy = -span.y; oy <= span.y; ++oy)
        {
          for (std::int32_t oz = -span.z; oz <= span.z; ++oz)
          {
            double3 delta = cell * double3(nearest.x + static_cast<double>(ox), nearest.y + static_cast<double>(oy),
                                           nearest.z + static_cast<double>(oz));
            if (delta.length() - radii[j] <= searchRadius) visit(j, delta);
          }
        }
      }
    }
  };

  // How far the window is followed. Beyond half the narrowest width of the cell a direction has left
  // the neighbourhood of this window for that of its own periodic images, and a window wide enough to
  // reach that far is not a window. The clearance sets the other end of the scale, since a window is of
  // the order of its own inscribed circle, and the larger of the two is taken so that neither a small
  // cell nor a wide passage cuts the measurement short unnecessarily.
  double cellReach = 0.5 * std::min({perpendicularWidths.x, perpendicularWidths.y, perpendicularWidths.z});

  // The clearance first, since it sets the scale of everything else. Half the cell is the natural place
  // to start looking but is not always far enough, a passage in a small cell being able to be wider than
  // that, so the search widens until it finds the surface it is looking for.
  double clearance = std::numeric_limits<double>::max();
  for (double searchRadius = cellReach;
       clearance == std::numeric_limits<double>::max() && searchRadius < 64.0 * cellReach; searchRadius *= 2.0)
  {
    forEachAtomImage(searchRadius, [&](std::size_t j, const double3& delta)
                     { clearance = std::min(clearance, delta.length() - radii[j]); });
  }
  if (clearance == std::numeric_limits<double>::max()) return window;
  window.freeRadius = clearance;

  // A bottleneck of no clearance is a passage that is shut, and a shut passage has no window: the plane
  // meets the atoms at the bottleneck itself and there is nothing free to measure.
  if (clearance <= 1.0e-9) return window;

  double reach = std::max(cellReach, 3.0 * clearance);

  std::vector<WindowDisc> discs;
  forEachAtomImage(reach,
                   [&](std::size_t j, const double3& delta)
                   {
                     double standoff = double3::dot(delta, axisNormal);
                     double squared = radii[j] * radii[j] - standoff * standoff;
                     if (squared <= 0.0) return;  // the atom does not reach the plane of the window

                     WindowDisc disc;
                     disc.first = double3::dot(delta, firstAxis);
                     disc.second = double3::dot(delta, secondAxis);
                     disc.radius = std::sqrt(squared);
                     disc.gap = disc.first * disc.first + disc.second * disc.second - squared;
                     if (disc.gap <= 0.0) return;  // cannot happen at positive clearance
                     if (std::sqrt(disc.first * disc.first + disc.second * disc.second) - disc.radius > reach) return;
                     discs.push_back(disc);
                   });

  // How far the free cross-section extends from the bottleneck in each direction of the plane, that is
  // the distance to the first atom met. What this traces out is the part of the window visible from the
  // bottleneck, which is a subset of the window whatever shape the window has.
  constexpr std::size_t directionCount = 360;  // even, so that opposite directions are sampled in pairs
  std::array<double, directionCount> extent{};
  std::array<double, directionCount> cosine{};
  std::array<double, directionCount> sine{};
  std::vector<char> bounds(discs.size(), 0);

  for (std::size_t i = 0; i < directionCount; ++i)
  {
    double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(directionCount);
    cosine[i] = std::cos(angle);
    sine[i] = std::sin(angle);

    double nearest = reach;
    std::size_t nearestDisc = discs.size();
    for (std::size_t d = 0; d < discs.size(); ++d)
    {
      const WindowDisc& disc = discs[d];
      double along = cosine[i] * disc.first + sine[i] * disc.second;
      if (along <= 0.0) continue;
      double discriminant = along * along - disc.gap;
      if (discriminant <= 0.0) continue;
      double distance = along - std::sqrt(discriminant);
      if (distance < nearest)
      {
        nearest = distance;
        nearestDisc = d;
      }
    }
    extent[i] = nearest;
    if (nearestDisc < discs.size())
      bounds[nearestDisc] = 1;
    else
      window.clipped = true;
  }

  // The ring is counted by image, since a window in a small cell can be ringed by several images of the
  // same atom, and each of them is one atom of the ring.
  for (char bounding : bounds) window.boundingAtoms += static_cast<std::size_t>(bounding);

  // The free chord through the bottleneck, in each direction and the one opposite it.
  window.smallestFreeWidth = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < directionCount / 2; ++i)
  {
    double width = extent[i] + extent[i + directionCount / 2];
    window.smallestFreeWidth = std::min(window.smallestFreeWidth, width);
    window.largestFreeWidth = std::max(window.largestFreeWidth, width);
  }

  // The largest-area ellipse centred on the bottleneck that fits inside what was traced out. An ellipse
  // of semi-axes (elongation·b) along its own first axis and b along the other reaches
  // b·elongation/sqrt(cos²+elongation²sin²) at the angle measured from that axis, so it fits exactly
  // when that is within the extent in every direction, which gives b in closed form for a given
  // orientation and elongation. Both of those are scanned, coarsely and then again around the best.
  std::array<double, directionCount> cosineSquared{};
  std::array<double, directionCount> sineSquared{};
  for (std::size_t i = 0; i < directionCount; ++i)
  {
    cosineSquared[i] = cosine[i] * cosine[i];
    sineSquared[i] = sine[i] * sine[i];
  }

  auto semiAxis = [&](std::size_t orientation, double elongation)
  {
    double smallest = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < directionCount; ++i)
    {
      std::size_t offset = (i + directionCount - orientation) % directionCount;
      double allowed = extent[i] * std::sqrt(cosineSquared[offset] / (elongation * elongation) + sineSquared[offset]);
      smallest = std::min(smallest, allowed);
    }
    return smallest;
  };

  std::size_t bestOrientation = 0;
  double bestElongation = 1.0;
  double bestSemiAxis = 0.0;
  double bestArea = 0.0;

  auto consider = [&](std::size_t orientation, double elongation)
  {
    double semi = semiAxis(orientation, elongation);
    double area = elongation * semi * semi;
    if (area > bestArea)
    {
      bestArea = area;
      bestSemiAxis = semi;
      bestElongation = elongation;
      bestOrientation = orientation;
    }
  };

  // The ellipse is symmetric under a half turn, so only half the orientations are distinct.
  constexpr std::size_t orientationStride = 3;
  constexpr std::size_t elongationSteps = 32;
  constexpr double largestElongation = 8.0;
  for (std::size_t orientation = 0; orientation < directionCount / 2; orientation += orientationStride)
  {
    for (std::size_t step = 0; step <= elongationSteps; ++step)
    {
      double elongation = std::pow(largestElongation, static_cast<double>(step) / static_cast<double>(elongationSteps));
      consider(orientation, elongation);
    }
  }

  std::size_t coarseOrientation = bestOrientation;
  double coarseElongation = bestElongation;
  for (std::size_t offset = 1; offset < orientationStride; ++offset)
  {
    consider((coarseOrientation + offset) % (directionCount / 2), coarseElongation);
    consider((coarseOrientation + directionCount / 2 - offset) % (directionCount / 2), coarseElongation);
  }
  double band = std::pow(largestElongation, 1.0 / static_cast<double>(elongationSteps));
  for (std::size_t step = 0; step <= 16; ++step)
  {
    double factor = std::pow(band, -1.0 + 2.0 * static_cast<double>(step) / 16.0);
    consider(bestOrientation, std::max(1.0, coarseElongation * factor));
  }

  window.minorAxis = 2.0 * bestSemiAxis;
  window.majorAxis = 2.0 * bestElongation * bestSemiAxis;
  window.majorAxisDirection = cosine[bestOrientation] * firstAxis + sine[bestOrientation] * secondAxis;
  window.measured = true;

  return window;
}

PoreWindow freeSphereWindow(const VoronoiNetwork& network)
{
  std::vector<std::size_t> allNodes(network.nodes.size());
  std::iota(allNodes.begin(), allNodes.end(), std::size_t{0});

  PercolatingPath path = widestPercolatingPath(network, allNodes);
  if (!path.percolates) return PoreWindow{};

  const VoronoiEdge& edge = network.edges[path.limitingEdge];
  if (!edge.hasBottleneckGeometry) return PoreWindow{};

  return PoreWindow::measure(network, edge.bottleneckPosition, edge.bottleneckDirection);
}

std::vector<ChannelWindow> channelWindows(const VoronoiNetwork& network, const ChannelAnalysis& channels)
{
  std::vector<ChannelWindow> windows;

  for (std::size_t poreIndex = 0; poreIndex < channels.pores.size(); ++poreIndex)
  {
    const VoronoiPore& pore = channels.pores[poreIndex];
    if (!pore.isChannel) continue;

    PercolatingPath path = widestPercolatingPath(network, pore.nodeIndices);
    if (!path.percolates) continue;

    ChannelWindow entry;
    entry.poreIndex = poreIndex;
    entry.dimensionality = pore.dimensionality;
    entry.freeSphereDiameter = 2.0 * path.radius;
    entry.limitingEdge = path.limitingEdge;

    const VoronoiEdge& edge = network.edges[path.limitingEdge];
    if (edge.hasBottleneckGeometry)
    {
      entry.window = PoreWindow::measure(network, edge.bottleneckPosition, edge.bottleneckDirection);
    }

    windows.push_back(std::move(entry));
  }

  return windows;
}
