module;

module phonon_kpath;

import std;

import double3;
import double3x3;
import simulationbox;

std::vector<PhononKPoint> buildPhononKPath(std::span<const PhononPathNode> nodes, std::size_t pointsPerSegment,
                                           const SimulationBox& box)
{
  std::vector<PhononKPoint> path;
  if (nodes.empty()) return path;

  const std::size_t steps = std::max<std::size_t>(1, pointsPerSegment);

  // Reciprocal-lattice basis (rows of the inverse cell), matching the dynamical-matrix convention.
  const double3x3 inverseCell = box.inverseCell;
  const double3 aStar(inverseCell.ax, inverseCell.bx, inverseCell.cx);
  const double3 bStar(inverseCell.ay, inverseCell.by, inverseCell.cy);
  const double3 cStar(inverseCell.az, inverseCell.bz, inverseCell.cz);
  constexpr double twoPi = 2.0 * std::numbers::pi;
  const auto cartesian = [&](const double3& kFractional)
  { return twoPi * (kFractional.x * aStar + kFractional.y * bStar + kFractional.z * cStar); };

  double cumulative = 0.0;
  path.push_back(PhononKPoint{.kFractional = nodes[0].kPoint, .pathCoordinate = 0.0, .label = nodes[0].label});
  double3 previousCartesian = cartesian(nodes[0].kPoint);

  for (std::size_t node = 1; node < nodes.size(); ++node)
  {
    if (nodes[node].startsNewSegment)
    {
      // Discontinuity: begin a new sub-path at the current coordinate without drawing a connecting segment.
      previousCartesian = cartesian(nodes[node].kPoint);
      path.push_back(
          PhononKPoint{.kFractional = nodes[node].kPoint, .pathCoordinate = cumulative, .label = nodes[node].label});
      continue;
    }

    const double3 start = nodes[node - 1].kPoint;
    const double3 end = nodes[node].kPoint;
    for (std::size_t step = 1; step <= steps; ++step)
    {
      const double fraction = static_cast<double>(step) / static_cast<double>(steps);
      const double3 kFractional = start + fraction * (end - start);
      const double3 kCartesian = cartesian(kFractional);
      cumulative += (kCartesian - previousCartesian).length();
      previousCartesian = kCartesian;

      const std::string label = (step == steps) ? nodes[node].label : std::string{};
      path.push_back(PhononKPoint{.kFractional = kFractional, .pathCoordinate = cumulative, .label = label});
    }
  }
  return path;
}
