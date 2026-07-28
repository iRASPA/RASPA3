module;

module exact_union_volume;

import std;

import double2;
import double3;
import voronoi_accessibility;
import exact_surface_patches;

// A plane bounding one atom's cell, as the direction of the neighbour that put it there and the signed
// distance from the atom's centre to it. The distance is negative when the plane lies behind the centre,
// which happens where a larger neighbour cuts deeply, and the sum below wants that sign.
struct CellPlane
{
  double3 axis;
  double distance{0.0};
};

// One of the lines a face is cut back by, in the two-dimensional frame of that face: the face keeps the
// side where the inner product with `normal` is at most `offset`.
struct CutLine
{
  double2 normal;
  double offset{0.0};
};

double cross(const double2& a, const double2& b) { return a.x * b.y - a.y * b.x; }
double2 scale(const double2& a, double factor) { return double2(a.x * factor, a.y * factor); }

// The sector of the disc between two directions, signed by the turn from the first to the second. Only
// the directions of the arguments are used, so either may lie outside the circle.
double sectorArea(double radius, const double2& from, const double2& to)
{
  return 0.5 * radius * radius * std::atan2(cross(from, to), double2::dot(from, to));
}

double triangleArea(const double2& from, const double2& to) { return 0.5 * cross(from, to); }

// The area swept from the origin as one edge of a polygon is traversed, counting only what is inside the
// disc: a triangle where the edge runs inside the circle, and a sector where it runs outside. Summed
// over a closed polygon the sweeps cancel except over the overlap, which is the area wanted.
double sweptArea(double radius, const double2& from, const double2& to)
{
  const double radiusSquared = radius * radius;
  const bool fromInside = double2::dot(from, from) <= radiusSquared;
  const bool toInside = double2::dot(to, to) <= radiusSquared;
  if (fromInside && toInside) return triangleArea(from, to);

  // Where the edge crosses the circle: the roots of |from + s (to - from)|^2 = radius^2.
  const double2 along = to - from;
  const double quadratic = double2::dot(along, along);
  if (quadratic < 1.0e-300) return 0.0;  // a repeated vertex sweeps nothing
  const double linear = 2.0 * double2::dot(from, along);
  const double constant = double2::dot(from, from) - radiusSquared;
  const double discriminant = linear * linear - 4.0 * quadratic * constant;

  if (fromInside || toInside)
  {
    // One end inside: there is exactly one crossing, and the other root lies outside [0, 1].
    const double root = std::sqrt(std::max(0.0, discriminant));
    const double crossing = fromInside ? (-linear + root) / (2.0 * quadratic) : (-linear - root) / (2.0 * quadratic);
    const double2 meeting = from + scale(along, std::clamp(crossing, 0.0, 1.0));
    return fromInside ? triangleArea(from, meeting) + sectorArea(radius, meeting, to)
                      : sectorArea(radius, from, meeting) + triangleArea(meeting, to);
  }

  // Both ends outside: either the edge misses the circle, and the sweep is one sector, or it cuts across
  // and the sweep is a sector, the chord, and a sector.
  if (discriminant <= 0.0) return sectorArea(radius, from, to);
  const double root = std::sqrt(discriminant);
  const double first = (-linear - root) / (2.0 * quadratic);
  const double second = (-linear + root) / (2.0 * quadratic);
  if (second <= 0.0 || first >= 1.0) return sectorArea(radius, from, to);

  const double2 entry = from + scale(along, std::max(0.0, first));
  const double2 exit = from + scale(along, std::min(1.0, second));
  return sectorArea(radius, from, entry) + triangleArea(entry, exit) + sectorArea(radius, exit, to);
}

// The polygon kept on one side of a line, by the usual convex clip: vertices on the near side are kept
// and an edge that changes side contributes the point where it crosses.
void clipToLine(std::vector<double2>& polygon, std::vector<double2>& scratch, const CutLine& line)
{
  scratch.clear();
  for (std::size_t i = 0; i < polygon.size(); ++i)
  {
    const double2& current = polygon[i];
    const double2& next = polygon[(i + 1) % polygon.size()];
    double here = double2::dot(line.normal, current) - line.offset;
    double there = double2::dot(line.normal, next) - line.offset;

    if (here <= 0.0) scratch.push_back(current);
    if ((here < 0.0 && there > 0.0) || (here > 0.0 && there < 0.0))
    {
      scratch.push_back(current + scale(next - current, here / (here - there)));
    }
  }
  polygon.swap(scratch);
}

// The area of one face: the disc in which the plane cuts the atom, cut back by the other planes of the
// cell. The disc is approached from outside, by cutting a square that contains it, so that a face no
// line touches comes out as the whole disc.
double faceArea(double discRadius, const std::vector<CutLine>& lines, std::vector<double2>& polygon,
                std::vector<double2>& scratch)
{
  const double corner = 2.0 * discRadius;
  polygon.assign({double2(-corner, -corner), double2(corner, -corner), double2(corner, corner),
                  double2(-corner, corner)});
  for (const CutLine& line : lines)
  {
    clipToLine(polygon, scratch, line);
    if (polygon.size() < 3) return 0.0;
  }

  double total = 0.0;
  for (std::size_t i = 0; i < polygon.size(); ++i)
  {
    total += sweptArea(discRadius, polygon[i], polygon[(i + 1) % polygon.size()]);
  }
  return std::abs(total);
}

double clippedDiscArea(double radius, const std::vector<double2>& normals, const std::vector<double>& offsets)
{
  if (radius <= 0.0) return 0.0;

  std::vector<CutLine> lines;
  lines.reserve(normals.size());
  for (std::size_t k = 0; k < normals.size(); ++k)
  {
    // A line whose normal has collapsed is not a line: it either leaves the disc alone or removes it.
    if (double2::dot(normals[k], normals[k]) < 1.0e-24)
    {
      if (offsets[k] < 0.0) return 0.0;
      continue;
    }
    lines.push_back({normals[k], offsets[k]});
  }

  std::vector<double2> polygon;
  std::vector<double2> scratch;
  return faceArea(radius, lines, polygon, scratch);
}

double unionOfBallsVolume(const VoronoiAccessibility& accessibility, std::size_t subdivisions)
{
  // The spheres. Each atom's exposed patch enters weighted by its own radius, which is what the field
  // x - p_i contributes there, and the classifier is not needed for a volume.
  return unionOfBallsVolume(accessibility, exactAccessibleSurfaceArea(accessibility, subdivisions, false));
}

double unionOfBallsVolume(const VoronoiAccessibility& accessibility, const ExactSurfaceAreaSample& patches)
{
  double total = patches.radiusWeightedArea;

  // The faces. Every plane that cuts an atom bounds that atom's cell there, and a plane that misses the
  // atom cannot cut a face either, since a face lies inside the atom, so the planes collected here are
  // all the ones any face is cut by.
  std::vector<CellPlane> planes;
  std::vector<CutLine> lines;
  std::vector<double2> polygon;
  std::vector<double2> scratch;

  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const double radius = accessibility.atomRadii[i];
    const double3 centre = accessibility.atomPositions[i];

    planes.clear();
    bool outsideOwnCell = false;
    for (const auto& [delta, neighbourRadius] :
         accessibility.neighbourAtoms(centre, radius + accessibility.maximumAtomRadius))
    {
      double distance = delta.length();
      if (distance < 1.0e-12)
      {
        // The atom itself, which the query returns along with the rest. A second atom on the very same
        // spot is a duplicate in the input rather than a geometry: a larger one buries this atom, and
        // one of the same size leaves a tie that no diagram can break, which is left alone here as it is
        // in the area.
        if (neighbourRadius > radius) outsideOwnCell = true;
        continue;
      }

      // The plane on which this atom and the neighbour are equally far into each other's territory.
      double planeDistance =
          (distance * distance + radius * radius - neighbourRadius * neighbourRadius) / (2.0 * distance);
      if (planeDistance >= radius) continue;  // the plane misses the atom and bounds nothing of it
      if (planeDistance <= -radius)
      {
        // The whole atom is on the far side of one of its own cell's planes, so its cell holds none of
        // it: the neighbour has swallowed it and carries the volume itself.
        outsideOwnCell = true;
        break;
      }
      planes.push_back({delta * (1.0 / distance), planeDistance});
    }
    if (outsideOwnCell) continue;

    for (std::size_t j = 0; j < planes.size(); ++j)
    {
      const CellPlane& face = planes[j];

      // A frame in the plane of the face. Which pair of perpendicular directions is used cannot matter,
      // an area being what is asked for.
      double3 helper = (std::abs(face.axis.x) < 0.9) ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
      double3 first = double3::cross(helper, face.axis);
      first = first * (1.0 / first.length());
      double3 second = double3::cross(face.axis, first);

      lines.clear();
      bool empty = false;
      for (std::size_t k = 0; k < planes.size() && !empty; ++k)
      {
        if (k == j) continue;
        const CellPlane& other = planes[k];
        CutLine line;
        line.normal = double2(double3::dot(first, other.axis), double3::dot(second, other.axis));
        line.offset = other.distance - face.distance * double3::dot(face.axis, other.axis);

        // A plane parallel to this face either leaves it alone or does away with it altogether.
        if (double2::dot(line.normal, line.normal) < 1.0e-24)
        {
          empty = line.offset < 0.0;
          continue;
        }
        lines.push_back(line);
      }
      if (empty) continue;

      double discRadius = std::sqrt(std::max(0.0, radius * radius - face.distance * face.distance));
      total += face.distance * faceArea(discRadius, lines, polygon, scratch);
    }
  }

  return total / 3.0;
}
