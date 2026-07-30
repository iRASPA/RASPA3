module;

module exact_union_volume;

import std;

import double2;
import double3;
import pore_accessibility;
import exact_boundary_components;
import exact_sphere_sweep;
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

// Whether a line has to be clipped against at all, given that only what falls inside the disc of radius
// `discRadius` about the origin is ever asked for.
//
// A line that leaves the whole disc on the side it keeps cannot change the area inside the disc, so it may
// be dropped outright; one that leaves the whole disc on the side it discards empties the face. Most of the
// planes bounding a cell do one or the other on any given face -- they were collected because they cut the
// atom, not because they cut this face of it -- and skipping them is what keeps the clip below over the
// handful that really cut it. A normal that has collapsed to nothing falls out as the first case or the
// second according to the sign of its offset, which is what a parallel plane should do anyway.
enum class LineEffect : std::uint8_t
{
  cuts,
  misses,
  empties
};

LineEffect lineEffect(const CutLine& line, double discRadius)
{
  const double reach = discRadius * std::sqrt(double2::dot(line.normal, line.normal));
  if (line.offset >= reach) return LineEffect::misses;
  if (line.offset <= -reach) return LineEffect::empties;
  return LineEffect::cuts;
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
    const CutLine line{normals[k], offsets[k]};
    switch (lineEffect(line, radius))
    {
      case LineEffect::misses:
        continue;
      case LineEffect::empties:
        return 0.0;
      case LineEffect::cuts:
        lines.push_back(line);
        break;
    }
  }

  std::vector<double2> polygon;
  std::vector<double2> scratch;
  return faceArea(radius, lines, polygon, scratch);
}

double unionOfBallsVolume(const PoreAccessibility& accessibility, std::size_t subdivisions)
{
  // The spheres. Each atom's exposed patch enters weighted by its own radius, which is what the field
  // x - p_i contributes there, and the classifier is not needed for a volume.
  return unionOfBallsVolume(accessibility, exactSurfaceArea(accessibility, subdivisions));
}

// The volume of the union, given something that says which planes bound each atom's cell.
//
// `cellPlanes(i, planes)` fills `planes` for atom i and returns false where the atom's cell holds none of it,
// which is the one thing the two routes below have to answer and the only thing they differ in.
template <typename Planes>
double sumOverCells(const PoreAccessibility& accessibility, const MeasuredPatches& patches,
                    Planes&& cellPlanes)
{
  double total = patches.radiusWeightedArea;

  std::vector<CellPlane> planes;
  std::vector<CutLine> lines;
  std::vector<double2> polygon;
  std::vector<double2> scratch;

  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const double radius = accessibility.atomRadii[i];
    if (!cellPlanes(i, planes)) continue;

    for (std::size_t j = 0; j < planes.size(); ++j)
    {
      const CellPlane& face = planes[j];

      // A frame in the plane of the face. Which pair of perpendicular directions is used cannot matter,
      // an area being what is asked for.
      double3 helper = (std::abs(face.axis.x) < 0.9) ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
      double3 first = double3::cross(helper, face.axis);
      first = first * (1.0 / first.length());
      double3 second = double3::cross(face.axis, first);

      const double discRadius = std::sqrt(std::max(0.0, radius * radius - face.distance * face.distance));

      lines.clear();
      bool empty = false;
      for (std::size_t k = 0; k < planes.size() && !empty; ++k)
      {
        if (k == j) continue;
        const CellPlane& other = planes[k];
        CutLine line;
        line.normal = double2(double3::dot(first, other.axis), double3::dot(second, other.axis));
        line.offset = other.distance - face.distance * double3::dot(face.axis, other.axis);

        switch (lineEffect(line, discRadius))
        {
          case LineEffect::misses:
            break;
          case LineEffect::empties:
            empty = true;
            break;
          case LineEffect::cuts:
            lines.push_back(line);
            break;
        }
      }
      if (empty) continue;

      total += face.distance * faceArea(discRadius, lines, polygon, scratch);
    }
  }

  return total / 3.0;
}

double unionOfBallsVolume(const PoreAccessibility& accessibility, const MeasuredPatches& patches)
{
  // Every plane that cuts an atom bounds that atom's cell there, and a plane that misses the atom cannot cut
  // a face either, since a face lies inside the atom, so the planes collected here are all the ones any face
  // is cut by.
  return sumOverCells(
      accessibility, patches,
      [&](std::size_t i, std::vector<CellPlane>& planes)
      {
        const double radius = accessibility.atomRadii[i];
        const double3 centre = accessibility.atomPositions[i];

        planes.clear();
        for (const auto& [delta, neighbourRadius] :
             accessibility.neighbourAtoms(centre, radius + accessibility.maximumAtomRadius))
        {
          double distance = delta.length();
          if (distance < 1.0e-12)
          {
            // The atom itself, which the query returns along with the rest. A second atom on the very same
            // spot is a duplicate in the input rather than a geometry: a larger one buries this atom, and one
            // of the same size leaves a tie that no diagram can break, which is left alone here as it is in
            // the area.
            if (neighbourRadius > radius) return false;
            continue;
          }

          // The plane on which this atom and the neighbour are equally far into each other's territory.
          double planeDistance =
              (distance * distance + radius * radius - neighbourRadius * neighbourRadius) / (2.0 * distance);
          if (planeDistance >= radius) continue;  // the plane misses the atom and bounds nothing of it

          // The whole atom is on the far side of one of its own cell's planes, so its cell holds none of it:
          // the neighbour has swallowed it and carries the volume itself.
          if (planeDistance <= -radius) return false;

          planes.push_back({delta * (1.0 / distance), planeDistance});
        }
        return true;
      });
}

double unionOfBallsVolume(const PoreAccessibility& accessibility, const MeasuredPatches& patches,
                          const BoundaryComponents& components)
{
  return sumOverCells(accessibility, patches,
                      [&](std::size_t i, std::vector<CellPlane>& planes)
                      {
                        const SphereBoundary& boundary = components.atoms[i];
                        planes.clear();
                        if (boundary.buried) return false;

                        const double radius = accessibility.atomRadii[i];
                        planes.reserve(boundary.circles.size());
                        for (const SphereCircle& circle : boundary.circles)
                        {
                          planes.push_back({circle.axis, radius * circle.cosineHalfAngle});
                        }
                        return true;
                      });
}
