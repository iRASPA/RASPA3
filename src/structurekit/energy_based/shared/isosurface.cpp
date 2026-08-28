module;

module energy_shared_isosurface;

import std;

import uint3;
import double3;
import double3x3;
import surface_curvature;

double largestPlausibleTriangleArea(const double3x3 &unitCell, uint3 gridSize)
{
  double3 voxelA = unitCell * double3(1.0 / static_cast<double>(gridSize.x), 0.0, 0.0);
  double3 voxelB = unitCell * double3(0.0, 1.0 / static_cast<double>(gridSize.y), 0.0);
  double3 voxelC = unitCell * double3(0.0, 0.0, 1.0 / static_cast<double>(gridSize.z));

  return 2.0 * std::max({double3::cross(voxelA, voxelB).length(), double3::cross(voxelB, voxelC).length(),
                         double3::cross(voxelC, voxelA).length()});
}

IsosurfaceArea accumulateTriangleAreas(const double3x3 &unitCell, uint3 gridSize, std::span<const double3> corners)
{
  return accumulateTriangleAreas(unitCell, gridSize, corners, std::span<const double3>{}, FieldSense::GrowsIntoVoid);
}


IsosurfaceArea accumulateTriangleAreas(const double3x3 &unitCell, uint3 gridSize, std::span<const double3> corners,
                                       std::span<const double3> gradients, FieldSense sense)
{
  const double largestPlausible = largestPlausibleTriangleArea(unitCell, gridSize);
  const bool withCurvature = gradients.size() >= corners.size();
  const double3x3 inverseCell = unitCell.inverse();

  IsosurfaceArea result;
  result.curvatureBands = curvatureBandsForGrid(unitCell, gridSize);

  for (std::size_t i = 0; i + 2 < corners.size(); i += 3)
  {
    double3 p1 = unitCell * corners[i];
    double3 p2 = unitCell * corners[i + 1];
    double3 p3 = unitCell * corners[i + 2];

    double area = 0.5 * double3::cross(p2 - p1, p3 - p1).length();
    if (std::isfinite(area) && area < largestPlausible)
    {
      result.area += area;
      ++result.numberOfTriangles;

      if (withCurvature)
      {
        std::array<double3, 3> triangleCorners{p1, p2, p3};
        std::array<double3, 3> triangleNormals{
            cartesianNormalOfGridGradient(inverseCell, gridSize, gradients[i], sense),
            cartesianNormalOfGridGradient(inverseCell, gridSize, gradients[i + 1], sense),
            cartesianNormalOfGridGradient(inverseCell, gridSize, gradients[i + 2], sense)};

        result.curvature.add(area, triangleCurvature(triangleCorners, triangleNormals), result.curvatureBands);
      }
    }
    else
    {
      ++result.numberOfRejectedTriangles;
    }
  }

  return result;
}


void writeIsosurfaceCurvature(std::ostream &stream, const IsosurfaceArea &surface, const char *name)
{
  const CurvatureAreas &areas = surface.curvature;
  if (areas.total() <= 0.0) return;

  std::print(stream, "#\n");
  std::print(stream, "# How that area divides by the shape of the surface.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Marching cubes places each vertex where the field crosses the iso-value along a cube edge\n");
  std::print(stream, "# and records the field's gradient there, which is the normal of the level set. Three of\n");
  std::print(stream, "# those to a triangle fix how the normal turns across it, and that is the two principal\n");
  std::print(stream, "# curvatures. The energy grows into the wall, so the outward normal is minus the gradient,\n");
  std::print(stream, "# and with the normal that way round a curvature is positive where the wall bulges out into\n");
  std::print(stream, "# the void: both positive is convex, both negative is the inside of a pocket and concave,\n");
  std::print(stream, "# one of each is a saddle.\n");
  std::print(stream, "#\n");
  std::print(stream, "# On an energy field these are properties of the surface and not of the grid. The field is a\n");
  std::print(stream, "# sum over the atoms, so it is smooth everywhere: the walls of neighbouring atoms blend into\n");
  std::print(stream, "# one another over a distance set by how fast the repulsion climbs, and the saddle between\n");
  std::print(stream, "# two of them is a genuine saddle with a curvature of its own. Refining the grid sharpens\n");
  std::print(stream, "# the estimate rather than moving the answer. The clearance-grid route is not like this and\n");
  std::print(stream, "# its own report says so: that field is a minimum over the atoms rather than a sum, so it\n");
  std::print(stream, "# creases wherever two are equally near, and there the split does move with the spacing.\n");
  std::print(stream, "#\n");
  std::print(stream, "# It is also not Richards's convex/saddle/concave split of the solvent-excluded surface. That\n");
  std::print(stream, "# one counts how many atoms a hard probe touches at once, on a surface built out of spheres\n");
  std::print(stream, "# and tori. This surface is a contour of a continuous energy and has no such construction\n");
  std::print(stream, "# behind it, so the two are not comparable even where both are trustworthy.\n");
  std::print(stream, "#\n");
  std::print(stream, "# A curvature below {:.5f} [1/Å], a radius of {:.1f} [Å] and longer, is reported as flat\n",
             surface.curvatureBands.flat, 1.0 / surface.curvatureBands.flat);
  std::print(stream, "# rather than given a sign it does not have. One above {:.5f} [1/Å] is a radius shorter than\n",
             surface.curvatureBands.sharpest);
  std::print(stream, "# a voxel, which the grid cannot carry, and is held out of the columns but not of the\n");
  std::print(stream, "# integrals. A large unresolved share means the grid is too coarse for this question.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The first four fractions are of the classified area and add to one; the last is the\n");
  std::print(stream, "# unresolved share of the whole.\n");
  std::print(stream, "#                              area [Å²]     convex     saddle    concave       flat unresolved\n");
  std::print(stream, "{:<26} {:13.5f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n", name, areas.classified(),
             areas.convexFraction(), areas.saddleFraction(), areas.concaveFraction(), areas.flatFraction(),
             areas.unresolvedFraction());
  std::print(stream, "# Triangles: {} convex, {} saddle, {} concave, {} flat, {} unresolved, {} of those\n",
             areas.numberOfConvex, areas.numberOfSaddle, areas.numberOfConcave, areas.numberOfFlat,
             areas.numberOfUnresolved, areas.numberOfDegenerate);
  std::print(stream, "# degenerate rather than too sharp.\n");
  std::print(stream, "Integral of the mean curvature:     {} [Å]\n", areas.integratedMeanCurvature);
  std::print(stream, "Integral of the Gaussian curvature: {} [-]\n", areas.integratedGaussianCurvature);
  std::print(stream, "# The second is 2 pi times the Euler characteristic over a closed surface, so it counts\n");
  std::print(stream, "# rather than measures: negative on a network of channels, and the more negative the more\n");
  std::print(stream, "# connected the network. It is the noisier of the two, the fit being a least-squares one\n");
  std::print(stream, "# over three normals rather than the angle deficit that satisfies discrete Gauss-Bonnet\n");
  std::print(stream, "# exactly, so read its sign and its size and not past the first digit or two.\n");
}
