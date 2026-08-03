module;

module energy_shared_isosurface;

import std;

import uint3;
import double3;
import double3x3;

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
  const double largestPlausible = largestPlausibleTriangleArea(unitCell, gridSize);

  IsosurfaceArea result;
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
    }
    else
    {
      ++result.numberOfRejectedTriangles;
    }
  }

  return result;
}
