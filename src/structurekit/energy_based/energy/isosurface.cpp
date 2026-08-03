module;

module energy_isosurface;

import std;

import uint3;
import double3;
import double3x3;
import marching_cubes;
import unit_cell;
import crystal;

import energy_shared_isosurface;

std::vector<double3> EnergyIsosurface::trianglesOfIsosurface(std::span<const float> field, uint3 gridSize,
                                                             double isoValue)
{
  std::size_t numberOfGridPoints = gridSize.x * gridSize.y * gridSize.z;
  if (field.size() < numberOfGridPoints)
  {
    throw std::runtime_error(std::format(
        "EnergyIsosurface: the field has {} values, too few for a {} x {} x {} grid\n", field.size(), gridSize.x,
        gridSize.y, gridSize.z));
  }

  // The cell is periodic, so the cube spanning its far face is the one that runs back round to the near face,
  // and leaving it out loses a whole boundary layer of surface. The GPU extractor reaches round with a modulo
  // on the sample index; this one cannot, so the field is given an extra plane on each side that repeats the
  // first. Either way every cube of the cell is counted exactly once.
  uint3 extended{gridSize.x + 1, gridSize.y + 1, gridSize.z + 1};

  MarchingCubes cube(static_cast<int>(extended.x), static_cast<int>(extended.y), static_cast<int>(extended.z));
  cube.init_all();

  // The field is laid out with x varying fastest, which is the order every field of this route is written in.
  for (std::size_t k = 0; k < extended.z; ++k)
  {
    for (std::size_t j = 0; j < extended.y; ++j)
    {
      for (std::size_t i = 0; i < extended.x; ++i)
      {
        std::size_t source =
            ((k % gridSize.z) * gridSize.y + (j % gridSize.y)) * gridSize.x + (i % gridSize.x);
        cube.set_data(static_cast<double>(field[source]), i, j, k);
      }
    }
  }

  cube.run(isoValue);

  // The repeated plane sits at index gridSize, which is fractional one, so the vertices scale by the original
  // grid rather than by the extended one.
  std::vector<double3> corners;
  corners.reserve(3 * cube.ntrigs());
  for (std::size_t i = 0; i < cube.ntrigs(); ++i)
  {
    const Triangle *triangle = cube.trig(static_cast<std::ptrdiff_t>(i));
    for (int index : {triangle->v1, triangle->v2, triangle->v3})
    {
      const Vertex *vertex = cube.vert(index);
      corners.emplace_back(vertex->x / static_cast<double>(gridSize.x), vertex->y / static_cast<double>(gridSize.y),
                           vertex->z / static_cast<double>(gridSize.z));
    }
  }

  return corners;
}


IsosurfaceArea EnergyIsosurface::areaOfIsosurface(const Crystal &framework, std::span<const float> field,
                                                  uint3 gridSize, double isoValue)
{
  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, isoValue);
  return accumulateTriangleAreas(framework.unitCell.cell, gridSize, corners);
}
