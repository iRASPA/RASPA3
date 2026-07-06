#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import skvoronoi;

// A single site in a cubic cell: the Voronoi cell is the cube itself.
TEST(voronoi, single_site_cubic_cell)
{
  double3x3 unitCell(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  SKVoronoi voronoi(unitCell, {double3(0.25, 0.5, 0.75)});

  SKVoronoiCell cell = voronoi.computeCell(0);

  EXPECT_NEAR(cell.volume, 1000.0, 1.0e-8);
  EXPECT_EQ(cell.faces.size(), 6);
  EXPECT_EQ(cell.verticesFractional.size(), 8);
  for (const SKVoronoiFace &face : cell.faces)
  {
    EXPECT_EQ(face.neighborIndex, 0);          // all faces come from periodic self-images
    EXPECT_NEAR(face.area, 100.0, 1.0e-8);
    EXPECT_NEAR(face.neighborDistance, 10.0, 1.0e-8);
  }
}

// A single site in an orthorhombic cell: the Voronoi cell is the full (non-cubic) box.
// A construction that ignored the metric would return a cell with equal face areas.
TEST(voronoi, single_site_orthorhombic_cell)
{
  double3x3 unitCell(double3(10.0, 0.0, 0.0), double3(0.0, 20.0, 0.0), double3(0.0, 0.0, 30.0));
  SKVoronoi voronoi(unitCell, {double3(0.1, 0.2, 0.3)});

  SKVoronoiCell cell = voronoi.computeCell(0);

  EXPECT_NEAR(cell.volume, 6000.0, 1.0e-8);
  EXPECT_EQ(cell.faces.size(), 6);

  std::multiset<int> distances;
  for (const SKVoronoiFace &face : cell.faces)
  {
    distances.insert(static_cast<int>(std::round(face.neighborDistance)));
  }
  EXPECT_EQ(distances, (std::multiset<int>{10, 10, 20, 20, 30, 30}));
}

// Two sites forming a bcc lattice: each Voronoi cell is the truncated octahedron with
// 8 hexagonal and 6 square faces and half the cell volume.
TEST(voronoi, bcc_truncated_octahedron)
{
  double a = 8.0;
  double3x3 unitCell(double3(a, 0.0, 0.0), double3(0.0, a, 0.0), double3(0.0, 0.0, a));
  SKVoronoi voronoi(unitCell, {double3(0.0, 0.0, 0.0), double3(0.5, 0.5, 0.5)});

  SKVoronoiCell cell = voronoi.computeCell(0);

  EXPECT_NEAR(cell.volume, 0.5 * a * a * a, 1.0e-8);
  EXPECT_EQ(cell.faces.size(), 14);

  std::size_t hexagons = 0;
  std::size_t squares = 0;
  for (const SKVoronoiFace &face : cell.faces)
  {
    if (face.vertexIndices.size() == 6) ++hexagons;
    if (face.vertexIndices.size() == 4) ++squares;
  }
  EXPECT_EQ(hexagons, 8);  // bisectors towards the 8 body-diagonal neighbors
  EXPECT_EQ(squares, 6);   // bisectors towards the 6 axial self-images
}

// The metric matters: in fractional coordinates the two sites are separated along x by 0.5,
// and the dividing face must lie at fractional x = 0.25 with area b*c.
TEST(voronoi, metric_aware_bisector_plane)
{
  double3x3 unitCell(double3(10.0, 0.0, 0.0), double3(0.0, 20.0, 0.0), double3(0.0, 0.0, 30.0));
  SKVoronoi voronoi(unitCell, {double3(0.0, 0.0, 0.0), double3(0.5, 0.0, 0.0)});

  SKVoronoiCell cell = voronoi.computeCell(0);
  EXPECT_NEAR(cell.volume, 0.5 * 6000.0, 1.0e-8);

  bool found = false;
  for (const SKVoronoiFace &face : cell.faces)
  {
    if (face.neighborIndex == 1 && face.neighborImage == int3(0, 0, 0))
    {
      found = true;
      EXPECT_NEAR(face.neighborDistance, 5.0, 1.0e-8);
      EXPECT_NEAR(face.area, 20.0 * 30.0, 1.0e-8);
      for (std::size_t index : face.vertexIndices)
      {
        EXPECT_NEAR(cell.verticesFractional[index].x, 0.25, 1.0e-8);
      }
    }
  }
  EXPECT_TRUE(found);
}

// The Voronoi cells of an arbitrary configuration tile the periodic cell: the volumes must
// sum to the cell volume. Uses a triclinic cell to exercise the full metric tensor.
TEST(voronoi, volumes_tile_triclinic_cell)
{
  double3x3 unitCell(double3(12.0, 0.0, 0.0), double3(3.0, 11.0, 0.0), double3(2.0, 1.5, 9.0));
  std::mt19937 generator(12345);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  std::vector<double3> positions;
  for (std::size_t i = 0; i < 24; ++i)
  {
    positions.push_back(double3(uniform(generator), uniform(generator), uniform(generator)));
  }

  SKVoronoi voronoi(unitCell, positions);
  std::vector<SKVoronoiCell> cells = voronoi.computeAllCells();

  double totalVolume = 0.0;
  for (const SKVoronoiCell &cell : cells)
  {
    EXPECT_GT(cell.volume, 0.0);
    totalVolume += cell.volume;
  }
  EXPECT_NEAR(totalVolume, voronoi.cellVolume(), 1.0e-6);
}

// Voronoi neighbor relations are symmetric: if j is a neighbor of i through image L, then
// i must be a neighbor of j through image -L, and the shared faces have equal area.
TEST(voronoi, faces_are_symmetric)
{
  double3x3 unitCell(double3(9.0, 0.0, 0.0), double3(0.0, 13.0, 0.0), double3(1.0, 2.0, 11.0));
  std::mt19937 generator(6789);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  std::vector<double3> positions;
  for (std::size_t i = 0; i < 12; ++i)
  {
    positions.push_back(double3(uniform(generator), uniform(generator), uniform(generator)));
  }

  SKVoronoi voronoi(unitCell, positions);
  std::vector<SKVoronoiCell> cells = voronoi.computeAllCells();

  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    for (const SKVoronoiFace &face : cells[i].faces)
    {
      bool found = false;
      for (const SKVoronoiFace &other : cells[face.neighborIndex].faces)
      {
        if (other.neighborIndex == i && other.neighborImage == -face.neighborImage)
        {
          EXPECT_NEAR(other.area, face.area, 1.0e-8);
          found = true;
        }
      }
      EXPECT_TRUE(found);
    }
  }
}
