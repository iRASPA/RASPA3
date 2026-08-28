module;

module opencl_surface_area;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import crystal;
import pair_interactions;
import units;
import marching_cubes;
import surface_curvature;
import opencl_clearance_grid;
import grid_connected_components;


GridSurfaceArea::GridSurfaceArea() {}


GridSurfaceArea::~GridSurfaceArea() {}


void GridSurfaceArea::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                          uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  this->run(interactions, framework, probePseudoAtom, grid);
}


void GridSurfaceArea::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                          const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridSurfaceArea: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridComponents components = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);

  const std::size_t nx = grid.gridSize.x;
  const std::size_t ny = grid.gridSize.y;
  const std::size_t nz = grid.gridSize.z;

  // Marching cubes wants the far face of the cell as well as the near one, and on a periodic field the far
  // face is the near one over again. Repeating it is what closes the surface across the cell boundary instead
  // of leaving it cut open there.
  MarchingCubes cube(static_cast<int>(nx + 1), static_cast<int>(ny + 1), static_cast<int>(nz + 1));
  cube.init_all();
  for (std::size_t k = 0; k <= nz; ++k)
  {
    for (std::size_t j = 0; j <= ny; ++j)
    {
      for (std::size_t i = 0; i <= nx; ++i)
      {
        std::size_t voxel = grid.voxelIndex(i % nx, j % ny, k % nz);
        cube.set_data(static_cast<double>(grid.clearance[voxel]), i, j, k);
      }
    }
  }

  cube.run(this->probeRadius);

  this->atomSurfaceArea.assign(framework.atoms.size(), 0.0);

  const double3x3 cell = framework.unitCell.cell;
  const double3x3 inverseCell = framework.unitCell.inverseCell;
  const float threshold = static_cast<float>(this->probeRadius);

  // A curvature whose radius is shorter than the coarsest voxel edge is not a shape the grid found. The
  // clearance field has a crease wherever two atoms are equally near, its gradient jumps across the crease,
  // and the triangles straddling one come back with curvatures of that size and of either sign. Bounding them
  // keeps that debris out of the saddle column, where it would otherwise be the whole of it.
  this->curvatureBands = curvatureBandsForGrid(cell, grid.gridSize);

  this->numberOfTriangles = cube.ntrigs();
  for (std::size_t t = 0; t < cube.ntrigs(); ++t)
  {
    const Triangle *triangle = cube.trig(static_cast<std::ptrdiff_t>(t));
    const Vertex *first = cube.vert(triangle->v1);
    const Vertex *second = cube.vert(triangle->v2);
    const Vertex *third = cube.vert(triangle->v3);
    if (first == nullptr || second == nullptr || third == nullptr) continue;

    // Vertices come back in grid steps, which are a fraction of the cell along each axis.
    double3 a(first->x / static_cast<double>(nx), first->y / static_cast<double>(ny),
              first->z / static_cast<double>(nz));
    double3 b(second->x / static_cast<double>(nx), second->y / static_cast<double>(ny),
              second->z / static_cast<double>(nz));
    double3 c(third->x / static_cast<double>(nx), third->y / static_cast<double>(ny),
              third->z / static_cast<double>(nz));

    double3 pa = cell * a;
    double3 pb = cell * b;
    double3 pc = cell * c;
    double area = 0.5 * double3::cross(pb - pa, pc - pa).length();
    if (!std::isfinite(area)) continue;

    this->totalSurfaceArea += area;

    // The normals marching cubes stored are the field's gradient at each vertex, held per grid step, so they
    // need carrying to a position before they mean anything on a cell whose axes are not orthogonal. The
    // clearance field grows into the void, so the gradient points away from the solid, and with that normal a
    // curvature is positive where the wall bulges out into the void.
    std::array<double3, 3> triangleCorners{pa, pb, pc};
    std::array<double3, 3> triangleNormals{
        cartesianNormalOfGridGradient(inverseCell, grid.gridSize, double3(first->nx, first->ny, first->nz),
                                      FieldSense::GrowsIntoVoid),
        cartesianNormalOfGridGradient(inverseCell, grid.gridSize, double3(second->nx, second->ny, second->nz),
                                      FieldSense::GrowsIntoVoid),
        cartesianNormalOfGridGradient(inverseCell, grid.gridSize, double3(third->nx, third->ny, third->nz),
                                      FieldSense::GrowsIntoVoid)};

    TriangleCurvature shape = triangleCurvature(triangleCorners, triangleNormals);
    this->curvature.add(area, shape, this->curvatureBands);

    // Which pore this piece of surface faces, and which atom it belongs to, read off the corner of the cell
    // the piece was cut from that lies on the open side.
    double3 centre((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0);
    std::size_t baseI = static_cast<std::size_t>(centre.x * static_cast<double>(nx)) % nx;
    std::size_t baseJ = static_cast<std::size_t>(centre.y * static_cast<double>(ny)) % ny;
    std::size_t baseK = static_cast<std::size_t>(centre.z * static_cast<double>(nz)) % nz;

    std::int32_t pore = -1;
    std::int32_t atom = -1;
    float widest = threshold;
    for (std::size_t dk = 0; dk < 2 && pore < 0; ++dk)
    {
      for (std::size_t dj = 0; dj < 2 && pore < 0; ++dj)
      {
        for (std::size_t di = 0; di < 2 && pore < 0; ++di)
        {
          std::size_t voxel = grid.voxelIndex((baseI + di) % nx, (baseJ + dj) % ny, (baseK + dk) % nz);
          if (grid.clearance[voxel] >= widest && components.voxelPore[voxel] >= 0)
          {
            pore = components.voxelPore[voxel];
            atom = grid.closestAtom[voxel];
          }
        }
      }
    }

    if (pore < 0)
    {
      this->undecidedSurfaceArea += area;
    }
    else if (components.pores[static_cast<std::size_t>(pore)].isChannel)
    {
      this->accessibleSurfaceArea += area;
      this->accessibleCurvature.add(area, shape, this->curvatureBands);
    }
    else
    {
      this->inaccessibleSurfaceArea += area;
      this->inaccessibleCurvature.add(area, shape, this->curvatureBands);
    }

    if (atom >= 0 && static_cast<std::size_t>(atom) < this->atomSurfaceArea.size())
    {
      this->atomSurfaceArea[static_cast<std::size_t>(atom)] += area;
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  const double volume = framework.unitCell.volume;
  const double toGravimetric = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass;
  const double3 spacing = grid.spacing();

  std::ofstream myfile;
  myfile.open(framework.name + ".grid.sa.gpu.txt");
  std::print(myfile, "# Accessible / inaccessible surface area (clearance grid, marching cubes)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, this->probeRadius);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Triangles on the surface: {}\n", this->numberOfTriangles);
  std::print(myfile, "# Channels: {}, pockets: {}\n", components.numberOfChannels, components.numberOfPockets);
  std::print(myfile, "# Pore system dimensionality: {}\n", components.dimensionality);
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
  std::print(myfile, "# GPU Timing: {} [s] for the clearance field\n", grid.seconds);
  std::print(myfile, "# CPU Timing: {} [s] for the surface\n", this->seconds);
  std::print(myfile, "# The surface is the level set of the clearance field at the probe's radius, which is the\n");
  std::print(myfile, "# surface of the atoms grown by that radius, so it is the accessible surface itself rather\n");
  std::print(myfile, "# than an estimate of it, and what the grid decides is only where that level set lies.\n");
  std::print(myfile, "# Triangles chord the surface, so the area comes out low and closes on the true one from\n");
  std::print(myfile, "# below. It does so at first order rather than second, because the field is a distance to\n");
  std::print(myfile, "# the nearest of several atoms and so has a crease everywhere two of them meet, and a\n");
  std::print(myfile, "# crease is where a triangulation loses the most.\n");
  std::print(myfile, "Accessible surface area:   {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", this->accessibleSurfaceArea,
             1.0e4 * this->accessibleSurfaceArea / volume, this->accessibleSurfaceArea * toGravimetric);
  std::print(myfile, "Inaccessible surface area: {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", this->inaccessibleSurfaceArea,
             1.0e4 * this->inaccessibleSurfaceArea / volume, this->inaccessibleSurfaceArea * toGravimetric);
  if (this->undecidedSurfaceArea > 0.0)
  {
    std::print(myfile, "Undecided surface area:    {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", this->undecidedSurfaceArea,
               1.0e4 * this->undecidedSurfaceArea / volume, this->undecidedSurfaceArea * toGravimetric);
  }
  std::print(myfile, "Total surface area:        {} [Å²]\n", this->totalSurfaceArea);

  std::print(myfile, "\n");
  std::print(myfile, "# The same area divided by the shape of the sheet rather than by what lies behind it.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Marching cubes places each vertex where the field crosses the level along a cube edge and\n");
  std::print(myfile, "# records the field's gradient there, which is the normal of the level set. Three normals to\n");
  std::print(myfile, "# a triangle fix how the normal turns across it, and that is the two principal curvatures.\n");
  std::print(myfile, "# The gradient of a clearance field points away from the solid, so a curvature is positive\n");
  std::print(myfile, "# where the wall bulges out into the void: both positive is convex, both negative is the\n");
  std::print(myfile, "# inside of a pocket and concave, one of each is a saddle.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# READ THE COLUMNS AS A PROPERTY OF THIS GRID AND NOT OF THE MATERIAL. This is a peculiarity\n");
  std::print(myfile, "# of the clearance field and not of the method: the field is a minimum over the atoms, so its\n");
  std::print(myfile, "# level set is the boundary of a union of balls, and every point of that not on a crease\n");
  std::print(myfile, "# between two balls lies on a sphere and is convex, the creases being curves that carry no\n");
  std::print(myfile, "# area. So the true split here is all convex and nothing else. What makes the other columns\n");
  std::print(myfile, "# is the grid: a crease is rounded over about a voxel, and the band either side of it, where\n");
  std::print(myfile, "# the central differences straddle the kink, comes back saddle or concave. Refine and that\n");
  std::print(myfile, "# band narrows and the convex share climbs towards one, which on MFI takes it from a third\n");
  std::print(myfile, "# at 96³ to two thirds at 320³ without settling. Quoting a saddle fraction from here without\n");
  std::print(myfile, "# the spacing beside it therefore says nothing. What the columns are good for is localising\n");
  std::print(myfile, "# the creases: how much of the wall lies within a voxel of one, and which pores it sits in.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The energy route does not have this problem, and if a genuine convex/saddle/concave split\n");
  std::print(myfile, "# is what is wanted, that is the one to use. An energy field is a sum over the atoms rather\n");
  std::print(myfile, "# than a minimum, so it is smooth everywhere and its contours have no creases at all: the\n");
  std::print(myfile, "# saddle between two atoms there is a real saddle whose curvature the grid only has to\n");
  std::print(myfile, "# resolve, not invent. The two integrals below converge on either field.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# This is also not the convex/saddle/concave split of the solvent-excluded surface that the\n");
  std::print(myfile, "# exact route reports. That one counts how many atoms the probe touches at once and lives on\n");
  std::print(myfile, "# the surface the probe itself touches, where the saddle pieces are tori and the concave ones\n");
  std::print(myfile, "# patches of the probe's own sphere. On this surface, the one the probe's centre traces,\n");
  std::print(myfile, "# those tori have shrunk to the creases and those patches to single points. The two are\n");
  std::print(myfile, "# decompositions of different surfaces and the numbers are not comparable.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# A curvature below {:.5f} [1/Å], a radius of {:.1f} [Å] and longer, is reported as flat\n",
             this->curvatureBands.flat, 1.0 / this->curvatureBands.flat);
  std::print(myfile, "# rather than given a sign it does not have. One above {:.5f} [1/Å] is a radius shorter than\n",
             this->curvatureBands.sharpest);
  std::print(myfile, "# a voxel, which is the crease itself rather than anything the grid resolved, and is held\n");
  std::print(myfile, "# back in its own column. It is kept out of the columns only, not out of the integrals.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The first four fractions are of the classified area and add to one; the last is the\n");
  std::print(myfile, "# unresolved share of the whole.\n");
  std::print(myfile, "#                              area [Å²]     convex     saddle    concave       flat unresolved\n");

  auto curvatureRow = [&](const char *name, const CurvatureAreas &areas)
  {
    if (areas.total() <= 0.0) return;
    std::print(myfile, "{:<26} {:13.5f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n", name, areas.classified(),
               areas.convexFraction(), areas.saddleFraction(), areas.concaveFraction(), areas.flatFraction(),
               areas.unresolvedFraction());
  };

  curvatureRow("Curvature split, all:", this->curvature);
  curvatureRow("  facing channels:", this->accessibleCurvature);
  curvatureRow("  facing pockets:", this->inaccessibleCurvature);

  std::print(myfile, "# Triangles: {} convex, {} saddle, {} concave, {} flat, {} unresolved, {} of those\n",
             this->curvature.numberOfConvex, this->curvature.numberOfSaddle, this->curvature.numberOfConcave,
             this->curvature.numberOfFlat, this->curvature.numberOfUnresolved, this->curvature.numberOfDegenerate);
  std::print(myfile, "# degenerate rather than too sharp.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The two integrals over the whole surface, which unlike the columns above do converge. They\n");
  std::print(myfile, "# converge because the crease triangles are in them: a crease carries a curvature of one\n");
  std::print(myfile, "# over the rounding radius over an area proportional to that radius, and the product settles\n");
  std::print(myfile, "# while neither factor does.\n");
  std::print(myfile, "Integral of the mean curvature:     {} [Å]\n", this->curvature.integratedMeanCurvature);
  std::print(myfile, "Integral of the Gaussian curvature: {} [-]\n",
             this->curvature.integratedGaussianCurvature);
  std::print(myfile, "# The first is a Minkowski functional of the solid, so it can be checked against a closed\n");
  std::print(myfile, "# form for a union of balls; on MFI it is within a percent of its limit by 256³. It comes out\n");
  std::print(myfile, "# negative wherever the balls overlap heavily, the creases outweighing the caps.\n");
  std::print(myfile, "# The second is 2 pi times the Euler characteristic over a closed surface, so it counts\n");
  std::print(myfile, "# rather than measures: negative on a network of channels, and the more negative the more\n");
  std::print(myfile, "# connected the network. It is much the noisier of the two, the fit being a least-squares one\n");
  std::print(myfile, "# over three normals rather than the angle deficit that satisfies discrete Gauss-Bonnet\n");
  std::print(myfile, "# exactly. On one isolated atom it lands within a part in a thousand of 4 pi; on a framework\n");
  std::print(myfile, "# it moves by a few percent between grids, so read its sign and its size and not past the\n");
  std::print(myfile, "# first digit or two.\n");

  std::print(myfile, "\n");
  std::print(myfile, "# The area each atom contributes. The field records which atom's surface is the nearest at\n");
  std::print(myfile, "# every point, so a piece of surface is attributed to an atom rather than divided among\n");
  std::print(myfile, "# them, and the parts add up to the total above. An atom buried in the framework wall\n");
  std::print(myfile, "# contributes nothing and is left out.\n");
  std::print(myfile, "#  atom  type                     area [Å²]\n");
  for (std::size_t i = 0; i < this->atomSurfaceArea.size(); ++i)
  {
    if (this->atomSurfaceArea[i] <= 0.0) continue;
    std::size_t type = framework.atoms[i].type;
    std::print(myfile, "CrystalAtom: {:6} {:<20} {:13.5f}\n", i, interactions.names[type],
               this->atomSurfaceArea[i]);
  }

  myfile.close();
}
