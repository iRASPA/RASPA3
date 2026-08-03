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
  const float threshold = static_cast<float>(this->probeRadius);

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
    }
    else
    {
      this->inaccessibleSurfaceArea += area;
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
  double3 spacing = grid.spacing();

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
