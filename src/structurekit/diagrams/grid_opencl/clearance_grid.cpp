module;

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module opencl_clearance_grid;

import std;

import opencl;
import int3;
import uint3;
import double2;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;


ClearanceGrid::ClearanceGrid() {}


ClearanceGrid::~ClearanceGrid() {}


double3 ClearanceGrid::fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const
{
  return double3(static_cast<double>(i) / static_cast<double>(this->gridSize.x),
                 static_cast<double>(j) / static_cast<double>(this->gridSize.y),
                 static_cast<double>(k) / static_cast<double>(this->gridSize.z));
}


double3 ClearanceGrid::cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const
{
  return this->unitCell.cell * this->fractionalPosition(i, j, k);
}


double ClearanceGrid::voxelVolume() const
{
  if (this->clearance.empty()) return 0.0;
  return this->unitCell.volume / static_cast<double>(this->clearance.size());
}


double ClearanceGrid::maximumClearance() const
{
  if (this->clearance.empty()) return 0.0;
  return static_cast<double>(*std::ranges::max_element(this->clearance));
}


void ClearanceGrid::writeTessellation(const Crystal &framework) const
{
  double3 spacing = this->spacing();
  const double voxel = this->voxelVolume();

  // How much of the cell falls to each atom. This is the useful part of a tessellation, and it is a count of
  // grid points, so it needs nothing but the labels themselves.
  std::vector<std::size_t> voxelsPerAtom(framework.atoms.size(), 0);
  for (std::int32_t atom : this->closestAtom)
  {
    if (atom >= 0 && static_cast<std::size_t>(atom) < voxelsPerAtom.size())
    {
      ++voxelsPerAtom[static_cast<std::size_t>(atom)];
    }
  }

  std::ofstream myfile;
  myfile.open(framework.name + ".grid.tessellation.gpu.txt");
  std::print(myfile, "# Tessellation of the cell by nearest atom surface (clearance grid)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", this->unitCell.volume);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", this->gridSize.x,
             this->gridSize.y, this->gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Volume one grid point stands for: {} [Å³]\n", voxel);
  std::print(myfile, "# GPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "#\n");
  std::print(myfile, "# A point belongs to the atom whose surface is nearest it, |x - p_i| - r_i, with r_i half\n");
  std::print(myfile, "# the atom's Lennard-Jones size parameter. These are the nearest-surface cells, the ones\n");
  std::print(myfile, "# the Apollonius diagram draws, and they are not the Voronoi cells: a large atom takes\n");
  std::print(myfile, "# ground from a small neighbour, and a boundary between two unequal atoms is curved. They\n");
  std::print(myfile, "# are also the cells this route's volume and area are divided along, so an atom's share\n");
  std::print(myfile, "# here is the share it is given there.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# column 1: atom\n");
  std::print(myfile, "# column 2: its radius [Å]\n");
  std::print(myfile, "# column 3: grid points nearer its surface than any other\n");
  std::print(myfile, "# column 4: the volume that comes to [Å³]\n");
  std::print(myfile, "# column 5: as a fraction of the cell\n");
  std::print(myfile, "#   atom   radius [Å]    points   volume [Å³]     fraction\n");
  for (std::size_t atom = 0; atom < voxelsPerAtom.size(); ++atom)
  {
    double share = static_cast<double>(voxelsPerAtom[atom]) * voxel;
    std::print(myfile, "CrystalAtom: {:6} {:11.5f} {:9} {:13.5f} {:12.8f}\n", atom,
               (atom < this->atomRadii.size()) ? this->atomRadii[atom] : 0.0, voxelsPerAtom[atom], share,
               share / this->unitCell.volume);
  }

  std::print(myfile, "\n");
  std::print(myfile, "# The grid itself: the atom nearest each point and the room there is at it, negative\n");
  std::print(myfile, "# inside an atom. Ordered with x varying fastest.\n");
  std::print(myfile, "#   i    j    k     atom  clearance [Å]\n");
  for (std::size_t k = 0; k < this->gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < this->gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < this->gridSize.x; ++i)
      {
        std::size_t voxelIndex = this->voxelIndex(i, j, k);
        std::print(myfile, "{:5} {:4} {:4} {:8} {:14.6f}\n", i, j, k, this->closestAtom[voxelIndex],
                   this->clearance[voxelIndex]);
      }
    }
  }

  myfile.close();
}


double3 ClearanceGrid::spacing() const
{
  double3 a = this->unitCell.cell[0];
  double3 b = this->unitCell.cell[1];
  double3 c = this->unitCell.cell[2];

  return double3(a.length() / static_cast<double>(this->gridSize.x), b.length() / static_cast<double>(this->gridSize.y),
                 c.length() / static_cast<double>(this->gridSize.z));
}


ClearanceGrid ClearanceGrid::compute(const PairInteractions &interactions, const Crystal &framework, uint3 gridSize)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("ClearanceGrid: no OpenCL device found, the grid route needs a GPU\n");
  }
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("ClearanceGrid: the grid must have at least one point along each axis\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  ClearanceGrid grid;
  grid.gridSize = gridSize;
  grid.unitCell = framework.unitCell;

  // Half the Lennard-Jones size parameter of the atom's own type, which is the radius the Voronoi and
  // Apollonius routes give an atom. No probe is folded in here: the field is the room for a probe centre and
  // a probe of any radius is applied afterwards by comparing against it.
  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  grid.atomRadii.reserve(framework.atoms.size());
  for (const CrystalAtom &atom : framework.atoms)
  {
    std::size_t type = atom.type;
    grid.atomRadii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.clearance.assign(numberOfVoxels, 0.0f);
  grid.closestAtom.assign(numberOfVoxels, -1);

  if (fractionalPositions.empty()) return grid;

  cl_int err;

  const char *source = ClearanceGrid::clearanceKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[2048];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    throw std::runtime_error(
        std::format("ClearanceGrid: OpenCL failed to build program at {} (error: {})\n", __LINE__,
                    std::string(buffer)));
  }

  cl_kernel kernel = clCreateKernel(program, "ClearanceField", &err);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  std::vector<cl_float4> positions(fractionalPositions.size());
  std::vector<cl_float> radii(fractionalPositions.size());
  for (std::size_t i = 0; i < fractionalPositions.size(); ++i)
  {
    positions[i] = {{cl_float(fractionalPositions[i].x), cl_float(fractionalPositions[i].y),
                     cl_float(fractionalPositions[i].z), 0.0f}};
    radii[i] = cl_float(grid.atomRadii[i]);
  }

  cl_int positionError, radiusError, clearanceError, atomError;
  cl_mem positionBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                         sizeof(cl_float4) * positions.size(), nullptr, &positionError);
  cl_mem radiusBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * radii.size(),
                                       nullptr, &radiusError);
  cl_mem clearanceBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY,
                                          sizeof(cl_float) * numberOfVoxels, nullptr, &clearanceError);
  cl_mem atomBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY, sizeof(cl_int) * numberOfVoxels,
                                     nullptr, &atomError);
  if (positionError != CL_SUCCESS || radiusError != CL_SUCCESS || clearanceError != CL_SUCCESS ||
      atomError != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clCreateBuffer failed at {}\n", __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), positionBuffer, CL_TRUE, 0,
                             sizeof(cl_float4) * positions.size(), positions.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), radiusBuffer, CL_TRUE, 0, sizeof(cl_float) * radii.size(),
                              radii.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
  }

  // The rows of the matrix that turns a fractional difference into a Cartesian one, so that the kernel can
  // build a distance with three dot products.
  double3x3 cell = framework.unitCell.cell;
  cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};

  // Wrapping a fractional difference into [-1/2, 1/2) picks the nearest image outright when the cell is
  // rectangular, because then each axis is independent. Once the cell is oblique it need not, so the
  // neighbouring images are visited too; one shell either way is enough for the reduced lattice vectors a
  // CIF supplies, and it costs nothing on the rectangular cells that do not need it.
  double longest = std::max({cell[0].length(), cell[1].length(), cell[2].length()});
  double skew = std::max({std::fabs(cell[0].y), std::fabs(cell[0].z), std::fabs(cell[1].x), std::fabs(cell[1].z),
                          std::fabs(cell[2].x), std::fabs(cell[2].y)});
  cl_int shell = (skew > 1.0e-8 * longest) ? 1 : 0;
  cl_int3 clShells = {{shell, shell, shell, 0}};

  // The least a fractional step along one axis can carry a point: the width of the cell perpendicular to the
  // plane of the other two axes. A fractional difference t is at least |t_k| w_k long in Cartesian terms
  // whatever the cell, because w_k is by definition the part of the axis that stands off that plane. The
  // kernel leans on this twice, once to dismiss a distant atom without a square root and once to settle the
  // nearest image without visiting the neighbouring cells; the fourth slot carries half the least of the
  // three, which is the threshold the second of those tests wants.
  double3 widths = framework.unitCell.perpendicularWidths();
  cl_float4 clWidths = {{cl_float(widths.x), cl_float(widths.y), cl_float(widths.z),
                         cl_float(0.5 * std::min({widths.x, widths.y, widths.z}))}};
  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  cl_int numberOfAtoms = static_cast<cl_int>(positions.size());

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &positionBuffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &radiusBuffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &clearanceBuffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &atomBuffer);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(kernel, 6, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_float4), &clWidths);
  err |= clSetKernelArg(kernel, 8, sizeof(cl_int), &numberOfAtoms);
  err |= clSetKernelArg(kernel, 9, sizeof(cl_int3), &clShells);
  err |= clSetKernelArg(kernel, 10, sizeof(cl_int3), &clGridSize);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                   static_cast<std::size_t>(gridSize.z)};
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 3, nullptr, globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), clearanceBuffer, CL_TRUE, 0,
                            sizeof(cl_float) * numberOfVoxels, grid.clearance.data(), 0, nullptr, nullptr);
  err |= clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), atomBuffer, CL_TRUE, 0, sizeof(cl_int) * numberOfVoxels,
                             grid.closestAtom.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ClearanceGrid: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(positionBuffer);
  clReleaseMemObject(radiusBuffer);
  clReleaseMemObject(clearanceBuffer);
  clReleaseMemObject(atomBuffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}


const char *ClearanceGrid::clearanceKernelSource = R"foo(
__kernel void ClearanceField(__global const float4 *position,
                             __global const float *radius,
                             __global float *clearance,
                             __global int *closestAtom,
                             const float4 cellRow0,
                             const float4 cellRow1,
                             const float4 cellRow2,
                             const float4 widths,
                             const int numberOfAtoms,
                             const int3 shells,
                             const int3 gridSize)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || iz >= gridSize.z) return;

  // Endpoint-exclusive sampling: fractional 0 and 1 are the same periodic point, so dividing by the grid
  // size rather than one less than it keeps the spacing uniform and every sample distinct.
  float4 s = (float4)((float)(ix) / (float)(gridSize.x),
                      (float)(iy) / (float)(gridSize.y),
                      (float)(iz) / (float)(gridSize.z),
                      0.0f);

  float best = MAXFLOAT;
  int bestAtom = -1;

  for(int iatom = 0; iatom < numberOfAtoms; iatom++)
  {
    float4 ds = s - position[iatom];
    ds -= rint(ds);
    ds.w = 0.0f;

    // Every image of this atom lies at least |ds_k| w_k away, for each of the three axes, since reducing into
    // [-1/2, 1/2) has already put the nearest integer at zero and stepping to another only lengthens the
    // fractional difference along that axis. An atom whose surface cannot come nearer than what has been
    // found already is passed over here, before the square root and before any image of it is built, and on a
    // framework of any size that is nearly all of them.
    float reach = fmax(fmax(fabs(ds.x) * widths.x, fabs(ds.y) * widths.y), fabs(ds.z) * widths.z);
    if(reach - radius[iatom] >= best) continue;

    float4 dr;
    dr.x = dot(cellRow0, ds);
    dr.y = dot(cellRow1, ds);
    dr.z = dot(cellRow2, ds);
    dr.w = 0.0f;
    float distance = sqrt(dot(dr, dr));

    float room = distance - radius[iatom];
    if(room < best)
    {
      best = room;
      bestAtom = iatom;
    }

    // Stepping to any other image moves some fractional component to half a cell or more, which by the same
    // bound costs at least half the least perpendicular width. So the reduced image is the nearest one
    // outright whenever it is nearer than that, and the neighbouring cells need not be looked at. Only an
    // atom most of a cell away fails this, and such an atom has almost always been dismissed above.
    if(distance < widths.w) continue;

    for(int a = -shells.x; a <= shells.x; a++)
    {
      for(int b = -shells.y; b <= shells.y; b++)
      {
        for(int c = -shells.z; c <= shells.z; c++)
        {
          if(a == 0 && b == 0 && c == 0) continue;

          float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);

          float4 image;
          image.x = dot(cellRow0, t);
          image.y = dot(cellRow1, t);
          image.z = dot(cellRow2, t);
          image.w = 0.0f;

          float other = sqrt(dot(image, image)) - radius[iatom];
          if(other < best)
          {
            best = other;
            bestAtom = iatom;
          }
        }
      }
    }
  }

  int index = (iz * gridSize.y + iy) * gridSize.x + ix;
  clearance[index] = best;
  closestAtom[index] = bestAtom;
}
)foo";
