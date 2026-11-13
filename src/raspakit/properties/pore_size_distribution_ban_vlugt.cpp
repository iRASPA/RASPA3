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

module pore_size_distribution_ban_vlugt;
    
import std;

import opencl;
import int3;
import uint3;
import double2;
import double3;
import double3x3;
import atom;
import framework;
import forcefield;

BanVlugtPoreSizeDistribution::BanVlugtPoreSizeDistribution(uint3 grid_size):
  grid_size(grid_size)
{
  OpenCL::initialize();
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* PSDShaderSourceCode = BanVlugtPoreSizeDistribution::PSDKernelSource;
    PSDProgram = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &PSDShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));}

    err = clBuildProgram(PSDProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(PSDProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      std::string message = std::format("Pore size distribution (Ban, Vlugt): OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__, __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    PSDKernel = clCreateKernel(PSDProgram, "PSD", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__ , __LINE__));}

    err = clGetKernelWorkGroupInfo(PSDKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(std::size_t), &PSDWorkGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));}
  }
}

BanVlugtPoreSizeDistribution::~BanVlugtPoreSizeDistribution()
{
  if (PSDKernel != nullptr) clReleaseKernel(PSDKernel);
  if (PSDProgram != nullptr) clReleaseProgram(PSDProgram);
}

void BanVlugtPoreSizeDistribution::run(const ForceField &forceField, const Framework &framework)
{
  const double3x3 unitCell = framework.simulationBox.cell;
  const double3 widths = framework.simulationBox.perpendicularWidths();
  const std::vector<double3> positions = framework.fractionalAtomPositionsUnitCell();
  const std::vector<double2> potentialParameters = framework.atomUnitCellLennardJonesPotentialParameters(forceField);
  std::chrono::steady_clock::time_point time_begin, time_end;
  cl_int err;

  const std::size_t nGrid = static_cast<std::size_t>(grid_size.x) * grid_size.y * grid_size.z;
  poreSizes.assign(nGrid, 0.0f);

  if(!OpenCL::clContext.has_value()) return;

  time_begin = std::chrono::steady_clock::now();

  if (positions.size() > 0 && nGrid > 0)
  {
    std::vector<cl_float4> pos(positions.size());
    std::vector<cl_float> radius(positions.size());
    for(size_t i = 0 ; i < positions.size(); ++i)
    {
      pos[i] = {{cl_float(positions[i].x), cl_float(positions[i].y), cl_float(positions[i].z), 0.0f}};

      // An atom is a hard sphere of its Lennard-Jones size across, so what is taken off a distance to its centre
      // is half of that. Taking the whole of it off left every pore the size of an atom too narrow.
      radius[i] = cl_float(0.5 * potentialParameters[i].y);
    }

    // The sizes are floats kept as their bits, which is what lets the device swap a wider one in atomically.
    std::vector<cl_int> pore_sizes_host(nGrid, std::bit_cast<cl_int>(0.0f));
    cl_mem pore_sizes = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                      sizeof(cl_int) * nGrid, pore_sizes_host.data(), &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputPos = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                    sizeof(cl_float4) * pos.size(), pos.data(), &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputRadius = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       sizeof(cl_float) * radius.size(), radius.data(), &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    // The rows of the unit cell, which carry a fractional displacement to a Cartesian one a component at a time,
    // and the distances between its faces, with half the least of them alongside: an image nearer than that is
    // the nearest outright.
    cl_float4 clRowX = {{cl_float(unitCell[0][0]), cl_float(unitCell[1][0]), cl_float(unitCell[2][0]), 0.0f}};
    cl_float4 clRowY = {{cl_float(unitCell[0][1]), cl_float(unitCell[1][1]), cl_float(unitCell[2][1]), 0.0f}};
    cl_float4 clRowZ = {{cl_float(unitCell[0][2]), cl_float(unitCell[1][2]), cl_float(unitCell[2][2]), 0.0f}};
    cl_float4 clWidths = {{cl_float(widths.x), cl_float(widths.y), cl_float(widths.z),
                           cl_float(0.5 * std::min({widths.x, widths.y, widths.z}))}};
    cl_int4 clGridSize = {{cl_int(grid_size.x), cl_int(grid_size.y), cl_int(grid_size.z), 0}};

    cl_int numberOfAtoms = static_cast<cl_int>(positions.size());
    err  = clSetKernelArg(PSDKernel, 0, sizeof(cl_mem), &inputPos);
    err |= clSetKernelArg(PSDKernel, 1, sizeof(cl_mem), &inputRadius);
    err |= clSetKernelArg(PSDKernel, 2, sizeof(cl_mem), &pore_sizes);
    err |= clSetKernelArg(PSDKernel, 3, sizeof(cl_float4), &clRowX);
    err |= clSetKernelArg(PSDKernel, 4, sizeof(cl_float4), &clRowY);
    err |= clSetKernelArg(PSDKernel, 5, sizeof(cl_float4), &clRowZ);
    err |= clSetKernelArg(PSDKernel, 6, sizeof(cl_float4), &clWidths);
    err |= clSetKernelArg(PSDKernel, 7, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(PSDKernel, 8, sizeof(cl_int4), &clGridSize);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clSetKernelArg failed {} : {}\n", __FILE__, __LINE__));}

    std::size_t global_work_size[3] = {static_cast<std::size_t>(grid_size.x), static_cast<std::size_t>(grid_size.y), static_cast<std::size_t>(grid_size.z)};
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), PSDKernel, 3, nullptr, global_work_size, nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel PSDKernel failed at {} line {}\n", __FILE__, __LINE__));}

    clFinish(OpenCL::clCommandQueue.value());

    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), pore_sizes, CL_TRUE, 0,
                              sizeof(cl_int) * nGrid, pore_sizes_host.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueReadBuffer failed {} : {}\n", __FILE__, __LINE__));}

    for (std::size_t i = 0; i < nGrid; ++i) poreSizes[i] = std::bit_cast<float>(pore_sizes_host[i]);

    clReleaseMemObject(inputPos);
    clReleaseMemObject(inputRadius);
    clReleaseMemObject(pore_sizes);
  }

  time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".PSD_ban_vlugt.gpu.txt");
  std::print(myfile, "# Pore size distribution using Ban, Vlugt method\n");
  std::print(myfile, "# Grid point i of n along an axis stands at the fraction i / n of the cell\n");
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  for (std::size_t k = 0; k < static_cast<std::size_t>(grid_size.z); ++k)
  {
    for (std::size_t j = 0; j < static_cast<std::size_t>(grid_size.y); ++j)
    {
      for (std::size_t i = 0; i < static_cast<std::size_t>(grid_size.x); ++i)
      {
        std::size_t idx = i + j * grid_size.x + k * grid_size.x * grid_size.y;
        std::print(myfile,"{} {} {} {}\n", i, j, k, poreSizes[idx]);
      }
    }
  }

  std::cout << "PSD completed in " << timing.count() << " seconds.\n";
  std::cout << "PSD results saved to " << framework.name << ".PSD_ban_vlugt.gpu.txt\n";

  myfile.close();

}

const char* BanVlugtPoreSizeDistribution::PSDKernelSource = R"foo(
// A displacement in fractions of the cell as a Cartesian one, through the rows of the cell matrix.
static float4 cartesianOf(float4 t, float4 cellRowX, float4 cellRowY, float4 cellRowZ)
{
    return (float4)(dot(cellRowX, t), dot(cellRowY, t), dot(cellRowZ, t), 0.0f);
}

__kernel void PSD(
    __global const float4 *position,   // fractional, in the unit cell
    __global const float *radius,      // of each atom as a hard sphere
    volatile __global int *pore_sizes, // 3D grid flattened: [z*Y*X + y*X + x], floats kept as their bits
    const float4 cellRowX,
    const float4 cellRowY,
    const float4 cellRowZ,
    const float4 widths,               // between the faces of the cell across each axis; .w half the least
    const int numberOfAtoms,
    const int4 grid_size)
{
    int ix = get_global_id(0);
    int iy = get_global_id(1);
    int iz = get_global_id(2);

    // The grid is of a periodic cell: the point after the last along an axis is the first again, so point i of n
    // stands at i / n and no point is on the far face.
    float4 grid_position = (float4)(
        (float)(ix) / (float)(grid_size.x),
        (float)(iy) / (float)(grid_size.y),
        (float)(iz) / (float)(grid_size.z),
        0.0f);

    // Step 1: the distance to the surface of the closest atom, by the closest copy of each.
    float closest_distance = 1e5f;
    for(int iatom = 0; iatom < numberOfAtoms; iatom++)
    {
        float4 t = grid_position - position[iatom];
        t -= rint(t);
        t.w = 0.0f;
        float size = radius[iatom];

        // Every copy of this atom is at least |t_k| times the width across axis k away, the nearest whole
        // number having been taken off already, so one that cannot come nearer than what is in hand is passed.
        float reach = fmax(fmax(fabs(t.x) * widths.x, fabs(t.y) * widths.y), fabs(t.z) * widths.z);
        if(reach - size >= closest_distance) continue;

        float distance = length(cartesianOf(t, cellRowX, cellRowY, cellRowZ));
        closest_distance = fmin(closest_distance, distance - size);

        // The copy nearest in fractions is the nearest in Angstrom if it is nearer than half the least width.
        // Past that, in a skewed cell or a narrow one, a copy a cell or more over can be nearer, and one nearer
        // than this one is no more than distance / width + 1/2 cells over along any axis.
        if(distance < widths.w) continue;
        int na = (int)floor(distance / widths.x + 0.5f);
        int nb = (int)floor(distance / widths.y + 0.5f);
        int nc = (int)floor(distance / widths.z + 0.5f);
        for(int a = -na; a <= na; a++)
        for(int b = -nb; b <= nb; b++)
        for(int c = -nc; c <= nc; c++)
        {
            if(a == 0 && b == 0 && c == 0) continue;
            float4 image = t + (float4)((float)(a), (float)(b), (float)(c), 0.0f);
            float other = length(cartesianOf(image, cellRowX, cellRowY, cellRowZ));
            closest_distance = fmin(closest_distance, other - size);
        }
    }

    // A point inside an atom is the centre of no sphere.
    if(!(closest_distance > 0.0f)) return;

    // Step 2: every grid point inside the sphere of that radius keeps it if it is the widest laid over it yet.
    //
    // A step of one grid point along an axis is at least the width across that axis over the number of points
    // along it in Angstrom, whatever the angles of the cell, so that many steps is as far as the sphere can reach.
    // The steps are of the unit cell the grid is of. Each is taken to Angstrom through the cell matrix, so that
    // what is compared with the radius is a distance and not a count of grid points, and the point it lands on is
    // found round the cell: a sphere at a face covers the points behind the opposite one.
    int num_points_x = (int)floor(closest_distance * (float)(grid_size.x) / widths.x);
    int num_points_y = (int)floor(closest_distance * (float)(grid_size.y) / widths.y);
    int num_points_z = (int)floor(closest_distance * (float)(grid_size.z) / widths.z);

    float radius_squared = closest_distance * closest_distance;
    int new_val = as_int(closest_distance);

    for (int dz_ = -num_points_z; dz_ <= num_points_z; dz_++) {
        int zz = ((iz + dz_) % grid_size.z + grid_size.z) % grid_size.z;
        for (int dy_ = -num_points_y; dy_ <= num_points_y; dy_++) {
            int yy = ((iy + dy_) % grid_size.y + grid_size.y) % grid_size.y;
            for (int dx_ = -num_points_x; dx_ <= num_points_x; dx_++) {
                int xx = ((ix + dx_) % grid_size.x + grid_size.x) % grid_size.x;

                float4 step = (float4)(
                    (float)(dx_) / (float)(grid_size.x),
                    (float)(dy_) / (float)(grid_size.y),
                    (float)(dz_) / (float)(grid_size.z),
                    0.0f);
                float4 offset = cartesianOf(step, cellRowX, cellRowY, cellRowZ);
                if (dot(offset, offset) > radius_squared) continue;

                int neighbor_idx = zz * grid_size.y * grid_size.x + yy * grid_size.x + xx;

                // Atomically update the pore size if this value is larger
                int seen = pore_sizes[neighbor_idx];
                while (closest_distance > as_float(seen)) {
                    int found = atomic_cmpxchg(&pore_sizes[neighbor_idx], seen, new_val);
                    if (found == seen) break;
                    seen = found;
                }
            }
        }
    }
}
)foo";
