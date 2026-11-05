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

module mc_opencl_pore_size_distribution;

import std;

import double3;
import float4;
import double3x3;
import randomnumbers;
import sampled_structure;
import opencl;

namespace
{
// A buffer is let go of however the scope is left, which an error thrown halfway through otherwise would not do.
struct ReleasedOnExit
{
  cl_mem memory;

  explicit ReleasedOnExit(cl_mem memory) : memory(memory) {}
  ReleasedOnExit(const ReleasedOnExit &) = delete;
  ReleasedOnExit &operator=(const ReleasedOnExit &) = delete;
  ~ReleasedOnExit() { clReleaseMemObject(memory); }
};
}  // namespace

MC_OpenCL_PoreSizeDistribution::MC_OpenCL_PoreSizeDistribution(std::size_t numberOfBins) : numberOfBins(numberOfBins),
             histogram(numberOfBins), histogram_cummulative(numberOfBins)
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* poreSizeDistributionShaderSourceCode = MC_OpenCL_PoreSizeDistribution::poreSizeDistributionKernelSource;
    poreSizeDistributionProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &poreSizeDistributionShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(poreSizeDistributionProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(poreSizeDistributionProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("MC_OpenCL_PoreSizeDistribution: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      clReleaseProgram(poreSizeDistributionProgram);
      throw std::runtime_error(message);
    }

    poreSizeDistributionKernel = clCreateKernel(poreSizeDistributionProgram, "ComputePoreSizeDistribution", &err);
    if (err != CL_SUCCESS)
    {
      clReleaseProgram(poreSizeDistributionProgram);
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(poreSizeDistributionKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &poreSizeDistributionWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      clReleaseKernel(poreSizeDistributionKernel);
      clReleaseProgram(poreSizeDistributionProgram);
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
};

MC_OpenCL_PoreSizeDistribution::~MC_OpenCL_PoreSizeDistribution()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    clReleaseKernel(poreSizeDistributionKernel);
    clReleaseProgram(poreSizeDistributionProgram);
  }
}


void MC_OpenCL_PoreSizeDistribution::run(const SampledStructure &structure,
                                         std::optional<std::size_t> numberOfIterations,
                                         std::optional<std::size_t> numberOfInnerSteps,
                                         std::optional<double> maximumRange)
{
  RandomNumber random{std::nullopt};
  cl_int err;
  std::chrono::steady_clock::time_point time_begin, time_end;

  std::size_t number_of_iterations = numberOfIterations.value_or(10000);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);

  // Up to a whole number of work groups, which need not be a power of two. The work items past the number asked for
  // find that out for themselves.
  size_t global_work_size = ((number_of_inner_steps + poreSizeDistributionWorkGroupSize - 1) / poreSizeDistributionWorkGroupSize) *
                            poreSizeDistributionWorkGroupSize;

  double3x3 unit_cell = structure.unitCell.cell;

  // The kernel is handed the cell's obtuse superbase and not the cell as given. It works in Cartesian positions
  // and uses a cell only to bring a displacement to its shortest image, which is found in that base; see
  // `UnitCell::shortestImage`. The two are one lattice. The points are still drawn in the cell as given.
  double3x3 image_cell = structure.unitCell.obtuseCell;
  double3x3 inverse_image_cell = structure.unitCell.obtuseInverse;

  std::size_t numberOfAtoms = structure.size();

  // There is no buffer of no size, so a box with nothing in it has one entry that nothing reads.
  std::vector<cl_float4> framework_positions(std::max<std::size_t>(numberOfAtoms, 1));
  std::vector<cl_float> sigma(std::max<std::size_t>(numberOfAtoms, 1));


  time_begin = std::chrono::steady_clock::now();

  for (size_t i = 0; i < numberOfAtoms; i++)
  {
    double3 position = structure.positions[i];
    framework_positions[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

    sigma[i] = cl_float(structure.radii[i]);
  }


  // upload framework position array
  cl_mem framework_positions_mem =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * framework_positions.size(), nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
  }
  ReleasedOnExit releaseFrameworkPositions(framework_positions_mem);

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), framework_positions_mem, CL_TRUE, 0, sizeof(float4) * framework_positions.size(),
                               framework_positions.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
  }


  // upload equilibrium size array
  cl_mem inputSigma =
      clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
  }
  ReleasedOnExit releaseSigma(inputSigma);

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(),
                             sigma.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
  }



  // setup unit cell matrix
  cl_float4 clCella = {{cl_float(image_cell[0][0]), cl_float(image_cell[1][0]), cl_float(image_cell[2][0]), cl_float(0.0)}};
  cl_float4 clCellb = {{cl_float(image_cell[0][1]), cl_float(image_cell[1][1]), cl_float(image_cell[2][1]), cl_float(0.0)}};
  cl_float4 clCellc = {{cl_float(image_cell[0][2]), cl_float(image_cell[1][2]), cl_float(image_cell[2][2]), cl_float(0.0)}};

  cl_float4 inverse_clCella = {{cl_float(inverse_image_cell[0][0]), cl_float(inverse_image_cell[1][0]), cl_float(inverse_image_cell[2][0]), cl_float(0.0)}};
  cl_float4 inverse_clCellb = {{cl_float(inverse_image_cell[0][1]), cl_float(inverse_image_cell[1][1]), cl_float(inverse_image_cell[2][1]), cl_float(0.0)}};
  cl_float4 inverse_clCellc = {{cl_float(inverse_image_cell[0][2]), cl_float(inverse_image_cell[1][2]), cl_float(inverse_image_cell[2][2]), cl_float(0.0)}};

  cl_int number_of_atoms = cl_int(numberOfAtoms);
  cl_int number_of_inner_steps_parameter = cl_int(number_of_inner_steps);

  double delta_r = maximumRange.value_or(10.0) / static_cast<double>(numberOfBins);

  // A second run is a run of its own and does not add to the first.
  std::ranges::fill(histogram, 0.0);
  std::ranges::fill(histogram_cummulative, 0.0);

  for(std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double3 fractional_position = {random.uniform(), random.uniform(), random.uniform()};
    double3 cartesian_position =  unit_cell * fractional_position;
    cl_float4 posA = {{cl_float(cartesian_position.x), cl_float(cartesian_position.y), cl_float(cartesian_position.z), 0.0f}};

    // A point of the void is covered by the sphere about itself if by no other, so that is where its size starts,
    // and a point in an atom is no part of the curve and needs no device to say so.
    std::optional<double> ownRadius = structure.freeRadius(cartesian_position);
    if (!ownRadius.has_value()) continue;
    const float seedRadius = static_cast<float>(std::min(ownRadius.value(), static_cast<double>(std::numeric_limits<float>::max())));

    float largest_radius[1] = {seedRadius};

    // With no centres to draw the sphere about the point itself is all there is.
    if (global_work_size > 0)
    {
      // upload random Cartesian positions
      std::size_t number_of_random_position{global_work_size};
      std::vector<cl_float4> random_cartesian_positions(number_of_random_position);
      for (size_t j = 0; j < number_of_random_position; j++)
      {
        fractional_position = {random.uniform(), random.uniform(), random.uniform()};
        cartesian_position =  unit_cell * fractional_position;
        random_cartesian_positions[j] = {{cl_float(cartesian_position.x), cl_float(cartesian_position.y), cl_float(cartesian_position.z), 0.0f}};
      }
      cl_mem random_cartesian_positions_mem =
            clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) *  random_cartesian_positions.size(), nullptr, &err);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
      }
      ReleasedOnExit releaseRandomPositions(random_cartesian_positions_mem);

      err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), random_cartesian_positions_mem, CL_TRUE, 0, sizeof(float4) * random_cartesian_positions.size(),
                                   random_cartesian_positions.data(), 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
      }

      // setup output array
      cl_mem largest_radius_mem =
            clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float), nullptr, &err);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
      }
      ReleasedOnExit releaseLargestRadius(largest_radius_mem);

      largest_radius[0] = seedRadius;
      err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), largest_radius_mem, CL_TRUE, 0,
                                 sizeof(cl_float), largest_radius, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
      }

      err = clSetKernelArg(poreSizeDistributionKernel,   0, sizeof(cl_mem), &framework_positions_mem);
      err |= clSetKernelArg(poreSizeDistributionKernel,  1, sizeof(cl_mem), &random_cartesian_positions_mem);
      err |= clSetKernelArg(poreSizeDistributionKernel,  2, sizeof(cl_mem), &inputSigma);
      err |= clSetKernelArg(poreSizeDistributionKernel,  3, sizeof(cl_mem), &largest_radius_mem);
      err |= clSetKernelArg(poreSizeDistributionKernel,  4, sizeof(cl_float4), &posA);
      err |= clSetKernelArg(poreSizeDistributionKernel,  5, sizeof(cl_int), &number_of_atoms);
      err |= clSetKernelArg(poreSizeDistributionKernel,  6, sizeof(cl_int), &number_of_inner_steps_parameter);
      err |= clSetKernelArg(poreSizeDistributionKernel,  7, sizeof(cl_float4), &clCella);
      err |= clSetKernelArg(poreSizeDistributionKernel,  8, sizeof(cl_float4), &clCellb);
      err |= clSetKernelArg(poreSizeDistributionKernel,  9, sizeof(cl_float4), &clCellc);
      err |= clSetKernelArg(poreSizeDistributionKernel, 10, sizeof(cl_float4), &inverse_clCella);
      err |= clSetKernelArg(poreSizeDistributionKernel, 11, sizeof(cl_float4), &inverse_clCellb);
      err |= clSetKernelArg(poreSizeDistributionKernel, 12, sizeof(cl_float4), &inverse_clCellc);

      err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), poreSizeDistributionKernel, 1, nullptr, &global_work_size,
                                    &poreSizeDistributionWorkGroupSize, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed (err: {}) {} : {}\n", err, __FILE__, __LINE__));
      }

      // Read the buffer back to the array
      err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), largest_radius_mem, CL_TRUE, 0, sizeof(float),
                                &largest_radius, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error("MC_OpenCL_PoreSizeDistribution: error in clEnqueueReadBuffer");
      }
    }

    if (!(delta_r > 0.0)) continue;

    // A size past the end of the range is drawn in no bin, and is at or above every one of them.
    const double bin = static_cast<double>(largest_radius[0]) / delta_r;
    const std::size_t index = bin < static_cast<double>(numberOfBins) ? static_cast<std::size_t>(bin) : numberOfBins;
    if(index < numberOfBins) histogram[index]++;
    for(std::size_t k = 0; k <= index && k < numberOfBins; ++k)
    {
      histogram_cummulative[k]++;
    }
  }
  clFinish(OpenCL::clCommandQueue.value());

  time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(structure.name + ".mc.psd.gpu.txt");
  std::print(myfile, "# Pore-size distribution using Mont Carlo-based method\n");
  structure.writeHeader(myfile);
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps (sample points per atom): {}\n", number_of_inner_steps);
  std::print(myfile, "# GPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "# column 1: diameter d [A]\n");
  std::print(myfile, "# column 2: pore size distribution, the share of the cell per unit of diameter [1/A]\n");
  std::print(myfile, "# column 3: cummulative pore volume\n");
  std::print(myfile, "# value at d=0 is related to the void-fraction\n");

  // The bins are of radius and the curve is against the diameter, so a bin is twice as wide as `delta_r` here.
  double normalization = 1.0 / static_cast<double>(number_of_iterations);
  double width = 2.0 * delta_r;
  for(std::size_t index = 0; index < numberOfBins; ++index)
  {
    std::print(myfile, "{} {} {}\n", width * (static_cast<double>(index) + 0.5),
               width > 0.0 ? histogram[index] * normalization / width : 0.0,
               histogram_cummulative[index] * normalization);
  }

  myfile.close();
}

const char* MC_OpenCL_PoreSizeDistribution::poreSizeDistributionKernelSource = R"foo(

// Function to perform the atomic max
// https://stackoverflow.com/questions/18950732/atomic-max-for-floats-in-opencl
//
inline void AtomicMax(volatile __global float *source, const float operand)
{
  union 
  {
    unsigned int intVal;
    float floatVal;
  } newVal;
  union
  {
    unsigned int intVal;
    float floatVal;
  } prevVal;
  do 
  {
     prevVal.floatVal = *source;
     newVal.floatVal = max(prevVal.floatVal,operand);
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

// The shortest image of a displacement, in Angstrom. The cell handed in is an obtuse superbase of the lattice:
// the image nearest in fractions of it to begin with, and then whichever of its three vectors, their sums in
// pairs and the sum of all three shortens it most is taken off, until none does. In a skewed cell the image
// nearest in fractions alone is not the shortest, and a sphere reaching more than half way across the cell was
// measured against the wrong copy of an atom.
float4 ShortestImage(float4 dr, float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc,
                     float4 cella, float4 cellb, float4 cellc)
{
  float4 shorteners[7];
  shorteners[0] = (float4)(cella.x, cellb.x, cellc.x, 0.0f);
  shorteners[1] = (float4)(cella.y, cellb.y, cellc.y, 0.0f);
  shorteners[2] = (float4)(cella.z, cellb.z, cellc.z, 0.0f);
  shorteners[3] = shorteners[0] + shorteners[1];
  shorteners[4] = shorteners[0] + shorteners[2];
  shorteners[5] = shorteners[1] + shorteners[2];
  shorteners[6] = shorteners[3] + shorteners[2];

  dr.w = 0.0f;
  dr -= rint(dot(inverse_cella, dr)) * shorteners[0] + rint(dot(inverse_cellb, dr)) * shorteners[1] +
        rint(dot(inverse_cellc, dr)) * shorteners[2];
  for(int step = 0; step < 16; step++)
  {
    int best = -1;
    float mostGained = 0.0f;
    float sign = 1.0f;
    for(int k = 0; k < 7; k++)
    {
      float halfSquare = 0.5f * dot(shorteners[k], shorteners[k]);
      float along = dot(dr, shorteners[k]);
      float gained = fabs(along) - halfSquare;
      if(gained > mostGained && gained > 1.0e-5f * halfSquare)
      {
        mostGained = gained;
        best = k;
        sign = along > 0.0f ? 1.0f : -1.0f;
      }
    }
    if(best < 0) break;
    dr -= sign * shorteners[best];
  }
  return dr;
}

// check overlap with other atoms
int CheckSphereOverLap(float4 posA, float *smallest_radius, int numberOfFrameworkAtoms, __global float4 *frameworkPositions, 
                       float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc,
                       float4 cella, float4 cellb, float4 cellc, __global float *radius_j)
{
  *smallest_radius = 1000.0;
  for(int jatom = 0; jatom < numberOfFrameworkAtoms; ++jatom)
  {
    float4 dr = ShortestImage(posA - frameworkPositions[jatom], inverse_cella, inverse_cellb, inverse_cellc,
                              cella, cellb, cellc);

    float rr = dot(dr, dr);

    if(rr < radius_j[jatom] * radius_j[jatom])
    {
      return true;
    }

    float radius = sqrt(rr) - radius_j[jatom];
    if(radius < *smallest_radius)
    {
      *smallest_radius = radius;
    }
  }
  return false;
}

__kernel void ComputePoreSizeDistribution(__global float4 *frameworkPositions,
                                          __global float4 *randomCartesianPositions,
                                          __global float *radius_j,
                                          __global float *largest_radius,
                                          const float4 posA,
                                          const int numberOfFrameworkAtoms,
                                          const int numberOfInnerSteps,
                                          const float4 cella,
                                          const float4 cellb,
                                          const float4 cellc,
                                          const float4 inverse_cella,
                                          const float4 inverse_cellb,
                                          const float4 inverse_cellc)
{
  int global_id = get_global_id(0);
  if(global_id >= numberOfInnerSteps) return;

  // The point itself is known to be in the void: the host looked before it came here.
  float radius = 0.0;
  float4 posB = randomCartesianPositions[global_id];

  if(CheckSphereOverLap(posB, &radius, numberOfFrameworkAtoms, frameworkPositions, 
                     inverse_cella, inverse_cellb, inverse_cellc,
                     cella, cellb, cellc, radius_j)) return;

  float4 dr = ShortestImage(posA - posB, inverse_cella, inverse_cellb, inverse_cellc, cella, cellb, cellc);

  float rr = dot(dr, dr);

  if(rr > radius * radius) return;

  AtomicMax(largest_radius, radius);
}
)foo";
