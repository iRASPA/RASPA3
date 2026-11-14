module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module mc_opencl_surface_area;

import std;

import opencl;
import double3;
import float4;
import double3x3;
import randomnumbers;
import sampled_structure;
import unit_cell;

MC_OpenCL_SurfaceArea::MC_OpenCL_SurfaceArea()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* surfaceAreaShaderSourceCode = MC_OpenCL_SurfaceArea::surfaceAreaKernelSource;
    surfaceAreaProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &surfaceAreaShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(surfaceAreaProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(surfaceAreaProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("MC_OpenCL_SurfaceArea: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    surfaceAreaKernel = clCreateKernel(surfaceAreaProgram, "ComputeSurfaceArea", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(surfaceAreaKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &surfaceAreaWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
}

MC_OpenCL_SurfaceArea::~MC_OpenCL_SurfaceArea()
{
  if (surfaceAreaKernel != nullptr) clReleaseKernel(surfaceAreaKernel);
  if (surfaceAreaProgram != nullptr) clReleaseProgram(surfaceAreaProgram);
}

void MC_OpenCL_SurfaceArea::run(const SampledStructure &structure, const SampledProbe &probe,
                                std::optional<std::size_t> numberOfIterations,
                                std::optional<std::size_t> numberOfInnerSteps)
{
  RandomNumber random{std::nullopt};
  // A route made where there was no device has no kernel to run.
  if (!OpenCL::clContext.has_value() || !OpenCL::clCommandQueue.has_value() || surfaceAreaKernel == nullptr ||
      surfaceAreaWorkGroupSize == 0)
  {
    throw std::runtime_error("MC_OpenCL_SurfaceArea: no OpenCL device found\n");
  }

  cl_int err;
  std::chrono::steady_clock::time_point time_begin, time_end;

  // The buffers of this run, let go of when the run is left, by its last line or by a failure before it. They
  // used to be let go of by neither, and every run left three of them on the device.
  OpenCL::Handles made;

  // the defaults of MC_SurfaceArea, so that the two routes sample the same number of points unless asked otherwise
  std::size_t number_of_iterations = numberOfIterations.value_or(100);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(1000);
  // the estimate divides by both, and an empty buffer cannot be made on the device
  if (number_of_iterations == 0 || number_of_inner_steps == 0)
  {
    throw std::invalid_argument(
        "MC_OpenCL_SurfaceArea: the number of iterations and of inner steps must be at least one\n");
  }

  double3x3 unit_cell = structure.unitCell.cell;
  double3x3 inverse_unit_cell = structure.unitCell.inverseCell;

  std::size_t numberOfAtoms = structure.size();
  size_t global_work_size =
      ((numberOfAtoms + surfaceAreaWorkGroupSize - 1) / surfaceAreaWorkGroupSize) * surfaceAreaWorkGroupSize;

  std::vector<cl_float4> pos(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);

  std::vector<cl_float> output(numberOfAtoms);

  time_begin = std::chrono::steady_clock::now();

  for (size_t i = 0; i < numberOfAtoms; i++)
  {
    double3 position = structure.positions[i];
    pos[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

    sigma[i] = cl_float(structure.radii[i]);
  }

  double sum{};
  double sum_of_squares{};

  // A structure of no atoms has no surface, and the device is not asked about it: it will not make a buffer of
  // no size, and there is nothing for it to work out. What follows the block is the same either way.
  if (numberOfAtoms > 0)
  {
    // upload position array
    cl_mem inputPos = made.keep(
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * pos.size(), nullptr, &err));
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputPos, CL_TRUE, 0, sizeof(float4) * pos.size(),
                                 pos.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    // upload equilibrium size array
    cl_mem inputSigma = made.keep(
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err));
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(),
                               sigma.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    // setup output array
    cl_mem outputArray = made.keep(
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * output.size(), nullptr, &err));
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputArray, CL_TRUE, 0,
                               sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
    }


    // setup unit cell matrix
    cl_float4 clCella = {{cl_float(unit_cell[0][0]), cl_float(unit_cell[1][0]), cl_float(unit_cell[2][0]), cl_float(0.0)}};
    cl_float4 clCellb = {{cl_float(unit_cell[0][1]), cl_float(unit_cell[1][1]), cl_float(unit_cell[2][1]), cl_float(0.0)}};
    cl_float4 clCellc = {{cl_float(unit_cell[0][2]), cl_float(unit_cell[1][2]), cl_float(unit_cell[2][2]), cl_float(0.0)}};

    cl_float4 inverse_clCella = {{cl_float(inverse_unit_cell[0][0]), cl_float(inverse_unit_cell[1][0]), cl_float(inverse_unit_cell[2][0]), cl_float(0.0)}};
    cl_float4 inverse_clCellb = {{cl_float(inverse_unit_cell[0][1]), cl_float(inverse_unit_cell[1][1]), cl_float(inverse_unit_cell[2][1]), cl_float(0.0)}};
    cl_float4 inverse_clCellc = {{cl_float(inverse_unit_cell[0][2]), cl_float(inverse_unit_cell[1][2]), cl_float(inverse_unit_cell[2][2]), cl_float(0.0)}};

    // The widths across the cell, which bound the images of an atom that can reach a point.
    double3 widths = structure.unitCell.perpendicularWidths();
    cl_float4 inverseWidths = {{cl_float(1.0 / widths.x), cl_float(1.0 / widths.y), cl_float(1.0 / widths.z), cl_float(0.0)}};

    std::size_t number_of_random_unit_vectors{number_of_inner_steps};
    std::vector<cl_float4> random_unit_vectors(number_of_random_unit_vectors);
    cl_mem random_unit_vectors_mem = made.keep(clCreateBuffer(
        OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * random_unit_vectors.size(), nullptr, &err));
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_int end_index = cl_int(numberOfAtoms);
    cl_int number_of_slices = cl_int(number_of_inner_steps);

    err = clSetKernelArg(surfaceAreaKernel, 0, sizeof(cl_mem), &inputPos);
    err |= clSetKernelArg(surfaceAreaKernel, 1, sizeof(cl_mem), &inputSigma);
    err |= clSetKernelArg(surfaceAreaKernel, 2, sizeof(cl_mem), &random_unit_vectors_mem);
    err |= clSetKernelArg(surfaceAreaKernel, 3, sizeof(cl_mem), &outputArray);
    err |= clSetKernelArg(surfaceAreaKernel, 4, sizeof(cl_int), &end_index);
    err |= clSetKernelArg(surfaceAreaKernel, 5, sizeof(cl_int), &number_of_slices);
    err |= clSetKernelArg(surfaceAreaKernel, 6, sizeof(cl_float4), &clCella);
    err |= clSetKernelArg(surfaceAreaKernel, 7, sizeof(cl_float4), &clCellb);
    err |= clSetKernelArg(surfaceAreaKernel, 8, sizeof(cl_float4), &clCellc);
    err |= clSetKernelArg(surfaceAreaKernel, 9, sizeof(cl_float4), &inverse_clCella);
    err |= clSetKernelArg(surfaceAreaKernel, 10, sizeof(cl_float4), &inverse_clCellb);
    err |= clSetKernelArg(surfaceAreaKernel, 11, sizeof(cl_float4), &inverse_clCellc);
    err |= clSetKernelArg(surfaceAreaKernel, 12, sizeof(cl_float4), &inverseWidths);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clSetKernelArg failed {} : {}\n", __FILE__, __LINE__));
    }


    for(std::size_t i = 0; i < number_of_iterations; ++i)
    {
      for (size_t j = 0; j < number_of_random_unit_vectors; j++)
      {
        double3 vec = random.randomVectorOnUnitSphere();
        random_unit_vectors[j] = {{cl_float(vec.x), cl_float(vec.y), cl_float(vec.z), 0.0f}};
      }

      // upload random Cartesian positions
      err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), random_unit_vectors_mem, CL_TRUE, 0, sizeof(float4) * random_unit_vectors.size(),
                                   random_unit_vectors.data(), 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
      }

      err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), surfaceAreaKernel, 1, nullptr, &global_work_size,
                                    &surfaceAreaWorkGroupSize, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed (err: {}) {} : {}\n", err, __FILE__, __LINE__));
      }

      // Read the buffer back to the array
      err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), outputArray, CL_TRUE, 0, output.size() * sizeof(float),
                                output.data(), 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error("MC_OpenCL_SurfaceArea: error in clEnqueueReadBuffer");
      }

      // One independent reading of the total area: keep sum and sum of squares for the mean and its error.
      double surface_area = std::accumulate(output.begin(), output.end(), 0.0);
      sum += surface_area;
      sum_of_squares += surface_area * surface_area;

    }

    clFinish(OpenCL::clCommandQueue.value());
  }


  time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  this->seconds = timing.count();
  double n = static_cast<double>(number_of_iterations);
  this->surfaceArea = n > 0.0 ? sum / n : 0.0;

  // Student-t critical values for a two-sided 95% interval (df = 1..20); beyond that the normal limit.
  constexpr std::array<double, 21> student_t_95{
      0.0,   12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
      2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086};
  this->surfaceAreaError = 0.0;
  if (number_of_iterations >= 3)
  {
    // Passes that all read the same have no scatter, and the two terms can then round to a difference below zero.
    double standard_error = std::sqrt(std::max(0.0, sum_of_squares - sum * sum / n) / (n * (n - 1.0)));
    std::size_t degrees_of_freedom = number_of_iterations - 1;
    double t = degrees_of_freedom < student_t_95.size() ? student_t_95[degrees_of_freedom] : 1.959963984540054;
    this->surfaceAreaError = t * standard_error;
  }

  std::ofstream myfile;
  myfile.open(structure.name + ".mc.sa.gpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  structure.writeHeader(myfile);
  probe.writeHeader(myfile);
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps (sample points per atom): {}\n", number_of_inner_steps);
  std::print(myfile, "# GPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "# The area is the mean over independent passes; the error beside it is the\n");
  std::print(myfile, "# half-width of the 95% confidence interval (Student-t times the standard\n");
  std::print(myfile, "# error of the mean from the sum and sum of squares of the passes).\n");
  std::print(myfile, "{} +/- {} [Å²]\n", this->surfaceArea, this->surfaceAreaError);
  std::print(myfile, "{} +/- {} [m²/cm³]\n", 1.0e4 * this->surfaceArea / structure.unitCell.volume,
             1.0e4 * this->surfaceAreaError / structure.unitCell.volume);
  std::print(myfile, "{} +/- {} [m²/g]\n", this->surfaceArea * structure.gravimetricFactor(),
             this->surfaceAreaError * structure.gravimetricFactor());
  myfile.close();
}

const char* MC_OpenCL_SurfaceArea::surfaceAreaKernelSource = R"foo(
// Whether a point lies inside an atom, other than the atom whose sphere it is on as it stands in this cell: images of
// that one count, since a sphere can be buried by its own repeats. An image only reaches the point if its fraction
// along an axis is within its radius over the width of the cell across that axis, a fraction f being at least f times
// that width long whatever the shape of the cell, and that is usually the one image rounding gives or none at all.
// The last term is allowance for rounding in these floats.
inline bool isBuried(float4 sample, int iatom, int numberOfAtoms, __global float4 *position, __global float *sigma,
                     float4 cella, float4 cellb, float4 cellc,
                     float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc, float4 inverseWidths)
{
  for(int jatom = 0; jatom < numberOfAtoms; ++jatom)
  {
    float4 dr = sample - position[jatom];

    float4 ds;
    ds.x = dot(inverse_cella, dr);
    ds.y = dot(inverse_cellb, dr);
    ds.z = dot(inverse_cellc, dr);
    ds.w = 0.0f;

    float4 whole = rint(ds);
    float4 t = ds - whole;

    float radius_sq = sigma[jatom] * sigma[jatom];
    float4 reach = sigma[jatom] * inverseWidths + 1.0e-4f;
    float4 first = ceil(-reach - t);
    float4 last = floor(reach - t);
    for(float nz = first.z; nz <= last.z; nz += 1.0f)
    {
      for(float ny = first.y; ny <= last.y; ny += 1.0f)
      {
        for(float nx = first.x; nx <= last.x; nx += 1.0f)
        {
          if(jatom == iatom && nx == whole.x && ny == whole.y && nz == whole.z) continue;

          float4 f = t + (float4)(nx, ny, nz, 0.0f);
          float4 image;
          image.x = dot(cella, f);
          image.y = dot(cellb, f);
          image.z = dot(cellc, f);
          image.w = 0.0f;
          if(dot(image, image) < radius_sq) return true;
        }
      }
    }
  }
  return false;
}

__kernel void ComputeSurfaceArea(__global float4 *position,
                                 __global float *sigma,
                                 __global float4 *randomCartesianPositions,
                                 __global float *output,
                                 const int numberOfAtoms,
                                 const int numberOfSlices,
                                 const float4 cella,
                                 const float4 cellb,
                                 const float4 cellc,
                                 const float4 inverse_cella,
                                 const float4 inverse_cellb,
                                 const float4 inverse_cellc,
                                 const float4 inverseWidths)
{
  int iatom = get_global_id(0);
  // integer counts: a float stops at 2^24 (16777216.0f + 1.0f == 16777216.0f), which froze the total while the
  // count went on, inflating the exposed fraction
  ulong counted = 0;
  ulong total = 0;

  if(iatom < numberOfAtoms)
  {
    float radius_i = sigma[iatom];
    float4 sphere_center = position[iatom];

    for(int slice = 0; slice < numberOfSlices; ++slice)
    {
      float4 unit_vector = randomCartesianPositions[slice];
      float4 sample = sphere_center + radius_i * unit_vector;

      bool overlap = isBuried(sample, iatom, numberOfAtoms, position, sigma,
                              cella, cellb, cellc, inverse_cella, inverse_cellb, inverse_cellc, inverseWidths);

      if(!overlap)
      {
        ++counted;
      }

      ++total;
    }

    output[ iatom ] = ((float)counted / (float)total) * 4.0f * M_PI_F * radius_i * radius_i;
  }
}
)foo";
