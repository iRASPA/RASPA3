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

void MC_OpenCL_SurfaceArea::run(const SampledStructure &structure, const SampledProbe &probe,
                                std::optional<std::size_t> numberOfIterations,
                                std::optional<std::size_t> numberOfInnerSteps)
{
  RandomNumber random{std::nullopt};
  cl_int err;
  std::chrono::steady_clock::time_point time_begin, time_end;

  std::size_t number_of_iterations = numberOfIterations.value_or(100);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);

  double3x3 unit_cell = structure.unitCell.cell;
  double3x3 inverse_unit_cell = structure.unitCell.inverseCell;

  std::size_t numberOfAtoms = structure.size();
  size_t global_work_size = (numberOfAtoms + surfaceAreaWorkGroupSize - 1) & ~(surfaceAreaWorkGroupSize - 1);

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

  // upload position array
  cl_mem inputPos =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * pos.size(), nullptr, &err);
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
  cl_mem inputSigma =
      clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
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
  cl_mem outputArray =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * output.size(), nullptr, &err);
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

  // Same image shell as SampledStructure::overlaps: MIC alone misses burials when 2r exceeds a cell edge.
  double max_radius = structure.radii.empty() ? 0.0 : *std::ranges::max_element(structure.radii);
  const double3x3 &cell = structure.unitCell.cell;
  double shortest_edge = std::min({double3(cell[0][0], cell[0][1], cell[0][2]).length(),
                                   double3(cell[1][0], cell[1][1], cell[1][2]).length(),
                                   double3(cell[2][0], cell[2][1], cell[2][2]).length()});
  cl_int use_minimum_image = (2.0 * max_radius <= shortest_edge) ? 1 : 0;

  double3 a(cell[0][0], cell[0][1], cell[0][2]);
  double3 b(cell[1][0], cell[1][1], cell[1][2]);
  double3 c(cell[2][0], cell[2][1], cell[2][2]);
  double spread = 0.5 * (a.length() + b.length() + c.length());
  double reach = 2.0 * spread + max_radius;
  double3 widths = structure.unitCell.perpendicularWidths();
  auto along = [&](double width)
  { return static_cast<cl_int>(std::clamp(std::ceil(reach / std::max(width, 1.0e-9)), 1.0, 8.0)); };
  cl_int shell_x = along(widths.x);
  cl_int shell_y = along(widths.y);
  cl_int shell_z = along(widths.z);

  std::size_t number_of_random_unit_vectors{number_of_inner_steps};
  std::vector<cl_float4> random_unit_vectors(number_of_random_unit_vectors);
  cl_mem random_unit_vectors_mem =
          clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) *  random_unit_vectors.size(), nullptr, &err);
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
  err |= clSetKernelArg(surfaceAreaKernel, 12, sizeof(cl_int), &use_minimum_image);
  err |= clSetKernelArg(surfaceAreaKernel, 13, sizeof(cl_int), &shell_x);
  err |= clSetKernelArg(surfaceAreaKernel, 14, sizeof(cl_int), &shell_y);
  err |= clSetKernelArg(surfaceAreaKernel, 15, sizeof(cl_int), &shell_z);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clSetKernelArg failed {} : {}\n", __FILE__, __LINE__));
  }


  double sum{};
  double sum_of_squares{};
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

  clReleaseMemObject(random_unit_vectors_mem);


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
    double standard_error = std::sqrt((sum_of_squares - sum * sum / n) / (n * (n - 1.0)));
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
                                 const int use_minimum_image,
                                 const int shell_x,
                                 const int shell_y,
                                 const int shell_z)
{
  int iatom = get_global_id(0);
  float counted = 0.0f;
  float total = 0.0f;

  if(iatom < numberOfAtoms)
  {
    float radius_i = sigma[iatom];
    float4 sphere_center = position[iatom];

    for(int slice = 0; slice < numberOfSlices; ++slice)
    {
      float4 unit_vector = randomCartesianPositions[slice];
      float4 sample = sphere_center + radius_i * unit_vector;

      bool overlap = false;
      if(use_minimum_image)
      {
        // Large cells: one minimum-image test per other atom is enough.
        for(int jatom = 0; jatom < numberOfAtoms; ++jatom)
        {
          if(jatom == iatom) continue;

          float4 dr = sample - position[jatom];
          float4 ds;
          ds.x = dot(inverse_cella, dr);
          ds.y = dot(inverse_cellb, dr);
          ds.z = dot(inverse_cellc, dr);
          ds.w = 0.0f;
          float4 t = ds - rint(ds);
          dr.x = dot(cella, t);
          dr.y = dot(cellb, t);
          dr.z = dot(cellc, t);
          dr.w = 0.0f;

          if(dot(dr, dr) < sigma[jatom] * sigma[jatom])
          {
            overlap = true;
            break;
          }
        }
      }
      else
      {
        // Small cells: search enough lattice images that a sphere with 2r > L is still caught,
        // including periodic copies of the atom whose sphere was sampled.
        for(int jatom = 0; jatom < numberOfAtoms && !overlap; ++jatom)
        {
          float radius_sq = sigma[jatom] * sigma[jatom];
          for(int nc = -shell_z; nc <= shell_z && !overlap; ++nc)
          {
            for(int nb = -shell_y; nb <= shell_y && !overlap; ++nb)
            {
              for(int na = -shell_x; na <= shell_x; ++na)
              {
                if(jatom == iatom && na == 0 && nb == 0 && nc == 0) continue;

                float4 lattice = (float4)((float)na, (float)nb, (float)nc, 0.0f);
                float4 translation;
                translation.x = dot(cella, lattice);
                translation.y = dot(cellb, lattice);
                translation.z = dot(cellc, lattice);
                translation.w = 0.0f;

                float4 dr = sample - (position[jatom] + translation);
                if(dot(dr, dr) < radius_sq)
                {
                  overlap = true;
                  break;
                }
              }
            }
          }
        }
      }

      if(!overlap)
      {
        counted += 1.0f;
      }

      total += 1.0f;
    }

    output[ iatom ] = (counted / total) * 4.0f * M_PI_F * radius_i * radius_i;
  }
}
)foo";
