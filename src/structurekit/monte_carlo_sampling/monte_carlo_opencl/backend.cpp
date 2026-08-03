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

module mc_opencl_backend;

import std;

import double3;
import double3x3;
import unit_cell;
import sampled_structure;
import sampling_backend;
import opencl;

namespace
{
const char *samplingKernelSource = R"foo(
// The void of a crystal, measured against the atoms of one unit cell under the minimum-image convention.
// Everything below is one of two things: the distance from a point to the nearest atom surface, or the same
// along a segment. The rest is which points and which segments.

inline float3 minimumImage(float3 dr, float4 cella, float4 cellb, float4 cellc,
                           float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc)
{
  float4 wide = (float4)(dr.x, dr.y, dr.z, 0.0f);

  float4 ds;
  ds.x = dot(inverse_cella, wide);
  ds.y = dot(inverse_cellb, wide);
  ds.z = dot(inverse_cellc, wide);
  ds.w = 0.0f;

  float4 t = ds - rint(ds);

  float3 out;
  out.x = dot(cella, t);
  out.y = dot(cellb, t);
  out.z = dot(cellc, t);
  return out;
}

// The signed clearance: how much room there is for a sphere centred here, or how deep inside the nearest
// atom it is when there is none. This is min over atoms of (distance to the centre less the atom's radius),
// which is what a pore size is measured with.
inline float clearanceAt(float3 position, int numberOfAtoms, __global const float4 *atoms,
                         float4 cella, float4 cellb, float4 cellc,
                         float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc)
{
  float smallest = MAXFLOAT;
  for (int j = 0; j < numberOfAtoms; ++j)
  {
    float4 atom = atoms[j];
    float3 dr = minimumImage(position - atom.xyz, cella, cellb, cellc,
                             inverse_cella, inverse_cellb, inverse_cellc);
    smallest = min(smallest, length(dr) - atom.w);
  }
  return smallest;
}

// The largest sphere that can travel in a straight line from `origin` to `origin + displacement`, and where
// along the way it is hemmed in most closely, returned as the fraction of the way along. Exact rather than
// sampled: the clearance along a segment is a minimum over atoms of a distance that is smallest somewhere
// definite, and that place is found in closed form for each atom.
inline float segmentBottleneck(float3 origin, float3 displacement, int numberOfAtoms,
                               __global const float4 *atoms, float *where,
                               float4 cella, float4 cellb, float4 cellc,
                               float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc)
{
  float3 midpoint = origin + 0.5f * displacement;
  float lengthSquared = dot(displacement, displacement);

  float smallest = MAXFLOAT;
  float at = 0.0f;

  for (int j = 0; j < numberOfAtoms; ++j)
  {
    float4 atom = atoms[j];

    // From the start of the segment to the nearest image of the atom, by way of the midpoint so that the
    // image chosen is the one nearest the segment as a whole rather than the one nearest an end of it.
    float3 dr = minimumImage(atom.xyz - midpoint, cella, cellb, cellc,
                             inverse_cella, inverse_cellb, inverse_cellc) + 0.5f * displacement;

    float t = lengthSquared > 0.0f ? clamp(dot(dr, displacement) / lengthSquared, 0.0f, 1.0f) : 0.0f;
    float3 closest = dr - t * displacement;

    float clearance = length(closest) - atom.w;
    if (clearance < smallest)
    {
      smallest = clearance;
      at = t;
    }
  }

  *where = at;
  return smallest;
}

// A stream of directions per work item, so that a walk does not depend on how the walks were shared out.
//
// This is PCG rather than one of the shift-register generators that are usual in kernels, and the reason is
// that the walks are searches. A walk uphill and a bend sideways both keep only what improves on what they
// had, so a generator whose directions are the least bit correlated does not make the answer noisier, it
// makes it smaller: the search covers less of the sphere and finds less room. Every work item here starts
// from a seed one apart from its neighbour's, which is exactly the case the cheap generators are worst at.
inline uint nextRandom(ulong *state)
{
  ulong previous = *state;
  *state = previous * 6364136223846793005UL + 1442695040888963407UL;

  uint xorshifted = (uint)(((previous >> 18) ^ previous) >> 27);
  uint rotation = (uint)(previous >> 59);
  return (xorshifted >> rotation) | (xorshifted << ((32u - rotation) & 31u));
}

inline float uniform01(ulong *state)
{
  return (float)(nextRandom(state) >> 8) * (1.0f / 16777216.0f);
}

// Marsaglia's method: a point of the unit disc lifted onto the sphere, which is uniform there.
inline float3 randomOnUnitSphere(ulong *state)
{
  for (int attempt = 0; attempt < 64; ++attempt)
  {
    float u = 2.0f * uniform01(state) - 1.0f;
    float v = 2.0f * uniform01(state) - 1.0f;
    float s = u * u + v * v;
    if (s >= 1.0f || s <= 0.0f) continue;

    float root = sqrt(1.0f - s);
    return (float3)(2.0f * u * root, 2.0f * v * root, 1.0f - 2.0f * s);
  }
  return (float3)(0.0f, 0.0f, 1.0f);
}

inline ulong seedOf(uint seed)
{
  ulong state = (ulong)seed * 6364136223846793005UL + 1442695040888963407UL;
  for (int i = 0; i < 4; ++i) nextRandom(&state);
  return state;
}

__kernel void ComputeClearances(__global const float4 *atoms,
                                __global const float4 *positions,
                                __global float *clearances,
                                const int numberOfAtoms,
                                const int numberOfCases,
                                const float4 cella, const float4 cellb, const float4 cellc,
                                const float4 inverse_cella, const float4 inverse_cellb,
                                const float4 inverse_cellc)
{
  int id = get_global_id(0);
  if (id >= numberOfCases) return;

  clearances[id] = clearanceAt(positions[id].xyz, numberOfAtoms, atoms, cella, cellb, cellc,
                               inverse_cella, inverse_cellb, inverse_cellc);
}

__kernel void ComputeStraightWays(__global const float4 *atoms,
                                  __global const float4 *origins,
                                  __global const float4 *displacements,
                                  __global float4 *ways,       // xyz the bottleneck, w its radius
                                  const int numberOfAtoms,
                                  const int numberOfCases,
                                  const float4 cella, const float4 cellb, const float4 cellc,
                                  const float4 inverse_cella, const float4 inverse_cellb,
                                  const float4 inverse_cellc)
{
  int id = get_global_id(0);
  if (id >= numberOfCases) return;

  float3 origin = origins[id].xyz;
  float3 displacement = displacements[id].xyz;

  float at = 0.0f;
  float radius = segmentBottleneck(origin, displacement, numberOfAtoms, atoms, &at,
                                   cella, cellb, cellc, inverse_cella, inverse_cellb, inverse_cellc);

  float3 position = origin + at * displacement;
  ways[id] = (float4)(position.x, position.y, position.z, radius);
}

__kernel void WalkUphill(__global const float4 *atoms,
                         __global const float4 *starts,   // xyz the place, w the room there
                         __global const uint *seeds,
                         __global float4 *ends,
                         const int numberOfAtoms,
                         const int numberOfCases,
                         const int numberOfSteps,
                         const float shrink,
                         const float4 cella, const float4 cellb, const float4 cellc,
                         const float4 inverse_cella, const float4 inverse_cellb,
                         const float4 inverse_cellc)
{
  int id = get_global_id(0);
  if (id >= numberOfCases) return;

  ulong state = seedOf(seeds[id]);

  float4 start = starts[id];
  float3 best = start.xyz;
  float bestRadius = start.w;
  float step = fmax(bestRadius, 0.1f);

  for (int attempt = 0; attempt < numberOfSteps; ++attempt)
  {
    float3 trial = best + step * randomOnUnitSphere(&state);
    float clearance = clearanceAt(trial, numberOfAtoms, atoms, cella, cellb, cellc,
                                  inverse_cella, inverse_cellb, inverse_cellc);

    if (clearance > bestRadius)
    {
      best = trial;
      bestRadius = clearance;
    }

    step *= shrink;
  }

  ends[id] = (float4)(best.x, best.y, best.z, bestRadius);
}

// The bent way, as an explicit stack rather than by recursion, there being none to be had here. The tree is
// the same one the processor walks: at each level the straight line's tightest point is walked out sideways
// and, where that finds room, the two halves are put through the same again.
#define MAXIMUM_BENDING_DEPTH 4

typedef struct
{
  float3 origin;
  float3 displacement;
  float3 straightPosition;
  float3 straightDirection;
  float3 lifted;
  float3 firstPosition;
  float3 firstDirection;
  float straightRadius;
  float firstRadius;
  int depth;
  int stage;
} BendFrame;

__kernel void ComputeWidestWays(__global const float4 *atoms,
                                __global const float4 *origins,
                                __global const float4 *displacements,
                                __global const uint *seeds,
                                __global float4 *ways,        // xyz the bottleneck, w its radius
                                __global float4 *directions,  // xyz the way's direction there
                                const int numberOfAtoms,
                                const int numberOfCases,
                                const int depth,
                                const int refinementSteps,
                                const float bendShrink,
                                const float4 cella, const float4 cellb, const float4 cellc,
                                const float4 inverse_cella, const float4 inverse_cellb,
                                const float4 inverse_cellc)
{
  int id = get_global_id(0);
  if (id >= numberOfCases) return;

  ulong state = seedOf(seeds[id]);

  BendFrame stack[MAXIMUM_BENDING_DEPTH + 1];
  int top = 0;

  stack[0].origin = origins[id].xyz;
  stack[0].displacement = displacements[id].xyz;
  stack[0].depth = min(depth, MAXIMUM_BENDING_DEPTH);
  stack[0].stage = 0;

  float returnRadius = 0.0f;
  float3 returnPosition = (float3)(0.0f, 0.0f, 0.0f);
  float3 returnDirection = (float3)(0.0f, 0.0f, 1.0f);

  while (top >= 0)
  {
    if (stack[top].stage == 0)
    {
      float3 origin = stack[top].origin;
      float3 displacement = stack[top].displacement;

      float at = 0.0f;
      float straightRadius = segmentBottleneck(origin, displacement, numberOfAtoms, atoms, &at,
                                               cella, cellb, cellc,
                                               inverse_cella, inverse_cellb, inverse_cellc);

      float len = length(displacement);
      float3 straightDirection = len > 0.0f ? displacement / len : (float3)(0.0f, 0.0f, 1.0f);
      float3 straightPosition = origin + at * displacement;

      stack[top].straightRadius = straightRadius;
      stack[top].straightPosition = straightPosition;
      stack[top].straightDirection = straightDirection;

      if (stack[top].depth == 0 || len <= 0.0f)
      {
        returnRadius = straightRadius;
        returnPosition = straightPosition;
        returnDirection = straightDirection;
        --top;
        continue;
      }

      // Walk the tightest point out across the line. Along it the room grows towards either end and says
      // nothing about the passage, so only the part perpendicular to it is kept.
      float3 lifted = straightPosition;
      float liftedRadius = straightRadius;
      float step = fmax(fabs(straightRadius), 0.1f);

      for (int attempt = 0; attempt < refinementSteps; ++attempt)
      {
        float3 across = randomOnUnitSphere(&state);
        across -= dot(across, straightDirection) * straightDirection;

        float norm = length(across);
        if (norm > 0.0f)
        {
          float3 trial = lifted + (step / norm) * across;
          float clearance = clearanceAt(trial, numberOfAtoms, atoms, cella, cellb, cellc,
                                        inverse_cella, inverse_cellb, inverse_cellc);

          if (clearance > liftedRadius)
          {
            lifted = trial;
            liftedRadius = clearance;
          }
        }

        step *= bendShrink;
      }

      if (liftedRadius <= straightRadius)
      {
        returnRadius = straightRadius;
        returnPosition = straightPosition;
        returnDirection = straightDirection;
        --top;
        continue;
      }

      stack[top].lifted = lifted;
      stack[top].stage = 1;

      int child = top + 1;
      stack[child].origin = origin;
      stack[child].displacement = lifted - origin;
      stack[child].depth = stack[top].depth - 1;
      stack[child].stage = 0;
      top = child;
      continue;
    }

    if (stack[top].stage == 1)
    {
      stack[top].firstRadius = returnRadius;
      stack[top].firstPosition = returnPosition;
      stack[top].firstDirection = returnDirection;
      stack[top].stage = 2;

      int child = top + 1;
      stack[child].origin = stack[top].lifted;
      stack[child].displacement = stack[top].origin + stack[top].displacement - stack[top].lifted;
      stack[child].depth = stack[top].depth - 1;
      stack[child].stage = 0;
      top = child;
      continue;
    }

    // The narrowest point of the bent way is the narrower of the two halves' own, and the bent way is taken
    // only where it beats the straight one.
    float bentRadius = returnRadius;
    float3 bentPosition = returnPosition;
    float3 bentDirection = returnDirection;

    if (stack[top].firstRadius <= bentRadius)
    {
      bentRadius = stack[top].firstRadius;
      bentPosition = stack[top].firstPosition;
      bentDirection = stack[top].firstDirection;
    }

    if (bentRadius > stack[top].straightRadius)
    {
      returnRadius = bentRadius;
      returnPosition = bentPosition;
      returnDirection = bentDirection;
    }
    else
    {
      returnRadius = stack[top].straightRadius;
      returnPosition = stack[top].straightPosition;
      returnDirection = stack[top].straightDirection;
    }
    --top;
  }

  ways[id] = (float4)(returnPosition.x, returnPosition.y, returnPosition.z, returnRadius);
  directions[id] = (float4)(returnDirection.x, returnDirection.y, returnDirection.z, 0.0f);
}
)foo";

// The compiled program and its four kernels, built once for the life of the process. Compiling costs more
// than any single batch does, and a roadmap asks for hundreds of them.
struct SamplingDevice
{
  cl_program program{nullptr};
  cl_kernel clearances{nullptr};
  cl_kernel straightWays{nullptr};
  cl_kernel walksUphill{nullptr};
  cl_kernel widestWays{nullptr};
  std::size_t workGroupSize{64};

  SamplingDevice()
  {
    cl_int err;

    program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &samplingKernelSource, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("mc_opencl_backend: clCreateProgramWithSource failed\n");
    }

    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t length = 0;
      std::string log(16384, '\0');
      clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, log.size(), log.data(),
                            &length);
      log.resize(length);
      throw std::runtime_error(std::format("mc_opencl_backend: failed to build the sampling kernels: {}\n", log));
    }

    clearances = create("ComputeClearances");
    straightWays = create("ComputeStraightWays");
    walksUphill = create("WalkUphill");
    widestWays = create("ComputeWidestWays");

    // The bent way carries a stack in private memory and so wants the smallest of the four work groups; one
    // size for all of them keeps the staging simple and costs the other three nothing measurable.
    workGroupSize = std::min({sizeOf(clearances), sizeOf(straightWays), sizeOf(walksUphill),
                              sizeOf(widestWays)});
    workGroupSize = std::max<std::size_t>(workGroupSize, 1);
  }

  ~SamplingDevice()
  {
    if (clearances != nullptr) clReleaseKernel(clearances);
    if (straightWays != nullptr) clReleaseKernel(straightWays);
    if (walksUphill != nullptr) clReleaseKernel(walksUphill);
    if (widestWays != nullptr) clReleaseKernel(widestWays);
    if (program != nullptr) clReleaseProgram(program);
  }

  SamplingDevice(const SamplingDevice &) = delete;
  SamplingDevice &operator=(const SamplingDevice &) = delete;

 private:
  cl_kernel create(const char *name)
  {
    cl_int err;
    cl_kernel kernel = clCreateKernel(program, name, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("mc_opencl_backend: clCreateKernel failed for '{}'\n", name));
    }
    return kernel;
  }

  static std::size_t sizeOf(cl_kernel kernel)
  {
    std::size_t size = 1;
    clGetKernelWorkGroupInfo(kernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(std::size_t),
                             &size, nullptr);
    return size;
  }
};

std::shared_ptr<SamplingDevice> device()
{
  static std::shared_ptr<SamplingDevice> shared = std::make_shared<SamplingDevice>();
  return shared;
}

// A device buffer that releases itself, so that a kernel that throws does not leak one.
class Buffer
{
 public:
  Buffer(cl_mem_flags flags, std::size_t bytes, const void *from = nullptr)
  {
    cl_int err;
    this->memory = clCreateBuffer(OpenCL::clContext.value(), flags, std::max<std::size_t>(bytes, 1), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("mc_opencl_backend: clCreateBuffer failed\n");
    }

    if (from != nullptr && bytes > 0)
    {
      err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), this->memory, CL_TRUE, 0, bytes, from, 0, nullptr,
                                 nullptr);
      if (err != CL_SUCCESS)
      {
        clReleaseMemObject(this->memory);
        throw std::runtime_error("mc_opencl_backend: clEnqueueWriteBuffer failed\n");
      }
    }
  }

  ~Buffer()
  {
    if (this->memory != nullptr) clReleaseMemObject(this->memory);
  }

  Buffer(const Buffer &) = delete;
  Buffer &operator=(const Buffer &) = delete;

  void read(std::size_t bytes, void *into) const
  {
    if (bytes == 0) return;

    cl_int err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), this->memory, CL_TRUE, 0, bytes, into, 0,
                                     nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("mc_opencl_backend: clEnqueueReadBuffer failed\n");
    }
  }

  const cl_mem *operator&() const { return &this->memory; }

 private:
  cl_mem memory{nullptr};
};

// The cell, as the kernels want it: the rows of the matrix and of its inverse, so that a minimum image is
// three dot products each way.
struct DeviceCell
{
  cl_float4 a, b, c;
  cl_float4 inverseA, inverseB, inverseC;

  explicit DeviceCell(const UnitCell &unitCell)
  {
    const double3x3 &cell = unitCell.cell;
    const double3x3 &inverse = unitCell.inverseCell;

    auto row = [](const double3x3 &m, std::size_t i) {
      return cl_float4{{static_cast<cl_float>(m[0][i]), static_cast<cl_float>(m[1][i]),
                        static_cast<cl_float>(m[2][i]), 0.0f}};
    };

    a = row(cell, 0);
    b = row(cell, 1);
    c = row(cell, 2);
    inverseA = row(inverse, 0);
    inverseB = row(inverse, 1);
    inverseC = row(inverse, 2);
  }
};

// The atoms, as a position and a contact radius packed together, which is one fetch per atom in the inner
// loop rather than two.
std::vector<cl_float4> deviceAtoms(const SampledStructure &structure)
{
  std::vector<cl_float4> atoms(structure.size());
  for (std::size_t i = 0; i < structure.size(); ++i)
  {
    atoms[i] = {{static_cast<cl_float>(structure.positions[i].x), static_cast<cl_float>(structure.positions[i].y),
                 static_cast<cl_float>(structure.positions[i].z), static_cast<cl_float>(structure.radii[i])}};
  }
  return atoms;
}

std::vector<cl_float4> devicePoints(std::span<const double3> points)
{
  std::vector<cl_float4> staged(points.size());
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    staged[i] = {{static_cast<cl_float>(points[i].x), static_cast<cl_float>(points[i].y),
                  static_cast<cl_float>(points[i].z), 0.0f}};
  }
  return staged;
}

void enqueue(cl_kernel kernel, std::size_t cases, std::size_t workGroupSize)
{
  if (cases == 0) return;

  std::size_t global = ((cases + workGroupSize - 1) / workGroupSize) * workGroupSize;

  cl_int err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &global, &workGroupSize,
                                      0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("mc_opencl_backend: clEnqueueNDRangeKernel failed (error {})\n", err));
  }
}
}  // namespace

bool samplingOpenCLAvailable()
{
  return OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value() && OpenCL::clCommandQueue.has_value();
}

SamplingBackend samplingBackendOpenCL()
{
  if (!samplingOpenCLAvailable())
  {
    throw std::runtime_error("The sampled routes were asked for on the GPU, and there is no device to run on\n");
  }

  SamplingBackend backend;
  backend.name = "gpu";

  backend.clearances = [](const SampledStructure &structure, std::span<const double3> positions,
                          std::span<double> into)
  {
    if (positions.empty()) return;

    std::shared_ptr<SamplingDevice> gpu = device();
    DeviceCell cell(structure.unitCell);

    std::vector<cl_float4> atoms = deviceAtoms(structure);
    std::vector<cl_float4> staged = devicePoints(positions);

    Buffer atomsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * atoms.size(), atoms.data());
    Buffer pointsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * staged.size(), staged.data());
    Buffer output(CL_MEM_WRITE_ONLY, sizeof(cl_float) * staged.size());

    cl_int numberOfAtoms = static_cast<cl_int>(atoms.size());
    cl_int numberOfCases = static_cast<cl_int>(staged.size());

    cl_int err = clSetKernelArg(gpu->clearances, 0, sizeof(cl_mem), &atomsBuffer);
    err |= clSetKernelArg(gpu->clearances, 1, sizeof(cl_mem), &pointsBuffer);
    err |= clSetKernelArg(gpu->clearances, 2, sizeof(cl_mem), &output);
    err |= clSetKernelArg(gpu->clearances, 3, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(gpu->clearances, 4, sizeof(cl_int), &numberOfCases);
    err |= clSetKernelArg(gpu->clearances, 5, sizeof(cl_float4), &cell.a);
    err |= clSetKernelArg(gpu->clearances, 6, sizeof(cl_float4), &cell.b);
    err |= clSetKernelArg(gpu->clearances, 7, sizeof(cl_float4), &cell.c);
    err |= clSetKernelArg(gpu->clearances, 8, sizeof(cl_float4), &cell.inverseA);
    err |= clSetKernelArg(gpu->clearances, 9, sizeof(cl_float4), &cell.inverseB);
    err |= clSetKernelArg(gpu->clearances, 10, sizeof(cl_float4), &cell.inverseC);
    if (err != CL_SUCCESS) throw std::runtime_error("mc_opencl_backend: clSetKernelArg failed (clearances)\n");

    enqueue(gpu->clearances, staged.size(), gpu->workGroupSize);

    std::vector<cl_float> answers(staged.size());
    output.read(sizeof(cl_float) * answers.size(), answers.data());

    for (std::size_t i = 0; i < answers.size(); ++i) into[i] = static_cast<double>(answers[i]);
  };

  backend.straightWays = [](const SampledStructure &structure, std::span<const double3> origins,
                            std::span<const double3> displacements, std::span<SampledWay> into)
  {
    if (origins.empty()) return;

    std::shared_ptr<SamplingDevice> gpu = device();
    DeviceCell cell(structure.unitCell);

    std::vector<cl_float4> atoms = deviceAtoms(structure);
    std::vector<cl_float4> stagedOrigins = devicePoints(origins);
    std::vector<cl_float4> stagedDisplacements = devicePoints(displacements);

    Buffer atomsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * atoms.size(), atoms.data());
    Buffer originsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * stagedOrigins.size(), stagedOrigins.data());
    Buffer displacementsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * stagedDisplacements.size(),
                               stagedDisplacements.data());
    Buffer output(CL_MEM_WRITE_ONLY, sizeof(cl_float4) * stagedOrigins.size());

    cl_int numberOfAtoms = static_cast<cl_int>(atoms.size());
    cl_int numberOfCases = static_cast<cl_int>(stagedOrigins.size());

    cl_int err = clSetKernelArg(gpu->straightWays, 0, sizeof(cl_mem), &atomsBuffer);
    err |= clSetKernelArg(gpu->straightWays, 1, sizeof(cl_mem), &originsBuffer);
    err |= clSetKernelArg(gpu->straightWays, 2, sizeof(cl_mem), &displacementsBuffer);
    err |= clSetKernelArg(gpu->straightWays, 3, sizeof(cl_mem), &output);
    err |= clSetKernelArg(gpu->straightWays, 4, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(gpu->straightWays, 5, sizeof(cl_int), &numberOfCases);
    err |= clSetKernelArg(gpu->straightWays, 6, sizeof(cl_float4), &cell.a);
    err |= clSetKernelArg(gpu->straightWays, 7, sizeof(cl_float4), &cell.b);
    err |= clSetKernelArg(gpu->straightWays, 8, sizeof(cl_float4), &cell.c);
    err |= clSetKernelArg(gpu->straightWays, 9, sizeof(cl_float4), &cell.inverseA);
    err |= clSetKernelArg(gpu->straightWays, 10, sizeof(cl_float4), &cell.inverseB);
    err |= clSetKernelArg(gpu->straightWays, 11, sizeof(cl_float4), &cell.inverseC);
    if (err != CL_SUCCESS) throw std::runtime_error("mc_opencl_backend: clSetKernelArg failed (straight ways)\n");

    enqueue(gpu->straightWays, stagedOrigins.size(), gpu->workGroupSize);

    std::vector<cl_float4> answers(stagedOrigins.size());
    output.read(sizeof(cl_float4) * answers.size(), answers.data());

    for (std::size_t i = 0; i < answers.size(); ++i)
    {
      double length = displacements[i].length();

      into[i] = SampledWay{
          .radius = static_cast<double>(answers[i].s[3]),
          .position = double3(static_cast<double>(answers[i].s[0]), static_cast<double>(answers[i].s[1]),
                              static_cast<double>(answers[i].s[2])),
          .direction = length > 0.0 ? (1.0 / length) * displacements[i] : double3(0.0, 0.0, 1.0)};
    }
  };

  backend.walksUphill = [](const SampledStructure &structure, std::span<const SampledPeak> starts,
                           std::span<const std::uint32_t> seeds, std::size_t steps, std::span<SampledPeak> into)
  {
    if (starts.empty()) return;

    std::shared_ptr<SamplingDevice> gpu = device();
    DeviceCell cell(structure.unitCell);

    std::vector<cl_float4> atoms = deviceAtoms(structure);

    std::vector<cl_float4> stagedStarts(starts.size());
    for (std::size_t i = 0; i < starts.size(); ++i)
    {
      stagedStarts[i] = {{static_cast<cl_float>(starts[i].position.x), static_cast<cl_float>(starts[i].position.y),
                          static_cast<cl_float>(starts[i].position.z), static_cast<cl_float>(starts[i].radius)}};
    }

    std::vector<cl_uint> stagedSeeds(seeds.begin(), seeds.end());

    Buffer atomsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * atoms.size(), atoms.data());
    Buffer startsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * stagedStarts.size(), stagedStarts.data());
    Buffer seedsBuffer(CL_MEM_READ_ONLY, sizeof(cl_uint) * stagedSeeds.size(), stagedSeeds.data());
    Buffer output(CL_MEM_WRITE_ONLY, sizeof(cl_float4) * stagedStarts.size());

    cl_int numberOfAtoms = static_cast<cl_int>(atoms.size());
    cl_int numberOfCases = static_cast<cl_int>(stagedStarts.size());
    cl_int numberOfSteps = static_cast<cl_int>(steps);

    const double finalStepFraction = 1.0e-4;
    cl_float shrink = static_cast<cl_float>(
        steps > 0 ? std::pow(finalStepFraction, 1.0 / static_cast<double>(steps)) : 1.0);

    cl_int err = clSetKernelArg(gpu->walksUphill, 0, sizeof(cl_mem), &atomsBuffer);
    err |= clSetKernelArg(gpu->walksUphill, 1, sizeof(cl_mem), &startsBuffer);
    err |= clSetKernelArg(gpu->walksUphill, 2, sizeof(cl_mem), &seedsBuffer);
    err |= clSetKernelArg(gpu->walksUphill, 3, sizeof(cl_mem), &output);
    err |= clSetKernelArg(gpu->walksUphill, 4, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(gpu->walksUphill, 5, sizeof(cl_int), &numberOfCases);
    err |= clSetKernelArg(gpu->walksUphill, 6, sizeof(cl_int), &numberOfSteps);
    err |= clSetKernelArg(gpu->walksUphill, 7, sizeof(cl_float), &shrink);
    err |= clSetKernelArg(gpu->walksUphill, 8, sizeof(cl_float4), &cell.a);
    err |= clSetKernelArg(gpu->walksUphill, 9, sizeof(cl_float4), &cell.b);
    err |= clSetKernelArg(gpu->walksUphill, 10, sizeof(cl_float4), &cell.c);
    err |= clSetKernelArg(gpu->walksUphill, 11, sizeof(cl_float4), &cell.inverseA);
    err |= clSetKernelArg(gpu->walksUphill, 12, sizeof(cl_float4), &cell.inverseB);
    err |= clSetKernelArg(gpu->walksUphill, 13, sizeof(cl_float4), &cell.inverseC);
    if (err != CL_SUCCESS) throw std::runtime_error("mc_opencl_backend: clSetKernelArg failed (walks)\n");

    enqueue(gpu->walksUphill, stagedStarts.size(), gpu->workGroupSize);

    std::vector<cl_float4> answers(stagedStarts.size());
    output.read(sizeof(cl_float4) * answers.size(), answers.data());

    for (std::size_t i = 0; i < answers.size(); ++i)
    {
      into[i] = SampledPeak{
          .radius = static_cast<double>(answers[i].s[3]),
          .position = double3(static_cast<double>(answers[i].s[0]), static_cast<double>(answers[i].s[1]),
                              static_cast<double>(answers[i].s[2]))};
    }
  };

  backend.widestWays = [](const SampledStructure &structure, std::span<const double3> origins,
                          std::span<const double3> displacements, std::span<const std::uint32_t> seeds,
                          std::size_t depth, std::span<SampledWay> into)
  {
    if (origins.empty()) return;

    std::shared_ptr<SamplingDevice> gpu = device();
    DeviceCell cell(structure.unitCell);

    std::vector<cl_float4> atoms = deviceAtoms(structure);
    std::vector<cl_float4> stagedOrigins = devicePoints(origins);
    std::vector<cl_float4> stagedDisplacements = devicePoints(displacements);
    std::vector<cl_uint> stagedSeeds(seeds.begin(), seeds.end());

    Buffer atomsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * atoms.size(), atoms.data());
    Buffer originsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * stagedOrigins.size(), stagedOrigins.data());
    Buffer displacementsBuffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * stagedDisplacements.size(),
                               stagedDisplacements.data());
    Buffer seedsBuffer(CL_MEM_READ_ONLY, sizeof(cl_uint) * stagedSeeds.size(), stagedSeeds.data());
    Buffer wayOutput(CL_MEM_WRITE_ONLY, sizeof(cl_float4) * stagedOrigins.size());
    Buffer directionOutput(CL_MEM_WRITE_ONLY, sizeof(cl_float4) * stagedOrigins.size());

    cl_int numberOfAtoms = static_cast<cl_int>(atoms.size());
    cl_int numberOfCases = static_cast<cl_int>(stagedOrigins.size());
    cl_int bendingDepth = static_cast<cl_int>(depth);
    cl_int refinementSteps = static_cast<cl_int>(samplingRefinementSteps);
    cl_float bendShrink =
        static_cast<cl_float>(std::pow(1.0e-3, 1.0 / static_cast<double>(samplingRefinementSteps)));

    cl_int err = clSetKernelArg(gpu->widestWays, 0, sizeof(cl_mem), &atomsBuffer);
    err |= clSetKernelArg(gpu->widestWays, 1, sizeof(cl_mem), &originsBuffer);
    err |= clSetKernelArg(gpu->widestWays, 2, sizeof(cl_mem), &displacementsBuffer);
    err |= clSetKernelArg(gpu->widestWays, 3, sizeof(cl_mem), &seedsBuffer);
    err |= clSetKernelArg(gpu->widestWays, 4, sizeof(cl_mem), &wayOutput);
    err |= clSetKernelArg(gpu->widestWays, 5, sizeof(cl_mem), &directionOutput);
    err |= clSetKernelArg(gpu->widestWays, 6, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(gpu->widestWays, 7, sizeof(cl_int), &numberOfCases);
    err |= clSetKernelArg(gpu->widestWays, 8, sizeof(cl_int), &bendingDepth);
    err |= clSetKernelArg(gpu->widestWays, 9, sizeof(cl_int), &refinementSteps);
    err |= clSetKernelArg(gpu->widestWays, 10, sizeof(cl_float), &bendShrink);
    err |= clSetKernelArg(gpu->widestWays, 11, sizeof(cl_float4), &cell.a);
    err |= clSetKernelArg(gpu->widestWays, 12, sizeof(cl_float4), &cell.b);
    err |= clSetKernelArg(gpu->widestWays, 13, sizeof(cl_float4), &cell.c);
    err |= clSetKernelArg(gpu->widestWays, 14, sizeof(cl_float4), &cell.inverseA);
    err |= clSetKernelArg(gpu->widestWays, 15, sizeof(cl_float4), &cell.inverseB);
    err |= clSetKernelArg(gpu->widestWays, 16, sizeof(cl_float4), &cell.inverseC);
    if (err != CL_SUCCESS) throw std::runtime_error("mc_opencl_backend: clSetKernelArg failed (widest ways)\n");

    enqueue(gpu->widestWays, stagedOrigins.size(), gpu->workGroupSize);

    std::vector<cl_float4> ways(stagedOrigins.size());
    std::vector<cl_float4> directions(stagedOrigins.size());
    wayOutput.read(sizeof(cl_float4) * ways.size(), ways.data());
    directionOutput.read(sizeof(cl_float4) * directions.size(), directions.data());

    for (std::size_t i = 0; i < ways.size(); ++i)
    {
      into[i] = SampledWay{
          .radius = static_cast<double>(ways[i].s[3]),
          .position = double3(static_cast<double>(ways[i].s[0]), static_cast<double>(ways[i].s[1]),
                              static_cast<double>(ways[i].s[2])),
          .direction = double3(static_cast<double>(directions[i].s[0]), static_cast<double>(directions[i].s[1]),
                               static_cast<double>(directions[i].s[2]))};
    }
  };

  return backend;
}
