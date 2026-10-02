module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module spatial_decomposition_opencl_mesh;

import std;

import double3;
import double3x3;
import int3;
import simulationbox;
import opencl;
import spatial_decomposition_pppm;
import spatial_decomposition_opencl_handles;

using OpenCLDevice::check;
using OpenCLDevice::roundUp;

namespace
{
constexpr std::size_t influenceGroup = 128;
constexpr std::size_t influencePartials = 14;
constexpr std::size_t influencePointsPerItem = 8;  // INFLUENCE_POINTS of the kernel source
constexpr std::size_t pointGroup = 64;
constexpr std::size_t maximumTile = 16;
constexpr double fixedPointScale = 16777216.0;  // 2^24: 6e-8 e resolution, +-128 e range per mesh point

std::size_t divideRoundUp(std::size_t a, std::size_t b) { return (a + b - 1) / b; }
}  // namespace

void OpenCLMesh::initialize(cl_context clContext, cl_device_id clDevice, std::size_t orderValue)
{
  context = clContext;
  device = clDevice;
  order = std::clamp<std::size_t>(orderValue, 3, 7);
  const std::string options = std::format("-cl-mad-enable -cl-no-signed-zeros -D MESH_ORDER={}", order);
  program = OpenCLDevice::buildProgram(context, device, openclMeshKernelSource, options.c_str(), "OpenCL mesh");
  spreadKernel = OpenCLDevice::createKernel(program.get(), "spreadCharges");
  realForwardKernel = OpenCLDevice::createKernel(program.get(), "fftRealForward");
  realBackwardKernel = OpenCLDevice::createKernel(program.get(), "fftRealBackward");
  fftKernel = OpenCLDevice::createKernel(program.get(), "fftLines");
  influenceKernel = OpenCLDevice::createKernel(program.get(), "applyInfluence");
  interpolateKernel = OpenCLDevice::createKernel(program.get(), "interpolateForces");

  cl_ulong local = 0;
  if (clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local), &local, nullptr) == CL_SUCCESS && local > 0)
  {
    localMemory = static_cast<std::size_t>(local);
  }
  std::size_t groupLimit = 0;
  if (clGetKernelWorkGroupInfo(fftKernel.get(), device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(groupLimit), &groupLimit,
                               nullptr) == CL_SUCCESS &&
      groupLimit > 0)
  {
    maxGroupSize = std::min<std::size_t>(256, groupLimit);
  }
  parameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(Parameters)));
}

void OpenCLMesh::planAxis(AxisPlan& plan, std::uint32_t N, std::uint32_t localLength, std::uint32_t axisStride,
                          std::uint32_t lineStride, std::uint32_t innerCount, std::uint32_t outerStride,
                          std::uint32_t outerCount)
{
  plan.N = N;
  plan.localLength = localLength;
  plan.axisStride = axisStride;
  plan.lineStride = lineStride;
  plan.innerCount = innerCount;
  plan.outerStride = outerStride;

  // radices: 4 while two factors 2 remain, then 2, 3, 5 (the mesh sizes have no other prime factors)
  std::uint32_t remainder = N;
  std::vector<std::uint32_t> radices;
  while (remainder % 4 == 0)
  {
    radices.push_back(4);
    remainder /= 4;
  }
  if (remainder % 2 == 0)
  {
    radices.push_back(2);
    remainder /= 2;
  }
  while (remainder % 3 == 0)
  {
    radices.push_back(3);
    remainder /= 3;
  }
  while (remainder % 5 == 0)
  {
    radices.push_back(5);
    remainder /= 5;
  }
  if (remainder != 1 || radices.size() > 16)
  {
    throw std::runtime_error(std::format("[OpenCL mesh]: mesh size {} is not a product of 2, 3 and 5\n", N));
  }
  plan.radixCode = 0;
  for (std::size_t s = 0; s < radices.size(); ++s)
  {
    const std::uint32_t code = radices[s] == 2 ? 0u : radices[s] == 3 ? 1u : radices[s] == 4 ? 2u : 3u;
    plan.radixCode |= code << (2 * s);
  }
  plan.stages = static_cast<std::uint32_t>(radices.size());

  // tile of adjacent lines per work-group: two local buffers of tile * localLength complex values
  const std::size_t budget = std::min<std::size_t>(localMemory, 65536) - 1024;
  const std::size_t perLine = 2 * static_cast<std::size_t>(localLength) * sizeof(float) * 2;
  if (perLine > budget)
  {
    throw std::runtime_error(
        std::format("[OpenCL mesh]: a mesh axis of {} points does not fit the local memory of the device ({} bytes)\n",
                    N, localMemory));
  }
  plan.tile = static_cast<std::uint32_t>(std::clamp<std::size_t>(budget / perLine, 1, maximumTile));
  plan.tile = std::min(plan.tile, innerCount);
  plan.tilesPerOuter = static_cast<std::uint32_t>(divideRoundUp(innerCount, plan.tile));
  plan.groups = static_cast<std::size_t>(plan.tilesPerOuter) * outerCount;
  plan.groupSize = std::clamp<std::size_t>(roundUp(static_cast<std::size_t>(N) * plan.tile / 2, 32), 32, maxGroupSize);

  std::vector<float> twiddle(2 * static_cast<std::size_t>(N));
  for (std::size_t i = 0; i < N; ++i)
  {
    const double angle = -2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(N);
    twiddle[2 * i] = static_cast<float>(std::cos(angle));
    twiddle[2 * i + 1] = static_cast<float>(std::sin(angle));
  }
  plan.twiddle.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, twiddle.size() * sizeof(float)));
  OpenCL::writeBuffer(plan.twiddle.get(), twiddle.size() * sizeof(float), twiddle.data());
}

void OpenCLMesh::setup(int3 meshSize, double alphaValue, double conversionFactorValue)
{
  mesh = meshSize;
  // the real-to-complex transform along z packs pairs of samples: Kz must be even
  while (mesh.z % 2 != 0)
    mesh.z = static_cast<std::int32_t>(PPPM::nextFFTFriendly(static_cast<std::size_t>(mesh.z) + 1));
  alpha = alphaValue;
  conversionFactor = conversionFactorValue;
  const std::size_t points = meshPoints();
  const std::size_t spectrum = spectrumPoints();

  meshBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, points * sizeof(std::int32_t)));
  dataBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, spectrum * 2 * sizeof(float)));
  potentialBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, points * sizeof(float)));
  influenceGroups = divideRoundUp(spectrum, influenceGroup * influencePointsPerItem);
  partialBuffer.reset(OpenCL::createBuffer(CL_MEM_WRITE_ONLY, influenceGroups * influencePartials * sizeof(float)));
  hostPartials.assign(influenceGroups * influencePartials, 0.0f);

  // the fixed-point mesh starts zeroed (the forward transform zeroes it again after every step)
  {
    const std::vector<std::int32_t> zeros(points, 0);
    OpenCL::writeBuffer(meshBuffer.get(), points * sizeof(std::int32_t), zeros.data());
  }

  auto uploadModuli = [&](OpenCLDevice::MemHandle& buffer, std::int32_t K)
  {
    const std::vector<double> moduli = PPPM::bsplineModuli(order, K);
    std::vector<float> values(moduli.size());
    for (std::size_t m = 0; m < moduli.size(); ++m) values[m] = static_cast<float>(moduli[m]);
    buffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, values.size() * sizeof(float)));
    OpenCL::writeBuffer(buffer.get(), values.size() * sizeof(float), values.data());
  };
  uploadModuli(moduliX, mesh.x);
  uploadModuli(moduliY, mesh.y);
  uploadModuli(moduliZ, mesh.z);

  const std::uint32_t Kx = static_cast<std::uint32_t>(mesh.x);
  const std::uint32_t Ky = static_cast<std::uint32_t>(mesh.y);
  const std::uint32_t Kz = static_cast<std::uint32_t>(mesh.z);
  // x-major mesh (ix * Ky + iy) * Kz + iz; half spectrum (ix * Ky + iy) * Hz + kz with Hz = Kz / 2 + 1
  const std::uint32_t M = Kz / 2;
  const std::uint32_t Hz = M + 1;
  planAxis(planZ, M, Hz, 1, Hz, Kx * Ky, 0, 1);
  planAxis(planY, Ky, Ky, Hz, 1, Hz, Ky * Hz, Kx);
  planAxis(planX, Kx, Kx, Ky * Hz, 1, Ky * Hz, 0, 1);
  {
    std::vector<float> twiddle(2 * static_cast<std::size_t>(Hz));
    for (std::size_t k = 0; k < Hz; ++k)
    {
      const double angle = -2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(Kz);
      twiddle[2 * k] = static_cast<float>(std::cos(angle));
      twiddle[2 * k + 1] = static_cast<float>(std::sin(angle));
    }
    halfTwiddle.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, twiddle.size() * sizeof(float)));
    OpenCL::writeBuffer(halfTwiddle.get(), twiddle.size() * sizeof(float), twiddle.data());
  }

  parameters.scale = static_cast<float>(fixedPointScale);
  parameters.inverseScale = static_cast<float>(1.0 / fixedPointScale);
  parameters.alphaFactor = static_cast<float>(-0.25 / (alpha * alpha));
  parameters.inverseFourAlphaSquared = static_cast<float>(0.25 / (alpha * alpha));
  parameters.meshX = Kx;
  parameters.meshY = Ky;
  parameters.meshZ = Kz;
  parametersChanged = true;

  // fixed kernel arguments
  const cl_mem spreadBuffers[1] = {parameterBuffer.get()};
  check(clSetKernelArg(spreadKernel.get(), 1, sizeof(cl_mem), &spreadBuffers[0]), "clSetKernelArg (spreadCharges)");
  cl_mem meshMem = meshBuffer.get();
  cl_mem potentialMem = potentialBuffer.get();
  check(clSetKernelArg(spreadKernel.get(), 2, sizeof(cl_mem), &meshMem), "clSetKernelArg (spreadCharges)");
  const cl_mem influenceBuffers[6] = {dataBuffer.get(), moduliX.get(),         moduliY.get(),
                                      moduliZ.get(),    parameterBuffer.get(), partialBuffer.get()};
  OpenCLDevice::setBufferArguments(influenceKernel.get(), influenceBuffers, "clSetKernelArg (applyInfluence)");
  check(clSetKernelArg(interpolateKernel.get(), 1, sizeof(cl_mem), &potentialMem),
        "clSetKernelArg (interpolateForces)");
  check(clSetKernelArg(interpolateKernel.get(), 2, sizeof(cl_mem), &spreadBuffers[0]),
        "clSetKernelArg (interpolateForces)");

  // the real z transforms: 0 mesh/spectrum in, 1 spectrum/potential out, 2 twiddle, 3 half twiddle, 4 M,
  // 5 radix code, 6 stages, 7 lines, 8 tile, (9 inverse scale), local buffers
  const std::size_t localBytes = static_cast<std::size_t>(planZ.localLength) * planZ.tile * 2 * sizeof(float);
  const cl_uint lines = Kx * Ky;
  auto setRealArguments = [&](cl_kernel kernel, cl_mem in, cl_mem out, const char* name)
  {
    cl_mem twiddle = planZ.twiddle.get();
    cl_mem half = halfTwiddle.get();
    check(clSetKernelArg(kernel, 0, sizeof(cl_mem), &in), name);
    check(clSetKernelArg(kernel, 1, sizeof(cl_mem), &out), name);
    check(clSetKernelArg(kernel, 2, sizeof(cl_mem), &twiddle), name);
    check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &half), name);
    check(clSetKernelArg(kernel, 4, sizeof(cl_uint), &planZ.N), name);
    check(clSetKernelArg(kernel, 5, sizeof(cl_uint), &planZ.radixCode), name);
    check(clSetKernelArg(kernel, 6, sizeof(cl_uint), &planZ.stages), name);
    check(clSetKernelArg(kernel, 7, sizeof(cl_uint), &lines), name);
    check(clSetKernelArg(kernel, 8, sizeof(cl_uint), &planZ.tile), name);
  };
  setRealArguments(realForwardKernel.get(), meshMem, dataBuffer.get(), "clSetKernelArg (fftRealForward)");
  check(clSetKernelArg(realForwardKernel.get(), 9, sizeof(float), &parameters.inverseScale),
        "clSetKernelArg (fftRealForward)");
  check(clSetKernelArg(realForwardKernel.get(), 10, localBytes, nullptr), "clSetKernelArg (fftRealForward, local A)");
  check(clSetKernelArg(realForwardKernel.get(), 11, localBytes, nullptr), "clSetKernelArg (fftRealForward, local B)");
  setRealArguments(realBackwardKernel.get(), dataBuffer.get(), potentialMem, "clSetKernelArg (fftRealBackward)");
  check(clSetKernelArg(realBackwardKernel.get(), 9, localBytes, nullptr), "clSetKernelArg (fftRealBackward, local A)");
  check(clSetKernelArg(realBackwardKernel.get(), 10, localBytes, nullptr), "clSetKernelArg (fftRealBackward, local B)");
}

void OpenCLMesh::setSlots(std::size_t slots)
{
  const std::uint32_t value = static_cast<std::uint32_t>(slots);
  if (parameters.numberOfSlots != value)
  {
    parameters.numberOfSlots = value;
    parametersChanged = true;
  }
}

void OpenCLMesh::setAlpha(double alphaValue)
{
  if (alphaValue == alpha) return;
  alpha = alphaValue;
  parameters.alphaFactor = static_cast<float>(-0.25 / (alpha * alpha));
  parameters.inverseFourAlphaSquared = static_cast<float>(0.25 / (alpha * alpha));
  parametersChanged = true;
}

void OpenCLMesh::updateBox(const SimulationBox& box)
{
  float inverse[9];
  OpenCLDevice::matrixToFloats(box.inverseCell, inverse);
  const float prefactor = static_cast<float>(conversionFactor * 2.0 * std::numbers::pi / box.volume);
  bool changed = parametersChanged || parameters.prefactor != prefactor;
  for (std::size_t k = 0; k < 9; ++k) changed = changed || parameters.inverseCell[k] != inverse[k];
  if (!changed) return;
  std::copy(std::begin(inverse), std::end(inverse), std::begin(parameters.inverseCell));
  parameters.prefactor = prefactor;
  parametersChanged = true;
}

void OpenCLMesh::enqueueTransform(cl_command_queue queue, const AxisPlan& plan, float sign)
{
  cl_kernel kernel = fftKernel.get();
  cl_mem data = dataBuffer.get();
  cl_mem twiddle = plan.twiddle.get();
  const std::size_t localBytes = static_cast<std::size_t>(plan.localLength) * plan.tile * 2 * sizeof(float);
  check(clSetKernelArg(kernel, 0, sizeof(cl_mem), &data), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 1, sizeof(cl_mem), &twiddle), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 2, sizeof(cl_uint), &plan.N), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 3, sizeof(cl_uint), &plan.radixCode), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 4, sizeof(cl_uint), &plan.stages), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 5, sizeof(cl_uint), &plan.axisStride), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 6, sizeof(cl_uint), &plan.lineStride), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 7, sizeof(cl_uint), &plan.innerCount), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 8, sizeof(cl_uint), &plan.outerStride), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 9, sizeof(cl_uint), &plan.tile), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 10, sizeof(cl_uint), &plan.tilesPerOuter), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 11, sizeof(float), &sign), "clSetKernelArg (fftLines)");
  check(clSetKernelArg(kernel, 12, localBytes, nullptr), "clSetKernelArg (fftLines, local A)");
  check(clSetKernelArg(kernel, 13, localBytes, nullptr), "clSetKernelArg (fftLines, local B)");
  const std::size_t global = plan.groups * plan.groupSize;
  const std::size_t local = plan.groupSize;
  check(clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (fftLines)");
}

void OpenCLMesh::enqueueRealTransform(cl_command_queue queue, bool forward)
{
  cl_kernel kernel = forward ? realForwardKernel.get() : realBackwardKernel.get();
  const std::size_t global = planZ.groups * planZ.groupSize;
  const std::size_t local = planZ.groupSize;
  check(clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, &local, 0, nullptr, nullptr),
        forward ? "clEnqueueNDRangeKernel (fftRealForward)" : "clEnqueueNDRangeKernel (fftRealBackward)");
}

void OpenCLMesh::enqueueSpread(cl_command_queue queue, cl_mem position)
{
  const std::size_t slots = parameters.numberOfSlots;
  const std::size_t slotGlobal = roundUp(std::max<std::size_t>(slots, 1), pointGroup);
  const std::size_t slotLocal = pointGroup;
  check(clSetKernelArg(spreadKernel.get(), 0, sizeof(cl_mem), &position), "clSetKernelArg (spreadCharges)");
  check(clEnqueueNDRangeKernel(queue, spreadKernel.get(), 1, nullptr, &slotGlobal, &slotLocal, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (spreadCharges)");
}

void OpenCLMesh::enqueueForwardTransforms(cl_command_queue queue)
{
  enqueueRealTransform(queue, true);
  enqueueTransform(queue, planY, -1.0f);
  enqueueTransform(queue, planX, -1.0f);
}

void OpenCLMesh::enqueueInfluence(cl_command_queue queue)
{
  const std::size_t influenceGlobal = influenceGroups * influenceGroup;
  const std::size_t influenceLocal = influenceGroup;
  check(clEnqueueNDRangeKernel(queue, influenceKernel.get(), 1, nullptr, &influenceGlobal, &influenceLocal, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel (applyInfluence)");
}

void OpenCLMesh::enqueueBackwardTransforms(cl_command_queue queue)
{
  enqueueTransform(queue, planX, 1.0f);
  enqueueTransform(queue, planY, 1.0f);
  enqueueRealTransform(queue, false);
}

void OpenCLMesh::enqueueInterpolate(cl_command_queue queue, cl_mem position, cl_mem force)
{
  const std::size_t slots = parameters.numberOfSlots;
  const std::size_t slotGlobal = roundUp(std::max<std::size_t>(slots, 1), pointGroup);
  const std::size_t slotLocal = pointGroup;
  check(clSetKernelArg(interpolateKernel.get(), 0, sizeof(cl_mem), &position), "clSetKernelArg (interpolateForces)");
  check(clSetKernelArg(interpolateKernel.get(), 3, sizeof(cl_mem), &force), "clSetKernelArg (interpolateForces)");
  check(
      clEnqueueNDRangeKernel(queue, interpolateKernel.get(), 1, nullptr, &slotGlobal, &slotLocal, 0, nullptr, nullptr),
      "clEnqueueNDRangeKernel (interpolateForces)");
}

void OpenCLMesh::enqueue(cl_command_queue queue, cl_mem position, cl_mem force, bool interpolate)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    check(clEnqueueWriteBuffer(queue, parameterBuffer.get(), CL_FALSE, 0, sizeof(Parameters), &parameters, 0, nullptr,
                               nullptr),
          "clEnqueueWriteBuffer (mesh parameters)");
  }
  enqueueSpread(queue, position);
  enqueueForwardTransforms(queue);
  enqueueInfluence(queue);
  enqueueBackwardTransforms(queue);
  if (interpolate) enqueueInterpolate(queue, position, force);
}

void OpenCLMesh::profile(cl_command_queue queue, cl_mem position, cl_mem force)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    check(clEnqueueWriteBuffer(queue, parameterBuffer.get(), CL_TRUE, 0, sizeof(Parameters), &parameters, 0, nullptr,
                               nullptr),
          "clEnqueueWriteBuffer (mesh parameters)");
  }
  // every stage runs on the state its predecessor leaves (the spreading is reset by the forward transforms, which
  // zero the mesh and are repeated untimed before every timed spreading)
  auto timed = [&](auto&& stage, auto&& reset) -> double
  {
    double best = std::numeric_limits<double>::max();
    for (std::size_t run = 0; run < 3; ++run)
    {
      reset();
      check(clFinish(queue), "clFinish");
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      stage();
      check(clFinish(queue), "clFinish");
      best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    }
    return best;
  };
  auto nothing = [] {};
  stageTimes[0] = timed([&] { enqueueSpread(queue, position); }, [&] { enqueueForwardTransforms(queue); });
  stageTimes[1] = timed([&] { enqueueForwardTransforms(queue); }, [&] { enqueueSpread(queue, position); });
  stageTimes[2] = timed([&] { enqueueInfluence(queue); }, nothing);
  stageTimes[3] = timed([&] { enqueueBackwardTransforms(queue); }, nothing);
  stageTimes[4] = timed([&] { enqueueInterpolate(queue, position, force); }, nothing);
  profiled = true;
}

void OpenCLMesh::enqueueRead(cl_command_queue queue, cl_event* event)
{
  check(clEnqueueReadBuffer(queue, partialBuffer.get(), CL_FALSE, 0, hostPartials.size() * sizeof(float),
                            hostPartials.data(), 0, nullptr, event),
        "clEnqueueReadBuffer (mesh partials)");
}

void OpenCLMesh::collect(double& energy, double3x3& strain, double& singleIonSum, double3x3& singleIonStrain) const
{
  double sums[influencePartials] = {};
  for (std::size_t g = 0; g < influenceGroups; ++g)
  {
    const float* partial = hostPartials.data() + g * influencePartials;
    for (std::size_t q = 0; q < influencePartials; ++q) sums[q] += static_cast<double>(partial[q]);
  }
  // symmetric tensors: xx xy xz yy yz zz
  auto tensor = [&](std::size_t first) -> double3x3
  {
    double3x3 t{};
    t.ax = sums[first];
    t.bx = t.ay = sums[first + 1];
    t.cx = t.az = sums[first + 2];
    t.by = sums[first + 3];
    t.cy = t.bz = sums[first + 4];
    t.cz = sums[first + 5];
    return t;
  };
  energy = sums[0];
  strain = tensor(1);
  singleIonSum = sums[7];
  singleIonStrain = tensor(8);
}

std::string OpenCLMesh::status() const
{
  std::string text = std::format(
      "    particle-mesh Ewald on the device: mesh {} x {} x {}, B-spline order {}, alpha {:.6f} A^-1 (single "
      "precision, fixed-point spreading, real-to-complex FFT of {} x {} x {} spectrum points, tiles {}/{}/{} lines)\n",
      mesh.x, mesh.y, mesh.z, order, alpha, mesh.x, mesh.y, mesh.z / 2 + 1, planX.tile, planY.tile, planZ.tile);
  if (profiled)
  {
    text += std::format(
        "    device mesh stages (sampled, ms): spreading {:.3f}, forward FFTs with the conversion {:.3f}, influence "
        "{:.3f}, backward FFTs {:.3f}, interpolation {:.3f}\n",
        1e3 * stageTimes[0], 1e3 * stageTimes[1], 1e3 * stageTimes[2], 1e3 * stageTimes[3], 1e3 * stageTimes[4]);
  }
  return text;
}
