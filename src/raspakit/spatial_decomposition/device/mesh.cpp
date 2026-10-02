module;

module spatial_decomposition_device_mesh;

import std;

import double3;
import double3x3;
import int3;
import simulationbox;
import spatial_decomposition_pppm;
import spatial_decomposition_device_context;
import spatial_decomposition_device_kernels;

namespace
{
constexpr std::size_t influenceGroup = 128;
constexpr std::size_t influencePartials = 14;
constexpr std::size_t influencePointsPerItem = 8;  // INFLUENCE_POINTS of the kernel source
constexpr std::size_t pointGroup = 64;
constexpr std::size_t maximumTile = 16;  // lines per FFT work-group (a power of two); see planAxis
constexpr double fixedPointScale = 16777216.0;  // 2^24: 6e-8 e resolution, +-128 e range per mesh point

std::size_t divideRoundUp(std::size_t a, std::size_t b) { return (a + b - 1) / b; }
std::size_t roundUp(std::size_t value, std::size_t multiple) { return ((value + multiple - 1) / multiple) * multiple; }

template <typename Matrix>
void matrixToFloats(const Matrix& m, float* out)
{
  out[0] = static_cast<float>(m.ax);
  out[1] = static_cast<float>(m.ay);
  out[2] = static_cast<float>(m.az);
  out[3] = static_cast<float>(m.bx);
  out[4] = static_cast<float>(m.by);
  out[5] = static_cast<float>(m.bz);
  out[6] = static_cast<float>(m.cx);
  out[7] = static_cast<float>(m.cy);
  out[8] = static_cast<float>(m.cz);
}

/// A device buffer holding a host table (blocking upload).
void uploadTable(DeviceContext& context, DeviceBufferOwner& buffer, const void* data, std::size_t bytes)
{
  buffer.allocate(context, bytes, DeviceMemory::Device);
  if (bytes > 0) context.write(buffer.get(), 0, bytes, data, true);
}
}  // namespace

void DeviceMesh::initialize(DeviceContext& deviceContext, std::size_t orderValue)
{
  context = &deviceContext;
  order = std::clamp<std::size_t>(orderValue, 3, 7);
  const std::string program = std::format("mesh{}", order);
  const std::string source = std::format("#define MESH_ORDER {}\n", order) + deviceKernelMeshSource;
  auto kernel = [&](const char* name)
  { return context->compileKernel(program, source.c_str(), DeviceMath::Relaxed, name); };
  spreadKernel = kernel("spreadCharges");
  realForwardKernel = kernel("fftRealForward");
  realBackwardKernel = kernel("fftRealBackward");
  fftKernel = kernel("fftLines");
  influenceKernel = kernel("applyInfluence");
  interpolateKernel = kernel("interpolateForces");

  if (context->localMemorySize() > 0) localMemory = context->localMemorySize();
  const std::size_t groupLimit = context->maxGroupSize(fftKernel);
  if (groupLimit > 0) maxGroupSize = std::min<std::size_t>(256, groupLimit);
  parameterBuffer.allocate(*context, sizeof(Parameters), DeviceMemory::Device);
}

void DeviceMesh::planAxis(AxisPlan& plan, std::uint32_t N, std::uint32_t localLength, std::uint32_t axisStride,
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
    throw std::runtime_error(std::format("[Device mesh]: mesh size {} is not a product of 2, 3 and 5\n", N));
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
        std::format("[Device mesh]: a mesh axis of {} points does not fit the local memory of the device ({} bytes)\n",
                    N, localMemory));
  }
  // a power of two (the kernels index the tile with shifts), at most maximumTile: smaller tiles let several
  // work-groups share a compute unit's local memory
  std::size_t tile = std::clamp<std::size_t>(budget / perLine, 1, maximumTile);
  tile = std::min<std::size_t>(tile, innerCount);
  plan.tileShift = 0;
  while ((std::size_t{2} << plan.tileShift) <= tile) ++plan.tileShift;
  plan.tile = std::uint32_t{1} << plan.tileShift;
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
  uploadTable(*context, plan.twiddle, twiddle.data(), twiddle.size() * sizeof(float));
}

void DeviceMesh::setup(int3 meshSize, double alphaValue, double conversionFactorValue)
{
  mesh = meshSize;
  // the real-to-complex transform along z packs pairs of samples: Kz must be even
  while (mesh.z % 2 != 0)
    mesh.z = static_cast<std::int32_t>(PPPM::nextFFTFriendly(static_cast<std::size_t>(mesh.z) + 1));
  alpha = alphaValue;
  conversionFactor = conversionFactorValue;
  const std::size_t points = meshPoints();
  const std::size_t spectrum = spectrumPoints();

  // the fixed-point mesh starts zeroed (the forward transform zeroes it again after every step)
  {
    const std::vector<std::int32_t> zeros(points, 0);
    uploadTable(*context, meshBuffer, zeros.data(), points * sizeof(std::int32_t));
  }
  dataBuffer.allocate(*context, spectrum * 2 * sizeof(float), DeviceMemory::Device);
  potentialBuffer.allocate(*context, points * sizeof(float), DeviceMemory::Device);
  influenceGroups = divideRoundUp(spectrum, influenceGroup * influencePointsPerItem);
  partialBuffer.allocate(*context, influenceGroups * influencePartials * sizeof(float), DeviceMemory::Device);
  hostPartials.assign(influenceGroups * influencePartials, 0.0f);

  auto uploadModuli = [&](DeviceBufferOwner& buffer, std::int32_t K)
  {
    const std::vector<double> moduli = PPPM::bsplineModuli(order, K);
    std::vector<float> values(moduli.size());
    for (std::size_t m = 0; m < moduli.size(); ++m) values[m] = static_cast<float>(moduli[m]);
    uploadTable(*context, buffer, values.data(), values.size() * sizeof(float));
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
  lines = Kx * Ky;
  {
    std::vector<float> twiddle(2 * static_cast<std::size_t>(Hz));
    for (std::size_t k = 0; k < Hz; ++k)
    {
      const double angle = -2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(Kz);
      twiddle[2 * k] = static_cast<float>(std::cos(angle));
      twiddle[2 * k + 1] = static_cast<float>(std::sin(angle));
    }
    uploadTable(*context, halfTwiddle, twiddle.data(), twiddle.size() * sizeof(float));
  }

  parameters.scale = static_cast<float>(fixedPointScale);
  parameters.inverseScale = static_cast<float>(1.0 / fixedPointScale);
  parameters.alphaFactor = static_cast<float>(-0.25 / (alpha * alpha));
  parameters.inverseFourAlphaSquared = static_cast<float>(0.25 / (alpha * alpha));
  parameters.meshX = Kx;
  parameters.meshY = Ky;
  parameters.meshZ = Kz;
  parametersChanged = true;
}

void DeviceMesh::setSlots(std::size_t slots)
{
  const std::uint32_t value = static_cast<std::uint32_t>(slots);
  if (parameters.numberOfSlots != value)
  {
    parameters.numberOfSlots = value;
    parametersChanged = true;
  }
}

void DeviceMesh::setAlpha(double alphaValue)
{
  if (alphaValue == alpha) return;
  alpha = alphaValue;
  parameters.alphaFactor = static_cast<float>(-0.25 / (alpha * alpha));
  parameters.inverseFourAlphaSquared = static_cast<float>(0.25 / (alpha * alpha));
  parametersChanged = true;
}

void DeviceMesh::updateBox(const SimulationBox& box)
{
  float inverse[9];
  matrixToFloats(box.inverseCell, inverse);
  const float prefactor = static_cast<float>(conversionFactor * 2.0 * std::numbers::pi / box.volume);
  bool changed = parametersChanged || parameters.prefactor != prefactor;
  for (std::size_t k = 0; k < 9; ++k) changed = changed || parameters.inverseCell[k] != inverse[k];
  if (!changed) return;
  std::copy(std::begin(inverse), std::end(inverse), std::begin(parameters.inverseCell));
  parameters.prefactor = prefactor;
  parametersChanged = true;
}

void DeviceMesh::writeParameters(bool blocking)
{
  if (!parametersChanged) return;
  parametersChanged = false;
  context->write(parameterBuffer.get(), 0, sizeof(Parameters), &parameters, blocking);
}

void DeviceMesh::enqueueTransform(const AxisPlan& plan, float sign)
{
  const std::size_t localBytes = static_cast<std::size_t>(plan.localLength) * plan.tile * 2 * sizeof(float);
  const DeviceArg arguments[] = {DeviceArg::of(dataBuffer.get()),      DeviceArg::of(plan.twiddle.get()),
                                 DeviceArg::value(plan.N),             DeviceArg::value(plan.radixCode),
                                 DeviceArg::value(plan.stages),        DeviceArg::value(plan.axisStride),
                                 DeviceArg::value(plan.lineStride),    DeviceArg::value(plan.innerCount),
                                 DeviceArg::value(plan.outerStride),   DeviceArg::value(plan.tileShift),
                                 DeviceArg::value(plan.tilesPerOuter), DeviceArg::value(sign),
                                 DeviceArg::local(localBytes),         DeviceArg::local(localBytes)};
  context->launch(fftKernel, arguments, plan.groups, plan.groupSize);
}

void DeviceMesh::enqueueRealTransform(bool forward)
{
  const std::size_t localBytes = static_cast<std::size_t>(planZ.localLength) * planZ.tile * 2 * sizeof(float);
  if (forward)
  {
    const DeviceArg arguments[] = {DeviceArg::of(meshBuffer.get()),     DeviceArg::of(dataBuffer.get()),
                                   DeviceArg::of(planZ.twiddle.get()),  DeviceArg::of(halfTwiddle.get()),
                                   DeviceArg::value(planZ.N),           DeviceArg::value(planZ.radixCode),
                                   DeviceArg::value(planZ.stages),      DeviceArg::value(lines),
                                   DeviceArg::value(planZ.tileShift),   DeviceArg::value(parameters.inverseScale),
                                   DeviceArg::local(localBytes),        DeviceArg::local(localBytes)};
    context->launch(realForwardKernel, arguments, planZ.groups, planZ.groupSize);
  }
  else
  {
    const DeviceArg arguments[] = {DeviceArg::of(dataBuffer.get()),    DeviceArg::of(potentialBuffer.get()),
                                   DeviceArg::of(planZ.twiddle.get()), DeviceArg::of(halfTwiddle.get()),
                                   DeviceArg::value(planZ.N),          DeviceArg::value(planZ.radixCode),
                                   DeviceArg::value(planZ.stages),     DeviceArg::value(lines),
                                   DeviceArg::value(planZ.tileShift),  DeviceArg::local(localBytes),
                                   DeviceArg::local(localBytes)};
    context->launch(realBackwardKernel, arguments, planZ.groups, planZ.groupSize);
  }
}

void DeviceMesh::enqueueSpread(DeviceBuffer position)
{
  const std::size_t slots = parameters.numberOfSlots;
  const std::size_t groups = roundUp(std::max<std::size_t>(slots, 1), pointGroup) / pointGroup;
  const DeviceArg arguments[] = {DeviceArg::of(position), DeviceArg::of(parameterBuffer.get()),
                                 DeviceArg::of(meshBuffer.get())};
  context->launch(spreadKernel, arguments, groups, pointGroup);
}

void DeviceMesh::enqueueForwardTransforms()
{
  enqueueRealTransform(true);
  enqueueTransform(planY, -1.0f);
  enqueueTransform(planX, -1.0f);
}

void DeviceMesh::enqueueInfluence()
{
  const DeviceArg arguments[] = {DeviceArg::of(dataBuffer.get()),      DeviceArg::of(moduliX.get()),
                                 DeviceArg::of(moduliY.get()),         DeviceArg::of(moduliZ.get()),
                                 DeviceArg::of(parameterBuffer.get()), DeviceArg::of(partialBuffer.get())};
  context->launch(influenceKernel, arguments, influenceGroups, influenceGroup);
}

void DeviceMesh::enqueueBackwardTransforms()
{
  enqueueTransform(planX, 1.0f);
  enqueueTransform(planY, 1.0f);
  enqueueRealTransform(false);
}

void DeviceMesh::enqueueInterpolate(DeviceBuffer position, DeviceBuffer force)
{
  const std::size_t slots = parameters.numberOfSlots;
  const std::size_t groups = roundUp(std::max<std::size_t>(slots, 1), pointGroup) / pointGroup;
  const DeviceArg arguments[] = {DeviceArg::of(position), DeviceArg::of(potentialBuffer.get()),
                                 DeviceArg::of(parameterBuffer.get()), DeviceArg::of(force)};
  context->launch(interpolateKernel, arguments, groups, pointGroup);
}

void DeviceMesh::enqueue(DeviceBuffer position, DeviceBuffer force, bool interpolate)
{
  writeParameters(false);
  enqueueSpread(position);
  enqueueForwardTransforms();
  enqueueInfluence();
  enqueueBackwardTransforms();
  if (interpolate) enqueueInterpolate(position, force);
}

void DeviceMesh::profile(DeviceBuffer position, DeviceBuffer force)
{
  writeParameters(true);
  // every stage runs on the state its predecessor leaves (the spreading is reset by the forward transforms, which
  // zero the mesh and are repeated untimed before every timed spreading)
  auto timed = [&](auto&& stage, auto&& reset) -> double
  {
    double best = std::numeric_limits<double>::max();
    for (std::size_t run = 0; run < 3; ++run)
    {
      reset();
      context->finish();
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      stage();
      context->finish();
      best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    }
    return best;
  };
  auto nothing = [] {};
  stageTimes[0] = timed([&] { enqueueSpread(position); }, [&] { enqueueForwardTransforms(); });
  stageTimes[1] = timed([&] { enqueueForwardTransforms(); }, [&] { enqueueSpread(position); });
  stageTimes[2] = timed([&] { enqueueInfluence(); }, nothing);
  stageTimes[3] = timed([&] { enqueueBackwardTransforms(); }, nothing);
  stageTimes[4] = timed([&] { enqueueInterpolate(position, force); }, nothing);
  profiled = true;
}

void DeviceMesh::enqueueRead()
{
  context->read(partialBuffer.get(), 0, hostPartials.size() * sizeof(float), hostPartials.data());
}

void DeviceMesh::collect(double& energy, double3x3& strain, double& singleIonSum, double3x3& singleIonStrain) const
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

std::string DeviceMesh::status() const
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
