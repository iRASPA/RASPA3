module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module spatial_decomposition_opencl_bonded;

import std;

import double3;
import double3x3;
import opencl;
import spatial_decomposition_opencl_handles;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_backend;
import spatial_decomposition_device_bonded_topology;

using OpenCLDevice::check;
using OpenCLDevice::roundUp;

namespace
{
constexpr std::size_t groupSize = 64;     // BONDED_GROUP and TERM_GROUP of the kernel source
constexpr std::size_t atomPartials = 20;  // ATOM_PARTIALS
constexpr std::size_t termPartials = 8;   // TERM_PARTIALS

template <typename Buffer>
void ensureCapacity(Buffer& buffer, std::size_t& capacity, std::size_t required, std::size_t elementBytes,
                    cl_mem_flags flags)
{
  if (required <= capacity) return;
  capacity = required + required / 4 + 1;
  buffer.reset(OpenCL::createBuffer(flags, capacity * elementBytes));
}

// blocking upload of a static table (the shared queue; complete before the step's queue uses it)
template <typename Value>
void uploadTable(OpenCLDevice::MemHandle& buffer, const std::vector<Value>& values)
{
  const std::size_t bytes = std::max<std::size_t>(1, values.size()) * sizeof(Value);
  buffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, bytes));
  if (!values.empty()) OpenCL::writeBuffer(buffer.get(), values.size() * sizeof(Value), values.data());
}
}  // namespace

void OpenCLBonded::initialize(cl_context clContext, cl_device_id clDevice)
{
  context = clContext;
  device = clDevice;
  program = OpenCLDevice::buildDeviceProgram(context, device, deviceKernelBondedSource,
                                             "-cl-mad-enable -cl-no-signed-zeros", "OpenCL bonded");
  termKernel = OpenCLDevice::createKernel(program.get(), "bondedTerms");
  atomKernel = OpenCLDevice::createKernel(program.get(), "bondedAtoms");
  parameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(Parameters)));
}

void OpenCLBonded::setTopology(const BondedTopology& topology)
{
  numberOfMolecules = topology.molecules.size();
  numberOfTerms = topology.terms.size();
  chargeSquaredSum = topology.chargeSquaredSum;
  parameters.numberOfInstances = static_cast<std::uint32_t>(topology.numberOfInstances);
  parametersChanged = true;
  termGroups = roundUp(std::max<std::size_t>(topology.numberOfInstances, 1), groupSize) / groupSize;

  uploadTable(termsBuffer, topology.terms);
  uploadTable(gradientOffsetBuffer, topology.gradientOffset);
  uploadTable(atomGradientStartBuffer, topology.atomGradientStart);
  uploadTable(atomGradientsBuffer, topology.atomGradients);
  uploadTable(instanceMoleculeBuffer, topology.instanceMolecule);
  uploadTable(moleculeInfoBuffer, topology.molecules);
  uploadTable(massBuffer, topology.massOfAtom);
  termGradientBuffer.reset(OpenCL::createBuffer(
      CL_MEM_READ_WRITE, std::max<std::size_t>(1, topology.numberOfGradients) * 4 * sizeof(float)));
}

void OpenCLBonded::setParameters(double alpha, double factor, bool useCharge)
{
  alphaValue = alpha;
  conversionFactor = factor;
  parameters.alpha = static_cast<float>(alpha);
  parameters.twoAlphaOverSqrtPi = static_cast<float>(2.0 * alpha * std::numbers::inv_sqrtpi);
  parameters.selfPrefactor = static_cast<float>(conversionFactor * alpha * std::numbers::inv_sqrtpi);
  parameters.coulombFactor = static_cast<float>(conversionFactor);
  parameters.useCharge = useCharge ? 1u : 0u;
  parametersChanged = true;
}

void OpenCLBonded::setAlpha(double alpha)
{
  alphaValue = alpha;
  const float value = static_cast<float>(alpha);
  if (value == parameters.alpha) return;
  parameters.alpha = value;
  parameters.twoAlphaOverSqrtPi = static_cast<float>(2.0 * alpha * std::numbers::inv_sqrtpi);
  parameters.selfPrefactor = static_cast<float>(conversionFactor * alpha * std::numbers::inv_sqrtpi);
  parametersChanged = true;
}

void OpenCLBonded::setLayout(cl_command_queue queue, std::span<const std::uint32_t> slotMolecule)
{
  const std::size_t slots = slotMolecule.size();
  atomGroups = roundUp(std::max<std::size_t>(slots, 1), groupSize) / groupSize;
  const std::uint32_t atomPartialOffset = static_cast<std::uint32_t>(termGroups * termPartials);
  if (parameters.numberOfSlots != static_cast<std::uint32_t>(slots) ||
      parameters.atomPartialOffset != atomPartialOffset)
  {
    parameters.numberOfSlots = static_cast<std::uint32_t>(slots);
    parameters.atomPartialOffset = atomPartialOffset;
    parametersChanged = true;
  }

  ensureCapacity(slotMoleculeBuffer, slotCapacity, slots, sizeof(std::uint32_t), CL_MEM_READ_ONLY);
  // one buffer (one read-back) for the partials of both kernels
  const std::size_t partials = termGroups * termPartials + atomGroups * atomPartials;
  if (partials > partialCapacity)
  {
    partialCapacity = partials + partials / 4 + 1;
    partialBuffer.reset(OpenCL::createBuffer(CL_MEM_WRITE_ONLY, partialCapacity * sizeof(float)));
  }
  hostPartials.assign(partials, 0.0f);

  check(clEnqueueWriteBuffer(queue, slotMoleculeBuffer.get(), CL_FALSE, 0, slots * sizeof(std::uint32_t),
                             slotMolecule.data(), 0, nullptr, nullptr),
        "clEnqueueWriteBuffer (slot molecules)");
}

void OpenCLBonded::enqueue(cl_command_queue queue, cl_mem relative, cl_mem force)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    check(clEnqueueWriteBuffer(queue, parameterBuffer.get(), CL_FALSE, 0, sizeof(Parameters), &parameters, 0, nullptr,
                               nullptr),
          "clEnqueueWriteBuffer (bonded parameters)");
  }
  const std::size_t local = groupSize;
  if (parameters.numberOfInstances > 0)
  {
    const cl_mem termBuffers[8] = {relative,
                                   instanceMoleculeBuffer.get(),
                                   moleculeInfoBuffer.get(),
                                   termsBuffer.get(),
                                   gradientOffsetBuffer.get(),
                                   parameterBuffer.get(),
                                   termGradientBuffer.get(),
                                   partialBuffer.get()};
    OpenCLDevice::setBufferArguments(termKernel.get(), termBuffers, "clSetKernelArg (bondedTerms)");
    const std::size_t termGlobal = termGroups * groupSize;
    check(clEnqueueNDRangeKernel(queue, termKernel.get(), 1, nullptr, &termGlobal, &local, 0, nullptr, nullptr),
          "clEnqueueNDRangeKernel (bondedTerms)");
  }
  const cl_mem atomBuffers[10] = {relative,
                                  slotMoleculeBuffer.get(),
                                  moleculeInfoBuffer.get(),
                                  massBuffer.get(),
                                  atomGradientStartBuffer.get(),
                                  atomGradientsBuffer.get(),
                                  termGradientBuffer.get(),
                                  parameterBuffer.get(),
                                  force,
                                  partialBuffer.get()};
  OpenCLDevice::setBufferArguments(atomKernel.get(), atomBuffers, "clSetKernelArg (bondedAtoms)");
  const std::size_t atomGlobal = atomGroups * groupSize;
  check(clEnqueueNDRangeKernel(queue, atomKernel.get(), 1, nullptr, &atomGlobal, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (bondedAtoms)");
}

void OpenCLBonded::enqueueRead(cl_command_queue queue, cl_event* event)
{
  check(clEnqueueReadBuffer(queue, partialBuffer.get(), CL_FALSE, 0, hostPartials.size() * sizeof(float),
                            hostPartials.data(), 0, nullptr, event),
        "clEnqueueReadBuffer (bonded partials)");
}

DeviceBondedResults OpenCLBonded::collect() const
{
  // the term partials (per kind) followed by the atom partials
  double termSums[termPartials] = {};
  if (parameters.numberOfInstances > 0)
  {
    for (std::size_t g = 0; g < termGroups; ++g)
    {
      const float* partial = hostPartials.data() + g * termPartials;
      for (std::size_t q = 0; q < termPartials; ++q) termSums[q] += static_cast<double>(partial[q]);
    }
  }
  double sums[atomPartials] = {};
  const float* atomBase = hostPartials.data() + parameters.atomPartialOffset;
  for (std::size_t g = 0; g < atomGroups; ++g)
  {
    const float* partial = atomBase + g * atomPartials;
    for (std::size_t q = 0; q < atomPartials; ++q) sums[q] += static_cast<double>(partial[q]);
  }
  auto tensor = [&](std::size_t first) -> double3x3
  {
    double3x3 t{};
    t.ax = sums[first];
    t.ay = sums[first + 1];
    t.az = sums[first + 2];
    t.bx = sums[first + 3];
    t.by = sums[first + 4];
    t.bz = sums[first + 5];
    t.cx = sums[first + 6];
    t.cy = sums[first + 7];
    t.cz = sums[first + 8];
    return t;
  };
  DeviceBondedResults results{};
  if (parameters.useCharge != 0)
  {
    results.self = -conversionFactor * alphaValue * std::numbers::inv_sqrtpi * chargeSquaredSum;
    results.exclusion = (sums[0] + sums[1]) - results.self;
  }
  results.exclusionStrain = tensor(2);
  results.correction = tensor(11);
  results.bond = termSums[0];
  results.bend = termSums[1];
  results.torsion = termSums[2];
  results.improperTorsion = termSums[3];
  results.intraVDW = termSums[4];
  results.intraCoulomb = termSums[5];
  return results;
}

std::string OpenCLBonded::status() const
{
  return std::format(
      "    bonded terms on the device: {} molecules, {} bonded and intramolecular pair terms over the components "
      "({} instances), self and exclusion corrections, virial correction (single precision, positions relative to "
      "the molecule)\n",
      numberOfMolecules, numberOfTerms, parameters.numberOfInstances);
}
