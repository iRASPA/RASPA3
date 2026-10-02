module;

module spatial_decomposition_device_bonded;

import std;

import double3;
import double3x3;
import spatial_decomposition_device_context;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_bonded_topology;

namespace
{
constexpr std::size_t groupSize = 64;     // BONDED_GROUP and TERM_GROUP of the kernel source
constexpr std::size_t atomPartials = 20;  // ATOM_PARTIALS
constexpr std::size_t termPartials = 8;   // TERM_PARTIALS

std::size_t roundUp(std::size_t value, std::size_t multiple) { return ((value + multiple - 1) / multiple) * multiple; }

// blocking upload of a static table (complete before the step uses it)
template <typename Value>
void uploadTable(DeviceContext& context, DeviceBufferOwner& buffer, const std::vector<Value>& values)
{
  buffer.allocate(context, values.size() * sizeof(Value), DeviceMemory::Device);
  if (!values.empty()) context.write(buffer.get(), 0, values.size() * sizeof(Value), values.data(), true);
}
}  // namespace

void DeviceBonded::initialize(DeviceContext& deviceContext)
{
  context = &deviceContext;
  termKernel = context->compileKernel("bonded", deviceKernelBondedSource, DeviceMath::Relaxed, "bondedTerms");
  atomKernel = context->compileKernel("bonded", deviceKernelBondedSource, DeviceMath::Relaxed, "bondedAtoms");
  parameterBuffer.allocate(*context, sizeof(Parameters), DeviceMemory::Device);
}

void DeviceBonded::setTopology(const BondedTopology& topology)
{
  numberOfMolecules = topology.molecules.size();
  numberOfTerms = topology.terms.size();
  chargeSquaredSum = topology.chargeSquaredSum;
  parameters.numberOfInstances = static_cast<std::uint32_t>(topology.numberOfInstances);
  parametersChanged = true;
  termGroups = roundUp(std::max<std::size_t>(topology.numberOfInstances, 1), groupSize) / groupSize;

  uploadTable(*context, termsBuffer, topology.terms);
  uploadTable(*context, gradientOffsetBuffer, topology.gradientOffset);
  uploadTable(*context, atomGradientStartBuffer, topology.atomGradientStart);
  uploadTable(*context, atomGradientsBuffer, topology.atomGradients);
  uploadTable(*context, instanceMoleculeBuffer, topology.instanceMolecule);
  uploadTable(*context, moleculeInfoBuffer, topology.molecules);
  uploadTable(*context, massBuffer, topology.massOfAtom);
  termGradientBuffer.allocate(*context, topology.numberOfGradients * 4 * sizeof(float), DeviceMemory::Device);
}

void DeviceBonded::setParameters(double alpha, double factor, bool useCharge)
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

void DeviceBonded::setAlpha(double alpha)
{
  alphaValue = alpha;
  const float value = static_cast<float>(alpha);
  if (value == parameters.alpha) return;
  parameters.alpha = value;
  parameters.twoAlphaOverSqrtPi = static_cast<float>(2.0 * alpha * std::numbers::inv_sqrtpi);
  parameters.selfPrefactor = static_cast<float>(conversionFactor * alpha * std::numbers::inv_sqrtpi);
  parametersChanged = true;
}

void DeviceBonded::setLayout(std::span<const std::uint32_t> slotMolecule)
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

  if (slots > slotCapacity)
  {
    slotCapacity = slots + slots / 4 + 1;
    slotMoleculeBuffer.allocate(*context, slotCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
  }
  // one buffer (one read-back) for the partials of both kernels
  const std::size_t partials = termGroups * termPartials + atomGroups * atomPartials;
  if (partials > partialCapacity)
  {
    partialCapacity = partials + partials / 4 + 1;
    partialBuffer.allocate(*context, partialCapacity * sizeof(float), DeviceMemory::Device);
  }
  hostPartials.assign(partials, 0.0f);

  if (slots > 0) context->write(slotMoleculeBuffer.get(), 0, slots * sizeof(std::uint32_t), slotMolecule.data(), false);
}

void DeviceBonded::enqueue(DeviceBuffer relative, DeviceBuffer force)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    context->write(parameterBuffer.get(), 0, sizeof(Parameters), &parameters, false);
  }
  if (parameters.numberOfInstances > 0)
  {
    const DeviceArg arguments[] = {DeviceArg::of(relative),
                                   DeviceArg::of(instanceMoleculeBuffer.get()),
                                   DeviceArg::of(moleculeInfoBuffer.get()),
                                   DeviceArg::of(termsBuffer.get()),
                                   DeviceArg::of(gradientOffsetBuffer.get()),
                                   DeviceArg::of(parameterBuffer.get()),
                                   DeviceArg::of(termGradientBuffer.get()),
                                   DeviceArg::of(partialBuffer.get())};
    context->launch(termKernel, arguments, termGroups, groupSize);
  }
  const DeviceArg arguments[] = {DeviceArg::of(relative),
                                 DeviceArg::of(slotMoleculeBuffer.get()),
                                 DeviceArg::of(moleculeInfoBuffer.get()),
                                 DeviceArg::of(massBuffer.get()),
                                 DeviceArg::of(atomGradientStartBuffer.get()),
                                 DeviceArg::of(atomGradientsBuffer.get()),
                                 DeviceArg::of(termGradientBuffer.get()),
                                 DeviceArg::of(parameterBuffer.get()),
                                 DeviceArg::of(force),
                                 DeviceArg::of(partialBuffer.get())};
  context->launch(atomKernel, arguments, atomGroups, groupSize);
}

void DeviceBonded::enqueueRead()
{
  context->read(partialBuffer.get(), 0, hostPartials.size() * sizeof(float), hostPartials.data());
}

DeviceBondedResults DeviceBonded::collect() const
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

std::string DeviceBonded::status() const
{
  return std::format(
      "    bonded terms on the device: {} molecules, {} bonded and intramolecular pair terms over the components "
      "({} instances), self and exclusion corrections, virial correction (single precision, positions relative to "
      "the molecule)\n",
      numberOfMolecules, numberOfTerms, parameters.numberOfInstances);
}
