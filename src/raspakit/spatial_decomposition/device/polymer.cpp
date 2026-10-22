module;

module spatial_decomposition_device_polymer;

import std;

import double3x3;
import spatial_decomposition_device_context;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_polymer_topology;

namespace
{
constexpr std::size_t groupSize = 64;     // POLYMER_GROUP of the kernel source
constexpr std::size_t termPartials = 17;  // POLYMER_TERM_PARTIALS
constexpr std::size_t atomPartials = 9;   // POLYMER_ATOM_PARTIALS

std::size_t groupsFor(std::size_t count) { return (std::max<std::size_t>(count, 1) + groupSize - 1) / groupSize; }

template <typename Value>
void uploadTable(DeviceContext& context, DeviceBufferOwner& buffer, const std::vector<Value>& values)
{
  buffer.allocate(context, values.size() * sizeof(Value), DeviceMemory::Device);
  if (!values.empty()) context.write(buffer.get(), 0, values.size() * sizeof(Value), values.data(), true);
}
}  // namespace

void DevicePolymerOverlay::initialize(DeviceContext& deviceContext)
{
  context = &deviceContext;
  // the overlay source reuses the term functions of the bonded source: one program of both
  const std::string source = std::string(deviceKernelBondedSource) + deviceKernelPolymerSource;
  termKernel = context->compileKernel("polymer", source.c_str(), DeviceMath::Relaxed, "polymerTerms");
  atomKernel = context->compileKernel("polymer", source.c_str(), DeviceMath::Relaxed, "polymerAtoms");
  parameterBuffer.allocate(*context, sizeof(Parameters), DeviceMemory::Device);
  lambdaBuffer.allocate(*context, sizeof(float), DeviceMemory::Device);
  lambdaCapacity = 1;
}

void DevicePolymerOverlay::setTopology(const PolymerOverlayTopology& topology)
{
  terms = topology.terms.size();
  atoms = topology.atoms.size();
  ramps = topology.rampSteps.size();
  termGroups = groupsFor(terms);
  atomGroups = groupsFor(atoms);
  parameters.numberOfTerms = static_cast<std::uint32_t>(terms);
  parameters.numberOfAtoms = static_cast<std::uint32_t>(atoms);
  parametersChanged = true;

  uploadTable(*context, termsBuffer, topology.terms);
  uploadTable(*context, atomsBuffer, topology.atoms);
  uploadTable(*context, refStartBuffer, topology.refStart);
  uploadTable(*context, refLinkEndBuffer, topology.refLinkEnd);
  uploadTable(*context, refsBuffer, topology.refs);
  termGradientBuffer.allocate(*context, std::max<std::size_t>(terms, 1) * 4 * 4 * sizeof(float), DeviceMemory::Device);
  if (ramps > lambdaCapacity)
  {
    lambdaCapacity = ramps + ramps / 4 + 1;
    lambdaBuffer.allocate(*context, lambdaCapacity * sizeof(float), DeviceMemory::Device);
  }
  // the two kernels reduce into their own buffers (both index the partials from their group zero)
  const std::size_t termFloats = termGroups * termPartials;
  const std::size_t atomFloats = atomGroups * atomPartials;
  if (termFloats > termPartialCapacity)
  {
    termPartialCapacity = termFloats + termFloats / 4 + 1;
    termPartialBuffer.allocate(*context, termPartialCapacity * sizeof(float), DeviceMemory::Device);
  }
  if (atomFloats > atomPartialCapacity)
  {
    atomPartialCapacity = atomFloats + atomFloats / 4 + 1;
    atomPartialBuffer.allocate(*context, atomPartialCapacity * sizeof(float), DeviceMemory::Device);
  }
  hostPartials.assign(termFloats + atomFloats, 0.0f);
}

void DevicePolymerOverlay::setLambdas(std::span<const float> lambdas)
{
  if (lambdas.empty()) return;
  // non-blocking (the step is waited for before the next write): two host copies, so that a write enqueued
  // right before this one keeps its memory
  std::vector<float>& host = lambdaHost[lambdaTurn];
  lambdaTurn ^= 1;
  host.assign(lambdas.begin(), lambdas.end());
  context->write(lambdaBuffer.get(), 0, host.size() * sizeof(float), host.data(), false);
}

void DevicePolymerOverlay::setPairSetup(const DevicePolymerPairSetup& pairSetup)
{
  if (pairSetup == setup && setupApplied) return;
  setupApplied = true;
  setup = pairSetup;
  std::copy(std::begin(setup.cell), std::end(setup.cell), std::begin(parameters.cell));
  std::copy(std::begin(setup.inverseCell), std::end(setup.inverseCell), std::begin(parameters.inverseCell));
  parameters.cutOffVDWSquared = setup.cutOffVDWSquared;
  parameters.cutOffChargeSquared = setup.cutOffChargeSquared;
  parameters.alpha = setup.alpha;
  parameters.alphaSquared = setup.alpha * setup.alpha;
  parameters.alphaOverSqrtPi = static_cast<float>(static_cast<double>(setup.alpha) * std::numbers::inv_sqrtpi);
  parameters.twoAlphaOverSqrtPi = static_cast<float>(2.0 * static_cast<double>(setup.alpha) * std::numbers::inv_sqrtpi);
  parameters.coulombFactor = setup.coulombFactor;
  parameters.useCharge = setup.useCharge;
  parameters.orthorhombic = setup.orthorhombic;
  parameters.switchMode = setup.switchMode;
  parameters.switchDistanceSquared = setup.switchDistanceSquared;
  parameters.switchDistance = setup.switchDistance;
  parameters.switchInverseWidth = setup.switchInverseWidth;
  parameters.switchInverseCutOff3 = setup.switchInverseCutOff3;
  parameters.switchA12 = setup.switchA12;
  parameters.switchA6 = setup.switchA6;
  parametersChanged = true;
}

void DevicePolymerOverlay::enqueue(DeviceBuffer position, DeviceBuffer slotOfAtom, DeviceBuffer relative,
                                   DeviceBuffer moleculeOfAtom, DeviceBuffer moleculeCenter, DeviceBuffer force)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    context->write(parameterBuffer.get(), 0, sizeof(Parameters), &parameters, false);
  }
  if (terms == 0) return;
  {
    const DeviceArg arguments[] = {DeviceArg::of(position),
                                   DeviceArg::of(slotOfAtom),
                                   DeviceArg::of(termsBuffer.get()),
                                   DeviceArg::of(lambdaBuffer.get()),
                                   DeviceArg::of(parameterBuffer.get()),
                                   DeviceArg::of(termGradientBuffer.get()),
                                   DeviceArg::of(termPartialBuffer.get())};
    context->launch(termKernel, arguments, termGroups, groupSize);
  }
  {
    const DeviceArg arguments[] = {DeviceArg::of(atomsBuffer.get()),
                                   DeviceArg::of(refStartBuffer.get()),
                                   DeviceArg::of(refLinkEndBuffer.get()),
                                   DeviceArg::of(refsBuffer.get()),
                                   DeviceArg::of(termGradientBuffer.get()),
                                   DeviceArg::of(slotOfAtom),
                                   DeviceArg::of(relative),
                                   DeviceArg::of(moleculeOfAtom),
                                   DeviceArg::of(moleculeCenter),
                                   DeviceArg::of(parameterBuffer.get()),
                                   DeviceArg::of(force),
                                   DeviceArg::of(atomPartialBuffer.get())};
    context->launch(atomKernel, arguments, atomGroups, groupSize);
  }
}

void DevicePolymerOverlay::enqueueRead()
{
  if (terms == 0) return;
  // the atom partials follow the term partials in the host copy
  const std::size_t termFloats = termGroups * termPartials;
  const std::size_t atomFloats = atomGroups * atomPartials;
  context->read(termPartialBuffer.get(), 0, termFloats * sizeof(float), hostPartials.data());
  context->read(atomPartialBuffer.get(), 0, atomFloats * sizeof(float), hostPartials.data() + termFloats);
}

DevicePolymerResults DevicePolymerOverlay::collect() const
{
  DevicePolymerResults results{};
  if (terms == 0) return results;
  double termSums[termPartials] = {};
  for (std::size_t g = 0; g < termGroups; ++g)
  {
    const float* partial = hostPartials.data() + g * termPartials;
    for (std::size_t q = 0; q < termPartials; ++q) termSums[q] += static_cast<double>(partial[q]);
  }
  double atomSums[atomPartials] = {};
  const float* atomBase = hostPartials.data() + termGroups * termPartials;
  for (std::size_t g = 0; g < atomGroups; ++g)
  {
    const float* partial = atomBase + g * atomPartials;
    for (std::size_t q = 0; q < atomPartials; ++q) atomSums[q] += static_cast<double>(partial[q]);
  }
  auto tensor = [](const double* values) -> double3x3
  {
    double3x3 t{};
    t.ax = values[0];
    t.ay = values[1];
    t.az = values[2];
    t.bx = values[3];
    t.by = values[4];
    t.bz = values[5];
    t.cx = values[6];
    t.cy = values[7];
    t.cz = values[8];
    return t;
  };
  results.crossLink = termSums[0];
  results.bond = termSums[1];
  results.bend = termSums[2];
  results.torsion = termSums[3];
  results.moleculeMoleculeVDW = termSums[4];
  results.moleculeMoleculeCharge = termSums[5];
  results.ewaldExclusion = termSums[6];
  results.intraVDW = termSums[7];
  results.strain = tensor(termSums + 8);
  results.correction = tensor(atomSums);
  return results;
}

std::string DevicePolymerOverlay::status() const
{
  return std::format(
      "    polymer overlay on the device: {} terms over {} atoms, {} ramps (single precision, wrapped slot "
      "positions imaged along the terms)\n",
      terms, atoms, ramps);
}
