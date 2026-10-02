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
import atom;
import molecule;
import component;
import forcefield;
import system;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import units;
import intra_molecular_potentials;
import opencl;
import spatial_decomposition_opencl_handles;

using OpenCLDevice::check;
using OpenCLDevice::roundUp;

namespace
{
constexpr std::size_t groupSize = 64;
constexpr std::size_t partialsPerGroup = 26;
constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
constexpr std::size_t maximumAtomsPerMolecule = 256;

template <typename Buffer>
void ensureCapacity(Buffer& buffer, std::size_t& capacity, std::size_t required, std::size_t elementBytes,
                    cl_mem_flags flags)
{
  if (required <= capacity) return;
  capacity = required + required / 4 + 1;
  buffer.reset(OpenCL::createBuffer(flags, capacity * elementBytes));
}
}  // namespace

bool OpenCLBonded::supports(const System& system, std::string& reason)
{
  for (const Component& component : system.components)
  {
    const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
    const char* unsupported = nullptr;
    if (!potentials.ureyBradleys.empty())
      unsupported = "Urey-Bradley";
    else if (!potentials.inversionBends.empty())
      unsupported = "inversion-bend";
    else if (!potentials.outOfPlaneBends.empty())
      unsupported = "out-of-plane-bend";
    else if (!potentials.bondBonds.empty())
      unsupported = "bond-bond";
    else if (!potentials.bondBends.empty())
      unsupported = "bond-bend";
    else if (!potentials.bondTorsions.empty())
      unsupported = "bond-torsion";
    else if (!potentials.bendBends.empty())
      unsupported = "bend-bend";
    else if (!potentials.bendTorsions.empty())
      unsupported = "bend-torsion";
    if (unsupported)
    {
      reason = std::format("{} terms (component '{}')", unsupported, component.name);
      return false;
    }
    if (component.atoms.size() > maximumAtomsPerMolecule)
    {
      reason =
          std::format("molecules with more than {} atoms (component '{}')", maximumAtomsPerMolecule, component.name);
      return false;
    }
  }
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    if (atom.groupId != 0)
    {
      reason = "dU/dlambda group atoms";
      return false;
    }
  }
  return true;
}

void OpenCLBonded::initialize(cl_context clContext, cl_device_id clDevice)
{
  context = clContext;
  device = clDevice;
  program = OpenCLDevice::buildProgram(context, device, openclBondedKernelSource, "-cl-mad-enable -cl-no-signed-zeros",
                                       "OpenCL bonded");
  kernel = OpenCLDevice::createKernel(program.get(), "bondedAtoms");
  parameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(Parameters)));
}

void OpenCLBonded::setTopology(const System& system)
{
  terms.clear();
  atomTermStart.clear();
  atomTerms.clear();
  molecules.clear();
  masses.clear();

  // the terms of every component, flattened, with the references per component atom (CSR)
  std::vector<std::uint32_t> componentAtomOffset(system.components.size(), 0);
  std::vector<std::vector<std::uint32_t>> referencesPerAtom;
  for (std::size_t c = 0; c < system.components.size(); ++c)
  {
    const Component& component = system.components[c];
    const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
    componentAtomOffset[c] = static_cast<std::uint32_t>(referencesPerAtom.size());
    const std::size_t atomsInComponent = component.atoms.size();
    referencesPerAtom.resize(referencesPerAtom.size() + atomsInComponent);
    auto addTerm = [&](std::uint32_t kind, std::size_t type, std::span<const std::size_t> identifiers,
                       std::span<const double> values)
    {
      Term term{};
      term.kind = kind;
      term.type = static_cast<std::uint32_t>(type);
      for (std::size_t k = 0; k < identifiers.size(); ++k) term.atoms[k] = static_cast<std::uint32_t>(identifiers[k]);
      for (std::size_t k = 0; k < std::min<std::size_t>(values.size(), 6); ++k)
      {
        term.parameters[k] = static_cast<float>(values[k]);
      }
      const std::uint32_t index = static_cast<std::uint32_t>(terms.size());
      terms.push_back(term);
      for (std::uint32_t role = 0; role < identifiers.size(); ++role)
      {
        referencesPerAtom[componentAtomOffset[c] + identifiers[role]].push_back((index << 2) | role);
      }
    };
    for (const BondPotential& bond : potentials.bonds)
    {
      addTerm(0, std::to_underlying(bond.type), bond.identifiers, bond.parameters);
    }
    for (const BendPotential& bend : potentials.bends)
    {
      addTerm(1, std::to_underlying(bend.type), bend.identifiers, bend.parameters);
    }
    for (const TorsionPotential& torsion : potentials.torsions)
    {
      addTerm(2, std::to_underlying(torsion.type), torsion.identifiers, torsion.parameters);
    }
    for (const TorsionPotential& torsion : potentials.improperTorsions)
    {
      addTerm(3, std::to_underlying(torsion.type), torsion.identifiers, torsion.parameters);
    }
    for (const VanDerWaalsPotential& pair : potentials.vanDerWaals)
    {
      const double values[2] = {pair.scaling * 4.0 * pair.parameters[0], pair.parameters[1] * pair.parameters[1]};
      addTerm(4, std::to_underlying(pair.type), pair.identifiers, values);
    }
    for (const CoulombPotential& pair : potentials.coulombs)
    {
      const double values[1] = {pair.scaling * Units::CoulombicConversionFactor * pair.chargeA * pair.chargeB};
      addTerm(5, std::to_underlying(pair.type), pair.identifiers, values);
    }
  }
  atomTermStart.reserve(referencesPerAtom.size() + 1);
  atomTermStart.push_back(0);
  for (const std::vector<std::uint32_t>& references : referencesPerAtom)
  {
    atomTerms.insert(atomTerms.end(), references.begin(), references.end());
    atomTermStart.push_back(static_cast<std::uint32_t>(atomTerms.size()));
  }

  molecules.reserve(system.moleculeData.size());
  for (const Molecule& molecule : system.moleculeData)
  {
    molecules.push_back(MoleculeInfo{static_cast<std::uint32_t>(molecule.atomIndex),
                                     static_cast<std::uint32_t>(molecule.numberOfAtoms),
                                     componentAtomOffset[molecule.componentId], 0u});
  }
  masses.reserve(system.forceField.pseudoAtoms.size());
  for (const auto& pseudoAtom : system.forceField.pseudoAtoms) masses.push_back(static_cast<float>(pseudoAtom.mass));
  chargeSquaredSum = 0.0;
  for (const Atom& atom : system.spanOfMoleculeAtoms()) chargeSquaredSum += atom.charge * atom.charge;

  // blocking uploads (the shared queue; complete before the step's queue uses them)
  auto upload = [&](OpenCLDevice::MemHandle& buffer, const auto& values)
  {
    using Value = typename std::remove_cvref_t<decltype(values)>::value_type;
    const std::size_t bytes = std::max<std::size_t>(1, values.size()) * sizeof(Value);
    buffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, bytes));
    if (!values.empty()) OpenCL::writeBuffer(buffer.get(), values.size() * sizeof(Value), values.data());
  };
  upload(termsBuffer, terms);
  upload(atomTermStartBuffer, atomTermStart);
  upload(atomTermsBuffer, atomTerms);
  upload(moleculeInfoBuffer, molecules);
  upload(massBuffer, masses);
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

void OpenCLBonded::setLayout(cl_command_queue queue, std::span<const std::uint32_t> slotOfSorted,
                             std::span<const std::uint32_t> originalToSorted, std::size_t slots)
{
  const std::size_t numberOfAtoms = originalToSorted.size();
  slotMolecule.assign(slots, noAtom);
  slotOfOriginal.assign(numberOfAtoms, noAtom);
  referenceOfSorted.assign(numberOfAtoms, 0);
  for (std::size_t m = 0; m < molecules.size(); ++m)
  {
    const MoleculeInfo& molecule = molecules[m];
    const std::uint32_t reference = originalToSorted[molecule.firstAtom];
    for (std::uint32_t b = 0; b < molecule.numberOfAtoms; ++b)
    {
      const std::uint32_t sorted = originalToSorted[molecule.firstAtom + b];
      const std::uint32_t slot = slotOfSorted[sorted];
      slotMolecule[slot] = (static_cast<std::uint32_t>(m) << 8) | b;
      slotOfOriginal[molecule.firstAtom + b] = slot;
      referenceOfSorted[sorted] = reference;
    }
  }
  if (parameters.numberOfSlots != static_cast<std::uint32_t>(slots))
  {
    parameters.numberOfSlots = static_cast<std::uint32_t>(slots);
    parametersChanged = true;
  }

  ensureCapacity(slotMoleculeBuffer, slotCapacity, slots, sizeof(std::uint32_t), CL_MEM_READ_ONLY);
  ensureCapacity(slotOfOriginalBuffer, atomCapacity, numberOfAtoms, sizeof(std::uint32_t), CL_MEM_READ_ONLY);
  groups = roundUp(std::max<std::size_t>(slots, 1), groupSize) / groupSize;
  if (groups > partialCapacity)
  {
    partialCapacity = groups + groups / 4 + 1;
    partialBuffer.reset(OpenCL::createBuffer(CL_MEM_WRITE_ONLY, partialCapacity * partialsPerGroup * sizeof(float)));
  }
  hostPartials.assign(groups * partialsPerGroup, 0.0f);

  check(clEnqueueWriteBuffer(queue, slotMoleculeBuffer.get(), CL_FALSE, 0, slots * sizeof(std::uint32_t),
                             slotMolecule.data(), 0, nullptr, nullptr),
        "clEnqueueWriteBuffer (slot molecules)");
  if (numberOfAtoms > 0)
  {
    check(clEnqueueWriteBuffer(queue, slotOfOriginalBuffer.get(), CL_FALSE, 0, numberOfAtoms * sizeof(std::uint32_t),
                               slotOfOriginal.data(), 0, nullptr, nullptr),
          "clEnqueueWriteBuffer (slot of atom)");
  }
}

void OpenCLBonded::enqueue(cl_command_queue queue, cl_mem position, cl_mem relative, cl_mem typeOf, cl_mem force)
{
  if (parametersChanged)
  {
    parametersChanged = false;
    check(clEnqueueWriteBuffer(queue, parameterBuffer.get(), CL_FALSE, 0, sizeof(Parameters), &parameters, 0, nullptr,
                               nullptr),
          "clEnqueueWriteBuffer (bonded parameters)");
  }
  const cl_mem buffers[13] = {position,
                              relative,
                              typeOf,
                              slotMoleculeBuffer.get(),
                              moleculeInfoBuffer.get(),
                              slotOfOriginalBuffer.get(),
                              atomTermStartBuffer.get(),
                              atomTermsBuffer.get(),
                              termsBuffer.get(),
                              massBuffer.get(),
                              parameterBuffer.get(),
                              force,
                              partialBuffer.get()};
  OpenCLDevice::setBufferArguments(kernel.get(), buffers, "clSetKernelArg (bondedAtoms)");
  const std::size_t global = groups * groupSize;
  const std::size_t local = groupSize;
  check(clEnqueueNDRangeKernel(queue, kernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (bondedAtoms)");
}

void OpenCLBonded::enqueueRead(cl_command_queue queue, cl_event* event)
{
  check(clEnqueueReadBuffer(queue, partialBuffer.get(), CL_FALSE, 0, hostPartials.size() * sizeof(float),
                            hostPartials.data(), 0, nullptr, event),
        "clEnqueueReadBuffer (bonded partials)");
}

OpenCLBonded::Results OpenCLBonded::collect() const
{
  double sums[partialsPerGroup] = {};
  for (std::size_t g = 0; g < groups; ++g)
  {
    const float* partial = hostPartials.data() + g * partialsPerGroup;
    for (std::size_t q = 0; q < partialsPerGroup; ++q) sums[q] += static_cast<double>(partial[q]);
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
  Results results{};
  if (parameters.useCharge != 0)
  {
    results.self = -conversionFactor * alphaValue * std::numbers::inv_sqrtpi * chargeSquaredSum;
    results.exclusion = (sums[0] + sums[1]) - results.self;
  }
  results.bond = sums[2];
  results.bend = sums[3];
  results.torsion = sums[4];
  results.improperTorsion = sums[5];
  results.exclusionStrain = tensor(6);
  results.correction = tensor(15);
  results.intraVDW = sums[24];
  results.intraCoulomb = sums[25];
  return results;
}

std::string OpenCLBonded::status() const
{
  return std::format(
      "    bonded terms on the device: {} molecules, {} bonded and intramolecular pair terms over the components, "
      "self and exclusion corrections, virial correction (single precision, positions relative to the molecule)\n",
      molecules.size(), terms.size());
}
