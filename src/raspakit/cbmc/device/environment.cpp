module;

module cbmc_device_environment;

import std;

import atom;
import double3x3;
import forcefield;
import vdwparameters;
import simulationbox;
import units;
import spatial_decomposition_settings;
import spatial_decomposition_device_context;
import spatial_decomposition_device_step;
import cbmc_device_kernels;

namespace
{
constexpr std::uint32_t frameworkBit = 0x80000000u;
// Two distinct sentinels: the molecule id of framework atoms (never excluded) and the "skip no molecule" id of
// the kernel parameters; they must differ, or the skip test would drop the framework.
constexpr std::uint32_t frameworkMolecule = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint32_t noMolecule = std::numeric_limits<std::uint32_t>::max() - 1;

void matrixToFloats(const double3x3& m, float* out)
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

std::uint32_t isOrthorhombic(const double3x3& cell)
{
  return (cell.ay == 0.0 && cell.az == 0.0 && cell.bx == 0.0 && cell.bz == 0.0 && cell.cx == 0.0 && cell.cy == 0.0)
             ? 1u
             : 0u;
}

std::size_t grown(std::size_t needed) { return needed + needed / 4 + 16; }
}  // namespace

bool CBMC::DeviceGrowEnvironment::available(PairDevice device) { return deviceAvailable(device); }

std::string CBMC::DeviceGrowEnvironment::unavailableReason(PairDevice device)
{
  return deviceUnavailableReason(device);
}

bool CBMC::DeviceGrowEnvironment::supportsForceField(const ForceField& forceField,
                                                     std::span<const std::uint8_t> presentTypes, std::string& reason)
{
  const std::size_t types = forceField.pseudoAtoms.size();
  for (std::size_t a = 0; a < types; ++a)
  {
    if (a >= presentTypes.size() || !presentTypes[a]) continue;
    for (std::size_t b = 0; b < types; ++b)
    {
      if (b >= presentTypes.size() || !presentTypes[b]) continue;
      const VDWParameters& parameters = forceField(a, b);
      if (parameters.type != VDWParameters::Type::LennardJones && parameters.type != VDWParameters::Type::None)
      {
        reason = std::format("the pair potential of '{}'-'{}' is '{}', the device grower evaluates plain 12-6 "
                             "Lennard-Jones only",
                             forceField.pseudoAtoms[a].name, forceField.pseudoAtoms[b].name,
                             VDWParameters::nameOfType(parameters.type));
        return false;
      }
    }
  }
  if (forceField.useCharge && forceField.chargeMethod != ForceField::ChargeMethod::Ewald)
  {
    reason = "the device grower evaluates Ewald electrostatics only";
    return false;
  }
  if (forceField.computePolarization)
  {
    reason = "the device grower does not evaluate polarization";
    return false;
  }
  return true;
}

CBMC::DeviceGrowEnvironment::DeviceGrowEnvironment(PairDevice device) : kind(device)
{
  context = createDeviceContext(device);
  kernel = context->compileKernel("cbmc_grow", cbmcDeviceGrowKernelSource, DeviceMath::Relaxed, "growEnergies");
  parameterBuffer.allocate(*context, sizeof(DeviceGrowParameters), DeviceMemory::Device);
}

CBMC::DeviceGrowEnvironment::~DeviceGrowEnvironment()
{
  if (context)
  {
    context->finish();
    parameterBuffer.reset();
    lennardJonesBuffer.reset();
    backgroundBuffer.reset();
    backgroundTypeBuffer.reset();
    backgroundMoleculeBuffer.reset();
    trialBuffer.reset();
    trialTypeBuffer.reset();
    partialBuffer.reset();
  }
}

std::string CBMC::DeviceGrowEnvironment::deviceName() const { return context ? context->deviceName() : ""; }

void CBMC::DeviceGrowEnvironment::setForceField(const ForceField& forceField, const SimulationBox& box)
{
  const std::size_t types = forceField.pseudoAtoms.size();
  if (uploadedForceField != &forceField || uploadedTypes != types)
  {
    // 4 epsilon, sigma^2, shift per type pair; a 'None' pair is all zeros (u = 0)
    std::vector<std::array<float, 4>> table(types * types, {0.0f, 0.0f, 0.0f, 0.0f});
    for (std::size_t a = 0; a < types; ++a)
    {
      for (std::size_t b = 0; b < types; ++b)
      {
        const VDWParameters& pair = forceField(a, b);
        if (pair.type != VDWParameters::Type::LennardJones) continue;
        const double sigma = pair.parameters.y;
        table[a * types + b] = {static_cast<float>(4.0 * pair.parameters.x), static_cast<float>(sigma * sigma),
                                static_cast<float>(pair.shift), 0.0f};
      }
    }
    lennardJonesBuffer.allocate(*context, table.size() * sizeof(std::array<float, 4>), DeviceMemory::Device);
    context->write(lennardJonesBuffer.get(), 0, table.size() * sizeof(std::array<float, 4>), table.data(), true);
    uploadedForceField = &forceField;
    uploadedTypes = types;
  }

  matrixToFloats(box.cell, parameters.cell);
  matrixToFloats(box.inverseCell, parameters.inverseCell);
  parameters.orthorhombic = isOrthorhombic(box.cell);
  parameters.cutOffFrameworkVDWSquared =
      static_cast<float>(forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW);
  parameters.cutOffMoleculeVDWSquared = static_cast<float>(forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW);
  parameters.cutOffChargeSquared = static_cast<float>(forceField.cutOffCoulomb * forceField.cutOffCoulomb);
  parameters.alpha = static_cast<float>(forceField.EwaldAlpha);
  parameters.coulombFactor = static_cast<float>(Units::CoulombicConversionFactor);
  parameters.overlapCriteria = static_cast<float>(forceField.energyOverlapCriteria);
  parameters.useCharge = forceField.useCharge ? 1u : 0u;
  parameters.numberOfTypes = static_cast<std::uint32_t>(types);
}

void CBMC::DeviceGrowEnvironment::ensureBackgroundCapacity(std::size_t atoms)
{
  if (atoms <= backgroundCapacity) return;
  backgroundCapacity = grown(atoms);
  backgroundBuffer.allocate(*context, backgroundCapacity * sizeof(std::array<float, 4>), DeviceMemory::Device);
  backgroundTypeBuffer.allocate(*context, backgroundCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
  backgroundMoleculeBuffer.allocate(*context, backgroundCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
  // a reallocation loses the cached framework
  uploadedFrameworkSize = std::numeric_limits<std::size_t>::max();
}

void CBMC::DeviceGrowEnvironment::setBackground(std::span<const Atom> frameworkAtoms,
                                                std::span<const Atom> moleculeAtoms)
{
  const std::size_t total = frameworkAtoms.size() + moleculeAtoms.size();
  ensureBackgroundCapacity(total);

  std::array<float, 9> cell{};
  std::copy(std::begin(parameters.cell), std::end(parameters.cell), cell.begin());
  const bool frameworkChanged = uploadedFrameworkData != frameworkAtoms.data() ||
                                uploadedFrameworkSize != frameworkAtoms.size() || uploadedCell != cell;

  auto pack = [&](std::span<const Atom> atoms, std::uint32_t typeBits, bool framework)
  {
    hostPositions.resize(atoms.size());
    hostTypes.resize(atoms.size());
    hostMolecules.resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      const Atom& atom = atoms[i];
      hostPositions[i] = {static_cast<float>(atom.position.x), static_cast<float>(atom.position.y),
                          static_cast<float>(atom.position.z), static_cast<float>(atom.charge * atom.scalingCoulomb)};
      hostTypes[i] = static_cast<std::uint32_t>(atom.type) | typeBits;
      hostMolecules[i] = framework ? frameworkMolecule : atom.moleculeId;
    }
  };
  auto upload = [&](std::size_t offset, std::size_t count)
  {
    if (count == 0) return;
    context->write(backgroundBuffer.get(), offset * sizeof(std::array<float, 4>), count * sizeof(std::array<float, 4>),
                   hostPositions.data(), true);
    context->write(backgroundTypeBuffer.get(), offset * sizeof(std::uint32_t), count * sizeof(std::uint32_t),
                   hostTypes.data(), true);
    context->write(backgroundMoleculeBuffer.get(), offset * sizeof(std::uint32_t), count * sizeof(std::uint32_t),
                   hostMolecules.data(), true);
  };

  if (frameworkChanged)
  {
    pack(frameworkAtoms, frameworkBit, true);
    upload(0, frameworkAtoms.size());
    uploadedFrameworkData = frameworkAtoms.data();
    uploadedFrameworkSize = frameworkAtoms.size();
    uploadedCell = cell;
  }
  numberOfFramework = frameworkAtoms.size();
  pack(moleculeAtoms, 0u, false);
  upload(numberOfFramework, moleculeAtoms.size());
  parameters.numberOfBackground = static_cast<std::uint32_t>(total);
}

void CBMC::DeviceGrowEnvironment::evaluate(std::span<const Atom> trialAtoms, std::size_t atomsPerSet,
                                           std::uint32_t trialMoleculeId, std::optional<std::size_t> skipMolecule,
                                           std::vector<DeviceTrialEnergy>& energies)
{
  energies.clear();
  if (atomsPerSet == 0 || atomsPerSet > maximumAtomsPerSet || trialAtoms.size() % atomsPerSet != 0)
  {
    throw std::invalid_argument(std::format(
        "CBMC device grower: a trial set of {} atoms is not supported (1..{} atoms, batch a multiple of the set)",
        atomsPerSet, maximumAtomsPerSet));
  }
  const std::size_t sets = trialAtoms.size() / atomsPerSet;
  if (sets == 0) return;

  if (trialAtoms.size() > trialCapacity)
  {
    trialCapacity = grown(trialAtoms.size());
    trialBuffer.allocate(*context, trialCapacity * sizeof(std::array<float, 4>), DeviceMemory::Device);
    trialTypeBuffer.allocate(*context, trialCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
  }
  if (sets * partials > partialCapacity)
  {
    partialCapacity = grown(sets * partials);
    partialBuffer.allocate(*context, partialCapacity * sizeof(float), DeviceMemory::Device);
  }

  hostTrial.resize(trialAtoms.size());
  hostTrialTypes.resize(trialAtoms.size());
  for (std::size_t i = 0; i < trialAtoms.size(); ++i)
  {
    const Atom& atom = trialAtoms[i];
    hostTrial[i] = {static_cast<float>(atom.position.x), static_cast<float>(atom.position.y),
                    static_cast<float>(atom.position.z), static_cast<float>(atom.charge * atom.scalingCoulomb)};
    hostTrialTypes[i] = static_cast<std::uint32_t>(atom.type);
  }
  context->write(trialBuffer.get(), 0, hostTrial.size() * sizeof(std::array<float, 4>), hostTrial.data(), false);
  context->write(trialTypeBuffer.get(), 0, hostTrialTypes.size() * sizeof(std::uint32_t), hostTrialTypes.data(),
                 false);

  parameters.atomsPerSet = static_cast<std::uint32_t>(atomsPerSet);
  parameters.trialMoleculeId = trialMoleculeId;
  parameters.skipMoleculeId = skipMolecule.has_value() ? static_cast<std::uint32_t>(skipMolecule.value()) : noMolecule;
  context->write(parameterBuffer.get(), 0, sizeof(DeviceGrowParameters), &parameters, false);

  const DeviceArg arguments[] = {DeviceArg::of(backgroundBuffer.get()),  DeviceArg::of(backgroundTypeBuffer.get()),
                                 DeviceArg::of(backgroundMoleculeBuffer.get()), DeviceArg::of(trialBuffer.get()),
                                 DeviceArg::of(trialTypeBuffer.get()),   DeviceArg::of(lennardJonesBuffer.get()),
                                 DeviceArg::of(parameterBuffer.get()),   DeviceArg::of(partialBuffer.get())};
  context->launch(kernel, arguments, sets, groupSize);

  hostPartials.resize(sets * partials);
  context->read(partialBuffer.get(), 0, sets * partials * sizeof(float), hostPartials.data());
  context->finish();

  energies.resize(sets);
  for (std::size_t s = 0; s < sets; ++s)
  {
    const float* partial = hostPartials.data() + s * partials;
    energies[s] = DeviceTrialEnergy{.vdw = static_cast<double>(partial[0]),
                                    .charge = static_cast<double>(partial[1]),
                                    .overlap = partial[2] > 0.5f};
  }
}
