module;

module spatial_decomposition_device_resident;

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import component;
import simulationbox;
import system;
import spatial_decomposition_cell_list;
import spatial_decomposition_device_context;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_step;

namespace
{
constexpr std::size_t atomGroup = 256;      // ATOM_GROUP of the kernel source
constexpr std::size_t moleculeGroup = 64;   // MOLECULE_GROUP of the kernel source

/// A double-float: the kernels' float2 (hi, lo); passed by value to the kernels.
struct Float2
{
  float x{0.0f};
  float y{0.0f};
};
static_assert(sizeof(Float2) == 8);

Float2 split(double value)
{
  const float hi = static_cast<float>(value);
  return Float2{hi, static_cast<float>(value - static_cast<double>(hi))};
}

double join(float hi, float lo) { return static_cast<double>(hi) + static_cast<double>(lo); }

std::size_t groupsFor(std::size_t count, std::size_t group) { return std::max<std::size_t>(1, (count + group - 1) / group); }

void put3(std::vector<float>& hi, std::vector<float>& lo, std::size_t index, const double3& value, float w = 0.0f)
{
  const Float2 x = split(value.x), y = split(value.y), z = split(value.z);
  hi[4 * index] = x.x;
  hi[4 * index + 1] = y.x;
  hi[4 * index + 2] = z.x;
  hi[4 * index + 3] = w;
  lo[4 * index] = x.y;
  lo[4 * index + 1] = y.y;
  lo[4 * index + 2] = z.y;
  lo[4 * index + 3] = 0.0f;
}

void put4(std::vector<float>& hi, std::vector<float>& lo, std::size_t index, const simd_quatd& value)
{
  const Float2 x = split(value.ix), y = split(value.iy), z = split(value.iz), r = split(value.r);
  hi[4 * index] = x.x;
  hi[4 * index + 1] = y.x;
  hi[4 * index + 2] = z.x;
  hi[4 * index + 3] = r.x;
  lo[4 * index] = x.y;
  lo[4 * index + 1] = y.y;
  lo[4 * index + 2] = z.y;
  lo[4 * index + 3] = r.y;
}

double3 get3(const std::vector<float>& hi, const std::vector<float>& lo, std::size_t index)
{
  return double3(join(hi[4 * index], lo[4 * index]), join(hi[4 * index + 1], lo[4 * index + 1]),
                 join(hi[4 * index + 2], lo[4 * index + 2]));
}

simd_quatd get4(const std::vector<float>& hi, const std::vector<float>& lo, std::size_t index)
{
  return simd_quatd(join(hi[4 * index], lo[4 * index]), join(hi[4 * index + 1], lo[4 * index + 1]),
                    join(hi[4 * index + 2], lo[4 * index + 2]), join(hi[4 * index + 3], lo[4 * index + 3]));
}
}  // namespace

bool DeviceResident::supports(const System& system, std::string& reason)
{
  if (system.thermobarostat.has_value())
  {
    reason = "a barostat (NPT ensembles)";
    return false;
  }
  for (const Component& component : system.components)
  {
    if (component.isSemiFlexible())
    {
      reason = std::format("semi-flexible molecules (component '{}')", component.name);
      return false;
    }
  }
  if (!system.spanOfGroupData().empty())
  {
    reason = "rigid-group state";
    return false;
  }
  return true;
}

void DeviceResident::initialize(DeviceContext& c)
{
  context = &c;
  const std::string_view program = "resident";
  atomsA = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentAtomsA");
  moleculesA = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentMoleculesA");
  pack = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentPack");
  torques = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentTorques");
  atomsB = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentAtomsB");
  moleculesB = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentMoleculesB");
  scaleAtoms = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentScaleAtoms");
  scaleMolecules = c.compileKernel(program, deviceKernelResidentSource, DeviceMath::Strict, "residentScaleMolecules");
  atoms = molecules = components = references = 0;
}

void DeviceResident::ensureBuffers(const System& system)
{
  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t numberOfComponents = system.components.size();
  std::size_t numberOfReferences = 0;
  for (const Component& component : system.components) numberOfReferences += component.atoms.size();
  if (numberOfAtoms == atoms && numberOfMolecules == molecules && numberOfComponents == components &&
      numberOfReferences == references)
  {
    return;
  }
  atoms = numberOfAtoms;
  molecules = numberOfMolecules;
  components = numberOfComponents;
  references = numberOfReferences;
  atomGroups = groupsFor(atoms, atomGroup);
  moleculeGroups = groupsFor(molecules, moleculeGroup);

  DeviceContext& c = *context;
  const std::size_t atomBytes = atoms * 4 * sizeof(float);
  const std::size_t moleculeBytes = molecules * 4 * sizeof(float);
  positionHi.allocate(c, atomBytes, DeviceMemory::Device);
  positionLo.allocate(c, atomBytes, DeviceMemory::Device);
  for (std::size_t k = 0; k < 2; ++k)
  {
    velocityHi[k].allocate(c, atomBytes, DeviceMemory::Device);
    velocityLo[k].allocate(c, atomBytes, DeviceMemory::Device);
    moleculeVelocityHi[k].allocate(c, moleculeBytes, DeviceMemory::Device);
    moleculeVelocityLo[k].allocate(c, moleculeBytes, DeviceMemory::Device);
    momentumHi[k].allocate(c, moleculeBytes, DeviceMemory::Device);
    momentumLo[k].allocate(c, moleculeBytes, DeviceMemory::Device);
  }
  atomMass.allocate(c, atomBytes, DeviceMemory::Device);
  atomInfo.allocate(c, atoms * 4 * sizeof(std::uint32_t), DeviceMemory::Device);
  slotOfAtom.allocate(c, atoms * sizeof(std::uint32_t), DeviceMemory::Device);
  translationHi.allocate(c, atomBytes, DeviceMemory::Device);
  translationLo.allocate(c, atomBytes, DeviceMemory::Device);
  comHi.allocate(c, moleculeBytes, DeviceMemory::Device);
  comLo.allocate(c, moleculeBytes, DeviceMemory::Device);
  orientationHi.allocate(c, moleculeBytes, DeviceMemory::Device);
  orientationLo.allocate(c, moleculeBytes, DeviceMemory::Device);
  moleculeGradient.allocate(c, moleculeBytes, DeviceMemory::Device);
  moleculeTorque.allocate(c, moleculeBytes, DeviceMemory::Device);
  moleculeMass.allocate(c, moleculeBytes, DeviceMemory::Device);
  moleculeInfo.allocate(c, molecules * 4 * sizeof(std::uint32_t), DeviceMemory::Device);
  componentInertia.allocate(c, components * 16 * sizeof(float), DeviceMemory::Device);
  componentReferenceOffset.allocate(c, components * sizeof(std::uint32_t), DeviceMemory::Device);
  referenceHi.allocate(c, references * 4 * sizeof(float), DeviceMemory::Device);
  referenceLo.allocate(c, references * 4 * sizeof(float), DeviceMemory::Device);
  packPartials.allocate(c, atomGroups * 2 * sizeof(float), DeviceMemory::Device);
  kineticPartials.allocate(c, atomGroups * 2 * sizeof(float), DeviceMemory::Device);
  moleculePartials.allocate(c, moleculeGroups * 4 * sizeof(float), DeviceMemory::Device);
  packHost.assign(atomGroups * 2, 0.0f);
  kineticHost.assign(atomGroups * 2 + moleculeGroups * 4, 0.0f);
  current = 0;
}

void DeviceResident::writeFloat4(DeviceBufferOwner& buffer, std::span<const float> values)
{
  if (values.empty()) return;
  context->write(buffer.get(), 0, values.size() * sizeof(float), values.data(), true);
}

void DeviceResident::readFloat4(const DeviceBufferOwner& buffer, std::size_t count, std::vector<float>& into)
{
  into.assign(4 * count, 0.0f);
  if (count == 0) return;
  context->read(buffer.get(), 0, into.size() * sizeof(float), into.data());
}

void DeviceResident::upload(const System& system, DeviceStep& step, const CellList& cells)
{
  ensureBuffers(system);
  const std::span<const Atom> atomData = system.spanOfMoleculeAtoms();
  const std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  const std::span<const Molecule> moleculeData = system.moleculeData;
  const std::vector<Component>& componentData = system.components;

  // component tables: inertia and the body-fixed reference positions
  std::vector<float> inertia(16 * components, 0.0f);
  std::vector<std::uint32_t> offsets(components, 0);
  std::vector<float> refHi(4 * references, 0.0f), refLo(4 * references, 0.0f);
  std::size_t offset = 0;
  for (std::size_t c = 0; c < components; ++c)
  {
    const Component& component = componentData[c];
    const double3 I = component.inertiaVector;
    const double3 invI = component.inverseInertiaVector;
    const Float2 ix = split(I.x), iy = split(I.y), iz = split(I.z);
    const Float2 jx = split(invI.x), jy = split(invI.y), jz = split(invI.z);
    float* table = inertia.data() + 16 * c;
    table[0] = ix.x;
    table[1] = iy.x;
    table[2] = iz.x;
    table[4] = ix.y;
    table[5] = iy.y;
    table[6] = iz.y;
    table[8] = jx.x;
    table[9] = jy.x;
    table[10] = jz.x;
    table[12] = jx.y;
    table[13] = jy.y;
    table[14] = jz.y;
    offsets[c] = static_cast<std::uint32_t>(offset);
    for (std::size_t b = 0; b < component.atoms.size(); ++b) put3(refHi, refLo, offset + b, component.atoms[b].position);
    offset += component.atoms.size();
  }
  writeFloat4(componentInertia, inertia);
  if (!offsets.empty())
  {
    context->write(componentReferenceOffset.get(), 0, offsets.size() * sizeof(std::uint32_t), offsets.data(), true);
  }
  writeFloat4(referenceHi, refHi);
  writeFloat4(referenceLo, refLo);

  // molecule records and the per-molecule gradients and torques of the current atom gradients (as
  // Integrators::updateCenterOfMassAndQuaternionGradients computes them)
  std::vector<float> hi(4 * molecules, 0.0f), lo(4 * molecules, 0.0f);
  std::vector<float> hi2(4 * molecules, 0.0f), lo2(4 * molecules, 0.0f);
  std::vector<std::uint32_t> info(4 * molecules, 0);
  std::vector<float> masses(4 * molecules, 0.0f);
  std::vector<float> gradients(4 * molecules, 0.0f), torqueValues(4 * molecules, 0.0f);
  std::vector<std::uint32_t> atomInfoValues(4 * atoms, 0);
  std::vector<float> atomMasses(4 * atoms, 0.0f);
  for (std::size_t m = 0; m < molecules; ++m)
  {
    const Molecule& molecule = moleculeData[m];
    const Component& component = componentData[molecule.componentId];
    info[4 * m] = static_cast<std::uint32_t>(molecule.atomIndex);
    info[4 * m + 1] = static_cast<std::uint32_t>(molecule.numberOfAtoms);
    info[4 * m + 2] = static_cast<std::uint32_t>(molecule.componentId);
    info[4 * m + 3] = component.rigid ? 1u : 0u;
    const Float2 mass = split(molecule.mass), invMass = split(molecule.invMass);
    masses[4 * m] = mass.x;
    masses[4 * m + 1] = mass.y;
    masses[4 * m + 2] = invMass.x;
    masses[4 * m + 3] = invMass.y;
    put3(hi, lo, m, molecule.centerOfMassPosition);
    put4(hi2, lo2, m, molecule.orientation);

    double3 gradient{};
    for (std::size_t b = 0; b < molecule.numberOfAtoms; ++b) gradient += dynamics[molecule.atomIndex + b].gradient;
    gradients[4 * m] = static_cast<float>(gradient.x);
    gradients[4 * m + 1] = static_cast<float>(gradient.y);
    gradients[4 * m + 2] = static_cast<float>(gradient.z);
    if (component.rigid)
    {
      const simd_quatd q = molecule.orientation;
      const double3x3 M = double3x3::buildRotationMatrix(q);
      double3 torque{};
      for (std::size_t b = 0; b < molecule.numberOfAtoms; ++b)
      {
        const double atomMassValue = component.definedAtoms[b].second;
        const double3 F = M * (dynamics[molecule.atomIndex + b].gradient - gradient * atomMassValue * molecule.invMass);
        torque += double3::cross(F, component.atoms[b].position);
      }
      const simd_quatd orientationGradient = -2.0 * q * simd_quatd(0.0, torque);
      torqueValues[4 * m] = static_cast<float>(orientationGradient.ix);
      torqueValues[4 * m + 1] = static_cast<float>(orientationGradient.iy);
      torqueValues[4 * m + 2] = static_cast<float>(orientationGradient.iz);
      torqueValues[4 * m + 3] = static_cast<float>(orientationGradient.r);
    }
    for (std::size_t b = 0; b < molecule.numberOfAtoms; ++b)
    {
      const std::size_t i = molecule.atomIndex + b;
      atomInfoValues[4 * i] = static_cast<std::uint32_t>(m);
      atomInfoValues[4 * i + 1] = static_cast<std::uint32_t>(b);
      atomInfoValues[4 * i + 2] = static_cast<std::uint32_t>(molecule.atomIndex);
      atomInfoValues[4 * i + 3] = component.rigid ? 0u : 1u;
      const double massValue = component.definedAtoms[b].second;
      const Float2 am = split(massValue), aim = split(1.0 / massValue);
      atomMasses[4 * i] = am.x;
      atomMasses[4 * i + 1] = am.y;
      atomMasses[4 * i + 2] = aim.x;
      atomMasses[4 * i + 3] = aim.y;
    }
  }
  writeFloat4(comHi, hi);
  writeFloat4(comLo, lo);
  writeFloat4(orientationHi, hi2);
  writeFloat4(orientationLo, lo2);
  writeFloat4(moleculeMass, masses);
  writeFloat4(moleculeGradient, gradients);
  writeFloat4(moleculeTorque, torqueValues);
  if (molecules > 0)
  {
    context->write(moleculeInfo.get(), 0, info.size() * sizeof(std::uint32_t), info.data(), true);
  }
  // velocities and momenta into both buffer sets (current and pending)
  for (std::size_t m = 0; m < molecules; ++m)
  {
    put3(hi, lo, m, moleculeData[m].velocity);
    put4(hi2, lo2, m, moleculeData[m].orientationMomentum);
  }
  for (std::size_t k = 0; k < 2; ++k)
  {
    writeFloat4(moleculeVelocityHi[k], hi);
    writeFloat4(moleculeVelocityLo[k], lo);
    writeFloat4(momentumHi[k], hi2);
    writeFloat4(momentumLo[k], lo2);
  }

  // atoms
  if (atoms > 0)
  {
    context->write(atomInfo.get(), 0, atomInfoValues.size() * sizeof(std::uint32_t), atomInfoValues.data(), true);
  }
  writeFloat4(atomMass, atomMasses);
  hi.assign(4 * atoms, 0.0f);
  lo.assign(4 * atoms, 0.0f);
  for (std::size_t i = 0; i < atoms; ++i) put3(hi, lo, i, atomData[i].position, static_cast<float>(atomData[i].charge));
  writeFloat4(positionHi, hi);
  writeFloat4(positionLo, lo);
  hi.assign(4 * atoms, 0.0f);
  lo.assign(4 * atoms, 0.0f);
  for (std::size_t i = 0; i < atoms; ++i)
  {
    if (atomInfoValues[4 * i + 3] != 0u) put3(hi, lo, i, dynamics[i].velocity);
  }
  for (std::size_t k = 0; k < 2; ++k)
  {
    writeFloat4(velocityHi[k], hi);
    writeFloat4(velocityLo[k], lo);
  }
  current = 0;

  // the gradients into the force buffer of the step (slot order; the dummy slots zero)
  setLayout(step, cells, system.simulationBox);
  const std::span<const std::uint32_t> slots = step.slotsOfSorted();
  stage.assign(4 * step.numberOfSlots(), 0.0f);
  for (std::size_t i = 0; i < atoms; ++i)
  {
    const std::size_t slot = slots[cells.originalToSorted[i]];
    const double3 g = dynamics[i].gradient;
    stage[4 * slot] = static_cast<float>(g.x);
    stage[4 * slot + 1] = static_cast<float>(g.y);
    stage[4 * slot + 2] = static_cast<float>(g.z);
  }
  if (!stage.empty()) context->write(targets.forces, 0, stage.size() * sizeof(float), stage.data(), true);
}

void DeviceResident::setLayout(DeviceStep& step, const CellList& cells, const SimulationBox& box)
{
  targets = step.residentTargets();
  const std::span<const std::uint32_t> slots = step.slotsOfSorted();
  stageU.assign(atoms, 0);
  std::vector<float> hi(4 * atoms, 0.0f), lo(4 * atoms, 0.0f);
  for (std::size_t i = 0; i < atoms; ++i)
  {
    const std::uint32_t sorted = cells.originalToSorted[i];
    stageU[i] = slots[sorted];
    const int3 w = cells.wrap[sorted];
    const double3 translation =
        box.cell * double3(static_cast<double>(w.x), static_cast<double>(w.y), static_cast<double>(w.z));
    put3(hi, lo, i, translation);
  }
  if (atoms > 0) context->write(slotOfAtom.get(), 0, stageU.size() * sizeof(std::uint32_t), stageU.data(), true);
  writeFloat4(translationHi, hi);
  writeFloat4(translationLo, lo);
}

DeviceEvent DeviceResident::enqueuePack()
{
  const std::uint32_t count = static_cast<std::uint32_t>(atoms);
  const std::uint32_t writeRelative = targets.relativeEnabled ? 1u : 0u;
  const DeviceArg arguments[] = {DeviceArg::of(positionHi.get()),
                                 DeviceArg::of(positionLo.get()),
                                 DeviceArg::of(translationHi.get()),
                                 DeviceArg::of(translationLo.get()),
                                 DeviceArg::of(slotOfAtom.get()),
                                 DeviceArg::of(atomInfo.get()),
                                 DeviceArg::of(targets.buildPositions),
                                 DeviceArg::of(targets.compactReference),
                                 DeviceArg::of(targets.positions),
                                 DeviceArg::of(targets.relative),
                                 DeviceArg::of(packPartials.get()),
                                 DeviceArg::value(count),
                                 DeviceArg::value(writeRelative)};
  context->launch(pack, arguments, atomGroups, atomGroup);
  context->read(packPartials.get(), 0, packHost.size() * sizeof(float), packHost.data());
  const DeviceEvent event = context->mark();
  context->flush();
  return event;
}

DeviceEvent DeviceResident::enqueueFirstHalf(const Scaling& scaling, double timeStep)
{
  const std::uint32_t atomCount = static_cast<std::uint32_t>(atoms);
  const std::uint32_t moleculeCount = static_cast<std::uint32_t>(molecules);
  const Float2 scaleT = split(scaling.translational);
  const Float2 scaleR = split(scaling.rotational);
  const Float2 halfDt = split(0.5 * timeStep);
  const Float2 dt = split(timeStep);
  const Float2 dtTenth = split(0.5 * timeStep / 5.0);
  const Float2 dtFifth = split(timeStep / 5.0);
  {
    const DeviceArg arguments[] = {DeviceArg::of(positionHi.get()),
                                   DeviceArg::of(positionLo.get()),
                                   DeviceArg::of(velocityHi[current].get()),
                                   DeviceArg::of(velocityLo[current].get()),
                                   DeviceArg::of(atomMass.get()),
                                   DeviceArg::of(atomInfo.get()),
                                   DeviceArg::of(slotOfAtom.get()),
                                   DeviceArg::of(targets.forces),
                                   DeviceArg::value(atomCount),
                                   DeviceArg::value(scaleT),
                                   DeviceArg::value(halfDt),
                                   DeviceArg::value(dt)};
    context->launch(atomsA, arguments, atomGroups, atomGroup);
  }
  {
    const DeviceArg arguments[] = {DeviceArg::of(positionHi.get()),
                                   DeviceArg::of(positionLo.get()),
                                   DeviceArg::of(atomMass.get()),
                                   DeviceArg::of(comHi.get()),
                                   DeviceArg::of(comLo.get()),
                                   DeviceArg::of(moleculeVelocityHi[current].get()),
                                   DeviceArg::of(moleculeVelocityLo[current].get()),
                                   DeviceArg::of(orientationHi.get()),
                                   DeviceArg::of(orientationLo.get()),
                                   DeviceArg::of(momentumHi[current].get()),
                                   DeviceArg::of(momentumLo[current].get()),
                                   DeviceArg::of(moleculeGradient.get()),
                                   DeviceArg::of(moleculeTorque.get()),
                                   DeviceArg::of(moleculeMass.get()),
                                   DeviceArg::of(moleculeInfo.get()),
                                   DeviceArg::of(componentInertia.get()),
                                   DeviceArg::of(componentReferenceOffset.get()),
                                   DeviceArg::of(referenceHi.get()),
                                   DeviceArg::of(referenceLo.get()),
                                   DeviceArg::value(moleculeCount),
                                   DeviceArg::value(scaleT),
                                   DeviceArg::value(scaleR),
                                   DeviceArg::value(halfDt),
                                   DeviceArg::value(dt),
                                   DeviceArg::value(dtTenth),
                                   DeviceArg::value(dtFifth)};
    context->launch(moleculesA, arguments, moleculeGroups, moleculeGroup);
  }
  return enqueuePack();
}

DeviceResident::Displacement DeviceResident::collectDisplacement() const
{
  Displacement result{};
  for (std::size_t g = 0; g < atomGroups; ++g)
  {
    result.sinceBuild = std::max(result.sinceBuild, static_cast<double>(packHost[2 * g]));
    result.sinceCompaction = std::max(result.sinceCompaction, static_cast<double>(packHost[2 * g + 1]));
  }
  return result;
}

DeviceEvent DeviceResident::enqueueSecondHalf(double timeStep)
{
  const std::uint32_t atomCount = static_cast<std::uint32_t>(atoms);
  const std::uint32_t moleculeCount = static_cast<std::uint32_t>(molecules);
  const Float2 halfDt = split(0.5 * timeStep);
  const std::size_t pending = 1 - current;
  {
    const DeviceArg arguments[] = {DeviceArg::of(targets.forces),
                                   DeviceArg::of(slotOfAtom.get()),
                                   DeviceArg::of(atomMass.get()),
                                   DeviceArg::of(orientationHi.get()),
                                   DeviceArg::of(moleculeMass.get()),
                                   DeviceArg::of(moleculeInfo.get()),
                                   DeviceArg::of(componentReferenceOffset.get()),
                                   DeviceArg::of(referenceHi.get()),
                                   DeviceArg::of(moleculeGradient.get()),
                                   DeviceArg::of(moleculeTorque.get()),
                                   DeviceArg::value(moleculeCount)};
    context->launch(torques, arguments, moleculeGroups, moleculeGroup);
  }
  {
    const DeviceArg arguments[] = {DeviceArg::of(velocityHi[current].get()),
                                   DeviceArg::of(velocityLo[current].get()),
                                   DeviceArg::of(velocityHi[pending].get()),
                                   DeviceArg::of(velocityLo[pending].get()),
                                   DeviceArg::of(atomMass.get()),
                                   DeviceArg::of(atomInfo.get()),
                                   DeviceArg::of(slotOfAtom.get()),
                                   DeviceArg::of(targets.forces),
                                   DeviceArg::of(kineticPartials.get()),
                                   DeviceArg::value(atomCount),
                                   DeviceArg::value(halfDt)};
    context->launch(atomsB, arguments, atomGroups, atomGroup);
  }
  {
    const DeviceArg arguments[] = {DeviceArg::of(moleculeVelocityHi[current].get()),
                                   DeviceArg::of(moleculeVelocityLo[current].get()),
                                   DeviceArg::of(moleculeVelocityHi[pending].get()),
                                   DeviceArg::of(moleculeVelocityLo[pending].get()),
                                   DeviceArg::of(momentumHi[current].get()),
                                   DeviceArg::of(momentumLo[current].get()),
                                   DeviceArg::of(momentumHi[pending].get()),
                                   DeviceArg::of(momentumLo[pending].get()),
                                   DeviceArg::of(orientationHi.get()),
                                   DeviceArg::of(orientationLo.get()),
                                   DeviceArg::of(moleculeGradient.get()),
                                   DeviceArg::of(moleculeTorque.get()),
                                   DeviceArg::of(moleculeMass.get()),
                                   DeviceArg::of(moleculeInfo.get()),
                                   DeviceArg::of(componentInertia.get()),
                                   DeviceArg::of(moleculePartials.get()),
                                   DeviceArg::value(moleculeCount),
                                   DeviceArg::value(halfDt)};
    context->launch(moleculesB, arguments, moleculeGroups, moleculeGroup);
  }
  context->read(kineticPartials.get(), 0, atomGroups * 2 * sizeof(float), kineticHost.data());
  context->read(moleculePartials.get(), 0, moleculeGroups * 4 * sizeof(float), kineticHost.data() + 2 * atomGroups);
  const DeviceEvent event = context->mark();
  context->flush();
  return event;
}

DeviceResident::Kinetic DeviceResident::collectKinetic() const
{
  Kinetic result{};
  for (std::size_t g = 0; g < atomGroups; ++g)
  {
    result.translational += join(kineticHost[2 * g], kineticHost[2 * g + 1]);
  }
  const float* moleculeSums = kineticHost.data() + 2 * atomGroups;
  for (std::size_t g = 0; g < moleculeGroups; ++g)
  {
    result.translational += join(moleculeSums[4 * g], moleculeSums[4 * g + 1]);
    result.rotational += join(moleculeSums[4 * g + 2], moleculeSums[4 * g + 3]);
  }
  return result;
}

void DeviceResident::acceptVelocities() { current = 1 - current; }

void DeviceResident::enqueueScale(const Scaling& scaling)
{
  if (scaling.translational == 1.0 && scaling.rotational == 1.0) return;
  const std::uint32_t atomCount = static_cast<std::uint32_t>(atoms);
  const std::uint32_t moleculeCount = static_cast<std::uint32_t>(molecules);
  const Float2 scaleT = split(scaling.translational);
  const Float2 scaleR = split(scaling.rotational);
  {
    const DeviceArg arguments[] = {DeviceArg::of(velocityHi[current].get()), DeviceArg::of(velocityLo[current].get()),
                                   DeviceArg::of(atomInfo.get()), DeviceArg::value(atomCount),
                                   DeviceArg::value(scaleT)};
    context->launch(scaleAtoms, arguments, atomGroups, atomGroup);
  }
  {
    const DeviceArg arguments[] = {DeviceArg::of(moleculeVelocityHi[current].get()),
                                   DeviceArg::of(moleculeVelocityLo[current].get()),
                                   DeviceArg::of(momentumHi[current].get()),
                                   DeviceArg::of(momentumLo[current].get()),
                                   DeviceArg::of(moleculeInfo.get()),
                                   DeviceArg::value(moleculeCount),
                                   DeviceArg::value(scaleT),
                                   DeviceArg::value(scaleR)};
    context->launch(scaleMolecules, arguments, moleculeGroups, moleculeGroup);
  }
}

void DeviceResident::downloadPositions(System& system)
{
  std::vector<float> hi, lo;
  readFloat4(positionHi, atoms, hi);
  readFloat4(positionLo, atoms, lo);
  context->finish();
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < atoms; ++i) atomData[i].position = get3(hi, lo, i);
}

void DeviceResident::download(System& system, const DeviceStep& step, const CellList& cells)
{
  std::vector<float> hi, lo, hi2, lo2, hi3, lo3, hi4, lo4, gradients, torqueValues, forces;
  readFloat4(positionHi, atoms, hi);
  readFloat4(positionLo, atoms, lo);
  readFloat4(velocityHi[current], atoms, hi2);
  readFloat4(velocityLo[current], atoms, lo2);
  readFloat4(comHi, molecules, hi3);
  readFloat4(comLo, molecules, lo3);
  readFloat4(moleculeVelocityHi[current], molecules, hi4);
  readFloat4(moleculeVelocityLo[current], molecules, lo4);
  readFloat4(moleculeGradient, molecules, gradients);
  readFloat4(moleculeTorque, molecules, torqueValues);
  std::vector<float> qHi, qLo, pHi, pLo;
  readFloat4(orientationHi, molecules, qHi);
  readFloat4(orientationLo, molecules, qLo);
  readFloat4(momentumHi[current], molecules, pHi);
  readFloat4(momentumLo[current], molecules, pLo);
  forces.assign(4 * step.numberOfSlots(), 0.0f);
  if (!forces.empty()) context->read(targets.forces, 0, forces.size() * sizeof(float), forces.data());
  context->finish();

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  const std::span<const std::uint32_t> slots = step.slotsOfSorted();
  for (std::size_t i = 0; i < atoms; ++i)
  {
    atomData[i].position = get3(hi, lo, i);
    const std::size_t slot = slots[cells.originalToSorted[i]];
    dynamics[i].gradient = double3(static_cast<double>(forces[4 * slot]), static_cast<double>(forces[4 * slot + 1]),
                                   static_cast<double>(forces[4 * slot + 2]));
  }
  for (std::size_t m = 0; m < molecules; ++m)
  {
    Molecule& molecule = system.moleculeData[m];
    const bool rigid = system.components[molecule.componentId].rigid;
    molecule.centerOfMassPosition = get3(hi3, lo3, m);
    molecule.velocity = get3(hi4, lo4, m);
    molecule.gradient = double3(static_cast<double>(gradients[4 * m]), static_cast<double>(gradients[4 * m + 1]),
                                static_cast<double>(gradients[4 * m + 2]));
    if (rigid)
    {
      molecule.orientation = get4(qHi, qLo, m);
      molecule.orientationMomentum = get4(pHi, pLo, m);
      molecule.orientationGradient =
          simd_quatd(static_cast<double>(torqueValues[4 * m]), static_cast<double>(torqueValues[4 * m + 1]),
                     static_cast<double>(torqueValues[4 * m + 2]), static_cast<double>(torqueValues[4 * m + 3]));
    }
    else
    {
      molecule.orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);
      for (std::size_t b = 0; b < molecule.numberOfAtoms; ++b)
      {
        const std::size_t i = molecule.atomIndex + b;
        dynamics[i].velocity = get3(hi2, lo2, i);
      }
    }
  }
}
