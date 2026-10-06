module;

module spatial_decomposition_force_engine;

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import molecule;
import component;
import vdwparameters;
import forcefield;
import simulationbox;
import running_energy;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_intra_pair;
import intra_molecular_exclusions;
import interactions_ewald;
import integrators_update;
import integrators_compute;
import thermostat;
import thermobarostat;
import system;
import spatial_decomposition_settings;
import spatial_decomposition_cell_list;
import spatial_decomposition_pppm;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_cluster_kernel;
import spatial_decomposition_device_context;
import spatial_decomposition_device_step;
import spatial_decomposition_device_resident;
import spatial_decomposition_worker_team;

SpatialDecompositionForceEngine::SpatialDecompositionForceEngine(const SpatialDecompositionSettings& s) : settings(s)
{
  settings.numberOfThreads = std::max<std::size_t>(1, settings.numberOfThreads);
}

SpatialDecompositionForceEngine::~SpatialDecompositionForceEngine() = default;
SpatialDecompositionForceEngine::SpatialDecompositionForceEngine(SpatialDecompositionForceEngine&&) noexcept = default;
SpatialDecompositionForceEngine& SpatialDecompositionForceEngine::operator=(
    SpatialDecompositionForceEngine&&) noexcept = default;

bool SpatialDecompositionForceEngine::supports(const System& system, std::string& reason)
{
  if (system.framework.has_value() || system.numberOfFrameworkAtoms > 0)
  {
    reason = "a framework";
    return false;
  }
  if (system.hasExternalField)
  {
    reason = "an external field";
    return false;
  }
  if (system.forceField.computePolarization)
  {
    reason = "polarization";
    return false;
  }
  if (!system.crossLinks.empty())
  {
    reason = "cross-link bonds";
    return false;
  }
  if (molecularDynamicsHasParticleExchange(system.molecularDynamicsEnsemble))
  {
    reason = std::format("the '{}' ensemble (particle exchange during the MD stages)",
                         molecularDynamicsEnsembleName(system.molecularDynamicsEnsemble));
    return false;
  }
  for (const Component& component : system.components)
  {
    if (component.hasFractionalMolecule)
    {
      reason = std::format("fractional (CFCMC) molecules (component '{}')", component.name);
      return false;
    }
  }
  if (system.forceField.omitInterInteractions)
  {
    reason = "'OmitInterInteractions'";
    return false;
  }
  if (system.forceField.settings.useDualCutOff)
  {
    reason = "'UseDualCutOff'";
    return false;
  }
  return true;
}

void SpatialDecompositionForceEngine::prepareBondedWork(const System& system)
{
  const std::size_t threads = settings.numberOfThreads;
  const std::size_t numberOfMolecules = system.moleculeData.size();
  // about eight chunks per thread: fine enough to balance the threads that are free next to the FFTs of thread 0,
  // coarse enough to keep the atomic counter out of the way
  bondedChunkSize = std::max<std::size_t>(1, numberOfMolecules / (8 * threads));
  const std::size_t chunks = (numberOfMolecules + bondedChunkSize - 1) / bondedChunkSize;
  chunkEnergies.assign(chunks, RunningEnergy{});
  chunkStrain.assign(chunks, double3x3{});
  chunkCorrection.assign(chunks, double3x3{});
  atomCenterOfMass.resize(system.spanOfMoleculeAtoms().size());
  partitionedMolecules = numberOfMolecules;
}

void SpatialDecompositionForceEngine::runOnTeam(const std::function<void(std::size_t, std::size_t)>& body)
{
  if (!team || team->size() == 1)
  {
    body(0, 1);
    return;
  }
  const std::size_t members = team->size();
  team->run([&](std::size_t member) { body(member, members); });
}

void SpatialDecompositionForceEngine::refreshStaticSums(const System& system)
{
  // position-independent sums over the atoms (the MD drivers keep the molecules fixed): recomputed only when
  // the number of atoms changes
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  staticNetCharge = 0.0;
  staticScaledCountPerType.assign(system.forceField.pseudoAtoms.size(), 0.0);
  for (const Atom& atom : atoms)
  {
    staticNetCharge += atom.scalingCoulomb * atom.charge;
    staticScaledCountPerType[atom.type] += atom.scalingVDW;
  }
  staticSumsAtoms = atoms.size();
}

void SpatialDecompositionForceEngine::initialize(System& system)
{
  std::string reason;
  if (!supports(system, reason))
  {
    throw std::runtime_error(std::format("[Spatial decomposition]: the force engine does not support {}\n", reason));
  }

  const ForceField& forceField = system.forceField;
  double cutoff = forceField.cutOffMoleculeVDW;
  if (forceField.useCharge) cutoff = std::max(cutoff, forceField.cutOffCoulomb);

  cellList.setup(system.simulationBox, cutoff, settings.verletSkin, settings.numberOfThreads, settings.domainGrid);
  team = std::make_unique<WorkerTeam>(settings.numberOfThreads);

  deviceKernel = settings.pairDevice != PairDevice::CPU;
  useMesh = forceField.usesEwaldFourier();
  deviceMesh = deviceKernel && useMesh && settings.deviceMesh;
  if (useMesh && !deviceMesh)
  {
    pppm.initialize(system.simulationBox, forceField.EwaldAlpha, settings.meshSpacing, settings.interpolationOrder,
                    settings.numberOfThreads, Units::CoulombicConversionFactor);
  }

  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  fx.assign(numberOfAtoms, 0.0);
  fy.assign(numberOfAtoms, 0.0);
  fz.assign(numberOfAtoms, 0.0);
  refreshStaticSums(system);
  localForce.assign(settings.numberOfThreads, {});
  clusterKernelDouble.resize(settings.numberOfThreads);
  clusterKernelMixed.resize(settings.numberOfThreads);
  deviceBonded = false;
  deviceBondedFallback.clear();
  if (deviceKernel)
  {
    if (!DeviceStep::available(settings.pairDevice))
    {
      throw std::runtime_error(std::format(
          "[Spatial decomposition]: 'PairDevice' is '{}' but no {} device is available on this machine ({})\n",
          pairDeviceName(settings.pairDevice), pairDeviceName(settings.pairDevice),
          DeviceStep::unavailableReason(settings.pairDevice)));
    }
    devicePairs.initialize(settings.pairDevice);
  }
  prepareKernel(system);
  if (deviceKernel && !fastKernel)
  {
    throw std::runtime_error(
        "[Spatial decomposition]: the device pair kernel covers plain Lennard-Jones (truncated or shifted) between "
        "fully coupled atoms with Ewald or no electrostatics; use 'PairDevice' : 'CPU' for this system\n");
  }
  if (deviceMesh)
  {
    devicePairs.enableMesh(PPPM::chooseMesh(system.simulationBox, settings.meshSpacing, settings.interpolationOrder),
                           settings.interpolationOrder, forceField.EwaldAlpha, Units::CoulombicConversionFactor);
  }
  if (deviceKernel && settings.deviceBonded)
  {
    if (DeviceStep::supportsBonded(system, deviceBondedFallback))
    {
      devicePairs.enableBonded(system, forceField.EwaldAlpha, Units::CoulombicConversionFactor,
                               forceField.useCharge && useMesh);
      deviceBonded = true;
    }
  }
  // the resident integrator needs the complete forces on the device
  residentEnabled = false;
  residentValid = false;
  residentHostCurrent = true;
  residentFallback.clear();
  resident = DeviceResident{};
  residentPendingScale = {};
  residentKinetic = {};
  residentSteps = 0;
  if (deviceKernel && settings.resident)
  {
    if (useMesh && !deviceMesh)
    {
      residentFallback = "the mesh on the host";
    }
    else if (!deviceBonded)
    {
      residentFallback = deviceBondedFallback.empty()
                             ? std::string("the bonded terms on the host")
                             : std::format("the bonded terms on the host (the device kernels do not cover {})",
                                           deviceBondedFallback);
    }
    else if (DeviceResident::supports(system, residentFallback))
    {
      resident.initialize(devicePairs.deviceContext());
      residentEnabled = true;
    }
  }
  rebuildRequested.assign(settings.numberOfThreads, 1);
  influenceUpdateRequested = 0;
  threadEnergies.assign(settings.numberOfThreads, RunningEnergy{});
  threadStrain.assign(settings.numberOfThreads, double3x3{});
  threadCorrection.assign(settings.numberOfThreads, double3x3{});
  prepareBondedWork(system);
  timing = Timings{};
  initializedFlag = true;
}

void SpatialDecompositionForceEngine::refreshCutoffs(System& system)
{
  // the cutoffs may have been re-derived (automatic cutoffs after a cell change): the lists follow them
  const ForceField& forceField = system.forceField;
  double cutoff = forceField.cutOffMoleculeVDW;
  if (forceField.useCharge) cutoff = std::max(cutoff, forceField.cutOffCoulomb);
  if (cutoff != cellList.cutoff)
  {
    cellList.setup(system.simulationBox, cutoff, settings.verletSkin, settings.numberOfThreads, settings.domainGrid);
  }
  if (fastCoulomb && (forceField.EwaldAlpha != ewaldTable.alpha || !ewaldTable.spans(forceField.cutOffCoulomb)))
  {
    // the kernel choice depends on whether the table can span the (new) cutoff; the cluster lists of the
    // specialised kernel follow the choice at the next build
    prepareKernel(system);
    cellList.numberOfBuilds = 0;
  }
}

RunningEnergy SpatialDecompositionForceEngine::computeGradients(System& system, bool withVirial)
{
  if (!initializedFlag) initialize(system);

  const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  virialRequested = withVirial;
  // the host state is the authoritative one again (the resident integrator uploads it at its next step)
  residentValid = false;
  residentHostCurrent = true;
  if (residentEnabled) devicePairs.setResident(false);

  refreshCutoffs(system);

  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  if (numberOfAtoms != fx.size())
  {
    fx.assign(numberOfAtoms, 0.0);
    fy.assign(numberOfAtoms, 0.0);
    fz.assign(numberOfAtoms, 0.0);
    cellList.numberOfBuilds = 0;  // forces a rebuild
    if (deviceKernel) devicePairs.setExclusions(system);
  }
  if (system.moleculeData.size() != partitionedMolecules || atomCenterOfMass.size() != numberOfAtoms)
  {
    prepareBondedWork(system);
  }
  team->resetWork();

  for (RunningEnergy& energy : threadEnergies) energy = RunningEnergy{};
  for (double3x3& strain : threadStrain) strain = double3x3{};
  for (double3x3& correction : threadCorrection) correction = double3x3{};
  if (deviceBonded)
  {
    // the chunk accumulators are not written when the device does the bonded work
    for (RunningEnergy& energy : chunkEnergies) energy = RunningEnergy{};
    for (double3x3& strain : chunkStrain) strain = double3x3{};
    for (double3x3& correction : chunkCorrection) correction = double3x3{};
  }
  reciprocalEnergy = 0.0;

  team->run([&](std::size_t thread) { step(thread, system); });

  RunningEnergy total{};
  for (const RunningEnergy& energy : threadEnergies) total += energy;
  for (const RunningEnergy& energy : chunkEnergies) total += energy;
  double3x3 strain{};
  double3x3 correction{};
  if (withVirial)
  {
    for (std::size_t t = 0; t < settings.numberOfThreads; ++t)
    {
      strain += threadStrain[t];
      correction += threadCorrection[t];
    }
    for (std::size_t c = 0; c < chunkStrain.size(); ++c)
    {
      strain += chunkStrain[c];
      correction += chunkCorrection[c];
    }
  }
  total = finishStep(system, total, strain, correction);

  ++timing.steps;
  timing.total += std::chrono::steady_clock::now() - begin;
  if (deviceKernel)
  {
    timing.device = devicePairs.deviceTime();
    timing.devicePrune = devicePairs.pruneTime();
    timing.deviceMesh = devicePairs.meshTime();
    timing.deviceBonded = devicePairs.bondedTime();
    timing.deviceBuild = devicePairs.buildTime();
  }
  return total;
}

// the mesh energy with the net-charge correction and, when requested, the molecular pressure tensor from the
// strain derivatives (pairs, exclusions, the mesh) and the atomic-to-molecular virial correction of the step
RunningEnergy SpatialDecompositionForceEngine::finishStep(const System& system, RunningEnergy total, double3x3 strain,
                                                           double3x3 correction)
{
  if (staticSumsAtoms != system.spanOfMoleculeAtoms().size()) refreshStaticSums(system);
  const double netCharge = staticNetCharge;
  if (useMesh)
  {
    total.ewald_fourier += reciprocalEnergy;

    // Net-charge correction (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)), position independent
    const double singleIonSum = deviceMesh ? deviceSingleIonSum : pppm.singleIonFourierSum();
    const double uIon =
        -(singleIonSum - Units::CoulombicConversionFactor * system.forceField.EwaldAlpha / std::sqrt(std::numbers::pi));
    total.ewald_fourier += uIon * netCharge * netCharge;
  }

  if (virialRequested)
  {
    // Assemble the molecular pressure tensor exactly like System::computeMolecularPressure: the strain
    // derivatives of the inter-molecular pairs, the reciprocal sum (including the net-charge correction) and the
    // exclusions, minus the tail correction on the diagonal, corrected from the atomic to the molecular
    // (center-of-mass) virial, negated and symmetrized.
    if (deviceMesh)
    {
      strain += -(deviceReciprocalStrain - (netCharge * netCharge) * deviceSingleIonStrain);
    }
    else if (useMesh)
    {
      strain += -(pppm.reciprocalStrainTensor() - (netCharge * netCharge) * pppm.singleIonStrainTensor());
    }

    // tail correction to the pressure virial, summed over pseudo-atom types instead of atom pairs
    const ForceField& forceField = system.forceField;
    const std::vector<double>& scaledCountPerType = staticScaledCountPerType;
    const double preFactor = -2.0 * std::numbers::pi / (3.0 * system.simulationBox.volume);
    double tailCorrection = 0.0;
    for (std::size_t typeA = 0; typeA < scaledCountPerType.size(); ++typeA)
    {
      if (scaledCountPerType[typeA] == 0.0) continue;
      for (std::size_t typeB = 0; typeB < scaledCountPerType.size(); ++typeB)
      {
        if (scaledCountPerType[typeB] == 0.0) continue;
        tailCorrection += scaledCountPerType[typeA] * scaledCountPerType[typeB] * preFactor *
                          forceField(typeA, typeB).tailCorrectionPressure;
      }
    }
    strain.ax -= tailCorrection;
    strain.by -= tailCorrection;
    strain.cz -= tailCorrection;

    double3x3 tensor = -(strain - correction);
    double temp = 0.5 * (tensor.ay + tensor.bx);
    tensor.ay = tensor.bx = temp;
    temp = 0.5 * (tensor.az + tensor.cx);
    tensor.az = tensor.cx = temp;
    temp = 0.5 * (tensor.bz + tensor.cy);
    tensor.bz = tensor.cy = temp;
    pressureTensor = tensor;
  }
  return total;
}

// list rebuild of the resident step from the host positions (the binning and the slot layout are host work):
// the device lists, the slot layout of the integrator and the slot positions of the new layout
void SpatialDecompositionForceEngine::residentRebuild(System& system)
{
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  const SimulationBox& box = system.simulationBox;
  const double3 widths = box.perpendicularWidths();
  const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
  if (cellList.listCutoff > halfWidth)
  {
    throw std::runtime_error(
        std::format("[Spatial decomposition]: cutoff + skin ({:.3f} A) exceeds half the smallest perpendicular "
                    "width of the box ({:.3f} A); use a smaller cutoff or 'VerletSkin'\n",
                    cellList.listCutoff, halfWidth));
  }
  cellList.updateCellGrid(box);
  cellList.bin(box, system.spanOfMoleculeAtoms(), system.components);
  ++timing.rebuilds;
  devicePairs.beginBuild(cellList, box, 1);
  resident.setLayout(devicePairs, cellList, box);
  resident.enqueuePack();
  devicePairs.finishBuild();
  timing.rebuild += std::chrono::steady_clock::now() - start;
}

RunningEnergy SpatialDecompositionForceEngine::residentVelocityVerlet(System& system)
{
  if (!initializedFlag) initialize(system);
  if (!residentEnabled)
  {
    throw std::runtime_error("[Spatial decomposition]: the resident integrator is not enabled for this system\n");
  }
  const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  virialRequested = true;
  DeviceContext& context = devicePairs.deviceContext();
  const double dt = system.timeStep;
  const SimulationBox& box = system.simulationBox;
  refreshCutoffs(system);

  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  if (numberOfAtoms != fx.size())
  {
    fx.assign(numberOfAtoms, 0.0);
    fy.assign(numberOfAtoms, 0.0);
    fz.assign(numberOfAtoms, 0.0);
    cellList.numberOfBuilds = 0;
    devicePairs.setExclusions(system);
  }
  if (system.moleculeData.size() != partitionedMolecules || atomCenterOfMass.size() != numberOfAtoms)
  {
    prepareBondedWork(system);
  }
  bool rebuild = cellList.numberOfBuilds == 0 || cellList.numberOfAtoms != numberOfAtoms;

  // the host state becomes the device state after a host-side evaluation (stage start, restart)
  if (!residentValid)
  {
    devicePairs.setResident(true);
    if (rebuild)
    {
      // bin the host positions and lay out the device lists; the slot positions follow with the upload
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      cellList.updateCellGrid(box);
      cellList.bin(box, system.spanOfMoleculeAtoms(), system.components);
      ++timing.rebuilds;
      devicePairs.beginBuild(cellList, box, 1);
      devicePairs.finishBuild();
      timing.rebuild += std::chrono::steady_clock::now() - start;
      rebuild = false;
    }
    resident.upload(system, devicePairs, cellList);
    context.wait(resident.enqueuePack());
    residentValid = true;
    residentHostCurrent = true;
    residentPendingScale = {};
    residentKinetic.translational = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    residentKinetic.rotational =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                    system.framework, system.spanOfFrameworkGroupData());
    // the kinetic virial of the barostat (molecular coupling): sum_k M_k |V_k|^2 of the centres of mass
    residentKinetic.virialTrace = 0.0;
    if (system.thermobarostat.has_value())
    {
      const std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
      for (const Molecule& molecule : system.moleculeData)
      {
        const Component& component = system.components[molecule.componentId];
        double3 momentum{};
        if (component.rigid)
        {
          momentum = molecule.mass * molecule.velocity;
        }
        else
        {
          for (std::size_t b = 0; b < molecule.numberOfAtoms; ++b)
            momentum += component.definedAtoms[b].second * dynamics[molecule.atomIndex + b].velocity;
        }
        residentKinetic.virialTrace += double3::dot(momentum, momentum) * molecule.invMass;
      }
    }
  }

  // first half: the barostat chain and the first half-kick of the cell velocity (from the virial of the current
  // configuration and the kinetic virial of the current velocities), the thermostat factor of the end of the
  // last step (deferred) times the one of this step, the barostat coupling of the velocities, kick, drift (of
  // the cell and the centres of mass), free rotor, cartesian positions, slot positions and the displacement check
  DeviceResident::Scaling scaling = residentPendingScale;
  residentPendingScale = {};
  DeviceResident::Coupling coupling{};
  Thermobarostat* barostat = system.thermobarostat.has_value() ? &*system.thermobarostat : nullptr;
  if (barostat != nullptr)
  {
    const double chainScale = barostat->chainStep(barostat->barostatKineticEnergy());
    barostat->logVolumeVelocity *= chainScale;
    barostat->cellVelocity = barostat->cellVelocity * chainScale;
  }
  if (system.thermostat.has_value())
  {
    const std::pair<double, double> factor =
        system.thermostat->NoseHooverNVT(residentKinetic.translational, residentKinetic.rotational);
    scaling.translational *= factor.first;
    scaling.rotational *= factor.second;
  }
  if (barostat != nullptr)
  {
    coupling = residentBarostatFirstHalf(system, pressureTensor.trace());
    resident.setBox(box);
  }
  std::chrono::steady_clock::time_point waitStart = std::chrono::steady_clock::now();
  context.wait(resident.enqueueFirstHalf(scaling, dt, coupling));
  timing.deviceWait += std::chrono::steady_clock::now() - waitStart;
  DeviceResident::Displacement displacement = resident.collectDisplacement();
  if (displacement.sinceBuild > 0.25 * cellList.skin * cellList.skin) rebuild = true;
  if (rebuild)
  {
    resident.downloadPositions(system);
    residentRebuild(system);
    displacement.sinceCompaction = 0.0;
  }

  // forces, then the second half: molecular gradients and torques, kick, barostat coupling of the velocities,
  // kinetic energies and the kinetic virial
  devicePairs.enqueueResident(box, displacement.sinceCompaction > devicePairs.pruneDisplacementThreshold());
  DeviceEvent secondHalf = resident.enqueueSecondHalf(dt, coupling);
  waitStart = std::chrono::steady_clock::now();
  const DeviceStep::Results results =
      devicePairs.wait([&] { secondHalf = resident.enqueueSecondHalf(dt, coupling); });
  context.wait(secondHalf);
  timing.deviceWait += std::chrono::steady_clock::now() - waitStart;
  const DeviceResident::Kinetic kinetic = resident.collectKinetic();
  resident.acceptVelocities();
  residentHostCurrent = false;

  // the thermostat at the end of the step: the factor is applied with the first half of the next step, the
  // reported kinetic energies (and the kinetic virial of the next step) are those of the scaled velocities
  if (system.thermostat.has_value())
  {
    const std::pair<double, double> factor = system.thermostat->NoseHooverNVT(kinetic.translational, kinetic.rotational);
    residentPendingScale.translational = factor.first;
    residentPendingScale.rotational = factor.second;
  }
  residentKinetic.translational =
      residentPendingScale.translational * residentPendingScale.translational * kinetic.translational;
  residentKinetic.rotational = residentPendingScale.rotational * residentPendingScale.rotational * kinetic.rotational;
  residentKinetic.virialTrace =
      residentPendingScale.translational * residentPendingScale.translational * kinetic.virialTrace;

  // energies and the pressure tensor of the step
  RunningEnergy total{};
  total.moleculeMoleculeVDW = results.energyVDW;
  total.moleculeMoleculeCharge = results.energyCharge;
  total.ewald_self = results.bonded.self;
  total.ewald_exclusion = results.bonded.exclusion;
  total.bond = results.bonded.bond;
  total.bend = results.bonded.bend;
  total.torsion = results.bonded.torsion;
  total.improperTorsion = results.bonded.improperTorsion;
  total.intraVDW = results.bonded.intraVDW;
  total.intraCoul = results.bonded.intraCoulomb;
  reciprocalEnergy = 0.0;
  if (deviceMesh)
  {
    reciprocalEnergy = results.reciprocalEnergy;
    deviceReciprocalStrain = results.reciprocalStrain;
    deviceSingleIonSum = results.singleIonSum;
    deviceSingleIonStrain = results.singleIonStrain;
  }
  total = finishStep(system, total, results.pairStrain + results.bonded.exclusionStrain, results.bonded.correction);
  total.translationalKineticEnergy = residentKinetic.translational;
  total.rotationalKineticEnergy = residentKinetic.rotational;
  if (system.thermostat.has_value()) total.NoseHooverEnergy = system.thermostat->getEnergy();

  // the barostat at the end of the step: the second half-kick of the cell velocity from the virial of the new
  // configuration and the kinetic virial of the propagated (unscaled) velocities, then the chain
  if (barostat != nullptr)
  {
    const double mtkFactor =
        1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat->translationalDegreesOfFreedom));
    const double scalarForce =
        3.0 *
        (pressureTensor.trace() + mtkFactor * kinetic.virialTrace - 3.0 * barostat->pressure * box.volume) /
        barostat->logVolumeMass;
    barostat->logVolumeVelocity += 0.5 * dt * scalarForce;
    const double finalScale = barostat->chainStep(barostat->barostatKineticEnergy());
    barostat->logVolumeVelocity *= finalScale;
    barostat->cellVelocity = barostat->cellVelocity * finalScale;
    total.thermobarostatEnergy = barostat->energy(box.volume);
  }

  ++timing.steps;
  ++residentSteps;
  timing.total += std::chrono::steady_clock::now() - begin;
  timing.device = devicePairs.deviceTime();
  timing.devicePrune = devicePairs.pruneTime();
  timing.deviceMesh = devicePairs.meshTime();
  timing.deviceBonded = devicePairs.bondedTime();
  timing.deviceBuild = devicePairs.buildTime();
  return total;
}

// The isotropic Martyna-Tobias-Klein barostat of the resident step (molecular coupling), the host part of the
// first half: x = ln V, xddot = 3 G_epsilon / W with G_epsilon = virial + alpha 2K - 3 P V, the velocity
// propagator S = exp(-dt/2 (xdot/3 + xdot/N_f)) of the centre-of-mass velocities, the cell propagation
// cell' = exp(dt xdot/3) cell with the drift factor dt phi_1(dt xdot/3) of the centres of mass, and the
// force-field parameters that follow the cell (the same sequence as the NPT step of the driver on the host).
DeviceResident::Coupling SpatialDecompositionForceEngine::residentBarostatFirstHalf(System& system,
                                                                                    double pressureVirialTrace)
{
  Thermobarostat& barostat = *system.thermobarostat;
  const double dt = system.timeStep;
  const double mtkFactor =
      1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
  const double scalarForce = 3.0 *
                             (pressureVirialTrace + mtkFactor * residentKinetic.virialTrace -
                              3.0 * barostat.pressure * system.simulationBox.volume) /
                             barostat.logVolumeMass;
  barostat.logVolumeVelocity += 0.5 * dt * scalarForce;
  const double rate = barostat.logVolumeVelocity / 3.0;
  const double3x3 cellRate(rate, rate, rate);

  DeviceResident::Coupling coupling{};
  coupling.enabled = true;
  coupling.propagator = velocityPropagator(cellRate, 0.5 * dt, barostat.translationalDegreesOfFreedom).ax;
  const double argument = dt * rate;
  coupling.cellFactor = std::exp(argument);
  coupling.driftFactor = dt * (argument != 0.0 ? std::expm1(argument) / argument : 1.0);

  const double3x3 cell = coupling.cellFactor * system.simulationBox.cell;
  if (!std::isfinite(cell.determinant()) || cell.determinant() <= 1.0e-10)
  {
    throw std::runtime_error("[Spatial decomposition]: the barostat produced an invalid or singular cell\n");
  }
  system.simulationBox = SimulationBox(cell);
  const ForceField& forceField = system.forceField;
  const double3 widths = system.simulationBox.perpendicularWidths();
  const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
  const double requiredWidth =
      2.0 * std::max({forceField.cutOffFrameworkVDWAutomatic ? 0.0 : forceField.cutOffFrameworkVDW,
                      forceField.cutOffMoleculeVDWAutomatic ? 0.0 : forceField.cutOffMoleculeVDW,
                      forceField.cutOffCoulombAutomatic ? 0.0 : forceField.cutOffCoulomb});
  if (2.0 * halfWidth <= requiredWidth)
  {
    throw std::runtime_error(
        std::format("[Spatial decomposition]: the barostat cell violates the minimum-image cutoff requirement "
                    "(widths: {}, {}, {}; required > {})\n",
                    widths.x, widths.y, widths.z, requiredWidth));
  }
  if (cellList.listCutoff > halfWidth)
  {
    throw std::runtime_error(
        std::format("[Spatial decomposition]: the box shrank so that cutoff + skin ({:.3f} A) exceeds half the "
                    "smallest perpendicular width ({:.3f} A); use a smaller cutoff or 'VerletSkin'\n",
                    cellList.listCutoff, halfWidth));
  }
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  barostat.logVolumePosition += dt * barostat.logVolumeVelocity;
  return coupling;
}

void SpatialDecompositionForceEngine::downloadResidentState(System& system)
{
  if (!residentEnabled || !residentValid || residentHostCurrent) return;
  // the deferred thermostat factor belongs to the state
  resident.enqueueScale(residentPendingScale);
  residentPendingScale = {};
  resident.download(system, devicePairs, cellList);
  residentHostCurrent = true;
}

void SpatialDecompositionForceEngine::step(std::size_t thread, System& system)
{
  WorkerTeam& workers = *team;
  const std::size_t threads = settings.numberOfThreads;
  const SimulationBox& box = system.simulationBox;
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const bool timed = (thread == 0);
  std::chrono::steady_clock::time_point mark = std::chrono::steady_clock::now();
  auto lap = [&](std::chrono::duration<double>& bucket)
  {
    if (!timed) return;
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    bucket += now - mark;
    mark = now;
  };

  // phase 0: refresh positions of the owned atoms, decide whether the lists must be rebuilt. A change of the box
  // (barostat) by itself does not invalidate the lists: the atoms move with the cell and the displacement check
  // covers that motion, and the pair distances are always taken with the minimum image of the current box. The
  // box must only stay large enough for the minimum-image lists. The mesh influence function does depend on the
  // cell: thread 0 decides here whether it must be recomputed (a flag read by all threads after the barrier).
  const bool boxChanged = cellList.boxChanged(box);
  workers.phase(
      [&]
      {
        if (thread == 0 && boxChanged)
        {
          const double3 widths = box.perpendicularWidths();
          const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
          if (cellList.listCutoff > halfWidth)
          {
            throw std::runtime_error(
                std::format("[Spatial decomposition]: the box shrank so that cutoff + skin ({:.3f} A) exceeds half "
                            "the smallest perpendicular width ({:.3f} A); use a smaller cutoff or 'VerletSkin'\n",
                            cellList.listCutoff, halfWidth));
          }
        }
        if (thread == 0)
        {
          influenceUpdateRequested =
              useMesh && !deviceMesh && pppm.influenceFunctionOutdated(box, system.forceField.EwaldAlpha) ? 1 : 0;
          if (influenceUpdateRequested != 0)
          {
            pppm.beginInfluenceFunction(box, system.forceField.EwaldAlpha, threads);
          }
        }
        bool request = cellList.numberOfBuilds == 0 || cellList.numberOfAtoms != atoms.size();
        if (!request) request = cellList.refreshPositionsAndCheck(thread, atoms);
        rebuildRequested[thread] = request ? 1 : 0;
      });

  // phase 0b (only after a cell change, NPT): the influence function G(m) over the half spectrum is expensive on
  // a fine mesh (~K^3 / 2 wave vectors with an exponential each), so its x-slabs are computed by all threads and
  // the single-ion sums are reduced by thread 0. Only thread 0 reads the reduced values later in the step.
  if (influenceUpdateRequested != 0)
  {
    lap(timing.rebuild);
    workers.phase([&] { pppm.computeInfluenceSlab(thread, threads); });
    if (thread == 0)
    {
      pppm.finishInfluenceFunction();
      ++timing.influenceUpdates;
    }
    lap(timing.influence);
  }

  bool rebuild = false;
  for (std::size_t t = 0; t < threads; ++t) rebuild = rebuild || (rebuildRequested[t] != 0);

  // phase 1 (only when needed): serial binning, then parallel list construction. With the device pairs, thread 0
  // lays out the device slots and starts the list build on the device instead.
  if (rebuild)
  {
    workers.phase(
        [&]
        {
          if (thread == 0)
          {
            cellList.updateCellGrid(box);
            cellList.bin(box, atoms, system.components);
            ++timing.rebuilds;
            if (deviceKernel) devicePairs.beginBuild(cellList, box, threads);
          }
        });
    if (!deviceKernel)
    {
      workers.phase(
          [&]
          {
            cellList.buildLists(thread, box);
            if (fastKernel && usesClusterKernel()) buildClusterLists(thread);
          });
    }
  }
  lap(timing.rebuild);
  // phase 1b (device pairs): thread 0 waits for a new list on the device while every thread stages the wrapped
  // positions of its owned atoms for the device
  if (deviceKernel)
  {
    workers.phase(
        [&]
        {
          if (thread == 0 && rebuild) devicePairs.finishBuild();
          devicePairs.packPositions(thread, cellList.domains[thread].ownedAtoms, cellList, box);
        });
    lap(timing.pack);
  }

  // phase 2: gather the compact positions (owned atoms and shifted ghost images), short-range pairs of the owned
  // atoms and, with the mesh, charge spreading into the private sub-box buffer of the thread
  RunningEnergy& energy = threadEnergies[thread];
  workers.phase(
      [&]
      {
        if (deviceKernel)
        {
          // the device evaluates all pairs (and, when enabled, the mesh and the bonded terms) while the threads
          // continue with whatever stayed on the host
          if (thread == 0) devicePairs.enqueue(box);
        }
        else
        {
          cellList.gatherPositions(thread, box);
          pairPhase(thread, system, energy);
        }
        if (useMesh && !deviceMesh)
        {
          const CellList::DomainLists& domain = cellList.domains[thread];
          pppm.spread(thread, domain.ownedAtoms, cellList.x.data(), cellList.y.data(), cellList.z.data(),
                      cellList.charge.data(), cellList.scalingCoulomb.data());
        }
      });
  lap(timing.pairs);

  // phase 3: the owners collect the ghost forces of the other threads; with the mesh, assemble the charge mesh
  // from the sub-box buffers (parallel x-slabs)
  if (threads > 1 && (!deviceKernel || (useMesh && !deviceMesh)))
  {
    workers.phase(
        [&]
        {
          if (!deviceKernel) collectGhostForces(thread);
          if (useMesh && !deviceMesh) pppm.reduceMeshes(thread, threads);
        });
  }
  lap(timing.mesh);

  // phase 4: the bonded terms and the charge self / exclusion corrections, by molecule, in chunks from a shared
  // counter. They start from zeroed gradients (the pair + mesh gradients are added in phase 5) and need nothing
  // from the mesh, so with the mesh they run on the free threads while thread 0 drives the FFTs: forward
  // transform || bonded, influence function on all threads, backward transform || bonded (thread 0 joins the
  // remaining chunks after each transform).
  // With the device pairs, thread 0 collects the device results after its share of the bonded work (the wait
  // is idle time only when the device is slower than the mesh and bonded phases together).
  auto collectDeviceResults = [&]
  {
    if (!deviceKernel || thread != 0) return;
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    const DeviceStep::Results results = devicePairs.wait();
    energy.moleculeMoleculeVDW += results.energyVDW;
    energy.moleculeMoleculeCharge += results.energyCharge;
    if (virialRequested) threadStrain[0] += results.pairStrain;
    if (deviceMesh)
    {
      reciprocalEnergy = results.reciprocalEnergy;
      deviceReciprocalStrain = results.reciprocalStrain;
      deviceSingleIonSum = results.singleIonSum;
      deviceSingleIonStrain = results.singleIonStrain;
    }
    if (deviceBonded)
    {
      energy.ewald_self += results.bonded.self;
      energy.ewald_exclusion += results.bonded.exclusion;
      energy.bond += results.bonded.bond;
      energy.bend += results.bonded.bend;
      energy.torsion += results.bonded.torsion;
      energy.improperTorsion += results.bonded.improperTorsion;
      energy.intraVDW += results.bonded.intraVDW;
      energy.intraCoul += results.bonded.intraCoulomb;
      if (virialRequested)
      {
        threadStrain[0] += results.bonded.exclusionStrain;
        threadCorrection[0] += results.bonded.correction;
      }
    }
    timing.deviceWait += std::chrono::steady_clock::now() - start;
  };
  if (useMesh && !deviceMesh)
  {
    workers.phase(
        [&]
        {
          if (thread == 0) pppm.forwardTransform();
          if (!deviceBonded) bondedWork(system);
        });
    workers.phase([&] { pppm.applyInfluence(thread, threads, virialRequested); });
    workers.phase(
        [&]
        {
          if (thread == 0) reciprocalEnergy = pppm.backwardTransform(threads);
          if (!deviceBonded) bondedWork(system);
          collectDeviceResults();
        });
  }
  else
  {
    workers.phase(
        [&]
        {
          if (!deviceBonded) bondedWork(system);
          collectDeviceResults();
        });
  }
  lap(timing.bonded);

  // phase 5: mesh gradients of the owned atoms, then add the pair + mesh gradients into the system
  workers.phase([&] { scatterPhase(thread, system); });
  lap(timing.mesh);
}

void SpatialDecompositionForceEngine::prepareKernel(const System& system)
{
  const ForceField& forceField = system.forceField;
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  numberOfPseudoAtomTypes = forceField.pseudoAtoms.size();

  // The specialised kernel covers plain 12-6 Lennard-Jones (truncated or shifted) between all pseudo-atom types
  // present, fully coupled atoms (no lambda scaling), and Ewald or no electrostatics.
  std::vector<std::uint8_t> present(numberOfPseudoAtomTypes, 0);
  bool unitScaling = true;
  for (const Atom& atom : atoms)
  {
    present[atom.type] = 1;
    if (atom.scalingVDW != 1.0 || atom.scalingCoulomb != 1.0) unitScaling = false;
  }
  bool lennardJonesOnly = true;
  lennardJones.assign(numberOfPseudoAtomTypes * numberOfPseudoAtomTypes, LennardJonesPair{});
  for (std::size_t a = 0; a < numberOfPseudoAtomTypes; ++a)
  {
    for (std::size_t b = 0; b < numberOfPseudoAtomTypes; ++b)
    {
      const VDWParameters& parameters = forceField(a, b);
      LennardJonesPair& pair = lennardJones[a * numberOfPseudoAtomTypes + b];
      if (parameters.type == VDWParameters::Type::LennardJones)
      {
        pair.epsilon4 = 4.0 * parameters.parameters.x;
        pair.inverseSigma2 = 1.0 / (parameters.parameters.y * parameters.parameters.y);
        pair.shift = parameters.shift;
      }
      else if (parameters.type == VDWParameters::Type::None)
      {
        pair = LennardJonesPair{};
      }
      else if (present[a] && present[b])
      {
        lennardJonesOnly = false;
      }
    }
  }

  fastCoulomb = false;
  if (forceField.useCharge)
  {
    if (forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      if (!ewaldTable.matches(forceField.EwaldAlpha, forceField.cutOffCoulomb))
      {
        ewaldTable.build(forceField.EwaldAlpha, forceField.cutOffCoulomb);
      }
      if (ewaldTable.spans(forceField.cutOffCoulomb))
      {
        fastCoulomb = true;
      }
      else
      {
        lennardJonesOnly = false;  // cutoff beyond the table cap: the generic kernel handles the tail exactly
      }
    }
    else
    {
      lennardJonesOnly = false;  // the real-space shifted schemes go through the generic kernel
    }
  }
  fastKernel = lennardJonesOnly && unitScaling;

  if (fastKernel && deviceKernel)
  {
    devicePairs.setParameters(lennardJones, numberOfPseudoAtomTypes, forceField.useCharge && fastCoulomb,
                              forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb, Units::CoulombicConversionFactor,
                              forceField.EwaldAlpha, settings.verletSkin, settings.pruneSkin);
    devicePairs.setEwaldAlpha(forceField.EwaldAlpha);
    devicePairs.setExclusions(system);
  }
  else if (fastKernel && usesClusterKernel())
  {
    const bool charge = forceField.useCharge && fastCoulomb;
    const EwaldRealSpaceTable* table = charge ? &ewaldTable : nullptr;
    if (mixedPrecision())
    {
      clusterKernelMixed.setParameters(lennardJones, numberOfPseudoAtomTypes, charge, forceField.cutOffMoleculeVDW,
                                       forceField.cutOffCoulomb, Units::CoulombicConversionFactor, table,
                                       forceField.EwaldAlpha, settings.verletSkin, settings.pruneSkin);
    }
    else
    {
      clusterKernelDouble.setParameters(lennardJones, numberOfPseudoAtomTypes, charge, forceField.cutOffMoleculeVDW,
                                        forceField.cutOffCoulomb, Units::CoulombicConversionFactor, table,
                                        forceField.EwaldAlpha, settings.verletSkin, settings.pruneSkin);
    }
  }
}

void SpatialDecompositionForceEngine::buildClusterLists(std::size_t thread)
{
  const CellList::DomainLists& domain = cellList.domains[thread];
  if (mixedPrecision())
  {
    clusterKernelMixed.buildLists(thread, domain);
  }
  else
  {
    clusterKernelDouble.buildLists(thread, domain);
  }
}

void SpatialDecompositionForceEngine::pairPhase(std::size_t thread, const System& system, RunningEnergy& energy)
{
  if (!fastKernel)
  {
    pairLoop<false>(thread, system, energy);
  }
  else if (usesClusterKernel())
  {
    clusterPairLoop(thread, energy);
  }
  else
  {
    pairLoop<true>(thread, system, energy);
  }
}

void SpatialDecompositionForceEngine::clusterPairLoop(std::size_t thread, RunningEnergy& energy)
{
  const CellList::DomainLists& domain = cellList.domains[thread];
  const std::size_t owned = domain.ownedAtoms.size();
  const bool withVirial = virialRequested;
  double energyVDW = 0.0;
  double energyCharge = 0.0;
  double3x3 strain{};
  std::vector<LocalForce>& forces = localForce[thread];

  if (mixedPrecision())
  {
    forces.assign(clusterKernelMixed.paddedLocalAtoms(thread), LocalForce{});
    clusterKernelMixed.refreshPositions(thread, domain);
    clusterKernelMixed.pruneIfNeeded(thread, domain);
    clusterKernelMixed.compute(thread, forces.data(), withVirial, energyVDW, energyCharge, strain);
  }
  else
  {
    forces.assign(clusterKernelDouble.paddedLocalAtoms(thread), LocalForce{});
    clusterKernelDouble.refreshPositions(thread, domain);
    clusterKernelDouble.pruneIfNeeded(thread, domain);
    clusterKernelDouble.compute(thread, forces.data(), withVirial, energyVDW, energyCharge, strain);
  }

  // owned forces into the shared arrays (this thread is the only writer of its atoms), plus this thread's own
  // images (an owned atom seen through a periodic shift)
  const LocalForce* force = forces.data();
  for (std::size_t k = 0; k < owned; ++k)
  {
    const std::uint32_t i = domain.ownedAtoms[k];
    fx[i] = force[k].x;
    fy[i] = force[k].y;
    fz[i] = force[k].z;
  }
  if (domain.imagesByOwner.size() > thread)
  {
    for (const std::uint32_t slot : domain.imagesByOwner[thread])
    {
      const std::uint32_t i = domain.imageAtom[slot];
      const LocalForce& f = force[owned + slot];
      fx[i] += f.x;
      fy[i] += f.y;
      fz[i] += f.z;
    }
  }

  energy.moleculeMoleculeVDW += energyVDW;
  energy.moleculeMoleculeCharge += energyCharge;
  if (withVirial) threadStrain[thread] += strain;
}

template <bool Fast>
void SpatialDecompositionForceEngine::pairLoop(std::size_t thread, const System& system, RunningEnergy& energy)
{
  const ForceField& forceField = system.forceField;
  const CellList::DomainLists& domain = cellList.domains[thread];
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double coulombFactor = Units::CoulombicConversionFactor;

  const std::size_t owned = domain.ownedAtoms.size();
  const std::size_t local = domain.numberOfLocalAtoms();
  const LocalAtom* positions = domain.positions.data();
  const std::uint16_t* type = domain.localType.data();
  const double* scalingVDW = domain.localScalingVDW.data();
  const double* scalingCoulomb = domain.localScalingCoulomb.data();
  const std::uint32_t* neighbourStart = domain.neighbourStart.data();
  const std::uint32_t* neighbourList = domain.neighbourList.data();
  const LennardJonesPair* lj = lennardJones.data();
  const std::size_t numberOfTypes = numberOfPseudoAtomTypes;
  const EwaldRealSpaceTable& table = ewaldTable;

  std::vector<LocalForce>& forces = localForce[thread];
  forces.assign(local, LocalForce{});
  LocalForce* force = forces.data();

  const bool withVirial = virialRequested;
  double sxx = 0.0, syx = 0.0, szx = 0.0, sxy = 0.0, syy = 0.0, szy = 0.0, sxz = 0.0, syz = 0.0, szz = 0.0;
  double energyVDW = 0.0;
  double energyCharge = 0.0;

  for (std::size_t k = 0; k < owned; ++k)
  {
    const LocalAtom& atomI = positions[k];
    const double xi = atomI.x;
    const double yi = atomI.y;
    const double zi = atomI.z;
    const double chargeI = atomI.charge;
    const std::size_t typeI = type[k];
    const LennardJonesPair* ljRow = lj + typeI * numberOfTypes;
    const double scalingVDWI = scalingVDW[k];
    const double scalingCoulombI = scalingCoulomb[k];
    double fxi = 0.0, fyi = 0.0, fzi = 0.0;

    const std::uint32_t begin = neighbourStart[k];
    const std::uint32_t end = neighbourStart[k + 1];
    for (std::uint32_t n = begin; n < end; ++n)
    {
      const std::uint32_t j = neighbourList[n];
      const LocalAtom& atomJ = positions[j];
      const double dx = xi - atomJ.x;
      const double dy = yi - atomJ.y;
      const double dz = zi - atomJ.z;
      const double rr = dx * dx + dy * dy + dz * dz;

      // pair energy and gradient factor (dU/dr / r)
      double factor = 0.0;
      bool interacting = false;
      if (rr < cutOffVDWSquared)
      {
        if constexpr (Fast)
        {
          const LennardJonesPair& p = ljRow[type[j]];
          const double s2 = rr * p.inverseSigma2;
          const double s6 = s2 * s2 * s2;
          const double rri3 = 1.0 / s6;
          const double rri6 = rri3 * rri3;
          energyVDW += p.epsilon4 * (rri6 - rri3) - p.shift;
          factor += 12.0 * p.epsilon4 * rri3 * (0.5 - rri3) / rr;
        }
        else
        {
          const Potentials::PairDerivatives<1> factors =
              Potentials::potentialVDW<1>(forceField, scalingVDWI, scalingVDW[j], rr, typeI, type[j]);
          energyVDW += factors.energy;
          factor += factors.firstDerivativeFactor;
        }
        interacting = true;
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        const double chargeProduct = chargeI * atomJ.charge;
        if (chargeProduct != 0.0)
        {
          if constexpr (Fast)
          {
            double u, dudrr;
            table.evaluate(rr, u, dudrr);
            const double prefactor = coulombFactor * chargeProduct;
            energyCharge += prefactor * u;
            factor += 2.0 * prefactor * dudrr;
          }
          else if (scalingCoulombI * scalingCoulomb[j] != 0.0)
          {
            // a Coulomb-decoupled pair (fractional molecule at lambda <= 0.5) contributes neither energy nor
            // force; skipping it keeps soft-core overlaps away from the 1/r singularity of the non-Ewald methods
            const double r = std::sqrt(rr);
            const Potentials::PairDerivatives<1> factors = Potentials::potentialCoulomb<1>(
                forceField, scalingCoulombI, scalingCoulomb[j], r, chargeI, atomJ.charge);
            energyCharge += factors.energy;
            factor += factors.firstDerivativeFactor;
          }
          interacting = true;
        }
      }
      if (!interacting) continue;

      const double gx = factor * dx;
      const double gy = factor * dy;
      const double gz = factor * dz;
      fxi += gx;
      fyi += gy;
      fzi += gz;
      LocalForce& fj = force[j];
      fj.x -= gx;
      fj.y -= gy;
      fj.z -= gz;
      if (withVirial)
      {
        sxx += gx * dx;
        syx += gy * dx;
        szx += gz * dx;
        sxy += gx * dy;
        syy += gy * dy;
        szy += gz * dy;
        sxz += gx * dz;
        syz += gy * dz;
        szz += gz * dz;
      }
    }
    force[k].x += fxi;
    force[k].y += fyi;
    force[k].z += fzi;
  }

  // owned forces into the shared arrays (this thread is the only writer of its atoms), plus this thread's own
  // images (an owned atom seen through a periodic shift)
  for (std::size_t k = 0; k < owned; ++k)
  {
    const std::uint32_t i = domain.ownedAtoms[k];
    fx[i] = force[k].x;
    fy[i] = force[k].y;
    fz[i] = force[k].z;
  }
  if (domain.imagesByOwner.size() > thread)
  {
    for (const std::uint32_t slot : domain.imagesByOwner[thread])
    {
      const std::uint32_t i = domain.imageAtom[slot];
      const LocalForce& f = force[owned + slot];
      fx[i] += f.x;
      fy[i] += f.y;
      fz[i] += f.z;
    }
  }

  energy.moleculeMoleculeVDW += energyVDW;
  energy.moleculeMoleculeCharge += energyCharge;
  if (withVirial)
  {
    double3x3& strain = threadStrain[thread];
    strain.ax += sxx;
    strain.bx += syx;
    strain.cx += szx;
    strain.ay += sxy;
    strain.by += syy;
    strain.cy += szy;
    strain.az += sxz;
    strain.bz += syz;
    strain.cz += szz;
  }
}

void SpatialDecompositionForceEngine::collectGhostForces(std::size_t thread)
{
  const std::size_t threads = settings.numberOfThreads;
  for (std::size_t t = 0; t < threads; ++t)
  {
    if (t == thread) continue;
    const CellList::DomainLists& other = cellList.domains[t];
    if (other.imagesByOwner.size() <= thread) continue;
    const std::size_t ownedByOther = other.ownedAtoms.size();
    const LocalForce* force = localForce[t].data();
    for (const std::uint32_t slot : other.imagesByOwner[thread])
    {
      const std::uint32_t j = other.imageAtom[slot];
      const LocalForce& f = force[ownedByOther + slot];
      fx[j] += f.x;
      fy[j] += f.y;
      fz[j] += f.z;
    }
  }
}

namespace
{
inline void addOuterProduct(double3x3& tensor, const double3& arm, const double3& gradient)
{
  tensor.ax += arm.x * gradient.x;
  tensor.ay += arm.x * gradient.y;
  tensor.az += arm.x * gradient.z;
  tensor.bx += arm.y * gradient.x;
  tensor.by += arm.y * gradient.y;
  tensor.bz += arm.y * gradient.z;
  tensor.cx += arm.z * gradient.x;
  tensor.cy += arm.z * gradient.y;
  tensor.cz += arm.z * gradient.z;
}

/// The scaled (1-4) pairs of one molecule with the intramolecular pair model (Potentials::intraMolecularVDW /
/// intraMolecularCoulomb at the pair's scaling). The cell lists leave these pairs out together with the excluded
/// ones; every other pair of the molecule is evaluated by the pair kernels like a pair of two molecules.
void addScaledPairGradient(RunningEnergy& energy, const ForceField& forceField, const SimulationBox& box,
                           const IntraMolecularExclusions& exclusions, std::span<const Atom> atoms,
                           std::span<AtomDynamics> dynamics)
{
  for (const IntraMolecularExclusions::ScaledPair& pair : exclusions.scaledPairs)
  {
    const Atom& atomA = atoms[pair.atomA];
    const Atom& atomB = atoms[pair.atomB];
    const double3 dr = box.applyPeriodicBoundaryConditions(atomA.position - atomB.position);
    const double rr = double3::dot(dr, dr);

    double factor = 0.0;
    const Potentials::PairDerivatives<1> vdw = Potentials::intraMolecularVDW<1>(
        forceField, pair.scalingVDW, rr, static_cast<std::size_t>(atomA.type), static_cast<std::size_t>(atomB.type));
    energy.intraVDW += vdw.energy;
    factor += vdw.firstDerivativeFactor;
    if (forceField.useCharge)
    {
      const Potentials::PairDerivatives<1> coulomb =
          Potentials::intraMolecularCoulomb<1>(forceField, pair.scalingCoulomb, atomA.scalingCoulomb,
                                               atomB.scalingCoulomb, std::sqrt(rr), atomA.charge, atomB.charge);
      energy.intraCoul += coulomb.energy;
      factor += coulomb.firstDerivativeFactor;
    }
    if (factor == 0.0) continue;
    const double3 gradient = factor * dr;
    dynamics[pair.atomA].gradient += gradient;
    dynamics[pair.atomB].gradient -= gradient;
  }
}
}  // namespace

void SpatialDecompositionForceEngine::bondedWork(System& system)
{
  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t chunk = bondedChunkSize;
  const bool withVirial = virialRequested;

  // every chunk is taken by exactly one thread per step; its sums go to the chunk's own slots
  WorkerTeam& workers = *team;
  for (std::size_t begin = workers.claimWork(chunk); begin < numberOfMolecules; begin = workers.claimWork(chunk))
  {
    const std::size_t end = std::min(begin + chunk, numberOfMolecules);
    RunningEnergy energy{};
    double3x3 strain{};
    double3x3 correction{};
    for (std::size_t m = begin; m < end; ++m)
    {
      const Molecule& molecule = system.moleculeData[m];
      const Component& component = system.components[molecule.componentId];
      std::span<const Atom> moleculeAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      std::span<AtomDynamics> moleculeDynamics = dynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms);

      // the gradients of the system are assembled here first (the pair + mesh gradients are added in the
      // scatter phase, after this work is complete)
      for (AtomDynamics& atomDynamics : moleculeDynamics) atomDynamics.gradient = double3(0.0, 0.0, 0.0);

      Interactions::addChargeSelfEnergy(energy, forceField, moleculeAtoms);
      Interactions::addIntraMolecularChargeExclusionGradient(energy, forceField, box, system.components, moleculeAtoms,
                                                             moleculeDynamics, withVirial ? &strain : nullptr);

      if (withVirial)
      {
        // atomic-to-molecular virial correction of the non-bonded gradients about the mass-weighted center of
        // mass: the exclusion part here, the pair + mesh part in the scatter phase (which reads the center of
        // mass stored per atom); the bonded gradients added below cancel in the molecular virial
        double totalMass = 0.0;
        double3 com(0.0, 0.0, 0.0);
        for (const Atom& atom : moleculeAtoms)
        {
          const double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
          com += mass * atom.position;
          totalMass += mass;
        }
        com = com / totalMass;
        for (std::size_t k = 0; k < moleculeAtoms.size(); ++k)
        {
          atomCenterOfMass[molecule.atomIndex + k] = com;
          addOuterProduct(correction, moleculeAtoms[k].position - com, moleculeDynamics[k].gradient);
        }
      }

      // the bonded terms and the scaled pairs (the non-excluded, unscaled pairs of the molecule are in the cell
      // lists and evaluated by the pair kernels with the molecule-molecule pairs)
      energy += component.intraMolecularPotentials.computeInternalBondedGradient(box, moleculeAtoms, moleculeDynamics);
      addScaledPairGradient(energy, forceField, box, component.intraMolecularPotentials.exclusions, moleculeAtoms,
                            moleculeDynamics);
    }
    const std::size_t index = begin / chunk;
    chunkEnergies[index] = energy;
    chunkStrain[index] = strain;
    chunkCorrection[index] = correction;
  }
}

void SpatialDecompositionForceEngine::scatterPhase(std::size_t thread, System& system)
{
  const CellList::DomainLists& domain = cellList.domains[thread];
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  if (deviceKernel && deviceBonded)
  {
    // the device forces are the complete gradients (pairs, mesh, exclusions, bonded terms) and the device did
    // the virial correction
    for (const std::uint32_t i : domain.ownedAtoms)
    {
      dynamics[cellList.sortedToOriginal[i]].gradient = devicePairs.force(i);
    }
    return;
  }
  if (deviceKernel)
  {
    // the pair gradients of the owned atoms from the device (the mesh gradients are added below)
    for (const std::uint32_t i : domain.ownedAtoms)
    {
      const double3 gradient = devicePairs.force(i);
      fx[i] = gradient.x;
      fy[i] = gradient.y;
      fz[i] = gradient.z;
    }
  }
  if (useMesh && !deviceMesh)
  {
    pppm.interpolate(domain.ownedAtoms, cellList.x.data(), cellList.y.data(), cellList.z.data(), cellList.charge.data(),
                     cellList.scalingCoulomb.data(), fx.data(), fy.data(), fz.data());
  }
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  if (virialRequested)
  {
    double3x3 correction{};
    for (const std::uint32_t i : domain.ownedAtoms)
    {
      const std::uint32_t original = cellList.sortedToOriginal[i];
      const double3 gradient(fx[i], fy[i], fz[i]);
      dynamics[original].gradient += gradient;
      addOuterProduct(correction, atoms[original].position - atomCenterOfMass[original], gradient);
    }
    threadCorrection[thread] += correction;
  }
  else
  {
    for (const std::uint32_t i : domain.ownedAtoms)
    {
      dynamics[cellList.sortedToOriginal[i]].gradient += double3(fx[i], fy[i], fz[i]);
    }
  }
}

SpatialDecompositionForceEngine::Validation SpatialDecompositionForceEngine::validate(System& system)
{
  Validation result{};

  // reference: the exact O(N^2) + direct Ewald code
  const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  RunningEnergy reference = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
  std::vector<double3> referenceGradients;
  referenceGradients.reserve(system.atomDynamics.size());
  for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) referenceGradients.push_back(dynamics.gradient);

  const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  RunningEnergy engine = computeGradients(system);
  const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  result.referenceSeconds = std::chrono::duration<double>(t1 - t0).count();
  result.engineSeconds = std::chrono::duration<double>(t2 - t1).count();

  result.engineEnergy = engine.potentialEnergy();
  result.referenceEnergy = reference.potentialEnergy();
  result.engineReciprocalEnergy = engine.ewald_fourier;
  result.referenceReciprocalEnergy = reference.ewald_fourier;

  double sumSquaredDifference = 0.0;
  double sumSquared = 0.0;
  double maximum = 0.0;
  std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  for (std::size_t i = 0; i < dynamics.size(); ++i)
  {
    const double3 difference = dynamics[i].gradient - referenceGradients[i];
    const double d2 = double3::dot(difference, difference);
    sumSquaredDifference += d2;
    sumSquared += double3::dot(referenceGradients[i], referenceGradients[i]);
    maximum = std::max(maximum, std::sqrt(d2));
  }
  const double n = static_cast<double>(std::max<std::size_t>(1, dynamics.size()));
  result.maximumGradientDifference = maximum;
  result.rmsGradientDifference = std::sqrt(sumSquaredDifference / n);
  result.rmsGradient = std::sqrt(sumSquared / n);
  return result;
}

std::string SpatialDecompositionForceEngine::writeStatus() const
{
  std::string result;
  result += std::format("Spatial-decomposition force engine\n");
  result += std::format(
      "================================================================================================================"
      "========\n");
  result += std::format("    threads: {}\n", settings.numberOfThreads);
  if (deviceKernel)
  {
    result += devicePairs.status();
  }
  else if (fastKernel && usesClusterKernel())
  {
    auto describe = [&](const auto& kernel)
    {
      using Kernel = std::remove_cvref_t<decltype(kernel)>;
      std::string text = std::format(
          "    pair kernel: cluster {} x {} ({} precision) Lennard-Jones{}\n", Kernel::clusterI, Kernel::clusterJ,
          mixedPrecision() ? "mixed" : "double",
          fastCoulomb ? (Kernel::analyticEwald ? " + analytic Ewald real space" : " + tabulated Ewald real space")
                      : "");
      if (kernel.pruning())
      {
        text += std::format("    blocks: {} in the Verlet list, {} after pruning (skin {:.3f} A, {} prunes)\n",
                            kernel.numberOfBlocks(), kernel.numberOfPrunedBlocks(), settings.pruneSkin,
                            kernel.numberOfPrunes());
      }
      else
      {
        text += std::format("    blocks: {} (no pruning)\n", kernel.numberOfBlocks());
      }
      return text;
    };
    result += mixedPrecision() ? describe(clusterKernelMixed) : describe(clusterKernelDouble);
  }
  else if (fastKernel)
  {
    result += std::format("    pair kernel: scalar (double precision) Lennard-Jones{}\n",
                          fastCoulomb ? " + tabulated Ewald real space" : "");
  }
  else
  {
    result += std::format("    pair kernel: generic (Potentials::potentialVDW / potentialCoulomb)\n");
  }
  result += cellList.status();
  if (useMesh && !deviceMesh)
  {
    result += pppm.status();
  }
  else if (!useMesh)
  {
    result += std::format("    no reciprocal-space sum (charge method without Ewald Fourier part)\n");
  }
  if (deviceKernel && !deviceBonded)
  {
    result += std::format("    bonded terms on the host{}\n",
                          deviceBondedFallback.empty()
                              ? std::string{}
                              : std::format(" (the device kernels do not cover {})", deviceBondedFallback));
  }
  if (residentEnabled)
  {
    result +=
        "    resident integrator: positions, velocities and molecule records on the device in double-float "
        "(emulated double), forces and torques in single precision; the host keeps the thermostat and barostat "
        "chains\n";
  }
  else if (deviceKernel && settings.resident)
  {
    result += std::format("    integration on the host ({})\n", residentFallback);
  }
  result += "\n";
  return result;
}

std::string SpatialDecompositionForceEngine::writeTimings() const
{
  std::string result;
  result += std::format("Spatial-decomposition force engine timings\n");
  result += std::format(
      "================================================================================================================"
      "========\n");
  result += std::format("    force evaluations:        {}\n", timing.steps);
  if (residentSteps > 0)
  {
    result += std::format("    resident MD steps:        {} (integrator on the device, two synchronizations per step)\n",
                          residentSteps);
  }
  result +=
      std::format("    neighbour-list rebuilds:  {} ({:.2f} steps per rebuild)\n", timing.rebuilds,
                  timing.rebuilds > 0 ? static_cast<double>(timing.steps) / static_cast<double>(timing.rebuilds) : 0.0);
  if (timing.influenceUpdates > 0)
  {
    result += std::format("    influence-function updates: {} (cell changes)\n", timing.influenceUpdates);
  }
  result += std::format("    total:                    {:14.4f} [s]\n", timing.total.count());
  result += std::format("    rebuilds:                 {:14.4f} [s]\n", timing.rebuild.count());
  if (timing.influenceUpdates > 0)
  {
    result += std::format("    influence function:       {:14.4f} [s]\n", timing.influence.count());
  }
  if (deviceKernel)
  {
    result += std::format("    device list builds:       {:14.4f} [s] (layout, upload and build; in 'rebuilds')\n",
                          timing.deviceBuild.count());
    result += std::format("    device pair kernel:       {:14.4f} [s] (kernel time, sampled after the list builds)\n",
                          timing.device.count());
    if (devicePairs.pruning())
    {
      result += std::format("    device list compaction:   {:14.4f} [s] (kernel time, sampled after the list builds)\n",
                            timing.devicePrune.count());
    }
    if (deviceMesh)
    {
      result += std::format("    device mesh (PPPM):       {:14.4f} [s] (spreading, FFTs, influence, interpolation)\n",
                            timing.deviceMesh.count());
    }
    if (deviceBonded)
    {
      result +=
          std::format("    device bonded + exclusions: {:12.4f} [s] (kernel time, sampled after the list builds)\n",
                      timing.deviceBonded.count());
    }
    result += std::format("    wait for the device:      {:14.4f} [s] (idle time of thread 0 after the host work)\n",
                          timing.deviceWait.count());
    result += std::format("    position staging:         {:14.4f} [s] (packing by all threads, wait for a new list)\n",
                          timing.pack.count());
    result += std::format("    {:<30}{:14.4f} [s] (thread 0; includes the timing samples)\n",
                          (useMesh && !deviceMesh) ? "device enqueue + spreading:" : "device enqueue:",
                          timing.pairs.count());
  }
  else
  {
    result += std::format("    pairs + spreading:        {:14.4f} [s]\n", timing.pairs.count());
  }
  result += std::format("    reduction, interpolation, scatter: {:5.4f} [s]\n", timing.mesh.count());
  if (useMesh && !deviceMesh)
  {
    result += std::format("    FFTs || bonded + exclusions: {:11.4f} [s]\n", timing.bonded.count());
  }
  else if (deviceBonded)
  {
    result += std::format("    device wait phase:        {:14.4f} [s]\n", timing.bonded.count());
  }
  else
  {
    result += std::format("    bonded + exclusions:      {:14.4f} [s]\n", timing.bonded.count());
  }
  if (timing.steps > 0)
  {
    result += std::format("    per force evaluation:     {:14.4f} [ms]\n",
                          1000.0 * timing.total.count() / static_cast<double>(timing.steps));
  }
  result += "\n";
  return result;
}
