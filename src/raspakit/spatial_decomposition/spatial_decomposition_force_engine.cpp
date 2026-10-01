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
import interactions_ewald;
import integrators_update;
import thermobarostat;
import system;
import spatial_decomposition_settings;
import spatial_decomposition_cell_list;
import spatial_decomposition_pppm;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_cluster_kernel;
import spatial_decomposition_worker_team;

SpatialDecompositionForceEngine::SpatialDecompositionForceEngine(const SpatialDecompositionSettings& s) : settings(s)
{
  settings.numberOfThreads = std::max<std::size_t>(1, settings.numberOfThreads);
}

SpatialDecompositionForceEngine::~SpatialDecompositionForceEngine() = default;

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
  if (system.forceField.useDualCutOff)
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
  atomCenterOfMass.resize(system.spanOfMoleculeAtoms().size());
  partitionedMolecules = numberOfMolecules;
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

  useMesh = forceField.usesEwaldFourier();
  if (useMesh)
  {
    pppm.initialize(system.simulationBox, forceField.EwaldAlpha, settings.meshSpacing, settings.interpolationOrder,
                    settings.numberOfThreads, Units::CoulombicConversionFactor);
  }

  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  fx.assign(numberOfAtoms, 0.0);
  fy.assign(numberOfAtoms, 0.0);
  fz.assign(numberOfAtoms, 0.0);
  localForce.assign(settings.numberOfThreads, {});
  clusterKernelDouble.resize(settings.numberOfThreads);
  clusterKernelMixed.resize(settings.numberOfThreads);
  prepareKernel(system);
  rebuildRequested.assign(settings.numberOfThreads, 1);
  influenceUpdateRequested = 0;
  threadEnergies.assign(settings.numberOfThreads, RunningEnergy{});
  threadStrain.assign(settings.numberOfThreads, double3x3{});
  threadCorrection.assign(settings.numberOfThreads, double3x3{});
  prepareBondedWork(system);
  timing = Timings{};
  initializedFlag = true;
}

RunningEnergy SpatialDecompositionForceEngine::computeGradients(System& system, bool withVirial)
{
  if (!initializedFlag) initialize(system);

  const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  virialRequested = withVirial;

  // the cutoffs may have been re-derived (automatic cutoffs after a cell change): the lists follow them
  {
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

  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  if (numberOfAtoms != fx.size())
  {
    fx.assign(numberOfAtoms, 0.0);
    fy.assign(numberOfAtoms, 0.0);
    fz.assign(numberOfAtoms, 0.0);
    cellList.numberOfBuilds = 0;  // forces a rebuild
  }
  if (system.moleculeData.size() != partitionedMolecules || atomCenterOfMass.size() != numberOfAtoms)
  {
    prepareBondedWork(system);
  }
  nextBondedChunk.store(0, std::memory_order_relaxed);

  for (RunningEnergy& energy : threadEnergies) energy = RunningEnergy{};
  for (double3x3& strain : threadStrain) strain = double3x3{};
  for (double3x3& correction : threadCorrection) correction = double3x3{};
  reciprocalEnergy = 0.0;

  team->run([&](std::size_t thread) { step(thread, system); });

  RunningEnergy total{};
  for (const RunningEnergy& energy : threadEnergies) total += energy;

  double netCharge = 0.0;
  if (useMesh)
  {
    total.ewald_fourier += reciprocalEnergy;

    // Net-charge correction (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)), position independent
    for (const Atom& atom : system.spanOfMoleculeAtoms()) netCharge += atom.scalingCoulomb * atom.charge;
    const double uIon = -(pppm.singleIonFourierSum() - Units::CoulombicConversionFactor * system.forceField.EwaldAlpha /
                                                           std::sqrt(std::numbers::pi));
    total.ewald_fourier += uIon * netCharge * netCharge;
  }

  if (withVirial)
  {
    // Assemble the molecular pressure tensor exactly like System::computeMolecularPressure: the strain
    // derivatives of the inter-molecular pairs, the reciprocal sum (including the net-charge correction) and the
    // exclusions, minus the tail correction on the diagonal, corrected from the atomic to the molecular
    // (center-of-mass) virial, negated and symmetrized.
    double3x3 strain{};
    double3x3 correction{};
    for (std::size_t t = 0; t < settings.numberOfThreads; ++t)
    {
      strain += threadStrain[t];
      correction += threadCorrection[t];
    }
    if (useMesh)
    {
      strain += -(pppm.reciprocalStrainTensor() - (netCharge * netCharge) * pppm.singleIonStrainTensor());
    }

    // tail correction to the pressure virial, summed over pseudo-atom types instead of atom pairs
    const ForceField& forceField = system.forceField;
    std::vector<double> scaledCountPerType(forceField.pseudoAtoms.size(), 0.0);
    for (const Atom& atom : system.spanOfMoleculeAtoms()) scaledCountPerType[atom.type] += atom.scalingVDW;
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

  ++timing.steps;
  timing.total += std::chrono::steady_clock::now() - begin;
  return total;
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
              useMesh && pppm.influenceFunctionOutdated(box, system.forceField.EwaldAlpha) ? 1 : 0;
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

  // phase 1 (only when needed): serial binning, then parallel list construction
  if (rebuild)
  {
    workers.phase(
        [&]
        {
          if (thread == 0)
          {
            cellList.updateCellGrid(box);
            cellList.bin(box, atoms);
            ++timing.rebuilds;
          }
        });
    workers.phase(
        [&]
        {
          cellList.buildLists(thread, box);
          if (fastKernel && usesClusterKernel()) buildClusterLists(thread);
        });
  }
  lap(timing.rebuild);

  // phase 2: gather the compact positions (owned atoms and shifted ghost images), short-range pairs of the owned
  // atoms and, with the mesh, charge spreading into the private sub-box buffer of the thread
  RunningEnergy& energy = threadEnergies[thread];
  workers.phase(
      [&]
      {
        cellList.gatherPositions(thread, box);
        pairPhase(thread, system, energy);
        if (useMesh)
        {
          const CellList::DomainLists& domain = cellList.domains[thread];
          pppm.spread(thread, domain.ownedAtoms, cellList.x.data(), cellList.y.data(), cellList.z.data(),
                      cellList.charge.data(), cellList.scalingCoulomb.data());
        }
      });
  lap(timing.pairs);

  // phase 3: the owners collect the ghost forces of the other threads; with the mesh, assemble the charge mesh
  // from the sub-box buffers (parallel x-slabs)
  if (threads > 1)
  {
    workers.phase(
        [&]
        {
          collectGhostForces(thread);
          if (useMesh) pppm.reduceMeshes(thread, threads);
        });
  }
  lap(timing.mesh);

  // phase 4: the bonded terms and the charge self / exclusion corrections, by molecule, in chunks from a shared
  // counter. They start from zeroed gradients (the pair + mesh gradients are added in phase 5) and need nothing
  // from the mesh, so with the mesh they run on the free threads while thread 0 drives the FFTs: forward
  // transform || bonded, influence function on all threads, backward transform || bonded (thread 0 joins the
  // remaining chunks after each transform).
  if (useMesh)
  {
    workers.phase(
        [&]
        {
          if (thread == 0) pppm.forwardTransform();
          bondedWork(thread, system, energy);
        });
    workers.phase([&] { pppm.applyInfluence(thread, threads, virialRequested); });
    workers.phase(
        [&]
        {
          if (thread == 0) reciprocalEnergy = pppm.backwardTransform(threads);
          bondedWork(thread, system, energy);
        });
  }
  else
  {
    workers.phase([&] { bondedWork(thread, system, energy); });
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

  if (fastKernel && usesClusterKernel())
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
}  // namespace

void SpatialDecompositionForceEngine::bondedWork(std::size_t thread, System& system, RunningEnergy& energy)
{
  const ForceField& forceField = system.forceField;
  const SimulationBox& box = system.simulationBox;
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t chunk = bondedChunkSize;

  const bool withVirial = virialRequested;
  double3x3 strain{};
  double3x3 correction{};

  for (std::size_t begin = nextBondedChunk.fetch_add(chunk, std::memory_order_relaxed); begin < numberOfMolecules;
       begin = nextBondedChunk.fetch_add(chunk, std::memory_order_relaxed))
  {
    const std::size_t end = std::min(begin + chunk, numberOfMolecules);
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
      Interactions::addIntraMolecularChargeExclusionGradient(energy, forceField, box, moleculeAtoms, moleculeDynamics,
                                                             withVirial ? &strain : nullptr);

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

      energy += component.intraMolecularPotentials.computeInternalGradient(moleculeAtoms, moleculeDynamics);
    }
  }

  if (withVirial)
  {
    threadStrain[thread] += strain;
    threadCorrection[thread] += correction;
  }
}

void SpatialDecompositionForceEngine::scatterPhase(std::size_t thread, System& system)
{
  const CellList::DomainLists& domain = cellList.domains[thread];
  if (useMesh)
  {
    pppm.interpolate(domain.ownedAtoms, cellList.x.data(), cellList.y.data(), cellList.z.data(), cellList.charge.data(),
                     cellList.scalingCoulomb.data(), fx.data(), fy.data(), fz.data());
  }
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
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
  if (fastKernel && usesClusterKernel())
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
  if (useMesh)
  {
    result += pppm.status();
  }
  else
  {
    result += std::format("    no reciprocal-space sum (charge method without Ewald Fourier part)\n");
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
  result += std::format("    pairs + spreading:        {:14.4f} [s]\n", timing.pairs.count());
  result += std::format("    reduction, interpolation, scatter: {:5.4f} [s]\n", timing.mesh.count());
  if (useMesh)
  {
    result += std::format("    FFTs || bonded + exclusions: {:11.4f} [s]\n", timing.bonded.count());
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
