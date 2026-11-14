module;

module mc_moves_hybridmc;

import std;

import ewald_structure_factor;

import running_energy;
import randomnumbers;
import system;
import atom;
import atom_dynamics;
import molecule;
import component;
import integrators;
import integrators_update;
import integrators_compute;
import thermostat;
import interactions_ewald;
import interactions_intermolecular;
import interactions_framework_molecule;
import mc_moves_move_types;
import mc_moves_volume;
import simulationbox;
import double3;
import spatial_decomposition_settings;
import force_engine;
import spatial_decomposition_integrator;

namespace
{
/// The force engine kept in 'System::hybridMCForceEngine' between moves, with what it was built for: a change of
/// the atom count (swap moves) or of the system it belongs to (a copied System shares the pointer) rebuilds it; a
/// changed box is handled by the engine itself.
struct HybridMCEngine
{
  ForceEngine engine;
  const System* owner{nullptr};
  std::size_t numberOfAtoms{0};
  std::size_t numberOfMolecules{0};
  /// The Verlet skin the engine was built with (clamped to the box of that moment, see 'settingsFor').
  double verletSkin{0.0};

  HybridMCEngine(const System& system, const SpatialDecompositionSettings& settings)
      : engine(settingsFor(system, settings)), verletSkin(settingsFor(system, settings).verletSkin)
  {
  }

  /// The move integrates on the host (the state is reset between the short trajectories and restored on
  /// rejection), the engine supplies the forces only. The Verlet skin is clamped to what the box allows (the MC
  /// cutoffs are often exactly half the box, which leaves no room for the MD default): the trajectories are short,
  /// so the extra list rebuilds of a thin skin cost little.
  static SpatialDecompositionSettings settingsFor(const System& system, SpatialDecompositionSettings settings)
  {
    settings.resident = false;
    double cutoff = system.forceField.cutOffMoleculeVDW;
    if (system.forceField.useCharge) cutoff = std::max(cutoff, system.forceField.cutOffCoulomb);
    const double3 widths = system.simulationBox.perpendicularWidths();
    const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
    settings.verletSkin = std::clamp(halfWidth - cutoff - 1.0e-9 * halfWidth, 0.0, settings.verletSkin);
    return settings;
  }
};

/// The engine for this system, (re)initialised when needed; nullptr when the engine does not cover the system (the
/// move then takes the exact path, and stays on it).
HybridMCEngine* hybridMCEngine(System& system)
{
  if (!system.hybridMCUseForceEngine || system.numberOfHybridMCSteps == 0) return nullptr;
  std::string reason;
  if (!ForceEngine::supports(system, reason))
  {
    system.hybridMCUseForceEngine = false;
    system.hybridMCForceEngine.reset();
    return nullptr;
  }

  std::shared_ptr<HybridMCEngine> holder = std::static_pointer_cast<HybridMCEngine>(system.hybridMCForceEngine);
  // The skin was clamped to the box the engine was built for. A volume move (or the automatic cut-offs re-derived
  // after one) can leave less room than that: cutoff + skin then exceeds half the box and the engine threw at its
  // next list update. An engine whose skin no longer fits is rebuilt for the box as it is now.
  const double fittingSkin = HybridMCEngine::settingsFor(system, system.forceEngineSettings).verletSkin;
  if (!holder || holder->owner != &system || holder->verletSkin > fittingSkin)
  {
    holder = std::make_shared<HybridMCEngine>(system, system.forceEngineSettings);
    system.hybridMCForceEngine = holder;
  }
  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  if (!holder->engine.initialized() || holder->owner != &system || holder->numberOfAtoms != numberOfAtoms ||
      holder->numberOfMolecules != system.moleculeData.size())
  {
    try
    {
      holder->engine.initialize(system);
    }
    catch (const std::exception& e)
    {
      // the engine refuses the system (e.g. the box is too small for the cutoff): the exact path from here on
      std::print(std::cerr, "[Hybrid MC]: force engine unavailable, using the exact integrator: {}", e.what());
      system.hybridMCUseForceEngine = false;
      system.hybridMCForceEngine.reset();
      return nullptr;
    }
    holder->owner = &system;
    holder->numberOfAtoms = numberOfAtoms;
    holder->numberOfMolecules = system.moleculeData.size();
  }
  return holder.get();
}

/// The configurations a trajectory may not end in: a molecule inside a blocking pocket of its component, or a
/// declared chiral centre that the trajectory inverted (a hard constraint, which every move that changes a
/// molecule's internal geometry enforces).
bool trajectoryEndsInRefusedConfiguration(const System& system, std::span<const Molecule> molecules,
                                          std::span<const Atom> atoms)
{
  std::size_t offset = 0;
  for (const Molecule& molecule : molecules)
  {
    const Component& component = system.components[molecule.componentId];
    const std::span<const Atom> moleculeAtoms = atoms.subspan(offset, molecule.numberOfAtoms);
    if (system.insideBlockedPockets(component, moleculeAtoms) ||
        !component.intraMolecularPotentials.preservesChiralCenters(moleculeAtoms))
    {
      return true;
    }
    offset += molecule.numberOfAtoms;
  }
  return false;
}

/**
 * \brief The hybrid MC move with the trajectory integrated by the spatial-decomposition force engine, on the
 * system's own arrays (restored on rejection).
 *
 * The acceptance is exact: velocity Verlet is time-reversible and volume-preserving for any force that depends on
 * the positions only, so the engine forces (cell lists, particle-mesh Ewald, optionally mixed precision) make a
 * valid proposal, and the Metropolis test uses the exact Hamiltonian, U_exact + K, with both U_exact(old) and
 * U_exact(new) a full recomputation (System::computeTotalEnergies, which also refreshes the stored Ewald
 * structure factors for the configuration it evaluates).
 */
std::optional<RunningEnergy> hybridMCMoveWithEngine(RandomNumber& random, System& system, ForceEngine& engine,
                                                    double dt)
{
  const Move::Types move = Move::Types::HybridMC;

  // the rigid-group states of the semi-flexible molecules follow the (authoritative) atom positions, which the
  // other moves change without updating them
  system.initializeGroupData();

  // the state restored on rejection (the move integrates in place)
  const std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  const std::span<AtomDynamics> dynamicsData = system.spanOfMoleculeDynamics();
  const std::span<GroupState> groupDataSpan = system.spanOfGroupData();
  const std::vector<Atom> savedAtoms(atomData.begin(), atomData.end());
  const std::vector<AtomDynamics> savedDynamics(dynamicsData.begin(), dynamicsData.end());
  const std::vector<Molecule> savedMolecules(system.moleculeData);
  const std::vector<GroupState> savedGroups(groupDataSpan.begin(), groupDataSpan.end());
  const double savedTimeStep = system.timeStep;
  // an NVE proposal: no thermostat along the trajectory
  std::optional<Thermostat> savedThermostat = std::exchange(system.thermostat, std::nullopt);
  system.timeStep = dt;

  // The acceptance uses the exact Hamiltonian. This path never sees a framework, an external field or
  // polarization (ForceEngine::supports), but U_old is still taken from a full recomputation so a drifted
  // running energy cannot leak into the Metropolis test.
  const std::chrono::steady_clock::time_point time_old_energy = std::chrono::steady_clock::now();
  const RunningEnergy oldEnergy = system.computeTotalEnergies();
  // 'computeTotalEnergies' refreshes the stored Ewald structure factors and the electric field for the configuration
  // it evaluates: the state to restore on rejection is the one it just left, not the one before it
  const EwaldStructureFactor savedStoredEik(system.storedEik);
  std::vector<double3> savedElectricField;
  if (system.forceField.computePolarization)
  {
    const std::span<const double3> electricField = system.spanOfMoleculeElectricField();
    savedElectricField.assign(electricField.begin(), electricField.end());
  }
  system.mc_moves_cputime[move][Move::Timing::Energy] += (std::chrono::steady_clock::now() - time_old_energy);

  // Boltzmann velocities (the reference energy carries their kinetic energy); a bulk fluid has a free lab frame
  Integrators::initializeVelocities(random, system.moleculeData, atomData, dynamicsData, system.components,
                                    system.temperature, std::nullopt, {}, {}, &system.forceField, groupDataSpan);
  // Total momentum is conserved along this trajectory (no rigid framework, no external field). Removing the
  // centre-of-mass velocity puts the draw on the P = 0 hyperplane, which the reverse move then shares.
  if (!system.hasExternalField)
  {
    Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, atomData, dynamicsData, system.components,
                                                 std::nullopt, {}, {}, &system.forceField, groupDataSpan);
  }
  const double initialKinetic =
      Integrators::computeTranslationalKineticEnergy(system.moleculeData, atomData, dynamicsData, system.components,
                                                     std::nullopt, {}, {}, &system.forceField, groupDataSpan) +
      Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, groupDataSpan);

  // the forces for the first half-kick, then the trajectory
  const std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
  engine.computeGradients(system, false);
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, atomData, dynamicsData, system.components,
                                                        groupDataSpan);
  RunningEnergy trajectoryEnergy{};
  for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    trajectoryEnergy = SpatialDecompositionIntegrator::engineVelocityVerlet(system, engine);
  }
  const std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);
  system.thermostat = std::move(savedThermostat);

  if (trajectoryEndsInRefusedConfiguration(system, system.moleculeData, atomData))
  {
    std::copy(savedAtoms.begin(), savedAtoms.end(), atomData.begin());
    std::copy(savedDynamics.begin(), savedDynamics.end(), dynamicsData.begin());
    std::copy(savedGroups.begin(), savedGroups.end(), groupDataSpan.begin());
    system.moleculeData = savedMolecules;
    system.storedEik = savedStoredEik;
    system.timeStep = savedTimeStep;
    return std::nullopt;
  }

  const double finalKinetic = trajectoryEnergy.translationalKineticEnergy + trajectoryEnergy.rotationalKineticEnergy;
  // the exact energy of the new configuration (this also refreshes the stored Ewald structure factors)
  const std::chrono::steady_clock::time_point time_energy = std::chrono::steady_clock::now();
  RunningEnergy finalEnergy = system.computeTotalEnergies();
  system.mc_moves_cputime[move][Move::Timing::Energy] += (std::chrono::steady_clock::now() - time_energy);
  const double deltaH = (finalEnergy.potentialEnergy() + finalKinetic) - (oldEnergy.potentialEnergy() + initialKinetic);

  // a configuration a CBMC move could not retrace (a pair, framework or external-field energy above the overlap
  // criterion: outside a confining wall, on an atom) is refused as the volume move refuses it
  if (!MC_Moves::volumeTrialOverlaps(system, finalEnergy) && random.uniform() < std::exp(-system.beta * deltaH))
  {
    system.mc_moves_statistics.addAccepted(move);
    // the adaptive HMC step size is a property of the move, not of the system (the MD drivers read 'timeStep')
    system.timeStep = savedTimeStep;
    finalEnergy.translationalKineticEnergy = trajectoryEnergy.translationalKineticEnergy;
    finalEnergy.rotationalKineticEnergy = trajectoryEnergy.rotationalKineticEnergy;
    return finalEnergy;
  }

  // rejected: the old configuration (the engine notices the displacement at its next evaluation)
  std::copy(savedAtoms.begin(), savedAtoms.end(), atomData.begin());
  std::copy(savedDynamics.begin(), savedDynamics.end(), dynamicsData.begin());
  std::copy(savedGroups.begin(), savedGroups.end(), groupDataSpan.begin());
  system.moleculeData = savedMolecules;
  system.storedEik = savedStoredEik;
  if (system.forceField.computePolarization)
  {
    std::copy(savedElectricField.begin(), savedElectricField.end(), system.spanOfMoleculeElectricField().begin());
  }
  system.timeStep = savedTimeStep;
  // the molecule cell list was rebuilt by 'computeTotalEnergies' for the rejected trial positions
  system.invalidateCellList();
  return std::nullopt;
}
}  // namespace

std::string MC_Moves::hybridMCForceEngineReport(const System& system)
{
  const std::shared_ptr<HybridMCEngine> holder = std::static_pointer_cast<HybridMCEngine>(system.hybridMCForceEngine);
  if (!holder || !holder->engine.initialized()) return {};
  return std::format("Hybrid MC: {}{}\n", holder->engine.writeStatus(), holder->engine.writeTimings());
}

std::optional<RunningEnergy> MC_Moves::hybridMCMove(RandomNumber& random, System& system)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::HybridMC;

  system.mc_moves_statistics.addTrial(move);

  const bool flexibleFramework = system.framework && !system.framework->rigid;
  const bool hasMovableFramework = flexibleFramework && system.framework->hasMobileAtoms();

  // The total momentum is conserved along the trajectory of a bulk fluid or of a flexible framework without Fixed
  // atoms (the lab frame is free), and the centre-of-mass velocity is then removed from the drawn velocities. A
  // rigid framework, the Fixed atoms of a mixed one, or an external field pin the lab frame: the momentum is not
  // conserved and is left intact.
  const bool momentumIsConserved =
      (system.numberOfFrameworkAtoms == 0 || (flexibleFramework && system.framework->numberOfFixedAtoms() == 0)) &&
      !system.hasExternalField;

  // Reject when there is nothing to propagate: no molecules without mobile framework atoms, or one single-atom
  // molecule whose only velocity is the centre-of-mass drift that is removed. A single polyatomic molecule has
  // rotational or internal degrees of freedom, and a single atom moves when the lab frame is pinned.
  if (!hasMovableFramework)
  {
    if (system.moleculeData.empty()) return std::nullopt;
    if (system.moleculeData.size() == 1 && momentumIsConserved && system.moleculeData.front().numberOfAtoms == 1)
    {
      return std::nullopt;
    }
  }

  if (HybridMCEngine* engineHolder = hybridMCEngine(system))
  {
    return hybridMCMoveWithEngine(random, system, engineHolder->engine, system.mc_moves_statistics.getMaxChange(move));
  }

  // all copied data: moleculeData, moleculeAtomPositions, moleculeDynamics,
  //                  frameworkAtomPositions/frameworkDynamics (flexible), thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents,
  //                 fixedFrameworkStoredEik (rigid framework only)
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomData.begin(), atomData.end());

  std::span<AtomDynamics> dynamicsData = system.spanOfMoleculeDynamics();
  std::vector<AtomDynamics> moleculeDynamics(dynamicsData.begin(), dynamicsData.end());

  std::vector<Molecule> moleculeData(system.moleculeData);

  // Rigid-body state for every rigid group of every semi-flexible molecule, derived from the current
  // (authoritative) atom positions. Integrated on this local copy and committed on acceptance.
  std::vector<GroupState> groupData;
  {
    std::size_t atomIndex{};
    for (const Molecule& molecule : moleculeData)
    {
      const Component& component = system.components[molecule.componentId];
      if (component.isSemiFlexible())
      {
        std::span<const Atom> span = std::span(&moleculeAtomPositions[atomIndex], molecule.numberOfAtoms);
        for (std::size_t g = 0; g != component.fragmentGraph.fragments.size(); ++g)
        {
          if (component.fragmentGraph.fragments[g].isRigidBody())
          {
            groupData.push_back(component.deriveFragmentState(g, span));
          }
        }
      }
      atomIndex += molecule.numberOfAtoms;
    }
  }
  std::span<GroupState> groupDataSpan(groupData);

  std::span<Atom> frameworkAtomData = system.spanOfFrameworkAtoms();
  std::vector<Atom> frameworkAtomPositions;
  std::span<AtomDynamics> frameworkDynamicsData = system.spanOfFrameworkDynamics();
  std::vector<AtomDynamics> frameworkDynamics;
  if (flexibleFramework)
  {
    frameworkAtomPositions.assign(frameworkAtomData.begin(), frameworkAtomData.end());
    frameworkDynamics.assign(frameworkDynamicsData.begin(), frameworkDynamicsData.end());
  }

  std::span<Atom> trialFrameworkAtoms = flexibleFramework ? std::span<Atom>(frameworkAtomPositions) : frameworkAtomData;
  std::span<AtomDynamics> trialFrameworkDynamics =
      flexibleFramework ? std::span<AtomDynamics>(frameworkDynamics) : std::span<AtomDynamics>{};

  // Rigid mobile groups of a mixed framework: the same per-group state the MD integrator uses.
  std::vector<GroupState> frameworkGroupData;
  if (flexibleFramework && system.framework->numberOfRigidGroups() > 0)
  {
    for (std::size_t g = 0; g != system.framework->groups.size(); ++g)
    {
      if (system.framework->groups[g].isRigidBody())
      {
        frameworkGroupData.push_back(system.framework->deriveGroupState(g, trialFrameworkAtoms));
      }
    }
  }
  std::span<GroupState> frameworkGroupDataSpan(frameworkGroupData);

  // Hybrid MC is an NVE proposal; do not carry over an MD thermostat.
  std::optional<Thermostat> thermostat = std::nullopt;

  // get Timestep from the max change
  double dt = system.mc_moves_statistics.getMaxChange(move);

  // The acceptance uses the exact Hamiltonian. The MD gradient routines carry the pair, Ewald, bonded, cross-link and
  // polymer-link terms only (no external field, no polarization); a trajectory driven by those forces is still a
  // time-reversible, volume-preserving proposal, so the Metropolis test on the exact energies samples the full
  // Boltzmann distribution. 'computeTotalEnergies' also refreshes the stored Ewald structure factors and the
  // electric field for the configuration it evaluates; the current ones are saved for the rejection path.
  time_begin = std::chrono::steady_clock::now();
  const RunningEnergy oldEnergy = system.computeTotalEnergies();
  const EwaldStructureFactor savedStoredEik(system.storedEik);
  std::vector<double3> savedElectricField;
  if (system.forceField.computePolarization)
  {
    const std::span<const double3> electricField = system.spanOfMoleculeElectricField();
    savedElectricField.assign(electricField.begin(), electricField.end());
  }
  system.mc_moves_cputime[move][Move::Timing::Energy] += (std::chrono::steady_clock::now() - time_begin);

  // initialize the velocities according to Boltzmann distribution
  // NOTE: it is important that the reference energy has the initial kinetic energies
  Integrators::initializeVelocities(random, moleculeData, moleculeAtomPositions, moleculeDynamics, system.components,
                                    system.temperature, system.framework, trialFrameworkAtoms, trialFrameworkDynamics,
                                    &system.forceField, groupDataSpan, frameworkGroupDataSpan);

  // Remove the centre-of-mass drift when the momentum is conserved: the P = 0 hyperplane the velocities are drawn
  // on is then invariant along the trajectory, so the reversed final momenta can be drawn by the reverse move.
  if (momentumIsConserved)
  {
    Integrators::removeCenterOfMassVelocityDrift(
        moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, system.framework, trialFrameworkAtoms,
        trialFrameworkDynamics, &system.forceField, groupDataSpan, frameworkGroupDataSpan);
  }

  // Gradients must live on the trial copies: Velocity Verlet's first half-kick uses them.
  RunningEnergy referenceEnergy = Integrators::updateGradients(
      moleculeData, moleculeAtomPositions, moleculeDynamics, trialFrameworkAtoms, system.forceField,
      system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      trialFrameworkDynamics, &system.crossLinks, &system.polymerLinks);
  Integrators::updateCenterOfMassAndQuaternionGradients(moleculeData, moleculeAtomPositions, moleculeDynamics,
                                                        system.components, groupDataSpan, system.framework,
                                                        trialFrameworkDynamics, frameworkGroupDataSpan);
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, system.framework, trialFrameworkAtoms,
      trialFrameworkDynamics, &system.forceField, groupDataSpan, frameworkGroupDataSpan);
  referenceEnergy.rotationalKineticEnergy = Integrators::computeRotationalKineticEnergy(
      moleculeData, system.components, groupDataSpan, system.framework, frameworkGroupDataSpan);
  RunningEnergy currentEnergy = referenceEnergy;

  // integrate for N steps
  time_begin = std::chrono::steady_clock::now();
  for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, dt, thermostat, trialFrameworkAtoms,
        system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
        system.numberOfMoleculesPerComponent, system.framework, trialFrameworkDynamics, groupDataSpan,
        frameworkGroupDataSpan, &system.crossLinks, &system.polymerLinks);
  }
  time_end = std::chrono::steady_clock::now();

  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);

  if (trajectoryEndsInRefusedConfiguration(system, moleculeData, moleculeAtomPositions))
  {
    return std::nullopt;
  }

  // The trial configuration goes into the system for its exact energy (the only energy routines are the
  // system's own); the current one is kept for the rejection path.
  const std::vector<Atom> savedAtoms(atomData.begin(), atomData.end());
  const std::vector<AtomDynamics> savedDynamics(dynamicsData.begin(), dynamicsData.end());
  const std::vector<Molecule> savedMolecules(std::move(system.moleculeData));
  const std::vector<GroupState> savedGroups(std::move(system.groupData));
  const std::vector<GroupState> savedFrameworkGroups(std::move(system.frameworkGroupData));
  std::vector<Atom> savedFrameworkAtoms;
  std::vector<AtomDynamics> savedFrameworkDynamics;
  if (flexibleFramework)
  {
    savedFrameworkAtoms.assign(frameworkAtomData.begin(), frameworkAtomData.end());
    savedFrameworkDynamics.assign(frameworkDynamicsData.begin(), frameworkDynamicsData.end());
  }

  system.moleculeData = std::move(moleculeData);
  system.groupData = std::move(groupData);
  system.frameworkGroupData = std::move(frameworkGroupData);
  std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomData.begin());
  std::copy(moleculeDynamics.begin(), moleculeDynamics.end(), dynamicsData.begin());
  if (flexibleFramework)
  {
    std::copy(frameworkAtomPositions.begin(), frameworkAtomPositions.end(), frameworkAtomData.begin());
    std::copy(frameworkDynamics.begin(), frameworkDynamics.end(), frameworkDynamicsData.begin());
  }
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                        system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                        system.spanOfFrameworkGroupData());

  // the exact energy of the new configuration (this also refreshes the stored Ewald structure factors and the
  // electric field)
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy newEnergy = system.computeTotalEnergies();
  system.mc_moves_cputime[move][Move::Timing::Energy] += (std::chrono::steady_clock::now() - time_begin);

  // the Metropolis rule on the exact Hamiltonian H = U + K: a trajectory that lowers H is always accepted (the
  // proposal is reversible and volume-preserving, so this is detailed balance for the Boltzmann distribution)
  const double initialKinetic = referenceEnergy.translationalKineticEnergy + referenceEnergy.rotationalKineticEnergy;
  const double finalKinetic = currentEnergy.translationalKineticEnergy + currentEnergy.rotationalKineticEnergy;
  const double deltaH = (newEnergy.potentialEnergy() + finalKinetic) - (oldEnergy.potentialEnergy() + initialKinetic);

  // accept or reject based on energy difference
  if (!MC_Moves::volumeTrialOverlaps(system, newEnergy) && random.uniform() < std::exp(-system.beta * deltaH))
  {
    system.mc_moves_statistics.addAccepted(move);

    // The stored Ewald structure factors and the electric field already describe the new configuration. The
    // returned energy is the exact one (tail corrections, external field and polarization included); the caller
    // overwrites system.runningEnergies with it.
    newEnergy.translationalKineticEnergy = currentEnergy.translationalKineticEnergy;
    newEnergy.rotationalKineticEnergy = currentEnergy.rotationalKineticEnergy;
    return newEnergy;
  }

  // rejected: the old configuration, structure factors and electric field
  std::copy(savedAtoms.begin(), savedAtoms.end(), atomData.begin());
  std::copy(savedDynamics.begin(), savedDynamics.end(), dynamicsData.begin());
  system.moleculeData = savedMolecules;
  system.groupData = savedGroups;
  system.frameworkGroupData = savedFrameworkGroups;
  if (flexibleFramework)
  {
    std::copy(savedFrameworkAtoms.begin(), savedFrameworkAtoms.end(), frameworkAtomData.begin());
    std::copy(savedFrameworkDynamics.begin(), savedFrameworkDynamics.end(), frameworkDynamicsData.begin());
  }
  system.storedEik = savedStoredEik;
  if (system.forceField.computePolarization)
  {
    std::copy(savedElectricField.begin(), savedElectricField.end(), system.spanOfMoleculeElectricField().begin());
  }
  // the molecule cell list was rebuilt by 'computeTotalEnergies' for the rejected trial positions
  system.invalidateCellList();
  return std::nullopt;
}
