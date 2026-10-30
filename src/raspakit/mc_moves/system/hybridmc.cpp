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

  HybridMCEngine(const System& system, const SpatialDecompositionSettings& settings)
      : engine(settingsFor(system, settings))
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
  if (!holder || holder->owner != &system)
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

/**
 * \brief The hybrid MC move with the trajectory integrated by the spatial-decomposition force engine, on the
 * system's own arrays (restored on rejection).
 *
 * The acceptance is exact: velocity Verlet is time-reversible and volume-preserving for any force that depends on
 * the positions only, so the engine forces (cell lists, particle-mesh Ewald, optionally mixed precision) make a
 * valid proposal, and the Metropolis test uses the exact Hamiltonian, U_exact + K, with U_exact(old) the running
 * energy of the system and U_exact(new) a full recomputation (System::computeTotalEnergies, which also refreshes
 * the stored Ewald structure factors for the new configuration).
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
  const EwaldStructureFactor savedStoredEik(system.storedEik);
  const double savedTimeStep = system.timeStep;
  // an NVE proposal: no thermostat along the trajectory
  std::optional<Thermostat> savedThermostat = std::exchange(system.thermostat, std::nullopt);
  system.timeStep = dt;

  // Boltzmann velocities (the reference energy carries their kinetic energy); a bulk fluid has a free lab frame
  Integrators::initializeVelocities(random, system.moleculeData, atomData, dynamicsData, system.components,
                                    system.temperature, std::nullopt, {}, {}, &system.forceField, groupDataSpan);
  Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, atomData, dynamicsData, system.components,
                                               std::nullopt, {}, {}, &system.forceField, groupDataSpan);
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

  const double finalKinetic = trajectoryEnergy.translationalKineticEnergy + trajectoryEnergy.rotationalKineticEnergy;
  const double initialPotential = system.runningEnergies.potentialEnergy();
  // the exact energy of the new configuration (this also refreshes the stored Ewald structure factors)
  const std::chrono::steady_clock::time_point time_energy = std::chrono::steady_clock::now();
  RunningEnergy finalEnergy = system.computeTotalEnergies();
  system.mc_moves_cputime[move][Move::Timing::Energy] += (std::chrono::steady_clock::now() - time_energy);
  const double deltaH = (finalEnergy.potentialEnergy() + finalKinetic) - (initialPotential + initialKinetic);

  if (random.uniform() < std::exp(-system.beta * deltaH))
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
  system.timeStep = savedTimeStep;
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
  const bool hasMovableFramework = flexibleFramework && system.numberOfFrameworkAtoms > 0;
  // Reject when there is nothing to propagate. For a flexible framework the atoms themselves
  // are degrees of freedom, so a framework-only (or single-molecule) system is allowed.
  if (system.moleculeData.size() <= 1 && !hasMovableFramework)
  {
    return std::nullopt;
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

  // Hybrid MC is an NVE proposal; do not carry over an MD thermostat.
  std::optional<Thermostat> thermostat = std::nullopt;

  // get Timestep from the max change
  double dt = system.mc_moves_statistics.getMaxChange(move);

  // initialize the velocities according to Boltzmann distribution
  // NOTE: it is important that the reference energy has the initial kinetic energies
  Integrators::initializeVelocities(random, moleculeData, moleculeAtomPositions, moleculeDynamics, system.components,
                                    system.temperature, system.framework, trialFrameworkAtoms, trialFrameworkDynamics,
                                    &system.forceField, groupDataSpan);

  // Remove COM drift for bulk fluids and for flexible frameworks (lab frame is free).
  // A rigid framework pins the lab frame, so adsorbate COM momentum is left intact.
  if (system.numberOfFrameworkAtoms == 0 || flexibleFramework)
  {
    Integrators::removeCenterOfMassVelocityDrift(moleculeData, moleculeAtomPositions, moleculeDynamics,
                                                 system.components, system.framework, trialFrameworkAtoms,
                                                 trialFrameworkDynamics, &system.forceField, groupDataSpan);
  }

  // Gradients must live on the trial copies: Velocity Verlet's first half-kick uses them.
  RunningEnergy referenceEnergy = Integrators::updateGradients(
      moleculeData, moleculeAtomPositions, moleculeDynamics, trialFrameworkAtoms, system.forceField,
      system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      trialFrameworkDynamics, &system.crossLinks);
  Integrators::updateCenterOfMassAndQuaternionGradients(moleculeData, moleculeAtomPositions, moleculeDynamics,
                                                        system.components, groupDataSpan);
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, system.framework, trialFrameworkAtoms,
      trialFrameworkDynamics, &system.forceField, groupDataSpan);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculeData, system.components, groupDataSpan);
  RunningEnergy currentEnergy = referenceEnergy;

  // integrate for N steps
  time_begin = std::chrono::steady_clock::now();
  for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, dt, thermostat, trialFrameworkAtoms,
        system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
        system.numberOfMoleculesPerComponent, system.framework, trialFrameworkDynamics, groupDataSpan, {},
        &system.crossLinks);
  }
  time_end = std::chrono::steady_clock::now();

  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);

  // the Metropolis rule on the Hamiltonian: a trajectory that lowers H is always accepted (the proposal is
  // reversible and volume-preserving, so this is detailed balance for the Boltzmann distribution)
  const double deltaH = currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy();

  // accept or reject based on energy difference
  if (random.uniform() < std::exp(-system.beta * deltaH))
  {
    system.mc_moves_statistics.addAccepted(move);

    system.moleculeData = moleculeData;
    system.groupData = groupData;

    std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomData.begin());
    std::copy(moleculeDynamics.begin(), moleculeDynamics.end(), dynamicsData.begin());
    if (flexibleFramework)
    {
      std::copy(frameworkAtomPositions.begin(), frameworkAtomPositions.end(), frameworkAtomData.begin());
      std::copy(frameworkDynamics.begin(), frameworkDynamics.end(), frameworkDynamicsData.begin());
    }

    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData());
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    // The MD gradient/energy routines do not include the LJ tail corrections. The tail term is
    // constant along the NVE trajectory (N, V and the scaling factors do not change), so it cancels
    // in the acceptance rule, but it must be restored in the returned energy: the caller overwrites
    // system.runningEnergies with it, and both the tail energy and its dU/dlambda contribution would
    // otherwise be lost.
    currentEnergy += system.computeTailCorrectionEnergies();
    return currentEnergy;
  }
  return std::nullopt;
}
