module;

module mc_moves_hybridmc;

import std;

import running_energy;
import randomnumbers;
import system;
import atom;
import atom_dynamics;
import molecule;
import integrators;
import integrators_update;
import integrators_compute;
import thermostat;
import interactions_ewald;
import mc_moves_move_types;

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

  std::span<Atom> frameworkAtomData = system.spanOfFrameworkAtoms();
  std::vector<Atom> frameworkAtomPositions;
  std::span<AtomDynamics> frameworkDynamicsData = system.spanOfFrameworkDynamics();
  std::vector<AtomDynamics> frameworkDynamics;
  if (flexibleFramework)
  {
    frameworkAtomPositions.assign(frameworkAtomData.begin(), frameworkAtomData.end());
    frameworkDynamics.assign(frameworkDynamicsData.begin(), frameworkDynamicsData.end());
  }

  std::span<Atom> trialFrameworkAtoms =
      flexibleFramework ? std::span<Atom>(frameworkAtomPositions) : frameworkAtomData;
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
                                    &system.forceField);

  // Remove COM drift for bulk fluids and for flexible frameworks (lab frame is free).
  // A rigid framework pins the lab frame, so adsorbate COM momentum is left intact.
  if (system.numberOfFrameworkAtoms == 0 || flexibleFramework)
  {
    Integrators::removeCenterOfMassVelocityDrift(moleculeData, moleculeAtomPositions, moleculeDynamics,
                                                 system.components, system.framework, trialFrameworkAtoms,
                                                 trialFrameworkDynamics, &system.forceField);
  }

  // Gradients must live on the trial copies: Velocity Verlet's first half-kick uses them.
  RunningEnergy referenceEnergy = Integrators::updateGradients(
      moleculeData, moleculeAtomPositions, moleculeDynamics, trialFrameworkAtoms, system.forceField,
      system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      trialFrameworkDynamics);
  Integrators::updateCenterOfMassAndQuaternionGradients(moleculeData, moleculeAtomPositions, moleculeDynamics,
                                                        system.components);
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, system.framework, trialFrameworkAtoms,
      trialFrameworkDynamics, &system.forceField);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
  RunningEnergy currentEnergy = referenceEnergy;

  // integrate for N steps
  time_begin = std::chrono::steady_clock::now();
  for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculeData, moleculeAtomPositions, moleculeDynamics, system.components, dt, thermostat, trialFrameworkAtoms,
        system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
        system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent,
        system.framework, trialFrameworkDynamics);
  }
  time_end = std::chrono::steady_clock::now();

  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);

  double drift = std::abs(currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy());

  // accept or reject based on energy difference
  if (random.uniform() < std::exp(-system.beta * drift))
  {
    system.mc_moves_statistics.addAccepted(move);

    system.moleculeData = moleculeData;
    system.timeStep = dt;

    std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomData.begin());
    std::copy(moleculeDynamics.begin(), moleculeDynamics.end(), dynamicsData.begin());
    if (flexibleFramework)
    {
      std::copy(frameworkAtomPositions.begin(), frameworkAtomPositions.end(), frameworkAtomData.begin());
      std::copy(frameworkDynamics.begin(), frameworkDynamics.end(), frameworkDynamicsData.begin());
    }

    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    return currentEnergy;
  }
  return std::nullopt;
}
