module;

module molecular_dynamics;

import std;

import stringutils;
import hardware_info;
import archive;
import system;
import framework;
import randomnumbers;
import input_reader;
import component;
import averages;
import property_loading;
import units;
import property_enthalpy;
import simulationbox;
import forcefield;
import sample_movies;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import atom;
import atom_dynamics;
import molecule;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import property_widom;
import property_simulationbox;
import property_energy;
import property_msd;
import property_vacf;
import mc_moves;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_cputime;
import mc_moves_statistics;
import property_pressure;
import transition_matrix;
import interactions_ewald;
import equation_of_states;
import integrators;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import interpolation_energy_grid;
import simulation_schedule;
import thermobarostat;
import minimization_cell_layout;
import elastic_constants;
import property_elastic_constants_fluctuation;

namespace
{
std::array<double, 6> stressVoigt(const double3x3& stress)
{
  return {stress.ax,
          stress.by,
          stress.cz,
          0.5 * (stress.bz + stress.cy),
          0.5 * (stress.az + stress.cx),
          0.5 * (stress.ay + stress.bx)};
}

void applyVelocityMatrix(System& system, const double3x3& matrix)
{
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<GroupState> groupData = system.spanOfGroupData();
  std::size_t atomIndex{};
  std::size_t groupIndex{};
  for (Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    molecule.velocity = matrix * molecule.velocity;
    if (component.isSemiFlexible())
    {
      // Semi-flexible molecule: the barostat couples to the rigid-group center-of-mass velocities and
      // to the flexible-atom velocities; orientation momenta are not coupled to the cell.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity = matrix * state.velocity;
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
          moleculeDynamics[atomIndex + i].velocity = matrix * moleculeDynamics[atomIndex + i].velocity;
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else if (!component.rigid)
    {
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        moleculeDynamics[atomIndex + i].velocity = matrix * moleculeDynamics[atomIndex + i].velocity;
    }
    atomIndex += molecule.numberOfAtoms;
  }
  if (system.framework && !system.framework->rigid)
    for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics()) dynamics.velocity = matrix * dynamics.velocity;
}

void propagateCell(System& system, const double3x3& cellVelocity)
{
  std::vector<double3> positions;
  std::vector<double3> velocities;
  std::vector<double3*> targets;
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<GroupState> groupData = system.spanOfGroupData();
  std::size_t atomIndex{};
  std::size_t groupIndex{};
  for (Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      positions.push_back(molecule.centerOfMassPosition);
      velocities.push_back(molecule.velocity);
      targets.push_back(&molecule.centerOfMassPosition);
    }
    else if (component.isSemiFlexible())
    {
      // Semi-flexible molecule: the cell drives the rigid-group centers of mass (preserving the
      // internal group geometry; the atoms are rebuilt from the GroupState afterwards) and the
      // flexible atoms individually.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          positions.push_back(state.centerOfMassPosition);
          velocities.push_back(state.velocity);
          targets.push_back(&state.centerOfMassPosition);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          positions.push_back(moleculeAtoms[atomIndex + i].position);
          velocities.push_back(moleculeDynamics[atomIndex + i].velocity);
          targets.push_back(&moleculeAtoms[atomIndex + i].position);
        }
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else
    {
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        positions.push_back(moleculeAtoms[atomIndex + i].position);
        velocities.push_back(moleculeDynamics[atomIndex + i].velocity);
        targets.push_back(&moleculeAtoms[atomIndex + i].position);
      }
    }
    atomIndex += molecule.numberOfAtoms;
  }
  if (system.framework && !system.framework->rigid)
  {
    std::span<Atom> atoms = system.spanOfFrameworkAtoms();
    std::span<AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      positions.push_back(atoms[i].position);
      velocities.push_back(dynamics[i].velocity);
      targets.push_back(&atoms[i].position);
    }
  }

  double3x3 cell = system.simulationBox.cell;
  const bool upper = system.thermobarostat->cellType == CellMinimizationType::RegularUpperTriangle ||
                     system.thermobarostat->cellType == CellMinimizationType::MonoclinicUpperTriangle;
  propagateCellAndPosition(cell, positions, velocities, cellVelocity, system.timeStep, upper);
  if (!std::isfinite(cell.determinant()) || cell.determinant() <= 1.0e-10)
    throw std::runtime_error("Thermobarostat produced an invalid or singular cell");
  for (std::size_t i = 0; i != targets.size(); ++i) *targets[i] = positions[i];
  system.simulationBox = SimulationBox(cell);
  const double3 widths = system.simulationBox.perpendicularWidths();
  const double requiredWidth =
      2.0 * std::max({system.forceField.cutOffFrameworkVDWAutomatic ? 0.0 : system.forceField.cutOffFrameworkVDW,
                      system.forceField.cutOffMoleculeVDWAutomatic ? 0.0 : system.forceField.cutOffMoleculeVDW,
                      system.forceField.cutOffCoulombAutomatic ? 0.0 : system.forceField.cutOffCoulomb});
  if (std::min({widths.x, widths.y, widths.z}) <= requiredWidth)
    throw std::runtime_error(std::format(
        "Thermobarostat cell violates the minimum-image cutoff requirement (widths: {}, {}, {}; required > {})",
        widths.x, widths.y, widths.z, requiredWidth));
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
}

// Offset of a molecule's rigid-group states inside System::groupData: the entries are ordered like
// the molecules (component-major), with one GroupState per rigid group of the component.
std::size_t groupDataOffsetOfMolecule(const System& system, std::size_t componentId, std::size_t moleculeId)
{
  std::size_t offset{};
  for (std::size_t c = 0; c < componentId; ++c)
  {
    offset += system.numberOfMoleculesPerComponent[c] * system.components[c].numberOfRigidGroups();
  }
  return offset + moleculeId * system.components[componentId].numberOfRigidGroups();
}

void refreshParticleNumberDependentState(RandomNumber& random, System& system, std::size_t selectedComponent,
                                         MC_Moves::ParticleExchangeResult result, std::size_t exchangedMolecule)
{
  if (result == MC_Moves::ParticleExchangeResult::Rejected) return;

  const Component& component = system.components[selectedComponent];

  if (result == MC_Moves::ParticleExchangeResult::Deleted && component.numberOfRigidGroups() > 0)
  {
    // drop the rigid-group states of the deleted molecule; the entries of the molecules behind it
    // shift down so that groupData stays aligned with the molecule layout
    const std::size_t offset = groupDataOffsetOfMolecule(system, selectedComponent, exchangedMolecule);
    const auto first = system.groupData.begin() + static_cast<std::ptrdiff_t>(offset);
    system.groupData.erase(first, first + static_cast<std::ptrdiff_t>(component.numberOfRigidGroups()));
  }

  if (result == MC_Moves::ParticleExchangeResult::Inserted)
  {
    const std::size_t selectedMolecule = exchangedMolecule;
    Molecule& molecule =
        system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];
    std::span<AtomDynamics> dynamics =
        system.spanOfMoleculeDynamics().subspan(molecule.atomIndex, molecule.numberOfAtoms);

    if (component.isSemiFlexible())
    {
      // Derive the rigid-group states of the freshly grown molecule from its atom positions and give
      // the groups thermal velocities; flexible atoms get per-atom thermal velocities and rigid-group
      // atoms carry no independent velocity (mirrors Integrators::initializeVelocities).
      std::span<Atom> atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
      std::vector<GroupState> states;
      states.reserve(component.numberOfRigidGroups());
      for (std::size_t g = 0; g != component.groups.size(); ++g)
      {
        if (component.groups[g].rigid)
        {
          GroupState state = component.deriveGroupState(g, atoms);
          Integrators::initializeGroupVelocity(random, state, component.groups[g], system.temperature);
          states.push_back(state);
        }
      }
      const std::size_t offset = groupDataOffsetOfMolecule(system, selectedComponent, selectedMolecule);
      system.groupData.insert(system.groupData.begin() + static_cast<std::ptrdiff_t>(offset), states.begin(),
                              states.end());

      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (component.rigidGroupContaining(i).has_value())
        {
          dynamics[i].velocity = double3(0.0, 0.0, 0.0);
          continue;
        }
        const double mass = component.definedAtoms[i].second;
        dynamics[i].velocity = double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                               std::sqrt(Units::KB * system.temperature / mass);
      }
    }
    else
    {
      Integrators::initializeMoleculeVelocity(random, molecule, dynamics, component, system.temperature);
    }
  }

  const bool flexibleFrameworkConstraint =
      system.framework && !system.framework->rigid &&
      system.numberOfFrameworkAtoms + system.spanOfMoleculeAtoms().size() > 1uz;
  const bool moleculeConstraint = !system.framework.has_value() && system.numberOfMolecules() > 1uz;
  system.translationalCenterOfMassConstraint =
      flexibleFrameworkConstraint || moleculeConstraint ? std::min<std::size_t>(3, system.translationalDegreesOfFreedom)
                                                       : 0;

  if (system.thermostat)
    system.thermostat->refreshDegreesOfFreedom(system.translationalDegreesOfFreedom,
                                                system.rotationalDegreesOfFreedom,
                                                system.translationalCenterOfMassConstraint);
  if (system.thermobarostat)
  {
    const std::size_t effectiveDegreesOfFreedom =
        system.translationalDegreesOfFreedom - system.translationalCenterOfMassConstraint;
    system.thermobarostat->refreshDegreesOfFreedom(random, effectiveDegreesOfFreedom, system.simulationBox.volume);
  }

  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components,
                                                        system.spanOfGroupData());
}

void attemptParticleExchange(RandomNumber& random, System& system)
{
  if (!molecularDynamicsHasParticleExchange(system.molecularDynamicsEnsemble) || system.components.empty()) return;

  const std::size_t selectedComponent =
      std::min(static_cast<std::size_t>(random.uniform() * static_cast<double>(system.components.size())),
               system.components.size() - 1);
  if (system.components[selectedComponent].mc_moves_probabilities.getProbability(Move::Types::SwapCBMC) <= 0.0)
    return;

  std::size_t exchangedMolecule{};
  const MC_Moves::ParticleExchangeResult result =
      MC_Moves::performMolecularDynamicsSwap(random, system, selectedComponent, exchangedMolecule);
  refreshParticleNumberDependentState(random, system, selectedComponent, result, exchangedMolecule);
}

RunningEnergy thermobarostatVelocityVerlet(System& system)
{
  Thermobarostat& barostat = *system.thermobarostat;
  const auto pressureBefore = system.computeMolecularPressure();
  const double3x3 kineticBefore = computeMolecularKineticVirial(system);

  const double barostatKinetic =
      molecularDynamicsUsesIsotropicBarostat(barostat.ensemble)
          ? 0.5 * barostat.logVolumeMass * barostat.logVolumeVelocity * barostat.logVolumeVelocity
          : 0.5 * barostat.cellMass *
                std::transform_reduce(&barostat.cellVelocity.m[0], &barostat.cellVelocity.m[16], 0.0, std::plus<>(),
                                      [](double value) { return value * value; });
  const double chainScale = barostat.chainStep(barostatKinetic);
  barostat.logVolumeVelocity *= chainScale;
  barostat.cellVelocity = barostat.cellVelocity * chainScale;

  if (system.thermostat)
  {
    const auto scaling = system.thermostat->NoseHooverNVT(
        Integrators::computeTranslationalKineticEnergy(
            system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
            system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
            system.spanOfGroupData()),
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components,
                                                    system.spanOfGroupData()));
    Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                 system.components, scaling, system.framework, system.spanOfFrameworkDynamics(),
                                 system.spanOfGroupData());
  }

  double3x3 cellRate{};
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
  {
    const double mtkFactor =
        1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    const double scalarForce = (pressureBefore.second.trace() + mtkFactor * kineticBefore.trace() -
                                3.0 * barostat.pressure * system.simulationBox.volume) /
                               barostat.logVolumeMass;
    barostat.logVolumeVelocity += 0.5 * system.timeStep * scalarForce;
    cellRate =
        double3x3(barostat.logVolumeVelocity / 3.0, barostat.logVolumeVelocity / 3.0, barostat.logVolumeVelocity / 3.0);
  }
  else
  {
    double3x3 correctedKinetic = kineticBefore;
    const double mtkCorrection =
        kineticBefore.trace() / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    correctedKinetic.ax += mtkCorrection;
    correctedKinetic.by += mtkCorrection;
    correctedKinetic.cz += mtkCorrection;
    const double3x3 acceleration =
        cellForce(pressureBefore.second, correctedKinetic, system.simulationBox.volume, barostat.pressure,
                  barostat.cellMass, barostat.cellType, barostat.monoclinicAngle);
    barostat.cellVelocity += 0.5 * system.timeStep * acceleration;
    barostat.cellVelocity = projectCellTensor(barostat.cellVelocity, barostat.cellType, barostat.monoclinicAngle);
    cellRate = barostat.cellVelocity;
  }

  applyVelocityMatrix(system,
                      velocityPropagator(cellRate, 0.5 * system.timeStep, system.translationalDegreesOfFreedom));
  Integrators::updateVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                system.components, system.timeStep, system.framework, system.spanOfFrameworkAtoms(),
                                system.spanOfFrameworkDynamics(), &system.forceField, system.spanOfGroupData());
  propagateCell(system, cellRate);
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
    barostat.logVolumePosition += system.timeStep * barostat.logVolumeVelocity;
  Integrators::noSquishFreeRotorOrderTwo(system.moleculeData, system.components, system.timeStep,
                                         system.spanOfGroupData());
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                        system.spanOfGroupData());
  RunningEnergy energies = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics());
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components,
                                                        system.spanOfGroupData());
  Integrators::updateVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                system.components, system.timeStep, system.framework, system.spanOfFrameworkAtoms(),
                                system.spanOfFrameworkDynamics(), &system.forceField, system.spanOfGroupData());
  applyVelocityMatrix(system,
                      velocityPropagator(cellRate, 0.5 * system.timeStep, system.translationalDegreesOfFreedom));

  energies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData());
  energies.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData());
  const auto pressureAfter = system.computeMolecularPressure();
  const double3x3 kineticAfter = computeMolecularKineticVirial(system);
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
  {
    const double mtkFactor =
        1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    const double scalarForce = (pressureAfter.second.trace() + mtkFactor * kineticAfter.trace() -
                                3.0 * barostat.pressure * system.simulationBox.volume) /
                               barostat.logVolumeMass;
    barostat.logVolumeVelocity += 0.5 * system.timeStep * scalarForce;
  }
  else
  {
    double3x3 correctedKinetic = kineticAfter;
    const double mtkCorrection =
        kineticAfter.trace() / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    correctedKinetic.ax += mtkCorrection;
    correctedKinetic.by += mtkCorrection;
    correctedKinetic.cz += mtkCorrection;
    barostat.cellVelocity +=
        0.5 * system.timeStep *
        cellForce(pressureAfter.second, correctedKinetic, system.simulationBox.volume, barostat.pressure,
                  barostat.cellMass, barostat.cellType, barostat.monoclinicAngle);
    barostat.cellVelocity = projectCellTensor(barostat.cellVelocity, barostat.cellType, barostat.monoclinicAngle);
  }

  if (system.thermostat)
  {
    const auto scaling =
        system.thermostat->NoseHooverNVT(energies.translationalKineticEnergy, energies.rotationalKineticEnergy);
    Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                 system.components, scaling, system.framework, system.spanOfFrameworkDynamics(),
                                 system.spanOfGroupData());
    energies.NoseHooverEnergy = system.thermostat->getEnergy();
  }
  const double finalBarostatKinetic =
      molecularDynamicsUsesIsotropicBarostat(barostat.ensemble)
          ? 0.5 * barostat.logVolumeMass * barostat.logVolumeVelocity * barostat.logVolumeVelocity
          : 0.5 * barostat.cellMass *
                std::transform_reduce(&barostat.cellVelocity.m[0], &barostat.cellVelocity.m[16], 0.0, std::plus<>(),
                                      [](double value) { return value * value; });
  const double finalScale = barostat.chainStep(finalBarostatKinetic);
  barostat.logVolumeVelocity *= finalScale;
  barostat.cellVelocity = barostat.cellVelocity * finalScale;
  energies.thermobarostatEnergy = barostat.energy(system.simulationBox.volume);
  return energies;
}
}  // namespace

RunningEnergy molecularDynamicsStep(System& system)
{
  if (system.thermobarostat) return thermobarostatVelocityVerlet(system);
  return Integrators::velocityVerlet(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      system.spanOfFrameworkDynamics(), system.spanOfGroupData());
}

MolecularDynamics::MolecularDynamics() : random(std::nullopt) {};

MolecularDynamics::MolecularDynamics(InputReader& reader) noexcept
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      systems(std::move(reader.systems)),
      outputJsons(systems.size()),
      estimation(reader.numberOfBlocks, reader.numberOfProductionCycles)
{
}

MolecularDynamics::MolecularDynamics(const SimulationSchedule& schedule, const std::vector<System>& systems,
                                     std::optional<std::size_t> randomSeed, std::size_t numberOfBlocks,
                                     bool outputToFiles)
    : outputToFiles(outputToFiles),
      random(RandomNumber(randomSeed)),
      numberOfProductionCycles(schedule.numberOfProductionCycles),
      numberOfPreInitializationCycles(schedule.numberOfPreInitializationCycles),
      numberOfInitializationCycles(schedule.numberOfInitializationCycles),
      numberOfEquilibrationCycles(schedule.numberOfEquilibrationCycles),
      printEvery(schedule.printEvery),
      writeRestartEvery(5000),
      writeBinaryRestartEvery(schedule.writeBinaryRestartEvery),
      rescaleWangLandauEvery(schedule.rescaleWangLandauEvery),
      optimizeMCMovesEvery(schedule.optimizeMCMovesEvery),
      systems(systems),
      outputJsons(systems.size()),
      estimation(numberOfBlocks, schedule.numberOfProductionCycles)
{
}

System& MolecularDynamics::randomSystem()
{
  return systems[std::size_t(random.uniform() * static_cast<double>(systems.size()))];
}

void MolecularDynamics::run()
{
  // Semi-flexible molecules (components with rigid/flexible 'Groups') are integrated with a per-group
  // rigid-body scheme. The fixed-cell (NVE / NVT) and constant-pressure (NPT / NPT-PR) ensembles couple
  // the barostat to the rigid-group centers of mass and the flexible atoms. In the grand-canonical MD
  // ensembles the particle-exchange step splices the per-group rigid-body states (GroupState) of the
  // inserted/deleted molecule in or out of System::groupData (see refreshParticleNumberDependentState).
  switch (simulationStage)
  {
    case SimulationStage::Uninitialized:
      setup();
      break;
    case SimulationStage::PreInitialization:
      goto continuePreInitializationStage;
    case SimulationStage::Initialization:
      goto continueInitializationStage;
    case SimulationStage::Equilibration:
      goto continueEquilibrationStage;
    case SimulationStage::Production:
      goto continueProductionStage;
    default:
      break;
  }

continuePreInitializationStage:
  preInitialize();
continueInitializationStage:
  initialize();
continueEquilibrationStage:
  equilibrate();
continueProductionStage:
  production();

  tearDown();
}

void MolecularDynamics::createOutputFiles()
{
  std::filesystem::create_directories("output");
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.txt", system.temperature, system.input_pressure, system_id);
    streams.emplace_back(fileNameString, std::ios::out);
    ++system_id;
  }
}

void MolecularDynamics::createInterpolationGrids()
{
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    system.createExternalFieldInterpolationGrid(stream, system_id);
    system.createFrameworkInterpolationGrids(stream);

    ++system_id;
  }
}

void MolecularDynamics::setup()
{
  for (std::size_t system_id{0}; System& system : systems)
  {
    if (system.propertyElasticConstantsFluctuation)
    {
      if (system.molecularDynamicsEnsemble != MolecularDynamicsEnsemble::NVT)
        throw std::runtime_error("Stress-fluctuation elastic constants require fixed-cell NVT molecular dynamics");
      if (system.elasticConstantsSampleEvery == 0)
        throw std::runtime_error("ElasticConstantsSampleEvery must be positive");
      if (system.hasExternalField || system.forceField.computePolarization)
        throw std::runtime_error("Stress-fluctuation elastic constants do not support external fields or polarization");
    }
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    // switch the fractional molecule on in the first system, and off in all others
    if (system_id == 0uz)
      system.containsTheFractionalMolecule = true;
    else
      system.containsTheFractionalMolecule = false;

    // inactive fractional molecules must not contribute dUdlambda: set their groupId to zero
    system.initializeGibbsSwapFractionalMoleculeGroupIds();

    // if the MC/MD hybrid move is on, make sure that interpolation-method include gradients
    if (system.forceField.interpolationScheme == ForceField::InterpolationScheme::Polynomial)
    {
      system.forceField.interpolationScheme = ForceField::InterpolationScheme::Tricubic;
    }

    ++system_id;
  }

  if (outputToFiles)
  {
    createOutputFiles();

    for (std::size_t system_id{0}; const System& system : systems)
    {
      std::ostream stream(streams[system_id].rdbuf());

      std::print(stream, "{}", system.writeOutputHeader());
      std::print(stream, "Random seed: {}\n\n", random.seed);
      std::print(stream, "{}\n", HardwareInfo::writeInfo());
      std::print(stream, "{}", Units::printStatus());
      std::print(stream, "{}", system.writeSystemStatus());
      std::print(stream, "{}", system.forceField.printPseudoAtomStatus());
      std::print(stream, "{}", system.forceField.printForceFieldStatus());
      std::print(stream, "{}", system.writeComponentStatus());
      std::print(stream, "{}", system.reactions.printStatus());

      ++system_id;
    }
  }

  createInterpolationGrids();

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.initializeGroupData();
    system.precomputeTotalRigidEnergy();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData());
    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData());

    if (outputToFiles)
    {
      std::ostream stream(streams[system_id].rdbuf());
      stream << system.runningEnergies.printMC("Recomputed from scratch");
      std::print(stream, "\n\n\n\n");
    }

    ++system_id;
  };
}

void MolecularDynamics::tearDown()
{
  if (outputToFiles)
  {
    output();
  }
}

void MolecularDynamics::preInitialize(std::function<void()> call_back_function, std::size_t callBackEvery)
{
  std::size_t totalNumberOfMolecules{0uz};
  std::size_t totalNumberOfComponents{0uz};
  std::size_t numberOfStepsPerCycle{0uz};

  if (simulationStage == SimulationStage::PreInitialization) goto continuePreInitializationStage;
  simulationStage = SimulationStage::PreInitialization;

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();
  }

  for (currentCycle = 0uz; currentCycle != numberOfPreInitializationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    totalNumberOfMolecules = std::transform_reduce(
        systems.begin(), systems.end(), 0uz, [](const std::size_t& acc, const std::size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules(); });
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

    for (std::size_t j = 0uz; j != numberOfStepsPerCycle; j++)
    {
      std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      std::size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMovePreInitialization(random, selectedSystem, selectSecondSystem, selectedComponent,
                                                   fractionalMoleculeSystem);

      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
        }
      }
    }

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system_id].rdbuf());

          std::print(stream, "{}",
                     system.writePreInitializationStatusReport(currentCycle, numberOfPreInitializationCycles));
          std::print(stream, "{}\n\n\n\n", system.runningEnergies.printMC(""));
          std::flush(stream);
        }

        ++system_id;
      }
    }

    if (currentCycle % callBackEvery == 0uz)
    {
      if (call_back_function)
      {
        call_back_function();
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
        Archive<std::ofstream> archive(ofile);
        archive << *this;
        ofile.close();
        if (ofile)
        {
          std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
        }
      }
    }

  continuePreInitializationStage:;
  }
}

void MolecularDynamics::initialize(std::function<void()> call_back_function, std::size_t callBackEvery)
{
  std::size_t totalNumberOfMolecules{0uz};
  std::size_t totalNumberOfComponents{0uz};
  std::size_t numberOfStepsPerCycle{0uz};

  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  for (System& system : systems)
  {
    system.precomputeTotalRigidEnergy();
    system.runningEnergies = system.computeTotalEnergies();
  }

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    totalNumberOfMolecules = std::transform_reduce(
        systems.begin(), systems.end(), 0uz, [](const std::size_t& acc, const std::size_t& b) { return acc + b; },
        [](const System& system) { return system.numberOfMolecules(); });
    totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

    numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

    for (std::size_t j = 0uz; j != numberOfStepsPerCycle; j++)
    {
      // move to 'slide' when implemented in llvm
      // [[maybe_unused]] auto s = std::ranges::views::iota(0uz, systems.size());
      // std::ranges::views::slide(s, 2uz);

      std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      std::size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem);

      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
        }
      }
    }

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system_id].rdbuf());

          std::print(stream, "{}", system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
          std::print(stream, "{}\n\n\n\n", system.runningEnergies.printMC(""));
          std::flush(stream);
        }

        ++system_id;
      }
    }

    if (currentCycle % callBackEvery == 0uz)
    {
      if (call_back_function)
      {
        call_back_function();
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
        Archive<std::ofstream> archive(ofile);
        archive << *this;
        ofile.close();
        if (ofile)
        {
          std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
        }
      }
    }

  continueInitializationStage:;
  }
}

void MolecularDynamics::equilibrate(std::function<void()> call_back_function, std::size_t callBackEvery)
{
  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.initializeGroupData();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData());
    Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                      system.spanOfMoleculeDynamics(), system.components, system.temperature,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(),
                                      &system.forceField, system.spanOfGroupData());

    Integrators::removeCenterOfMassVelocityDrift(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData());
    if (system.thermostat.has_value())
    {
      const bool flexibleFrameworkConstraint =
          system.framework && !system.framework->rigid &&
          system.numberOfFrameworkAtoms + system.spanOfMoleculeAtoms().size() > 1uz;
      if (flexibleFrameworkConstraint || (!system.framework.has_value() && system.numberOfMolecules() > 1uz))
      {
        system.translationalCenterOfMassConstraint = 3;
        system.thermostat->translationalCenterOfMassConstraint = 3;
      }
      system.thermostat->initialize(random);
    }
    if (system.thermobarostat.has_value())
    {
      system.thermobarostat->translationalDegreesOfFreedom =
          system.translationalDegreesOfFreedom - system.translationalCenterOfMassConstraint;
      system.thermobarostat->logVolumePosition = std::log(system.simulationBox.volume);
      system.thermobarostat->initialize(random);
    }

    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData());
    if (system.thermostat.has_value())
    {
      system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
    }
    if (system.thermobarostat.has_value())
      system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);
    system.referenceEnergy = system.runningEnergies.conservedEnergy();

    if (outputToFiles)
    {
      std::ostream stream(streams[system_id].rdbuf());
      stream << system.runningEnergies.printMD("Recomputed from scratch", system.referenceEnergy);
      std::print(stream, "\n\n\n\n");
    }

    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }

    ++system_id;
  };

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    for (System& system : systems)
    {
      if (currentCycle % 3uz == 0uz) attemptParticleExchange(random, system);
      system.runningEnergies = molecularDynamicsStep(system);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift +=
          std::abs((system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);

      if (system.propertyElasticConstantsFluctuation && currentCycle % system.elasticConstantsSampleEvery == 0)
      {
        System samplingSystem = system;
        const auto molecularPressure = samplingSystem.computeMolecularPressure();
        system.currentEnergyStatus = molecularPressure.first;
        system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
        const double3x3 kineticPressure = computeMolecularKineticVirial(system) / system.simulationBox.volume;
        const double effectiveKineticEntities =
            static_cast<double>(system.translationalDegreesOfFreedom - system.translationalCenterOfMassConstraint) /
            3.0;
        const ElasticFluctuationTerms terms(stressVoigt(system.currentExcessPressureTensor),
                                            stressVoigt(kineticPressure), computeAffineBornTensor(system), system.beta,
                                            system.simulationBox.volume, effectiveKineticEntities);
        system.propertyElasticConstantsFluctuation->addSample(estimation.currentBin, terms, 1.0);
      }
    }

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        system.loadings =
            LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);

        if (outputToFiles)
        {
          std::ostream stream(streams[system_id].rdbuf());

          std::print(stream, "{}", system.writeEquilibrationStatusReportMD(currentCycle, numberOfEquilibrationCycles));
          std::flush(stream);
        }

        ++system_id;
      }
    }

    if (currentCycle % callBackEvery == 0uz)
    {
      if (call_back_function)
      {
        call_back_function();
      }
    }

    if (currentCycle % rescaleWangLandauEvery == 0uz)
    {
      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.WangLandauIteration(
              PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
              system.containsTheFractionalMolecule);
        }
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % printEvery == 0uz)
    {
      // write restart
      if (outputToFiles)
      {
        std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
        Archive<std::ofstream> archive(ofile);
        archive << *this;
        ofile.close();
        if (ofile)
        {
          std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
        }
      }
    }
  continueEquilibrationStage:;
  }
}

void MolecularDynamics::production(std::function<void()> call_back_function, std::size_t callBackEvery)
{
  std::chrono::steady_clock::time_point t1, t2;
  double minBias{0.0};

  if (simulationStage == SimulationStage::Production) goto continueProductionStage;
  simulationStage = SimulationStage::Production;

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.initializeGroupData();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData());
    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData());
    if (system.thermostat.has_value())
    {
      system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
    }
    if (system.thermobarostat.has_value())
      system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);
    system.referenceEnergy = system.runningEnergies.conservedEnergy();

    if (outputToFiles)
    {
      std::ostream stream(streams[system_id].rdbuf());
      stream << system.runningEnergies.printMD("Recomputed from scratch", system.referenceEnergy);
      std::print(stream, "\n");
    }

    system.mc_moves_statistics.clearMoveStatistics();
    system.mc_moves_cputime.clearTimingStatistics();
    integratorsCPUTime.clearTimingStatistics();

    system.accumulatedDrift = 0.0;

    for (Component& component : system.components)
    {
      component.mc_moves_statistics.clearMoveStatistics();
      component.mc_moves_cputime.clearTimingStatistics();

      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
    }

    ++system_id;
  };

  minBias = std::numeric_limits<double>::max();
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      double currentMinBias =
          *std::min_element(component.lambdaGC.biasFactor.cbegin(), component.lambdaGC.biasFactor.cend());
      minBias = currentMinBias < minBias ? currentMinBias : minBias;
    }
  }
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      component.lambdaGC.normalize(minBias);
    }
  }

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfProductionCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    estimation.setCurrentSample(currentCycle);

    for ([[maybe_unused]] System& system : systems)
    {
      // add the sample energy to the averages
      if (currentCycle % 10uz == 0uz || currentCycle % printEvery == 0uz)
      {
        std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
        std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
        system.currentEnergyStatus = molecularPressure.first;
        system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
        std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
        system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
        system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
      }
    }

    for (System& system : systems)
    {
      if (currentCycle % 3uz == 0uz) attemptParticleExchange(random, system);
      system.runningEnergies = molecularDynamicsStep(system);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift +=
          std::abs((system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
    }

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    // sample properties
    for (std::size_t system_id{0}; System& system : systems)
    {
      system.sampleProperties(system_id, estimation.currentBin, currentCycle);
      // MD forces from the integrator step above are current; reuse them for the force-based RDF.
      if (system.forceBasedRDFSampleDue(currentCycle))
      {
        system.sampleForceBasedRDFFromCurrentGradients(currentCycle, estimation.currentBin);
      }

      ++system_id;
    }

    if (currentCycle % printEvery == 0uz)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        if (outputToFiles)
        {
          std::ostream stream(streams[system_id].rdbuf());
          std::print(stream, "{}", system.writeProductionStatusReportMD(currentCycle, numberOfProductionCycles));
          std::flush(stream);
        }

        ++system_id;
      }
    }

    if (currentCycle % callBackEvery == 0uz)
    {
      if (call_back_function)
      {
        call_back_function();
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    // output properties to files
    for (std::size_t system_id{0}; System& system : systems)
    {
      if (system.propertyConventionalRadialDistributionFunction.has_value())
      {
        system.propertyConventionalRadialDistributionFunction->writeOutput(
            system.forceField, system_id, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, currentCycle);
      }

      if (system.propertyRadialDistributionFunction.has_value())
      {
        system.propertyRadialDistributionFunction->writeOutput(
            system.forceField, system_id, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, currentCycle);
      }
      if (system.propertyDensityGrid.has_value())
      {
        system.propertyDensityGrid->writeOutput(system_id, system.simulationBox, system.forceField, system.framework,
                                                system.components, currentCycle);
      }

      if (system.propertyMSD.has_value())
      {
        system.propertyMSD->writeOutput(system_id, system.components, currentCycle);
      }

      if (system.propertyVACF.has_value())
      {
        system.propertyVACF->writeOutput(system_id, system.components, currentCycle);
      }

      if (system.propertyMoleculeProperties.has_value())
      {
        system.propertyMoleculeProperties->writeOutput(system_id, system.components, currentCycle);
      }

      ++system_id;
    }

    // write binary-restart file
    if (currentCycle % printEvery == 0uz)
    {
      if (outputToFiles)
      {
        std::ofstream ofile("restart_data.bin_temp", std::ios::binary);
        Archive<std::ofstream> archive(ofile);
        archive << *this;
        ofile.close();
        if (ofile)
        {
          std::filesystem::rename("restart_data.bin_temp", "restart_data.bin");
        }
      }
    }
    t2 = std::chrono::steady_clock::now();
    totalSimulationTime += (t2 - t1);
  continueProductionStage:;
  }

  // output properties to files
  for (std::size_t system_id{0}; System& system : systems)
  {
    if (system.propertyConventionalRadialDistributionFunction.has_value())
    {
      system.propertyConventionalRadialDistributionFunction->writeOutput(
          system.forceField, system_id, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, currentCycle);
    }

    if (system.propertyRadialDistributionFunction.has_value())
    {
      system.propertyRadialDistributionFunction->writeOutput(system.forceField, system_id, system.simulationBox.volume,
                                                             system.totalNumberOfPseudoAtoms, currentCycle);
    }
    if (system.propertyDensityGrid.has_value())
    {
      system.propertyDensityGrid->writeOutput(system_id, system.simulationBox, system.forceField, system.framework,
                                              system.components, currentCycle);
    }

    if (system.propertyMSD.has_value())
    {
      system.propertyMSD->writeOutput(system_id, system.components, currentCycle);
    }

    if (system.propertyVACF.has_value())
    {
      system.propertyVACF->writeOutput(system_id, system.components, currentCycle);
    }

    if (system.propertyMoleculeProperties.has_value())
    {
      system.propertyMoleculeProperties->writeOutput(system_id, system.components, currentCycle);
    }

    ++system_id;
  }
}

void MolecularDynamics::output()
{
  MCMoveCpuTime total;
  MCMoveStatistics countTotal;
  for (const System& system : systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_statistics;
    for (const Component& component : system.components)
    {
      countTotal += component.mc_moves_statistics;
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    std::print(stream, "\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "                             Simulation finished!\n");
    std::print(stream, "===============================================================================\n");
    std::print(stream, "\n");

    std::string status_line{std::format("Final state after {} cycles\n\n", numberOfProductionCycles)};

    std::print(stream, "Production run CPU timings of the MD simulation\n");
    std::print(stream, "===============================================================================\n\n");

    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(stream, "{}", integratorsCPUTime.writeIntegratorsCPUTimeStatistics(totalSimulationTime));
    std::print(stream, "\n\n");

    std::print(
        stream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));

    std::print(stream, "Temperature averages and statistics:\n");
    std::print(stream, "===============================================================================\n\n");
    std::print(stream, "{}", system.averageTemperature.writeAveragesStatistics("Total"));
    std::print(stream, "{}", system.averageTranslationalTemperature.writeAveragesStatistics("Translational"));
    std::print(stream, "{}", system.averageRotationalTemperature.writeAveragesStatistics("Rotational"));

    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    if (system.propertyElasticConstantsFluctuation)
    {
      std::print(stream, "{}", system.propertyElasticConstantsFluctuation->writeAveragesStatistics());
      outputJsons[system_id]["properties"]["elasticConstantsFromStressFluctuations"] =
          system.propertyElasticConstantsFluctuation->jsonAveragesStatistics();
    }

    std::print(
        stream, "{}",
        system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swappableComponents, system.components));
    std::print(
        stream, "{}",
        system.averagePartialMolarProperties.writeAveragesStatistics(system.swappableComponents, system.components));
    std::print(stream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));

    if (system.propertyElasticConstantsFluctuation)
    {
      const std::string jsonFileName =
          std::format("output/output_{}_{}.s{}.json", system.temperature, system.input_pressure, system_id);
      std::ofstream json(jsonFileName);
      json << outputJsons[system_id].dump(4);
    }

    ++system_id;
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MolecularDynamics& mc)
{
  archive << mc.versionNumber;

  archive << mc.outputToFiles;
  archive << mc.random;

  archive << mc.numberOfProductionCycles;
  archive << mc.numberOfSteps;
  archive << mc.numberOfPreInitializationCycles;
  archive << mc.numberOfInitializationCycles;
  archive << mc.numberOfEquilibrationCycles;

  archive << mc.printEvery;
  archive << mc.writeRestartEvery;
  archive << mc.writeBinaryRestartEvery;
  archive << mc.rescaleWangLandauEvery;
  archive << mc.optimizeMCMovesEvery;

  archive << mc.currentCycle;
  archive << mc.absoluteCurrentCycle;
  archive << mc.simulationStage;

  archive << mc.systems;
  archive << mc.fractionalMoleculeSystem;

  archive << mc.estimation;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MolecularDynamics& mc)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > mc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MolecularDynamics' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> mc.outputToFiles;
  archive >> mc.random;

  archive >> mc.numberOfProductionCycles;
  archive >> mc.numberOfSteps;
  archive >> mc.numberOfPreInitializationCycles;
  archive >> mc.numberOfInitializationCycles;
  archive >> mc.numberOfEquilibrationCycles;

  archive >> mc.printEvery;
  archive >> mc.writeRestartEvery;
  archive >> mc.writeBinaryRestartEvery;
  archive >> mc.rescaleWangLandauEvery;
  archive >> mc.optimizeMCMovesEvery;

  archive >> mc.currentCycle;
  archive >> mc.absoluteCurrentCycle;
  archive >> mc.simulationStage;

  archive >> mc.systems;
  archive >> mc.fractionalMoleculeSystem;

  archive >> mc.estimation;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
  }
  std::cout << std::format("Magic number read correctly: {} vs {}\n", magicNumber,
                           static_cast<std::uint64_t>(0x6f6b6179));
  return archive;
}
