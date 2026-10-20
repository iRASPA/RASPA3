module;

module dctst;

import std;

import stringutils;
import hardware_info;
import system;
import framework;
import randomnumbers;
import input_reader;
import component;
import averages;
import units;
import simulationbox;
import forcefield;
import running_energy;
import atom;
import atom_dynamics;
import molecule;
import fragment;
import int3;
import double3;
import double3x3;
import simd_quatd;
import cbmc;
import cbmc_constants;
import mc_moves;
import mc_moves_move_types;
import mc_moves_translation;
import mc_moves_random_translation;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_partial_reinsertion;
import integrators;
import integrators_compute;
import integrators_update;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import thermostat;
import json;
import graceful_shutdown;
import dctst_settings;
import dctst_barrier;

// =================================================================================================
// Accumulators
// =================================================================================================

DcTSTTransmissionAccumulator::DcTSTTransmissionAccumulator(std::size_t numberOfBlocks, std::size_t numberOfShells,
                                                           std::size_t numberOfComponents, std::size_t numberOfTimes)
    : numberOfBlocks(numberOfBlocks),
      numberOfShells(numberOfShells),
      numberOfComponents(numberOfComponents),
      numberOfTimes(numberOfTimes),
      speed(numberOfBlocks, 0.0),
      self(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      selfTwoState(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      selfCommitted(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      collective(numberOfBlocks, std::vector<double>(numberOfShells * numberOfComponents * numberOfTimes, 0.0)),
      numberOfTrajectoryPairs(numberOfBlocks, 0uz),
      selfEvent(numberOfBlocks, 0.0),
      collectiveEvent(numberOfBlocks, 0.0),
      speedEventTerminated(numberOfBlocks, 0.0),
      selfEventConditioned(numberOfBlocks, 0.0),
      collectiveEventConditioned(numberOfBlocks, 0.0),
      collectiveEventBeforeCommitment(numberOfBlocks, 0.0),
      collectiveEventAfterCommitment(numberOfBlocks, 0.0)
{
}

void DcTSTConstrainedMoveStatistics::optimize(double targetAcceptance)
{
  if (translationTrials > 0)
  {
    const double acceptance = static_cast<double>(translationAccepted) / static_cast<double>(translationTrials);
    const double scaling = std::clamp(acceptance / targetAcceptance, 0.5, 1.5);
    maximumTranslation = std::clamp(maximumTranslation * scaling, 0.01, 10.0);
  }
  translationTotalTrials += translationTrials;
  translationTotalAccepted += translationAccepted;
  translationTrials = 0;
  translationAccepted = 0;

  reinsertionTotalTrials += reinsertionTrials;
  reinsertionTotalAccepted += reinsertionAccepted;
  reinsertionTrials = 0;
  reinsertionAccepted = 0;
}

// =================================================================================================
// Helpers
// =================================================================================================

namespace
{
/// The velocity of atom 'bead' of a molecule: for a rigid molecule the center-of-mass velocity plus the
/// rotational part omega x r (omega from the orientation momentum, in the same convention as the kinetic
/// energy of the integrators), for a flexible molecule the atomic velocity.
double3 beadVelocity(const Component& component, const Molecule& molecule, std::span<const AtomDynamics> dynamics,
                     std::size_t bead)
{
  if (component.rigid)
  {
    if (molecule.numberOfAtoms == 1) return molecule.velocity;

    simd_quatd p = molecule.orientationMomentum;
    p.r = -p.r;
    const simd_quatd pq = p * molecule.orientation;
    const double3 omegaBody = 0.5 * double3(pq.ix, pq.iy, pq.iz) * component.inverseInertiaVector;
    const double3 rBody = component.atoms[bead].position;
    return molecule.velocity + molecule.orientation * double3::cross(omegaBody, rBody);
  }
  return dynamics[molecule.atomIndex + bead].velocity;
}

/// Reverses all momenta (time reversal of the trajectory).
void reverseVelocities(std::span<Molecule> moleculeData, std::span<AtomDynamics> dynamics,
                       std::span<GroupState> groupData)
{
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity = -molecule.velocity;
    molecule.orientationMomentum = simd_quatd(-molecule.orientationMomentum.r, -molecule.orientationMomentum.ix,
                                              -molecule.orientationMomentum.iy, -molecule.orientationMomentum.iz);
  }
  for (AtomDynamics& dynamic : dynamics)
  {
    dynamic.velocity = -dynamic.velocity;
  }
  for (GroupState& state : groupData)
  {
    state.velocity = -state.velocity;
    state.orientationMomentum = simd_quatd(-state.orientationMomentum.r, -state.orientationMomentum.ix,
                                           -state.orientationMomentum.iy, -state.orientationMomentum.iz);
  }
}

/// Body-frame angular velocity from the orientation momentum (convention of the integrators).
double3 bodyAngularVelocity(const simd_quatd& orientation, const simd_quatd& orientationMomentum,
                            const double3& inverseInertia)
{
  simd_quatd p = orientationMomentum;
  p.r = -p.r;
  const simd_quatd pq = p * orientation;
  return 0.5 * double3(pq.ix, pq.iy, pq.iz) * inverseInertia;
}

/// Orientation momentum from the body-frame angular velocity (inverse of 'bodyAngularVelocity', the
/// same expression the integrators use to draw initial velocities).
simd_quatd orientationMomentumFromBodyAngularVelocity(const simd_quatd& q, const double3& I, const double3& w)
{
  simd_quatd p(0.0, 0.0, 0.0, 0.0);
  p.ix = 2.0 * (q.r * (I.x * w.x) - q.iz * (I.y * w.y) + q.iy * (I.z * w.z));
  p.iy = 2.0 * (q.iz * (I.x * w.x) + q.r * (I.y * w.y) - q.ix * (I.z * w.z));
  p.iz = 2.0 * (-q.iy * (I.x * w.x) + q.ix * (I.y * w.y) + q.r * (I.z * w.z));
  p.r = 2.0 * (-q.ix * (I.x * w.x) - q.iy * (I.y * w.y) - q.iz * (I.z * w.z));
  return p;
}

/**
 * \brief One step of the stochastic bath on all molecules (Andersen or Langevin).
 *
 * Andersen: every molecule is, with probability 1 - exp(-rate dt), given a fresh Maxwell-Boltzmann
 * velocity (center of mass and angular momentum of rigid molecules and rigid fragments, atomic velocities
 * of flexible atoms). Between collisions the dynamics is exactly Hamiltonian.
 * Langevin: an exact Ornstein-Uhlenbeck step v -> c v + sqrt(1 - c^2) sqrt(kT/m) xi with c = exp(-rate dt)
 * on the same degrees of freedom (angular velocities in the body frame with the principal inertia).
 */
void applyBath(RandomNumber& rng, DcTSTSettings::TrajectoryThermostat type, double rate, double timeStep,
               double temperature, std::span<Molecule> moleculeData, std::span<AtomDynamics> dynamics,
               std::span<GroupState> groupData, const std::vector<Component>& components)
{
  const double kT = Units::KB * temperature;
  const double c1 = std::exp(-rate * timeStep);
  const double c2 = std::sqrt(std::max(0.0, 1.0 - c1 * c1));
  const double collisionProbability = 1.0 - c1;

  auto gaussian3 = [&]() { return double3(rng.Gaussian(), rng.Gaussian(), rng.Gaussian()); };

  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    std::span<AtomDynamics> moleculeDynamics = dynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms);
    const std::size_t numberOfRigidFragments = component.isSemiFlexible() ? component.numberOfRigidFragments() : 0uz;

    if (type == DcTSTSettings::TrajectoryThermostat::Andersen)
    {
      if (rng.uniform() < collisionProbability)
      {
        if (component.isSemiFlexible())
        {
          std::size_t rigidRank{};
          for (const Fragment& group : component.fragmentGraph.fragments)
          {
            if (!group.isRigidBody()) continue;
            Integrators::initializeGroupVelocity(rng, groupData[groupIndex + rigidRank], group, temperature);
            ++rigidRank;
          }
          for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
          {
            if (component.rigidFragmentContaining(i).has_value()) continue;
            moleculeDynamics[i].velocity = gaussian3() * std::sqrt(kT / component.definedAtoms[i].second);
          }
        }
        else
        {
          Integrators::initializeMoleculeVelocity(rng, molecule, moleculeDynamics, component, temperature);
        }
      }
    }
    else if (type == DcTSTSettings::TrajectoryThermostat::Langevin)
    {
      if (component.isSemiFlexible())
      {
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (!group.isRigidBody()) continue;
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity = c1 * state.velocity + c2 * std::sqrt(kT / group.mass) * gaussian3();
          const double3 omega =
              c1 * bodyAngularVelocity(state.orientation, state.orientationMomentum, group.inverseInertiaVector) +
              c2 * sqrt(group.inverseInertiaVector * kT) * gaussian3();
          state.orientationMomentum =
              orientationMomentumFromBodyAngularVelocity(state.orientation, group.inertiaVector, omega);
          ++rigidRank;
        }
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          if (component.rigidFragmentContaining(i).has_value()) continue;
          moleculeDynamics[i].velocity =
              c1 * moleculeDynamics[i].velocity + c2 * std::sqrt(kT / component.definedAtoms[i].second) * gaussian3();
        }
      }
      else if (component.rigid)
      {
        molecule.velocity = c1 * molecule.velocity + c2 * std::sqrt(kT / molecule.mass) * gaussian3();
        if (molecule.numberOfAtoms > 1)
        {
          const double3 omega = c1 * bodyAngularVelocity(molecule.orientation, molecule.orientationMomentum,
                                                         component.inverseInertiaVector) +
                                c2 * sqrt(component.inverseInertiaVector * kT) * gaussian3();
          molecule.orientationMomentum =
              orientationMomentumFromBodyAngularVelocity(molecule.orientation, component.inertiaVector, omega);
        }
      }
      else
      {
        double3 momentum{};
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          const double mass = component.definedAtoms[i].second;
          moleculeDynamics[i].velocity = c1 * moleculeDynamics[i].velocity + c2 * std::sqrt(kT / mass) * gaussian3();
          momentum += mass * moleculeDynamics[i].velocity;
        }
        molecule.velocity = momentum * molecule.invMass;
      }
    }
    groupIndex += numberOfRigidFragments;
  }
}

/// Mean and standard error of the mean of a set of block values.
std::pair<double, double> meanAndError(std::span<const double> values)
{
  if (values.empty()) return {0.0, 0.0};
  const double n = static_cast<double>(values.size());
  const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
  if (values.size() < 2) return {mean, 0.0};
  double sumSquares{0.0};
  for (double v : values) sumSquares += (v - mean) * (v - mean);
  // 95% confidence interval of the mean (Student t for the usual 5 blocks is ~2.78; use 2 for simplicity
  // consistent with the other RASPA block averages which report 2 sigma)
  return {mean, 2.0 * std::sqrt(sumSquares / (n - 1.0) / n)};
}
}  // namespace

// =================================================================================================
// Construction and setup
// =================================================================================================

DynamicallyCorrectedTST::DynamicallyCorrectedTST(InputReader& reader)
    : random(reader.randomSeed),
      settings(reader.dctstSettings),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      numberOfBlocks(reader.numberOfBlocks),
      numberOfReplicas(std::max(1uz, reader.dctstSettings.numberOfReplicas))
{
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  barrier = DcTSTBarrier(templateSystem.simulationBox, settings.barrierPosition, settings.barrierNormal,
                         settings.latticeDistance);

  // the cage lattice of the occupancy-resolved and event-terminated modes: one cage per unit cell. Its origin
  // is 'CageCenter' (unit-cell fractional) when given; otherwise the barrier point shifted half a lattice
  // distance along the normal (the center of cage B), which is only correct when the in-plane position of
  // 'BarrierPosition' (a supercell fraction) falls on a cage center
  double3 lattice_origin = barrier.point + 0.5 * barrier.lambda * barrier.normal;
  if (settings.cageCenter.has_value())
  {
    if (!templateSystem.framework.has_value())
    {
      throw std::runtime_error("dcTST: 'CageCenter' requires a framework (it is a unit-cell fractional position)\n");
    }
    lattice_origin = templateSystem.framework->simulationBox.cell * settings.cageCenter.value();
    // the cage center must be midway between two barrier planes
    const double offset = std::abs(barrier.reducedCoordinate(barrier.projection(templateSystem.simulationBox, lattice_origin)));
    if (std::abs(offset - 0.5 * barrier.lambda) > 1e-3 * barrier.lambda)
    {
      throw std::runtime_error(
          std::format("dcTST: 'CageCenter' must lie half a 'LatticeDistance' from the barrier plane along the normal "
                      "(it lies {:.4f} [A] from the nearest plane, expected {:.4f})\n",
                      offset, 0.5 * barrier.lambda));
    }
  }
  if (templateSystem.framework.has_value())
  {
    cageLattice =
        DcTSTCageLattice(templateSystem.framework->simulationBox, templateSystem.framework->numberOfUnitCells, lattice_origin);
  }
  else
  {
    cageLattice = DcTSTCageLattice(templateSystem.simulationBox, int3(1, 1, 1), lattice_origin);
  }

  // the event-terminated estimator classifies the hops by lattice axis: the barrier normal must be a
  // lattice direction and the column of cages A', A, B, B' needs at least three cells along it
  if (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event)
  {
    const DcTSTCageLattice::AxisAlignment alignment = cageLattice.axisAlongDirection(barrier.normal);
    if (alignment.cosine < 0.999)
    {
      throw std::runtime_error(
          std::format("dcTST: 'TrajectoryTermination' \"Event\" requires the barrier normal along a lattice vector "
                      "of the framework unit cell (best alignment |cos| = {:.4f})\n",
                      alignment.cosine));
    }
    if (cageLattice.numberOfUnitCells[alignment.axis] < 3)
    {
      throw std::runtime_error(
          "dcTST: 'TrajectoryTermination' \"Event\" needs at least three unit cells along the barrier normal (the "
          "column of cages A', A, B, B')\n");
    }
    if (std::abs(cageLattice.axisLength(alignment.axis) - barrier.lambda) > 1e-3 * barrier.lambda)
    {
      throw std::runtime_error(
          std::format("dcTST: 'TrajectoryTermination' \"Event\" requires 'LatticeDistance' ({}) equal to the unit-cell "
                      "length along the normal ({})\n",
                      barrier.lambda, cageLattice.axisLength(alignment.axis)));
    }
  }

  // the trajectories are integrated with the system time step; the transmission functions are stored
  // every 'sampleTrajectoryEvery' steps
  settings.sampleTrajectoryEvery = std::max(1uz, settings.sampleTrajectoryEvery);
  numberOfTrajectorySteps =
      std::max(1uz, static_cast<std::size_t>(std::llround(settings.trajectoryTime / templateSystem.timeStep)));
  numberOfTimes = numberOfTrajectorySteps / settings.sampleTrajectoryEvery + 1uz;

  // the trajectories are NVE: no thermostat is used by this driver
  templateSystem.thermostat = std::nullopt;
  templateSystem.thermobarostat = std::nullopt;

  systems.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  randoms.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    randoms.emplace_back(random.seed + replicaId + 1);
  }

  const std::size_t numberOfShells = settings.maximumCollectiveShell + 2uz;
  const std::size_t numberOfComponents = systems.front().components.size();
  transmission.assign(numberOfReplicas,
                      DcTSTTransmissionAccumulator(numberOfBlocks, numberOfShells, numberOfComponents, numberOfTimes));
  constrainedMoves.assign(numberOfReplicas, DcTSTConstrainedMoveStatistics{});
  for (DcTSTConstrainedMoveStatistics& statistics : constrainedMoves)
  {
    statistics.maximumTranslation = settings.maximumTranslationInPlane;
  }
  occupancyStates.assign(numberOfReplicas, DcTSTOccupancyState{});
  completedTransmissionCycles.assign(numberOfReplicas, 0uz);
  lastSnapshotTime.assign(numberOfReplicas, std::chrono::steady_clock::time_point{});
}

void DynamicallyCorrectedTST::run()
{
  setup();
  loadFreeEnergyInput();
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Transmission,
           settings.numberOfConstrainedEquilibrationCycles + settings.numberOfTransmissionCycles);
  output();
}

void DynamicallyCorrectedTST::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    // the integrator needs gradients: a polynomial interpolation grid does not provide them
    if (system.forceField.settings.interpolationScheme == ForceFieldSettings::InterpolationScheme::Polynomial)
    {
      system.forceField.settings.interpolationScheme = ForceFieldSettings::InterpolationScheme::Tricubic;
    }
  }

  std::filesystem::create_directories("output");
  std::filesystem::create_directories("dctst");

  const System& front = systems.front();
  stream.open(std::format("output/output_{}_{}.dctst.txt", front.temperature, front.input_pressure), std::ios::out);
  outputJsonFileName = std::format("output/output_{}_{}.dctst.json", front.temperature, front.input_pressure);

  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  const Component& tagged = front.components[settings.taggedComponent];
  std::print(stream, "Dynamically corrected transition state theory (dcTST)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Tagged component:                        {} ({})\n", tagged.name, settings.taggedComponent);
  std::print(stream, "Reaction bead (starting bead):           {}\n", tagged.startingBead);
  std::print(stream, "Barrier position (fractional):           {} {} {}\n", settings.barrierPosition.x,
             settings.barrierPosition.y, settings.barrierPosition.z);
  std::print(stream, "Barrier point (Cartesian):               {:.5f} {:.5f} {:.5f} [A]\n", barrier.point.x,
             barrier.point.y, barrier.point.z);
  std::print(stream, "Barrier normal (Cartesian, unit):        {:.5f} {:.5f} {:.5f}\n", barrier.normal.x,
             barrier.normal.y, barrier.normal.z);
  std::print(stream, "Lattice (cage-to-cage) distance lambda:  {} [A]\n", barrier.lambda);
  std::print(stream, "Cages along the normal in the box:       {:.3f}\n",
             double3::dot(double3::abs(barrier.normal), front.simulationBox.perpendicularWidths()) / barrier.lambda);
  std::print(stream, "Number of replicas (threads):            {}\n", numberOfReplicas);
  std::print(stream, "Number of blocks per replica:            {}\n", numberOfBlocks);
  std::print(stream, "Initialization cycles (constrained):     {}\n", numberOfInitializationCycles);
  std::print(stream, "Constrained equilibration cycles:        {}\n", settings.numberOfConstrainedEquilibrationCycles);
  std::print(stream, "Transmission cycles (constrained):       {}\n", settings.numberOfTransmissionCycles);
  std::print(stream, "Trajectory pairs launched every:         {} cycles\n", settings.launchTrajectoryEvery);
  std::print(stream, "Trajectory pairs per configuration:      {}\n", settings.trajectoriesPerConfiguration);
  std::print(stream, "Trajectory time:                         {} [ps] ({} steps of {} [ps])\n",
             settings.trajectoryTime, numberOfTrajectorySteps, front.timeStep);
  std::print(stream, "Transmission functions stored every:     {} steps ({} points)\n", settings.sampleTrajectoryEvery,
             numberOfTimes);
  std::print(stream, "Collective shells accumulated:           0..{} and all\n", settings.maximumCollectiveShell);
  switch (settings.trajectoryThermostat)
  {
    case DcTSTSettings::TrajectoryThermostat::None:
      std::print(stream, "Trajectory bath:                         none (NVE; diagnostic mode at low loading)\n");
      break;
    case DcTSTSettings::TrajectoryThermostat::Andersen:
      std::print(stream, "Trajectory bath:                         Andersen, collision rate {} [1/ps] per molecule\n",
                 settings.bathRate);
      break;
    case DcTSTSettings::TrajectoryThermostat::Langevin:
      std::print(stream, "Trajectory bath:                         Langevin, friction {} [1/ps]\n", settings.bathRate);
      break;
  }
  std::print(stream, "Commit distance (committed kappa_S):     {} [A]\n",
             (settings.commitDistance > 0.0) ? settings.commitDistance : 0.5 * barrier.lambda);
  if (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event)
  {
    std::print(stream,
               "Trajectory termination:                  event-terminated, at the {}-th local axial hop event after "
               "commitment (cap {} [ps])\n",
               settings.terminationEvents, settings.trajectoryTime);
  }
  else
  {
    std::print(stream, "Trajectory termination:                  fixed time\n");
  }
  std::print(stream, "Maximum in-plane translation (initial):  {} [A]\n", settings.maximumTranslationInPlane);
  std::print(stream, "In-plane reinsertion probability:        {}\n", settings.reinsertionInPlaneProbability);
  std::print(stream, "Per-pair records:                        {}\n",
             settings.writePairRecords ? "dctst/pairs_{T}_{P}.r{k}.txt (one line per completed pair)" : "off");
  if (settings.snapshotInterval > 0.0)
    std::print(stream, "Accumulator snapshots:                   dctst/accumulator_{{T}}_{{P}}.r{{k}}.json every {} [s]\n",
               settings.snapshotInterval);
  else
    std::print(stream, "Accumulator snapshots:                   off\n");
  std::print(stream, "Termination signal (TERM/INT/USR1):      stops at the next cycle boundary (an unfinished pair is "
                     "discarded) and writes the results of the completed pairs\n");
  std::print(
      stream,
      "Cage lattice (occupancies):              {} x {} x {} cages, center of cage B at {:.4f} {:.4f} {:.4f} [A]\n",
      cageLattice.numberOfUnitCells.x, cageLattice.numberOfUnitCells.y, cageLattice.numberOfUnitCells.z,
      cageLattice.origin.x, cageLattice.origin.y, cageLattice.origin.z);
  {
    // the same in unit-cell fractional coordinates, reduced to [0, 1): the user can check that this is a cage
    // center of the framework (for LTA with the sodalite cage at the origin: 0.5 0.5 0.5)
    const double3 f = cageLattice.inverseCell * cageLattice.origin;
    const double3 reduced(f.x - std::floor(f.x), f.y - std::floor(f.y), f.z - std::floor(f.z));
    std::print(stream, "Cage center (unit-cell fractional):      {:.4f} {:.4f} {:.4f}{}\n", reduced.x, reduced.y,
               reduced.z,
               settings.cageCenter.has_value() ? " ('CageCenter')"
                                               : " (derived from 'BarrierPosition'; set 'CageCenter' to be explicit)");
  }
  if (settings.cageOccupancies.has_value())
  {
    std::print(stream,
               "Occupancy-resolved mode:                 n_A = {} (tagged molecule included), n_B = {}; moves of other "
               "molecules that change either occupancy are rejected\n",
               settings.cageOccupancies->first, settings.cageOccupancies->second);
  }
  std::print(stream, "\n");
  std::flush(stream);

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["initialConditions"] = front.jsonSystemStatus();
  outputJson["initialization"]["components"] = front.jsonComponentStatus();
  outputJson["dcTST"]["taggedComponent"] = tagged.name;
  outputJson["dcTST"]["reactionBead"] = tagged.startingBead;
  outputJson["dcTST"]["barrierPosition"] = {settings.barrierPosition.x, settings.barrierPosition.y,
                                            settings.barrierPosition.z};
  outputJson["dcTST"]["barrierNormal"] = {barrier.normal.x, barrier.normal.y, barrier.normal.z};
  outputJson["dcTST"]["latticeDistance"] = barrier.lambda;
  outputJson["dcTST"]["numberOfReplicas"] = numberOfReplicas;
  outputJson["dcTST"]["trajectoryTime"] = settings.trajectoryTime;
  outputJson["dcTST"]["timeStep"] = front.timeStep;
  if (settings.cageOccupancies.has_value())
  {
    outputJson["dcTST"]["cageOccupancies"] = {settings.cageOccupancies->first, settings.cageOccupancies->second};
  }
  outputJson["dcTST"]["trajectoryTermination"] =
      (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event) ? "Event" : "FixedTime";
  outputJson["dcTST"]["terminationEvents"] = settings.terminationEvents;

  // per-replica output files
  replicaStreams.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(
        std::format("output/output_{}_{}.dctst.r{}.txt", system.temperature, system.input_pressure, replicaId),
        std::ios::out);
    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "dcTST: replica {} of {}\n", replicaId, numberOfReplicas);
    std::print(replicaStream, "Random seed of this replica: {}\n\n", randoms[replicaId].seed);
    std::print(replicaStream, "{}", system.writeSystemStatus());
    std::print(replicaStream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(replicaStream, "{}", system.forceField.printForceFieldStatus());
    std::print(replicaStream, "{}", system.writeComponentStatus());
    std::print(replicaStream, "{}", system.writeNumberOfPseudoAtoms());
    std::flush(replicaStream);
  }

  plateauIndices = computePlateauIndices();

  // per-replica per-pair records
  if (settings.writePairRecords)
  {
    pairStreams.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      const System& system = systems[replicaId];
      pairStreams.emplace_back(
          std::format("dctst/pairs_{}_{}.r{}.txt", system.temperature, system.input_pressure, replicaId), std::ios::out);
      writePairRecordHeader(replicaId);
    }
  }

  // interpolation grids are computed once and shared (copied) between the replicas
  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t replicaId = 1; replicaId < systems.size(); ++replicaId)
  {
    systems[replicaId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[replicaId].interpolationGrids = systems.front().interpolationGrids;
  }
}

void DynamicallyCorrectedTST::loadFreeEnergyInput()
{
  const System& front = systems.front();
  const Component& tagged = front.components[settings.taggedComponent];

  if (settings.barrierProbability.has_value())
  {
    barrierProbability = settings.barrierProbability;
    barrierProbabilityError = settings.barrierProbabilityError;
    barrierProbabilitySource = "'BarrierProbability' from the input";
  }
  else if (settings.freeEnergyInput.has_value())
  {
    const std::string& path = settings.freeEnergyInput.value();
    std::ifstream file(path);
    if (!file.is_open())
    {
      throw std::runtime_error(std::format("dcTST: cannot open 'FreeEnergyInput' '{}'\n", path));
    }
    nlohmann::json json;
    try
    {
      json = nlohmann::json::parse(file);
    }
    catch (const std::exception& e)
    {
      throw std::runtime_error(std::format("dcTST: 'FreeEnergyInput' '{}' is not valid JSON: {}\n", path, e.what()));
    }
    if (!json.contains("dcTST") || !json.contains("results"))
    {
      throw std::runtime_error(std::format(
          "dcTST: 'FreeEnergyInput' '{}' is not the JSON output of the free-energy driver (no 'dcTST'/'results')\n",
          path));
    }
    const nlohmann::json& header = json["dcTST"];
    const nlohmann::json& results = json["results"];

    // the free-energy run must describe the same barrier
    if (header.contains("taggedComponent") && header["taggedComponent"].get<std::string>() != tagged.name)
    {
      throw std::runtime_error(std::format("dcTST: 'FreeEnergyInput' was computed for component '{}', this run tags '{}'\n",
                                           header["taggedComponent"].get<std::string>(), tagged.name));
    }
    if (header.contains("latticeDistance") &&
        std::abs(header["latticeDistance"].get<double>() - barrier.lambda) > 1e-6 * barrier.lambda)
    {
      throw std::runtime_error(std::format("dcTST: 'FreeEnergyInput' has 'LatticeDistance' {}, this run {}\n",
                                           header["latticeDistance"].get<double>(), barrier.lambda));
    }
    if (header.contains("barrierNormal") && header["barrierNormal"].is_array() && header["barrierNormal"].size() == 3)
    {
      const double3 n(header["barrierNormal"][0].get<double>(), header["barrierNormal"][1].get<double>(),
                      header["barrierNormal"][2].get<double>());
      if (std::abs(std::abs(double3::dot(n, barrier.normal)) - 1.0) > 1e-6)
      {
        throw std::runtime_error("dcTST: 'FreeEnergyInput' was computed for a different barrier normal\n");
      }
    }
    if (header.contains("temperature") && std::abs(header["temperature"].get<double>() - front.temperature) > 1e-6)
    {
      throw std::runtime_error(std::format("dcTST: 'FreeEnergyInput' was computed at {} K, this run is at {} K\n",
                                           header["temperature"].get<double>(), front.temperature));
    }
    if (header.contains("numberOfMolecules") && header["numberOfMolecules"].get<std::size_t>() != front.numberOfMolecules())
    {
      std::print(stream,
                 "WARNING: 'FreeEnergyInput' was computed with {} molecules, this run has {} (different box or "
                 "loading)\n\n",
                 header["numberOfMolecules"].get<std::size_t>(), front.numberOfMolecules());
    }

    if (settings.cageOccupancies.has_value())
    {
      if (!results.contains("conditioned") || !results["conditioned"].contains("barrierProbability"))
      {
        throw std::runtime_error(
            "dcTST: 'CageOccupancies' needs a free-energy run with the same 'CageOccupancies' ('results.conditioned' "
            "missing in 'FreeEnergyInput')\n");
      }
      const nlohmann::json& conditioned = results["conditioned"];
      if (conditioned.contains("cageOccupancies") &&
          (conditioned["cageOccupancies"][0].get<std::size_t>() != settings.cageOccupancies->first ||
           conditioned["cageOccupancies"][1].get<std::size_t>() != settings.cageOccupancies->second))
      {
        throw std::runtime_error(std::format(
            "dcTST: 'FreeEnergyInput' is conditioned on occupancies [{}, {}], this run on [{}, {}]\n",
            conditioned["cageOccupancies"][0].get<std::size_t>(), conditioned["cageOccupancies"][1].get<std::size_t>(),
            settings.cageOccupancies->first, settings.cageOccupancies->second));
      }
      barrierProbability = conditioned["barrierProbability"].get<double>();
      barrierProbabilityError = conditioned.value("barrierProbabilityError", 0.0);
      barrierProbabilitySource = std::format("free-energy input '{}', conditioned", path);
    }
    else
    {
      if (results.contains("conditioned"))
      {
        throw std::runtime_error(
            "dcTST: 'FreeEnergyInput' is an occupancy-conditioned run; set the same 'CageOccupancies' here\n");
      }
      if (!results.contains("barrierProbability"))
      {
        throw std::runtime_error("dcTST: 'FreeEnergyInput' has no 'results.barrierProbability'\n");
      }
      barrierProbability = results["barrierProbability"].get<double>();
      barrierProbabilityError = results.value("barrierProbabilityError", 0.0);
      barrierProbabilitySource = std::format("free-energy input '{}'", path);
    }
    if (barrierProbability.value() <= 0.0)
    {
      throw std::runtime_error("dcTST: 'FreeEnergyInput' has a non-positive P(q*)\n");
    }
  }

  if (barrierProbability.has_value())
  {
    std::print(stream, "Barrier probability density P(q*):       {:.6e} +/- {:.3e} [1/A] ({})\n\n",
               barrierProbability.value(), barrierProbabilityError, barrierProbabilitySource);
    outputJson["dcTST"]["barrierProbability"] = barrierProbability.value();
    outputJson["dcTST"]["barrierProbabilitySource"] = barrierProbabilitySource;
  }
  else
  {
    std::print(stream, "Barrier probability density P(q*):       not given; rates and D are not reported\n\n");
  }
  std::flush(stream);
}

std::size_t DynamicallyCorrectedTST::taggedMoleculeIndex(const System& system) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < settings.taggedComponent; ++i)
  {
    index += system.numberOfMoleculesPerComponent[i];
  }
  return index;
}

// =================================================================================================
// Stages
// =================================================================================================

void DynamicallyCorrectedTST::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  const std::string_view stageName = (stage == SimulationStage::Initialization) ? "Initialization" : "Transmission";
  std::print(stream, "\n{} stage: {} MC cycles on {} replicas/threads\n", stageName, numberOfCycles, numberOfReplicas);
  std::flush(stream);

  const std::size_t equilibrationCycles =
      (stage == SimulationStage::Transmission) ? settings.numberOfConstrainedEquilibrationCycles : 0uz;
  const std::size_t productionCycles = numberOfCycles - std::min(numberOfCycles, equilibrationCycles);

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, equilibrationCycles, productionCycles]()
          {
            System& system = systems[replicaId];
            std::ostream replicaStream(replicaStreams[replicaId].rdbuf());

            if (stage == SimulationStage::Initialization)
            {
              createConstrainedConfiguration(replicaId);
            }
            if (stage == SimulationStage::Transmission)
            {
              if (settings.cageOccupancies.has_value())
              {
                prepareOccupancies(replicaId);
              }
            }
            if (stage != SimulationStage::Initialization)
            {
              system.mc_moves_statistics.clearMoveStatistics();
              system.mc_moves_cputime.clearTimingStatistics();
              for (Component& component : system.components)
              {
                component.mc_moves_statistics.clearMoveStatistics();
                component.mc_moves_cputime.clearTimingStatistics();
              }
              // the in-plane move counters restart with the transmission stage; the adapted maximum
              // displacement is kept
              DcTSTConstrainedMoveStatistics& moves = constrainedMoves[replicaId];
              const double maximumTranslation = moves.maximumTranslation;
              moves = DcTSTConstrainedMoveStatistics{};
              moves.maximumTranslation = maximumTranslation;
            }

            system.precomputeTotalRigidEnergy();
            system.runningEnergies = system.computeTotalEnergies();
            std::print(replicaStream, "\n{} stage\n",
                       (stage == SimulationStage::Initialization) ? "Initialization" : "Transmission");
            std::print(replicaStream, "{}", system.runningEnergies.printMC("Recomputed from scratch"));
            std::print(replicaStream, "\n\n");
            std::flush(replicaStream);

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, productionCycles));
            const bool snapshots = (stage == SimulationStage::Transmission) && settings.snapshotInterval > 0.0;
            lastSnapshotTime[replicaId] = std::chrono::steady_clock::now();

            for (std::size_t cycle = 0; cycle != numberOfCycles; ++cycle)
            {
              // a termination signal (e.g. the scheduler's warning before the wall-time kill): stop at this
              // cycle boundary; the analysis then runs over the completed trajectory pairs
              if (GracefulShutdown::requested())
              {
                stoppedEarly.store(true, std::memory_order_relaxed);
                std::print(replicaStream, "Termination signal received: stopping at cycle {} of {} ({})\n\n", cycle,
                           numberOfCycles, (cycle >= equilibrationCycles) ? "production" : "equilibration");
                std::flush(replicaStream);
                break;
              }

              const bool production = cycle >= equilibrationCycles;
              if (production)
              {
                estimation.setCurrentSample(cycle - equilibrationCycles);
              }
              const std::size_t block = estimation.currentBin;

              performMonteCarloCycle(replicaId);

              if (stage == SimulationStage::Transmission && production && settings.launchTrajectoryEvery > 0uz &&
                  (cycle - equilibrationCycles) % settings.launchTrajectoryEvery == 0uz)
              {
                launchTrajectories(replicaId, block, cycle - equilibrationCycles);
                if (GracefulShutdown::requested()) continue;  // the pair was abandoned; the loop head stops the run
              }
              if (stage == SimulationStage::Transmission && production)
              {
                completedTransmissionCycles[replicaId] = cycle - equilibrationCycles + 1uz;
              }
              if (snapshots && production)
              {
                const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
                if (std::chrono::duration<double>(now - lastSnapshotTime[replicaId]).count() >= settings.snapshotInterval)
                {
                  writeAccumulatorSnapshot(replicaId, false);
                  lastSnapshotTime[replicaId] = now;
                }
              }

              if (cycle % optimizeMCMovesEvery == 0uz)
              {
                system.optimizeMCMoves();
                constrainedMoves[replicaId].optimize();
              }

              if (cycle % printEvery == 0uz)
              {
                std::print(replicaStream, "Cycle {} of {} ({})\n", cycle, numberOfCycles,
                           production ? "production" : "equilibration");
                std::print(replicaStream, "{}", system.runningEnergies.printMC(""));
                {
                  const DcTSTConstrainedMoveStatistics& moves = constrainedMoves[replicaId];
                  const std::size_t translationTrials = moves.translationTotalTrials + moves.translationTrials;
                  const std::size_t translationAccepted = moves.translationTotalAccepted + moves.translationAccepted;
                  const std::size_t reinsertionTrials = moves.reinsertionTotalTrials + moves.reinsertionTrials;
                  const std::size_t reinsertionAccepted = moves.reinsertionTotalAccepted + moves.reinsertionAccepted;
                  std::print(replicaStream,
                             "    in-plane translation: {}/{} accepted ({:.2f}%), maximum displacement {:.4f} [A]\n",
                             translationAccepted, translationTrials,
                             100.0 * static_cast<double>(translationAccepted) /
                                 static_cast<double>(std::max(1uz, translationTrials)),
                             moves.maximumTranslation);
                  std::print(replicaStream, "    in-plane reinsertion: {}/{} accepted ({:.2f}%)\n", reinsertionAccepted,
                             reinsertionTrials,
                             100.0 * static_cast<double>(reinsertionAccepted) /
                                 static_cast<double>(std::max(1uz, reinsertionTrials)));
                }
                if (stage == SimulationStage::Transmission)
                {
                  const DcTSTTransmissionAccumulator& accumulator = transmission[replicaId];
                  if (settings.cageOccupancies.has_value())
                  {
                    const DcTSTOccupancyState& state = occupancyStates[replicaId];
                    std::print(
                        replicaStream,
                        "    occupancy constraint: cage A {} (n = {}), cage B {} (n = {}); {} cage changes of other "
                        "molecules of which {} rejected, {} in-plane moves rejected\n",
                        state.cageA, state.occupancy[state.cageA] + 1uz, state.cageB, state.occupancy[state.cageB],
                        state.cageChanges, state.rejectedMoves, state.rejectedTaggedMoves);
                  }
                  const std::size_t pairs = std::accumulate(accumulator.numberOfTrajectoryPairs.begin(),
                                                            accumulator.numberOfTrajectoryPairs.end(), 0uz);
                  double speed{0.0};
                  double selfEnd{0.0};
                  double collectiveEnd{0.0};
                  for (std::size_t b = 0; b < numberOfBlocks; ++b)
                  {
                    speed += accumulator.speed[b];
                    selfEnd += accumulator.self[b][numberOfTimes - 1];
                    // all shells and components (the analysis sums the shells cumulatively)
                    for (std::size_t s = 0; s < accumulator.numberOfShells; ++s)
                    {
                      for (std::size_t c = 0; c < accumulator.numberOfComponents; ++c)
                      {
                        collectiveEnd +=
                            accumulator.collective[b][accumulator.collectiveIndex(s, c, numberOfTimes - 1)];
                      }
                    }
                  }
                  if (pairs > 0uz && speed > 0.0)
                  {
                    std::print(replicaStream,
                               "    trajectory pairs: {}, <|qdot|> = {:.5f} [A/ps], kappa_S(t_max) = {:.5f}, "
                               "kappa_C(t_max) = {:.5f}, mean |energy drift| = {:.3e}\n",
                               pairs, speed / static_cast<double>(pairs), selfEnd / speed,
                               (selfEnd + collectiveEnd) / speed,
                               accumulator.sumRelativeEnergyDrift /
                                   static_cast<double>(std::max(1uz, accumulator.numberOfTrajectories)));
                  }
                }
                std::print(replicaStream, "\n");
                std::flush(replicaStream);
              }
            }

            if (snapshots)
            {
              writeAccumulatorSnapshot(replicaId, true);
            }

            // consistency check of the running energies at the end of the stage
            RunningEnergy recomputed = system.computeTotalEnergies();
            std::print(replicaStream, "\n{}", recomputed.printMCDiff(system.runningEnergies));
            if (stage == SimulationStage::Transmission)
            {
              std::print(replicaStream, "\n{}\n", system.writeMCMoveStatistics());
            }
            std::flush(replicaStream);
          });
    }
  }

  if (stoppedEarly.load(std::memory_order_relaxed))
  {
    if (stage == SimulationStage::Transmission)
    {
      const auto [minimum, maximum] =
          std::minmax_element(completedTransmissionCycles.begin(), completedTransmissionCycles.end());
      std::print(stream,
                 "Termination signal received: the transmission stage was stopped after {}..{} of {} production cycles "
                 "per replica; the results are over the completed trajectory pairs\n",
                 *minimum, *maximum, productionCycles);
    }
    else
    {
      std::print(stream, "Termination signal received during the {} stage: stopping\n", stageName);
    }
    std::flush(stream);
  }

  const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  switch (stage)
  {
    case SimulationStage::Initialization:
      totalInitializationTime += (t2 - t1);
      break;
    case SimulationStage::Transmission:
      totalTransmissionTime += (t2 - t1);
      break;
    default:
      break;
  }
  totalSimulationTime += (t2 - t1);
}

// =================================================================================================
// Monte Carlo
// =================================================================================================

void DynamicallyCorrectedTST::performMonteCarloCycle(std::size_t replicaId)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];

  const std::size_t numberOfSteps = std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();
  // the occupancy constraint exists from the start of the transmission stage ('prepareOccupancies')
  const bool occupancyMode = occupancyStates[replicaId].active;

  for (std::size_t step = 0; step != numberOfSteps; ++step)
  {
    const std::size_t selectedComponent = system.randomComponent(rng);

    // both stages sample the constrained ensemble: the tagged molecule is restricted to its in-plane
    // moves, all other molecules move with the single-molecule moves of their component
    const std::size_t numberOfMolecules = system.numberOfMoleculesPerComponent[selectedComponent];
    if (numberOfMolecules == 0uz) continue;
    const std::size_t selectedMolecule = std::min(
        numberOfMolecules - 1uz, static_cast<std::size_t>(rng.uniform() * static_cast<double>(numberOfMolecules)));

    std::optional<RunningEnergy> energyDifference;
    if (selectedComponent == settings.taggedComponent && selectedMolecule == 0uz)
    {
      energyDifference = (rng.uniform() < settings.reinsertionInPlaneProbability) ? reinsertionInPlaneMove(replicaId)
                                                                                  : translationInPlaneMove(replicaId);
    }
    else
    {
      Component& component = system.components[selectedComponent];
      const Move::Types moveType = component.mc_moves_probabilities.sample(rng);
      component.mc_moves_statistics.addAllCounts(moveType);

      // occupancy-resolved mode: the old configuration of the molecule, to undo an accepted move that changes
      // the occupancy of cage A or B (the constraint is a hard wall in configuration space: rejecting such
      // moves keeps detailed balance since the allowed set is the same seen from either side)
      const std::size_t moleculeIndex =
          occupancyMode ? system.moleculeIndexOfComponent(selectedComponent, selectedMolecule) : 0uz;
      std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
      std::vector<Atom> oldAtoms;
      Molecule oldMolecule;
      if (occupancyMode)
      {
        oldAtoms.assign(moleculeAtoms.begin(), moleculeAtoms.end());
        oldMolecule = system.moleculeData[moleculeIndex];
      }

      switch (moveType)
      {
        case Move::Types::RandomTranslation:
          energyDifference = MC_Moves::randomTranslationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::Rotation:
          energyDifference = MC_Moves::rotationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::RandomRotation:
          energyDifference = MC_Moves::randomRotationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::ReinsertionCBMC:
          energyDifference = MC_Moves::reinsertionMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::PartialReinsertionCBMC:
          energyDifference = MC_Moves::partialReinsertionMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::Translation:
        default:
          energyDifference = MC_Moves::translationMove(rng, system, selectedComponent, selectedMolecule);
          break;
      }

      if (occupancyMode && energyDifference.has_value())
      {
        DcTSTOccupancyState& state = occupancyStates[replicaId];
        const std::size_t oldCage = state.cageOfMolecule[moleculeIndex];
        const std::size_t newCage = cageLattice.cageOf(moleculeAtoms[component.startingBead].position);
        if (newCage != oldCage)
        {
          ++state.cageChanges;
          if (oldCage == state.cageA || oldCage == state.cageB || newCage == state.cageA || newCage == state.cageB)
          {
            // undo: positions, molecule record, cell list and the Ewald structure factor
            if (system.forceField.useCharge && system.forceField.usesEwaldFourier())
            {
              Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, system.components, oldAtoms, moleculeAtoms);
              Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
            }
            std::copy(oldAtoms.cbegin(), oldAtoms.cend(), moleculeAtoms.begin());
            system.cellListAtomsMoved(moleculeAtoms);
            system.moleculeData[moleculeIndex] = oldMolecule;
            energyDifference = std::nullopt;
            ++state.rejectedMoves;
          }
          else
          {
            --state.occupancy[oldCage];
            ++state.occupancy[newCage];
            state.cageOfMolecule[moleculeIndex] = newCage;
          }
        }
      }
    }
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
    }
  }
}

std::optional<RunningEnergy> DynamicallyCorrectedTST::translationInPlaneMove(std::size_t replicaId)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTConstrainedMoveStatistics& statistics = constrainedMoves[replicaId];

  const std::size_t selectedComponent = settings.taggedComponent;
  Component& component = system.components[selectedComponent];
  std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, 0uz);
  Molecule& molecule = system.moleculeData[taggedMoleculeIndex(system)];

  ++statistics.translationTrials;

  const double du = statistics.maximumTranslation * 2.0 * (rng.uniform() - 0.5);
  const double dv = statistics.maximumTranslation * 2.0 * (rng.uniform() - 0.5);
  const double3 displacement = du * barrier.inPlaneU + dv * barrier.inPlaneV;

  std::pair<Molecule, std::vector<Atom>> trialMolecule = component.translate(molecule, moleculeAtoms, displacement);

  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  // occupancy-resolved mode: the window the bead moves to must have the requested occupancies
  std::pair<std::size_t, std::size_t> trialCages{};
  DcTSTOccupancyState& state = occupancyStates[replicaId];
  if (state.active)
  {
    trialCages = cagesOfWindow(trialMolecule.second[component.startingBead].position);
    if (!occupanciesSatisfied(state, trialCages.first, trialCages.second))
    {
      ++state.rejectedTaggedMoves;
      return std::nullopt;
    }
  }

  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialMolecule.second, moleculeAtoms);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, system.frameworkCellList(),
      system.spanOfFrameworkAtoms(), trialMolecule.second, moleculeAtoms);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(), trialMolecule.second,
      moleculeAtoms);
  if (!interMolecule.has_value()) return std::nullopt;

  const RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, trialMolecule.second, moleculeAtoms);

  const RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  if (rng.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    ++statistics.translationAccepted;
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), moleculeAtoms.begin());
    system.cellListAtomsMoved(moleculeAtoms);
    molecule = trialMolecule.first;
    if (state.active)
    {
      state.cageA = trialCages.first;
      state.cageB = trialCages.second;
    }
    return energyDifference;
  }
  return std::nullopt;
}

std::optional<RunningEnergy> DynamicallyCorrectedTST::reinsertionInPlaneMove(std::size_t replicaId)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTConstrainedMoveStatistics& statistics = constrainedMoves[replicaId];

  const std::size_t selectedComponent = settings.taggedComponent;
  const Component& component = system.components[selectedComponent];
  std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, 0uz);
  Molecule& molecule = system.moleculeData[taggedMoleculeIndex(system)];

  ++statistics.reinsertionTrials;

  // the reaction bead is pinned at a uniformly random point of the barrier plane; the remainder of the
  // molecule (orientation of a rigid molecule, the chain of a flexible one) is grown with CBMC. The trial
  // molecule carries the id of the old molecule, which excludes the old copy from the background.
  const double3 point = barrier.randomPointOnPlane(rng, system.simulationBox);

  // occupancy-resolved mode: the window of the new point must have the requested occupancies
  std::pair<std::size_t, std::size_t> trialCages{};
  DcTSTOccupancyState& state = occupancyStates[replicaId];
  if (state.active)
  {
    trialCages = cagesOfWindow(point);
    if (!occupanciesSatisfied(state, trialCages.first, trialCages.second))
    {
      ++state.rejectedTaggedMoves;
      return std::nullopt;
    }
  }

  const CBMC::GrowContext context = system.makeGrowContext();
  std::optional<CBMC::GrowResult> growData = CBMC::growNewMolecule(
      rng, context, component,
      CBMC::NewMoleculeIdentity{.componentId = selectedComponent,
                                .moleculeId = static_cast<std::size_t>(moleculeAtoms.front().moleculeId)},
      CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = point});
  if (!growData) return std::nullopt;

  if (system.insideBlockedPockets(component, growData->atoms))
  {
    return std::nullopt;
  }

  const CBMC::RetraceResult retraceData =
      CBMC::retraceMolecule(rng, context, component, moleculeAtoms, {.firstBead = CBMC::FirstBeadScheme::Pinned});

  const RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, growData->atoms, moleculeAtoms);

  const double logAcceptance = -system.beta * energyFourierDifference.potentialEnergy() +
                               growData->logRosenbluthWeight - retraceData.logRosenbluthWeight;

  if (rng.uniform() < std::exp(logAcceptance))
  {
    ++statistics.reinsertionAccepted;
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
    {
      moleculeAtoms[i].position = growData->atoms[i].position;
    }
    system.cellListAtomsMoved(moleculeAtoms);

    Molecule record = growData->molecule;
    record.atomIndex = molecule.atomIndex;
    record.numberOfAtoms = molecule.numberOfAtoms;
    molecule = record;
    if (state.active)
    {
      state.cageA = trialCages.first;
      state.cageB = trialCages.second;
    }

    return (growData->energies - retraceData.energies) + energyFourierDifference;
  }
  return std::nullopt;
}

void DynamicallyCorrectedTST::createConstrainedConfiguration(std::size_t replicaId)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  std::ostream replicaStream(replicaStreams[replicaId].rdbuf());

  const std::size_t taggedComponent = settings.taggedComponent;
  const Component& tagged = system.components[taggedComponent];

  // the number of molecules of every component the system was created with ('CreateNumberOfMolecules')
  const std::vector<std::size_t> targetNumberOfMolecules = system.numberOfMoleculesPerComponent;
  if (targetNumberOfMolecules[taggedComponent] == 0uz)
  {
    throw std::runtime_error(
        std::format("dcTST: the tagged component '{}' has no molecules (replica {})\n", tagged.name, replicaId));
  }

  // Order of creation (as in RASPA2): the system is emptied, the tagged molecule is grown first with its
  // reaction bead on the barrier plane, and only then are the other molecules grown around it. Placing the
  // tagged molecule into an already filled system fails at high loading -- every window may hold a
  // molecule -- whereas in the empty framework the growth on a window centre always succeeds. The
  // molecules created with the system are discarded (they are a random configuration anyway).
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    while (system.numberOfMoleculesPerComponent[componentId] > 0uz)
    {
      const std::size_t last = system.numberOfMoleculesPerComponent[componentId] - 1uz;
      system.deleteMolecule(componentId, last, system.spanOfMolecule(componentId, last));
    }
  }

  // a grown molecule is valid when it exists and does not overlap with its environment
  const auto isValid = [&](const std::optional<CBMC::GrowResult>& growData)
  { return growData.has_value() && growData->energies.potentialEnergy() <= system.forceField.energyOverlapCriteria; };

  // the context is rebuilt per growth: every inserted molecule changes the background it spans
  const auto creationContext = [&]()
  { return system.makeGrowContext(CBMC::CutOffMode::Full).withChainScheme(CBMC::ChainScheme::ConfigurationalBias); };

  // --- the tagged molecule: reaction bead pinned on the barrier plane ---------------------------------
  const auto growTaggedAt = [&](double3 point) -> std::optional<CBMC::GrowResult>
  {
    std::optional<CBMC::GrowResult> growData = CBMC::growNewMolecule(
        rng, creationContext(), tagged,
        CBMC::NewMoleculeIdentity{.componentId = taggedComponent, .moleculeId = system.numberOfMolecules()},
        CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = point});
    if (!isValid(growData) || system.insideBlockedPockets(tagged, growData->atoms)) return std::nullopt;
    return growData;
  };

  // the window centres of this plane (midpoints of neighbouring cage centres along the normal), in random
  // order so that the replicas start in different windows
  std::vector<double3> windowPoints;
  const DcTSTCageLattice::AxisAlignment alignment = cageLattice.axisAlongDirection(barrier.normal);
  for (std::size_t c = 0; c < cageLattice.numberOfCages(); ++c)
  {
    const double3 mid =
        0.5 * (cageLattice.center(c) + cageLattice.center(cageLattice.neighborAlongAxis(c, alignment.axis, alignment.sign)));
    if (std::abs(barrier.reducedCoordinate(barrier.projection(system.simulationBox, mid))) < 0.1 * barrier.lambda)
    {
      windowPoints.push_back(barrier.projectOntoPlane(system.simulationBox, mid, 0.0));
    }
  }
  for (std::size_t i = windowPoints.size(); i > 1uz; --i)
  {
    const std::size_t j = std::min(i - 1uz, static_cast<std::size_t>(rng.uniform() * static_cast<double>(i)));
    std::swap(windowPoints[i - 1uz], windowPoints[j]);
  }

  std::optional<CBMC::GrowResult> taggedGrow{};
  std::string placement{};
  std::size_t attempts{0};
  constexpr std::size_t jitterPerWindow = 200uz;
  constexpr std::size_t uniformAttempts = 100000uz;
  for (const double3& window : windowPoints)
  {
    ++attempts;
    taggedGrow = growTaggedAt(window);
    if (taggedGrow.has_value())
    {
      placement = "a window centre";
      break;
    }
    for (std::size_t j = 0; j < jitterPerWindow && !taggedGrow.has_value(); ++j)
    {
      ++attempts;
      const double3 point = window + (rng.uniform() - 0.5) * 2.0 * barrier.inPlaneU +
                            (rng.uniform() - 0.5) * 2.0 * barrier.inPlaneV;
      taggedGrow = growTaggedAt(barrier.projectOntoPlane(system.simulationBox, point, 0.0));
      placement = "a point near a window centre";
    }
    if (taggedGrow.has_value()) break;
  }
  // no cage lattice that matches the plane (bulk system, or a barrier through a different cross-section):
  // uniform samples of the plane
  for (std::size_t attempt = 0; attempt < uniformAttempts && !taggedGrow.has_value(); ++attempt)
  {
    ++attempts;
    taggedGrow = growTaggedAt(barrier.randomPointOnPlane(rng, system.simulationBox));
    placement = "a uniformly sampled point of the plane";
  }
  if (!taggedGrow.has_value())
  {
    throw std::runtime_error(
        std::format("dcTST: the tagged molecule could not be grown with its reaction bead on the barrier plane of "
                    "the empty framework after {} attempts (replica {}); check 'BarrierPosition', 'BarrierNormal' "
                    "and the blocking pockets\n",
                    attempts, replicaId));
  }
  system.insertMolecule(taggedComponent, taggedGrow->molecule, taggedGrow->atoms);

  // --- the other molecules: grown around the tagged one -------------------------------------------------
  // (the same creation scheme as 'System::createInitialMolecules': the first valid growth starts a short
  // Markov chain of further growths accepted with the Rosenbluth weights, which corrects the conformational
  // bias of a single growth of a flexible molecule; a rigid molecule is grown once)
  const auto growCreationCandidate = [&](std::size_t componentId) -> CBMC::GrowResult
  {
    const Component& component = system.components[componentId];
    const CBMC::NewMoleculeIdentity identity{.componentId = componentId, .moleculeId = system.numberOfMolecules()};

    std::optional<CBMC::GrowResult> current = std::nullopt;
    do
    {
      current = CBMC::growNewMolecule(rng, creationContext(), component, identity);
    } while (!isValid(current));

    const std::size_t numberOfGrowths = component.rigid ? 1uz : CBMC::Constants::creationCandidateGrowths;
    for (std::size_t attempt = 1; attempt < numberOfGrowths; ++attempt)
    {
      std::optional<CBMC::GrowResult> candidate = CBMC::growNewMolecule(rng, creationContext(), component, identity);
      if (!isValid(candidate)) continue;

      const double logAcceptance = candidate->logRosenbluthWeight - current->logRosenbluthWeight;
      if (logAcceptance >= 0.0 || rng.uniform() < std::exp(logAcceptance))
      {
        current = std::move(candidate);
      }
    }
    return std::move(current.value());
  };

  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    while (system.numberOfMoleculesPerComponent[componentId] < targetNumberOfMolecules[componentId])
    {
      const CBMC::GrowResult growData = growCreationCandidate(componentId);
      system.insertMolecule(componentId, growData.molecule, growData.atoms);
    }
  }
  system.invalidateCellList();

  const std::span<const Atom> taggedAtoms = system.spanOfMolecule(taggedComponent, 0uz);
  std::print(replicaStream, "Constrained starting configuration\n");
  std::print(replicaStream, "    tagged molecule grown first on {} after {} attempt(s) (reduced q = {:.3e} [A])\n",
             placement, attempts,
             barrier.reducedCoordinate(
                 barrier.projection(system.simulationBox, taggedAtoms[tagged.startingBead].position)));
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    std::print(replicaStream, "    component {} ({}): {} molecules\n", componentId, system.components[componentId].name,
               system.numberOfMoleculesPerComponent[componentId]);
  }
  std::print(replicaStream, "\n");
  std::flush(replicaStream);
}

// =================================================================================================
// Occupancy-resolved mode
// =================================================================================================

std::pair<std::size_t, std::size_t> DynamicallyCorrectedTST::cagesOfWindow(double3 beadPosition) const
{
  // the bead is on the plane, equidistant from the centers of A and B: step a quarter cage to either side
  const double3 step = 0.25 * barrier.lambda * barrier.normal;
  return {cageLattice.cageOf(beadPosition - step), cageLattice.cageOf(beadPosition + step)};
}

bool DynamicallyCorrectedTST::occupanciesSatisfied(const DcTSTOccupancyState& state, std::size_t cageA,
                                                   std::size_t cageB) const
{
  if (cageA == cageB) return false;
  const auto& [nA, nB] = settings.cageOccupancies.value();
  return state.occupancy[cageA] + 1uz == nA && state.occupancy[cageB] == nB;
}

std::vector<std::size_t> DynamicallyCorrectedTST::countOccupancies(const System& system) const
{
  std::vector<std::size_t> occupancy(cageLattice.numberOfCages(), 0uz);
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  for (const Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    ++occupancy[cageLattice.cageOf(atoms[molecule.atomIndex + component.startingBead].position)];
  }
  return occupancy;
}

void DynamicallyCorrectedTST::prepareOccupancies(std::size_t replicaId)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTOccupancyState& state = occupancyStates[replicaId];
  std::ostream replicaStream(replicaStreams[replicaId].rdbuf());

  const auto [nA, nB] = settings.cageOccupancies.value();
  const std::size_t taggedIndex = taggedMoleculeIndex(system);
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t numberOfCages = cageLattice.numberOfCages();

  if (nA - 1uz + nB > numberOfMolecules - 1uz)
  {
    throw std::runtime_error(
        std::format("dcTST: 'CageOccupancies' [{}, {}] needs {} molecules besides the tagged one, the system has {}\n",
                    nA, nB, nA - 1uz + nB, numberOfMolecules - 1uz));
  }

  // bookkeeping from the current configuration (the tagged molecule is excluded from the occupancies)
  auto rebuild = [&]()
  {
    const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    state.cageOfMolecule.assign(numberOfMolecules, 0uz);
    state.occupancy.assign(numberOfCages, 0uz);
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      const Molecule& molecule = system.moleculeData[j];
      const Component& component = system.components[molecule.componentId];
      state.cageOfMolecule[j] = cageLattice.cageOf(atoms[molecule.atomIndex + component.startingBead].position);
      if (j != taggedIndex) ++state.occupancy[state.cageOfMolecule[j]];
    }
  };
  rebuild();
  state.cageChanges = 0;
  state.rejectedMoves = 0;
  state.rejectedTaggedMoves = 0;

  {
    const Molecule& tagged = system.moleculeData[taggedIndex];
    const Component& component = system.components[tagged.componentId];
    const double3 bead = system.spanOfMoleculeAtoms()[tagged.atomIndex + component.startingBead].position;
    std::tie(state.cageA, state.cageB) = cagesOfWindow(bead);
  }
  if (state.cageA == state.cageB)
  {
    throw std::runtime_error(
        "dcTST: cages A and B coincide; the occupancy-resolved mode needs at least two unit cells along the barrier "
        "normal\n");
  }
  const bool hasOtherCages = numberOfCages > 2uz;

  // molecules are moved between cages by a pinned CBMC growth at a random point of the target cage; the
  // first non-overlapping growth is taken (the constrained equilibration relaxes the configuration)
  const CBMC::GrowContext context = system.makeGrowContext(CBMC::CutOffMode::Full);
  auto relocate = [&](std::size_t j, std::size_t targetCage)
  {
    Molecule& molecule = system.moleculeData[j];
    const Component& component = system.components[molecule.componentId];
    std::size_t offset{0};
    for (std::size_t c = 0; c < molecule.componentId; ++c) offset += system.numberOfMoleculesPerComponent[c];
    std::span<Atom> moleculeAtoms = system.spanOfMolecule(molecule.componentId, j - offset);

    for (std::size_t attempt = 0; attempt < 100000uz; ++attempt)
    {
      const double3 point = cageLattice.randomPointInCage(rng, targetCage);
      std::optional<CBMC::GrowResult> growData = CBMC::growNewMolecule(
          rng, context, component,
          CBMC::NewMoleculeIdentity{.componentId = molecule.componentId,
                                    .moleculeId = static_cast<std::size_t>(moleculeAtoms.front().moleculeId)},
          CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = point});
      if (!growData) continue;
      if (system.insideBlockedPockets(component, growData->atoms)) continue;
      if (cageLattice.cageOf(growData->atoms[component.startingBead].position) != targetCage) continue;

      for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
      {
        moleculeAtoms[i].position = growData->atoms[i].position;
      }
      Molecule record = growData->molecule;
      record.atomIndex = molecule.atomIndex;
      record.numberOfAtoms = molecule.numberOfAtoms;
      molecule = record;
      system.invalidateCellList();

      --state.occupancy[state.cageOfMolecule[j]];
      ++state.occupancy[targetCage];
      state.cageOfMolecule[j] = targetCage;
      return;
    }
    throw std::runtime_error(
        std::format("dcTST: could not relocate a molecule into cage {} (replica {})\n", targetCage, replicaId));
  };

  auto randomMoleculeWhere = [&](auto&& predicate) -> std::optional<std::size_t>
  {
    std::vector<std::size_t> candidates;
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      if (j != taggedIndex && predicate(state.cageOfMolecule[j])) candidates.push_back(j);
    }
    if (candidates.empty()) return std::nullopt;
    return candidates[std::min(candidates.size() - 1uz,
                               static_cast<std::size_t>(rng.uniform() * static_cast<double>(candidates.size())))];
  };
  auto randomOtherCage = [&]()
  {
    for (;;)
    {
      const std::size_t cage =
          std::min(numberOfCages - 1uz, static_cast<std::size_t>(rng.uniform() * static_cast<double>(numberOfCages)));
      if (cage != state.cageA && cage != state.cageB) return cage;
    }
  };

  const std::size_t initialA = state.occupancy[state.cageA];
  const std::size_t initialB = state.occupancy[state.cageB];
  std::size_t relocations{0};
  for (const auto [cage, target] : {std::pair{state.cageA, nA - 1uz}, std::pair{state.cageB, nB}})
  {
    while (state.occupancy[cage] > target)
    {
      if (!hasOtherCages)
      {
        throw std::runtime_error("dcTST: no cage besides A and B to move the excess molecules to\n");
      }
      relocate(randomMoleculeWhere([&](std::size_t c) { return c == cage; }).value(), randomOtherCage());
      ++relocations;
    }
  }
  for (const auto [cage, target] : {std::pair{state.cageA, nA - 1uz}, std::pair{state.cageB, nB}})
  {
    while (state.occupancy[cage] < target)
    {
      const std::optional<std::size_t> j =
          randomMoleculeWhere([&](std::size_t c) { return c != state.cageA && c != state.cageB; });
      if (!j.has_value())
      {
        throw std::runtime_error("dcTST: not enough molecules outside cages A and B to fill them\n");
      }
      relocate(j.value(), cage);
      ++relocations;
    }
  }

  std::print(replicaStream,
             "Occupancy constraint: cage A = {} (had {} other molecules, now {}), cage B = {} (had {}, now {}); {} "
             "molecules relocated\n",
             state.cageA, initialA, state.occupancy[state.cageA], state.cageB, initialB, state.occupancy[state.cageB],
             relocations);
  std::flush(replicaStream);

  // from here on the moves enforce the constraint
  state.active = true;
}

// =================================================================================================
// Trajectories
// =================================================================================================

void DynamicallyCorrectedTST::launchTrajectories(std::size_t replicaId, std::size_t block, std::size_t productionCycle)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTTransmissionAccumulator& accumulator = transmission[replicaId];

  const std::size_t taggedIndex = taggedMoleculeIndex(system);
  const Component& taggedComponent = system.components[settings.taggedComponent];
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t numberOfShells = accumulator.numberOfShells;
  const std::size_t numberOfComponents = accumulator.numberOfComponents;
  const std::size_t sampleEvery = settings.sampleTrajectoryEvery;

  // the trajectories run on copies of the dynamical state; the Monte Carlo state of the system is untouched
  const std::span<const Atom> systemAtoms = system.spanOfMoleculeAtoms();
  const std::span<const AtomDynamics> systemDynamics = system.spanOfMoleculeDynamics();
  const std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::optional<Thermostat> noThermostat = std::nullopt;

  // rigid-group states of semi-flexible molecules, derived from the current atom positions
  std::vector<GroupState> initialGroupData;
  {
    std::size_t atomIndex{};
    for (const Molecule& molecule : system.moleculeData)
    {
      const Component& component = system.components[molecule.componentId];
      if (component.isSemiFlexible())
      {
        std::span<const Atom> span = std::span(&systemAtoms[atomIndex], molecule.numberOfAtoms);
        for (std::size_t g = 0; g != component.fragmentGraph.fragments.size(); ++g)
        {
          if (component.fragmentGraph.fragments[g].isRigidBody())
          {
            initialGroupData.push_back(component.deriveFragmentState(g, span));
          }
        }
      }
      atomIndex += molecule.numberOfAtoms;
    }
  }

  // the initial reaction coordinates and cage indices of all molecules (q = 0 for the tagged bead)
  std::vector<std::int64_t> initialCage(numberOfMolecules);
  std::vector<std::size_t> bucketOfMolecule(numberOfMolecules);
  for (std::size_t j = 0; j < numberOfMolecules; ++j)
  {
    const Molecule& molecule = system.moleculeData[j];
    const Component& component = system.components[molecule.componentId];
    const double q =
        barrier.projection(system.simulationBox, systemAtoms[molecule.atomIndex + component.startingBead].position);
    initialCage[j] = barrier.cageIndex(q);
    bucketOfMolecule[j] = std::min(DcTSTBarrier::shellOfCage(initialCage[j]), settings.maximumCollectiveShell + 1uz);
  }

  // working copies and the per-trajectory observables
  std::vector<Molecule> moleculeData;
  std::vector<Atom> atoms;
  std::vector<AtomDynamics> dynamics;
  std::vector<GroupState> groupData;
  std::vector<double> unwrappedQ(numberOfMolecules);
  std::vector<double3> previousPosition(numberOfMolecules);

  std::vector<double> deltaSelf(numberOfTimes);
  std::vector<double> deltaSelfTwoState(numberOfTimes);
  std::vector<double> deltaSelfCommitted(numberOfTimes);
  std::vector<double> deltaCollective(numberOfShells * numberOfComponents * numberOfTimes);

  const DcTSTSettings::TrajectoryThermostat bathType = settings.trajectoryThermostat;
  const double commitDistance = (settings.commitDistance > 0.0) ? settings.commitDistance : 0.5 * barrier.lambda;

  // --- event-terminated estimator: the column of cages A', A, B, B' of the tagged molecule's window -------------
  using EventChannel = DcTSTTransmissionAccumulator::EventChannel;
  const bool eventMode = settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event;
  constexpr std::size_t noCage = std::numeric_limits<std::size_t>::max();
  DcTSTCageLattice::AxisAlignment axis{};
  std::size_t cageA{noCage}, cageB{noCage}, cageAPrime{noCage}, cageC{noCage};
  double3 hopThreshold{};  // fractional distance from the committed cage center at which a hop is registered
  std::vector<std::size_t> initialCommittedCage(numberOfMolecules, noCage);
  std::vector<std::size_t> committedCage(numberOfMolecules, noCage);
  if (eventMode)
  {
    axis = cageLattice.axisAlongDirection(barrier.normal);
    const double3 taggedBead =
        systemAtoms[system.moleculeData[taggedIndex].atomIndex + taggedComponent.startingBead].position;
    std::tie(cageA, cageB) = cagesOfWindow(taggedBead);
    cageAPrime = cageLattice.neighborAlongAxis(cageA, axis.axis, -axis.sign);
    cageC = cageLattice.neighborAlongAxis(cageB, axis.axis, axis.sign);
    // a molecule has hopped to the neighboring cage along a lattice axis when its reaction bead is more than
    // half a cell plus the commit distance away from the center of the cage it was committed to (the same
    // convention as the commitment of the tagged molecule: measured from the window it crossed). This
    // hysteresis keeps a bead oscillating about a window from generating spurious events.
    for (std::size_t i = 0; i < 3; ++i)
    {
      hopThreshold[i] = 0.5 + commitDistance / cageLattice.axisLength(static_cast<int>(i));
    }
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      if (j == taggedIndex) continue;  // committed to A or B only once it has reached the commit distance
      const Molecule& molecule = system.moleculeData[j];
      const Component& component = system.components[molecule.componentId];
      initialCommittedCage[j] = cageLattice.cageOf(systemAtoms[molecule.atomIndex + component.startingBead].position);
    }
  }
  // the window between two cages of the column along the normal: 0 = W(A'|A), 1 = W(A|B), 2 = W(B|B'); the
  // cages A' and B' coincide when the box is three cells long along the normal, in which case the three
  // windows are all the windows of the column
  auto columnWindow = [&](std::size_t from, std::size_t to) -> std::optional<std::size_t>
  {
    auto pair = [&](std::size_t c1, std::size_t c2) { return (from == c1 && to == c2) || (from == c2 && to == c1); };
    if (pair(cageA, cageB)) return 1uz;
    if (pair(cageAPrime, cageA)) return 0uz;
    if (pair(cageB, cageC)) return 2uz;
    return std::nullopt;
  };

  struct TrajectoryOutcome
  {
    double meanDrift{0.0};
    double maxDrift{0.0};
    bool aborted{false};  // interrupted by a termination signal; the pair is discarded
    bool longJump{false};
    bool committed{false};
    bool eventTerminated{false};
    std::size_t terminationStep{0};  // the last integrated step
    std::optional<EventChannel> terminatingChannel{};
    std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)> events{};  // after commitment
    double finalSelf{0.0};    // net cage count of the tagged molecule at termination
    double finalLocalQ{0.0};  // Q_local: signed axial crossings of the other molecules in the column
    // decomposition in the frame of the tagged hop, [phase][channel] with phase 0 before and 1 after the
    // commitment of the tagged molecule: signed contribution to Q_local (others only), number of crossings
    // (tagged included) and the sum of their times (since launch, resp. since commitment) [ps]
    std::size_t commitStep{0};
    std::array<std::array<double, static_cast<std::size_t>(EventChannel::Count)>, 2> frameQ{};
    std::array<std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)>, 2> eventCount{};
    std::array<std::array<double, static_cast<std::size_t>(EventChannel::Count)>, 2> eventTime{};
  };

  // integrates one trajectory from the current working copies and fills the delta arrays; reports the
  // energy drift (meaningful for NVE only), whether the tagged molecule left the pair of cages A/B
  // ('long jump') and whether it reached the commit distance; in event mode the trajectory ends at the
  // n-th local axial hop event after commitment and the stored values are held from then on
  auto integrateTrajectory = [&](double direction, TrajectoryOutcome& outcome)
  {
    std::span<GroupState> groupSpan(groupData);
    const std::span<Molecule> moleculeSpan(moleculeData);
    const std::span<Atom> atomSpan(atoms);
    const std::span<AtomDynamics> dynamicsSpan(dynamics);

    // gradients on the working copies (the first half-kick of velocity Verlet uses them) and the reference
    // conserved energy
    RunningEnergy referenceEnergy = Integrators::updateGradients(
        moleculeSpan, atomSpan, dynamicsSpan, frameworkAtoms, system.forceField, system.simulationBox,
        system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
        system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent,
        system.framework, {}, &system.crossLinks);
    Integrators::updateCenterOfMassAndQuaternionGradients(moleculeSpan, atomSpan, dynamicsSpan, system.components,
                                                          groupSpan);
    referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        moleculeSpan, atomSpan, dynamicsSpan, system.components, system.framework, frameworkAtoms, {},
        &system.forceField, groupSpan);
    referenceEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeSpan, system.components, groupSpan);
    const double conservedReference = referenceEnergy.conservedEnergy();

    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      const Molecule& molecule = moleculeData[j];
      const Component& component = system.components[molecule.componentId];
      previousPosition[j] = atoms[molecule.atomIndex + component.startingBead].position;
      unwrappedQ[j] = (j == taggedIndex) ? 0.0 : barrier.projection(system.simulationBox, previousPosition[j]);
    }

    // t = 0: the tagged bead is on the surface, moving into B (direction +1) or A (direction -1)
    std::fill(deltaCollective.begin(), deltaCollective.end(), 0.0);
    deltaSelf[0] = 0.5 * direction;
    deltaSelfTwoState[0] = 0.5 * direction;
    deltaSelfCommitted[0] = 0.5 * direction;
    outcome = TrajectoryOutcome{};
    bool& committed = outcome.committed;
    double committedValue{0.0};

    // event mode: the frame of the tagged hop (mirrored for the time-reversed trajectory)
    committedCage = initialCommittedCage;
    double localQ{0.0};
    std::size_t eventsAfterCommitment{0};
    auto channelOf = [&](std::size_t j, std::size_t window, int alongNormal) -> EventChannel
    {
      // in the frame of the tagged hop: direction +1 is the direction of the tagged molecule's initial
      // velocity, window 0 is W(A'|A), 1 is W(A|B), 2 is W(B|B') with A the cage it started from
      const int frameDirection = (direction > 0.0) ? alongNormal : -alongNormal;
      const std::size_t frameWindow = (direction > 0.0) ? window : 2uz - window;
      if (j == taggedIndex)
      {
        return (frameDirection > 0) ? EventChannel::TaggedForwardHop : EventChannel::TaggedBackHop;
      }
      switch (frameWindow)
      {
        case 0:
          return (frameDirection > 0) ? EventChannel::OtherAPrimeToA : EventChannel::OtherAtoAPrime;
        case 1:
          return (frameDirection > 0) ? EventChannel::OtherAtoB : EventChannel::OtherBtoA;
        default:
          return (frameDirection > 0) ? EventChannel::OtherBtoBp : EventChannel::OtherBptoB;
      }
    };

    // stores the observables of the current state at time index t
    auto storeSample = [&](std::size_t t)
    {
      const double qTagged = unwrappedQ[taggedIndex];
      // the tagged molecule: net cage count from the dividing surface (B = +1/2, A = -1/2, further cages
      // +-3/2, ...), the two-state (side of the plane) variant of Bennett-Chandler, and the committed
      // variant (side of the plane, frozen once the commit distance has been reached)
      deltaSelf[t] = static_cast<double>(barrier.cageIndex(qTagged)) + 0.5;
      deltaSelfTwoState[t] = (qTagged > 0.0) ? 0.5 : -0.5;
      deltaSelfCommitted[t] = committed ? committedValue : deltaSelfTwoState[t];
      for (std::size_t j = 0; j < numberOfMolecules; ++j)
      {
        if (j == taggedIndex) continue;
        const double delta = static_cast<double>(barrier.cageIndex(unwrappedQ[j]) - initialCage[j]);
        if (delta != 0.0)
        {
          deltaCollective[accumulator.collectiveIndex(bucketOfMolecule[j], moleculeData[j].componentId, t)] += delta;
        }
      }
    };

    double driftSum{0.0};
    std::size_t step{1};
    for (; step <= numberOfTrajectorySteps; ++step)
    {
      // a trajectory can take hours at a long cap: poll the termination signal inside the integration too
      if ((step & 255uz) == 0uz && GracefulShutdown::requested())
      {
        outcome.aborted = true;
        return;
      }

      const RunningEnergy energies = Integrators::velocityVerlet(
          moleculeSpan, atomSpan, dynamicsSpan, system.components, system.timeStep, noThermostat, frameworkAtoms,
          system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
          system.numberOfMoleculesPerComponent, system.framework, {}, groupSpan, {}, &system.crossLinks);

      if (bathType != DcTSTSettings::TrajectoryThermostat::None)
      {
        applyBath(rng, bathType, settings.bathRate, system.timeStep, system.temperature, moleculeSpan, dynamicsSpan,
                  groupSpan, system.components);
      }
      else
      {
        const double drift = std::abs((energies.conservedEnergy() - conservedReference) / conservedReference);
        driftSum += drift;
        outcome.maxDrift = std::max(outcome.maxDrift, drift);
      }

      for (std::size_t j = 0; j < numberOfMolecules; ++j)
      {
        const Molecule& molecule = moleculeData[j];
        const Component& component = system.components[molecule.componentId];
        const double3 position = atoms[molecule.atomIndex + component.startingBead].position;
        unwrappedQ[j] += barrier.displacement(system.simulationBox, previousPosition[j], position);
        previousPosition[j] = position;
      }

      // commitment and long-jump detection of the tagged molecule (checked every step)
      const double qTagged = unwrappedQ[taggedIndex];
      if (!committed && std::abs(qTagged) >= commitDistance)
      {
        committed = true;
        committedValue = (qTagged > 0.0) ? 0.5 : -0.5;
        outcome.commitStep = step;
        if (eventMode) committedCage[taggedIndex] = (qTagged > 0.0) ? cageB : cageA;
      }
      if (qTagged >= barrier.lambda || qTagged < -barrier.lambda) outcome.longJump = true;

      // event mode: hop detection of every molecule with a committed cage; a hop along the normal between
      // two cages of the column is a local axial event
      bool terminate{false};
      if (eventMode)
      {
        for (std::size_t j = 0; j < numberOfMolecules && !terminate; ++j)
        {
          const std::size_t cage = committedCage[j];
          if (cage == noCage) continue;
          const double3 f = cageLattice.fractionalOffset(system.simulationBox, previousPosition[j], cage);
          for (std::size_t i = 0; i < 3; ++i)
          {
            const int hop = (f[i] >= hopThreshold[i]) ? 1 : (f[i] <= -hopThreshold[i]) ? -1 : 0;
            if (hop == 0) continue;
            const std::size_t newCage = cageLattice.neighborAlongAxis(cage, static_cast<int>(i), hop);
            if (static_cast<int>(i) == axis.axis)
            {
              const std::optional<std::size_t> window = columnWindow(cage, newCage);
              if (window.has_value())
              {
                const int alongNormal = hop * axis.sign;
                const EventChannel channel = channelOf(j, window.value(), alongNormal);
                const std::size_t channelIndex = static_cast<std::size_t>(channel);
                const std::size_t phase = committed ? 1uz : 0uz;
                const int frameDirection = (direction > 0.0) ? alongNormal : -alongNormal;
                if (j != taggedIndex)
                {
                  localQ += static_cast<double>(alongNormal);
                  outcome.frameQ[phase][channelIndex] += static_cast<double>(frameDirection);
                }
                ++outcome.eventCount[phase][channelIndex];
                outcome.eventTime[phase][channelIndex] +=
                    static_cast<double>(step - (committed ? outcome.commitStep : 0uz)) * system.timeStep;
                if (committed)
                {
                  ++outcome.events[static_cast<std::size_t>(channel)];
                  if (++eventsAfterCommitment >= settings.terminationEvents)
                  {
                    terminate = true;
                    outcome.terminatingChannel = channel;
                  }
                }
              }
            }
            committedCage[j] = newCage;
            break;  // at most one hop per molecule per step
          }
        }
      }

      if (step % sampleEvery == 0uz) storeSample(step / sampleEvery);
      if (terminate) break;
    }

    const std::size_t lastStep = std::min(step, numberOfTrajectorySteps);
    outcome.terminationStep = lastStep;
    outcome.eventTerminated = (step <= numberOfTrajectorySteps);
    outcome.meanDrift = driftSum / static_cast<double>(lastStep);
    outcome.finalSelf = static_cast<double>(barrier.cageIndex(unwrappedQ[taggedIndex])) + 0.5;
    outcome.finalLocalQ = localQ;
    // hold the final values for the remaining stored time points (the index lastStep / sampleEvery has been
    // stored inside the loop when lastStep is a multiple of sampleEvery, and otherwise holds the sample of the
    // last stored time before termination)
    for (std::size_t t = lastStep / sampleEvery + 1uz; t < numberOfTimes; ++t)
    {
      storeSample(t);
    }
  };

  std::vector<double> forwardSelf(numberOfTimes);
  std::vector<double> forwardSelfTwoState(numberOfTimes);
  std::vector<double> forwardSelfCommitted(numberOfTimes);
  std::vector<double> forwardCollective(deltaCollective.size());

  for (std::size_t k = 0; k < settings.trajectoriesPerConfiguration; ++k)
  {
    // fresh copies of the Monte Carlo configuration with Maxwell-Boltzmann velocities
    moleculeData.assign(system.moleculeData.begin(), system.moleculeData.end());
    atoms.assign(systemAtoms.begin(), systemAtoms.end());
    dynamics.assign(systemDynamics.begin(), systemDynamics.end());
    groupData = initialGroupData;
    std::span<GroupState> groupSpan(groupData);

    Integrators::createCartesianPositions(moleculeData, atoms, system.components, groupSpan);
    Integrators::initializeVelocities(rng, moleculeData, atoms, dynamics, system.components, system.temperature,
                                      system.framework, frameworkAtoms, {}, &system.forceField, groupSpan);
    if (system.numberOfFrameworkAtoms == 0uz)
    {
      Integrators::removeCenterOfMassVelocityDrift(moleculeData, atoms, dynamics, system.components, system.framework,
                                                   frameworkAtoms, {}, &system.forceField, groupSpan);
    }

    // orient the initial velocity of the reaction bead towards B
    const double3 velocity =
        beadVelocity(taggedComponent, moleculeData[taggedIndex], dynamics, taggedComponent.startingBead);
    double qdot = double3::dot(velocity, barrier.normal);
    if (qdot < 0.0)
    {
      reverseVelocities(moleculeData, dynamics, groupSpan);
      qdot = -qdot;
    }

    // the initial state (positions and momenta) of the pair
    const std::vector<Molecule> initialMoleculeData = moleculeData;
    const std::vector<Atom> initialAtoms = atoms;
    const std::vector<AtomDynamics> initialDynamics = dynamics;
    const std::vector<GroupState> initialGroups = groupData;

    // forward trajectory
    TrajectoryOutcome forward;
    integrateTrajectory(+1.0, forward);
    if (forward.aborted) return;
    forwardSelf = deltaSelf;
    forwardSelfTwoState = deltaSelfTwoState;
    forwardSelfCommitted = deltaSelfCommitted;
    forwardCollective = deltaCollective;

    // time-reversed trajectory from the same phase point
    moleculeData = initialMoleculeData;
    atoms = initialAtoms;
    dynamics = initialDynamics;
    groupData = initialGroups;
    groupSpan = std::span<GroupState>(groupData);
    reverseVelocities(moleculeData, dynamics, groupSpan);
    TrajectoryOutcome backward;
    integrateTrajectory(-1.0, backward);
    if (backward.aborted) return;

    // accumulate |qdot| (Delta^+ - Delta^-)
    accumulator.speed[block] += qdot;
    accumulator.sumSpeedSquared += qdot * qdot;
    ++accumulator.numberOfTrajectoryPairs[block];
    accumulator.numberOfTrajectories += 2uz;
    accumulator.sumRelativeEnergyDrift += forward.meanDrift + backward.meanDrift;
    accumulator.maximumRelativeEnergyDrift =
        std::max({accumulator.maximumRelativeEnergyDrift, forward.maxDrift, backward.maxDrift});
    accumulator.numberOfLongJumpTrajectories += (forward.longJump ? 1uz : 0uz) + (backward.longJump ? 1uz : 0uz);
    accumulator.numberOfCommittedTrajectories += (forward.committed ? 1uz : 0uz) + (backward.committed ? 1uz : 0uz);

    std::vector<double>& self = accumulator.self[block];
    std::vector<double>& collective = accumulator.collective[block];
    std::vector<double>& selfTwoState = accumulator.selfTwoState[block];
    std::vector<double>& selfCommitted = accumulator.selfCommitted[block];
    for (std::size_t t = 0; t < numberOfTimes; ++t)
    {
      self[t] += qdot * (forwardSelf[t] - deltaSelf[t]);
      selfTwoState[t] += qdot * (forwardSelfTwoState[t] - deltaSelfTwoState[t]);
      selfCommitted[t] += qdot * (forwardSelfCommitted[t] - deltaSelfCommitted[t]);
    }
    for (std::size_t i = 0; i < collective.size(); ++i)
    {
      collective[i] += qdot * (forwardCollective[i] - deltaCollective[i]);
    }

    // event-terminated estimator: values at termination
    accumulator.selfEvent[block] += qdot * (forward.finalSelf - backward.finalSelf);
    accumulator.collectiveEvent[block] +=
        qdot * ((forward.finalSelf + forward.finalLocalQ) - (backward.finalSelf + backward.finalLocalQ));
    for (std::size_t d = 0; const TrajectoryOutcome* outcome : {&forward, &backward})
    {
      accumulator.integratedTimeSum += static_cast<double>(outcome->terminationStep) * system.timeStep;
      if (outcome->eventTerminated)
      {
        ++accumulator.numberOfEventTerminatedTrajectories;
        const double time = static_cast<double>(outcome->terminationStep) * system.timeStep;
        accumulator.terminationTimeSum += time;
        accumulator.terminationTimeMaximum = std::max(accumulator.terminationTimeMaximum, time);
        ++accumulator.terminatingEvents[d][static_cast<std::size_t>(outcome->terminatingChannel.value())];

        // conditioned estimator: this member alone, with the sign of its launch direction
        const double sign = (d == 0uz) ? 1.0 : -1.0;
        accumulator.speedEventTerminated[block] += qdot;
        accumulator.selfEventConditioned[block] += qdot * sign * outcome->finalSelf;
        accumulator.collectiveEventConditioned[block] += qdot * sign * (outcome->finalSelf + outcome->finalLocalQ);
      }
      else if (outcome->committed)
      {
        ++accumulator.numberOfCappedTrajectories;
      }
      else
      {
        ++accumulator.numberOfUncommittedTrajectories;
      }
      for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
      {
        accumulator.allEvents[d][c] += outcome->events[c];
      }
      // decomposition of Q_local by phase and channel (frame of the tagged hop: the forward and the
      // time-reversed member add with the same sign)
      for (std::size_t phase = 0; phase < 2; ++phase)
      {
        double phaseQ{0.0};
        for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
        {
          phaseQ += outcome->frameQ[phase][c];
          accumulator.localQChannel[phase][c] += qdot * outcome->frameQ[phase][c];
          accumulator.channelEvents[phase][c] += outcome->eventCount[phase][c];
          accumulator.channelTimeSum[phase][c] += outcome->eventTime[phase][c];
        }
        if (phase == 0uz)
          accumulator.collectiveEventBeforeCommitment[block] += qdot * phaseQ;
        else
          accumulator.collectiveEventAfterCommitment[block] += qdot * phaseQ;
      }
      ++d;
    }

    // per-pair record (flushed, so that a killed run keeps every completed pair)
    if (settings.writePairRecords && replicaId < pairStreams.size())
    {
      auto plateauAverage = [&](const std::vector<double>& forwardValues, const std::vector<double>& backwardValues)
      {
        double sum{0.0};
        for (std::size_t t : plateauIndices) sum += forwardValues[t] - backwardValues[t];
        return sum / static_cast<double>(std::max(1uz, plateauIndices.size()));
      };
      double collectivePlateau{0.0};
      for (std::size_t t : plateauIndices)
      {
        for (std::size_t s = 0; s < numberOfShells; ++s)
        {
          for (std::size_t c = 0; c < numberOfComponents; ++c)
          {
            const std::size_t index = accumulator.collectiveIndex(s, c, t);
            collectivePlateau += forwardCollective[index] - deltaCollective[index];
          }
        }
      }
      collectivePlateau /= static_cast<double>(std::max(1uz, plateauIndices.size()));

      auto status = [](const TrajectoryOutcome& outcome) -> int
      { return outcome.eventTerminated ? 0 : (outcome.committed ? 1 : 2); };
      auto channel = [](const TrajectoryOutcome& outcome) -> int
      { return outcome.terminatingChannel.has_value() ? static_cast<int>(outcome.terminatingChannel.value()) : -1; };
      auto phaseSum = [](const TrajectoryOutcome& outcome, std::size_t phase)
      {
        double sum{0.0};
        for (double q : outcome.frameQ[phase]) sum += q;
        return sum;
      };
      auto time = [&](std::size_t steps) { return static_cast<double>(steps) * system.timeStep; };

      std::ostream pairStream(pairStreams[replicaId].rdbuf());
      std::print(pairStream,
                 "{} {} {:.8g} {:g} {:g} {:g} {:g} {} {} {:.6g} {:.6g} {} {} {:.6g} {:.6g} {:g} {:g} {:g} {:g} {} {} "
                 "{:.4e} {:.4e} {:.8g} {:.8g} {:.8g} {:.8g} {:.8g}\n",
                 productionCycle, block, qdot, forward.finalSelf, backward.finalSelf, forward.finalLocalQ,
                 backward.finalLocalQ, status(forward), status(backward), time(forward.terminationStep),
                 time(backward.terminationStep), channel(forward), channel(backward),
                 forward.committed ? time(forward.commitStep) : 0.0, backward.committed ? time(backward.commitStep) : 0.0,
                 phaseSum(forward, 0), phaseSum(forward, 1), phaseSum(backward, 0), phaseSum(backward, 1),
                 forward.longJump ? 1 : 0, backward.longJump ? 1 : 0, forward.meanDrift, backward.meanDrift,
                 plateauAverage(forwardSelf, deltaSelf), plateauAverage(forwardSelfTwoState, deltaSelfTwoState),
                 plateauAverage(forwardSelfCommitted, deltaSelfCommitted),
                 forwardSelfCommitted.back() - deltaSelfCommitted.back(), collectivePlateau);
      std::flush(pairStream);
    }
  }
}

// =================================================================================================
// Plateau window, per-pair records and accumulator snapshots
// =================================================================================================

std::vector<std::size_t> DynamicallyCorrectedTST::computePlateauIndices() const
{
  const double timeStep = systems.front().timeStep;
  double plateauBegin = 0.5 * settings.trajectoryTime;
  double plateauEnd = settings.trajectoryTime;
  if (settings.plateauWindow.has_value())
  {
    plateauBegin = settings.plateauWindow->first;
    plateauEnd = settings.plateauWindow->second;
  }
  std::vector<std::size_t> indices;
  for (std::size_t t = 0; t < numberOfTimes; ++t)
  {
    const double time = static_cast<double>(t * settings.sampleTrajectoryEvery) * timeStep;
    if (time >= plateauBegin - 1e-12 && time <= plateauEnd + 1e-12) indices.push_back(t);
  }
  if (indices.empty()) indices.push_back(numberOfTimes - 1uz);
  return indices;
}

void DynamicallyCorrectedTST::writePairRecordHeader(std::size_t replicaId)
{
  const System& system = systems[replicaId];
  const Component& tagged = system.components[settings.taggedComponent];
  const double commitDistance = (settings.commitDistance > 0.0) ? settings.commitDistance : 0.5 * barrier.lambda;
  const double plateauBegin = static_cast<double>(plateauIndices.front() * settings.sampleTrajectoryEvery) * system.timeStep;
  const double plateauEnd = static_cast<double>(plateauIndices.back() * settings.sampleTrajectoryEvery) * system.timeStep;

  std::ostream file(pairStreams[replicaId].rdbuf());
  std::print(file, "# dcTST per-pair records of component '{}', replica {} of {} (seed {})\n", tagged.name, replicaId,
             numberOfReplicas, randoms[replicaId].seed);
  std::print(file, "# one line per completed trajectory pair (forward '+' and time-reversed '-'), written after the pair;\n");
  std::print(file, "# every estimator is a ratio of sums over the pairs, so files of different replicas and runs can be "
                   "concatenated\n");
  std::print(file, "# trajectory cap {} ps, time step {} ps, commit distance {:.4f} A, plateau window [{}, {}] ps, "
                   "termination: {}\n",
             settings.trajectoryTime, system.timeStep, commitDistance, plateauBegin, plateauEnd,
             (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event)
                 ? std::format("event ({} local axial event(s) after commitment)", settings.terminationEvents)
                 : std::string("fixed time"));
  std::print(file, "# column 1: production cycle of the transmission stage\n");
  std::print(file, "# column 2: block\n");
  std::print(file, "# column 3: |qdot| [A/ps], initial speed of the reaction bead along the normal (the weight)\n");
  std::print(file, "# column 4: Delta_S+ (final), net cage count of the tagged molecule at termination (+1/2 in B, -1/2 in A)\n");
  std::print(file, "# column 5: Delta_S- (final)\n");
  std::print(file, "# column 6: Q_local+, signed axial window crossings of the other molecules in the column (lab frame)\n");
  std::print(file, "# column 7: Q_local-\n");
  std::print(file, "# column 8: status+ (0 ended by an event, 1 capped after commitment, 2 never committed)\n");
  std::print(file, "# column 9: status-\n");
  std::print(file, "# column 10: termination time+ [ps]\n");
  std::print(file, "# column 11: termination time- [ps]\n");
  std::print(file, "# column 12: terminating channel+ (-1 none; 0 tagged B->A, 1 tagged B->B', 2 other A->B, 3 other B->A, "
                   "4 other A->A', 5 other A'->A, 6 other B->B', 7 other B'->B; frame of the tagged hop)\n");
  std::print(file, "# column 13: terminating channel-\n");
  std::print(file, "# column 14: commit time+ [ps] (0 if never committed)\n");
  std::print(file, "# column 15: commit time-\n");
  std::print(file, "# column 16: Q_before+, crossings of the other molecules before the commitment (frame of the tagged hop)\n");
  std::print(file, "# column 17: Q_after+, after the commitment\n");
  std::print(file, "# column 18: Q_before-\n");
  std::print(file, "# column 19: Q_after-\n");
  std::print(file, "# column 20: long jump+ (tagged molecule left cages A/B)\n");
  std::print(file, "# column 21: long jump-\n");
  std::print(file, "# column 22: mean |relative energy drift|+\n");
  std::print(file, "# column 23: mean |relative energy drift|-\n");
  std::print(file, "# column 24: plateau average of Delta_S+ - Delta_S- (net cage count)\n");
  std::print(file, "# column 25: plateau average of Delta_S+ - Delta_S- (two-state, side of the plane)\n");
  std::print(file, "# column 26: plateau average of Delta_S+ - Delta_S- (committed)\n");
  std::print(file, "# column 27: Delta_S+ - Delta_S- (committed) at the last stored time\n");
  std::print(file, "# column 28: plateau average of the sum over all other molecules of Delta_j+ - Delta_j- (kappa_C - kappa_S, all molecules)\n");
  std::print(file, "# estimators (w = column 3, sums over the pairs):\n");
  std::print(file, "#   kappa_S^event = sum w (c4 - c5) / sum w;   kappa_C^event = sum w ((c4 + c6) - (c5 + c7)) / sum w\n");
  std::print(file, "#   kappa_S (plateau) = sum w c24 / sum w;   kappa_S committed (plateau) = sum w c26 / sum w;   "
                   "kappa_C (plateau, all) = sum w (c24 + c28) / sum w\n");
  std::print(file, "#   conditioned (event-ended members only, s = +1 for '+', -1 for '-'): kappa^cond = 2 sum w s Delta / sum w\n");
  std::print(file, "#   epsilon_before = sum w (c16 + c18) / sum w;   epsilon_after = sum w (c17 + c19) / sum w\n");
  std::flush(file);
}

void DynamicallyCorrectedTST::writeAccumulatorSnapshot(std::size_t replicaId, bool final) const
{
  const System& system = systems[replicaId];
  const DcTSTTransmissionAccumulator& accumulator = transmission[replicaId];

  nlohmann::json json;
  json["replica"] = replicaId;
  json["numberOfReplicas"] = numberOfReplicas;
  json["seed"] = randoms[replicaId].seed;
  json["final"] = final;
  json["transmissionCyclesCompleted"] = completedTransmissionCycles[replicaId];
  json["numberOfTransmissionCycles"] = settings.numberOfTransmissionCycles;
  json["numberOfBlocks"] = accumulator.numberOfBlocks;
  json["numberOfShells"] = accumulator.numberOfShells;
  json["numberOfComponents"] = accumulator.numberOfComponents;
  json["numberOfTimes"] = accumulator.numberOfTimes;
  json["sampleTrajectoryEvery"] = settings.sampleTrajectoryEvery;
  json["timeStep"] = system.timeStep;
  json["trajectoryTime"] = settings.trajectoryTime;
  json["plateauIndices"] = {plateauIndices.front(), plateauIndices.back()};
  json["collectiveIndex"] = "(shell * numberOfComponents + component) * numberOfTimes + time";

  json["speed"] = accumulator.speed;
  json["numberOfTrajectoryPairs"] = accumulator.numberOfTrajectoryPairs;
  json["self"] = accumulator.self;
  json["selfTwoState"] = accumulator.selfTwoState;
  json["selfCommitted"] = accumulator.selfCommitted;
  json["collective"] = accumulator.collective;
  json["sumSpeedSquared"] = accumulator.sumSpeedSquared;
  json["sumRelativeEnergyDrift"] = accumulator.sumRelativeEnergyDrift;
  json["maximumRelativeEnergyDrift"] = accumulator.maximumRelativeEnergyDrift;
  json["numberOfTrajectories"] = accumulator.numberOfTrajectories;
  json["numberOfLongJumpTrajectories"] = accumulator.numberOfLongJumpTrajectories;
  json["numberOfCommittedTrajectories"] = accumulator.numberOfCommittedTrajectories;
  json["selfEvent"] = accumulator.selfEvent;
  json["collectiveEvent"] = accumulator.collectiveEvent;
  json["speedEventTerminated"] = accumulator.speedEventTerminated;
  json["selfEventConditioned"] = accumulator.selfEventConditioned;
  json["collectiveEventConditioned"] = accumulator.collectiveEventConditioned;
  json["collectiveEventBeforeCommitment"] = accumulator.collectiveEventBeforeCommitment;
  json["collectiveEventAfterCommitment"] = accumulator.collectiveEventAfterCommitment;
  json["localQChannel"] = accumulator.localQChannel;
  json["channelTimeSum"] = accumulator.channelTimeSum;
  json["channelEvents"] = accumulator.channelEvents;
  json["terminatingEvents"] = accumulator.terminatingEvents;
  json["allEvents"] = accumulator.allEvents;
  json["numberOfEventTerminatedTrajectories"] = accumulator.numberOfEventTerminatedTrajectories;
  json["numberOfCappedTrajectories"] = accumulator.numberOfCappedTrajectories;
  json["numberOfUncommittedTrajectories"] = accumulator.numberOfUncommittedTrajectories;
  json["terminationTimeSum"] = accumulator.terminationTimeSum;
  json["terminationTimeMaximum"] = accumulator.terminationTimeMaximum;
  json["integratedTimeSum"] = accumulator.integratedTimeSum;

  // written to a temporary file and renamed, so that the snapshot is always complete
  const std::string fileName =
      std::format("dctst/accumulator_{}_{}.r{}.json", system.temperature, system.input_pressure, replicaId);
  const std::string temporaryName = fileName + ".tmp";
  {
    std::ofstream file(temporaryName, std::ios::out | std::ios::trunc);
    file << json.dump();
  }
  std::error_code error;
  std::filesystem::rename(temporaryName, fileName, error);
}

// =================================================================================================
// Analysis
// =================================================================================================

DynamicallyCorrectedTST::TransmissionResult DynamicallyCorrectedTST::analyzeTransmission() const
{
  TransmissionResult result;
  const DcTSTTransmissionAccumulator& front = transmission.front();
  const std::size_t numberOfShells = front.numberOfShells;
  const std::size_t numberOfComponents = front.numberOfComponents;
  const System& system = systems.front();

  result.time.resize(numberOfTimes);
  for (std::size_t t = 0; t < numberOfTimes; ++t)
  {
    result.time[t] = static_cast<double>(t * settings.sampleTrajectoryEvery) * system.timeStep;
  }

  // plateau window
  double plateauBegin = 0.5 * settings.trajectoryTime;
  double plateauEnd = settings.trajectoryTime;
  if (settings.plateauWindow.has_value())
  {
    plateauBegin = settings.plateauWindow->first;
    plateauEnd = settings.plateauWindow->second;
  }
  // plateauIndices (member) was set in setup()

  // pooled sums and per-block (replica x block) transmission functions
  double speed{0.0};
  std::vector<double> self(numberOfTimes, 0.0);
  std::vector<double> selfTwoState(numberOfTimes, 0.0);
  std::vector<double> selfCommitted(numberOfTimes, 0.0);
  std::vector<double> collective(numberOfShells * numberOfComponents * numberOfTimes, 0.0);
  std::vector<double> blockSpeeds;
  std::vector<std::vector<double>> blockSelf;
  std::vector<std::vector<double>> blockSelfCommitted;
  std::vector<std::vector<std::vector<double>>> blockCollective;  // [block][shell][time], cumulative, kappa_C
  std::size_t numberOfTrajectories{0};
  std::size_t numberOfLongJumps{0};
  std::size_t numberOfCommitted{0};
  double selfEvent{0.0};
  double collectiveEvent{0.0};
  std::vector<double> blockSelfEvent;
  std::vector<double> blockCollectiveEvent;
  std::vector<double> blockEventExcess;           // (kappa_C^event - kappa_S^event) per block
  std::vector<double> blockEventBackCorrelation;  // (kappa_S^event - kappa_S^committed) per block
  double speedEventTerminated{0.0};
  double selfEventConditioned{0.0};
  double collectiveEventConditioned{0.0};
  std::vector<double> blockSelfEventConditioned;
  std::vector<double> blockCollectiveEventConditioned;
  double collectiveEventBefore{0.0};
  double collectiveEventAfter{0.0};
  std::vector<double> blockCollectiveEventBefore;
  std::vector<double> blockCollectiveEventAfter;
  std::array<std::array<double, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
      localQChannel{};
  std::array<std::array<double, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
      channelTimeSum{};
  double terminationTimeSum{0.0};
  double integratedTimeSum{0.0};
  for (const DcTSTTransmissionAccumulator& accumulator : transmission)
  {
    result.numberOfTrajectoryPairs +=
        std::accumulate(accumulator.numberOfTrajectoryPairs.begin(), accumulator.numberOfTrajectoryPairs.end(), 0uz);
    result.meanRelativeEnergyDrift += accumulator.sumRelativeEnergyDrift;
    result.maximumRelativeEnergyDrift =
        std::max(result.maximumRelativeEnergyDrift, accumulator.maximumRelativeEnergyDrift);
    numberOfTrajectories += accumulator.numberOfTrajectories;
    numberOfLongJumps += accumulator.numberOfLongJumpTrajectories;
    numberOfCommitted += accumulator.numberOfCommittedTrajectories;
    result.numberOfEventTerminatedTrajectories += accumulator.numberOfEventTerminatedTrajectories;
    result.numberOfCappedTrajectories += accumulator.numberOfCappedTrajectories;
    result.numberOfUncommittedTrajectories += accumulator.numberOfUncommittedTrajectories;
    terminationTimeSum += accumulator.terminationTimeSum;
    integratedTimeSum += accumulator.integratedTimeSum;
    result.maximumTerminationTime = std::max(result.maximumTerminationTime, accumulator.terminationTimeMaximum);
    for (std::size_t d = 0; d < 2; ++d)
    {
      for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
      {
        result.terminatingEvents[d][c] += accumulator.terminatingEvents[d][c];
        result.allEvents[d][c] += accumulator.allEvents[d][c];
        // d doubles as the phase index here (both run over 0, 1)
        localQChannel[d][c] += accumulator.localQChannel[d][c];
        channelTimeSum[d][c] += accumulator.channelTimeSum[d][c];
        result.channelEvents[d][c] += accumulator.channelEvents[d][c];
      }
    }
    for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
    {
      speed += accumulator.speed[b];
      selfEvent += accumulator.selfEvent[b];
      collectiveEvent += accumulator.collectiveEvent[b];
      if (accumulator.numberOfTrajectoryPairs[b] > 0uz && accumulator.speed[b] > 0.0)
      {
        blockSelfEvent.push_back(accumulator.selfEvent[b] / accumulator.speed[b]);
        blockCollectiveEvent.push_back(accumulator.collectiveEvent[b] / accumulator.speed[b]);
        // the differences are far less noisy than the kappas themselves (same trajectories)
        blockEventExcess.push_back((accumulator.collectiveEvent[b] - accumulator.selfEvent[b]) /
                                   accumulator.speed[b]);
        blockEventBackCorrelation.push_back((accumulator.selfEvent[b] - accumulator.selfCommitted[b].back()) /
                                            accumulator.speed[b]);
        blockCollectiveEventBefore.push_back(accumulator.collectiveEventBeforeCommitment[b] / accumulator.speed[b]);
        blockCollectiveEventAfter.push_back(accumulator.collectiveEventAfterCommitment[b] / accumulator.speed[b]);
      }
      collectiveEventBefore += accumulator.collectiveEventBeforeCommitment[b];
      collectiveEventAfter += accumulator.collectiveEventAfterCommitment[b];
      speedEventTerminated += accumulator.speedEventTerminated[b];
      selfEventConditioned += accumulator.selfEventConditioned[b];
      collectiveEventConditioned += accumulator.collectiveEventConditioned[b];
      if (accumulator.speedEventTerminated[b] > 0.0)
      {
        // per trajectory member: kappa = 2 sum |qdot| s Delta / sum |qdot| (see the accumulator)
        blockSelfEventConditioned.push_back(2.0 * accumulator.selfEventConditioned[b] /
                                            accumulator.speedEventTerminated[b]);
        blockCollectiveEventConditioned.push_back(2.0 * accumulator.collectiveEventConditioned[b] /
                                                  accumulator.speedEventTerminated[b]);
      }
      for (std::size_t t = 0; t < numberOfTimes; ++t)
      {
        self[t] += accumulator.self[b][t];
        selfTwoState[t] += accumulator.selfTwoState[b][t];
        selfCommitted[t] += accumulator.selfCommitted[b][t];
      }
      for (std::size_t i = 0; i < collective.size(); ++i) collective[i] += accumulator.collective[b][i];

      if (accumulator.numberOfTrajectoryPairs[b] == 0uz || accumulator.speed[b] <= 0.0) continue;
      blockSpeeds.push_back(accumulator.speed[b] / static_cast<double>(accumulator.numberOfTrajectoryPairs[b]));
      std::vector<double> kappaSelf(numberOfTimes);
      std::vector<double> kappaCommitted(numberOfTimes);
      for (std::size_t t = 0; t < numberOfTimes; ++t)
      {
        kappaSelf[t] = accumulator.self[b][t] / accumulator.speed[b];
        kappaCommitted[t] = accumulator.selfCommitted[b][t] / accumulator.speed[b];
      }
      blockSelf.push_back(kappaSelf);
      blockSelfCommitted.push_back(kappaCommitted);

      std::vector<std::vector<double>> kappaCollective(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
      std::vector<double> cumulative(numberOfTimes, 0.0);
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        for (std::size_t c = 0; c < numberOfComponents; ++c)
        {
          for (std::size_t t = 0; t < numberOfTimes; ++t)
          {
            cumulative[t] += accumulator.collective[b][accumulator.collectiveIndex(s, c, t)];
          }
        }
        for (std::size_t t = 0; t < numberOfTimes; ++t)
        {
          kappaCollective[s][t] = kappaSelf[t] + cumulative[t] / accumulator.speed[b];
        }
      }
      blockCollective.push_back(kappaCollective);
    }
  }
  result.meanRelativeEnergyDrift /= static_cast<double>(std::max(2uz, 2uz * result.numberOfTrajectoryPairs));
  if (numberOfTrajectories > 0uz)
  {
    result.longJumpFraction = static_cast<double>(numberOfLongJumps) / static_cast<double>(numberOfTrajectories);
    result.committedFraction = static_cast<double>(numberOfCommitted) / static_cast<double>(numberOfTrajectories);
  }

  result.meanSpeed =
      (result.numberOfTrajectoryPairs > 0uz) ? speed / static_cast<double>(result.numberOfTrajectoryPairs) : 0.0;
  result.meanSpeedError = meanAndError(blockSpeeds).second;
  if (result.numberOfEventTerminatedTrajectories > 0uz)
  {
    result.meanTerminationTime = terminationTimeSum / static_cast<double>(result.numberOfEventTerminatedTrajectories);
  }
  if (numberOfTrajectories > 0uz)
  {
    result.meanIntegratedTime = integratedTimeSum / static_cast<double>(numberOfTrajectories);
  }
  if (speed > 0.0)
  {
    result.selfEvent = selfEvent / speed;
    result.collectiveEvent = collectiveEvent / speed;
    result.selfEventError = meanAndError(blockSelfEvent).second;
    result.collectiveEventError = meanAndError(blockCollectiveEvent).second;
  }
  if (speedEventTerminated > 0.0)
  {
    result.selfEventConditioned = 2.0 * selfEventConditioned / speedEventTerminated;
    result.collectiveEventConditioned = 2.0 * collectiveEventConditioned / speedEventTerminated;
    result.selfEventConditionedError = meanAndError(blockSelfEventConditioned).second;
    result.collectiveEventConditionedError = meanAndError(blockCollectiveEventConditioned).second;
  }
  if (numberOfTrajectories > 0uz)
  {
    result.cappedFraction =
        static_cast<double>(result.numberOfCappedTrajectories + result.numberOfUncommittedTrajectories) /
        static_cast<double>(numberOfTrajectories);
  }

  result.self.assign(numberOfTimes, 0.0);
  result.selfTwoState.assign(numberOfTimes, 0.0);
  result.selfCommitted.assign(numberOfTimes, 0.0);
  result.selfError.assign(numberOfTimes, 0.0);
  result.collective.assign(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
  result.collectiveError.assign(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
  result.collectivePerComponent.assign(
      numberOfShells, std::vector<std::vector<double>>(numberOfComponents, std::vector<double>(numberOfTimes, 0.0)));
  result.collectivePlateau.assign(numberOfShells, 0.0);
  result.collectivePlateauError.assign(numberOfShells, 0.0);
  if (speed <= 0.0) return result;

  for (std::size_t t = 0; t < numberOfTimes; ++t)
  {
    result.self[t] = self[t] / speed;
    result.selfTwoState[t] = selfTwoState[t] / speed;
    result.selfCommitted[t] = selfCommitted[t] / speed;
    std::vector<double> values;
    for (const std::vector<double>& k : blockSelf) values.push_back(k[t]);
    result.selfError[t] = meanAndError(values).second;
  }
  // cumulative shells (the per-component parts are the kappa_X contributions of each shell, cumulative)
  std::vector<double> cumulativePerComponent(numberOfComponents * numberOfTimes, 0.0);
  for (std::size_t s = 0; s < numberOfShells; ++s)
  {
    for (std::size_t c = 0; c < numberOfComponents; ++c)
    {
      for (std::size_t t = 0; t < numberOfTimes; ++t)
      {
        cumulativePerComponent[c * numberOfTimes + t] += collective[front.collectiveIndex(s, c, t)];
        result.collectivePerComponent[s][c][t] = cumulativePerComponent[c * numberOfTimes + t] / speed;
      }
    }
    for (std::size_t t = 0; t < numberOfTimes; ++t)
    {
      double sum{0.0};
      for (std::size_t c = 0; c < numberOfComponents; ++c) sum += result.collectivePerComponent[s][c][t];
      result.collective[s][t] = result.self[t] + sum;
      std::vector<double> values;
      for (const std::vector<std::vector<double>>& k : blockCollective) values.push_back(k[s][t]);
      result.collectiveError[s][t] = meanAndError(values).second;
    }
  }

  // plateau values: time average over the window, error from the block values
  auto windowAverage = [&](const std::vector<double>& values)
  {
    double sum{0.0};
    for (std::size_t t : plateauIndices) sum += values[t];
    return sum / static_cast<double>(plateauIndices.size());
  };
  result.selfPlateau = windowAverage(result.self);
  result.selfTwoStatePlateau = windowAverage(result.selfTwoState);
  result.selfCommittedPlateau = windowAverage(result.selfCommitted);
  {
    std::vector<double> values;
    for (const std::vector<double>& k : blockSelf) values.push_back(windowAverage(k));
    result.selfPlateauError = meanAndError(values).second;
    values.clear();
    for (const std::vector<double>& k : blockSelfCommitted) values.push_back(windowAverage(k));
    result.selfCommittedPlateauError = meanAndError(values).second;
  }
  for (std::size_t s = 0; s < numberOfShells; ++s)
  {
    result.collectivePlateau[s] = windowAverage(result.collective[s]);
    std::vector<double> values;
    for (const std::vector<std::vector<double>>& k : blockCollective) values.push_back(windowAverage(k[s]));
    result.collectivePlateauError[s] = meanAndError(values).second;
  }
  result.plateauBegin = plateauBegin;
  result.plateauEnd = plateauEnd;
  // first-order back-correlation factor and collective excess of the event-terminated estimator, relative to
  // the committed kappa_S at the end of the trajectories (its values are held after termination)
  const double selfCommittedFinal = result.selfCommitted.back();
  if (std::abs(selfCommittedFinal) > 1e-12)
  {
    result.backCorrelationFactor = result.selfEvent / selfCommittedFinal;
    result.collectiveExcess = result.collectiveEvent / selfCommittedFinal - 1.0;
    // block errors of the differences; epsilon = (kappa_C^event - kappa_S^event)/kappa_S^committed + (f - 1)
    const double backCorrelationError = meanAndError(blockEventBackCorrelation).second;
    const double excessError = meanAndError(blockEventExcess).second;
    result.backCorrelationFactorError = backCorrelationError / std::abs(selfCommittedFinal);
    result.collectiveExcessError =
        std::sqrt(excessError * excessError + backCorrelationError * backCorrelationError) /
        std::abs(selfCommittedFinal);
    result.backCorrelationFactorConditioned = result.selfEventConditioned / selfCommittedFinal;
    result.collectiveExcessConditioned = result.collectiveEventConditioned / selfCommittedFinal - 1.0;
    // a capped trajectory contributes Delta = +-1 and Q_local = 0, i.e. exactly kappa_S^committed; the
    // deviations f - 1 and epsilon therefore scale with the event-terminated fraction 1 - p_cap
    const double eventFraction = 1.0 - result.cappedFraction;
    if (eventFraction > 1e-12)
    {
      result.backCorrelationFactorCapCorrected = 1.0 + (result.backCorrelationFactor - 1.0) / eventFraction;
      result.collectiveExcessCapCorrected = result.collectiveExcess / eventFraction;
    }
    // decomposition of epsilon - (f - 1) = (kappa_C^event - kappa_S^event) / kappa_S^committed into the
    // crossings before and after commitment, and per channel
    if (speed > 0.0)
    {
      result.collectiveExcessBeforeCommitment = collectiveEventBefore / speed / selfCommittedFinal;
      result.collectiveExcessAfterCommitment = collectiveEventAfter / speed / selfCommittedFinal;
      result.collectiveExcessBeforeCommitmentError =
          meanAndError(blockCollectiveEventBefore).second / std::abs(selfCommittedFinal);
      result.collectiveExcessAfterCommitmentError =
          meanAndError(blockCollectiveEventAfter).second / std::abs(selfCommittedFinal);
      for (std::size_t phase = 0; phase < 2; ++phase)
      {
        for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
        {
          result.collectiveExcessChannel[phase][c] = localQChannel[phase][c] / speed / selfCommittedFinal;
          result.channelMeanTime[phase][c] =
              (result.channelEvents[phase][c] > 0uz)
                  ? channelTimeSum[phase][c] / static_cast<double>(result.channelEvents[phase][c])
                  : 0.0;
        }
      }
    }
  }
  return result;
}

// =================================================================================================
// Output
// =================================================================================================

void DynamicallyCorrectedTST::output()
{
  const System& system = systems.front();
  const Component& tagged = system.components[settings.taggedComponent];
  const std::size_t numberOfShells = settings.maximumCollectiveShell + 2uz;
  const std::size_t numberOfComponents = system.components.size();

  const TransmissionResult kappa = analyzeTransmission();
  const bool occupancyMode = settings.cageOccupancies.has_value();
  const bool eventMode = settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event;

  // --- transmission-function file ---------------------------------------------------------------
  if (kappa.numberOfTrajectoryPairs > 0uz)
  {
    std::ofstream file(std::format("dctst/transmission_function_{}_{}.txt", system.temperature, system.input_pressure));
    std::print(file, "# dcTST transmission functions of component '{}' ({} trajectory pairs)\n", tagged.name,
               kappa.numberOfTrajectoryPairs);
    if (eventMode)
    {
      std::print(file,
                 "# event-terminated trajectories: each trajectory is stopped after {} local axial event(s) following "
                 "commitment (cap {} ps);\n# its final values are held for the later time points\n",
                 settings.terminationEvents, settings.trajectoryTime);
    }
    std::print(file, "# column 1: t [ps]\n# column 2: kappa_S(t) (net cage count)\n# column 3: error\n");
    std::print(file, "# column 4: kappa_S(t) two-state (side of the plane, Bennett-Chandler)\n");
    std::print(file, "# column 5: kappa_S(t) committed (side of the plane, frozen once |q| >= {:.4f} A)\n",
               (settings.commitDistance > 0.0) ? settings.commitDistance : 0.5 * barrier.lambda);
    std::size_t column = 6;
    for (std::size_t s = 0; s < numberOfShells; ++s)
    {
      if (s + 1 == numberOfShells)
        std::print(file, "# column {}: kappa_C(t) all molecules\n# column {}: error\n", column, column + 1);
      else
        std::print(file, "# column {}: kappa_C(t) molecules within {} cage(s) of A/B\n# column {}: error\n", column, s,
                   column + 1);
      column += 2;
    }
    for (std::size_t s = 0; s < numberOfShells; ++s)
    {
      for (std::size_t c = 0; c < numberOfComponents; ++c)
      {
        std::print(file, "# column {}: kappa_X(t) contribution of component '{}', shell {}\n", column,
                   system.components[c].name, (s + 1 == numberOfShells) ? std::string("all") : std::to_string(s));
        ++column;
      }
    }
    for (std::size_t t = 0; t < numberOfTimes; ++t)
    {
      std::print(file, "{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}", kappa.time[t], kappa.self[t], kappa.selfError[t],
                 kappa.selfTwoState[t], kappa.selfCommitted[t]);
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        std::print(file, " {:12.6f} {:12.6f}", kappa.collective[s][t], kappa.collectiveError[s][t]);
      }
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        for (std::size_t c = 0; c < numberOfComponents; ++c)
        {
          std::print(file, " {:12.6f}", kappa.collectivePerComponent[s][c][t]);
        }
      }
      std::print(file, "\n");
    }
  }

  // --- summary ----------------------------------------------------------------------------------
  std::print(stream, "\n\nResults\n");
  std::print(stream, "===============================================================================\n\n");

  const double kBT_over_m = Units::KB * system.temperature / tagged.totalMass;
  const double thermalSpeed = std::sqrt(2.0 * kBT_over_m / std::numbers::pi);  // <|v_n|> of a free molecule

  if (barrierProbability.has_value())
  {
    std::print(stream, "Barrier probability density P(q*) ({}):\n", barrierProbabilitySource);
    std::print(stream, "    P(q*) = exp(-beta F*)/int exp(-beta F): {:.6e} +/- {:.3e} [1/A]\n", barrierProbability.value(),
               barrierProbabilityError);
  }
  else
  {
    std::print(stream,
               "Barrier probability density P(q*): not given (set 'FreeEnergyInput' to the JSON output of the "
               "free-energy driver or 'BarrierProbability'); only the transmission coefficients are reported\n");
  }
  std::print(stream, "\n");

  if (occupancyMode)
  {
    std::print(stream, "Occupancy-resolved mode: n_A = {} (tagged molecule included), n_B = {}\n",
               settings.cageOccupancies->first, settings.cageOccupancies->second);
    std::size_t cageChanges{0};
    std::size_t rejectedMoves{0};
    std::size_t rejectedTaggedMoves{0};
    for (const DcTSTOccupancyState& state : occupancyStates)
    {
      cageChanges += state.cageChanges;
      rejectedMoves += state.rejectedMoves;
      rejectedTaggedMoves += state.rejectedTaggedMoves;
    }
    std::print(stream, "    accepted moves of other molecules changing their cage: {}\n", cageChanges);
    std::print(stream, "    of which rejected by the constraint (A or B involved):  {}\n", rejectedMoves);
    std::print(stream, "    in-plane moves of the tagged molecule rejected:       {}\n\n", rejectedTaggedMoves);
  }

  if (kappa.numberOfTrajectoryPairs > 0uz)
  {
    std::print(stream, "Transmission functions ({} trajectory pairs, {} replicas x {} blocks):\n",
               kappa.numberOfTrajectoryPairs, numberOfReplicas, numberOfBlocks);
    if (stoppedEarly.load(std::memory_order_relaxed))
    {
      const auto [minimum, maximum] =
          std::minmax_element(completedTransmissionCycles.begin(), completedTransmissionCycles.end());
      std::print(stream,
                 "    STOPPED EARLY by a termination signal:  {}..{} of {} production cycles completed per replica\n",
                 *minimum, *maximum, settings.numberOfTransmissionCycles);
    }
    std::print(stream, "    <|qdot|> at the barrier:             {:.6f} +/- {:.6f} [A/ps]\n", kappa.meanSpeed,
               kappa.meanSpeedError);
    std::print(stream, "    sqrt(2 kT/(pi m)) of the molecule:   {:.6f} [A/ps]\n", thermalSpeed);
    if (settings.trajectoryThermostat == DcTSTSettings::TrajectoryThermostat::None)
    {
      std::print(stream, "    bath:                                none (NVE)\n");
      std::print(stream, "    mean |relative energy drift|:        {:.3e}\n", kappa.meanRelativeEnergyDrift);
      std::print(stream, "    maximum |relative energy drift|:     {:.3e}\n", kappa.maximumRelativeEnergyDrift);
    }
    else
    {
      std::print(
          stream, "    bath:                                {} with rate {} [1/ps]\n",
          settings.trajectoryThermostat == DcTSTSettings::TrajectoryThermostat::Andersen ? "Andersen" : "Langevin",
          settings.bathRate);
    }
    std::print(stream, "    long-jump fraction:                  {:.4f} (tagged molecule left cages A/B)\n",
               kappa.longJumpFraction);
    std::print(stream, "    committed fraction:                  {:.4f} (tagged molecule reached |q| >= {:.3f} A)\n",
               kappa.committedFraction,
               (settings.commitDistance > 0.0) ? settings.commitDistance : 0.5 * barrier.lambda);
    std::print(stream, "    plateau window:                      [{}, {}] [ps]\n", kappa.plateauBegin,
               kappa.plateauEnd);
    std::print(stream, "    kappa_S (plateau):                   {:.5f} +/- {:.5f}\n", kappa.selfPlateau,
               kappa.selfPlateauError);
    std::print(stream, "    kappa_S two-state (plateau):         {:.5f}\n", kappa.selfTwoStatePlateau);
    std::print(stream, "    kappa_S committed (plateau):         {:.5f} +/- {:.5f}\n", kappa.selfCommittedPlateau,
               kappa.selfCommittedPlateauError);
    for (std::size_t s = 0; s < numberOfShells; ++s)
    {
      if (s + 1 == numberOfShells)
        std::print(stream, "    kappa_C (plateau, all molecules):    {:.5f} +/- {:.5f}\n", kappa.collectivePlateau[s],
                   kappa.collectivePlateauError[s]);
      else
        std::print(stream, "    kappa_C (plateau, shell <= {}):       {:.5f} +/- {:.5f}\n", s,
                   kappa.collectivePlateau[s], kappa.collectivePlateauError[s]);
    }
    std::print(stream, "\n");

    if (eventMode)
    {
      using EventChannel = DcTSTTransmissionAccumulator::EventChannel;
      const std::size_t total = 2uz * kappa.numberOfTrajectoryPairs;
      auto fraction = [&](std::size_t n) { return static_cast<double>(n) / static_cast<double>(std::max(1uz, total)); };
      std::print(stream, "Event-terminated estimator ({} local axial event(s) after commitment, cap {} [ps]):\n",
                 settings.terminationEvents, settings.trajectoryTime);
      std::print(stream, "    trajectories ended by an event:      {} ({:.4f})\n",
                 kappa.numberOfEventTerminatedTrajectories, fraction(kappa.numberOfEventTerminatedTrajectories));
      std::print(stream, "    trajectories capped after commitment: {} ({:.4f})\n", kappa.numberOfCappedTrajectories,
                 fraction(kappa.numberOfCappedTrajectories));
      std::print(stream, "    trajectories never committed:        {} ({:.4f})\n",
                 kappa.numberOfUncommittedTrajectories, fraction(kappa.numberOfUncommittedTrajectories));
      std::print(stream, "    termination time (event-ended):      mean {:.4f}, maximum {:.4f} [ps]\n",
                 kappa.meanTerminationTime, kappa.maximumTerminationTime);
      std::print(stream, "    integrated time per trajectory:      mean {:.4f} [ps]\n", kappa.meanIntegratedTime);
      std::print(stream, "    kappa_S^event (net cage count):      {:.5f} +/- {:.5f}\n", kappa.selfEvent,
                 kappa.selfEventError);
      std::print(stream, "    kappa_S^committed (final):           {:.5f}\n", kappa.selfCommitted.back());
      std::print(stream, "    kappa_C^event (self + Q_local):      {:.5f} +/- {:.5f}\n", kappa.collectiveEvent,
                 kappa.collectiveEventError);
      std::print(stream, "    f = kappa_S^event / kappa_S^committed:        {:.5f} +/- {:.5f}\n",
                 kappa.backCorrelationFactor, kappa.backCorrelationFactorError);
      std::print(stream, "    epsilon = kappa_C^event / kappa_S^committed - 1: {:.5f} +/- {:.5f}\n",
                 kappa.collectiveExcess, kappa.collectiveExcessError);
      std::print(stream, "    (errors from the block values of the differences kappa_C^event - kappa_S^event and\n"
                         "     kappa_S^event - kappa_S^committed, which share the trajectories and are much less noisy\n"
                         "     than the kappas themselves; the error of k^TST is not included)\n");
      std::print(stream, "  Cap-corrected (deviations divided by the event-ended fraction 1 - p_cap, p_cap = {:.4f}):\n",
                 kappa.cappedFraction);
      std::print(stream, "    f_corr = 1 + (f - 1) / (1 - p_cap):              {:.5f}\n",
                 kappa.backCorrelationFactorCapCorrected);
      std::print(stream, "    epsilon_corr = epsilon / (1 - p_cap):            {:.5f}\n",
                 kappa.collectiveExcessCapCorrected);
      std::print(stream, "  Conditioned on the event-ended trajectories only (per trajectory, not per pair):\n");
      std::print(stream, "    kappa_S^event|cond:                  {:.5f} +/- {:.5f}\n", kappa.selfEventConditioned,
                 kappa.selfEventConditionedError);
      std::print(stream, "    kappa_C^event|cond:                  {:.5f} +/- {:.5f}\n",
                 kappa.collectiveEventConditioned, kappa.collectiveEventConditionedError);
      std::print(stream, "    f_cond = kappa_S^event|cond / kappa_S^committed:        {:.5f}\n",
                 kappa.backCorrelationFactorConditioned);
      std::print(stream, "    epsilon_cond = kappa_C^event|cond / kappa_S^committed - 1: {:.5f}\n",
                 kappa.collectiveExcessConditioned);
      std::print(stream,
                 "    (f_corr/epsilon_corr and f_cond/epsilon_cond agree when the capped trajectories are an "
                 "unbiased sample;\n"
                 "    a cap T >~ 5 tau_event keeps p_cap below a few percent and makes all three estimates "
                 "coincide)\n");
      std::print(stream,
                 "    channels in the frame of the tagged hop (A = cage of departure, B = arrival, A' behind A, B' "
                 "beyond B;\n"
                 "    the time-reversed trajectory is mirrored so that its own hop direction is 'forward'):\n");
      std::print(stream, "    {:<26} {:>10} {:>10} {:>10}   {:>10} {:>10} {:>10}\n", "channel", "term.fwd", "term.rev",
                 "term.all", "all.fwd", "all.rev", "all.all");
      static constexpr std::array<const char*, static_cast<std::size_t>(EventChannel::Count)> channelNames = {
          "tagged back-hop B->A", "tagged forward hop B->B'",
          "other A->B",           "other B->A",
          "other A->A'",          "other A'->A",
          "other B->B'",          "other B'->B"};
      std::size_t terminatingSum{0};
      for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
      {
        std::print(stream, "    {:<26} {:>10} {:>10} {:>10}   {:>10} {:>10} {:>10}\n", channelNames[c],
                   kappa.terminatingEvents[0][c], kappa.terminatingEvents[1][c],
                   kappa.terminatingEvents[0][c] + kappa.terminatingEvents[1][c], kappa.allEvents[0][c],
                   kappa.allEvents[1][c], kappa.allEvents[0][c] + kappa.allEvents[1][c]);
        terminatingSum += kappa.terminatingEvents[0][c] + kappa.terminatingEvents[1][c];
      }
      std::print(stream, "    {:<26} {:>32}\n", "terminating events total", terminatingSum);
      std::print(stream,
                 "    (Q_local counts the axial crossings of the other molecules in the column from t = 0; events\n"
                 "     before the commitment of the tagged molecule do not count towards termination)\n\n");

      // decomposition of epsilon: tagged part (f - 1) and the Q_local crossings before/after commitment, per
      // channel
      std::print(stream, "  Decomposition of the collective excess, epsilon = (f - 1) + epsilon_before + epsilon_after:\n");
      std::print(stream, "    (f - 1), tagged molecule after commitment:         {:>9.5f} +/- {:.5f}\n",
                 kappa.backCorrelationFactor - 1.0, kappa.backCorrelationFactorError);
      std::print(stream, "    epsilon_before, others before tagged commitment:   {:>9.5f} +/- {:.5f}\n",
                 kappa.collectiveExcessBeforeCommitment, kappa.collectiveExcessBeforeCommitmentError);
      std::print(stream, "    epsilon_after, others after tagged commitment:     {:>9.5f} +/- {:.5f}\n",
                 kappa.collectiveExcessAfterCommitment, kappa.collectiveExcessAfterCommitmentError);
      std::print(stream, "    sum:                                               {:>9.5f}   (epsilon = {:.5f})\n",
                 (kappa.backCorrelationFactor - 1.0) + kappa.collectiveExcessBeforeCommitment +
                     kappa.collectiveExcessAfterCommitment,
                 kappa.collectiveExcess);
      std::print(stream,
                 "    per channel (frame of the tagged hop; 'events' are all crossings of both members of the pair,\n"
                 "    <t> the mean time since launch (before) or since commitment (after), 'epsilon' the |qdot|-weighted\n"
                 "    signed contribution to epsilon; the tagged channels act through f and carry no epsilon):\n");
      std::print(stream, "    {:<26} {:>8} {:>9} {:>10}   {:>8} {:>9} {:>10}\n", "channel", "events", "<t> [ps]",
                 "epsilon", "events", "<t> [ps]", "epsilon");
      std::print(stream, "    {:<26} {:^29}   {:^29}\n", "", "before commitment", "after commitment");
      for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
      {
        std::print(stream, "    {:<26} {:>8} {:>9.3f} {:>10.5f}   {:>8} {:>9.3f} {:>10.5f}\n", channelNames[c],
                   kappa.channelEvents[0][c], kappa.channelMeanTime[0][c], kappa.collectiveExcessChannel[0][c],
                   kappa.channelEvents[1][c], kappa.channelMeanTime[1][c], kappa.collectiveExcessChannel[1][c]);
      }
      std::print(stream, "\n");
    }

    // rates and diffusion coefficients: k^TST = <|qdot|>/2 P(q*); D = kappa k^TST lambda^2
    auto report = [&](std::string_view label, double importedProbability, double importedProbabilityError)
    {
      const double kTST = 0.5 * kappa.meanSpeed * importedProbability;
      const double kTSTError =
          kTST * std::sqrt(std::pow(kappa.meanSpeedError / std::max(kappa.meanSpeed, 1e-300), 2) +
                           std::pow(importedProbabilityError / std::max(importedProbability, 1e-300), 2));
      const double lambdaSquared = barrier.lambda * barrier.lambda;
      const double conversion = Units::DiffusionConversionFactor;
      std::print(stream, "Rates and diffusion coefficients ({}):\n", label);
      std::print(stream, "    k^TST:                               {:.6e} +/- {:.2e} [1/ps]\n", kTST, kTSTError);
      const double kS = kappa.selfPlateau * kTST;
      const double kSError = kS * std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                                            std::pow(kappa.selfPlateauError / std::max(kappa.selfPlateau, 1e-300), 2));
      std::print(stream, "    k_(A->B) = kappa_S k^TST:            {:.6e} +/- {:.2e} [1/ps]\n", kS, kSError);
      std::print(stream, "    D_S = kappa_S k^TST lambda^2:        {:.6e} +/- {:.2e} [m^2/s]\n",
                 kS * lambdaSquared * conversion, kSError * lambdaSquared * conversion);
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        const double kC = kappa.collectivePlateau[s] * kTST;
        const double kCError =
            kC *
            std::sqrt(
                std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                std::pow(kappa.collectivePlateauError[s] / std::max(std::abs(kappa.collectivePlateau[s]), 1e-300), 2));
        if (s + 1 == numberOfShells)
          std::print(stream, "    D_C = kappa_C k^TST lambda^2 (all):  {:.6e} +/- {:.2e} [m^2/s]\n",
                     kC * lambdaSquared * conversion, kCError * lambdaSquared * conversion);
        else
          std::print(stream, "    D_C (shell <= {}):                    {:.6e} +/- {:.2e} [m^2/s]\n", s,
                     kC * lambdaSquared * conversion, kCError * lambdaSquared * conversion);
      }
      double kSEvent{0.0}, kSEventError{0.0}, kCEvent{0.0}, kCEventError{0.0};
      if (eventMode)
      {
        kSEvent = kappa.selfEvent * kTST;
        kSEventError = std::abs(kSEvent) *
                       std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                                 std::pow(kappa.selfEventError / std::max(std::abs(kappa.selfEvent), 1e-300), 2));
        kCEvent = kappa.collectiveEvent * kTST;
        kCEventError =
            std::abs(kCEvent) *
            std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                      std::pow(kappa.collectiveEventError / std::max(std::abs(kappa.collectiveEvent), 1e-300), 2));
        std::print(stream, "    event-terminated estimator:\n");
        std::print(stream, "    k_(A->B) = kappa_S^event k^TST:      {:.6e} +/- {:.2e} [1/ps]\n", kSEvent,
                   kSEventError);
        std::print(stream, "    D_S = kappa_S^event k^TST lambda^2:  {:.6e} +/- {:.2e} [m^2/s]\n",
                   kSEvent * lambdaSquared * conversion, kSEventError * lambdaSquared * conversion);
        std::print(stream, "    D_C = kappa_C^event k^TST lambda^2:  {:.6e} +/- {:.2e} [m^2/s]\n",
                   kCEvent * lambdaSquared * conversion, kCEventError * lambdaSquared * conversion);
      }
      if (occupancyMode)
      {
        std::print(stream,
                   "    (occupancy-resolved: k^TST and kappa are conditioned on n_A = {}, n_B = {}; the D values are\n"
                   "     the hypothetical diffusion coefficients of a lattice in which every window had these\n"
                   "     occupancies; recombine the rates k(n_A -> n_B) with P(n_A, n_B) over the sampled pairs)\n\n",
                   settings.cageOccupancies->first, settings.cageOccupancies->second);
      }
      else
      {
        std::print(stream,
                   "    (D assumes a lattice of equivalent cages with one window per lattice direction;\n"
                   "     k_(A->B) is the rate to one specific neighboring cage)\n\n");
      }

      nlohmann::json& json = outputJson["results"][std::string(label)];
      json["kTST"] = kTST;
      json["kTSTError"] = kTSTError;
      json["kSelf"] = kS;
      json["DSelf"] = kS * lambdaSquared * conversion;
      json["DSelfError"] = kSError * lambdaSquared * conversion;
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        const std::string key = (s + 1 == numberOfShells) ? "all" : std::to_string(s);
        json["DCollective"][key] = kappa.collectivePlateau[s] * kTST * lambdaSquared * conversion;
        json["DCollectiveError"][key] = kappa.collectivePlateauError[s] * kTST * lambdaSquared * conversion;
      }
      if (eventMode)
      {
        json["kSelfEvent"] = kSEvent;
        json["kSelfEventError"] = kSEventError;
        json["DSelfEvent"] = kSEvent * lambdaSquared * conversion;
        json["DSelfEventError"] = kSEventError * lambdaSquared * conversion;
        json["DCollectiveEvent"] = kCEvent * lambdaSquared * conversion;
        json["DCollectiveEventError"] = kCEventError * lambdaSquared * conversion;
      }
    };
    if (barrierProbability.has_value())
    {
      report(barrierProbabilitySource, barrierProbability.value(), barrierProbabilityError);
    }
    else
    {
      std::print(stream, "Rates and diffusion coefficients: not available without P(q*)\n\n");
    }
  }
  else
  {
    std::print(stream, "Transmission functions: no trajectories were launched\n\n");
  }

  std::print(stream, "Timings\n");
  std::print(stream, "    initialization:                      {:.3f} [s]\n", totalInitializationTime.count());
  std::print(stream, "    transmission:                        {:.3f} [s]\n", totalTransmissionTime.count());
  std::print(stream, "    total:                               {:.3f} [s]\n", totalSimulationTime.count());
  std::flush(stream);

  outputJson["results"]["numberOfTrajectoryPairs"] = kappa.numberOfTrajectoryPairs;
  outputJson["results"]["stoppedEarly"] = stoppedEarly.load(std::memory_order_relaxed);
  outputJson["results"]["transmissionCyclesCompleted"] = completedTransmissionCycles;
  outputJson["results"]["transmissionCyclesRequested"] = settings.numberOfTransmissionCycles;
  outputJson["results"]["meanSpeed"] = kappa.meanSpeed;
  outputJson["results"]["kappaSelf"] = kappa.selfPlateau;
  outputJson["results"]["kappaSelfError"] = kappa.selfPlateauError;
  outputJson["results"]["kappaSelfTwoState"] = kappa.selfTwoStatePlateau;
  outputJson["results"]["kappaSelfCommitted"] = kappa.selfCommittedPlateau;
  outputJson["results"]["kappaSelfCommittedError"] = kappa.selfCommittedPlateauError;
  outputJson["results"]["longJumpFraction"] = kappa.longJumpFraction;
  outputJson["results"]["committedFraction"] = kappa.committedFraction;
  for (std::size_t s = 0; s < numberOfShells; ++s)
  {
    const std::string key = (s + 1 == numberOfShells) ? "all" : std::to_string(s);
    outputJson["results"]["kappaCollective"][key] = kappa.collectivePlateau[s];
    outputJson["results"]["kappaCollectiveError"][key] = kappa.collectivePlateauError[s];
  }
  if (eventMode)
  {
    using EventChannel = DcTSTTransmissionAccumulator::EventChannel;
    nlohmann::json& json = outputJson["results"]["eventTerminated"];
    json["kappaSelfEvent"] = kappa.selfEvent;
    json["kappaSelfEventError"] = kappa.selfEventError;
    json["kappaSelfCommittedFinal"] = kappa.selfCommitted.back();
    json["kappaCollectiveEvent"] = kappa.collectiveEvent;
    json["kappaCollectiveEventError"] = kappa.collectiveEventError;
    json["backCorrelationFactor"] = kappa.backCorrelationFactor;
    json["backCorrelationFactorError"] = kappa.backCorrelationFactorError;
    json["collectiveExcess"] = kappa.collectiveExcess;
    json["collectiveExcessError"] = kappa.collectiveExcessError;
    json["cappedFraction"] = kappa.cappedFraction;
    json["backCorrelationFactorCapCorrected"] = kappa.backCorrelationFactorCapCorrected;
    json["collectiveExcessCapCorrected"] = kappa.collectiveExcessCapCorrected;
    json["kappaSelfEventConditioned"] = kappa.selfEventConditioned;
    json["kappaSelfEventConditionedError"] = kappa.selfEventConditionedError;
    json["kappaCollectiveEventConditioned"] = kappa.collectiveEventConditioned;
    json["kappaCollectiveEventConditionedError"] = kappa.collectiveEventConditionedError;
    json["backCorrelationFactorConditioned"] = kappa.backCorrelationFactorConditioned;
    json["collectiveExcessConditioned"] = kappa.collectiveExcessConditioned;
    json["numberOfEventTerminatedTrajectories"] = kappa.numberOfEventTerminatedTrajectories;
    json["numberOfCappedTrajectories"] = kappa.numberOfCappedTrajectories;
    json["numberOfUncommittedTrajectories"] = kappa.numberOfUncommittedTrajectories;
    json["meanTerminationTime"] = kappa.meanTerminationTime;
    json["maximumTerminationTime"] = kappa.maximumTerminationTime;
    json["meanIntegratedTime"] = kappa.meanIntegratedTime;
    static constexpr std::array<const char*, static_cast<std::size_t>(EventChannel::Count)> channelKeys = {
        "taggedBackHop",  "taggedForwardHop", "otherAtoB", "otherBtoA",
        "otherAtoAPrime", "otherAPrimeToA",   "otherBtoC", "otherCtoB"};
    for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
    {
      json["terminatingEvents"]["forward"][channelKeys[c]] = kappa.terminatingEvents[0][c];
      json["terminatingEvents"]["reverse"][channelKeys[c]] = kappa.terminatingEvents[1][c];
      json["allEvents"]["forward"][channelKeys[c]] = kappa.allEvents[0][c];
      json["allEvents"]["reverse"][channelKeys[c]] = kappa.allEvents[1][c];
      for (std::size_t phase = 0; phase < 2; ++phase)
      {
        nlohmann::json& entry =
            json["decomposition"]["channels"][channelKeys[c]][(phase == 0uz) ? "beforeCommitment" : "afterCommitment"];
        entry["events"] = kappa.channelEvents[phase][c];
        entry["meanTime"] = kappa.channelMeanTime[phase][c];
        entry["collectiveExcess"] = kappa.collectiveExcessChannel[phase][c];
      }
    }
    json["decomposition"]["backCorrelation"] = kappa.backCorrelationFactor - 1.0;
    json["decomposition"]["backCorrelationError"] = kappa.backCorrelationFactorError;
    json["decomposition"]["collectiveExcessBeforeCommitment"] = kappa.collectiveExcessBeforeCommitment;
    json["decomposition"]["collectiveExcessBeforeCommitmentError"] = kappa.collectiveExcessBeforeCommitmentError;
    json["decomposition"]["collectiveExcessAfterCommitment"] = kappa.collectiveExcessAfterCommitment;
    json["decomposition"]["collectiveExcessAfterCommitmentError"] = kappa.collectiveExcessAfterCommitmentError;
  }
  outputJson["results"]["meanRelativeEnergyDrift"] = kappa.meanRelativeEnergyDrift;
  if (barrierProbability.has_value())
  {
    outputJson["results"]["barrierProbability"] = barrierProbability.value();
    outputJson["results"]["barrierProbabilityError"] = barrierProbabilityError;
    outputJson["results"]["barrierProbabilitySource"] = barrierProbabilitySource;
  }
  outputJson["timings"]["total"] = totalSimulationTime.count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}
