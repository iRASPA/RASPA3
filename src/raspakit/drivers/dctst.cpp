module;

module dctst;

import std;

import mc_moves_acceptance;

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
import move_statistics;
import mc_moves_statistics;
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
import interactions_polarization;
import thermostat;
import json;
import graceful_shutdown;
import dctst_settings;
import dctst_barrier;
import dctst_moves;

// =================================================================================================
// Accumulators
// =================================================================================================

DcTSTTransmissionAccumulator::DcTSTTransmissionAccumulator(std::size_t numberOfBlocks, std::size_t numberOfShells,
                                                           std::size_t numberOfComponents, std::size_t numberOfTimes,
                                                           std::size_t numberOfHandovers,
                                                           std::size_t numberOfEventOrders)
    : numberOfBlocks(numberOfBlocks),
      numberOfShells(numberOfShells),
      numberOfComponents(numberOfComponents),
      numberOfTimes(numberOfTimes),
      speed(numberOfBlocks, 0.0),
      self(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      selfTwoState(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      selfCommitted(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      selfSmooth(numberOfBlocks, std::vector<double>(numberOfTimes, 0.0)),
      collective(numberOfBlocks, std::vector<double>(numberOfShells * numberOfComponents * numberOfTimes, 0.0)),
      numberOfTrajectoryPairs(numberOfBlocks, 0uz),
      byEvent(std::max(1uz, numberOfEventOrders)),
      collectiveEventBeforeCommitment(numberOfBlocks, 0.0),
      hybridTagged(numberOfHandovers, std::vector<double>(numberOfBlocks, 0.0)),
      hybridOthers(numberOfHandovers, std::vector<double>(numberOfBlocks, 0.0)),
      hybridOccupancy(numberOfHandovers),
      hybridNeighbourhood(numberOfHandovers),
      hybridOccupancySigned(numberOfHandovers),
      hybridNeighbourhoodSigned(numberOfHandovers),
      hybridEventBefore(numberOfHandovers, 0uz),
      hybridTruncated(numberOfHandovers, 0uz)
{
  for (EventOrderSums& sums : byEvent)
  {
    sums.self.assign(numberOfBlocks, 0.0);
    sums.collective.assign(numberOfBlocks, 0.0);
    sums.speedTerminated.assign(numberOfBlocks, 0.0);
    sums.selfConditioned.assign(numberOfBlocks, 0.0);
    sums.collectiveConditioned.assign(numberOfBlocks, 0.0);
    sums.collectiveAfterCommitment.assign(numberOfBlocks, 0.0);
  }
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
                     std::span<const GroupState> groupData, std::size_t groupIndex, std::size_t bead)
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
  // a bead inside a rigid fragment of a semi-flexible molecule carries no velocity of its own: it moves with the
  // fragment's rigid-body state (the states follow the rigid fragments in ascending fragment order)
  if (const std::optional<std::size_t> fragmentIndex = component.rigidFragmentContaining(bead))
  {
    const Fragment& fragment = component.fragmentGraph.fragments[*fragmentIndex];
    std::size_t rigidRank{};
    for (std::size_t f = 0; f < *fragmentIndex; ++f)
    {
      if (component.fragmentGraph.fragments[f].isRigidBody()) ++rigidRank;
    }
    const GroupState& state = groupData[groupIndex + rigidRank];
    const std::size_t k = static_cast<std::size_t>(std::ranges::find(fragment.atoms, bead) - fragment.atoms.begin());

    simd_quatd p = state.orientationMomentum;
    p.r = -p.r;
    const simd_quatd pq = p * state.orientation;
    const double3 omegaBody = 0.5 * double3(pq.ix, pq.iy, pq.iz) * fragment.inverseInertiaVector;
    return state.velocity + state.orientation * double3::cross(omegaBody, fragment.bodyFixedPositions[k]);
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
    molecule.orientationMomentum = simd_quatd(-molecule.orientationMomentum.ix, -molecule.orientationMomentum.iy,
                                              -molecule.orientationMomentum.iz, -molecule.orientationMomentum.r);
  }
  for (AtomDynamics& dynamic : dynamics)
  {
    dynamic.velocity = -dynamic.velocity;
  }
  for (GroupState& state : groupData)
  {
    state.velocity = -state.velocity;
    state.orientationMomentum = simd_quatd(-state.orientationMomentum.ix, -state.orientationMomentum.iy,
                                           -state.orientationMomentum.iz, -state.orientationMomentum.r);
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
               double temperature, std::span<Molecule> moleculeData, std::span<const Atom> atoms,
               std::span<AtomDynamics> dynamics, std::span<GroupState> groupData,
               const std::vector<Component>& components)
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
    std::span<const Atom> moleculeAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
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
            moleculeDynamics[i].velocity = gaussian3() * std::sqrt(kT / component.atomMass(i, moleculeAtoms[i]));
          }
          Integrators::constrainMoleculeVelocities(component, moleculeAtoms, moleculeDynamics);
          Integrators::syncNonRigidMoleculeCenterOfMass(molecule, moleculeAtoms, moleculeDynamics, component, groupData,
                                                        groupIndex);
        }
        else
        {
          Integrators::initializeMoleculeVelocity(rng, molecule, moleculeAtoms, moleculeDynamics, component,
                                                  temperature);
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
          moleculeDynamics[i].velocity = c1 * moleculeDynamics[i].velocity +
                                         c2 * std::sqrt(kT / component.atomMass(i, moleculeAtoms[i])) * gaussian3();
        }
        // the random kick is drawn for free atoms: keep the part the constraints allow
        Integrators::constrainMoleculeVelocities(component, moleculeAtoms, moleculeDynamics);
        Integrators::syncNonRigidMoleculeCenterOfMass(molecule, moleculeAtoms, moleculeDynamics, component, groupData,
                                                      groupIndex);
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
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          const double mass = component.atomMass(i, moleculeAtoms[i]);
          moleculeDynamics[i].velocity = c1 * moleculeDynamics[i].velocity + c2 * std::sqrt(kT / mass) * gaussian3();
        }
        // the random kick is drawn for free atoms: keep the part the constraints allow
        Integrators::constrainMoleculeVelocities(component, moleculeAtoms, moleculeDynamics);
        double3 momentum{};
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          momentum += component.atomMass(i, moleculeAtoms[i]) * moleculeDynamics[i].velocity;
        }
        molecule.velocity = momentum * molecule.invMass;
      }
    }
    groupIndex += numberOfRigidFragments;
  }
}

/**
 * \brief Takes the momenta and the forces from the molecules of the components that are held in place.
 *
 * A component without a Monte Carlo move (extra-framework cations on their sites) is not sampled by the Monte
 * Carlo stages: its molecules are part of the fixed environment of the constrained ensemble. Called after the
 * forces have been computed, this leaves such a molecule nothing to be kicked with in the next step of velocity
 * Verlet: its position does not change, and it carries no kinetic energy.
 */
void holdInPlace(std::span<const char> heldComponent, std::span<Molecule> moleculeData, std::span<AtomDynamics> dynamics,
                 std::span<GroupState> groupData, const std::vector<Component>& components)
{
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    const std::size_t numberOfRigidFragments = component.isSemiFlexible() ? component.numberOfRigidFragments() : 0uz;
    if (heldComponent[molecule.componentId])
    {
      molecule.velocity = double3(0.0, 0.0, 0.0);
      molecule.gradient = double3(0.0, 0.0, 0.0);
      molecule.orientationMomentum = simd_quatd(0.0, 0.0, 0.0, 0.0);
      molecule.orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);
      for (AtomDynamics& dynamic : dynamics.subspan(molecule.atomIndex, molecule.numberOfAtoms))
      {
        dynamic.velocity = double3(0.0, 0.0, 0.0);
        dynamic.gradient = double3(0.0, 0.0, 0.0);
      }
      for (std::size_t g = groupIndex; g < groupIndex + numberOfRigidFragments && g < groupData.size(); ++g)
      {
        groupData[g].velocity = double3(0.0, 0.0, 0.0);
        groupData[g].gradient = double3(0.0, 0.0, 0.0);
        groupData[g].orientationMomentum = simd_quatd(0.0, 0.0, 0.0, 0.0);
        groupData[g].orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);
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

namespace
{
// Runs 'body' on a replica thread and keeps the exception it throws (an exception that leaves a std::jthread calls
// std::terminate); the caller rethrows it after the threads have joined.
template <typename F>
auto keepingException(std::exception_ptr& error, F body)
{
  return [&error, body = std::move(body)]() mutable
  {
    try
    {
      body();
    }
    catch (...)
    {
      error = std::current_exception();
    }
  };
}

void rethrowWorkerException(const std::vector<std::exception_ptr>& errors)
{
  for (const std::exception_ptr& error : errors)
  {
    if (error) std::rethrow_exception(error);
  }
}

// The physics the drivers support, checked here so that programmatic construction is held to the same rules as the
// input reader. Both drivers: a rigid framework (the barrier and the cages are fixed in the laboratory frame). The
// transmission driver integrates trajectories with the molecular-dynamics integrator, which evaluates neither
// polarization nor external-field forces (the Monte Carlo stages would include them) and does not give the
// reaction-bead velocity of a semi-flexible tagged molecule.
void validateDcTSTSystem(const System& system, std::size_t taggedComponent, std::string_view driverName,
                         bool integratesTrajectories)
{
  if (taggedComponent >= system.components.size())
  {
    throw std::runtime_error(std::format("dcTST: '{}': the tagged component {} does not exist\n", driverName,
                                         taggedComponent));
  }
  if (system.framework.has_value() && system.framework->hasMobileAtoms())
  {
    throw std::runtime_error(std::format("dcTST: '{}' supports rigid frameworks only\n", driverName));
  }
  // the bonds between molecules tie the tagged molecule to others (its moves and the occupancy constraint move one
  // molecule) and change the topology the free-energy profile is taken in
  if (!system.crossLinks.empty() || !system.polymerLinks.empty())
  {
    throw std::runtime_error(std::format(
        "dcTST: '{}' does not support bonds between molecules (cross-links or polymer bonds)\n", driverName));
  }
  if (!integratesTrajectories) return;
  if (system.forceField.computePolarization)
  {
    throw std::runtime_error(std::format("dcTST: '{}' does not support polarization (the trajectories are "
                                         "integrated without polarization forces)\n",
                                         driverName));
  }
  if (system.components[taggedComponent].isSemiFlexible())
  {
    throw std::runtime_error(std::format("dcTST: '{}': the tagged component '{}' is semi-flexible, which is not "
                                         "supported (rigid or fully flexible molecules only)\n",
                                         driverName, system.components[taggedComponent].name));
  }
  if (system.hasExternalField)
  {
    throw std::runtime_error(std::format("dcTST: '{}' does not support an external field (the trajectories are "
                                         "integrated without its forces)\n",
                                         driverName));
  }
}
}  // namespace

DynamicallyCorrectedTST::DynamicallyCorrectedTST(InputReader& reader)
    : random(reader.randomSeed),
      settings(reader.dctstSettings),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      printEvery(std::max(1uz, reader.printEvery)),
      optimizeMCMovesEvery(std::max(1uz, reader.optimizeMCMovesEvery)),
      numberOfBlocks(reader.numberOfBlocks),
      numberOfReplicas(std::max(1uz, reader.dctstSettings.numberOfReplicas))
{
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  barrier = DcTSTBarrier(templateSystem.simulationBox, settings.barrierPosition, settings.barrierNormal,
                         settings.latticeDistance);
  validateDcTSTSystem(templateSystem, settings.taggedComponent, "DynamicallyCorrectedTST", true);
  // the moves the drivers perform, and a tagged molecule that can change its orientation and conformation
  DcTSTMoves::validate(templateSystem, settings.taggedComponent, "DynamicallyCorrectedTST",
                       settings.reinsertionInPlaneProbability, settings.cageOccupancies.has_value());


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
    cageLattice = DcTSTCageLattice(templateSystem.framework->simulationBox, templateSystem.framework->numberOfUnitCells,
                                   lattice_origin, settings.cagesPerUnitCell);
  }
  else
  {
    cageLattice = DcTSTCageLattice(templateSystem.simulationBox, int3(1, 1, 1), lattice_origin);
  }

  // the cages must partition space along the barrier: the normal perpendicular to the cage planes of one axis
  // (the reciprocal vector, not the lattice vector, in a non-orthogonal cell), 'LatticeDistance' their spacing
  // (one cage per unit cell unless 'CagesPerUnitCell' says otherwise) and the barrier midway between centers.
  // The occupancy-resolved and event-terminated modes are built on the cages; otherwise the occupancy classes
  // are left out. The event-terminated estimator classifies the hops by cage axis; with two cells along the
  // normal the column A', A, B, B' closes on itself (A' = B, B' = A), with three A' = B'.
  cageLatticeInconsistency = cageLattice.inconsistencyWith(templateSystem.simulationBox, barrier);
  if (!templateSystem.framework.has_value() && !cageLatticeInconsistency.has_value())
  {
    cageLatticeInconsistency = "no framework";
  }
  if (cageLatticeInconsistency.has_value() &&
      (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event ||
       settings.cageOccupancies.has_value()))
  {
    throw std::runtime_error(std::format("dcTST: {} requires a cage lattice consistent with the barrier: {}\n",
                                         settings.cageOccupancies.has_value() ? "'CageOccupancies'"
                                                                              : "'TrajectoryTermination' \"Event\"",
                                         cageLatticeInconsistency.value()));
  }

  // the settings a programmatic construction can get wrong (the input reader checks the keywords it reads)
  if (!(settings.trajectoryTime > 0.0) || !(settings.trajectoryTime < 1.0e300))
  {
    throw std::runtime_error("dcTST: 'TrajectoryTime' must be positive and finite [ps]\n");
  }
  // The tagged molecule is committed to cage A or B, by the side of the plane it is on, once it is the commit
  // distance from the plane. Beyond half a lattice distance it has passed the center of that cage, and at a
  // lattice distance it is in the next one: the cage it was assigned, and the offset from its center the hops
  // are measured by, were then those of the wrong cage.
  if (!(settings.commitDistance >= 0.0) || !(settings.commitDistance <= 0.5 * barrier.lambda * (1.0 + 1e-12)))
  {
    throw std::runtime_error(std::format(
        "dcTST: 'CommitDistance' ({} [A]) must lie between 0 (the default, half a lattice distance) and half the "
        "lattice distance ({} [A])\n",
        settings.commitDistance, 0.5 * barrier.lambda));
  }
  if (settings.plateauWindow.has_value())
  {
    const auto& [begin, end] = settings.plateauWindow.value();
    if (!(begin >= 0.0) || !(end >= begin) || !(end <= settings.trajectoryTime * (1.0 + 1e-12)))
    {
      throw std::runtime_error(std::format(
          "dcTST: 'PlateauWindow' [{}, {}] needs 0 <= begin <= end <= 'TrajectoryTime' ({} [ps]); a window outside "
          "the trajectories was replaced by their last point without a word\n",
          begin, end, settings.trajectoryTime));
    }
  }
  // n_A counts the tagged molecule: zero wrapped n_A - 1 around to the largest unsigned integer
  if (settings.cageOccupancies.has_value() && settings.cageOccupancies->first == 0uz)
  {
    throw std::runtime_error(
        "dcTST: 'CageOccupancies' n_A must be at least 1 (the tagged molecule itself is in cage A)\n");
  }

  // The handover times belong to the event-terminated trajectories. The input reader drops the default list for
  // "FixedTime"; settings filled in by hand kept it, and the hybrid estimator was then switched on with
  // accumulators that are never filled (a rate of zero and a concerted fraction of -1 in the output).
  if (settings.trajectoryTermination != DcTSTSettings::TrajectoryTermination::Event)
  {
    settings.handoverTimes.clear();
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
  // every event order k = 1..n is accumulated from the same trajectories in event mode
  const bool eventTermination = settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event;
  const std::size_t numberOfEventOrders = eventTermination ? settings.terminationEvents : 1uz;
  transmission.assign(numberOfReplicas,
                      DcTSTTransmissionAccumulator(numberOfBlocks, numberOfShells, numberOfComponents, numberOfTimes,
                                                   settings.handoverTimes.size(), numberOfEventOrders));
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
  loadOccupancyRates();
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  // a termination signal during the initialization skips the transmission stage (its constrained equilibration
  // and the set-up of the barrier states included); the analysis of what was completed is still written
  if (!GracefulShutdown::requested())
  {
    runStage(SimulationStage::Transmission,
             settings.numberOfConstrainedEquilibrationCycles + settings.numberOfTransmissionCycles);
  }
  output();
}

void DynamicallyCorrectedTST::setup()
{
  for (System& system : systems)
  {
    // FIXED terms are held with SHAKE/RATTLE: refuse what cannot be held before the first step
    system.requireSupportedDynamicsConstraints();
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

    // the free-energy profile beta F(q) for the committor (the smooth kappa_S and the hybrid estimator)
    if (results.contains("profile") && results["profile"].contains("q") && results["profile"].contains("freeEnergy"))
    {
      std::vector<double> q, betaF;
      const nlohmann::json& profile = results["profile"];
      for (std::size_t i = 0; i < std::min(profile["q"].size(), profile["freeEnergy"].size()); ++i)
      {
        const double f = profile["freeEnergy"][i].get<double>();
        if (f < 0.0) continue;  // undefined bin (-1.0)
        q.push_back(profile["q"][i].get<double>());
        betaF.push_back(f);
      }
      // a conditioned profile runs over cage A, its q measured from the centre of A with the window at lambda/2
      const double barrierOffset = settings.cageOccupancies.has_value() ? 0.5 * barrier.lambda : 0.0;
      buildCommittorTable(q, betaF, barrierOffset);
    }

    // the occupancy distributions P(n_A, n_B) (pilot) and P(n_A, n_B | q*) (window stage at the plane): their
    // ratio is the conditioned k^TST(n_A -> n_B) / k^TST of the occupancy classes
    if (results.contains("cageOccupancies"))
    {
      const nlohmann::json& occupancies = results["cageOccupancies"];
      auto readMap = [&](const char* key, std::map<std::pair<std::size_t, std::size_t>, double>& target)
      {
        if (!occupancies.contains(key) || !occupancies[key].is_object()) return;
        for (const auto& [pairKey, value] : occupancies[key].items())
        {
          const std::size_t comma = pairKey.find(',');
          if (comma == std::string::npos) continue;
          target[{std::stoul(pairKey.substr(0, comma)), std::stoul(pairKey.substr(comma + 1))}] = value.get<double>();
        }
      };
      readMap("joint", occupancyWeights);
      readMap("jointAtBarrier", barrierOccupancyWeights);
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

  // the committor (smooth kappa_S; hybrid estimator)
  if (!committorTable.empty())
  {
    auto quantile = [&](double p)
    {
      const double width = barrier.lambda / static_cast<double>(committorTable.size());
      for (std::size_t i = 0; i < committorTable.size(); ++i)
      {
        if (committorTable[i] >= p) return -0.5 * barrier.lambda + (static_cast<double>(i) + 0.5) * width;
      }
      return 0.5 * barrier.lambda;
    };
    std::print(stream, "Committor p_B(q) from the free-energy profile ({} points over [-lambda/2, lambda/2)):\n",
               committorTable.size());
    std::print(stream, "    p_B = 0.1, 0.5, 0.9 at q =           {:.3f}, {:.3f}, {:.3f} [A] (barrier at 0)\n",
               quantile(0.1), quantile(0.5), quantile(0.9));
    std::print(stream, "    (used for the committor-smoothed kappa_S and the hybrid collective excess)\n\n");
  }
  else
  {
    std::print(stream, "Committor p_B(q):                        no free-energy profile; the smooth kappa_S and the "
                       "hybrid estimator are not available\n\n");
  }
  if (!occupancyWeights.empty())
  {
    std::print(stream, "Occupancy distributions of the free-energy input: P(n_A, n_B) over {} pairs", occupancyWeights.size());
    if (!barrierOccupancyWeights.empty())
    {
      std::print(stream, ", P(n_A, n_B | q*) over {} pairs\n    (their ratio is k^TST(n_A -> n_B)/k^TST, combined with "
                         "the occupancy classes of this run into the rates k(n -> n'))\n\n",
                 barrierOccupancyWeights.size());
    }
    else
    {
      std::print(stream, "; no P(n_A, n_B | q*) (older free-energy run): the occupancy classes give kappa(n, n') only\n\n");
    }
  }

  // hybrid estimator: handover times in steps; needs the committor, i.e. a free-energy profile
  handoverSteps.clear();
  if (!settings.handoverTimes.empty())
  {
    if (committorTable.empty())
    {
      std::print(stream,
                 "WARNING: 'HandoverTimes' given but no free-energy profile is available for the committor "
                 "('FreeEnergyInput' with 'results.profile' is required); the hybrid estimator is disabled\n\n");
    }
    else if (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event &&
             cageLattice.numberOfUnitCells[cageLattice.axisAlongDirection(barrier.normal).axis] < 4)
    {
      // the Markov model of the hybrid estimator resolves the column A', A, B, B' as four distinct cages; on a
      // ring of two (A' = B, B' = A) or three (A' = B') cages along the normal its occupancies alias
      std::print(stream,
                 "WARNING: 'HandoverTimes' given but the cage lattice has fewer than four cages along the barrier "
                 "normal (the column A', A, B, B' would alias); the hybrid estimator is disabled\n\n");
    }
    else if (settings.trajectoryTermination == DcTSTSettings::TrajectoryTermination::Event &&
             (cageLattice.numberOfUnitCells[(cageLattice.axisAlongDirection(barrier.normal).axis + 1) % 3] < 2 ||
              cageLattice.numberOfUnitCells[(cageLattice.axisAlongDirection(barrier.normal).axis + 2) % 3] < 2))
    {
      // a cage that is alone along a lateral axis is its own lateral neighbour: the neighbourhood state held n_A
      // as a lateral occupancy of A, and the Markov closure counted hops from A into A as channels that end the
      // state
      std::print(stream,
                 "WARNING: 'HandoverTimes' given but the cage lattice has a single cage along an axis beside the "
                 "barrier normal (a cage would be its own lateral neighbour); the hybrid estimator is disabled\n\n");
    }
    else
    {
      const double timeStep = front.timeStep;
      for (double tau : settings.handoverTimes)
      {
        const std::size_t steps = std::max(1uz, static_cast<std::size_t>(std::llround(tau / timeStep)));
        if (steps < numberOfTrajectorySteps) handoverSteps.push_back(steps);
      }
      std::sort(handoverSteps.begin(), handoverSteps.end());
      handoverSteps.erase(std::unique(handoverSteps.begin(), handoverSteps.end()), handoverSteps.end());
      if (handoverSteps.size() != settings.handoverTimes.size())
      {
        // the accumulators were sized with settings.handoverTimes; keep the two lists aligned
        settings.handoverTimes.clear();
        for (std::size_t steps : handoverSteps) settings.handoverTimes.push_back(static_cast<double>(steps) * timeStep);
        for (DcTSTTransmissionAccumulator& accumulator : transmission)
        {
          accumulator.hybridTagged.assign(handoverSteps.size(), std::vector<double>(numberOfBlocks, 0.0));
          accumulator.hybridOthers.assign(handoverSteps.size(), std::vector<double>(numberOfBlocks, 0.0));
          accumulator.hybridOccupancy.assign(handoverSteps.size(), {});
          accumulator.hybridNeighbourhood.assign(handoverSteps.size(), {});
          accumulator.hybridOccupancySigned.assign(handoverSteps.size(), {});
          accumulator.hybridNeighbourhoodSigned.assign(handoverSteps.size(), {});
          accumulator.hybridEventBefore.assign(handoverSteps.size(), 0uz);
          accumulator.hybridTruncated.assign(handoverSteps.size(), 0uz);
        }
      }
      std::print(stream, "Hybrid estimator of the collective excess (committor-smoothed displacements at the handover "
                         "times, Markov closure of the race):\n");
      std::print(stream, "    handover times after commitment:     ");
      for (std::size_t steps : handoverSteps) std::print(stream, "{:g} ", static_cast<double>(steps) * timeStep);
      std::print(stream, "[ps]\n\n");
      outputJson["dcTST"]["handoverTimes"] = settings.handoverTimes;
    }
  }
  std::flush(stream);
}

void DynamicallyCorrectedTST::buildCommittorTable(const std::vector<double>& qGrid, const std::vector<double>& betaF,
                                                  double barrierOffset)
{
  committorTable.clear();
  if (qGrid.size() < 4uz) return;

  // beta F on a fine uniform grid over [-lambda/2, lambda/2) around the barrier (the range of the profile is not
  // wrapped, so a conditioned profile over cage A and a bit of B keeps its asymmetry): at every bin center inside
  // the range of the profile, the linear interpolation of its points. The points were dropped into the bins they
  // fall in instead; those of the free-energy driver lie on bin edges, which moved the profile, and the committor
  // with it, by half a bin. Bins that neither the profile nor its mirror image reaches take the lowest beta F of
  // the profile (the cage interiors, where e^{beta F} is negligible against the barrier region anyway).
  constexpr std::size_t numberOfPoints = 400;
  const double width = barrier.lambda / static_cast<double>(numberOfPoints);
  std::vector<std::pair<double, double>> profile;
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < qGrid.size(); ++i)
  {
    const double q = qGrid[i] - barrierOffset;
    if (q < -0.5 * barrier.lambda || q >= 0.5 * barrier.lambda) continue;
    profile.emplace_back(q, betaF[i]);
    minimum = std::min(minimum, betaF[i]);
  }
  if (profile.size() < 2uz) return;
  std::sort(profile.begin(), profile.end());

  // a profile that does not cover the whole period around the barrier (the conditioned profile over cage A
  // ends a quarter cage into B) takes the missing bins from their mirror images about the barrier, so that the
  // committor does not drop to the minimum discontinuously where the profile ends
  std::vector<double> f(numberOfPoints, minimum);
  std::vector<std::size_t> count(numberOfPoints, 0uz);
  {
    // (a bin center within half a bin of the first or the last point takes the value of that point)
    std::size_t segment{0};
    for (std::size_t i = 0; i < numberOfPoints; ++i)
    {
      const double x = -0.5 * barrier.lambda + (static_cast<double>(i) + 0.5) * width;
      if (x < profile.front().first - 0.5 * width || x > profile.back().first + 0.5 * width) continue;
      while (segment + 2uz < profile.size() && profile[segment + 1uz].first < x) ++segment;
      const auto& [qa, fa] = profile[segment];
      const auto& [qb, fb] = profile[segment + 1uz];
      const double t = (qb - qa > 1e-12) ? std::clamp((x - qa) / (qb - qa), 0.0, 1.0) : 0.0;
      f[i] = (1.0 - t) * fa + t * fb;
      count[i] = 1uz;
    }
  }
  // linear interpolation between the filled bins
  const auto interpolate = [&]() -> bool
  {
    std::vector<std::size_t> filled;
    for (std::size_t i = 0; i < numberOfPoints; ++i)
    {
      if (count[i] > 0uz) filled.push_back(i);
    }
    if (filled.size() < 2uz) return false;
    for (std::size_t k = 0; k + 1 < filled.size(); ++k)
    {
      const std::size_t a = filled[k], b = filled[k + 1];
      for (std::size_t i = a + 1; i < b; ++i)
      {
        const double t = static_cast<double>(i - a) / static_cast<double>(b - a);
        f[i] = (1.0 - t) * f[a] + t * f[b];
        count[i] = 1uz;
      }
    }
    return true;
  };
  // The mirror images only for the bins beyond the range of the profile. They were taken for every bin without
  // a point of its own: a profile of 100 points fills every fourth of the 400 bins, the bin before each of
  // them has a measured mirror image and took beta F(-q), and the two bins between were interpolated from
  // beta F(q) to beta F(-q). Under the exponential the higher shoulder of an asymmetric barrier then counted
  // on both sides and p_B(0) went to 1/2.
  const std::vector<std::size_t> known = count;
  for (std::size_t i = 0; i < numberOfPoints; ++i)
  {
    const std::size_t mirror = numberOfPoints - 1uz - i;
    if (known[i] == 0uz && known[mirror] > 0uz)
    {
      f[i] = f[mirror];
      count[i] = 1uz;
    }
  }
  // (what neither the profile nor its mirror image covers, between the two)
  if (!interpolate()) return;

  // p_B(q) = int_{-lambda/2}^{q} e^{beta F} / int e^{beta F}, evaluated at the bin centres (midpoint rule)
  const double reference = *std::max_element(f.begin(), f.end());
  committorTable.assign(numberOfPoints, 0.0);
  double cumulative{0.0};
  for (std::size_t i = 0; i < numberOfPoints; ++i)
  {
    const double weight = std::exp(f[i] - reference);
    committorTable[i] = cumulative + 0.5 * weight;
    cumulative += weight;
  }
  for (double& p : committorTable) p /= cumulative;
}

double DynamicallyCorrectedTST::smoothCageCoordinate(double q) const
{
  // u = q / lambda - 1/2 is zero at the centre of cage B and -1 at the centre of cage A; m = floor(u) selects
  // the window between cages m and m + 1 (at q = (m + 1) lambda), and the committor across that window
  // interpolates the step of the net cage count floor(q / lambda) + 1/2 from m + 1/2 to m + 3/2
  const double u = q / barrier.lambda - 0.5;
  const double m = std::floor(u);
  const double reduced = (u - m - 0.5) * barrier.lambda;  // in [-lambda/2, lambda/2), zero at the window
  const double width = barrier.lambda / static_cast<double>(committorTable.size());
  const std::size_t bin = std::min(committorTable.size() - 1uz,
                                   static_cast<std::size_t>(std::floor((reduced + 0.5 * barrier.lambda) / width)));
  return m + 0.5 + committorTable[bin];
}

double DynamicallyCorrectedTST::OccupancyRates::rate(std::size_t n, std::size_t nPrime, bool tableFirst) const
{
  if (tableFirst)
  {
    auto it = measured.find({n, nPrime});
    if (it != measured.end()) return it->second.rate;
  }
  if (!fitValid) return 0.0;
  return std::exp(fitConstant + fitDeparture * static_cast<double>(n) - fitArrival * static_cast<double>(nPrime));
}

void DynamicallyCorrectedTST::loadOccupancyRates()
{
  externalRates.clear();
  if (settings.occupancyRateInputs.empty()) return;

  std::print(stream, "Occupancy-resolved rates of conditioned runs ('OccupancyRateInputs'; they override the classes of "
                     "this run):\n");
  for (const std::string& path : settings.occupancyRateInputs)
  {
    std::ifstream file(path);
    if (!file.is_open())
    {
      throw std::runtime_error(std::format("dcTST: cannot open 'OccupancyRateInputs' entry '{}'\n", path));
    }
    nlohmann::json json;
    try
    {
      json = nlohmann::json::parse(file);
    }
    catch (const std::exception& e)
    {
      throw std::runtime_error(std::format("dcTST: 'OccupancyRateInputs' entry '{}' is not valid JSON: {}\n", path, e.what()));
    }
    if (!json.contains("dcTST") || !json["dcTST"].contains("cageOccupancies") || !json.contains("results"))
    {
      throw std::runtime_error(std::format(
          "dcTST: 'OccupancyRateInputs' entry '{}' is not the JSON output of an occupancy-conditioned transmission run "
          "('dcTST.cageOccupancies' missing)\n",
          path));
    }
    const std::size_t n = json["dcTST"]["cageOccupancies"][0].get<std::size_t>();
    const std::size_t nPrime = json["dcTST"]["cageOccupancies"][1].get<std::size_t>();
    // the rate k_(A->B) = kappa_S k^TST is in results[<P(q*) source>]["kSelfSmooth"/"kSelf"]
    std::optional<double> k;
    for (const auto& [key, value] : json["results"].items())
    {
      if (!value.is_object()) continue;
      if (value.contains("kSelfSmooth") && value["kSelfSmooth"].is_number() && value["kSelfSmooth"].get<double>() > 0.0)
      {
        k = value["kSelfSmooth"].get<double>();
        break;
      }
      if (value.contains("kSelf") && value["kSelf"].is_number())
      {
        k = value["kSelf"].get<double>();
        break;
      }
    }
    if (!k.has_value() || k.value() <= 0.0)
    {
      throw std::runtime_error(std::format(
          "dcTST: 'OccupancyRateInputs' entry '{}' has no positive rate k_(A->B) ('results.<source>.kSelf'; the run "
          "needs a 'FreeEnergyInput')\n",
          path));
    }
    OccupancyRates::Entry entry;
    entry.rate = k.value();
    entry.weight = 2.0 * static_cast<double>(json["results"].value("numberOfTrajectoryPairs", 500uz));
    entry.source = path;
    if (externalRates.contains({n, nPrime}))
    {
      std::print(stream, "    WARNING: k({} -> {}) given more than once; the last one is used\n", n, nPrime);
    }
    externalRates[{n, nPrime}] = entry;
    std::print(stream, "    k({:>2} -> {:>2}) = {:.6e} [1/ps]   ({})\n", n, nPrime, k.value(), path);
  }
  std::print(stream, "\n");
  std::flush(stream);
}

namespace
{
/// Weighted least-squares fit log k = c + a n - b n' through the measured rates (three parameters; needs three
/// pairs that are not collinear).
void fitOccupancyRates(DynamicallyCorrectedTST::OccupancyRates& rates)
{
  rates.fitValid = false;
  if (rates.measured.size() < 3uz) return;
  double S[3][3]{}, r[3]{};
  for (const auto& [pair, entry] : rates.measured)
  {
    const double x[3] = {1.0, static_cast<double>(pair.first), -static_cast<double>(pair.second)};
    const double y = std::log(entry.rate);
    const double w = std::max(entry.weight, 1e-12);
    for (std::size_t i = 0; i < 3; ++i)
    {
      r[i] += w * x[i] * y;
      for (std::size_t j = 0; j < 3; ++j) S[i][j] += w * x[i] * x[j];
    }
  }
  auto determinant = [](const double M[3][3])
  {
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  };
  const double det = determinant(S);
  const double scale = std::abs(S[0][0] * S[1][1] * S[2][2]);
  if (std::abs(det) <= 1e-10 * std::max(scale, 1e-300)) return;
  auto solve = [&](std::size_t column)
  {
    double M[3][3];
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j) M[i][j] = (j == column) ? r[i] : S[i][j];
    return determinant(M) / det;
  };
  rates.fitConstant = solve(0);
  rates.fitDeparture = solve(1);
  rates.fitArrival = solve(2);
  rates.fitValid = true;
  double residual{0.0}, weightSum{0.0};
  for (const auto& [pair, entry] : rates.measured)
  {
    const double predicted = rates.fitConstant + rates.fitDeparture * static_cast<double>(pair.first) -
                             rates.fitArrival * static_cast<double>(pair.second);
    residual += entry.weight * std::pow(std::log(entry.rate) - predicted, 2);
    weightSum += entry.weight;
  }
  rates.fitRms = std::sqrt(residual / std::max(weightSum, 1e-300));
}

// An accepted move of the tagged molecule commits its new electric field and, with molecule-molecule polarization,
// the change of the field on every other atom into the stored field (the same commit as the standard moves).
void commitElectricField(System& system, std::size_t selectedComponent, std::size_t selectedMolecule,
                         std::span<const double3> electricFieldMoleculeNew,
                         std::span<const double3> electricFieldNeighborDelta)
{
  if (!system.forceField.omitInterPolarization)
  {
    std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
    for (std::size_t i = 0; i < storedElectricField.size(); ++i)
    {
      storedElectricField[i] += electricFieldNeighborDelta[i];
    }
  }
  std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
  std::copy(electricFieldMoleculeNew.begin(), electricFieldMoleculeNew.end(), electricFieldMolecule.begin());
}
}  // namespace


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

  std::vector<std::exception_ptr> workerErrors(numberOfReplicas);
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          keepingException(workerErrors[replicaId], [this, replicaId, stage, numberOfCycles, equilibrationCycles, productionCycles]()
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

              // the step sizes are tuned before the sampled cycles only: an adaptive kernel while the
              // barrier ensemble is sampled breaks detailed balance
              if (!(stage == SimulationStage::Transmission && production) && cycle % optimizeMCMovesEvery == 0uz)
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
                        if (c != settings.taggedComponent) continue;  // kappa_C is of the tagged species
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
          }));
    }
  }
  rethrowWorkerException(workerErrors);

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
      // a component without any move (all probabilities zero, e.g. immobile cations) samples 'Count'
      if (moveType == Move::Types::Count) continue;
      component.mc_moves_statistics.addAllCounts(moveType);

      // occupancy-resolved mode: the old configuration of the molecule, to undo an accepted move that changes
      // the occupancy of cage A or B (the constraint is a hard wall in configuration space: rejecting such
      // moves keeps detailed balance since the allowed set is the same seen from either side)
      const std::size_t moleculeIndex =
          occupancyMode ? system.moleculeIndexOfComponent(selectedComponent, selectedMolecule) : 0uz;
      std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
      std::vector<Atom> oldAtoms;
      Molecule oldMolecule;
      std::vector<double3> oldElectricField;
      std::variant<MoveStatistics<double>, MoveStatistics<double3>> oldMoveStatistics;
      SubMoveStatistics oldSubMoveStatistics;
      if (occupancyMode)
      {
        oldMoveStatistics = component.mc_moves_statistics.stats[std::to_underlying(moveType)];
        oldSubMoveStatistics = component.mc_moves_statistics.subMoveStatistics(moveType);
        oldAtoms.assign(moleculeAtoms.begin(), moleculeAtoms.end());
        oldMolecule = system.moleculeData[moleculeIndex];
        // an accepted move commits the electric field (of the molecule and, with molecule-molecule polarization,
        // of every other molecule) into the stored field: keep a copy so that the undo below can restore it
        if (system.forceField.computePolarization)
        {
          oldElectricField = system.electricField;
        }
      }

      // every single-molecule move of the component, the conformational ones of flexible molecules included
      // ('DcTSTMoves::validate' refuses the rest up front; there is no fall-back to a translation)
      energyDifference = DcTSTMoves::perform(rng, system, selectedComponent, selectedMolecule, moveType);

      // (the constraint is on the molecules of the tagged component; the others move freely between the cages)
      if (occupancyMode && energyDifference.has_value() && selectedComponent == settings.taggedComponent)
      {
        DcTSTOccupancyState& state = occupancyStates[replicaId];
        const std::size_t oldCage = state.cageOfMolecule[moleculeIndex];
        const std::size_t newCage = cageLattice.cageOf(moleculeAtoms[component.startingBead].position);
        if (newCage != oldCage)
        {
          ++state.cageChanges;
          if (oldCage == state.cageA || oldCage == state.cageB || newCage == state.cageA || newCage == state.cageB)
          {
            // undo: positions, molecule record, cell list, the Ewald structure factor and the stored electric field
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
            if (system.forceField.computePolarization)
            {
              system.electricField = oldElectricField;
            }
            energyDifference = std::nullopt;
            component.mc_moves_statistics.revokeAccepted(moveType, oldMoveStatistics, oldSubMoveStatistics);
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

  // Polarization follows the bookkeeping of MC_Moves::translationMove (and DynamicallyCorrectedTSTFreeEnergy): the
  // field on the moved molecule (new and old) from the framework, the reciprocal space and, with molecule-molecule
  // polarization, the other molecules, plus the change of the field on every other atom; an accepted move commits
  // both into the stored field so that the running polarization energy stays consistent with a full recomputation.
  const bool polarization = system.forceField.computePolarization;
  const bool interPolarization = polarization && !system.forceField.omitInterPolarization;
  std::vector<double3> electricFieldMoleculeNew(moleculeAtoms.size());
  std::vector<double3> electricFieldMoleculeOld(moleculeAtoms.size());
  std::vector<double3> electricFieldNeighborDelta;

  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialMolecule.second, moleculeAtoms);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkMolecule =
      polarization
          ? Interactions::computeFrameworkMoleculeEnergyDifference(
                system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second,
                moleculeAtoms)
          : Interactions::computeFrameworkMoleculeEnergyDifference(
                system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                system.frameworkCellList(), system.spanOfFrameworkAtoms(), trialMolecule.second, moleculeAtoms);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interMolecule;
  if (interPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), trialMolecule.second, moleculeAtoms);
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(system.forceField, system.simulationBox,
                                                                        system.cellList(), system.spanOfMoleculeAtoms(),
                                                                        trialMolecule.second, moleculeAtoms);
  }
  if (!interMolecule.has_value()) return std::nullopt;

  const RunningEnergy ewaldFourierEnergy =
      polarization
          ? Interactions::energyDifferenceEwaldFourier(
                system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
                system.storedEik, system.trialEik, system.forceField, system.simulationBox, system.components,
                electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, moleculeAtoms)
          : Interactions::energyDifferenceEwaldFourier(
                system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                system.forceField, system.simulationBox, system.components, trialMolecule.second, moleculeAtoms);

  RunningEnergy polarizationDifference;
  if (polarization)
  {
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, moleculeAtoms);
    if (interPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  const RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() +
                                         interMolecule.value() + ewaldFourierEnergy + polarizationDifference;

  if (rng.uniform() < MC_Moves::boundedExp(-system.beta * energyDifference.potentialEnergy()))
  {
    ++statistics.translationAccepted;
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), moleculeAtoms.begin());
    system.cellListAtomsMoved(moleculeAtoms);
    molecule = trialMolecule.first;
    if (polarization)
    {
      commitElectricField(system, selectedComponent, 0uz, electricFieldMoleculeNew, electricFieldNeighborDelta);
    }
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

  // Polarization follows the bookkeeping of MC_Moves::reinsertionMove (and DynamicallyCorrectedTSTFreeEnergy): the
  // field on the reinserted molecule (new and old) from the framework, the reciprocal space and, with
  // molecule-molecule polarization, the other molecules (whose inter-molecular energy is in the Rosenbluth weights
  // already, so that energy is discarded), plus the change of the field on every other atom; an accepted move
  // commits both into the stored field.
  const bool polarization = system.forceField.computePolarization;
  const bool interPolarization = polarization && !system.forceField.omitInterPolarization;
  std::vector<double3> electricFieldMoleculeNew(moleculeAtoms.size());
  std::vector<double3> electricFieldMoleculeOld(moleculeAtoms.size());
  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (polarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(
        system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), electricFieldMoleculeNew,
        electricFieldMoleculeOld, growData->atoms, moleculeAtoms);
    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, electricFieldMoleculeOld,
        growData->atoms, moleculeAtoms);
    if (interPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] const std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
              electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), growData->atoms, moleculeAtoms);
    }
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, growData->atoms, moleculeAtoms);
    if (interPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  const double logAcceptance =
      -system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy()) +
      growData->logRosenbluthWeight - retraceData.logRosenbluthWeight;

  if (rng.uniform() < MC_Moves::boundedExp(logAcceptance))
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
    if (polarization)
    {
      commitElectricField(system, selectedComponent, 0uz, electricFieldMoleculeNew, electricFieldNeighborDelta);
    }
    if (state.active)
    {
      state.cageA = trialCages.first;
      state.cageB = trialCages.second;
    }

    return (growData->energies - retraceData.energies) + energyFourierDifference + polarizationDifference;
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
  // Not every component: extra-framework cations, and any component without a move, stay where the system has
  // them. Such molecules were deleted and regrown at random positions with the rest, and a component without
  // moves then stayed there for the whole run (a cation off its site, in a window as soon as not).
  const auto keptInPlace = [&](std::size_t componentId)
  {
    if (componentId == taggedComponent) return false;
    const Component& component = system.components[componentId];
    if (component.type == Component::Type::Cation) return true;
    for (std::size_t index = 0; index < std::to_underlying(Move::Types::Count); ++index)
    {
      if (component.mc_moves_probabilities.getProbability(static_cast<Move::Types>(index)) > 0.0) return false;
    }
    return true;
  };
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (keptInPlace(componentId)) continue;
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
  {
    return system.makeGrowContext(CBMC::CutOffMode::Full)
        .withChainScheme(CBMC::ChainScheme::ConfigurationalBias)
        .forPlacement();
  };

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
    // Halfway along the minimum image of the step to the neighbour, and on this plane only. The mean of the two
    // centers is the middle of the supercell for the pair across its boundary (the windows of the plane q = 0
    // where that is the seam, all of them refused), and the midpoints of the other planes q = k lambda were
    // moved here along the normal: k lattice distances, which in an oblique cell is not a lattice translation
    // and put the window point into the wall.
    const double3 centerOfCage = cageLattice.center(c);
    const double3 mid = centerOfCage + 0.5 * system.simulationBox.applyPeriodicBoundaryConditions(
                                                 cageLattice.center(cageLattice.neighborAlongAxis(
                                                     c, alignment.axis, alignment.sign)) -
                                                 centerOfCage);
    if (std::abs(barrier.projection(system.simulationBox, mid)) < 0.1 * barrier.lambda)
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

    // (bounded, and interruptible: at a loading the framework cannot hold the loop never ended)
    constexpr std::size_t maximumGrowthAttempts = 1000000uz;
    std::optional<CBMC::GrowResult> current = std::nullopt;
    for (std::size_t attempt = 0; !isValid(current); ++attempt)
    {
      if (attempt >= maximumGrowthAttempts)
      {
        throw std::runtime_error(std::format(
            "dcTST: molecule {} of component '{}' could not be grown after {} attempts (replica {}); the loading "
            "does not fit next to the tagged molecule on the barrier\n",
            system.numberOfMoleculesPerComponent[componentId] + 1uz, component.name, maximumGrowthAttempts,
            replicaId));
      }
      if ((attempt & 255uz) == 255uz && GracefulShutdown::requested())
      {
        throw std::runtime_error("dcTST: interrupted by a termination signal while creating the configuration\n");
      }
      current = CBMC::growNewMolecule(rng, creationContext(), component, identity);
    }

    const std::size_t numberOfGrowths = component.rigid ? 1uz : CBMC::Constants::creationCandidateGrowths;
    for (std::size_t attempt = 1; attempt < numberOfGrowths; ++attempt)
    {
      std::optional<CBMC::GrowResult> candidate = CBMC::growNewMolecule(rng, creationContext(), component, identity);
      if (!isValid(candidate)) continue;

      const double logAcceptance = candidate->logRosenbluthWeight - current->logRosenbluthWeight;
      if (logAcceptance >= 0.0 || rng.uniform() < MC_Moves::boundedExp(logAcceptance))
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
  // the bead is on the plane, the boundary between A and B: the two cages of its column on either side of it
  // (a quarter-cage step along the normal to either side left the column in an oblique cell, for a bead near a
  // lateral cage wall)
  return cageLattice.cagesAcrossNearestPlane(beadPosition, barrier.normal);
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
    // the occupancies are those of the tagged species: (n_A, n_B) label the rates of its hops, and a cation or
    // a co-adsorbate in the cage is part of the environment, not a molecule that hops with that rate
    if (molecule.componentId != settings.taggedComponent) continue;
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

  // the molecules of the tagged component besides the tagged one: only they are counted and moved (every
  // molecule was, so that cations were taken off their sites to fill a cage and counted as its occupancy)
  const std::size_t numberOfOthers = system.numberOfMoleculesPerComponent[settings.taggedComponent] - 1uz;
  if (nA - 1uz + nB > numberOfOthers)
  {
    throw std::runtime_error(std::format(
        "dcTST: 'CageOccupancies' [{}, {}] needs {} molecules of the tagged component besides the tagged one, the "
        "system has {}\n",
        nA, nB, nA - 1uz + nB, numberOfOthers));
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
      if (j != taggedIndex && molecule.componentId == settings.taggedComponent)
      {
        ++state.occupancy[state.cageOfMolecule[j]];
      }
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
  const CBMC::GrowContext context = system.makeGrowContext(CBMC::CutOffMode::Full).forPlacement();
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
      // (as every other placement of a molecule has it: the growth refuses a bead above the overlap criterion, and
      // the beads of a grown molecule can still add up to more than it)
      if (growData->energies.potentialEnergy() > system.forceField.energyOverlapCriteria) continue;
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
      if (j != taggedIndex && system.moleculeData[j].componentId == settings.taggedComponent &&
          predicate(state.cageOfMolecule[j]))
      {
        candidates.push_back(j);
      }
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

  // The components the Monte Carlo stages do not move (no move with a probability: immobile cations) are held in
  // place in the trajectories too. They were given thermal velocities and integrated with the rest: degrees of
  // freedom the constrained ensemble does not sample, released from positions no move had relaxed (a point-charge
  // cation without a repulsion that holds it off the framework oxygens falls into them).
  std::vector<char> heldComponent(system.components.size(), 0);
  bool anyHeld{false};
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (componentId == settings.taggedComponent) continue;
    bool moves{false};
    for (std::size_t index = 0; index < std::to_underlying(Move::Types::Count); ++index)
    {
      if (system.components[componentId].mc_moves_probabilities.getProbability(static_cast<Move::Types>(index)) > 0.0)
      {
        moves = true;
      }
    }
    heldComponent[componentId] = moves ? 0 : 1;
    anyHeld = anyHeld || (!moves && system.numberOfMoleculesPerComponent[componentId] > 0uz);
  }

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
  // event mode: fractional offset (unit-cell basis) of the reaction bead from the center of its committed cage,
  // accumulated from the minimum-image step displacements so that it is not folded back by the supercell;
  // the hop threshold (half a cell plus the commit distance) can then exceed half the box along an axis
  std::vector<double3> cageOffset(numberOfMolecules);

  std::vector<double> deltaSelf(numberOfTimes);
  std::vector<double> deltaSelfTwoState(numberOfTimes);
  std::vector<double> deltaSelfCommitted(numberOfTimes);
  std::vector<double> deltaSelfSmooth(numberOfTimes);
  std::vector<double> deltaCollective(numberOfShells * numberOfComponents * numberOfTimes);
  const bool smoothMode = !committorTable.empty();

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
  // occupancy classes: the numbers of other molecules in the two cages sharing the tagged window at launch
  std::size_t othersInA{0}, othersInB{0};
  {
    const double3 taggedBead =
        systemAtoms[system.moleculeData[taggedIndex].atomIndex + taggedComponent.startingBead].position;
    const auto [launchA, launchB] = cagesOfWindow(taggedBead);
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      if (j == taggedIndex) continue;
      const Molecule& molecule = system.moleculeData[j];
      // of the tagged species only (see countOccupancies)
      if (molecule.componentId != settings.taggedComponent) continue;
      const Component& component = system.components[molecule.componentId];
      const std::size_t cage = cageLattice.cageOf(systemAtoms[molecule.atomIndex + component.startingBead].position);
      if (cage == launchA) ++othersInA;
      else if (cage == launchB) ++othersInB;
    }
  }
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
      hopThreshold[i] = 0.5 + commitDistance / cageLattice.planeSpacing(static_cast<int>(i));
    }
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      if (j == taggedIndex) continue;  // committed to A or B only once it has reached the commit distance
      const Molecule& molecule = system.moleculeData[j];
      // Only the molecules of the tagged species have a cage: the events, Q_local, the column occupancies and
      // the Markov closure are those of one species hopping with one set of rates k(n -> n'). A cation or a
      // co-adsorbate was counted into (n_A', n_A, n_B, n_B') and its count multiplied by the rate of the tagged
      // species; its window crossings were events. It stays in the per-component collective sums.
      if (molecule.componentId != settings.taggedComponent) continue;
      const Component& component = system.components[molecule.componentId];
      initialCommittedCage[j] = cageLattice.cageOf(systemAtoms[molecule.atomIndex + component.startingBead].position);
    }
  }
  // the window between two cages of the column along the normal: 0 = W(A'|A), 1 = W(A|B), 2 = W(B|B'). The
  // window is identified by the cage the molecule leaves and the direction of the hop along the normal, so
  // that it stays well defined when cages of the column coincide (three cells: A' = B'; two cells: A' = B and
  // B' = A, where windows 0 and 2 are the same far window of the two-cage ring)
  auto columnWindow = [&](std::size_t from, int alongNormal) -> std::optional<std::size_t>
  {
    if (from == cageA) return (alongNormal > 0) ? std::optional<std::size_t>{1uz} : std::optional<std::size_t>{0uz};
    if (from == cageB) return (alongNormal > 0) ? std::optional<std::size_t>{2uz} : std::optional<std::size_t>{1uz};
    if (from == cageAPrime && alongNormal > 0) return 0uz;
    if (from == cageC && alongNormal < 0) return 2uz;
    return std::nullopt;
  };

  // --- hybrid estimator: committor-smoothed displacements and column occupancies at the handover times ----------
  using ColumnState = DcTSTTransmissionAccumulator::ColumnState;
  using NeighbourhoodState = DcTSTTransmissionAccumulator::NeighbourhoodState;
  const bool hybridMode = eventMode && hybridEnabled();
  const std::size_t numberOfHandovers = hybridMode ? handoverSteps.size() : 0uz;
  // the smooth cage coordinate is clamped to the column (A' = -3/2 ... B' = +3/2): displacements beyond it are
  // not column crossings, exactly as in Q_local
  auto smoothCoordinate = [&](double q) { return std::clamp(smoothCageCoordinate(q), -1.5, 1.5); };
  std::vector<double> smoothQ(numberOfMolecules, 0.0);
  struct Handover
  {
    bool reached{false};      // the handover time was reached before the end of the trajectory
    bool eventBefore{false};  // the first post-commitment event occurred before the handover time
    double tagged{0.0};       // tagged cage count: smooth C at tau, or the net count at the first event
    double others{0.0};       // others: smoothed column displacement at tau, or Q_local at the first event
    int orientation{1};       // the side the tagged molecule committed to (+1: along the normal), the frame of 'state'
    ColumnState state{};      // occupancies of the column at tau (frame of the hop), when no event occurred before
    NeighbourhoodState neighbourhood{};  // column plus the four lateral neighbours of A and of B
  };

  // one local axial event after the commitment of the tagged molecule; the state of the trajectory at its k-th
  // event is the k-th order event-terminated estimate, so every event up to the n-th is recorded
  struct PostCommitmentEvent
  {
    EventChannel channel{};
    double frameQ{0.0};  // signed contribution to Q_local in the frame of the tagged hop (zero for the tagged channels)
    double time{0.0};    // since commitment [ps]
    std::size_t step{0};
    double self{0.0};    // net cage count of the tagged molecule just after the event
    double localQ{0.0};  // Q_local (lab frame) just after the event
  };

  struct TrajectoryOutcome
  {
    double meanDrift{0.0};
    double maxDrift{0.0};
    bool aborted{false};  // interrupted by a termination signal; the pair is discarded
    bool longJump{false};
    bool committed{false};
    bool eventTerminated{false};     // reached the n-th event before the cap
    std::size_t terminationStep{0};  // the step of the n-th event, or the last integrated step
    std::size_t integratedSteps{0};  // steps integrated (beyond the termination in hybrid mode)
    std::vector<Handover> handover;  // [handover time]
    // the local axial events after commitment, in order (at most n)
    std::vector<PostCommitmentEvent> events;
    double finalSelf{0.0};    // net cage count of the tagged molecule at termination
    double finalLocalQ{0.0};  // Q_local: signed axial crossings of the other molecules in the column
    // the crossings before the commitment of the tagged molecule, per channel in the frame of the tagged hop:
    // signed contribution to Q_local (others only), number of crossings (tagged included) and the sum of their
    // times since launch [ps]
    std::size_t commitStep{0};
    std::array<double, static_cast<std::size_t>(EventChannel::Count)> frameQBefore{};
    std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)> eventCountBefore{};
    std::array<double, static_cast<std::size_t>(EventChannel::Count)> eventTimeBefore{};

    // the k-th order view of the member: its state at the k-th event when it reached it, else at termination
    struct OrderView
    {
      bool ended{false};  // reached the k-th event before the cap
      double self{0.0};
      double localQ{0.0};
      std::size_t step{0};
      std::optional<EventChannel> channel{};
    };
    [[nodiscard]] OrderView at(std::size_t order) const
    {
      if (order >= 1uz && events.size() >= order)
      {
        const PostCommitmentEvent& event = events[order - 1uz];
        return {true, event.self, event.localQ, event.step, event.channel};
      }
      return {false, finalSelf, finalLocalQ, terminationStep, std::nullopt};
    }
  };

  // integrates one trajectory from the current working copies and fills the delta arrays; reports the
  // energy drift (meaningful for NVE only), whether the tagged molecule left the pair of cages A/B
  // ('long jump') and whether it reached the commit distance; in event mode the trajectory ends at the
  // n-th local axial hop event after commitment (every event up to it is recorded) and the stored values are
  // held from then on
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
    if (anyHeld) holdInPlace(heldComponent, moleculeSpan, dynamicsSpan, groupSpan, system.components);
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
      cageOffset[j] = (eventMode && initialCommittedCage[j] != noCage)
                          ? cageLattice.fractionalOffset(system.simulationBox, previousPosition[j],
                                                         initialCommittedCage[j])
                          : double3{};
    }

    // t = 0: the tagged bead is on the surface, moving into B (direction +1) or A (direction -1)
    std::fill(deltaCollective.begin(), deltaCollective.end(), 0.0);
    deltaSelf[0] = 0.5 * direction;
    deltaSelfTwoState[0] = 0.5 * direction;
    deltaSelfCommitted[0] = 0.5 * direction;
    deltaSelfSmooth[0] = 0.5 * direction;
    outcome = TrajectoryOutcome{};
    outcome.handover.assign(numberOfHandovers, Handover{});
    bool& committed = outcome.committed;
    double committedValue{0.0};
    // the side of the plane the tagged molecule committed to (+1: the cage along the normal, lab B)
    int committedOrientation{0};

    // event mode: the frame of the tagged hop (mirrored for the time-reversed trajectory)
    committedCage = initialCommittedCage;
    double localQ{0.0};
    std::size_t eventsAfterCommitment{0};
    bool terminated{false};  // the terminating event has occurred (the first-event quantities are frozen)
    // hybrid mode: the smoothed column displacement of the other molecules, the values at the first event after
    // commitment, and the next handover time to record
    double smoothOthers{0.0};
    std::optional<std::size_t> firstEventStep{};
    double firstEventTagged{0.0}, firstEventOthers{0.0};
    std::size_t nextHandover{0};
    // hybrid mode: the column index of a cage (A' = -2, A = -1, B = 0, B' = +1 in the lab frame; unique with four
    // or more cages along the normal), the smooth coordinate of the center of such a cage, and the smooth
    // coordinate of another molecule from its displacement along the normal from the center of its committed cage
    auto columnIndex = [&](std::size_t cage) -> std::optional<int>
    {
      if (cage == noCage) return std::nullopt;
      if (cage == cageAPrime) return -2;
      if (cage == cageA) return -1;
      if (cage == cageB) return 0;
      if (cage == cageC) return 1;
      return std::nullopt;
    };
    auto columnCenter = [](int index) { return static_cast<double>(index) + 0.5; };
    auto smoothOtherCoordinate = [&](std::size_t j, int index)
    {
      const double x = double3::dot(cageLattice.cell * cageOffset[j], barrier.normal);
      return smoothCageCoordinate((static_cast<double>(index) + 0.5) * barrier.lambda + x);
    };
    std::vector<char> inColumnBefore(numberOfMolecules, 0);
    std::vector<int> columnIndexBefore(numberOfMolecules, 0);
    if (hybridMode)
    {
      for (std::size_t j = 0; j < numberOfMolecules; ++j)
      {
        if (j == taggedIndex)
        {
          smoothQ[j] = smoothCoordinate(unwrappedQ[j]);
          continue;
        }
        const std::optional<int> index = columnIndex(committedCage[j]);
        inColumnBefore[j] = index.has_value();
        if (index.has_value())
        {
          columnIndexBefore[j] = index.value();
          smoothQ[j] = smoothOtherCoordinate(j, index.value());
        }
      }
    }
    auto frameOrientation = [&]() -> int
    { return (committedOrientation != 0) ? committedOrientation : ((direction > 0.0) ? 1 : -1); };
    auto columnState = [&]() -> ColumnState
    {
      std::size_t count[4]{};  // lab frame: A', A, B, B'
      for (std::size_t j = 0; j < numberOfMolecules; ++j)
      {
        const std::size_t cage = committedCage[j];
        if (cage == noCage) continue;
        if (cage == cageAPrime) ++count[0];
        else if (cage == cageA) ++count[1];
        else if (cage == cageB) ++count[2];
        else if (cage == cageC) ++count[3];
      }
      // The frame of the hop: the tagged molecule is in B. That is the cage it committed to, not the one it was
      // launched towards: a member that turns back before the commit distance sits in the cage it came from,
      // and in the frame of its launch direction the state had the tagged molecule in A (n_B could be zero).
      ColumnState state;
      if (frameOrientation() > 0)
      {
        state.nAPrime = static_cast<std::uint16_t>(count[0]);
        state.nA = static_cast<std::uint16_t>(count[1]);
        state.nB = static_cast<std::uint16_t>(count[2]);
        state.nBPrime = static_cast<std::uint16_t>(count[3]);
      }
      else
      {
        state.nAPrime = static_cast<std::uint16_t>(count[3]);
        state.nA = static_cast<std::uint16_t>(count[2]);
        state.nB = static_cast<std::uint16_t>(count[1]);
        state.nBPrime = static_cast<std::uint16_t>(count[0]);
      }
      return state;
    };
    auto neighbourhoodState = [&]() -> NeighbourhoodState
    {
      NeighbourhoodState state;
      state.column = columnState();
      const int a0 = (axis.axis == 0) ? 1 : 0;
      const int a1 = (axis.axis == 2) ? 1 : 2;
      auto fill = [&](std::size_t cage, std::array<std::uint16_t, 4>& dest)
      {
        const std::array<std::size_t, 4> neighbours{cageLattice.neighborAlongAxis(cage, a0, +1),
                                                    cageLattice.neighborAlongAxis(cage, a0, -1),
                                                    cageLattice.neighborAlongAxis(cage, a1, +1),
                                                    cageLattice.neighborAlongAxis(cage, a1, -1)};
        for (std::size_t i = 0; i < 4; ++i)
        {
          std::size_t count{0};
          for (std::size_t j = 0; j < numberOfMolecules; ++j)
          {
            if (committedCage[j] == neighbours[i]) ++count;
          }
          dest[i] = static_cast<std::uint16_t>(count);
        }
        std::sort(dest.begin(), dest.end());
      };
      // frame of the tagged hop (see columnState): A is the cage the tagged molecule is not in
      if (frameOrientation() > 0)
      {
        fill(cageA, state.lateralA);
        fill(cageB, state.lateralB);
      }
      else
      {
        fill(cageB, state.lateralA);
        fill(cageA, state.lateralB);
      }
      return state;
    };
    auto channelOf = [&](std::size_t j, std::size_t window, int alongNormal) -> EventChannel
    {
      // in the frame of the tagged hop: direction +1 is the direction of the tagged molecule's initial
      // velocity, window 0 is W(A'|A), 1 is W(A|B), 2 is W(B|B') with A the cage it started from
      const int frameDirection = (direction > 0.0) ? alongNormal : -alongNormal;
      const std::size_t frameWindow = (direction > 0.0) ? window : 2uz - window;
      if (j == taggedIndex)
      {
        // relative to its commitment: on from the cage it committed to, or back through the window it left. By
        // the launch direction alone a member that had turned back and committed to A had its recrossing A -> B
        // counted as running on from B to B', and its hop A -> A' as the return from B to A.
        const int commitDirection = (committedOrientation != 0) ? alongNormal * committedOrientation : frameDirection;
        return (commitDirection > 0) ? EventChannel::TaggedForwardHop : EventChannel::TaggedBackHop;
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
      // the committor-smoothed variant: the step of the net cage count replaced by p_B across each window
      if (smoothMode) deltaSelfSmooth[t] = smoothCoordinate(qTagged);
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

      RunningEnergy energies = Integrators::velocityVerlet(
          moleculeSpan, atomSpan, dynamicsSpan, system.components, system.timeStep, noThermostat, frameworkAtoms,
          system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
          system.numberOfMoleculesPerComponent, system.framework, {}, groupSpan, {}, &system.crossLinks);
      if (anyHeld)
      {
        // the second half-kick of the step gave the held molecules a velocity: taken back, with the forces, before
        // the next step moves them, and left out of the kinetic energy of the conserved quantity
        holdInPlace(heldComponent, moleculeSpan, dynamicsSpan, groupSpan, system.components);
        energies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
            moleculeSpan, atomSpan, dynamicsSpan, system.components, system.framework, frameworkAtoms, {},
            &system.forceField, groupSpan);
        energies.rotationalKineticEnergy =
            Integrators::computeRotationalKineticEnergy(moleculeSpan, system.components, groupSpan);
      }

      if (bathType != DcTSTSettings::TrajectoryThermostat::None)
      {
        applyBath(rng, bathType, settings.bathRate, system.timeStep, system.temperature, moleculeSpan, atomSpan,
                  dynamicsSpan, groupSpan, system.components);
        if (anyHeld) holdInPlace(heldComponent, moleculeSpan, dynamicsSpan, groupSpan, system.components);
      }
      else
      {
        // relative to the conserved energy at launch, or absolute when that is zero (a free particle at rest on
        // the barrier): the quotient was 0/0 or x/0 there, which -ffast-math leaves undefined and JSON writes as
        // null
        const double driftScale = (std::abs(conservedReference) > 1.0e-12) ? std::abs(conservedReference) : 1.0;
        const double drift = std::abs(energies.conservedEnergy() - conservedReference) / driftScale;
        driftSum += drift;
        outcome.maxDrift = std::max(outcome.maxDrift, drift);
      }

      for (std::size_t j = 0; j < numberOfMolecules; ++j)
      {
        const Molecule& molecule = moleculeData[j];
        const Component& component = system.components[molecule.componentId];
        const double3 position = atoms[molecule.atomIndex + component.startingBead].position;
        const double3 stepDisplacement = system.simulationBox.applyPeriodicBoundaryConditions(position - previousPosition[j]);
        unwrappedQ[j] += double3::dot(stepDisplacement, barrier.normal);
        if (eventMode) cageOffset[j] += cageLattice.inverseCell * stepDisplacement;
        previousPosition[j] = position;
      }

      // commitment and long-jump detection of the tagged molecule (checked every step)
      const double qTagged = unwrappedQ[taggedIndex];
      if (!committed && std::abs(qTagged) >= commitDistance)
      {
        committed = true;
        committedValue = (qTagged > 0.0) ? 0.5 : -0.5;
        outcome.commitStep = step;
        committedOrientation = (qTagged > 0.0) ? 1 : -1;
        if (eventMode)
        {
          committedCage[taggedIndex] = (qTagged > 0.0) ? cageB : cageA;
          // the tagged bead is within the commit distance of the window, i.e. inside the cage it commits to
          cageOffset[taggedIndex] = cageLattice.fractionalOffset(system.simulationBox, previousPosition[taggedIndex],
                                                                 committedCage[taggedIndex]);
        }
      }
      if (!terminated && (qTagged >= barrier.lambda || qTagged < -barrier.lambda)) outcome.longJump = true;

      // event mode: hop detection of every molecule with a committed cage; a hop along the normal between
      // two cages of the column is a local axial event. After the terminating event (hybrid mode integrates
      // on to the last handover time) the cages are still tracked but no events are counted.
      const bool terminatedBefore = terminated;
      if (eventMode)
      {
        for (std::size_t j = 0; j < numberOfMolecules; ++j)
        {
          if (committedCage[j] == noCage) continue;
          // the axial direction first, then the lateral ones, and every axis crossed in this step (a diagonal
          // crossing used to drop the axial hop when a lateral axis came first)
          const std::array<std::size_t, 3> axisOrder{
              static_cast<std::size_t>(axis.axis), static_cast<std::size_t>(axis.axis == 0 ? 1 : 0),
              static_cast<std::size_t>(axis.axis == 2 ? 1 : 2)};
          for (const std::size_t i : axisOrder)
          {
            const std::size_t cage = committedCage[j];
            const double3& f = cageOffset[j];
            const int hop = (f[i] >= hopThreshold[i]) ? 1 : (f[i] <= -hopThreshold[i]) ? -1 : 0;
            if (hop == 0) continue;
            const std::size_t newCage = cageLattice.neighborAlongAxis(cage, static_cast<int>(i), hop);
            if (static_cast<int>(i) == axis.axis && !terminated)
            {
              const int alongNormal = hop * axis.sign;
              const std::optional<std::size_t> window = columnWindow(cage, alongNormal);
              if (window.has_value())
              {
                const EventChannel channel = channelOf(j, window.value(), alongNormal);
                const std::size_t channelIndex = static_cast<std::size_t>(channel);
                const int frameDirection = (direction > 0.0) ? alongNormal : -alongNormal;
                const double frameQ = (j != taggedIndex) ? static_cast<double>(frameDirection) : 0.0;
                if (j != taggedIndex) localQ += static_cast<double>(alongNormal);
                if (!committed)
                {
                  outcome.frameQBefore[channelIndex] += frameQ;
                  ++outcome.eventCountBefore[channelIndex];
                  outcome.eventTimeBefore[channelIndex] += static_cast<double>(step) * system.timeStep;
                }
                else
                {
                  // the state just after the event (the net cage count of the tagged molecule has already
                  // changed when its own hop is registered)
                  const double self = static_cast<double>(barrier.cageIndex(unwrappedQ[taggedIndex])) + 0.5;
                  outcome.events.push_back(PostCommitmentEvent{
                      channel, frameQ, static_cast<double>(step - outcome.commitStep) * system.timeStep, step, self,
                      localQ});
                  if (!firstEventStep.has_value())
                  {
                    // the first event after commitment, as the hybrid estimator uses it
                    firstEventStep = step;
                    firstEventTagged = self;
                    firstEventOthers = localQ;
                  }
                  if (++eventsAfterCommitment >= settings.terminationEvents)
                  {
                    terminated = true;
                    outcome.terminationStep = step;
                    outcome.finalSelf = self;
                    outcome.finalLocalQ = localQ;
                  }
                }
              }
            }
            committedCage[j] = newCage;
            cageOffset[j][i] -= static_cast<double>(hop);  // now measured from the center of the new cage
          }
        }
      }

      // hybrid mode: the smoothed axial displacement of the other molecules while they are in the column,
      // evaluated from the cage each molecule is committed to (so that it does not depend on where in the
      // periodic box the molecule started, and stays unique on a ring of four or more cages). A molecule that
      // leaves the column (laterally, or along the normal beyond A' or B') takes its contribution back to the
      // center of the cage it left, and one that enters starts from the center of the cage it enters: like
      // Q_local, only hops between cages of the column count, not the position inside a cage at the entry or
      // the exit.
      if (hybridMode)
      {
        smoothQ[taggedIndex] = smoothCoordinate(unwrappedQ[taggedIndex]);
        for (std::size_t j = 0; j < numberOfMolecules; ++j)
        {
          if (j == taggedIndex) continue;
          const std::optional<int> index = columnIndex(committedCage[j]);
          if (index.has_value())
          {
            const double c = smoothOtherCoordinate(j, index.value());
            smoothOthers += c - (inColumnBefore[j] ? smoothQ[j] : columnCenter(index.value()));
            smoothQ[j] = c;
          }
          else if (inColumnBefore[j])
          {
            smoothOthers += columnCenter(columnIndexBefore[j]) - smoothQ[j];
          }
          inColumnBefore[j] = index.has_value();
          if (index.has_value()) columnIndexBefore[j] = index.value();
        }
      }

      // hybrid mode: record the state at the handover times after commitment
      if (hybridMode && committed)
      {
        while (nextHandover < numberOfHandovers && step == outcome.commitStep + handoverSteps[nextHandover])
        {
          Handover& handover = outcome.handover[nextHandover];
          handover.reached = true;
          if (firstEventStep.has_value())
          {
            handover.eventBefore = true;
            handover.tagged = firstEventTagged;
            handover.others = firstEventOthers;
          }
          else
          {
            handover.tagged = smoothQ[taggedIndex];
            handover.others = smoothOthers;
            handover.state = columnState();
            handover.neighbourhood = neighbourhoodState();
            handover.orientation = frameOrientation();
          }
          ++nextHandover;
        }
      }

      // the stored transmission functions hold their values from the terminating event on
      if (!terminatedBefore && step % sampleEvery == 0uz) storeSample(step / sampleEvery);
      if (terminated && !terminatedBefore)
      {
        for (std::size_t t = step / sampleEvery + 1uz; t < numberOfTimes; ++t) storeSample(t);
      }
      // after the terminating event every remaining handover takes the values of the first event (termination
      // needs at least one event after commitment), so the hybrid estimator needs no further integration
      if (terminated) break;
    }

    const std::size_t lastStep = std::min(step, numberOfTrajectorySteps);
    outcome.integratedSteps = lastStep;
    outcome.eventTerminated = terminated;
    outcome.meanDrift = driftSum / static_cast<double>(lastStep);
    if (!terminated)
    {
      outcome.terminationStep = lastStep;
      outcome.finalSelf = static_cast<double>(barrier.cageIndex(unwrappedQ[taggedIndex])) + 0.5;
      outcome.finalLocalQ = localQ;
      // hold the final values for the remaining stored time points (the index lastStep / sampleEvery has been
      // stored inside the loop when lastStep is a multiple of sampleEvery)
      for (std::size_t t = lastStep / sampleEvery + 1uz; t < numberOfTimes; ++t)
      {
        storeSample(t);
      }
    }
    // hybrid mode: handover times that were not reached (capped or never committed) take the values at the end
    for (std::size_t k = nextHandover; k < numberOfHandovers; ++k)
    {
      Handover& handover = outcome.handover[k];
      handover.reached = false;
      if (firstEventStep.has_value())
      {
        handover.eventBefore = true;
        handover.tagged = firstEventTagged;
        handover.others = firstEventOthers;
      }
      else
      {
        handover.tagged = smoothQ[taggedIndex];
        handover.others = smoothOthers;
      }
    }
  };

  std::vector<double> forwardSelf(numberOfTimes);
  std::vector<double> forwardSelfTwoState(numberOfTimes);
  std::vector<double> forwardSelfCommitted(numberOfTimes);
  std::vector<double> forwardSelfSmooth(numberOfTimes);
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
    if (anyHeld) holdInPlace(heldComponent, moleculeData, dynamics, groupSpan, system.components);
    // the guest center-of-mass velocity is not removed, not even without a framework: the dividing plane is fixed
    // in the laboratory frame, and subtracting the drift would scale the velocity of the tagged bead along the
    // normal by sqrt(1 - 1/N) (zero for a single molecule), and with it <|qdot|>, k^TST and kappa

    // orient the initial velocity of the reaction bead towards B
    const double3 velocity =
        beadVelocity(taggedComponent, moleculeData[taggedIndex], dynamics, groupSpan,
                     system.groupDataOffsetOfMolecule(settings.taggedComponent, 0), taggedComponent.startingBead);
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
    forwardSelfSmooth = deltaSelfSmooth;
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
    std::vector<double>& selfSmooth = accumulator.selfSmooth[block];
    for (std::size_t t = 0; t < numberOfTimes; ++t)
    {
      self[t] += qdot * (forwardSelf[t] - deltaSelf[t]);
      selfTwoState[t] += qdot * (forwardSelfTwoState[t] - deltaSelfTwoState[t]);
      selfCommitted[t] += qdot * (forwardSelfCommitted[t] - deltaSelfCommitted[t]);
      if (smoothMode) selfSmooth[t] += qdot * (forwardSelfSmooth[t] - deltaSelfSmooth[t]);
    }
    for (std::size_t i = 0; i < collective.size(); ++i)
    {
      collective[i] += qdot * (forwardCollective[i] - deltaCollective[i]);
    }

    // occupancy classes: the pair sits on the boundary between (m_A + 1, m_B) (the forward member hops A -> B)
    // and (m_B + 1, m_A) (the time-reversed member hops B -> A); by microscopic reversibility both share the
    // recrossing factor <|qdot| (Delta^+ - Delta^-)> / <|qdot|> of the pair, so both classes take the pair
    // difference. Splitting the members (Delta^+ into one class, -Delta^- into the other) gave 2 p^+ - 1 =
    // kappa + P(B -> B) - P(A -> A): every non-reactive bounce added to one class and subtracted from the other.
    {
      auto plateau = [&](const std::vector<double>& forwardValues, const std::vector<double>& backwardValues)
      {
        double sum{0.0};
        for (std::size_t t : plateauIndices) sum += forwardValues[t] - backwardValues[t];
        return sum / static_cast<double>(std::max(1uz, plateauIndices.size()));
      };
      const double pairSelf = plateau(forwardSelf, deltaSelf);
      const double pairSmooth = smoothMode ? plateau(forwardSelfSmooth, deltaSelfSmooth) : 0.0;
      const double pairCommitted = plateau(forwardSelfCommitted, deltaSelfCommitted);
      auto addClass = [&](std::size_t nFrom, std::size_t nTo)
      {
        DcTSTTransmissionAccumulator::ClassSums& sums = accumulator.classes[{nFrom, nTo}];
        if (sums.speed.empty())
        {
          sums.speed.assign(accumulator.numberOfBlocks, 0.0);
          sums.self.assign(accumulator.numberOfBlocks, 0.0);
          sums.selfSmooth.assign(accumulator.numberOfBlocks, 0.0);
          sums.selfCommitted.assign(accumulator.numberOfBlocks, 0.0);
        }
        sums.speed[block] += qdot;
        sums.self[block] += qdot * pairSelf;
        if (smoothMode) sums.selfSmooth[block] += qdot * pairSmooth;
        sums.selfCommitted[block] += qdot * pairCommitted;
        ++sums.members;
      };
      if (!cageLatticeInconsistency.has_value())
      {
        addClass(othersInA + 1uz, othersInB);
        addClass(othersInB + 1uz, othersInA);
      }
    }

    // event-terminated estimator at every order k = 1..n: the values of the members at their k-th event (at
    // termination when they had fewer events), the channel of that event and the events up to it
    for (std::size_t o = 0; o < accumulator.byEvent.size(); ++o)
    {
      DcTSTTransmissionAccumulator::EventOrderSums& sums = accumulator.byEvent[o];
      const std::size_t order = o + 1uz;
      const TrajectoryOutcome::OrderView forwardView = forward.at(order);
      const TrajectoryOutcome::OrderView backwardView = backward.at(order);
      sums.self[block] += qdot * (forwardView.self - backwardView.self);
      sums.collective[block] +=
          qdot * ((forwardView.self + forwardView.localQ) - (backwardView.self + backwardView.localQ));
      for (std::size_t d = 0; const TrajectoryOutcome* outcome : {&forward, &backward})
      {
        const TrajectoryOutcome::OrderView view = outcome->at(order);
        if (view.ended)
        {
          ++sums.numberOfEventTerminatedTrajectories;
          const double time = static_cast<double>(view.step) * system.timeStep;
          sums.terminationTimeSum += time;
          sums.terminationTimeMaximum = std::max(sums.terminationTimeMaximum, time);
          ++sums.terminatingEvents[d][static_cast<std::size_t>(view.channel.value())];

          // conditioned estimator: this member alone, with the sign of its launch direction
          const double sign = (d == 0uz) ? 1.0 : -1.0;
          sums.speedTerminated[block] += qdot;
          sums.selfConditioned[block] += qdot * sign * view.self;
          sums.collectiveConditioned[block] += qdot * sign * (view.self + view.localQ);
        }
        else if (outcome->committed)
        {
          ++sums.numberOfCappedTrajectories;
        }
        // the events after commitment up to the k-th (frame of the tagged hop: the forward and the time-reversed
        // member add with the same sign)
        double phaseQ{0.0};
        const std::size_t count = std::min(order, outcome->events.size());
        for (std::size_t i = 0; i < count; ++i)
        {
          const PostCommitmentEvent& event = outcome->events[i];
          const std::size_t c = static_cast<std::size_t>(event.channel);
          ++sums.allEvents[d][c];
          phaseQ += event.frameQ;
          sums.localQChannelAfter[c] += qdot * event.frameQ;
          ++sums.channelEventsAfter[c];
          sums.channelTimeSumAfter[c] += event.time;
        }
        sums.collectiveAfterCommitment[block] += qdot * phaseQ;
        ++d;
      }
    }
    // hybrid estimator: pair differences at the handover times and the occupancy histogram
    for (std::size_t h = 0; h < numberOfHandovers; ++h)
    {
      accumulator.hybridTagged[h][block] += qdot * (forward.handover[h].tagged - backward.handover[h].tagged);
      accumulator.hybridOthers[h][block] += qdot * (forward.handover[h].others - backward.handover[h].others);
      for (double launch = 1.0; const TrajectoryOutcome* outcome : {&forward, &backward})
      {
        const Handover& handover = outcome->handover[h];
        // The expected sign of the first event is evaluated in the frame of the commitment, E(state); in the lab
        // frame it is orientation * E, and the pair difference takes the forward member with +1 and the
        // time-reversed one with -1: the sign is +1 for a member that committed the way it was launched and -1
        // for one that turned back first. Both entered with +1, so that the two members of a pair that bounced
        // into the same cage, in the same configuration, added their tails instead of cancelling them.
        const double memberSign = launch * static_cast<double>(handover.orientation);
        launch = -1.0;
        if (handover.eventBefore)
        {
          // also when the handover time lies beyond the end of the trajectory: the event came first
          ++accumulator.hybridEventBefore[h];
        }
        else if (!handover.reached)
        {
          ++accumulator.hybridTruncated[h];
        }
        else
        {
          std::pair<double, std::size_t>& entry = accumulator.hybridOccupancy[h][handover.state];
          entry.first += qdot;
          ++entry.second;
          std::pair<double, std::size_t>& neighbourhood = accumulator.hybridNeighbourhood[h][handover.neighbourhood];
          neighbourhood.first += qdot;
          ++neighbourhood.second;
          accumulator.hybridOccupancySigned[h][handover.state] += memberSign * qdot;
          accumulator.hybridNeighbourhoodSigned[h][handover.neighbourhood] += memberSign * qdot;
        }
      }
    }

    for (const TrajectoryOutcome* outcome : {&forward, &backward})
    {
      // the flux of the members that reached the cap without a single event after commitment (or uncommitted)
      if (eventMode && outcome->events.empty() && !outcome->eventTerminated) accumulator.eventFreeCappedSpeed += qdot;
      accumulator.integratedTimeSum += static_cast<double>(outcome->integratedSteps) * system.timeStep;
      if (!outcome->committed) ++accumulator.numberOfUncommittedTrajectories;
      // the crossings of the other molecules before the commitment of the tagged molecule, per channel (frame of
      // the tagged hop: the forward and the time-reversed member add with the same sign)
      double phaseQ{0.0};
      for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
      {
        phaseQ += outcome->frameQBefore[c];
        accumulator.localQChannelBefore[c] += qdot * outcome->frameQBefore[c];
        accumulator.channelEventsBefore[c] += outcome->eventCountBefore[c];
        accumulator.channelTimeSumBefore[c] += outcome->eventTimeBefore[c];
      }
      accumulator.collectiveEventBeforeCommitment[block] += qdot * phaseQ;
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
            if (c != settings.taggedComponent) continue;  // kappa_C is of the tagged species
            const std::size_t index = accumulator.collectiveIndex(s, c, t);
            collectivePlateau += forwardCollective[index] - deltaCollective[index];
          }
        }
      }
      collectivePlateau /= static_cast<double>(std::max(1uz, plateauIndices.size()));

      auto status = [](const TrajectoryOutcome& outcome) -> int
      { return outcome.eventTerminated ? 0 : (outcome.committed ? 1 : 2); };
      auto channel = [](const TrajectoryOutcome& outcome) -> int
      { return outcome.eventTerminated ? static_cast<int>(outcome.events.back().channel) : -1; };
      auto beforeSum = [](const TrajectoryOutcome& outcome)
      {
        double sum{0.0};
        for (double q : outcome.frameQBefore) sum += q;
        return sum;
      };
      auto afterSum = [](const TrajectoryOutcome& outcome)
      {
        double sum{0.0};
        for (const PostCommitmentEvent& event : outcome.events) sum += event.frameQ;
        return sum;
      };
      auto time = [&](std::size_t steps) { return static_cast<double>(steps) * system.timeStep; };

      std::ostream pairStream(pairStreams[replicaId].rdbuf());
      std::print(pairStream,
                 "{} {} {:.8g} {:g} {:g} {:g} {:g} {} {} {:.6g} {:.6g} {} {} {:.6g} {:.6g} {:g} {:g} {:g} {:g} {} {} "
                 "{:.4e} {:.4e} {:.8g} {:.8g} {:.8g} {:.8g} {:.8g}",
                 productionCycle, block, qdot, forward.finalSelf, backward.finalSelf, forward.finalLocalQ,
                 backward.finalLocalQ, status(forward), status(backward), time(forward.terminationStep),
                 time(backward.terminationStep), channel(forward), channel(backward),
                 forward.committed ? time(forward.commitStep) : 0.0, backward.committed ? time(backward.commitStep) : 0.0,
                 beforeSum(forward), afterSum(forward), beforeSum(backward), afterSum(backward),
                 forward.longJump ? 1 : 0, backward.longJump ? 1 : 0, forward.meanDrift, backward.meanDrift,
                 plateauAverage(forwardSelf, deltaSelf), plateauAverage(forwardSelfTwoState, deltaSelfTwoState),
                 plateauAverage(forwardSelfCommitted, deltaSelfCommitted),
                 forwardSelfCommitted.back() - deltaSelfCommitted.back(), collectivePlateau);
      // the lower event orders k = 1..n-1 (the n-th is in columns 4-13): ten columns per order
      for (std::size_t order = 1; order < accumulator.byEvent.size(); ++order)
      {
        const TrajectoryOutcome::OrderView f = forward.at(order);
        const TrajectoryOutcome::OrderView b = backward.at(order);
        auto orderStatus = [](const TrajectoryOutcome& outcome, const TrajectoryOutcome::OrderView& view) -> int
        { return view.ended ? 0 : (outcome.committed ? 1 : 2); };
        auto orderChannel = [](const TrajectoryOutcome::OrderView& view) -> int
        { return view.channel.has_value() ? static_cast<int>(view.channel.value()) : -1; };
        std::print(pairStream, " {:g} {:g} {:g} {:g} {} {} {:.6g} {:.6g} {} {}", f.self, b.self, f.localQ, b.localQ,
                   orderStatus(forward, f), orderStatus(backward, b), time(f.step), time(b.step), orderChannel(f),
                   orderChannel(b));
      }
      std::print(pairStream, "\n");
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
  const std::size_t numberOfEventOrders = transmission[replicaId].byEvent.size();
  if (numberOfEventOrders > 1uz)
  {
    std::print(file,
               "# columns 29..{}: the lower event orders k = 1..{} (columns 4-13 hold the {}-th), ten columns per order "
               "starting at 29 + 10 (k - 1):\n",
               28uz + 10uz * (numberOfEventOrders - 1uz), numberOfEventOrders - 1uz, numberOfEventOrders);
    std::print(file,
               "#   Delta_S+ and Delta_S- at the k-th event, Q_local+ and Q_local-, status+ and status- (0 reached the "
               "k-th event, 1 capped before it after commitment, 2 never committed),\n"
               "#   time+ and time- of the k-th event [ps], channel+ and channel- of the k-th event (members with fewer "
               "events carry their values at termination and channel -1)\n");
  }
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
  json["selfSmooth"] = accumulator.selfSmooth;
  json["collective"] = accumulator.collective;
  {
    nlohmann::json classes = nlohmann::json::object();
    for (const auto& [pair, sums] : accumulator.classes)
    {
      nlohmann::json entry;
      entry["speed"] = sums.speed;
      entry["self"] = sums.self;
      entry["selfSmooth"] = sums.selfSmooth;
      entry["selfCommitted"] = sums.selfCommitted;
      entry["members"] = sums.members;
      classes[std::format("{},{}", pair.first, pair.second)] = std::move(entry);
    }
    json["classes"] = classes;
  }
  json["sumSpeedSquared"] = accumulator.sumSpeedSquared;
  json["sumRelativeEnergyDrift"] = accumulator.sumRelativeEnergyDrift;
  json["maximumRelativeEnergyDrift"] = accumulator.maximumRelativeEnergyDrift;
  json["numberOfTrajectories"] = accumulator.numberOfTrajectories;
  json["numberOfLongJumpTrajectories"] = accumulator.numberOfLongJumpTrajectories;
  json["numberOfCommittedTrajectories"] = accumulator.numberOfCommittedTrajectories;
  // the event-terminated sums: the n-th order under the top-level keys, every order k = 1..n under 'byEvent'
  auto writeOrder = [](nlohmann::json& target, const DcTSTTransmissionAccumulator::EventOrderSums& sums)
  {
    target["selfEvent"] = sums.self;
    target["collectiveEvent"] = sums.collective;
    target["speedEventTerminated"] = sums.speedTerminated;
    target["selfEventConditioned"] = sums.selfConditioned;
    target["collectiveEventConditioned"] = sums.collectiveConditioned;
    target["collectiveEventAfterCommitment"] = sums.collectiveAfterCommitment;
    target["localQChannelAfter"] = sums.localQChannelAfter;
    target["channelTimeSumAfter"] = sums.channelTimeSumAfter;
    target["channelEventsAfter"] = sums.channelEventsAfter;
    target["terminatingEvents"] = sums.terminatingEvents;
    target["allEvents"] = sums.allEvents;
    target["numberOfEventTerminatedTrajectories"] = sums.numberOfEventTerminatedTrajectories;
    target["numberOfCappedTrajectories"] = sums.numberOfCappedTrajectories;
    target["terminationTimeSum"] = sums.terminationTimeSum;
    target["terminationTimeMaximum"] = sums.terminationTimeMaximum;
  };
  writeOrder(json, accumulator.byEvent.back());
  json["terminationEvents"] = accumulator.byEvent.size();
  json["byEvent"] = nlohmann::json::array();
  for (std::size_t o = 0; o < accumulator.byEvent.size(); ++o)
  {
    nlohmann::json entry;
    entry["order"] = o + 1uz;
    writeOrder(entry, accumulator.byEvent[o]);
    json["byEvent"].push_back(std::move(entry));
  }
  json["collectiveEventBeforeCommitment"] = accumulator.collectiveEventBeforeCommitment;
  json["localQChannelBefore"] = accumulator.localQChannelBefore;
  json["channelTimeSumBefore"] = accumulator.channelTimeSumBefore;
  json["channelEventsBefore"] = accumulator.channelEventsBefore;
  json["numberOfUncommittedTrajectories"] = accumulator.numberOfUncommittedTrajectories;
  json["integratedTimeSum"] = accumulator.integratedTimeSum;
  json["eventFreeCappedSpeed"] = accumulator.eventFreeCappedSpeed;
  if (hybridEnabled())
  {
    nlohmann::json& hybrid = json["hybrid"];
    hybrid["handoverTimes"] = settings.handoverTimes;
    hybrid["tagged"] = accumulator.hybridTagged;
    hybrid["others"] = accumulator.hybridOthers;
    hybrid["eventBefore"] = accumulator.hybridEventBefore;
    hybrid["truncated"] = accumulator.hybridTruncated;
    hybrid["occupancyKey"] = "nA',nA,nB,nB' (frame of the tagged hop, tagged molecule counted in B): [sum |qdot|, members]";
    hybrid["neighbourhoodKey"] =
        "nA',nA,nB,nB'|a1,a2,a3,a4|b1,b2,b3,b4 (column plus sorted lateral occupancies of A and of B)";
    for (std::size_t k = 0; k < accumulator.hybridOccupancy.size(); ++k)
    {
      nlohmann::json entry = nlohmann::json::object();
      for (const auto& [state, weight] : accumulator.hybridOccupancy[k])
      {
        entry[std::format("{},{},{},{}", state.nAPrime, state.nA, state.nB, state.nBPrime)] = {weight.first, weight.second};
      }
      hybrid["occupancy"].push_back(entry);
    }
    for (std::size_t k = 0; k < accumulator.hybridNeighbourhood.size(); ++k)
    {
      nlohmann::json entry = nlohmann::json::object();
      for (const auto& [state, weight] : accumulator.hybridNeighbourhood[k])
      {
        entry[std::format("{},{},{},{}|{},{},{},{}|{},{},{},{}", state.column.nAPrime, state.column.nA, state.column.nB,
                          state.column.nBPrime, state.lateralA[0], state.lateralA[1], state.lateralA[2],
                          state.lateralA[3], state.lateralB[0], state.lateralB[1], state.lateralB[2],
                          state.lateralB[3])] = {weight.first, weight.second};
      }
      hybrid["neighbourhood"].push_back(entry);
    }
    // the weights signed by the commitment of the member, which the Markov tail is summed with (the snapshot held
    // the unsigned ones only, and the hybrid estimate could not be rebuilt from it)
    hybrid["signedKey"] = "the same keys: sum of sign * |qdot|";
    for (std::size_t k = 0; k < accumulator.hybridOccupancySigned.size(); ++k)
    {
      nlohmann::json entry = nlohmann::json::object();
      for (const auto& [state, weight] : accumulator.hybridOccupancySigned[k])
      {
        entry[std::format("{},{},{},{}", state.nAPrime, state.nA, state.nB, state.nBPrime)] = weight;
      }
      hybrid["occupancySigned"].push_back(entry);
    }
    for (std::size_t k = 0; k < accumulator.hybridNeighbourhoodSigned.size(); ++k)
    {
      nlohmann::json entry = nlohmann::json::object();
      for (const auto& [state, weight] : accumulator.hybridNeighbourhoodSigned[k])
      {
        entry[std::format("{},{},{},{}|{},{},{},{}|{},{},{},{}", state.column.nAPrime, state.column.nA, state.column.nB,
                          state.column.nBPrime, state.lateralA[0], state.lateralA[1], state.lateralA[2],
                          state.lateralA[3], state.lateralB[0], state.lateralB[1], state.lateralB[2],
                          state.lateralB[3])] = weight;
      }
      hybrid["neighbourhoodSigned"].push_back(entry);
    }
  }

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
  std::vector<double> selfSmooth(numberOfTimes, 0.0);
  std::vector<double> collective(numberOfShells * numberOfComponents * numberOfTimes, 0.0);
  std::vector<double> blockSpeeds;
  std::vector<std::vector<double>> blockSelf;
  std::vector<std::vector<double>> blockSelfCommitted;
  std::vector<std::vector<double>> blockSelfSmooth;
  std::vector<std::vector<std::vector<double>>> blockCollective;  // [block][shell][time], cumulative, kappa_C
  std::size_t numberOfTrajectories{0};
  std::size_t numberOfLongJumps{0};
  std::size_t numberOfCommitted{0};
  double collectiveEventBefore{0.0};
  std::vector<double> blockCollectiveEventBefore;
  using ChannelArray = std::array<double, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>;
  ChannelArray localQChannelBefore{};
  ChannelArray channelTimeSumBefore{};
  double integratedTimeSum{0.0};
  double eventFreeCappedSpeed{0.0};
  for (const DcTSTTransmissionAccumulator& accumulator : transmission)
  {
    eventFreeCappedSpeed += accumulator.eventFreeCappedSpeed;
    result.numberOfTrajectoryPairs +=
        std::accumulate(accumulator.numberOfTrajectoryPairs.begin(), accumulator.numberOfTrajectoryPairs.end(), 0uz);
    result.meanRelativeEnergyDrift += accumulator.sumRelativeEnergyDrift;
    result.maximumRelativeEnergyDrift =
        std::max(result.maximumRelativeEnergyDrift, accumulator.maximumRelativeEnergyDrift);
    numberOfTrajectories += accumulator.numberOfTrajectories;
    numberOfLongJumps += accumulator.numberOfLongJumpTrajectories;
    numberOfCommitted += accumulator.numberOfCommittedTrajectories;
    result.numberOfUncommittedTrajectories += accumulator.numberOfUncommittedTrajectories;
    integratedTimeSum += accumulator.integratedTimeSum;
    for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
    {
      localQChannelBefore[c] += accumulator.localQChannelBefore[c];
      channelTimeSumBefore[c] += accumulator.channelTimeSumBefore[c];
      result.channelEventsBefore[c] += accumulator.channelEventsBefore[c];
    }
    for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
    {
      speed += accumulator.speed[b];
      if (accumulator.numberOfTrajectoryPairs[b] > 0uz && accumulator.speed[b] > 0.0)
      {
        blockCollectiveEventBefore.push_back(accumulator.collectiveEventBeforeCommitment[b] / accumulator.speed[b]);
      }
      collectiveEventBefore += accumulator.collectiveEventBeforeCommitment[b];
      for (std::size_t t = 0; t < numberOfTimes; ++t)
      {
        self[t] += accumulator.self[b][t];
        selfTwoState[t] += accumulator.selfTwoState[b][t];
        selfCommitted[t] += accumulator.selfCommitted[b][t];
        selfSmooth[t] += accumulator.selfSmooth[b][t];
      }
      for (std::size_t i = 0; i < collective.size(); ++i) collective[i] += accumulator.collective[b][i];

      if (accumulator.numberOfTrajectoryPairs[b] == 0uz || accumulator.speed[b] <= 0.0) continue;
      blockSpeeds.push_back(accumulator.speed[b] / static_cast<double>(accumulator.numberOfTrajectoryPairs[b]));
      std::vector<double> kappaSelf(numberOfTimes);
      std::vector<double> kappaCommitted(numberOfTimes);
      std::vector<double> kappaSmooth(numberOfTimes);
      for (std::size_t t = 0; t < numberOfTimes; ++t)
      {
        kappaSelf[t] = accumulator.self[b][t] / accumulator.speed[b];
        kappaCommitted[t] = accumulator.selfCommitted[b][t] / accumulator.speed[b];
        kappaSmooth[t] = accumulator.selfSmooth[b][t] / accumulator.speed[b];
      }
      blockSelf.push_back(kappaSelf);
      blockSelfCommitted.push_back(kappaCommitted);
      if (!committorTable.empty()) blockSelfSmooth.push_back(kappaSmooth);

      std::vector<std::vector<double>> kappaCollective(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
      std::vector<double> cumulative(numberOfTimes, 0.0);
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        for (std::size_t c = 0; c < numberOfComponents; ++c)
        {
          if (c != settings.taggedComponent) continue;  // kappa_C is of the tagged species
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
  if (numberOfTrajectories > 0uz)
  {
    result.meanIntegratedTime = integratedTimeSum / static_cast<double>(numberOfTrajectories);
  }

  // the event-terminated estimator at every order k = 1..n: pooled sums and block values over the replicas
  // (the parts relative to kappa_S^committed follow once the transmission functions are known)
  const std::size_t numberOfEventOrders = front.byEvent.size();
  result.eventOrders.assign(numberOfEventOrders, TransmissionResult::EventOrder{});
  // per order and block: (kappa_C^event - kappa_S^event), (kappa_S^event - kappa_S^committed), epsilon_after sums
  std::vector<std::vector<double>> blockEventExcess(numberOfEventOrders);
  std::vector<std::vector<double>> blockEventBackCorrelation(numberOfEventOrders);
  std::vector<std::vector<double>> blockCollectiveEventAfter(numberOfEventOrders);
  std::vector<double> collectiveEventAfter(numberOfEventOrders, 0.0);
  std::vector<ChannelArray> localQChannelAfter(numberOfEventOrders);
  std::vector<ChannelArray> channelTimeSumAfter(numberOfEventOrders);
  for (std::size_t o = 0; o < numberOfEventOrders; ++o)
  {
    TransmissionResult::EventOrder& order = result.eventOrders[o];
    order.order = o + 1uz;
    localQChannelAfter[o].fill(0.0);
    channelTimeSumAfter[o].fill(0.0);
    double selfEvent{0.0}, collectiveEvent{0.0};
    double speedEventTerminated{0.0}, selfEventConditioned{0.0}, collectiveEventConditioned{0.0};
    double terminationTimeSum{0.0};
    std::vector<double> blockSelfEvent, blockCollectiveEvent;
    std::vector<double> blockSelfEventConditioned, blockCollectiveEventConditioned;
    for (const DcTSTTransmissionAccumulator& accumulator : transmission)
    {
      if (o >= accumulator.byEvent.size()) continue;
      const DcTSTTransmissionAccumulator::EventOrderSums& sums = accumulator.byEvent[o];
      order.numberOfEventTerminatedTrajectories += sums.numberOfEventTerminatedTrajectories;
      order.numberOfCappedTrajectories += sums.numberOfCappedTrajectories;
      terminationTimeSum += sums.terminationTimeSum;
      order.maximumTerminationTime = std::max(order.maximumTerminationTime, sums.terminationTimeMaximum);
      for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
      {
        for (std::size_t d = 0; d < 2; ++d)
        {
          order.terminatingEvents[d][c] += sums.terminatingEvents[d][c];
          order.allEvents[d][c] += sums.allEvents[d][c];
        }
        localQChannelAfter[o][c] += sums.localQChannelAfter[c];
        channelTimeSumAfter[o][c] += sums.channelTimeSumAfter[c];
        order.channelEventsAfter[c] += sums.channelEventsAfter[c];
      }
      for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
      {
        selfEvent += sums.self[b];
        collectiveEvent += sums.collective[b];
        collectiveEventAfter[o] += sums.collectiveAfterCommitment[b];
        if (accumulator.numberOfTrajectoryPairs[b] > 0uz && accumulator.speed[b] > 0.0)
        {
          blockSelfEvent.push_back(sums.self[b] / accumulator.speed[b]);
          blockCollectiveEvent.push_back(sums.collective[b] / accumulator.speed[b]);
          // the differences are far less noisy than the kappas themselves (same trajectories)
          blockEventExcess[o].push_back((sums.collective[b] - sums.self[b]) / accumulator.speed[b]);
          blockEventBackCorrelation[o].push_back((sums.self[b] - accumulator.selfCommitted[b].back()) /
                                                 accumulator.speed[b]);
          blockCollectiveEventAfter[o].push_back(sums.collectiveAfterCommitment[b] / accumulator.speed[b]);
        }
        speedEventTerminated += sums.speedTerminated[b];
        selfEventConditioned += sums.selfConditioned[b];
        collectiveEventConditioned += sums.collectiveConditioned[b];
        if (sums.speedTerminated[b] > 0.0)
        {
          // per trajectory member: kappa = 2 sum |qdot| s Delta / sum |qdot| (see the accumulator)
          blockSelfEventConditioned.push_back(2.0 * sums.selfConditioned[b] / sums.speedTerminated[b]);
          blockCollectiveEventConditioned.push_back(2.0 * sums.collectiveConditioned[b] / sums.speedTerminated[b]);
        }
      }
    }
    if (order.numberOfEventTerminatedTrajectories > 0uz)
    {
      order.meanTerminationTime = terminationTimeSum / static_cast<double>(order.numberOfEventTerminatedTrajectories);
    }
    if (speed > 0.0)
    {
      order.selfEvent = selfEvent / speed;
      order.collectiveEvent = collectiveEvent / speed;
      order.selfEventError = meanAndError(blockSelfEvent).second;
      order.collectiveEventError = meanAndError(blockCollectiveEvent).second;
    }
    if (speedEventTerminated > 0.0)
    {
      order.selfEventConditioned = 2.0 * selfEventConditioned / speedEventTerminated;
      order.collectiveEventConditioned = 2.0 * collectiveEventConditioned / speedEventTerminated;
      order.selfEventConditionedError = meanAndError(blockSelfEventConditioned).second;
      order.collectiveEventConditionedError = meanAndError(blockCollectiveEventConditioned).second;
    }
    if (numberOfTrajectories > 0uz)
    {
      order.cappedFraction =
          static_cast<double>(order.numberOfCappedTrajectories + result.numberOfUncommittedTrajectories) /
          static_cast<double>(numberOfTrajectories);
    }
  }

  result.self.assign(numberOfTimes, 0.0);
  result.selfTwoState.assign(numberOfTimes, 0.0);
  result.selfCommitted.assign(numberOfTimes, 0.0);
  result.selfSmooth.assign(numberOfTimes, 0.0);
  result.selfSmoothError.assign(numberOfTimes, 0.0);
  result.selfError.assign(numberOfTimes, 0.0);
  result.smoothAvailable = !committorTable.empty();
  result.collective.assign(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
  result.collectiveError.assign(numberOfShells, std::vector<double>(numberOfTimes, 0.0));
  result.collectivePerComponent.assign(
      numberOfShells, std::vector<std::vector<double>>(numberOfComponents, std::vector<double>(numberOfTimes, 0.0)));
  result.collectivePlateau.assign(numberOfShells, 0.0);
  result.collectivePlateauError.assign(numberOfShells, 0.0);
  if (speed <= 0.0) return result;
  // every pair has two members of the same |qdot|
  result.eventFreeCappedFraction = eventFreeCappedSpeed / (2.0 * speed);

  for (std::size_t t = 0; t < numberOfTimes; ++t)
  {
    result.self[t] = self[t] / speed;
    result.selfTwoState[t] = selfTwoState[t] / speed;
    result.selfCommitted[t] = selfCommitted[t] / speed;
    std::vector<double> values;
    for (const std::vector<double>& k : blockSelf) values.push_back(k[t]);
    result.selfError[t] = meanAndError(values).second;
    if (result.smoothAvailable)
    {
      result.selfSmooth[t] = selfSmooth[t] / speed;
      values.clear();
      for (const std::vector<double>& k : blockSelfSmooth) values.push_back(k[t]);
      result.selfSmoothError[t] = meanAndError(values).second;
    }
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
      // kappa_C of the tagged species: its own molecules only. D_C = kappa_C k^TST lambda^2 is then
      // (k_B T / N_0) L_00, the collective diffusivity of that species. The sum over every component, as was
      // taken, is (k_B T / N_0) sum_c L_0c: neither that nor the collective diffusivity of the mixture, and with
      // the cage spacing of the tagged species applied to all. The cross terms stay in the per-component columns.
      result.collective[s][t] = result.self[t] + result.collectivePerComponent[s][settings.taggedComponent][t];
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
  if (result.smoothAvailable)
  {
    result.selfSmoothPlateau = windowAverage(result.selfSmooth);
    std::vector<double> values;
    for (const std::vector<double>& k : blockSelfSmooth) values.push_back(windowAverage(k));
    result.selfSmoothPlateauError = meanAndError(values).second;
  }
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

  // occupancy classes of the trajectory members (post-stratification) and the rates k(n -> n')
  {
    struct PooledClass
    {
      double speed{0.0};
      double self{0.0};
      double selfSmooth{0.0};
      double selfCommitted{0.0};
      std::size_t members{0};
      std::vector<double> blockKappa;
      std::vector<double> blockKappaSmooth;
    };
    std::map<std::pair<std::size_t, std::size_t>, PooledClass> pooled;
    for (const DcTSTTransmissionAccumulator& accumulator : transmission)
    {
      for (const auto& [pair, sums] : accumulator.classes)
      {
        PooledClass& entry = pooled[pair];
        entry.members += sums.members;
        for (std::size_t b = 0; b < sums.speed.size(); ++b)
        {
          entry.speed += sums.speed[b];
          entry.self += sums.self[b];
          entry.selfSmooth += sums.selfSmooth[b];
          entry.selfCommitted += sums.selfCommitted[b];
          if (sums.speed[b] > 0.0)
          {
            entry.blockKappa.push_back(sums.self[b] / sums.speed[b]);
            if (result.smoothAvailable) entry.blockKappaSmooth.push_back(sums.selfSmooth[b] / sums.speed[b]);
          }
        }
      }
    }
    const double memberSpeed = 2.0 * speed;
    for (const auto& [pair, entry] : pooled)
    {
      TransmissionResult::OccupancyClass item;
      item.nFrom = pair.first;
      item.nTo = pair.second;
      item.members = entry.members;
      item.weightFraction = (memberSpeed > 0.0) ? entry.speed / memberSpeed : 0.0;
      item.kappa = (entry.speed > 0.0) ? entry.self / entry.speed : 0.0;
      item.kappaError = meanAndError(entry.blockKappa).second;
      item.kappaSmooth = (entry.speed > 0.0) ? entry.selfSmooth / entry.speed : 0.0;
      item.kappaSmoothError = meanAndError(entry.blockKappaSmooth).second;
      item.kappaCommitted = (entry.speed > 0.0) ? entry.selfCommitted / entry.speed : 0.0;
      if (auto it = occupancyWeights.find(pair); it != occupancyWeights.end()) item.pilotWeight = it->second;
      if (auto itB = barrierOccupancyWeights.find(pair);
          itB != barrierOccupancyWeights.end() && item.pilotWeight > 0.0)
      {
        item.tstRatio = itB->second / item.pilotWeight;
      }
      result.classes.push_back(item);
    }
    std::sort(result.classes.begin(), result.classes.end(),
              [](const TransmissionResult::OccupancyClass& a, const TransmissionResult::OccupancyClass& b)
              { return a.weightFraction > b.weightFraction; });

    OccupancyRates rates;
    const bool haveKTST = barrierProbability.has_value() && result.meanSpeed > 0.0;
    const double kTST = haveKTST ? 0.5 * result.meanSpeed * barrierProbability.value() : 0.0;
    for (const TransmissionResult::OccupancyClass& item : result.classes)
    {
      const double kappaUsed = result.smoothAvailable ? item.kappaSmooth : item.kappa;
      if (item.pilotWeight > 0.0 && item.tstRatio > 0.0)
      {
        result.classesReweightedKappa += item.pilotWeight * kappaUsed * item.tstRatio;
        result.classesReweightedCoverage += item.pilotWeight;
      }
      if (item.members < minimumClassMembers || item.tstRatio <= 0.0 || !haveKTST) continue;
      OccupancyRates::Entry entry;
      entry.kappa = kappaUsed;
      entry.tstRatio = item.tstRatio;
      entry.rate = kappaUsed * item.tstRatio * kTST;
      entry.weight = static_cast<double>(item.members);
      entry.source = "classes";
      if (entry.rate > 0.0) rates.measured[{item.nFrom, item.nTo}] = entry;
    }
    for (const auto& [pair, entry] : externalRates) rates.measured[pair] = entry;
    fitOccupancyRates(rates);
    if (!rates.measured.empty() || rates.fitValid) result.rates = std::move(rates);
  }

  // first-order back-correlation factor and collective excess of the event-terminated estimator, relative to
  // the committed kappa_S at the end of the trajectories (its values are held after termination)
  const double selfCommittedFinal = result.selfCommitted.back();
  if (std::abs(selfCommittedFinal) > 1e-12)
  {
    for (std::size_t o = 0; o < numberOfEventOrders; ++o)
    {
      TransmissionResult::EventOrder& order = result.eventOrders[o];
      order.backCorrelationFactor = order.selfEvent / selfCommittedFinal;
      order.collectiveExcess = order.collectiveEvent / selfCommittedFinal - 1.0;
      // block errors of the differences; epsilon = (kappa_C^event - kappa_S^event)/kappa_S^committed + (f - 1)
      const double backCorrelationError = meanAndError(blockEventBackCorrelation[o]).second;
      order.backCorrelationFactorError = backCorrelationError / std::abs(selfCommittedFinal);
      // the two differences are strongly anti-correlated (the same kappa_S^event enters both with opposite signs):
      // the error of their sum, kappa_C^event - kappa_S^committed, from the block values of the sum itself
      std::vector<double> blockSum(blockEventExcess[o].size());
      for (std::size_t b = 0; b < blockSum.size(); ++b)
      {
        blockSum[b] = blockEventExcess[o][b] + blockEventBackCorrelation[o][b];
      }
      order.collectiveExcessError = meanAndError(blockSum).second / std::abs(selfCommittedFinal);
      order.backCorrelationFactorConditioned = order.selfEventConditioned / selfCommittedFinal;
      order.collectiveExcessConditioned = order.collectiveEventConditioned / selfCommittedFinal - 1.0;
      // A member that reached the cap without any event contributes Delta = +-1 and no crossing after its
      // commitment, i.e. exactly kappa_S^committed: the deviation f - 1 scales with the flux fraction 1 - p_0 of
      // the members that had an event. The divisor was the fraction that reached the k-th event, counted per
      // trajectory: for k >= 2 a member capped after 1..k-1 events was taken as having contributed nothing while
      // its events (the first one or two carry nearly all of f - 1 and epsilon) were in the sums, and the
      // correction was inflated by 1 / p(k events), to a negative f when many were capped. (epsilon below.)
      const double eventFraction = 1.0 - result.eventFreeCappedFraction;
      if (eventFraction > 1e-12)
      {
        order.backCorrelationFactorCapCorrected = 1.0 + (order.backCorrelationFactor - 1.0) / eventFraction;
      }
      // epsilon_after of the decomposition epsilon - (f - 1) = epsilon_before + epsilon_after, and per channel
      if (speed > 0.0)
      {
        order.collectiveExcessAfterCommitment = collectiveEventAfter[o] / speed / selfCommittedFinal;
        order.collectiveExcessAfterCommitmentError =
            meanAndError(blockCollectiveEventAfter[o]).second / std::abs(selfCommittedFinal);
        for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
        {
          order.collectiveExcessChannelAfter[c] = localQChannelAfter[o][c] / speed / selfCommittedFinal;
          order.channelMeanTimeAfter[c] =
              (order.channelEventsAfter[c] > 0uz)
                  ? channelTimeSumAfter[o][c] / static_cast<double>(order.channelEventsAfter[c])
                  : 0.0;
        }
      }
    }
    // epsilon_before: the crossings before commitment (independent of the event order), and per channel
    if (speed > 0.0)
    {
      result.collectiveExcessBeforeCommitment = collectiveEventBefore / speed / selfCommittedFinal;
      result.collectiveExcessBeforeCommitmentError =
          meanAndError(blockCollectiveEventBefore).second / std::abs(selfCommittedFinal);
      for (std::size_t c = 0; c < static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count); ++c)
      {
        result.collectiveExcessChannelBefore[c] = localQChannelBefore[c] / speed / selfCommittedFinal;
        result.channelMeanTimeBefore[c] =
            (result.channelEventsBefore[c] > 0uz)
                ? channelTimeSumBefore[c] / static_cast<double>(result.channelEventsBefore[c])
                : 0.0;
      }
    }

    // epsilon_corr: the crossings before the commitment are those of every member, capped or not, and are not
    // rescaled
    {
      const double eventFraction = 1.0 - result.eventFreeCappedFraction;
      for (TransmissionResult::EventOrder& order : result.eventOrders)
      {
        if (eventFraction > 1e-12)
        {
          order.collectiveExcessCapCorrected =
              result.collectiveExcessBeforeCommitment +
              (order.collectiveExcess - result.collectiveExcessBeforeCommitment) / eventFraction;
        }
      }
    }

    // hybrid estimator at the handover times
    if (hybridEnabled() && speed > 0.0)
    {
      using ColumnState = DcTSTTransmissionAccumulator::ColumnState;
      using NeighbourhoodState = DcTSTTransmissionAccumulator::NeighbourhoodState;
      const std::size_t numberOfHandovers = handoverSteps.size();
      const std::size_t numberOfMembers = 2uz * result.numberOfTrajectoryPairs;
      // (of the tagged component, as the occupancies of the states are: the cations and co-adsorbates were counted
      // in, and the lateral cages of the column states held them as molecules that hop)
      const double meanOccupancy = static_cast<double>(system.numberOfMoleculesPerComponent[settings.taggedComponent]) /
                                   static_cast<double>(std::max(1uz, cageLattice.numberOfCages()));

      // expected sign of the first axial event from the full neighbourhood (first-order Markov closure):
      // the six axial channels with s = +-1, and (variant 'lateral') the refill of A from each of its four
      // lateral neighbours n_L k(n_L -> n_A) and the escape from B to each lateral n_B k(n_B -> n_L'), s = 0
      auto expectedSign = [&](const NeighbourhoodState& state, bool tableFirst, bool lateral, bool& usedFit) -> double
      {
        if (!result.rates.has_value()) return 0.0;
        const OccupancyRates& rates = result.rates.value();
        auto channelRate = [&](std::size_t nFrom, std::size_t nTo) -> double
        {
          if (nFrom == 0uz) return 0.0;
          if (tableFirst && !rates.isMeasured(nFrom, nTo)) usedFit = true;
          return static_cast<double>(nFrom) * rates.rate(nFrom, nTo, tableFirst);
        };
        const ColumnState& column = state.column;
        const double forward =
            channelRate(column.nAPrime, column.nA) + channelRate(column.nA, column.nB) + channelRate(column.nB, column.nBPrime);
        const double backward =
            channelRate(column.nA, column.nAPrime) + channelRate(column.nB, column.nA) + channelRate(column.nBPrime, column.nB);
        double total = forward + backward;
        if (lateral)
        {
          // every lateral hop that changes n_A or n_B ends the state: into A and out of A, into B and out of B,
          // through the four lateral windows of each (sixteen channels). Only the refill of A and the escape
          // from B were counted, which left the total too small and the expected sign too large.
          for (std::uint16_t nL : state.lateralA)
          {
            total += channelRate(nL, column.nA) + channelRate(column.nA, nL);
          }
          for (std::uint16_t nL : state.lateralB)
          {
            total += channelRate(column.nB, nL) + channelRate(nL, column.nB);
          }
        }
        return (total > 0.0) ? (forward - backward) / total : 0.0;
      };
      auto wrapColumn = [&](const ColumnState& column) -> NeighbourhoodState
      {
        NeighbourhoodState state;
        state.column = column;
        const auto mean = static_cast<std::uint16_t>(std::max(0L, std::lround(meanOccupancy)));
        state.lateralA.fill(mean);
        state.lateralB.fill(mean);
        return state;
      };

      for (std::size_t k = 0; k < numberOfHandovers; ++k)
      {
        TransmissionResult::Hybrid hybrid;
        hybrid.time = static_cast<double>(handoverSteps[k]) * system.timeStep;
        double taggedSum{0.0}, othersSum{0.0};
        std::vector<double> blockExcess;
        std::size_t eventBefore{0}, truncated{0};
        std::map<ColumnState, std::pair<double, std::size_t>> occupancy;
        std::map<NeighbourhoodState, std::pair<double, std::size_t>> neighbourhood;
        // the weights of the Markov tail: signed by the commitment of the member (see hybridOccupancySigned)
        std::map<ColumnState, double> occupancySigned;
        std::map<NeighbourhoodState, double> neighbourhoodSigned;
        for (const DcTSTTransmissionAccumulator& accumulator : transmission)
        {
          if (k >= accumulator.hybridTagged.size()) continue;
          if (k < accumulator.hybridOccupancySigned.size())
          {
            for (const auto& [state, weight] : accumulator.hybridOccupancySigned[k]) occupancySigned[state] += weight;
          }
          if (k < accumulator.hybridNeighbourhoodSigned.size())
          {
            for (const auto& [state, weight] : accumulator.hybridNeighbourhoodSigned[k])
            {
              neighbourhoodSigned[state] += weight;
            }
          }
          for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
          {
            taggedSum += accumulator.hybridTagged[k][b];
            othersSum += accumulator.hybridOthers[k][b];
            if (accumulator.numberOfTrajectoryPairs[b] > 0uz && accumulator.speed[b] > 0.0)
            {
              blockExcess.push_back((accumulator.hybridTagged[k][b] + accumulator.hybridOthers[k][b] -
                                     accumulator.selfCommitted[b].back()) /
                                    accumulator.speed[b]);
            }
          }
          eventBefore += accumulator.hybridEventBefore[k];
          truncated += accumulator.hybridTruncated[k];
          for (const auto& [state, weight] : accumulator.hybridOccupancy[k])
          {
            std::pair<double, std::size_t>& entry = occupancy[state];
            entry.first += weight.first;
            entry.second += weight.second;
          }
          if (k < accumulator.hybridNeighbourhood.size())
          {
            for (const auto& [state, weight] : accumulator.hybridNeighbourhood[k])
            {
              std::pair<double, std::size_t>& entry = neighbourhood[state];
              entry.first += weight.first;
              entry.second += weight.second;
            }
          }
        }
        hybrid.concertedTagged = taggedSum / speed / selfCommittedFinal;
        hybrid.concertedOthers = othersSum / speed / selfCommittedFinal;
        hybrid.concerted = (taggedSum + othersSum) / speed / selfCommittedFinal - 1.0;
        hybrid.concertedError = meanAndError(blockExcess).second / std::abs(selfCommittedFinal);
        hybrid.fractionEventBefore = static_cast<double>(eventBefore) / static_cast<double>(std::max(1uz, numberOfMembers));
        hybrid.fractionTruncated = static_cast<double>(truncated) / static_cast<double>(std::max(1uz, numberOfMembers));
        hybrid.numberOfStates = occupancy.size();
        // the histogram, in descending weight (the weights are fractions of the total |qdot| of all members, so
        // that they sum to the fraction of members that had no event before tau)
        for (const auto& [state, weight] : occupancy)
        {
          hybrid.occupancy.emplace_back(state, weight.first / (2.0 * speed));
          hybrid.occupancyCounts.push_back(weight.second);
        }
        {
          std::vector<std::size_t> order(hybrid.occupancy.size());
          std::iota(order.begin(), order.end(), 0uz);
          std::sort(order.begin(), order.end(),
                    [&](std::size_t a, std::size_t b) { return hybrid.occupancy[a].second > hybrid.occupancy[b].second; });
          std::vector<std::pair<ColumnState, double>> sortedOccupancy;
          std::vector<std::size_t> sortedCounts;
          for (std::size_t i : order)
          {
            sortedOccupancy.push_back(hybrid.occupancy[i]);
            sortedCounts.push_back(hybrid.occupancyCounts[i]);
          }
          hybrid.occupancy = std::move(sortedOccupancy);
          hybrid.occupancyCounts = std::move(sortedCounts);
        }
        hybrid.numberOfNeighbourhoodStates = neighbourhood.size();
        if (result.rates.has_value())
        {
          hybrid.markovAvailable = true;
          double table{0.0}, tableAxial{0.0}, fit{0.0}, fitAxial{0.0}, unmeasured{0.0};
          auto accumulate = [&](const NeighbourhoodState& state, double weight)
          {
            bool usedFit{false};
            bool dummy{false};
            table += weight * expectedSign(state, true, true, usedFit);
            tableAxial += weight * expectedSign(state, true, false, dummy);
            fit += weight * expectedSign(state, false, true, dummy);
            fitAxial += weight * expectedSign(state, false, false, dummy);
            if (usedFit) unmeasured += weight;
          };
          if (!neighbourhood.empty())
          {
            for (const auto& [state, weight] : neighbourhoodSigned) accumulate(state, weight);
          }
          else
          {
            for (const auto& [state, weight] : occupancySigned) accumulate(wrapColumn(state), weight);
          }
          // per member: sum_members |qdot| E[s] / sum_pairs |qdot|, relative to kappa_S^committed
          hybrid.markovTable = table / speed / selfCommittedFinal;
          hybrid.markovTableAxial = tableAxial / speed / selfCommittedFinal;
          hybrid.markovFit = fit / speed / selfCommittedFinal;
          hybrid.markovFitAxial = fitAxial / speed / selfCommittedFinal;
          hybrid.unmeasuredWeight = unmeasured / (2.0 * speed);
          hybrid.hybrid = hybrid.concerted + hybrid.markovTable;
          hybrid.hybridFit = hybrid.concerted + hybrid.markovFit;
        }
        result.hybrid.push_back(std::move(hybrid));
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
    std::print(file, "# column 6: kappa_S(t) committor-smoothed (Ruiz-Montero; 0 if no free-energy profile)\n");
    std::print(file, "# column 7: error of the committor-smoothed kappa_S\n");
    std::size_t column = 8;
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
      std::print(file, "{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}", kappa.time[t], kappa.self[t],
                 kappa.selfError[t], kappa.selfTwoState[t], kappa.selfCommitted[t],
                 kappa.smoothAvailable ? kappa.selfSmooth[t] : 0.0,
                 kappa.smoothAvailable ? kappa.selfSmoothError[t] : 0.0);
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
    if (kappa.smoothAvailable)
    {
      std::print(stream, "    kappa_S smooth (plateau):            {:.5f} +/- {:.5f}\n", kappa.selfSmoothPlateau,
                 kappa.selfSmoothPlateauError);
    }
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

    if (!kappa.classes.empty())
    {
      std::print(stream, "Occupancy classes of the trajectory members (post-stratification; n_from includes the tagged "
                         "molecule, n_to is the arrival cage at launch):\n");
      std::print(stream, "    {:>6} {:>6} {:>8} {:>8} {:>10} {:>9} {:>10} {:>9} {:>8} {:>8}\n", "n_from", "n_to",
                 "members", "weight", "kappa", "+/-", "kappa_sm", "+/-", "P(n,n')", "kTST rel");
      for (const TransmissionResult::OccupancyClass& item : kappa.classes)
      {
        if (item.weightFraction < 0.005 && item.members < minimumClassMembers) continue;
        std::print(stream, "    {:>6} {:>6} {:>8} {:>8.4f} {:>10.5f} {:>9.5f} {:>10.5f} {:>9.5f} {:>8.4f} {:>8.4f}\n",
                   item.nFrom, item.nTo, item.members, item.weightFraction, item.kappa, item.kappaError,
                   item.kappaSmooth, item.kappaSmoothError, item.pilotWeight, item.tstRatio);
      }
      if (kappa.classesReweightedCoverage > 0.0)
      {
        std::print(stream,
                   "    reweighted kappa = sum P(n,n') kappa(n,n') [P(n,n'|q*)/P(n,n')] over the classes: {:.5f} "
                   "(coverage {:.3f} of P(n,n'))\n",
                   kappa.classesReweightedKappa, kappa.classesReweightedCoverage);
        std::print(stream,
                   "    (equals kappa_S of this run when the classes cover the occupancy distribution and use the "
                   "same kappa; the consistency of the rates is sum P k vs kappa k^TST)\n");
      }
      std::print(stream, "\n");
    }
    if (kappa.rates.has_value())
    {
      const OccupancyRates& rates = kappa.rates.value();
      std::print(stream, "Occupancy-resolved rates k(n -> n') [1/ps] (post-stratified classes, overridden by "
                         "'OccupancyRateInputs'):\n");
      std::print(stream, "    {:>6} {:>6} {:>14} {:>8} {:>10} {:>10} {}\n", "n", "n'", "k", "weight", "kappa",
                 "kTST rel", "source");
      for (const auto& [pair, entry] : rates.measured)
      {
        std::print(stream, "    {:>6} {:>6} {:>14.6e} {:>8.0f} {:>10.5f} {:>10.4f} {}\n", pair.first, pair.second,
                   entry.rate, entry.weight, entry.kappa, entry.tstRatio, entry.source);
      }
      if (rates.fitValid)
      {
        std::print(stream, "    fit log k = c + a n - b n': c = {:.4f}, a = {:.4f}, b = {:.4f} (weighted rms {:.4f})\n",
                   rates.fitConstant, rates.fitDeparture, rates.fitArrival, rates.fitRms);
      }
      else if (rates.measured.size() < 3uz)
      {
        std::print(stream, "    (fewer than three measured pairs: the fit is not available)\n");
      }
      if (barrierProbability.has_value() && kappa.classesReweightedCoverage > 0.0)
      {
        const double kTST = 0.5 * kappa.meanSpeed * barrierProbability.value();
        const double kappaUsed = kappa.smoothAvailable ? kappa.selfSmoothPlateau : kappa.selfPlateau;
        std::print(stream,
                   "    consistency: sum P(n,n') k(n->n') = {:.6e} [1/ps],  kappa k^TST = {:.6e} [1/ps]\n",
                   kappa.classesReweightedKappa * kTST, kappaUsed * kTST);
      }
      std::print(stream, "\n");
    }

    if (eventMode)
    {
      using EventChannel = DcTSTTransmissionAccumulator::EventChannel;
      using EventOrder = TransmissionResult::EventOrder;
      const EventOrder& last = kappa.lastEventOrder();
      const std::size_t total = 2uz * kappa.numberOfTrajectoryPairs;
      auto fraction = [&](std::size_t n) { return static_cast<double>(n) / static_cast<double>(std::max(1uz, total)); };
      std::print(stream, "Event-terminated estimator ({} local axial event(s) after commitment, cap {} [ps]):\n",
                 settings.terminationEvents, settings.trajectoryTime);
      std::print(stream, "    trajectories ended by an event:      {} ({:.4f})\n",
                 last.numberOfEventTerminatedTrajectories, fraction(last.numberOfEventTerminatedTrajectories));
      std::print(stream, "    trajectories capped after commitment: {} ({:.4f})\n", last.numberOfCappedTrajectories,
                 fraction(last.numberOfCappedTrajectories));
      std::print(stream, "    trajectories never committed:        {} ({:.4f})\n",
                 kappa.numberOfUncommittedTrajectories, fraction(kappa.numberOfUncommittedTrajectories));
      std::print(stream, "    termination time (event-ended):      mean {:.4f}, maximum {:.4f} [ps]\n",
                 last.meanTerminationTime, last.maximumTerminationTime);
      std::print(stream, "    integrated time per trajectory:      mean {:.4f} [ps]\n", kappa.meanIntegratedTime);
      std::print(stream, "    kappa_S^event (net cage count):      {:.5f} +/- {:.5f}\n", last.selfEvent,
                 last.selfEventError);
      std::print(stream, "    kappa_S^committed (final):           {:.5f}\n", kappa.selfCommitted.back());
      std::print(stream, "    kappa_C^event (self + Q_local):      {:.5f} +/- {:.5f}\n", last.collectiveEvent,
                 last.collectiveEventError);
      std::print(stream, "    f = kappa_S^event / kappa_S^committed:        {:.5f} +/- {:.5f}\n",
                 last.backCorrelationFactor, last.backCorrelationFactorError);
      std::print(stream, "    epsilon = kappa_C^event / kappa_S^committed - 1: {:.5f} +/- {:.5f}\n",
                 last.collectiveExcess, last.collectiveExcessError);
      std::print(stream, "    (errors from the block values of the differences kappa_C^event - kappa_S^event and\n"
                         "     kappa_S^event - kappa_S^committed, which share the trajectories and are much less noisy\n"
                         "     than the kappas themselves; the error of k^TST is not included)\n");
      std::print(stream,
                 "  Cap-corrected (deviations divided by 1 - p_0, p_0 = {:.4f} the flux fraction of the members that\n"
                 "  reached the cap without any event; {:.4f} of the trajectories did not reach the last event):\n",
                 kappa.eventFreeCappedFraction, last.cappedFraction);
      std::print(stream, "    f_corr = 1 + (f - 1) / (1 - p_0):                {:.5f}\n",
                 last.backCorrelationFactorCapCorrected);
      std::print(stream, "    epsilon_corr = eps_before + (eps - eps_before) / (1 - p_0): {:.5f}\n",
                 last.collectiveExcessCapCorrected);
      std::print(stream, "  Conditioned on the event-ended trajectories only (per trajectory, not per pair):\n");
      std::print(stream, "    kappa_S^event|cond:                  {:.5f} +/- {:.5f}\n", last.selfEventConditioned,
                 last.selfEventConditionedError);
      std::print(stream, "    kappa_C^event|cond:                  {:.5f} +/- {:.5f}\n",
                 last.collectiveEventConditioned, last.collectiveEventConditionedError);
      std::print(stream, "    f_cond = kappa_S^event|cond / kappa_S^committed:        {:.5f}\n",
                 last.backCorrelationFactorConditioned);
      std::print(stream, "    epsilon_cond = kappa_C^event|cond / kappa_S^committed - 1: {:.5f}\n",
                 last.collectiveExcessConditioned);
      std::print(stream,
                 "    (f_corr/epsilon_corr and f_cond/epsilon_cond agree when the capped trajectories are an "
                 "unbiased sample;\n"
                 "    a cap T >~ 5 tau_event keeps p_cap below a few percent and makes all three estimates "
                 "coincide)\n");

      if (kappa.eventOrders.size() > 1uz)
      {
        // the lower orders from the same trajectories: a trajectory passes through its events 1..n-1 on the way to
        // the n-th, and its state at the k-th event does not depend on where the integration stops
        std::print(stream,
                   "  Every event order k = 1..{} from the same trajectories (the state at the k-th event after "
                   "commitment,\n"
                   "  identical to a run with 'TerminationEvents' = k and the same cap; f and epsilon should be "
                   "independent of k\n"
                   "  once the hopping beyond the k-th event is uncorrelated with the tagged hop):\n",
                   kappa.eventOrders.size());
        std::print(stream, "    {:>3} {:>9} {:>8} {:>9} {:>8} {:>9} {:>8} {:>9} {:>8} {:>9} {:>9} {:>7} {:>9} {:>9}\n",
                   "k", "kappa_S", "+/-", "kappa_C", "+/-", "f", "+/-", "epsilon", "+/-", "f_cond", "eps_cond", "p_cap",
                   "eps_after", "<t_k>[ps]");
        for (const EventOrder& order : kappa.eventOrders)
        {
          std::print(
              stream,
              "    {:>3} {:>9.5f} {:>8.5f} {:>9.5f} {:>8.5f} {:>9.5f} {:>8.5f} {:>9.5f} {:>8.5f} {:>9.5f} {:>9.5f} "
              "{:>7.4f} {:>9.5f} {:>9.4f}\n",
              order.order, order.selfEvent, order.selfEventError, order.collectiveEvent, order.collectiveEventError,
              order.backCorrelationFactor, order.backCorrelationFactorError, order.collectiveExcess,
              order.collectiveExcessError, order.backCorrelationFactorConditioned, order.collectiveExcessConditioned,
              order.cappedFraction, order.collectiveExcessAfterCommitment, order.meanTerminationTime);
        }
        std::print(stream,
                   "    (p_cap at order k counts the trajectories capped before their k-th event; eps_after is the "
                   "part of epsilon\n"
                   "     from the other molecules after commitment up to the k-th event, eps_before = {:.5f} is "
                   "common to all k;\n"
                   "     the channel tables below are for k = {}, those of every k are in the JSON output "
                   "'results.eventTerminated.byEvent')\n",
                   kappa.collectiveExcessBeforeCommitment, kappa.eventOrders.size());
      }

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
                   last.terminatingEvents[0][c], last.terminatingEvents[1][c],
                   last.terminatingEvents[0][c] + last.terminatingEvents[1][c], last.allEvents[0][c],
                   last.allEvents[1][c], last.allEvents[0][c] + last.allEvents[1][c]);
        terminatingSum += last.terminatingEvents[0][c] + last.terminatingEvents[1][c];
      }
      std::print(stream, "    {:<26} {:>32}\n", "terminating events total", terminatingSum);
      std::print(stream,
                 "    (Q_local counts the axial crossings of the other molecules in the column from t = 0; events\n"
                 "     before the commitment of the tagged molecule do not count towards termination)\n\n");

      // decomposition of epsilon: tagged part (f - 1) and the Q_local crossings before/after commitment, per
      // channel
      std::print(stream, "  Decomposition of the collective excess, epsilon = (f - 1) + epsilon_before + epsilon_after:\n");
      std::print(stream, "    (f - 1), tagged molecule after commitment:         {:>9.5f} +/- {:.5f}\n",
                 last.backCorrelationFactor - 1.0, last.backCorrelationFactorError);
      std::print(stream, "    epsilon_before, others before tagged commitment:   {:>9.5f} +/- {:.5f}\n",
                 kappa.collectiveExcessBeforeCommitment, kappa.collectiveExcessBeforeCommitmentError);
      std::print(stream, "    epsilon_after, others after tagged commitment:     {:>9.5f} +/- {:.5f}\n",
                 last.collectiveExcessAfterCommitment, last.collectiveExcessAfterCommitmentError);
      std::print(stream, "    sum:                                               {:>9.5f}   (epsilon = {:.5f})\n",
                 (last.backCorrelationFactor - 1.0) + kappa.collectiveExcessBeforeCommitment +
                     last.collectiveExcessAfterCommitment,
                 last.collectiveExcess);
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
                   kappa.channelEventsBefore[c], kappa.channelMeanTimeBefore[c], kappa.collectiveExcessChannelBefore[c],
                   last.channelEventsAfter[c], last.channelMeanTimeAfter[c], last.collectiveExcessChannelAfter[c]);
      }
      std::print(stream, "\n");

      if (!kappa.hybrid.empty())
      {
        std::print(stream, "  Hybrid estimator of the collective excess at the handover times tau after commitment:\n");
        std::print(stream,
                   "    epsilon_conc(tau): the committor-smoothed cage count of the tagged molecule and the smoothed\n"
                   "    axial displacement of the other molecules of the column at tau (members with an event before\n"
                   "    tau enter with their first-event values), relative to kappa_S^committed, minus one;\n"
                   "    epsilon_Markov(tau): the expected sign of the first event from the full neighbourhood at tau\n"
                   "    (column plus the four lateral neighbours of A and of B) with the occupancy-resolved rates\n"
                   "    k(n -> n') of this run ('lateral': refill of A from each neighbour and escape from B to each\n"
                   "    neighbour end the race with s = 0; 'axial': the six axial channels only);\n"
                   "    epsilon_hyb = epsilon_conc + epsilon_Markov should not depend on tau once the hopping beyond\n"
                   "    tau is Markovian, and tends to the first-event epsilon = {:.5f} +/- {:.5f} for tau -> infinity.\n",
                   kappa.eventOrders.front().collectiveExcess, kappa.eventOrders.front().collectiveExcessError);
        std::print(stream, "    {:>8} {:>10} {:>9} {:>10} {:>10} {:>8} {:>8} {:>6}", "tau [ps]", "eps_conc", "+/-",
                   "tagged", "others", "p_event", "p_trunc", "states");
        if (kappa.hybrid.front().markovAvailable)
        {
          std::print(stream, " {:>10} {:>10} {:>10} {:>10} {:>8}", "Markov.lat", "Markov.ax", "eps_hyb", "eps_hyb.ax",
                     "p_fit");
        }
        std::print(stream, "\n");
        for (const TransmissionResult::Hybrid& hybrid : kappa.hybrid)
        {
          std::print(stream, "    {:>8.3f} {:>10.5f} {:>9.5f} {:>10.5f} {:>10.5f} {:>8.4f} {:>8.4f} {:>6}", hybrid.time,
                     hybrid.concerted, hybrid.concertedError, hybrid.concertedTagged - 1.0, hybrid.concertedOthers,
                     hybrid.fractionEventBefore, hybrid.fractionTruncated, hybrid.numberOfStates);
          if (hybrid.markovAvailable)
          {
            std::print(stream, " {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>8.4f}", hybrid.markovTable,
                       hybrid.markovTableAxial, hybrid.hybrid, hybrid.concerted + hybrid.markovTableAxial,
                       hybrid.unmeasuredWeight);
          }
          std::print(stream, "\n");
        }
        std::print(stream,
                   "    ('tagged' is the smooth analogue of f - 1, 'others' that of epsilon_before + epsilon_after;\n"
                   "     p_event: fraction of members whose first event preceded tau; p_trunc: fraction that did not\n"
                   "     reach tau (capped or never committed); states: distinct occupancy states at tau)\n");
        if (!kappa.hybrid.front().markovAvailable)
        {
          std::print(stream,
                     "    (no rates k(n -> n'): epsilon_Markov is not evaluated; the occupancy and neighbourhood\n"
                     "     histograms at every tau are in the JSON output, 'results.eventTerminated.hybrid')\n");
        }
        else
        {
          const OccupancyRates& rates = kappa.rates.value();
          std::print(stream, "    (rates: {} measured pairs", rates.measured.size());
          if (rates.fitValid)
          {
            std::print(stream, ", fit log k = c + a n - b n' with c = {:.4f}, a = {:.4f}, b = {:.4f} (rms {:.4f}); p_fit: weight of the states "
                               "that needed the fit;\n     with the fit for every pair: Markov.lat {:.5f}, Markov.ax "
                               "{:.5f}, eps_hyb {:.5f} at the last tau",
                       rates.fitConstant, rates.fitDeparture, rates.fitArrival, rates.fitRms,
                       kappa.hybrid.back().markovFit, kappa.hybrid.back().markovFitAxial, kappa.hybrid.back().hybridFit);
          }
          std::print(stream, ")\n");
        }
        // the dominant occupancy states at the last handover time
        const TransmissionResult::Hybrid& lastHandover = kappa.hybrid.back();
        std::print(stream, "    occupancy states (n_A', n_A, n_B, n_B') at tau = {:g} ps, weight fraction (members):\n",
                   lastHandover.time);
        for (std::size_t i = 0; i < std::min(12uz, lastHandover.occupancy.size()); ++i)
        {
          const auto& [state, weight] = lastHandover.occupancy[i];
          std::print(stream, "      ({:>2}, {:>2}, {:>2}, {:>2})  {:.4f} ({})\n", state.nAPrime, state.nA, state.nB,
                     state.nBPrime, weight, lastHandover.occupancyCounts[i]);
        }
        std::print(stream, "\n");
      }
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
      // relative errors from the magnitudes (a noisy kappa_S can be negative: max(kappa, 1e-300) then blew the
      // error up by 1e300, and a negative prefactor printed a negative error)
      const double kSError =
          std::abs(kS) * std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                                   std::pow(kappa.selfPlateauError / std::max(std::abs(kappa.selfPlateau), 1e-300), 2));
      std::print(stream, "    k_(A->B) = kappa_S k^TST:            {:.6e} +/- {:.2e} [1/ps]\n", kS, kSError);
      std::print(stream, "    D_S = kappa_S k^TST lambda^2:        {:.6e} +/- {:.2e} [m^2/s]\n",
                 kS * lambdaSquared * conversion, kSError * lambdaSquared * conversion);
      double kSSmooth{0.0}, kSSmoothError{0.0};
      if (kappa.smoothAvailable)
      {
        kSSmooth = kappa.selfSmoothPlateau * kTST;
        kSSmoothError =
            std::abs(kSSmooth) *
            std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                      std::pow(kappa.selfSmoothPlateauError / std::max(std::abs(kappa.selfSmoothPlateau), 1e-300), 2));
        std::print(stream, "    k_(A->B) = kappa_S^smooth k^TST:     {:.6e} +/- {:.2e} [1/ps]\n", kSSmooth, kSSmoothError);
        std::print(stream, "    D_S = kappa_S^smooth k^TST lambda^2: {:.6e} +/- {:.2e} [m^2/s]\n",
                   kSSmooth * lambdaSquared * conversion, kSSmoothError * lambdaSquared * conversion);
      }
      std::vector<double> kCErrors(numberOfShells, 0.0);
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        const double kC = kappa.collectivePlateau[s] * kTST;
        const double kCError =
            std::abs(kC) *
            std::sqrt(
                std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                std::pow(kappa.collectivePlateauError[s] / std::max(std::abs(kappa.collectivePlateau[s]), 1e-300), 2));
        kCErrors[s] = kCError;
        if (s + 1 == numberOfShells)
          std::print(stream, "    D_C = kappa_C k^TST lambda^2 (all):  {:.6e} +/- {:.2e} [m^2/s]\n"
                             "    (of the tagged component; the cross terms with other components are the\n"
                             "     per-component columns of the transmission file)\n",
                     kC * lambdaSquared * conversion, kCError * lambdaSquared * conversion);
        else
          std::print(stream, "    D_C (shell <= {}):                    {:.6e} +/- {:.2e} [m^2/s]\n", s,
                     kC * lambdaSquared * conversion, kCError * lambdaSquared * conversion);
      }
      // the event-terminated rates at every order k (the n-th is the headline value)
      struct EventRates
      {
        double kSelf{0.0}, kSelfError{0.0}, kCollective{0.0}, kCollectiveError{0.0};
      };
      auto eventRates = [&](const TransmissionResult::EventOrder& order) -> EventRates
      {
        EventRates rates;
        rates.kSelf = order.selfEvent * kTST;
        rates.kSelfError = std::abs(rates.kSelf) *
                           std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                                     std::pow(order.selfEventError / std::max(std::abs(order.selfEvent), 1e-300), 2));
        rates.kCollective = order.collectiveEvent * kTST;
        rates.kCollectiveError =
            std::abs(rates.kCollective) *
            std::sqrt(std::pow(kTSTError / std::max(kTST, 1e-300), 2) +
                      std::pow(order.collectiveEventError / std::max(std::abs(order.collectiveEvent), 1e-300), 2));
        return rates;
      };
      EventRates lastRates;
      if (eventMode)
      {
        lastRates = eventRates(kappa.lastEventOrder());
        std::print(stream, "    event-terminated estimator:\n");
        std::print(stream, "    k_(A->B) = kappa_S^event k^TST:      {:.6e} +/- {:.2e} [1/ps]\n", lastRates.kSelf,
                   lastRates.kSelfError);
        std::print(stream, "    D_S = kappa_S^event k^TST lambda^2:  {:.6e} +/- {:.2e} [m^2/s]\n",
                   lastRates.kSelf * lambdaSquared * conversion, lastRates.kSelfError * lambdaSquared * conversion);
        std::print(stream, "    D_C = kappa_C^event k^TST lambda^2:  {:.6e} +/- {:.2e} [m^2/s]\n",
                   lastRates.kCollective * lambdaSquared * conversion,
                   lastRates.kCollectiveError * lambdaSquared * conversion);
        if (kappa.eventOrders.size() > 1uz)
        {
          std::print(stream, "    per event order k:  {:>3} {:>14} {:>14} {:>14} {:>14}\n", "k", "k_(A->B) [1/ps]",
                     "+/-", "D_S [m^2/s]", "D_C [m^2/s]");
          for (const TransmissionResult::EventOrder& order : kappa.eventOrders)
          {
            const EventRates rates = eventRates(order);
            std::print(stream, "                        {:>3} {:>14.6e} {:>14.2e} {:>14.6e} {:>14.6e}\n", order.order,
                       rates.kSelf, rates.kSelfError, rates.kSelf * lambdaSquared * conversion,
                       rates.kCollective * lambdaSquared * conversion);
          }
        }
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
      if (kappa.smoothAvailable)
      {
        json["kSelfSmooth"] = kSSmooth;
        json["kSelfSmoothError"] = kSSmoothError;
        json["DSelfSmooth"] = kSSmooth * lambdaSquared * conversion;
        json["DSelfSmoothError"] = kSSmoothError * lambdaSquared * conversion;
      }
      for (std::size_t s = 0; s < numberOfShells; ++s)
      {
        const std::string key = (s + 1 == numberOfShells) ? "all" : std::to_string(s);
        json["DCollective"][key] = kappa.collectivePlateau[s] * kTST * lambdaSquared * conversion;
        // the same error as the text output (the k^TST error included)
        json["DCollectiveError"][key] = kCErrors[s] * lambdaSquared * conversion;
      }
      if (eventMode)
      {
        auto writeRates = [&](nlohmann::json& target, const EventRates& rates)
        {
          target["kSelfEvent"] = rates.kSelf;
          target["kSelfEventError"] = rates.kSelfError;
          target["DSelfEvent"] = rates.kSelf * lambdaSquared * conversion;
          target["DSelfEventError"] = rates.kSelfError * lambdaSquared * conversion;
          target["DCollectiveEvent"] = rates.kCollective * lambdaSquared * conversion;
          target["DCollectiveEventError"] = rates.kCollectiveError * lambdaSquared * conversion;
        };
        writeRates(json, lastRates);
        json["byEvent"] = nlohmann::json::array();
        for (const TransmissionResult::EventOrder& order : kappa.eventOrders)
        {
          nlohmann::json entry;
          entry["order"] = order.order;
          writeRates(entry, eventRates(order));
          json["byEvent"].push_back(std::move(entry));
        }
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
  if (kappa.smoothAvailable)
  {
    outputJson["results"]["kappaSelfSmooth"] = kappa.selfSmoothPlateau;
    outputJson["results"]["kappaSelfSmoothError"] = kappa.selfSmoothPlateauError;
  }
  if (!kappa.classes.empty())
  {
    nlohmann::json& classes = outputJson["results"]["occupancyClasses"];
    classes["reweightedKappa"] = kappa.classesReweightedKappa;
    classes["reweightedCoverage"] = kappa.classesReweightedCoverage;
    for (const TransmissionResult::OccupancyClass& item : kappa.classes)
    {
      nlohmann::json entry;
      entry["nFrom"] = item.nFrom;
      entry["nTo"] = item.nTo;
      entry["members"] = item.members;
      entry["weightFraction"] = item.weightFraction;
      entry["kappa"] = item.kappa;
      entry["kappaError"] = item.kappaError;
      entry["kappaSmooth"] = item.kappaSmooth;
      entry["kappaSmoothError"] = item.kappaSmoothError;
      entry["kappaCommitted"] = item.kappaCommitted;
      entry["pilotWeight"] = item.pilotWeight;
      entry["tstRatio"] = item.tstRatio;
      classes["classes"].push_back(std::move(entry));
    }
  }
  if (kappa.rates.has_value())
  {
    const OccupancyRates& rates = kappa.rates.value();
    nlohmann::json& json = outputJson["results"]["occupancyRates"];
    for (const auto& [pair, entry] : rates.measured)
    {
      nlohmann::json item;
      item["rate"] = entry.rate;
      item["weight"] = entry.weight;
      item["kappa"] = entry.kappa;
      item["tstRatio"] = entry.tstRatio;
      item["source"] = entry.source;
      json["measured"][std::format("{},{}", pair.first, pair.second)] = std::move(item);
    }
    json["fitValid"] = rates.fitValid;
    json["fitConstant"] = rates.fitConstant;
    json["fitDeparture"] = rates.fitDeparture;
    json["fitArrival"] = rates.fitArrival;
    json["fitRms"] = rates.fitRms;
  }
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
    static constexpr std::array<const char*, static_cast<std::size_t>(EventChannel::Count)> channelKeys = {
        "taggedBackHop",  "taggedForwardHop", "otherAtoB", "otherBtoA",
        "otherAtoAPrime", "otherAPrimeToA",   "otherBtoC", "otherCtoB"};
    // the quantities of one event order; the n-th order is written under the top-level keys (the run's headline
    // estimate), and every order k = 1..n under 'byEvent'
    auto writeOrder = [&](nlohmann::json& target, const TransmissionResult::EventOrder& order)
    {
      target["order"] = order.order;
      target["kappaSelfEvent"] = order.selfEvent;
      target["kappaSelfEventError"] = order.selfEventError;
      target["kappaCollectiveEvent"] = order.collectiveEvent;
      target["kappaCollectiveEventError"] = order.collectiveEventError;
      target["backCorrelationFactor"] = order.backCorrelationFactor;
      target["backCorrelationFactorError"] = order.backCorrelationFactorError;
      target["collectiveExcess"] = order.collectiveExcess;
      target["collectiveExcessError"] = order.collectiveExcessError;
      target["cappedFraction"] = order.cappedFraction;
      target["eventFreeCappedFraction"] = kappa.eventFreeCappedFraction;
      target["backCorrelationFactorCapCorrected"] = order.backCorrelationFactorCapCorrected;
      target["collectiveExcessCapCorrected"] = order.collectiveExcessCapCorrected;
      target["kappaSelfEventConditioned"] = order.selfEventConditioned;
      target["kappaSelfEventConditionedError"] = order.selfEventConditionedError;
      target["kappaCollectiveEventConditioned"] = order.collectiveEventConditioned;
      target["kappaCollectiveEventConditionedError"] = order.collectiveEventConditionedError;
      target["backCorrelationFactorConditioned"] = order.backCorrelationFactorConditioned;
      target["collectiveExcessConditioned"] = order.collectiveExcessConditioned;
      target["numberOfEventTerminatedTrajectories"] = order.numberOfEventTerminatedTrajectories;
      target["numberOfCappedTrajectories"] = order.numberOfCappedTrajectories;
      target["meanTerminationTime"] = order.meanTerminationTime;
      target["maximumTerminationTime"] = order.maximumTerminationTime;
      for (std::size_t c = 0; c < static_cast<std::size_t>(EventChannel::Count); ++c)
      {
        target["terminatingEvents"]["forward"][channelKeys[c]] = order.terminatingEvents[0][c];
        target["terminatingEvents"]["reverse"][channelKeys[c]] = order.terminatingEvents[1][c];
        target["allEvents"]["forward"][channelKeys[c]] = order.allEvents[0][c];
        target["allEvents"]["reverse"][channelKeys[c]] = order.allEvents[1][c];
        nlohmann::json& before = target["decomposition"]["channels"][channelKeys[c]]["beforeCommitment"];
        before["events"] = kappa.channelEventsBefore[c];
        before["meanTime"] = kappa.channelMeanTimeBefore[c];
        before["collectiveExcess"] = kappa.collectiveExcessChannelBefore[c];
        nlohmann::json& after = target["decomposition"]["channels"][channelKeys[c]]["afterCommitment"];
        after["events"] = order.channelEventsAfter[c];
        after["meanTime"] = order.channelMeanTimeAfter[c];
        after["collectiveExcess"] = order.collectiveExcessChannelAfter[c];
      }
      target["decomposition"]["backCorrelation"] = order.backCorrelationFactor - 1.0;
      target["decomposition"]["backCorrelationError"] = order.backCorrelationFactorError;
      target["decomposition"]["collectiveExcessBeforeCommitment"] = kappa.collectiveExcessBeforeCommitment;
      target["decomposition"]["collectiveExcessBeforeCommitmentError"] = kappa.collectiveExcessBeforeCommitmentError;
      target["decomposition"]["collectiveExcessAfterCommitment"] = order.collectiveExcessAfterCommitment;
      target["decomposition"]["collectiveExcessAfterCommitmentError"] = order.collectiveExcessAfterCommitmentError;
    };
    writeOrder(json, kappa.lastEventOrder());
    json["terminationEvents"] = settings.terminationEvents;
    json["kappaSelfCommittedFinal"] = kappa.selfCommitted.back();
    json["numberOfUncommittedTrajectories"] = kappa.numberOfUncommittedTrajectories;
    json["meanIntegratedTime"] = kappa.meanIntegratedTime;
    json["byEvent"] = nlohmann::json::array();
    for (const TransmissionResult::EventOrder& order : kappa.eventOrders)
    {
      nlohmann::json entry;
      writeOrder(entry, order);
      json["byEvent"].push_back(std::move(entry));
    }
    if (!kappa.hybrid.empty())
    {
      nlohmann::json& hybridJson = json["hybrid"];
      hybridJson["occupancyKey"] =
          "nA',nA,nB,nB' (frame of the tagged hop, tagged molecule counted in B): weight fraction of |qdot| over all "
          "members; the fractions sum to 1 - fractionEventBefore - fractionTruncated";
      hybridJson["meanOccupancy"] =
          static_cast<double>(system.numberOfMoleculesPerComponent[settings.taggedComponent]) /
          static_cast<double>(std::max(1uz, cageLattice.numberOfCages()));
      if (kappa.rates.has_value())
      {
        const OccupancyRates& rates = kappa.rates.value();
        for (const auto& [pair, entry] : rates.measured)
        {
          nlohmann::json item;
          item["rate"] = entry.rate;
          item["weight"] = entry.weight;
          item["source"] = entry.source;
          hybridJson["rates"]["measured"][std::format("{},{}", pair.first, pair.second)] = std::move(item);
        }
        hybridJson["rates"]["fitValid"] = rates.fitValid;
        hybridJson["rates"]["fitConstant"] = rates.fitConstant;
        hybridJson["rates"]["fitDeparture"] = rates.fitDeparture;
        hybridJson["rates"]["fitArrival"] = rates.fitArrival;
        hybridJson["rates"]["fitRms"] = rates.fitRms;
      }
      for (const TransmissionResult::Hybrid& hybrid : kappa.hybrid)
      {
        nlohmann::json entry;
        entry["time"] = hybrid.time;
        entry["collectiveExcessConcerted"] = hybrid.concerted;
        entry["collectiveExcessConcertedError"] = hybrid.concertedError;
        entry["concertedTagged"] = hybrid.concertedTagged - 1.0;
        entry["concertedOthers"] = hybrid.concertedOthers;
        entry["fractionEventBefore"] = hybrid.fractionEventBefore;
        entry["fractionTruncated"] = hybrid.fractionTruncated;
        entry["numberOfStates"] = hybrid.numberOfStates;
        entry["numberOfNeighbourhoodStates"] = hybrid.numberOfNeighbourhoodStates;
        if (hybrid.markovAvailable)
        {
          entry["collectiveExcessMarkov"] = hybrid.markovTable;
          entry["collectiveExcessMarkovAxial"] = hybrid.markovTableAxial;
          entry["collectiveExcessMarkovFit"] = hybrid.markovFit;
          entry["collectiveExcessMarkovFitAxial"] = hybrid.markovFitAxial;
          entry["unmeasuredWeight"] = hybrid.unmeasuredWeight;
          entry["collectiveExcessHybrid"] = hybrid.hybrid;
          entry["collectiveExcessHybridAxial"] = hybrid.concerted + hybrid.markovTableAxial;
          entry["collectiveExcessHybridFit"] = hybrid.hybridFit;
        }
        nlohmann::json occupancy = nlohmann::json::object();
        for (std::size_t i = 0; i < hybrid.occupancy.size(); ++i)
        {
          const auto& [state, weight] = hybrid.occupancy[i];
          occupancy[std::format("{},{},{},{}", state.nAPrime, state.nA, state.nB, state.nBPrime)] = {
              weight, hybrid.occupancyCounts[i]};
        }
        entry["occupancy"] = occupancy;
        hybridJson["handover"].push_back(entry);
      }
    }
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
