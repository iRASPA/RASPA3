module;

module integrators_compute;

import std;

import molecule;
import atom;
import atom_dynamics;
import double3;
import component;
import framework;
import forcefield;
import simd_quatd;
import double3x3;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import running_energy;
import integrators_cputime;

double Integrators::computeTranslationalKineticEnergy(std::span<const Molecule> moleculeData,
                                                      [[maybe_unused]] std::span<const Atom> moleculeAtomPositions,
                                                      std::span<const AtomDynamics> moleculeDynamics,
                                                      const std::vector<Component>& components,
                                                      const std::optional<Framework>& framework,
                                                      std::span<const Atom> frameworkAtomPositions,
                                                      std::span<const AtomDynamics> frameworkDynamics,
                                                      const ForceField* forceField,
                                                      std::span<const GroupState> groupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double energy{};
  std::size_t index{};
  std::size_t groupIndex{};
  for (const Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    if (component.isSemiFlexible() && !groupData.empty())
    {
      // Center-of-mass kinetic energy of each rigid group plus per-atom kinetic energy of the
      // flexible atoms.
      std::span<const AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          const GroupState& state = groupData[groupIndex + rigidRank];
          energy += 0.5 * group.mass * double3::dot(state.velocity, state.velocity);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          double mass = component.definedAtoms[i].second;
          energy += 0.5 * mass * double3::dot(span[i].velocity, span[i].velocity);
        }
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else if (component.rigid)
    {
      // Accumulate kinetic energy: 0.5 * mass * velocity squared
      energy += 0.5 * molecule.mass * double3::dot(molecule.velocity, molecule.velocity);
    }
    else
    {
      // Flexible molecule: all kinetic energy is carried by the atoms
      std::span<const AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); i++)
      {
        double mass = component.definedAtoms[i].second;
        energy += 0.5 * mass * double3::dot(span[i].velocity, span[i].velocity);
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && !framework->rigid)
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework kinetic energy requires force-field and dynamics data");
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      energy += 0.5 * mass * double3::dot(frameworkDynamics[i].velocity, frameworkDynamics[i].velocity);
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeTranslationalKineticEnergy += end - begin;
  return energy;
}

double Integrators::computeRotationalKineticEnergy(std::span<const Molecule> moleculeData,
                                                   const std::vector<Component> components,
                                                   std::span<const GroupState> groupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double3 ang_vel;
  double energy{};

  auto rotationalEnergy = [&](const simd_quatd& orientationMomentum, const simd_quatd& orientation,
                              const double3& inertiaVector, const double3& inverseInertiaVector)
  {
    simd_quatd p = orientationMomentum;
    p.r = -p.r;
    simd_quatd pq = p * orientation;
    double3 omega = 0.5 * double3(pq.ix, pq.iy, pq.iz) * inverseInertiaVector;
    return 0.5 * double3::dot(inertiaVector, omega * omega);
  };

  std::size_t groupIndex{};
  for (const Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];

    if (component.isSemiFlexible())
    {
      if (!groupData.empty())
      {
        std::size_t rigidRank{};
        for (const MoleculeGroup& group : component.groups)
        {
          if (group.rigid)
          {
            const GroupState& state = groupData[groupIndex + rigidRank];
            energy += rotationalEnergy(state.orientationMomentum, state.orientation, group.inertiaVector,
                                       group.inverseInertiaVector);
            ++rigidRank;
          }
        }
        groupIndex += component.numberOfRigidGroups();
      }
      continue;
    }

    // Flexible molecules carry no rigid-body orientation momentum
    if (!component.rigid) continue;

    // Get inertia and inverse inertia vectors for the molecule's component
    double3 inertiaVector = component.inertiaVector;
    double3 inverseInertiaVector = component.inverseInertiaVector;

    energy += rotationalEnergy(molecule.orientationMomentum, molecule.orientation, inertiaVector, inverseInertiaVector);
    (void)ang_vel;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeRotationalKineticEnergy += end - begin;
  return energy;
}

double3 Integrators::computeCenterOfMass(std::span<const Molecule> moleculeData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double3 com{};
  double totalMass{};

  for (const Molecule& molecule : moleculeData)
  {
    // Accumulate total mass and weighted positions
    totalMass += molecule.mass;
    com += molecule.mass * molecule.centerOfMassPosition;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeCenterOfMass += end - begin;

  // Return the center of mass position
  return com / totalMass;
}

// The velocity of the center of mass is the average velocity of all objects in the system weighted by their masses
double3 Integrators::computeCenterOfMassVelocity(std::span<const Molecule> moleculeData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double3 com_velocity{};
  double totalMass{};

  for (const Molecule& molecule : moleculeData)
  {
    // Accumulate total mass and weighted velocities
    totalMass += molecule.mass;
    com_velocity += molecule.mass * molecule.velocity;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeCenterOfMassVelocity += end - begin;
  // Return the center of mass velocity
  return com_velocity / totalMass;
}

double3 Integrators::computeLinearMomentum(std::span<const Molecule> moleculeData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  double3 com_momentum{};
  for (const Molecule& molecule : moleculeData)
  {
    // Accumulate linear momentum: mass * velocity
    com_momentum += molecule.mass * molecule.velocity;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeLinearMomentum += end - begin;
  return com_momentum;
}
