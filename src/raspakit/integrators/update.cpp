module;

module integrators_update;

import std;

import molecule;
import double3;
import atom;
import atom_dynamics;
import component;
import simd_quatd;
import double3x3;
import rigid;
import forcefield;
import simulationbox;
import running_energy;
import forcefield;
import interactions_ewald;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_internal;
import integrators_cputime;
import integrators_compute;
import randomnumbers;
import units;
import interpolation_energy_grid;
import framework;

void Integrators::scaleVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                  std::span<AtomDynamics> moleculeDynamics, const std::vector<Component>& components,
                                  std::pair<double, double> scaling, const std::optional<Framework>& framework,
                                  std::span<AtomDynamics> frameworkDynamics, std::span<GroupState> groupData,
                                  std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  // Scale velocities and orientation momenta of each molecule
  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    molecule.velocity *= scaling.first;
    molecule.orientationMomentum *= scaling.second;

    if (component.isSemiFlexible())
    {
      // Semi-flexible molecule: scale each rigid group's rigid-body state and the flexible atoms.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity *= scaling.first;
          state.orientationMomentum *= scaling.second;
          ++rigidRank;
        }
      }
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          span[i].velocity *= scaling.first;
        }
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else if (!component.rigid)
    {
      // For flexible molecules all degrees of freedom are translational: scale the atomic velocities
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (AtomDynamics& dynamics : span)
      {
        dynamics.velocity *= scaling.first;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && framework->hasMobileAtoms())
  {
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      frameworkGroupData[rigidRank].velocity *= scaling.first;
      frameworkGroupData[rigidRank].orientationMomentum *= scaling.second;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkDynamics.size(); ++i)
    {
      if (framework->isFlexibleAtom(i)) frameworkDynamics[i].velocity *= scaling.first;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.scaleVelocities += end - begin;
}

void Integrators::removeCenterOfMassVelocityDrift(std::span<Molecule> moleculeData,
                                                  std::span<Atom> moleculeAtomPositions,
                                                  std::span<AtomDynamics> moleculeDynamics,
                                                  const std::vector<Component>& components,
                                                  const std::optional<Framework>& framework,
                                                  std::span<const Atom> frameworkAtomPositions,
                                                  std::span<AtomDynamics> frameworkDynamics,
                                                  const ForceField* forceField, std::span<GroupState> groupData,
                                                  std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  const bool flexibleFramework = framework && framework->hasMobileAtoms();
  const bool hasSemiFlexible = std::ranges::any_of(
      moleculeData, [&](const Molecule& m) { return components[m.componentId].isSemiFlexible(); });

  double3 totalVelocity;
  if (flexibleFramework || hasSemiFlexible)
  {
    if (flexibleFramework && (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size()))
    {
      throw std::runtime_error("Flexible-framework center-of-mass removal requires force-field and dynamics data");
    }

    // General momentum accounting: rigid groups and flexible atoms of semi-flexible molecules
    // contribute per group / per atom; other molecules through their center-of-mass velocity.
    double3 totalMomentum{};
    double totalMass{};
    std::size_t index{};
    std::size_t groupIndex{};
    for (const Molecule& molecule : moleculeData)
    {
      const Component& component = components[molecule.componentId];
      if (component.isSemiFlexible())
      {
        std::span<const AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
        std::size_t rigidRank{};
        for (const MoleculeGroup& group : component.groups)
        {
          if (group.rigid)
          {
            const GroupState& state = groupData[groupIndex + rigidRank];
            totalMass += group.mass;
            totalMomentum += group.mass * state.velocity;
            ++rigidRank;
          }
        }
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          if (!component.rigidGroupContaining(i).has_value())
          {
            double mass = component.definedAtoms[i].second;
            totalMass += mass;
            totalMomentum += mass * span[i].velocity;
          }
        }
        groupIndex += component.numberOfRigidGroups();
      }
      else
      {
        totalMass += molecule.mass;
        totalMomentum += molecule.mass * molecule.velocity;
      }
      index += molecule.numberOfAtoms;
    }
    if (flexibleFramework)
    {
      std::size_t rigidRank{};
      for (const FrameworkGroup& group : framework->groups)
      {
        if (!group.isRigidBody()) continue;
        totalMass += group.mass;
        totalMomentum += group.mass * frameworkGroupData[rigidRank].velocity;
        ++rigidRank;
      }
      for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
      {
        if (!framework->isFlexibleAtom(i)) continue;
        double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
        totalMass += mass;
        totalMomentum += mass * frameworkDynamics[i].velocity;
      }
    }
    totalVelocity = totalMomentum / totalMass;
  }
  else
  {
    // Preserve the established rigid-framework/molecule-only arithmetic path.
    totalVelocity = computeCenterOfMassVelocity(moleculeData);
  }
  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    molecule.velocity -= totalVelocity;

    if (component.isSemiFlexible())
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          groupData[groupIndex + rigidRank].velocity -= totalVelocity;
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          span[i].velocity -= totalVelocity;
        }
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else if (!component.rigid)
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (AtomDynamics& dynamics : span)
      {
        dynamics.velocity -= totalVelocity;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (flexibleFramework)
  {
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      frameworkGroupData[rigidRank].velocity -= totalVelocity;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkDynamics.size(); ++i)
    {
      if (framework->isFlexibleAtom(i)) frameworkDynamics[i].velocity -= totalVelocity;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.removeCenterOfMassVelocity += end - begin;
}

void Integrators::updatePositions(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                  std::span<const AtomDynamics> moleculeDynamics,
                                  const std::vector<Component>& components, double dt,
                                  const std::optional<Framework>& framework,
                                  std::span<Atom> frameworkAtomPositions,
                                  std::span<const AtomDynamics> frameworkDynamics, std::span<GroupState> groupData,
                                  std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  // Update the center of mass positions for each molecule
  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];

    if (component.isSemiFlexible())
    {
      // Drift the center of mass of each rigid group; the atoms are rebuilt later from the group
      // state (after the free-rotor step) in 'createCartesianPositions'.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.centerOfMassPosition += dt * state.velocity;
          ++rigidRank;
        }
      }
      // Flexible atoms are integrated directly.
      std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
      std::span<const AtomDynamics> spanDynamics = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          span[i].position += dt * spanDynamics[i].velocity;
        }
      }
      groupIndex += component.numberOfRigidGroups();
      index += molecule.numberOfAtoms;
      continue;
    }

    molecule.centerOfMassPosition += dt * molecule.velocity;

    // For flexible molecules the atomic positions are the integration variables
    if (!component.rigid)
    {
      std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
      std::span<const AtomDynamics> spanDynamics = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); ++i)
      {
        span[i].position += dt * spanDynamics[i].velocity;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && framework->hasMobileAtoms())
  {
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      frameworkGroupData[rigidRank].centerOfMassPosition += dt * frameworkGroupData[rigidRank].velocity;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      if (framework->isFlexibleAtom(i))
      {
        frameworkAtomPositions[i].position += dt * frameworkDynamics[i].velocity;
      }
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updatePositions += end - begin;
}

void Integrators::updateVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                   std::span<AtomDynamics> moleculeDynamics, const std::vector<Component>& components,
                                   double dt, const std::optional<Framework>& framework,
                                   std::span<const Atom> frameworkAtomPositions,
                                   std::span<AtomDynamics> frameworkDynamics, const ForceField* forceField,
                                   std::span<GroupState> groupData, std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Update velocities and orientation momenta based on gradients
  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];

    if (component.isSemiFlexible())
    {
      // Half-kick each rigid group's center-of-mass velocity and orientation momentum, and the
      // flexible atoms' velocities.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity -= 0.5 * dt * state.gradient / group.mass;
          state.orientationMomentum -= 0.5 * dt * state.orientationGradient;
          ++rigidRank;
        }
      }
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
        {
          double mass = component.definedAtoms[i].second;
          span[i].velocity -= 0.5 * dt * span[i].gradient / mass;
        }
      }
      groupIndex += component.numberOfRigidGroups();
      index += molecule.numberOfAtoms;
      continue;
    }

    molecule.velocity -= 0.5 * dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum -= 0.5 * dt * molecule.orientationGradient;

    // For flexible molecules the atomic velocities are the integration variables
    if (!component.rigid)
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); i++)
      {
        double mass = component.definedAtoms[i].second;
        span[i].velocity -= 0.5 * dt * span[i].gradient / mass;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && framework->hasMobileAtoms())
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework velocity update requires force-field and dynamics data");
    }
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      GroupState& state = frameworkGroupData[rigidRank];
      state.velocity -= 0.5 * dt * state.gradient / group.mass;
      state.orientationMomentum -= 0.5 * dt * state.orientationGradient;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      if (!framework->isFlexibleAtom(i)) continue;
      double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      frameworkDynamics[i].velocity -= 0.5 * dt * frameworkDynamics[i].gradient / mass;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateVelocities += end - begin;
}

void Integrators::initializeMoleculeVelocity(RandomNumber& random, Molecule& molecule,
                                             std::span<AtomDynamics> moleculeDynamics, const Component& component,
                                             double temperature)
{
  if (component.rigid)
  {
    molecule.velocity = double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                        std::sqrt(Units::KB * temperature / molecule.mass);

    const double3 I = component.inertiaVector;
    const double3 invI = component.inverseInertiaVector;
    const double3 angularVelocity =
        double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * sqrt(invI * Units::KB * temperature);

    const simd_quatd q = molecule.orientation;
    molecule.orientationMomentum.ix =
        2.0 * (q.r * (I.x * angularVelocity.x) - q.iz * (I.y * angularVelocity.y) + q.iy * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iy =
        2.0 * (q.iz * (I.x * angularVelocity.x) + q.r * (I.y * angularVelocity.y) - q.ix * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iz =
        2.0 * (-q.iy * (I.x * angularVelocity.x) + q.ix * (I.y * angularVelocity.y) + q.r * (I.z * angularVelocity.z));
    molecule.orientationMomentum.r =
        2.0 * (-q.ix * (I.x * angularVelocity.x) - q.iy * (I.y * angularVelocity.y) - q.iz * (I.z * angularVelocity.z));
  }
  else
  {
    if (moleculeDynamics.size() != molecule.numberOfAtoms)
      throw std::runtime_error("Flexible-molecule dynamics span has the wrong size");
    double3 momentum{};
    for (std::size_t i = 0; i != moleculeDynamics.size(); i++)
    {
      const double mass = component.definedAtoms[i].second;
      moleculeDynamics[i].velocity =
          double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * std::sqrt(Units::KB * temperature / mass);
      momentum += mass * moleculeDynamics[i].velocity;
    }
    molecule.velocity = momentum * molecule.invMass;
    molecule.orientationMomentum = simd_quatd(0.0, 0.0, 0.0, 0.0);
  }
}

void Integrators::initializeGroupVelocity(RandomNumber& random, GroupState& state, const MoleculeGroup& group,
                                          double temperature)
{
  state.velocity = double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                   std::sqrt(Units::KB * temperature / group.mass);

  const double3 I = group.inertiaVector;
  const double3 invI = group.inverseInertiaVector;
  const double3 angularVelocity =
      double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * sqrt(invI * Units::KB * temperature);

  const simd_quatd q = state.orientation;
  state.orientationMomentum.ix =
      2.0 * (q.r * (I.x * angularVelocity.x) - q.iz * (I.y * angularVelocity.y) + q.iy * (I.z * angularVelocity.z));
  state.orientationMomentum.iy =
      2.0 * (q.iz * (I.x * angularVelocity.x) + q.r * (I.y * angularVelocity.y) - q.ix * (I.z * angularVelocity.z));
  state.orientationMomentum.iz =
      2.0 * (-q.iy * (I.x * angularVelocity.x) + q.ix * (I.y * angularVelocity.y) + q.r * (I.z * angularVelocity.z));
  state.orientationMomentum.r =
      2.0 * (-q.ix * (I.x * angularVelocity.x) - q.iy * (I.y * angularVelocity.y) - q.iz * (I.z * angularVelocity.z));
}

void Integrators::initializeFrameworkGroupVelocity(RandomNumber& random, GroupState& state,
                                                   const FrameworkGroup& group, double temperature)
{
  state.velocity = double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                   std::sqrt(Units::KB * temperature / group.mass);

  const double3 I = group.inertiaVector;
  const double3 invI = group.inverseInertiaVector;
  const double3 angularVelocity =
      double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * sqrt(invI * Units::KB * temperature);

  const simd_quatd q = state.orientation;
  state.orientationMomentum.ix =
      2.0 * (q.r * (I.x * angularVelocity.x) - q.iz * (I.y * angularVelocity.y) + q.iy * (I.z * angularVelocity.z));
  state.orientationMomentum.iy =
      2.0 * (q.iz * (I.x * angularVelocity.x) + q.r * (I.y * angularVelocity.y) - q.ix * (I.z * angularVelocity.z));
  state.orientationMomentum.iz =
      2.0 * (-q.iy * (I.x * angularVelocity.x) + q.ix * (I.y * angularVelocity.y) + q.r * (I.z * angularVelocity.z));
  state.orientationMomentum.r =
      2.0 * (-q.ix * (I.x * angularVelocity.x) - q.iy * (I.y * angularVelocity.y) - q.iz * (I.z * angularVelocity.z));
}

void Integrators::initializeVelocities(RandomNumber& random, std::span<Molecule> moleculeData,
                                       [[maybe_unused]] std::span<Atom> moleculeAtomPositions,
                                       std::span<AtomDynamics> moleculeDynamics,
                                       const std::vector<Component> components, double temperature,
                                       const std::optional<Framework>& framework,
                                       std::span<const Atom> frameworkAtomPositions,
                                       std::span<AtomDynamics> frameworkDynamics, const ForceField* forceField,
                                       std::span<GroupState> groupData, std::span<GroupState> frameworkGroupData)
{
  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    std::span<AtomDynamics> dynamics = moleculeDynamics.subspan(index, molecule.numberOfAtoms);

    if (component.isSemiFlexible())
    {
      // Rigid groups get a center-of-mass velocity and an orientation momentum; flexible atoms get
      // per-atom velocities. Rigid-group atoms carry no independent velocity.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          initializeGroupVelocity(random, groupData[groupIndex + rigidRank], group, temperature);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (component.rigidGroupContaining(i).has_value())
        {
          dynamics[i].velocity = double3(0.0, 0.0, 0.0);
          continue;
        }
        double mass = component.definedAtoms[i].second;
        dynamics[i].velocity =
            double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * std::sqrt(Units::KB * temperature / mass);
      }
      groupIndex += component.numberOfRigidGroups();
      index += molecule.numberOfAtoms;
      continue;
    }

    initializeMoleculeVelocity(random, molecule, dynamics, component, temperature);
    index += molecule.numberOfAtoms;
  }
  if (framework && framework->hasMobileAtoms())
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework velocity initialization requires force-field and dynamics data");
    }
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      initializeFrameworkGroupVelocity(random, frameworkGroupData[rigidRank], group, temperature);
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      if (framework->isFixedAtom(i) || framework->rigidGroupContaining(i).has_value())
      {
        frameworkDynamics[i].velocity = double3(0.0, 0.0, 0.0);
        continue;
      }
      double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      frameworkDynamics[i].velocity =
          double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * std::sqrt(Units::KB * temperature / mass);
    }
  }
}

void Integrators::createCartesianPositions(std::span<Molecule> moleculeData,
                                           std::span<Atom> moleculeAtomPositions, std::vector<Component> components,
                                           std::span<const GroupState> groupData,
                                           const std::optional<Framework>& framework,
                                           std::span<Atom> frameworkAtomPositions,
                                           std::span<const GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  std::size_t groupIndex{};
  // Convert molecule positions and orientations to atom positions
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);

    if (component.isSemiFlexible())
    {
      // Rebuild the rigid-group atoms from their rigid-body state; the flexible atoms are
      // authoritative. Synchronize the molecule center-of-mass record from all atoms. When no group
      // state is supplied (calls outside molecular dynamics) the atom positions are already
      // authoritative, so only the center-of-mass record is refreshed.
      if (!groupData.empty())
      {
        std::size_t rigidRank{};
        for (std::size_t g = 0; g != component.groups.size(); ++g)
        {
          if (component.groups[g].rigid)
          {
            component.regenerateGroupAtoms(groupData[groupIndex + rigidRank], g, span);
            ++rigidRank;
          }
        }
      }
      double3 com{};
      for (std::size_t i = 0; i != span.size(); i++)
      {
        com += component.definedAtoms[i].second * span[i].position;
      }
      molecule.centerOfMassPosition = com * molecule.invMass;
      groupIndex += component.numberOfRigidGroups();
      index += molecule.numberOfAtoms;
      continue;
    }

    if (component.rigid)
    {
      simd_quatd q = molecule.orientation;

      // Calculate positions of each atom in the molecule
      for (std::size_t i = 0; i != span.size(); i++)
      {
        span[i].position = molecule.centerOfMassPosition + q * component.atoms[i].position;
      }
    }
    else
    {
      // Flexible molecule: the atomic positions are authoritative,
      // synchronize the molecule record by recomputing the center of mass
      double3 com{};
      for (std::size_t i = 0; i != span.size(); i++)
      {
        double mass = component.definedAtoms[i].second;
        com += mass * span[i].position;
      }
      molecule.centerOfMassPosition = com * molecule.invMass;
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && !frameworkGroupData.empty())
  {
    std::size_t rigidRank{};
    for (std::size_t g = 0; g != framework->groups.size(); ++g)
    {
      if (!framework->groups[g].isRigidBody()) continue;
      framework->regenerateGroupAtoms(frameworkGroupData[rigidRank], g, frameworkAtomPositions);
      ++rigidRank;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.createCartesianPositions += end - begin;
}

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculeData, const std::vector<Component> components,
                                            double dt, std::span<GroupState> groupData,
                                            const std::optional<Framework>& framework,
                                            std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];

    if (component.isSemiFlexible())
    {
      // Free-rotor step for each rigid group about its own principal-axis inertia.
      std::size_t rigidRank{};
      for (const MoleculeGroup& group : component.groups)
      {
        if (!group.rigid) continue;
        GroupState& state = groupData[groupIndex + rigidRank];
        std::pair<simd_quatd, simd_quatd> pq = std::make_pair(state.orientationMomentum, state.orientation);
        pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, group.inverseInertiaVector);
        state.orientationMomentum = pq.first;
        state.orientation = pq.second;
        ++rigidRank;
      }
      groupIndex += component.numberOfRigidGroups();
      continue;
    }

    // Flexible molecules carry no rigid-body orientation
    if (!component.rigid) continue;

    double3 inverseInertiaVector = component.inverseInertiaVector;
    std::pair<simd_quatd, simd_quatd> pq = std::make_pair(molecule.orientationMomentum, molecule.orientation);
    pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, inverseInertiaVector);
    molecule.orientationMomentum = pq.first;
    molecule.orientation = pq.second;
  }
  if (framework && !frameworkGroupData.empty())
  {
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      GroupState& state = frameworkGroupData[rigidRank];
      std::pair<simd_quatd, simd_quatd> pq = std::make_pair(state.orientationMomentum, state.orientation);
      pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, group.inverseInertiaVector);
      state.orientationMomentum = pq.first;
      state.orientation = pq.second;
      ++rigidRank;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.noSquishFreeRotorOrderTwo += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculeData,
                                                            std::span<Atom> moleculeAtomPositions,
                                                            std::span<const AtomDynamics> moleculeDynamics,
                                                            std::vector<Component> components,
                                                            [[maybe_unused]] std::span<GroupState> groupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  // Update center of mass velocities and orientation momenta for each molecule
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    std::span<const AtomDynamics> spanDynamics = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);

    double3 com{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double mass = components[molecule.componentId].definedAtoms[i].second;
      com += mass * span[i].position;
    }
    com = com * molecule.invMass;

    double3 com_velocity{};
    double3 angularMomentum{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double mass = components[molecule.componentId].definedAtoms[i].second;
      com_velocity += mass * spanDynamics[i].velocity;

      double3 dr = span[i].position - com;
      angularMomentum += mass * double3::cross(dr, spanDynamics[i].velocity);
    }
    molecule.velocity = com_velocity * molecule.invMass;

    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    simd_quatd omega(0.0, M * angularMomentum);
    molecule.orientationMomentum = 2.0 * q * omega;

    index += molecule.numberOfAtoms;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionVelocities += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculeData,
                                                           std::span<Atom> moleculeAtomPositions,
                                                           std::span<const AtomDynamics> moleculeDynamics,
                                                           std::vector<Component> components,
                                                           std::span<GroupState> groupData,
                                                           const std::optional<Framework>& framework,
                                                           std::span<const AtomDynamics> frameworkDynamics,
                                                           std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  std::size_t groupIndex{};

  // Update gradients of center of mass and orientation for each molecule
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    std::span<const AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);

    if (component.isSemiFlexible())
    {
      // Callers outside molecular dynamics (e.g. smart MC gradient moves) carry no rigid-group
      // state; the per-atom gradients are already complete, so there is nothing to gather.
      if (groupData.empty())
      {
        index += molecule.numberOfAtoms;
        continue;
      }
      // Gather the net force and torque on each rigid group from its atoms' gradients. The flexible
      // atoms keep their per-atom gradients (used directly by 'updateVelocities').
      std::size_t rigidRank{};
      for (std::size_t g = 0; g != component.groups.size(); ++g)
      {
        const MoleculeGroup& group = component.groups[g];
        if (!group.rigid) continue;

        GroupState& state = groupData[groupIndex + rigidRank];
        double3 com_gradient{};
        for (std::size_t atom : group.atoms)
        {
          com_gradient += span[atom].gradient;
        }
        state.gradient = com_gradient;

        double3 torque{};
        simd_quatd q = state.orientation;
        double3x3 M = double3x3::buildRotationMatrix(q);
        for (std::size_t k = 0; k != group.atoms.size(); ++k)
        {
          std::size_t atom = group.atoms[k];
          double atomMass = component.definedAtoms[atom].second;
          double3 F = M * (span[atom].gradient - com_gradient * atomMass / group.mass);
          torque += double3::cross(F, group.bodyFixedPositions[k]);
        }
        state.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);
        ++rigidRank;
      }
      groupIndex += component.numberOfRigidGroups();
      index += molecule.numberOfAtoms;
      continue;
    }

    double3 com_gradient{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      com_gradient += span[i].gradient;
    }
    molecule.gradient = com_gradient;

    // Flexible molecules carry no rigid-body orientation: no torque on the quaternion
    if (!component.rigid)
    {
      molecule.orientationGradient = simd_quatd(0.0, 0.0, 0.0, 0.0);
      index += molecule.numberOfAtoms;
      continue;
    }

    double3 torque{};
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double atomMass = component.definedAtoms[i].second;
      double3 F = M * (span[i].gradient - com_gradient * atomMass * molecule.invMass);
      double3 dr = component.atoms[i].position;
      torque += double3::cross(F, dr);
    }

    // Compute orientation gradient based on torque
    molecule.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);

    index += molecule.numberOfAtoms;
  }

  if (framework && !frameworkGroupData.empty() && frameworkDynamics.size() == framework->atoms.size())
  {
    std::size_t rigidRank{};
    for (std::size_t g = 0; g != framework->groups.size(); ++g)
    {
      const FrameworkGroup& group = framework->groups[g];
      if (!group.isRigidBody()) continue;

      GroupState& state = frameworkGroupData[rigidRank];
      double3 com_gradient{};
      for (std::size_t atom : group.atoms)
      {
        com_gradient += frameworkDynamics[atom].gradient;
      }
      state.gradient = com_gradient;

      double3 torque{};
      const simd_quatd q = state.orientation;
      const double3x3 M = double3x3::buildRotationMatrix(q);
      for (std::size_t k = 0; k != group.atoms.size(); ++k)
      {
        const std::size_t atom = group.atoms[k];
        const double atomMass = group.atomMasses[k];
        const double3 F = M * (frameworkDynamics[atom].gradient - com_gradient * atomMass / group.mass);
        torque += double3::cross(F, group.bodyFixedPositions[k]);
      }
      state.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);
      ++rigidRank;
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionGradients += end - begin;
}

RunningEnergy Integrators::updateGradients(
    std::span<const Molecule> moleculeData, std::span<const Atom> moleculeAtomPositions,
    std::span<AtomDynamics> moleculeDynamics, std::span<const Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent, const std::optional<Framework>& framework,
    std::span<AtomDynamics> frameworkDynamics)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Initialize gradients to zero
  for (AtomDynamics& dynamics : moleculeDynamics)
  {
    dynamics.gradient = double3(0.0, 0.0, 0.0);
  }
  const bool flexibleFramework =
      framework && framework->hasMobileAtoms() && frameworkDynamics.size() == frameworkAtomPositions.size();
  if (flexibleFramework)
  {
    for (AtomDynamics& dynamics : frameworkDynamics)
    {
      dynamics.gradient = double3(0.0, 0.0, 0.0);
    }
  }

  // Compute gradients and energies due to interactions
  RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeGradient(
      forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, moleculeDynamics, interpolationGrids,
      framework, frameworkDynamics);
  RunningEnergy intermolecularEnergy = Interactions::computeInterMolecularGradient(
      forceField, simulationBox, moleculeAtomPositions, moleculeDynamics);
  double netChargeFramework = 0.0;
  for (const Atom& atom : frameworkAtomPositions)
  {
    netChargeFramework += atom.charge;
  }
  RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierGradient(
      eik_x, eik_y, eik_z, eik_xy, trialEik, fixedFrameworkStoredEik, forceField, simulationBox, components,
      numberOfMoleculesPerComponent, moleculeAtomPositions, moleculeDynamics, netChargeFramework, framework,
      frameworkAtomPositions, frameworkDynamics);

  RunningEnergy internal_energies{};
  if (flexibleFramework)
  {
    internal_energies += Interactions::computeFrameworkIntraMolecularGradient(
        forceField, *framework, simulationBox, frameworkAtomPositions, frameworkDynamics);
  }
  std::size_t molecule_index = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (numberOfMoleculesPerComponent[i] > 0)
    {
      std::span<const Molecule> span_molecules = {&moleculeData[molecule_index], numberOfMoleculesPerComponent[i]};
      internal_energies += Interactions::computeIntraMolecularGradient(
          components[i].intraMolecularPotentials, span_molecules, moleculeAtomPositions, moleculeDynamics);
    }
    molecule_index += numberOfMoleculesPerComponent[i];
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateGradients += end - begin;
  return frameworkMoleculeEnergy + intermolecularEnergy + ewaldEnergy + internal_energies;
}
