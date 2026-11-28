module;

module integrators_update;

import std;

import finite_value;

import ewald_structure_factor;

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
import cross_links;
import interactions_cross_link;
import polymer_links;
import interactions_polymer_link;
import thermostat;

// molecule.velocity / centerOfMassPosition are what VACF and MSD read. Rigid molecules integrate
// those fields directly. Semi-flexible and flexible molecules integrate groups and atoms; the
// molecule record has to be the mass-weighted centre of mass of those degrees of freedom.
void Integrators::syncNonRigidMoleculeCenterOfMass(Molecule& molecule, std::span<const Atom> atoms,
                                                   std::span<const AtomDynamics> dynamics, const Component& component,
                                                   std::span<const GroupState> groupData, std::size_t groupIndex)
{
  if (component.rigid) return;

  double3 com{};
  double3 momentum{};
  if (component.isSemiFlexible() && !groupData.empty())
  {
    std::size_t rigidRank{};
    for (const Fragment& group : component.fragmentGraph.fragments)
    {
      if (!group.isRigidBody()) continue;
      const GroupState& state = groupData[groupIndex + rigidRank];
      com += group.mass * state.centerOfMassPosition;
      momentum += group.mass * state.velocity;
      ++rigidRank;
    }
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      if (component.rigidFragmentContaining(i).has_value()) continue;
      const double mass = component.atomMass(i, atoms[i]);
      com += mass * atoms[i].position;
      momentum += mass * dynamics[i].velocity;
    }
  }
  else
  {
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      const double mass = component.atomMass(i, atoms[i]);
      com += mass * atoms[i].position;
      momentum += mass * dynamics[i].velocity;
    }
  }
  molecule.centerOfMassPosition = com * molecule.invMass;
  molecule.velocity = momentum * molecule.invMass;
}

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
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
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
        if (!component.rigidFragmentContaining(i).has_value())
        {
          span[i].velocity *= scaling.first;
        }
      }
      groupIndex += component.numberOfRigidFragments();
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
        std::span<const Atom> atoms = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (group.isRigidBody())
          {
            const GroupState& state = groupData[groupIndex + rigidRank];
            totalMass += group.mass;
            totalMomentum += group.mass * state.velocity;
            ++rigidRank;
          }
        }
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          if (!component.rigidFragmentContaining(i).has_value())
          {
            double mass = component.atomMass(i, atoms[i]);
            totalMass += mass;
            totalMomentum += mass * span[i].velocity;
          }
        }
        groupIndex += component.numberOfRigidFragments();
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
    totalVelocity = totalMass > 0.0 ? totalMomentum / totalMass : double3{};
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
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
        {
          groupData[groupIndex + rigidRank].velocity -= totalVelocity;
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidFragmentContaining(i).has_value())
        {
          span[i].velocity -= totalVelocity;
        }
      }
      groupIndex += component.numberOfRigidFragments();
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
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
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
        if (!component.rigidFragmentContaining(i).has_value())
        {
          span[i].position += dt * spanDynamics[i].velocity;
        }
      }
      syncNonRigidMoleculeCenterOfMass(molecule, span, spanDynamics, component, groupData, groupIndex);
      groupIndex += component.numberOfRigidFragments();
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
      syncNonRigidMoleculeCenterOfMass(molecule, span, spanDynamics, component, {}, 0);
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
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity -= 0.5 * dt * state.gradient / group.mass;
          state.orientationMomentum -= 0.5 * dt * state.orientationGradient;
          ++rigidRank;
        }
      }
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      std::span<const Atom> atoms = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidFragmentContaining(i).has_value())
        {
          // a massless (virtual) site takes no kick: its constraints carry it
          double mass = component.atomMass(i, atoms[i]);
          if (mass > 0.0) span[i].velocity -= 0.5 * dt * span[i].gradient / mass;
        }
      }
      syncNonRigidMoleculeCenterOfMass(molecule, atoms, span, component, groupData, groupIndex);
      groupIndex += component.numberOfRigidFragments();
      index += molecule.numberOfAtoms;
      continue;
    }

    molecule.velocity -= 0.5 * dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum -= 0.5 * dt * molecule.orientationGradient;

    // For flexible molecules the atomic velocities are the integration variables
    if (!component.rigid)
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      std::span<const Atom> atoms = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); i++)
      {
        // a massless (virtual) site takes no kick: its constraints carry it
        double mass = component.atomMass(i, atoms[i]);
        if (mass > 0.0) span[i].velocity -= 0.5 * dt * span[i].gradient / mass;
      }
      syncNonRigidMoleculeCenterOfMass(molecule, atoms, span, component, {}, 0);
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
      if (mass > 0.0) frameworkDynamics[i].velocity -= 0.5 * dt * frameworkDynamics[i].gradient / mass;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateVelocities += end - begin;
}

namespace
{
// The inverse mass of a bead; a massless site takes no part in the correction of a constraint.
double inverseMassOf(const Component& component, std::size_t bead, const Atom& atom)
{
  const double mass = component.atomMass(bead, atom);
  return mass > 0.0 ? 1.0 / mass : 0.0;
}

void constrainMoleculePositions(const Component& component, std::span<Atom> atoms, std::span<AtomDynamics> dynamics,
                                double dt)
{
  const std::vector<Component::DistanceConstraint>& constraints = component.distanceConstraints;

  // the bond vectors before the drift: the directions of the constraint forces of this step
  static thread_local std::vector<double3> reference;
  reference.resize(constraints.size());
  for (std::size_t k = 0; k != constraints.size(); ++k)
  {
    const std::size_t a = constraints[k].a, b = constraints[k].b;
    reference[k] = (atoms[a].position - atoms[b].position) - dt * (dynamics[a].velocity - dynamics[b].velocity);
  }

  const double inverseTimeStep = dt != 0.0 ? 1.0 / dt : 0.0;
  for (std::size_t iteration = 0; iteration != component.constraintMaximumIterations; ++iteration)
  {
    bool converged = true;
    for (std::size_t k = 0; k != constraints.size(); ++k)
    {
      const std::size_t a = constraints[k].a, b = constraints[k].b;
      const double lengthSquared = constraints[k].length * constraints[k].length;
      const double3 separation = atoms[a].position - atoms[b].position;
      const double difference = lengthSquared - double3::dot(separation, separation);
      // (a NaN is no converged constraint: the comparison is not defined to fail on one under -ffast-math, and a
      // diverged molecule went on as converged)
      if (!isFiniteValue(difference))
      {
        throw std::runtime_error(std::format(
            "[Molecular dynamics]: SHAKE met a position that is not a number in a molecule of component '{}'; "
            "reduce the time step.\n",
            component.name));
      }
      if (std::abs(difference) <= 2.0 * component.constraintTolerance * lengthSquared) continue;
      converged = false;

      const double inverseMassA = inverseMassOf(component, a, atoms[a]);
      const double inverseMassB = inverseMassOf(component, b, atoms[b]);
      const double projection = double3::dot(separation, reference[k]);
      // the bond turned (nearly) perpendicular to where it was a step ago: no correction along the old bond
      // restores its length
      if (!(std::abs(projection) > 1.0e-6 * lengthSquared) || !(inverseMassA + inverseMassB > 0.0))
      {
        throw std::runtime_error(std::format(
            "[Molecular dynamics]: SHAKE cannot restore the constrained distance between atoms {} and {} of a "
            "molecule of component '{}' (the bond rotated too far in one step, or both atoms are massless); "
            "reduce the time step.\n",
            a, b, component.name));
      }
      const double multiplier = difference / (2.0 * (inverseMassA + inverseMassB) * projection);
      atoms[a].position += (multiplier * inverseMassA) * reference[k];
      atoms[b].position -= (multiplier * inverseMassB) * reference[k];
      dynamics[a].velocity += (multiplier * inverseMassA * inverseTimeStep) * reference[k];
      dynamics[b].velocity -= (multiplier * inverseMassB * inverseTimeStep) * reference[k];
    }
    if (converged) return;
  }
  throw std::runtime_error(std::format(
      "[Molecular dynamics]: SHAKE did not converge to a relative tolerance of {:g} in {} iterations for a molecule "
      "of component '{}'; reduce the time step, or raise 'ConstraintMaximumIterations' or 'ConstraintTolerance'.\n",
      component.constraintTolerance, component.constraintMaximumIterations, component.name));
}
}  // namespace

void Integrators::constrainMoleculeVelocities(const Component& component, std::span<const Atom> atoms,
                                              std::span<AtomDynamics> dynamics)
{
  const std::vector<Component::DistanceConstraint>& constraints = component.distanceConstraints;
  if (constraints.empty()) return;

  // The scale the residual is measured against is that of the velocities coming in: the velocity along a
  // constraint is driven to zero, so a criterion relative to what is left of it would never be met.
  double scaleSquared{};
  for (const Component::DistanceConstraint& constraint : constraints)
  {
    scaleSquared = std::max({scaleSquared, double3::dot(dynamics[constraint.a].velocity, dynamics[constraint.a].velocity),
                             double3::dot(dynamics[constraint.b].velocity, dynamics[constraint.b].velocity)});
  }
  const double scale = std::sqrt(scaleSquared);

  for (std::size_t iteration = 0; iteration != component.constraintMaximumIterations; ++iteration)
  {
    bool converged = true;
    for (const Component::DistanceConstraint& constraint : constraints)
    {
      const std::size_t a = constraint.a, b = constraint.b;
      const double3 separation = atoms[a].position - atoms[b].position;
      const double3 relativeVelocity = dynamics[a].velocity - dynamics[b].velocity;
      const double rate = double3::dot(separation, relativeVelocity);
      const double separationSquared = double3::dot(separation, separation);
      if (!isFiniteValue(rate))
      {
        throw std::runtime_error(std::format(
            "[Molecular dynamics]: RATTLE met a velocity that is not a number in a molecule of component '{}'; "
            "reduce the time step.\n",
            component.name));
      }
      if (std::abs(rate) <= component.constraintTolerance * std::sqrt(separationSquared) * scale) continue;
      converged = false;

      const double inverseMassA = inverseMassOf(component, a, atoms[a]);
      const double inverseMassB = inverseMassOf(component, b, atoms[b]);
      if (!(inverseMassA + inverseMassB > 0.0) || !(separationSquared > 0.0))
      {
        throw std::runtime_error(std::format(
            "[Molecular dynamics]: RATTLE cannot constrain the velocities of atoms {} and {} of a molecule of "
            "component '{}' (the atoms coincide, or both are massless).\n",
            a, b, component.name));
      }
      const double multiplier = rate / ((inverseMassA + inverseMassB) * separationSquared);
      dynamics[a].velocity -= (multiplier * inverseMassA) * separation;
      dynamics[b].velocity += (multiplier * inverseMassB) * separation;
    }
    if (converged) return;
  }
  throw std::runtime_error(std::format(
      "[Molecular dynamics]: RATTLE did not converge to a relative tolerance of {:g} in {} iterations for a "
      "molecule of component '{}'; raise 'ConstraintMaximumIterations' or 'ConstraintTolerance'.\n",
      component.constraintTolerance, component.constraintMaximumIterations, component.name));
}

void Integrators::applyPositionConstraints(std::span<const Molecule> moleculeData,
                                           std::span<Atom> moleculeAtomPositions,
                                           std::span<AtomDynamics> moleculeDynamics,
                                           const std::vector<Component>& components, double dt)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  for (const Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    // FIXED terms that are no distance constraint would be left free: refuse the dynamics of such a molecule
    if (!component.unsupportedDynamicsConstraints.empty()) component.requireSupportedDynamicsConstraints();
    if (!component.distanceConstraints.empty())
    {
      constrainMoleculePositions(component, moleculeAtomPositions.subspan(index, molecule.numberOfAtoms),
                                 moleculeDynamics.subspan(index, molecule.numberOfAtoms), dt);
    }
    index += molecule.numberOfAtoms;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updatePositions += end - begin;
}

void Integrators::applyVelocityConstraints(std::span<const Molecule> moleculeData,
                                           std::span<const Atom> moleculeAtomPositions,
                                           std::span<AtomDynamics> moleculeDynamics,
                                           const std::vector<Component>& components)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  for (const Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    if (!component.distanceConstraints.empty())
    {
      constrainMoleculeVelocities(component, moleculeAtomPositions.subspan(index, molecule.numberOfAtoms),
                                  moleculeDynamics.subspan(index, molecule.numberOfAtoms));
    }
    index += molecule.numberOfAtoms;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateVelocities += end - begin;
}

namespace
{
// Body-frame angular velocity from the orientation momentum (the convention of the free rotor).
double3 bodyAngularVelocity(const simd_quatd& orientation, const simd_quatd& orientationMomentum,
                            const double3& inverseInertia)
{
  simd_quatd p = orientationMomentum;
  p.r = -p.r;
  const simd_quatd pq = p * orientation;
  return 0.5 * double3(pq.ix, pq.iy, pq.iz) * inverseInertia;
}

// Orientation momentum from the body-frame angular velocity (the inverse of 'bodyAngularVelocity', as the
// velocities are drawn).
simd_quatd orientationMomentumFromBodyAngularVelocity(const simd_quatd& q, const double3& I, const double3& w)
{
  simd_quatd p(0.0, 0.0, 0.0, 0.0);
  p.ix = 2.0 * (q.r * (I.x * w.x) - q.iz * (I.y * w.y) + q.iy * (I.z * w.z));
  p.iy = 2.0 * (q.iz * (I.x * w.x) + q.r * (I.y * w.y) - q.ix * (I.z * w.z));
  p.iz = 2.0 * (-q.iy * (I.x * w.x) + q.ix * (I.y * w.y) + q.r * (I.z * w.z));
  p.r = 2.0 * (-q.ix * (I.x * w.x) - q.iy * (I.y * w.y) - q.iz * (I.z * w.z));
  return p;
}
}  // namespace

void Integrators::setMoleculeVelocities(Molecule& molecule, std::span<const Atom> atoms,
                                        std::span<AtomDynamics> dynamics, const Component& component,
                                        std::span<GroupState> groupData, std::size_t groupIndex,
                                        std::span<const double3> velocities)
{
  if (velocities.size() != molecule.numberOfAtoms || atoms.size() != molecule.numberOfAtoms ||
      dynamics.size() != molecule.numberOfAtoms)
  {
    throw std::runtime_error("Setting the velocities of a molecule requires one velocity per atom");
  }

  // The rigid motion closest to the atomic velocities: the momentum and the angular momentum about the centre of
  // mass are kept, so velocities that are those of a rigid body are reproduced exactly.
  const auto setBody = [&](std::span<const std::size_t> members, const auto& massOf, double3& velocity,
                           simd_quatd& orientationMomentum, const simd_quatd& orientation, const double3& inertia,
                           const double3& inverseInertia)
  {
    double3 com{}, momentum{};
    double mass{};
    for (std::size_t k = 0; k != members.size(); ++k)
    {
      const double m = massOf(k);
      com += m * atoms[members[k]].position;
      momentum += m * velocities[members[k]];
      mass += m;
    }
    if (!(mass > 0.0)) return;
    com = com / mass;
    velocity = momentum / mass;
    if (members.size() < 2) return;
    double3 angularMomentum{};
    for (std::size_t k = 0; k != members.size(); ++k)
    {
      angularMomentum += massOf(k) * double3::cross(atoms[members[k]].position - com, velocities[members[k]] - velocity);
    }
    // space-fixed to body-fixed (principal axes), where the inertia is diagonal
    const double3 omega = (double3x3::buildRotationMatrix(orientation) * angularMomentum) * inverseInertia;
    orientationMomentum = orientationMomentumFromBodyAngularVelocity(orientation, inertia, omega);
  };

  if (component.rigid)
  {
    std::vector<std::size_t> members(molecule.numberOfAtoms);
    std::iota(members.begin(), members.end(), 0uz);
    setBody(members, [&](std::size_t k) { return component.atomMass(k, atoms[k]); }, molecule.velocity,
            molecule.orientationMomentum, molecule.orientation, component.inertiaVector,
            component.inverseInertiaVector);
    return;
  }

  const bool semiFlexible = component.isSemiFlexible() && !groupData.empty();
  if (semiFlexible)
  {
    std::size_t rigidRank{};
    for (const Fragment& group : component.fragmentGraph.fragments)
    {
      if (!group.isRigidBody()) continue;
      GroupState& state = groupData[groupIndex + rigidRank];
      setBody(group.atoms, [&](std::size_t k) { return group.atomMasses[k]; }, state.velocity,
              state.orientationMomentum, state.orientation, group.inertiaVector, group.inverseInertiaVector);
      ++rigidRank;
    }
  }
  for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
  {
    dynamics[i].velocity =
        semiFlexible && component.rigidFragmentContaining(i).has_value() ? double3(0.0, 0.0, 0.0) : velocities[i];
  }
  // velocities along a constrained bond are not those of the model: keep the part the constraints allow
  constrainMoleculeVelocities(component, atoms, dynamics);
  syncNonRigidMoleculeCenterOfMass(molecule, atoms, dynamics, component, groupData, groupIndex);
}

void Integrators::langevinVelocityStep(Thermostat& thermostat, double dt, std::span<Molecule> moleculeData,
                                       std::span<const Atom> moleculeAtomPositions,
                                       std::span<AtomDynamics> moleculeDynamics,
                                       const std::vector<Component>& components,
                                       const std::optional<Framework>& framework,
                                       std::span<const Atom> frameworkAtomPositions,
                                       std::span<AtomDynamics> frameworkDynamics, const ForceField* forceField,
                                       std::span<GroupState> groupData, std::span<GroupState> frameworkGroupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  const auto kineticEnergy = [&]()
  {
    return computeTranslationalKineticEnergy(moleculeData, moleculeAtomPositions, moleculeDynamics, components,
                                             framework, frameworkAtomPositions, frameworkDynamics, forceField,
                                             groupData, frameworkGroupData) +
           computeRotationalKineticEnergy(moleculeData, components, groupData, framework, frameworkGroupData);
  };
  const double kineticBefore = kineticEnergy();

  RandomNumber& random = thermostat.langevinRandom;
  const double kT = Units::KB * thermostat.temperature;
  const double c1 = std::exp(-thermostat.langevinFriction * dt);
  const double c2 = std::sqrt(std::max(0.0, 1.0 - c1 * c1));
  const auto gaussian3 = [&]() { return double3(random.Gaussian(), random.Gaussian(), random.Gaussian()); };
  // a rigid body: its centre of mass, and its angular velocity in the body frame with the principal inertia
  const auto kickBody = [&](double3& velocity, simd_quatd& orientationMomentum, const simd_quatd& orientation,
                            double mass, const double3& inertia, const double3& inverseInertia, bool rotates)
  {
    velocity = c1 * velocity + c2 * std::sqrt(kT / mass) * gaussian3();
    if (!rotates) return;
    const double3 omega = c1 * bodyAngularVelocity(orientation, orientationMomentum, inverseInertia) +
                          c2 * sqrt(inverseInertia * kT) * gaussian3();
    orientationMomentum = orientationMomentumFromBodyAngularVelocity(orientation, inertia, omega);
  };

  std::size_t index{};
  std::size_t groupIndex{};
  for (Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    std::span<AtomDynamics> dynamics = moleculeDynamics.subspan(index, molecule.numberOfAtoms);
    std::span<const Atom> atoms = moleculeAtomPositions.subspan(index, molecule.numberOfAtoms);

    if (component.rigid)
    {
      kickBody(molecule.velocity, molecule.orientationMomentum, molecule.orientation, molecule.mass,
               component.inertiaVector, component.inverseInertiaVector, molecule.numberOfAtoms > 1);
    }
    else
    {
      const bool semiFlexible = component.isSemiFlexible();
      if (semiFlexible)
      {
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (!group.isRigidBody()) continue;
          GroupState& state = groupData[groupIndex + rigidRank];
          kickBody(state.velocity, state.orientationMomentum, state.orientation, group.mass, group.inertiaVector,
                   group.inverseInertiaVector, true);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (semiFlexible && component.rigidFragmentContaining(i).has_value()) continue;
        const double mass = component.atomMass(i, atoms[i]);
        if (!(mass > 0.0)) continue;
        dynamics[i].velocity = c1 * dynamics[i].velocity + c2 * std::sqrt(kT / mass) * gaussian3();
      }
      // the kicks are drawn for free atoms: keep the part the constraints allow
      constrainMoleculeVelocities(component, atoms, dynamics);
      syncNonRigidMoleculeCenterOfMass(molecule, atoms, dynamics, component, groupData, groupIndex);
      if (semiFlexible) groupIndex += component.numberOfRigidFragments();
    }
    index += molecule.numberOfAtoms;
  }

  if (framework && framework->hasMobileAtoms())
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework Langevin thermostat requires force-field and dynamics data");
    }
    std::size_t rigidRank{};
    for (const FrameworkGroup& group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      GroupState& state = frameworkGroupData[rigidRank];
      kickBody(state.velocity, state.orientationMomentum, state.orientation, group.mass, group.inertiaVector,
               group.inverseInertiaVector, true);
      ++rigidRank;
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      if (!framework->isFlexibleAtom(i)) continue;
      const double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      if (!(mass > 0.0)) continue;
      frameworkDynamics[i].velocity = c1 * frameworkDynamics[i].velocity + c2 * std::sqrt(kT / mass) * gaussian3();
    }
  }

  thermostat.langevinHeat += kineticEnergy() - kineticBefore;

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.scaleVelocities += end - begin;
}

void Integrators::initializeMoleculeVelocity(RandomNumber& random, Molecule& molecule,
                                             std::span<const Atom> moleculeAtoms,
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
    if (moleculeDynamics.size() != molecule.numberOfAtoms || moleculeAtoms.size() != molecule.numberOfAtoms)
      throw std::runtime_error("Flexible-molecule dynamics span has the wrong size");
    for (std::size_t i = 0; i != moleculeDynamics.size(); i++)
    {
      // a massless (virtual) site has no thermal velocity (sqrt(kT / 0) is an infinity)
      const double mass = component.atomMass(i, moleculeAtoms[i]);
      moleculeDynamics[i].velocity = mass > 0.0 ? double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                                                      std::sqrt(Units::KB * temperature / mass)
                                                : double3(0.0, 0.0, 0.0);
    }
    // The Maxwell-Boltzmann distribution of a constrained molecule is that of its free atoms projected on the
    // velocities the constraints allow.
    constrainMoleculeVelocities(component, moleculeAtoms, moleculeDynamics);
    double3 momentum{};
    for (std::size_t i = 0; i != moleculeDynamics.size(); i++)
    {
      momentum += component.atomMass(i, moleculeAtoms[i]) * moleculeDynamics[i].velocity;
    }
    molecule.velocity = momentum * molecule.invMass;
    molecule.orientationMomentum = simd_quatd(0.0, 0.0, 0.0, 0.0);
  }
}

void Integrators::initializeGroupVelocity(RandomNumber& random, GroupState& state, const Fragment& group,
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
                                       std::span<Atom> moleculeAtomPositions, std::span<AtomDynamics> moleculeDynamics,
                                       const std::vector<Component>& components, double temperature,
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
    std::span<const Atom> atoms = moleculeAtomPositions.subspan(index, molecule.numberOfAtoms);

    if (component.isSemiFlexible())
    {
      // Rigid groups get a center-of-mass velocity and an orientation momentum; flexible atoms get
      // per-atom velocities. Rigid-group atoms carry no independent velocity.
      std::size_t rigidRank{};
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
        {
          initializeGroupVelocity(random, groupData[groupIndex + rigidRank], group, temperature);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (component.rigidFragmentContaining(i).has_value())
        {
          dynamics[i].velocity = double3(0.0, 0.0, 0.0);
          continue;
        }
        double mass = component.atomMass(i, atoms[i]);
        dynamics[i].velocity = mass > 0.0 ? double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                                                std::sqrt(Units::KB * temperature / mass)
                                          : double3(0.0, 0.0, 0.0);
      }
      constrainMoleculeVelocities(component, atoms, dynamics);
      syncNonRigidMoleculeCenterOfMass(molecule, atoms, dynamics, component, groupData, groupIndex);
      groupIndex += component.numberOfRigidFragments();
      index += molecule.numberOfAtoms;
      continue;
    }

    initializeMoleculeVelocity(random, molecule, atoms, dynamics, component, temperature);
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
      frameworkDynamics[i].velocity = mass > 0.0 ? double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                                                       std::sqrt(Units::KB * temperature / mass)
                                                 : double3(0.0, 0.0, 0.0);
    }
  }
}

void Integrators::createCartesianPositions(std::span<Molecule> moleculeData,
                                           std::span<Atom> moleculeAtomPositions, const std::vector<Component>& components,
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
        for (std::size_t g = 0; g != component.fragmentGraph.fragments.size(); ++g)
        {
          if (component.fragmentGraph.fragments[g].isRigidBody())
          {
            component.regenerateFragmentAtoms(groupData[groupIndex + rigidRank], g, span);
            ++rigidRank;
          }
        }
      }
      double3 com{};
      for (std::size_t i = 0; i != span.size(); i++)
      {
        com += component.atomMass(i, span[i]) * span[i].position;
      }
      molecule.centerOfMassPosition = com * molecule.invMass;
      groupIndex += component.numberOfRigidFragments();
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
        double mass = component.atomMass(i, span[i]);
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

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculeData, const std::vector<Component>& components,
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
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (!group.isRigidBody()) continue;
        GroupState& state = groupData[groupIndex + rigidRank];
        std::pair<simd_quatd, simd_quatd> pq = std::make_pair(state.orientationMomentum, state.orientation);
        pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, group.inverseInertiaVector);
        state.orientationMomentum = pq.first;
        state.orientation = pq.second;
        ++rigidRank;
      }
      groupIndex += component.numberOfRigidFragments();
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
                                                            const std::vector<Component>& components,
                                                            std::span<GroupState> groupData)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  std::size_t groupIndex{};
  // Update center of mass velocities and orientation momenta for each molecule
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    std::span<const AtomDynamics> spanDynamics = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
    const Component& component = components[molecule.componentId];

    if (!component.rigid)
    {
      // A flexible molecule has no orientation of its own (no orientation momentum); its velocity is that of its
      // momentum, of which the atoms of a rigid fragment of a semi-flexible molecule hold none themselves: their
      // group state does.
      const bool semiFlexible = component.isSemiFlexible() && !groupData.empty();
      double3 momentum{};
      for (std::size_t i = 0; i != span.size(); i++)
      {
        if (semiFlexible && component.rigidFragmentContaining(i).has_value()) continue;
        momentum += component.atomMass(i, span[i]) * spanDynamics[i].velocity;
      }
      if (semiFlexible)
      {
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (!group.isRigidBody()) continue;
          momentum += group.mass * groupData[groupIndex + rigidRank].velocity;
          ++rigidRank;
        }
      }
      molecule.velocity = momentum * molecule.invMass;
      molecule.orientationMomentum = simd_quatd{};
      if (component.isSemiFlexible()) groupIndex += component.numberOfRigidFragments();
      index += molecule.numberOfAtoms;
      continue;
    }

    double3 com{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double mass = component.atomMass(i, span[i]);
      com += mass * span[i].position;
    }
    com = com * molecule.invMass;

    double3 com_velocity{};
    double3 angularMomentum{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double mass = component.atomMass(i, span[i]);
      com_velocity += mass * spanDynamics[i].velocity;

      double3 dr = span[i].position - com;
      angularMomentum += mass * double3::cross(dr, spanDynamics[i].velocity);
    }
    molecule.velocity = com_velocity * molecule.invMass;

    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    // in the body frame, without the axes of no moment of inertia (a linear or point molecule), which the
    // integrator does not rotate about
    double3 bodyAngularMomentum = M * angularMomentum;
    if (component.inverseInertiaVector.x == 0.0) bodyAngularMomentum.x = 0.0;
    if (component.inverseInertiaVector.y == 0.0) bodyAngularMomentum.y = 0.0;
    if (component.inverseInertiaVector.z == 0.0) bodyAngularMomentum.z = 0.0;
    simd_quatd omega(0.0, bodyAngularMomentum);
    molecule.orientationMomentum = 2.0 * q * omega;

    index += molecule.numberOfAtoms;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionVelocities += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculeData,
                                                           std::span<Atom> moleculeAtomPositions,
                                                           std::span<const AtomDynamics> moleculeDynamics,
                                                           const std::vector<Component>& components,
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
      for (std::size_t g = 0; g != component.fragmentGraph.fragments.size(); ++g)
      {
        const Fragment& group = component.fragmentGraph.fragments[g];
        if (!group.isRigidBody()) continue;

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
      groupIndex += component.numberOfRigidFragments();
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
    const SimulationBox& simulationBox, const std::vector<Component>& components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    EwaldStructureFactor& trialEik,
    EwaldStructureFactor& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent, const std::optional<Framework>& framework,
    std::span<AtomDynamics> frameworkDynamics, const CrossLinkTable* crossLinks, PolymerLinkTable* polymerLinks)
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
      internal_energies += Interactions::computeIntraMolecularGradient(forceField, simulationBox,
                                                                       components[i].intraMolecularPotentials,
                                                                       span_molecules, moleculeAtomPositions,
                                                                       moleculeDynamics);
    }
    molecule_index += numberOfMoleculesPerComponent[i];
  }
  if (crossLinks && !crossLinks->empty())
  {
    internal_energies += Interactions::computeCrossLinkGradient(forceField, simulationBox, components,
                                                                numberOfMoleculesPerComponent, moleculeAtomPositions,
                                                                moleculeDynamics, *crossLinks);
  }
  if (polymerLinks && !polymerLinks->empty())
  {
    internal_energies += Interactions::computePolymerLinkGradient(forceField, simulationBox, components,
                                                                  numberOfMoleculesPerComponent, moleculeAtomPositions,
                                                                  moleculeDynamics, *polymerLinks, nullptr, nullptr,
                                                                  true);
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updateGradients += end - begin;
  return frameworkMoleculeEnergy + intermolecularEnergy + ewaldEnergy + internal_energies;
}
