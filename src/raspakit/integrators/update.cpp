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
                                  std::span<AtomDynamics> frameworkDynamics)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  // Scale velocities and orientation momenta of each molecule
  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity *= scaling.first;
    molecule.orientationMomentum *= scaling.second;

    // For flexible molecules all degrees of freedom are translational: scale the atomic velocities
    if (!components[molecule.componentId].rigid)
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (AtomDynamics& dynamics : span)
      {
        dynamics.velocity *= scaling.first;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && !framework->rigid)
  {
    for (AtomDynamics& dynamics : frameworkDynamics)
    {
      dynamics.velocity *= scaling.first;
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
                                                  const ForceField* forceField)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  double3 totalVelocity;
  const bool flexibleFramework = framework && !framework->rigid;
  if (flexibleFramework)
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework center-of-mass removal requires force-field and dynamics data");
    }

    double3 totalMomentum{};
    double totalMass{};
    for (const Molecule& molecule : moleculeData)
    {
      totalMass += molecule.mass;
      totalMomentum += molecule.mass * molecule.velocity;
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      totalMass += mass;
      totalMomentum += mass * frameworkDynamics[i].velocity;
    }
    totalVelocity = totalMomentum / totalMass;
  }
  else
  {
    // Preserve the established rigid-framework/molecule-only arithmetic path.
    totalVelocity = computeCenterOfMassVelocity(moleculeData);
  }
  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity -= totalVelocity;

    if (!components[molecule.componentId].rigid)
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
    for (AtomDynamics& dynamics : frameworkDynamics)
    {
      dynamics.velocity -= totalVelocity;
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
                                  std::span<const AtomDynamics> frameworkDynamics)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  // Update the center of mass positions for each molecule
  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    molecule.centerOfMassPosition += dt * molecule.velocity;

    // For flexible molecules the atomic positions are the integration variables
    if (!components[molecule.componentId].rigid)
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
  if (framework && !framework->rigid)
  {
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      frameworkAtomPositions[i].position += dt * frameworkDynamics[i].velocity;
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.updatePositions += end - begin;
}

void Integrators::updateVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                   std::span<AtomDynamics> moleculeDynamics, const std::vector<Component>& components,
                                   double dt, const std::optional<Framework>& framework,
                                   std::span<const Atom> frameworkAtomPositions,
                                   std::span<AtomDynamics> frameworkDynamics, const ForceField* forceField)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Update velocities and orientation momenta based on gradients
  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity -= 0.5 * dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum -= 0.5 * dt * molecule.orientationGradient;

    // For flexible molecules the atomic velocities are the integration variables
    if (!components[molecule.componentId].rigid)
    {
      std::span<AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
      for (std::size_t i = 0; i != span.size(); i++)
      {
        double mass = components[molecule.componentId].definedAtoms[i].second;
        span[i].velocity -= 0.5 * dt * span[i].gradient / mass;
      }
    }
    index += molecule.numberOfAtoms;
  }
  if (framework && !framework->rigid)
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework velocity update requires force-field and dynamics data");
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
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

void Integrators::initializeVelocities(RandomNumber& random, std::span<Molecule> moleculeData,
                                       [[maybe_unused]] std::span<Atom> moleculeAtomPositions,
                                       std::span<AtomDynamics> moleculeDynamics,
                                       const std::vector<Component> components, double temperature,
                                       const std::optional<Framework>& framework,
                                       std::span<const Atom> frameworkAtomPositions,
                                       std::span<AtomDynamics> frameworkDynamics, const ForceField* forceField)
{
  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    std::span<AtomDynamics> dynamics = moleculeDynamics.subspan(index, molecule.numberOfAtoms);
    initializeMoleculeVelocity(random, molecule, dynamics, components[molecule.componentId], temperature);
    index += molecule.numberOfAtoms;
  }
  if (framework && !framework->rigid)
  {
    if (forceField == nullptr || frameworkDynamics.size() != frameworkAtomPositions.size())
    {
      throw std::runtime_error("Flexible-framework velocity initialization requires force-field and dynamics data");
    }
    for (std::size_t i = 0; i != frameworkAtomPositions.size(); ++i)
    {
      double mass = forceField->pseudoAtoms[frameworkAtomPositions[i].type].mass;
      frameworkDynamics[i].velocity =
          double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * std::sqrt(Units::KB * temperature / mass);
    }
  }
}

void Integrators::createCartesianPositions(std::span<Molecule> moleculeData,
                                           std::span<Atom> moleculeAtomPositions, std::vector<Component> components)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};
  // Convert molecule positions and orientations to atom positions
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);

    if (components[molecule.componentId].rigid)
    {
      simd_quatd q = molecule.orientation;

      // Calculate positions of each atom in the molecule
      for (std::size_t i = 0; i != span.size(); i++)
      {
        span[i].position = molecule.centerOfMassPosition + q * components[molecule.componentId].atoms[i].position;
      }
    }
    else
    {
      // Flexible molecule: the atomic positions are authoritative,
      // synchronize the molecule record by recomputing the center of mass
      double3 com{};
      for (std::size_t i = 0; i != span.size(); i++)
      {
        double mass = components[molecule.componentId].definedAtoms[i].second;
        com += mass * span[i].position;
      }
      molecule.centerOfMassPosition = com * molecule.invMass;
    }
    index += molecule.numberOfAtoms;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.createCartesianPositions += end - begin;
}

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculeData, const std::vector<Component> components,
                                            double dt)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (Molecule& molecule : moleculeData)
  {
    // Flexible molecules carry no rigid-body orientation
    if (!components[molecule.componentId].rigid) continue;

    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    std::pair<simd_quatd, simd_quatd> pq = std::make_pair(molecule.orientationMomentum, molecule.orientation);
    pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, inverseInertiaVector);
    molecule.orientationMomentum = pq.first;
    molecule.orientation = pq.second;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.noSquishFreeRotorOrderTwo += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculeData,
                                                            std::span<Atom> moleculeAtomPositions,
                                                            std::span<const AtomDynamics> moleculeDynamics,
                                                            std::vector<Component> components)
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
                                                           std::vector<Component> components)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::size_t index{};

  // Update gradients of center of mass and orientation for each molecule
  for (Molecule& molecule : moleculeData)
  {
    std::span<const AtomDynamics> span = std::span(&moleculeDynamics[index], molecule.numberOfAtoms);
    double3 com_gradient{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      com_gradient += span[i].gradient;
    }
    molecule.gradient = com_gradient;

    // Flexible molecules carry no rigid-body orientation: no torque on the quaternion
    if (!components[molecule.componentId].rigid)
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
      double atomMass = components[molecule.componentId].definedAtoms[i].second;
      double3 F = M * (span[i].gradient - com_gradient * atomMass * molecule.invMass);
      double3 dr = components[molecule.componentId].atoms[i].position;
      torque += double3::cross(F, dr);
    }

    // Compute orientation gradient based on torque
    molecule.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);

    index += molecule.numberOfAtoms;
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
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& totalEik,
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
      framework && !framework->rigid && frameworkDynamics.size() == frameworkAtomPositions.size();
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
      eik_x, eik_y, eik_z, eik_xy, totalEik, fixedFrameworkStoredEik, forceField, simulationBox, components,
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
