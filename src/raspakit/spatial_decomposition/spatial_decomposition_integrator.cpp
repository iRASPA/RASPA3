module;

module spatial_decomposition_integrator;

import std;

import double3;
import double3x3;
import running_energy;
import atom;
import atom_dynamics;
import molecule;
import component;
import forcefield;
import system;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import force_engine;

namespace SpatialDecompositionIntegrator
{
RangeViews viewsOf(System& system, const MoleculeRange& range)
{
  const std::span<GroupState> groupData = system.spanOfGroupData();
  return {
      std::span<Molecule>(system.moleculeData).subspan(range.moleculeBegin, range.moleculeEnd - range.moleculeBegin),
      system.spanOfMoleculeAtoms().subspan(range.atomBegin, range.atomEnd - range.atomBegin),
      system.spanOfMoleculeDynamics().subspan(range.atomBegin, range.atomEnd - range.atomBegin),
      groupData.empty() ? std::span<GroupState>{}
                        : groupData.subspan(range.groupBegin, range.groupEnd - range.groupBegin)};
}

std::size_t rigidGroupsOf(const Component& component)
{
  return component.isSemiFlexible() ? component.numberOfRigidFragments() : 0;
}

/// Cuts the molecules into `parts` contiguous ranges balanced on the atom count. The atoms of the molecules are
/// stored consecutively in molecule order, the rigid-group states likewise (one per rigid fragment of a
/// semi-flexible component), which is what the Integrators functions assume as well.
std::vector<MoleculeRange> moleculeRanges(const System& system, std::size_t parts)
{
  const std::vector<Molecule>& molecules = system.moleculeData;
  const std::size_t numberOfAtoms = system.spanOfMoleculeAtoms().size();
  std::vector<MoleculeRange> ranges(std::max<std::size_t>(1, parts));
  std::size_t molecule{}, atom{}, group{};
  for (std::size_t part = 0; part < ranges.size(); ++part)
  {
    MoleculeRange& range = ranges[part];
    range.moleculeBegin = molecule;
    range.atomBegin = atom;
    range.groupBegin = group;
    const bool last = part + 1 == ranges.size();
    const std::size_t atomTarget = ((part + 1) * numberOfAtoms) / ranges.size();
    while (molecule < molecules.size() && (last || atom < atomTarget))
    {
      const Component& component = system.components[molecules[molecule].componentId];
      atom += molecules[molecule].numberOfAtoms;
      group += rigidGroupsOf(component);
      ++molecule;
    }
    range.moleculeEnd = molecule;
    range.atomEnd = atom;
    range.groupEnd = group;
  }
  return ranges;
}

/// Fills the table for the molecules of the range from the current atoms (the same sums as moleculeCenterOfMass).
void centersOfMass(const System& system, const MoleculeRange& range, bool withPositions, CenterOfMassTable& table)
{
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  const std::span<const GroupState> groupData = system.spanOfGroupData();
  std::size_t atomIndex = range.atomBegin;
  std::size_t groupIndex = range.groupBegin;
  for (std::size_t m = range.moleculeBegin; m < range.moleculeEnd; ++m)
  {
    const Molecule& molecule = system.moleculeData[m];
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      table.mass[m] = molecule.mass;
      table.position[m] = molecule.centerOfMassPosition;
      table.velocity[m] = molecule.velocity;
    }
    else
    {
      // positions are stored unwrapped per molecule, so no minimum image is needed
      const bool semiFlexible = component.isSemiFlexible() && !groupData.empty();
      double mass{};
      double3 weightedPosition{};
      double3 momentum{};
      for (std::size_t i = 0; i < molecule.numberOfAtoms; ++i)
      {
        const double atomMass = component.definedAtoms[i].second;
        mass += atomMass;
        if (withPositions) weightedPosition += atomMass * atoms[atomIndex + i].position;
        if (!semiFlexible || !component.rigidFragmentContaining(i).has_value())
          momentum += atomMass * dynamics[atomIndex + i].velocity;
      }
      if (semiFlexible)
      {
        // the atoms of a rigid group carry no velocities of their own; the group state does
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (!group.isRigidBody()) continue;
          momentum += group.mass * groupData[groupIndex + rigidRank].velocity;
          ++rigidRank;
        }
      }
      table.mass[m] = mass;
      table.position[m] = withPositions ? weightedPosition / mass : double3{};
      table.velocity[m] = momentum / mass;
    }
    atomIndex += molecule.numberOfAtoms;
    groupIndex += rigidGroupsOf(component);
  }
}

/// Kinetic virial sum_k M_k V_k V_k^T of the molecules of the range from the table (molecular coupling).
double3x3 tableKineticVirial(const MoleculeRange& range, const CenterOfMassTable& table)
{
  double3x3 stress{};
  for (std::size_t m = range.moleculeBegin; m < range.moleculeEnd; ++m)
  {
    const double mass = table.mass[m];
    const double3& v = table.velocity[m];
    stress.ax += mass * v.x * v.x;
    stress.ay += mass * v.x * v.y;
    stress.az += mass * v.x * v.z;
    stress.bx += mass * v.y * v.x;
    stress.by += mass * v.y * v.y;
    stress.bz += mass * v.y * v.z;
    stress.cx += mass * v.z * v.x;
    stress.cy += mass * v.z * v.y;
    stress.cz += mass * v.z * v.z;
  }
  return stress;
}

/// Translational and rotational kinetic energy of the molecules of a range.
std::pair<double, double> rangeKineticEnergies(System& system, const MoleculeRange& range)
{
  const RangeViews views = viewsOf(system, range);
  const double translational =
      Integrators::computeTranslationalKineticEnergy(views.molecules, views.atoms, views.dynamics, system.components,
                                                     std::nullopt, {}, {}, &system.forceField, views.groups, {});
  const double rotational =
      Integrators::computeRotationalKineticEnergy(views.molecules, system.components, views.groups, std::nullopt, {});
  return {translational, rotational};
}

// Velocity Verlet with the engine forces (Integrators::velocityVerlet with updateGradients replaced)
RunningEnergy engineVelocityVerlet(System& system, ForceEngine& engine)
{
  if (engine.usesResident())
  {
    // the whole step on the device (the host state is refreshed by the driver when it needs it)
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    RunningEnergy energies = engine.residentVelocityVerlet(system);
    integratorsCPUTime.velocityVerlet += std::chrono::steady_clock::now() - start;
    return energies;
  }

  // a local reference: a thread_local is not captured by the lambdas below, each worker thread would see its own
  static thread_local IntegratorWorkspace workspace;
  IntegratorWorkspace& w = workspace;
  w.prepare(system, engine.numberOfThreads());
  const std::vector<Component>& components = system.components;
  const ForceField* forceField = &system.forceField;
  const double dt = system.timeStep;
  const bool thermostat = system.thermostat.has_value();

  std::pair<double, double> scaling{1.0, 1.0};
  if (thermostat)
  {
    forEachRange(engine, w, [&](std::size_t r, const MoleculeRange& range)
                 { w.kinetic[r] = rangeKineticEnergies(system, range); });
    const auto [translational, rotational] = w.totalKinetic();
    scaling = system.thermostat->NoseHooverNVT(translational, rotational);
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  forEachRange(engine, w,
               [&](std::size_t, const MoleculeRange& range)
               {
                 const RangeViews v = viewsOf(system, range);
                 if (thermostat)
                 {
                   Integrators::scaleVelocities(v.molecules, v.atoms, v.dynamics, components, scaling, std::nullopt, {},
                                                v.groups, {});
                 }
                 Integrators::updateVelocities(v.molecules, v.atoms, v.dynamics, components, dt, std::nullopt, {}, {},
                                               forceField, v.groups, {});
                 Integrators::updatePositions(v.molecules, v.atoms, v.dynamics, components, dt, std::nullopt, {}, {},
                                              v.groups, {});
                 Integrators::noSquishFreeRotorOrderTwo(v.molecules, components, dt, v.groups, std::nullopt, {});
                 Integrators::createCartesianPositions(v.molecules, v.atoms, components, v.groups, std::nullopt, {},
                                                       {});
               });

  RunningEnergy runningEnergies = engine.computeGradients(system, true);

  forEachRange(engine, w,
               [&](std::size_t r, const MoleculeRange& range)
               {
                 const RangeViews v = viewsOf(system, range);
                 Integrators::updateCenterOfMassAndQuaternionGradients(v.molecules, v.atoms, v.dynamics, components,
                                                                       v.groups, std::nullopt, {}, {});
                 Integrators::updateVelocities(v.molecules, v.atoms, v.dynamics, components, dt, std::nullopt, {}, {},
                                               forceField, v.groups, {});
                 w.kinetic[r] = rangeKineticEnergies(system, range);
               });
  auto [translational, rotational] = w.totalKinetic();

  if (thermostat)
  {
    scaling = system.thermostat->NoseHooverNVT(translational, rotational);
    forEachRange(engine, w,
                 [&](std::size_t, const MoleculeRange& range)
                 {
                   const RangeViews v = viewsOf(system, range);
                   Integrators::scaleVelocities(v.molecules, v.atoms, v.dynamics, components, scaling, std::nullopt, {},
                                                v.groups, {});
                 });
    // the scaling multiplies every translational velocity by scaling.first and every orientation momentum by
    // scaling.second: the kinetic energies after it follow without another pass over the atoms
    translational *= scaling.first * scaling.first;
    rotational *= scaling.second * scaling.second;
    runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
  }
  runningEnergies.translationalKineticEnergy = translational;
  runningEnergies.rotationalKineticEnergy = rotational;

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.velocityVerlet += end - begin;
  return runningEnergies;
}
}  // namespace SpatialDecompositionIntegrator
