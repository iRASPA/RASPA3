module;

module molecular_dynamics_spatial_decomposition;

import std;

import hardware_info;
import archive;
import graceful_shutdown;
import system;
import randomnumbers;
import input_reader;
import component;
import averages;
import property_loading;
import units;
import simulationbox;
import forcefield;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import energy_dudlambda;
import running_energy;
import atom;
import atom_dynamics;
import molecule;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import mc_moves;
import mc_moves_cputime;
import mc_moves_statistics;
import integrators;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import thermobarostat;
import minimization_cell_layout;
import elastic_constants;
import spatial_decomposition_settings;
import force_engine;
import interactions_intermolecular;
import interactions_framework_molecule;

namespace
{
// ---------------------------------------------------------------------------------------------------------------
// Host-side integrator passes over molecule ranges
//
// Every pass of the velocity-Verlet step over the molecules (velocity and position updates, kinetic energies, the
// barostat coupling of the velocities and positions) is independent per molecule, so the driver cuts the molecules
// into contiguous ranges balanced on the atom count and runs the passes on the engine's worker team, one range per
// thread. The Integrators functions index the atoms and rigid-group states relative to the spans they are given,
// so they run unchanged on a range. Consecutive per-molecule passes between two global reductions are fused into
// one team task; the partial sums of the ranges are reduced in range order, so the results do not depend on the
// scheduling. The framework terms of the Integrators functions are not needed (the engine rejects frameworks).
// ---------------------------------------------------------------------------------------------------------------

/// A contiguous range of molecules with the matching ranges of atoms and rigid-group states.
struct MoleculeRange
{
  std::size_t moleculeBegin{}, moleculeEnd{};
  std::size_t atomBegin{}, atomEnd{};
  std::size_t groupBegin{}, groupEnd{};
};

/// The spans of a range (the rigid-group span is empty when the system keeps no group state).
struct RangeViews
{
  std::span<Molecule> molecules;
  std::span<Atom> atoms;
  std::span<AtomDynamics> dynamics;
  std::span<GroupState> groups;
};

RangeViews viewsOf(System& system, const MoleculeRange& range)
{
  const std::span<GroupState> groupData = system.spanOfGroupData();
  return {std::span<Molecule>(system.moleculeData).subspan(range.moleculeBegin, range.moleculeEnd - range.moleculeBegin),
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

/// Mass, centre of mass and centre-of-mass velocity per molecule (the coupled points of the molecular barostat
/// coupling), filled range by range and kept up to date through the velocity scalings of the step, so that the
/// atoms are read once per pass instead of once per use.
struct CenterOfMassTable
{
  std::vector<double> mass;
  std::vector<double3> position;
  std::vector<double3> velocity;
  void resize(std::size_t n)
  {
    mass.resize(n);
    position.resize(n);
    velocity.resize(n);
  }
};

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
  const double translational = Integrators::computeTranslationalKineticEnergy(
      views.molecules, views.atoms, views.dynamics, system.components, std::nullopt, {}, {}, &system.forceField,
      views.groups, {});
  const double rotational =
      Integrators::computeRotationalKineticEnergy(views.molecules, system.components, views.groups, std::nullopt, {});
  return {translational, rotational};
}

/// The per-step state of the integrator passes: the ranges, the centre-of-mass table and the partial sums per
/// range (persisting between steps so that nothing is reallocated).
struct IntegratorWorkspace
{
  std::vector<MoleculeRange> ranges;
  CenterOfMassTable table;
  std::vector<double3> comDisplacement;
  std::vector<std::pair<double, double>> kinetic;  ///< translational, rotational kinetic energy per range
  std::vector<double3x3> virial;                   ///< kinetic virial per range

  void prepare(const System& system, std::size_t parts)
  {
    ranges = moleculeRanges(system, parts);
    table.resize(system.moleculeData.size());
    kinetic.assign(ranges.size(), {});
    virial.assign(ranges.size(), {});
  }
  std::pair<double, double> totalKinetic() const
  {
    std::pair<double, double> total{};
    for (const auto& [translational, rotational] : kinetic)
    {
      total.first += translational;
      total.second += rotational;
    }
    return total;
  }
  double3x3 totalVirial() const
  {
    double3x3 total{};
    for (const double3x3& part : virial) total += part;
    return total;
  }
};

/// Runs body(r, range) for every range on the engine's worker team (range r on member r).
template <typename Body>
void forEachRange(ForceEngine& engine, IntegratorWorkspace& w, Body&& body)
{
  engine.runOnTeam(
      [&](std::size_t member, std::size_t members)
      {
        for (std::size_t r = member; r < w.ranges.size(); r += members) body(r, w.ranges[r]);
      });
}

// Barostat coupling of the velocities (molecules only; the driver rejects frameworks). With molecular coupling the
// cell acts on the centre-of-mass velocity V of every molecule only: each atom (or rigid group) of a non-rigid
// molecule receives the same increment (S - 1) V, so the velocities relative to the centre of mass are untouched;
// the table holds V of every molecule on entry and S V on return. With atomic coupling every flexible atom is
// scaled individually (the table is not used).
void applyVelocityMatrix(System& system, const MoleculeRange& range, const double3x3& matrix,
                         BarostatCoupling coupling, CenterOfMassTable& table)
{
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<GroupState> groupData = system.spanOfGroupData();
  std::size_t atomIndex = range.atomBegin;
  std::size_t groupIndex = range.groupBegin;
  for (std::size_t m = range.moleculeBegin; m < range.moleculeEnd; ++m)
  {
    Molecule& molecule = system.moleculeData[m];
    const Component& component = system.components[molecule.componentId];
    molecule.velocity = matrix * molecule.velocity;
    if (!component.rigid && coupling == BarostatCoupling::Molecular)
    {
      const double3 comVelocity = table.velocity[m];
      const double3 increment = matrix * comVelocity - comVelocity;
      table.velocity[m] = matrix * comVelocity;
      if (component.isSemiFlexible())
      {
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (group.isRigidBody())
          {
            groupData[groupIndex + rigidRank].velocity += increment;
            ++rigidRank;
          }
        }
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          if (!component.rigidFragmentContaining(i).has_value()) moleculeDynamics[atomIndex + i].velocity += increment;
        }
        groupIndex += component.numberOfRigidFragments();
      }
      else
      {
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i) moleculeDynamics[atomIndex + i].velocity += increment;
      }
    }
    else if (component.isSemiFlexible())
    {
      std::size_t rigidRank{};
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          state.velocity = matrix * state.velocity;
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidFragmentContaining(i).has_value())
          moleculeDynamics[atomIndex + i].velocity = matrix * moleculeDynamics[atomIndex + i].velocity;
      }
      groupIndex += component.numberOfRigidFragments();
    }
    else if (!component.rigid)
    {
      for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        moleculeDynamics[atomIndex + i].velocity = matrix * moleculeDynamics[atomIndex + i].velocity;
    }
    if (component.rigid) table.velocity[m] = molecule.velocity;
    atomIndex += molecule.numberOfAtoms;
  }
}

// Cell and position propagation. The coupled points (centres of mass of rigid molecules; with molecular coupling the
// centres of mass of all molecules, with atomic coupling the rigid groups and the flexible atoms individually) are
// propagated with the cell by 'propagateCellAndPosition'. With molecular coupling the atoms and rigid groups of a
// non-rigid molecule then follow their centre of mass, R_com' - R_com, plus the plain drift of their velocity
// relative to the centre of mass, dt (v - V): the internal geometry is not strained by the cell.
//
// With molecular coupling the coupled points are taken from the table (filled with positions by centersOfMass
// after the velocity half-kick); on return the table holds the propagated centres of mass and `comDisplacement`
// the shift R_com' - R_com of every molecule, which distributeCenterOfMassDisplacement applies to the atoms and
// rigid groups range by range. With atomic coupling the points are updated in place here.
void propagateCell(System& system, const double3x3& cellVelocity, BarostatCoupling coupling,
                   CenterOfMassTable& table, std::vector<double3>& comDisplacement)
{
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<GroupState> groupData = system.spanOfGroupData();

  // the coupled points are at most one per atom (plus the rigid groups); the buffers persist between steps so
  // that the NPT step does not reallocate them every time
  static thread_local std::vector<double3> positions;
  static thread_local std::vector<double3> velocities;
  static thread_local std::vector<double3*> targets;
  const std::size_t numberOfMolecules = system.moleculeData.size();
  positions.clear();
  velocities.clear();
  targets.clear();
  std::size_t atomIndex{};
  std::size_t groupIndex{};
  const bool molecular = coupling == BarostatCoupling::Molecular;
  if (molecular)
  {
    positions.assign(table.position.begin(), table.position.begin() + static_cast<std::ptrdiff_t>(numberOfMolecules));
    velocities.assign(table.velocity.begin(), table.velocity.begin() + static_cast<std::ptrdiff_t>(numberOfMolecules));
  }
  const std::size_t capacity = moleculeAtoms.size() + groupData.size();
  if (!molecular)
  {
    positions.reserve(capacity);
    velocities.reserve(capacity);
    targets.reserve(capacity);
  }
  for (Molecule& molecule : system.moleculeData)
  {
    if (molecular) break;
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      positions.push_back(molecule.centerOfMassPosition);
      velocities.push_back(molecule.velocity);
      targets.push_back(&molecule.centerOfMassPosition);
    }
    else if (component.isSemiFlexible())
    {
      std::size_t rigidRank{};
      for (const Fragment& group : component.fragmentGraph.fragments)
      {
        if (group.isRigidBody())
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
        if (!component.rigidFragmentContaining(i).has_value())
        {
          positions.push_back(moleculeAtoms[atomIndex + i].position);
          velocities.push_back(moleculeDynamics[atomIndex + i].velocity);
          targets.push_back(&moleculeAtoms[atomIndex + i].position);
        }
      }
      groupIndex += component.numberOfRigidFragments();
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

  double3x3 cell = system.simulationBox.cell;
  const bool upper = system.thermobarostat->cellType == CellMinimizationType::RegularUpperTriangle ||
                     system.thermobarostat->cellType == CellMinimizationType::MonoclinicUpperTriangle;
  propagateCellAndPosition(cell, positions, velocities, cellVelocity, system.timeStep, upper);
  if (!std::isfinite(cell.determinant()) || cell.determinant() <= 1.0e-10)
    throw std::runtime_error("Thermobarostat produced an invalid or singular cell");
  if (molecular)
  {
    comDisplacement.resize(numberOfMolecules);
    for (std::size_t m = 0; m < numberOfMolecules; ++m)
    {
      Molecule& molecule = system.moleculeData[m];
      if (system.components[molecule.componentId].rigid)
      {
        molecule.centerOfMassPosition = positions[m];
        comDisplacement[m] = double3{};
      }
      else
      {
        comDisplacement[m] = positions[m] - table.position[m];
      }
      table.position[m] = positions[m];
    }
  }
  else
  {
    for (std::size_t i = 0; i != targets.size(); ++i)
    {
      if (targets[i] != nullptr) *targets[i] = positions[i];
    }
  }
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

// Molecular coupling: the atoms and rigid groups of the non-rigid molecules of the range follow their centre of
// mass, R_com' - R_com, plus the plain drift of their velocity relative to the centre of mass, dt (v - V).
void distributeCenterOfMassDisplacement(System& system, const MoleculeRange& range, const CenterOfMassTable& table,
                                        const std::vector<double3>& comDisplacement)
{
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<GroupState> groupData = system.spanOfGroupData();
  const double dt = system.timeStep;
  std::size_t atomIndex = range.atomBegin;
  std::size_t groupIndex = range.groupBegin;
  for (std::size_t m = range.moleculeBegin; m < range.moleculeEnd; ++m)
  {
    const Molecule& molecule = system.moleculeData[m];
    const Component& component = system.components[molecule.componentId];
    if (!component.rigid)
    {
      const double3 displacement = comDisplacement[m];
      const double3 comVelocity = table.velocity[m];
      if (component.isSemiFlexible())
      {
        std::size_t rigidRank{};
        for (const Fragment& group : component.fragmentGraph.fragments)
        {
          if (group.isRigidBody())
          {
            GroupState& state = groupData[groupIndex + rigidRank];
            state.centerOfMassPosition += displacement + dt * (state.velocity - comVelocity);
            ++rigidRank;
          }
        }
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
        {
          if (!component.rigidFragmentContaining(i).has_value())
            moleculeAtoms[atomIndex + i].position +=
                displacement + dt * (moleculeDynamics[atomIndex + i].velocity - comVelocity);
        }
      }
      else
      {
        for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
          moleculeAtoms[atomIndex + i].position +=
              displacement + dt * (moleculeDynamics[atomIndex + i].velocity - comVelocity);
      }
    }
    atomIndex += molecule.numberOfAtoms;
    groupIndex += rigidGroupsOf(component);
  }
}

// Whether the per-component energy decomposition can be taken from the engine's running energies instead of the
// exact O(N^2) code. The engine rejects frameworks, external fields, polarization and cross-links, so the status is
// the intra-molecular terms per component plus the inter-molecular VDW, real-space Coulomb and reciprocal terms;
// with a single component those all belong to the (0, 0) pair. The reciprocal (mesh) sum cannot be attributed to
// component pairs, so mixtures keep the exact decomposition.
bool engineProvidesEnergyStatus(const System& system) { return system.components.size() == 1; }

EnergyStatus energyStatusFromRunningEnergies(const System& system)
{
  const RunningEnergy& e = system.runningEnergies;
  EnergyStatus status(1, system.framework.has_value() ? 1 : 0, system.components.size());

  EnergyIntra& intra = status.intraComponentEnergies[0];
  intra.bond = e.bond;
  intra.ureyBradley = e.ureyBradley;
  intra.bend = e.bend;
  intra.inversionBend = e.inversionBend;
  intra.outOfPlaneBend = e.outOfPlaneBend;
  intra.torsion = e.torsion;
  intra.improperTorsion = e.improperTorsion;
  intra.bondBond = e.bondBond;
  intra.bondBend = e.bondBend;
  intra.bondTorsion = e.bondTorsion;
  intra.bendBend = e.bendBend;
  intra.bendTorsion = e.bendTorsion;
  intra.vanDerWaals = e.intraVDW;
  intra.coulomb = e.intraCoul;

  // same convention as the exact code: self and exclusion terms are part of the Fourier entry
  EnergyInter& inter = status.componentEnergy(0, 0);
  inter.VanDerWaals = EnergyDuDlambda(e.moleculeMoleculeVDW, 0.0);
  inter.VanDerWaalsTailCorrection = EnergyDuDlambda(e.tail, 0.0);
  inter.CoulombicReal = EnergyDuDlambda(e.moleculeMoleculeCharge, 0.0);
  inter.CoulombicFourier = EnergyDuDlambda(e.ewald_fourier + e.ewald_self + e.ewald_exclusion, 0.0);

  status.translationalKineticEnergy = e.translationalKineticEnergy;
  status.rotationalKineticEnergy = e.rotationalKineticEnergy;
  status.noseHooverEnergy = e.NoseHooverEnergy;
  status.sumTotal();
  return status;
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
                 Integrators::createCartesianPositions(v.molecules, v.atoms, components, v.groups, std::nullopt, {}, {});
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

// Thermobarostat step (NPT / NPT-PR) with the engine forces. The barostat is driven by the virial of the
// points it couples to (molecular coupling: the centres of mass of all molecules; atomic coupling: centers of
// rigid molecules and rigid groups, flexible atoms individually), obtained from the engine's molecular pressure
// tensor with 'computeBarostatVirial'; its kinetic partner is 'computeMolecularKineticVirial'. The reported
// pressure is the estimator of the same coupling (see 'barostatPressureTensor'), so its average is the set point.
RunningEnergy engineThermobarostatVelocityVerlet(System& system, ForceEngine& engine)
{
  if (engine.usesResident())
  {
    // the whole step on the device (isotropic barostat with molecular coupling; the engine refuses the resident
    // integrator for the other barostats), the chains on the host from the device reductions
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    RunningEnergy energies = engine.residentVelocityVerlet(system);
    integratorsCPUTime.velocityVerlet += std::chrono::steady_clock::now() - start;
    return energies;
  }

  Thermobarostat& barostat = *system.thermobarostat;
  // a local reference: a thread_local is not captured by the lambdas below, each worker thread would see its own
  static thread_local IntegratorWorkspace workspace;
  IntegratorWorkspace& w = workspace;
  w.prepare(system, engine.numberOfThreads());
  const std::vector<Component>& components = system.components;
  const ForceField* forceField = &system.forceField;
  const double dt = system.timeStep;
  const BarostatCoupling coupling = barostat.coupling;
  const bool molecular = coupling == BarostatCoupling::Molecular;
  const bool thermostat = system.thermostat.has_value();

  const double3x3 pressureBefore = computeBarostatVirial(system, engine.molecularPressureTensor(), coupling);
  // one pass: the centre-of-mass velocities of the molecules and their kinetic virial (molecular coupling), the
  // kinetic energies for the thermostat
  forEachRange(engine, w,
               [&](std::size_t r, const MoleculeRange& range)
               {
                 if (molecular)
                 {
                   centersOfMass(system, range, false, w.table);
                   w.virial[r] = tableKineticVirial(range, w.table);
                 }
                 if (thermostat) w.kinetic[r] = rangeKineticEnergies(system, range);
               });
  const double3x3 kineticBefore = molecular ? w.totalVirial() : computeMolecularKineticVirial(system, coupling);

  const double barostatKinetic = barostat.barostatKineticEnergy();
  const double chainScale = barostat.chainStep(barostatKinetic);
  barostat.logVolumeVelocity *= chainScale;
  barostat.cellVelocity = barostat.cellVelocity * chainScale;

  // the thermostat scaling is applied in the next pass, together with the barostat coupling
  std::pair<double, double> scaling{1.0, 1.0};
  if (thermostat)
  {
    const auto [translational, rotational] = w.totalKinetic();
    scaling = system.thermostat->NoseHooverNVT(translational, rotational);
  }

  double3x3 cellRate{};
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
  {
    const double mtkFactor =
        1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    // x = ln(V) = 3 epsilon: xddot = 3 G_epsilon / W with G_epsilon = virial + alpha 2K - 3 P V (MTK)
    const double scalarForce = 3.0 *
                               (pressureBefore.trace() + mtkFactor * kineticBefore.trace() -
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
        cellForce(pressureBefore, correctedKinetic, system.simulationBox.volume, barostat.pressure, barostat.cellMass,
                  barostat.cellType, barostat.monoclinicAngle);
    barostat.cellVelocity += 0.5 * system.timeStep * acceleration;
    barostat.cellVelocity = projectCellTensor(barostat.cellVelocity, barostat.cellType, barostat.monoclinicAngle);
    cellRate = barostat.cellVelocity;
  }

  const double3x3 propagator = velocityPropagator(cellRate, 0.5 * dt, barostat.translationalDegreesOfFreedom);
  // one pass: thermostat scaling, barostat coupling of the velocities, the half-kick, and (molecular coupling) the
  // centres of mass with their velocities as the coupled points of the cell propagation
  forEachRange(engine, w,
               [&](std::size_t, const MoleculeRange& range)
               {
                 const RangeViews v = viewsOf(system, range);
                 if (thermostat)
                 {
                   Integrators::scaleVelocities(v.molecules, v.atoms, v.dynamics, components, scaling, std::nullopt, {},
                                                v.groups, {});
                   // every translational velocity is scaled by scaling.first, so the centres of mass are too
                   if (molecular)
                   {
                     for (std::size_t m = range.moleculeBegin; m < range.moleculeEnd; ++m)
                       w.table.velocity[m] *= scaling.first;
                   }
                 }
                 applyVelocityMatrix(system, range, propagator, coupling, w.table);
                 Integrators::updateVelocities(v.molecules, v.atoms, v.dynamics, components, dt, std::nullopt, {}, {},
                                               forceField, v.groups, {});
                 if (molecular) centersOfMass(system, range, true, w.table);
               });
  propagateCell(system, cellRate, coupling, w.table, w.comDisplacement);
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
    barostat.logVolumePosition += dt * barostat.logVolumeVelocity;
  forEachRange(engine, w,
               [&](std::size_t, const MoleculeRange& range)
               {
                 const RangeViews v = viewsOf(system, range);
                 if (molecular) distributeCenterOfMassDisplacement(system, range, w.table, w.comDisplacement);
                 Integrators::noSquishFreeRotorOrderTwo(v.molecules, components, dt, v.groups, std::nullopt, {});
                 Integrators::createCartesianPositions(v.molecules, v.atoms, components, v.groups, std::nullopt, {}, {});
               });

  RunningEnergy energies = engine.computeGradients(system, true);

  // one pass: the second half-kick, the barostat coupling of the velocities, the centre-of-mass velocities with
  // their kinetic virial (molecular coupling) and the kinetic energies
  forEachRange(engine, w,
               [&](std::size_t r, const MoleculeRange& range)
               {
                 const RangeViews v = viewsOf(system, range);
                 Integrators::updateCenterOfMassAndQuaternionGradients(v.molecules, v.atoms, v.dynamics, components,
                                                                       v.groups, std::nullopt, {}, {});
                 Integrators::updateVelocities(v.molecules, v.atoms, v.dynamics, components, dt, std::nullopt, {}, {},
                                               forceField, v.groups, {});
                 if (molecular) centersOfMass(system, range, false, w.table);
                 applyVelocityMatrix(system, range, propagator, coupling, w.table);
                 if (molecular) w.virial[r] = tableKineticVirial(range, w.table);
                 w.kinetic[r] = rangeKineticEnergies(system, range);
               });

  {
    const auto [translational, rotational] = w.totalKinetic();
    energies.translationalKineticEnergy = translational;
    energies.rotationalKineticEnergy = rotational;
  }
  const double3x3 pressureAfter = computeBarostatVirial(system, engine.molecularPressureTensor(), coupling);
  const double3x3 kineticAfter = molecular ? w.totalVirial() : computeMolecularKineticVirial(system, coupling);
  if (molecularDynamicsUsesIsotropicBarostat(barostat.ensemble))
  {
    const double mtkFactor =
        1.0 + 3.0 / static_cast<double>(std::max<std::size_t>(1, barostat.translationalDegreesOfFreedom));
    // x = ln(V) = 3 epsilon: xddot = 3 G_epsilon / W with G_epsilon = virial + alpha 2K - 3 P V (MTK)
    const double scalarForce = 3.0 *
                               (pressureAfter.trace() + mtkFactor * kineticAfter.trace() -
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
    barostat.cellVelocity += 0.5 * system.timeStep *
                             cellForce(pressureAfter, correctedKinetic, system.simulationBox.volume, barostat.pressure,
                                       barostat.cellMass, barostat.cellType, barostat.monoclinicAngle);
    barostat.cellVelocity = projectCellTensor(barostat.cellVelocity, barostat.cellType, barostat.monoclinicAngle);
  }

  if (thermostat)
  {
    scaling = system.thermostat->NoseHooverNVT(energies.translationalKineticEnergy, energies.rotationalKineticEnergy);
    forEachRange(engine, w,
                 [&](std::size_t, const MoleculeRange& range)
                 {
                   const RangeViews v = viewsOf(system, range);
                   Integrators::scaleVelocities(v.molecules, v.atoms, v.dynamics, components, scaling, std::nullopt, {},
                                                v.groups, {});
                 });
    // the reported kinetic energies are those of the scaled velocities, the state the chain energy belongs to
    // (as in the NVT step)
    energies.translationalKineticEnergy *= scaling.first * scaling.first;
    energies.rotationalKineticEnergy *= scaling.second * scaling.second;
    energies.NoseHooverEnergy = system.thermostat->getEnergy();
  }
  const double finalBarostatKinetic = barostat.barostatKineticEnergy();
  const double finalScale = barostat.chainStep(finalBarostatKinetic);
  barostat.logVolumeVelocity *= finalScale;
  barostat.cellVelocity = barostat.cellVelocity * finalScale;
  energies.thermobarostatEnergy = barostat.energy(system.simulationBox.volume);
  return energies;
}
}  // namespace

MolecularDynamicsSpatialDecomposition::MolecularDynamicsSpatialDecomposition() : random(std::nullopt) {};

MolecularDynamicsSpatialDecomposition::MolecularDynamicsSpatialDecomposition(InputReader& reader) noexcept
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
      engineSettings(reader.spatialDecompositionSettings),
      outputJsons(systems.size()),
      estimation(reader.numberOfBlocks, reader.numberOfProductionCycles)
{
}

void MolecularDynamicsSpatialDecomposition::run()
{
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
      throw std::runtime_error(
          "MolecularDynamicsSpatialDecomposition::run(): no resume dispatch for the checkpointed simulation stage");
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

void MolecularDynamicsSpatialDecomposition::createOutputFiles()
{
  const std::ios::openmode mode = (simulationStage != SimulationStage::Uninitialized) ? std::ios::app : std::ios::out;

  std::filesystem::create_directories("output");
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::string fileNameString =
        std::format("output/output_{}_{}.s{}.txt", system.temperature, system.input_pressure, system_id);
    streams.emplace_back(fileNameString, mode);
    ++system_id;
  }
}

void MolecularDynamicsSpatialDecomposition::checkpointIfDue(std::size_t cycle)
{
  if (cycle % writeBinaryRestartEvery == 0uz && outputToFiles)
  {
    writeBinaryRestartFile(*this);
  }

  if (GracefulShutdown::requested())
  {
    writeBinaryRestartFile(*this);
    for (std::ofstream& outputStream : streams) std::flush(outputStream);
    GracefulShutdown::exitAfterCheckpoint();
  }
}

void MolecularDynamicsSpatialDecomposition::setup()
{
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::string reason;
    if (!ForceEngine::supports(system, reason))
    {
      throw std::runtime_error(
          std::format("MolecularDynamicsSpatialDecomposition: system {} uses {}, which the spatial-decomposition "
                      "force engine does not support; use 'SimulationType': 'MolecularDynamics'\n",
                      system_id, reason));
    }
    if (system.propertyElasticConstantsFluctuation)
    {
      throw std::runtime_error(
          "MolecularDynamicsSpatialDecomposition: stress-fluctuation elastic constants are not available with the "
          "spatial-decomposition force engine; use 'SimulationType': 'MolecularDynamics'\n");
    }

    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);

    if (system_id == 0uz)
      system.containsTheFractionalMolecule = true;
    else
      system.containsTheFractionalMolecule = false;
    system.initializeGibbsSwapFractionalMoleculeGroupIds();

    if (system.forceField.interpolationScheme == ForceField::InterpolationScheme::Polynomial)
    {
      system.forceField.interpolationScheme = ForceField::InterpolationScheme::Tricubic;
    }

    ++system_id;
  }

  // Build the engines now (cell grid, sub-domain layout, FFT plans) so that a box that is too small for the
  // requested thread count or Verlet skin is reported before any MC cycles are spent.
  if (engines.size() != systems.size())
  {
    engines.clear();
    for (System& system : systems)
    {
      engines.emplace_back(engineSettings);
      engines.back().initialize(system);
    }
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

      std::print(stream, "Spatial-decomposition settings\n");
      std::print(stream, "========================================================================================================================\n");
      std::print(stream, "    number of threads (sub-domains): {}\n", engineSettings.numberOfThreads);
      std::print(stream, "    Verlet skin:                     {} [A]\n", engineSettings.verletSkin);
      std::print(stream, "    PPPM mesh spacing:               {} [A]\n", engineSettings.meshSpacing);
      std::print(stream, "    PPPM interpolation order:        {}\n", engineSettings.interpolationOrder);
      std::print(stream, "    pair kernel device:              {}\n", pairDeviceName(engineSettings.pairDevice));
      if (engineSettings.pairDevice == PairDevice::CPU)
      {
        std::print(stream, "    pair kernel precision:           {}\n",
                   pairPrecisionName(engineSettings.pairPrecision));
      }
      if (isDevice(engineSettings.pairDevice) || engineSettings.pairPrecision == PairPrecision::Mixed)
      {
        std::print(stream, "    prune skin:                      {} [A]\n", engineSettings.pruneSkin);
      }
      if (engineSettings.domainGrid.has_value())
      {
        std::print(stream, "    domain grid:                     {} x {} x {}\n", engineSettings.domainGrid->x,
                   engineSettings.domainGrid->y, engineSettings.domainGrid->z);
      }
      std::print(stream, "    MD ensemble:                     {}\n",
                 molecularDynamicsEnsembleName(system.molecularDynamicsEnsemble));
      if (system.thermobarostat.has_value())
      {
        std::print(stream, "    barostat coupling:               {} (the reported pressure is the {} estimator)\n",
                   barostatCouplingName(system.thermobarostat->coupling),
                   system.thermobarostat->coupling == BarostatCoupling::Molecular ? "molecular" : "atomic");
      }
      std::print(stream, "\n\n");

      ++system_id;
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.initializeGroupData();
    system.initializeFrameworkGroupData();
    system.precomputeTotalRigidEnergy();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                          system.spanOfFrameworkGroupData());
    system.precomputeTotalGradients();
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                    system.framework, system.spanOfFrameworkGroupData());

    if (outputToFiles)
    {
      std::ostream stream(streams[system_id].rdbuf());
      stream << system.runningEnergies.printMC("Recomputed from scratch");
      std::print(stream, "\n\n\n\n");
    }

    ++system_id;
  };
}

void MolecularDynamicsSpatialDecomposition::tearDown()
{
  if (outputToFiles)
  {
    output();
  }
}

void MolecularDynamicsSpatialDecomposition::preInitialize()
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

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    checkpointIfDue(currentCycle);

  continuePreInitializationStage:;
  }
}

void MolecularDynamicsSpatialDecomposition::initialize()
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
      std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
      System& selectedSystem = systems[selectedSystemPair.first];
      System& selectSecondSystem = systems[selectedSystemPair.second];

      std::size_t selectedComponent = selectedSystem.randomComponent(random);
      MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem);
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

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    checkpointIfDue(currentCycle);

  continueInitializationStage:;
  }
}

void MolecularDynamicsSpatialDecomposition::startEngines(std::string_view stageName)
{
  if (engines.size() != systems.size())
  {
    engines.clear();
    for (std::size_t i = 0; i < systems.size(); ++i)
    {
      engines.emplace_back(engineSettings);
    }
  }

  for (std::size_t system_id{0}; System& system : systems)
  {
    ForceEngine& engine = engines[system_id];
    if (!engine.initialized()) engine.initialize(system);

    // the exact code as the reference for this configuration; the engine leaves its own forces in the system
    const ForceEngine::Validation validation = engine.validate(system);
    system.runningEnergies = engine.computeGradients(system, true) + system.computeTailCorrectionEnergies();
    updateReportedPressure(system_id, false);

    if (outputToFiles)
    {
      std::ostream stream(streams[system_id].rdbuf());
      std::print(stream, "{}", engine.writeStatus());
      std::print(stream, "Spatial-decomposition force engine check at the start of the {} stage\n", stageName);
      std::print(stream, "========================================================================================================================\n");
      std::print(stream, "    potential energy      engine {:20.10f} exact {:20.10f} [K]  difference {:.3e}\n",
                 validation.engineEnergy, validation.referenceEnergy,
                 validation.engineEnergy - validation.referenceEnergy);
      std::print(stream, "    Ewald reciprocal      engine {:20.10f} exact {:20.10f} [K]  difference {:.3e}\n",
                 validation.engineReciprocalEnergy, validation.referenceReciprocalEnergy,
                 validation.engineReciprocalEnergy - validation.referenceReciprocalEnergy);
      std::print(
          stream,
          "    gradients             rms difference {:.3e}, maximum difference {:.3e}, rms gradient {:.3e} [K/A]\n",
          validation.rmsGradientDifference, validation.maximumGradientDifference, validation.rmsGradient);
      std::print(stream,
                 "    wall time             engine {:.4f} [s] (including the first neighbour-list build), exact "
                 "{:.4f} [s]\n",
                 validation.engineSeconds, validation.referenceSeconds);
      std::print(stream, "\n\n");
      std::flush(stream);
    }

    ++system_id;
  }
}

void MolecularDynamicsSpatialDecomposition::ensureEngines()
{
  if (engines.size() == systems.size() &&
      std::all_of(engines.begin(), engines.end(), [](const ForceEngine& engine) { return engine.initialized(); }))
  {
    return;
  }
  engines.clear();
  for (std::size_t system_id{0}; system_id < systems.size(); ++system_id)
  {
    engines.emplace_back(engineSettings);
    engines.back().initialize(systems[system_id]);
    recomputeGradients(system_id);

    if (outputToFiles && system_id < streams.size())
    {
      std::ostream stream(streams[system_id].rdbuf());
      std::print(stream, "{}    (rebuilt after the binary restart)\n\n", engines.back().writeStatus());
      std::flush(stream);
    }
  }
}

void MolecularDynamicsSpatialDecomposition::recomputeGradients(std::size_t systemId)
{
  System& system = systems[systemId];
  ForceEngine& engine = engines[systemId];
  system.runningEnergies = engine.computeGradients(system, true) + system.computeTailCorrectionEnergies();
  Integrators::updateCenterOfMassAndQuaternionGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.spanOfGroupData(), system.framework, system.spanOfFrameworkDynamics(), system.spanOfFrameworkGroupData());
  updateReportedPressure(systemId, false);
}

void MolecularDynamicsSpatialDecomposition::refreshHostState(std::size_t systemId, bool production, bool always)
{
  System& system = systems[systemId];
  ForceEngine& engine = engines[systemId];
  if (!engine.usesResident()) return;
  bool needed = always || (currentCycle % printEvery == 0uz);
  if (!needed && production)
  {
    needed = system.samplePDBMovie.has_value() || system.writeLammpsData.has_value() ||
             system.propertyConventionalRadialDistributionFunction.has_value() ||
             system.propertyMoleculeProperties.has_value() || system.propertyMoleculeShape.has_value() ||
             system.propertyMoleculeBackbone.has_value() || system.propertyMSD.has_value() ||
             system.propertyVACF.has_value() || system.propertyEndToEndACF.has_value() ||
             system.propertyDensityGrid.has_value() || system.forceBasedRDFSampleDue(currentCycle);
  }
  if (needed) engine.downloadResidentState(system);
}

void MolecularDynamicsSpatialDecomposition::updateReportedPressure(std::size_t systemId, bool accumulate)
{
  System& system = systems[systemId];
  ForceEngine& engine = engines[systemId];
  const double volume = system.simulationBox.volume;
  if (!system.thermobarostat.has_value())
  {
    system.currentExcessPressureTensor = engine.molecularPressureTensor() / volume;
    return;
  }

  // The reported pressure is the estimator the barostat drives to the external pressure. 'sampleProperties' adds
  // the molecular ideal-gas part N_molecules k T / V to the excess tensor, so the excess tensor stored here is the
  // barostat tensor minus that part: with molecular coupling simply the molecular virial over the volume, with
  // atomic coupling the atomic virial plus the ideal-gas part of the extra coupled atoms.
  const Thermobarostat& barostat = *system.thermobarostat;
  const double3x3 pressure = barostatPressureTensor(system, engine.molecularPressureTensor(), barostat.coupling,
                                                    barostat.translationalDegreesOfFreedom);
  const double molecularIdeal = static_cast<double>(system.numberOfMolecules()) / (system.beta * volume);
  system.currentExcessPressureTensor = pressure;
  system.currentExcessPressureTensor.ax -= molecularIdeal;
  system.currentExcessPressureTensor.by -= molecularIdeal;
  system.currentExcessPressureTensor.cz -= molecularIdeal;

  if (accumulate)
  {
    if (barostatPressureWindowSum.size() != systems.size())
    {
      barostatPressureWindowSum.assign(systems.size(), double3x3{});
      barostatPressureWindowCount.assign(systems.size(), 0uz);
    }
    barostatPressureWindowSum[systemId] += pressure;
    ++barostatPressureWindowCount[systemId];
  }
}

std::string MolecularDynamicsSpatialDecomposition::writeBarostatPressureWindow(std::size_t systemId)
{
  const System& system = systems[systemId];
  if (!system.thermobarostat.has_value() || systemId >= barostatPressureWindowCount.size() ||
      barostatPressureWindowCount[systemId] == 0)
    return {};

  const Thermobarostat& barostat = *system.thermobarostat;
  const double conversion = 1e-5 * Units::PressureConversionFactor;
  const double3x3 tensor =
      conversion * barostatPressureWindowSum[systemId] / static_cast<double>(barostatPressureWindowCount[systemId]);
  std::ostringstream stream;
  std::print(stream,
             "Barostat pressure tensor ({} coupling, average over the last {} steps; the estimator whose "
             "average is the external pressure):\n",
             barostatCouplingName(barostat.coupling), barostatPressureWindowCount[systemId]);
  std::print(stream, "------------------------------------------------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} [bar]\n", tensor.ax, tensor.bx, tensor.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} [bar]\n", tensor.ay, tensor.by, tensor.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} [bar]\n", tensor.az, tensor.bz, tensor.cz);
  std::print(stream, "Barostat pressure:   {: .6e} [bar]   (external pressure {: .6e} [bar])\n\n", tensor.trace() / 3.0,
             conversion * barostat.pressure);

  barostatPressureWindowSum[systemId] = double3x3{};
  barostatPressureWindowCount[systemId] = 0;
  return stream.str();
}

RunningEnergy MolecularDynamicsSpatialDecomposition::tailCorrectionEnergies(std::size_t systemId)
{
  const System& system = systems[systemId];
  if (tailCorrectionCounts.size() != systems.size()) tailCorrectionCounts.resize(systems.size());
  TailCorrectionCounts& counts = tailCorrectionCounts[systemId];
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  if (counts.numberOfAtoms != atoms.size())
  {
    const std::size_t numberOfPseudoAtomTypes = system.forceField.numberOfPseudoAtoms;
    counts.effectiveTypeCounts.assign(numberOfPseudoAtomTypes, 0.0);
    for (std::vector<double>& group : counts.groupCounts) group.assign(numberOfPseudoAtomTypes, 0.0);
    Interactions::updateEffectiveTypeCounts(counts.effectiveTypeCounts, counts.groupCounts, atoms, {});
    counts.numberOfAtoms = atoms.size();
  }
  return Interactions::computeFrameworkMoleculeTailEnergyAggregated(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(),
                                                                    counts.effectiveTypeCounts, counts.groupCounts) +
         Interactions::computeInterMolecularTailEnergyAggregated(system.forceField, system.simulationBox,
                                                                 counts.effectiveTypeCounts, counts.groupCounts);
}

RunningEnergy MolecularDynamicsSpatialDecomposition::molecularDynamicsStep(std::size_t systemId)
{
  System& system = systems[systemId];
  ForceEngine& engine = engines[systemId];
  RunningEnergy energies =
      system.thermobarostat ? engineThermobarostatVelocityVerlet(system, engine) : engineVelocityVerlet(system, engine);
  updateReportedPressure(systemId, true);
  // the engine returns the gradient-based energies; the tail corrections are added for the updated volume
  return energies + tailCorrectionEnergies(systemId);
}

void MolecularDynamicsSpatialDecomposition::equilibrate()
{
  if (simulationStage == SimulationStage::Equilibration)
  {
    ensureEngines();
    goto continueEquilibrationStage;
  }
  simulationStage = SimulationStage::Equilibration;

  for (std::size_t system_id{0}; System& system : systems)
  {
    system.initializeGroupData();
    system.initializeFrameworkGroupData();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                          system.spanOfFrameworkGroupData());
    Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                      system.spanOfMoleculeDynamics(), system.components, system.temperature,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(),
                                      &system.forceField, system.spanOfGroupData(), system.spanOfFrameworkGroupData());

    Integrators::removeCenterOfMassVelocityDrift(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    if (system.thermostat.has_value())
    {
      if (system.numberOfMolecules() > 1uz)
      {
        system.translationalCenterOfMassConstraint = 3;
        system.thermostat->translationalCenterOfMassConstraint = 3;
      }
      system.thermostat->initialize(random);
    }
    if (system.thermobarostat.has_value())
    {
      system.thermobarostat->translationalDegreesOfFreedom =
          barostatTranslationalDegreesOfFreedom(system, system.thermobarostat->coupling);
      system.thermobarostat->logVolumePosition = std::log(system.simulationBox.volume);
      system.thermobarostat->initialize(random);
    }
    ++system_id;
  }

  startEngines("equilibration");

  for (std::size_t system_id{0}; System& system : systems)
  {
    recomputeGradients(system_id);
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                    system.framework, system.spanOfFrameworkGroupData());
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
    for (std::size_t system_id{0}; System& system : systems)
    {
      system.runningEnergies = molecularDynamicsStep(system_id);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift += std::abs((system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
      refreshHostState(system_id, false);
      ++system_id;
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
          std::print(stream, "{}", writeBarostatPressureWindow(system_id));
          std::flush(stream);
        }

        ++system_id;
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
      if (outputToFiles)
      {
        writeBinaryRestartFile(*this);
      }
    }
  continueEquilibrationStage:;
  }
  // the host state of the stage end (the production stage starts from it)
  for (std::size_t system_id{0}; system_id < systems.size(); ++system_id) refreshHostState(system_id, false, true);
}

void MolecularDynamicsSpatialDecomposition::production()
{
  std::chrono::steady_clock::time_point t1, t2;

  if (simulationStage == SimulationStage::Production)
  {
    ensureEngines();
    goto continueProductionStage;
  }
  simulationStage = SimulationStage::Production;

  for (System& system : systems)
  {
    system.initializeGroupData();
    system.initializeFrameworkGroupData();
    Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                          system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                          system.spanOfFrameworkGroupData());
  }

  startEngines("production");

  for (std::size_t system_id{0}; System& system : systems)
  {
    recomputeGradients(system_id);
    system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    system.runningEnergies.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                    system.framework, system.spanOfFrameworkGroupData());
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

  numberOfSteps = 0uz;
  for (currentCycle = 0uz; currentCycle != numberOfProductionCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    t1 = std::chrono::steady_clock::now();

    estimation.setCurrentSample(currentCycle);

    for (std::size_t system_id{0}; System& system : systems)
    {
      system.runningEnergies = molecularDynamicsStep(system_id);

      system.conservedEnergy = system.runningEnergies.conservedEnergy();
      system.accumulatedDrift += std::abs((system.conservedEnergy - system.referenceEnergy) / system.referenceEnergy);
      refreshHostState(system_id, true);
      ++system_id;
    }

    // Energy decomposition for the averages. The engine's running energies are exact totals and, for a single
    // component, the complete decomposition: sampled every cycle at no cost. Only a mixture needs the exact
    // O(N^2) code for the per-pair split; that is sampled every 'PrintEvery' cycles. The pressure tensor comes
    // from the engine's virial every cycle in both cases (see molecularDynamicsStep).
    for (System& system : systems)
    {
      if (engineProvidesEnergyStatus(system))
      {
        system.currentEnergyStatus = energyStatusFromRunningEnergies(system);
        system.averageEnergies.addSample(estimation.currentBin, system.currentEnergyStatus, system.weight());
      }
      else if (currentCycle % printEvery == 0uz)
      {
        std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
        system.currentEnergyStatus = system.computeMolecularPressure().first;
        std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
        system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
        system.averageEnergies.addSample(estimation.currentBin, system.currentEnergyStatus, system.weight());
      }
    }

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    for (std::size_t system_id{0}; System& system : systems)
    {
      system.sampleProperties(system_id, estimation.currentBin, currentCycle);
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
          std::print(stream, "{}", writeBarostatPressureWindow(system_id));
          std::flush(stream);
        }
        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

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
      if (system.propertyEndToEndACF.has_value())
      {
        system.propertyEndToEndACF->writeOutput(system_id, system.components, currentCycle);
      }
      if (system.propertyMoleculeProperties.has_value())
      {
        system.propertyMoleculeProperties->writeOutput(system_id, system.components, currentCycle);
      }
      if (system.propertyMoleculeShape.has_value())
      {
        system.propertyMoleculeShape->writeOutput(system_id, system.components, currentCycle);
      }
      if (system.propertyMoleculeBackbone.has_value())
      {
        system.propertyMoleculeBackbone->writeOutput(system_id, system.components, currentCycle);
      }
      ++system_id;
    }

    if (currentCycle % printEvery == 0uz)
    {
      if (outputToFiles)
      {
        writeBinaryRestartFile(*this);
      }
    }
    t2 = std::chrono::steady_clock::now();
    totalSimulationTime += (t2 - t1);
  continueProductionStage:;
  }
  for (std::size_t system_id{0}; system_id < systems.size(); ++system_id) refreshHostState(system_id, true, true);

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
    if (system.propertyEndToEndACF.has_value())
    {
      system.propertyEndToEndACF->writeOutput(system_id, system.components, currentCycle);
    }
    if (system.propertyMoleculeProperties.has_value())
    {
      system.propertyMoleculeProperties->writeOutput(system_id, system.components, currentCycle);
    }
    if (system.propertyMoleculeShape.has_value())
    {
      system.propertyMoleculeShape->writeOutput(system_id, system.components, currentCycle);
    }
    if (system.propertyMoleculeBackbone.has_value())
    {
      system.propertyMoleculeBackbone->writeOutput(system_id, system.components, currentCycle);
    }
    ++system_id;
  }
}

void MolecularDynamicsSpatialDecomposition::output()
{
  for (std::size_t system_id{0}; System& system : systems)
  {
    std::ostream stream(streams[system_id].rdbuf());

    std::print(stream, "\n");
    std::print(stream, "========================================================================================================================\n");
    std::print(stream, "                             Simulation finished!\n");
    std::print(stream, "========================================================================================================================\n");
    std::print(stream, "\n");

    std::print(stream, "Production run CPU timings of the MD simulation\n");
    std::print(stream, "========================================================================================================================\n\n");

    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(stream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(stream, "{}", integratorsCPUTime.writeIntegratorsCPUTimeStatistics(totalSimulationTime));
    std::print(stream, "\n");
    if (system_id < engines.size())
    {
      std::print(stream, "{}", engines[system_id].writeTimings());
    }
    std::print(stream, "\n");

    std::print(
        stream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));

    std::print(stream, "Temperature averages and statistics:\n");
    std::print(stream, "========================================================================================================================\n\n");
    std::print(stream, "{}", system.averageTemperature.writeAveragesStatistics("Total"));
    std::print(stream, "{}", system.averageTranslationalTemperature.writeAveragesStatistics("Translational"));
    std::print(stream, "{}", system.averageRotationalTemperature.writeAveragesStatistics("Rotational"));

    std::print(stream, "{}", system.averagePressure.writeAveragesStatistics());

    std::print(
        stream, "{}",
        system.averageEnthalpiesOfAdsorption.writeAveragesStatistics(system.swappableComponents, system.components));
    std::print(
        stream, "{}",
        system.averagePartialMolarProperties.writeAveragesStatistics(system.swappableComponents, system.components));
    std::print(stream, "{}",
               system.averageLoadings.writeAveragesStatistics(system.components, system.frameworkMass(), std::nullopt));

    ++system_id;
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MolecularDynamicsSpatialDecomposition& md)
{
  archive << md.versionNumber;

  archive << md.outputToFiles;
  archive << md.random;

  archive << md.numberOfProductionCycles;
  archive << md.numberOfSteps;
  archive << md.numberOfPreInitializationCycles;
  archive << md.numberOfInitializationCycles;
  archive << md.numberOfEquilibrationCycles;

  archive << md.printEvery;
  archive << md.writeRestartEvery;
  archive << md.writeBinaryRestartEvery;
  archive << md.rescaleWangLandauEvery;
  archive << md.optimizeMCMovesEvery;

  archive << md.currentCycle;
  archive << md.absoluteCurrentCycle;
  archive << md.simulationStage;

  archive << md.systems;
  archive << md.fractionalMoleculeSystem;
  archive << md.engineSettings;

  archive << md.estimation;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MolecularDynamicsSpatialDecomposition& md)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > md.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'MolecularDynamicsSpatialDecomposition' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> md.outputToFiles;
  archive >> md.random;

  archive >> md.numberOfProductionCycles;
  archive >> md.numberOfSteps;
  archive >> md.numberOfPreInitializationCycles;
  archive >> md.numberOfInitializationCycles;
  archive >> md.numberOfEquilibrationCycles;

  archive >> md.printEvery;
  archive >> md.writeRestartEvery;
  archive >> md.writeBinaryRestartEvery;
  archive >> md.rescaleWangLandauEvery;
  archive >> md.optimizeMCMovesEvery;

  archive >> md.currentCycle;
  archive >> md.absoluteCurrentCycle;
  archive >> md.simulationStage;

  archive >> md.systems;
  archive >> md.fractionalMoleculeSystem;

  // The engine settings (thread count, sub-domain grid, skin, mesh) describe the hardware layout and the numerical
  // approximation, not the state of the simulation: a restarted run keeps the values of its own input file.
  SpatialDecompositionSettings archivedSettings;
  archive >> archivedSettings;

  archive >> md.estimation;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(
        std::format("MolecularDynamicsSpatialDecomposition: Invalid magic number {} at the end of the restart data\n",
                    magicNumber));
  }

  // the engines are rebuilt lazily from the settings when the interrupted MD stage resumes
  md.engines.clear();
  return archive;
}
