module;

export module spatial_decomposition_integrator;

import std;

import double3;
import double3x3;
import running_energy;
import atom;
import atom_dynamics;
import molecule;
import component;
import system;
import force_engine;

/**
 * \brief Host-side integrator passes over molecule ranges, and the velocity-Verlet step with the engine forces.
 *
 * Every pass of the velocity-Verlet step over the molecules (velocity and position updates, kinetic energies, the
 * barostat coupling of the velocities and positions) is independent per molecule, so the molecules are cut into
 * contiguous ranges balanced on the atom count and the passes run on the engine's worker team, one range per
 * thread. The Integrators functions index the atoms and rigid-group states relative to the spans they are given,
 * so they run unchanged on a range. Consecutive per-molecule passes between two global reductions are fused into
 * one team task; the partial sums of the ranges are reduced in range order, so the results do not depend on the
 * scheduling. The framework terms of the Integrators functions are not needed (the engine rejects frameworks).
 *
 * The NVE / NVT step ('engineVelocityVerlet') is shared by the spatial-decomposition MD driver and the hybrid
 * Monte Carlo move; the thermobarostat step lives with the driver.
 */
export namespace SpatialDecompositionIntegrator
{
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

RangeViews viewsOf(System& system, const MoleculeRange& range);

std::size_t rigidGroupsOf(const Component& component);

/// Cuts the molecules into `parts` contiguous ranges balanced on the atom count. The atoms of the molecules are
/// stored consecutively in molecule order, the rigid-group states likewise (one per rigid fragment of a
/// semi-flexible component), which is what the Integrators functions assume as well.
std::vector<MoleculeRange> moleculeRanges(const System& system, std::size_t parts);

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
void centersOfMass(const System& system, const MoleculeRange& range, bool withPositions, CenterOfMassTable& table);

/// Kinetic virial sum_k M_k V_k V_k^T of the molecules of the range from the table (molecular coupling).
double3x3 tableKineticVirial(const MoleculeRange& range, const CenterOfMassTable& table);

/// Translational and rotational kinetic energy of the molecules of a range.
std::pair<double, double> rangeKineticEnergies(System& system, const MoleculeRange& range);

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

/// Velocity Verlet with the engine forces (Integrators::velocityVerlet with updateGradients replaced): one step
/// of `system.timeStep`, with the Nose-Hoover thermostat of the system when it has one (NVE otherwise). The
/// returned energies are those of the engine (pair, reciprocal, bonded terms; no tail correction) plus the kinetic
/// energies after the step.
RunningEnergy engineVelocityVerlet(System& system, ForceEngine& engine);
}  // namespace SpatialDecompositionIntegrator
