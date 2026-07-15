module;

export module cbmc_interactions_intermolecular;

import std;

import energy_status;
import potential_correction_vdw;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_status_inter;
import running_energy;
import units;
import threadpool;

export namespace CBMC
{
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;
}
