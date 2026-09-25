module;

export module cbmc_intermolecular;

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
/// Inter-molecular energy of the trial atoms 'atoms' against the background 'moleculeAtoms' (pairs
/// within the same molecule id are skipped, as are all background atoms of 'skipBackgroundMolecule'
/// when given -- the molecule being regrown). std::nullopt on a hard overlap.
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;
}
