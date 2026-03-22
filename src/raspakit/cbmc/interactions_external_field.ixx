module;

export module cbmc_interactions_external_field;

import std;

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
import interpolation_energy_grid;

export namespace CBMC
{
[[nodiscard]] std::optional<RunningEnergy> computeExternalFieldEnergy(bool hasExternalField,
       const ForceField &forceField,
       const SimulationBox &simulationBox,
       const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
       double cutOffVDW, double cutOffCoulomb,
       std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) noexcept;
}
