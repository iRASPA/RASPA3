module;

module cbmc_interactions_external_field;

import <numbers>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>; 
import <cmath>;
import <optional>;
import <thread>;
import <future>;
import <type_traits>;

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
import force_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;


[[nodiscard]] std::optional<RunningEnergy> 
CBMC::computeExternalFieldEnergy([[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox, 
                                 [[maybe_unused]] std::span<const Atom> moleculeAtoms, [[maybe_unused]] double cutOffVDW, [[maybe_unused]] double cutOffCoulomb, 
                                 [[maybe_unused]] std::span<Atom> atoms, [[maybe_unused]] std::make_signed_t<std::size_t> skip) noexcept
{
  RunningEnergy energySum;

  return energySum;
}
