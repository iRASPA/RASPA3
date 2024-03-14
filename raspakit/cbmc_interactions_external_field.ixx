export module cbmc_interactions_external_field;

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
    

export namespace CBMC                                                                                                   
{ 
  [[nodiscard]] std::optional<RunningEnergy> 
  computeExternalFieldEnergy(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                             double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms) noexcept;
}
