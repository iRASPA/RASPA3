module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>
#endif

export module cbmc_interactions_framework_molecule;

#ifndef USE_LEGACY_HEADERS
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
#endif

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
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) noexcept;
}
