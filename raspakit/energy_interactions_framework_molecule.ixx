module;

export module energy_interactions_framework_molecule;

import energy_status;
import potential_energy_vdw;
import potential_energy_coulomb;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import energy_status_inter;
import energy_status_intra;
import units;

import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <tuple>;
import <thread>;

export std::optional<EnergyStatus> computeFrameworkMoleculeEnergyPar(size_t numberOfComponents, const ForceField &forceField, const SimulationBox &simulationBox, std::span<Atom> &atoms, std::make_signed_t<std::size_t> skip, std::span<const Atom> &frameworkAtoms) noexcept;

//export std::optional<EnergyStatus> computeFrameworkMoleculeEnergySerial(size_t numberOfComponents, const ForceField &forceField, const SimulationBox &simulationBox, std::span<Atom> &atoms, std::make_signed_t<std::size_t> skip, std::span<const Atom> &frameworkAtoms) noexcept;
