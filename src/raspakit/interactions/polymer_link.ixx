module;

export module interactions_polymer_link;

import std;

import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import simulationbox;
import forcefield;
import component;
import polymer_links;

/**
 * Energies, gradients and strain derivatives of the polymer bonds of a system (see PolymerLinkTable).
 *
 * The overlay consists of
 *  - the bonded terms across the polymer bonds (bond, bends, torsions), booked in 'crossLink';
 *  - the corrections of the pairs whose non-bonded treatment the new topology changes (value(wanted) -
 *    value(counted): the pair kernels keep counting these pairs as ordinary pairs and the correction removes
 *    exactly what they added and adds the excluded / scaled form), booked in the pair slots;
 *  - the overrides of the bonded terms of the monomers that contain a retyped atom ('replacement - original'),
 *    booked in the bonded slots.
 *
 * Every term is multiplied by the ramp factor of the reaction that created it. The derived term lists of the
 * table are rebuilt here when they are out of date (after a reaction or a restart), which is why the table is
 * taken by non-const reference.
 */
export namespace Interactions
{
/// Total overlay energy of the current configuration.
RunningEnergy computePolymerLinkEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                       const std::vector<Component> &components,
                                       const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                       std::span<const Atom> moleculeAtoms, PolymerLinkTable &table);

/// Total overlay energy with all gradients added to 'moleculeDynamics'. When 'strain' and 'molecularCorrection'
/// are given, the strain derivative (sum of gradient (x) separation, inter-molecular convention) and the
/// atomic-to-molecular virial correction sum (r_i - R_com) (x) g_i of the terms that act between molecules are
/// accumulated as well (the intramolecular overrides cancel in the molecular virial and are left out).
RunningEnergy computePolymerLinkGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                         const std::vector<Component> &components,
                                         const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                         std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
                                         PolymerLinkTable &table, double3x3 *strain = nullptr,
                                         double3x3 *molecularCorrection = nullptr);

/// The energy and the strain derivative for the molecular pressure (System::computeMolecularPressure): only the
/// gradients of the terms that act between molecules are added to 'pressureDynamics' (the intramolecular
/// overrides cancel in the molecular virial).
std::pair<RunningEnergy, double3x3> computePolymerLinkEnergyStrainDerivative(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    std::span<AtomDynamics> pressureDynamics, PolymerLinkTable &table);
}  // namespace Interactions
