module;

export module interactions_cross_link;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import simulationbox;
import forcefield;
import component;
import cross_links;

/**
 * Energies, gradients and strain derivatives of the cross-links of a system (see CrossLinkTable).
 *
 * A cross-link is an inter-molecular bonded interaction between two reactive sites. Its energy is the
 * bond potential (plus the constant formation energy and the optional junction bends) minus the
 * non-bonded pair interactions that the inter-molecular sum counted for the linked pair and the 1-3
 * pairs across the link: the inter-molecular kernels are left untouched (they still sum over these
 * pairs, because the atoms belong to different molecules) and the correction here removes exactly
 * what they added and adds what an intramolecular pair would have received (the Ewald exclusion
 * -erf(alpha r)/r, or the shifted-potential completion of the finite-cutoff charge methods). The
 * net effect is that a linked pair is treated exactly like a bonded pair inside one molecule.
 *
 * The corrections are booked in the same RunningEnergy slots as the terms they cancel
 * (moleculeMoleculeVDW, moleculeMoleculeCharge, ewald_exclusion); the bonded part goes to
 * 'crossLink'. All routines return zero immediately when the table holds no link.
 */
export namespace Interactions
{
/**
 * \brief Resolves a CrossLinkSite to an index in the molecule-atom span.
 *
 * Molecules of a component are stored contiguously per component in component order, each with the
 * fixed number of atoms of its component.
 */
struct CrossLinkAtomLayout
{
  std::vector<std::size_t> componentOffset{};
  std::vector<std::size_t> atomsPerMolecule{};

  static CrossLinkAtomLayout make(const std::vector<Component> &components,
                                  const std::vector<std::size_t> &numberOfMoleculesPerComponent);

  std::size_t moleculeStart(std::size_t componentId, std::size_t moleculeIndex) const
  {
    return componentOffset[componentId] + moleculeIndex * atomsPerMolecule[componentId];
  }
  std::size_t index(const CrossLinkSite &site) const
  {
    return moleculeStart(site.componentId, site.moleculeIndex) + site.atomIndex;
  }
};

/// Total cross-link energy of the current configuration.
RunningEnergy computeCrossLinkEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                     const std::vector<Component> &components,
                                     const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                     std::span<const Atom> moleculeAtoms, const CrossLinkTable &table);

/// Energy of the given links (which need not be present in the table; their bondTypeId refers to the
/// table's bond types) in the current configuration. Used by the topology moves.
RunningEnergy computeCrossLinkEnergyOfLinks(const ForceField &forceField, const SimulationBox &simulationBox,
                                            const std::vector<Component> &components,
                                            const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                            std::span<const Atom> moleculeAtoms, const CrossLinkTable &table,
                                            std::span<const CrossLink> links);

/// Change of the cross-link energy when the atoms of molecule 'moleculeIndex' of 'componentId' move
/// from 'oldAtoms' to 'newAtoms' (both spans in the molecule's local atom order); all other atoms are
/// taken from 'moleculeAtoms'.
RunningEnergy computeCrossLinkEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                               const std::vector<Component> &components,
                                               const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                               std::span<const Atom> moleculeAtoms, const CrossLinkTable &table,
                                               std::size_t componentId, std::size_t moleculeIndex,
                                               std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms);

/// The tethers of molecule 'moleculeIndex' of 'componentId': one per link of the molecule, with the
/// partner side frozen at its current positions (see CrossLinkTether). Empty for an unlinked molecule.
std::vector<CrossLinkTether> makeCrossLinkTethers(const std::vector<Component> &components,
                                                  const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                                  std::span<const Atom> moleculeAtoms, const CrossLinkTable &table,
                                                  std::size_t componentId, std::size_t moleculeIndex);

/// The energy of the tethers 'selected' (indices into 'tethers') for the trial conformation
/// 'moleculeAtoms' of the regrown molecule (local atom order). The same terms as the table-based
/// energy of the link (bond, formation energy, junction bends, and the non-bonded exclusions across
/// the link), so the energies of a regrowth are consistent with the system's accounting.
RunningEnergy computeCrossLinkTetherEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                           const Component &component, std::span<const Atom> moleculeAtoms,
                                           std::span<const CrossLinkTether> tethers,
                                           std::span<const std::size_t> selected);

/// Total cross-link energy with the gradients added to 'moleculeDynamics' (indexed like 'moleculeAtoms').
RunningEnergy computeCrossLinkGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                       const std::vector<Component> &components,
                                       const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                       std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
                                       const CrossLinkTable &table);

/// Total cross-link energy, gradients and the strain-derivative tensor (same convention as the
/// inter-molecular strain routines: sum of gradient (x) separation).
std::pair<RunningEnergy, double3x3> computeCrossLinkEnergyStrainDerivative(
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    std::span<AtomDynamics> moleculeDynamics, const CrossLinkTable &table);
}  // namespace Interactions
