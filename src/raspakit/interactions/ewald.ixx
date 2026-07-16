module;

export module interactions_ewald;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import energy_status;
import simulationbox;
import forcefield;
import framework;
import component;

export namespace Interactions
{
/**
 * \brief Computes the Ewald Fourier energy contribution for a single ion.
 *
 * Calculates the Fourier-space part of the Ewald summation for a single ion,
 * given its position and charge. This function is useful for computing
 * corrections due to net charges in the system or for inserting ions.
 *
 * References:
 * - An Exact Ewald Summation Method in Theory and Practice, S. Stenberg and B. Stenqvist,
 *   J. Phys. Chem. A 2020, 124, 3943−3946; https://doi.org/10.1021/acs.jpca.0c01684
 * - Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections
 *   to the Ewald potential, Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks,
 *   J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param position The position of the ion.
 * \param charge The charge of the ion.
 * \return The Fourier energy contribution for the ion.
 */
double computeEwaldFourierEnergySingleIon(std::vector<std::complex<double>> &eik_x,
                                          std::vector<std::complex<double>> &eik_y,
                                          std::vector<std::complex<double>> &eik_z,
                                          std::vector<std::complex<double>> &eik_xy, const ForceField &forceField,
                                          const SimulationBox &simulationBox, double3 position, double charge);

/**
 * \brief Precomputes the Ewald Fourier terms for a rigid framework.
 *
 * Calculates and stores the Fourier-space components for a rigid framework,
 * which can be reused in subsequent Ewald summation calculations to improve efficiency.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Storage for the precomputed Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param frameworkAtoms The atoms of the rigid framework.
 */
void precomputeEwaldFourierRigid(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms);

/**
 * \brief Computes the Ewald Fourier energy for the system.
 *
 * Calculates the Fourier-space part of the Ewald summation for the entire system,
 * including interactions between molecules and between molecules and the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Storage for the total Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param moleculeAtoms Positions and properties of the molecules' atoms.
 * \param netChargeFramework Net charge of the framework, used for the net-charge correction
 *                           (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)).
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierEnergy(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
    double netChargeFramework = 0.0);

/**
 * \brief Computes the energy difference due to atom position changes in the Ewald Fourier summation.
 *
 * Calculates the change in Fourier-space Ewald energy when atoms are moved from old positions to new positions.
 * Useful for Monte Carlo moves or molecular dynamics steps.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param storedEik Previously stored Fourier components of the system.
 * \param trialEik Updated Fourier components after the move.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param newatoms The new positions and properties of the atoms.
 * \param oldatoms The old positions and properties of the atoms.
 * \param netCharge The current total net charge of the system (framework plus adsorbates) before the
 *                  move; used for the net-charge correction when the move changes the net charge
 *                  (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)).
 * \param netChargeDerivativeExternal The per-group summed charge of dU/dlambda group-tagged atoms outside
 *                  'newatoms'/'oldatoms' (e.g. the partner fractional molecule in chained pair
 *                  moves); needed for the dU/dlambda of the net-charge correction.
 * \return The running energy containing the Ewald Fourier energy difference.
 */
RunningEnergy energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<const Atom> newatoms, std::span<const Atom> oldatoms,
    double netCharge = 0.0, const std::array<double, maximumNumberOfDUDlambdaGroups> &netChargeDerivativeExternal = {});

RunningEnergy energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms, double netCharge = 0.0,
    const std::array<double, maximumNumberOfDUDlambdaGroups> &netChargeDerivativeExternal = {});

/**
 * \brief Computes the difference in electric field due to atom position changes in the Ewald Fourier summation.
 *
 * Calculates the change in electric field when atoms are moved from old positions to new positions,
 * using the Fourier-space Ewald summation.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Previously stored Fourier components of the system.
 * \param trialEik Updated Fourier components after the move.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param electricField Output array to store the computed electric fields difference.
 * \param newatoms The new positions and properties of the atoms.
 * \param oldatoms The old positions and properties of the atoms.
 * \return The running energy containing the Ewald Fourier energy difference.
 */
RunningEnergy eletricFieldEwaldFourierEnergyDifference(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms);

void computeEwaldFourierElectricFieldDifference(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldNew, std::span<double3> electricFieldOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms);

/**
 * \brief Computes the Ewald Fourier energy and its gradient (forces) on atoms.
 *
 * Calculates the Fourier-space part of the Ewald summation and computes the forces acting on each atom.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomData Positions and properties of the atoms (updated with computed gradients).
 * \param netChargeFramework Net charge of the framework, used for the net-charge correction
 *                           (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)).
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierGradient(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> atomData,
    std::span<AtomDynamics> atomDynamics, double netChargeFramework = 0.0,
    const std::optional<Framework>& framework = std::nullopt, std::span<const Atom> frameworkAtoms = {},
    std::span<AtomDynamics> frameworkDynamics = {});

/**
 * \brief Computes the reciprocal-space (Fourier) Ewald force acting on the atoms of a single molecule.
 *
 * The reciprocal Ewald energy is a functional of the total structure factor S(k) = sum_j q_j s_j exp(i k.r_j)
 * of the whole system (rigid framework plus all molecules). The force on atom a of the selected molecule is
 *   F_a = -dU/dr_a = -prefactor sum_k A(k) 2 q_a s_a Im[ conj(S(k)) exp(i k.r_a) ] k,
 * which depends on the *complete* structure factor but only on the positions of the selected molecule's atoms.
 * This lets the force on one molecule be evaluated in O(N_k * n_selected) time by reusing a precomputed total
 * structure factor, instead of the O(N_k * N_atoms) full-system gradient.
 *
 * For a single-molecule move the caller supplies:
 *   - \p storedEik = S_old (the maintained structure factor) together with the current atom positions to obtain
 *     the force at the old configuration, or
 *   - \p storedEik = S_new (as produced by \ref energyDifferenceEwaldFourier in \c trialEik) together with the
 *     trial atom positions to obtain the force at the new configuration.
 *
 * Only the accumulated gradient stored in \p atomDynamics is updated; intramolecular self/exclusion corrections
 * are deliberately omitted because they cancel in the net (center-of-mass) force of a rigid translation.
 *
 * \param eik_x,eik_y,eik_z,eik_xy Scratch buffers for the exp(i k.r) recurrences (resized as needed).
 * \param storedEik The total structure factor S(k) of the whole system (framework + molecules).
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box.
 * \param atoms The atoms of the selected molecule at the configuration whose force is requested.
 * \param atomDynamics Per-atom gradient accumulator (size = atoms.size()).
 */
void computeEwaldFourierGradientSingleMolecule(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> atoms,
    std::span<AtomDynamics> atomDynamics);

/**
 * \brief Computes the Ewald Fourier energy and its strain derivative.
 *
 * Calculates the Fourier-space part of the Ewald summation, including the strain derivative, which is necessary
 * for pressure and stress calculations in simulations involving variable cell shapes.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Stored Fourier components of the system (unused).
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param frameworkComponents The framework components in the system.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomData Positions and properties of the atoms (updated with computed gradients).
 * \param netChargeFramework Net charge of the framework.
 * \param netChargePerComponent Net charges per component.
 * \return A pair containing the energy status and the strain derivative tensor.
 *
 * The net-charge correction (Bogusz et al., J. Chem. Phys. 108, 7070 (1998)) is computed
 * internally from the wave-vector sum, so that it remains consistent with the current
 * simulation box; its contribution to the strain derivative is included.
 */
std::pair<EnergyStatus, double3x3> computeEwaldFourierEnergyStrainDerivative(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::optional<Framework> &framework,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> atomData, std::span<AtomDynamics> atomDynamics, double netChargeFramework,
    std::vector<double> netChargePerComponent) noexcept;

/**
 * \brief Accepts a move by updating the stored Ewald Fourier components.
 *
 * Updates the stored Fourier components after a move is accepted, ensuring that future energy calculations
 * are based on the new configuration.
 *
 * \param forceField The force field parameters.
 * \param storedEik The stored Fourier components to be updated.
 * \param trialEik The new Fourier components after the move.
 */
void acceptEwaldMove(const ForceField &forceField,
                     std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
                     std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &trialEik);

/**
 * \brief Computes the electric potential at molecule positions using the Ewald Fourier summation.
 *
 * Calculates the electrostatic potential at each atom's position due to the Fourier-space Ewald summation,
 * including contributions from the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param electricPotentialMolecules Output array to store the computed electric potentials for each molecule.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param moleculeAtomPositions Positions and properties of the molecules' atoms.
 */
void computeEwaldFourierElectrostaticPotential(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    [[maybe_unused]] std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik,
    std::span<double> electricPotentialMolecules, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> moleculeAtomPositions);

/**
 * \brief Computes the electric field at molecule positions using the Ewald Fourier summation.
 *
 * Calculates the electrostatic field at each atom's position due to the Fourier-space Ewald summation,
 * including contributions from the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Storage for the total Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param electricFieldMolecules Output array to store the computed electric fields for each molecule.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomData Positions and properties of the atoms (updated with computed fields).
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierElectricField(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldMolecules,
    const std::vector<Component> &components, const std::vector<std::size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> atomData);

/**
 * \brief Computes the periodic Coulomb potential matrix used in charge equilibration.
 *
 * Fills the N x N matrix V (row-major) with the Ewald lattice sum of 1/|r_i - r_j + L| over all
 * periodic images L (with neutralizing background). The diagonal contains the interaction of a
 * site with its own periodic images (self/Wigner potential). The matrix is in units of
 * [1/Angstrom]; no Coulomb conversion factor is applied.
 *
 * The Ewald parameters (alpha and the number of wave vectors) are derived internally from the
 * simulation box, because charge equilibration runs during CIF-reading, before the force-field
 * Ewald parameters have been initialized. The real-space part uses the minimum-image convention
 * with a cutoff of half the smallest perpendicular box width; the Fourier part uses the same
 * exp(ik.r) recursion as the other Ewald routines and accumulates each wave vector as a
 * symmetric rank-2 update of the matrix.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param simulationBox The simulation box parameters.
 * \param atoms The atoms (only positions are used).
 * \param potentialMatrix Output buffer of size atoms.size() * atoms.size(), row-major.
 */
void computeEwaldFourierChargeEquilibrationPotentialMatrix(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    const SimulationBox &simulationBox, std::span<const Atom> atoms, std::span<double> potentialMatrix);

}  // namespace Interactions
