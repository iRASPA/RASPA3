module;

export module integrators_update;

import std;

import molecule;
import molecule;
import atom;
import atom_dynamics;
import component;
import running_energy;
import simulationbox;
import forcefield;
import randomnumbers;
import interpolation_energy_grid;
import framework;

export namespace Integrators
{
/**
 * \brief Scales the velocities and orientation momenta of molecules.
 *
 * For rigid molecules the center-of-mass velocity and orientation momentum are scaled;
 * for flexible molecules the atomic velocities are scaled as well.
 *
 * \param moleculeData Span of molecules whose velocities are to be scaled.
 * \param moleculeAtomPositions Span of atoms corresponding to the molecules.
 * \param components Vector of component definitions.
 * \param scaling Pair of scaling factors for velocity and orientation momentum.
 * \param framework Optional framework definition.
 * \param frameworkDynamics Per-framework-atom velocities to scale for a flexible framework.
 */
void scaleVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                     std::span<AtomDynamics> moleculeDynamics, const std::vector<Component>& components,
                     std::pair<double, double> scaling, const std::optional<Framework>& framework = std::nullopt,
                     std::span<AtomDynamics> frameworkDynamics = {});

/**
 * \brief Removes total movable-system center-of-mass velocity.
 *
 * \param moleculeData Span of molecules whose velocities are to be scaled.
 * \param moleculeAtomPositions Span of atoms corresponding to the molecules.
 * \param components Vector of component definitions.
 * \param framework Optional framework definition.
 * \param frameworkAtomPositions Framework atoms used to obtain pseudo-atom masses.
 * \param frameworkDynamics Per-framework-atom velocities.
 * \param forceField Force field containing framework pseudo-atom masses.
 */
void removeCenterOfMassVelocityDrift(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                     std::span<AtomDynamics> moleculeDynamics,
                                     const std::vector<Component>& components,
                                     const std::optional<Framework>& framework = std::nullopt,
                                     std::span<const Atom> frameworkAtomPositions = {},
                                     std::span<AtomDynamics> frameworkDynamics = {},
                                     const ForceField* forceField = nullptr);
/**
 * \brief Updates the positions of molecules based on their velocities.
 *
 * For rigid molecules the center-of-mass position is propagated; for flexible molecules
 * the atomic positions are the integration variables and are propagated directly.
 *
 * \param moleculeData Span of molecules whose positions are to be updated.
 * \param moleculeAtomPositions Span of atoms corresponding to the molecules.
 * \param components Vector of component definitions.
 * \param dt Time step for the update.
 * \param framework Optional framework definition.
 * \param frameworkAtomPositions Framework atoms to propagate when flexible.
 * \param frameworkDynamics Per-framework-atom velocities.
 */
void updatePositions(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                     std::span<const AtomDynamics> moleculeDynamics, const std::vector<Component>& components,
                     double dt, const std::optional<Framework>& framework = std::nullopt,
                     std::span<Atom> frameworkAtomPositions = {},
                     std::span<const AtomDynamics> frameworkDynamics = {});

/**
 * \brief Updates the velocities and orientation momenta of molecules based on gradients.
 *
 * For rigid molecules the center-of-mass velocity and orientation momentum are updated;
 * for flexible molecules the atomic velocities are updated from the atomic gradients.
 *
 * \param moleculeData Span of molecules whose velocities are to be updated.
 * \param moleculeAtomPositions Span of atoms corresponding to the molecules.
 * \param components Vector of component definitions.
 * \param dt Time step for the update.
 * \param framework Optional framework definition.
 * \param frameworkAtomPositions Framework atoms used to obtain pseudo-atom masses.
 * \param frameworkDynamics Per-framework-atom velocities and gradients.
 * \param forceField Force field containing framework pseudo-atom masses.
 */
void updateVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                      std::span<AtomDynamics> moleculeDynamics, const std::vector<Component>& components, double dt,
                      const std::optional<Framework>& framework = std::nullopt,
                      std::span<const Atom> frameworkAtomPositions = {},
                      std::span<AtomDynamics> frameworkDynamics = {}, const ForceField* forceField = nullptr);

/**
 * \brief Initializes the velocities according to the Boltzmann distribution
 *
 * Rigid molecules get a center-of-mass velocity and an orientation momentum; flexible
 * molecules get per-atom velocities (their molecule record is synchronized to the
 * center-of-mass velocity).
 *
 * \param random A random number generator
 * \param moleculeData Span of molecules to initialize velocity for.
 * \param moleculeAtomPositions Span of atoms corresponding to the molecules.
 * \param components Components to get the inertiaVector
 * \param temperature Temperature to set velocities to.
 * \param framework Optional framework definition.
 * \param frameworkAtomPositions Framework atoms used to obtain pseudo-atom masses.
 * \param frameworkDynamics Per-framework-atom velocities to initialize.
 * \param forceField Force field containing framework pseudo-atom masses.
 */
void initializeVelocities(RandomNumber& random, std::span<Molecule> moleculeData,
                          std::span<Atom> moleculeAtomPositions, std::span<AtomDynamics> moleculeDynamics,
                          const std::vector<Component> components, double temperature,
                          const std::optional<Framework>& framework = std::nullopt,
                          std::span<const Atom> frameworkAtomPositions = {},
                          std::span<AtomDynamics> frameworkDynamics = {}, const ForceField* forceField = nullptr);

/**
 * \brief Converts molecule positions and orientations into Cartesian atom positions.
 *
 * For rigid molecules the atom positions are reconstructed from the center of mass and
 * orientation. For flexible molecules the atomic positions are authoritative and the
 * molecule center-of-mass record is synchronized from them instead.
 *
 * \param moleculeData Span of molecules.
 * \param moleculeAtomPositions Span to store the resulting atom positions.
 * \param components Vector of component definitions.
 */
void createCartesianPositions(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                              std::vector<Component> components);

/**
 * \brief Performs second-order NoSquish free rotor integration.
 *
 * \param moleculeData Span of molecules to be updated.
 * \param components Vector of component definitions.
 * \param dt Time step for the integration.
 */
void noSquishFreeRotorOrderTwo(std::span<Molecule> moleculeData, const std::vector<Component> components, double dt);

/**
 * \brief Updates center of mass velocities and orientation momenta of molecules.
 *
 * \param moleculeData Span of molecules to be updated.
 * \param moleculeAtomPositions Span of atom positions corresponding to the molecules.
 * \param components Vector of component definitions.
 */
void updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                               std::span<const AtomDynamics> moleculeDynamics,
                                               std::vector<Component> components);

/**
 * \brief Updates gradients of center of mass and orientation for molecules.
 *
 * \param moleculeData Span of molecules to be updated.
 * \param moleculeAtomPositions Span of atom positions corresponding to the molecules.
 * \param components Vector of component definitions.
 */
void updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions,
                                              std::span<const AtomDynamics> moleculeDynamics,
                                              std::vector<Component> components);

/**
 * \brief Updates gradients and computes energies due to interactions.
 *
 * \param moleculeAtomPositions Span of molecule atom positions.
 * \param frameworkAtomPositions Span of framework atom positions.
 * \param forceField Force field parameters.
 * \param simulationBox Simulation box parameters.
 * \param components Vector of component definitions.
 * \param eik_x Vector of complex exponentials in x-direction.
 * \param eik_y Vector of complex exponentials in y-direction.
 * \param eik_z Vector of complex exponentials in z-direction.
 * \param eik_xy Vector of complex exponentials in xy-plane.
 * \param fixedFrameworkStoredEik Stored complex exponentials for the fixed framework.
 * \param numberOfMoleculesPerComponent Vector of molecule counts per component.
 * \return The total running energy computed from interactions.
 */
RunningEnergy updateGradients(
    std::span<const Molecule> moleculeData, std::span<const Atom> moleculeAtomPositions,
    std::span<AtomDynamics> moleculeDynamics, std::span<const Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& totalEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent,
    const std::optional<Framework>& framework = std::nullopt,
    std::span<AtomDynamics> frameworkDynamics = {});
}  // namespace Integrators
