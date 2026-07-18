module;

export module cbmc_generate_trialorientations_mc;

import std;

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import atom;
import simulationbox;
import cbmc_chain_data;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import framework;
import component;
import interpolation_energy_grid;
import bond_potential;
import intra_molecular_potentials;

export namespace CBMC
{
enum class MoveType : std::size_t
{
  BondLengthChange = 0,
  BendAngleChange = 1,
  ConeChange = 2
};

std::vector<Atom> generateTrialOrientationsMonteCarloScheme(
    RandomNumber &random,  std::size_t numberOfTrialMovesPerOpenBead, double beta, 
    const Component &component, const std::vector<Atom> chain_atoms,
    std::size_t previousBead, std::size_t currentBead, std::vector<std::size_t> nextBeads,
    const Potentials::IntraMolecularPotentials &intraMolecularInteractions);

/// Positions of the next-beads selected by the coupled-decoupled torsion step, together with the
/// torsion Rosenbluth weight accumulated over the torsion trial directions.
struct TorsionOrientation
{
  std::vector<Atom> positions;
  double rosenbluthWeight;
};

/**
 * \brief Coupled-decoupled base orientation of a rigid group hinged on the placed anchor.
 *
 * Rigid-body counterpart of 'generateTrialOrientationsMonteCarloScheme'. The group atoms 'nextBeads'
 * keep their body-fixed geometry, so the only degrees of freedom are the rotations of the body about
 * the anchor 'currentBead'. These decompose exactly: rotations about the previous-anchor bond leave
 * every junction bend invariant (each bend is an angle to that axis), so the bends are determined by
 * the remaining 'tilt' and the crossing torsions by the spin about the bond.
 *
 * This routine samples the tilt from the Boltzmann distribution of all junction bends with a small
 * Metropolis MC over rigid-body rotations about the anchor (initialized on a cone drawn from the
 * primary junction bend potential). Like the flexible bend MC it contributes no Rosenbluth weight.
 * The caller applies the shared torsion step ('selectTorsionOrientation') afterwards to bias the
 * spin about the junction bond by the crossing torsions.
 *
 * Without a previous bead (the group is the growth seed) or without junction bends, all tilts are
 * equally likely and a uniformly random rotation is returned.
 */
std::vector<Atom> generateRigidUnitOrientationMonteCarloScheme(
    RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta, const Component &component,
    const std::vector<Atom> &chainAtoms, std::optional<std::size_t> previousBead, std::size_t currentBead,
    const std::vector<std::size_t> &nextBeads, const Potentials::IntraMolecularPotentials &intra);

/**
 * \brief Coupled-decoupled torsion step for a single trial direction.
 *
 * Given a base orientation of the next-beads (bond lengths and bend angles already sampled), this
 * rotates it by `numberOfTorsionTrials` random angles around the last bond vector, weights each
 * rotation by its torsion Boltzmann factor, and selects one. Shared by CBMC insertion/deletion and
 * recoil growth so the torsion bias is generated in exactly one place.
 *
 * \param pinFirstToBase When true the first torsion trial is the (un-rotated) base orientation and
 *        it is the one selected -- used by the retrace, where the old configuration must be trial 0.
 */
TorsionOrientation selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials, double beta,
                                            const std::vector<Atom> &chainAtoms,
                                            const std::vector<Atom> &baseOrientation, std::size_t previousBead,
                                            std::size_t currentBead, const std::vector<std::size_t> &nextBeads,
                                            double3 lastBondVector, const Potentials::IntraMolecularPotentials &intra,
                                            bool pinFirstToBase);
}  // namespace CBMC
