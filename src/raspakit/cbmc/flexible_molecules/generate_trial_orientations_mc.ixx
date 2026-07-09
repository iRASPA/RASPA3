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
