module;

module cbmc_recoil_growth;

import std;

import randomnumbers;
import units;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import forcefield;
import running_energy;
import framework;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_growth_context;
import cbmc_generate_trialorientations_mc;
import cbmc_segments;
import cbmc_ring_closure;

// A single growth step ('segment') of the chain. Shared with CBMC insertion/deletion: a set of
// 'nextBeads' grown from 'currentBead' (a whole rigid group when 'rigidUnit'), optionally with a
// 'previousBead' defining the bend/torsion reference. The growth sequence is deterministic (it only
// depends on the connectivity table and the group definition), so it is precomputed once.
using Segment = CBMC::GrowSegment;

// A generated trial direction: candidate positions of the 'nextBeads' plus the torsion Rosenbluth
// weight accumulated while selecting the torsion rotation.
struct Trial
{
  std::vector<Atom> positions;
  double torsionWeight{1.0};
};

// Bundles the (mostly constant) data needed to evaluate energies and grow feelers. The shared
// simulation environment (force field, box, grids, cut-offs, beta, ...) lives in 'env'; the recoil
// context adds the recoil-specific parameters and the precomputed growth sequence.
struct RecoilContext
{
  const CBMC::GrowContext &env;
  const Component &component;
  std::make_signed_t<std::size_t> skipBackgroundMolecule;
  std::size_t numberOfTrialDirections;  // k
  std::size_t recoilLength;             // l
  std::vector<Segment> segments;
};

// Energy of a trial placement, split into the external (non-bonded) part and the intramolecular
// van-der-Waals and Coulomb contributions with the already-placed beads. The sum of the two is
// used for opening a direction; only the external part is accumulated into the returned energies
// (the intramolecular terms are added once at the end via 'computeInternalEnergies').
struct TrialEnergy
{
  RunningEnergy external{};
  RunningEnergy intra{};

  double potentialEnergy() const { return external.potentialEnergy() + intra.potentialEnergy(); }
};

// Compute the energy of placing 'trialPositions' for the segment's next-beads.
// Returns std::nullopt when the placement has a hard overlap (a 'closed' direction).
std::optional<TrialEnergy> computeTrialEnergy(const RecoilContext &ctx, const Segment &segment,
                                              const std::vector<Atom> &contextAtoms,
                                              const std::vector<Atom> &trialPositions)
{
  std::vector<std::vector<Atom>> trialPositionSets{trialPositions};

  std::vector<CBMC::ChainTrial> external = CBMC::computeExternalNonOverlappingEnergies(
      ctx.env, ctx.component, trialPositionSets, -1, ctx.skipBackgroundMolecule);

  if (external.empty()) return std::nullopt;

  std::vector<Atom> candidate(contextAtoms.begin(), contextAtoms.end());
  for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
  {
    candidate[segment.nextBeads[k]] = trialPositions[k];
  }
  RunningEnergy intra = segment.intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(candidate);

  return TrialEnergy{external.front().energy, intra};
}

// Generate 'numberOfTrials' trial directions for a segment, sampling bond lengths and bend angles from
// their Boltzmann distributions and selecting a torsion rotation with the coupled-decoupled scheme.
std::vector<Trial> generateSegmentTrials(RandomNumber &random, const RecoilContext &ctx, const Segment &segment,
                                         const std::vector<Atom> &contextAtoms, std::size_t numberOfTrials)
{
  std::vector<Trial> trials;
  trials.reserve(numberOfTrials);

  if (segment.rigidUnit)
  {
    // Rigid group grown as a single rigid body hinged on the (already-placed) anchor, with the same
    // coupled-decoupled scheme as flexible beads: the junction bends are sampled by a small
    // Metropolis MC over the body tilt (no weight), the spin about the junction bond is biased by
    // the crossing torsions (torsion weight, shared 'selectTorsionOrientation').
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      std::vector<Atom> base_orientation = CBMC::generateRigidUnitOrientationMonteCarloScheme(
          random, ctx.env.forceField.numberOfTrialMovesPerOpenBead, ctx.env.beta, ctx.component, contextAtoms,
          segment.previousBead, segment.currentBead, segment.nextBeads, segment.intra);

      if (segment.previousBead.has_value())
      {
        double3 last_bond_vector =
            (contextAtoms[segment.previousBead.value()].position - contextAtoms[segment.currentBead].position)
                .normalized();

        CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
            random, ctx.env.forceField.numberOfTorsionTrialDirections, ctx.env.beta, contextAtoms, base_orientation,
            segment.previousBead.value(), segment.currentBead, segment.nextBeads, last_bond_vector, segment.intra,
            false);

        trials.push_back({torsion.positions, torsion.rosenbluthWeight});
      }
      else
      {
        // Seed group: no junction terms exist yet, every orientation is equally likely.
        trials.push_back({base_orientation, 1.0});
      }
    }
    return trials;
  }

  if (segment.ringUnit)
  {
    // Flexible ring grown with ring-closure CBMC, mirroring the CBMC insertion: the internal ring
    // conformation is sampled by an internal Monte-Carlo (no weight, kept closed by the closure bond)
    // and the spin about the junction bond is biased by the junction-crossing torsions and bends.
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      std::vector<Atom> base_conformation = CBMC::generateRingConformationMonteCarloScheme(
          random, ctx.env.forceField.numberOfTrialMovesPerOpenBead, ctx.env.beta, ctx.component, contextAtoms,
          segment.previousBead, segment.currentBead, segment.nextBeads, segment.intra);

      if (segment.previousBead.has_value())
      {
        Potentials::IntraMolecularPotentials spinPotentials =
            CBMC::ringSpinPotentials(segment.intra, segment.currentBead, segment.nextBeads);

        double3 last_bond_vector =
            (contextAtoms[segment.previousBead.value()].position - contextAtoms[segment.currentBead].position)
                .normalized();

        CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
            random, ctx.env.forceField.numberOfTorsionTrialDirections, ctx.env.beta, contextAtoms, base_conformation,
            segment.previousBead.value(), segment.currentBead, segment.nextBeads, last_bond_vector, spinPotentials,
            false);

        trials.push_back({torsion.positions, torsion.rosenbluthWeight});
      }
      else
      {
        // Seed ring: no junction terms exist yet, every orientation is equally likely.
        trials.push_back({base_conformation, 1.0});
      }
    }
    return trials;
  }

  if (!segment.previousBead.has_value())
  {
    // Growing a single bond with no previous bead (e.g. dimer, or starting in the middle of a chain).
    // Fallback C-C bond length (Angstrom) used only when the segment has no bond potential defined.
    constexpr double defaultBondLength = 1.54;
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      double bond_length = segment.intra.bonds.empty()
                               ? defaultBondLength
                               : segment.intra.bonds.front().generateBondLength(random, ctx.env.beta);
      double3 unit_vector = random.randomVectorOnUnitSphere();

      Atom trial_atom = contextAtoms[segment.nextBeads[0]];
      trial_atom.position = contextAtoms[segment.currentBead].position + bond_length * unit_vector;

      trials.push_back({{trial_atom}, 1.0});
    }
    return trials;
  }

  std::size_t previous_bead = segment.previousBead.value();
  std::size_t current_bead = segment.currentBead;

  double3 last_bond_vector =
      (contextAtoms[previous_bead].position - contextAtoms[current_bead].position).normalized();

  for (std::size_t i = 0; i != numberOfTrials; ++i)
  {
    // Generate the base orientation (bond lengths, bend angles and, at branch points, the coupled
    // branch-branch arrangement) with exactly the same internal Monte-Carlo scheme as CBMC, so that
    // the intramolecular Boltzmann distribution of the trial directions - and hence the ideal-gas
    // Rosenbluth weight - is identical to CBMC for linear and branched molecules alike.
    std::vector<Atom> base_orientation = CBMC::generateTrialOrientationsMonteCarloScheme(
        random, ctx.env.forceField.numberOfTrialMovesPerOpenBead, ctx.env.beta, ctx.component, contextAtoms,
        previous_bead, current_bead, segment.nextBeads, segment.intra);

    // Select a torsion rotation with the coupled-decoupled scheme (shared with CBMC insertion).
    CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
        random, ctx.env.forceField.numberOfTorsionTrialDirections, ctx.env.beta, contextAtoms, base_orientation,
        previous_bead, current_bead, segment.nextBeads, last_bond_vector, segment.intra, false);

    trials.push_back({torsion.positions, torsion.rosenbluthWeight});
  }

  return trials;
}

// Test whether an open pathway ('feeler') of 'depth' more segments can be grown starting at segment
// 'seg'. Trial positions are generated on the fly and directions are opened stochastically. Returns
// true as soon as a surviving pathway is found, or when the end of the chain is reached.
bool feelerExists(RandomNumber &random, const RecoilContext &ctx, std::size_t seg, std::size_t depth,
                  const std::vector<Atom> &contextAtoms)
{
  if (seg >= ctx.segments.size()) return true;
  if (depth == 0) return true;

  const Segment &segment = ctx.segments[seg];

  // Trial directions are generated lazily, one at a time: the feeler descends into the first open
  // direction, so directions beyond it never need to be generated (the trials are i.i.d., which
  // makes lazy generation statistically identical to generating all 'k' up front).
  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = generateSegmentTrials(random, ctx, segment, contextAtoms, 1).front();

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, segment, contextAtoms, trial.positions);
    if (!energy.has_value()) continue;  // hard overlap: closed

    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
    if (random.uniform() < open_probability)
    {
      std::vector<Atom> next_atoms(contextAtoms.begin(), contextAtoms.end());
      for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
      {
        next_atoms[segment.nextBeads[k]] = trial.positions[k];
      }

      if (feelerExists(random, ctx, seg + 1, depth - 1, next_atoms)) return true;
    }
  }

  return false;
}

// Result of the recursive recoil-growth search.
enum class GrowResult { Complete, DeadEnd, Discard };

// Bookkeeping of the successful growth step at each segment, needed afterwards for the weight.
struct GrowRecord
{
  Trial selected{};
  TrialEnergy energy{};
  double openProbability{1.0};
  std::size_t triedCount{0};  ///< Directions explored at this segment, including the selected one.
};

// Sequential recoil-growth search (Consta et al., Mol. Phys. 97 (1999) 1243, section 2.1): at each
// segment the trial directions are generated one at a time, and the chain tentatively advances into
// the first direction that is open (stochastically, with probability min(1, exp(-beta u))). When
// all k directions of a segment fail, the chain recoils one segment and explores the unused
// directions there. Recoiling is only allowed while the segment is not frozen: a segment freezes
// as soon as the chain has attained a length of 'recoilLength' segments counted from and including
// it (i.e. the chain head reached 'recoilLength - 1' segments beyond it). Recoiling into frozen
// territory discards the entire chain. Note that the viability of the selected direction is
// established by the actual continuation of the chain itself; an independent feeler check here
// would violate detailed balance.
GrowResult growRecursive(RandomNumber &random, const RecoilContext &ctx, std::size_t seg,
                         std::vector<Atom> &atoms, std::size_t &maxHead, std::vector<GrowRecord> &records)
{
  if (seg == ctx.segments.size()) return GrowResult::Complete;

  const Segment &segment = ctx.segments[seg];

  for (std::size_t j = 0; j != ctx.numberOfTrialDirections; ++j)
  {
    Trial trial = generateSegmentTrials(random, ctx, segment, atoms, 1).front();

    std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, segment, atoms, trial.positions);
    if (!energy.has_value()) continue;  // hard overlap: closed

    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
    if (random.uniform() >= open_probability) continue;  // stochastically closed

    std::vector<Atom> saved(segment.nextBeads.size());
    for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
    {
      saved[k] = atoms[segment.nextBeads[k]];
      atoms[segment.nextBeads[k]] = trial.positions[k];
    }
    maxHead = std::max(maxHead, seg);

    GrowResult result = growRecursive(random, ctx, seg + 1, atoms, maxHead, records);
    if (result == GrowResult::Complete)
    {
      records[seg] = {trial, energy.value(), open_probability, j + 1};
      return GrowResult::Complete;
    }

    for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
    {
      atoms[segment.nextBeads[k]] = saved[k];
    }

    if (result == GrowResult::Discard) return GrowResult::Discard;

    // Dead end below this direction: trying the next direction means recoiling to this segment,
    // which is forbidden once this segment is frozen.
    if (maxHead + 1 >= seg + ctx.recoilLength) return GrowResult::Discard;
  }

  return GrowResult::DeadEnd;
}

// Torsion Rosenbluth weight of the existing (old) orientation of a segment, mirroring the CBMC
// deletion scheme: the real orientation is trial 0 and (numberOfTorsionTrials - 1) random rotations
// around the last bond vector complete the torsion Rosenbluth sum.
double computeOldTorsionWeight(RandomNumber &random, const RecoilContext &ctx, const Segment &segment,
                               const std::vector<Atom> &oldAtoms)
{
  if (!segment.previousBead.has_value()) return 1.0;

  std::size_t previous_bead = segment.previousBead.value();
  std::size_t current_bead = segment.currentBead;

  // Rigid units need no special case here: their bend MC carries no weight, so the old-orientation
  // weight is the same pinned torsion selection as for flexible segments. Ring units must restrict
  // the torsion selection to the junction-crossing terms (as on growth), since their internal ring
  // torsions are sampled by the conformational Monte-Carlo and must not be weighted here.
  double3 last_bond_vector = (oldAtoms[previous_bead].position - oldAtoms[current_bead].position).normalized();

  std::vector<Atom> old_orientation(segment.nextBeads.size());
  for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
  {
    old_orientation[k] = oldAtoms[segment.nextBeads[k]];
  }

  Potentials::IntraMolecularPotentials torsionPotentials =
      segment.ringUnit ? CBMC::ringSpinPotentials(segment.intra, current_bead, segment.nextBeads) : segment.intra;

  // Mirror the CBMC deletion torsion scheme: the real orientation is trial 0 (angle 0) and the
  // remaining torsion trials are random rotations. Only the Rosenbluth weight is needed here.
  CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
      random, ctx.env.forceField.numberOfTorsionTrialDirections, ctx.env.beta, oldAtoms, old_orientation,
      previous_bead, current_bead, segment.nextBeads, last_bond_vector, torsionPotentials, true);

  return torsion.rosenbluthWeight;
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growRecoilGrowthMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule)
{
  const ForceField &forceField = context.forceField;

  RecoilContext ctx{context,
                    component,
                    skipBackgroundMolecule,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    CBMC::buildGrowSegments(component, beadsAlreadyPlaced)};

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<GrowRecord> records(ctx.segments.size());
  std::size_t maxHead = 0;

  // Grow the chain with the sequential recoil search; a dead end or a recoil into frozen territory
  // discards the chain (the trial move is rejected).
  if (growRecursive(random, ctx, 0, chain_atoms, maxHead, records) != GrowResult::Complete) return std::nullopt;

  // Weight of the grown chain, eq. (8) of Consta et al.: per segment W_i = (m_i / k) *
  // exp(-beta u_i) / p_i * w_torsion_i. Directions explored during the growth and found closed or
  // dead-ended count as unavailable (their randomness is already part of the proposal); only the
  // never-tried directions are probed with fresh feelers. The '1/k' normalization is the paper's
  // 'k^(N-1)' denominator; it makes the ideal-gas weight identical to CBMC (there m_i == k), so
  // the same IdealGasRosenbluthWeight reference is valid for both schemes.
  double chain_rosenbluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  for (std::size_t seg = 0; seg != ctx.segments.size(); ++seg)
  {
    const Segment &segment = ctx.segments[seg];
    const GrowRecord &record = records[seg];

    std::size_t numberOfFeelers = 1;
    for (std::size_t j = record.triedCount; j < ctx.numberOfTrialDirections; ++j)
    {
      Trial alternative = generateSegmentTrials(random, ctx, segment, chain_atoms, 1).front();

      std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, segment, chain_atoms, alternative.positions);
      if (!energy.has_value()) continue;  // hard overlap: closed

      double open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
      if (random.uniform() >= open_probability) continue;  // stochastically closed

      std::vector<Atom> feeler_atoms(chain_atoms.begin(), chain_atoms.end());
      for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
      {
        feeler_atoms[segment.nextBeads[k]] = alternative.positions[k];
      }

      // The feeler length 'l' includes the trial bead itself, so probe 'l - 1' segments beyond it.
      if (feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, feeler_atoms)) ++numberOfFeelers;
    }

    chain_rosenbluth_weight *= static_cast<double>(numberOfFeelers) /
                                static_cast<double>(ctx.numberOfTrialDirections) *
                                std::exp(-ctx.env.beta * record.energy.potentialEnergy()) / record.openProbability *
                                record.selected.torsionWeight;

    // Only the external part is accumulated here; the intramolecular van der Waals and Coulomb
    // terms are added once at the end through 'computeInternalEnergies'.
    chain_external_energies += record.energy.external;

    if (chain_rosenbluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  // Only bond, bend, and torsion are taken into account during the growth (intra van der Waals
  // and Coulomb enter through the selection of the beads). Correct the Rosenbluth weight with
  // the Boltzmann factor of the remaining internal interactions (Urey-Bradley, inversion-bend,
  // out-of-plane-bend, improper torsion, and the cross-terms); the returned energies already
  // contain these terms through 'computeInternalEnergies'.
  RunningEnergy unsampled_internal_energies =
      component.intraMolecularPotentials.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
  chain_rosenbluth_weight *= std::exp(-ctx.env.beta * unsampled_internal_energies.potentialEnergy());

  // Copy this configuration so that it can be used as a starting point (parity with CBMC insertion).
  component.grownAtoms = chain_atoms;

  // Build a valid molecule record (center-of-mass, mass) for the grown flexible chain
  double3 com{};
  for (std::size_t i = 0; i != chain_atoms.size(); ++i)
  {
    com += component.definedAtoms[i].second * chain_atoms[i].position;
  }
  com = com / component.totalMass;
  Molecule molecule = Molecule(com, simd_quatd(0.0, 0.0, 0.0, 1.0), component.totalMass,
                               static_cast<std::size_t>(chain_atoms.front().componentId), chain_atoms.size());

  return ChainGrowData(molecule, chain_atoms, chain_external_energies + internal_energies, chain_rosenbluth_weight,
                       0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceRecoilGrowthMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  const ForceField &forceField = context.forceField;

  RecoilContext ctx{context,
                    component,
                    -1,
                    std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
                    std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength),
                    CBMC::buildGrowSegments(component, beadsAlreadyPlaced)};

  std::vector<Atom> old_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_rosenbluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  for (std::size_t seg = 0; seg != ctx.segments.size(); ++seg)
  {
    const Segment &segment = ctx.segments[seg];

    // Energy and torsion weight of the existing (old) orientation of this segment.
    std::vector<Atom> old_positions(segment.nextBeads.size());
    for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
    {
      old_positions[k] = old_atoms[segment.nextBeads[k]];
    }

    std::optional<TrialEnergy> old_energy = computeTrialEnergy(ctx, segment, old_atoms, old_positions);
    TrialEnergy selected_energy = old_energy.value_or(TrialEnergy{});
    double selected_potential = selected_energy.potentialEnergy();
    double open_probability = std::min(1.0, std::exp(-ctx.env.beta * selected_potential));

    double torsion_weight = computeOldTorsionWeight(random, ctx, segment, old_atoms);

    // The old backbone direction is always counted as available (the old chain itself is its
    // feeler; its openness is accounted for by the p_i factors), plus any of the (k-1) additional
    // trial directions that are open and have a surviving feeler of length 'l' (which includes the
    // trial bead itself, hence 'l - 1' segments beyond it).
    std::size_t numberOfFeelers = 1;
    if (ctx.numberOfTrialDirections > 1)
    {
      std::vector<Trial> alternatives =
          generateSegmentTrials(random, ctx, segment, old_atoms, ctx.numberOfTrialDirections - 1);

      for (const Trial &trial : alternatives)
      {
        std::optional<TrialEnergy> energy = computeTrialEnergy(ctx, segment, old_atoms, trial.positions);
        if (!energy.has_value()) continue;

        double alternative_open_probability = std::min(1.0, std::exp(-ctx.env.beta * energy->potentialEnergy()));
        if (random.uniform() >= alternative_open_probability) continue;

        std::vector<Atom> next_atoms(old_atoms.begin(), old_atoms.end());
        for (std::size_t k = 0; k != segment.nextBeads.size(); ++k)
        {
          next_atoms[segment.nextBeads[k]] = trial.positions[k];
        }

        if (feelerExists(random, ctx, seg + 1, ctx.recoilLength - 1, next_atoms)) ++numberOfFeelers;
      }
    }

    // Same (m_i / k) normalization as the growth routine so that insertion and deletion weights
    // are directly comparable (and comparable to the CBMC ideal-gas reference).
    chain_rosenbluth_weight *= static_cast<double>(numberOfFeelers) /
                                static_cast<double>(ctx.numberOfTrialDirections) *
                                std::exp(-ctx.env.beta * selected_potential) / open_probability * torsion_weight;

    // Only the external part is accumulated here; the intramolecular van der Waals and Coulomb
    // terms are added once at the end through 'computeInternalEnergies'.
    chain_external_energies += selected_energy.external;
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(old_atoms);

  // Same correction as the growth routine: the Rosenbluth weight of the old chain picks up the
  // Boltzmann factor of the internal interactions that are not sampled during the retracing.
  RunningEnergy unsampled_internal_energies =
      component.intraMolecularPotentials.computeInternalEnergiesNotSampledDuringGrowth(old_atoms);
  chain_rosenbluth_weight *= std::exp(-ctx.env.beta * unsampled_internal_energies.potentialEnergy());

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosenbluth_weight, 0.0);
}
