module;

module cbmc_operators;

import std;

import randomnumbers;
import units;
import forcefield;
import component;
import atom;
import double3;
import double3x3;
import cbmc_util;
import cbmc_growth_plan;
import cbmc_move_statistics;
import move_statistics;
import intra_molecular_potentials;
import chiral_center;
import bond_potential;
import bend_potential;
import torsion_potential;

// ---------------------------------------------------------------------------------------------------
// Internal Monte-Carlo move types of the flexible-bead conformation sampler; the per-bead statistics
// of these moves live in 'Component::cbmc_moves_statistics'.
// ---------------------------------------------------------------------------------------------------
namespace
{
// Sentinel bound by the base-conformation seed selection when the reservoir is empty (only while the
// reservoir is being built), which triggers the cold-start branch.
const std::vector<Atom> emptyConformation{};

enum class MoveType : std::size_t
{
  BondLengthChange = 0,
  BendAngleChange = 1,
  ConeChange = 2
};

// A bend is 'spin-variant' when it involves a placed atom other than the previous or current bead.
// Rotating the next-beads about the previous-current axis (the torsion spin) changes such bends,
// while bends among {previous, current, next-beads} transform rigidly (previous and current lie on
// the axis). Spin-variant bends are excluded from the internal bend MC (no Rosenbluth weight) and
// weighted exactly in the torsion selection instead.
bool isSpinVariantBend(const BendPotential &bend, std::optional<std::size_t> previousBead, std::size_t currentBead,
                       const std::vector<std::size_t> &nextBeads)
{
  for (std::size_t id : bend.identifiers)
  {
    if (id == currentBead) continue;
    if (previousBead.has_value() && id == previousBead.value()) continue;
    if (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end()) continue;
    return true;
  }
  return false;
}

// A torsion transforms rigidly (is spin-invariant) when all four atoms belong to the ring body
// ('currentBead' plus the 'nextBeads'); it is spin-variant when it also involves an outside atom.
bool isSpinVariantRingTorsion(const std::array<std::size_t, 4> &identifiers, std::size_t currentBead,
                              const std::vector<std::size_t> &nextBeads)
{
  bool hasRingAtom = false;
  bool hasOutsideAtom = false;
  for (std::size_t id : identifiers)
  {
    bool inRingBody = (id == currentBead) || (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end());
    bool moves = std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end();
    if (moves) hasRingAtom = true;
    if (!inRingBody) hasOutsideAtom = true;
  }
  return hasRingAtom && hasOutsideAtom;
}

// Potentials restricted to the spin-selected (junction-crossing) terms of a ring-closure step: the
// spin-variant torsions (the torsion selection filters the spin-variant bends itself, so all bends
// are kept). The internal ring torsions are sampled by the conformational MC and must be excluded
// here to avoid double counting.
Potentials::IntraMolecularPotentials ringSpinPotentials(const Potentials::IntraMolecularPotentials &intra,
                                                        std::size_t currentBead,
                                                        const std::vector<std::size_t> &nextBeads)
{
  Potentials::IntraMolecularPotentials spin{};
  spin.bends = intra.bends;
  spin.torsions.reserve(intra.torsions.size());
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads)) spin.torsions.push_back(torsion);
  }
  return spin;
}

// Signed volume of the tetrahedron of a chiral center; its sign is the center's parity.
double chiralSignedVolume(const std::array<std::size_t, 4> &ids, const std::vector<Atom> &atoms)
{
  double3 p0 = atoms[ids[0]].position;
  double3 d1 = atoms[ids[1]].position - p0;
  double3 d2 = atoms[ids[2]].position - p0;
  double3 d3 = atoms[ids[3]].position - p0;
  return double3::dot(d1, double3::cross(d2, d3));
}

// Whether swapping the positions of two branch beads would flip the parity of any chiral center that
// involves one of the swapped beads (evaluated against the pre-swap geometry, which is valid).
bool branchSwapFlipsChirality(const Component &component, const std::vector<Atom> &beforeSwap,
                              const std::vector<Atom> &afterSwap, std::size_t beadA, std::size_t beadB)
{
  for (const ChiralCenter &center : component.intraMolecularPotentials.chiralCenters)
  {
    bool involved = false;
    for (std::size_t id : center.ids)
    {
      if (id == beadA || id == beadB) involved = true;
    }
    if (!involved) continue;

    double before = chiralSignedVolume(center.ids, beforeSwap);
    double after = chiralSignedVolume(center.ids, afterSwap);
    if (before * after < 0.0) return true;
  }
  return false;
}
}  // namespace

// ---------------------------------------------------------------------------------------------------
// Coupled-decoupled torsion (spin) step: shared by all operators and both directions.
// ---------------------------------------------------------------------------------------------------
struct TorsionOrientation
{
  std::vector<Atom> positions;
  double rosenbluthWeight;
};

static TorsionOrientation selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials, double beta,
                                                   const std::vector<Atom> &chainAtoms,
                                                   const std::vector<Atom> &baseOrientation, std::size_t previousBead,
                                                   std::size_t currentBead, const std::vector<std::size_t> &nextBeads,
                                                   double3 lastBondVector,
                                                   const Potentials::IntraMolecularPotentials &intra,
                                                   bool pinFirstToBase)
{
  std::vector<std::pair<std::vector<Atom>, double>> torsion_orientations(numberOfTorsionTrials);
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  // Bends to placed atoms other than the previous bead are not invariant under the spin about the
  // previous-current axis and must enter the selection energy here (Rosenbluth-weighted identically
  // on growth and retrace).
  std::vector<const BendPotential *> spin_variant_bends{};
  for (const BendPotential &bend : intra.bends)
  {
    if (isSpinVariantBend(bend, previousBead, currentBead, nextBeads)) spin_variant_bends.push_back(&bend);
  }

  for (std::size_t j = 0; j != numberOfTorsionTrials; ++j)
  {
    double random_angle = (pinFirstToBase && j == 0) ? 0.0 : (2.0 * random.uniform() - 1.0) * std::numbers::pi;

    std::vector<Atom> rotated_atoms = baseOrientation;
    for (std::size_t k = 0; k != rotated_atoms.size(); ++k)
    {
      rotated_atoms[k].position =
          chainAtoms[currentBead].position +
          lastBondVector.rotateAroundAxis(baseOrientation[k].position - chainAtoms[currentBead].position, random_angle);
    }

    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      chain_atoms[nextBeads[k]] = rotated_atoms[k];
    }

    double torsion_energy = intra.calculateTorsionEnergies(chain_atoms);
    for (const BendPotential *bend : spin_variant_bends)
    {
      torsion_energy += bend->calculateEnergy(chain_atoms[bend->identifiers[0]].position,
                                              chain_atoms[bend->identifiers[1]].position,
                                              chain_atoms[bend->identifiers[2]].position, std::nullopt);
    }
    torsion_orientations[j] = {rotated_atoms, torsion_energy};
  }

  std::vector<double> logTorsionBoltzmannFactors{};
  logTorsionBoltzmannFactors.reserve(numberOfTorsionTrials);
  std::transform(torsion_orientations.begin(), torsion_orientations.end(),
                 std::back_inserter(logTorsionBoltzmannFactors),
                 [&](const std::pair<std::vector<Atom>, double> &v) { return -beta * std::get<1>(v); });

  double rosenbluth_weight_torsion =
      std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(), 0.0,
                      [](const double &acc, const double &logFactor) { return acc + std::exp(logFactor); });

  std::size_t selected_torsion = pinFirstToBase ? 0 : CBMC::selectTrialPosition(random, logTorsionBoltzmannFactors);

  return {torsion_orientations[selected_torsion].first,
          rosenbluth_weight_torsion / static_cast<double>(numberOfTorsionTrials)};
}

// ---------------------------------------------------------------------------------------------------
// Flexible-bead base conformation: bond lengths, bend angles, and (at branch points) the coupled
// branch arrangement sampled by an internal Metropolis MC. Carries no Rosenbluth weight. Ported from
// the former 'generateTrialOrientationsMonteCarloScheme', with chirality-protected branch swaps.
// ---------------------------------------------------------------------------------------------------
static std::vector<Atom> generateFlexibleBaseConformation(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead,
                                                          double beta, const Component &component,
                                                          const std::vector<Atom> &molecule_atoms,
                                                          std::size_t previousBead, std::size_t currentBead,
                                                          const std::vector<std::size_t> &nextBeads,
                                                          const Potentials::IntraMolecularPotentials &intra)
{
  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double3 last_bond_vector = (chain_atoms[previousBead].position - chain_atoms[currentBead].position).normalized();

  // Seed the base conformation. A cold start (empty reservoir, i.e. only while the reservoir is being
  // built) samples each next-bead's bond length from its bond potential and its direction on the bend
  // cone; otherwise an independent, well-mixed ideal-gas conformation is drawn from the component's
  // reservoir and rigidly rotated so its previous-current bond aligns with the current one. The base
  // conformation carries no Rosenbluth weight -- it is only the starting point for the internal
  // Metropolis MC below -- but the seed matters in practice: the internal MC only sees the bonded
  // energy, so a self-overlapping cold seed of a floppy molecule cannot be relaxed here and would be
  // rejected downstream by every trial direction. Drawing an independent member per grow (rather than
  // reusing one shared warm start) also decorrelates successive grows.
  const std::vector<Atom> &warmStart =
      component.conformationReservoir.empty()
          ? emptyConformation
          : component.conformationReservoir[random.uniform_integer(0, component.conformationReservoir.size() - 1)];
  if (warmStart.empty())
  {
    for (std::size_t next_bead : nextBeads)
    {
      std::optional<BondPotential> bond = intra.findBondPotential(currentBead, next_bead);
      double bond_length =
          bond.transform([&](const BondPotential &b) { return b.generateBondLength(random, beta); }).value_or(1.54);
      double angle = intra.bends.empty() ? 120.0 * Units::DegreesToRadians
                                         : intra.bends.front().generateBendAngle(random, beta);
      double3 vec = random.randomVectorOnCone(last_bond_vector, angle);
      chain_atoms[next_bead].position = chain_atoms[currentBead].position + bond_length * vec;
    }
  }
  else
  {
    double3 v = warmStart[previousBead].position - warmStart[currentBead].position;
    double3x3 rotation_matrix = double3x3::computeRotationMatrix(v, last_bond_vector);
    for (std::size_t next_bead : nextBeads)
    {
      double3 va = warmStart[next_bead].position - warmStart[currentBead].position;
      chain_atoms[next_bead].position = chain_atoms[currentBead].position + rotation_matrix * va;
    }
  }

  std::size_t number_of_trials = 2 * numberOfTrialMovesPerOpenBead * nextBeads.size();
  if (warmStart.empty()) number_of_trials *= 2;

  // Swap the positions of two branch beads (improves branch-arrangement sampling), but only when the
  // swap preserves the parity of every chiral center that involves a swapped bead.
  if (nextBeads.size() >= 2 && random.uniform() < 0.5)
  {
    std::size_t A = nextBeads[0];
    std::size_t B = nextBeads[1];
    double3 dr_A = chain_atoms[A].position - chain_atoms[currentBead].position;
    double3 dr_B = chain_atoms[B].position - chain_atoms[currentBead].position;
    double bond_length_A = dr_A.length();
    double bond_length_B = dr_B.length();
    double ratio_A = bond_length_A / bond_length_B;
    double ratio_B = bond_length_B / bond_length_A;

    std::vector<Atom> swapped(chain_atoms.begin(), chain_atoms.end());
    swapped[A].position = chain_atoms[currentBead].position + ratio_A * dr_B;
    swapped[B].position = chain_atoms[currentBead].position + ratio_B * dr_A;

    if (!branchSwapFlipsChirality(component, chain_atoms, swapped, A, B))
    {
      chain_atoms[A].position = swapped[A].position;
      chain_atoms[B].position = swapped[B].position;
    }
  }

  // With more than one bend the azimuth is constrained; Boltzmann-select it on a fine grid first.
  if (intra.bends.size() > 1)
  {
    constexpr std::size_t numberOfAzimuthAngles = 72;
    std::vector<double> logAzimuthBoltzmannFactors(numberOfAzimuthAngles);

    for (std::size_t next_bead : nextBeads)
    {
      double3 bond_vector = chain_atoms[next_bead].position - chain_atoms[currentBead].position;
      for (std::size_t r = 0; r != numberOfAzimuthAngles; ++r)
      {
        double azimuth = 2.0 * std::numbers::pi * static_cast<double>(r) / numberOfAzimuthAngles;
        chain_atoms[next_bead].position =
            chain_atoms[currentBead].position + last_bond_vector.rotateAroundAxis(bond_vector, azimuth);
        logAzimuthBoltzmannFactors[r] = -beta * intra.calculateBendSmallMCEnergies(chain_atoms);
      }
      std::size_t selected = CBMC::selectTrialPosition(random, logAzimuthBoltzmannFactors);
      double azimuth = 2.0 * std::numbers::pi * static_cast<double>(selected) / numberOfAzimuthAngles;
      chain_atoms[next_bead].position =
          chain_atoms[currentBead].position + last_bond_vector.rotateAroundAxis(bond_vector, azimuth);
    }
  }

  double current_bond_energy = intra.calculateBondSmallMCEnergies(chain_atoms);
  double current_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    std::vector<double> move_probabilities{};
    if (intra.bends.size() == 1)
    {
      move_probabilities = {0.3, 0.7, 0.0};
    }
    else
    {
      move_probabilities = trial < 20 ? std::vector<double>{0.0, 0.0, 1.0} : std::vector<double>{0.2, 0.4, 0.4};
    }

    MoveType move_type = MoveType(random.categoricalDistribution(move_probabilities));

    switch (move_type)
    {
      case MoveType::BondLengthChange:
      {
        component.cbmc_moves_statistics[currentBead].bondLengthChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bondLengthChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];
        std::optional<BondPotential> bond = intra.findBondPotential(currentBead, selected_next_bead);
        if (!bond.has_value())
        {
          throw std::runtime_error(std::format("CBMC: bond-potential can not be found (internal error)\n"));
        }

        double3 current_bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;
        double current_bond_length = current_bond_vector.length();
        double3 saved_current_position = chain_atoms[selected_next_bead].position;

        if (bond->type != BondType::Fixed)
        {
          double max_change = component.cbmc_moves_statistics[currentBead].bondLengthChange.maxChange;
          double new_bond_length = current_bond_length + (2.0 * random.uniform() - 1.0) * max_change;
          double ratio = new_bond_length / current_bond_length;
          chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + ratio * current_bond_vector;

          double new_bond_energy = intra.calculateBondSmallMCEnergies(chain_atoms);
          double energy_difference = new_bond_energy - current_bond_energy;

          component.cbmc_moves_statistics[currentBead].bondLengthChange.constructed += 1;
          component.cbmc_moves_statistics[currentBead].bondLengthChange.totalConstructed += 1;

          if (random.uniform() < ratio * ratio * std::exp(-beta * energy_difference))
          {
            current_bond_energy += energy_difference;
            component.cbmc_moves_statistics[currentBead].bondLengthChange.accepted += 1;
            component.cbmc_moves_statistics[currentBead].bondLengthChange.totalAccepted += 1;
          }
          else
          {
            chain_atoms[selected_next_bead].position = saved_current_position;
          }
        }
        break;
      }
      case MoveType::BendAngleChange:
      {
        component.cbmc_moves_statistics[currentBead].bendAngleChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];

        double current_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                              chain_atoms[selected_next_bead].position);
        double current_sinus = std::sin(current_angle);

        double3 dr = chain_atoms[currentBead].position - chain_atoms[selected_next_bead].position;
        double3 perpendicular_vector = double3::perpendicular(dr, last_bond_vector);

        double max_change = component.cbmc_moves_statistics[currentBead].bendAngleChange.maxChange;
        double random_angle = (2.0 * random.uniform() - 1.0) * max_change;

        double3 bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;
        double3 new_vec = perpendicular_vector.rotateAroundAxis(bond_vector, random_angle);
        double3 saved_current_position = chain_atoms[selected_next_bead].position;

        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + new_vec;

        double new_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                          chain_atoms[selected_next_bead].position);
        double new_sinus = std::sin(new_angle);

        double new_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);
        double energy_difference = new_bend_energy - current_bend_energy;

        component.cbmc_moves_statistics[currentBead].bendAngleChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalConstructed += 1;

        if (random.uniform() < (new_sinus / current_sinus) * std::exp(-beta * energy_difference))
        {
          current_bend_energy = new_bend_energy;
          component.cbmc_moves_statistics[currentBead].bendAngleChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].bendAngleChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
        break;
      }
      case MoveType::ConeChange:
      {
        component.cbmc_moves_statistics[currentBead].conePositionChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];
        double3 bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;

        double max_change = component.cbmc_moves_statistics[currentBead].conePositionChange.maxChange;
        double random_angle = (2.0 * random.uniform() - 1.0) * max_change;
        double3 new_vec = last_bond_vector.rotateAroundAxis(bond_vector, random_angle);
        double3 saved_current_position = chain_atoms[selected_next_bead].position;
        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + new_vec;

        double new_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);
        double energy_difference = new_bend_energy - current_bend_energy;

        component.cbmc_moves_statistics[currentBead].conePositionChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalConstructed += 1;

        if (random.uniform() < std::exp(-beta * energy_difference))
        {
          current_bend_energy = new_bend_energy;
          component.cbmc_moves_statistics[currentBead].conePositionChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].conePositionChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
        break;
      }
      default:
        std::unreachable();
    }
  }

  std::vector<Atom> next_bead_atoms(nextBeads.size());
  for (std::size_t i = 0; i != nextBeads.size(); ++i) next_bead_atoms[i] = chain_atoms[nextBeads[i]];
  return next_bead_atoms;
}

// ---------------------------------------------------------------------------------------------------
// Rigid-body tilt: samples the junction-bend tilt of a rigid fragment hinged on the anchor with a
// rigid-rotation Metropolis MC. Carries no Rosenbluth weight. Ported from the former
// 'generateRigidUnitOrientationMonteCarloScheme'.
// ---------------------------------------------------------------------------------------------------
static std::vector<Atom> generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta,
                                           const Component &component, const std::vector<Atom> &chainAtoms,
                                           std::optional<std::size_t> previousBead, std::size_t currentBead,
                                           const std::vector<std::size_t> &nextBeads,
                                           const Potentials::IntraMolecularPotentials &intra)
{
  double3 anchor_reference = component.atoms[currentBead].position;
  double3 anchor_position = chainAtoms[currentBead].position;
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  auto placeWithRotation = [&](const double3x3 &rotation)
  {
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      double3 offset = component.atoms[nextBeads[k]].position - anchor_reference;
      chain_atoms[nextBeads[k]].position = anchor_position + rotation * offset;
    }
  };

  // Seed group or no junction bends: every tilt equally likely, return a uniform orientation.
  if (!previousBead.has_value() || intra.bends.empty())
  {
    placeWithRotation(random.randomRotationMatrix());
    std::vector<Atom> result(nextBeads.size());
    for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
    return result;
  }

  double3 last_bond_vector = (chainAtoms[previousBead.value()].position - anchor_position).normalized();

  std::size_t inner = nextBeads[0];
  for (std::size_t atom : nextBeads)
  {
    if (component.connectivityTable[atom, currentBead])
    {
      inner = atom;
      break;
    }
  }

  double bend_angle = 120.0 * Units::DegreesToRadians;
  for (const BendPotential &bend : intra.bends)
  {
    if (bend.identifiers[1] != currentBead) continue;
    if ((bend.identifiers[0] == previousBead.value() && bend.identifiers[2] == inner) ||
        (bend.identifiers[0] == inner && bend.identifiers[2] == previousBead.value()))
    {
      bend_angle = bend.generateBendAngle(random, beta);
      break;
    }
  }

  double3 target_direction = random.randomVectorOnCone(last_bond_vector, bend_angle);
  double3 body_inner = (component.atoms[inner].position - anchor_reference).normalized();
  double3x3 alignment = double3x3::computeRotationMatrix(body_inner, target_direction);

  std::vector<double3> aligned_offsets(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    aligned_offsets[k] = alignment * (component.atoms[nextBeads[k]].position - anchor_reference);
  }

  constexpr std::size_t numberOfRollAngles = 72;
  std::vector<double> logRollBoltzmannFactors(numberOfRollAngles);
  double roll_offset = 2.0 * std::numbers::pi * random.uniform();
  for (std::size_t r = 0; r != numberOfRollAngles; ++r)
  {
    double roll_angle = roll_offset + 2.0 * std::numbers::pi * static_cast<double>(r) / numberOfRollAngles;
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      chain_atoms[nextBeads[k]].position =
          anchor_position + target_direction.rotateAroundAxis(aligned_offsets[k], roll_angle);
    }
    logRollBoltzmannFactors[r] = -beta * intra.calculateBendSmallMCEnergies(chain_atoms);
  }

  std::size_t selected_roll = CBMC::selectTrialPosition(random, logRollBoltzmannFactors);
  double roll_angle = roll_offset + 2.0 * std::numbers::pi * static_cast<double>(selected_roll) / numberOfRollAngles;
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    chain_atoms[nextBeads[k]].position =
        anchor_position + target_direction.rotateAroundAxis(aligned_offsets[k], roll_angle);
  }

  double current_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);

  constexpr double maximumRotationAngle = 0.15;
  std::size_t number_of_trials = 2 * numberOfTrialMovesPerOpenBead * nextBeads.size();
  std::vector<double3> saved_positions(nextBeads.size());

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    double3 axis = random.randomVectorOnUnitSphere();
    double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      saved_positions[k] = chain_atoms[nextBeads[k]].position;
      chain_atoms[nextBeads[k]].position =
          anchor_position + axis.rotateAroundAxis(saved_positions[k] - anchor_position, angle);
    }
    double trial_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);
    if (random.uniform() < std::exp(-beta * (trial_bend_energy - current_bend_energy)))
    {
      current_bend_energy = trial_bend_energy;
    }
    else
    {
      for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]].position = saved_positions[k];
    }
  }

  std::vector<Atom> result(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
  return result;
}

// ---------------------------------------------------------------------------------------------------
// Ring-closure conformation: samples the internal conformation of a cyclic cluster (kept closed by
// its closure bonds) plus the junction tilt. Carries no Rosenbluth weight. Ported from the former
// 'generateRingConformationMonteCarloScheme'. A cyclic cluster may pass through rigid-body
// fragments (e.g. a macrocycle of rigid rings linked by flexible bridges); those fragments are
// moved as whole rigid units so their internal geometry stays exact.
// ---------------------------------------------------------------------------------------------------
static bool bendIsSpinVariantRing(const std::array<std::size_t, 3> &identifiers, std::optional<std::size_t> previousBead,
                                  std::size_t currentBead, const std::vector<std::size_t> &nextBeads)
{
  for (std::size_t id : identifiers)
  {
    if (id == currentBead) continue;
    if (previousBead.has_value() && id == previousBead.value()) continue;
    if (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end()) continue;
    return true;
  }
  return false;
}

static std::vector<Atom> randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                            double3 anchorPosition)
{
  double3x3 rotation = random.randomRotationMatrix();
  std::vector<Atom> result = ringAtoms;
  for (Atom &atom : result) atom.position = anchorPosition + rotation * (atom.position - anchorPosition);
  return result;
}

static std::vector<Atom> generateRingConformation(RandomNumber &random, const ForceField &forceField, double beta,
                                                  const Component &component, const std::vector<Atom> &chainAtoms,
                                                  std::optional<std::size_t> previousBead, std::size_t currentBead,
                                                  const std::vector<std::size_t> &nextBeads,
                                                  const Potentials::IntraMolecularPotentials &intra)
{
  // Seed geometry: an independent, well-mixed ideal-gas conformation from the reservoir (so different
  // grows start from different ring puckers), or the reference geometry while the reservoir is still
  // being built. A reservoir member keeps the exact internal geometry of any rigid sub-fragment, so
  // rigid parts of the ring stay rigid. The declared chirality reference below stays 'component.atoms'.
  const std::vector<Atom> &seed =
      component.conformationReservoir.empty()
          ? component.atoms
          : component.conformationReservoir[random.uniform_integer(0, component.conformationReservoir.size() - 1)];
  double3 anchor_reference = seed[currentBead].position;
  double3 anchor_position = chainAtoms[currentBead].position;
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  auto placeWithRotation = [&](const double3x3 &rotation)
  {
    for (std::size_t atom : nextBeads)
    {
      double3 offset = seed[atom].position - anchor_reference;
      chain_atoms[atom].position = anchor_position + rotation * offset;
    }
  };

  // Chiral centers whose four atoms all have known positions during this step (the ring body, the
  // anchor, and the junction's placed neighbor) keep the parity of the reference geometry: the
  // bond/bend/torsion model is achiral (a mirror image has the same energy), so without this guard
  // the internal MC could invert a declared stereocenter, e.g. flip a cis ring fusion to trans.
  std::vector<std::pair<const std::array<std::size_t, 4> *, double>> monitored_centers{};
  {
    auto isKnown = [&](std::size_t id)
    {
      if (id == currentBead) return true;
      if (previousBead.has_value() && id == previousBead.value()) return true;
      return std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end();
    };
    for (const ChiralCenter &center : component.intraMolecularPotentials.chiralCenters)
    {
      if (std::all_of(center.ids.begin(), center.ids.end(), isKnown))
      {
        monitored_centers.push_back({&center.ids, chiralSignedVolume(center.ids, component.atoms)});
      }
    }
  }
  auto parityPreserved = [&]()
  {
    for (const auto &[ids, referenceVolume] : monitored_centers)
    {
      if (chiralSignedVolume(*ids, chain_atoms) * referenceVolume < 0.0) return false;
    }
    return true;
  };

  // A proper rotation of the reference geometry preserves the parity of every chiral center that
  // lies entirely inside the rigidly placed ring body (its signed volume keeps the reference sign),
  // so only a center that also involves the junction's placed neighbour can come out with the wrong
  // parity. Reflecting the whole placed body through a plane that contains the previous-current axis
  // flips exactly those junction-involving centers while preserving all internal distances and the
  // bend angles to the anchor -- the bonded model is achiral, so the reflected seed has identical
  // energy. Because a reflection is improper it also flips any body-internal center, so it is only a
  // valid fix when it restores every monitored parity at once; we therefore try it, keep it only when
  // 'parityPreserved()' then holds, and otherwise fall back to re-rolling random orientations. This
  // replaces an unconditional (up to 1000) random re-roll that could also fail outright, and it makes
  // the dominant single-junction-stereocentre case an O(1), deterministic correction.
  placeWithRotation(random.randomRotationMatrix());
  if (!parityPreserved() && previousBead.has_value())
  {
    double3 axis = (chain_atoms[previousBead.value()].position - anchor_position).normalized();
    double3 normal = double3::perpendicular(axis, random.randomVectorOnUnitSphere());
    std::vector<double3> beforeReflection(nextBeads.size());
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      beforeReflection[k] = chain_atoms[nextBeads[k]].position;
      double3 relative = beforeReflection[k] - anchor_position;
      chain_atoms[nextBeads[k]].position =
          anchor_position + (relative - 2.0 * double3::dot(relative, normal) * normal);
    }
    if (!parityPreserved())
    {
      for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]].position = beforeReflection[k];
    }
  }
  for (std::size_t attempt = 0; attempt != 1000 && !parityPreserved(); ++attempt)
  {
    placeWithRotation(random.randomRotationMatrix());
  }

  // Move units of the conformational MC: a single-atom (flexible) fragment moves by per-atom
  // displacement, a rigid-body fragment moves as one unit (translation or rotation about its
  // center) so its internal geometry is preserved exactly. The growth plan places a rigid fragment
  // either entirely inside this step or entirely before it, never partially.
  const FragmentGraph &graph = component.fragmentGraph;
  std::vector<std::vector<std::size_t>> moveUnits{};
  {
    std::map<std::size_t, std::size_t> rigidFragmentUnits{};
    for (std::size_t atom : nextBeads)
    {
      std::size_t fragmentIndex = graph.atomFragmentIds[atom];
      if (!graph.fragments[fragmentIndex].isRigidBody())
      {
        moveUnits.push_back({atom});
        continue;
      }
      auto [it, inserted] = rigidFragmentUnits.insert({fragmentIndex, moveUnits.size()});
      if (inserted) moveUnits.push_back({});
      moveUnits[it->second].push_back(atom);
    }
  }

  // Fixed-bond neighbours of each atom. A Fixed bond is a holonomic distance constraint that carries
  // no energy, so a free Cartesian displacement of a ring atom would stretch it with no penalty and be
  // accepted. Any move of a flexible ring atom is therefore built as a rotation about its fixed
  // neighbour(s), which preserves those bond lengths exactly (a rotation about an axis through a point
  // preserves the distance to it). 'intra.bonds' only involves atoms placed in this step, all of which
  // have positions, so every listed endpoint is usable as a pivot.
  std::vector<std::vector<std::size_t>> fixedNeighbors(chain_atoms.size());
  for (const BondPotential &bond : intra.bonds)
  {
    if (bond.type != BondType::Fixed) continue;
    fixedNeighbors[bond.identifiers[0]].push_back(bond.identifiers[1]);
    fixedNeighbors[bond.identifiers[1]].push_back(bond.identifiers[0]);
  }

  // Conformer-hopping (crankshaft) candidates: a flexible (single-atom fragment) ring atom rotated by
  // a large angle about the line through two of its positioned bonded neighbours. This preserves both
  // of those bond lengths exactly while flipping the local pucker (chair <-> twist-boat) and
  // axial <-> equatorial placement -- the barrier crossing that the small adaptive moves almost never
  // make on their own. The axis is chosen to include every Fixed-bond neighbour of the atom (fixed
  // neighbours first), so the crankshaft can never break a Fixed bond; an atom with three or more
  // Fixed bonds is over-constrained (its position is pinned, e.g. a fused-ring junction with fixed
  // bond lengths) and is left to a concerted move, not offered here. The proposal is symmetric
  // (uniform +/- angle), so plain Metropolis on the base energy keeps the sampled distribution exact.
  struct CrankshaftCandidate
  {
    std::size_t atom;
    std::size_t axisA;
    std::size_t axisB;
  };
  std::vector<CrankshaftCandidate> crankshafts{};
  {
    std::vector<bool> positioned(chain_atoms.size(), false);
    positioned[currentBead] = true;
    if (previousBead.has_value()) positioned[previousBead.value()] = true;
    for (std::size_t atom : nextBeads) positioned[atom] = true;
    for (std::size_t atom : nextBeads)
    {
      if (graph.fragments[graph.atomFragmentIds[atom]].isRigidBody()) continue;
      if (fixedNeighbors[atom].size() >= 3) continue;

      // Fixed neighbours first (they must lie on the axis), then any other positioned bonded
      // neighbours; the first two form the crankshaft axis.
      std::vector<std::size_t> axisNeighbors = fixedNeighbors[atom];
      for (std::size_t other = 0; other != chain_atoms.size(); ++other)
      {
        if (other == atom || !positioned[other]) continue;
        if (!component.connectivityTable[atom, other]) continue;
        if (std::find(fixedNeighbors[atom].begin(), fixedNeighbors[atom].end(), other) != fixedNeighbors[atom].end())
          continue;
        axisNeighbors.push_back(other);
      }
      if (axisNeighbors.size() >= 2) crankshafts.push_back({atom, axisNeighbors[0], axisNeighbors[1]});
    }
  }

  Potentials::IntraMolecularPotentials baseIntra{};
  baseIntra.bonds = intra.bonds;
  for (const BendPotential &bend : intra.bends)
  {
    if (!bendIsSpinVariantRing(bend.identifiers, previousBead, currentBead, nextBeads)) baseIntra.bends.push_back(bend);
  }
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (!isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads)) baseIntra.torsions.push_back(torsion);
  }

  auto baseEnergy = [&](std::vector<Atom> &atoms)
  {
    return baseIntra.calculateBondSmallMCEnergies(atoms) + baseIntra.calculateBendSmallMCEnergies(atoms) +
           baseIntra.calculateTorsionEnergies(atoms);
  };

  double current_energy = baseEnergy(chain_atoms);

  // Adaptive internal-MC step sizes: the maximum displacement and rotation angle are read from the
  // anchor bead's CBMC statistics and adapted towards the target acceptance ratio between sweeps by
  // 'System::optimizeMCMoves'. These moves carry no Rosenbluth weight, so their step size affects
  // only sampling efficiency, not detailed balance.
  MoveStatistics<double> &displacementStats = component.cbmc_moves_statistics[currentBead].ringDisplacementChange;
  MoveStatistics<double> &rotationStats = component.cbmc_moves_statistics[currentBead].ringRotationChange;
  MoveStatistics<double> &crankshaftStats = component.cbmc_moves_statistics[currentBead].ringCrankshaftMove;
  const double maximumDisplacement = displacementStats.maxChange;
  const double maximumRotationAngle = rotationStats.maxChange;
  const bool haveJunction = previousBead.has_value();
  std::size_t number_of_trials = 2 * forceField.numberOfTrialMovesPerOpenBead * nextBeads.size();
  std::vector<double3> saved(nextBeads.size());

  auto acceptOrReject = [&](MoveStatistics<double> &stats, std::size_t unitSize, auto restore)
  {
    stats.counts += 1.0;
    stats.totalCounts += 1.0;
    stats.constructed += 1.0;
    stats.totalConstructed += 1.0;
    double trial_energy = baseEnergy(chain_atoms);
    if (parityPreserved() && random.uniform() < std::exp(-beta * (trial_energy - current_energy)))
    {
      current_energy = trial_energy;
      stats.accepted += 1.0;
      stats.totalAccepted += 1.0;
    }
    else
    {
      restore(unitSize);
    }
  };

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    // Conformer-hopping crankshaft: a large-angle rotation of one flexible ring atom about the line
    // through two of its neighbours (see 'crankshafts' above). Attempted a fraction of the time
    // ('CBMCRingCrankshaftProbability') so local relaxation still dominates; it supplies the barrier
    // crossings between ring conformers. Tracked in its own statistics: the angle is deliberately
    // full-range and never adapted, and pooling its acceptances into the adaptive rotation statistics
    // would distort that step-size optimization.
    if (!crankshafts.empty() && random.uniform() < forceField.cbmcRingCrankshaftProbability)
    {
      const CrankshaftCandidate &c = crankshafts[random.uniform_integer(0, crankshafts.size() - 1)];
      double3 pivot = chain_atoms[c.axisA].position;
      double3 axis = (chain_atoms[c.axisB].position - pivot).normalized();
      double angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
      double3 savedPosition = chain_atoms[c.atom].position;
      chain_atoms[c.atom].position = pivot + axis.rotateAroundAxis(savedPosition - pivot, angle);
      acceptOrReject(crankshaftStats, 1, [&](std::size_t) { chain_atoms[c.atom].position = savedPosition; });
      continue;
    }

    if (haveJunction && random.uniform() < forceField.cbmcRingTiltProbability)
    {
      double3 axis = random.randomVectorOnUnitSphere();
      double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
      for (std::size_t k = 0; k != nextBeads.size(); ++k)
      {
        saved[k] = chain_atoms[nextBeads[k]].position;
        chain_atoms[nextBeads[k]].position = anchor_position + axis.rotateAroundAxis(saved[k] - anchor_position, angle);
      }
      acceptOrReject(rotationStats, nextBeads.size(),
                     [&](std::size_t n)
                     {
                       for (std::size_t k = 0; k != n; ++k) chain_atoms[nextBeads[k]].position = saved[k];
                     });
    }
    else
    {
      const std::vector<std::size_t> &unit = moveUnits[random.uniform_integer(0, moveUnits.size() - 1)];
      for (std::size_t k = 0; k != unit.size(); ++k) saved[k] = chain_atoms[unit[k]].position;
      auto restore = [&](std::size_t n)
      { for (std::size_t k = 0; k != n; ++k) chain_atoms[unit[k]].position = saved[k]; };

      if (unit.size() == 1)
      {
        // Local single-atom move, constraint-preserving. With no Fixed bond the atom is displaced
        // freely; with Fixed bonds it is rotated about its fixed neighbour(s) so those exact lengths
        // are kept: one fixed neighbour leaves a full sphere (rotate about a random axis through it),
        // two leave a circle (rotate about the line through both), and three or more pin the atom, so
        // it moves only through the whole-cluster rotation or a concerted move.
        const std::vector<std::size_t> &fixed = fixedNeighbors[unit[0]];
        if (fixed.empty())
        {
          double3 displacement{(2.0 * random.uniform() - 1.0) * maximumDisplacement,
                               (2.0 * random.uniform() - 1.0) * maximumDisplacement,
                               (2.0 * random.uniform() - 1.0) * maximumDisplacement};
          chain_atoms[unit[0]].position += displacement;
          acceptOrReject(displacementStats, unit.size(), restore);
        }
        else if (fixed.size() <= 2)
        {
          double3 pivot = chain_atoms[fixed[0]].position;
          double3 axis = fixed.size() == 1 ? random.randomVectorOnUnitSphere()
                                           : (chain_atoms[fixed[1]].position - pivot).normalized();
          double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
          chain_atoms[unit[0]].position = pivot + axis.rotateAroundAxis(chain_atoms[unit[0]].position - pivot, angle);
          acceptOrReject(rotationStats, unit.size(), restore);
        }
      }
      else if (random.uniform() < 0.5)
      {
        // Rigid fragment: symmetric whole-unit translation.
        double3 displacement{(2.0 * random.uniform() - 1.0) * maximumDisplacement,
                             (2.0 * random.uniform() - 1.0) * maximumDisplacement,
                             (2.0 * random.uniform() - 1.0) * maximumDisplacement};
        for (std::size_t atom : unit) chain_atoms[atom].position += displacement;
        acceptOrReject(displacementStats, unit.size(), restore);
      }
      else
      {
        // Rigid fragment: symmetric rotation about the unit center.
        double3 center{};
        for (std::size_t atom : unit) center += chain_atoms[atom].position;
        center = center / static_cast<double>(unit.size());
        double3 axis = random.randomVectorOnUnitSphere();
        double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
        for (std::size_t atom : unit)
        {
          chain_atoms[atom].position = center + axis.rotateAroundAxis(chain_atoms[atom].position - center, angle);
        }
        acceptOrReject(rotationStats, unit.size(), restore);
      }
    }
  }

  std::vector<Atom> result(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
  return result;
}

// ---------------------------------------------------------------------------------------------------
// Base-conformation dispatch shared by grow / retrace / recoil.
// ---------------------------------------------------------------------------------------------------
static std::vector<Atom> sampleBaseConformation(RandomNumber &random, const ForceField &forceField, double beta,
                                                const Component &component, const std::vector<Atom> &chainAtoms,
                                                const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return generateRingConformation(random, forceField, beta, component, chainAtoms,
                                    step.previousBead, step.currentBead, step.nextBeads, step.intra);
  }
  if (step.rigidBody)
  {
    return generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                             step.previousBead, step.currentBead, step.nextBeads, step.intra);
  }
  return generateFlexibleBaseConformation(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                          step.previousBead.value(), step.currentBead, step.nextBeads, step.intra);
}

// The intramolecular potentials used in the torsion (spin) selection: the junction-crossing subset
// for a ring-closure step, the full step potentials otherwise.
static Potentials::IntraMolecularPotentials torsionSelectionPotentials(const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return ringSpinPotentials(step.intra, step.currentBead, step.nextBeads);
  }
  return step.intra;
}

// ---------------------------------------------------------------------------------------------------
// Public operators.
// ---------------------------------------------------------------------------------------------------
std::vector<CBMC::StepTrial> CBMC::generateGrowTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                                      const Component &component, const std::vector<Atom> &chainAtoms,
                                                      const GrowStep &step, std::size_t numberOfTrialDirections)
{
  std::vector<StepTrial> trials(numberOfTrialDirections);

  // Seed step: no orientational reference exists yet, every orientation is equally likely.
  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      // Rigid seed: each direction is an independent uniform orientation of the body about the anchor.
      for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
                     1.0};
      }
      return trials;
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      // Ring seed: sample one internal conformation and rigidly rotate it for the other directions.
      std::vector<Atom> base =
          generateRingConformation(random, forceField, beta, component, chainAtoms,
                                   std::nullopt, step.currentBead, step.nextBeads, step.intra);
      trials[0] = {base, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, base, chainAtoms[step.currentBead].position), 1.0};
      }
      return trials;
    }

    // Flexible single-bond seed: Boltzmann bond length in a uniformly random direction.
    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
    {
      double bond_length = bond ? bond->generateBondLength(random, beta) : 1.54;
      double3 unit_vector = random.randomVectorOnUnitSphere();
      Atom trial_atom = chainAtoms[step.nextBeads[0]];
      trial_atom.position = chainAtoms[step.currentBead].position + bond_length * unit_vector;
      trials[i] = {{trial_atom}, 1.0};
    }
    return trials;
  }

  // Attach / ring-closure with a junction: one shared base conformation, one torsion spin per
  // direction about the junction bond.
  std::vector<Atom> base = sampleBaseConformation(random, forceField, beta, component, chainAtoms, step);
  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  double3 last_bond_vector =
      (chainAtoms[step.previousBead.value()].position - chainAtoms[step.currentBead].position).normalized();

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion =
        selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, chainAtoms, base,
                                 step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                                 torsionIntra, false);
    trials[i] = {torsion.positions, torsion.rosenbluthWeight};
  }
  return trials;
}

std::vector<CBMC::StepTrial> CBMC::generateRetraceTrials(RandomNumber &random, const ForceField &forceField,
                                                         double beta, const Component &component,
                                                         const std::vector<Atom> &chainAtoms, const GrowStep &step,
                                                         std::size_t numberOfTrialDirections)
{
  std::vector<StepTrial> trials(numberOfTrialDirections);

  // The old positions of the step's next-beads (trial direction 0).
  std::vector<Atom> old_orientation(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) old_orientation[k] = chainAtoms[step.nextBeads[k]];

  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      trials[0] = {old_orientation, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
                     1.0};
      }
      return trials;
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      trials[0] = {old_orientation, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, old_orientation, chainAtoms[step.currentBead].position), 1.0};
      }
      return trials;
    }

    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    trials[0] = {old_orientation, 1.0};
    for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
    {
      double bond_length = bond ? bond->generateBondLength(random, beta) : 1.54;
      double3 unit_vector = random.randomVectorOnUnitSphere();
      Atom trial_atom = chainAtoms[step.nextBeads[0]];
      trial_atom.position = chainAtoms[step.currentBead].position + bond_length * unit_vector;
      trials[i] = {{trial_atom}, 1.0};
    }
    return trials;
  }

  // Attach / ring-closure with a junction. The old orientation is the shared torsion base for every
  // trial direction, pinned as torsion trial 0 of the first trial direction.
  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  double3 last_bond_vector =
      (chainAtoms[step.previousBead.value()].position - chainAtoms[step.currentBead].position).normalized();

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion =
        selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, chainAtoms, old_orientation,
                                 step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                                 torsionIntra, i == 0);
    trials[i] = {i == 0 ? old_orientation : torsion.positions, torsion.rosenbluthWeight};
  }
  return trials;
}

CBMC::StepTrial CBMC::generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, const std::vector<Atom> &contextAtoms,
                                          const GrowStep &step)
{
  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      return {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, contextAtoms,
                                std::nullopt, step.currentBead, step.nextBeads, step.intra),
              1.0};
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      return {generateRingConformation(random, forceField, beta, component, contextAtoms,
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
              1.0};
    }
    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    double bond_length = bond ? bond->generateBondLength(random, beta) : 1.54;
    double3 unit_vector = random.randomVectorOnUnitSphere();
    Atom trial_atom = contextAtoms[step.nextBeads[0]];
    trial_atom.position = contextAtoms[step.currentBead].position + bond_length * unit_vector;
    return {{trial_atom}, 1.0};
  }

  std::vector<Atom> base = sampleBaseConformation(random, forceField, beta, component, contextAtoms, step);
  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  double3 last_bond_vector =
      (contextAtoms[step.previousBead.value()].position - contextAtoms[step.currentBead].position).normalized();

  TorsionOrientation torsion =
      selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, contextAtoms, base,
                               step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                               torsionIntra, false);
  return {torsion.positions, torsion.rosenbluthWeight};
}

double CBMC::oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                           const std::vector<Atom> &oldAtoms, const GrowStep &step)
{
  if (!step.previousBead.has_value()) return 1.0;

  std::vector<Atom> old_orientation(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) old_orientation[k] = oldAtoms[step.nextBeads[k]];

  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  double3 last_bond_vector =
      (oldAtoms[step.previousBead.value()].position - oldAtoms[step.currentBead].position).normalized();

  TorsionOrientation torsion =
      selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, oldAtoms, old_orientation,
                               step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                               torsionIntra, true);
  return torsion.rosenbluthWeight;
}
