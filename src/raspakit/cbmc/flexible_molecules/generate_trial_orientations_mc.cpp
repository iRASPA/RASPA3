module;

module cbmc_generate_trialorientations_mc;

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
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import forcefield;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import bond_potential;
import bend_potential;
import intra_molecular_potentials;
import cbmc_move_statistics;
import connectivity_table;

// A bend is 'spin-variant' when it involves a placed atom other than the previous or current bead.
// Rotating the next-beads about the previous-current axis (the torsion spin) changes such bends --
// e.g. the bend to the second ring neighbor when growing off an aromatic ring atom -- while bends
// among {previous, current, next-beads} transform rigidly (previous and current lie on the axis).
// Spin-variant bends are therefore excluded from the internal bend MC (which carries no Rosenbluth
// weight) and weighted exactly in the torsion selection instead.
static bool isSpinVariantBend(const BendPotential &bend, std::optional<std::size_t> previousBead,
                              std::size_t currentBead, const std::vector<std::size_t> &nextBeads)
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

std::vector<Atom> CBMC::generateTrialOrientationsMonteCarloScheme(
    RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta, 
    const Component &component, const std::vector<Atom> molecule_atoms,
    std::size_t previousBead, std::size_t currentBead, std::vector<std::size_t> nextBeads,
    const Potentials::IntraMolecularPotentials &intraMolecularInteractions)
{
  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double3 last_bond_vector = (chain_atoms[previousBead].position - chain_atoms[currentBead].position).normalized();

  // Calculate initial positions for the next beads
  if (component.grownAtoms.empty())
  {
    // Generate positions from scratch
    for (std::size_t next_bead : nextBeads)
    {
      std::optional<BondPotential> bond = intraMolecularInteractions.findBondPotential(currentBead, next_bead);

      double bond_length =
          bond.transform([&random, &beta](const BondPotential &bond) { return bond.generateBondLength(random, beta); })
              .value_or(1.54);
      double angle = intraMolecularInteractions.bends.empty()
                         ? 120.0 * Units::DegreesToRadians
                         : intraMolecularInteractions.bends.front().generateBendAngle(random, beta);

      double3 vec = random.randomVectorOnCone(last_bond_vector, angle);

      chain_atoms[next_bead].position = chain_atoms[currentBead].position + bond_length * vec;
    }
  }
  else
  {
    // Restore from previous configuration stored in component
    double3 v = component.grownAtoms[previousBead].position - component.grownAtoms[currentBead].position;

    double3x3 rotation_matrix = double3x3::computeRotationMatrix(v, last_bond_vector);

    for (std::size_t next_bead : nextBeads)
    {
      double3 va = component.grownAtoms[next_bead].position - component.grownAtoms[currentBead].position;

      double3 vec = rotation_matrix * va;

      chain_atoms[next_bead].position = chain_atoms[currentBead].position + vec;
    }
  }

  std::size_t number_of_trials = 2 * numberOfTrialMovesPerOpenBead * nextBeads.size();
  if(component.grownAtoms.empty())
  {
    number_of_trials *= 2;
  }

  // Randomly swap the positions of two beads (only when there is no chirality)
  // This is needed especially when starting from an previous stored configuration
  // because otherwise the arrangement of the branches would hardly change.
  if(nextBeads.size() >= 2)
  {
    if(random.uniform() < 0.5)
    {
      std::size_t A = nextBeads[0];
      std::size_t B = nextBeads[1];
      double3 dr_A = chain_atoms[A].position - chain_atoms[currentBead].position;
      double3 dr_B = chain_atoms[B].position - chain_atoms[currentBead].position;
      double bond_length_A = dr_A.length();
      double bond_length_B = dr_B.length();
      double ratio_A = bond_length_A / bond_length_B;
      double ratio_B = bond_length_B / bond_length_A;

      chain_atoms[A].position = chain_atoms[currentBead].position + ratio_A * dr_B;
      chain_atoms[B].position = chain_atoms[currentBead].position + ratio_B * dr_A;
    }
  }

  // With more than one bend the azimuth of a bead on its cone around the last bond vector is
  // constrained as well: e.g. a CH2 growing off an aromatic ring atom has bends to both ring
  // neighbors, which together force it towards the ring plane. Neither initialization above fixes
  // the azimuth: the from-scratch cone draw samples it uniformly, and the restore from the stored
  // configuration preserves only the bend to the previous bead (the positions of the other placed
  // neighbors differ from the stored ones). An azimuth that starts far off against stiff bends
  // cannot reliably be recovered by the local Metropolis refinement below, so Boltzmann-select each
  // bead's azimuth from the total bend energy on a fine grid first (a Gibbs sweep over the beads);
  // the Metropolis walk afterwards only needs to sample local fluctuations.
  if (intraMolecularInteractions.bends.size() > 1)
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
        logAzimuthBoltzmannFactors[r] = -beta * intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);
      }

      std::size_t selected = CBMC::selectTrialPosition(random, logAzimuthBoltzmannFactors);
      double azimuth = 2.0 * std::numbers::pi * static_cast<double>(selected) / numberOfAzimuthAngles;
      chain_atoms[next_bead].position =
          chain_atoms[currentBead].position + last_bond_vector.rotateAroundAxis(bond_vector, azimuth);
    }
  }

  double current_bond_energy = intraMolecularInteractions.calculateBondSmallMCEnergies(chain_atoms);
  double current_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    std::vector<double> move_probabilities{};
    if (intraMolecularInteractions.bends.size() == 1)
    {
      move_probabilities = {0.3, 0.7, 0.0};
    }
    else
    {
      if(trial < 20)
      {
        move_probabilities = {0.0, 0.0, 1.0};
      }
      else
      {
        move_probabilities = {0.2, 0.4, 0.4};
      }
    }

    // pick the move accoording to the prescribed weights
    MoveType move_type = MoveType(random.categoricalDistribution(move_probabilities));

    switch (move_type)
    {
      case MoveType::BondLengthChange:
      {
        // Move: change bond-length
        // bond-energies are affected by only bond-terms
        // bond-bond needs to be handled after growing the molecule and correcting for it

        component.cbmc_moves_statistics[currentBead].bondLengthChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bondLengthChange.totalCounts += 1;

        std::size_t selected_next_bead_index = random.uniform_integer(0, nextBeads.size() - 1);
        std::size_t selected_next_bead = nextBeads[selected_next_bead_index];

        std::optional<BondPotential> bond = intraMolecularInteractions.findBondPotential(currentBead, selected_next_bead);

        if(!bond.has_value())
        {
          throw std::runtime_error(std::format("CBMC: bond-potential can not be found (internal error)\n"));
        }

        double3 current_bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;
        double current_bond_length = current_bond_vector.length();

        double3 saved_current_position = chain_atoms[selected_next_bead].position;

        if(bond->type != BondType::Fixed)
        {
          double max_change = component.cbmc_moves_statistics[currentBead].bondLengthChange.maxChange;
          double new_bond_length = current_bond_length + (2.0 * random.uniform() - 1.0) * max_change;

          double ratio = new_bond_length / current_bond_length;
          chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + ratio * current_bond_vector;

          double new_bond_energy = intraMolecularInteractions.calculateBondSmallMCEnergies(chain_atoms);

          double energy_difference = new_bond_energy - current_bond_energy;

          component.cbmc_moves_statistics[currentBead].bondLengthChange.constructed += 1;
          component.cbmc_moves_statistics[currentBead].bondLengthChange.totalConstructed += 1;

          if (random.uniform() < ratio * ratio * std::exp(-beta * energy_difference))
          {
            // update component statistics
            current_bond_energy += energy_difference;

            component.cbmc_moves_statistics[currentBead].bondLengthChange.accepted += 1;
            component.cbmc_moves_statistics[currentBead].bondLengthChange.totalAccepted += 1;
          }
          else
          {
            chain_atoms[selected_next_bead].position = saved_current_position;
          }
        }
        
#if defined(DEBUG)
        double old_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                          saved_current_position);
        double new_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                          chain_atoms[selected_next_bead].position);
        if (std::fabs(new_angle - old_angle) > 1e-5)
        {
          throw std::runtime_error(std::format("CBMC: bond-angle change in 'MoveType::BondLengthChange' ({} vs {})\n",
                                               new_angle * Units::RadiansToDegrees,
                                               old_angle * Units::RadiansToDegrees));
        }
#endif
        break;
      }
      case MoveType::BendAngleChange:
      {
        // Move: change bend angle of one chosen bead
        // bend-energies are affected by: 1) bend-terms, 2) urey-bradley 3) bond-bend, 4) bend-bend, 5) inversion-bend,
        // 6) out-of-plane bend

        component.cbmc_moves_statistics[currentBead].bendAngleChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalCounts += 1;

        std::size_t selected_next_bead_index = random.uniform_integer(0, nextBeads.size() - 1);
        std::size_t selected_next_bead = nextBeads[selected_next_bead_index];

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

        double new_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

        double energy_difference = new_bend_energy - current_bend_energy;

        component.cbmc_moves_statistics[currentBead].bendAngleChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalConstructed += 1;

        // compute energy
        if (random.uniform() < (new_sinus / current_sinus) * std::exp(-beta * energy_difference))
        {
          // update component statistics
          current_bend_energy = new_bend_energy;

          component.cbmc_moves_statistics[currentBead].bendAngleChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].bendAngleChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
#if defined(DEBUG)
        double old_bond = (chain_atoms[currentBead].position - saved_current_position).length();
        double new_bond = (chain_atoms[currentBead].position - chain_atoms[selected_next_bead].position).length();
        if (std::fabs(new_bond - old_bond) > 1e-5)
        {
          throw std::runtime_error(
              std::format("CBMC: bond-length change in 'MoveType::BendAngleChange' ({} vs {})\n", new_bond, old_bond));
        }
#endif
        break;
      }
      case MoveType::ConeChange:
      {
        // Move: rotate one position on the cone while keeping the cone angle fixed
        // bend-energies are affected by: 1) bend-terms, 2) bond-bend, 3) bend-bend, 4) inversion-bend, 5) out-of-plane
        // bend

        component.cbmc_moves_statistics[currentBead].conePositionChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalCounts += 1;

        std::size_t selected_next_bead_index = random.uniform_integer(0, nextBeads.size() - 1);
        std::size_t selected_next_bead = nextBeads[selected_next_bead_index];

        double3 bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;

        double max_change = component.cbmc_moves_statistics[currentBead].conePositionChange.maxChange;
        double random_angle = (2.0 * random.uniform() - 1.0) * max_change;
        double3 new_vec = last_bond_vector.rotateAroundAxis(bond_vector, random_angle);

        double3 saved_current_position = chain_atoms[selected_next_bead].position;
        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + new_vec;

        double new_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

        double energy_difference = new_bend_energy - current_bend_energy;

        component.cbmc_moves_statistics[currentBead].conePositionChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalConstructed += 1;

        // compute energy
        if (random.uniform() < std::exp(-beta * energy_difference))
        {
          // update component statistics
          current_bend_energy = new_bend_energy;

          component.cbmc_moves_statistics[currentBead].conePositionChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].conePositionChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
#if defined(DEBUG)
        double old_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                          saved_current_position);
        double new_angle = double3::angle(chain_atoms[previousBead].position, chain_atoms[currentBead].position,
                                          chain_atoms[selected_next_bead].position);
        if (std::fabs(new_angle - old_angle) > 1e-5)
        {
          throw std::runtime_error(std::format("CBMC: bend-angle change in 'MoveType::ConeChange' ({} vs {})\n",
                                               new_angle * Units::RadiansToDegrees,
                                               old_angle * Units::RadiansToDegrees));
        }
#endif
        break;
      }
      default:
        std::unreachable();
    }
  }

  double recomputed_bond_energy = intraMolecularInteractions.calculateBondSmallMCEnergies(chain_atoms);
  double recomputed_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

  double bond_drift = std::fabs(recomputed_bond_energy - current_bond_energy);
  double bend_drift = std::fabs(recomputed_bend_energy - current_bend_energy);

  if (bond_drift > 1e-4 || bend_drift > 1e-4)
  {
    throw std::runtime_error(std::format(
        "CBMC: internal drift in the trial-orientation Monte-Carlo scheme (bond {}, bend {})\n", bond_drift,
        bend_drift));
  }

  std::vector<Atom> next_bead_atoms(nextBeads.size());
  for (std::size_t i = 0; i != nextBeads.size(); ++i)
  {
    next_bead_atoms[i] = chain_atoms[nextBeads[i]];
  }

  return next_bead_atoms;
}

std::vector<Atom> CBMC::generateRigidUnitOrientationMonteCarloScheme(
    RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta, const Component &component,
    const std::vector<Atom> &chainAtoms, std::optional<std::size_t> previousBead, std::size_t currentBead,
    const std::vector<std::size_t> &nextBeads, const Potentials::IntraMolecularPotentials &intra)
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

  // Without a previous bead (seed group) or without junction bends every tilt is equally likely:
  // return a uniformly random orientation.
  if (!previousBead.has_value() || intra.bends.empty())
  {
    placeWithRotation(random.randomRotationMatrix());

    std::vector<Atom> result(nextBeads.size());
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      result[k] = chain_atoms[nextBeads[k]];
    }
    return result;
  }

  double3 last_bond_vector = (chainAtoms[previousBead.value()].position - anchor_position).normalized();

  // Initial orientation: align the body-fixed anchor->inner direction (inner: the group atom bonded
  // to the anchor) on a cone around the previous-anchor bond whose opening angle is drawn from the
  // primary junction bend potential, with a uniform roll. This starts the Metropolis walk close to
  // the bend minimum; the walk then equilibrates all junction bends simultaneously.
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

  // The roll about the anchor->inner axis leaves the primary bend (previous-anchor-inner) unchanged
  // but controls every other junction bend (e.g. the second ipso bend of an aromatic ring, which
  // forces the exocyclic bond into the ring plane). A uniform roll can start up to ~60 degrees off
  // against stiff bends, which the local Metropolis refinement cannot reliably recover. Instead,
  // Boltzmann-select the initial roll from the full junction bend energy on a fine grid; the
  // Metropolis walk afterwards only needs to sample local fluctuations.
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
  double roll_angle =
      roll_offset + 2.0 * std::numbers::pi * static_cast<double>(selected_roll) / numberOfRollAngles;
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    chain_atoms[nextBeads[k]].position =
        anchor_position + target_direction.rotateAroundAxis(aligned_offsets[k], roll_angle);
  }

  double current_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);

  // Metropolis MC over rigid-body rotations about the anchor. Proposals rotate the whole body by a
  // small angle about a uniformly random axis through the anchor; this is symmetric with respect to
  // the Haar (uniform rotation) measure, so plain Metropolis acceptance on the junction bend
  // energies samples the tilt from its Boltzmann distribution. No Rosenbluth weight is accumulated:
  // as with the flexible bend MC, the bend bias cancels between growth and retrace.
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
      for (std::size_t k = 0; k != nextBeads.size(); ++k)
      {
        chain_atoms[nextBeads[k]].position = saved_positions[k];
      }
    }
  }

  std::vector<Atom> result(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    result[k] = chain_atoms[nextBeads[k]];
  }
  return result;
}

CBMC::TorsionOrientation CBMC::selectTorsionOrientation(
    RandomNumber &random, std::size_t numberOfTorsionTrials, double beta, const std::vector<Atom> &chainAtoms,
    const std::vector<Atom> &baseOrientation, std::size_t previousBead, std::size_t currentBead,
    const std::vector<std::size_t> &nextBeads, double3 lastBondVector,
    const Potentials::IntraMolecularPotentials &intra, bool pinFirstToBase)
{
  std::vector<std::pair<std::vector<Atom>, double>> torsion_orientations(numberOfTorsionTrials);
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  // Bends to placed atoms other than the previous bead (e.g. the second ring neighbor when growing
  // off an aromatic ring atom) are NOT invariant under the spin about the previous-current axis:
  // that atom does not lie on the rotation axis. Such spin-variant bends must enter the selection
  // energy here, where they are Rosenbluth-weighted identically on growth and retrace; otherwise
  // the spin would silently randomize a bend the internal bend MC just equilibrated.
  std::vector<const BendPotential *> spin_variant_bends{};
  for (const BendPotential &bend : intra.bends)
  {
    if (isSpinVariantBend(bend, previousBead, currentBead, nextBeads))
    {
      spin_variant_bends.push_back(&bend);
    }
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

#if defined(DEBUG)
    // the torsion-rotation must leave the bend angles and bond lengths unchanged
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      double old_angle = double3::angle(chainAtoms[previousBead].position, chainAtoms[currentBead].position,
                                        baseOrientation[k].position);
      double new_angle = double3::angle(chainAtoms[previousBead].position, chainAtoms[currentBead].position,
                                        rotated_atoms[k].position);
      if (std::fabs(new_angle - old_angle) > 1e-5)
      {
        throw std::runtime_error(std::format("CBMC: bend-angle change in torsion-rotation ({} vs {})\n",
                                             new_angle * Units::RadiansToDegrees,
                                             old_angle * Units::RadiansToDegrees));
      }

      double old_bond_length = (chainAtoms[currentBead].position - baseOrientation[k].position).length();
      double new_bond_length = (chainAtoms[currentBead].position - rotated_atoms[k].position).length();
      if (std::fabs(new_bond_length - old_bond_length) > 1e-5)
      {
        throw std::runtime_error(std::format("CBMC: bond-distance change in torsion-rotation ({} vs {})\n",
                                             new_bond_length, old_bond_length));
      }
    }
#endif
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
