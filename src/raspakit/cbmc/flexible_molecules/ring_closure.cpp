module;

module cbmc_ring_closure;

import std;

import randomnumbers;
import units;
import component;
import atom;
import double3;
import double3x3;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;

// A bend is spin-variant (changes under a rigid spin of the ring about the previous-anchor bond) when
// it involves a placed atom other than the previous bead, the anchor, or the ring atoms. This mirrors
// 'isSpinVariantBend' in 'generate_trial_orientations_mc.cpp': a bend on the axis
// (previous-anchor-ring) is spin-invariant, so the primary junction (tilt) bends are sampled by the
// ring conformational Monte-Carlo, while off-axis junction bends are weighted by the torsion step.
static bool bendIsSpinVariant(const std::array<std::size_t, 3> &identifiers, std::optional<std::size_t> previousBead,
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

bool CBMC::isSpinVariantRingTorsion(const std::array<std::size_t, 4> &identifiers, std::size_t currentBead,
                                    const std::vector<std::size_t> &nextBeads)
{
  bool hasRingAtom = false;
  bool hasOutsideAtom = false;
  for (std::size_t id : identifiers)
  {
    bool inRingBody =
        (id == currentBead) || (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end());
    bool moves = std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end();
    if (moves) hasRingAtom = true;
    if (!inRingBody) hasOutsideAtom = true;
  }
  return hasRingAtom && hasOutsideAtom;
}

std::vector<Atom> CBMC::randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                           double3 anchorPosition)
{
  double3x3 rotation = random.randomRotationMatrix();
  std::vector<Atom> result = ringAtoms;
  for (Atom &atom : result)
  {
    atom.position = anchorPosition + rotation * (atom.position - anchorPosition);
  }
  return result;
}

Potentials::IntraMolecularPotentials CBMC::ringSpinPotentials(const Potentials::IntraMolecularPotentials &intra,
                                                             std::size_t currentBead,
                                                             const std::vector<std::size_t> &nextBeads)
{
  Potentials::IntraMolecularPotentials spin{};

  // 'selectTorsionOrientation' filters the spin-variant bends itself, so keep all bends.
  spin.bends = intra.bends;

  spin.torsions.reserve(intra.torsions.size());
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads))
    {
      spin.torsions.push_back(torsion);
    }
  }

  return spin;
}

std::vector<Atom> CBMC::generateRingConformationMonteCarloScheme(
    RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta, const Component &component,
    const std::vector<Atom> &chainAtoms, std::optional<std::size_t> previousBead, std::size_t currentBead,
    const std::vector<std::size_t> &nextBeads, const Potentials::IntraMolecularPotentials &intra)
{
  double3 anchor_reference = component.atoms[currentBead].position;
  double3 anchor_position = chainAtoms[currentBead].position;

  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  // Initial conformation: the component reference ring geometry (a valid closed ring) placed with a
  // random overall rotation about the anchor. The internal Monte-Carlo below then samples the
  // conformational fluctuations, and (with a junction) the rigid rotations equilibrate the tilt.
  auto placeWithRotation = [&](const double3x3 &rotation)
  {
    for (std::size_t atom : nextBeads)
    {
      double3 offset = component.atoms[atom].position - anchor_reference;
      chain_atoms[atom].position = anchor_position + rotation * offset;
    }
  };
  placeWithRotation(random.randomRotationMatrix());

  // Base-conformation potentials: everything sampled here, i.e. the internal ring bonds/bends/torsions
  // and the junction (tilt) bends, but not the spin-selected junction-crossing bends and torsions.
  Potentials::IntraMolecularPotentials baseIntra{};
  baseIntra.bonds = intra.bonds;
  for (const BendPotential &bend : intra.bends)
  {
    if (!bendIsSpinVariant(bend.identifiers, previousBead, currentBead, nextBeads))
    {
      baseIntra.bends.push_back(bend);
    }
  }
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (!isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads))
    {
      baseIntra.torsions.push_back(torsion);
    }
  }

  auto baseEnergy = [&](std::vector<Atom> &atoms)
  {
    return baseIntra.calculateBondSmallMCEnergies(atoms) + baseIntra.calculateBendSmallMCEnergies(atoms) +
           baseIntra.calculateTorsionEnergies(atoms);
  };

  double current_energy = baseEnergy(chain_atoms);

  // Metropolis Monte-Carlo sampling the internal ring conformation from its Boltzmann distribution.
  // Single-atom Cartesian displacements (symmetric proposal in Cartesian space, so they carry the
  // r^2 sin(theta) volume elements automatically) sample the internal bonds/bends/torsions; the
  // closure bond keeps the ring closed. Small rigid rotations of the whole ring about the anchor
  // sample the tilt (junction bends). As with the flexible bend Monte-Carlo, no Rosenbluth weight is
  // accumulated: the internal energy is absorbed into the sampling and the junction-crossing terms are
  // weighted by the caller's torsion step.
  constexpr double maximumDisplacement = 0.1;    // Angstrom
  constexpr double maximumRotationAngle = 0.15;  // radian
  const bool haveJunction = previousBead.has_value();

  std::size_t number_of_trials = numberOfTrialMovesPerOpenBead * nextBeads.size();

  std::vector<double3> saved(nextBeads.size());

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    if (haveJunction && random.uniform() < 0.25)
    {
      // Rigid rotation of the whole ring about the anchor: samples the junction (tilt) bends. The
      // proposal is symmetric with respect to the uniform rotation measure.
      double3 axis = random.randomVectorOnUnitSphere();
      double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;

      for (std::size_t k = 0; k != nextBeads.size(); ++k)
      {
        saved[k] = chain_atoms[nextBeads[k]].position;
        chain_atoms[nextBeads[k]].position = anchor_position + axis.rotateAroundAxis(saved[k] - anchor_position, angle);
      }

      double trial_energy = baseEnergy(chain_atoms);
      if (random.uniform() < std::exp(-beta * (trial_energy - current_energy)))
      {
        current_energy = trial_energy;
      }
      else
      {
        for (std::size_t k = 0; k != nextBeads.size(); ++k)
        {
          chain_atoms[nextBeads[k]].position = saved[k];
        }
      }
    }
    else
    {
      // Single-atom Cartesian displacement: samples the internal ring conformation.
      std::size_t idx = random.uniform_integer(0, nextBeads.size() - 1);
      std::size_t atom = nextBeads[idx];

      double3 displacement{(2.0 * random.uniform() - 1.0) * maximumDisplacement,
                           (2.0 * random.uniform() - 1.0) * maximumDisplacement,
                           (2.0 * random.uniform() - 1.0) * maximumDisplacement};

      double3 saved_position = chain_atoms[atom].position;
      chain_atoms[atom].position = saved_position + displacement;

      double trial_energy = baseEnergy(chain_atoms);
      if (random.uniform() < std::exp(-beta * (trial_energy - current_energy)))
      {
        current_energy = trial_energy;
      }
      else
      {
        chain_atoms[atom].position = saved_position;
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
