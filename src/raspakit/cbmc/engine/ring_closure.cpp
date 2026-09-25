module;

module cbmc_ring_closure;

import std;

import atom;
import double3;
import double3x3;
import randomnumbers;
import cbmc_growth_context;
import component;
import move_statistics;
import chiral_center;
import bond_potential;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_constants;
import cbmc_statistics;
import cbmc_growth_plan;

namespace Constants = CBMC::Constants;

// Ring-closure conformation: samples the internal conformation of a cyclic cluster (kept closed by
// its closure bonds) plus the junction tilt. Carries no Rosenbluth weight. Ported from the former
// 'generateRingConformationMonteCarloScheme'. A cyclic cluster may pass through rigid-body
// fragments (e.g. a macrocycle of rigid rings linked by flexible bridges); those fragments are
// moved as whole rigid units so their internal geometry stays exact.
std::vector<Atom> CBMC::randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                           double3 anchorPosition)
{
  double3x3 rotation = random.randomRotationMatrix();
  std::vector<Atom> result = ringAtoms;
  for (Atom &atom : result) atom.position = anchorPosition + rotation * (atom.position - anchorPosition);
  return result;
}

std::vector<Atom> CBMC::generateRingConformation(RandomNumber &random, const GrowthSettings &settings, double beta,
                                                 const Component &component, const std::vector<Atom> &chainAtoms,
                                                 const GrowStep &step)
{
  const std::optional<std::size_t> previousBead = step.previousBead;
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;

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
  // anchor, and the junction's placed neighbor; listed in the plan) keep the parity of the reference
  // geometry: the bond/bend/torsion model is achiral (a mirror image has the same energy), so without
  // this guard the internal MC could invert a declared stereocenter, e.g. flip a cis ring fusion to
  // trans.
  std::vector<std::pair<const std::array<std::size_t, 4> *, double>> monitored_centers{};
  monitored_centers.reserve(step.ring.monitoredChiralCenters.size());
  for (const ChiralCenter &center : step.ring.monitoredChiralCenters)
  {
    monitored_centers.push_back({&center.ids, chiralSignedVolume(center.ids, component.atoms)});
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
  for (std::size_t attempt = 0; attempt != Constants::ringParityMaximumReorientations && !parityPreserved();
       ++attempt)
  {
    placeWithRotation(random.randomRotationMatrix());
  }

  // The step-constant data of the internal MC, derived once when the plan was built: the move units
  // (a flexible ring atom alone, a rigid sub-fragment as one unit so its internal geometry is exact),
  // the Fixed-bond neighbours every move of an atom must pivot about (a Fixed bond is a holonomic
  // constraint without energy, so a free displacement would stretch it unpenalized), and the
  // conformer-hopping crankshaft candidates: a flexible ring atom rotated by a large angle about the
  // line through two of its positioned bonded neighbours, which preserves both bond lengths exactly
  // while flipping the local pucker (chair <-> twist-boat) and axial <-> equatorial placement -- the
  // barrier crossing the small adaptive moves almost never make. The proposal is symmetric (uniform
  // +/- angle), so plain Metropolis on the base energy keeps the sampled distribution exact.
  const std::vector<std::vector<std::size_t>> &moveUnits = step.ring.moveUnits;
  const std::vector<std::vector<std::size_t>> &fixedNeighbors = step.ring.fixedNeighbors;
  const std::vector<CBMC::GrowStep::RingCrankshaft> &crankshafts = step.ring.crankshafts;

  // The bonded terms the internal MC samples (bonds, spin-invariant bends and torsions), precomputed
  // in the plan; the spin-variant terms are weighted in the torsion selection.
  const Potentials::IntraMolecularPotentials &baseIntra = step.ringInternalPotentials;

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
  std::size_t number_of_trials = 2 * settings.numberOfTrialMovesPerOpenBead * nextBeads.size();
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
    if (!crankshafts.empty() && random.uniform() < settings.ringCrankshaftProbability)
    {
      const CBMC::GrowStep::RingCrankshaft &c = crankshafts[random.uniform_integer(0, crankshafts.size() - 1)];
      double3 pivot = chain_atoms[c.axisA].position;
      double3 axis = (chain_atoms[c.axisB].position - pivot).normalized();
      double angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
      double3 savedPosition = chain_atoms[c.atom].position;
      chain_atoms[c.atom].position = pivot + axis.rotateAroundAxis(savedPosition - pivot, angle);
      acceptOrReject(crankshaftStats, 1, [&](std::size_t) { chain_atoms[c.atom].position = savedPosition; });
      continue;
    }

    if (haveJunction && random.uniform() < settings.ringTiltProbability)
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
