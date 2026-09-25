module;

module cbmc_ring_closure;

import std;

import atom;
import double3;
import double3x3;
import randomnumbers;
import forcefield;
import component;
import fragment;
import fragment_graph;
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

std::vector<Atom> CBMC::generateRingConformation(RandomNumber &random, const ForceField &forceField, double beta,
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
  for (std::size_t attempt = 0; attempt != Constants::ringParityMaximumReorientations && !parityPreserved();
       ++attempt)
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
  // preserves the distance to it). 'step.intra.bonds' only involves atoms placed in this step, all of
  // which have positions, so every listed endpoint is usable as a pivot.
  std::vector<std::vector<std::size_t>> fixedNeighbors(chain_atoms.size());
  for (const BondPotential &bond : step.intra.bonds)
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
