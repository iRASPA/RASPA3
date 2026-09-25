module;

module cbmc_ring_closure;

import std;

import atom;
import double3;
import double3x3;
import randomnumbers;
import cbmc_grow_context;
import component;
import move_statistics;
import chiral_center;
import bond_potential;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_constants;
import cbmc_statistics;
import cbmc_grow_step;

namespace
{
namespace Constants = CBMC::Constants;

// Samples the internal conformation of a cyclic cluster (kept closed by its closure bonds) plus the
// junction tilt. Carries no Rosenbluth weight. The step-constant data (move units, Fixed-bond pivots,
// crankshaft candidates, monitored chiral centres, the ring-internal bonded terms) was derived once
// when the plan was built ('GrowStep::ring'); the sampler binds it to one chain and runs:
//
//   1. seed:      an independent ring conformation, rigidly oriented about the anchor with the parity
//                 of every monitored chiral centre preserved;
//   2. relax:     an internal Metropolis Monte-Carlo on the ring-internal bonded terms, one of four
//                 symmetric moves per trial (crankshaft, whole-ring tilt, single-atom move, rigid-unit
//                 move), all constraint-preserving, all parity-guarded;
//   3. result:    the positions of the step's next-beads.
//
// A cyclic cluster may pass through rigid-body fragments (e.g. a macrocycle of rigid rings linked by
// flexible bridges); those fragments are moved as whole rigid units so their internal geometry stays
// exact.
class RingSampler
{
 public:
  RingSampler(RandomNumber &random, const CBMC::GrowthSettings &settings, double beta, const Component &component,
              std::vector<Atom> &chainAtoms, const CBMC::GrowStep &step)
      : random_(random),
        settings_(settings),
        beta_(beta),
        component_(component),
        step_(step),
        ring_(step.ring),
        nextBeads_(step.nextBeads),
        chain_(chainAtoms),
        scratch_(chainAtoms, step),
        anchorPosition_(chainAtoms[step.currentBead].position),
        displacementStats_(component.cbmcMoveStatistics[step.currentBead].ringDisplacementChange),
        rotationStats_(component.cbmcMoveStatistics[step.currentBead].ringRotationChange),
        crankshaftStats_(component.cbmcMoveStatistics[step.currentBead].ringCrankshaftMove),
        saved_(step.nextBeads.size())
  {
    // Reference parities of the chiral centres whose four atoms all have positions during this step:
    // the bond/bend/torsion model is achiral (a mirror image has the same energy), so without this
    // guard the internal MC could invert a declared stereocentre, e.g. flip a cis ring fusion to trans.
    monitoredCenters_.reserve(ring_.monitoredChiralCenters.size());
    for (const ChiralCenter &center : ring_.monitoredChiralCenters)
    {
      monitoredCenters_.push_back({&center.ids, CBMC::chiralSignedVolume(center.ids, component.atoms)});
    }
  }

  std::vector<Atom> sample()
  {
    seed();
    relax();

    std::vector<Atom> result(nextBeads_.size());
    for (std::size_t k = 0; k != nextBeads_.size(); ++k) result[k] = chain_[nextBeads_[k]];
    return result;
  }

 private:
  // ----------------------------------------------------------------------------------------------
  // Seeding.
  // ----------------------------------------------------------------------------------------------

  // Seed geometry: an independent, well-mixed ideal-gas conformation from the reservoir (so different
  // grows start from different ring puckers), or the reference geometry while the reservoir is still
  // being built. A reservoir member keeps the exact internal geometry of any rigid sub-fragment, so
  // rigid parts of the ring stay rigid. The declared chirality reference stays 'component.atoms'.
  const std::vector<Atom> &seedGeometry()
  {
    const std::vector<std::vector<Atom>> &reservoir = component_.conformationReservoir;
    if (reservoir.empty()) return component_.atoms;
    return reservoir[random_.uniform_integer(0, reservoir.size() - 1)];
  }

  void placeSeedWithRotation(const std::vector<Atom> &seedAtoms, const double3x3 &rotation)
  {
    const double3 anchorReference = seedAtoms[step_.currentBead].position;
    for (std::size_t atom : nextBeads_)
    {
      chain_[atom].position = anchorPosition_ + rotation * (seedAtoms[atom].position - anchorReference);
    }
  }

  // Reflects the placed ring body through a plane containing the previous-current axis.
  //
  // A proper rotation of the reference geometry preserves the parity of every chiral centre that lies
  // entirely inside the rigidly placed ring body, so only a centre that also involves the junction's
  // placed neighbour can come out with the wrong parity. Reflecting the body through a plane that
  // contains the axis flips exactly those junction-involving centres while preserving all internal
  // distances and the bend angles to the anchor -- the bonded model is achiral, so the reflected seed
  // has identical energy. Because a reflection is improper it also flips any body-internal centre, so
  // it is only a valid fix when it restores every monitored parity at once: the caller keeps it only
  // when 'parityPreserved()' then holds. This makes the dominant single-junction-stereocentre case an
  // O(1), deterministic correction instead of a random re-roll.
  void reflectThroughJunctionPlane()
  {
    const double3 axis = (chain_[step_.previousBead.value()].position - anchorPosition_).normalized();
    const double3 normal = double3::perpendicular(axis, random_.randomVectorOnUnitSphere());
    for (std::size_t atom : nextBeads_)
    {
      const double3 relative = chain_[atom].position - anchorPosition_;
      chain_[atom].position = anchorPosition_ + (relative - 2.0 * double3::dot(relative, normal) * normal);
    }
  }

  void seed()
  {
    const std::vector<Atom> &seedAtoms = seedGeometry();

    placeSeedWithRotation(seedAtoms, random_.randomRotationMatrix());
    if (!parityPreserved() && step_.previousBead.has_value())
    {
      saveNextBeads();
      reflectThroughJunctionPlane();
      if (!parityPreserved()) restoreNextBeads();
    }
    for (std::size_t attempt = 0; attempt != Constants::ringParityMaximumReorientations && !parityPreserved();
         ++attempt)
    {
      placeSeedWithRotation(seedAtoms, random_.randomRotationMatrix());
    }
  }

  // ----------------------------------------------------------------------------------------------
  // The internal Monte-Carlo.
  // ----------------------------------------------------------------------------------------------

  // The bonded terms the internal MC samples (bonds, spin-invariant bends and torsions); the
  // spin-variant terms are weighted in the torsion selection.
  double baseEnergy()
  {
    const Potentials::IntraMolecularPotentials &intra = ring_.internalPotentials;
    return intra.calculateBondSmallMCEnergies(chain_) + intra.calculateBendSmallMCEnergies(chain_) +
           intra.calculateTorsionEnergies(chain_);
  }

  bool parityPreserved() const
  {
    for (const auto &[ids, referenceVolume] : monitoredCenters_)
    {
      if (CBMC::chiralSignedVolume(*ids, chain_) * referenceVolume < 0.0) return false;
    }
    return true;
  }

  // Adaptive internal-MC step sizes: the maximum displacement and rotation angle are read from the
  // anchor bead's CBMC statistics and adapted towards the target acceptance ratio between sweeps by
  // 'System::optimizeMCMoves'. These moves carry no Rosenbluth weight, so their step size affects
  // only sampling efficiency, not detailed balance.
  void relax()
  {
    currentEnergy_ = baseEnergy();
    const bool haveJunction = step_.previousBead.has_value();
    const std::size_t numberOfTrials = 2 * settings_.numberOfTrialMovesPerOpenBead * nextBeads_.size();

    for (std::size_t trial = 0; trial != numberOfTrials; ++trial)
    {
      if (!ring_.crankshafts.empty() && random_.uniform() < settings_.ringCrankshaftProbability)
      {
        crankshaftMove();
      }
      else if (haveJunction && random_.uniform() < settings_.ringTiltProbability)
      {
        tiltMove();
      }
      else
      {
        const std::vector<std::size_t> &unit = ring_.moveUnits[random_.uniform_integer(0, ring_.moveUnits.size() - 1)];
        if (unit.size() == 1)
        {
          singleAtomMove(unit[0]);
        }
        else
        {
          rigidUnitMove(unit);
        }
      }
    }
  }

  // Metropolis on the base energy for a symmetric proposal already written into 'chain_'; 'restore'
  // undoes it on rejection (or on a parity violation).
  template <typename Restore>
  void acceptOrReject(MoveStatistics<double> &stats, Restore restore)
  {
    stats.counts += 1.0;
    stats.totalCounts += 1.0;
    stats.constructed += 1.0;
    stats.totalConstructed += 1.0;
    const double trialEnergy = baseEnergy();
    if (parityPreserved() && random_.uniform() < std::exp(-beta_ * (trialEnergy - currentEnergy_)))
    {
      currentEnergy_ = trialEnergy;
      stats.accepted += 1.0;
      stats.totalAccepted += 1.0;
    }
    else
    {
      restore();
    }
  }

  double3 randomDisplacement(double maximum)
  {
    return {(2.0 * random_.uniform() - 1.0) * maximum, (2.0 * random_.uniform() - 1.0) * maximum,
            (2.0 * random_.uniform() - 1.0) * maximum};
  }

  double randomAngle(double maximum) { return (2.0 * random_.uniform() - 1.0) * maximum; }

  // Conformer-hopping crankshaft: a large-angle rotation of one flexible ring atom about the line
  // through two of its positioned bonded neighbours, which preserves both bond lengths exactly while
  // flipping the local pucker (chair <-> twist-boat) and axial <-> equatorial placement -- the barrier
  // crossing the small adaptive moves almost never make. Attempted a fraction of the time
  // ('CBMCRingCrankshaftProbability') so local relaxation still dominates. Tracked in its own
  // statistics: the angle is deliberately full-range and never adapted, and pooling its acceptances
  // into the adaptive rotation statistics would distort that step-size optimization.
  void crankshaftMove()
  {
    const CBMC::GrowStep::RingCrankshaft &c = ring_.crankshafts[random_.uniform_integer(0, ring_.crankshafts.size() - 1)];
    const double3 pivot = chain_[c.axisA].position;
    const double3 axis = (chain_[c.axisB].position - pivot).normalized();
    const double angle = randomAngle(std::numbers::pi);
    const double3 savedPosition = chain_[c.atom].position;
    chain_[c.atom].position = pivot + axis.rotateAroundAxis(savedPosition - pivot, angle);
    acceptOrReject(crankshaftStats_, [&] { chain_[c.atom].position = savedPosition; });
  }

  // Whole-ring tilt about the anchor (junction steps only): relaxes the junction bends.
  void tiltMove()
  {
    const double3 axis = random_.randomVectorOnUnitSphere();
    const double angle = randomAngle(rotationStats_.maxChange);
    saveNextBeads();
    for (std::size_t atom : nextBeads_)
    {
      chain_[atom].position = anchorPosition_ + axis.rotateAroundAxis(chain_[atom].position - anchorPosition_, angle);
    }
    acceptOrReject(rotationStats_, [&] { restoreNextBeads(); });
  }

  // Local single-atom move, constraint-preserving. A Fixed bond is a holonomic constraint without
  // energy, so a free displacement would stretch it unpenalized: with no Fixed bond the atom is
  // displaced freely; with Fixed bonds it is rotated about its fixed neighbour(s) so those exact
  // lengths are kept -- one fixed neighbour leaves a full sphere (rotate about a random axis through
  // it), two leave a circle (rotate about the line through both), and three or more pin the atom, so
  // it moves only through the whole-cluster rotation or a concerted move.
  void singleAtomMove(std::size_t atom)
  {
    const std::vector<std::size_t> &fixed = ring_.fixedNeighbors[atom];
    if (fixed.size() > 2) return;

    const double3 savedPosition = chain_[atom].position;
    auto restore = [&] { chain_[atom].position = savedPosition; };
    if (fixed.empty())
    {
      chain_[atom].position += randomDisplacement(displacementStats_.maxChange);
      acceptOrReject(displacementStats_, restore);
      return;
    }
    const double3 pivot = chain_[fixed[0]].position;
    const double3 axis =
        fixed.size() == 1 ? random_.randomVectorOnUnitSphere() : (chain_[fixed[1]].position - pivot).normalized();
    const double angle = randomAngle(rotationStats_.maxChange);
    chain_[atom].position = pivot + axis.rotateAroundAxis(savedPosition - pivot, angle);
    acceptOrReject(rotationStats_, restore);
  }

  // Rigid sub-fragment: a symmetric whole-unit translation or rotation about the unit centre, so its
  // internal geometry stays exact.
  void rigidUnitMove(const std::vector<std::size_t> &unit)
  {
    for (std::size_t k = 0; k != unit.size(); ++k) saved_[k] = chain_[unit[k]].position;
    auto restore = [&]
    {
      for (std::size_t k = 0; k != unit.size(); ++k) chain_[unit[k]].position = saved_[k];
    };

    if (random_.uniform() < 0.5)
    {
      const double3 displacement = randomDisplacement(displacementStats_.maxChange);
      for (std::size_t atom : unit) chain_[atom].position += displacement;
      acceptOrReject(displacementStats_, restore);
      return;
    }

    double3 center{};
    for (std::size_t atom : unit) center += chain_[atom].position;
    center = center / static_cast<double>(unit.size());
    const double3 axis = random_.randomVectorOnUnitSphere();
    const double angle = randomAngle(rotationStats_.maxChange);
    for (std::size_t atom : unit)
    {
      chain_[atom].position = center + axis.rotateAroundAxis(chain_[atom].position - center, angle);
    }
    acceptOrReject(rotationStats_, restore);
  }

  void saveNextBeads()
  {
    for (std::size_t k = 0; k != nextBeads_.size(); ++k) saved_[k] = chain_[nextBeads_[k]].position;
  }
  void restoreNextBeads()
  {
    for (std::size_t k = 0; k != nextBeads_.size(); ++k) chain_[nextBeads_[k]].position = saved_[k];
  }

  RandomNumber &random_;
  const CBMC::GrowthSettings &settings_;
  const double beta_;
  const Component &component_;
  const CBMC::GrowStep &step_;
  const CBMC::GrowStep::RingSamplerData &ring_;
  const std::vector<std::size_t> &nextBeads_;

  /// The chain the step grows in; the ring body's beads are edited in place and restored by 'scratch_'
  /// when the sampler is destroyed (see the scratch contract in cbmc_operators).
  std::vector<Atom> &chain_;
  const CBMC::ScratchBeads scratch_;
  const double3 anchorPosition_;
  std::vector<std::pair<const std::array<std::size_t, 4> *, double>> monitoredCenters_{};
  double currentEnergy_{0.0};

  MoveStatistics<double> &displacementStats_;
  MoveStatistics<double> &rotationStats_;
  MoveStatistics<double> &crankshaftStats_;
  std::vector<double3> saved_;  ///< Scratch for undoing a multi-atom move (sized to the next-beads).
};
}  // namespace

std::vector<Atom> CBMC::randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                           double3 anchorPosition)
{
  double3x3 rotation = random.randomRotationMatrix();
  std::vector<Atom> result = ringAtoms;
  for (Atom &atom : result) atom.position = anchorPosition + rotation * (atom.position - anchorPosition);
  return result;
}

std::vector<Atom> CBMC::generateRingConformation(RandomNumber &random, const GrowthSettings &settings, double beta,
                                                 const Component &component, std::vector<Atom> &chainAtoms,
                                                 const GrowStep &step)
{
  return RingSampler(random, settings, beta, component, chainAtoms, step).sample();
}
