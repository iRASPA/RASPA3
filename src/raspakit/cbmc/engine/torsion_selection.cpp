module;

module cbmc_torsion_selection;

import std;

import atom;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import running_energy;
import intra_molecular_potentials;
import bend_potential;
import cbmc_util;
import cbmc_grow_step;
import cbmc_closure_guide;
import cbmc_lookahead_guide;

CBMC::TorsionOrientation CBMC::selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials,
                                                        double beta, const ForceField &forceField,
                                                        const SimulationBox &simulationBox,
                                                        std::vector<Atom> &chainAtoms,
                                                        const std::vector<Atom> &baseOrientation,
                                                        const GrowStep &step, double3 lastBondVector,
                                                        bool pinFirstToBase)
{
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const Potentials::IntraMolecularPotentials &intra = step.spin.potentials;
  const double3 anchor = chainAtoms[currentBead].position;

  // The next-beads of the chain are scratch for the trial spins and are restored on return.
  const ScratchBeads scratch(chainAtoms, nextBeads);

  // Writes the base orientation spun by 'angle' about the previous-current axis into the chain.
  auto placeSpin = [&](double angle)
  {
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      Atom &bead = chainAtoms[nextBeads[k]];
      bead = baseOrientation[k];
      bead.position = anchor + lastBondVector.rotateAroundAxis(baseOrientation[k].position - anchor, angle);
    }
  };

  // Closure guides (fixed-endpoint regrowth, see cbmc_closure_guide): a bias log g(D) on the distance
  // between a guided next bead and its closure target, added to every trial's log factor and divided
  // out of the weight for the selected spin again -- the sampled distribution is then exact for any g,
  // which only steers the spin towards geometries the segment can still close from. The tables need
  // the temperature and are filled by 'Component::prepareGrowthPlans'; an unprepared guide is a setup
  // error, not something the hot path repairs.
  const std::vector<GrowStep::SpinSelectionData::ClosureGuide> &guides = step.spin.guides;
  for (const GrowStep::SpinSelectionData::ClosureGuide &guide : guides)
  {
    if (!guide.table)
    {
      throw std::logic_error(std::format(
          "CBMC: growth step at anchor bead {} carries a closure guide (bead {} -> {}) without a prepared table; "
          "call Component::prepareGrowthPlans(beta) before growing with it\n",
          currentBead, nextBeads[guide.nextBeadIndex], guide.targetBead));
    }
  }
  // The lookahead guide (cbmc_lookahead_guide): the same construction, a bias log g on the dihedral of
  // the first grown bead against the reference neighbour of the previous bead, averaging the terms that
  // couple the spin to the beads of the following steps.
  const std::optional<GrowStep::SpinSelectionData::LookaheadGuide> &lookahead = step.spin.lookahead;
  if (lookahead.has_value() && !lookahead->table)
  {
    throw std::logic_error(std::format(
        "CBMC: growth step at anchor bead {} carries a lookahead guide without a prepared table; "
        "call Component::prepareGrowthPlans(beta) before growing with it\n",
        currentBead));
  }
  const bool hasGuides = step.spin.hasGuides();
  auto logGuideOfPlacedSpin = [&]() -> double
  {
    double logGuide = 0.0;
    for (const GrowStep::SpinSelectionData::ClosureGuide &guide : guides)
    {
      const double distance =
          (chainAtoms[nextBeads[guide.nextBeadIndex]].position - chainAtoms[guide.targetBead].position).length();
      logGuide += guide.table->logGuideAt(distance);
    }
    if (lookahead.has_value())
    {
      const double dihedral =
          dihedralAngle(chainAtoms[lookahead->referenceBead].position, chainAtoms[step.previousBead.value()].position,
                        anchor, chainAtoms[nextBeads[0]].position);
      logGuide += lookahead->table->logGuideAt(dihedral);
    }
    return logGuide;
  };

  // Only the spin angle and the log factor of each trial are kept; the positions of the selected spin
  // are regenerated afterwards (one placement instead of one vector per trial).
  std::vector<double> angles(numberOfTorsionTrials);
  std::vector<double> logTorsionBoltzmannFactors(numberOfTorsionTrials);
  std::vector<double> logGuides(hasGuides ? numberOfTorsionTrials : 0, 0.0);

  for (std::size_t j = 0; j != numberOfTorsionTrials; ++j)
  {
    const double angle = (pinFirstToBase && j == 0) ? 0.0 : (2.0 * random.uniform() - 1.0) * std::numbers::pi;
    placeSpin(angle);

    double torsion_energy = intra.calculateTorsionEnergies(chainAtoms);
    for (const BendPotential &bend : step.spin.variantBends)
    {
      torsion_energy += bend.calculateEnergy(chainAtoms[bend.identifiers[0]].position,
                                             chainAtoms[bend.identifiers[1]].position,
                                             chainAtoms[bend.identifiers[2]].position, std::nullopt);
    }
    // The step's spin-routed unsampled terms (see the growth plan's term classification): couplings
    // to placed geometry or to the spin angle, steering the spin choice through the selection.
    if (step.spin.hasUnsampledTerms)
    {
      torsion_energy += intra.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
    }
    // The spin-routed intramolecular van der Waals / Coulomb pairs (grown bead to a placed bead a few
    // bonds away, see 'GrowStep::NonBondedData'): their distances are set by the spin, so a 1-5 clash
    // is resolved here rather than after the spin has been chosen on the bare torsion.
    if (step.spin.hasNonBondedTerms)
    {
      torsion_energy += intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(forceField, simulationBox, chainAtoms).potentialEnergy();
    }

    angles[j] = angle;
    logTorsionBoltzmannFactors[j] = -beta * torsion_energy;
    if (hasGuides)
    {
      logGuides[j] = logGuideOfPlacedSpin();
      logTorsionBoltzmannFactors[j] += logGuides[j];
    }
  }

  const std::size_t selected_torsion =
      pinFirstToBase ? 0 : CBMC::selectTrialPosition(random, logTorsionBoltzmannFactors);

  placeSpin(angles[selected_torsion]);
  std::vector<Atom> positions(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) positions[k] = chainAtoms[nextBeads[k]];

  double rosenbluth_weight_torsion;
  if (!hasGuides)
  {
    rosenbluth_weight_torsion =
        std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(), 0.0,
                        [](const double &acc, const double &logFactor) { return acc + std::exp(logFactor); });
  }
  else
  {
    // Guided: log-sum-exp (a guide can shift the log factors by tens of units), and the guide of the
    // selected spin divided out (see above).
    const double maxLogFactor =
        *std::max_element(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end());
    const double logRosenbluthSum =
        maxLogFactor + std::log(std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(),
                                                0.0, [&](const double &acc, const double &logFactor)
                                                { return acc + std::exp(logFactor - maxLogFactor); }));
    rosenbluth_weight_torsion = std::exp(logRosenbluthSum - logGuides[selected_torsion]);
  }

  return {std::move(positions), rosenbluth_weight_torsion / static_cast<double>(numberOfTorsionTrials)};
}

void CBMC::prepareLookaheadGuides(double beta, std::span<GrowStep> plan,
                                  std::map<std::string, std::shared_ptr<const LookaheadGuideTable>> &memo)
{
  for (GrowStep &step : plan)
  {
    if (!step.spin.lookahead.has_value()) continue;
    GrowStep::SpinSelectionData::LookaheadGuide &guide = step.spin.lookahead.value();
    auto it = memo.find(guide.signature);
    if (it == memo.end())
    {
      it = memo.emplace(guide.signature,
                        std::make_shared<const LookaheadGuideTable>(buildLookaheadGuideTable(beta, guide.model)))
               .first;
    }
    guide.table = it->second;
  }
}
