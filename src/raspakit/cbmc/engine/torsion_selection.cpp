module;

module cbmc_torsion_selection;

import std;

import atom;
import double3;
import randomnumbers;
import running_energy;
import intra_molecular_potentials;
import bend_potential;
import cbmc_util;
import cbmc_growth_plan;

CBMC::TorsionOrientation CBMC::selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials,
                                                        double beta, std::vector<Atom> &chainAtoms,
                                                        const std::vector<Atom> &baseOrientation,
                                                        const GrowStep &step, double3 lastBondVector,
                                                        bool pinFirstToBase)
{
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const Potentials::IntraMolecularPotentials &intra = step.torsionSelectionPotentials;
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

  // Only the spin angle and the log Boltzmann factor of each trial are kept; the positions of the
  // selected spin are regenerated afterwards (one placement instead of one vector per trial).
  std::vector<double> angles(numberOfTorsionTrials);
  std::vector<double> logTorsionBoltzmannFactors(numberOfTorsionTrials);

  for (std::size_t j = 0; j != numberOfTorsionTrials; ++j)
  {
    const double angle = (pinFirstToBase && j == 0) ? 0.0 : (2.0 * random.uniform() - 1.0) * std::numbers::pi;
    placeSpin(angle);

    double torsion_energy = intra.calculateTorsionEnergies(chainAtoms);
    for (const BendPotential &bend : step.spinVariantBends)
    {
      torsion_energy += bend.calculateEnergy(chainAtoms[bend.identifiers[0]].position,
                                             chainAtoms[bend.identifiers[1]].position,
                                             chainAtoms[bend.identifiers[2]].position, std::nullopt);
    }
    // The step's spin-routed unsampled terms (see the growth plan's term classification): couplings
    // to placed geometry or to the spin angle, steering the spin choice through the selection.
    if (step.torsionSelectionHasUnsampledTerms)
    {
      torsion_energy += intra.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
    }

    angles[j] = angle;
    logTorsionBoltzmannFactors[j] = -beta * torsion_energy;
  }

  const double rosenbluth_weight_torsion =
      std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(), 0.0,
                      [](const double &acc, const double &logFactor) { return acc + std::exp(logFactor); });

  const std::size_t selected_torsion =
      pinFirstToBase ? 0 : CBMC::selectTrialPosition(random, logTorsionBoltzmannFactors);

  placeSpin(angles[selected_torsion]);
  std::vector<Atom> positions(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) positions[k] = chainAtoms[nextBeads[k]];

  return {std::move(positions), rosenbluth_weight_torsion / static_cast<double>(numberOfTorsionTrials)};
}
