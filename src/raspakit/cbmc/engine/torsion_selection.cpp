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
                                                        double beta, const std::vector<Atom> &chainAtoms,
                                                        const std::vector<Atom> &baseOrientation,
                                                        const GrowStep &step, double3 lastBondVector,
                                                        bool pinFirstToBase)
{
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const Potentials::IntraMolecularPotentials &intra = step.torsionSelectionPotentials;

  std::vector<std::pair<std::vector<Atom>, double>> torsion_orientations(numberOfTorsionTrials);
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

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
    for (const BendPotential &bend : step.spinVariantBends)
    {
      torsion_energy += bend.calculateEnergy(chain_atoms[bend.identifiers[0]].position,
                                             chain_atoms[bend.identifiers[1]].position,
                                             chain_atoms[bend.identifiers[2]].position, std::nullopt);
    }
    // The step's spin-routed unsampled terms (see the growth plan's term classification): couplings
    // to placed geometry or to the spin angle, steering the spin choice through the selection.
    if (step.torsionSelectionHasUnsampledTerms)
    {
      torsion_energy += intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms).potentialEnergy();
    }
    torsion_orientations[j] = {std::move(rotated_atoms), torsion_energy};
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

  return {std::move(torsion_orientations[selected_torsion].first),
          rosenbluth_weight_torsion / static_cast<double>(numberOfTorsionTrials)};
}
