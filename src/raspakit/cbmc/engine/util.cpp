module;

module cbmc_util;

import std;

import atom;
import double3;
import randomnumbers;
import cbmc_grow_step;

void CBMC::placeStepBeads(std::vector<Atom> &chainAtoms, const GrowStep &step, std::span<const Atom> positions)
{
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) chainAtoms[step.nextBeads[k]] = positions[k];
}

std::vector<Atom> CBMC::stepBeadPositions(std::span<const Atom> chainAtoms, const GrowStep &step)
{
  std::vector<Atom> positions(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) positions[k] = chainAtoms[step.nextBeads[k]];
  return positions;
}

// logBoltzmannFactors are (-beta U)
std::size_t CBMC::selectTrialPosition(RandomNumber &random, std::span<const double> logBoltzmannFactors)
{
  // Energies are always bounded from below [-U_max, infinity>
  // Find the lowest energy value, i.e. the largest value of (-beta U)
  auto match = std::max_element(logBoltzmannFactors.begin(), logBoltzmannFactors.end());

  if (match == logBoltzmannFactors.end())
  {
    throw std::runtime_error("[cbmc-utils]: no maximum value found\n");
  }
  const double largest_value = *match;

  // Standard trick: shift the Boltzmann factors down to avoid numerical problems
  // The largest shifted factor is 1 (the lowest energy).
  double sumShiftedBoltzmannFactors = 0.0;
  for (double logFactor : logBoltzmannFactors)
  {
    sumShiftedBoltzmannFactors += std::exp(logFactor - largest_value);
  }

  // select the Boltzmann factor (the shifted factors are recomputed on the fly; the exp is cheaper than
  // the heap allocation this used to make per selection)
  std::size_t selected = 0;
  double cumw = std::exp(logBoltzmannFactors[0] - largest_value);
  const double ws = random.uniform() * sumShiftedBoltzmannFactors;
  while (selected + 1 < logBoltzmannFactors.size() && cumw < ws)
  {
    cumw += std::exp(logBoltzmannFactors[++selected] - largest_value);
  }

  return selected;
}

double CBMC::chiralSignedVolume(const std::array<std::size_t, 4> &ids, std::span<const Atom> atoms)
{
  double3 p0 = atoms[ids[0]].position;
  double3 d1 = atoms[ids[1]].position - p0;
  double3 d2 = atoms[ids[2]].position - p0;
  double3 d3 = atoms[ids[3]].position - p0;
  return double3::dot(d1, double3::cross(d2, d3));
}
