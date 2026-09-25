module;

module cbmc_util;

import std;

import atom;
import double3;
import randomnumbers;

// LogBoltzmannFactors are (-Beta U)
std::size_t CBMC::selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors)
{
  std::vector<double> ShiftedBoltzmannFactors(LogBoltzmannFactors.size());

  // Energies are always bounded from below [-U_max, infinity>
  // Find the lowest energy value, i.e. the largest value of (-Beta U)
  std::vector<double>::iterator match = std::max_element(LogBoltzmannFactors.begin(), LogBoltzmannFactors.end());

  if (match == LogBoltzmannFactors.end())
  {
    throw std::runtime_error("[cbmc-utils]: no maximum value found\n");
  }
  double largest_value = *match;

  // Standard trick: shift the Boltzmann factors down to avoid numerical problems
  // The largest value of 'ShiftedBoltzmannFactors' will be 1 (which corresponds to the lowest energy).
  double SumShiftedBoltzmannFactors = 0.0;
  for (std::size_t i = 0; i < LogBoltzmannFactors.size(); ++i)
  {
    ShiftedBoltzmannFactors[i] = std::exp(LogBoltzmannFactors[i] - largest_value);
    SumShiftedBoltzmannFactors += ShiftedBoltzmannFactors[i];
  }

  // select the Boltzmann factor
  std::size_t selected = 0;
  double cumw = ShiftedBoltzmannFactors[0];
  double ws = random.uniform() * SumShiftedBoltzmannFactors;
  while (selected + 1 < ShiftedBoltzmannFactors.size() && cumw < ws)
  {
    cumw += ShiftedBoltzmannFactors[++selected];
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
