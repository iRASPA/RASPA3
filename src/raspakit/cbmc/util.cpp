module;

module cbmc_util;

import std;

import atom;
import double3x3;
import double3;
import simd_quatd;
import randomnumbers;
import stringutils;

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
