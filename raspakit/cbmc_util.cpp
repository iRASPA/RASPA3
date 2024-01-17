module cbmc_util;

import <vector>;

import atom;
import double3x3;
import double3;
import randomnumbers;


std::vector<Atom> CBMC::rotateRandomlyAround(RandomNumber &random, std::vector<Atom> atoms, size_t startingBead)
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> randomlyRotatedAtoms{};
  for (size_t i = 0; i < atoms.size(); ++i)
  {
    Atom b = atoms[i];
    b.position = atoms[startingBead].position + randomRotationMatrix * (b.position - atoms[startingBead].position);
    randomlyRotatedAtoms.push_back(b);
  }
  return randomlyRotatedAtoms;
}


// LogBoltzmannFactors are (-Beta U)
size_t CBMC::selectTrialPosition(RandomNumber &random, std::vector <double> LogBoltzmannFactors) noexcept
{
  std::vector<double> ShiftedBoltzmannFactors(LogBoltzmannFactors.size());

  // Energies are always bounded from below [-U_max, infinity>
  // Find the lowest energy value, i.e. the largest value of (-Beta U)
  double largest_value = *std::max_element(LogBoltzmannFactors.begin(), LogBoltzmannFactors.end());

  // Standard trick: shift the Boltzmann factors down to avoid numerical problems
  // The largest value of 'ShiftedBoltzmannFactors' will be 1 (which corresponds to the lowest energy).
  double SumShiftedBoltzmannFactors = 0.0;
  for (size_t i = 0; i < LogBoltzmannFactors.size(); ++i)
  {
    ShiftedBoltzmannFactors[i] = exp(LogBoltzmannFactors[i] - largest_value);
    SumShiftedBoltzmannFactors += ShiftedBoltzmannFactors[i];
  }

  // select the Boltzmann factor
  size_t selected = 0;
  double cumw = ShiftedBoltzmannFactors[0];
  double ws = random.uniform() * SumShiftedBoltzmannFactors;
  while (cumw < ws)
    cumw += ShiftedBoltzmannFactors[++selected];

  return selected;
}
