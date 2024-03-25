export module cbmc_util;

import <vector>;
import <cmath>;
import <algorithm>;

import atom;
import double3x3;
import double3;
import randomnumbers;


export namespace CBMC
{
  std::vector<Atom> rotateRandomlyAround(RandomNumber &random, std::vector<Atom> atoms, size_t startingBead);

  // LogBoltzmannFactors are (-Beta U)
  size_t selectTrialPosition(RandomNumber &random, std::vector <double> LogBoltzmannFactors);
}
