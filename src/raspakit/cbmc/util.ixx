module;

export module cbmc_util;

import std;

import atom;
import double3x3;
import double3;
import simd_quatd;
import randomnumbers;

export namespace CBMC
{
// LogBoltzmannFactors are (-Beta U)
std::size_t selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors);
}  // namespace CBMC
