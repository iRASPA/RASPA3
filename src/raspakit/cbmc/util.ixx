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
std::vector<Atom> rotateRandomlyAround(RandomNumber &random, std::span<Atom> atoms, std::size_t startingBead);
std::vector<Atom> rotateRandomlyAround(simd_quatd &q, std::span<Atom> atoms, std::size_t startingBead);

// LogBoltzmannFactors are (-Beta U)
std::size_t selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors);
}  // namespace CBMC
