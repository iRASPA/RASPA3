module;

export module cbmc_util;

import std;

import atom;
import randomnumbers;

export namespace CBMC
{
/// Rosenbluth selection among trial directions given their log Boltzmann factors (-beta U): returns
/// the index of the selected trial, drawn with probability proportional to exp(logBoltzmannFactor).
std::size_t selectTrialPosition(RandomNumber &random, std::vector<double> LogBoltzmannFactors);

/// Signed volume of the tetrahedron spanned by the four atoms of a chiral center; its sign is the
/// center's parity.
double chiralSignedVolume(const std::array<std::size_t, 4> &ids, std::span<const Atom> atoms);
}  // namespace CBMC
