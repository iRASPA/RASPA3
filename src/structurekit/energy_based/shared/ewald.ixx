module;

export module energy_shared_ewald;

import std;

import double3;
import crystal;

// The pieces of an Ewald sum that neither backend owns: where to put the split, which wave vectors to keep,
// and what a charged cell costs. Both builders take them from here, because two fields built with different
// splits or different wave vector sets are not comparable no matter how well each is computed.

// A wave vector, held by the three whole numbers that index it rather than by its length. Written that way
// the phase at a point is just 2*pi times the dot product with the fractional position, which saves the
// kernel from ever having to leave fractional coordinates.
export struct WaveVector
{
  std::int32_t h, k, l;
  double weightedReal;       // the framework's structure factor at this wave vector, real part, times its weight
  double weightedImaginary;  // and the imaginary part
};

// Where to put the split, and how far out in wave vectors to go. Both follow from asking that each half be
// left with no more than the requested relative error at the edge of where it is summed. The near part dies
// like erfc(alpha r) and the far part like exp(-k^2/4 alpha^2), so the same logarithm sets both.
export struct EwaldSplit
{
  double alpha{0.0};
  double largestWaveVector{0.0};
};

export EwaldSplit ewaldSplit(double cutOff, double relativePrecision, std::optional<double> alphaOverride = {});

// The wave vectors worth keeping, each with the framework's structure factor already folded in. Only half of
// them are kept: a wave vector and its opposite contribute equally, so one stands for both and is counted
// twice at the end.
export std::vector<WaveVector> waveVectors(const Crystal &framework,
                                           const std::vector<double3> &fractionalPositions, double alpha,
                                           double largestWaveVector);

// A cell that does not add up to zero is neutralised by spreading the opposite charge evenly through it,
// which shifts the potential everywhere by the same amount. For a guest that is itself neutral the shift
// cancels between its sites and none of this matters; for one that is not, it does.
export double neutralisingBackground(double netCharge, double alpha, double volume);
