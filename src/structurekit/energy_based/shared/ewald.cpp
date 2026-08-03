module;

module energy_shared_ewald;

import std;

import double3;
import double3x3;
import unit_cell;
import crystal;
import units;

EwaldSplit ewaldSplit(double cutOff, double relativePrecision, std::optional<double> alphaOverride)
{
  double decay = std::sqrt(-std::log(relativePrecision));

  EwaldSplit split;
  split.alpha = alphaOverride.value_or(decay / cutOff);
  split.largestWaveVector = 2.0 * split.alpha * decay;

  return split;
}

std::vector<WaveVector> waveVectors(const Crystal &framework, const std::vector<double3> &fractionalPositions,
                                    double alpha, double largestWaveVector)
{
  std::vector<WaveVector> vectors;

  double3x3 cell = framework.unitCell.cell;
  double3x3 reciprocal = framework.unitCell.inverseCell.transpose();  // rows are the reciprocal vectors over 2*pi

  // A wave vector's index along an axis is bounded by its length times that axis's length, over 2*pi, whatever
  // the shape of the cell, because the dot product of the two is exactly 2*pi times the index.
  auto bound = [&](std::size_t axis)
  {
    return static_cast<std::int32_t>(std::floor(largestWaveVector * cell[axis].length() / (2.0 * std::numbers::pi)));
  };
  std::int32_t hMax = bound(0), kMax = bound(1), lMax = bound(2);

  const double largestSquared = largestWaveVector * largestWaveVector;
  const double alphaFactor = 1.0 / (4.0 * alpha * alpha);
  const double volume = framework.unitCell.volume;

  for (std::int32_t h = 0; h <= hMax; ++h)
  {
    for (std::int32_t k = (h == 0) ? 0 : -kMax; k <= kMax; ++k)
    {
      // Half of reciprocal space, taking the first non-zero index to be positive, so that no vector and its
      // opposite are both present and the origin is left out.
      std::int32_t lStart = (h == 0 && k == 0) ? 1 : -lMax;
      for (std::int32_t l = lStart; l <= lMax; ++l)
      {
        double3 waveVector = 2.0 * std::numbers::pi *
                             (static_cast<double>(h) * reciprocal[0] + static_cast<double>(k) * reciprocal[1] +
                              static_cast<double>(l) * reciprocal[2]);

        double lengthSquared = waveVector.length_squared();
        if (lengthSquared > largestSquared || lengthSquared <= 0.0) continue;

        // The weight each wave carries: how much of it survives the split, divided by its length squared.
        // The factor of two is the opposite vector this one stands in for.
        double weight =
            2.0 * (4.0 * std::numbers::pi / volume) * std::exp(-lengthSquared * alphaFactor) / lengthSquared;

        double real = 0.0;
        double imaginary = 0.0;
        for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
        {
          double charge = framework.atoms[j].charge;
          if (charge == 0.0) continue;

          double phase = 2.0 * std::numbers::pi *
                         (static_cast<double>(h) * fractionalPositions[j].x +
                          static_cast<double>(k) * fractionalPositions[j].y +
                          static_cast<double>(l) * fractionalPositions[j].z);
          real += charge * std::cos(phase);
          imaginary -= charge * std::sin(phase);
        }

        vectors.push_back(WaveVector{h, k, l, weight * real, weight * imaginary});
      }
    }
  }

  return vectors;
}

double neutralisingBackground(double netCharge, double alpha, double volume)
{
  return -Units::CoulombicConversionFactor * std::numbers::pi * netCharge / (alpha * alpha * volume);
}
