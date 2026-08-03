module;

module energy_shared_electrostatic_potential_grid;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import units;

import energy_shared_ewald;

ElectrostaticPotentialGrid::ElectrostaticPotentialGrid() {}

ElectrostaticPotentialGrid::~ElectrostaticPotentialGrid() {}

double ElectrostaticPotentialGrid::exactPotentialAt(const Crystal &framework, double3 fractionalPosition) const
{
  std::vector<double3> fractionalPositions = framework.fractionalPositions;

  double3x3 cell = framework.unitCell.cell;
  double3 widths = framework.unitCell.perpendicularWidths();
  int3 shells(static_cast<std::int32_t>(std::floor(this->cutOff / widths.x + 0.5)),
              static_cast<std::int32_t>(std::floor(this->cutOff / widths.y + 0.5)),
              static_cast<std::int32_t>(std::floor(this->cutOff / widths.z + 0.5)));

  double near = 0.0;
  for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
  {
    double charge = framework.atoms[j].charge;
    if (charge == 0.0) continue;

    double3 ds = fractionalPosition - fractionalPositions[j];
    ds.x -= std::rint(ds.x);
    ds.y -= std::rint(ds.y);
    ds.z -= std::rint(ds.z);

    for (std::int32_t a = -shells.x; a <= shells.x; ++a)
    {
      for (std::int32_t b = -shells.y; b <= shells.y; ++b)
      {
        for (std::int32_t c = -shells.z; c <= shells.z; ++c)
        {
          double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));
          double3 dr = cell * t;
          double r = dr.length();
          if (r >= this->cutOff || r < 1.0e-8) continue;

          near += charge * std::erfc(this->alpha * r) / r;
        }
      }
    }
  }
  near *= Units::CoulombicConversionFactor;

  std::vector<WaveVector> vectors = waveVectors(framework, fractionalPositions, this->alpha, this->largestWaveVector);

  double far = 0.0;
  for (const WaveVector &vector : vectors)
  {
    double phase = 2.0 * std::numbers::pi *
                   (static_cast<double>(vector.h) * fractionalPosition.x +
                    static_cast<double>(vector.k) * fractionalPosition.y +
                    static_cast<double>(vector.l) * fractionalPosition.z);
    far += vector.weightedReal * std::cos(phase) - vector.weightedImaginary * std::sin(phase);
  }
  far *= Units::CoulombicConversionFactor;

  double background = neutralisingBackground(this->netCharge, this->alpha, framework.unitCell.volume);

  return near + far + background;
}

std::string ElectrostaticPotentialGrid::splitIndependenceCheck(const PairInteractions &interactions, const Crystal &framework)
{
  std::string report;

  // Points scattered through the cell on no particular lattice, so that nothing is tested at a place where
  // some symmetry of the grid might rescue a mistake.
  std::vector<double3> samples{{0.1234, 0.2345, 0.3456}, {0.5000, 0.5000, 0.5000}, {0.7777, 0.1111, 0.9999},
                               {0.3333, 0.6667, 0.2500}, {0.0500, 0.9500, 0.5500}};

  // The split is put in five different places. Every one of them divides the same lattice sum differently,
  // and every one must give the same total back.
  //
  // The choice is not free at both ends. Put the split too far towards the far half and the near half has
  // not died away by the cutoff, so summing it only to there throws away a tail that matters, and the total
  // comes out short. That is a limit of the cutoff rather than a fault in either half, so the splits tried
  // here are all ones the cutoff can carry: alpha times the cutoff from 3.5 upwards leaves less than 1e-7 of
  // the near half beyond it.
  std::vector<double> alphaTimesCutOff{4.0, 4.5, 5.0, 5.5};
  std::vector<double> alphas;
  for (double product : alphaTimesCutOff) alphas.push_back(product / interactions.cutOffCoulomb);

  report += std::format("# Ewald split-independence check for {}\n", framework.name);
  report += "#\n";
  report += "# Where the sum is cut between its near and far halves is a free choice, so the potential must\n";
  report += "# not depend on it. Each column below is the same potential computed with the split somewhere\n";
  report += "# else. They agree only if the two halves are each right, since an error in one is not matched\n";
  report += "# by the same error in the other as alpha moves.\n";
  report += "#\n";
  report += std::format("# Coulomb cutoff: {} [Å], potential in [K] per unit charge\n", interactions.cutOffCoulomb);
  report += std::format("# Splits tried: alpha times the cutoff from {} to {}\n", alphaTimesCutOff.front(),
                        alphaTimesCutOff.back());
  report += "\n";

  report += std::format("{:>34}", "fractional position");
  for (double alpha : alphas) report += std::format("{:>16}", std::format("alpha={:.4f}", alpha));
  report += std::format("{:>16}\n", "spread");

  double worstSpread = 0.0;
  double largestMagnitude = 0.0;

  for (const double3 &sample : samples)
  {
    report += std::format("{:>10.4f}{:>12.4f}{:>12.4f}", sample.x, sample.y, sample.z);

    std::vector<double> values;
    for (double alpha : alphas)
    {
      ElectrostaticPotentialGrid grid;
      grid.cutOff = interactions.cutOffCoulomb;
      grid.alpha = alpha;
      // The far half has to reach far enough that what is dropped from it is no larger than what is dropped
      // from the near half, and where that is depends on where the split was put.
      grid.largestWaveVector = 2.0 * alpha * std::sqrt(-std::log(1.0e-9));
      for (const CrystalAtom &atom : framework.atoms) grid.netCharge += atom.charge;

      double value = grid.exactPotentialAt(framework, sample) * Units::EnergyToKelvin;
      values.push_back(value);
      report += std::format("{:>16.6f}", value);
    }

    auto [low, high] = std::ranges::minmax(values);
    double spread = high - low;
    worstSpread = std::max(worstSpread, spread);
    for (double value : values) largestMagnitude = std::max(largestMagnitude, std::abs(value));

    report += std::format("{:>16.2e}\n", spread);
  }

  report += "\n";
  report += std::format("Largest spread over the splits: {:.4e} [K] per unit charge\n", worstSpread);
  report += std::format("Relative to the largest value:  {:.4e}\n",
                        largestMagnitude > 0.0 ? worstSpread / largestMagnitude : 0.0);
  report += std::format("Net charge in the cell:         {:+.6e}\n", [&]
                        {
                          double total = 0.0;
                          for (const CrystalAtom &atom : framework.atoms) total += atom.charge;
                          return total;
                        }());
  report += "\n";
  report += "# A spread of a part in ten million or better means the two halves add up. Anything larger means\n";
  report += "# one of them is wrong, or that the far half is not being carried far enough for the split used.\n";

  return report;
}
