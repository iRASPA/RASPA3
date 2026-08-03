module;

module energy_electrostatic_potential_grid;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import units;

import energy_shared_electrostatic_potential_grid;
import energy_shared_ewald;

ElectrostaticPotentialGrid ElectrostaticPotentialGridCPU::compute(const PairInteractions &interactions,
                                                                  const Crystal &framework, uint3 gridSize,
                                                                  double relativePrecision,
                                                                  std::optional<double> alphaOverride)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  ElectrostaticPotentialGrid grid;
  grid.gridSize = gridSize;
  grid.backend = "cpu";
  grid.unitCell = framework.unitCell;
  grid.cutOff = interactions.cutOffCoulomb;
  grid.relativePrecision = relativePrecision;

  EwaldSplit split = ewaldSplit(grid.cutOff, relativePrecision, alphaOverride);
  grid.alpha = split.alpha;
  grid.largestWaveVector = split.largestWaveVector;

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.smoothPotential.assign(numberOfVoxels, 0.0f);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  for (const CrystalAtom &atom : framework.atoms)
  {
    grid.netCharge += atom.charge;
    grid.largestFrameworkCharge = std::max(grid.largestFrameworkCharge, std::abs(atom.charge));
  }

  std::vector<WaveVector> vectors = waveVectors(framework, fractionalPositions, grid.alpha, grid.largestWaveVector);
  grid.numberOfWaveVectors = vectors.size();

  double background = neutralisingBackground(grid.netCharge, grid.alpha, framework.unitCell.volume);

  const std::size_t numberOfWaves = vectors.size();

#pragma omp parallel for schedule(static)
  for (std::int64_t iz = 0; iz < static_cast<std::int64_t>(gridSize.z); ++iz)
  {
    for (std::size_t iy = 0; iy < gridSize.y; ++iy)
    {
      for (std::size_t ix = 0; ix < gridSize.x; ++ix)
      {
        double3 s(static_cast<double>(ix) / static_cast<double>(gridSize.x),
                  static_cast<double>(iy) / static_cast<double>(gridSize.y),
                  static_cast<double>(iz) / static_cast<double>(gridSize.z));

        // Holding the wave vectors as whole numbers of reciprocal cells means the phase is just their dot
        // product with the fractional position, so nothing here needs to know the shape of the cell.
        double total = 0.0;
        for (std::size_t i = 0; i < numberOfWaves; ++i)
        {
          double phase = 2.0 * std::numbers::pi *
                         (static_cast<double>(vectors[i].h) * s.x + static_cast<double>(vectors[i].k) * s.y +
                          static_cast<double>(vectors[i].l) * s.z);
          total += Units::CoulombicConversionFactor *
                   (vectors[i].weightedReal * std::cos(phase) - vectors[i].weightedImaginary * std::sin(phase));
        }

        grid.smoothPotential[(static_cast<std::size_t>(iz) * gridSize.y + iy) * gridSize.x + ix] =
            static_cast<float>(total + background);
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}
