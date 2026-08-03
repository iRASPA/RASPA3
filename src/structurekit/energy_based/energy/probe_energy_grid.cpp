module;

module energy_probe_energy_grid;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import units;

import energy_shared_probe_energy_grid;

ProbeEnergyGrid ProbeEnergyGridCPU::compute(const PairInteractions &interactions, const Crystal &framework,
                                            std::string probePseudoAtom, uint3 gridSize)
{
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("ProbeEnergyGridCPU: the grid must have at least one point along each axis\n");
  }

  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("ProbeEnergyGridCPU: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  ProbeEnergyGrid grid;
  grid.gridSize = gridSize;
  grid.backend = "cpu";
  grid.unitCell = framework.unitCell;
  grid.probeName = probePseudoAtom;
  grid.probeEpsilon = interactions[probeType.value()].strengthParameter;
  grid.probeSigma = interactions[probeType.value()].sizeParameter;
  grid.cutOff = interactions.cutOffVDW;
  grid.ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.energy.assign(numberOfVoxels, 0.0f);
  grid.strongestAtom.assign(numberOfVoxels, -1);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  const std::size_t numberOfAtoms = fractionalPositions.size();

  // The mixing is the force field's own, taken for the pair rather than assembled here, so that the field
  // agrees with what a simulation with this probe would feel.
  std::vector<double> epsilonTimesFour(numberOfAtoms);
  std::vector<double> sigmaSquared(numberOfAtoms);
  std::vector<double> shiftValue(numberOfAtoms);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    std::size_t atomType = framework.atoms[i].type;
    double sigma = interactions(probeType.value(), atomType).sizeParameter;
    double epsilon = interactions(probeType.value(), atomType).strengthParameter;

    epsilonTimesFour[i] = 4.0 * epsilon;
    sigmaSquared[i] = sigma * sigma;

    // A truncated potential is usually shifted so that it reaches the cutoff at zero rather than stepping
    // there, and the force field says whether it is. Leaving the step in puts the energy out by however many
    // neighbours are in range, which varies from place to place and so bends the landscape rather than
    // raising it.
    shiftValue[i] = interactions(probeType.value(), atomType).shift;
  }

  double3x3 cell = framework.unitCell.cell;

  // A fractional difference t is at least |t_k| w_k long in Cartesian terms, w_k being the width of the cell
  // perpendicular to the plane of the other two axes. Reducing into [-1/2, 1/2) leaves |t_k| at least
  // |n_k| - 1/2 for an image n_k cells out, so an image can only come within the cutoff while
  // |n_k| <= cutoff / w_k + 1/2, and that settles how far to look without any reference to the cell's shape.
  double3 widths = framework.unitCell.perpendicularWidths();
  grid.numberOfImageShells = int3(static_cast<std::int32_t>(std::floor(grid.cutOff / widths.x + 0.5)),
                                  static_cast<std::int32_t>(std::floor(grid.cutOff / widths.y + 0.5)),
                                  static_cast<std::int32_t>(std::floor(grid.cutOff / widths.z + 0.5)));

  const int3 shells = grid.numberOfImageShells;
  const double cutOff = grid.cutOff;
  const double cutOffSquared = cutOff * cutOff;
  const double ceiling = grid.ceiling;

#pragma omp parallel for schedule(dynamic)
  for (std::int64_t iz = 0; iz < static_cast<std::int64_t>(gridSize.z); ++iz)
  {
    for (std::size_t iy = 0; iy < gridSize.y; ++iy)
    {
      for (std::size_t ix = 0; ix < gridSize.x; ++ix)
      {
        // Endpoint-exclusive sampling, as in the clearance field: fractional 0 and 1 are the same periodic
        // point, so dividing by the grid size rather than one less keeps the spacing uniform and every
        // sample distinct. The two fields have to be sampled alike for anything to be read off both.
        double3 s(static_cast<double>(ix) / static_cast<double>(gridSize.x),
                  static_cast<double>(iy) / static_cast<double>(gridSize.y),
                  static_cast<double>(iz) / static_cast<double>(gridSize.z));

        double total = 0.0;

        // The atom pulling hardest, and the nearest one in case none of them pulls. An atom's claim is its
        // own contribution to the sum, over all of its images within the cutoff, so that an atom is not
        // credited twice for being near in a small cell.
        //
        // Inside a wall every term is held at the ceiling, which makes the atoms there exactly equal and the
        // strongest of them a matter of which was looked at first. Distance separates them and energy does
        // not, so that region falls back on the nearest atom, which is also what the geometric route would
        // say about it.
        double strongestPull = 0.0;
        double nearestDistanceSquared = std::numeric_limits<double>::max();
        std::int32_t pullingAtom = -1;
        std::int32_t nearestAtom = -1;

        for (std::size_t iatom = 0; iatom < numberOfAtoms; ++iatom)
        {
          double contribution = 0.0;
          double closest = std::numeric_limits<double>::max();
          double3 ds = s - fractionalPositions[iatom];
          ds.x -= std::rint(ds.x);
          ds.y -= std::rint(ds.y);
          ds.z -= std::rint(ds.z);

          const double epsilon4 = epsilonTimesFour[iatom];
          const double sigma2 = sigmaSquared[iatom];
          const double shift = shiftValue[iatom];

          for (std::int32_t a = -shells.x; a <= shells.x; ++a)
          {
            for (std::int32_t b = -shells.y; b <= shells.y; ++b)
            {
              for (std::int32_t c = -shells.z; c <= shells.z; ++c)
              {
                double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));

                // The same inequality that sizes the shells above, used again to throw out an image before it
                // is built.
                double reach = std::max({std::abs(t.x) * widths.x, std::abs(t.y) * widths.y,
                                         std::abs(t.z) * widths.z});
                if (reach > cutOff) continue;

                double3 dr = cell * t;
                double rr = double3::dot(dr, dr);
                if (rr >= cutOffSquared) continue;

                // A grid point can land on an atom's centre, where the pair energy has no value. The
                // separation is held off zero so that the sum stays a number; the point is buried far above
                // the ceiling either way, so nothing read off the field can turn on the figure chosen here.
                rr = std::max(rr, 1.0e-6);
                closest = std::min(closest, rr);

                double ratio = sigma2 / rr;
                double ratio3 = ratio * ratio * ratio;

                // Each term is held down before it is added rather than the sum afterwards, so that no term
                // ever overflows to something a later addition would turn into a value that is not a number.
                contribution += std::min(epsilon4 * ratio3 * (ratio3 - 1.0) - shift, ceiling);
              }
            }
          }

          total += contribution;

          if (contribution < strongestPull)
          {
            strongestPull = contribution;
            pullingAtom = static_cast<std::int32_t>(iatom);
          }
          if (closest < nearestDistanceSquared)
          {
            nearestDistanceSquared = closest;
            nearestAtom = static_cast<std::int32_t>(iatom);
          }
        }

        const std::size_t voxel = (static_cast<std::size_t>(iz) * gridSize.y + iy) * gridSize.x + ix;
        grid.energy[voxel] = static_cast<float>(std::min(total, ceiling));
        grid.strongestAtom[voxel] = (pullingAtom >= 0) ? pullingAtom : nearestAtom;
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}
