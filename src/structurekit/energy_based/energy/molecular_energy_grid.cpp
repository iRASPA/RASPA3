module;

module energy_molecular_energy_grid;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import units;

import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

namespace
{
// The far half of the electrostatic sum, read off between grid points. It is built from waves no shorter than
// a few Ångström, so it bends very little over one spacing and reading it off straight lines costs almost
// nothing. The near half, which has the 1/r in it and bends a great deal, is never read off but summed
// exactly below.
double smoothPotentialAt(std::span<const float> smoothPotential, double3 s, uint3 gridSize)
{
  s.x -= std::floor(s.x);
  s.y -= std::floor(s.y);
  s.z -= std::floor(s.z);

  double3 scaled(s.x * static_cast<double>(gridSize.x), s.y * static_cast<double>(gridSize.y),
                 s.z * static_cast<double>(gridSize.z));

  std::size_t lowX = static_cast<std::size_t>(std::floor(scaled.x)) % gridSize.x;
  std::size_t lowY = static_cast<std::size_t>(std::floor(scaled.y)) % gridSize.y;
  std::size_t lowZ = static_cast<std::size_t>(std::floor(scaled.z)) % gridSize.z;

  double fracX = scaled.x - std::floor(scaled.x);
  double fracY = scaled.y - std::floor(scaled.y);
  double fracZ = scaled.z - std::floor(scaled.z);

  std::size_t highX = (lowX + 1) % gridSize.x;
  std::size_t highY = (lowY + 1) % gridSize.y;
  std::size_t highZ = (lowZ + 1) % gridSize.z;

  auto at = [&](std::size_t i, std::size_t j, std::size_t k)
  { return static_cast<double>(smoothPotential[(k * gridSize.y + j) * gridSize.x + i]); };
  auto mix = [](double a, double b, double t) { return a + t * (b - a); };

  double c00 = mix(at(lowX, lowY, lowZ), at(highX, lowY, lowZ), fracX);
  double c10 = mix(at(lowX, highY, lowZ), at(highX, highY, lowZ), fracX);
  double c01 = mix(at(lowX, lowY, highZ), at(highX, lowY, highZ), fracX);
  double c11 = mix(at(lowX, highY, highZ), at(highX, highY, highZ), fracX);

  return mix(mix(c00, c10, fracY), mix(c01, c11, fracY), fracZ);
}
}  // namespace

MolecularEnergyGrid MolecularEnergyGridCPU::compute(const PairInteractions &interactions, const Crystal &framework,
                                                    const LinearProbe &probe, uint3 gridSize,
                                                    std::size_t numberOfOrientations, double temperature,
                                                    const ElectrostaticPotentialGrid *potential)
{
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: the grid must have at least one point along each axis\n");
  }
  if (probe.sites.empty())
  {
    throw std::runtime_error("MolecularEnergyGridCPU: the probe has no sites\n");
  }
  if (numberOfOrientations == 0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: at least one orientation is needed\n");
  }
  if (temperature <= 0.0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: the free energy needs a temperature above zero\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  MolecularEnergyGrid grid;
  grid.gridSize = gridSize;
  grid.backend = "cpu";
  grid.unitCell = framework.unitCell;
  grid.probe = probe;
  grid.overHemisphere = probe.headTailSymmetric;
  grid.temperature = temperature;
  grid.cutOff = interactions.cutOffVDW;
  grid.ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  // Electrostatics are acted on only when the molecule has charges to act on and a potential was supplied to
  // act with. A charged molecule without one is the case worth being loud about, so it is recorded rather
  // than quietly treated as neutral.
  grid.chargesIncluded = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;
  grid.chargesIgnored = probe.isCharged() && !grid.chargesIncluded;
  if (grid.chargesIncluded)
  {
    grid.ewaldAlpha = potential->alpha;
    grid.numberOfWaveVectors = potential->numberOfWaveVectors;
    if (potential->gridSize.x != gridSize.x || potential->gridSize.y != gridSize.y ||
        potential->gridSize.z != gridSize.z)
    {
      throw std::runtime_error("MolecularEnergyGridCPU: the potential is on a different grid than this one\n");
    }
  }

  // A molecule that is the same end for end gives the same energy for a direction and its opposite, so half
  // the sphere holds every distinct orientation and sampling the whole of it would only duplicate.
  std::vector<double3> directions = orientationSet(numberOfOrientations, probe.headTailSymmetric);
  grid.numberOfOrientations = directions.size();

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.freeEnergy.assign(numberOfVoxels, 0.0f);
  grid.minimumEnergy.assign(numberOfVoxels, 0.0f);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  const std::size_t numberOfAtoms = fractionalPositions.size();
  const std::size_t numberOfSites = probe.sites.size();
  const std::size_t numberOfDirections = directions.size();

  double3x3 cell = framework.unitCell.cell;
  double3x3 inverseCell = framework.unitCell.inverseCell;

  // Where each site sits relative to the centre of mass, as a fractional displacement. It depends on the
  // orientation and the site but not on the grid point, so it is built once here.
  std::vector<double3> siteOffsets(numberOfDirections * numberOfSites);
  for (std::size_t o = 0; o < numberOfDirections; ++o)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      siteOffsets[o * numberOfSites + s] = inverseCell * (probe.sites[s].offset * directions[o]);
    }
  }

  // The mixed size and strength for every site against every framework atom, taken from the force field so
  // that its own rule applies rather than one repeated here.
  std::vector<double> epsilonTimesFour(numberOfSites * numberOfAtoms);
  std::vector<double> sigmaSquared(numberOfSites * numberOfAtoms);
  std::vector<double> chargeProduct(numberOfSites * numberOfAtoms, 0.0);
  std::vector<double> shiftValue(numberOfSites * numberOfAtoms, 0.0);
  std::vector<double> siteCharge(numberOfSites);
  for (std::size_t s = 0; s < numberOfSites; ++s) siteCharge[s] = probe.sites[s].charge;
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    std::size_t atomType = framework.atoms[i].type;
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double sigma = interactions(probe.sites[s].type, atomType).sizeParameter;
      double epsilon = interactions(probe.sites[s].type, atomType).strengthParameter;
      epsilonTimesFour[s * numberOfAtoms + i] = 4.0 * epsilon;
      sigmaSquared[s * numberOfAtoms + i] = sigma * sigma;

      // The near half of the lattice sum is done pair by pair alongside the dispersion, since the walk over
      // neighbours is already being made and an erfc is all that is added to it.
      chargeProduct[s * numberOfAtoms + i] =
          Units::CoulombicConversionFactor * probe.sites[s].charge * framework.atoms[i].charge;

      // The force field's own truncation convention, applied here as everywhere else in this route.
      shiftValue[s * numberOfAtoms + i] = interactions(probe.sites[s].type, atomType).shift;
    }
  }

  // The molecule reaches beyond its own centre, so an atom that is out of range of the centre may still be
  // within range of a site. The reach is widened by the half-length of the molecule before the image range is
  // settled, which keeps the same inequality valid for every site at once.
  double coulombCutOff = grid.chargesIncluded ? potential->cutOff : 0.0;
  double longestCutOff = std::max(grid.cutOff, coulombCutOff);

  double3 widths = framework.unitCell.perpendicularWidths();
  double reach = longestCutOff + 0.5 * probe.length();
  grid.numberOfImageShells = int3(static_cast<std::int32_t>(std::floor(reach / widths.x + 0.5)),
                                  static_cast<std::int32_t>(std::floor(reach / widths.y + 0.5)),
                                  static_cast<std::int32_t>(std::floor(reach / widths.z + 0.5)));

  const int3 shells = grid.numberOfImageShells;
  const double cutOff = grid.cutOff;
  const double cutOffSquared = cutOff * cutOff;
  const double coulombCutOffSquared = coulombCutOff * coulombCutOff;
  const double longestCutOffSquared = longestCutOff * longestCutOff;
  const double ceiling = grid.ceiling;
  const double kT = Units::KB * temperature;
  const double beta = 1.0 / kT;
  const double alpha = grid.chargesIncluded ? potential->alpha : 0.0;
  const bool useCharges = grid.chargesIncluded;
  std::span<const float> smoothPotential =
      useCharges ? std::span<const float>(potential->smoothPotential) : std::span<const float>{};

#pragma omp parallel for schedule(dynamic)
  for (std::int64_t iz = 0; iz < static_cast<std::int64_t>(gridSize.z); ++iz)
  {
    for (std::size_t iy = 0; iy < gridSize.y; ++iy)
    {
      for (std::size_t ix = 0; ix < gridSize.x; ++ix)
      {
        // Endpoint-exclusive sampling, as in the other fields of this route, so that what is read off one may
        // be read off another. The point is the molecule's centre of mass; its sites hang off it.
        double3 centre(static_cast<double>(ix) / static_cast<double>(gridSize.x),
                       static_cast<double>(iy) / static_cast<double>(gridSize.y),
                       static_cast<double>(iz) / static_cast<double>(gridSize.z));

        // The average of exp(-U/kT) over orientations is accumulated as it goes, rescaled whenever a deeper
        // orientation turns up. Taken plainly the exponential would overflow at a strong site long before the
        // physics did; carried this way every term is at most one and the least is harmless.
        double least = std::numeric_limits<double>::max();
        double sum = 0.0;

        for (std::size_t o = 0; o < numberOfDirections; ++o)
        {
          double total = 0.0;

          for (std::size_t site = 0; site < numberOfSites; ++site)
          {
            double3 s = centre + siteOffsets[o * numberOfSites + site];

            const double *epsilonForSite = epsilonTimesFour.data() + site * numberOfAtoms;
            const double *sigmaForSite = sigmaSquared.data() + site * numberOfAtoms;
            const double *chargeForSite = chargeProduct.data() + site * numberOfAtoms;
            const double *shiftForSite = shiftValue.data() + site * numberOfAtoms;

            if (useCharges) total += siteCharge[site] * smoothPotentialAt(smoothPotential, s, gridSize);

            for (std::size_t iatom = 0; iatom < numberOfAtoms; ++iatom)
            {
              double3 ds = s - fractionalPositions[iatom];
              ds.x -= std::rint(ds.x);
              ds.y -= std::rint(ds.y);
              ds.z -= std::rint(ds.z);

              const double epsilon4 = epsilonForSite[iatom];
              const double sigma2 = sigmaForSite[iatom];
              const double charges = chargeForSite[iatom];

              for (std::int32_t a = -shells.x; a <= shells.x; ++a)
              {
                for (std::int32_t b = -shells.y; b <= shells.y; ++b)
                {
                  for (std::int32_t c = -shells.z; c <= shells.z; ++c)
                  {
                    double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));

                    // A fractional difference is at least |t_k| w_k long whatever the shape of the cell, so an
                    // image beyond the cutoff is thrown out here before it is built.
                    double far = std::max({std::abs(t.x) * widths.x, std::abs(t.y) * widths.y,
                                           std::abs(t.z) * widths.z});
                    if (far > cutOff) continue;

                    double3 dr = cell * t;
                    double rr = double3::dot(dr, dr);
                    if (rr >= longestCutOffSquared) continue;

                    // A site can land on an atom's centre, where the pair energy has no value. Holding the
                    // separation off zero keeps the sum a number; such an orientation is buried far above the
                    // ceiling and contributes nothing to the average either way.
                    rr = std::max(rr, 1.0e-6);

                    if (rr < cutOffSquared)
                    {
                      double ratio = sigma2 / rr;
                      double ratio3 = ratio * ratio * ratio;
                      total += std::min(epsilon4 * ratio3 * (ratio3 - 1.0) - shiftForSite[iatom], ceiling);
                    }

                    // The near half of the electrostatic sum. The dispersion has already walked to this
                    // neighbour and worked out how far away it is, so all this costs is the erfc.
                    if (useCharges && rr < coulombCutOffSquared && charges != 0.0)
                    {
                      double r = std::sqrt(rr);
                      total += std::min(charges * std::erfc(alpha * r) / r, ceiling);
                    }
                  }
                }
              }
            }
          }

          total = std::min(total, ceiling);

          if (total < least)
          {
            // A deeper orientation resets what the sum is measured from, so everything gathered so far is
            // scaled down to the new floor before this one is added at its full weight.
            sum = sum * std::exp(-beta * (least - total)) + 1.0;
            least = total;
          }
          else
          {
            sum += std::exp(-beta * (total - least));
          }
        }

        std::size_t index = (static_cast<std::size_t>(iz) * gridSize.y + iy) * gridSize.x + ix;
        grid.minimumEnergy[index] = static_cast<float>(least);

        // -kT ln <exp(-U/kT)>, written from the floor upwards. The average can only be at most one when
        // measured from the least, so the logarithm is never positive and the free energy never falls below
        // the minimum: the gap between them is what the molecule pays for having to be turned a particular
        // way.
        grid.freeEnergy[index] = static_cast<float>(
            std::min(least - kT * std::log(sum / static_cast<double>(numberOfDirections)), ceiling));
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}

std::vector<std::int32_t> MolecularEnergyGridCPU::strongestAtoms(const PairInteractions &interactions,
                                                                 const Crystal &framework, const LinearProbe &probe,
                                                                 uint3 gridSize, std::size_t numberOfOrientations,
                                                                 double temperature,
                                                                 const ElectrostaticPotentialGrid *potential)
{
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: the grid must have at least one point along each axis\n");
  }
  if (probe.sites.empty()) throw std::runtime_error("MolecularEnergyGridCPU: the probe has no sites\n");
  if (numberOfOrientations == 0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: at least one orientation is needed\n");
  }
  if (temperature <= 0.0)
  {
    throw std::runtime_error("MolecularEnergyGridCPU: the weighting over orientations needs a temperature\n");
  }

  const std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  std::vector<std::int32_t> labels(numberOfVoxels, -1);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return labels;

  const bool useCharges = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;

  std::vector<double3> directions = orientationSet(numberOfOrientations, probe.headTailSymmetric);

  const std::size_t numberOfAtoms = fractionalPositions.size();
  const std::size_t numberOfSites = probe.sites.size();
  const std::size_t numberOfDirections = directions.size();

  double3x3 cell = framework.unitCell.cell;
  double3x3 inverseCell = framework.unitCell.inverseCell;

  std::vector<double3> siteOffsets(numberOfDirections * numberOfSites);
  for (std::size_t o = 0; o < numberOfDirections; ++o)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      siteOffsets[o * numberOfSites + s] = inverseCell * (probe.sites[s].offset * directions[o]);
    }
  }

  std::vector<double> epsilonTimesFour(numberOfSites * numberOfAtoms);
  std::vector<double> sigmaSquared(numberOfSites * numberOfAtoms);
  std::vector<double> chargeProduct(numberOfSites * numberOfAtoms, 0.0);
  std::vector<double> shiftValue(numberOfSites * numberOfAtoms, 0.0);
  std::vector<double> siteCharge(numberOfSites);
  for (std::size_t s = 0; s < numberOfSites; ++s) siteCharge[s] = probe.sites[s].charge;
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    std::size_t atomType = framework.atoms[i].type;
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double sigma = interactions(probe.sites[s].type, atomType).sizeParameter;
      double epsilon = interactions(probe.sites[s].type, atomType).strengthParameter;
      epsilonTimesFour[s * numberOfAtoms + i] = 4.0 * epsilon;
      sigmaSquared[s * numberOfAtoms + i] = sigma * sigma;
      chargeProduct[s * numberOfAtoms + i] =
          Units::CoulombicConversionFactor * probe.sites[s].charge * framework.atoms[i].charge;
      shiftValue[s * numberOfAtoms + i] = interactions(probe.sites[s].type, atomType).shift;
    }
  }

  const double coulombCutOff = useCharges ? potential->cutOff : 0.0;
  const double cutOff = interactions.cutOffVDW;
  const double longestCutOff = std::max(cutOff, coulombCutOff);

  double3 widths = framework.unitCell.perpendicularWidths();
  double reach = longestCutOff + 0.5 * probe.length();
  const int3 shells(static_cast<std::int32_t>(std::floor(reach / widths.x + 0.5)),
                    static_cast<std::int32_t>(std::floor(reach / widths.y + 0.5)),
                    static_cast<std::int32_t>(std::floor(reach / widths.z + 0.5)));

  const double cutOffSquared = cutOff * cutOff;
  const double coulombCutOffSquared = coulombCutOff * coulombCutOff;
  const double longestCutOffSquared = longestCutOff * longestCutOff;
  const double ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;
  const double beta = 1.0 / (Units::KB * temperature);
  const double alpha = useCharges ? potential->alpha : 0.0;
  std::span<const float> smoothPotential =
      useCharges ? std::span<const float>(potential->smoothPotential) : std::span<const float>{};

  // What one framework atom, over all of its images in range, does to the molecule held one way. It is the
  // part of the energy that can be laid at one atom's door, and the label is decided on it.
  auto contributionOf = [&](double3 centre, std::size_t orientation, std::size_t iatom, double &closestSquared)
  {
    double contribution = 0.0;
    closestSquared = std::numeric_limits<double>::max();

    for (std::size_t site = 0; site < numberOfSites; ++site)
    {
      double3 s = centre + siteOffsets[orientation * numberOfSites + site];

      double3 ds = s - fractionalPositions[iatom];
      ds.x -= std::rint(ds.x);
      ds.y -= std::rint(ds.y);
      ds.z -= std::rint(ds.z);

      const double epsilon4 = epsilonTimesFour[site * numberOfAtoms + iatom];
      const double sigma2 = sigmaSquared[site * numberOfAtoms + iatom];
      const double charges = chargeProduct[site * numberOfAtoms + iatom];
      const double shift = shiftValue[site * numberOfAtoms + iatom];

      for (std::int32_t a = -shells.x; a <= shells.x; ++a)
      {
        for (std::int32_t b = -shells.y; b <= shells.y; ++b)
        {
          for (std::int32_t c = -shells.z; c <= shells.z; ++c)
          {
            double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));

            double far =
                std::max({std::abs(t.x) * widths.x, std::abs(t.y) * widths.y, std::abs(t.z) * widths.z});
            if (far > cutOff) continue;

            double3 dr = cell * t;
            double rr = double3::dot(dr, dr);
            if (rr >= longestCutOffSquared) continue;

            rr = std::max(rr, 1.0e-6);
            closestSquared = std::min(closestSquared, rr);

            if (rr < cutOffSquared)
            {
              double ratio = sigma2 / rr;
              double ratio3 = ratio * ratio * ratio;
              contribution += std::min(epsilon4 * ratio3 * (ratio3 - 1.0) - shift, ceiling);
            }

            if (useCharges && rr < coulombCutOffSquared && charges != 0.0)
            {
              double r = std::sqrt(rr);
              contribution += std::min(charges * std::erfc(alpha * r) / r, ceiling);
            }
          }
        }
      }
    }

    return contribution;
  };

#pragma omp parallel
  {
    // How the molecule is turned at this point, weighed. It is the one thing that has to be known before the
    // atoms can be gone through, so it is gathered first and the atoms are visited afterwards.
    std::vector<double> weight(numberOfDirections, 0.0);

#pragma omp for schedule(dynamic)
    for (std::int64_t iz = 0; iz < static_cast<std::int64_t>(gridSize.z); ++iz)
    {
      for (std::size_t iy = 0; iy < gridSize.y; ++iy)
      {
        for (std::size_t ix = 0; ix < gridSize.x; ++ix)
        {
          double3 centre(static_cast<double>(ix) / static_cast<double>(gridSize.x),
                         static_cast<double>(iy) / static_cast<double>(gridSize.y),
                         static_cast<double>(iz) / static_cast<double>(gridSize.z));

          double least = std::numeric_limits<double>::max();
          for (std::size_t o = 0; o < numberOfDirections; ++o)
          {
            double total = 0.0;
            double ignored = 0.0;
            for (std::size_t iatom = 0; iatom < numberOfAtoms; ++iatom)
            {
              total += contributionOf(centre, o, iatom, ignored);
            }

            // The reciprocal half of the electrostatic sum belongs to no one atom, so it takes no part in
            // the label below. It does belong in the weight, which is about how the molecule sits.
            if (useCharges)
            {
              for (std::size_t site = 0; site < numberOfSites; ++site)
              {
                double3 s = centre + siteOffsets[o * numberOfSites + site];
                total += siteCharge[site] * smoothPotentialAt(smoothPotential, s, gridSize);
              }
            }

            weight[o] = std::min(total, ceiling);
            least = std::min(least, weight[o]);
          }

          // Measured from the deepest orientation, so the largest weight is one whatever the site is worth
          // and nothing overflows. Inside a wall every orientation is at the ceiling and the weights come
          // out equal, which is the right answer there: no orientation is preferred because none is possible.
          double normalisation = 0.0;
          for (std::size_t o = 0; o < numberOfDirections; ++o)
          {
            weight[o] = std::exp(-beta * (weight[o] - least));
            normalisation += weight[o];
          }
          for (std::size_t o = 0; o < numberOfDirections; ++o) weight[o] /= normalisation;

          double strongestPull = 0.0;
          double nearestDistanceSquared = std::numeric_limits<double>::max();
          std::int32_t pullingAtom = -1;
          std::int32_t nearestAtom = -1;

          for (std::size_t iatom = 0; iatom < numberOfAtoms; ++iatom)
          {
            double pull = 0.0;
            double closestOverOrientations = std::numeric_limits<double>::max();
            for (std::size_t o = 0; o < numberOfDirections; ++o)
            {
              double closest = 0.0;
              pull += weight[o] * contributionOf(centre, o, iatom, closest);
              closestOverOrientations = std::min(closestOverOrientations, closest);
            }

            if (pull < strongestPull)
            {
              strongestPull = pull;
              pullingAtom = static_cast<std::int32_t>(iatom);
            }
            if (closestOverOrientations < nearestDistanceSquared)
            {
              nearestDistanceSquared = closestOverOrientations;
              nearestAtom = static_cast<std::int32_t>(iatom);
            }
          }

          const std::size_t voxel = (static_cast<std::size_t>(iz) * gridSize.y + iy) * gridSize.x + ix;
          labels[voxel] = (pullingAtom >= 0) ? pullingAtom : nearestAtom;
        }
      }
    }
  }

  return labels;
}
