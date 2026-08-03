module;

module energy_shared_molecular_surface_area;

import std;

import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import grid_area_curve;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_isosurface;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

MolecularSurfaceArea::MolecularSurfaceArea() {}

MolecularSurfaceArea::~MolecularSurfaceArea() {}

MolecularSurfaceArea MolecularSurfaceArea::fromGrid(const EnergyBackend &backend, const MolecularEnergyGrid &grid,
                                                    const Crystal &framework, double isoValue)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  MolecularSurfaceArea result;
  result.grid = grid;
  result.isoValue = isoValue;

  if (grid.numberOfVoxels() == 0) return result;

  // Both fields are laid out with x varying fastest, which is the order the iso-surface extraction reads, so
  // they can be handed over as they are.
  IsosurfaceArea free = backend.isosurfaceArea(framework, grid.freeEnergy, grid.gridSize, isoValue);
  IsosurfaceArea least = backend.isosurfaceArea(framework, grid.minimumEnergy, grid.gridSize, isoValue);

  result.freeEnergyArea = free.area;
  result.minimumEnergyArea = least.area;
  result.numberOfFreeEnergyTriangles = free.numberOfTriangles;
  result.numberOfMinimumEnergyTriangles = least.numberOfTriangles;
  result.numberOfRejectedTriangles = free.numberOfRejectedTriangles + least.numberOfRejectedTriangles;

  if (framework.mass > 0.0)
  {
    double perGram = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass;
    result.gravimetricArea = free.area * perGram;
    result.minimumEnergyGravimetricArea = least.area * perGram;
  }
  if (framework.unitCell.volume > 0.0)
  {
    result.volumetricArea = 1.0e4 * free.area / framework.unitCell.volume;
    result.minimumEnergyVolumetricArea = 1.0e4 * least.area / framework.unitCell.volume;
  }

  std::chrono::duration<double> surfaces = std::chrono::steady_clock::now() - time_begin;
  result.seconds = surfaces.count();
  std::chrono::steady_clock::time_point time_curve = std::chrono::steady_clock::now();

  // The whole family of areas, in one pass over the field, by the coarea formula. The range is set by the
  // landscape rather than chosen: the deepest point the molecule can reach, and as far above zero as that is
  // below it, which covers the wells at one end and the barriers at the other with zero in the middle.
  double deepest = 0.0;
  for (const float value : grid.freeEnergy) deepest = std::min(deepest, static_cast<double>(value));

  if (deepest < 0.0)
  {
    result.curve = areaAgainstLevel(grid.gridSize, framework.unitCell, grid.freeEnergy, deepest, -deepest, 200);

    // Marching cubes at a handful of the same levels. The two are separate constructions, one gathering the
    // gradient over the cell and the other triangulating a surface, so agreement between them is evidence
    // about the grid rather than a tautology.
    for (const double fraction : {-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75})
    {
      const double level = fraction * (-deepest);
      IsosurfaceArea spot = backend.isosurfaceArea(framework, grid.freeEnergy, grid.gridSize, level);
      result.crossCheck.push_back({level, result.curve.areaAt(level), spot.area});
    }

    // The scale, from the one level that was going to be triangulated anyway.
    const double atIsoValue = result.curve.areaAt(isoValue);
    if (atIsoValue > 0.0 && free.area > 0.0) result.anchorScale = free.area / atIsoValue;
  }

  std::chrono::duration<double> checking = std::chrono::steady_clock::now() - time_curve;
  result.crossCheckSeconds = checking.count() - result.curve.seconds;

  return result;
}

void MolecularSurfaceArea::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                               double isoValue, uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                               bool useElectrostatics, double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);

  *this = MolecularSurfaceArea::fromGrid(backend, field.grid, framework, isoValue);
  this->potential = field.potential;

  if (this->grid.numberOfVoxels() == 0) return;

  double3 spacing = this->grid.spacing();
  double density = 1e-3 * framework.mass /
                   (framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                    Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.sa." + this->grid.backend + ".txt");
  std::print(myfile, "# Energy-based surface area for a rigid linear molecule\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", density);
  std::print(myfile, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  std::print(myfile, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end" : "whole sphere");
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             this->grid.gridSize.x, this->grid.gridSize.y, this->grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->grid.cutOff);
  std::print(myfile, "# Iso-value: {} [internal] {:.4f} [K], as passed to the single-site route unchanged\n",
             this->isoValue, this->isoValue * Units::EnergyToKelvin);
  std::print(myfile, "# Triangles: {} free-energy, {} minimum-energy (discarded as implausibly large: {})\n",
             this->numberOfFreeEnergyTriangles, this->numberOfMinimumEnergyTriangles,
             this->numberOfRejectedTriangles);

  if (this->grid.chargesIncluded)
  {
    std::print(myfile, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors\n",
               this->potential.alpha, this->potential.numberOfWaveVectors);
  }
  if (this->grid.chargesIgnored)
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# WARNING: this molecule carries partial charges and they have not been acted on.\n");
  }

  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# Timing ({}): {} [s] for the two iso-surfaces\n", this->grid.backend, this->seconds);
  std::print(myfile, "#\n");
  std::print(myfile, "# Two surfaces of the same landscape. The free-energy one is the surface the molecule\n");
  std::print(myfile, "# meets at this temperature, its orientations averaged over; the minimum-energy one is\n");
  std::print(myfile, "# the surface it would meet if it could always be turned the best way at no cost. The\n");
  std::print(myfile, "# second encloses the first, always, and the gap is the reach a molecule gives up to\n");
  std::print(myfile, "# the freedom it would have to surrender to get there.\n");
  std::print(myfile, "\n");

  std::print(myfile, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "Free energy surface:", this->freeEnergyArea, this->gravimetricArea, this->volumetricArea);
  std::print(myfile, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "Minimum energy surface:", this->minimumEnergyArea, this->minimumEnergyGravimetricArea,
             this->minimumEnergyVolumetricArea);
  std::print(myfile, "{:34} {:14.4f} [-]\n", "Orientational excess:", this->orientationalExcess());
  std::print(myfile, "\n");

  std::print(myfile, "# The area of a level set converges more slowly than anything else read off this\n");
  std::print(myfile, "# landscape, because it measures how crumpled a surface is and so counts noise rather\n");
  std::print(myfile, "# than averaging it away. Both areas should be checked against the grid spacing and the\n");
  std::print(myfile, "# number of orientations before either is quoted. The minimum-energy area is the worse\n");
  std::print(myfile, "# behaved of the two: a least over a discrete set of orientations is kinked wherever the\n");
  std::print(myfile, "# best orientation changes, and those kinks do not smooth away with more of them.\n");

  if (!this->curve.points.empty())
  {
    const double toKelvin = Units::EnergyToKelvin;

    std::print(myfile, "\n");
    std::print(myfile, "# THE AREA AT EVERY LEVEL\n");
    std::print(myfile, "#\n");
    std::print(myfile, "# The two numbers above are areas at one level, and there is no level that is the right\n");
    std::print(myfile, "# one. Everything else on this landscape describes a region, and a region can be pinned\n");
    std::print(myfile, "# to a contour that means something; an area describes a surface, and the surface is the\n");
    std::print(myfile, "# level. Nor is there a Boltzmann weight to be had: the energy is the same on every part\n");
    std::print(myfile, "# of one surface by construction, so weighting the surface by exp(-A/kT) multiplies it\n");
    std::print(myfile, "# by a constant and changes nothing, and the only variation there is across the\n");
    std::print(myfile, "# triangles is interpolation error. So the whole family is given here instead.\n");
    std::print(myfile, "#\n");
    std::print(myfile, "# The temperature does reach the landscape, through the orientational average, and the\n");
    std::print(myfile, "# effect is not small: the free-energy area of a molecule falls steadily as the\n");
    std::print(myfile, "# temperature rises, while the minimum-energy area, which has no temperature in it,\n");
    std::print(myfile, "# does not move at all.\n");
    std::print(myfile, "#\n");
    std::print(myfile, "# This is not marching cubes run two hundred times. It is the coarea formula, which\n");
    std::print(myfile, "# says that gathering the magnitude of the gradient into bins of the field gives the\n");
    std::print(myfile, "# area of every level set at once, in one pass. Each grid point is spread over the\n");
    std::print(myfile, "# range of values it covers rather than dropped into one bin, since a voxel is not a\n");
    std::print(myfile, "# point and treating it as one makes the curve alias.\n");
    std::print(myfile, "#\n");
    std::print(myfile, "# Levels: {} bins of {:.4f} [K] from {:.2f} to {:.2f} [K]\n", this->curve.points.size(),
               this->curve.binWidth * toKelvin, this->curve.lowestLevel * toKelvin,
               this->curve.highestLevel * toKelvin);
    std::print(myfile, "# The field moves {:.4f} [K] across one voxel on average, which is the resolution of\n",
               this->curve.meanSpanPerVoxel * toKelvin);
    std::print(myfile, "# this curve: bins finer than that are correlated with their neighbours, not\n");
    std::print(myfile, "# independent measurements.\n");
    std::print(myfile, "# Outside the range: {:.6f} of the cell below, {:.6f} above\n",
               this->curve.fractionBelowRange, this->curve.fractionAboveRange);
    const double perSurface = 0.5 * this->seconds;
    std::print(myfile, "# Timing ({}): {} [s] for the curve, all {} levels of it\n", this->grid.backend,
               this->curve.seconds, this->curve.points.size());
    std::print(myfile, "# Timing ({}): {} [s] for one iso-surface, so {:.2f} [s] to triangulate them all\n",
               this->grid.backend, perSurface, perSurface * static_cast<double>(this->curve.points.size()));
    std::print(myfile, "# Timing ({}): {} [s] for the {} spot surfaces checked against below\n",
               this->grid.backend, this->crossCheckSeconds, this->crossCheck.size());

    if (!this->crossCheck.empty())
    {
      std::print(myfile, "#\n");
      std::print(myfile, "# Against marching cubes at the same levels. The two are separate constructions, one\n");
      std::print(myfile, "# gathering the gradient over the cell and the other triangulating a surface, and on\n");
      std::print(myfile, "# any grid you would actually use they do not agree. The curve comes out high, and it\n");
      std::print(myfile, "# comes out high for a reason: a central difference taken across a wall that climbs\n");
      std::print(myfile, "# as the twelfth power overstates the gradient whenever the voxel is wide compared to\n");
      std::print(myfile, "# the distance over which the field turns. Refining the grid closes it, but slowly.\n");
      std::print(myfile, "# On MFI with argon the curve runs 32, 18, 7, 4 and 2 percent high at 64³, 128³, 192³,\n");
      std::print(myfile, "# 256³ and 320³, while marching cubes moves by one percent over that whole range.\n");
      std::print(myfile, "#\n");
      std::print(myfile, "#     level [K]      coarea [Å²]    marching [Å²]     difference\n");
      for (const std::array<double, 3> &row : this->crossCheck)
      {
        const double relative = (row[2] > 0.0) ? (row[1] - row[2]) / row[2] : 0.0;
        std::print(myfile, "# {:13.4f} {:16.4f} {:16.4f} {:14.5f}\n", row[0] * toKelvin, row[1], row[2], relative);
      }
      std::print(myfile, "#\n");
      std::print(myfile, "# What makes the curve usable well before it has converged is that the error above is\n");
      std::print(myfile, "# nearly the same at every level rather than wandering with it. The shape is therefore\n");
      std::print(myfile, "# right where the scale is not, and the scale is fixed by the one level that had to be\n");
      std::print(myfile, "# triangulated anyway. Column 3 below is the curve multiplied by {:.5f}, which brings\n",
                 this->anchorScale);
      std::print(myfile, "# it onto the marching-cubes area at the iso-value, and is the column to read.\n");
    }

    std::print(myfile, "#\n");
    std::print(myfile, "# column 1: level [K]\n");
    std::print(myfile, "# column 2: area of that level set as the coarea pass gives it [Å²]\n");
    std::print(myfile, "# column 3: the same anchored to marching cubes at the iso-value [Å²]\n");
    std::print(myfile, "# column 4: column 3 per unit mass [m²/g]\n");
    std::print(myfile, "# column 5: share of the cell below the level\n");
    std::print(myfile, "#  level [K]       coarea [Å²]     anchored [Å²]       [m²/g]      below\n");

    const double perGram = (framework.mass > 0.0)
                               ? Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass
                               : 0.0;
    for (const AreaCurvePoint &point : this->curve.points)
    {
      const double anchored = point.area * this->anchorScale;
      std::print(myfile, "{:11.4f} {:17.4f} {:17.4f} {:12.4f} {:10.6f}\n", point.level * toKelvin, point.area,
                 anchored, anchored * perGram, point.fractionBelow);
    }
  }

  myfile.close();
}
