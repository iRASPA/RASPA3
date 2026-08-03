module;

export module energy_shared_molecular_surface_area;

import std;

import uint3;
import crystal;
import pair_interactions;
import grid_area_curve;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// The surface a rigid linear molecule sees in a framework, taken as an iso-surface of its energy landscape.
//
// For a single atom there is one energy at a point and so one surface. A molecule has an energy for every way
// it can be turned, and the two reductions of that give two different surfaces, both of which are reported
// because the gap between them is the interesting part.
//
// The surface of `freeEnergy` is the one to quote. It bounds the region the molecule can occupy at temperature
// T once the turning has been averaged over, so a narrow window that the molecule can only thread in one
// orientation counts for little, as it should.
//
// The surface of `minimumEnergy` bounds the region the molecule could occupy if it were always turned the best
// way, free of charge. It is the zero-temperature limit, and the region it bounds always contains the
// free-energy one, since the least over orientations never exceeds the average. Quoting it alone would
// overstate what is open to the molecule. Note that the containment is of the regions and not of their areas,
// which nesting alone does not order, though in practice the larger region is also the rougher one and comes
// out with the larger area.
//
// The area of a level set is a far more delicate thing to converge than a barrier or a spatial average. A
// barrier reads one number off the field and an average smooths the whole of it, but an area measures how
// crumpled a surface is, so any noise in the field is counted rather than cancelled. Both the grid and the
// orientation quadrature therefore have to be taken further here than elsewhere, and the minimum-energy
// surface is the worse of the two: a least over a discrete set of orientations has kinks that never smooth
// out, and its area creeps upward with the number of orientations instead of settling.
export struct MolecularSurfaceArea
{
  // Areas of the two iso-surfaces, per unit cell, in Å².
  double freeEnergyArea{0.0};
  double minimumEnergyArea{0.0};

  // The same areas per unit mass [m²/g] and per unit volume [m²/cm³], for the free-energy surface.
  double gravimetricArea{0.0};
  double volumetricArea{0.0};

  // And for the minimum-energy surface.
  double minimumEnergyGravimetricArea{0.0};
  double minimumEnergyVolumetricArea{0.0};

  std::size_t numberOfFreeEnergyTriangles{0};
  std::size_t numberOfMinimumEnergyTriangles{0};
  std::size_t numberOfRejectedTriangles{0};

  // The area at every level rather than at one, on the free-energy landscape.
  //
  // The area is the one quantity of the energy route that cannot be freed of the level it is read at. A
  // diameter, a volume and a pore size all describe a region, and a region can be pinned to a contour that
  // means something; an area describes a surface, and the surface is the level. Nor can a Boltzmann weight
  // rescue it. Weighting the surface by exp(-A/kT) does nothing whatever, since A is the same on every part
  // of it by construction, and the weight is one constant that divides straight back out; what little
  // variation there is across the triangles is interpolation error, so weighting by it would be weighting by
  // the discretization.
  //
  // What the temperature does reach is the landscape itself, through the orientational average, and there it
  // has a large effect: carbon dioxide in MFI loses a fifth of its free-energy area between 100 and 1000 K,
  // while the minimum-energy area, which has no temperature in it, does not move at all.
  //
  // So rather than one area, or an average of areas under some invented weight, the whole family is reported
  // and the choice of where to read it left to whoever knows what the number is for.
  AreaCurve curve;

  // The curve read at a few levels beside marching cubes run at those same levels. The two are different
  // constructions and do not agree on a coarse grid: the curve comes out high, because a central difference
  // taken across a wall that climbs as the twelfth power overstates the gradient badly when the voxel is wide.
  // The gap closes as the grid is refined, but slowly, and much more slowly than marching cubes settles. On
  // MFI with argon the curve is 32%, 18%, 7%, 4% and 2% high at 64³, 128³, 192³, 256³ and 320³, while
  // marching cubes moves by only 1% over the same range and is within 0.3% of its limit by 128³.
  std::vector<std::array<double, 3>> crossCheck;

  // What saves the curve at a usable grid is that the error is nearly the same at every level: at 128³ it runs
  // between 17% and 18% across the whole range rather than wandering. The shape is therefore trustworthy where
  // the scale is not, and one marching-cubes run at one level supplies the scale. This is the factor that
  // does it, and the anchored column of the report is the curve multiplied by it.
  double anchorScale{1.0};

  double isoValue{0.0};
  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  // The two headline iso-surfaces, and the spot surfaces the curve is checked against, kept apart so that the
  // cost of the curve can be set against the cost of the thing it replaces.
  double seconds{0.0};
  double crossCheckSeconds{0.0};

  MolecularSurfaceArea();
  ~MolecularSurfaceArea();

  // How much of the zero-temperature surface falls away once the cost of orientation is charged for. It is a
  // measure rather than a share, since the two areas are not ordered by anything stronger than habit, and it
  // may come out negative for a molecule whose free-energy surface is the more crumpled of the two.
  double orientationalExcess() const
  {
    if (this->minimumEnergyArea <= 0.0) return 0.0;
    return (this->minimumEnergyArea - this->freeEnergyArea) / this->minimumEnergyArea;
  }

  static MolecularSurfaceArea fromGrid(const EnergyBackend &backend, const MolecularEnergyGrid &grid,
                                       const Crystal &framework, double isoValue);

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe, double isoValue,
           uint3 gridSize, std::size_t numberOfOrientations, double temperature, bool useElectrostatics = true,
           double relativePrecision = 1e-6);
};
