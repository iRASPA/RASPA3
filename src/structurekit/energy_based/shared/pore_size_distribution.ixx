module;

export module energy_shared_pore_size_distribution;

import std;

import uint3;
import crystal;
import pair_interactions;
import grid_pore_size;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// The pore-size distribution a molecule sees, measured on its energy landscape rather than on the geometry.
//
// Gelb and Gubbins ask, at each point of the void, for the diameter of the largest sphere that fits in the
// void and covers the point. That needs the void to have a boundary, which is the one thing an energy field
// does not come with: there is no surface, only a level set that moves as the level is changed. So the level
// has to be named, and it is `isoValue`. The void is where the molecule's energy is below it, everything
// downstream is the ordinary Gelb-Gubbins construction on that region, and the curve is a curve at one
// energy in the way that the geometric one is a curve at one set of radii.
//
// Choosing the level is the whole of what this method asks of its user, and it is not a nuisance parameter
// to be got rid of. At `isoValue` zero the void is where the molecule is not repelled, which is the closest
// thing to the geometric answer: for a hard sphere the energy is infinite inside and zero outside, the level
// set is the contact surface whatever finite level is picked, and the curve becomes the geometric one
// exactly. Above zero the void takes in the windows a molecule crosses by paying a few kT, and the cavities
// begin to join up into channels; below zero it shrinks onto the adsorption sites themselves. Running a few
// levels says more than any single one of them.
//
// Two curves are reported, from the two reductions of the landscape over orientations, and they mean what
// they mean in the surface area next door. The free-energy curve is the one to quote, since it charges the
// molecule for the freedom it gives up in a narrow place. The minimum-energy curve grants the best
// orientation everywhere at no cost, so its void is always the larger of the two and its pores the wider.
//
// A third curve weights each point by how much of its time the molecule spends there, and it is the one that
// answers the question a measurement asks. The two above count room: a wide cavity the molecule dislikes and
// a snug one it never leaves contribute in proportion to their volumes. The weighted one counts molecules,
// so the snug pore counts for as much as its Boltzmann factor says, and the answer is about where the
// adsorbate is rather than about where there is space.
//
// It is also the one curve here with no iso-value in it. Weighting the counting is not by itself enough to
// get rid of the level, which was the first thing tried and does not work: the size at a point is the width
// of the largest sphere that fits in the void, so raising the level widens every sphere and slides the whole
// axis, and re-weighting the points does nothing about an axis that has moved. On MFI with argon, moving the
// level over the 900 K either side of zero moves Di from 3.05 to 3.28 Å and carries the weighted curve along
// with the volumetric one, the two shifting by the same tenth in their means.
//
// So the weighted curve is measured against a contour that does not move: A = 0, where the molecule is as
// well off as it would be in the gas. That is the energetic reading of touching, and it is not an arbitrary
// level in the way the others are: for a single Lennard-Jones pair it is the surface r = sigma exactly, so
// the contour sits on the same contact surface the geometric route builds its void from, and the two curves
// are measuring the same thing. What is left in the level's place is the temperature, which is a property of
// the experiment rather than of the analysis, and which the curve ought to depend on.
//
// The contour has to close somewhere, and a little of the molecule's time is spent on the far side of it,
// where no size is defined and which the curve therefore leaves out. How much is reported.
export struct EnergyPoreSizeDistribution
{
  // The energy below which the molecule is taken to be in the void, in internal units.
  double isoValue{0.0};

  // The temperature the weighting is done at, in Kelvin.
  double temperature{0.0};

  // How many kT a piece must cost the molecule to leave before its width is left out of the accessible
  // curve. A piece that runs through the crystal is accessible whatever this is; the threshold decides only
  // the closed ones, which on a landscape are not all shut.
  double thresholdInKT{30.0};

  // How far a sphere may reach past its radius and still count as covering a grid point.
  double slack{0.0};

  // The curve on the free-energy landscape, which is the one to quote.
  std::vector<PoreSizeCurvePoint> points;
  double largestDiameter{0.0};
  double voidFraction{0.0};
  double accessibleVoidFraction{0.0};

  // The same landscape, on the fixed A = 0 contour, weighted by exp(-A/kT) instead of by volume.
  std::vector<PoreSizeCurvePoint> occupancyPoints;

  // The mean pore a unit of volume sits in, and the mean pore a molecule sits in, both over that contour so
  // that they are answers about the same region. Which of the two is larger is worth looking at. A spherical
  // probe goes for the tighter places, where the walls close in from both sides and the well is deepest, and
  // the second comes out below the first: argon in MFI sits in 2.27 Å against 2.37 Å of room at 300 K, and
  // the gap closes as the temperature rises until at 1000 K the weights are flat and the two agree. A
  // molecule with a shape of its own can go the other way, because a narrow place costs it the freedom to
  // turn, and above about 500 K that cost sends CO2 in MFI into the wider pores instead.
  double volumetricMeanDiameter{0.0};
  double occupancyMeanDiameter{0.0};

  // The share of the molecule's time that is spent inside the reference contour, and so is spoken for by the
  // weighted curve. The rest is spent pressed into the wall, where the curve has no size to offer it. It
  // falls as the temperature rises, and a value far below one is the sign that the curve is being asked
  // about a temperature at which the notion of sitting in a pore is itself getting thin.
  double occupancyInsideReference{0.0};

  // The share of the adsorbed molecules that sit in void reachable from outside the crystal. This is the
  // quantity that ought to be used to correct an uptake for blocked pockets, and it is not the same as the
  // share of the volume that is reachable, since a sealed pocket is usually a deep one.
  double reachableOccupancyFraction{0.0};

  // How far the void at this level runs: zero for cavities sealed from one another, and one, two or three
  // for a channel, a layer and a network. It is a property of the level and not of the framework, and
  // watching it change as the level is raised is the point of being able to raise it.
  int dimensionality{0};

  // The same curve on the minimum-energy landscape.
  std::vector<PoreSizeCurvePoint> minimumEnergyPoints;
  double minimumEnergyLargestDiameter{0.0};
  double minimumEnergyVoidFraction{0.0};

  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;

  double sweepSeconds{0.0};
  double seconds{0.0};

  EnergyPoreSizeDistribution();
  ~EnergyPoreSizeDistribution();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double isoValue, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, double thresholdInKT = 30.0, bool useElectrostatics = true,
           double relativePrecision = 1e-6, std::optional<double> maximumDiameter = std::nullopt,
           std::optional<std::size_t> numberOfBins = std::nullopt);
};
