module;

module energy_shared_well_surface;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import surface_curvature;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// Everything that has to know about Kelvin, kept apart from the arithmetic that does not.
//
// The measurement next door works in the units the field is held in and so needs nothing from the unit system;
// a temperature in Kelvin is a thing a *report* and a *driver* deal in. Splitting them is not tidiness: the
// conversion factors are mutable statics belonging to the engine, and a translation unit that touches them
// drags that dependency along with it. Keeping them out of the measurement is what lets the measurement be
// linked, and tested, against the structural library alone.

void writeWellSurface(std::ostream &stream, const WellSurface &surface)
{
  const double toKelvin = Units::EnergyToKelvin;

  std::print(stream, "#\n");
  std::print(stream, "# THE WELL SURFACE: WHERE THE MOLECULE ACTUALLY SITS\n");
  std::print(stream, "#\n");
  std::print(stream, "# The iso-surface above is the inner turning point. Zero energy is where the repulsion\n");
  std::print(stream, "# balances the attraction on the way in, so it is the closest a molecule with nothing to\n");
  std::print(stream, "# spare gets before it is thrown back; it is not where it sits. It sits at the bottom of\n");
  std::print(stream, "# the well, a little further out, and that is the surface measured here.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Construction: from every vertex of the zero surface, step out along the wall normal and\n");
  std::print(stream, "# stop where the energy stops falling. The normal is held fixed along the ray, and that\n");
  std::print(stream, "# matters: if it were taken to be the local gradient direction at each step then\n");
  std::print(stream, "# dU/dn = -|grad U|, which vanishes only at the critical points of the landscape, and the\n");
  std::print(stream, "# whole surface would collapse onto the handful of adsorption sites. With the normal fixed\n");
  std::print(stream, "# at the wall, dU/dn = 0 is a condition on a surface, and that surface is the locus of\n");
  std::print(stream, "# wells.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Reading the compression. A normal offset of t multiplies area by (1 + t k1)(1 + t k2),\n");
  std::print(stream, "# so where the wall is convex the well surface is the larger of the two --- on an isolated\n");
  std::print(stream, "# atom of radius R it is (1 + t/R)^2 times the zero surface --- and where it is concave,\n");
  std::print(stream, "# inside a pocket or in the corner where three atoms meet, it is the smaller.\n");
  std::print(stream, "#\n");
  std::print(stream, "# NARROW PORES GO TWO WAYS AND BOTH ARE ACCOUNTED FOR HERE. Once a pore is narrow enough\n");
  std::print(stream, "# that the wells of the walls facing each other have merged into a single minimum, every\n");
  std::print(stream, "# wall around it maps onto that one well. Where the wall curves round on itself --- a\n");
  std::print(stream, "# cylindrical channel, a spherical cage --- the offset contracts as it merges: the surface\n");
  std::print(stream, "# is carried in towards the axis or the centre and its area falls away to a curve or to a\n");
  std::print(stream, "# point. That is the force field saying the region is volume and not extra surface, and no\n");
  std::print(stream, "# area at a fixed level can say it. But where the two walls are flat and parallel there is\n");
  std::print(stream, "# no curvature to do the contracting: both map onto the mid-plane with their area intact,\n");
  std::print(stream, "# and the pair of them then report twice the one surface that is really there. So a well\n");
  std::print(stream, "# is tested for whose it is, by carrying the ray on past the minimum, and the area on\n");
  std::print(stream, "# shared wells is halved. `Well surface (shared halved)` is the number to quote.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Rays walked in strides of {:.5f} [Å], up to {:.2f} [Å]\n", surface.step, surface.longestWalk);
  std::print(stream, "# Triangles: {} mapped, {} of them on a shared well, {} folded over\n",
             surface.numberOfTriangles, surface.numberOfSharedTriangles, surface.numberOfFoldedTriangles);
  std::print(stream, "# {} had no well on their normal; {} were discarded as implausibly large before the walk\n",
             surface.numberOfUnmappedTriangles, surface.numberOfRejectedTriangles);
  std::print(stream, "# Vertices already rising at the wall: {}\n", surface.numberOfVerticesAtWall);
  std::print(stream, "# Timing: {} [s] for the map\n", surface.seconds);
  std::print(stream, "\n");

  std::print(stream, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "Well surface (shared halved):", surface.deduplicatedArea(), surface.gravimetricArea,
             surface.volumetricArea);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Well surface as the rays gave it:", surface.area);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Of that, on shared wells:", surface.sharedArea);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Zero surface mapped from:", surface.zeroArea);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Of that, unmapped:", surface.unmappedZeroArea);
  std::print(stream, "{:34} {:14.6f} [-]\n", "Compression:", surface.compression());
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Folded over:", surface.foldedArea);
  std::print(stream, "{:34} {:14.5f} [Å]\n", "Mean walk out:", surface.meanWalk);
  std::print(stream, "{:34} {:14.4f} [K]\n", "Mean well depth:", surface.meanDepth * toKelvin);
  std::print(stream, "{:34} {:14.4f} [K]\n", "Deepest well:", surface.deepestWell * toKelvin);
  std::print(stream, "\n");

  std::print(stream, "# The Boltzmann-weighted area is the integral of exp(-U_min/kT) over the well surface at\n");
  std::print(stream, "# {:.2f} [K]. On a level set there is no such quantity to be had --- the energy is the same\n",
             surface.thermalEnergy * toKelvin);
  std::print(stream, "# on all of it by construction, so the weight is one constant that divides straight back\n");
  std::print(stream, "# out. Here the depth varies over the surface, and the variation is physical rather than\n");
  std::print(stream, "# discretization: a shallow well on a convex cap facing a wide pore counts for little, a\n");
  std::print(stream, "# deep one in a corner touched on several sides at once counts for a great deal. It is the\n");
  std::print(stream, "# quantity that behaves like a Henry coefficient per unit area, and it is not an area: it\n");
  std::print(stream, "# is an area times a dimensionless weight that has no upper bound.\n");
  std::print(stream, "\n");
  std::print(stream, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g]\n",
             "Boltzmann-weighted area:", surface.deduplicatedWeightedArea(), surface.gravimetricWeightedArea);
  std::print(stream, "{:34} {:14.5f} [-]\n", "Mean weight:", surface.enhancement());
  std::print(stream, "\n");

  const CurvatureAreas &bare = surface.curvature;
  const CurvatureAreas &weighted = surface.weightedCurvature;
  if (bare.total() <= 0.0) return;

  std::print(stream, "# How the well surface divides by shape, and how that division changes once each patch is\n");
  std::print(stream, "# weighted by how much a molecule wants to be on it. The normals are the wall normals the\n");
  std::print(stream, "# rays were fired along, which are the normals of the mapped surface exactly where the\n");
  std::print(stream, "# walk is the same length everywhere and closely otherwise.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The folded area is left out of both rows. An offset carries a curvature k to k/(1 + tk),\n");
  std::print(stream, "# which keeps its sign only while 1 + tk stays positive, and losing that is what folding\n");
  std::print(stream, "# is: a folded patch comes back with both curvatures reversed, so a concave corner the\n");
  std::print(stream, "# offset has overshot would be counted here as a convex cap. Past its focal point the\n");
  std::print(stream, "# offset is not a surface a molecule sits on. So these two rows cover {:.4f} [Å²] of the\n",
             bare.total());
  std::print(stream, "# {:.4f} [Å²] the rays gave, the rest being fold.\n", surface.area);
  std::print(stream, "#\n");
  std::print(stream, "# The two rows are the same surface counted two ways: the bare one says how much surface\n");
  std::print(stream, "# there is of each shape, the weighted one how much of the adsorption happens on it.\n");
  std::print(stream, "#\n");
  std::print(stream, "# READ THE SHIFT BETWEEN THEM WITH CARE, and read it as suggestive rather than measured.\n");
  std::print(stream, "# The normals here are the wall normals, which are the mapped surface's own only where the\n");
  std::print(stream, "# walk is the same length from one vertex to the next; where it is not, the mapped surface\n");
  std::print(stream, "# tilts away from them by an amount set by how fast the walk length varies along it. That\n");
  std::print(stream, "# is worst in the narrowest places, and the narrowest places are where the wells are\n");
  std::print(stream, "# deepest, which is where the weight puts its emphasis. So the weighted row leans hardest\n");
  std::print(stream, "# on the part of the surface whose normals are least trustworthy. Doing better wants the\n");
  std::print(stream, "# mapped mesh's own vertex normals, which means welding the triangle soup into a mesh.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Fractions of the classified area, so the first four add to one; the last is the\n");
  std::print(stream, "# unresolved share of the whole.\n");
  std::print(stream, "#                              area [Å²]     convex     saddle    concave       flat unresolved\n");
  std::print(stream, "{:<26} {:13.5f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n", "Well shape, by area:",
             bare.classified(), bare.convexFraction(), bare.saddleFraction(), bare.concaveFraction(),
             bare.flatFraction(), bare.unresolvedFraction());
  std::print(stream, "{:<26} {:13.5f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n", "Well shape, by weight:",
             weighted.classified(), weighted.convexFraction(), weighted.saddleFraction(),
             weighted.concaveFraction(), weighted.flatFraction(), weighted.unresolvedFraction());
  std::print(stream, "Integral of the mean curvature:     {} [Å]\n", bare.integratedMeanCurvature);
  std::print(stream, "Integral of the Gaussian curvature: {} [-]\n", surface.integratedGaussianCurvature);
  std::print(stream, "# The two are taken over different parts of the surface, and that is deliberate. The\n");
  std::print(stream, "# Gaussian one covers all of it, folds included, because it is 2 pi times the Euler\n");
  std::print(stream, "# characteristic and a count needs the surface closed. The mean one cannot be had that way\n");
  std::print(stream, "# --- a folded patch returns H with its sign reversed and there is no recovering it ---\n");
  std::print(stream, "# so it covers the unfolded part, as the two rows above do.\n");
  std::print(stream, "#\n");
  std::print(stream, "# USE THE GAUSSIAN ONE AS A CHECK ON THE MAP. Set it against the same integral over the\n");
  std::print(stream, "# zero surface in the level-set report. A normal offset is a continuous deformation and\n");
  std::print(stream, "# cannot change the topology, so on a map that is behaving the two agree, and on MFI with\n");
  std::print(stream, "# a third of the area folded they agree to about six percent, which is the accuracy of the\n");
  std::print(stream, "# estimator rather than of the map. Where they do not agree the map is straining: a folded\n");
  std::print(stream, "# patch keeps the sign of K only if both of its principal directions turned over, and one\n");
  std::print(stream, "# of them turning over alone reverses it. On a narrow one-dimensional channel such as ABW,\n");
  std::print(stream, "# with getting on for half the area folded, this comes out nowhere near the zero surface's\n");
  std::print(stream, "# value, and that is the honest signal that the well surface there is barely a surface.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The weighted row's own integrals carry the weight through as well, so they are neither\n");
  std::print(stream, "# Minkowski functionals nor multiples of the Euler characteristic, and are not reported.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The bare area here converges more slowly with the grid than the level-set area does. The\n");
  std::print(stream, "# well of a Lennard-Jones pair lies only 2^(1/6) - 1 of sigma outside its zero crossing, a\n");
  std::print(stream, "# third of an Ångström for a typical probe, which is a few voxels on any usable grid. The\n");
  std::print(stream, "# depth of a minimum is second-order accurate in how well its position was found, so the\n");
  std::print(stream, "# depths and the weighted area settle quickly; the position is first-order, and the area is\n");
  std::print(stream, "# built out of positions. Check it against the spacing before quoting it.\n");
}

MolecularWellSurface::MolecularWellSurface() {}

MolecularWellSurface::~MolecularWellSurface() {}

void MolecularWellSurface::run(const EnergyBackend &backend, const PairInteractions &interactions,
                               const Crystal &framework, const LinearProbe &probe, double level, uint3 gridSize,
                               std::size_t numberOfOrientations, double temperature, double longestWalk,
                               bool useElectrostatics, double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);

  this->grid = field.grid;
  this->potential = field.potential;
  this->isoValue = level;

  if (this->grid.numberOfVoxels() == 0) return;

  // The free-energy landscape rather than the minimum-energy one. The walk asks where a molecule sits at this
  // temperature, and that is a question about the landscape whose orientations have been averaged over; the
  // minimum-energy landscape answers where it would sit if it could be turned the best way for free.
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
  std::vector<double3> corners = backend.isosurfaceTriangles(this->grid.freeEnergy, this->grid.gridSize, level);
  std::chrono::duration<double> extraction = std::chrono::steady_clock::now() - time_begin;

  this->surface = wellSurfaceOfField(framework, this->grid.freeEnergy, this->grid.gridSize, corners, level,
                                     temperature * Units::KelvinToEnergy, longestWalk);

  double3 spacing = this->grid.spacing();
  double density = 1e-3 * framework.mass / (framework.unitCell.volume * Units::Angstrom * Units::Angstrom *
                                            Units::Angstrom * Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.wells." + this->grid.backend + ".txt");
  std::print(myfile, "# The adsorption surface: the locus of energy wells over the framework\n");
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
  std::print(myfile, "# Iso-value the walk started from: {} [internal] {:.4f} [K]\n", this->isoValue,
             this->isoValue * Units::EnergyToKelvin);
  std::print(myfile, "# Landscape: the free-energy one, its orientations averaged over\n");

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
  std::print(myfile, "# Timing ({}): {} [s] for the zero surface\n", this->grid.backend, extraction.count());

  writeWellSurface(myfile, this->surface);
  myfile.close();
}
