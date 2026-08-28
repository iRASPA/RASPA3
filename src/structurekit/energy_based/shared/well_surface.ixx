module;

export module energy_shared_well_surface;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import surface_curvature;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

// The surface a molecule actually sits on, rather than the one it is forbidden to enter.
//
// Every area on an energy landscape elsewhere in this library is the area of a level set, and the level has to
// be chosen. Choosing zero is conventional and it is not arbitrary --- it is where the repulsion balances the
// attraction --- but it is the wrong place to look for an adsorbed molecule. Zero energy on the way in is the
// *inner turning point*: the closest a molecule with no kinetic energy gets before it is thrown back. It sits
// further out than that, at the bottom of the well, and the whole of adsorption happens there.
//
// So this route does not read the landscape at a different level. It maps one surface onto another. From every
// vertex of the zero surface it steps out into the void along the wall normal and stops where the energy stops
// falling, and the locus of those stopping points is a second surface: the locus of wells. Its area is what
// the force field has to say about how much room there is to adsorb on, and unlike the area of a level set it
// is not a number with a free parameter in it.
//
// WHY THE NORMAL IS HELD FIXED. "Walk down -grad U" and "walk until dU/dn = 0" are the same instruction only
// if n is fixed at the start. If n is taken to be the local gradient direction at each step, then
// dU/dn = grad U . n = -|grad U|, which vanishes only where the gradient itself does, and gradient descent
// runs to the isolated critical points of the landscape: the whole surface would collapse onto the handful of
// adsorption sites and there would be no area left to measure. Fixing n at the wall makes dU/dn = 0 a
// condition on a two-dimensional locus, which is the surface wanted. The walk is therefore a straight ray
// along the outward normal of the zero surface, and the stopping point is the first minimum of the energy
// along that ray.
//
// WHAT THE MAP DOES TO AREA, which is the physical content of the result.
//
// Where the wall is convex --- the outer face of an atom bulging into a wide pore --- the well lies a short way
// out and the offset *expands* the surface, by (1 + t k1)(1 + t k2) for a walk of t at curvatures k1, k2. On a
// sphere of radius R that is (1 + t/R)^2, so the well surface of an isolated atom is larger than its zero
// surface. Where the wall is concave --- the inside of a pocket, a corner where three atoms meet --- the same
// offset *contracts* it, and where it is a saddle the two factors fight.
//
// Where the pore is narrow enough that the wells of opposite walls have merged into one, there is a single
// minimum in the middle of the channel and the rays from both walls end on it. The map is then many-to-one: a
// patch of the zero surface maps onto a curve or onto a point, and its contribution to the area goes to
// nothing. That is not a failure of the construction, it is the force field saying that such a region is
// volume and not extra surface, which is exactly the thing an area at a fixed level cannot say. It shows up
// here as `foldedArea`, the mapped triangles that turned inside out on the way, and as a `compression` well
// below one.
//
// THE BOLTZMANN WEIGHT. A level set has one energy on all of it by construction, so weighting it by
// exp(-U/kT) multiplies it by a constant and says nothing --- the report of the level-set route says so at
// length. The well surface is different: the depth of the well varies over it, and that variation is real
// rather than discretization. A shallow well on a convex cap sticking into a big void gets a weight near one,
// a deep well in a corner where the molecule is touched on several sides at once gets a large one. So
// `weightedArea` is the integral of exp(-U_min/kT) over the well surface, and it is the quantity that behaves
// like a Henry coefficient per unit area: it counts the surface in proportion to how much a molecule wants to
// be on it. The curvature split is reported both ways for the same reason, and the shift between the two
// columns is the whole story of where adsorption happens.
//
// WHAT TO DISTRUST. The zero surface comes from marching cubes on a grid and so does the walk, which reads the
// field between its nodes. The well of a Lennard-Jones pair lies only about a tenth of sigma outside its zero
// crossing --- 2^(1/6) sigma against sigma, so a third of an Ångström for a typical probe --- which is a
// handful of voxels on any grid you would use. The *depth* of a minimum is second-order accurate in how well
// its position was found, so `meanDepth` and `weightedArea` settle quickly; the *position* is first-order, and
// the area of the mapped surface is built out of positions. So the bare area here converges more slowly than
// the level-set area does and should be checked against the grid before it is quoted.

// Samples a periodic scalar field held on a grid, at points that are not its nodes.
//
// Two things are interpolated and they are interpolated differently, on purpose. The value is trilinear in the
// eight nodes around the point, which is the cheapest thing that is continuous. The gradient is *not* the
// derivative of that trilinear form, which would be piecewise constant along each axis and would jump across
// every voxel face; it is the trilinear interpolation of the central differences held at the eight nodes. That
// is continuous, it is the same estimate of the gradient that marching cubes puts at its vertices, and it is
// second-order where the derivative of the trilinear form is first-order.
//
// The consequence worth knowing is that `gradient` is not exactly the derivative of `value`. For finding where
// a directional derivative vanishes that is the better bargain: the criterion is then a continuous function of
// position and can be bracketed and bisected, where the derivative of the trilinear form would let the walk
// stop on a voxel face that is an artefact of the interpolation rather than a feature of the field.
export struct PeriodicFieldSampler
{
  uint3 gridSize{0, 0, 0};
  std::span<const float> field;

  bool usable() const;

  // The stored value at a node, the index wrapping into the cell however far outside it is.
  double nodeValue(std::int64_t i, std::int64_t j, std::int64_t k) const;

  // The central difference at a node, held per grid step, which is the convention
  // `cartesianGradientOfGridGradient` expects.
  double3 nodeGradient(std::int64_t i, std::int64_t j, std::int64_t k) const;

  double value(double3 fractional) const;
  double3 gradient(double3 fractional) const;
};

// Where one vertex of the zero surface ended up.
export struct WellWalk
{
  double distance{0.0};  // Å along the outward normal
  double depth{0.0};     // the energy there, in internal units

  // False when the ray ran the whole way to `longestWalk` without the energy turning back up. That happens
  // where the pore is wide enough that the far wall is out of reach and the field has flattened to nothing
  // before any second wall pulls it back, and it means there is no well on that normal to map to.
  bool found{false};

  // True where the energy was already rising at the wall, so the well is the wall. Rare, and worth counting
  // rather than hiding: it is the signature of a place too tight for the zero surface to have a void side.
  bool atWall{false};

  // True where this well is not this wall's own but one it shares with a wall facing it, the two having merged
  // into a single minimum in the middle of the gap. Such a well is reached by rays from both walls, so the area
  // around it is measured twice and has to be halved to be counted once.
  //
  // Told apart from a wall's own well by carrying on past the minimum. If the energy climbs the whole way back
  // up through the iso-value then there is solid immediately beyond and the minimum sits between two walls. If
  // instead the climb crests and the energy turns back down, there is a separate well over there belonging to
  // the far wall, and this one belongs to this wall alone.
  bool shared{false};
};

// Follows one ray. `start` and `normal` are Cartesian, the normal a unit vector pointing away from the solid.
// `isoValue` is the level the surface was taken at, needed only to decide whether the well is shared.
//
// Exported because it is the whole of the physics and is worth being able to test on a field whose wells are
// known analytically, without a framework or a marching-cubes mesh in the way.
export WellWalk walkToWell(const PeriodicFieldSampler &sampler, const UnitCell &unitCell, double3 start,
                           double3 normal, double step, double longestWalk, double isoValue);

// The well surface, and how it differs from the zero surface it was mapped from.
export struct WellSurface
{
  double zeroArea{0.0};          // Å², the surface walked from, per unit cell
  double unmappedZeroArea{0.0};  // Å² of that whose rays found no well
  double area{0.0};              // Å², the well surface as the rays found it, shared wells counted twice

  // The same per unit mass and per unit volume, with the shared wells counted once: these are of
  // `deduplicatedArea` and not of `area`, that being the number to quote.
  double gravimetricArea{0.0};  // m²/g
  double volumetricArea{0.0};   // m²/cm³

  // Å², the integral of exp(-U_min/kT) over the well surface, and the same per unit mass with shared wells
  // counted once. Larger than the bare area whenever the typical well is deeper than kT, which for any real
  // framework at room temperature it is.
  double weightedArea{0.0};
  double gravimetricWeightedArea{0.0};

  // Å² of mapped area whose triangle turned inside out under the map, which is where the rays have crossed
  // each other. Counted in `area` as well, there being no honest way to subtract it, but reported so that a
  // result which is mostly fold can be recognised as one.
  double foldedArea{0.0};

  // Å² of mapped area sitting on wells that are shared with a wall facing it, and the same weighted. Two walls
  // whose wells have merged both map onto the one well between them, so this area has been measured twice.
  //
  // There are two ways a narrow pore can go and both happen. Where the wall curves round on itself --- a
  // cylindrical channel, a spherical cage --- the map contracts as well as merging: the offset carries the
  // whole wall towards the axis or the centre and the area falls away to a curve or a point, which is the
  // force field saying the region is volume and not surface. Where the two walls are flat and parallel there is
  // no curvature to contract, both walls map onto the mid-plane with their area intact, and the honest total is
  // half of what the two of them give. `foldedArea` catches the first, this catches the second.
  double sharedArea{0.0};
  double sharedWeightedArea{0.0};

  std::size_t numberOfTriangles{0};
  std::size_t numberOfFoldedTriangles{0};
  std::size_t numberOfSharedTriangles{0};
  std::size_t numberOfUnmappedTriangles{0};
  std::size_t numberOfRejectedTriangles{0};
  std::size_t numberOfVerticesAtWall{0};

  double meanWalk{0.0};     // Å, averaged over the well surface by area
  double meanDepth{0.0};    // internal units, averaged over the well surface by area
  double deepestWell{0.0};  // internal units
  double longestWalk{0.0};  // Å, the cap the rays were given
  double step{0.0};         // Å, the stride they were walked in
  double isoValue{0.0};     // the level the zero surface was taken at

  // kT of the Boltzmann weight, in the same internal units the field is held in rather than in Kelvin. The
  // field is in internal units, so this is the form in which the weight is dimensionally honest; turning a
  // temperature into it is the caller's business and is where the choice of unit system belongs.
  double thermalEnergy{0.0};

  // The shape of the well surface, by area and by Boltzmann weight, over the unfolded part of it. A folded
  // patch has both of its curvatures reversed, the offset having carried it past a focal point, so it is left
  // out of both rather than counted as the opposite of what it is.
  //
  // READ THE SHIFT BETWEEN THE TWO ROWS WITH CARE. The normals used are the wall normals the rays were fired
  // along. For a normal offset of constant length those are exactly the normals of the mapped surface; where
  // the walk varies in length from vertex to vertex they are not, the mapped surface tilting away from them by
  // an amount set by how fast the walk length changes along the surface. That is worst in the narrowest places,
  // where the walk length changes fastest --- and the narrowest places are also where the wells are deepest and
  // so where the Boltzmann weight puts most of its emphasis. So the weighted row leans hardest on the part of
  // the surface whose normals are least trustworthy, and a shift between the two rows is suggestive rather than
  // measured. Getting it properly would want the mapped mesh's own vertex normals, which means welding the
  // triangle soup into a mesh with connectivity first.
  //
  // Only the area columns of the weighted one mean anything. Its curvature integrals carry the weight through
  // as well, so they are neither Minkowski functionals nor multiples of the Euler characteristic, and nothing
  // reads them.
  CurvatureAreas curvature;
  CurvatureAreas weightedCurvature;
  CurvatureBands curvatureBands;

  // The integral of the Gaussian curvature over the whole mapped surface, folds and all, unlike everything in
  // the two above.
  //
  // It is kept apart because it is the one quantity here that needs the surface closed, and dropping the folded
  // patches would open the surface up and lose that. It is 2 pi times the Euler characteristic, and a normal
  // offset is a continuous deformation and so cannot change the topology, which makes it a check on the map:
  // compare it with the same integral over the zero surface. On MFI, with a third of the area folded, the two
  // land within about six percent of one another, which is the accuracy of the estimator rather than of the map.
  //
  // Where they disagree badly the map is straining rather than the arithmetic being wrong. A folded patch keeps
  // the sign of K only when both of its principal directions turned over; one of them turning over alone
  // reverses it. On a narrow one-dimensional channel such as ABW, where getting on for half the area folds,
  // this comes out nowhere near the zero surface's value, and that is worth reading as the signal it is.
  //
  // The mean curvature cannot be had the same way. A folded patch comes back with H reversed and there is no
  // recovering it without knowing which of its two principal directions did the folding, so that integral,
  // inside `curvature` above, covers the unfolded part only.
  double integratedGaussianCurvature{0.0};

  double seconds{0.0};

  // The well surface with each shared well counted once instead of twice. This is the number to quote.
  double deduplicatedArea() const { return this->area - 0.5 * this->sharedArea; }
  double deduplicatedWeightedArea() const { return this->weightedArea - 0.5 * this->sharedWeightedArea; }

  // How much of the zero surface survived the map, the shared wells counted once. Above one where the wall is
  // mostly convex and the offset expanded it, below one where the pores are narrow enough for opposite wells to
  // have merged and the wall curves round on itself.
  double compression() const
  {
    const double mapped = this->zeroArea - this->unmappedZeroArea;
    return (mapped > 0.0) ? this->deduplicatedArea() / mapped : 0.0;
  }

  // The area-averaged Boltzmann weight: how much more a molecule wants to be on this surface than on one whose
  // wells are all of depth kT.
  double enhancement() const { return (this->area > 0.0) ? this->weightedArea / this->area : 0.0; }
};

// How far a ray may be walked before the search is given up, in Å. The well of a pair potential is a fraction
// of an Ångström outside its zero crossing, and a ray crossing a channel to a well in the middle of it has a
// few to travel; past this there is no wall within reach and nothing to find.
export inline constexpr double defaultLongestWalk = 6.0;

// Maps the zero surface onto the locus of wells and measures it.
//
// `corners` is the zero surface as marching cubes hands it back, three fractional corners to a triangle, and
// `field` is the landscape it was extracted from. The two have to be the same field at the same iso-value or
// the rays will start off the surface. `thermalEnergy` is kT of the Boltzmann weight in the same units the
// field is held in, not in Kelvin; pass zero to leave the weighting out, in which case every weight is one.
export WellSurface wellSurfaceOfField(const Crystal &framework, std::span<const float> field, uint3 gridSize,
                                      std::span<const double3> corners, double isoValue, double thermalEnergy,
                                      double longestWalk = defaultLongestWalk);

// Writes the well surface with the commentary that says how to read it.
export void writeWellSurface(std::ostream &stream, const WellSurface &surface);

// The driver: builds the landscape, extracts the zero surface, maps it, and writes the report. Takes the same
// arguments as the level-set surface area next door so the two can be asked for in one run off one landscape.
export struct MolecularWellSurface
{
  WellSurface surface;
  MolecularEnergyGrid grid;
  ElectrostaticPotentialGrid potential;
  double isoValue{0.0};

  MolecularWellSurface();
  ~MolecularWellSurface();

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           const LinearProbe &probe, double level, uint3 gridSize, std::size_t numberOfOrientations,
           double temperature, double longestWalk = defaultLongestWalk, bool useElectrostatics = true,
           double relativePrecision = 1.0e-6);
};
