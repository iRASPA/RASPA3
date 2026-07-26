module;

export module skapolloniusdiagram;

import std;

import int3;
import double3;
import double3x3;

// The Apollonius diagram (additively weighted Voronoi diagram) of a periodic set of spheres.
//
// Space is divided by the additively weighted distance, the clearance
//
//     d_i(x) = |x - x_i| - r_i,
//
// so the cell of site i is the set of points that a probe can reach more freely from i than from
// any other site. This differs fundamentally from the radical (power) diagram in module skvoronoi,
// which uses |x - x_i|² - r_i². The power distance is what makes the radical diagram cheap: its
// bisectors are planes, so cells are convex polyhedra. The clearance bisector between two spheres
// of unequal radius is instead one sheet of a hyperboloid of revolution about the line joining
// them, degenerating to a plane only when the radii are equal. Apollonius cells are therefore
// bounded by curved patches, are star shaped about their site but not convex, and may be empty
// when one sphere is engulfed by another.
//
// The clearance is the quantity that actually governs whether a probe fits, so this diagram, and
// not the radical one, is what pore geometry should be read from. The radical diagram agrees with
// it exactly when all radii are equal and drifts from it as the radii spread.
//
// Representation. Vertices are points equidistant in clearance from four sites, which is to say
// centres of spheres tangent to and outside all four; this is the classical Apollonius problem and
// has a closed-form solution. Edges and faces are curved and are stored implicitly, by the sites
// whose sheets meet along them, rather than as sampled geometry: an edge is an arc of the trisector
// curve of three sites, and a face is a patch of the bisector sheet of two sites. That keeps the
// representation exact. Sampled positions along an edge can be recovered from its defining sites.
//
// Degenerate input. Four sites is the general case but not the only one: a sphere may touch five or
// more at once, and then all of them meet at that point. It is one vertex of the diagram, of degree
// higher than four, and it is held as one, carrying the whole cotangent set; see SKApolloniusVertex.
// Since the tangency of each quadruple drawn from such a set is satisfied, the vertex is constructed
// once per quadruple, and the copies are gathered afterwards. This costs nothing on a generic input,
// where no two vertices coincide and every set has four members, and it means the complex closes on a
// degenerate one instead of dissolving into a cluster of coincident vertices that pair with nothing.
//
// Robustness. Every vertex reported is certified: its tangent sphere is checked against an exhaustive
// neighbourhood search and contains no site, so no reported vertex is spurious and no reported radius
// overstates the room available. Vertex geometry is therefore sound, and quantities read from vertices
// alone, such as the largest sphere that fits, are trustworthy.
//
// Two regions can be asked for, and the difference matters; see SKApolloniusRegion. A tangency
// solution of negative radius is a point equidistant in clearance from four sites that lies inside the
// spheres. It is a genuine vertex of the diagram as a partition of space, and it is no place a probe
// can be. Over `FreeSpace` such solutions are dropped, which clips the complex against the boundary of
// the union of the spheres: an arc running into that boundary stops there rather than at a vertex, so
// Euler's formula and full valence everywhere do not hold and their failure is not a defect.
// `verification` keeps the two apart, counting clipped arcs as `truncatedTriples` and genuine failures
// to pair as `unpairedTriples`. Over `EntireSpace` nothing is clipped and the complex closes exactly.
//
// Vertices are located completely in either case. Candidate quadruples are drawn from radical-diagram
// adjacency, which is a heuristic, and are then backed by an exhaustive subdivision sweep in the manner
// of Wang et al., which refines the cell until a box admits only four nearest sites and so cannot pass
// over a vertex.
//
// Ordering the vertices along a trisector, which is what decides the edges, is done from the shape of the
// trisector itself, and a trisector is always a plane conic. Subtracting the tangency conditions pairwise
// gives one equation per pair, linear in the position and in the clearance alike, so unless all three radii
// are alike the clearance can be eliminated between the two, leaving the plane the curve lies in. Within
// that plane the clearance is an affine function of position, so what is left of the tangency is a quadratic,
// and laying the first axis of the plane along the gradient of the clearance leaves it already in principal
// position, an ellipse where that gradient is smaller than one and the branch of a hyperbola where it is
// larger. Its own parameter then orders the vertices, and the direction an edge leaves a vertex in follows
// from the tangent there, both in closed form.
//
// Doing it this way is what keeps degenerate arrangements from being separate cases with thresholds between
// them, where an answer would turn on which side of a threshold an input fell:
//
//   - Collinear centres bisect in coaxial surfaces of revolution and meet in a circle standing at one
//     clearance the whole way round, which nothing ordered by clearance can order. Here it is just a
//     vanishing gradient, so the ellipse is a circle, ordered round by the same formula as any other.
//   - Nearly collinear centres are a small gradient and a nearly round ellipse, ordered by that same formula,
//     so how nearly collinear they are does not enter into it.
//   - Radii nearly alike are a large gradient and a nearly straight branch, again the same formula, tending
//     to the straight line it should.
//   - Radii exactly alike are the one arrangement standing apart, and only because there is then no
//     clearance to eliminate: both equations are planes in their own right and the trisector is the line they
//     cut in, which is the ordinary Voronoi edge. This is no exotic input, being the case in which the whole
//     diagram reduces to the ordinary Voronoi diagram, and it is an exact condition on the radii rather than
//     a judgement about how nearly equal they are.
//
// Collinear centres also defeat the tangency solve, because two of the solve's rows then have parallel
// coefficients of the centre and the centre cannot be written as a function of the radius. The solve takes
// the line of solutions from the null space of its rows instead, which exists either way and needs no such
// function, so that arrangement is not special there either.
//
// The cell may be triclinic. Nothing here assumes the lattice vectors are orthogonal: the neighbour search
// sizes its reach from the perpendicular widths of the cell rather than from lattice vector lengths, and
// the sweep measures a box by the sphere that encloses it, which is what a sheared box needs. A test
// rebuilds the same crystal in a sheared basis of the same lattice and requires the diagram to come out
// unchanged. Skewed cells do cost more, since the sweep subdivides boxes that are parallelepipeds in
// fractional space and a sheared one reaches much further than its volume suggests, so more of them carry
// enough sites to need splitting.
//
// Trisectors carrying no vertex. When no fourth site ever intrudes on a trisector, the curve is closed
// and has no vertex on it anywhere, so nothing that assembles a diagram from its vertices can reach it.
// Wang et al. discuss these in their Sections 4.2 and 5.5 and close the gap by replicating the offending
// site into four perturbed copies, which reintroduces vertices artificially. They are instead found here
// directly: the subdivision sweep already visits the neighbourhood of every trisector, so it records the
// triples it passes, and a triple left without vertices is then tested by sweeping the clearance along
// its trisector. That decides whether the curve is bounded, which is exactly the condition for it to be
// closed, and where it is bounded the same sweep traverses it and certifies that no site intrudes. Such
// a curve is reported as an edge with `isLoop` set and counted in `verification.vertexlessLoops`.
//
// A ring of this kind is not a topological defect to be removed: it is a ring-shaped channel, and its
// narrowest point is a bottleneck like any other edge's. What a ring does need is to be put on the right
// face. A ring is usually not a region of its own but a hole punched in a neighbouring one, where a third
// site bites into the middle of a bisector patch without reaching its rim, and recording it as a separate
// face instead both invents a face and leaves Euler's formula one over. Which of the two it is follows
// from the band of clearances the ring spans, since its own curve is the only place on the patch where the
// third site is equally near: see `ringBoundsOwnFace` in the implementation.
//
// There is no such thing as a bisector patch carrying no edge at all. A patch lies on one branch of a
// hyperboloid of revolution and so is topologically a plane, other sites dominate it far from the pair, so
// every region of it is bounded, and the boundary of a bounded region consists of points equally near some
// third site, which is an edge. A region that looks edgeless is one bounded solely by a ring, which is the
// case above. This is not left as an argument: a test sweeps every bisector patch in clearance and angle,
// counts the regions on the surface itself, and requires the diagram to hold exactly one face for each.

// One edge leaving a vertex: the triple of the vertex's sites whose trisector it runs along, named by
// position in the vertex's site list, and the direction it sets off in.
//
// The trisector runs through the vertex both ways but only one half need be in the diagram, since on
// the other half a site outside the triple may be the nearer. `direction` is a half that is, following
// Wang et al., "Robust Computation of 3D Apollonius Diagrams" (Computer Graphics Forum 39(7), 2020),
// whose vertex record likewise carries a tangent per outgoing edge.
export struct SKApolloniusBranch
{
  std::array<std::size_t, 3> sites;  // positions within the vertex's site list
  double3 direction;                 // unit outgoing direction along the trisector
};

// A vertex of the diagram: the centre of a sphere tangent to, and outside, four or more sites. Sites
// are named by index together with the periodic image that participates in the tangency.
//
// Four sites is the general case, and then four edges leave the vertex, one per triple. More sites is a
// degenerate configuration, where a sphere happens to touch five or more at once. That is a single
// vertex of the diagram belonging to all of them, not several vertices at one point, and it is held as
// one here: `siteIndices` carries the whole cotangent set. It is then no longer true that every triple
// of the sites carries an edge, so the edges are listed explicitly in `branches` rather than being
// implied by the sites. Where the tangency points are in general position on the tangent sphere their
// convex hull is simplicial and a set of n sites yields 2n - 4 branches, which is 4 when n is 4.
export struct SKApolloniusVertex
{
  double3 position;  // Cartesian, wrapped into the home unit cell [Å]
  double radius;     // radius of the tangent sphere, equal to the clearance here [Å]
  std::vector<std::size_t> siteIndices;
  std::vector<int3> siteImages;
  std::vector<SKApolloniusBranch> branches;  // one per edge leaving the vertex
};

// An edge: the arc of the trisector curve of three sites running between two vertices. Along it
// the clearance is stationary with respect to those three sites, so its narrowest point is the
// bottleneck a probe must pass to travel between the two vertices.
//
// An edge need not have endpoints. When no fourth site ever intrudes on a trisector, the curve is
// closed and carries no vertex, and the whole of it is one edge: a ring-shaped channel. Then `isLoop`
// is true and `from`, `to` and `toImage` are meaningless, while the sites, bottleneck and length
// describe the ring as for any other edge.
//
// Where the bottleneck sits is kept along with how wide it is. It is the narrowest cross-section of
// the passage, so the plane through `bottleneckPosition` perpendicular to `bottleneckDirection`, the
// tangent of the arc there, cuts the window a probe has to get through. Both come from the samples
// that measure the bottleneck in the first place and cost nothing beyond them, and both are in the
// frame in which `from` sits in the home cell.
export struct SKApolloniusEdge
{
  std::size_t from;
  std::size_t to;
  int3 toImage;  // lattice shift applied to `to` to place it at the far end of this arc
  std::array<std::size_t, 3> siteIndices;
  std::array<int3, 3> siteImages;  // images as seen from `from`
  double bottleneckRadius;         // smallest tangent-sphere radius along the arc [Å]
  double3 bottleneckPosition;      // where along the arc that narrowest point sits [Å]
  double3 bottleneckDirection;     // unit tangent of the arc there, the direction the passage runs in
  double length;                   // arc length [Å]
  bool isLoop;                     // closed trisector carrying no vertex; `from` and `to` unused
};

// A face: a connected region of the clearance bisector sheet between two sites, one face per site so
// that each cell owns its own, as with SKVoronoiFace. `edgeIndices` is the unordered set of bounding
// edges. `isClosed` records whether the boundary satisfies the cycle condition, that every vertex on
// the region is met by exactly two of its edges; the edges are not put into cyclic order, which for a
// curved patch would require sampling its geometry.
//
// One bisector sheet may carry more than one such region, and they are separate faces: unlike a Voronoi
// or power diagram, where the sheet between two sites is a single convex polygon, an Apollonius sheet
// can be cut by other sheets into several disconnected pieces (Wang et al., Section 5.4). Counting one
// face per sheet instead of one per region is what makes Euler's formula miss.
//
// A region need not be a disc either. Its boundary can have several components: an outer rim of arcs
// together with holes bitten out of its middle, each hole a closed trisector carrying no vertex. All of
// them appear in `edgeIndices`, which is why that is a set and not a cycle.
export struct SKApolloniusFace
{
  std::size_t site1;
  std::size_t site2;
  int3 site2Image;
  std::vector<std::size_t> edgeIndices;
  bool isClosed;
};

// A cell: everything incident on one site. An engulfed sphere has no cell at all, which is a
// genuine feature of the additively weighted diagram rather than a failure.
export struct SKApolloniusCell
{
  std::size_t siteIndex;
  std::vector<std::size_t> faceIndices;
  std::vector<std::size_t> vertexIndices;
  bool isEmpty;
};

// Which part of space the diagram is built over.
//
// A point equidistant in clearance from four sites has a signed clearance, and where that value is
// negative the point lies inside the spheres. Both choices below are the same construction; they
// differ only in whether such points are admitted.
export enum class SKApolloniusRegion
{
  // Only non-negative clearance: the diagram of the space a probe can actually occupy. The complex is
  // then clipped against the boundary of the union of the spheres, so arcs running inside the spheres
  // stop at that boundary and the counts of a closed complex do not apply. This is what pore geometry
  // wants, and it is cheaper, since the interiors of the spheres are never explored.
  FreeSpace,

  // Every clearance, negative included: the diagram as a partition of all of space, in the sense of
  // Wang et al. Equation 15, where the common weighted distance r may be of either sign. Nothing is
  // clipped, so every vertex is four-valent, every triple pairs, and Euler's formula holds.
  EntireSpace
};

export struct SKApolloniusDiagram
{
  // Whether the complex came out combinatorially consistent. In the diagram over `EntireSpace` every
  // vertex carries an edge along each of its branches and every triple supporting a vertex supports
  // exactly one more, since an arc has two ends. Those identities are weakened over `FreeSpace`, which
  // clips arcs that run inside the spheres; `truncatedTriples` counts the clipped ones, and each costs
  // its surviving vertex one unit of valence. A shortfall beyond that is a real failure of construction,
  // and these counters are the honest measure of it.
  //
  // They hold on a degenerate input as well. A vertex where more than four sites are cotangent has more
  // than four branches, so valence is compared against a vertex's own branch count rather than against
  // four, and the identities are then the same ones. `degenerateVertices` reports that the input is
  // degenerate; it does not mean anything is wrong.
  struct Verification
  {
    std::size_t vertexCount{0};
    std::size_t verticesOfFullValence{0};  // vertices carrying an edge along every branch
    std::size_t unpairedTriples{0};        // triples that could not be paired along their trisector
    std::size_t truncatedTriples{0};       // arcs clipped by the boundary of the free region
    std::size_t overpairedTriples{0};      // triples seen more than twice, i.e. degenerate
    std::size_t coincidentVertices{0};     // distinct vertices left sharing a position, which is a defect
    std::size_t unclosedFaces{0};
    std::size_t vertexlessLoops{0};       // closed trisectors carrying no vertex, recovered separately
    std::size_t ringsOfUncertainFace{0};  // rings that are a hole in one of several regions of a sheet
    // Places where the sweep refined as far as it goes, a thousandth of an Ångström, and still had more
    // sites able to be nearest than it can take in quadruples cheaply, which needs seventeen of them
    // cotangent to within that. Those boxes are solved over the sites they have, so this is not a defect
    // either, but it is the one case the sweep resolves neither by refining nor by finishing early, and it
    // is reported so that it is never silent. An ordinary degeneracy does not reach it: the sweep finishes a
    // box as soon as its sites are few or as soon as splitting stops reducing them, and the second of those
    // is what a degeneracy does, which is what makes it no harder than the general case.
    std::size_t degenerateSweepBoxes{0};
    std::size_t degenerateVertices{0};    // vertices where more than four sites are cotangent
    std::size_t ambiguousBranches{0};     // branches whose direction could not be decided

    // Arcs along which a sample point could not be placed, so their narrowest point is not proven.
    // A bottleneck read from the samples that did land can only come out too wide, never too narrow,
    // and a bottleneck too wide is a passage reported open that is shut.
    std::size_t unsampledArcs{0};

    // True when the complex closes as far as the free region allows: every triple is paired except
    // those clipped by the boundary, every vertex carries an edge along each of its branches except
    // where a clipped arc leaves it one short, and no face is left open except around a clipped arc. An
    // arc lies on the trisector of three sites and so borders three bisector patches, one for each pair
    // among them, which bounds how many faces one clipped arc can leave unclosed.
    //
    // Two vertices at one point would be a cotangent set that failed to be gathered into the single
    // vertex it is, and a branch of undecided direction is an edge that may have been sent the wrong way,
    // so both are disqualifying.
    bool isComplete() const
    {
      return vertexCount > 0 && coincidentVertices == 0 && ambiguousBranches == 0 && unpairedTriples == 0 &&
             overpairedTriples == 0 && unsampledArcs == 0 &&
             verticesOfFullValence + truncatedTriples >= vertexCount && unclosedFaces <= 3 * truncatedTriples;
    }
  };

  std::vector<SKApolloniusVertex> vertices;
  std::vector<SKApolloniusEdge> edges;
  std::vector<SKApolloniusFace> faces;
  std::vector<SKApolloniusCell> cells;
  Verification verification;

  // Radius of the largest sphere that fits anywhere, and where it sits. This is the diagram's
  // deepest vertex, and it is exact rather than approximated from a discretisation.
  double largestEmptySphereRadius() const;
  double3 largestEmptySpherePosition() const;

  // `neighbourRings` controls how far candidate site quadruples are drawn from radical-diagram
  // adjacency: 1 uses direct neighbours, 2 also uses neighbours of neighbours. Raising it widens
  // the search, at cubic cost in the neighbour count, and is the lever to pull when
  // `verification.isComplete()` is false.
  static SKApolloniusDiagram create(const double3x3& unitCell, const std::vector<double3>& fractionalPositions,
                                    const std::vector<double>& radii, std::size_t neighbourRings = 1,
                                    SKApolloniusRegion region = SKApolloniusRegion::FreeSpace);
};

// The spheres tangent to all four given spheres, the classical Apollonius problem in three dimensions.
// Four tangency conditions fix the centre and radius up to a quadratic, so there are at most two solutions;
// fewer are returned when the quadratic has no real root, when `allowNegativeRadius` is false and no root is
// non-negative, or when the four spheres admit no isolated tangent sphere in the first place, as four
// coplanar spheres of equal radii do not: they admit either none or a whole family, never one or two.
//
// Coplanar centres are otherwise no obstacle, nor are three collinear ones, which arise together and which
// no arrangement of the four spheres avoids once three of their centres line up.
//
// A negative radius is a real solution of the same equations: the point is equidistant in clearance
// from the four sites with that clearance negative, so it lies inside them. Such a solution is a
// vertex of the diagram over all of space but not of the diagram of the free space.
export struct SKApolloniusTangentSphere
{
  double3 centre;
  double radius;
};

export std::vector<SKApolloniusTangentSphere> skApolloniusTangentSpheres(const std::array<double3, 4>& centres,
                                                                        const std::array<double, 4>& radii,
                                                                        bool allowNegativeRadius = false);

// The point on the trisector curve of three sites at parameter `t` in [0,1] between two vertex
// positions, obtained by projecting the straight chord onto the curve. Used to sample edges for
// bottleneck and length calculations, and available to callers that need edge geometry.
//
// An arc between two vertices of the free space can still dip inside the sites along the way: the
// window between two cages is such an arc, open at both ends and closed in the middle. Measuring
// its bottleneck means following it there, so `allowNegativeRadius` admits the stretch of the curve
// that lies inside the sites. Without it the dip returns nothing and a passage no probe can pass is
// left looking as wide as its ends.
export std::optional<SKApolloniusTangentSphere> skApolloniusTrisectorPoint(const std::array<double3, 3>& centres,
                                                                          const std::array<double, 3>& radii,
                                                                          const double3& from, const double3& to,
                                                                          double t, bool allowNegativeRadius = false);
