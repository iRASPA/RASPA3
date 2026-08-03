module;

export module exact_void_split;

import std;

import double3;
import pore_accessibility;
import exact_surface_patches;
import exact_boundary_components;

// The void of a structure divided between the channels and the pockets, measured rather than sampled.
//
// How much void there is comes from the volume of the union of the inflated atoms, which is a finite sum
// of closed-form pieces; see exact_union_volume. How it divides was for a long time the part that had to
// be sampled, on the grounds that the natural piece of a volume is the void part of one cell, a region
// that can be disconnected, so that no single label applies to it. That argument is about cells, and the
// division does not have to go through them.
//
// A pore that does not percolate is bounded, and its boundary lies wholly on the spheres of the inflated
// atoms: it is exactly the arcs the surface-area sweep labels as facing that pore. So the divergence
// theorem gives its volume from those arcs alone, with no cell and no connected-component search,
//
//   3 |pore|  =  - sum_a [ R_i |a| + (p_i + L_a - c) . m_a ] ,
//
// summed over its arcs, with L_a carrying each arc into the frame the pore's nodes are assembled in and c
// any fixed point. The channels then take what is left of the void, which needs no argument of its own:
// it is a subtraction. Both halves are therefore as exact as the total, and nothing here is drawn at
// random.
//
// What the split still rests on is the classifier, exactly as the surface area does. The geometry fixes
// the totals and the network labels them.
// One pocket, as the surface around it gives it: not only how much room it holds but where the room is.
//
// The same boundary that gives the volume gives the first moment of the region as well, so the centroid costs
// one more sum per arc and no further geometry. What it is good for is a blocking sphere, which a simulation
// wants written where the pocket is: `freeRadius` is the largest ball about the centroid that holds no atom,
// and such a ball lies wholly in the pocket, being free space and connected and so unable to cross the neck
// that sealed the pocket off. That makes the pair below a blocking sphere that is safe by construction, with
// nothing sampled anywhere in it.
//
// `coveringRadius` is the other end of it: the farthest the pocket's own surface gets from the centroid, so a
// ball of that radius holds the whole pocket and a simulation that rejects everything inside it rejects
// everything the pocket could have held. It too is measured rather than sampled --- on each patch the farthest
// point is either the direction straight away from the centre, if the patch holds it, or on the patch's own
// edges, where the distance along a bounding circle is a single sinusoid --- and it is an upper bound on
// nothing: it is attained. What it does not promise is that the ball reaches no channel. Reaching into the
// framework is harmless, a molecule having no business inside an atom, but a ball wide enough to poke through
// the neck would block part of a pore, and only the pore can say whether it does.
//
// `equivalentRadius`, the radius of the ball of the pocket's own volume, sits between the two and says how
// round the pocket is: all three agree for a spherical cavity, and `coveringRadius` runs away from the others
// for one drawn out along a channel, which is where one ball will not do.
//
// `channelRadius` is what settles the one thing the pocket alone cannot: how far a ball may be taken before it
// reaches somewhere a molecule belongs. The channels are open regions whose closure is bounded by the patches
// that face them, and the centre is not in one, so the nearest point of the accessible void is on one of those
// patches and the distance to it is a minimum of the same kind as the reach. A ball of `blockingRadius` is
// therefore both inside the pocket's own reach and clear of every channel, which is a blocking sphere with
// nothing sampled and nothing assumed; and where the reach is the smaller of the two it covers the pocket
// entirely, which `coversPocket` reports.
export struct PocketGeometry
{
  double3 centre{0.0, 0.0, 0.0};             // Cartesian, in the home cell
  double3 centreFractional{0.0, 0.0, 0.0};
  double volume{0.0};                        // Å³
  double area{0.0};                          // Å²
  double equivalentRadius{0.0};              // Å, of the ball of the same volume
  double freeRadius{0.0};                    // Å, of the largest ball about the centre that holds no atom
  double coveringRadius{0.0};                // Å, of the smallest ball about the centre that holds the pocket

  // Å, to the nearest point of the boundary of the accessible void. Left at the largest double where a
  // structure has no accessible void at all, its whole void being sealed cages, and then nothing caps the
  // reach. Not an infinity: the build turns those off, so a finite marker is the one that survives.
  double channelRadius{std::numeric_limits<double>::max()};

  bool hasChannel() const { return channelRadius < 0.5 * std::numeric_limits<double>::max(); }

  // Whether the centroid lies in a channel rather than in the pocket it is the centroid of. A pocket bent
  // round a corner can have its centroid outside itself: in the framework, which costs nothing, since the ball
  // about it is capped by the pocket's own reach either way; or in the void beyond, where the reasoning behind
  // both radii is measured from a point on the wrong side of the wall and the ball about it would block a pore.
  // Only the second is dangerous, and it is the one the classifier can prove, so it is refused rather than
  // capped and the structure falls back to sampling.
  bool centreInChannel{false};

  // The sphere to block this pocket with: as far as the pocket goes, or as far as it may go without reaching a
  // channel, whichever is the shorter. It is never smaller than `freeRadius`, the nearest wall to the centre
  // being the pocket's own.
  double blockingRadius() const { return std::min(coveringRadius, channelRadius); }

  // Whether that sphere holds the whole pocket, so that nothing is left for a second one.
  bool coversPocket() const { return coveringRadius <= channelRadius; }
};

export struct ExactVoidSplit
{
  double voidVolume{0.0};          // Å³, the whole of it
  double accessibleVolume{0.0};    // Å³, the channels
  double inaccessibleVolume{0.0};  // Å³, the pockets

  std::size_t numberOfPockets{0};  // bounded surfaces enclosing void

  // One entry per pocket, in the order they were counted. Filled by the route that goes through the connected
  // surfaces, a pocket having a centre only where its boundary is a closed surface of its own.
  std::vector<PocketGeometry> pockets;

  // Bounded surfaces enclosing solid rather than void, whose volume was subtracted: a cluster of atoms
  // standing inside a pocket takes up room that is not part of it. None occurs in a framework, whose atoms
  // are all bonded into one network, but the sum does not need to be told so.
  std::size_t numberOfEnclosedSolids{0};

  // Filled by the route that goes through the connected surfaces: how many there were, and how many of them
  // had their side settled by a free ball reached from the surface rather than inferred from proximity.
  std::size_t numberOfSurfaces{0};
  std::size_t provedSurfaces{0};
  std::size_t provedPockets{0};

  // Bounded surfaces where the percolation analysis disagrees with the geometry about which side the void is
  // on, and what that is worth.
  //
  // The geometry decides it. A bounded surface has the void on one side and the solid on the other, and which
  // way round that is fixes the sign of the volume it encloses: the sphere normals point into the region when
  // the void is inside it and away from it when the solid is. Void enclosed by a bounded surface can reach
  // nothing outside it, so a positive volume is itself a proof of a pocket, owing nothing to the network. The
  // network is consulted only for a surface enclosing solid, where the question is whether the void
  // \emph{around} the cluster is sealed.
  //
  // Where the two differ it is the network that is short of information rather than the geometry: a pocket too
  // small to hold a node with room for the probe has no node to be labelled, and the attribution of a point
  // on its wall then falls to whichever node is nearest in line of sight, which is in the neighbouring
  // channel. So these two numbers measure the network, and are reported for that.
  std::size_t signDisagreements{0};
  double signDisagreementVolume{0.0};  // Å³

  // Area the classifier could not place. It is surface belonging to no pore, so a pocket behind it is
  // missing from the sum and its volume has gone silently to the channels. Zero is the only value at
  // which the split above is what it claims to be, which is why it is reported rather than absorbed.
  double undecidedArea{0.0};  // Å²

  // The largest closure defect over the pockets: the length of the integral of the normal over a pocket's
  // boundary, divided by its area. That integral vanishes over any closed surface, so this measures whether
  // each pocket's boundary was in fact assembled whole and in one frame, and it is the check that the
  // periodic bookkeeping is right rather than merely plausible. Where the labelling is consistent it is a
  // residue of the quadrature, parts in 1e11 of the area or below.
  //
  // It is also the only check there is on the classifier's labelling itself. Arcs given to the wrong pore
  // leave both boundaries open and show up here at once; the same mistake in the surface area is invisible,
  // moving area between two pockets without altering either total.
  double closureDefect{0.0};

  // Whether the two above leave the split trustworthy. A caller that wants a number regardless should
  // sample instead of reinterpreting this one.
  bool reliable{true};

  // Why not, in a few words, for the report. Empty when the split stands.
  std::string rejection;
};

// The split, from a sweep that classified its arcs. `patches` must come from such a sweep, or there is
// nothing here to divide.
//
// The reference route, as the arc-by-arc sweep it is built on is: nothing outside the tests calls it, and
// what it is for is to reach the same numbers as the route below from an argument with nothing in common
// with it. Where the two agree the split is as good as the geometry; where they do not, the closure defect
// each reports says which of them failed and on which pocket.
export ExactVoidSplit exactVoidSplit(const PoreAccessibility& accessibility,
                                     const MeasuredPatches& patches, double cellVolume);

// The same, sweeping the patches itself.
export ExactVoidSplit exactVoidSplit(const PoreAccessibility& accessibility, double cellVolume,
                                     std::size_t subdivisions = 1);

// The split taken over the connected surfaces of the boundary instead of over the classifier's pores.
//
// The two differ in what has to be got right. Arc by arc, a pocket's boundary is whatever set of arcs the
// classifier assigned to that pore, and it is a closed surface only if every one of those decisions agreed;
// where they did not, the boundary is open, the divergence theorem does not apply to it, and the volume it
// returns is not the pocket's. Surface by surface, the boundary is a connected surface found by geometry,
// closed or not on its own account, and the classifier says only which side of it the void is on -- a single
// decision for the whole of it, which cannot disagree with itself.
//
// So the closure defect here is a quadrature residue and nothing else, and the classifier's remaining job is
// small enough to be done properly: it is asked at the most open point of each surface, where a free ball
// can usually be reached and the answer proved.
//
// A bounded surface facing a pocket contributes the volume it encloses. The atoms of a cluster floating
// inside a pocket carry such a surface too, with the void outside it rather than inside, and its
// contribution comes out negative -- which is right, that volume being solid and not part of the pocket. A
// bounded surface facing a channel is the outside of a cluster in open void and encloses no void at all, so
// it is left out; the channels take what the subtraction leaves, as before.
export ExactVoidSplit exactVoidSplitByComponents(const PoreAccessibility& accessibility,
                                                 const BoundaryComponents& components,
                                                 const std::vector<ComponentVerdict>& verdicts,
                                                 const MeasuredPatches& patches, double cellVolume);

// The same, decomposing the boundary and sweeping it itself.
export ExactVoidSplit exactVoidSplitByComponents(const PoreAccessibility& accessibility, double cellVolume,
                                                 std::size_t subdivisions = 1);
