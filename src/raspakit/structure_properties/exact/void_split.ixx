module;

export module exact_void_split;

import std;

import voronoi_accessibility;
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
export struct ExactVoidSplit
{
  double voidVolume{0.0};          // Å³, the whole of it
  double accessibleVolume{0.0};    // Å³, the channels
  double inaccessibleVolume{0.0};  // Å³, the pockets

  std::size_t numberOfPockets{0};  // bounded surfaces enclosing void

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

// The split, from a classified sweep of the patches. `patches` must come from a sweep with the arcs
// classified, or there is nothing here to divide.
export ExactVoidSplit exactVoidSplit(const VoronoiAccessibility& accessibility,
                                     const ExactSurfaceAreaSample& patches, double cellVolume);

// The same, sweeping the patches itself.
export ExactVoidSplit exactVoidSplit(const VoronoiAccessibility& accessibility, double cellVolume,
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
export ExactVoidSplit exactVoidSplitByComponents(const VoronoiAccessibility& accessibility,
                                                 const BoundaryComponents& components,
                                                 const std::vector<ComponentLabel>& labels,
                                                 const ExactSurfaceAreaSample& patches, double cellVolume);

// The same, decomposing the boundary and sweeping it itself.
export ExactVoidSplit exactVoidSplitByComponents(const VoronoiAccessibility& accessibility, double cellVolume,
                                                 std::size_t subdivisions = 1);
