module;

export module voronoi_blocking_spheres;

import std;

import double3;
import unit_cell;
import crystal;
import pair_interactions;
import pore_accessibility;
import exact_void_split;

// What a blocking sphere is and how it is written out is the same wherever the spheres came from.
// Re-exported because a caller asking this route for them is asking for that type back.
export import blocking_spheres;

// Spheres covering the pockets a probe cannot reach, for a simulation to reject insertions in.
//
// A blocking sphere has two jobs and they pull against each other: it has to hold the pocket, or a molecule
// still gets into the part it misses, and it may not reach any channel, or the simulation loses part of its own
// pore. Both are settled by the surface, pocket by pocket, and neither is sampled.
//
// The boundary of a pocket is a closed surface, so the divergence theorem gives its volume and, with one more
// moment, its centroid; and about that centroid two distances are extrema over the same patches, each of them a
// sinusoid along a bounding circle and so available in closed form. The first is how far the pocket's own wall
// gets from the centre, which is the radius that holds the pocket exactly. The second is how near the boundary
// of the accessible void comes, which is the radius past which a sphere would block a pore. The lesser of the
// two is the sphere, and it needs no allowance for the resolution of anything, there being no resolution.
//
// Two things follow that the sampled route below could not have. Every pocket gets a sphere, however small it
// is, since a pocket is found from its boundary rather than by a point landing inside it. And a pocket gets one
// sphere: over the zeolite atlas, at a helium probe and at a nitrogen one, the covering radius is the lesser of
// the two everywhere --- the channel is never nearer to a pocket's centre than the pocket's own wall is, by
// 0.33 Å at the tightest and about 4 Å in the median --- so a single sphere holds the whole pocket and there is
// nothing left over to cover.
//
// The one case the argument does not reach is a pocket bent far enough round the framework for its centroid to
// lie in the void beyond it, where both radii are measured from the wrong side of a wall. That is refused
// rather than capped, and the structure falls back to sampling; it does not arise on the atlas.
// One sphere per pocket of a measured split: at the pocket's centre, of its blocking radius.
export std::vector<BlockingSphere> exactBlockingSpheres(const ExactVoidSplit& split);

// Why the measured spheres cannot be written for a structure, and empty where they can. It is the split's own
// rejection, plus the one thing that is a reason here and nowhere else: a pocket whose centroid the classifier
// puts in a channel, which is a centre no radius makes safe.
export std::string measuredSpheresRefused(const ExactVoidSplit& split);

// The sampled route, kept for the structures a measured split has to refuse and for the comparison.
//
// Points of the void that a probe cannot reach are grouped by pocket and each pocket covered by spheres taking
// the widest first. A sphere reaches as far as the pocket does and no further than the nearest position the
// probe can occupy from a channel. Both are measured between sample points, which are positions of the probe's
// *centre*, so neither needs a further allowance for the probe's size -- an earlier version added one to the
// first and subtracted one from the second, and since the two limits are typically only a few Å apart that left
// nothing in between: the radii collapsed and a pocket came back covered by hundreds of spheres too small to
// block anything, or was abandoned altogether.
//
// Below both limits sits a floor that holds regardless: a sphere of the clearance at its centre holds no atom,
// so it lies in the free space, and being connected it cannot cross the neck that made the pocket a pocket.
// Such a sphere is inside the pocket by construction. That floor is the sampled route's version of the argument
// the measured one makes throughout.
//
// It runs over whatever accessibility classifier is handed to it, so that the same spheres can be found from a
// network taken from the radical diagram or from the Apollonius diagram. The probe radius is not a parameter: it
// is already in the inflated radii the classifier carries, which is what the clearance is measured against.
export std::vector<BlockingSphere> computeBlockingSpheres(const PoreAccessibility& accessibility,
                                                          std::size_t numberOfSamples);

// Writes the `.block` file a simulation reads back, which is a count and one line of fractional centre and
// radius per sphere and nothing else, and beside it an account of where the spheres came from.
export void writeBlockingSpheres(const std::string& frameworkName, const std::string& diagramName,
                                 const std::string& probeName, double probeRadius,
                                 const std::vector<BlockingSphere>& spheres, const std::vector<PocketGeometry>& pockets,
                                 const std::string& fallbackReason);

export struct VoronoiBlockingSpheres
{
  std::vector<BlockingSphere> spheres;

  // Whether the spheres came from the surface or from sampling, and why not from the surface when they did not.
  bool measured{false};
  std::string fallbackReason;

  // The pockets the spheres stand for, when they were measured.
  std::vector<PocketGeometry> pockets;

  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
