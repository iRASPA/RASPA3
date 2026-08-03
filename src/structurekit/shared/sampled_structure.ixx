module;

export module sampled_structure;

import std;

import double3;
import unit_cell;

// The narrowest point of a straight way through the structure: how much room there is there, and where.
export struct SegmentBottleneck
{
  double radius{0.0};  // Å, negative when the way is blocked and this is how deep into an atom it goes
  double3 position{};  // Cartesian, the point of the segment where that is measured
};

// A crystal as the samplers see it: a cell, the centres of the atoms of one unit cell, and the radius of
// the sphere around each of them that the thing being sampled may not enter.
//
// The routines that throw points at a structure need nothing else. They never ask what element an atom is,
// what it is bonded to or what it interacts with, only where it is and how far away it has to be kept, so
// taking a framework and a force field would be taking two objects to read four numbers off them. Splitting
// them here is what keeps the samplers testable against a hand-written cell and a pair of spheres, and what
// lets the caller decide what "in contact" means: half an atom's own diameter for a point probe, the mixed
// diameter of the pair for a probe with a size of its own.
//
// The metadata below the geometry is not used in any calculation. It is what the reports are headed with,
// and it travels with the structure so that a sampler can write its own file without knowing where the
// structure came from.
export struct SampledStructure
{
  std::string name;                     // also the stem of the file a sampler writes
  std::size_t spaceGroupHallNumber{1};

  UnitCell unitCell;
  std::vector<double3> positions;  // Cartesian, one unit cell
  std::vector<double> radii;       // Å, the distance from each centre at which contact is made

  double mass{0.0};  // g/mol, one unit cell

  std::size_t size() const { return positions.size(); }

  // kg/m³.
  double density() const;

  // What an area in Å² is multiplied by to become one in m²/g.
  double gravimetricFactor() const;

  // What a volume in Å³ is multiplied by to become one in cm³/g.
  double gravimetricVolumeFactor() const;

  // True when `position` lies within the contact radius of some atom other than `skip`, under the minimum
  // image convention. Pass an index at or past the end to test against every atom.
  bool overlaps(const double3 &position, std::size_t skip) const;

  // The largest sphere centred on `position` that reaches no atom, or nothing when the position is inside
  // one. This is the clearance min_j(|x - x_j| - r_j), which is what a pore size is measured with.
  std::optional<double> freeRadius(const double3 &position) const;

  // The same clearance, signed and always returned: how deep inside the nearest atom `position` is when it
  // is inside one. What that buys over `freeRadius` is a value to compare when there is nothing to measure
  // yet, so that a walk that starts inside an atom can still be told which way is out.
  double clearance(const double3 &position) const;

  // The largest sphere that can travel in a straight line from `position` to `position + displacement`, and
  // where along the way it is hemmed in most closely. The radius is negative when no sphere can pass.
  //
  // The clearance along a segment is a minimum of distances to points, each of which is smallest somewhere
  // definite, so this is exact rather than sampled: min over atoms of the distance from the atom to the
  // segment, less that atom's radius.
  //
  // `displacement` is a Cartesian vector and is taken as given, not wrapped, so a caller crossing a boundary
  // passes the minimum image of the two endpoints' difference and gets the segment through the boundary. The
  // nearest image of an atom is taken to the midpoint and held for the whole segment, which can only pick
  // the wrong image for an atom half a cell away, and such an atom is never the nearest.
  SegmentBottleneck segmentBottleneck(const double3 &position, const double3 &displacement) const;

  // Its radius alone, for the callers that have nothing to do with where it is.
  double segmentClearance(const double3 &position, const double3 &displacement) const;

  // The comment block every report here begins with: what was measured and how big it is.
  void writeHeader(std::ostream &stream) const;
};

// The probe a surface area is measured with, for the header of the report. What it is kept apart from the
// per-atom contact radii for is that those have the mixing rule already applied and cannot be read back.
export struct SampledProbe
{
  std::string name;
  double sizeParameter{0.0};    // Å
  double wellDepthFactor{1.0};  // what the contact radii were scaled by

  void writeHeader(std::ostream &stream) const;
};
