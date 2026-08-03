module;

export module grid_percolation;

import std;

import int3;
import uint3;

// A monotone sweep of a scalar field over the grid, which answers the question of when the field's
// super-level sets begin to run through the crystal rather than closing inside one cell.
//
// The field is read as an openness: larger means more open, and the sweep switches points on from the most
// open downwards, joining each to the neighbours already on. The moment some region closes on a periodic
// image of itself, the sweep has found the least open point that still has to be crossed to get through,
// and that value is what percolation costs.
//
// Both fields of the grid route are of this kind once the sign is right. The clearance is an openness
// already, and a probe passes where it is large. An energy is an openness once negated, and a molecule
// passes where it is small; the barrier to getting through is then the negative of the value reported here.
// Sweeping is what both have in common, so it lives here rather than in either of them.
export struct GridPercolation
{
  // Whether any region ever closed on a periodic image of itself.
  bool percolates{false};

  // The most open point anywhere in the field.
  double highestOpenness{0.0};

  // The openness at which a region first ran through the crystal: the least open point on the most open
  // route through. Everything more open than this is confined to one cell.
  double percolationOpenness{0.0};

  // The most open point of the region that first ran through, which is reachable without ever passing a
  // point less open than `percolationOpenness`.
  double highestOpennessOnPath{0.0};

  // The openness at which the pore system last ran in at least one, two, and three directions. Lowest
  // double where it never runs that far, so that a caller which negates gets a sensible symbol back.
  std::array<double, 3> opennessByDimension{};

  // How many independent lattice translations the first percolating region had, at the moment it closed.
  int dimensionalityAtThreshold{0};

  // How many points took part, which is those at or above the floor given to the sweep.
  std::size_t numberOfVoxels{0};

  // For each point, the openness of the tightest place on the most open route out of the crystal: the level
  // to which the field has to be allowed to fall before the point is joined to a region that runs through.
  // Empty unless asked for.
  //
  // It is the same question `percolationOpenness` answers, asked of one point rather than of the field. The
  // two agree at the point where the sweep first broke through, and everywhere else this is the smaller,
  // some other point having to be crossed first. Its use is that the difference between a point's own
  // openness and this is what has to be given up to leave, which on a clearance field is how much narrower
  // the tightest neck on the way out is than the room here, and on an energy field is the barrier to escape.
  //
  // Lowest double where no route out exists at any level, which happens only when nothing percolates at all.
  std::vector<float> escapeOpenness;

  double seconds{0.0};
};

// Sweeps `openness` over a grid of `gridSize` points, x varying fastest, taking part only in the points
// whose value is at least `lowestOpenness`. Neighbours are the six face neighbours, wrapped periodically.
//
// `recordEscape` fills `escapeOpenness` as well. It costs one more pass' worth of bookkeeping, the members
// of the regions that have not yet broken through being carried along in lists as the sweep goes, so it is
// off unless a caller wants it.
export GridPercolation sweepPercolation(uint3 gridSize, std::span<const float> openness, float lowestOpenness,
                                        bool recordEscape = false);
