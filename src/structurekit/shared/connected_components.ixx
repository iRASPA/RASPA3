module;

export module grid_connected_components;

import std;

import int3;
import uint3;
import double3;

// How many independent directions a set of lattice translations spans, in integers throughout. These are
// lattice translations, so the answer is a property of the integers themselves and there is no tolerance to
// choose: a vector is off a line when its cross product with that line does not vanish, and off a plane when
// its product with the plane's normal does not.
export int latticeVectorRank(const std::vector<int3> &vectors);

// One connected piece of the region left open by a scalar field held at or above some level.
//
// The field is read as an openness, in the same sense as the percolation sweep next door: larger means more
// open, and a point takes part when its value is at least the level given. The clearance is an openness
// already and the level is then a probe radius; an energy is one once negated, and the level is then the
// negative of the highest energy a molecule is allowed to sit at.
//
// The pore is followed until it closes on itself. Every journey that returns it to a periodic image of its
// starting point contributes the lattice translation it came back by, and the number of independent
// directions among those translations is the dimensionality: zero for a pocket sealed off from the rest of
// the crystal, and one, two or three for a channel, a layer, and a network running through it. This is the
// same quantity, and the same definition of it, that the Voronoi and Apollonius routes report.
export struct GridPore
{
  bool isChannel{false};
  int dimensionality{0};
  std::size_t numberOfVoxels{0};

  // Enough of the translations to fix the rank; a pore closing on itself in one direction has one.
  std::vector<int3> windingVectors;

  // The most open point of this pore, which for a clearance field is the largest sphere a channel can hold.
  double largestOpenness{0.0};
  int3 widestVoxel{0, 0, 0};

  // The centre of the pore, in fractional coordinates. It is taken over the pore followed out of the cell
  // rather than over the wrapped points, so a pocket lying across a cell boundary gets its own centre and not
  // a point halfway between its two halves. It means nothing for a channel, which has no centre to find.
  double3 centroidFractional{0.0, 0.0, 0.0};

  // A voxel of the pore, in grid coordinates, so the pore can be pointed at in a report.
  int3 exampleVoxel{0, 0, 0};
};

// The pores of an openness field at one level: the connected pieces of the set of points at or above it,
// together with what each of them turns out to be.
export struct GridComponents
{
  double lowestOpenness{0.0};

  // Which pore each grid point belongs to, or -1 where the field falls below the level.
  std::vector<std::int32_t> voxelPore;

  // Ordered channels first, then by how much room each pore holds, so the report reads from the pores that
  // matter down to the ones that barely exist.
  std::vector<GridPore> pores;
  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};

  // Pockets holding a single grid point. These are places where the probe only just fits and which hold no
  // volume worth the name, and a finer grid finds more of them rather than fewer, because in the limit they
  // are the isolated points where the probe fits exactly. They are counted as pockets above, and reported
  // apart from the rest so that a pocket count can be read either way.
  std::size_t numberOfSinglePointPockets{0};

  // The dimensionality of the pore system as a whole: the widest-running channel in it.
  int dimensionality{0};

  std::size_t numberOfOpenVoxels{0};
  double seconds{0.0};

  GridComponents();
  ~GridComponents();

  // The volume fraction of the cell that is open at this level, and the part of it in channels.
  double occupiedFraction() const;
  double channelFraction() const;

  static GridComponents compute(uint3 gridSize, std::span<const float> openness, double lowestOpenness);
};
