module;

export module brute_force_structure;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;

// A crystal as the brute-force checks see it, and the handful of geometric questions they ask of it.
//
// Everything in this directory exists to disagree with the rest of the library. A check that shares an
// argument with the thing it is checking cannot find a mistake in that argument, so nothing here is built
// on the Apollonius diagram, the Voronoi network, the pore classifier or the sampled roadmap; the only
// thing taken from outside is the cell, the atom centres and the radii, which is what was read out of the
// file and is not in dispute.
//
// The methods below are written the slow and obvious way on purpose. Where the exact routes solve for an
// answer, these loop over every atom of every image and take a minimum, and where the exact routes prune
// they do not. That is what makes them worth running: they are wrong in different ways.
export struct BruteForceStructure
{
  std::string name;

  UnitCell unitCell;
  std::vector<double3> positions;  // Cartesian, one unit cell
  std::vector<double> radii;       // Å, the distance from each centre that is out of bounds

  // How many cells out the image search runs, per axis. Chosen so that no image outside the shell can be
  // nearer than one inside it, which for a badly reduced cell takes more than the twenty-seven neighbours
  // the minimum-image convention looks at. It is set in `prepare`.
  int3 shell{1, 1, 1};

  // One image of the cell that the search covers, held in order of how far off it is.
  //
  // The order is what makes the search affordable without making it an assumption. An atom is first brought
  // to its nearest image in the ordinary way, at some displacement d from the point being asked about; the
  // same atom in the image translated by L is then at least |L| - |d| away, by the triangle inequality. So
  // once the search has found some atom at distance s, every image with |L| - |d| - r already past s holds
  // nothing that can beat it, and since the images are in order of |L| the rest of them go together.
  //
  // What is skipped is therefore skipped because it was shown not to matter, which is the difference
  // between this and the minimum-image convention: that asserts one cell is enough, and for a cell far from
  // reduced it is not. The shell searched is wide enough that no image outside it can hold the nearest atom
  // at all, and the bound above is what keeps walking it from costing anything.
  struct Image
  {
    double3 translation;
    double distance;  // Å, |translation|
  };

  std::size_t size() const { return this->positions.size(); }

  static BruteForceStructure create(std::string name, const UnitCell &unitCell,
                                    std::vector<double3> positions, std::vector<double> radii);

  // Works out the image shell and caches the fractional atom positions. Call once, before anything else.
  void prepare();

  // The clearance: how much room there is for a sphere centred at `position`, or how far inside the nearest
  // atom it is when there is none. Minimised over every atom of every image in the shell, with no attempt
  // to stop early.
  double clearance(const double3 &position) const;

  // The same, but stopping as soon as some atom is nearer than `bound`. For the many places that only need
  // to know whether there is room, not how much.
  bool hasRoomFor(const double3 &position, double bound) const;

  // The nearest image displacement from `from` to `to`, searched over the shell rather than assumed.
  double3 nearestImage(const double3 &from, const double3 &to) const;

  // The smallest sphere about `position` that reaches no atom is zero when the position is inside one, so
  // this says which of the two it is without the caller having to compare against a tolerance.
  bool insideAnAtom(const double3 &position) const { return this->clearance(position) < 0.0; }

  // Along the ray from `origin` in unit direction `direction`, the first distance at which the ray is
  // outside every atom, or nothing when it never is within `limit`. Every atom of every image contributes
  // an interval and the intervals are merged; this is what a brute force uses in place of an argument
  // about which atom matters.
  std::optional<double> firstExit(const double3 &origin, const double3 &direction, double limit) const;

  // The largest sphere that can travel in a straight line from `origin` to `origin + displacement`: the
  // minimum over every atom of every image of the distance from that atom to the segment, less its radius.
  double segmentClearance(const double3 &origin, const double3 &displacement) const;

  // The fractional coordinate of a Cartesian point, wrapped into the home cell.
  double3 wrappedFractional(const double3 &position) const;

  const std::vector<double3> &fractional() const { return this->fractionalPositions; }
  const std::vector<Image> &images() const { return this->imageList; }

  double largestRadius() const { return this->maximumRadius; }

 private:
  // The displacement from a point, given in fractional coordinates, to the nearest ordinary image of atom
  // `j`. Every image of that atom is this plus one of the translations above.
  double3 baseDisplacement(const double3 &fractional, std::size_t j) const;

  std::vector<double3> fractionalPositions;
  std::vector<Image> imageList;  // in order of how far off they are, the home cell first
  double maximumRadius{0.0};
  double spread{0.0};  // Å, the longest an ordinary nearest-image displacement can be
};
