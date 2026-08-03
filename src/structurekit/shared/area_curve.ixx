module;

export module grid_area_curve;

import std;

import uint3;
import unit_cell;

// The area of a level set, against the level, for every level at once.
//
// An area is the one quantity here that cannot be freed of the level it was measured at. A diameter, a volume
// and a pore size are all properties of a region, and a region can be pinned; an area is a property of a
// surface, and the surface *is* the level. So the honest thing to report is not one area but the whole family
// of them, and to leave the choice of where to read it off to whoever knows what the number is for.
//
// This is cheap, which is the reason it is worth doing. The obvious way to obtain the family would be to run
// marching cubes once per level, and at a hundred levels that is a hundred passes. The coarea formula gives
// all of them in one:
//
//     \int_\Omega g(A(x)) |grad A(x)| dx  =  \int g(a) S(a) da,
//
// for any weight g, where S(a) is the area of {A = a}. Taking g to be the indicator of a bin turns the left
// side into a sum over the grid points whose value falls in that bin, so gathering |grad A| into bins of A and
// dividing by the bin width gives S(a) over the whole range in a single pass over the field.
//
// The result is an estimator and not the same construction as marching cubes, so the two do not agree to the
// last digit, and the difference between them is a useful thing rather than a nuisance: it is a direct measure
// of how well the grid resolves the surface, and it closes as the grid is refined.
export struct AreaCurvePoint
{
  // The level, in whatever units the field carries.
  double level{0.0};

  // The area of that level set over one unit cell, in Å².
  double area{0.0};

  // The share of the cell lying below the level. It comes free with the pass and it is what says whether a
  // level is measuring the wall of a channel or the inside of a pocket nobody can reach.
  double fractionBelow{0.0};
};

export struct AreaCurve
{
  std::vector<AreaCurvePoint> points;

  double lowestLevel{0.0};
  double highestLevel{0.0};
  double binWidth{0.0};

  // The share of the field falling outside the range asked for, below and above. The high one is never small
  // on an energy field, most of the cell being wall, and it is not a fault; the low one being anything but
  // zero means the range has cut into the part of the landscape that was of interest.
  double fractionBelowRange{0.0};
  double fractionAboveRange{0.0};

  // How far the field moves across one voxel, averaged over the range. It is the resolution of the curve, and
  // asking for bins much finer than it buys nothing: the weight of each voxel is spread over the values it
  // covers, so finer bins are correlated with their neighbors rather than independent of them.
  double meanSpanPerVoxel{0.0};

  double seconds{0.0};

  // The area at a level, interpolated between the two bins that straddle it. Outside the range it is zero.
  double areaAt(double level) const;
};

// Gathers |grad A| into bins of A. The field is taken as periodic, its gradient by central differences on the
// grid, and the whole of the cell is visited. `numberOfBins` bins span [lowestLevel, highestLevel].
export AreaCurve areaAgainstLevel(uint3 gridSize, const UnitCell &unitCell, std::span<const float> field,
                                  double lowestLevel, double highestLevel, std::size_t numberOfBins);
