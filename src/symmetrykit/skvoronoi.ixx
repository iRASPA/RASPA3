module;

export module skvoronoi;

import std;

import int3;
import double3;
import double3x3;

// Metric-aware Voronoi construction with a native fractional-coordinate interface.
//
// The Voronoi tessellation of a periodic system is not invariant under the anisotropic
// scaling that maps a unit cell onto the fractional unit cube, so fractional positions
// cannot simply be tessellated in the unit cube. This class takes positions in fractional
// coordinates together with the unit-cell matrix h (lattice vectors as columns), bins and
// searches neighbors in fractional space, and accounts for the cell geometry exactly:
// candidate separations are mapped through h, so every bisector plane is the true
// Cartesian one (equivalent to cutting in fractional space with the metric tensor
// M = hᵀh, but cheaper). Results are reported in both fractional and Cartesian form.
//
// The construction is standalone and works for arbitrary (orthorhombic or triclinic)
// periodic cells with det(h) != 0.

// A single face of a Voronoi cell; the face is the metric perpendicular bisector between
// the cell's site and the periodic image `neighborImage` of site `neighborIndex`.
export struct SKVoronoiFace
{
  std::size_t neighborIndex;
  int3 neighborImage;
  double neighborDistance;  // Cartesian distance between the two sites [Å]
  std::vector<std::size_t> vertexIndices;  // ordered right-handed around the outward normal
  double area;              // Cartesian face area [Å²]
  double3 normalCartesian;  // outward unit normal in Cartesian space
};

// A Voronoi cell; vertex positions are stored relative to the site position, both in
// fractional coordinates and mapped to Cartesian coordinates.
export struct SKVoronoiCell
{
  std::size_t siteIndex;
  double3 sitePositionFractional;
  std::vector<double3> verticesFractional;  // relative to the site, fractional
  std::vector<double3> verticesCartesian;   // relative to the site, Cartesian [Å]
  std::vector<SKVoronoiFace> faces;
  double volume;             // Cartesian volume [Å³]
  double3 centroidCartesian; // relative to the site, Cartesian [Å]
};

export class SKVoronoi
{
 public:
  // When `radii` is empty the construction is the ordinary (unweighted) Voronoi diagram.
  // When per-site radii are given, the radical (power) diagram is built instead: the plane
  // between sites i and j is offset towards the smaller atom by ½(rᵢ² − rⱼ²)/|Δ|, which is
  // the tessellation used by zeo++/voro++ for radii-dependent pore analysis.
  SKVoronoi(const double3x3 &unitCell, const std::vector<double3> &fractionalPositions,
            const std::vector<double> &radii = {});

  SKVoronoiCell computeCell(std::size_t siteIndex) const;
  std::vector<SKVoronoiCell> computeAllCells() const;

  const double3x3 &unitCell() const { return _unitCell; }
  const double3x3 &metricTensor() const { return _metricTensor; }
  double cellVolume() const { return _cellVolume; }
  const std::vector<double3> &positions() const { return _positions; }

 private:
  double3x3 _unitCell;
  double3x3 _inverseUnitCell;
  double3x3 _metricTensor;  // M = hᵀh
  double _cellVolume;
  double3 _perpendicularWidths;
  std::vector<double3> _positions;           // fractional, wrapped into [0,1)
  std::vector<double3> _cartesianPositions;  // h·s for the wrapped fractional positions [Å]
  std::vector<double> _radiiSquared;         // per-site rᵢ² (all zero for the unweighted case)
  double _maximumRadiusSquared{0.0};         // max rᵢ², used to bound the weighted plane search

  // Cell list over the fractional unit cube; bin shells are walked outward with
  // metric-tensor distance bounds, enumerating periodic images on the fly.
  int3 _gridSize;
  std::vector<std::vector<std::size_t>> _bins;
  double _minimumBinWidth;    // smallest Cartesian perpendicular width of a bin [Å]
  double _initialHalfWidth;   // Cartesian half-width of the initial bounding cube [Å]
};
