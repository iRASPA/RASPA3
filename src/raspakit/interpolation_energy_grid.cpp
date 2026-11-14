module;

module interpolation_energy_grid;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import stringutils;
import units;
import polint;
import atom;
import framework;
import simulationbox;
import interactions_framework_molecule_grid;
import interactions_external_field_grid;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

// For a framework that is kept rigid it is effecient to precompute the energy and forces.
// The amount of points is compute from
//   (a) option 'SpacingVDWGrid 0.15' for Van der Waal grids
//   (b) option 'SpacingCoulombGrid 0.15' for the electrostatic grid
// 0.15 Angstrom per point seems to be a good value
// the sizes of the grids can be large, so due to memory-restrictions the amount and type of pseudo-atoms to
// be precomputed needs to chosen with care.

// for 2D bicubic interpolation, see Numerical Recipes 3.6 Interpolation in more dimensions
// the code here is the extension to 3D using the scheme of Lekien and Marsden
//   F. Lekien and J. Marsden, "Tricubic interpolation in three dimensions"
//   International Journal for Numerical Methods in Engineering, 63, 455-471, 2005
// and its quintic extension
//   H.A. Boateng and K. Bradach, "Triquintic interpolation in three dimensions"
//   Journal of Computational and Applied Mathematics, 430, 115254, 2023
//
// Both are local schemes: the interpolant on a cell is the unique polynomial (degree 3, resp. 5, in each variable)
// that matches the stored data (energy and its mixed derivatives) at the eight corners. Lekien and Marsden write
// it through a 64x64 matrix (216x216 for the quintic scheme) that maps the corner data to the 64 (216) monomial
// coefficients. That matrix-vector product is by far the most expensive part of an evaluation. The same
// interpolant is the tensor product of one-dimensional Hermite interpolants, so it is evaluated here directly as
//   U(t) = sum over corners c and stored derivatives d of  X[d][c] * h_dx(t_x; c_x) * h_dy(t_y; c_y) * h_dz(t_z; c_z)
// with h the cubic (quintic) Hermite basis functions on [0, 1]: 64 (216) multiply-adds instead of a few thousand.
// The gradient and the Hessian follow from the derivatives of the basis functions. The two formulations agree to
// round-off, and the interpolant keeps its C^1 (C^2) continuity across cells.

namespace
{
namespace Hermite
{
// One-dimensional Hermite bases on [0, 1]. value[d][c] is the polynomial whose d-th derivative is 1 at corner c
// (c = 0 at t = 0, c = 1 at t = 1) and whose other (derivative, corner) data are zero; first and second are its
// derivatives with respect to t. Order 2 is the cubic basis (f, f'), order 3 the quintic one (f, f', f'').
template <std::size_t Order>
struct Basis1D
{
  std::array<std::array<double, 2>, Order> value;
  std::array<std::array<double, 2>, Order> first;
  std::array<std::array<double, 2>, Order> second;
};

[[clang::always_inline]] inline Basis1D<2> cubic(double t)
{
  const double t2 = t * t, t3 = t2 * t;
  Basis1D<2> b;
  b.value[0][0] = 2.0 * t3 - 3.0 * t2 + 1.0;
  b.value[0][1] = -2.0 * t3 + 3.0 * t2;
  b.value[1][0] = t3 - 2.0 * t2 + t;
  b.value[1][1] = t3 - t2;
  b.first[0][0] = 6.0 * t2 - 6.0 * t;
  b.first[0][1] = -6.0 * t2 + 6.0 * t;
  b.first[1][0] = 3.0 * t2 - 4.0 * t + 1.0;
  b.first[1][1] = 3.0 * t2 - 2.0 * t;
  b.second[0][0] = 12.0 * t - 6.0;
  b.second[0][1] = -12.0 * t + 6.0;
  b.second[1][0] = 6.0 * t - 4.0;
  b.second[1][1] = 6.0 * t - 2.0;
  return b;
}

[[clang::always_inline]] inline Basis1D<3> quintic(double t)
{
  const double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
  Basis1D<3> b;
  b.value[0][0] = 1.0 - 10.0 * t3 + 15.0 * t4 - 6.0 * t5;
  b.value[0][1] = 10.0 * t3 - 15.0 * t4 + 6.0 * t5;
  b.value[1][0] = t - 6.0 * t3 + 8.0 * t4 - 3.0 * t5;
  b.value[1][1] = -4.0 * t3 + 7.0 * t4 - 3.0 * t5;
  b.value[2][0] = 0.5 * (t2 - 3.0 * t3 + 3.0 * t4 - t5);
  b.value[2][1] = 0.5 * (t3 - 2.0 * t4 + t5);
  b.first[0][0] = -30.0 * t2 + 60.0 * t3 - 30.0 * t4;
  b.first[0][1] = 30.0 * t2 - 60.0 * t3 + 30.0 * t4;
  b.first[1][0] = 1.0 - 18.0 * t2 + 32.0 * t3 - 15.0 * t4;
  b.first[1][1] = -12.0 * t2 + 28.0 * t3 - 15.0 * t4;
  b.first[2][0] = 0.5 * (2.0 * t - 9.0 * t2 + 12.0 * t3 - 5.0 * t4);
  b.first[2][1] = 0.5 * (3.0 * t2 - 8.0 * t3 + 5.0 * t4);
  b.second[0][0] = -60.0 * t + 180.0 * t2 - 120.0 * t3;
  b.second[0][1] = 60.0 * t - 180.0 * t2 + 120.0 * t3;
  b.second[1][0] = -36.0 * t + 96.0 * t2 - 60.0 * t3;
  b.second[1][1] = -24.0 * t + 84.0 * t2 - 60.0 * t3;
  b.second[2][0] = 0.5 * (2.0 - 18.0 * t + 36.0 * t2 - 20.0 * t3);
  b.second[2][1] = 0.5 * (6.0 * t - 24.0 * t2 + 20.0 * t3);
  return b;
}

// The (x, y, z) derivative orders of the i-th quantity stored at every node, i.e. of X[8 * i + corner] below.
// This is the order in which 'makeFrameworkInterpolationGrid' stores them.
constexpr std::array<std::array<std::uint8_t, 3>, 8> tricubicDerivatives = {
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};
constexpr std::array<std::array<std::uint8_t, 3>, 27> triquinticDerivatives = {
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1},
     {0, 0, 2}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 1, 1}, {0, 2, 1}, {1, 0, 2}, {0, 1, 2}, {2, 2, 0},
     {2, 0, 2}, {0, 2, 2}, {2, 1, 1}, {1, 2, 1}, {1, 1, 2}, {2, 2, 1}, {2, 1, 2}, {1, 2, 2}, {2, 2, 2}}};

struct Result
{
  double value{};
  double3 gradient{};
  double3x3 hessian{};
};

// The quantity index of the derivative pattern (dx, dy, dz) of a node (the inverse of the tables above), at
// (dx Order + dy) Order + dz.
template <std::size_t Order, std::size_t N>
constexpr std::array<std::uint8_t, Order * Order * Order> quantityIndex(
    const std::array<std::array<std::uint8_t, 3>, N>& derivatives)
{
  std::array<std::uint8_t, Order * Order * Order> index{};
  for (std::size_t i = 0; i < N; ++i)
  {
    index[(derivatives[i][0] * Order + derivatives[i][1]) * Order + derivatives[i][2]] = static_cast<std::uint8_t>(i);
  }
  return index;
}
constexpr std::array<std::uint8_t, 8> tricubicQuantity = quantityIndex<2>(tricubicDerivatives);
constexpr std::array<std::uint8_t, 27> triquinticQuantity = quantityIndex<3>(triquinticDerivatives);

// The tensor-product Hermite interpolant of the corner data X at the cell coordinate whose one-dimensional bases
// are bx, by, bz. Level 0: the value; 1: also the gradient; 2: also the Hessian (both with respect to the cell
// coordinate t in [0, 1]^3).
//
// With a = (dx, cx), b = (dy, cy), c = (dz, cz) the corner data form the tensor C[a][b][c] =
// X[8 quantity(dx, dy, dz) + cx + 2 cy + 4 cz] and the interpolant is sum_abc C[a][b][c] Bx[a] By[b] Bz[c]. It is
// contracted one axis at a time, x then y then z: M^3 + M^2 + M products per basis combination (M = 2 Order)
// instead of 3 M^3 for the flat sum over all terms, and the derivative bases reuse the partial contractions.
template <int Level, std::size_t Order>
[[clang::always_inline]] inline Result evaluate(const std::array<double, 8 * Order * Order * Order>& X,
                                                const std::array<std::uint8_t, Order * Order * Order>& quantity,
                                                const Basis1D<Order>& bx, const Basis1D<Order>& by,
                                                const Basis1D<Order>& bz)
{
  constexpr std::size_t M = 2 * Order;

  // contraction over x with Bx, Bx' and Bx'': every entry a short reduction over a = (dx, cx) held in registers
  // (the trip counts are compile-time constants: fully unrolled, every index is a constant and every partial sum a
  // register)
  std::array<std::array<double, M>, M> cv{}, cx{}, cxx{};
#pragma clang loop unroll(full)
  for (std::size_t oy = 0; oy < Order; ++oy)
  {
#pragma clang loop unroll(full)
    for (std::size_t cornerY = 0; cornerY < 2; ++cornerY)
    {
#pragma clang loop unroll(full)
      for (std::size_t oz = 0; oz < Order; ++oz)
      {
#pragma clang loop unroll(full)
        for (std::size_t cornerZ = 0; cornerZ < 2; ++cornerZ)
        {
          const std::size_t offset = 2 * cornerY + 4 * cornerZ;
          double sv = 0.0, sx = 0.0, sxx = 0.0;
#pragma clang loop unroll(full)
          for (std::size_t ox = 0; ox < Order; ++ox)
          {
            const std::size_t base = 8 * quantity[(ox * Order + oy) * Order + oz] + offset;
            const double v0 = X[base], v1 = X[base + 1];
            sv += v0 * bx.value[ox][0] + v1 * bx.value[ox][1];
            if constexpr (Level >= 1) sx += v0 * bx.first[ox][0] + v1 * bx.first[ox][1];
            if constexpr (Level >= 2) sxx += v0 * bx.second[ox][0] + v1 * bx.second[ox][1];
          }
          const std::size_t b = 2 * oy + cornerY, cc = 2 * oz + cornerZ;
          cv[b][cc] = sv;
          if constexpr (Level >= 1) cx[b][cc] = sx;
          if constexpr (Level >= 2) cxx[b][cc] = sxx;
        }
      }
    }
  }

  // contraction over y: the value, the x and y derivatives, and the xx, xy and yy second derivatives
  std::array<double, M> v{}, gx{}, gy{}, hxx{}, hxy{}, hyy{};
#pragma clang loop unroll(full)
  for (std::size_t b = 0; b < M; ++b)
  {
    const double hy = by.value[b / 2][b % 2];
    const double dy = by.first[b / 2][b % 2];
    const double ddy = by.second[b / 2][b % 2];
#pragma clang loop unroll(full)
    for (std::size_t cc = 0; cc < M; ++cc)
    {
      v[cc] += cv[b][cc] * hy;
      if constexpr (Level >= 1)
      {
        gx[cc] += cx[b][cc] * hy;
        gy[cc] += cv[b][cc] * dy;
      }
      if constexpr (Level >= 2)
      {
        hxx[cc] += cxx[b][cc] * hy;
        hxy[cc] += cx[b][cc] * dy;
        hyy[cc] += cv[b][cc] * ddy;
      }
    }
  }

  // contraction over z
  Result r{};
#pragma clang loop unroll(full)
  for (std::size_t cc = 0; cc < M; ++cc)
  {
    const double hz = bz.value[cc / 2][cc % 2];
    const double dz = bz.first[cc / 2][cc % 2];
    const double ddz = bz.second[cc / 2][cc % 2];
    r.value += v[cc] * hz;
    if constexpr (Level >= 1)
    {
      r.gradient.x += gx[cc] * hz;
      r.gradient.y += gy[cc] * hz;
      r.gradient.z += v[cc] * dz;
    }
    if constexpr (Level >= 2)
    {
      r.hessian.ax += hxx[cc] * hz;
      r.hessian.ay += hxy[cc] * hz;
      r.hessian.az += gx[cc] * dz;
      r.hessian.by += hyy[cc] * hz;
      r.hessian.bz += gy[cc] * dz;
      r.hessian.cz += v[cc] * ddz;
    }
  }
  return r;
}
}  // namespace Hermite

// The data of the eight corners of the cell with lower corner (x0, y0, z0): N quantities per node, stored
// node-contiguous in the order X[8 * quantity + corner], corners ordered x fastest.
template <std::size_t N>
[[clang::always_inline]] inline void gatherCorners(const double* data, const uint3& numberOfGridPoints, std::size_t x0,
                                                   std::size_t y0, std::size_t z0, std::array<double, 8 * N>& X)
{
  const std::size_t strideY = N * numberOfGridPoints.x;
  const std::size_t strideZ = strideY * numberOfGridPoints.y;
  const double* node = data + N * x0 + strideY * y0 + strideZ * z0;
  const double* corner[8] = {
      node,           node + N,           node + strideY,           node + strideY + N,
      node + strideZ, node + strideZ + N, node + strideZ + strideY, node + strideZ + strideY + N};
  for (std::size_t i = 0; i < N; ++i)
  {
    for (std::size_t c = 0; c < 8; ++c)
    {
      X[8 * i + c] = corner[c][i];
    }
  }
}

inline void prefetchLines(const double* begin, std::size_t count)
{
#if defined(__clang__) || defined(__GNUC__)
  for (std::size_t offset = 0; offset < count; offset += 8)
  {
    __builtin_prefetch(begin + offset);
  }
  __builtin_prefetch(begin + count - 1);
#else
  (void)begin;
  (void)count;
#endif
}

template <std::size_t N>
inline void prefetchCorners(const double* data, const uint3& numberOfGridPoints, std::size_t x0, std::size_t y0,
                            std::size_t z0)
{
  const std::size_t strideY = N * numberOfGridPoints.x;
  const std::size_t strideZ = strideY * numberOfGridPoints.y;
  const double* node = data + N * x0 + strideY * y0 + strideZ * z0;
  // two x-neighbours are adjacent in memory: one range of 2N values per (y, z) pair
  prefetchLines(node, 2 * N);
  prefetchLines(node + strideY, 2 * N);
  prefetchLines(node + strideZ, 2 * N);
  prefetchLines(node + strideZ + strideY, 2 * N);
}
}  // namespace

InterpolationEnergyGrid::InterpolationEnergyGrid(const SimulationBox unitCellBox, double3 origin, uint3 numberOfGridPoints, ForceFieldSettings::InterpolationScheme order)
      : unitCellBox(unitCellBox),
        origin(origin),
        numberOfGridPoints(numberOfGridPoints),
        numberOfCells(numberOfGridPoints.x - 1, numberOfGridPoints.y - 1, numberOfGridPoints.z -1),
        order(order),
        data(std::to_underlying(order) *
             (numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
{
}

const uint3 InterpolationEnergyGrid::parseExternalFieldGridDimensions(const std::string& filename)
{
  int3 gridDims;
  double ax, ay, az, bx, by, bz, cx, cy, cz;  // dummy variables currently

  std::ifstream infile(filename);

  if (!infile) {
      throw std::runtime_error("Cannot open cube file: " + filename);
  }

  std::string line;

  // skip comment lines
  std::getline(infile, line);
  std::getline(infile, line);

  if (!std::getline(infile, line))
  {
      throw std::runtime_error("Invalid cube file (missing origin line): " + filename);
  }

  // read X basis vector coordinates
  if (!std::getline(infile, line))
      throw std::runtime_error("Invalid cube file (missing X grid line): " + filename);
  {
      std::istringstream iss(line);
      if (!(iss >> gridDims.x >> ax >> ay >> az))
          throw std::runtime_error("Invalid X grid line in cube: " + filename);
  }

  // Y vector
  if (!std::getline(infile, line))
      throw std::runtime_error("Invalid cube file (missing Y grid line): " + filename);
  {
      std::istringstream iss(line);
      if (!(iss >> gridDims.y >> bx >> by >> bz))
          throw std::runtime_error("Invalid Y grid line in cube: " + filename);
  }

  // Z vector
  if (!std::getline(infile, line))
      throw std::runtime_error("Invalid cube file (missing Z grid line): " + filename);
  {
      std::istringstream iss(line);
      if (!(iss >> gridDims.z >> cx >> cy >> cz))
          throw std::runtime_error("Invalid Z grid line in cube: " + filename);
  }

  // The sign of a count gives the unit of its step vector, Bohr when positive (what Gaussian, CP2K, ORCA, Multiwfn
  // and Psi4 write) and Angstrom when negative; the steps are not used here (the grid spans the simulation cell), so
  // either is accepted. Only a mix of the two is malformed.
  if (gridDims.x < 0 && gridDims.y < 0 && gridDims.z < 0)
  {
    return {static_cast<std::size_t>(-1 * gridDims.x), static_cast<std::size_t>(-1 * gridDims.y), static_cast<std::size_t>(-1 * gridDims.z)};
  }
  if (gridDims.x > 0 && gridDims.y > 0 && gridDims.z > 0)
  {
    return {static_cast<std::size_t>(gridDims.x), static_cast<std::size_t>(gridDims.y),
            static_cast<std::size_t>(gridDims.z)};
  }
  throw std::runtime_error("All grid dimensions in the CUBE file must be only negative or only positive, not mixed");

  std::unreachable();
}

/**
 * \brief Parses external field grid based on the provided CUBE file.
 *
 * This function extracts parameters specific to external field grids,
 * allowing to use this data in construction of external field interpolation grid in simulations.
 *
 * \param filename The name of the CUBE file containing the external field data (in Kelvin). As in the CUBE format,
 * the voxels are listed with x as the slowest and z as the fastest index.
 * \return The energies of the grid points stored with x as the fastest index (std::layout_left), like the grid.
 */
const std::vector<double> InterpolationEnergyGrid::parseExternalFieldGridCube(const std::string& filename)
{
  uint3 gridDims = parseExternalFieldGridDimensions(filename);

  std::ifstream infile(filename);

  // skip the two comment lines, read the number of atoms from the origin line, and skip the three grid lines
  std::string line;
  std::getline(infile, line);
  std::getline(infile, line);
  std::getline(infile, line);
  std::int64_t numberOfAtoms{};
  if (!(std::istringstream(line) >> numberOfAtoms))
  {
    throw std::runtime_error("Invalid cube file (missing number of atoms): " + filename);
  }
  for (int i = 0; i < 3; ++i)
  {
    std::getline(infile, line);
  }

  // skip the atom lines; a negative number of atoms announces one more line (the data-set identifiers)
  const std::int64_t numberOfAtomLines = numberOfAtoms < 0 ? 1 - numberOfAtoms : numberOfAtoms;
  for (std::int64_t i = 0; i < numberOfAtomLines; ++i)
  {
    if (!std::getline(infile, line))
    {
      throw std::runtime_error("Invalid cube file (missing atom lines): " + filename);
    }
  }

  std::vector<double> gridFlat(gridDims.x * gridDims.y * gridDims.z);
  double value{};

  for (std::size_t ix = 0; ix < gridDims.x; ix++)
  {
    for (std::size_t iy = 0; iy < gridDims.y; iy++)
    {
      for (std::size_t iz = 0; iz < gridDims.z; iz++)
      {
        if (!(infile >> value)){
          throw std::runtime_error("Unexpected end of cube data: " + filename);
        }
        std::size_t idx = ix + gridDims.x * (iy + gridDims.y * iz);

        gridFlat[idx] = value * Units::KelvinToEnergy;
      }
    }
  }

  return gridFlat;
}

void InterpolationEnergyGrid::makeExternalFieldInterpolationGrid(std::ostream& stream, const ForceField& forceField, const SimulationBox &simulationBox)
{
  if (forceField.potentialEnergySurfaceType == ForceField::PotentialEnergySurfaceType::GridFile)
  {
    // the file holds only the energies, while the Hermite schemes store the derivatives at every node as well
    if (order != ForceFieldSettings::InterpolationScheme::Polynomial)
    {
      throw std::runtime_error(
          "External field 'GridFile' holds energies only and requires 'InterpolationScheme' 'Polynomial' (no "
          "molecular dynamics or gradient-based Monte Carlo moves)\n");
    }
    data = InterpolationEnergyGrid::parseExternalFieldGridCube(forceField.externalFieldGridFileName);
  }
  else
  {
    double3 delta = double3(1.0 / static_cast<double>(numberOfCells.x), 1.0 / static_cast<double>(numberOfCells.y),
                            1.0 / static_cast<double>(numberOfCells.z));
    std::mdspan<double, std::dextents<std::size_t, 4>, std::layout_left> data_cell(
        data.data(), std::to_underlying(order), numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

    double percent = 100.0 / static_cast<double>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z);

    double counter{};
    std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
    for (std::size_t k = 0; k < static_cast<std::size_t>(numberOfGridPoints.z); ++k)
    {
      for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfGridPoints.y); ++j)
      {
        for (std::size_t i = 0; i < static_cast<std::size_t>(numberOfGridPoints.x); ++i)
        {
          counter += 1.0;
          double3 s = double3(static_cast<double>(i) * delta.x, static_cast<double>(j) * delta.y,
                              static_cast<double>(k) * delta.z);
          double3 pos = simulationBox.cell * s;

          switch (order)
          {
            case ForceFieldSettings::InterpolationScheme::Polynomial:
            {
              std::array<double, 8> values = Interactions::calculateTricubicFractionalAtPositionExternalField(
                  forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

              data_cell[0, i, j, k] =
                  (values[0] > forceField.energyOverlapCriteria) ? 2.0 * forceField.energyOverlapCriteria : values[0];
            }
            break;
            case ForceFieldSettings::InterpolationScheme::Tricubic:
            {
              std::array<double, 8> values = Interactions::calculateTricubicFractionalAtPositionExternalField(
                  forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

              values[1] /= static_cast<double>(numberOfCells.x);
              values[2] /= static_cast<double>(numberOfCells.y);
              values[3] /= static_cast<double>(numberOfCells.z);
              values[4] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.y);
              values[5] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.z);
              values[6] /= static_cast<double>(numberOfCells.y) * static_cast<double>(numberOfCells.z);
              values[7] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.y) *
                          static_cast<double>(numberOfCells.z);

              if (values[0] > forceField.energyOverlapCriteria)
              {
                // if overlap then make values finite
                data_cell[0, i, j, k] = 2.0 * forceField.energyOverlapCriteria;
                data_cell[1, i, j, k] =
                    std::clamp(values[1], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                data_cell[2, i, j, k] =
                    std::clamp(values[2], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                data_cell[3, i, j, k] =
                    std::clamp(values[3], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                for (std::size_t l = 4; l < 8; ++l)
                {
                  data_cell[l, i, j, k] = 0.0;
                }
              }
              else
              {
                for (std::size_t l = 0; l < 8; ++l)
                {
                  data_cell[l, i, j, k] = values[l];
                }
              }
            }
            break;
            case ForceFieldSettings::InterpolationScheme::Triquintic:
              std::array<double, 27> values = Interactions::calculateTriquinticFractionalAtPositionExternalField(
                  forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

              values[1] /= static_cast<double>(numberOfCells.x);
              values[2] /= static_cast<double>(numberOfCells.y);
              values[3] /= static_cast<double>(numberOfCells.z);

              values[4] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]);
              values[5] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]);
              values[6] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[2]);
              values[7] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]);
              values[8] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
              values[9] /= static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

              values[10] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]);
              values[11] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[2]);
              values[12] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[1]);
              values[13] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]);
              values[14] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]);
              values[15] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[2]) *
                            static_cast<double>(numberOfCells[2]);
              values[16] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                            static_cast<double>(numberOfCells[2]);

              values[17] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]);
              values[18] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);
              values[19] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);
              values[20] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
              values[21] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
              values[22] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

              values[23] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]);
              values[24] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                            static_cast<double>(numberOfCells[2]);
              values[25] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                            static_cast<double>(numberOfCells[2]);

              values[26] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                            static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                            static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

              if (values[0] > forceField.energyOverlapCriteria)
              {
                // if overlap then make values finite
                data_cell[0, i, j, k] = 2.0 * forceField.energyOverlapCriteria;
                data_cell[1, i, j, k] =
                    std::clamp(values[1], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                data_cell[2, i, j, k] =
                    std::clamp(values[2], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                data_cell[3, i, j, k] =
                    std::clamp(values[3], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
                for (std::size_t l = 4; l < 27; ++l)
                {
                  data_cell[l, i, j, k] = 0.0;
                }
              }
              else
              {
                for (std::size_t l = 0; l < 27; ++l)
                {
                  data_cell[l, i, j, k] = values[l];
                }
              }
              break;
          }

          if (static_cast<std::size_t>(counter * percent) > static_cast<std::size_t>((counter - 1.0) * percent))
          {
            std::print(stream, "Percentage finished: {}\n", static_cast<std::size_t>(counter * percent));
            std::flush(stream);
          }
        }
      }
    }
    std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();

    std::chrono::duration<double> time_duration = time_end - time_begin;
    std::print(stream, "Grid done... ({:14f} [s])\n", time_duration.count());
    std::print(stream, "\n");
  }
}

// The grid files are stored row-order (std::layout_right)
// The grid is arranged with the x axis as the outer loop and the z axis as the inner loop
// calculateEnergyAtPosition using ForceFieldSettings::InterpolationGridType::EwaldReal calcalulates energy for a unit-charge
void InterpolationEnergyGrid::makeFrameworkInterpolationGrid(std::ostream& stream,
                                                    ForceFieldSettings::InterpolationGridType interpolationGridType,
                                                    const ForceField& forceField, const Framework& framework,
                                                    double cutOff, std::size_t pseudo_atom_index)
{
  double3 delta = double3(1.0 / static_cast<double>(numberOfCells.x), 1.0 / static_cast<double>(numberOfCells.y),
                          1.0 / static_cast<double>(numberOfCells.z));
  std::mdspan<double, std::dextents<std::size_t, 4>, std::layout_left> data_cell(
      data.data(), std::to_underlying(order), numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  double percent = 100.0 / static_cast<double>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z);

  const SimulationBox unit_cell_box = framework.simulationBox;
  const int3 number_of_unit_cells = unit_cell_box.smallestNumberOfUnitCellsForMinimumImagesConvention(cutOff);
  const std::vector<Atom> framework_atoms = framework.makeSuperCell(number_of_unit_cells);
  const SimulationBox super_cell_box = unit_cell_box.scaled(number_of_unit_cells);

  std::print(stream, "(Using {} {} {} periodic cells to create grid)\n", number_of_unit_cells.x, number_of_unit_cells.y,
             number_of_unit_cells.z);

  double counter{};
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
  for (std::size_t k = 0; k < static_cast<std::size_t>(numberOfGridPoints.z); ++k)
  {
    for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfGridPoints.y); ++j)
    {
      for (std::size_t i = 0; i < static_cast<std::size_t>(numberOfGridPoints.x); ++i)
      {
        counter += 1.0;
        double3 s = double3(static_cast<double>(i) * delta.x, static_cast<double>(j) * delta.y,
                            static_cast<double>(k) * delta.z);
        double3 pos = unit_cell_box.cell * s;

        switch (order)
        {
          case ForceFieldSettings::InterpolationScheme::Polynomial:
          {
            double value = Interactions::calculateEnergyAtPosition(interpolationGridType, forceField, super_cell_box,
                                                                   pos, pseudo_atom_index, framework_atoms);
            data_cell[0, i, j, k] =
                (value > forceField.energyOverlapCriteria) ? 2.0 * forceField.energyOverlapCriteria : value;
          }
          break;
          case ForceFieldSettings::InterpolationScheme::Tricubic:
          {
            std::array<double, 8> values = Interactions::calculateTricubicFractionalAtPosition(
                interpolationGridType, forceField, super_cell_box, pos, pseudo_atom_index, unit_cell_box,
                framework_atoms);

            values[1] /= static_cast<double>(numberOfCells.x);
            values[2] /= static_cast<double>(numberOfCells.y);
            values[3] /= static_cast<double>(numberOfCells.z);
            values[4] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.y);
            values[5] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.z);
            values[6] /= static_cast<double>(numberOfCells.y) * static_cast<double>(numberOfCells.z);
            values[7] /= static_cast<double>(numberOfCells.x) * static_cast<double>(numberOfCells.y) *
                         static_cast<double>(numberOfCells.z);

            if (values[0] > forceField.energyOverlapCriteria)
            {
              // if overlap then make values finite
              data_cell[0, i, j, k] = 2.0 * forceField.energyOverlapCriteria;
              data_cell[1, i, j, k] =
                  std::clamp(values[1], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              data_cell[2, i, j, k] =
                  std::clamp(values[2], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              data_cell[3, i, j, k] =
                  std::clamp(values[3], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              for (std::size_t l = 4; l < 8; ++l)
              {
                data_cell[l, i, j, k] = 0.0;
              }
            }
            else
            {
              for (std::size_t l = 0; l < 8; ++l)
              {
                data_cell[l, i, j, k] = values[l];
              }
            }
          }
          break;
          case ForceFieldSettings::InterpolationScheme::Triquintic:
            std::array<double, 27> values = Interactions::calculateTriquinticFractionalAtPosition(
                interpolationGridType, forceField, super_cell_box, pos, pseudo_atom_index, unit_cell_box,
                framework_atoms);

            values[1] /= static_cast<double>(numberOfCells.x);
            values[2] /= static_cast<double>(numberOfCells.y);
            values[3] /= static_cast<double>(numberOfCells.z);

            values[4] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]);
            values[5] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]);
            values[6] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[2]);
            values[7] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]);
            values[8] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
            values[9] /= static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

            values[10] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]);
            values[11] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[2]);
            values[12] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[1]);
            values[13] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]);
            values[14] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]);
            values[15] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[2]) *
                          static_cast<double>(numberOfCells[2]);
            values[16] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                          static_cast<double>(numberOfCells[2]);

            values[17] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]);
            values[18] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);
            values[19] /= static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);
            values[20] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
            values[21] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]);
            values[22] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

            values[23] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]);
            values[24] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                          static_cast<double>(numberOfCells[2]);
            values[25] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[2]) *
                          static_cast<double>(numberOfCells[2]);

            values[26] /= static_cast<double>(numberOfCells[0]) * static_cast<double>(numberOfCells[0]) *
                          static_cast<double>(numberOfCells[1]) * static_cast<double>(numberOfCells[1]) *
                          static_cast<double>(numberOfCells[2]) * static_cast<double>(numberOfCells[2]);

            if (values[0] > forceField.energyOverlapCriteria)
            {
              // if overlap then make values finite
              data_cell[0, i, j, k] = 2.0 * forceField.energyOverlapCriteria;
              data_cell[1, i, j, k] =
                  std::clamp(values[1], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              data_cell[2, i, j, k] =
                  std::clamp(values[2], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              data_cell[3, i, j, k] =
                  std::clamp(values[3], -forceField.energyOverlapCriteria, forceField.energyOverlapCriteria);
              for (std::size_t l = 4; l < 27; ++l)
              {
                data_cell[l, i, j, k] = 0.0;
              }
            }
            else
            {
              for (std::size_t l = 0; l < 27; ++l)
              {
                data_cell[l, i, j, k] = values[l];
              }
            }
            break;
        }

        if (static_cast<std::size_t>(counter * percent) > static_cast<std::size_t>((counter - 1.0) * percent))
        {
          std::print(stream, "Percentage finished: {}\n", static_cast<std::size_t>(counter * percent));
          std::flush(stream);
        }
      }
    }
  }
  std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> time_duration = time_end - time_begin;
  std::print(stream, "Grid done... ({:14f} [s])\n", time_duration.count());
  std::print(stream, "\n");
}

namespace
{
// The separable polynomial scheme: successive one-dimensional Lagrange interpolations (Numerical Recipes 3.6),
// energy only.
template <std::make_signed_t<std::size_t> NumPoints>
double polynomialInterpolate(const std::vector<double>& data, const uint3& numberOfGridPoints,
                             const uint3& numberOfCells, double3 s)
{
  using Index = std::make_signed_t<std::size_t>;

  double energy, dummy;
  std::array<double, NumPoints> yjtmp, yktmp, yltmp;
  std::array<double, NumPoints> xt, yt, zt;

  // Determine lower boundary
  // Note: the case s==1.0 will be handled by the last cell
  const Index x0 =
      std::min(static_cast<Index>(s.x * static_cast<double>(numberOfCells.x)), static_cast<Index>(numberOfCells.x - 1));
  const Index y0 =
      std::min(static_cast<Index>(s.y * static_cast<double>(numberOfCells.y)), static_cast<Index>(numberOfCells.y - 1));
  const Index z0 =
      std::min(static_cast<Index>(s.z * static_cast<double>(numberOfCells.z)), static_cast<Index>(numberOfCells.z - 1));

  // find the corresponding position within that cell (between 0.0 and 1.0)
  s.x = (s.x * static_cast<double>(numberOfCells.x)) - static_cast<double>(x0);
  s.y = (s.y * static_cast<double>(numberOfCells.y)) - static_cast<double>(y0);
  s.z = (s.z * static_cast<double>(numberOfCells.z)) - static_cast<double>(z0);

  const std::mdspan<const double, std::dextents<std::size_t, 4>, std::layout_left> data_cell(
      data.data(), 1, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  const Index mend = NumPoints / 2;
  const Index mstart = mend + 1 - NumPoints;  // for 4: -1, for 5: -2, for 6: -2

  // The stencil reaches across the cell boundaries. The grid is periodic with a period of numberOfCells
  // planes: the last stored plane (index numberOfCells) is the copy of the first one (s = 0 and s = 1), so the
  // wrap is modulo numberOfCells, not modulo numberOfGridPoints.

  for (Index l0 = 0; l0 < NumPoints; ++l0)
  {
    zt[static_cast<std::size_t>(l0)] = static_cast<double>(l0 + mstart);

    Index lp = z0 + l0 + mstart;
    if (lp < 0) lp += static_cast<Index>(numberOfCells.z);
    if (lp >= static_cast<Index>(numberOfCells.z)) lp -= static_cast<Index>(numberOfCells.z);

    for (Index k0 = 0; k0 < NumPoints; ++k0)
    {
      yt[static_cast<std::size_t>(k0)] = static_cast<double>(k0 + mstart);

      Index kp = y0 + k0 + mstart;
      if (kp < 0) kp += static_cast<Index>(numberOfCells.y);
      if (kp >= static_cast<Index>(numberOfCells.y)) kp -= static_cast<Index>(numberOfCells.y);

      for (Index j0 = 0; j0 < NumPoints; ++j0)
      {
        xt[static_cast<std::size_t>(j0)] = static_cast<double>(j0 + mstart);

        Index jp = x0 + j0 + mstart;
        if (jp < 0) jp += static_cast<Index>(numberOfCells.x);
        if (jp >= static_cast<Index>(numberOfCells.x)) jp -= static_cast<Index>(numberOfCells.x);

        yjtmp[static_cast<std::size_t>(j0)] =
            data_cell[0, static_cast<std::size_t>(jp), static_cast<std::size_t>(kp), static_cast<std::size_t>(lp)];
      }
      Interpolation::polint<NumPoints>(xt, yjtmp, s.x, &yktmp[static_cast<std::size_t>(k0)], &dummy);
    }
    Interpolation::polint<NumPoints>(yt, yktmp, s.y, &yltmp[static_cast<std::size_t>(l0)], &dummy);
  }
  Interpolation::polint<NumPoints>(zt, yltmp, s.z, &energy, &dummy);

  return energy;
}
}  // namespace

InterpolationEnergyGrid::Cell InterpolationEnergyGrid::locateCell(double3 pos) const
{
  const double3 s = (unitCellBox.inverseCell * pos).fract();

  // the scaled coordinate; the case s == 1.0 is handled by the last cell
  const double sx = s.x * static_cast<double>(numberOfCells.x);
  const double sy = s.y * static_cast<double>(numberOfCells.y);
  const double sz = s.z * static_cast<double>(numberOfCells.z);

  Cell cell;
  cell.x0 = std::min(static_cast<std::size_t>(sx), numberOfCells.x - 1);
  cell.y0 = std::min(static_cast<std::size_t>(sy), numberOfCells.y - 1);
  cell.z0 = std::min(static_cast<std::size_t>(sz), numberOfCells.z - 1);

  // the position within that cell (between 0.0 and 1.0)
  cell.t =
      double3(sx - static_cast<double>(cell.x0), sy - static_cast<double>(cell.y0), sz - static_cast<double>(cell.z0));
  return cell;
}

void InterpolationEnergyGrid::prefetch(double3 pos) const
{
  const Cell cell = locateCell(pos);
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Polynomial:
      break;
    case ForceFieldSettings::InterpolationScheme::Tricubic:
      prefetchCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      break;
    case ForceFieldSettings::InterpolationScheme::Triquintic:
      prefetchCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      break;
  }
}

// The gradient and the Hessian of the interpolant with respect to the cell coordinate converted to Cartesian
// space: first to the fractional coordinate (a factor of the number of cells per derivative), then with the
// inverse cell matrix.
double3 InterpolationEnergyGrid::cartesianGradient(const double3& cellGradient) const
{
  const double3 fractionalGradient(cellGradient.x * static_cast<double>(numberOfCells.x),
                                   cellGradient.y * static_cast<double>(numberOfCells.y),
                                   cellGradient.z * static_cast<double>(numberOfCells.z));
  return unitCellBox.inverseCell.transpose() * fractionalGradient;
}

double3x3 InterpolationEnergyGrid::cartesianHessian(const double3x3& cellHessian) const
{
  const double nx = static_cast<double>(numberOfCells.x);
  const double ny = static_cast<double>(numberOfCells.y);
  const double nz = static_cast<double>(numberOfCells.z);
  double3x3 hessian{};
  hessian.ax = cellHessian.ax * nx * nx;
  hessian.ay = cellHessian.ay * nx * ny;
  hessian.az = cellHessian.az * nx * nz;
  hessian.by = cellHessian.by * ny * ny;
  hessian.bz = cellHessian.bz * ny * nz;
  hessian.cz = cellHessian.cz * nz * nz;
  hessian.bx = hessian.ay;
  hessian.cx = hessian.az;
  hessian.cy = hessian.bz;
  return unitCellBox.inverseCell.transpose() * hessian * unitCellBox.inverseCell;
}

void InterpolationEnergyGrid::scaleEnergy(double factor)
{
  for (double& value : data) value *= factor;
}

double InterpolationEnergyGrid::interpolate(double3 pos) const
{
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Polynomial:
      return polynomialInterpolate<num_points_interpolation>(data, numberOfGridPoints, numberOfCells,
                                                             (unitCellBox.inverseCell * pos).fract());
    case ForceFieldSettings::InterpolationScheme::Tricubic:
    {
      const Cell cell = locateCell(pos);
      std::array<double, 64> X;
      gatherCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      return Hermite::evaluate<0>(X, Hermite::tricubicQuantity, Hermite::cubic(cell.t.x), Hermite::cubic(cell.t.y),
                                  Hermite::cubic(cell.t.z))
          .value;
    }
    case ForceFieldSettings::InterpolationScheme::Triquintic:
    {
      const Cell cell = locateCell(pos);
      std::array<double, 216> X;
      gatherCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      return Hermite::evaluate<0>(X, Hermite::triquinticQuantity, Hermite::quintic(cell.t.x),
                                  Hermite::quintic(cell.t.y), Hermite::quintic(cell.t.z))
          .value;
    }
  }

  return 0.0;
}

std::pair<double, double3> InterpolationEnergyGrid::interpolateGradient(double3 pos) const
{
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Polynomial:
      return {polynomialInterpolate<num_points_interpolation>(data, numberOfGridPoints, numberOfCells,
                                                              (unitCellBox.inverseCell * pos).fract()),
              double3(0.0, 0.0, 0.0)};
    case ForceFieldSettings::InterpolationScheme::Tricubic:
    {
      const Cell cell = locateCell(pos);
      std::array<double, 64> X;
      gatherCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      const Hermite::Result r = Hermite::evaluate<1>(X, Hermite::tricubicQuantity, Hermite::cubic(cell.t.x),
                                                     Hermite::cubic(cell.t.y), Hermite::cubic(cell.t.z));
      return {r.value, cartesianGradient(r.gradient)};
    }
    case ForceFieldSettings::InterpolationScheme::Triquintic:
    {
      const Cell cell = locateCell(pos);
      std::array<double, 216> X;
      gatherCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      const Hermite::Result r = Hermite::evaluate<1>(X, Hermite::triquinticQuantity, Hermite::quintic(cell.t.x),
                                                     Hermite::quintic(cell.t.y), Hermite::quintic(cell.t.z));
      return {r.value, cartesianGradient(r.gradient)};
    }
  }

  return {0.0, double3(0.0, 0.0, 0.0)};
}

std::tuple<double, double3, double3x3> InterpolationEnergyGrid::interpolateHessian(double3 pos) const
{
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Polynomial:
      return {polynomialInterpolate<num_points_interpolation>(data, numberOfGridPoints, numberOfCells,
                                                              (unitCellBox.inverseCell * pos).fract()),
              double3(0.0, 0.0, 0.0), double3x3{}};
    case ForceFieldSettings::InterpolationScheme::Tricubic:
    {
      // The tricubic interpolant is C^1 only; its Hessian is not continuous across cells and is not reported.
      const Cell cell = locateCell(pos);
      std::array<double, 64> X;
      gatherCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      const Hermite::Result r = Hermite::evaluate<1>(X, Hermite::tricubicQuantity, Hermite::cubic(cell.t.x),
                                                     Hermite::cubic(cell.t.y), Hermite::cubic(cell.t.z));
      return {r.value, cartesianGradient(r.gradient), double3x3{}};
    }
    case ForceFieldSettings::InterpolationScheme::Triquintic:
    {
      const Cell cell = locateCell(pos);
      std::array<double, 216> X;
      gatherCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      const Hermite::Result r = Hermite::evaluate<2>(X, Hermite::triquinticQuantity, Hermite::quintic(cell.t.x),
                                                     Hermite::quintic(cell.t.y), Hermite::quintic(cell.t.z));
      return {r.value, cartesianGradient(r.gradient), cartesianHessian(r.hessian)};
    }
  }

  return {0.0, double3(0.0, 0.0, 0.0), double3x3{}};
}

std::pair<double, double> InterpolationEnergyGrid::interpolate(const InterpolationEnergyGrid& other, double3 pos) const
{
  // The cell and the basis functions are shared only when the two grids coincide.
  if (other.order != order || !(other.numberOfGridPoints == numberOfGridPoints) ||
      order == ForceFieldSettings::InterpolationScheme::Polynomial)
  {
    other.prefetch(pos);
    return {interpolate(pos), other.interpolate(pos)};
  }

  const Cell cell = locateCell(pos);
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Tricubic:
    {
      prefetchCorners<8>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      std::array<double, 64> X, Y;
      gatherCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      gatherCorners<8>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, Y);
      const Hermite::Basis1D<2> bx = Hermite::cubic(cell.t.x), by = Hermite::cubic(cell.t.y),
                                bz = Hermite::cubic(cell.t.z);
      return {Hermite::evaluate<0>(X, Hermite::tricubicQuantity, bx, by, bz).value,
              Hermite::evaluate<0>(Y, Hermite::tricubicQuantity, bx, by, bz).value};
    }
    case ForceFieldSettings::InterpolationScheme::Triquintic:
    {
      prefetchCorners<27>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      std::array<double, 216> X, Y;
      gatherCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      gatherCorners<27>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, Y);
      const Hermite::Basis1D<3> bx = Hermite::quintic(cell.t.x), by = Hermite::quintic(cell.t.y),
                                bz = Hermite::quintic(cell.t.z);
      return {Hermite::evaluate<0>(X, Hermite::triquinticQuantity, bx, by, bz).value,
              Hermite::evaluate<0>(Y, Hermite::triquinticQuantity, bx, by, bz).value};
    }
    default:
      std::unreachable();
  }
}

std::pair<std::pair<double, double3>, std::pair<double, double3>> InterpolationEnergyGrid::interpolateGradient(
    const InterpolationEnergyGrid& other, double3 pos) const
{
  if (other.order != order || !(other.numberOfGridPoints == numberOfGridPoints) ||
      order == ForceFieldSettings::InterpolationScheme::Polynomial)
  {
    other.prefetch(pos);
    return {interpolateGradient(pos), other.interpolateGradient(pos)};
  }

  const Cell cell = locateCell(pos);
  switch (order)
  {
    case ForceFieldSettings::InterpolationScheme::Tricubic:
    {
      prefetchCorners<8>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      std::array<double, 64> X, Y;
      gatherCorners<8>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      gatherCorners<8>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, Y);
      const Hermite::Basis1D<2> bx = Hermite::cubic(cell.t.x), by = Hermite::cubic(cell.t.y),
                                bz = Hermite::cubic(cell.t.z);
      const Hermite::Result rx = Hermite::evaluate<1>(X, Hermite::tricubicQuantity, bx, by, bz);
      const Hermite::Result ry = Hermite::evaluate<1>(Y, Hermite::tricubicQuantity, bx, by, bz);
      return {{rx.value, cartesianGradient(rx.gradient)}, {ry.value, cartesianGradient(ry.gradient)}};
    }
    case ForceFieldSettings::InterpolationScheme::Triquintic:
    {
      prefetchCorners<27>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0);
      std::array<double, 216> X, Y;
      gatherCorners<27>(data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, X);
      gatherCorners<27>(other.data.data(), numberOfGridPoints, cell.x0, cell.y0, cell.z0, Y);
      const Hermite::Basis1D<3> bx = Hermite::quintic(cell.t.x), by = Hermite::quintic(cell.t.y),
                                bz = Hermite::quintic(cell.t.z);
      const Hermite::Result rx = Hermite::evaluate<1>(X, Hermite::triquinticQuantity, bx, by, bz);
      const Hermite::Result ry = Hermite::evaluate<1>(Y, Hermite::triquinticQuantity, bx, by, bz);
      return {{rx.value, cartesianGradient(rx.gradient)}, {ry.value, cartesianGradient(ry.gradient)}};
    }
    default:
      std::unreachable();
  }
}

void InterpolationEnergyGrid::writeOutput(std::size_t systemId, const SimulationBox &simulationBox,
                                      [[maybe_unused]] const ForceField &forceField)
{
  std::filesystem::create_directory("interpolation_grids");

  // only the energy of each vertex: the Hermite schemes store it together with its derivatives
  std::mdspan<const double, std::dextents<std::size_t, 4>, std::layout_left> data_cell(
      data.data(), std::to_underlying(order), numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  std::ofstream ostream(std::format("interpolation_grids/grid_external_field.s{}.cube", systemId));
  const double3x3 cell = simulationBox.cell;

  std::print(ostream, "Cube interpolation file\n");
  std::print(ostream, "Written by RASPA-3\n");
  std::print(ostream, "{} {} {} {}\n", 1, origin.x, origin.y, origin.z);

  // negative numbers of grid points: the vectors are in Angstrom
  std::print(ostream, "{} {} {} {}\n", -static_cast<std::int64_t>(numberOfGridPoints.x), cell.ax / static_cast<double>(numberOfCells.x), 
      cell.ay / static_cast<double>(numberOfCells.x), cell.az / static_cast<double>(numberOfCells.x));
  std::print(ostream, "{} {} {} {}\n", -static_cast<std::int64_t>(numberOfGridPoints.y), cell.bx / static_cast<double>(numberOfCells.y), 
      cell.by / static_cast<double>(numberOfCells.y), cell.bz / static_cast<double>(numberOfCells.y));
  std::print(ostream, "{} {} {} {}\n", -static_cast<std::int64_t>(numberOfGridPoints.z), cell.cx / static_cast<double>(numberOfCells.z), 
      cell.cy / static_cast<double>(numberOfCells.z), cell.cz / static_cast<double>(numberOfCells.z));

  // work around: needs to be 1 atom for iRASPA
  std::print(ostream, "1 0.0 0.0 0.0 0.0\n");

  // the cube order: x is the slowest and z the fastest index
  for (std::size_t i = 0; i < static_cast<std::size_t>(numberOfGridPoints.x); ++i)
  {
    for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfGridPoints.y); ++j)
    {
      for (std::size_t k = 0; k < static_cast<std::size_t>(numberOfGridPoints.z); ++k)
      {
        std::print(ostream, "{}\n", data_cell[0, i, j, k] * Units::EnergyToKelvin);
      }
    }
  }
}


Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const InterpolationEnergyGrid& s)
{
  archive << s.versionNumber;

  archive << s.unitCellBox;
  archive << s.origin;
  archive << s.numberOfGridPoints;
  archive << s.numberOfCells;
  archive << s.order;
  archive << s.data;

  // the end marker after the fields, where the reader looks for it
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, InterpolationEnergyGrid& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'System' at line {} in file {}\n", location.line(), location.file_name()));
  }

  archive >> s.unitCellBox;
  archive >> s.origin;
  archive >> s.numberOfGridPoints;
  archive >> s.numberOfCells;
  archive >> s.order;
  archive >> s.data;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("System: Error in binary restart\n"));
  }
#endif

  return archive;
}
