module;

module surface_curvature;

import std;

import uint3;
import double3;
import double3x3;

namespace
{

// A 3 x 3 system by elimination with partial pivoting. The matrix is the normal-equations matrix of the fit
// below, so it is symmetric and positive semi-definite, and the only way it is singular is that the triangle
// gave fewer than three independent equations. The pivots are tested against the size of the matrix rather
// than against an absolute number, the entries scaling as the square of an edge length and so with the grid.
bool solveThreeByThree(const double matrix[3][3], const double rightHandSide[3], double solution[3])
{
  double scale = 0.0;
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      scale = std::max(scale, std::abs(matrix[row][column]));
    }
  }
  if (!(scale > 0.0)) return false;

  double augmented[3][4];
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column) augmented[row][column] = matrix[row][column];
    augmented[row][3] = rightHandSide[row];
  }

  for (std::size_t column = 0; column < 3; ++column)
  {
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < 3; ++row)
    {
      if (std::abs(augmented[row][column]) > std::abs(augmented[pivot][column])) pivot = row;
    }
    if (std::abs(augmented[pivot][column]) <= 1.0e-12 * scale) return false;
    if (pivot != column)
    {
      for (std::size_t entry = 0; entry < 4; ++entry) std::swap(augmented[column][entry], augmented[pivot][entry]);
    }

    for (std::size_t row = column + 1; row < 3; ++row)
    {
      double factor = augmented[row][column] / augmented[column][column];
      for (std::size_t entry = column; entry < 4; ++entry) augmented[row][entry] -= factor * augmented[column][entry];
    }
  }

  for (std::size_t step = 0; step < 3; ++step)
  {
    std::size_t row = 2 - step;
    double sum = augmented[row][3];
    for (std::size_t column = row + 1; column < 3; ++column) sum -= augmented[row][column] * solution[column];
    solution[row] = sum / augmented[row][row];
  }

  return std::isfinite(solution[0]) && std::isfinite(solution[1]) && std::isfinite(solution[2]);
}

}  // namespace

CurvatureKind classifyCurvature(const TriangleCurvature &curvature, const CurvatureBands &bands)
{
  if (!curvature.resolved) return CurvatureKind::Unresolved;
  if (!std::isfinite(curvature.kappa1) || !std::isfinite(curvature.kappa2)) return CurvatureKind::Unresolved;
  if (bands.sharpest > 0.0 && curvature.sharpest() > bands.sharpest) return CurvatureKind::Unresolved;

  if (curvature.sharpest() <= bands.flat) return CurvatureKind::Flat;

  // kappa1 is the larger of the two, so the only way the pair straddles zero is the one test below.
  if (curvature.kappa1 > bands.flat && curvature.kappa2 < -bands.flat) return CurvatureKind::Saddle;

  // What is left has both curvatures on one side of zero, give or take one of them being too small to have a
  // side, which is a cylinder: a cylindrical pore wall belongs with the concave and a rod with the convex, so
  // the sign is taken from whichever of the two is the larger.
  double dominant = (std::abs(curvature.kappa1) >= std::abs(curvature.kappa2)) ? curvature.kappa1 : curvature.kappa2;
  return (dominant > 0.0) ? CurvatureKind::Convex : CurvatureKind::Concave;
}

void CurvatureAreas::add(double area, const TriangleCurvature &curvature, const CurvatureBands &bands)
{
  if (!std::isfinite(area) || area <= 0.0) return;

  // The two integrals take every triangle whose fit stood up, the band above them or not. They are the part of
  // this that converges as the grid is refined, and they converge precisely because the crease triangles are in
  // them: a crease contributes a curvature of one over the rounding radius over an area proportional to that
  // radius, and the product settles down while neither factor does. Holding those triangles back would leave
  // both integrals drifting with the grid, which is what the columns below do.
  if (curvature.resolved && std::isfinite(curvature.kappa1) && std::isfinite(curvature.kappa2))
  {
    this->integratedMeanCurvature += area * curvature.mean();
    this->integratedGaussianCurvature += area * curvature.gaussian();
  }

  switch (classifyCurvature(curvature, bands))
  {
    case CurvatureKind::Convex:
      this->convex += area;
      ++this->numberOfConvex;
      break;
    case CurvatureKind::Saddle:
      this->saddle += area;
      ++this->numberOfSaddle;
      break;
    case CurvatureKind::Concave:
      this->concave += area;
      ++this->numberOfConcave;
      break;
    case CurvatureKind::Flat:
      this->flat += area;
      ++this->numberOfFlat;
      break;
    case CurvatureKind::Unresolved:
      this->unresolved += area;
      ++this->numberOfUnresolved;
      if (!curvature.resolved) ++this->numberOfDegenerate;
      break;
  }
}

CurvatureBands curvatureBandsForGrid(const double3x3 &unitCell, uint3 gridSize)
{
  double3 alongA = unitCell * double3(1.0 / static_cast<double>(gridSize.x), 0.0, 0.0);
  double3 alongB = unitCell * double3(0.0, 1.0 / static_cast<double>(gridSize.y), 0.0);
  double3 alongC = unitCell * double3(0.0, 0.0, 1.0 / static_cast<double>(gridSize.z));

  double coarsest = std::max({alongA.length(), alongB.length(), alongC.length()});

  CurvatureBands bands;
  bands.sharpest = (coarsest > 0.0) ? 1.0 / coarsest : 0.0;

  return bands;
}


double3 cartesianGradientOfGridGradient(const double3x3 &inverseCell, uint3 gridSize, double3 gradient)
{
  double3 fractional(gradient.x * static_cast<double>(gridSize.x), gradient.y * static_cast<double>(gridSize.y),
                     gradient.z * static_cast<double>(gridSize.z));

  return transposedMultiply(inverseCell, fractional);
}

double3 cartesianNormalOfGridGradient(const double3x3 &inverseCell, uint3 gridSize, double3 gradient,
                                      FieldSense sense)
{
  double3 cartesian = cartesianGradientOfGridGradient(inverseCell, gridSize, gradient);
  if (sense == FieldSense::GrowsIntoSolid) cartesian = -cartesian;

  double length = cartesian.length();
  if (!(length > 0.0) || !std::isfinite(length)) return double3(0.0, 0.0, 0.0);

  return cartesian * (1.0 / length);
}

TriangleCurvature triangleCurvature(const std::array<double3, 3> &corners, const std::array<double3, 3> &normals)
{
  TriangleCurvature result;

  double3 face = double3::cross(corners[1] - corners[0], corners[2] - corners[0]);
  double twiceArea = face.length();
  if (!(twiceArea > 0.0) || !std::isfinite(twiceArea)) return result;
  face = face * (1.0 / twiceArea);

  // The frame is the triangle's own, but signed to agree with the field's normals, so that the sign of a
  // curvature means the same thing on every triangle whichever way marching cubes happened to wind this one.
  double3 average = normals[0] + normals[1] + normals[2];
  if (double3::dot(face, average) < 0.0) face = -face;

  double3 along = corners[1] - corners[0];
  along = along - face * double3::dot(face, along);
  double alongLength = along.length();
  if (!(alongLength > 0.0)) return result;

  double3 u = along * (1.0 / alongLength);
  double3 v = double3::cross(face, u);

  // The normal equations of the six-by-three least-squares fit, accumulated a row at a time. A row reads
  // (row . (a, b, c)) = value, with S = [[a, b], [b, c]] in the frame (u, v).
  double matrix[3][3] = {};
  double rightHandSide[3] = {};

  auto addRow = [&](double first, double second, double third, double value)
  {
    const double row[3] = {first, second, third};
    for (std::size_t a = 0; a < 3; ++a)
    {
      for (std::size_t b = 0; b < 3; ++b) matrix[a][b] += row[a] * row[b];
      rightHandSide[a] += row[a] * value;
    }
  };

  for (std::size_t corner = 0; corner < 3; ++corner)
  {
    std::size_t next = (corner + 1) % 3;

    double3 step = corners[next] - corners[corner];
    double3 turn = normals[next] - normals[corner];

    double stepU = double3::dot(step, u);
    double stepV = double3::dot(step, v);

    addRow(stepU, stepV, 0.0, double3::dot(turn, u));
    addRow(0.0, stepU, stepV, double3::dot(turn, v));
  }

  double shape[3] = {};
  if (!solveThreeByThree(matrix, rightHandSide, shape)) return result;

  // The eigenvalues of a symmetric two-by-two, written so that the square root cannot go negative on a matrix
  // that is symmetric to within rounding.
  double halfTrace = 0.5 * (shape[0] + shape[2]);
  double determinant = shape[0] * shape[2] - shape[1] * shape[1];
  double spread = std::sqrt(std::max(0.0, halfTrace * halfTrace - determinant));

  result.kappa1 = halfTrace + spread;
  result.kappa2 = halfTrace - spread;
  result.resolved = std::isfinite(result.kappa1) && std::isfinite(result.kappa2);

  return result;
}
