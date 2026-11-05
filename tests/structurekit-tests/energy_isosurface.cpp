#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import energy_isosurface;

// The processor extractor's vertex normals. The field is periodic and the extractor repeats the first plane
// at the far end, so the normals on the two planes that meet there have to be central differences read round
// the join: a one-sided difference stopping at the join tilts every normal on a cell face.

namespace
{

constexpr std::size_t n = 32;
const uint3 gridSize{n, n, n};

// cos(2πx) + cos(2πy)/2, which has a maximum on the plane x = 0 and so no gradient along x there. Its level
// set at one crosses that plane, so there are vertices on it to ask about.
std::vector<float> ridge()
{
  std::vector<float> field(n * n * n);
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        double x = static_cast<double>(i) / static_cast<double>(n);
        double y = static_cast<double>(j) / static_cast<double>(n);
        field[(k * n + j) * n + i] = static_cast<float>(std::cos(2.0 * std::numbers::pi * x) +
                                                        0.5 * std::cos(2.0 * std::numbers::pi * y));
      }
    }
  }
  return field;
}

}  // namespace


TEST(energy_isosurface, normals_on_the_cell_face_read_round_the_join)
{
  std::vector<float> field = ridge();
  std::vector<double3> gradients;
  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, 1.0, &gradients);
  ASSERT_EQ(gradients.size(), corners.size());

  std::size_t onTheFace = 0;
  double worstTilt = 0.0;
  for (std::size_t v = 0; v < corners.size(); ++v)
  {
    const double x = corners[v].x;
    if (!(x < 1.0e-9 || x > 1.0 - 1.0e-9)) continue;
    ++onTheFace;
    const double length = gradients[v].length();
    ASSERT_GT(length, 0.0);
    worstTilt = std::max(worstTilt, std::abs(gradients[v].x) / length);
  }

  ASSERT_GT(onTheFace, 0u);
  // The field is even in x about the face, so the exact gradient there has no x at all.
  EXPECT_LT(worstTilt, 1.0e-3);
}


namespace
{

// A field with nothing regular about it: two dozen plane waves of up to five periods across the cell, under an
// envelope that makes it several times steeper on one side of the cell than the other. Its level set at zero is
// folded tightly enough that some cubes are cut the ambiguous ways, which are the ones given a vertex of their own
// inside the cube, and the gradient is far from the same length on the edges of such a cube.
constexpr std::size_t m = 24;

double roughAt(double i, double j, double k)
{
  struct Wave
  {
    int a, b, c;
    double phase, amplitude;
  };
  static const std::vector<Wave> waves = []
  {
    std::vector<Wave> made;
    std::uint64_t state = 0x9E3779B97F4A7C15ull;
    auto next = [&]()
    {
      state = state * 6364136223846793005ull + 1442695040888963407ull;
      return static_cast<double>(state >> 11) / 9007199254740992.0;
    };
    for (int w = 0; w < 24; ++w)
    {
      const int a = static_cast<int>(next() * 11.0) - 5;
      const int b = static_cast<int>(next() * 11.0) - 5;
      const int c = static_cast<int>(next() * 11.0) - 5;
      const double phase = 2.0 * std::numbers::pi * next();
      made.push_back(Wave{a, b, c, phase, 0.2 + next()});
    }
    return made;
  }();

  double value = 0.0;
  for (const Wave &wave : waves)
  {
    value += wave.amplitude *
             std::cos(2.0 * std::numbers::pi * (wave.a * i + wave.b * j + wave.c * k) / static_cast<double>(m) +
                      wave.phase);
  }
  return value * std::exp(1.5 * std::sin(2.0 * std::numbers::pi * i / static_cast<double>(m)));
}

std::vector<float> rough()
{
  std::vector<float> field(m * m * m);
  for (std::size_t k = 0; k < m; ++k)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t i = 0; i < m; ++i)
        field[(k * m + j) * m + i] =
            static_cast<float>(roughAt(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)));
  return field;
}

}  // namespace


// The vertex a cube is given inside itself stands at the mean of the points where the surface crosses its edges, and
// its normal is the field's gradient there: the mean of the gradients at those points, each with its length, since the
// gradient is close to linear across a cube. It is not the mean of the unit normals at them, which is what it was, and
// which lets an edge where the field is nearly flat turn the normal as far as one where it is steep. The two are
// worked out here from the field alone, edge by edge, and are tens of degrees apart on a field like this one.
TEST(energy_isosurface, the_normal_inside_a_cube_is_the_mean_of_the_gradients_on_its_edges)
{
  const uint3 size{m, m, m};
  std::vector<float> field = rough();
  std::vector<double3> gradients;
  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(field, size, 0.0, &gradients);
  ASSERT_EQ(gradients.size(), corners.size());

  const int period = static_cast<int>(m);
  auto valueAt = [&](int i, int j, int k)
  {
    auto wrapped = [&](int index) { return static_cast<std::size_t>(((index % period) + period) % period); };
    double value = static_cast<double>(field[(wrapped(k) * m + wrapped(j)) * m + wrapped(i)]);
    if (std::abs(value) < std::numeric_limits<double>::epsilon()) value = std::numeric_limits<double>::epsilon();
    return value;
  };
  auto gradientAt = [&](int i, int j, int k)
  {
    return double3(0.5 * (valueAt(i + 1, j, k) - valueAt(i - 1, j, k)),
                   0.5 * (valueAt(i, j + 1, k) - valueAt(i, j - 1, k)),
                   0.5 * (valueAt(i, j, k + 1) - valueAt(i, j, k - 1)));
  };

  std::size_t inside = 0;
  std::size_t toldApart = 0;
  for (std::size_t v = 0; v < corners.size(); ++v)
  {
    const double3 at(corners[v].x * static_cast<double>(m), corners[v].y * static_cast<double>(m),
                     corners[v].z * static_cast<double>(m));
    int onGridPlanes = 0;
    for (double coordinate : {at.x, at.y, at.z})
    {
      if (std::abs(coordinate - std::round(coordinate)) < 1.0e-9) ++onGridPlanes;
    }
    if (onGridPlanes >= 2) continue;  // on an edge of the grid, where there is one gradient and nothing to average
    ASSERT_EQ(onGridPlanes, 0);
    ++inside;

    const int i = static_cast<int>(std::floor(at.x));
    const int j = static_cast<int>(std::floor(at.y));
    const int k = static_cast<int>(std::floor(at.z));
    double3 meanPosition(0.0, 0.0, 0.0);
    double3 meanGradient(0.0, 0.0, 0.0);
    double3 meanNormal(0.0, 0.0, 0.0);
    double crossings = 0.0;
    auto edge = [&](int a, int b, int c, int da, int db, int dc)
    {
      const double from = valueAt(a, b, c);
      const double to = valueAt(a + da, b + db, c + dc);
      if ((from < 0.0) == (to < 0.0)) return;
      const double u = from / (from - to);
      const double3 gradient = (1.0 - u) * gradientAt(a, b, c) + u * gradientAt(a + da, b + db, c + dc);
      meanPosition += double3(a + da * u, b + db * u, c + dc * u);
      meanGradient += gradient;
      meanNormal += gradient * (1.0 / gradient.length());
      crossings += 1.0;
    };
    for (int t : {0, 1})
    {
      for (int s : {0, 1})
      {
        edge(i, j + s, k + t, 1, 0, 0);
        edge(i + t, j, k + s, 0, 1, 0);
        edge(i + s, j + t, k, 0, 0, 1);
      }
    }
    ASSERT_GT(crossings, 0.0);
    meanPosition = meanPosition * (1.0 / crossings);
    ASSERT_LT((meanPosition - at).length(), 1.0e-9)
        << "the vertex is not the one inside the cube it was looked for in";

    const double3 expected = meanGradient * (1.0 / meanGradient.length());
    const double3 ofUnitNormals = meanNormal * (1.0 / meanNormal.length());
    const double3 given = gradients[v] * (1.0 / gradients[v].length());
    EXPECT_GT(double3::dot(given, expected), 1.0 - 1.0e-12)
        << "inside the cube at (" << i << ", " << j << ", " << k << ")";
    if (double3::dot(ofUnitNormals, expected) < std::cos(0.05)) ++toldApart;
  }

  EXPECT_GT(inside, 50uz) << "the surface has to have such vertices for this to have asked anything";
  EXPECT_GT(toldApart, inside / 2) << "and the two means have to differ at them";
}
