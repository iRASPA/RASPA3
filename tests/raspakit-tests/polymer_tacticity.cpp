#include <gtest/gtest.h>

import std;

import double3;
import polymer_links;

// The normalised chiral volume of a stereocentre (see ChiralVolume): its gradient against finite differences, its
// range, and the sign rule that classifies the dyads of a vinyl chain (meso: equal signs of the chain-ordered
// volumes of consecutive centres; racemo: opposite signs).

namespace
{
double volumeOf(const std::array<double3, 4> &r)
{
  return ChiralVolume::compute(r[1] - r[0], r[2] - r[0], r[3] - r[0]).value;
}
}  // namespace

TEST(Polymer_tacticity, chiral_volume_gradient_matches_finite_differences)
{
  std::mt19937 generator(12345);
  std::uniform_real_distribution<double> uniform(-2.0, 2.0);
  const double delta = 1e-6;
  const double tolerance = 1e-7;

  for (std::size_t trial = 0; trial < 50; ++trial)
  {
    std::array<double3, 4> r{};
    for (double3 &position : r) position = double3(uniform(generator), uniform(generator), uniform(generator));
    const ChiralVolume volume = ChiralVolume::compute(r[1] - r[0], r[2] - r[0], r[3] - r[0]);
    EXPECT_LE(std::abs(volume.value), 1.0);

    // the gradients sum to zero (the volume is translation invariant)
    const double3 total = volume.gradient[0] + volume.gradient[1] + volume.gradient[2] + volume.gradient[3];
    EXPECT_NEAR(total.x, 0.0, 1e-12);
    EXPECT_NEAR(total.y, 0.0, 1e-12);
    EXPECT_NEAR(total.z, 0.0, 1e-12);

    for (std::size_t i = 0; i < 4; ++i)
    {
      for (std::size_t k = 0; k < 3; ++k)
      {
        std::array<double3, 4> plus = r;
        std::array<double3, 4> minus = r;
        plus[i][k] += delta;
        minus[i][k] -= delta;
        const double numerical = (volumeOf(plus) - volumeOf(minus)) / (2.0 * delta);
        EXPECT_NEAR(volume.gradient[i][k], numerical, tolerance);
      }
    }
  }
}

TEST(Polymer_tacticity, chiral_volume_reference_values)
{
  // a tetrahedral centre: |V| = sqrt(1 - 3 cos^2 t + 2 cos^3 t) with cos t = -1/3, i.e. 0.7698
  const double c = 1.0 / std::sqrt(3.0);
  const double3 centre(0.0, 0.0, 0.0);
  const double3 p(c, c, c), n(-c, -c, c), z(-c, c, -c);
  const ChiralVolume tetrahedral = ChiralVolume::compute(p - centre, n - centre, z - centre);
  EXPECT_NEAR(std::abs(tetrahedral.value), 0.7698, 1e-3);
  // mirroring one substituent through the plane of the other two flips the sign
  const ChiralVolume mirrored = ChiralVolume::compute(p - centre, n - centre, double3(c, -c, c) - centre);
  EXPECT_NEAR(mirrored.value, -tetrahedral.value, 1e-12);
  // a planar centre has no chiral volume, and the normalisation makes the value independent of the bond lengths
  const ChiralVolume planar =
      ChiralVolume::compute(double3(1.0, 0.0, 0.0), double3(-0.5, 0.8, 0.0), double3(-0.5, -0.8, 0.0));
  EXPECT_NEAR(planar.value, 0.0, 1e-12);
  const ChiralVolume scaled = ChiralVolume::compute(1.54 * (p - centre), 1.33 * (n - centre), 1.52 * (z - centre));
  EXPECT_NEAR(scaled.value, tetrahedral.value, 1e-12);
}

TEST(Polymer_tacticity, dyad_sign_rule_on_a_planar_zigzag_chain)
{
  // backbone in the xy-plane: the CH centres are the odd atoms (all upper vertices of the zigzag), the CH2 the even
  // ones; a side group pointing to +z or -z; chain-ordered volume with p the previous and n the next backbone atom
  auto backbone = [](std::size_t i) { return double3(1.25 * static_cast<double>(i), (i % 2 == 1) ? 0.9 : 0.0, 0.0); };
  auto sideGroup = [&](std::size_t i, bool up) { return backbone(i) + double3(0.0, 0.6, up ? 1.3 : -1.3); };
  auto sign = [&](std::size_t i, bool up)
  {
    const double3 c = backbone(i);
    return ChiralVolume::compute(backbone(i - 1) - c, backbone(i + 1) - c, sideGroup(i, up) - c).value > 0.0 ? 1 : -1;
  };
  // isotactic (all side groups on the same side): consecutive centres have the same sign (meso dyads)
  EXPECT_EQ(sign(1, true), sign(3, true));
  EXPECT_EQ(sign(3, true), sign(5, true));
  EXPECT_EQ(sign(1, false), sign(3, false));
  // syndiotactic (alternating sides): consecutive centres have opposite signs (racemo dyads)
  EXPECT_EQ(sign(1, true), -sign(3, false));
  EXPECT_EQ(sign(3, false), -sign(5, true));
}
