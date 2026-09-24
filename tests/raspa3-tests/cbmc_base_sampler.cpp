#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import units;
import forcefield;
import component;
import randomnumbers;
import mc_moves_probabilities;
import cbmc_growth_plan;
import cbmc_operators;

// Distribution tests for the exact CBMC base-conformation sampler at FINITE potential strength: the
// trial base conformation carries no Rosenbluth weight, so its distribution must be the exact bonded
// Boltzmann distribution of the step. Soft potentials (k_bend = 700 K/rad^2 at 300 K) are used so
// that the solid-angle Jacobian visibly shapes the distributions: a wrong power of sin(theta) (the
// historical sin^2 bug) or a missing r^2 bond Jacobian fails these tests immediately.

namespace
{

// Two placed beads and one grown bead: the base sampler must draw the bond length from
// p(r) ~ r^2 exp(-beta u_bond) and the bend angle from p(theta) ~ sin(theta) exp(-beta u_bend).
constexpr std::string_view kThreeBeadJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.17, 1.41, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [700.0, 114.0]]
  ],
  "VanDerWaals" : "auto"
}
)";

// A branch step: two beads grown from the same center in one step. Their sibling-sibling bend
// (2-1-3) couples the two directions; the coupling is carried by the trial's Rosenbluth weight.
constexpr std::string_view kBranchJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH3", [1.54, 0.0, 0.0]],
      ["CH", [0.0, 0.0, 0.0]],
      ["CH3", [-0.51, 1.45, 0.0]],
      ["CH3", [-0.51, -0.72, -1.26]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [1, 3]
  ],
  "Bonds" : [
    [["CH", "CH3"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH", "CH3"], "HARMONIC", [700.0, 114.0]]
  ],
  "VanDerWaals" : "auto"
}
)";

// A quaternary center: THREE beads grown from the same center in one step, coupled by three
// sibling-sibling bends (2-1-3, 2-1-4, 3-1-4) on top of the three anchor bends. The same CH3-CH-CH3
// bend type covers all six triplets. Exercises the sampler's generality beyond a single sibling pair.
constexpr std::string_view kThreeBranchJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH3", [1.54, 0.0, 0.0]],
      ["CH", [0.0, 0.0, 0.0]],
      ["CH3", [-0.51, 1.45, 0.0]],
      ["CH3", [-0.51, -0.72, -1.26]],
      ["CH3", [-0.51, -0.72, 1.26]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [1, 3],
    [1, 4]
  ],
  "Bonds" : [
    [["CH", "CH3"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH", "CH3"], "HARMONIC", [700.0, 114.0]]
  ],
  "VanDerWaals" : "auto"
}
)";

// The same branch molecule with a declared stereocenter (bead 1; the reference geometry above is
// non-planar with signed volume < 0, i.e. S parity). Every sampled branch arrangement must keep it.
constexpr std::string_view kChiralBranchJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH3", [1.54, 0.0, 0.0]],
      ["CH", [0.0, 0.0, 0.0]],
      ["CH3", [-0.51, 1.45, 0.0]],
      ["CH3", [-0.51, -0.72, -1.26]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [1, 3]
  ],
  "Bonds" : [
    [["CH", "CH3"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH", "CH3"], "HARMONIC", [700.0, 114.0]]
  ],
  "VanDerWaals" : "auto",
  "ChiralCenters" : [[1, 0, 2, 3]]
}
)";

constexpr double kTemperature = 300.0;
constexpr double kBendK = 700.0;                                  // [K/rad^2]
constexpr double kBendTheta0 = 114.0 * std::numbers::pi / 180.0;  // [rad]
constexpr double kBondK = 96500.0;                                // [K/A^2]
constexpr double kBondR0 = 1.54;                                  // [A]

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH", false, 13.02, 0.0, 0.0, 6, false}},
                    {{0.0, 3.95}, {0.0, 3.75}, {0.0, 4.68}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                    12.0, true, true, false);
}

double bendEnergy(double theta) { return 0.5 * kBendK * (theta - kBendTheta0) * (theta - kBendTheta0); }
double bondEnergy(double r) { return 0.5 * kBondK * (r - kBondR0) * (r - kBondR0); }

double angleBetween(const double3 &a, const double3 &center, const double3 &c)
{
  double3 va = (a - center).normalized();
  double3 vc = (c - center).normalized();
  return std::acos(std::clamp(double3::dot(va, vc), -1.0, 1.0));
}

// Chi-squared of observed counts against expected bin masses of an unnormalized density, integrated
// per bin on a fine subgrid (midpoint rule) and normalized over the binned range.
double chiSquaredAgainstDensity(const std::vector<std::size_t> &counts, double lo, double hi,
                                const std::function<double(double)> &density)
{
  std::size_t bins = counts.size();
  std::size_t total = std::accumulate(counts.begin(), counts.end(), std::size_t{0});
  std::vector<double> masses(bins, 0.0);
  double norm = 0.0;
  constexpr std::size_t sub = 512;
  for (std::size_t b = 0; b != bins; ++b)
  {
    double a = lo + (hi - lo) * static_cast<double>(b) / static_cast<double>(bins);
    double w = (hi - lo) / static_cast<double>(bins);
    for (std::size_t s = 0; s != sub; ++s)
    {
      masses[b] += density(a + w * (static_cast<double>(s) + 0.5) / static_cast<double>(sub)) * w /
                   static_cast<double>(sub);
    }
    norm += masses[b];
  }
  double chi2 = 0.0;
  for (std::size_t b = 0; b != bins; ++b)
  {
    double expected = static_cast<double>(total) * masses[b] / norm;
    double diff = static_cast<double>(counts[b]) - expected;
    chi2 += diff * diff / std::max(expected, 1e-12);
  }
  return chi2;
}

// Two-sample chi-squared for equal-size samples.
double chiSquaredTwoSample(const std::vector<std::size_t> &a, const std::vector<std::size_t> &b)
{
  double chi2 = 0.0;
  for (std::size_t i = 0; i != a.size(); ++i)
  {
    double n1 = static_cast<double>(a[i]);
    double n2 = static_cast<double>(b[i]);
    if (n1 + n2 < 1.0) continue;
    chi2 += (n1 - n2) * (n1 - n2) / (n1 + n2);
  }
  return chi2;
}

}  // namespace

TEST(CBMC_BASE_SAMPLER, single_bead_bond_and_bend_marginals)
{
  ForceField forceField = makeZeroForceField();
  TemporaryFile file("three-bead-sampler.json", kThreeBeadJson);
  Component component(Component::Type::Adsorbate, 0, forceField, "three-bead-sampler", file.stemPath().string(), 5, 21,
                      MCMoveProbabilities(), std::nullopt, false);

  std::vector<Atom> chainAtoms = component.atoms;
  chainAtoms[0].position = double3(0.0, 0.0, 0.0);
  chainAtoms[1].position = double3(1.54, 0.0, 0.0);

  const std::vector<CBMC::GrowStep> &plan = component.growthPlan({0, 1});
  ASSERT_EQ(plan.size(), 1);
  ASSERT_EQ(plan.front().nextBeads.size(), 1);
  ASSERT_EQ(plan.front().nextBeads.front(), 2);

  RandomNumber random(9241);
  double beta = 1.0 / (Units::KB * kTemperature);

  constexpr std::size_t samples = 200'000;
  constexpr std::size_t bins = 20;
  double sigmaR = std::sqrt(kTemperature / kBondK);
  double rLo = kBondR0 - 5.0 * sigmaR;
  double rHi = kBondR0 + 5.0 * sigmaR;

  std::vector<std::size_t> cosThetaCounts(bins, 0), bondCounts(bins, 0);
  std::size_t outsideR = 0;
  for (std::size_t i = 0; i != samples; ++i)
  {
    std::vector<CBMC::StepTrial> trials =
        CBMC::generateGrowTrials(random, forceField, beta, component, chainAtoms, plan.front(), 1);
    const double3 &grown = trials.front().positions.front().position;

    double r = (grown - chainAtoms[1].position).length();
    double cosTheta = std::cos(angleBetween(chainAtoms[0].position, chainAtoms[1].position, grown));

    std::size_t cosBin = std::min(bins - 1, static_cast<std::size_t>((cosTheta + 1.0) / 2.0 * static_cast<double>(bins)));
    cosThetaCounts[cosBin] += 1;
    if (r < rLo || r >= rHi)
    {
      outsideR += 1;
    }
    else
    {
      bondCounts[static_cast<std::size_t>((r - rLo) / (rHi - rLo) * static_cast<double>(bins))] += 1;
    }
  }

  // p(cos theta) ~ exp(-u(theta)/T): the sin(theta) Jacobian is absorbed by the cos-theta measure.
  double chi2Bend = chiSquaredAgainstDensity(cosThetaCounts, -1.0, 1.0, [&](double c)
                                             { return std::exp(-bendEnergy(std::acos(c)) / kTemperature); });
  double chi2Bond = chiSquaredAgainstDensity(bondCounts, rLo, rHi, [&](double r)
                                             { return r * r * std::exp(-bondEnergy(r) / kTemperature); });

  EXPECT_LT(chi2Bend, 60.0) << "bend-angle distribution is not sin(theta) exp(-beta u)";
  EXPECT_LT(chi2Bond, 60.0) << "bond-length distribution is not r^2 exp(-beta u)";
  EXPECT_LT(outsideR, 10u);
}

TEST(CBMC_BASE_SAMPLER, branch_step_sibling_bend_coupling)
{
  ForceField forceField = makeZeroForceField();
  TemporaryFile file("branch-sampler.json", kBranchJson);
  Component component(Component::Type::Adsorbate, 0, forceField, "branch-sampler", file.stemPath().string(), 5, 21,
                      MCMoveProbabilities(), std::nullopt, false);

  std::vector<Atom> chainAtoms = component.atoms;
  chainAtoms[0].position = double3(1.54, 0.0, 0.0);
  chainAtoms[1].position = double3(0.0, 0.0, 0.0);

  const std::vector<CBMC::GrowStep> &plan = component.growthPlan({0, 1});
  ASSERT_EQ(plan.size(), 1);
  ASSERT_EQ(plan.front().nextBeads.size(), 2);

  RandomNumber random(70311);
  double beta = 1.0 / (Units::KB * kTemperature);

  constexpr std::size_t samples = 150'000;
  constexpr std::size_t bins = 20;

  // Sampler histograms: sibling angle (2-1-3) and one anchor angle (0-1-2), as cosines. The sibling
  // coupling is imposed on the base by rejection, so the raw positions follow the fully coupled
  // bonded Boltzmann density; the trial weight must be exactly one here (no torsions, no
  // spin-variant bends -- the coupling's normalization is a base property handled by
  // logBaseSamplerNormalization, not a Rosenbluth weight).
  std::vector<std::size_t> siblingCounts(bins, 0), anchorCounts(bins, 0);
  auto binOf = [&](double cosValue)
  { return std::min(bins - 1, static_cast<std::size_t>((cosValue + 1.0) / 2.0 * static_cast<double>(bins))); };
  for (std::size_t i = 0; i != samples; ++i)
  {
    std::vector<CBMC::StepTrial> trials =
        CBMC::generateGrowTrials(random, forceField, beta, component, chainAtoms, plan.front(), 1);
    const double3 &grownA = trials.front().positions[0].position;
    const double3 &grownB = trials.front().positions[1].position;
    ASSERT_NEAR(trials.front().torsionWeight, 1.0, 1e-12);
    siblingCounts[binOf(std::cos(angleBetween(grownA, chainAtoms[1].position, grownB)))] += 1;
    anchorCounts[binOf(std::cos(angleBetween(chainAtoms[0].position, chainAtoms[1].position, grownA)))] += 1;
  }

  // Reference: brute-force exact rejection sampling of the joint angular density
  //   p(theta_a, theta_b, dphi) ~ sin(theta_a) sin(theta_b) exp(-beta (u_a + u_b + u_ab)),
  // drawn in the cosine measure with acceptance exp(-u/T) (valid: harmonic bends are >= 0).
  std::vector<std::size_t> siblingReference(bins, 0), anchorReference(bins, 0);
  for (std::size_t i = 0; i != samples;)
  {
    double cosA = 2.0 * random.uniform() - 1.0;
    if (random.uniform() > std::exp(-bendEnergy(std::acos(cosA)) / kTemperature)) continue;
    double cosB = 2.0 * random.uniform() - 1.0;
    if (random.uniform() > std::exp(-bendEnergy(std::acos(cosB)) / kTemperature)) continue;
    double deltaPhi = 2.0 * std::numbers::pi * random.uniform();
    double sinA = std::sqrt(std::max(0.0, 1.0 - cosA * cosA));
    double sinB = std::sqrt(std::max(0.0, 1.0 - cosB * cosB));
    double cosSibling = cosA * cosB + sinA * sinB * std::cos(deltaPhi);
    if (random.uniform() > std::exp(-bendEnergy(std::acos(std::clamp(cosSibling, -1.0, 1.0))) / kTemperature)) continue;
    siblingReference[binOf(cosSibling)] += 1;
    anchorReference[binOf(cosA)] += 1;
    ++i;
  }

  EXPECT_LT(chiSquaredTwoSample(siblingCounts, siblingReference), 60.0)
      << "sibling-sibling bend is not Boltzmann-distributed";
  EXPECT_LT(chiSquaredTwoSample(anchorCounts, anchorReference), 60.0)
      << "anchor bend marginal distorted by the sibling coupling";
}

TEST(CBMC_BASE_SAMPLER, three_branch_step_three_sibling_bends)
{
  ForceField forceField = makeZeroForceField();
  TemporaryFile file("three-branch-sampler.json", kThreeBranchJson);
  Component component(Component::Type::Adsorbate, 0, forceField, "three-branch-sampler", file.stemPath().string(), 5,
                      21, MCMoveProbabilities(), std::nullopt, false);

  std::vector<Atom> chainAtoms = component.atoms;
  chainAtoms[0].position = double3(1.54, 0.0, 0.0);
  chainAtoms[1].position = double3(0.0, 0.0, 0.0);

  const std::vector<CBMC::GrowStep> &plan = component.growthPlan({0, 1});
  ASSERT_EQ(plan.size(), 1);
  ASSERT_EQ(plan.front().nextBeads.size(), 3);

  RandomNumber random(41927);
  double beta = 1.0 / (Units::KB * kTemperature);

  constexpr std::size_t samples = 150'000;
  constexpr std::size_t bins = 20;

  // Sampler histograms: one sibling angle (2-1-3) and one anchor angle (0-1-2), as cosines. All
  // three sibling bends are imposed on the base by rejection; the weight must still be exactly one.
  std::vector<std::size_t> siblingCounts(bins, 0), anchorCounts(bins, 0);
  auto binOf = [&](double cosValue)
  { return std::min(bins - 1, static_cast<std::size_t>((cosValue + 1.0) / 2.0 * static_cast<double>(bins))); };
  for (std::size_t i = 0; i != samples; ++i)
  {
    std::vector<CBMC::StepTrial> trials =
        CBMC::generateGrowTrials(random, forceField, beta, component, chainAtoms, plan.front(), 1);
    const double3 &grownA = trials.front().positions[0].position;
    const double3 &grownB = trials.front().positions[1].position;
    ASSERT_NEAR(trials.front().torsionWeight, 1.0, 1e-12);
    siblingCounts[binOf(std::cos(angleBetween(grownA, chainAtoms[1].position, grownB)))] += 1;
    anchorCounts[binOf(std::cos(angleBetween(chainAtoms[0].position, chainAtoms[1].position, grownA)))] += 1;
  }

  // Reference: brute-force exact rejection sampling of the joint density of three unit vectors,
  //   p ~ [prod_i sin(theta_i) exp(-beta u_anchor(theta_i))] exp(-beta sum_{i<j} u_sib(theta_ij)),
  // drawn in the cosine measure with acceptance exp(-u/T) (valid: harmonic bends are >= 0).
  std::vector<std::size_t> siblingReference(bins, 0), anchorReference(bins, 0);
  for (std::size_t i = 0; i != samples;)
  {
    std::array<double3, 3> v{};
    for (std::size_t k = 0; k != 3; ++k)
    {
      double cosT;
      do cosT = 2.0 * random.uniform() - 1.0;
      while (random.uniform() > std::exp(-bendEnergy(std::acos(cosT)) / kTemperature));
      double phi = 2.0 * std::numbers::pi * random.uniform();
      double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
      v[k] = double3(sinT * std::cos(phi), sinT * std::sin(phi), cosT);
    }
    double couplingEnergy = 0.0;
    for (std::size_t a = 0; a != 3; ++a)
    {
      for (std::size_t b = a + 1; b != 3; ++b)
      {
        couplingEnergy += bendEnergy(std::acos(std::clamp(double3::dot(v[a], v[b]), -1.0, 1.0)));
      }
    }
    if (random.uniform() > std::exp(-couplingEnergy / kTemperature)) continue;
    siblingReference[binOf(std::clamp(double3::dot(v[0], v[1]), -1.0, 1.0))] += 1;
    anchorReference[binOf(v[0].z)] += 1;  // the anchor axis is z in this frame
    ++i;
  }

  EXPECT_LT(chiSquaredTwoSample(siblingCounts, siblingReference), 60.0)
      << "sibling-sibling bend (3 coupled siblings) is not Boltzmann-distributed";
  EXPECT_LT(chiSquaredTwoSample(anchorCounts, anchorReference), 60.0)
      << "anchor bend marginal distorted by the three-sibling coupling";
}

// The base-sampler normalization must equal the brute-force integral of the step's base density:
// prod_i [ int r^2 e^{-beta u_bond} dr x 2 pi int sin(theta) e^{-beta u_anchor} dtheta ] x
// < e^{-beta sum u_sib} > over the independent per-bead draws. Checked for one sibling bend (two
// branches) and three sibling bends (three branches).
TEST(CBMC_BASE_SAMPLER, base_normalization_matches_brute_force)
{
  ForceField forceField = makeZeroForceField();
  double beta = 1.0 / (Units::KB * kTemperature);
  RandomNumber random(88711);

  // Brute-force one-dimensional factors (fine midpoint grids, generous ranges).
  double bondIntegral = 0.0;
  {
    constexpr std::size_t n = 200'000;
    double lo = kBondR0 - 10.0 * std::sqrt(kTemperature / kBondK), hi = kBondR0 + 10.0 * std::sqrt(kTemperature / kBondK);
    for (std::size_t i = 0; i != n; ++i)
    {
      double r = lo + (hi - lo) * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
      bondIntegral += r * r * std::exp(-bondEnergy(r) / kTemperature) * (hi - lo) / static_cast<double>(n);
    }
  }
  double coneIntegral = 0.0;
  {
    constexpr std::size_t n = 200'000;
    for (std::size_t i = 0; i != n; ++i)
    {
      double t = std::numbers::pi * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
      coneIntegral += 2.0 * std::numbers::pi * std::sin(t) * std::exp(-bendEnergy(t) / kTemperature) *
                      std::numbers::pi / static_cast<double>(n);
    }
  }

  // Mean sibling Boltzmann factor over the independent anchor-cone draws, for m grown beads.
  auto meanSiblingBoltzmann = [&](std::size_t m)
  {
    constexpr std::size_t n = 2'000'000;
    double sum = 0.0;
    std::vector<double3> v(m);
    for (std::size_t s = 0; s != n; ++s)
    {
      for (std::size_t k = 0; k != m; ++k)
      {
        double cosT;
        do cosT = 2.0 * random.uniform() - 1.0;
        while (random.uniform() > std::exp(-bendEnergy(std::acos(cosT)) / kTemperature));
        double phi = 2.0 * std::numbers::pi * random.uniform();
        double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
        v[k] = double3(sinT * std::cos(phi), sinT * std::sin(phi), cosT);
      }
      double u = 0.0;
      for (std::size_t a = 0; a != m; ++a)
        for (std::size_t b = a + 1; b != m; ++b)
          u += bendEnergy(std::acos(std::clamp(double3::dot(v[a], v[b]), -1.0, 1.0)));
      sum += std::exp(-u / kTemperature);
    }
    return sum / static_cast<double>(n);
  };

  struct Case
  {
    std::string_view json;
    std::string name;
    std::size_t branches;
  };
  for (const Case &c : {Case{kBranchJson, "branch-norm", 2}, Case{kThreeBranchJson, "three-branch-norm", 3}})
  {
    TemporaryFile file(c.name + ".json", c.json);
    Component component(Component::Type::Adsorbate, 0, forceField, c.name, file.stemPath().string(), 5, 21,
                        MCMoveProbabilities(), std::nullopt, false);
    const std::vector<CBMC::GrowStep> &plan = component.growthPlan({0, 1});
    ASSERT_EQ(plan.front().nextBeads.size(), c.branches);

    double expected = static_cast<double>(c.branches) * (std::log(bondIntegral) + std::log(coneIntegral)) +
                      std::log(meanSiblingBoltzmann(c.branches));
    double actual = CBMC::logBaseSamplerNormalization(beta, component, plan);
    EXPECT_NEAR(actual, expected, 0.02) << c.branches << " branches: base normalization mismatch";
  }
}

TEST(CBMC_BASE_SAMPLER, declared_chirality_preserved)
{
  ForceField forceField = makeZeroForceField();
  TemporaryFile file("chiral-branch-sampler.json", kChiralBranchJson);
  Component component(Component::Type::Adsorbate, 0, forceField, "chiral-branch-sampler", file.stemPath().string(), 5,
                      21, MCMoveProbabilities(), std::nullopt, false);

  std::vector<Atom> chainAtoms = component.atoms;
  chainAtoms[0].position = double3(1.54, 0.0, 0.0);
  chainAtoms[1].position = double3(0.0, 0.0, 0.0);

  const std::vector<CBMC::GrowStep> &plan = component.growthPlan({0, 1});
  ASSERT_EQ(plan.size(), 1);
  ASSERT_EQ(plan.front().nextBeads.size(), 2);

  RandomNumber random(551);
  double beta = 1.0 / (Units::KB * kTemperature);

  // The reference geometry declares S parity (signed volume < 0) for center 1 with neighbors 0,2,3.
  constexpr std::size_t samples = 20'000;
  for (std::size_t i = 0; i != samples; ++i)
  {
    std::vector<CBMC::StepTrial> trials =
        CBMC::generateGrowTrials(random, forceField, beta, component, chainAtoms, plan.front(), 1);
    double3 d1 = chainAtoms[0].position - chainAtoms[1].position;
    double3 d2 = trials.front().positions[0].position - chainAtoms[1].position;
    double3 d3 = trials.front().positions[1].position - chainAtoms[1].position;
    double signedVolume = double3::dot(d1, double3::cross(d2, d3));
    ASSERT_LT(signedVolume, 0.0) << "sample " << i << " flipped the declared stereocenter";
  }
}
