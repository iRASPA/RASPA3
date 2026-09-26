#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import molecule;
import units;
import forcefield;
import component;
import framework;
import simulationbox;
import interpolation_energy_grid;
import randomnumbers;
import mc_moves_probabilities;
import cbmc;
import cbmc_growth_plan;

// Fixed-endpoint (bridging) CBMC regrowth: a 'Partial-reinsertion' fixed set that lies on both sides
// of the regrown part. The interior segment is grown from one fixed side and closed onto the other
// (GrowStep::Kind::CloseBridge, with closure guides on the steps before it). The weights are exact,
// so a Metropolis chain of regrow/retrace moves accepted with W_new/W_old must reproduce the bonded
// Boltzmann distribution of the interior beads; the reference is an independent single-bead
// displacement Metropolis chain on the same energy function.

namespace
{

constexpr double kTemperature = 300.0;
constexpr double kBondK = 96500.0;  // [K/A^2]
constexpr double kBondR0 = 1.54;    // [A]

// Five beads, fixed {0, 1, 3, 4}: a single bead closes the bridge (one CloseBridge step, no guides).
constexpr std::string_view kPentaneJson =
R"({
  "CriticalTemperature" : 469.7,
  "CriticalPressure" : 3370000.0,
  "AcentricFactor" : 0.251,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.17, 1.41, 0.0]],
      ["CH2", [3.71, 1.41, 0.0]],
      ["CH2", [4.34, 2.82, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [700.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto",
  "Partial-reinsertion" : [
    [0, 1, 3, 4]
  ]
}
)";

// Six beads, fixed {0, 1, 4, 5}: bead 2 is attached with a closure guide toward bead 4, bead 3
// closes the bridge.
constexpr std::string_view kHexaneJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.17, 1.41, 0.0]],
      ["CH2", [3.71, 1.41, 0.0]],
      ["CH2", [4.34, 2.82, 0.0]],
      ["CH2", [5.88, 2.82, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [700.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto",
  "Partial-reinsertion" : [
    [0, 1, 4, 5]
  ]
}
)";

// The same six-bead chain with holonomic (FIXED) bonds.
constexpr std::string_view kHexaneFixedBondsJson =
R"({
  "CriticalTemperature" : 507.6,
  "CriticalPressure" : 3025000.0,
  "AcentricFactor" : 0.301,
  "pseudoAtoms" :
    [
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.17, 1.41, 0.0]],
      ["CH2", [3.71, 1.41, 0.0]],
      ["CH2", [4.34, 2.82, 0.0]],
      ["CH2", [5.88, 2.82, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "Bonds" : [
    [["CH2", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH2", "CH2", "CH2"], "HARMONIC", [62500.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto",
  "Partial-reinsertion" : [
    [0, 1, 4, 5],
    [0, 4, 5]
  ]
}
)";

// A branch: fixing {0, 2, 3} would require bead 1 to close onto three placed beads at once.
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

ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH", false, 13.02, 0.0, 0.0, 6, false}},
                    {{0.0, 3.95}, {0.0, 3.75}, {0.0, 4.68}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                    12.0, true, true, false);
}

struct Fixture
{
  ForceField forceField{makeZeroForceField()};
  TemporaryFile file;
  Component component;
  SimulationBox box{30.0, 30.0, 30.0};
  double beta{1.0 / (Units::KB * kTemperature)};
  std::optional<Framework> noFramework{};
  std::vector<std::optional<InterpolationEnergyGrid>> noGrids{forceField.pseudoAtoms.size() + 1};
  std::optional<InterpolationEnergyGrid> noExternalFieldGrid{};

  Fixture(std::string name, std::string_view json)
      : file(name + ".json", json),
        component(Component::Type::Adsorbate, 0, forceField, name, file.stemPath().string(), 5, 21,
                  MCMoveProbabilities(), std::nullopt, false)
  {
    component.prepareGrowthPlans(beta);
  }

  CBMC::GrowContext context() const
  {
    return CBMC::GrowContext(false, forceField, box, noGrids, noExternalFieldGrid, noFramework,
                             std::span<const Atom>{}, std::span<const Atom>{}, beta, CBMC::CutOffMode::Full);
  }

  CBMC::GrowResult initialMolecule(RandomNumber &random) const
  {
    const CBMC::GrowContext empty = context();
    std::optional<CBMC::GrowResult> grown;
    while (!grown) grown = CBMC::growNewMolecule(random, empty, component, {.componentId = 0, .moleculeId = 0});
    return *grown;
  }
};

double angleBetween(const double3 &a, const double3 &center, const double3 &c)
{
  double3 va = (a - center).normalized();
  double3 vc = (c - center).normalized();
  return std::acos(std::clamp(double3::dot(va, vc), -1.0, 1.0));
}

double dihedral(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  double3 b1 = b - a, b2 = c - b, b3 = d - c;
  double3 n1 = double3::cross(b1, b2);
  double3 n2 = double3::cross(b2, b3);
  double3 m1 = double3::cross(n1, b2.normalized());
  return std::atan2(double3::dot(m1, n2), double3::dot(n1, n2));
}

// Two-sample chi-squared for equal-size samples.
double chiSquaredTwoSample(const std::vector<double> &a, const std::vector<double> &b)
{
  double chi2 = 0.0;
  for (std::size_t i = 0; i != a.size(); ++i)
  {
    double n1 = a[i], n2 = b[i];
    if (n1 + n2 < 1.0) continue;
    chi2 += (n1 - n2) * (n1 - n2) / (n1 + n2);
  }
  return chi2;
}

struct Histogram
{
  double lo, hi;
  std::vector<double> counts;
  Histogram(double lo, double hi, std::size_t bins) : lo(lo), hi(hi), counts(bins, 0.0) {}
  void add(double x)
  {
    double f = std::clamp((x - lo) / (hi - lo), 0.0, 1.0 - 1e-12);
    counts[static_cast<std::size_t>(f * static_cast<double>(counts.size()))] += 1.0;
  }
};

// Observables of one configuration binned into a set of histograms.
using Observer = std::function<void(std::span<const Atom>, std::vector<Histogram> &)>;

// A Markov chain of fixed-endpoint regrow moves: the fixed beads keep their positions, the others
// are regrown; accept with W_new / W_old. Returns the acceptance fraction.
double runRegrowChain(RandomNumber &random, const Fixture &f, std::vector<Atom> &atoms, const Molecule &molecule,
                      const std::vector<std::size_t> &fixed, std::size_t iterations, std::size_t stride,
                      const Observer &observe, std::vector<Histogram> &histograms)
{
  const CBMC::GrowContext empty = f.context();
  std::size_t accepted = 0;
  for (std::size_t i = 0; i != iterations; ++i)
  {
    std::optional<CBMC::GrowResult> grown = CBMC::regrowMolecule(
        random, empty, f.component, molecule, atoms,
        {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = fixed});
    if (grown)
    {
      CBMC::RetraceResult old = CBMC::retraceMolecule(
          random, empty, f.component, atoms,
          {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = fixed});
      if (random.uniform() < std::exp(grown->logRosenbluthWeight - old.logRosenbluthWeight))
      {
        atoms = grown->atoms;
        ++accepted;
      }
    }
    if (i % stride == stride - 1) observe(atoms, histograms);
  }
  return static_cast<double>(accepted) / static_cast<double>(iterations);
}

// The reference: Metropolis on the full bonded energy of a contiguous movable segment of the chain
// with symmetric, volume-preserving proposals: a random displacement of one bead, a rotation of one
// bead about the axis through its two neighbours (crankshaft), and a rotation of the whole segment
// about the axis through the fixed beads on either side. The rotations cross the torsional barriers
// that displacements alone cross only slowly.
void runDisplacementChain(RandomNumber &random, const Fixture &f, std::vector<Atom> atoms,
                          const std::vector<std::size_t> &movable, std::size_t iterations, std::size_t stride,
                          const Observer &observe, std::vector<Histogram> &histograms)
{
  constexpr double kStep = 0.25;  // [A]
  auto rotateAbout = [&](std::size_t bead, const double3 &from, const double3 &to, double angle)
  {
    double3 axis = (to - from).normalized();
    atoms[bead].position = from + axis.rotateAroundAxis(atoms[bead].position - from, angle);
  };
  std::vector<double3> oldPositions(atoms.size());
  double energy = f.component.intraMolecularPotentials.computeInternalEnergies(atoms).potentialEnergy();
  for (std::size_t i = 0; i != iterations; ++i)
  {
    for (std::size_t k : movable) oldPositions[k] = atoms[k].position;
    std::size_t bead = movable[static_cast<std::size_t>(random.uniform() * static_cast<double>(movable.size())) %
                               movable.size()];
    double choice = random.uniform();
    if (choice < 1.0 / 3.0)
    {
      atoms[bead].position += kStep * double3(2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0,
                                              2.0 * random.uniform() - 1.0);
    }
    else if (choice < 2.0 / 3.0)
    {
      rotateAbout(bead, atoms[bead - 1].position, atoms[bead + 1].position,
                  (2.0 * random.uniform() - 1.0) * std::numbers::pi);
    }
    else
    {
      double angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
      double3 from = atoms[movable.front() - 1].position, to = atoms[movable.back() + 1].position;
      for (std::size_t k : movable) rotateAbout(k, from, to, angle);
    }
    double newEnergy = f.component.intraMolecularPotentials.computeInternalEnergies(atoms).potentialEnergy();
    if (random.uniform() < std::exp(-f.beta * (newEnergy - energy)))
    {
      energy = newEnergy;
    }
    else
    {
      for (std::size_t k : movable) atoms[k].position = oldPositions[k];
    }
    if (i % stride == stride - 1) observe(atoms, histograms);
  }
}

}  // namespace

TEST(CBMC_FIXED_ENDPOINT, plan_closes_the_bridge_and_guides_the_beads_before_it)
{
  Fixture f("fixed-endpoint-hexane-plan", kHexaneJson);

  // Fixed {0, 1, 4, 5}: bead 2 attached from bead 1 (guided toward bead 4), bead 3 closes onto 4.
  const std::vector<CBMC::GrowStep> &plan = f.component.growthPlan({0, 1, 4, 5});
  ASSERT_EQ(plan.size(), 2);

  EXPECT_EQ(plan[0].kind, CBMC::GrowStep::Kind::AttachFragment);
  EXPECT_EQ(plan[0].currentBead, 1);
  ASSERT_EQ(plan[0].nextBeads.size(), 1);
  EXPECT_EQ(plan[0].nextBeads[0], 2);
  ASSERT_EQ(plan[0].spin.guides.size(), 1);
  EXPECT_EQ(plan[0].spin.guides[0].nextBeadIndex, 0);
  EXPECT_EQ(plan[0].spin.guides[0].targetBead, 4);
  EXPECT_EQ(plan[0].spin.guides[0].path.numberOfBonds(), 2);
  ASSERT_TRUE(plan[0].spin.guides[0].table != nullptr) << "guide table not prepared for the plan's beta";
  // The guide prefers reachable distances: two 1.54 A bonds at 114 degrees span about 2.58 A, and
  // a distance beyond the stretched sub-chain is (almost) unreachable.
  EXPECT_GT(plan[0].spin.guides[0].table->logGuideAt(2.58), plan[0].spin.guides[0].table->logGuideAt(0.8));
  EXPECT_GT(plan[0].spin.guides[0].table->logGuideAt(2.58), plan[0].spin.guides[0].table->logGuideAt(3.5));

  // Bead 3 is bonded to two placed beads (2 and 4); either is the anchor, the other the closure.
  EXPECT_EQ(plan[1].kind, CBMC::GrowStep::Kind::CloseBridge);
  ASSERT_EQ(plan[1].nextBeads.size(), 1);
  EXPECT_EQ(plan[1].nextBeads[0], 3);
  ASSERT_TRUE(plan[1].previousBead.has_value());
  ASSERT_TRUE(plan[1].closureBead.has_value());
  EXPECT_TRUE((plan[1].currentBead == 2 && plan[1].closureBead.value() == 4 && plan[1].previousBead.value() == 1) ||
              (plan[1].currentBead == 4 && plan[1].closureBead.value() == 2 && plan[1].previousBead.value() == 5));
  EXPECT_TRUE(plan[1].bridge.anchorBond.has_value());
  EXPECT_TRUE(plan[1].bridge.closureBond.has_value());
  EXPECT_TRUE(plan[1].bridge.midBend.has_value());
  EXPECT_TRUE(plan[1].spin.guides.empty());

  // A single interior bead: one CloseBridge step, nothing to guide.
  const std::vector<CBMC::GrowStep> &single = f.component.growthPlan({0, 1, 3, 4, 5});
  ASSERT_EQ(single.size(), 1);
  EXPECT_EQ(single[0].kind, CBMC::GrowStep::Kind::CloseBridge);
  EXPECT_EQ(single[0].nextBeads[0], 2);
  EXPECT_TRUE((single[0].currentBead == 1 && single[0].closureBead.value() == 3) ||
              (single[0].currentBead == 3 && single[0].closureBead.value() == 1));

  // A fixed chain end without a placed neighbour: the seed step carries the guide toward bead 4.
  const std::vector<CBMC::GrowStep> &fromEnd = f.component.growthPlan({0, 4, 5});
  ASSERT_EQ(fromEnd.size(), 3);
  EXPECT_EQ(fromEnd[0].kind, CBMC::GrowStep::Kind::PlaceSeedFragment);
  EXPECT_EQ(fromEnd[0].nextBeads[0], 1);
  ASSERT_EQ(fromEnd[0].spin.guides.size(), 1);
  EXPECT_EQ(fromEnd[0].spin.guides[0].targetBead, 4);
  EXPECT_EQ(fromEnd[0].spin.guides[0].path.numberOfBonds(), 3);
  EXPECT_EQ(fromEnd[1].kind, CBMC::GrowStep::Kind::AttachFragment);
  ASSERT_EQ(fromEnd[1].spin.guides.size(), 1);
  EXPECT_EQ(fromEnd[1].spin.guides[0].path.numberOfBonds(), 2);
  EXPECT_EQ(fromEnd[2].kind, CBMC::GrowStep::Kind::CloseBridge);

  // Only the two chain ends fixed: the whole interior is bridged. The frontier scan visits the fixed
  // beads in order, so a guided seed grows from each end before the middle is attached and closed.
  const std::vector<CBMC::GrowStep> &interior = f.component.growthPlan({0, 5});
  ASSERT_EQ(interior.size(), 4);
  EXPECT_EQ(interior[0].kind, CBMC::GrowStep::Kind::PlaceSeedFragment);
  EXPECT_EQ(interior[0].nextBeads[0], 1);
  ASSERT_EQ(interior[0].spin.guides.size(), 1);
  EXPECT_EQ(interior[0].spin.guides[0].path.numberOfBonds(), 4);
  EXPECT_EQ(interior[1].kind, CBMC::GrowStep::Kind::PlaceSeedFragment);
  EXPECT_EQ(interior[1].nextBeads[0], 4);
  ASSERT_EQ(interior[1].spin.guides.size(), 1);
  EXPECT_EQ(interior[1].spin.guides[0].targetBead, 1);
  EXPECT_EQ(interior[1].spin.guides[0].path.numberOfBonds(), 3);
  EXPECT_EQ(interior[2].kind, CBMC::GrowStep::Kind::AttachFragment);
  ASSERT_EQ(interior[2].spin.guides.size(), 1);
  EXPECT_EQ(interior[2].spin.guides[0].path.numberOfBonds(), 2);
  EXPECT_EQ(interior[3].kind, CBMC::GrowStep::Kind::CloseBridge);
  EXPECT_EQ(interior[3].nextBeads[0], 3);

  // A connected fixed set is unchanged: plain attach steps, no closure, no guides.
  const std::vector<CBMC::GrowStep> &connected = f.component.growthPlan({0, 1, 2});
  ASSERT_EQ(connected.size(), 3);
  for (const CBMC::GrowStep &step : connected)
  {
    EXPECT_EQ(step.kind, CBMC::GrowStep::Kind::AttachFragment);
    EXPECT_FALSE(step.closureBead.has_value());
    EXPECT_TRUE(step.spin.guides.empty());
  }
}

TEST(CBMC_FIXED_ENDPOINT, over_determined_closure_is_rejected_at_plan_time)
{
  Fixture f("fixed-endpoint-branch", kBranchJson);
  EXPECT_THROW((void)f.component.growthPlan({0, 2, 3}), std::runtime_error);
}

TEST(CBMC_FIXED_ENDPOINT, single_bead_bridge_matches_metropolis_reference)
{
  Fixture f("fixed-endpoint-pentane", kPentaneJson);
  RandomNumber random(7411);
  CBMC::GrowResult initial = f.initialMolecule(random);

  const std::vector<std::size_t> fixed{0, 1, 3, 4};
  const std::vector<std::size_t> movable{2};
  constexpr std::size_t bins = 20;
  const double sigmaR = std::sqrt(kTemperature / kBondK);
  auto makeHistograms = [&]
  {
    return std::vector<Histogram>{Histogram(kBondR0 - 5.0 * sigmaR, kBondR0 + 5.0 * sigmaR, bins),
                                  Histogram(-1.0, 1.0, bins), Histogram(-std::numbers::pi, std::numbers::pi, bins)};
  };
  Observer observe = [](std::span<const Atom> atoms, std::vector<Histogram> &h)
  {
    h[0].add((atoms[2].position - atoms[1].position).length());
    h[1].add(std::cos(angleBetween(atoms[1].position, atoms[2].position, atoms[3].position)));
    h[2].add(dihedral(atoms[0].position, atoms[1].position, atoms[2].position, atoms[3].position));
  };

  constexpr std::size_t samples = 30'000;
  std::vector<Histogram> cbmc = makeHistograms();
  std::vector<Atom> atoms = initial.atoms;
  double acceptance = runRegrowChain(random, f, atoms, initial.molecule, fixed, 4 * samples, 4, observe, cbmc);
  EXPECT_GT(acceptance, 0.2) << "fixed-endpoint regrowth acceptance collapsed";

  std::vector<Histogram> reference = makeHistograms();
  runDisplacementChain(random, f, initial.atoms, movable, 40 * samples, 40, observe, reference);

  // Correlated samples inflate the chi-squared beyond its Poisson scale (19 degrees of freedom).
  EXPECT_LT(chiSquaredTwoSample(cbmc[0].counts, reference[0].counts), 100.0) << "bond length 1-2";
  EXPECT_LT(chiSquaredTwoSample(cbmc[1].counts, reference[1].counts), 100.0) << "bend angle 1-2-3";
  EXPECT_LT(chiSquaredTwoSample(cbmc[2].counts, reference[2].counts), 100.0) << "torsion 0-1-2-3";
}

TEST(CBMC_FIXED_ENDPOINT, guided_two_bead_bridge_matches_metropolis_reference)
{
  Fixture f("fixed-endpoint-hexane", kHexaneJson);
  RandomNumber random(9127);
  CBMC::GrowResult initial = f.initialMolecule(random);

  const std::vector<std::size_t> fixed{0, 1, 4, 5};
  const std::vector<std::size_t> movable{2, 3};
  constexpr std::size_t bins = 20;
  const double sigmaR = std::sqrt(kTemperature / kBondK);
  auto makeHistograms = [&]
  {
    return std::vector<Histogram>{Histogram(kBondR0 - 5.0 * sigmaR, kBondR0 + 5.0 * sigmaR, bins),
                                  Histogram(-1.0, 1.0, bins), Histogram(-1.0, 1.0, bins),
                                  Histogram(-std::numbers::pi, std::numbers::pi, bins),
                                  Histogram(-std::numbers::pi, std::numbers::pi, bins), Histogram(0.0, 4.0, bins),
                                  Histogram(kBondR0 - 5.0 * sigmaR, kBondR0 + 5.0 * sigmaR, bins)};
  };
  Observer observe = [](std::span<const Atom> atoms, std::vector<Histogram> &h)
  {
    h[0].add((atoms[3].position - atoms[2].position).length());
    h[1].add(std::cos(angleBetween(atoms[1].position, atoms[2].position, atoms[3].position)));
    h[2].add(std::cos(angleBetween(atoms[2].position, atoms[3].position, atoms[4].position)));
    h[3].add(dihedral(atoms[0].position, atoms[1].position, atoms[2].position, atoms[3].position));
    h[4].add(dihedral(atoms[1].position, atoms[2].position, atoms[3].position, atoms[4].position));
    h[5].add((atoms[4].position - atoms[2].position).length());  // the closure distance D
    h[6].add((atoms[4].position - atoms[3].position).length());
  };

  constexpr std::size_t samples = 30'000;
  std::vector<Histogram> cbmc = makeHistograms();
  std::vector<Atom> atoms = initial.atoms;
  double acceptance = runRegrowChain(random, f, atoms, initial.molecule, fixed, 4 * samples, 4, observe, cbmc);
  EXPECT_GT(acceptance, 0.1) << "fixed-endpoint regrowth acceptance collapsed";

  std::vector<Histogram> reference = makeHistograms();
  runDisplacementChain(random, f, initial.atoms, movable, 80 * samples, 80, observe, reference);

  EXPECT_LT(chiSquaredTwoSample(cbmc[0].counts, reference[0].counts), 100.0) << "bond length 2-3";
  EXPECT_LT(chiSquaredTwoSample(cbmc[1].counts, reference[1].counts), 100.0) << "bend angle 1-2-3";
  EXPECT_LT(chiSquaredTwoSample(cbmc[2].counts, reference[2].counts), 100.0) << "bend angle 2-3-4";
  EXPECT_LT(chiSquaredTwoSample(cbmc[3].counts, reference[3].counts), 100.0) << "torsion 0-1-2-3";
  EXPECT_LT(chiSquaredTwoSample(cbmc[4].counts, reference[4].counts), 100.0) << "torsion 1-2-3-4";
  EXPECT_LT(chiSquaredTwoSample(cbmc[5].counts, reference[5].counts), 100.0) << "closure distance 2-4";
  EXPECT_LT(chiSquaredTwoSample(cbmc[6].counts, reference[6].counts), 100.0) << "bond length 3-4";
}

TEST(CBMC_FIXED_ENDPOINT, guided_seed_from_fixed_chain_end_matches_metropolis_reference)
{
  // Fixed {0, 4, 5}: bead 1 is a guided seed (bead 0 has no placed neighbour), bead 2 a guided attach,
  // bead 3 closes onto 4. Exercises the three-bond (sampled) guide table and the guided seed operator.
  Fixture f("fixed-endpoint-hexane-end", kHexaneJson);
  RandomNumber random(5519);
  CBMC::GrowResult initial = f.initialMolecule(random);

  const std::vector<std::size_t> fixed{0, 4, 5};
  const std::vector<std::size_t> movable{1, 2, 3};
  constexpr std::size_t bins = 20;
  const double sigmaR = std::sqrt(kTemperature / kBondK);
  auto makeHistograms = [&]
  {
    return std::vector<Histogram>{Histogram(kBondR0 - 5.0 * sigmaR, kBondR0 + 5.0 * sigmaR, bins),
                                  Histogram(-1.0, 1.0, bins), Histogram(-1.0, 1.0, bins),
                                  Histogram(-std::numbers::pi, std::numbers::pi, bins), Histogram(0.0, 4.0, bins),
                                  Histogram(0.0, 5.0, bins)};
  };
  Observer observe = [](std::span<const Atom> atoms, std::vector<Histogram> &h)
  {
    h[0].add((atoms[1].position - atoms[0].position).length());
    h[1].add(std::cos(angleBetween(atoms[0].position, atoms[1].position, atoms[2].position)));
    h[2].add(std::cos(angleBetween(atoms[2].position, atoms[3].position, atoms[4].position)));
    h[3].add(dihedral(atoms[0].position, atoms[1].position, atoms[2].position, atoms[3].position));
    h[4].add((atoms[4].position - atoms[2].position).length());  // the closure distance D
    h[5].add((atoms[4].position - atoms[1].position).length());  // the seed's guide distance
  };

  constexpr std::size_t samples = 20'000;
  std::vector<Histogram> cbmc = makeHistograms();
  std::vector<Atom> atoms = initial.atoms;
  double acceptance = runRegrowChain(random, f, atoms, initial.molecule, fixed, 4 * samples, 4, observe, cbmc);
  EXPECT_GT(acceptance, 0.05) << "fixed-endpoint regrowth acceptance collapsed";

  std::vector<Histogram> reference = makeHistograms();
  runDisplacementChain(random, f, initial.atoms, movable, 100 * samples, 100, observe, reference);

  EXPECT_LT(chiSquaredTwoSample(cbmc[0].counts, reference[0].counts), 100.0) << "bond length 0-1";
  EXPECT_LT(chiSquaredTwoSample(cbmc[1].counts, reference[1].counts), 100.0) << "bend angle 0-1-2";
  EXPECT_LT(chiSquaredTwoSample(cbmc[2].counts, reference[2].counts), 100.0) << "bend angle 2-3-4";
  EXPECT_LT(chiSquaredTwoSample(cbmc[3].counts, reference[3].counts), 100.0) << "torsion 0-1-2-3";
  EXPECT_LT(chiSquaredTwoSample(cbmc[4].counts, reference[4].counts), 100.0) << "closure distance 2-4";
  EXPECT_LT(chiSquaredTwoSample(cbmc[5].counts, reference[5].counts), 100.0) << "seed guide distance 1-4";
}

TEST(CBMC_FIXED_ENDPOINT, fixed_bonds_are_closed_exactly)
{
  Fixture f("fixed-endpoint-hexane-fixed", kHexaneFixedBondsJson);
  RandomNumber random(31);
  CBMC::GrowResult initial = f.initialMolecule(random);
  const CBMC::GrowContext empty = f.context();

  for (const std::vector<std::size_t> &fixed : {std::vector<std::size_t>{0, 1, 4, 5}, std::vector<std::size_t>{0, 4, 5}})
  {
    std::vector<Atom> atoms = initial.atoms;
    std::size_t accepted = 0, generated = 0;
    for (std::size_t i = 0; i != 2000; ++i)
    {
      std::optional<CBMC::GrowResult> grown = CBMC::regrowMolecule(
          random, empty, f.component, initial.molecule, atoms,
          {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = fixed});
      if (!grown) continue;
      ++generated;
      for (std::size_t bond = 0; bond + 1 != grown->atoms.size(); ++bond)
      {
        EXPECT_NEAR((grown->atoms[bond + 1].position - grown->atoms[bond].position).length(), 1.54, 1e-9)
            << "bond " << bond << " after fixed-endpoint regrowth";
      }
      for (std::size_t bead : fixed)
      {
        EXPECT_EQ(grown->atoms[bead].position.x, atoms[bead].position.x);
        EXPECT_EQ(grown->atoms[bead].position.y, atoms[bead].position.y);
        EXPECT_EQ(grown->atoms[bead].position.z, atoms[bead].position.z);
      }
      CBMC::RetraceResult old = CBMC::retraceMolecule(
          random, empty, f.component, atoms,
          {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = fixed});
      EXPECT_TRUE(std::isfinite(old.logRosenbluthWeight));
      if (random.uniform() < std::exp(grown->logRosenbluthWeight - old.logRosenbluthWeight))
      {
        atoms = grown->atoms;
        ++accepted;
      }
    }
    // A stiff (62500 K/rad^2) three-bead segment between fixed ends is hard to close well: the
    // acceptance is a few percent; the single-bead-plus-closure set accepts far more often.
    EXPECT_GT(generated, 1000u) << "fixed set size " << fixed.size();
    EXPECT_GT(accepted, 10u) << "fixed set size " << fixed.size();
  }
}
