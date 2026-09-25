#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import double3;
import atom;
import units;
import forcefield;
import component;
import framework;
import simulationbox;
import interpolation_energy_grid;
import randomnumbers;
import mc_moves_probabilities;
import cbmc;
import cbmc_chain_data;
import cbmc_growth_context;
import cbmc_recoil_growth;

// Detailed-balance tests of recoil growth (RG) as a chain-growth scheme. RG and configurational-bias
// (CBMC) growth are two proposal schemes for the same regrow move; with their own Rosenbluth-like
// weights in the acceptance rule both must sample the same Boltzmann distribution.

namespace
{

// Pseudo-atoms of the probe chain and its environment: 'B' bulky chain beads, 'S' a small chain bead,
// 'M' a medium chain bead, 'W' wall beads.
constexpr std::uint16_t kWallType = 3;

struct ProbeParameters
{
  double epsilonWall;
  double sigmaBulky;
  double sigmaSmall;
  double sigmaMedium;
};

ForceField makeProbeForceField(const ProbeParameters &p)
{
  return ForceField({{"B", false, 14.0, 0.0, 0.0, 6, false},
                     {"S", false, 14.0, 0.0, 0.0, 6, false},
                     {"M", false, 15.0, 0.0, 0.0, 6, false},
                     {"W", false, 12.0, 0.0, 0.0, 6, false}},
                    {{40.0, p.sigmaBulky}, {40.0, p.sigmaSmall}, {40.0, p.sigmaMedium}, {p.epsilonWall, 3.6}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

// A four-bead chain B-B-S-M grown from bead 0: a seed step (bead 1), a bend step (bead 2), and a
// torsion step (bead 3). Recoil growth decides the availability of a bend-step direction by whether a
// torsion-step bead can be placed after it -- through the growth attempt for the directions it tries,
// through an explicit feeler for the others -- so the torsion-step bead is where the two probes must
// agree.
constexpr std::string_view kProbeChainJson =
R"({
  "CriticalTemperature" : 425.125,
  "CriticalPressure" : 3796000.0,
  "AcentricFactor" : 0.201,
  "StartingBead" : 0,
  "pseudoAtoms" :
    [
      ["B", [0.0, 0.0, 0.0]],
      ["B", [1.54, 0.0, 0.0]],
      ["S", [2.17, 1.41, 0.0]],
      ["M", [3.71, 1.41, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3]
  ],
  "Bonds" : [
    [["B", "B"], "FIXED", [1.54]],
    [["B", "S"], "FIXED", [1.54]],
    [["S", "M"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["B", "B", "S"], "HARMONIC", [62500.0, 114]],
    [["B", "S", "M"], "HARMONIC", [62500.0, 114]]
  ],
  "Torsions" : [
    [["B", "B", "S", "M"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ]
}
)";

// A narrow tube of frozen W beads (their own 'molecule', id 1) along y through the center of the
// x-z plane: rings of 'beadsPerRing' beads of radius 'tubeRadius', every 'ringSpacing' along y. The
// rest of the box is empty bulk.
//
// Inside the tube the bulky first two beads only fit near the axis, which aligns the first bond with
// the axis; the small bead 2 fits anywhere on its bend cone, so the plane of beads 0-1-2 takes a
// random azimuth about the axis; and the medium bead 3 fits at the radius of bead 2 but not further
// out. By the cylindrical symmetry about the aligned first bond, a trans spin of bead 3 keeps it at
// the radius of bead 2 (open) for EVERY azimuth of bead 2, whereas a gauche spin swings it outward
// into the wall (closed). The torsion-selected spin of the growth is therefore nearly always open in
// the tube while a uniform spin mostly is not; in the bulk both are always open. The Boltzmann
// partition of the chain between tube and bulk is thus sensitive to any mismatch in how the recoil
// feelers probe openness on grow versus retrace, which a homogeneous environment would hide (a
// configuration-independent mismatch cancels in W_new / W_old).
std::vector<Atom> makeTube(double3 boxLengths, double tubeRadius, std::size_t beadsPerRing, double ringSpacing)
{
  std::vector<Atom> walls{};
  const std::size_t rings = static_cast<std::size_t>(std::llround(boxLengths.y / ringSpacing));
  for (std::size_t ring = 0; ring != rings; ++ring)
  {
    for (std::size_t j = 0; j != beadsPerRing; ++j)
    {
      double angle = 2.0 * std::numbers::pi * static_cast<double>(j) / static_cast<double>(beadsPerRing);
      double3 position(0.5 * boxLengths.x + tubeRadius * std::cos(angle), ringSpacing * static_cast<double>(ring),
                       0.5 * boxLengths.z + tubeRadius * std::sin(angle));
      walls.emplace_back(position, 0.0, 1.0, 1.0, std::uint32_t{1}, kWallType, std::uint8_t{1}, std::uint8_t{0},
                         std::uint8_t{0});
    }
  }
  return walls;
}

double dihedralCosine(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  double3 b1 = b - a;
  double3 b2 = c - b;
  double3 b3 = d - c;
  double3 n1 = double3::cross(b1, b2).normalized();
  double3 n2 = double3::cross(b2, b3).normalized();
  return std::clamp(double3::dot(n1, n2), -1.0, 1.0);
}

// Two-sample chi-squared for equal-size samples.
double chiSquaredTwoSample(const std::vector<double> &a, const std::vector<double> &b)
{
  double chi2 = 0.0;
  for (std::size_t i = 0; i != a.size(); ++i)
  {
    double n1 = a[i];
    double n2 = b[i];
    if (n1 + n2 < 1.0) continue;
    chi2 += (n1 - n2) * (n1 - n2) / (n1 + n2);
  }
  return chi2;
}

// Distance of the molecule's center from the tube axis (in the x-z plane, minimum image).
double distanceToTubeAxis(std::span<const Atom> atoms, double3 boxLengths)
{
  double3 center{};
  for (const Atom &atom : atoms) center += atom.position;
  center = center / static_cast<double>(atoms.size());
  double dx = center.x - 0.5 * boxLengths.x;
  double dz = center.z - 0.5 * boxLengths.z;
  dx -= boxLengths.x * std::round(dx / boxLengths.x);
  dz -= boxLengths.z * std::round(dz / boxLengths.z);
  return std::sqrt(dx * dx + dz * dz);
}

struct RegrowChainResult
{
  std::vector<double> torsionCounts;
  std::vector<double> radiusCounts;  ///< Distance of the molecule center from the tube axis.
  double tubeFraction;
  std::size_t accepted;
  std::size_t recorded;
};

// A one-molecule Markov chain of full regrow moves: grow a fresh molecule with a random first bead,
// retrace the current one, accept with min(1, W_new / W_old). The stationary distribution of the
// molecule is Boltzmann for either growth scheme if and only if grow and retrace are mutually
// consistent, which for recoil growth requires the tried-and-failed directions of the growth and the
// explicitly probed alternatives to be judged 'available' by the same random experiment.
RegrowChainResult runRegrowChain(RandomNumber &random, const ForceField &forceField, Component &component,
                                 const SimulationBox &box, std::span<const Atom> obstacles, double beta,
                                 double tubeRadius, std::size_t iterations, std::size_t thinning, std::size_t bins)
{
  const std::optional<Framework> noFramework{};
  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids(forceField.pseudoAtoms.size() + 1);
  const std::optional<InterpolationEnergyGrid> noExternalFieldGrid{};

  const CBMC::GrowContext context{false,
                                  forceField,
                                  box,
                                  noGrids,
                                  noExternalFieldGrid,
                                  noFramework,
                                  std::span<const Atom>{},
                                  obstacles,
                                  beta,
                                  forceField.cutOffFrameworkVDW,
                                  forceField.cutOffMoleculeVDW,
                                  forceField.cutOffCoulomb};

  const double3 boxLengths(box.lengthA, box.lengthB, box.lengthC);
  const double maximumRadius = 0.5 * std::min(boxLengths.x, boxLengths.z);

  // Initial state: the first completed grow.
  std::vector<Atom> state{};
  for (;;)
  {
    std::optional<ChainGrowData> grown =
        CBMC::growMoleculeSwapInsertion(random, context, component, 0, 0, 1.0, std::uint8_t{0}, false);
    if (grown.has_value())
    {
      state = grown->atoms;
      break;
    }
  }

  RegrowChainResult result{std::vector<double>(bins, 0.0), std::vector<double>(bins, 0.0), 0.0, 0, 0};
  auto binOf = [&](double x, double lo, double hi)
  {
    return std::min(bins - 1, static_cast<std::size_t>(std::clamp((x - lo) / (hi - lo), 0.0, 1.0 - 1e-12) *
                                                       static_cast<double>(bins)));
  };

  for (std::size_t i = 0; i != iterations; ++i)
  {
    std::optional<ChainGrowData> grown =
        CBMC::growMoleculeSwapInsertion(random, context, component, 0, 0, 1.0, std::uint8_t{0}, false);
    if (grown.has_value())
    {
      ChainRetraceData retraced =
          CBMC::retraceMoleculeSwapDeletion(random, context, component, std::span<Atom>(state));
      if (random.uniform() < std::exp(grown->logRosenbluthWeight - retraced.logRosenbluthWeight))
      {
        state = grown->atoms;
        ++result.accepted;
      }
    }

    if (i % thinning == thinning - 1)
    {
      double cosPhi = dihedralCosine(state[0].position, state[1].position, state[2].position, state[3].position);
      double radius = distanceToTubeAxis(state, boxLengths);
      result.torsionCounts[binOf(cosPhi, -1.0, 1.0)] += 1.0;
      result.radiusCounts[binOf(radius, 0.0, maximumRadius)] += 1.0;
      if (radius < tubeRadius) result.tubeFraction += 1.0;
      ++result.recorded;
    }
  }
  result.tubeFraction /= static_cast<double>(result.recorded);
  return result;
}

}  // namespace

// Recoil growth versus configurational-bias growth as proposal schemes of the same regrow move on a
// four-bead chain partitioning between a narrow tube and an empty bulk region. Both Markov chains must
// reach the same stationary distribution: the tube occupation and the radial-position and torsion
// histograms are compared against the CBMC reference. A recoil scheme whose feelers probe openness
// with a differently distributed spin than its growth attempts counts the available directions m_i
// inconsistently between grow and retrace; inside the tube (where openness depends on the spin) that
// changes the recoil weights relative to the bulk (where it does not), and the tube occupation drifts
// away from the CBMC reference.
//
// Calibration (250k regrow moves per chain): with consistent feelers the recoil chain reproduces the
// CBMC tube occupation (0.32) to within 0.01 and the radial chi-squared stays below ~100; with feelers
// spun uniformly instead of torsion-selected the tube occupation drops to 0.27 and the radial
// chi-squared exceeds 500.
TEST(CBMC_RECOIL_GROWTH, regrow_markov_chain_matches_cbmc_in_tube)
{
  constexpr double kTemperature = 300.0;                // [K]
  const double3 kBoxLengths(20.0, 16.0, 20.0);          // [A]
  constexpr double kTubeRadius = 5.0;                   // [A]
  constexpr std::size_t kBeadsPerRing = 16;
  constexpr double kRingSpacing = 2.0;                  // [A]

  // Wall attraction strong enough for the tube to hold about a third of the population; bead sizes
  // such that only the bulky beads are confined to the axis (sigma_BW = 5.0 = tube radius).
  const ProbeParameters parameters{15.0, 6.4, 2.6, 4.0};

  ForceField cbmcForceField = makeProbeForceField(parameters);
  cbmcForceField.numberOfTrialDirections = 8;
  cbmcForceField.numberOfTorsionTrialDirections = 10;
  cbmcForceField.numberOfFirstBeadPositions = 25;

  ForceField recoilForceField = cbmcForceField;
  recoilForceField.useRecoilGrowth = true;
  recoilForceField.recoilGrowthNumberOfTrialDirections = 3;
  recoilForceField.recoilGrowthMaximumRecoilLength = 2;

  TemporaryFile file("recoil-probe-chain.json", kProbeChainJson);
  Component cbmcChain(Component::Type::Adsorbate, 0, cbmcForceField, "recoil-probe-chain", file.stemPath().string(),
                      5, 21, MCMoveProbabilities(), std::nullopt, false);
  Component recoilChain(Component::Type::Adsorbate, 0, recoilForceField, "recoil-probe-chain",
                        file.stemPath().string(), 5, 21, MCMoveProbabilities(), std::nullopt, false);

  const SimulationBox box(kBoxLengths.x, kBoxLengths.y, kBoxLengths.z);
  const std::vector<Atom> tube = makeTube(kBoxLengths, kTubeRadius, kBeadsPerRing, kRingSpacing);
  const double beta = 1.0 / (Units::KB * kTemperature);

  constexpr std::size_t iterations = 250'000;
  constexpr std::size_t thinning = 5;
  constexpr std::size_t bins = 20;

  RandomNumber randomCBMC(4711);
  RegrowChainResult cbmc = runRegrowChain(randomCBMC, cbmcForceField, cbmcChain, box, tube, beta, kTubeRadius,
                                          iterations, thinning, bins);
  RandomNumber randomRecoil(1999);
  RegrowChainResult recoil = runRegrowChain(randomRecoil, recoilForceField, recoilChain, box, tube, beta,
                                            kTubeRadius, iterations, thinning, bins);

  ASSERT_EQ(cbmc.recorded, recoil.recorded);
  EXPECT_GT(cbmc.accepted, iterations / 20) << "CBMC regrow acceptance collapsed";
  EXPECT_GT(recoil.accepted, iterations / 20) << "recoil-growth regrow acceptance collapsed";

  // Both chains must actually visit tube and bulk for the partition to be a meaningful observable.
  EXPECT_GT(cbmc.tubeFraction, 0.15) << "the tube is hardly populated; the test has lost its sensitivity";
  EXPECT_LT(cbmc.tubeFraction, 0.60) << "the bulk is hardly populated; the test has lost its sensitivity";

  const double torsionChi2 = chiSquaredTwoSample(cbmc.torsionCounts, recoil.torsionCounts);
  const double radiusChi2 = chiSquaredTwoSample(cbmc.radiusCounts, recoil.radiusCounts);

  // The molecule exchanges between tube and bulk only every ~100 moves, so the thinned samples are
  // strongly correlated and the chi-squared sits well above its independent-sample scale (~bins - 1);
  // the thresholds leave room for that while remaining far below the inconsistent-feeler values.
  EXPECT_LT(torsionChi2, 250.0) << "recoil-growth torsion distribution deviates from the CBMC reference";
  EXPECT_LT(radiusChi2, 250.0) << "recoil-growth radial distribution deviates from the CBMC reference";
  EXPECT_NEAR(recoil.tubeFraction, cbmc.tubeFraction, 0.05)
      << "recoil growth partitions the chain between tube and bulk differently from CBMC";
}

// The retrace of an existing molecule that overlaps with its surroundings has no defined recoil
// weight (the openness probability of the old direction is zero). An accepted state never overlaps,
// so this is an inconsistent simulation state and must be reported, not silently weighted as if the
// overlapping bead had zero energy.
TEST(CBMC_RECOIL_GROWTH, retrace_of_overlapping_old_configuration_throws)
{
  ForceField forceField = makeProbeForceField(ProbeParameters{15.0, 6.4, 2.6, 4.0});
  forceField.useRecoilGrowth = true;
  forceField.recoilGrowthNumberOfTrialDirections = 3;
  forceField.recoilGrowthMaximumRecoilLength = 2;

  TemporaryFile file("recoil-probe-chain.json", kProbeChainJson);
  Component chain(Component::Type::Adsorbate, 0, forceField, "recoil-probe-chain", file.stemPath().string(), 5, 21,
                  MCMoveProbabilities(), std::nullopt, false);

  const SimulationBox box(30.0, 30.0, 30.0);
  const double beta = 1.0 / (Units::KB * 300.0);
  const std::optional<Framework> noFramework{};
  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids(forceField.pseudoAtoms.size() + 1);
  const std::optional<InterpolationEnergyGrid> noExternalFieldGrid{};

  auto makeContext = [&](std::span<const Atom> obstacles)
  {
    return CBMC::GrowContext{false,
                             forceField,
                             box,
                             noGrids,
                             noExternalFieldGrid,
                             noFramework,
                             std::span<const Atom>{},
                             obstacles,
                             beta,
                             forceField.cutOffFrameworkVDW,
                             forceField.cutOffMoleculeVDW,
                             forceField.cutOffCoulomb};
  };
  auto makeWallBead = [](double3 position)
  {
    return Atom(position, 0.0, 1.0, 1.0, std::uint32_t{1}, kWallType, std::uint8_t{1}, std::uint8_t{0},
                std::uint8_t{0});
  };

  // Grow a molecule in an empty box.
  RandomNumber random(12345);
  std::vector<Atom> molecule{};
  {
    const CBMC::GrowContext empty = makeContext(std::span<const Atom>{});
    for (;;)
    {
      std::optional<ChainGrowData> grown =
          CBMC::growMoleculeSwapInsertion(random, empty, chain, 0, 0, 1.0, std::uint8_t{0}, false);
      if (grown.has_value())
      {
        molecule = grown->atoms;
        break;
      }
    }
  }

  // The chain is retraced from its starting bead (already placed).
  const std::vector<std::size_t> placed{chain.startingBead};

  // A single frozen wall bead (its own molecule, id 1) placed onto the last-grown bead of the molecule:
  // that bead now overlaps and the retrace has to refuse.
  const std::vector<Atom> obstacle{makeWallBead(molecule[3].position + double3(0.3, 0.0, 0.0))};
  const CBMC::GrowContext overlapping = makeContext(obstacle);
  try
  {
    (void)CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, overlapping, chain, std::span<Atom>(molecule),
                                                         placed);
    FAIL() << "the retrace of an overlapping configuration returned a weight";
  }
  catch (const std::runtime_error &error)
  {
    EXPECT_NE(std::string_view(error.what()).find("overlaps at growth step"), std::string_view::npos) << error.what();
  }

  // The same error must propagate through the public CBMC entry point (the move-level API), so that
  // the driver can report it; a 'noexcept' anywhere on that path would turn it into std::terminate.
  try
  {
    (void)CBMC::retraceMoleculeSwapDeletion(random, overlapping, chain, std::span<Atom>(molecule));
    FAIL() << "the public retrace of an overlapping configuration returned a weight";
  }
  catch (const std::runtime_error &error)
  {
    EXPECT_NE(std::string_view(error.what()).find("overlaps at growth step"), std::string_view::npos) << error.what();
  }

  // The same wall bead moved well away from the molecule is a valid environment: no throw.
  const std::vector<Atom> distant{makeWallBead(molecule[3].position + double3(10.0, 0.0, 0.0))};
  const CBMC::GrowContext valid = makeContext(distant);
  EXPECT_NO_THROW(
      (void)CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, valid, chain, std::span<Atom>(molecule), placed));
}
