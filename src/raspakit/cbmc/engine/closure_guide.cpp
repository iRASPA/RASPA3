module;

module cbmc_closure_guide;

import std;

import double3;
import randomnumbers;
import bond_potential;
import bend_potential;
import torsion_potential;
import cbmc_constants;

namespace
{
namespace Constants = CBMC::Constants;

double sampleBondLength(RandomNumber &random, double beta, const std::optional<BondPotential> &bond)
{
  return bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
}

// The maximum reach of the path: the sum over its bonds of the largest length the bond sampler
// produces in a fixed number of draws (exact for Fixed bonds, a ~4 sigma bound for harmonic ones),
// plus a margin. Distances beyond it are floored in the table.
double estimateReach(RandomNumber &random, double beta, const CBMC::ClosureGuidePath &path)
{
  double reach = 0.0;
  for (const std::optional<BondPotential> &bond : path.bonds)
  {
    double longest = 0.0;
    for (std::size_t s = 0; s != Constants::closureGuideReachSamples; ++s)
    {
      longest = std::max(longest, sampleBondLength(random, beta, bond));
    }
    reach += longest;
  }
  return reach + Constants::closureGuideReachMargin;
}

// Minimum of a torsion potential over the dihedral (a fine grid; the rejection sampler of the ideal
// sub-chain needs an envelope). Evaluated on a canonical geometry with dihedral phi.
double torsionMinimumEnergy(const TorsionPotential &torsion)
{
  const double3 posB{0.0, 0.0, 0.0};
  const double3 posC{0.0, 0.0, 1.0};
  const double3 posA{1.0, 0.0, 0.0};
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i != Constants::torsionMinimumGridPoints; ++i)
  {
    const double phi = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(Constants::torsionMinimumGridPoints);
    const double3 posD = posC + double3{std::cos(phi), std::sin(phi), 0.0};
    minimum = std::min(minimum, torsion.calculateEnergy(posA, posB, posC, posD));
  }
  return minimum;
}

// Converts a table of guide values to log space with the relative floor (see ClosureGuideTable).
CBMC::ClosureGuideTable finishTable(std::vector<double> guide, double spacing)
{
  const double maximum = *std::max_element(guide.begin(), guide.end());
  CBMC::ClosureGuideTable table{};
  table.spacing = spacing;
  table.logGuide.resize(guide.size());
  if (!(maximum > 0.0))
  {
    // Degenerate (no feasible geometry anywhere): a flat guide, i.e. no bias.
    std::fill(table.logGuide.begin(), table.logGuide.end(), 0.0);
    return table;
  }
  const double logMaximum = std::log(maximum);
  const double logFloor = logMaximum - Constants::closureGuideLogFloor;
  for (std::size_t i = 0; i != guide.size(); ++i)
  {
    table.logGuide[i] = guide[i] > 0.0 ? std::max(std::log(guide[i]), logFloor) : logFloor;
  }
  return table;
}

// Two-bond path: the guide is the closure-step partition function as a function of the distance D
// between the grown bead (a_0) and the target (a_2),
//
//   g(D) = < exp(-beta u_bend(theta(r1, r2, D))) / (r1 r2 D) >_{r1, r2}
//
// with the average over the two bond-length densities r^2 exp(-beta u_bond) (the bipolar volume element
// r1 r2 / D of the closure bead's position divided by the r1^2 r2^2 of the samplers; see
// cbmc_bridge_closure), and zero for a pair (r1, r2) that cannot span D. Evaluated by direct
// averaging over a fixed set of pairs on every grid point: smooth, exact in the tails, and it
// resolves the narrow peak of a stiff bend that a histogram would not.
CBMC::ClosureGuideTable buildTwoBondTable(RandomNumber &random, double beta, const CBMC::ClosureGuidePath &path,
                                          double reach)
{
  const double spacing = Constants::closureGuideSpacing;
  const std::size_t numberOfPoints = static_cast<std::size_t>(std::ceil(reach / spacing)) + 2;

  std::vector<std::pair<double, double>> pairs(Constants::closureGuideBondPairs);
  for (auto &[r1, r2] : pairs)
  {
    r1 = sampleBondLength(random, beta, path.bonds[0]);
    r2 = sampleBondLength(random, beta, path.bonds[1]);
  }

  const std::optional<BendPotential> &bend = path.bends[0];
  const double3 posB{0.0, 0.0, 0.0};

  std::vector<double> guide(numberOfPoints, 0.0);
  for (std::size_t k = 1; k != numberOfPoints; ++k)
  {
    const double D = static_cast<double>(k) * spacing;
    double sum = 0.0;
    for (const auto &[r1, r2] : pairs)
    {
      if (D > r1 + r2 || D < std::abs(r1 - r2)) continue;
      double energy = 0.0;
      if (bend.has_value())
      {
        const double cosTheta = std::clamp((r1 * r1 + r2 * r2 - D * D) / (2.0 * r1 * r2), -1.0, 1.0);
        const double theta = std::acos(cosTheta);
        const double3 posA{r1, 0.0, 0.0};
        const double3 posC{r2 * cosTheta, r2 * std::sin(theta), 0.0};
        energy = bend->calculateEnergy(posA, posB, posC, std::nullopt);
      }
      sum += std::exp(-beta * energy) / (r1 * r2 * D);
    }
    guide[k] = sum / static_cast<double>(pairs.size());
  }
  return finishTable(std::move(guide), spacing);
}

// Longer path: the radial density of the target around the grown bead, P(D) / (4 pi D^2), from a
// histogram of the end-to-end distance of the ideal sub-chain (bond lengths and bend angles from
// their Boltzmann samplers, torsions by rejection against their minimum). The first bend (to the
// bead's own anchor) and the torsions that reach outside the path are not part of the path and are
// left free; the guide is an approximation for efficiency only.
CBMC::ClosureGuideTable buildSampledTable(RandomNumber &random, double beta, const CBMC::ClosureGuidePath &path,
                                          double reach)
{
  const double spacing = Constants::closureGuideSpacing;
  const std::size_t numberOfPoints = static_cast<std::size_t>(std::ceil(reach / spacing)) + 2;
  const std::size_t numberOfBonds = path.numberOfBonds();

  std::vector<double> torsionMinima(path.torsions.size(), 0.0);
  for (std::size_t t = 0; t != path.torsions.size(); ++t)
  {
    if (path.torsions[t].has_value()) torsionMinima[t] = torsionMinimumEnergy(path.torsions[t].value());
  }

  std::vector<double> counts(numberOfPoints, 0.0);
  std::vector<double3> positions(numberOfBonds + 1);
  for (std::size_t sample = 0; sample != Constants::closureGuideChainSamples; ++sample)
  {
    positions[0] = double3{0.0, 0.0, 0.0};
    positions[1] = double3{0.0, 0.0, sampleBondLength(random, beta, path.bonds[0])};
    for (std::size_t i = 1; i != numberOfBonds; ++i)
    {
      // Place a_{i+1} from the centre a_i with the reference a_{i-1}.
      const double length = sampleBondLength(random, beta, path.bonds[i]);
      const std::optional<BendPotential> &bend = path.bends[i - 1];
      const double3 axis = (positions[i - 1] - positions[i]).normalized();
      const std::optional<TorsionPotential> &torsion = i >= 2 ? path.torsions[i - 2] : std::nullopt;

      double3 candidate{};
      for (std::size_t attempt = 0; attempt != Constants::closureGuideTorsionAttempts; ++attempt)
      {
        const double3 direction = bend.has_value() ? random.randomVectorOnCone(axis, bend->generateBendAngle(random, beta))
                                                   : random.randomVectorOnUnitSphere();
        candidate = positions[i] + length * direction;
        if (!torsion.has_value()) break;
        const double energy = torsion->calculateEnergy(positions[i - 2], positions[i - 1], positions[i], candidate);
        if (random.uniform() < std::exp(-beta * (energy - torsionMinima[i - 2]))) break;
      }
      positions[i + 1] = candidate;
    }

    const double distance = (positions[numberOfBonds] - positions[0]).length();
    const std::size_t bin = static_cast<std::size_t>(distance / spacing);
    if (bin < numberOfPoints) counts[bin] += 1.0;
  }

  // Bin i covers [i h, (i+1) h); its density estimate is assigned to the grid point at its centre by
  // attributing the count to the point i (the interpolation smooths the half-bin offset).
  std::vector<double> guide(numberOfPoints, 0.0);
  for (std::size_t i = 1; i != numberOfPoints; ++i)
  {
    const double D = (static_cast<double>(i) + 0.5) * spacing;
    guide[i] = counts[i] / (4.0 * std::numbers::pi * D * D * spacing * static_cast<double>(Constants::closureGuideChainSamples));
  }
  return finishTable(std::move(guide), spacing);
}
}  // namespace

std::string CBMC::ClosureGuidePath::signature() const
{
  std::string key{};
  for (const std::optional<BondPotential> &bond : bonds)
  {
    key += bond.has_value() ? std::format(";bond{}", static_cast<std::size_t>(bond->type)) : ";bond-none";
    if (bond.has_value())
      for (double p : bond->parameters) key += std::format(",{:.12e}", p);
  }
  for (const std::optional<BendPotential> &bend : bends)
  {
    key += bend.has_value() ? std::format(";bend{}", static_cast<std::size_t>(bend->type)) : ";bend-none";
    if (bend.has_value())
      for (double p : bend->parameters) key += std::format(",{:.12e}", p);
  }
  for (const std::optional<TorsionPotential> &torsion : torsions)
  {
    key += torsion.has_value() ? std::format(";torsion{}", static_cast<std::size_t>(torsion->type)) : ";torsion-none";
    if (torsion.has_value())
      for (double p : torsion->parameters) key += std::format(",{:.12e}", p);
  }
  return key;
}

double CBMC::ClosureGuideTable::logGuideAt(double distance) const
{
  if (logGuide.empty()) return 0.0;
  if (!(distance > 0.0)) return logGuide.front();
  const double x = distance / spacing;
  const std::size_t i = static_cast<std::size_t>(x);
  if (i + 1 >= logGuide.size()) return logGuide.back();
  const double f = x - static_cast<double>(i);
  return (1.0 - f) * logGuide[i] + f * logGuide[i + 1];
}

CBMC::ClosureGuideTable CBMC::buildClosureGuideTable(double beta, const ClosureGuidePath &path)
{
  if (path.numberOfBonds() < 2 || path.bends.size() + 1 != path.numberOfBonds() ||
      path.torsions.size() + 2 != path.numberOfBonds())
  {
    throw std::logic_error(std::format(
        "CBMC: a closure-guide path needs m >= 2 bonds, m-1 bends and m-2 torsions (got {}, {}, {})\n",
        path.numberOfBonds(), path.bends.size(), path.torsions.size()));
  }

  RandomNumber random(Constants::closureGuideSeed);  // frozen: one table per path signature
  const double reach = estimateReach(random, beta, path);
  return path.numberOfBonds() == 2 ? buildTwoBondTable(random, beta, path, reach)
                                   : buildSampledTable(random, beta, path, reach);
}
