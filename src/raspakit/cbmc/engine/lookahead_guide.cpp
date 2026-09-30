module;

module cbmc_lookahead_guide;

import std;

import double3;
import randomnumbers;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import intra_molecular_potentials;
import cbmc_constants;

namespace
{
namespace Constants = CBMC::Constants;

double sampleBondLength(RandomNumber &random, double beta, const std::optional<BondPotential> &bond)
{
  return bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
}

// Minimum of a torsion potential over the dihedral on a canonical geometry: the envelope of the
// rejection sampler of a modelled bead.
double torsionMinimumEnergy(const TorsionPotential &torsion)
{
  const double3 posB{0.0, 0.0, 0.0};
  const double3 posC{0.0, 0.0, 1.0};
  const double3 posA{1.0, 0.0, 0.0};
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i != Constants::torsionMinimumGridPoints; ++i)
  {
    const double phi =
        2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(Constants::torsionMinimumGridPoints);
    const double3 posD = posC + double3{std::cos(phi), std::sin(phi), 0.0};
    minimum = std::min(minimum, torsion.calculateEnergy(posA, posB, posC, posD));
  }
  return minimum;
}

// Minimum of a bend potential over [0, pi] (zero for the harmonic forms with an equilibrium angle in
// range; a fine grid lowered by a small slack otherwise, the grid can only overestimate the minimum).
double bendMinimumEnergy(const BendPotential &bend)
{
  if ((bend.type == BendType::Harmonic || bend.type == BendType::CoreShell) && bend.parameters[1] >= 0.0 &&
      bend.parameters[1] <= std::numbers::pi)
  {
    return 0.0;
  }
  const double3 posA{1.0, 0.0, 0.0};
  const double3 posB{0.0, 0.0, 0.0};
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i != Constants::bendMinimumGridPoints; ++i)
  {
    const double theta =
        std::numbers::pi * static_cast<double>(i) / static_cast<double>(Constants::bendMinimumGridPoints - 1);
    const double3 posC{std::cos(theta), std::sin(theta), 0.0};
    minimum = std::min(minimum, bend.calculateEnergy(posA, posB, posC, std::nullopt));
  }
  return minimum - Constants::bendMinimumSlack;
}

// The envelope of a bead's rejection: the minima of its torsion and sibling bends (a lower bound of
// their sum, so the acceptance exp(-beta (u - u_min)) never exceeds one).
double rejectionEnvelope(const CBMC::LookaheadBead &bead)
{
  double envelope = 0.0;
  if (bead.torsion.has_value()) envelope += torsionMinimumEnergy(bead.torsion.value());
  for (const BendPotential &bend : bead.siblingBends) envelope += bendMinimumEnergy(bend);
  return envelope;
}

// Places a modelled bead from its parent (see LookaheadBead). The rejection imposes the torsion and
// the sibling bends; an exhausted budget keeps the last draw (a negligible bias of the guide only).
void placeBead(RandomNumber &random, double beta, const CBMC::LookaheadBead &bead, double envelope,
               std::vector<double3> &positions)
{
  const double3 parent = positions[bead.parent];
  const double length = sampleBondLength(random, beta, bead.bond);
  const bool hasCone = bead.bend.has_value() && bead.bendReference.has_value();
  const double3 axis = hasCone ? (positions[bead.bendReference.value()] - parent).normalized() : double3{0.0, 0.0, 1.0};
  const bool hasRejection = bead.torsion.has_value() || !bead.siblingBends.empty();

  for (std::size_t attempt = 0; attempt != Constants::lookaheadGuideRejectionAttempts; ++attempt)
  {
    const double3 direction = hasCone ? random.randomVectorOnCone(axis, bead.bend->generateBendAngle(random, beta))
                                      : random.randomVectorOnUnitSphere();
    positions[bead.bead] = parent + length * direction;
    if (!hasRejection) return;

    double energy = 0.0;
    if (bead.torsion.has_value())
    {
      const auto &ids = bead.torsion->identifiers;
      energy += bead.torsion->calculateEnergy(positions[ids[0]], positions[ids[1]], positions[ids[2]], positions[ids[3]]);
    }
    for (const BendPotential &bend : bead.siblingBends)
    {
      const auto &ids = bend.identifiers;
      energy += bend.calculateEnergy(positions[ids[0]], positions[ids[1]], positions[ids[2]], std::nullopt);
    }
    if (random.uniform() < std::exp(-beta * (energy - envelope))) return;
  }
}

// The energy of the cross terms of a model on a set of positions.
double crossEnergy(const Potentials::IntraMolecularPotentials &terms, const std::vector<double3> &positions)
{
  double energy = 0.0;
  for (const BendPotential &bend : terms.bends)
  {
    const auto &ids = bend.identifiers;
    energy += bend.calculateEnergy(positions[ids[0]], positions[ids[1]], positions[ids[2]], std::nullopt);
  }
  for (const TorsionPotential &torsion : terms.torsions)
  {
    const auto &ids = torsion.identifiers;
    energy += torsion.calculateEnergy(positions[ids[0]], positions[ids[1]], positions[ids[2]], positions[ids[3]]);
  }
  for (const VanDerWaalsPotential &pair : terms.vanDerWaals)
  {
    energy += pair.calculateEnergy(positions[pair.identifiers[0]], positions[pair.identifiers[1]]);
  }
  for (const CoulombPotential &pair : terms.coulombs)
  {
    energy += pair.calculateEnergy(positions[pair.identifiers[0]], positions[pair.identifiers[1]]);
  }
  return energy;
}

double wrapAngle(double angle)
{
  while (angle > std::numbers::pi) angle -= 2.0 * std::numbers::pi;
  while (angle < -std::numbers::pi) angle += 2.0 * std::numbers::pi;
  return angle;
}

// Converts a table of guide values to log space with the relative floor (see LookaheadGuideTable).
CBMC::LookaheadGuideTable finishTable(std::vector<double> guide)
{
  const double maximum = *std::max_element(guide.begin(), guide.end());
  CBMC::LookaheadGuideTable table{};
  table.logGuide.resize(guide.size());
  if (!(maximum > 0.0))
  {
    std::fill(table.logGuide.begin(), table.logGuide.end(), 0.0);  // degenerate: no bias
    return table;
  }
  const double logMaximum = std::log(maximum);
  const double logFloor = logMaximum - Constants::lookaheadGuideLogFloor;
  for (std::size_t i = 0; i != guide.size(); ++i)
  {
    table.logGuide[i] = guide[i] > 0.0 ? std::max(std::log(guide[i]), logFloor) : logFloor;
  }
  return table;
}
}  // namespace

double CBMC::dihedralAngle(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  const double3 b1 = b - a;
  const double3 b2 = c - b;
  const double3 b3 = d - c;
  const double3 n1 = double3::cross(b1, b2);
  const double3 n2 = double3::cross(b2, b3);
  const double y = b2.length() * double3::dot(b1, n2);
  const double x = double3::dot(n1, n2);
  return std::atan2(y, x);
}

double CBMC::LookaheadGuideTable::logGuideAt(double dihedral) const
{
  const std::size_t size = logGuide.size();
  if (size == 0) return 0.0;
  const double x = (wrapAngle(dihedral) + std::numbers::pi) / (2.0 * std::numbers::pi) * static_cast<double>(size);
  const double floorX = std::floor(x);
  const std::size_t i = static_cast<std::size_t>(floorX) % size;
  const double f = x - floorX;
  return (1.0 - f) * logGuide[i] + f * logGuide[(i + 1) % size];
}

std::string CBMC::LookaheadGuideModel::signature() const
{
  std::map<std::size_t, std::string> roles{};
  roles[previousBead] = "p";
  roles[currentBead] = "c";
  for (std::size_t i = 0; i != backwardBeads.size(); ++i) roles[backwardBeads[i].bead] = std::format("b{}", i);
  for (std::size_t i = 0; i != nextBeads.size(); ++i) roles[nextBeads[i].bead] = std::format("n{}", i);
  for (std::size_t i = 0; i != futureBeads.size(); ++i) roles[futureBeads[i].bead] = std::format("f{}", i);
  const auto roleOf = [&](std::size_t id) -> std::string
  {
    auto it = roles.find(id);
    return it == roles.end() ? std::format("?{}", id) : it->second;
  };
  const auto appendParameters = [](std::string &key, const auto &parameters)
  {
    for (double p : parameters) key += std::format(",{:.12e}", p);
  };
  const auto appendTerm = [&](std::string &key, std::string_view tag, std::size_t type, const auto &identifiers,
                              const auto &parameters)
  {
    key += std::format(";{}{}", tag, type);
    for (std::size_t id : identifiers) key += ":" + roleOf(id);
    appendParameters(key, parameters);
  };

  std::string key{};
  key += std::format("depth{};grid{}", Constants::lookaheadGuideDepth, Constants::lookaheadGuideGridPoints);
  key += axisBond.has_value() ? std::format(";axis{}", static_cast<std::size_t>(axisBond->type)) : ";axis-none";
  if (axisBond.has_value()) appendParameters(key, axisBond->parameters);
  key += ";ref:" + roleOf(referenceBead) + ";first:" + roleOf(firstNextBead);

  const auto appendBead = [&](const LookaheadBead &bead)
  {
    key += ";[" + roleOf(bead.bead) + "<" + roleOf(bead.parent);
    key += bead.bendReference.has_value() ? "<" + roleOf(bead.bendReference.value()) : "<-";
    key += bead.torsionReference.has_value() ? "<" + roleOf(bead.torsionReference.value()) : "<-";
    key += bead.bond.has_value() ? std::format(";bond{}", static_cast<std::size_t>(bead.bond->type)) : ";bond-none";
    if (bead.bond.has_value()) appendParameters(key, bead.bond->parameters);
    key += bead.bend.has_value() ? std::format(";bend{}", static_cast<std::size_t>(bead.bend->type)) : ";bend-none";
    if (bead.bend.has_value()) appendParameters(key, bead.bend->parameters);
    if (bead.torsion.has_value())
      appendTerm(key, "tor", static_cast<std::size_t>(bead.torsion->type), bead.torsion->identifiers,
                 bead.torsion->parameters);
    for (const BendPotential &bend : bead.siblingBends)
      appendTerm(key, "sib", static_cast<std::size_t>(bend.type), bend.identifiers, bend.parameters);
    key += "]";
  };
  for (const LookaheadBead &bead : backwardBeads) appendBead(bead);
  for (const LookaheadBead &bead : nextBeads) appendBead(bead);
  for (const LookaheadBead &bead : futureBeads) appendBead(bead);

  const auto appendTerms = [&](std::string_view group, const Potentials::IntraMolecularPotentials &terms)
  {
    key += std::format(";{}", group);
    for (const BendPotential &t : terms.bends)
      appendTerm(key, "B", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
    for (const TorsionPotential &t : terms.torsions)
      appendTerm(key, "T", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
    for (const VanDerWaalsPotential &t : terms.vanDerWaals)
    {
      appendTerm(key, "V", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
      key += std::format(",s{:.12e}", t.scaling);
    }
    for (const CoulombPotential &t : terms.coulombs)
    {
      appendTerm(key, "C", static_cast<std::size_t>(t.type), t.identifiers,
                 std::array<double, 3>{t.chargeA, t.chargeB, t.scaling});
    }
  };
  appendTerms("inv", invariantTerms);
  appendTerms("spin", spinTerms);
  return key;
}

CBMC::LookaheadGuideTable CBMC::buildLookaheadGuideTable(double beta, const LookaheadGuideModel &model)
{
  RandomNumber random(Constants::lookaheadGuideSeed);  // frozen: one table per model signature

  const std::size_t gridPoints = Constants::lookaheadGuideGridPoints;
  const double spacing = 2.0 * std::numbers::pi / static_cast<double>(gridPoints);

  // Rejection envelopes of the modelled beads (setup only).
  const auto envelopesOf = [](const std::vector<LookaheadBead> &beads)
  {
    std::vector<double> envelopes(beads.size());
    for (std::size_t i = 0; i != beads.size(); ++i) envelopes[i] = rejectionEnvelope(beads[i]);
    return envelopes;
  };
  const std::vector<double> backwardEnvelopes = envelopesOf(model.backwardBeads);
  const std::vector<double> nextEnvelopes = envelopesOf(model.nextBeads);
  const std::vector<double> futureEnvelopes = envelopesOf(model.futureBeads);

  // The beads that spin with the grown beads about the axis: the grown and the future beads.
  std::vector<std::size_t> spinning{};
  for (const LookaheadBead &bead : model.nextBeads) spinning.push_back(bead.bead);
  for (const LookaheadBead &bead : model.futureBeads) spinning.push_back(bead.bead);

  // One sub-molecule draw in the step frame: current bead at the origin, previous bead on +z.
  std::vector<double3> positions(model.numberOfAtoms, double3{0.0, 0.0, 0.0});
  const auto drawSubMolecule = [&]()
  {
    positions[model.currentBead] = double3{0.0, 0.0, 0.0};
    positions[model.previousBead] = double3{0.0, 0.0, sampleBondLength(random, beta, model.axisBond)};
    for (std::size_t i = 0; i != model.backwardBeads.size(); ++i)
      placeBead(random, beta, model.backwardBeads[i], backwardEnvelopes[i], positions);
    for (std::size_t i = 0; i != model.nextBeads.size(); ++i)
      placeBead(random, beta, model.nextBeads[i], nextEnvelopes[i], positions);
    for (std::size_t i = 0; i != model.futureBeads.size(); ++i)
      placeBead(random, beta, model.futureBeads[i], futureEnvelopes[i], positions);
  };
  const auto referenceDihedral = [&](const std::vector<double3> &p)
  {
    return dihedralAngle(p[model.referenceBead], p[model.previousBead], p[model.currentBead], p[model.firstNextBead]);
  };

  // The rotation about the axis (previous - current, through the current bead at the origin) and the
  // sense in which it advances the reference dihedral, determined once on a test draw.
  const double3 axis{0.0, 0.0, 1.0};
  std::vector<double3> spun(model.numberOfAtoms);
  const auto spinTo = [&](double angle)
  {
    spun = positions;
    for (std::size_t bead : spinning) spun[bead] = axis.rotateAroundAxis(positions[bead], angle);
  };
  drawSubMolecule();
  const double testDihedral = referenceDihedral(positions);
  spinTo(0.3);
  const double sense = wrapAngle(referenceDihedral(spun) - testDihedral) > 0.0 ? 1.0 : -1.0;

  std::vector<double> guide(gridPoints, 0.0);
  for (std::size_t sample = 0; sample != Constants::lookaheadGuideSamples; ++sample)
  {
    drawSubMolecule();
    const double invariantEnergy = crossEnergy(model.invariantTerms, positions);
    const double drawnDihedral = referenceDihedral(positions);
    for (std::size_t k = 0; k != gridPoints; ++k)
    {
      const double target = -std::numbers::pi + static_cast<double>(k) * spacing;
      spinTo(sense * wrapAngle(target - drawnDihedral));
      guide[k] += std::exp(-beta * (invariantEnergy + crossEnergy(model.spinTerms, spun)));
    }
  }
  for (double &g : guide) g /= static_cast<double>(Constants::lookaheadGuideSamples);
  return finishTable(std::move(guide));
}
