module;

module cbmc_bridge_closure;

import std;

import atom;
import double3;
import randomnumbers;
import bond_potential;
import bend_potential;
import cbmc_constants;
import cbmc_grow_step;
import cbmc_closure_guide;

namespace
{
namespace Constants = CBMC::Constants;

double bondLength(RandomNumber &random, double beta, const std::optional<BondPotential> &bond)
{
  return bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
}

double midBendEnergy(const CBMC::GrowStep &step, const double3 &anchor, const double3 &next, const double3 &closure)
{
  if (!step.bridge.midBend.has_value()) return 0.0;
  return step.bridge.midBend->calculateEnergy(anchor, next, closure, std::nullopt);
}

// The bipolar base weight exp(-beta u_mid) / (r1 r2 D) of a closing-bead position; see the module
// description. A degenerate geometry (coincident beads) has no defined weight and returns zero.
double baseWeight(double beta, const CBMC::GrowStep &step, const double3 &anchor, const double3 &next,
                  const double3 &closure)
{
  const double r1 = (next - anchor).length();
  const double r2 = (next - closure).length();
  const double D = (closure - anchor).length();
  if (r1 < Constants::degenerateAxisLength || r2 < Constants::degenerateAxisLength || D < Constants::degenerateAxisLength)
  {
    return 0.0;
  }
  return std::exp(-beta * midBendEnergy(step, anchor, next, closure)) / (r1 * r2 * D);
}
}  // namespace

double3 CBMC::bridgeClosureAxis(std::span<const Atom> chainAtoms, const GrowStep &step)
{
  double3 axis = chainAtoms[step.closureBead.value()].position - chainAtoms[step.currentBead].position;
  if (axis.length() < Constants::degenerateAxisLength) axis = double3{0.0, 0.0, 1.0};
  return axis.normalized();
}

std::optional<CBMC::BridgeClosureBase> CBMC::sampleBridgeClosureBase(RandomNumber &random, double beta,
                                                                     std::span<const Atom> chainAtoms,
                                                                     const GrowStep &step)
{
  const std::size_t nextBead = step.nextBeads.front();
  const double3 anchor = chainAtoms[step.currentBead].position;
  const double3 closure = chainAtoms[step.closureBead.value()].position;
  const double D = (closure - anchor).length();
  if (D < Constants::degenerateAxisLength) return std::nullopt;

  const double r1 = bondLength(random, beta, step.bridge.anchorBond);
  const double r2 = bondLength(random, beta, step.bridge.closureBond);
  if (D > r1 + r2 || D < std::abs(r1 - r2)) return std::nullopt;

  // The circle: at angle psi from the axis as seen from the anchor, radius r1, uniform azimuth.
  const double3 axis = (closure - anchor) / D;
  const double cosPsi = std::clamp((r1 * r1 + D * D - r2 * r2) / (2.0 * r1 * D), -1.0, 1.0);

  Atom atom = chainAtoms[nextBead];
  atom.position = anchor + r1 * random.randomVectorOnCone(axis, std::acos(cosPsi));

  return BridgeClosureBase{atom, baseWeight(beta, step, anchor, atom.position, closure)};
}

double CBMC::bridgeClosureBaseWeight(double beta, std::span<const Atom> chainAtoms, const GrowStep &step)
{
  return baseWeight(beta, step, chainAtoms[step.currentBead].position, chainAtoms[step.nextBeads.front()].position,
                    chainAtoms[step.closureBead.value()].position);
}

void CBMC::prepareClosureGuides(double beta, std::span<GrowStep> plan,
                                std::map<std::string, std::shared_ptr<const ClosureGuideTable>> &memo)
{
  for (GrowStep &step : plan)
  {
    for (GrowStep::SpinSelectionData::ClosureGuide &guide : step.spin.guides)
    {
      auto it = memo.find(guide.signature);
      if (it == memo.end())
      {
        it = memo.emplace(guide.signature, std::make_shared<const ClosureGuideTable>(buildClosureGuideTable(beta, guide.path)))
                 .first;
      }
      guide.table = it->second;
    }
  }
}
