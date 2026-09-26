module;

module mc_moves_bridging_common;

import std;

import double3;
import double3x3;
import atom;
import randomnumbers;
import component;
import mc_moves_concerted_rotation_geometry;

namespace
{
// Two trimers are the same closure when every atom agrees to this distance [Angstrom]; the solver
// converges to ~1e-12, the tolerance only has to separate distinct solutions.
constexpr double sameTrimerTolerance = 1e-6;

bool sameTrimer(const ConcertedRotation::Trimer& a, const ConcertedRotation::Trimer& b)
{
  return (a.a3 - b.a3).length() < sameTrimerTolerance && (a.a4 - b.a4).length() < sameTrimerTolerance &&
         (a.a5 - b.a5).length() < sameTrimerTolerance;
}

// The seven backbone positions a1 ... a7 of a closure.
std::array<double3, 7> window(const Bridging::Anchors& anchors, const ConcertedRotation::Trimer& trimer)
{
  return {anchors.a1, anchors.a2, trimer.a3, trimer.a4, trimer.a5, anchors.a6, anchors.a7};
}
}  // namespace

std::optional<Bridging::Closure> Bridging::rebridge(RandomNumber& random, const Anchors& oldAnchors,
                                                    const ConcertedRotation::Trimer& oldTrimer,
                                                    const Anchors& newAnchors)
{
  const std::array<double3, 7> oldWindow = window(oldAnchors, oldTrimer);
  const ConcertedRotation::BackboneGeometry geometry =
      ConcertedRotation::BackboneGeometry::fromPositions(std::span<const double3, 7>(oldWindow));

  const std::vector<ConcertedRotation::Trimer> forwardSolutions =
      ConcertedRotation::rebridge(newAnchors.a1, newAnchors.a2, newAnchors.a6, newAnchors.a7, geometry);
  if (forwardSolutions.empty()) return std::nullopt;
  const ConcertedRotation::Trimer& chosen =
      forwardSolutions[static_cast<std::size_t>(random.uniform() * static_cast<double>(forwardSolutions.size()))];

  // The reverse move solves the closure between the old anchors with the same geometry (measured
  // from the new trimer, which carries it exactly) and must find the old trimer among its solutions.
  const std::vector<ConcertedRotation::Trimer> backwardSolutions =
      ConcertedRotation::rebridge(oldAnchors.a1, oldAnchors.a2, oldAnchors.a6, oldAnchors.a7, geometry);
  if (!std::ranges::any_of(backwardSolutions,
                           [&](const ConcertedRotation::Trimer& trimer) { return sameTrimer(trimer, oldTrimer); }))
  {
    return std::nullopt;
  }

  const std::array<double3, 7> newWindow = window(newAnchors, chosen);
  const double jacobianOld = ConcertedRotation::closureJacobian(std::span<const double3, 7>(oldWindow), oldAnchors.a8);
  const double jacobianNew = ConcertedRotation::closureJacobian(std::span<const double3, 7>(newWindow), newAnchors.a8);
  if (!(jacobianOld > 0.0) || !(jacobianNew > 0.0)) return std::nullopt;

  return Closure{.trimer = chosen,
                 .logWeight = std::log(static_cast<double>(forwardSolutions.size())) -
                              std::log(static_cast<double>(backwardSolutions.size())) + std::log(jacobianOld) -
                              std::log(jacobianNew)};
}

double Bridging::chartVolumeElement(const double3& a6, const double3& a7, const std::optional<double3>& a8)
{
  double3 u = a7 - a6;
  const double lengthSquared = double3::dot(u, u);
  if (!a8.has_value()) return lengthSquared;
  return lengthSquared * double3::cross(u.normalized(), a8.value() - a7).length();
}

void Bridging::placeTrimer(std::span<Atom> trialAtoms, std::span<const Atom> sourceAtoms,
                           const std::vector<Component::BridgingTopology::Unit>& units, std::size_t site,
                           const Anchors& oldAnchors, const ConcertedRotation::Trimer& oldTrimer,
                           const Anchors& newAnchors, const ConcertedRotation::Trimer& newTrimer)
{
  const std::array<double3, 5> oldPositions{oldAnchors.a2, oldTrimer.a3, oldTrimer.a4, oldTrimer.a5, oldAnchors.a6};
  const std::array<double3, 5> newPositions{newAnchors.a2, newTrimer.a3, newTrimer.a4, newTrimer.a5, newAnchors.a6};
  for (std::size_t m = 1; m <= 3; ++m)
  {
    const Component::BridgingTopology::Unit& unit = units[site + m];
    trialAtoms[unit.backboneAtom].position = newPositions[m];
    if (unit.sideAtoms.empty()) continue;
    const double3x3 oldFrame = ConcertedRotation::localFrame(oldPositions[m], oldPositions[m - 1], oldPositions[m + 1]);
    const double3x3 newFrame = ConcertedRotation::localFrame(newPositions[m], newPositions[m - 1], newPositions[m + 1]);
    for (std::size_t atom : unit.sideAtoms)
    {
      trialAtoms[atom].position = ConcertedRotation::transformRigidly(oldPositions[m], oldFrame, newPositions[m],
                                                                      newFrame, sourceAtoms[atom].position);
    }
  }
}
