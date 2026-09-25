module;

export module cbmc;

import std;

import atom;
import molecule;
import double3;
import randomnumbers;
import running_energy;
import component;
export import cbmc_results;
export import cbmc_growth_context;

// The entry points of the CBMC / recoil-growth machinery: grow a new molecule, regrow (part of) an
// existing one, retrace an existing one, and correct a result from the inner to the full cut-offs.
//
// A grow or retrace has two stages: the first bead is placed with the 'FirstBeadScheme' of the
// request, the remaining beads are grown fragment by fragment with the operator engine using the
// chain scheme of the context ('GrowContext::settings.chainScheme'). The two Rosenbluth weights
// multiply (their logarithms add) and the energies add.
//
// Error contract:
//  - A grow returns std::nullopt when the trial molecule can not be constructed (every trial of some
//    step overlaps or falls below 'minimumRosenbluthFactor', a recoil-growth dead end, ...). This is
//    an ordinary outcome of the move -- the caller counts it as a rejection.
//  - A retrace has no such outcome: the old configuration is an accepted state and always has a
//    weight.
//  - The functions throw std::runtime_error when the simulation state itself is inconsistent: no
//    growth step can be built for the requested placed set, the exact base sampler exhausts its
//    rejection budget, or the existing molecule overlaps with its environment. No weight is defined
//    in these cases and silently continuing would corrupt the acceptance rule, so the error
//    propagates to the driver, which reports the message and stops. The entry points are therefore
//    deliberately NOT noexcept.
export namespace CBMC
{
/// Identity and scaling attributes stamped on every atom of a freshly grown molecule.
struct NewMoleculeIdentity
{
  std::size_t componentId;
  std::size_t moleculeId;
  double scaling{1.0};
  std::uint8_t groupId{0};
  bool isFractional{false};
};

/// How the first bead of a molecule is placed (grow) or weighted (retrace).
enum class FirstBeadScheme : std::size_t
{
  /// 'numberOfFirstBeadPositions' uniformly random positions in the box; weight = sum of Boltzmann
  /// factors / number of positions. Insertion and deletion.
  MultipleFirstBead = 0,
  /// The multiple-first-bead reinsertion of Esselink et al.: as above against a background without
  /// the molecule itself; the grow retains the partial weight 'GrowResult::firstBeadStoredR', which
  /// the retrace of the old configuration needs ('RetraceRequest::storedR').
  Reinsertion = 1,
  /// A single trial at 'GrowRequest::firstBeadPosition', its Boltzmann factor as weight. Identity
  /// change: the new molecule takes the position of the old one.
  Pinned = 2,
  /// A single trial at 'GrowRequest::firstBeadPosition' with weight one: the caller sampled the
  /// position and accounts for its bias (distance-biased pair and group insertion).
  Fixed = 3,
  /// No first-bead stage: the beads in 'GrowRequest::beadsAlreadyPlaced' keep their positions from
  /// the given molecule and only the remaining beads are (re)grown. Partial reinsertion.
  AlreadyPlaced = 4,
};

/// What to grow. Designated initializers keep call sites readable:
///   CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = old.position}
struct GrowRequest
{
  FirstBeadScheme firstBead{FirstBeadScheme::MultipleFirstBead};
  /// Position of the first bead; required for 'Pinned' and 'Fixed'.
  std::optional<double3> firstBeadPosition{};
  /// Indices (into the component's atoms) of the beads that keep their positions; required for
  /// 'AlreadyPlaced'. Must be a valid placed set of the component's fragment graph.
  std::span<const std::size_t> beadsAlreadyPlaced{};
  /// Molecule id whose atoms in the context's background are ignored: the molecule being regrown,
  /// which is still present in the background but must not interact with its own trial positions.
  /// 'regrowMolecule' sets this to the molecule itself when left empty.
  std::optional<std::size_t> skipBackgroundMolecule{};
};

/// What to retrace; mirrors 'GrowRequest'.
struct RetraceRequest
{
  FirstBeadScheme firstBead{FirstBeadScheme::MultipleFirstBead};
  /// The retained partial first-bead weight of the matching grow ('GrowResult::firstBeadStoredR');
  /// 'Reinsertion' only.
  double storedR{0.0};
  /// See 'GrowRequest::beadsAlreadyPlaced'; 'AlreadyPlaced' only.
  std::span<const std::size_t> beadsAlreadyPlaced{};
  /// See 'GrowRequest::skipBackgroundMolecule'.
  std::optional<std::size_t> skipBackgroundMolecule{};
};

/// Grows a new molecule of 'component' from its reference geometry with the attributes of 'identity'.
/// First-bead schemes: MultipleFirstBead, Pinned, Fixed.
[[nodiscard]] std::optional<GrowResult> growNewMolecule(RandomNumber &random, const GrowContext &context,
                                                        const Component &component,
                                                        const NewMoleculeIdentity &identity,
                                                        const GrowRequest &request = {});

/// Regrows an existing molecule: the identity, charge, and scaling attributes are those of
/// 'moleculeAtoms' and the molecule record ('atomIndex', 'numberOfAtoms') is that of 'molecule'.
/// First-bead schemes: Reinsertion (the whole molecule at a new position), AlreadyPlaced (part of
/// the molecule, the placed beads keep their positions).
[[nodiscard]] std::optional<GrowResult> regrowMolecule(RandomNumber &random, const GrowContext &context,
                                                       const Component &component, const Molecule &molecule,
                                                       std::span<const Atom> moleculeAtoms,
                                                       const GrowRequest &request = {.firstBead =
                                                                                         FirstBeadScheme::Reinsertion});

/// The Rosenbluth weight and energies of an existing molecule, with the first bead weighted by the
/// scheme of the request. The trial directions of the remaining beads are drawn from 'random' for
/// every scheme, so the move stays reproducible under a fixed seed.
[[nodiscard]] RetraceResult retraceMolecule(RandomNumber &random, const GrowContext &context,
                                            const Component &component, std::span<const Atom> moleculeAtoms,
                                            const RetraceRequest &request = {});

/// Dual cut-off correction of a grown molecule: when the force field enables the dual cut-off
/// scheme, the external energy of the result is corrected from the inner cut-off it was grown with
/// to the full cut-offs and the Rosenbluth weight is multiplied by exp(-beta dU), so the result
/// behaves as if grown at the full cut-offs. Returns false when the molecule overlaps at the full
/// cut-offs (the caller rejects the move). A no-op returning true when the scheme is off. The
/// context supplies the background; 'skipBackgroundMolecule' as for the grow.
[[nodiscard]] bool applyDualCutOffCorrection(const GrowContext &context, const Component &component,
                                             GrowResult &result,
                                             std::optional<std::size_t> skipBackgroundMolecule = std::nullopt);

/// The same for a retraced molecule, whose atoms are 'moleculeAtoms'.
[[nodiscard]] bool applyDualCutOffCorrection(const GrowContext &context, const Component &component,
                                             std::span<const Atom> moleculeAtoms, RetraceResult &result,
                                             std::optional<std::size_t> skipBackgroundMolecule = std::nullopt);
}  // namespace CBMC
