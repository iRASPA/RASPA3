module;

module cbmc;

import std;

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import running_energy;
import cbmc_results;
import cbmc_grow_context;
import cbmc_first_bead;
import cbmc_chain_cbmc;
import cbmc_chain_recoil;
import cbmc_external_energy;
import cbmc_flexible_base;
import cbmc_grow_step;

// The entry points share one shape: place (or retrace) the first bead with the scheme of the request,
// then grow (or retrace) the remaining beads with the fragment-at-a-time operator engine, and combine
// the two Rosenbluth weights and energies.
//
// Every multi-atom molecule is grown with the operator engine: a fully rigid molecule is a single
// rigid seed fragment (placed with uniform random orientations), a flexible or semi-flexible molecule
// is grown fragment by fragment. Single-atom molecules are complete after the first bead. The chain
// scheme is the context's ('GrowContext::settings.chainScheme'), so a caller that needs a specific
// scheme (Widom, the ideal-gas reference grows) selects it there rather than by editing the force
// field.

namespace
{
using namespace CBMC;

std::optional<GrowResult> growChain(RandomNumber &random, const GrowContext &context, const Component &component,
                                    std::span<const Atom> moleculeAtoms,
                                    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  return context.settings.chainScheme == ChainScheme::RecoilGrowth
             ? growChainRecoil(random, context, component, moleculeAtoms, beadsAlreadyPlaced)
             : growChainCBMC(random, context, component, moleculeAtoms, beadsAlreadyPlaced);
}

RetraceResult retraceChain(RandomNumber &random, const GrowContext &context, const Component &component,
                           std::span<const Atom> moleculeAtoms, const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  return context.settings.chainScheme == ChainScheme::RecoilGrowth
             ? retraceChainRecoil(random, context, component, moleculeAtoms, beadsAlreadyPlaced)
             : retraceChainCBMC(random, context, component, moleculeAtoms, beadsAlreadyPlaced);
}

std::vector<std::size_t> placedSetOf(std::span<const std::size_t> beadsAlreadyPlaced, const char *entryPoint)
{
  if (beadsAlreadyPlaced.empty())
  {
    throw std::invalid_argument(
        std::format("CBMC::{}: FirstBeadScheme::AlreadyPlaced requires a non-empty 'beadsAlreadyPlaced'", entryPoint));
  }
  return std::vector<std::size_t>(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());
}

double3 requiredPosition(const std::optional<double3> &position, const char *entryPoint)
{
  if (!position.has_value())
  {
    throw std::invalid_argument(
        std::format("CBMC::{}: FirstBeadScheme::Pinned and ::Fixed require 'firstBeadPosition'", entryPoint));
  }
  return position.value();
}

/// First-bead stage of a grow; 'firstBead' carries the identity and scaling attributes of the molecule
/// and, for the pinned and fixed schemes, the position.
std::optional<FirstBeadData> growFirstBead(RandomNumber &random, const GrowContext &context, const Component &component,
                                           const Atom &firstBead, FirstBeadScheme scheme)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return growMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return growMultipleFirstBeadReinsertion(random, context, component, firstBead);
    case FirstBeadScheme::Pinned:
      return growPinnedFirstBead(context, component, firstBead);
    case FirstBeadScheme::Fixed:
      return growFixedFirstBead(context, component, firstBead);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
}

FirstBeadData retraceFirstBead(RandomNumber &random, const GrowContext &context, const Component &component,
                               const Atom &firstBead, const RetraceRequest &request)
{
  switch (request.firstBead)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return retraceMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return retraceMultipleFirstBeadReinsertion(context, component, firstBead, request.storedR);
    case FirstBeadScheme::Pinned:
      return retracePinnedFirstBead(context, component, firstBead);
    case FirstBeadScheme::Fixed:
      return retraceFixedFirstBead(context, component, firstBead);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
}

/// The name of a scheme for the error messages of the entry points.
const char *schemeName(FirstBeadScheme scheme)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return "MultipleFirstBead";
    case FirstBeadScheme::Reinsertion:
      return "Reinsertion";
    case FirstBeadScheme::Pinned:
      return "Pinned";
    case FirstBeadScheme::Fixed:
      return "Fixed";
    case FirstBeadScheme::AlreadyPlaced:
      return "AlreadyPlaced";
  }
  return "?";
}

/// Rejects a first-bead scheme an entry point does not support.
void requireScheme(FirstBeadScheme scheme, std::initializer_list<FirstBeadScheme> allowed, const char *entryPoint,
                   const char *allowedText)
{
  if (std::find(allowed.begin(), allowed.end(), scheme) != allowed.end()) return;
  throw std::invalid_argument(std::format("CBMC::{}: FirstBeadScheme::{} is not supported here (allowed: {})",
                                          entryPoint, schemeName(scheme), allowedText));
}

/// Combined result of the first-bead stage and the chain stage: energies add, Rosenbluth weights
/// multiply (both stages carry their weight as a logarithm, so the logs add).
GrowResult combine(const FirstBeadData &firstBeadData, GrowResult chainResult)
{
  chainResult.energies += firstBeadData.energies;
  chainResult.logRosenbluthWeight += firstBeadData.logRosenbluthWeight;
  chainResult.firstBeadStoredR = firstBeadData.storedR;
  return chainResult;
}

/// Grows the remainder of a molecule after its first bead was placed. 'templateAtoms' supplies the
/// identity, charge, and scaling attributes of every atom (the reference atoms stamped with the new
/// identity for a new molecule, the old atoms for a regrow); the reference geometry is translated so
/// the starting bead sits at the sampled first-bead position.
std::optional<GrowResult> growAfterFirstBead(RandomNumber &random, const GrowContext &context,
                                             const Component &component, std::span<const Atom> templateAtoms,
                                             const FirstBeadData &firstBeadData)
{
  if (component.atoms.size() == 1)
  {
    const std::vector<Atom> atoms{firstBeadData.atom};
    return GrowResult(component.createMoleculeRecord(atoms), atoms, firstBeadData.energies,
                      firstBeadData.logRosenbluthWeight, firstBeadData.storedR);
  }

  const double3 shift = firstBeadData.atom.position - component.atoms[component.startingBead].position;
  std::vector<Atom> atoms(component.atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    atoms[i] = templateAtoms[i];
    atoms[i].position = component.atoms[i].position + shift;
  }

  std::optional<GrowResult> chainResult = growChain(random, context, component, atoms, {component.startingBead});
  if (!chainResult) return std::nullopt;

  return combine(firstBeadData, std::move(*chainResult));
}

RetraceResult retraceAfterFirstBead(RandomNumber &random, const GrowContext &context, const Component &component,
                                    std::span<const Atom> moleculeAtoms, const FirstBeadData &firstBeadData)
{
  if (moleculeAtoms.size() == 1)
  {
    return RetraceResult(firstBeadData.energies, firstBeadData.logRosenbluthWeight);
  }

  RetraceResult chainResult = retraceChain(random, context, component, moleculeAtoms, {component.startingBead});
  chainResult.energies += firstBeadData.energies;
  chainResult.logRosenbluthWeight += firstBeadData.logRosenbluthWeight;
  return chainResult;
}
}  // namespace

std::optional<CBMC::GrowResult> CBMC::growNewMolecule(RandomNumber &random, const GrowContext &context,
                                                      const Component &component,
                                                      const NewMoleculeIdentity &identity, const GrowRequest &request)
{
  // A new molecule has no placed beads (AlreadyPlaced) and no old copy to reinsert (Reinsertion).
  requireScheme(request.firstBead, {FirstBeadScheme::MultipleFirstBead, FirstBeadScheme::Pinned, FirstBeadScheme::Fixed},
                "growNewMolecule", "MultipleFirstBead, Pinned, Fixed");

  // The reference atoms stamped with the identity of the new molecule.
  std::vector<Atom> templateAtoms = component.atoms;
  for (Atom &atom : templateAtoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(identity.moleculeId);
    atom.groupId = identity.groupId;
    atom.isFractional = identity.isFractional;
    atom.setScaling(identity.scaling);
  }

  Atom firstBead = templateAtoms[component.startingBead];
  if (request.firstBead == FirstBeadScheme::Pinned || request.firstBead == FirstBeadScheme::Fixed)
  {
    firstBead.position = requiredPosition(request.firstBeadPosition, "growNewMolecule");
  }

  const std::optional<FirstBeadData> firstBeadData =
      growFirstBead(random, context, component, firstBead, request.firstBead);
  if (!firstBeadData) return std::nullopt;

  return growAfterFirstBead(random, context, component, templateAtoms, *firstBeadData);
}

std::optional<CBMC::GrowResult> CBMC::regrowMolecule(RandomNumber &random, const GrowContext &context,
                                                     const Component &component, const Molecule &molecule,
                                                     std::span<const Atom> moleculeAtoms, const GrowRequest &request)
{
  // The old copy of the molecule is still in the background; it is excluded through its molecule id,
  // which the regrown atoms share (see GrowContext).
  requireScheme(request.firstBead, {FirstBeadScheme::Reinsertion, FirstBeadScheme::AlreadyPlaced}, "regrowMolecule",
                "Reinsertion, AlreadyPlaced");

  std::optional<GrowResult> result;
  if (request.firstBead == FirstBeadScheme::AlreadyPlaced)
  {
    result = growChain(random, context, component, moleculeAtoms,
                       placedSetOf(request.beadsAlreadyPlaced, "regrowMolecule"));
  }
  else
  {
    const std::optional<FirstBeadData> firstBeadData =
        growFirstBead(random, context, component, moleculeAtoms[component.startingBead], request.firstBead);
    if (!firstBeadData) return std::nullopt;

    result = growAfterFirstBead(random, context, component, moleculeAtoms, *firstBeadData);
  }
  if (!result) return std::nullopt;

  // The same molecule is reinserted: it keeps its record in the system.
  result->molecule.atomIndex = molecule.atomIndex;
  result->molecule.numberOfAtoms = molecule.numberOfAtoms;
  return result;
}

CBMC::RetraceResult CBMC::retraceMolecule(RandomNumber &random, const GrowContext &context, const Component &component,
                                          std::span<const Atom> moleculeAtoms, const RetraceRequest &request)
{
  if (request.firstBead == FirstBeadScheme::AlreadyPlaced)
  {
    return retraceChain(random, context, component, moleculeAtoms,
                        placedSetOf(request.beadsAlreadyPlaced, "retraceMolecule"));
  }

  const FirstBeadData firstBeadData =
      retraceFirstBead(random, context, component, moleculeAtoms[component.startingBead], request);

  return retraceAfterFirstBead(random, context, component, moleculeAtoms, firstBeadData);
}

namespace
{
/// The dual cut-off correction folded into a result's energies and log weight (grow and retrace
/// results carry both under the same names).
template <typename Result>
bool applyCorrection(const CBMC::GrowContext &context, const Component &component, std::span<const Atom> atoms,
                     Result &result)
{
  if (!context.forceField.useDualCutOff) return true;

  const std::optional<RunningEnergy> correction = CBMC::computeDualCutOffCorrection(context, component, atoms);
  if (!correction.has_value()) return false;

  result.energies += correction.value();
  result.multiplyRosenbluthWeight(-context.beta * correction->potentialEnergy());
  return true;
}
}  // namespace

bool CBMC::applyDualCutOffCorrection(const GrowContext &context, const Component &component, GrowResult &result)
{
  return applyCorrection(context, component, result.atoms, result);
}

bool CBMC::applyDualCutOffCorrection(const GrowContext &context, const Component &component,
                                     std::span<const Atom> moleculeAtoms, RetraceResult &result)
{
  return applyCorrection(context, component, moleculeAtoms, result);
}

double CBMC::logBaseSamplerNormalization(double beta, const Component &component, const std::vector<GrowStep> &plan)
{
  return logFlexibleBaseNormalization(beta, component, plan);
}
