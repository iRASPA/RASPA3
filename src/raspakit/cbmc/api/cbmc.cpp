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
                                    const std::vector<std::size_t> &beadsAlreadyPlaced,
                                    std::optional<std::size_t> skipBackgroundMolecule)
{
  return context.settings.chainScheme == ChainScheme::RecoilGrowth
             ? growChainRecoil(random, context, component, moleculeAtoms, beadsAlreadyPlaced,
                                                      skipBackgroundMolecule)
             : growChainCBMC(random, context, component, moleculeAtoms, beadsAlreadyPlaced,
                                                  skipBackgroundMolecule);
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
                                           const Atom &firstBead, FirstBeadScheme scheme,
                                           std::optional<std::size_t> skipBackgroundMolecule)
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return growMultipleFirstBead(random, context, component, firstBead);
    case FirstBeadScheme::Reinsertion:
      return growMultipleFirstBeadReinsertion(random, context, component, firstBead, skipBackgroundMolecule);
    case FirstBeadScheme::Pinned:
      return growPinnedFirstBead(context, component, firstBead, skipBackgroundMolecule);
    case FirstBeadScheme::Fixed:
      return growFixedFirstBead(context, component, firstBead, skipBackgroundMolecule);
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
      return retraceMultipleFirstBeadReinsertion(context, component, firstBead, request.storedR,
                                                 request.skipBackgroundMolecule);
    case FirstBeadScheme::Pinned:
      return retracePinnedFirstBead(context, component, firstBead, request.skipBackgroundMolecule);
    case FirstBeadScheme::Fixed:
      return retraceFixedFirstBead(context, component, firstBead, request.skipBackgroundMolecule);
    case FirstBeadScheme::AlreadyPlaced:
      break;
  }
  throw std::invalid_argument("CBMC: no first-bead stage for FirstBeadScheme::AlreadyPlaced");
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
                                             const FirstBeadData &firstBeadData,
                                             std::optional<std::size_t> skipBackgroundMolecule)
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

  std::optional<GrowResult> chainResult =
      growChain(random, context, component, atoms, {component.startingBead}, skipBackgroundMolecule);
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
  if (request.firstBead == FirstBeadScheme::AlreadyPlaced)
  {
    throw std::invalid_argument(
        "CBMC::growNewMolecule: a new molecule has no placed beads; use regrowMolecule for AlreadyPlaced");
  }

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
      growFirstBead(random, context, component, firstBead, request.firstBead, request.skipBackgroundMolecule);
  if (!firstBeadData) return std::nullopt;

  return growAfterFirstBead(random, context, component, templateAtoms, *firstBeadData,
                            request.skipBackgroundMolecule);
}

std::optional<CBMC::GrowResult> CBMC::regrowMolecule(RandomNumber &random, const GrowContext &context,
                                                     const Component &component, const Molecule &molecule,
                                                     std::span<const Atom> moleculeAtoms, const GrowRequest &request)
{
  // The molecule is regrown against a background that still contains its old copy: skip it.
  const std::optional<std::size_t> skipBackgroundMolecule =
      request.skipBackgroundMolecule.has_value()
          ? request.skipBackgroundMolecule
          : std::optional<std::size_t>{moleculeAtoms[component.startingBead].moleculeId};

  std::optional<GrowResult> result;
  if (request.firstBead == FirstBeadScheme::AlreadyPlaced)
  {
    result = growChain(random, context, component, moleculeAtoms,
                       placedSetOf(request.beadsAlreadyPlaced, "regrowMolecule"), skipBackgroundMolecule);
  }
  else
  {
    Atom firstBead = moleculeAtoms[component.startingBead];
    if (request.firstBead == FirstBeadScheme::Pinned || request.firstBead == FirstBeadScheme::Fixed)
    {
      firstBead.position = requiredPosition(request.firstBeadPosition, "regrowMolecule");
    }

    const std::optional<FirstBeadData> firstBeadData =
        growFirstBead(random, context, component, firstBead, request.firstBead, skipBackgroundMolecule);
    if (!firstBeadData) return std::nullopt;

    result = growAfterFirstBead(random, context, component, moleculeAtoms, *firstBeadData, skipBackgroundMolecule);
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

bool CBMC::applyDualCutOffCorrection(const GrowContext &context, const Component &component, GrowResult &result,
                                     std::optional<std::size_t> skipBackgroundMolecule)
{
  if (!context.forceField.useDualCutOff) return true;

  const std::optional<RunningEnergy> correction =
      computeDualCutOffCorrection(context, component, result.atoms, skipBackgroundMolecule);
  if (!correction.has_value()) return false;

  result.energies += correction.value();
  result.multiplyRosenbluthWeight(-context.beta * correction->potentialEnergy());
  return true;
}

bool CBMC::applyDualCutOffCorrection(const GrowContext &context, const Component &component,
                                     std::span<const Atom> moleculeAtoms, RetraceResult &result,
                                     std::optional<std::size_t> skipBackgroundMolecule)
{
  if (!context.forceField.useDualCutOff) return true;

  const std::optional<RunningEnergy> correction =
      computeDualCutOffCorrection(context, component, moleculeAtoms, skipBackgroundMolecule);
  if (!correction.has_value()) return false;

  result.energies += correction.value();
  result.multiplyRosenbluthWeight(-context.beta * correction->potentialEnergy());
  return true;
}

double CBMC::logBaseSamplerNormalization(double beta, const Component &component, const std::vector<GrowStep> &plan)
{
  return logFlexibleBaseNormalization(beta, component, plan);
}
