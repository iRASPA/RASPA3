module;

module cbmc;

import std;

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import energy_status;
import forcefield;
import cbmc_results;
import cbmc_first_bead;
import cbmc_chain_cbmc;
import cbmc_chain_recoil;
import framework;
import interpolation_energy_grid;

// The CBMC entry points share one shape: sample (or retrace) the first bead with the scheme of the
// move, then grow (or retrace) the remaining beads with the fragment-at-a-time operator engine, and
// combine the two Rosenbluth weights and energies. The helpers below hold that shared shape; each
// entry point only supplies its first-bead scheme.

// Every multi-atom molecule is grown with the fragment-at-a-time operator engine: a fully rigid
// molecule is a single rigid seed fragment (placed with uniform random orientations), a flexible or
// semi-flexible molecule is grown fragment by fragment. Single-atom molecules are handled before
// dispatch. The chain scheme is the context's ('GrowContext::chainScheme'), so a caller that needs a
// specific scheme (Widom, the ideal-gas reference grows) selects it there rather than by editing the
// force field.
static std::optional<ChainGrowData> growChain(RandomNumber &random, const CBMC::GrowContext &context,
                                              Component &component, std::span<Atom> molecule_atoms,
                                              const std::vector<std::size_t> &beadsAlreadyPlaced,
                                              std::optional<std::size_t> skipBackgroundMolecule = std::nullopt)
{
  return context.chainScheme == CBMC::ChainScheme::RecoilGrowth
             ? CBMC::growRecoilGrowthMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                            beadsAlreadyPlaced, skipBackgroundMolecule)
             : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                        beadsAlreadyPlaced, skipBackgroundMolecule);
}

static ChainRetraceData retraceChain(RandomNumber &random, const CBMC::GrowContext &context,
                                     const Component &component, std::span<Atom> molecule_atoms,
                                     const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  return context.chainScheme == CBMC::ChainScheme::RecoilGrowth
             ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                              beadsAlreadyPlaced)
             : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                          beadsAlreadyPlaced);
}

// First bead of a freshly inserted molecule: the reference starting bead carrying the identity and
// scaling attributes of the new molecule, optionally pinned at a given position (identity change,
// distance-biased pair insertion).
static Atom makeFirstBead(const Component &component, std::size_t selectedMolecule, double scaling,
                          std::uint8_t groupId, bool isFractional, std::optional<double3> position = std::nullopt)
{
  Atom firstBead = component.atoms[component.startingBead];
  if (position.has_value()) firstBead.position = position.value();
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);
  return firstBead;
}

// Combined result of the first-bead stage and the chain stage: energies add, Rosenbluth weights
// multiply (both stages carry their weight as a logarithm, so the logs add).
static ChainGrowData combineGrowData(const FirstBeadData &firstBeadData, const ChainGrowData &chainData,
                                     double storedR)
{
  return ChainGrowData(chainData.molecule, chainData.atoms, firstBeadData.energies + chainData.energies,
                       firstBeadData.logRosenbluthWeight + chainData.logRosenbluthWeight, storedR);
}

// Grows the remainder of a freshly inserted molecule after its first bead was sampled: a single-atom
// molecule is already complete; otherwise the reference geometry is translated so the starting bead
// sits at the sampled first-bead position, the identity and scaling attributes are applied to every
// atom, and the remaining beads are grown with the operator engine.
static std::optional<ChainGrowData> growNewMoleculeAtFirstBead(
    RandomNumber &random, const CBMC::GrowContext &context, Component &component, std::size_t selectedComponent,
    std::size_t selectedMolecule, double scaling, std::uint8_t groupId, bool isFractional,
    const FirstBeadData &firstBeadData, std::optional<std::size_t> skipBackgroundMolecule = std::nullopt)
{
  if (component.atoms.size() == 1)
  {
    return ChainGrowData(Molecule(double3(firstBeadData.atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass, selectedComponent, component.definedAtoms.size()),
                         {firstBeadData.atom}, firstBeadData.energies, firstBeadData.logRosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData.atom.position'
  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData.atom.position - component.atoms[component.startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData =
      growChain(random, context, component, molecule_atoms, {component.startingBead}, skipBackgroundMolecule);
  if (!chainData) return std::nullopt;

  return combineGrowData(firstBeadData, *chainData, 0.0);
}

// Retraces the remainder of an existing molecule after its first bead was retraced: a single-atom
// molecule is already complete; otherwise the remaining beads are retraced with the operator engine.
static ChainRetraceData retraceAfterFirstBead(RandomNumber &random, const CBMC::GrowContext &context,
                                              const Component &component, std::span<Atom> molecule_atoms,
                                              const FirstBeadData &firstBeadData)
{
  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData.energies, firstBeadData.logRosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData = retraceChain(random, context, component, molecule_atoms, {component.startingBead});

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.logRosenbluthWeight + chainData.logRosenbluthWeight, 0.0);
}

// Insertion:
// Insertion means growing a new molecule, and therefore the attributes of the atoms are
// taken from 'component.atoms'. The parameters 'scaling', 'groupId',
// and 'isFractional' are passed to the function and set on the atoms.

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeSwapInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    std::size_t selectedMolecule, double scaling, std::uint8_t groupId, bool isFractional)
{
  Atom firstBead = makeFirstBead(component, selectedMolecule, scaling, groupId, isFractional);

  std::optional<FirstBeadData> const firstBeadData =
      CBMC::growMoleculeMultipleFirstBeadSwapInsertion(random, context, component, firstBead);
  if (!firstBeadData) return std::nullopt;

  return growNewMoleculeAtFirstBead(random, context, component, selectedComponent, selectedMolecule, scaling, groupId,
                                    isFractional, *firstBeadData);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeSwapDeletion(RandomNumber &random, const GrowContext &context,
                                                                 const Component &component,
                                                                 std::span<Atom> molecule_atoms)
{
  const FirstBeadData firstBeadData = CBMC::retraceMultipleFirstBeadSwapDeletion(
      random, context, component, molecule_atoms[component.startingBead]);

  return retraceAfterFirstBead(random, context, component, molecule_atoms, firstBeadData);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeReinsertion(RandomNumber &random,
                                                                         const GrowContext &context,
                                                                         Component &component,
                                                                         std::size_t selectedComponent,
                                                                         Molecule &molecule,
                                                                         std::span<Atom> molecule_atoms)
{
  std::size_t startingBead = component.startingBead;
  // The molecule is regrown against a background that still contains its old copy: skip it.
  const std::optional<std::size_t> skipBackgroundMolecule = molecule_atoms[startingBead].moleculeId;

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMultipleFirstBeadReinsertion(
      random, context, component, molecule_atoms[startingBead], skipBackgroundMolecule);

  if (!firstBeadData) return std::nullopt;

  if (molecule_atoms.size() == 1)
  {
    Molecule firstBeadMolecule = Molecule(firstBeadData->atom.position, simd_quatd(), component.totalMass,
                                          selectedComponent, component.definedAtoms.size());
    firstBeadMolecule.atomIndex = molecule.atomIndex;
    firstBeadMolecule.numberOfAtoms = molecule.numberOfAtoms;

    return ChainGrowData(firstBeadMolecule, {firstBeadData->atom}, firstBeadData->energies,
                         firstBeadData->logRosenbluthWeight, firstBeadData->storedR);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'; the
  // identity and scaling attributes are copied from the old molecule, since we reinsert the same
  // molecule.
  std::vector<Atom> atoms = component.atoms;
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    atoms[i].position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atoms[i].charge = molecule_atoms[i].charge;
    atoms[i].scalingVDW = molecule_atoms[i].scalingVDW;
    atoms[i].scalingCoulomb = molecule_atoms[i].scalingCoulomb;
    atoms[i].moleculeId = molecule_atoms[i].moleculeId;
    atoms[i].componentId = molecule_atoms[i].componentId;
    atoms[i].groupId = molecule_atoms[i].groupId;
    atoms[i].isFractional = molecule_atoms[i].isFractional;
  }

  std::optional<ChainGrowData> chainData =
      growChain(random, context, component, atoms, {component.startingBead}, skipBackgroundMolecule);
  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return combineGrowData(*firstBeadData, *chainData, firstBeadData->storedR);
}

[[nodiscard]] std::optional<ChainRetraceData> CBMC::retraceMoleculeReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, [[maybe_unused]] Molecule &molecule,
    std::span<Atom> molecule_atoms, double storedR)
{
  const std::optional<FirstBeadData> firstBeadData = CBMC::retraceMultipleFirstBeadReinsertion(
      random, context, component, molecule_atoms[component.startingBead], storedR,
      std::optional<std::size_t>{molecule_atoms[component.startingBead].moleculeId});
  if (!firstBeadData) return std::nullopt;

  return retraceAfterFirstBead(random, context, component, molecule_atoms, *firstBeadData);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculePartialReinsertion(
    RandomNumber &random, const GrowContext &context, Component &component,
    [[maybe_unused]] std::size_t selectedComponent, Molecule &molecule, std::span<Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  const std::optional<std::size_t> skipBackgroundMolecule = moleculeAtoms.front().moleculeId;

  std::optional<ChainGrowData> chainData =
      growChain(random, context, component, moleculeAtoms, beadsAlreadyPlaced, skipBackgroundMolecule);

  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return ChainGrowData(chainData->molecule, chainData->atoms, chainData->energies, chainData->logRosenbluthWeight,
                       0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculePartialReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, [[maybe_unused]] Molecule &molecule,
    std::span<Atom> moleculeAtoms, const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  ChainRetraceData chainData = retraceChain(random, context, component, moleculeAtoms, beadsAlreadyPlaced);

  return ChainRetraceData(chainData.energies, chainData.logRosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeIdentityChangeInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    std::size_t selectedMolecule, const Atom &oldStartingBead, double scaling, std::uint8_t groupId, bool isFractional,
    std::optional<std::size_t> skipBackgroundMolecule)
{
  Atom firstBead =
      makeFirstBead(component, selectedMolecule, scaling, groupId, isFractional, oldStartingBead.position);

  std::optional<FirstBeadData> const firstBeadData =
      CBMC::growMultipleFirstBeadPartialInsertion(context, component, firstBead, skipBackgroundMolecule);
  if (!firstBeadData) return std::nullopt;

  return growNewMoleculeAtFirstBead(random, context, component, selectedComponent, selectedMolecule, scaling, groupId,
                                    isFractional, *firstBeadData, skipBackgroundMolecule);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeIdentityChangeDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component,
    std::span<Atom> molecule_atoms)
{
  const FirstBeadData firstBeadData =
      CBMC::retraceMultipleFirstBeadPartialDeletion(context, component, molecule_atoms[component.startingBead]);

  return retraceAfterFirstBead(random, context, component, molecule_atoms, firstBeadData);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculePairSecondSwapInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    std::size_t selectedMolecule, double3 fixedFirstBeadPosition, double scaling, std::uint8_t groupId,
    bool isFractional)
{
  Atom firstBead =
      makeFirstBead(component, selectedMolecule, scaling, groupId, isFractional, fixedFirstBeadPosition);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growFirstBeadAtFixedPosition(context, component, firstBead);
  if (!firstBeadData) return std::nullopt;

  return growNewMoleculeAtFirstBead(random, context, component, selectedComponent, selectedMolecule, scaling, groupId,
                                    isFractional, *firstBeadData);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculePairSecondSwapDeletion(RandomNumber &random,
                                                                           const GrowContext &context,
                                                                           const Component &component,
                                                                           std::span<Atom> molecule_atoms)
{
  // The first bead is pinned (weight one), but the chain retrace draws the trial directions of the
  // remaining beads from the caller's generator, like every other retrace: the move stays
  // reproducible under a fixed seed.
  const FirstBeadData firstBeadData =
      CBMC::retraceFirstBeadAtFixedPosition(context, component, molecule_atoms[component.startingBead]);

  return retraceAfterFirstBead(random, context, component, molecule_atoms, firstBeadData);
}
