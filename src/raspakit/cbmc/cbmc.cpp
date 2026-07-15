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
import cbmc_first_bead_data;
import cbmc_multiple_first_bead;
import cbmc_rigid_insertion;
import cbmc_rigid_deletion;
import cbmc_flexible_insertion;
import cbmc_flexible_deletion;
import cbmc_recoil_growth;
import framework;
import component;
import cbmc_chain_data;
import interpolation_energy_grid;

// Insertion:
// Insertion means growing a new molecule, and therefore the atrributes of the atoms are
// taken from 'component.atoms'. The parameters 'scaling', 'groupId',
// and 'isFractional' are passed to the function and set on the atoms.

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeSwapInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    Component::GrowType growType, std::size_t selectedMolecule, double scaling, std::uint8_t groupId,
    bool isFractional) noexcept
{
  std::size_t startingBead = component.startingBead;
  Atom firstBead = component.atoms[startingBead];
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData =
      CBMC::growMoleculeMultipleFirstBeadSwapInsertion(random, context, component, firstBead);

  if (!firstBeadData) return std::nullopt;

  if (component.atoms.size() == 1)
  {
    return ChainGrowData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass, selectedComponent, component.definedAtoms.size()),
                         {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData =
          CBMC::growRigidMoleculeChainInsertion(random, context, component, selectedComponent, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::growRecoilGrowthMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                     {component.startingBead})
                      : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                 {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;

  return ChainGrowData(chainData->molecule, chainData->atoms, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeSwapDeletion(RandomNumber &random, const GrowContext &context,
                                                                 const Component &component,
                                                                 Component::GrowType growType,
                                                                 std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData =
      CBMC::retraceMultipleFirstBeadSwapDeletion(random, context, component, molecule_atoms[startingBead]);

  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, context, component, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                       {component.startingBead})
                      : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                   {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeReinsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    Component::GrowType growType, Molecule &molecule, std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;
  const std::make_signed_t<std::size_t> skipBackgroundMolecule =
      static_cast<std::make_signed_t<std::size_t>>(molecule_atoms[startingBead].moleculeId);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMultipleFirstBeadReinsertion(
      random, context, component, molecule_atoms[startingBead], skipBackgroundMolecule);

  if (!firstBeadData) return std::nullopt;

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
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

  if (molecule_atoms.size() == 1)
  {
    Molecule firstBeadMolecule = Molecule(firstBeadData->atom.position, simd_quatd(), component.totalMass,
                                          selectedComponent, component.definedAtoms.size());
    firstBeadMolecule.atomIndex = molecule.atomIndex;
    firstBeadMolecule.numberOfAtoms = molecule.numberOfAtoms;

    return ChainGrowData(firstBeadMolecule, {firstBeadData->atom}, firstBeadData->energies,
                         firstBeadData->RosenbluthWeight, firstBeadData->storedR);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, context, component, selectedComponent, atoms,
                                                        skipBackgroundMolecule);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::growRecoilGrowthMoleculeChainInsertion(
                            random, context, component, atoms, {component.startingBead}, skipBackgroundMolecule)
                      : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, atoms,
                                                                 {component.startingBead}, skipBackgroundMolecule);
      break;
    default:
      std::unreachable();
  }
  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return ChainGrowData(chainData->molecule, chainData->atoms, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, firstBeadData->storedR);
}

[[nodiscard]] std::optional<ChainRetraceData> CBMC::retraceMoleculeReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, Component::GrowType growType,
    [[maybe_unused]] Molecule &molecule, std::span<Atom> molecule_atoms, double storedR) noexcept
{
  std::size_t startingBead = component.startingBead;

  const std::optional<FirstBeadData> firstBeadData = CBMC::retraceMultipleFirstBeadReinsertion(
      random, context, component, molecule_atoms[startingBead], storedR,
      static_cast<std::make_signed_t<std::size_t>>(molecule_atoms[startingBead].moleculeId));
  if (!firstBeadData)
  {
    return std::nullopt;
  }

  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, context, component, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                       {component.startingBead})
                      : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                   {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData->energies + chainData.energies,
                          firstBeadData->RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculePartialReinsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    Component::GrowType growType, Molecule &molecule, std::span<Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept
{
  const std::make_signed_t<std::size_t> skipBackgroundMolecule =
      static_cast<std::make_signed_t<std::size_t>>(moleculeAtoms.front().moleculeId);

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, context, component, selectedComponent, moleculeAtoms,
                                                        skipBackgroundMolecule);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::growRecoilGrowthMoleculeChainInsertion(
                            random, context, component, moleculeAtoms, beadsAlreadyPlaced, skipBackgroundMolecule)
                      : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, moleculeAtoms,
                                                                 beadsAlreadyPlaced, skipBackgroundMolecule);
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return ChainGrowData(chainData->molecule, chainData->atoms, chainData->energies, chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculePartialReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, Component::GrowType growType,
    [[maybe_unused]] Molecule &molecule, std::span<Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept
{
  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, context, component, moleculeAtoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, moleculeAtoms,
                                                                       beadsAlreadyPlaced)
                      : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, moleculeAtoms,
                                                                   beadsAlreadyPlaced);
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(chainData.energies, chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeIdentityChangeInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    Component::GrowType growType, std::size_t selectedMolecule, const Atom &oldStartingBead, double scaling,
    std::uint8_t groupId, bool isFractional, std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::size_t startingBead = component.startingBead;
  Atom firstBead = component.atoms[startingBead];
  firstBead.position = oldStartingBead.position;
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData =
      CBMC::growMultipleFirstBeadPartialInsertion(context, component, firstBead, skipBackgroundMolecule);

  if (!firstBeadData) return std::nullopt;

  if (component.atoms.size() == 1)
  {
    // update atom index
    return ChainGrowData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass,selectedComponent, component.definedAtoms.size()),
                         {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, context, component, selectedComponent, molecule_atoms,
                                                        skipBackgroundMolecule);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::growRecoilGrowthMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                     {component.startingBead}, skipBackgroundMolecule)
                      : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                 {component.startingBead}, skipBackgroundMolecule);
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;
  // update atom index

  return ChainGrowData(chainData->molecule, chainData->atoms, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeIdentityChangeDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, Component::GrowType growType,
    std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData =
      CBMC::retraceMultipleFirstBeadPartialDeletion(context, component, molecule_atoms[startingBead]);

  if (molecule_atoms.size() == 1)
  {
    // update atom index
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, context, component, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                       {component.startingBead})
                      : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                   {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculePairSecondSwapInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
    Component::GrowType growType, std::size_t selectedMolecule, double3 fixedFirstBeadPosition, double scaling,
    std::uint8_t groupId, bool isFractional) noexcept
{
  std::size_t startingBead = component.startingBead;
  Atom firstBead = component.atoms[startingBead];
  firstBead.position = fixedFirstBeadPosition;
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData =
      CBMC::growFirstBeadAtFixedPosition(context, component, firstBead);

  if (!firstBeadData) return std::nullopt;

  if (component.atoms.size() == 1)
  {
    return ChainGrowData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass, selectedComponent, component.definedAtoms.size()),
                         {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData =
          CBMC::growRigidMoleculeChainInsertion(random, context, component, selectedComponent, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::growRecoilGrowthMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                     {component.startingBead})
                      : CBMC::growFlexibleMoleculeChainInsertion(random, context, component, molecule_atoms,
                                                                 {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;

  return ChainGrowData(chainData->molecule, chainData->atoms, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculePairSecondSwapDeletion(const GrowContext &context,
                                                                           const Component &component,
                                                                           Component::GrowType growType,
                                                                           std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData =
      CBMC::retraceFirstBeadAtFixedPosition(context, component, molecule_atoms[startingBead]);

  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  RandomNumber random{std::nullopt};
  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, context, component, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = context.forceField.useRecoilGrowth
                      ? CBMC::retraceRecoilGrowthMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                       {component.startingBead})
                      : CBMC::retraceFlexibleMoleculeChainDeletion(random, context, component, molecule_atoms,
                                                                   {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}
