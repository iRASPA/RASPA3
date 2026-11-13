module;

module system;

import std;

import randomnumbers;
import simd_quatd;
import double3;
import double3x3;
import double4;
import atom;
import atom_dynamics;
import component;
import blocking_pockets;
import molecule;
import framework;
import forcefield;
import simulationbox;
import units;
import cbmc;
import cbmc_constants;
import cbmc_external_energy;
import interpolation_energy_grid;

// System molecules: insertion/deletion, initialization, and geometry helpers.

std::optional<double> System::frameworkMass() const
{
  if (!framework.has_value()) return std::nullopt;

  double mass = framework->mass;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].type == Component::Type::Cation)
    {
      mass += components[i].totalMass * static_cast<double>(numberOfIntegerMoleculesPerComponent[i]);
    }
  }
  return mass;
}

void System::insertFractionalMoleculeAtIndex(std::size_t selectedComponent, std::size_t moleculeIndex,
                                             [[maybe_unused]] const Molecule& molecule, std::vector<Atom> atoms)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  // Keep the cross-link sites of the molecules that shift up by one index consistent.
  crossLinks.moleculeInserted(selectedComponent, moleculeIndex);

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, moleculeIndex);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});
  insertAtomTags(static_cast<std::size_t>(atomOffset), atoms.size());

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, moleculeIndex);
  moleculeData.insert(moleculeIterator, molecule);

  // Keep electricField aligned with atomData: insert at the molecule position, not at the end.
  std::vector<double3>::const_iterator iterator_electric_field =
      iteratorForElectricField(selectedComponent, moleculeIndex);
  electricField.insert(iterator_electric_field, atoms.size(), double3(0.0, 0.0, 0.0));

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfFractionalMoleculesPerComponent[selectedComponent] += 1;

  // a fractional molecule has velocities like any other: every one of them counts in the temperature
  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // the Coulomb-scaled charge of the inserted atoms (a fractional molecule counts with its scaled charge)
  const double insertedCharge = scaledNetCharge(atoms);
  netCharge += insertedCharge;
  netChargeAdsorbates += insertedCharge;
  netChargePerComponent[selectedComponent] += insertedCharge;

  for (Atom& atom : atoms)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
    addAtomToTailCorrectionCounts(atom);
  }

  updateMoleculeAtomInformation(selectedComponent);
}

void System::insertFractionalMolecule(std::size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                                      std::vector<Atom> atoms, std::size_t moleculeIndex)
{
  const double l = 0.0;
  const std::uint8_t groupId = fractionalSlotDUdlambdaGroupId(selectedComponent, moleculeIndex);
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(moleculeIndex);
    atom.setScalingToFractional(l, groupId);
  }
  insertFractionalMoleculeAtIndex(selectedComponent, moleculeIndex, molecule, std::move(atoms));
}

void System::insertReactionFractionalMolecule(std::size_t selectedComponent, std::size_t moleculeIndex,
                                              [[maybe_unused]] const Molecule& molecule, std::vector<Atom> atoms,
                                              bool isReactant, double lambda, std::uint8_t dUdlambdaGroupId)
{
  // staged schedule from scaling.ixx; reactants are coupled with (1 - lambda)
  const double effectiveLambda = isReactant ? (1.0 - lambda) : lambda;
  for (Atom& atom : atoms)
  {
    atom.setScalingToFractional(effectiveLambda, dUdlambdaGroupId);
  }

  insertFractionalMoleculeAtIndex(selectedComponent, moleculeIndex, molecule, std::move(atoms));
}

/// Inserts a molecule into the vector of atoms.
///
/// Note: updates the numberOfMoleculesPerComponent, numberOfIntegerMoleculesPerComponent,
///       numberOfPseudoAtoms, totalNumberOfPseudoAtoms.
/// - Parameters:
///   - selectedComponent: the index of the component
///   - atoms: vector of atoms to be inserted
/// - returns:
void System::insertMolecule(std::size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                            std::vector<Atom> atoms)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  std::vector<Atom>::const_iterator iterator =
      iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});
  insertAtomTags(static_cast<std::size_t>(atomOffset), atoms.size());

  std::vector<Molecule>::iterator moleculeIterator =
      indexForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  moleculeData.insert(moleculeIterator, molecule);

  // Keep electricField aligned with atomData: insert at the molecule position, not at the end.
  std::vector<double3>::const_iterator iterator_electric_field =
      iteratorForElectricField(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  electricField.insert(iterator_electric_field, atoms.size(), double3(0.0, 0.0, 0.0));

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

  // the Coulomb-scaled charge of the inserted atoms (a fractional molecule counts with its scaled charge)
  const double insertedCharge = scaledNetCharge(atoms);
  netCharge += insertedCharge;
  netChargeAdsorbates += insertedCharge;
  netChargePerComponent[selectedComponent] += insertedCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(numberOfMoleculesPerComponent[selectedComponent]);
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
    addAtomToTailCorrectionCounts(atom);
  }

  updateMoleculeAtomInformation(selectedComponent);
}

void System::insertMoleculePolarization(std::size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                                        std::vector<Atom> atoms, std::span<double3> electric_field)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  std::vector<Atom>::const_iterator iterator =
      iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});
  insertAtomTags(static_cast<std::size_t>(atomOffset), atoms.size());

  std::vector<double3>::const_iterator iterator_electric_field =
      iteratorForElectricField(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  electricField.insert(iterator_electric_field, electric_field.begin(), electric_field.end());

  std::vector<Molecule>::iterator moleculeIterator =
      indexForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  moleculeData.insert(moleculeIterator, molecule);

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

  // the Coulomb-scaled charge of the inserted atoms (a fractional molecule counts with its scaled charge)
  const double insertedCharge = scaledNetCharge(atoms);
  netCharge += insertedCharge;
  netChargeAdsorbates += insertedCharge;
  netChargePerComponent[selectedComponent] += insertedCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(numberOfMoleculesPerComponent[selectedComponent]);
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
    addAtomToTailCorrectionCounts(atom);
  }

  updateMoleculeAtomInformation(selectedComponent);
}

void System::insertSerialReactionFractionalMolecule(std::size_t selectedComponent, std::size_t moleculeIndex,
                                                    [[maybe_unused]] const Molecule& molecule, std::vector<Atom> atoms,
                                                    double lambda, std::uint8_t dUdlambdaGroupId)
{
  // staged schedule from scaling.ixx: VDW switches on for lambda in [0, 0.5], Coulomb in [0.5, 1]
  for (Atom& atom : atoms)
  {
    atom.setScalingToFractional(lambda, dUdlambdaGroupId);
  }

  insertFractionalMoleculeAtIndex(selectedComponent, moleculeIndex, molecule, std::move(atoms));
}

void System::deleteFractionalMolecule(std::size_t selectedComponent, std::size_t selectedMolecule,
                                      const std::span<Atom> molecule)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  for (const Atom& atom : molecule)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] -= 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] -= 1;
    removeAtomFromTailCorrectionCounts(atom);
  }
  // the Coulomb-scaled charge of the removed atoms (a fractional molecule counts with its scaled charge); taken
  // before the erase below invalidates the span
  const double removedCharge = scaledNetCharge(molecule);

  // A linked molecule can not be removed (the moves guard against it); renumber the sites above it.
  crossLinks.moleculeDeleted(selectedComponent, selectedMolecule);

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  eraseAtomTags(static_cast<std::size_t>(atomOffset), molecule.size());
  atomData.erase(iterator, iterator + static_cast<std::vector<Atom>::difference_type>(molecule.size()));
  atomDynamics.erase(atomDynamics.cbegin() + atomOffset,
                     atomDynamics.cbegin() + atomOffset +
                         static_cast<std::vector<AtomDynamics>::difference_type>(molecule.size()));

  std::vector<double3>::const_iterator iterator_electric_field =
      iteratorForElectricField(selectedComponent, selectedMolecule);
  electricField.erase(iterator_electric_field,
                      iterator_electric_field + static_cast<std::vector<double3>::difference_type>(molecule.size()));

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, selectedMolecule);
  moleculeData.erase(moleculeIterator, moleculeIterator + 1);

  electricPotential.resize(electricPotential.size() - molecule.size());
  electricFieldNew.resize(electricFieldNew.size() - molecule.size());

  numberOfMoleculesPerComponent[selectedComponent] -= 1;
  numberOfFractionalMoleculesPerComponent[selectedComponent] -= 1;

  translationalDegreesOfFreedom -= components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom -= components[selectedComponent].rotationalDegreesOfFreedom;

  netCharge -= removedCharge;
  netChargeAdsorbates -= removedCharge;
  netChargePerComponent[selectedComponent] -= removedCharge;

  updateMoleculeAtomInformation(selectedComponent);
}

void System::deleteMolecule(std::size_t selectedComponent, std::size_t selectedMolecule, const std::span<Atom> molecule)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (const Atom& atom : molecule)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] -= 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] -= 1;
    removeAtomFromTailCorrectionCounts(atom);
  }
  // the Coulomb-scaled charge of the removed atoms (a fractional molecule counts with its scaled charge); taken
  // before the erase below invalidates the span
  const double removedCharge = scaledNetCharge(molecule);

  // A linked molecule can not be removed (the moves guard against it); renumber the sites above it.
  crossLinks.moleculeDeleted(selectedComponent, selectedMolecule);

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  eraseAtomTags(static_cast<std::size_t>(atomOffset), molecule.size());
  atomData.erase(iterator, iterator + static_cast<std::vector<Atom>::difference_type>(molecule.size()));
  atomDynamics.erase(atomDynamics.cbegin() + atomOffset,
                     atomDynamics.cbegin() + atomOffset +
                         static_cast<std::vector<AtomDynamics>::difference_type>(molecule.size()));

  std::vector<double3>::const_iterator iterator_electric_field =
      iteratorForElectricField(selectedComponent, selectedMolecule);
  electricField.erase(iterator_electric_field,
                      iterator_electric_field + static_cast<std::vector<double3>::difference_type>(molecule.size()));

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, selectedMolecule);
  moleculeData.erase(moleculeIterator, moleculeIterator + 1);

  electricPotential.resize(electricPotential.size() - molecule.size());
  electricFieldNew.resize(electricFieldNew.size() - molecule.size());

  numberOfMoleculesPerComponent[selectedComponent] -= 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] -= 1;

  netCharge -= removedCharge;
  netChargeAdsorbates -= removedCharge;
  netChargePerComponent[selectedComponent] -= removedCharge;

  translationalDegreesOfFreedom -= components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom -= components[selectedComponent].rotationalDegreesOfFreedom;

  updateMoleculeAtomInformation(selectedComponent);
}

void System::swapMolecules(std::size_t selectedComponent, std::size_t moleculeA, std::size_t moleculeB)
{
  exchangeMolecule(*this, selectedComponent, moleculeA, moleculeB);
}

void System::makeFractionalMoleculeWhole(std::size_t selectedComponent, std::size_t moleculeIndex)
{
  moveMolecule(selectedComponent, moleculeIndex, numberOfMoleculesPerComponent[selectedComponent] - 1);
  numberOfFractionalMoleculesPerComponent[selectedComponent] -= 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;
}

void System::makeWholeMoleculeFractional(std::size_t selectedComponent, std::size_t moleculeIndex,
                                         std::size_t slotIndex)
{
  moveMolecule(selectedComponent, moleculeIndex, slotIndex);
  numberOfFractionalMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] -= 1;
}

void System::exchangeAtomTags(std::size_t atomA, std::size_t atomB)
{
  std::swap(atomTags[atomA], atomTags[atomB]);
  atomIndexOfTag[atomTags[atomA]] = static_cast<std::uint32_t>(atomA);
  atomIndexOfTag[atomTags[atomB]] = static_cast<std::uint32_t>(atomB);
  std::swap(atomDynamics[atomA], atomDynamics[atomB]);
  std::swap(electricField[atomA], electricField[atomB]);
}

void System::checkPolymerLayoutIsFixed(std::string_view operation) const
{
  if (polymerLinks.activated)
  {
    throw std::runtime_error(std::format(
        "[Polymerization]: {} after the polymerization was activated; its group bookkeeping is laid out over the "
        "molecules at activation and follows no insertion, deletion or move\n",
        operation));
  }
}

void System::moveMolecule(std::size_t selectedComponent, std::size_t moleculeIndex, std::size_t newIndex)
{
  checkPolymerLayoutIsFixed("a molecule was inserted, deleted or moved");
  if (moleculeIndex == newIndex) return;

  // the per-atom arrays and the records of the molecules from the lower to the higher index rotate by one
  // molecule: a molecule moving up takes the first block to the back, one moving down the last block to the front
  const std::ptrdiff_t numberOfAtoms = static_cast<std::ptrdiff_t>(components[selectedComponent].atoms.size());
  const std::size_t lowIndex = std::min(moleculeIndex, newIndex);
  const std::size_t highIndex = std::max(moleculeIndex, newIndex);
  const std::ptrdiff_t first = iteratorForMolecule(selectedComponent, lowIndex) - atomData.begin();
  const std::ptrdiff_t last = first + static_cast<std::ptrdiff_t>(highIndex - lowIndex + 1) * numberOfAtoms;
  const std::ptrdiff_t middle = (moleculeIndex < newIndex) ? first + numberOfAtoms : last - numberOfAtoms;
  std::rotate(atomData.begin() + first, atomData.begin() + middle, atomData.begin() + last);
  std::rotate(atomDynamics.begin() + first, atomDynamics.begin() + middle, atomDynamics.begin() + last);
  std::rotate(electricField.begin() + first, electricField.begin() + middle, electricField.begin() + last);
  std::rotate(atomTags.begin() + first, atomTags.begin() + middle, atomTags.begin() + last);
  for (std::ptrdiff_t i = first; i < last; ++i)
  {
    atomIndexOfTag[atomTags[static_cast<std::size_t>(i)]] = static_cast<std::uint32_t>(i);
  }

  const std::ptrdiff_t firstRecord = indexForMolecule(selectedComponent, lowIndex) - moleculeData.begin();
  const std::ptrdiff_t lastRecord = firstRecord + static_cast<std::ptrdiff_t>(highIndex - lowIndex + 1);
  const std::ptrdiff_t middleRecord = (moleculeIndex < newIndex) ? firstRecord + 1 : lastRecord - 1;
  std::rotate(moleculeData.begin() + firstRecord, moleculeData.begin() + middleRecord,
              moleculeData.begin() + lastRecord);

  // the sites of the moved molecule and of the molecules in between are renumbered
  crossLinks.moleculeMoved(selectedComponent, moleculeIndex, newIndex);

  // molecule ids, atom indices of the records and the cell list
  updateMoleculeAtomInformation(selectedComponent);
}

void System::exchangeMolecule(System& other, std::size_t selectedComponent, std::size_t molecule,
                              std::size_t otherMolecule)
{
  checkPolymerLayoutIsFixed("two molecules were exchanged");
  other.checkPolymerLayoutIsFixed("two molecules were exchanged");
  std::span<Atom> atoms = spanOfMolecule(selectedComponent, molecule);
  std::span<Atom> otherAtoms = other.spanOfMolecule(selectedComponent, otherMolecule);
  const std::ptrdiff_t offset = atoms.data() - atomData.data();
  const std::ptrdiff_t otherOffset = otherAtoms.data() - other.atomData.data();
  const std::ptrdiff_t count = static_cast<std::ptrdiff_t>(atoms.size());

  std::swap_ranges(atoms.begin(), atoms.end(), otherAtoms.begin());
  std::swap_ranges(atomDynamics.begin() + offset, atomDynamics.begin() + offset + count,
                   other.atomDynamics.begin() + otherOffset);
  std::swap_ranges(electricField.begin() + offset, electricField.begin() + offset + count,
                   other.electricField.begin() + otherOffset);

  // A tag names an atom within its system: it follows the atom in an exchange within the system and stays with
  // the slot in an exchange between systems.
  if (&other == this)
  {
    std::swap_ranges(atomTags.begin() + offset, atomTags.begin() + offset + count, atomTags.begin() + otherOffset);
    for (std::ptrdiff_t i = 0; i < count; ++i)
    {
      atomIndexOfTag[atomTags[static_cast<std::size_t>(offset + i)]] = static_cast<std::uint32_t>(offset + i);
      atomIndexOfTag[atomTags[static_cast<std::size_t>(otherOffset + i)]] = static_cast<std::uint32_t>(otherOffset + i);
    }
  }

  // The molecule id of an atom and the atom index of a record name the slot, so they stay behind (the component
  // id is the same on both sides).
  const std::size_t index = moleculeIndexOfComponent(selectedComponent, molecule);
  const std::size_t otherIndex = other.moleculeIndexOfComponent(selectedComponent, otherMolecule);
  for (Atom& atom : atoms) atom.moleculeId = static_cast<std::uint32_t>(index);
  for (Atom& atom : otherAtoms) atom.moleculeId = static_cast<std::uint32_t>(otherIndex);
  std::swap(moleculeData[index], other.moleculeData[otherIndex]);
  std::swap(moleculeData[index].atomIndex, other.moleculeData[otherIndex].atomIndex);

  // the cross-link sites follow their molecules; a link can not reach into another system
  if (&other == this)
  {
    crossLinks.moleculesExchanged(selectedComponent, molecule, otherMolecule);
  }
  else if (crossLinks.moleculeIsLinked(selectedComponent, molecule) ||
           other.crossLinks.moleculeIsLinked(selectedComponent, otherMolecule))
  {
    throw std::runtime_error(std::format(
        "[System]: molecule {} of component {} was exchanged between systems while carrying a cross-link\n",
        crossLinks.moleculeIsLinked(selectedComponent, molecule) ? molecule : otherMolecule, selectedComponent));
  }

  // the cell lists describe the two slots as they are now
  cellListAtomsMoved(atoms);
  other.cellListAtomsMoved(otherAtoms);
}

namespace
{
// Exchanges the parts of two per-atom arrays beyond their framework atoms.
template <typename T>
void swapMobileTail(std::vector<T>& dataA, std::size_t fixedSizeA, std::vector<T>& dataB, std::size_t fixedSizeB)
{
  std::vector<T> mobileA(std::make_move_iterator(dataA.begin() + static_cast<std::ptrdiff_t>(fixedSizeA)),
                         std::make_move_iterator(dataA.end()));
  std::vector<T> mobileB(std::make_move_iterator(dataB.begin() + static_cast<std::ptrdiff_t>(fixedSizeB)),
                         std::make_move_iterator(dataB.end()));
  dataA.erase(dataA.begin() + static_cast<std::ptrdiff_t>(fixedSizeA), dataA.end());
  dataB.erase(dataB.begin() + static_cast<std::ptrdiff_t>(fixedSizeB), dataB.end());
  dataA.insert(dataA.end(), std::make_move_iterator(mobileB.begin()), std::make_move_iterator(mobileB.end()));
  dataB.insert(dataB.end(), std::make_move_iterator(mobileA.begin()), std::make_move_iterator(mobileA.end()));
}
}  // namespace

void System::exchangeConfiguration(System& other)
{
  checkPolymerLayoutIsFixed("the configurations of two systems were exchanged");
  other.checkPolymerLayoutIsFixed("the configurations of two systems were exchanged");
  swapMobileTail(atomData, numberOfFrameworkAtoms, other.atomData, other.numberOfFrameworkAtoms);
  swapMobileTail(atomDynamics, numberOfFrameworkAtoms, other.atomDynamics, other.numberOfFrameworkAtoms);
  // the tags travel with their atoms; replicas share the framework and so its tags, and then each inverse table
  // describes the other's atoms (should the framework tags differ, the atoms are numbered afresh)
  if (std::equal(atomTags.begin(), atomTags.begin() + static_cast<std::ptrdiff_t>(numberOfFrameworkAtoms),
                 other.atomTags.begin(),
                 other.atomTags.begin() + static_cast<std::ptrdiff_t>(other.numberOfFrameworkAtoms)))
  {
    swapMobileTail(atomTags, numberOfFrameworkAtoms, other.atomTags, other.numberOfFrameworkAtoms);
    std::swap(atomIndexOfTag, other.atomIndexOfTag);
    std::swap(lowestFreeTagBound, other.lowestFreeTagBound);
  }
  else
  {
    resetAtomTags();
    other.resetAtomTags();
  }
  std::swap(moleculeData, other.moleculeData);
  // rigid-body state of semi-flexible molecules: derived from the positions, so it moves with them
  std::swap(groupData, other.groupData);
  // the Monte Carlo cell lists describe the previous configurations; rebuilt lazily on first use
  invalidateCellList();
  other.invalidateCellList();
  if (!framework.has_value())
  {
    std::swap(simulationBox, other.simulationBox);
  }
  std::swap(numberOfMoleculesPerComponent, other.numberOfMoleculesPerComponent);
  std::swap(numberOfIntegerMoleculesPerComponent, other.numberOfIntegerMoleculesPerComponent);
  std::swap(numberOfFractionalMoleculesPerComponent, other.numberOfFractionalMoleculesPerComponent);
  swapMobileTail(electricPotential, numberOfFrameworkAtoms, other.electricPotential, other.numberOfFrameworkAtoms);
  swapMobileTail(electricField, numberOfFrameworkAtoms, other.electricField, other.numberOfFrameworkAtoms);
  swapMobileTail(electricFieldNew, numberOfFrameworkAtoms, other.electricFieldNew, other.numberOfFrameworkAtoms);
  std::swap(netChargeAdsorbates, other.netChargeAdsorbates);
  std::swap(netChargePerComponent, other.netChargePerComponent);
  std::swap(translationalCenterOfMassConstraint, other.translationalCenterOfMassConstraint);
  std::swap(translationalDegreesOfFreedom, other.translationalDegreesOfFreedom);
  std::swap(rotationalDegreesOfFreedom, other.rotationalDegreesOfFreedom);
  std::swap(containsTheFractionalMolecule, other.containsTheFractionalMolecule);
  // The cross-link topology belongs to the configuration; the bond types (Hamiltonian) stay.
  crossLinks.swapTopology(other.crossLinks);
  // the time-correlation properties belong to the trajectory: their time origins and sums follow the walker, so
  // every correlation is taken over one continuous trajectory (a replica writes those of the walker it holds)
  using std::swap;
  swap(propertyMSD, other.propertyMSD);
  swap(propertyVACF, other.propertyVACF);
  swap(propertyEndToEndACF, other.propertyEndToEndACF);
  // the lambda coordinates describe the fractional molecules of the configuration: every one of them moves along
  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    Component& mine = components[componentId];
    Component& theirs = other.components[componentId];
    std::swap(mine.lambdaGC.currentBin, theirs.lambdaGC.currentBin);
    std::swap(mine.lambdaGibbs.currentBin, theirs.lambdaGibbs.currentBin);
    std::swap(mine.lambdaPairSwap.currentBin, theirs.lambdaPairSwap.currentBin);
    std::swap(mine.lambdaPairSwapCB.currentBin, theirs.lambdaPairSwapCB.currentBin);
    std::swap(mine.lambdaGroupSwap.currentBin, theirs.lambdaGroupSwap.currentBin);
    std::swap(mine.lambdaGroupSwapCB.currentBin, theirs.lambdaGroupSwapCB.currentBin);
  }
  for (std::size_t reaction = 0; reaction < std::min(reactions.list.size(), other.reactions.list.size()); ++reaction)
  {
    std::swap(reactions.list[reaction].lambda.currentBin, other.reactions.list[reaction].lambda.currentBin);
    std::swap(reactions.list[reaction].lambdaProductSide.currentBin,
              other.reactions.list[reaction].lambdaProductSide.currentBin);
  }
}

void System::insertAtomTags(std::size_t offset, std::size_t count)
{
  // free tags exist exactly when the table is longer than the atom count: the lowest ones go first, searched
  // from the lower bound of the lowest free tag (every slot below it is in use)
  std::size_t freeTag = atomIndexOfTag.size() > atomTags.size() ? std::min(lowestFreeTagBound, atomIndexOfTag.size())
                                                                : atomIndexOfTag.size();
  atomTags.insert(atomTags.begin() + static_cast<std::ptrdiff_t>(offset), count, noAtom);
  for (std::size_t i = offset; i < offset + count; ++i)
  {
    while (freeTag < atomIndexOfTag.size() && atomIndexOfTag[freeTag] != noAtom) ++freeTag;
    if (freeTag == atomIndexOfTag.size()) atomIndexOfTag.push_back(noAtom);
    atomTags[i] = static_cast<std::uint32_t>(freeTag);
    atomIndexOfTag[freeTag] = static_cast<std::uint32_t>(i);
  }
  lowestFreeTagBound = freeTag;
  for (std::size_t i = offset + count; i < atomTags.size(); ++i)
  {
    atomIndexOfTag[atomTags[i]] = static_cast<std::uint32_t>(i);
  }
}

void System::eraseAtomTags(std::size_t offset, std::size_t count)
{
  for (std::size_t i = offset; i < offset + count; ++i)
  {
    atomIndexOfTag[atomTags[i]] = noAtom;
    lowestFreeTagBound = std::min<std::size_t>(lowestFreeTagBound, atomTags[i]);
  }
  atomTags.erase(atomTags.begin() + static_cast<std::ptrdiff_t>(offset),
                 atomTags.begin() + static_cast<std::ptrdiff_t>(offset + count));
  while (!atomIndexOfTag.empty() && atomIndexOfTag.back() == noAtom) atomIndexOfTag.pop_back();
  for (std::size_t i = offset; i < atomTags.size(); ++i)
  {
    atomIndexOfTag[atomTags[i]] = static_cast<std::uint32_t>(i);
  }
}

void System::resetAtomTags()
{
  atomTags.resize(atomData.size());
  std::iota(atomTags.begin(), atomTags.end(), std::uint32_t{0});
  atomIndexOfTag = atomTags;
  lowestFreeTagBound = 0;
}

void System::rebuildAtomIndexOfTag()
{
  // the free-tag marker is no tag: a table sized past it would overflow, so a restart that holds it is refused
  // before anything is written
  for (std::size_t i = 0; i < atomTags.size(); ++i)
  {
    if (atomTags[i] == noAtom)
    {
      throw std::runtime_error(std::format("Atom tags: atom {} holds the free-tag marker as its tag\n", i));
    }
  }
  // a tag names an atom that was alive at some time, so the table is as long as the largest atom count the system
  // ever held; a corrupt tag can ask for a table the machine does not have, which is reported as such
  const std::size_t highest =
      atomTags.empty() ? 0 : static_cast<std::size_t>(*std::max_element(atomTags.begin(), atomTags.end())) + 1;
  try
  {
    atomIndexOfTag.assign(highest, noAtom);
    lowestFreeTagBound = 0;
  }
  catch (const std::bad_alloc&)
  {
    throw std::runtime_error(
        std::format("Atom tags: the highest tag {} needs a table of {} entries, which can not be allocated\n",
                    highest - 1, highest));
  }
  for (std::size_t i = 0; i < atomTags.size(); ++i)
  {
    atomIndexOfTag[atomTags[i]] = static_cast<std::uint32_t>(i);
  }
}

void System::checkAtomTags() const
{
  if (atomTags.size() != atomData.size())
  {
    throw std::runtime_error(std::format("Atom tags: {} tags for {} atoms\n", atomTags.size(), atomData.size()));
  }
  for (std::size_t i = 0; i < atomTags.size(); ++i)
  {
    if (atomTags[i] >= atomIndexOfTag.size() || atomIndexOfTag[atomTags[i]] != i)
    {
      throw std::runtime_error(std::format("Atom tags: atom {} holds tag {} whose inverse entry is {}\n", i,
                                           atomTags[i],
                                           atomTags[i] < atomIndexOfTag.size() ? atomIndexOfTag[atomTags[i]] : 0));
    }
  }
  const std::size_t live = static_cast<std::size_t>(
      std::count_if(atomIndexOfTag.begin(), atomIndexOfTag.end(), [](std::uint32_t index) { return index != noAtom; }));
  if (live != atomTags.size() || (!atomIndexOfTag.empty() && atomIndexOfTag.back() == noAtom))
  {
    throw std::runtime_error(std::format("Atom tags: {} live entries in a table of {} for {} atoms\n", live,
                                         atomIndexOfTag.size(), atomTags.size()));
  }
}

void System::updateMoleculeAtomInformation(std::size_t firstComponent)
{
  // The atom array changed size or order: the cell list indices and molecule ids are stale.
  moleculeCellList.invalidate();

  std::size_t atom_index = numberOfFrameworkAtoms;
  std::size_t molecule_index{};
  for (std::size_t k = 0; k < firstComponent; k++)
  {
    atom_index += numberOfMoleculesPerComponent[k] * components[k].atoms.size();
    molecule_index += numberOfMoleculesPerComponent[k];
  }

  for (std::size_t k = firstComponent; k < components.size(); k++)
  {
    std::size_t numberOfAtoms = components[k].atoms.size();

    for (std::size_t i = 0; i < numberOfMoleculesPerComponent[k]; ++i)
    {
      moleculeData[molecule_index].atomIndex = atom_index - numberOfFrameworkAtoms;
      moleculeData[molecule_index].numberOfAtoms = numberOfAtoms;

      for (std::size_t j = 0; j < numberOfAtoms; ++j)
      {
        atomData[atom_index].moleculeId = static_cast<std::uint32_t>(molecule_index);
        atomData[atom_index].componentId = static_cast<std::uint8_t>(k);
        ++atom_index;
      }
      ++molecule_index;
    }
  }
}

void System::checkMoleculeIds()
{
  std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  std::size_t index = 0;
  std::size_t molecule_index = 0;
  for (std::size_t k = 0; k < components.size(); k++)
  {
    for (std::size_t i = 0; i < numberOfMoleculesPerComponent[k]; ++i)
    {
      for (std::size_t j = 0; j < components[k].atoms.size(); ++j)
      {
        if (moleculeAtoms[index].moleculeId != static_cast<std::uint32_t>(molecule_index))
        {
          throw std::runtime_error(std::format("Wrong molecule-id detected {} for global molecule {}\n",
                                               moleculeAtoms[index].moleculeId, molecule_index));
        }
        if (moleculeAtoms[index].componentId != static_cast<std::uint8_t>(k))
        {
          throw std::runtime_error(std::format("Wrong component-id detected {} for component {} molecule {}\n",
                                               moleculeAtoms[index].componentId, k, i));
        }
        ++index;
      }
      ++molecule_index;
    }
  }
  checkAtomTags();

  for (std::size_t componentId = 0; componentId < components.size(); componentId++)
  {
    if (numberOfGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
    {
      std::size_t indexFractionalMolecule = indexOfGCFractionalMoleculesPerComponent_CFCMC(componentId);
      std::span<Atom> fractionalMolecule = spanOfMolecule(componentId, indexFractionalMolecule);

      const std::uint8_t expectedGroupId = components[componentId].lambdaGC.dUdlambdaGroupId;
      for (const Atom& atom : fractionalMolecule)
      {
        if (atom.groupId != expectedGroupId)
        {
          throw std::runtime_error(std::format("Wrong group-id detected! ({} where it should be {})\n",
                                               static_cast<std::size_t>(atom.groupId),
                                               static_cast<std::size_t>(expectedGroupId)));
        }
      }
    }
  }
}


void System::createInitialMolecules(const std::vector<std::vector<double3>>& initialPositions)
{
  // keep a fixed seed
  RandomNumber random(1200);

  // Initial configurations are grown with CBMC even when the production moves use recoil growth: CBMC
  // always completes (a bad step only shrinks the Rosenbluth weight), whereas the absolute open/closed
  // test of recoil growth restarts a long chain indefinitely -- per-step attrition compounds over
  // hundreds of beads, and the recoil window aborts the attempt long before the end.
  // (The context is rebuilt per grow: every inserted molecule changes the background it spans.)
  auto creationContext = [&]()
  { return makeGrowContext(CBMC::CutOffMode::Full).withChainScheme(CBMC::ChainScheme::ConfigurationalBias); };

  // A grown molecule is valid when it exists and does not overlap with its environment.
  auto isValid = [&](const std::optional<CBMC::GrowResult>& growData)
  { return growData.has_value() && growData->energies.potentialEnergy() <= forceField.energyOverlapCriteria; };

  // The configuration a component is created with. A growth alone samples a conformation with
  // probability prod_i exp(-beta u_i) / w_i, not the Boltzmann distribution: the missing factor is its
  // Rosenbluth weight W = prod_i w_i. The first valid growth is therefore taken as the start of a short
  // Markov chain of 'creationCandidateGrowths' - 1 further growths, each accepted over the current one
  // with min(1, W_new / W_old) -- the reinsertion acceptance with the stored weights, whose stationary
  // distribution is exactly Boltzmann (detailed balance holds in the joint space of conformation and
  // trial set). Without it, every conformation whose growth was judged by terms that only enter at a
  // later step is over-represented: the ester torsion of a diacrylate chosen before its 1-5 partners
  // exist leaves a twisted or E ester that the subsequent Markov chain (or, worse, molecular dynamics)
  // may take microseconds to undo. Such a conformation carries a W orders of magnitude below the rest,
  // so the next ordinary candidate replaces it with probability close to one, while an ordinary
  // conformation is kept unless a better-weighted one comes along. (A single draw with probability
  // proportional to W among the candidates is not used: the spread of log W over the growths of a
  // flexible chain spans several units, which turns that draw into 'the largest weight wins' and
  // over-corrects towards the hardest-to-generate conformations.) A rigid molecule is grown once: its
  // one placement step has no conformation to bias. The Markov chain of the simulation still
  // equilibrates away from the created state; the creation chain just starts it far closer.
  auto growCreationCandidate = [&](const CBMC::NewMoleculeIdentity& identity) -> CBMC::GrowResult
  {
    const Component& component = components[identity.componentId];

    // The grow filters every trial against the blocking pockets, so the returned molecule needs no
    // further pocket check.
    std::optional<CBMC::GrowResult> current = std::nullopt;
    do
    {
      current = CBMC::growNewMolecule(random, creationContext(), component, identity);
    } while (!isValid(current));

    // An overlapping proposal has no weight and counts as rejected. In a dense background that is the
    // fate of most growths, and there the weight of a valid candidate is dominated by its fit into the
    // fluid rather than by its conformation, so the chain is deliberately not extended until it has a
    // fixed number of valid proposals (measured: 7x the creation time for a marginal gain). A dense
    // system of flexible molecules still relies on its initialization cycles -- partial reinsertion
    // with fixed endpoints regrows such a defect in place -- to finish the job.
    const std::size_t numberOfGrowths = component.rigid ? 1 : CBMC::Constants::creationCandidateGrowths;
    for (std::size_t attempt = 1; attempt < numberOfGrowths; ++attempt)
    {
      std::optional<CBMC::GrowResult> candidate = CBMC::growNewMolecule(random, creationContext(), component, identity);
      if (!isValid(candidate)) continue;

      const double logAcceptance = candidate->logRosenbluthWeight - current->logRosenbluthWeight;
      if (logAcceptance >= 0.0 || random.uniform() < std::exp(logAcceptance))
      {
        current = std::move(candidate);
      }
    }
    return std::move(current.value());
  };

  for (std::size_t componentId = 0; const Component& component : components)
  {
    if (component.hasFractionalMolecule)
    {
      numberOfMoleculesPerComponent[componentId] = 0;

      auto growFractionalMolecule = [&](std::uint8_t groupId) -> std::optional<CBMC::GrowResult>
      {
        return growCreationCandidate({.componentId = componentId,
                                      .moleculeId = numberOfMolecules(),
                                      .scaling = 0.0,
                                      .groupId = groupId,
                                      .isFractional = true});
      };

      if (numberOfGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfGCFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairGCFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      // group-swap slots: a component can hold several (e.g. two Cl slots of a CaCl2 group)
      for (std::size_t k = 0; k < numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[componentId]; ++k)
      {
        const std::size_t slot = indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(componentId) + k;
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      for (std::size_t k = 0; k < numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[componentId]; ++k)
      {
        const std::size_t slot = indexOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(componentId) + k;
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<CBMC::GrowResult> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfGibbsFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        // The conventional Gibbs fractional molecule starts at lambda = 0.5 (full VDW scaling), not switched off:
        // grow it at that scaling so that the placement sees the framework and the blocking pockets (a grow at
        // scaling zero has no energy and no pocket radius, and may land inside a wall or a pocket).
        const std::uint8_t groupId = components[componentId].lambdaGibbs.dUdlambdaGroupId;
        const std::optional<CBMC::GrowResult> growData = growCreationCandidate({.componentId = componentId,
                                                                                .moleculeId = numberOfMolecules(),
                                                                                .scaling = 0.5,
                                                                                .groupId = groupId,
                                                                                .isFractional = true});
        std::vector<Atom> atoms = growData->atoms;
        for (Atom& atom : atoms)
        {
          atom.setScalingToFractional(0.5, groupId);
        }
        insertFractionalMoleculeAtIndex(
            componentId, indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(componentId),
            growData->molecule, std::move(atoms));
      }
    }

    // Broensted protons: one single-site proton per candidate group, placed on the group's first
    // site (stored as fractional coordinates in the simulation box). The proton count is fully
    // determined by the tether groups, so no random growth is done for this component below.
    if (!component.tetheredProtonHopSiteGroups.empty())
    {
      for (const std::vector<double3>& siteGroup : component.tetheredProtonHopSiteGroups)
      {
        if (siteGroup.empty()) continue;

        const double3 position = simulationBox.cell * siteGroup.front();

        Molecule molecule(position, simd_quatd{}, component.totalMass, componentId, component.atoms.size());

        std::vector<Atom> molecule_atoms = component.atoms;
        for (Atom& atom : molecule_atoms)
        {
          atom.position = position;
        }

        // insertMolecule already adds the translational and rotational degrees of freedom of the proton
        insertMolecule(componentId, molecule, molecule_atoms);
      }

      componentId++;
      continue;
    }

    // Add atom position passed to 'createInitialMolecules'
    if (componentId < initialPositions.size())
    {
      std::vector<double3> positions = initialPositions[componentId];
      for (std::size_t i = 0; i < positions.size(); i += component.atoms.size())
      {
        std::span<double3> position_view = std::span<double3>(&positions[i], component.atoms.size());

        // the center of mass, for flexible molecules too
        double3 com(0.0, 0.0, 0.0);
        double totalMass = 0.0;
        for (std::size_t k = 0; k < position_view.size(); ++k)
        {
          double mass = forceField.pseudoAtoms[static_cast<std::size_t>(components[componentId].atoms[k].type)].mass;
          com += mass * position_view[k];
          totalMass += mass;
        }
        com /= totalMass;

        simd_quatd orientation{};

        if (components[componentId].rigid)
        {
          double3 reference_com = double3{0.0, 0.0, 0.0};

          std::vector<double3> reference_positions =
              components[componentId].atoms | std::views::transform(&Atom::position) | std::ranges::to<std::vector>();

          // Compute rotation matrix that describes going from the space-fixed frame to the body-fixed frame
          double3x3 rotation_matrix =
              double3x3::computeRotationMatrix(com, position_view, reference_com, reference_positions);

          // Get orientation: buildRotationMatrix(q) is the space-to-body rotation, so q is the quaternion of its
          // transpose
          orientation = rotation_matrix.transpose().quaternion();
        }

        Molecule molecule =
            Molecule(com, orientation, component.totalMass, componentId, components[componentId].atoms.size());

        std::vector<Atom> molecule_atoms = components[componentId].atoms;
        for (std::size_t k = 0; k < position_view.size(); ++k)
        {
          molecule_atoms[k].position = position_view[k];
        }

        insertMolecule(componentId, molecule, molecule_atoms);
      }
    }

    for (std::size_t i = 0; i < initialNumberOfMolecules[componentId]; ++i)
    {
      const CBMC::GrowResult growData =
          growCreationCandidate({.componentId = componentId, .moleculeId = numberOfMolecules()});

      insertMolecule(componentId, growData.molecule, growData.atoms);
    }

    componentId++;
  }
}

void System::computeAutomaticBlockingPockets()
{
  // The pockets depend on the framework and on the probe, neither of which is a property of a component, so
  // every component that asked for them is answered with the same spheres and the framework is measured once.
  std::optional<std::vector<double4>> computed{};

  for (Component& component : components)
  {
    if (!component.automaticBlockingPockets) continue;

    if (!framework.has_value())
    {
      throw std::runtime_error(
          std::format("[System]: component '{}' asks for its blocking pockets to be computed, but the system has "
                      "no framework to compute them from; a pocket is a cavity of a framework\n",
                      component.name));
    }

    if (!computed.has_value())
    {
      computed = BlockingPockets::compute(framework.value(), forceField);
    }

    component.blockingPockets.insert(component.blockingPockets.end(), computed->begin(), computed->end());

    // The spheres are listed like any others now, and a system built again from these components must not
    // measure the framework a second time and list them twice.
    component.automaticBlockingPockets = false;
  }
}

bool System::insideBlockedPockets(const Component& component, std::span<const Atom> molecule_atoms) const
{
  // Same test the CBMC grow applies to every trial set; a molecule returned by a CBMC grow therefore
  // never needs this check again. It remains for placements that bypass CBMC and for fractional
  // molecules whose 'scalingVDW' (which scales the pocket radius) is changed in place.
  return CBMC::insideBlockedPockets(framework, component, molecule_atoms);
}
std::vector<Atom> System::randomConfiguration(RandomNumber& random, std::size_t selectedComponent,
                                              const std::span<const Atom> molecule)
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  double3 position = simulationBox.randomPosition(random);
  std::size_t startingBead = components[selectedComponent].startingBead;
  for (std::size_t i = 0; i != molecule.size(); ++i)
  {
    copied_atoms[i].position =
        position + randomRotationMatrix * (molecule[i].position - molecule[startingBead].position);
  }
  return copied_atoms;
}

std::vector<Atom> System::equilibratedIdealGasConformation(RandomNumber& random, std::size_t selectedComponent)
{
  Component& component = components[selectedComponent];

  // Persistent scratch conformation stored on the Component: initialised once to the reference geometry,
  // then advanced in place so the reinsertion Markov chain stays warm between calls.
  std::vector<Atom>& scratchAtoms = component.grownIdealGasAtoms;
  const bool fresh = scratchAtoms.size() != component.atoms.size();
  if (fresh) scratchAtoms = component.atoms;

  Molecule scratchMolecule(double3(0.0, 0.0, 0.0), simd_quatd(0.0, 0.0, 0.0, 1.0), component.totalMass,
                           selectedComponent, component.atoms.size());

  if (component.atoms.size() < 2) return scratchAtoms;

  // Isolated (ideal-gas) growth context: CBMC growth feels only the intra-molecular potential and hence
  // samples exp(-beta * U_intra). Always configurational bias (see 'makeIdealGasGrowContext'): these
  // conformations seed the conformation reservoir and define the recoil-growth openness reference, and
  // recoil growth can not run before its own reference exists.
  const CBMC::GrowContext context = makeIdealGasGrowContext();

  // The reference geometry is a topology with placeholder coordinates (its atoms may coincide), not a
  // conformation: a retrace of it has no weight. The chain starts from a grown conformation instead; the
  // reinsertion moves below turn the Rosenbluth-weighted growth into the Boltzmann distribution.
  if (fresh)
  {
    for (std::size_t attempt = 0; attempt != 100; ++attempt)
    {
      std::optional<CBMC::GrowResult> growData = CBMC::regrowMolecule(
          random, context, component, scratchMolecule, scratchAtoms, {.firstBead = CBMC::FirstBeadScheme::Reinsertion});
      if (!growData) continue;
      std::copy(growData->atoms.begin(), growData->atoms.end(), scratchAtoms.begin());
      scratchMolecule = growData->molecule;
      break;
    }
  }

  // A handful of full-molecule reinsertion moves decorrelates the conformation from the starting geometry;
  // each accepted move regrows the whole chain from the ideal-gas Boltzmann distribution.
  constexpr std::size_t numberOfReinsertionMoves = 20;
  for (std::size_t move = 0; move != numberOfReinsertionMoves; ++move)
  {
    std::optional<CBMC::GrowResult> growData = CBMC::regrowMolecule(
        random, context, component, scratchMolecule, scratchAtoms, {.firstBead = CBMC::FirstBeadScheme::Reinsertion});
    if (!growData) continue;

    const CBMC::RetraceResult retraceData = CBMC::retraceMolecule(
        random, context, component, scratchAtoms,
        {.firstBead = CBMC::FirstBeadScheme::Reinsertion, .logStoredR = growData->firstBeadLogStoredR});

    // Metropolis acceptance for the reinsertion move in the isolated system (no Ewald/polarization/tail
    // corrections apply): accept with min(1, W_new / W_old), with the ratio evaluated in log space so it
    // stays exact for long chains whose raw Rosenbluth weights underflow to zero.
    if (random.uniform() < std::exp(growData->logRosenbluthWeight - retraceData.logRosenbluthWeight))
    {
      std::copy(growData->atoms.begin(), growData->atoms.end(), scratchAtoms.begin());
      scratchMolecule = growData->molecule;
    }
  }

  // Re-center on the mass center so the caller can place the molecule by its center of mass.
  double totalMass = 0.0;
  double3 com(0.0, 0.0, 0.0);
  for (const Atom& atom : scratchAtoms)
  {
    double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
    com += mass * atom.position;
    totalMass += mass;
  }
  if (totalMass > 0.0) com /= totalMass;
  for (Atom& atom : scratchAtoms) atom.position -= com;

  return scratchAtoms;
}

void System::prepareGrowthPlans()
{
  for (Component& component : components)
  {
    component.prepareGrowthPlans(beta);
  }
}

void System::buildConformationReservoirs()
{
  // The reservoirs are Boltzmann samples at the current temperature, and so are the base-coupling
  // constants of the growth plans that grow them: bring those up to date first (a no-op when the
  // temperature is unchanged).
  prepareGrowthPlans();

  // Number of independent seed conformations kept per component. A few hundred is enough to cover the
  // ring conformers and decorrelate successive grows; each is a cheap isolated-molecule regrow, so the
  // whole build is well under a second even for the largest adsorbates and is a one-time cost.
  constexpr std::size_t reservoirSize = 256;

  // A fixed local generator keeps the reservoir reproducible and independent of the simulation's random
  // stream (so enabling it does not shift any downstream move sequence).
  RandomNumber random(1400);

  for (std::size_t componentId = 0; componentId != components.size(); ++componentId)
  {
    Component& component = components[componentId];

    // Fully rigid components have no internal conformation to seed (the rigid tilt uses uniform random
    // orientations); single-atom components have nothing to grow.
    if (component.atoms.size() < 2 || component.rigid) continue;

    // The reservoir is only read by CBMC grows (ring closure seeds). A component no configured move can
    // grow (an MD-only protein read from a topology file) never consults it, and the 256 full-molecule
    // regrows with their lookahead tables cost minutes for a ringed peptide: skip the build.
    if (!component.mayUseGrowthPlans()) continue;

    // The reservoir carries the long-lived internal state a fresh grow cannot rebuild from potentials:
    // ring conformers ('generateRingConformation' seeds each ring-closure step from an independent,
    // well-mixed whole-molecule conformation so different grows start from different puckers) and the
    // exact internal geometry of rigid sub-fragments (fixed by construction, so no bond/bend potential
    // exists to sample it from -- a cold start would place those atoms at guessed geometry and every
    // grown configuration would be rejected downstream). A ring-free, fully flexible component has no
    // such state: every flexible step freshly samples its bond lengths, bend angles, and torsion spin,
    // and the cold-start path of 'generateFlexibleBaseConformation' (which doubles the internal-MC
    // length) relaxes the step-local geometry just as well. Skipping the build for those components
    // avoids its dominant cost: for a polymer chain the 256 x 20 full-molecule regrows take minutes
    // and dwarf the simulation itself.
    if (component.fragmentGraph.closureBonds.empty() && component.numberOfRigidFragments() == 0)
    {
      continue;
    }

    // Build into a temporary so the reservoir stays empty while it is being filled: the ideal-gas grows
    // that fill it then cold-start (their own Metropolis still targets exp(-beta * U_intra) regardless
    // of the seed), avoiding any feedback from a half-filled reservoir.
    //
    // The reservoir is only a warm start for the bond/bend Monte Carlo. It carries no Rosenbluth
    // weight; 'equilibratedIdealGasConformation' grows it with configurational bias even when the
    // production moves use recoil growth.
    std::vector<std::vector<Atom>> reservoir{};
    reservoir.reserve(reservoirSize);
    for (std::size_t i = 0; i != reservoirSize; ++i)
    {
      reservoir.push_back(equilibratedIdealGasConformation(random, componentId));
    }
    component.conformationReservoir = std::move(reservoir);
  }
}

void System::buildRecoilReferenceConformations()
{
  if (!forceField.settings.useRecoilGrowth) return;

  // The reference conformations are Boltzmann samples at the current temperature, grown by plans whose
  // base-coupling constants must be those of this temperature (a no-op when already prepared).
  prepareGrowthPlans();

  // The openness reference takes, per growth step, the maximum energy over these conformations (any
  // strain a valid equilibrated chain exhibits at a step is acceptable there), so more conformations
  // widen the accepted strain range toward its Boltzmann-typical upper edge. Count and seed are the
  // CBMC constants (see 'CBMC::Constants::recoilReferenceConformations' for why the build is fixed).
  constexpr std::size_t numberOfConformations = CBMC::Constants::recoilReferenceConformations;
  RandomNumber random(CBMC::Constants::recoilReferenceSeed);

  for (std::size_t componentId = 0; componentId != components.size(); ++componentId)
  {
    Component& component = components[componentId];
    if (component.atoms.size() < 2 || component.rigid) continue;

    // 'equilibratedIdealGasConformation' grows with configurational bias (which always completes)
    // rather than recoil growth itself: recoil growth cannot run before its own reference exists, and
    // CBMC samples the same ideal-gas Boltzmann distribution.
    std::vector<std::vector<Atom>> conformations{};
    conformations.reserve(numberOfConformations);
    for (std::size_t i = 0; i != numberOfConformations; ++i)
    {
      conformations.push_back(equilibratedIdealGasConformation(random, componentId));
    }
    component.setRecoilReferenceConformations(std::move(conformations));
  }
}

std::pair<Molecule, std::vector<Atom>> System::equilibratedIdealGasMoleculeRandomInBox(RandomNumber& random,
                                                                                      std::size_t selectedComponent)
{
  const Component& component = components[selectedComponent];

  // Body-frame atoms: rigid components keep their fixed reference geometry; flexible components draw an
  // internal conformation from the ideal-gas Boltzmann distribution exp(-beta * U_intra). In both cases the
  // body is placed with a uniformly random overall orientation at a uniformly random center of mass.
  const std::vector<Atom> bodyAtoms =
      component.rigid ? component.atoms : equilibratedIdealGasConformation(random, selectedComponent);

  simd_quatd q = random.randomSimdQuatd();
  double3x3 M = double3x3::buildRotationMatrixInverse(q);
  double3 com = simulationBox.randomPosition(random);

  std::vector<Atom> trial_atoms(bodyAtoms);
  for (std::size_t i = 0; i != bodyAtoms.size(); ++i)
  {
    trial_atoms[i].position = com + M * bodyAtoms[i].position;
  }

  return {{com, q, component.totalMass, selectedComponent, component.atoms.size()}, trial_atoms};
}


std::pair<std::vector<Molecule>, std::vector<Atom>> System::scaledCenterOfMassPositions(double scale) const
{
  std::vector<Molecule> scaledMolecules(moleculeData);
  // Molecule atoms only (framework atoms are not scaled by center-of-mass volume moves); this matches the
  // local Molecule::atomIndex convention and the molecule-atom spans expected by the energy routines.
  const std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
  std::vector<Atom> scaledAtoms(moleculeAtoms.begin(), moleculeAtoms.end());

  for (Molecule& molecule : scaledMolecules)
  {
    const std::span<Atom> span = {&scaledAtoms[molecule.atomIndex], molecule.numberOfAtoms};

    double totalMass = 0.0;
    double3 com(0.0, 0.0, 0.0);
    for (const Atom& atom : span)
    {
      double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
      com += mass * atom.position;
      totalMass += mass;
    }
    com /= totalMass;

    molecule.centerOfMassPosition = com * scale;

    double3 d = com * (scale - 1.0);
    for (Atom& atom : span)
    {
      atom.position += d;
    }
  }

  return {scaledMolecules, scaledAtoms};
}

std::pair<std::vector<Molecule>, std::vector<Atom>> System::scaledCenterOfMassPositions(const SimulationBox& oldBox,
                                                                                        const SimulationBox& newBox) const
{
  std::vector<Molecule> scaledMolecules(moleculeData);
  // Molecule atoms only (framework atoms are not scaled by center-of-mass volume moves); this matches the
  // local Molecule::atomIndex convention and the molecule-atom spans expected by the energy routines.
  const std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
  std::vector<Atom> scaledAtoms(moleculeAtoms.begin(), moleculeAtoms.end());

  for (Molecule& molecule : scaledMolecules)
  {
    const std::span<Atom> span = {&scaledAtoms[molecule.atomIndex], molecule.numberOfAtoms};

    double totalMass = 0.0;
    double3 com(0.0, 0.0, 0.0);
    for (const Atom& atom : span)
    {
      double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
      com += mass * atom.position;
      totalMass += mass;
    }
    com /= totalMass;

    // Scale only the center of mass through the box transformation, keeping the fractional center-of-mass
    // position fixed. Both rigid and flexible molecules are translated as a rigid body by the resulting
    // center-of-mass displacement, so internal geometry (and thus all intramolecular energies) is preserved.
    double3 newCom = newBox.cell * (oldBox.inverseCell * com);
    molecule.centerOfMassPosition = newCom;

    double3 d = newCom - com;
    for (Atom& atom : span)
    {
      atom.position += d;
    }
  }

  return {scaledMolecules, scaledAtoms};
}
