module;

module system;

import std;

import randomnumbers;
import simd_quatd;
import double3;
import double3x3;
import atom;
import atom_dynamics;
import component;
import molecule;
import framework;
import forcefield;
import simulationbox;
import units;
import cbmc;
import cbmc_chain_data;
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
  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, moleculeIndex);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, moleculeIndex);
  moleculeData.insert(moleculeIterator, molecule);

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricField.resize(electricField.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfFractionalMoleculesPerComponent[selectedComponent] += 1;

  netCharge += components[selectedComponent].netCharge;
  netChargeAdsorbates += components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] += components[selectedComponent].netCharge;

  for (Atom& atom : atoms)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
  }

  updateMoleculeAtomInformation();
}

void System::insertFractionalMolecule(std::size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                                      std::vector<Atom> atoms, std::size_t moleculeIndex)
{
  const double l = 0.0;
  const std::uint8_t groupId = fractionalSlotDUdlambdaGroupId(selectedComponent, moleculeIndex);
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint16_t>(moleculeIndex);
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
  std::vector<Atom>::const_iterator iterator =
      iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});

  std::vector<Molecule>::iterator moleculeIterator =
      indexForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  moleculeData.insert(moleculeIterator, molecule);

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricField.resize(electricField.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

  netCharge += components[selectedComponent].netCharge;
  netChargeAdsorbates += components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] += components[selectedComponent].netCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint16_t>(numberOfMoleculesPerComponent[selectedComponent]);
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
  }

  updateMoleculeAtomInformation();
}

void System::insertMoleculePolarization(std::size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                                        std::vector<Atom> atoms, std::span<double3> electric_field)
{
  std::vector<Atom>::const_iterator iterator =
      iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
  atomData.insert(iterator, atoms.begin(), atoms.end());
  atomDynamics.insert(atomDynamics.cbegin() + atomOffset, atoms.size(), AtomDynamics{});

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

  netCharge += components[selectedComponent].netCharge;
  netChargeAdsorbates += components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] += components[selectedComponent].netCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<std::uint16_t>(numberOfMoleculesPerComponent[selectedComponent]);
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] += 1;
  }

  updateMoleculeAtomInformation();
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
  for (const Atom& atom : molecule)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] -= 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] -= 1;
  }

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
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

  netCharge -= components[selectedComponent].netCharge;
  netChargeAdsorbates -= components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] -= components[selectedComponent].netCharge;

  updateMoleculeAtomInformation();
}

void System::deleteMolecule(std::size_t selectedComponent, std::size_t selectedMolecule, const std::span<Atom> molecule)
{
  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (const Atom& atom : molecule)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<std::size_t>(atom.type)] -= 1;
    totalNumberOfPseudoAtoms[static_cast<std::size_t>(atom.type)] -= 1;
  }

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
  std::vector<Atom>::difference_type atomOffset = iterator - atomData.cbegin();
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

  netCharge -= components[selectedComponent].netCharge;
  netChargeAdsorbates -= components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] -= components[selectedComponent].netCharge;

  translationalDegreesOfFreedom -= components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom -= components[selectedComponent].rotationalDegreesOfFreedom;

  updateMoleculeAtomInformation();
}

void System::updateMoleculeAtomInformation()
{
  std::size_t atom_index = numberOfFrameworkAtoms;
  std::size_t molecule_index{};

  for (std::size_t k = 0; k < components.size(); k++)
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
  RandomNumber random(std::time(0l));

  for (std::size_t componentId = 0; const Component& component : components)
  {
    if (component.hasFractionalMolecule)
    {
      numberOfMoleculesPerComponent[componentId] = 0;

      auto growFractionalMolecule = [&](std::uint8_t groupId) -> std::optional<ChainGrowData>
      {
        std::optional<ChainGrowData> growData = std::nullopt;
        do
        {
          const Component::GrowType growType = components[componentId].growType;
          growData = CBMC::growMoleculeSwapInsertion(
              random,
              CBMC::GrowContext{hasExternalField, forceField, simulationBox, interpolationGrids,
                                externalFieldInterpolationGrid, framework, spanOfFrameworkAtoms(),
                                spanOfMoleculeAtoms(), beta, forceField.cutOffFrameworkVDW,
                                forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb},
              components[componentId], componentId, growType, numberOfMolecules(), 0.0, groupId, true);
        } while (!growData || growData->energies.potentialEnergy() > forceField.energyOverlapCriteria);
        return growData;
      };

      if (numberOfGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfGCFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<ChainGrowData> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairGCFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<ChainGrowData> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<ChainGrowData> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<ChainGrowData> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::size_t slot = indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(componentId);
        const std::optional<ChainGrowData> growData =
            growFractionalMolecule(fractionalSlotDUdlambdaGroupId(componentId, slot));
        insertFractionalMolecule(componentId, growData->molecule, growData->atoms, slot);
      }

      if (numberOfGibbsFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
      {
        const std::uint8_t groupId = components[componentId].lambdaGibbs.dUdlambdaGroupId;
        const std::optional<ChainGrowData> growData = growFractionalMolecule(groupId);
        std::vector<Atom> atoms = growData->atoms;
        for (Atom& atom : atoms)
        {
          atom.setScalingToFractional(0.5, groupId);
        }
        insertFractionalMoleculeAtIndex(
            componentId, indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(componentId),
            growData->molecule, std::move(atoms));
      }

      translationalDegreesOfFreedom += components[componentId].translationalDegreesOfFreedom;
      rotationalDegreesOfFreedom += components[componentId].rotationalDegreesOfFreedom;
    }

    // Add atom position passed to 'createInitialMolecules'
    if (componentId < initialPositions.size())
    {
      std::vector<double3> positions = initialPositions[componentId];
      for (std::size_t i = 0; i < positions.size(); i += component.atoms.size())
      {
        std::span<double3> position_view = std::span<double3>(&positions[i], component.atoms.size());

        double3 com(0.0, 0.0, 0.0);
        simd_quatd orientation{};

        if (components[componentId].rigid)
        {
          // Compute center of mass
          double totalMass = 0.0;
          for (std::size_t k = 0; k < position_view.size(); ++k)
          {
            double mass = forceField.pseudoAtoms[static_cast<std::size_t>(components[componentId].atoms[k].type)].mass;
            com += mass * position_view[k];
            totalMass += mass;
          }
          com /= totalMass;

          double3 reference_com = double3{0.0, 0.0, 0.0};

          std::vector<double3> reference_positions =
              components[componentId].atoms | std::views::transform(&Atom::position) | std::ranges::to<std::vector>();

          // Compute rotation matrix that describes going from the space-fixed frame to the body-fixed frame
          double3x3 rotation_matrix =
              double3x3::computeRotationMatrix(com, position_view, reference_com, reference_positions);

          // Get orientation
          orientation = rotation_matrix.quaternion();
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
      std::optional<ChainGrowData> growData = std::nullopt;
      bool inside_blocked_pocket{false};
      do
      {
        do
        {
          Component::GrowType growType = components[componentId].growType;
          growData = CBMC::growMoleculeSwapInsertion(
              random,
              CBMC::GrowContext{hasExternalField, forceField, simulationBox, interpolationGrids,
                                externalFieldInterpolationGrid, framework, spanOfFrameworkAtoms(),
                                spanOfMoleculeAtoms(), beta, forceField.cutOffFrameworkVDW,
                                forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb},
              components[componentId], componentId, growType, numberOfMolecules(), 1.0, false, false);

        } while (!growData || growData->energies.potentialEnergy() > forceField.energyOverlapCriteria);

        std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());
        inside_blocked_pocket = insideBlockedPockets(components[componentId], newMolecule);

      } while (inside_blocked_pocket);

      insertMolecule(componentId, growData->molecule, growData->atoms);
    }

    componentId++;
  }
}

bool System::insideBlockedPockets(const Component& component, std::span<const Atom> molecule_atoms) const
{
  if (framework.has_value())
  {
    for (std::size_t i = 0; i != component.blockingPockets.size(); ++i)
    {
      double radius_squared = component.blockingPockets[i].w * component.blockingPockets[i].w;
      double3 pos =
          framework->simulationBox.cell *
          double3(component.blockingPockets[i].x, component.blockingPockets[i].y, component.blockingPockets[i].z);
      for (const Atom& atom : molecule_atoms)
      {
        double vdwScaling = atom.scalingVDW;
        double3 dr = atom.position - pos;
        dr = framework->simulationBox.applyPeriodicBoundaryConditions(dr);
        if (dr.length_squared() < vdwScaling * radius_squared)
        {
          return true;
        }
      }
    }
  }
  return false;
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
std::pair<std::vector<Molecule>, std::vector<Atom>> System::scaledCenterOfMassPositions(double scale) const
{
  std::vector<Molecule> scaledMolecules(moleculeData);
  std::vector<Atom> scaledAtoms(atomData);

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
  std::vector<Atom> scaledAtoms(atomData);

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
