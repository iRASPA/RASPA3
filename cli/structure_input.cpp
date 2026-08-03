module;

module structure_input;

import std;

import double3;
import simulationbox;
import framework;
import forcefield;
import vdwparameters;
import atom;

import unit_cell;
import crystal;
import pair_interactions;

namespace StructureInput
{

UnitCell makeUnitCell(const SimulationBox& simulationBox)
{
  UnitCell::Type type =
      simulationBox.type == SimulationBox::Type::Rectangular ? UnitCell::Type::Rectangular : UnitCell::Type::Triclinic;

  return UnitCell(simulationBox.cell, type);
}

Crystal makeCrystal(const Framework& framework)
{
  Crystal crystal{};

  crystal.name = framework.name;
  crystal.spaceGroupHallNumber = framework.spaceGroupHallNumber;
  crystal.unitCell = makeUnitCell(framework.simulationBox);
  crystal.mass = framework.unitCellMass;

  crystal.atoms.reserve(framework.unitCellAtoms.size());
  for (const Atom& atom : framework.unitCellAtoms)
  {
    crystal.atoms.push_back(
        CrystalAtom{.position = atom.position, .charge = atom.charge, .type = static_cast<std::size_t>(atom.type)});
  }

  crystal.fractionalPositions = framework.fractionalAtomPositionsUnitCell();

  return crystal;
}

PairInteractions makeInteractions(const ForceField& forceField)
{
  PairInteractions interactions{};

  std::size_t numberOfTypes = forceField.numberOfPseudoAtoms;
  interactions.numberOfTypes = numberOfTypes;

  interactions.names.reserve(numberOfTypes);
  interactions.charges.reserve(numberOfTypes);
  for (std::size_t type = 0; type < numberOfTypes; ++type)
  {
    interactions.names.push_back(forceField.pseudoAtoms[type].name);
    interactions.charges.push_back(forceField.pseudoAtoms[type].charge);
  }

  interactions.parameters.resize(numberOfTypes * numberOfTypes);
  for (std::size_t row = 0; row < numberOfTypes; ++row)
  {
    for (std::size_t column = 0; column < numberOfTypes; ++column)
    {
      const VDWParameters& parameters = forceField(row, column);
      interactions.parameters[row * numberOfTypes + column] =
          PairParameters{.sizeParameter = parameters.sizeParameter(),
                         .strengthParameter = parameters.strengthParameter(),
                         .shift = parameters.shift};
    }
  }

  interactions.cutOffVDW = forceField.cutOffFrameworkVDW;
  interactions.cutOffCoulomb = forceField.cutOffCoulomb;

  return interactions;
}

}  // namespace StructureInput
