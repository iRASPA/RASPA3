module;

export module mc_moves_rotation;

import std;

import randomnumbers;
import running_energy;
import atom;
import molecule;
import component;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> rotationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                          std::size_t selectedMolecule, const std::vector<Component> &components,
                                          Molecule &molecule, std::span<Atom> molecule_atoms);
}
