module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module mc_moves;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <fstream>;
#endif

import archive;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;
import running_energy;
import mc_moves_translation;

export namespace MC_Moves
{
void performRandomMove(RandomNumber& random, System& selectedSystem, System& selectedSecondystem,
                       size_t selectedComponent, size_t& fractionalMoleculeSystem);

void performRandomMoveProduction(RandomNumber& random, System& selectedSystem, System& selectedSecondSystem,
                                 size_t selectedComponent, size_t& fractionalMoleculeSystem, size_t currentBlock);
};  // namespace MC_Moves
