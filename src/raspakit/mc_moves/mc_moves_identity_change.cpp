module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cmath>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_identity_change;

#ifndef USE_LEGACY_HEADERS
import <optional>;
import <span>;
import <chrono>;
import <vector>;
import <cmath>;
import <tuple>;
#endif

import randomnumbers;
import running_energy;
import system;
import atom;
import cbmc;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import move_statistics;
import component;
import mc_moves_probabilities_particles;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;

std::optional<RunningEnergy> MC_Moves::identityChangeMove([[maybe_unused]] RandomNumber& random,
                                                          [[maybe_unused]] System& system,
                                                          [[maybe_unused]] size_t selectedComponent,
                                                          [[maybe_unused]] size_t selectedMolecule,
                                                          [[maybe_unused]] std::span<Atom> atoms)
{
  return std::nullopt;
}
