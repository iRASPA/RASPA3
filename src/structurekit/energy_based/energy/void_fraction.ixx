module;

export module energy_void_fraction;

import std;

import int3;
import double2;
import double3;
import double3x3;
import pair_interactions;
import crystal;

export struct EnergyVoidFraction
{
  EnergyVoidFraction() {};

  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnersteps);
};
