module;

export module mc_void_fraction;

import std;

import framework;
import forcefield;

export struct MC_VoidFraction
{
  std::vector<double> data;

  MC_VoidFraction() {};

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor, std::string probePseudoAtom, std::optional<std::size_t> numberOfIterations);
};
