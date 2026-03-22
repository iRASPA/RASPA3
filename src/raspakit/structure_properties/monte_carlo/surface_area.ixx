module;

export module mc_surface_area;

import std;

import framework;
import forcefield;

export struct MC_SurfaceArea
{
  std::vector<double> data;

  MC_SurfaceArea() {};

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::string probePseudoAtom, std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps) const;
};
