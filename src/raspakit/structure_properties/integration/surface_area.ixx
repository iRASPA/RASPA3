module;

export module integration_surface_area;

import std;

import framework;
import forcefield;

export struct Integration_SurfaceArea
{
  std::vector<double> data;

  Integration_SurfaceArea() {};

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::string probePseudoAtom, std::optional<std::size_t> numberOfSlices) const;
};
