module;

export module energy_surface_area;

import std;

import int3;
import double2;
import double3;
import double3x3;
import forcefield;
import framework;

export struct EnergySurfaceArea
{
  EnergySurfaceArea();

  void run(const ForceField &forceField, const Framework &framework, double isoValue,
                         std::string probePseudoAtom);
};
