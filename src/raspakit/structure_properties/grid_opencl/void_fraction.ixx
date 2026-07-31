module;

export module opencl_void_fraction;

import std;

import uint3;
import double3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;

// The void a probe has to itself, and how much of it the probe can actually get into.
//
// The void is the space where the probe's centre fits at all, which on the clearance field is the points at
// least the probe's radius from every atom surface. It divides in two by connectivity alone: the part of it
// lying in a channel is the volume a molecule reaches by walking in from outside the crystal, and the part
// lying in a pocket is sealed off behind a neck too narrow for it and holds nothing at equilibrium however
// much room is in there.
//
// Both parts come from counting grid points, so the answer is a volume fraction directly and its error is the
// error of the boundary between them, roughly one grid spacing over the width of the pore.
export struct GridVoidFraction
{
  double probeRadius{0.0};

  double voidFraction{0.0};
  double accessibleVolumeFraction{0.0};
  double inaccessibleVolumeFraction{0.0};

  double voidVolume{0.0};
  double accessibleVolume{0.0};
  double inaccessibleVolume{0.0};

  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};

  GridVoidFraction();
  ~GridVoidFraction();

  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
