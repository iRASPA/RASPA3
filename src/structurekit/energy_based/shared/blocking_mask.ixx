module;

export module energy_shared_blocking_mask;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_energy_backend;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;

// Blocking pockets taken from the exact surface-area route, applied to an energy landscape.
//
// The pockets themselves come from `VoronoiBlockingSpheres`: the surface-area sweep finds each inaccessible
// cage, the divergence theorem gives its volume and centroid, and a sphere is written that holds the cage
// without reaching a channel. Sampling is the fallback the same route already takes when a measured split
// has to be refused. Nothing here decides what a pocket is; it only asks that route and then writes the
// answer onto a grid.
//
// A point inside a sphere is pore the probe is not allowed into. Filling it with one flat overlap energy
// leaves marching cubes nothing to interpolate along an edge that is inside at both ends, and the rim of
// the sphere comes out as a staircase on the grid planes. A ramp, the energy rising with depth at
// `energyPerAngstrom`, gives the seam a gradient of the same order as the framework walls, which is what
// makes the level set follow the sphere. The iRASPA well-surface construction uses the same convention, at
// a thousand kelvin per angstrom, so that the energy half of a pocket and its distance half stay the same
// function of depth.
//
// The arithmetic here is a post-pass on a field that has already been built, which is why it can sit on
// either backend: the GPU and the processor both return a `ProbeEnergyGrid` or a `MolecularEnergyGrid`, and
// the spheres are written onto those. The well-field builder next door bakes the same spheres in while it
// sums, because it has a distance channel the energy grids do not.

export inline constexpr double blockedEnergyPerAngstromInKelvin = 1000.0;

// How far a point is, in angstrom, from the surface of the nearest blocking sphere: negative inside one,
// and a large positive number when there are no spheres. The spheres repeat with the unit cell.
export double blockingSphereDistance(double3 fractional, const UnitCell &unitCell,
                                     std::span<const BlockingSphere> spheres);

export void applyBlockingSpheresToEnergy(std::span<float> energy, uint3 gridSize, const UnitCell &unitCell,
                                         std::span<const BlockingSphere> spheres, double energyPerAngstrom,
                                         double ceiling);

// The spheres the exact surface-area route writes for this probe, measured where it can and sampled where
// it cannot. Empty when the probe is not in the force field or the structure has no inaccessible cages.
export std::vector<BlockingSphere> analysisBlockingSpheres(const PairInteractions &interactions,
                                                           const Crystal &framework, const std::string &probePseudoAtom);

// The same backend, with the spheres written onto every energy field it returns. An empty list of spheres
// is a no-op, so a caller that has already decided not to block can still hand the result of this down.
export EnergyBackend withBlockedPockets(EnergyBackend backend, std::span<const BlockingSphere> spheres,
                                        double energyPerAngstrom, double ceiling);
