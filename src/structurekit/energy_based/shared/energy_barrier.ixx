module;

export module energy_shared_energy_barrier;

import std;

import int3;
import uint3;
import crystal;
import pair_interactions;
import energy_shared_probe_energy_grid;
import energy_shared_energy_backend;

// The energetic counterpart of the pore diameters, read off the probe's energy field.
//
// The geometric route asks how wide a passage is and answers with the largest sphere that fits through it,
// which is Df. That question has a sharp edge to it: a sphere either fits or it does not, and a window a
// hair too narrow is as closed as a wall. A molecule does not see it that way. It crosses a window whose
// barrier is a few kT many times a second and never crosses one of eighty kT, and the difference between
// those two is nowhere in the geometry.
//
// So the same sweep is run on the energy instead. Points are switched on from the deepest upwards and joined
// to the neighbours already on, and the moment a region closes on a periodic image of itself the sweep has
// found the lowest energy at which a connected path runs through the crystal. That is the barrier: the
// saddle of the landscape, the least the probe must be raised to get from a cell to its neighbour.
//
// It is a global answer, which is what makes it worth having. A transition-state search or an elastic band
// returns the saddle on the path it was started near, and finding the lowest one means guessing every path
// worth trying. The sweep passes every point of the field in one pass and cannot miss the lowest, subject
// only to the grid resolving the window it goes through.
export struct EnergyBarrier
{
  // The lowest energy at which a path runs through the crystal: what the probe must be raised to in order to
  // get from one cell to the next. In the force field's internal units.
  double percolationBarrier{0.0};

  // The deepest point of the landscape anywhere, which is the strongest site the grid resolves.
  double deepestWell{0.0};

  // The deepest point of the region that first ran through, so the strongest site on the percolating
  // network rather than one in a pocket off it.
  double deepestWellOnPath{0.0};

  // The barrier to running in at least one, two, and three directions. Highest double where the pore system
  // never runs that far, which cannot happen once every point takes part.
  std::array<double, 3> barrierByDimension{};

  bool percolates{false};
  int dimensionalityAtBarrier{0};

  // The count above is read off at exactly the barrier, and a framework whose channels are alike in more
  // than one direction opens the second of them a hair above the first: DDR takes a further eight ten-
  // thousandths of a Kelvin to run in two directions rather than one. Quoting the count alone then makes a
  // two-dimensional network look one-dimensional. This gives the number of directions that are open within
  // the energy handed in, which at any sane temperature is the dimensionality anyone means, together with
  // what the last of them cost above the barrier.
  std::pair<int, double> dimensionalityWithin(double margin) const;

  // The temperature the barrier is reported against. It enters nothing that is computed, only the reading of
  // it: a barrier is worth knowing in units of kT.
  double temperature{0.0};

  std::size_t numberOfVoxels{0};
  double sweepSeconds{0.0};

  ProbeEnergyGrid grid;

  EnergyBarrier();
  ~EnergyBarrier();

  // Sweeps an energy field and reads the barrier off it. The field is any energy over the grid, so a route
  // that arrives at one differently, by averaging over a molecule's orientations for instance, gets its
  // barrier from here too rather than from a sweep of its own.
  static EnergyBarrier fromField(uint3 gridSize, std::span<const float> energy, double temperature);

  void run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
           std::string probePseudoAtom, uint3 gridSize, double temperature);
  void run(const PairInteractions &interactions, const Crystal &framework, const ProbeEnergyGrid &energyGrid,
           double temperature);
};
