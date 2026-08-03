module;

module mc_pore_volume;

import std;

import sampled_structure;
import sampled_roadmap;
import mc_backend;

void MC_PoreVolume::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                        std::optional<std::size_t> numberOfInnerSteps)
{
  run(structure, SampledRoadmap::build(structure, samplingBackendCPU(), numberOfIterations, numberOfInnerSteps));
}

void MC_PoreVolume::run(const SampledStructure &structure, const SampledRoadmap &roadmap)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->roadmap = roadmap;

  // The fraction, and how well the sample knows it. Nothing beyond the throwing stage went into either.
  this->voidFraction = roadmap.voidFraction;
  this->voidFractionError =
      roadmap.numberOfSamples > 0
          ? std::sqrt(this->voidFraction * (1.0 - this->voidFraction) / static_cast<double>(roadmap.numberOfSamples))
          : 0.0;

  this->numberOfChannels = roadmap.numberOfChannels;
  this->numberOfPockets = roadmap.numberOfPockets;

  // Only the sample points stand for volume. The pocket centres are nodes at places the sample already
  // covers, put there to carry the long hops, and counting them would count that volume twice over.
  //
  // A point in a piece too thinly sampled to call either way is counted in neither share, which is why the
  // two need not add up to the void fraction. That is the honest reading: it is room whose connection to
  // the rest of the crystal this sample did not settle.
  std::size_t reachable = 0;
  std::size_t shutIn = 0;
  std::size_t volumeNodes = roadmap.numberOfVolumeNodes();
  for (std::size_t node = 0; node < volumeNodes; ++node)
  {
    if (roadmap.isReachable(node))
      ++reachable;
    else if (roadmap.isShutIn(node))
      ++shutIn;
  }

  double perNode = volumeNodes > 0 ? this->voidFraction / static_cast<double>(volumeNodes) : 0.0;
  this->accessibleVolumeFraction = perNode * static_cast<double>(reachable);
  this->inaccessibleVolumeFraction = perNode * static_cast<double>(shutIn);
  this->unsettledVolumeFraction =
      this->voidFraction - this->accessibleVolumeFraction - this->inaccessibleVolumeFraction;

  double volume = structure.unitCell.volume;
  this->voidVolume = this->voidFraction * volume;
  this->accessibleVolume = this->accessibleVolumeFraction * volume;
  this->inaccessibleVolume = this->inaccessibleVolumeFraction * volume;

  double toGravimetric = structure.gravimetricVolumeFactor();
  this->gravimetricVoidVolume = this->voidVolume * toGravimetric;
  this->gravimetricAccessibleVolume = this->accessibleVolume * toGravimetric;
  this->gravimetricInaccessibleVolume = this->inaccessibleVolume * toGravimetric;

  this->density = structure.density();

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(std::format("{}.mc.av.{}.txt", structure.name, roadmap.backend));
  std::print(myfile, "# Void fraction and pore volume by sampling the void\n");
  structure.writeHeader(myfile);
  roadmap.writeHeader(myfile);
  std::print(myfile, "# Timing (volume, on the processor either way): {} [s]\n", this->seconds);
  std::print(myfile, "# The fraction is a binomial proportion over the points thrown and is unbiased; the\n");
  std::print(myfile, "# error beside it is its standard error. The split into reachable and shut-away is not\n");
  std::print(myfile, "# unbiased: a channel the sample failed to join up looks like a run of pockets, so the\n");
  std::print(myfile, "# accessible share is a lower bound and the inaccessible one an upper bound. Many more\n");
  std::print(myfile, "# pockets than the structure has is what too thin a sample looks like. The unsettled\n");
  std::print(myfile, "# share is room in pieces too thinly sampled to place either way, and it belongs to\n");
  std::print(myfile, "# neither of the other two rather than being guessed into one of them.\n");
  std::print(myfile, "Void fraction:                {} +/- {}\n", this->voidFraction, this->voidFractionError);
  std::print(myfile, "  of which accessible:        {}\n", this->accessibleVolumeFraction);
  std::print(myfile, "  of which inaccessible:      {}\n", this->inaccessibleVolumeFraction);
  std::print(myfile, "  of which unsettled:         {}\n", this->unsettledVolumeFraction);
  std::print(myfile, "Void volume:                  {} [Å³]  {} [cm³/g]\n", this->voidVolume,
             this->gravimetricVoidVolume);
  std::print(myfile, "  accessible volume:          {} [Å³]  {} [cm³/g]\n", this->accessibleVolume,
             this->gravimetricAccessibleVolume);
  std::print(myfile, "  inaccessible volume:        {} [Å³]  {} [cm³/g]\n", this->inaccessibleVolume,
             this->gravimetricInaccessibleVolume);
  std::print(myfile, "Crystal density:            {} [kg/m³]\n", this->density);
  std::print(myfile, "Number of channels: {}\n", this->numberOfChannels);
  std::print(myfile, "Number of pockets:  {}\n", this->numberOfPockets);
  myfile.close();
}
