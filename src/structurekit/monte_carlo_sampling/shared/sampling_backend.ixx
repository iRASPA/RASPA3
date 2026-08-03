module;

export module sampling_backend;

import std;

import double3;
import sampled_structure;

// The four things a sampler asks of a structure, in batches large enough to be worth handing to a device.
//
// Building a roadmap of the void is an algorithm about graphs and a great deal of arithmetic about spheres,
// and only the second of those is worth moving anywhere. The algorithm draws points, joins the ones a
// sphere can travel between, walks the roomiest of them uphill and then looks for ways between the pockets
// it finds -- and every one of those stages reduces to one of the four calls below, each of which is a loop
// over the atoms of the cell repeated some hundreds of thousands of times over.
//
// So the algorithm is written once, in `SampledRoadmap::build`, and what changes between the processor and
// the GPU is only which of these four is bound. That is what makes the two comparable: they are not two
// implementations of the same idea that have to be kept in step, they are the same implementation.
//
// Each call is given every case at once rather than one at a time. A device that is handed a hundred
// thousand independent segments has something to do; handed one, it spends longer being asked than
// answering.

// The narrowest point of a way from one place to another: how much room there is there, where it is, and
// which direction the way runs in there. The last two are what a window can be measured across.
export struct SampledWay
{
  double radius{0.0};  // Å, negative when no sphere can pass
  double3 position{};  // Cartesian
  double3 direction{};  // unit
};

// A place in the void and the room around it: where a walk uphill started, or where it ended.
export struct SampledPeak
{
  double radius{0.0};  // Å
  double3 position{};  // Cartesian
};

export struct SamplingBackend
{
  // How this backend is named in the files written from the roadmaps it builds, so that running both leaves
  // two sets of results side by side rather than one on top of the other.
  std::string name;

  // The signed clearance at each position: how much room there is for a sphere centred there, or how deep
  // inside the nearest atom it is when there is none.
  std::function<void(const SampledStructure &, std::span<const double3> positions, std::span<double> into)>
      clearances;

  // For each segment, the largest sphere that can travel in a straight line from `origins[i]` to
  // `origins[i] + displacements[i]`, and where along the way it is hemmed in most closely. Exact, not
  // sampled: the clearance along a segment is a minimum over atoms of a distance that is smallest somewhere
  // definite.
  std::function<void(const SampledStructure &, std::span<const double3> origins,
                     std::span<const double3> displacements, std::span<SampledWay> into)>
      straightWays;

  // From each starting place, a walk uphill of `steps` steps: a direction is tried at the current step
  // length and taken when it finds more room, and the step shrinks geometrically over the walk whether or
  // not it does, so the walk covers every scale from the free ball it starts in down to a resolution it
  // reaches on its last step. Each walk takes its directions from its own seed, so that the answer does not
  // depend on how the walks were shared out among whatever is running them.
  std::function<void(const SampledStructure &, std::span<const SampledPeak> starts,
                     std::span<const std::uint32_t> seeds, std::size_t steps, std::span<SampledPeak> into)>
      walksUphill;

  // The same segments as `straightWays`, but with the line allowed to bend. The point where the straight
  // line is hemmed in most closely is walked out sideways -- across the line rather than along it, since
  // along it the room grows towards either end and says nothing about the passage -- and if that finds a
  // roomier point, the two halves of the bent way are each put through the same again, `depth` times over.
  //
  // What comes back is a width some sphere was shown to manage, never one it was assumed to: every piece of
  // every way is measured exactly, and a way that bends more than this can find only makes the answer too
  // small.
  std::function<void(const SampledStructure &, std::span<const double3> origins,
                     std::span<const double3> displacements, std::span<const std::uint32_t> seeds,
                     std::size_t depth, std::span<SampledWay> into)>
      widestWays;
};

// How many cases are handed over at a time. Large enough that a device is busy and the cost of the call
// does not show, small enough that the arrays staged for it are megabytes rather than gigabytes.
export inline constexpr std::size_t samplingBatchSize = 1u << 18;

// How many times a bent way is bent again on its own halves, and how many directions are tried at each
// bend. Shared so that the two backends bend alike.
export inline constexpr std::size_t samplingRefinementDepth = 3;
export inline constexpr std::size_t samplingRefinementSteps = 24;
