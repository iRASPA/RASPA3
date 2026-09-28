module;

export module replica_round_trips;

import std;

import archive;
import json;

/**
 * \brief Round-trip diagnostic for one-dimensional replica exchange (temperature ladders).
 *
 * The per-pair swap acceptance only shows whether neighboring replicas can exchange; it does not
 * show whether configurations actually travel through the whole ladder. This tracker follows the
 * identity of every configuration ('walker'; walker k starts at temperature index k) as accepted
 * swaps move it up and down the ladder and records
 *
 *   - the number of round trips: a walker completes a round trip when it returns to the coldest
 *     replica after having visited the hottest one since it last left the coldest (Katzgraber,
 *     Trebst, Huse & Troyer, J. Stat. Mech. (2006) P03018),
 *   - the mean round-trip time (in swap sweeps) and the round trips of every individual walker,
 *   - the fraction f(T) of walkers at every temperature that are 'up-moving' (last extreme visited
 *     was the coldest replica). For an optimal ladder f(T) decreases linearly from 1 at the coldest
 *     to 0 at the hottest temperature; a plateau marks a bottleneck.
 *
 * Call 'recordAcceptedSwap' for every accepted exchange and 'endOfSweep' once per swap sweep.
 */
export struct ReplicaRoundTrips
{
  enum class Direction : std::size_t
  {
    Unlabeled = 0,  ///< Walker has not visited either end of the ladder yet.
    Up = 1,         ///< Last extreme visited was the coldest replica (moving up).
    Down = 2        ///< Last extreme visited was the hottest replica (moving down).
  };

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::size_t numberOfReplicas{0};  ///< Number of temperatures (replicas) in the ladder.
  std::size_t sweeps{0};            ///< Number of swap sweeps recorded.

  std::vector<std::size_t> walkerAtReplica;  ///< Walker (configuration identity) currently held by replica k.
  std::vector<std::size_t> direction;        ///< Direction label of every walker (Direction as size_t).
  std::vector<std::size_t> leftBottomSweep;  ///< Sweep at which every walker was last at the coldest replica.

  std::size_t roundTrips{0};                       ///< Total number of completed round trips.
  std::size_t roundTripSweepsTotal{0};             ///< Summed duration (in sweeps) of the completed round trips.
  std::vector<std::size_t> roundTripsPerWalker;    ///< Completed round trips of every walker.
  std::vector<std::size_t> upVisitsPerReplica;     ///< Sweeps in which the walker at replica k was labeled Up.
  std::vector<std::size_t> downVisitsPerReplica;   ///< Sweeps in which the walker at replica k was labeled Down.

  /// Resets the tracker for a ladder of 'n' replicas; walker k starts at replica k.
  void initialize(std::size_t n);

  /// Records an accepted exchange of the configurations held by replicas 'replicaA' and 'replicaB'.
  void recordAcceptedSwap(std::size_t replicaA, std::size_t replicaB);

  /// Updates the direction labels and the f(T) counters after a complete swap sweep.
  void endOfSweep();

  /// Fraction of labeled walkers at replica k that are moving up, or nullopt when no label was seen yet.
  std::optional<double> upFraction(std::size_t replica) const;

  /// Mean number of sweeps per completed round trip (0 when none completed).
  double meanRoundTripSweeps() const;

  /// Human-readable statistics block for the combined output file.
  std::string writeStatistics(const std::vector<double>& temperatures, std::size_t swapEvery) const;

  /// JSON statistics block.
  nlohmann::json jsonStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ReplicaRoundTrips& r);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ReplicaRoundTrips& r);
};
