module;

module replica_round_trips;

import std;

import archive;
import json;

void ReplicaRoundTrips::initialize(std::size_t n)
{
  numberOfReplicas = n;
  sweeps = 0uz;

  walkerAtReplica.resize(n);
  std::iota(walkerAtReplica.begin(), walkerAtReplica.end(), 0uz);
  direction.assign(n, static_cast<std::size_t>(Direction::Unlabeled));
  leftBottomSweep.assign(n, 0uz);

  roundTrips = 0uz;
  roundTripSweepsTotal = 0uz;
  roundTripsPerWalker.assign(n, 0uz);
  upVisitsPerReplica.assign(n, 0uz);
  downVisitsPerReplica.assign(n, 0uz);
}

void ReplicaRoundTrips::recordAcceptedSwap(std::size_t replicaA, std::size_t replicaB)
{
  std::swap(walkerAtReplica[replicaA], walkerAtReplica[replicaB]);
}

void ReplicaRoundTrips::endOfSweep()
{
  if (numberOfReplicas < 2uz) return;

  ++sweeps;

  // label the walkers at the two ends of the ladder; a walker arriving at the coldest replica while
  // labeled Down (it has visited the hottest replica since it last left the coldest) closes a round trip
  {
    const std::size_t walker = walkerAtReplica.front();
    if (direction[walker] == static_cast<std::size_t>(Direction::Down))
    {
      ++roundTrips;
      ++roundTripsPerWalker[walker];
      roundTripSweepsTotal += sweeps - leftBottomSweep[walker];
    }
    direction[walker] = static_cast<std::size_t>(Direction::Up);
    leftBottomSweep[walker] = sweeps;
  }
  {
    const std::size_t walker = walkerAtReplica.back();
    direction[walker] = static_cast<std::size_t>(Direction::Down);
  }

  // f(T) histogram: at every temperature, count whether the resident walker moves up or down
  for (std::size_t replica = 0; replica < numberOfReplicas; ++replica)
  {
    switch (static_cast<Direction>(direction[walkerAtReplica[replica]]))
    {
      case Direction::Up:
        ++upVisitsPerReplica[replica];
        break;
      case Direction::Down:
        ++downVisitsPerReplica[replica];
        break;
      case Direction::Unlabeled:
        break;
    }
  }
}

std::optional<double> ReplicaRoundTrips::upFraction(std::size_t replica) const
{
  const std::size_t labeled = upVisitsPerReplica[replica] + downVisitsPerReplica[replica];
  if (labeled == 0uz) return std::nullopt;
  return static_cast<double>(upVisitsPerReplica[replica]) / static_cast<double>(labeled);
}

double ReplicaRoundTrips::meanRoundTripSweeps() const
{
  if (roundTrips == 0uz) return 0.0;
  return static_cast<double>(roundTripSweepsTotal) / static_cast<double>(roundTrips);
}

std::string ReplicaRoundTrips::writeStatistics(const std::vector<double>& temperatures, std::size_t swapEvery) const
{
  std::ostringstream stream;

  std::print(stream, "Replica round-trip statistics (configuration diffusion through the ladder)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    a round trip is completed when a configuration returns to the coldest replica\n");
  std::print(stream, "    after having visited the hottest one; f(T) is the fraction of configurations at T\n");
  std::print(stream, "    that last visited the coldest replica (ideal ladder: linear from 1 to 0)\n\n");

  std::print(stream, "    sweeps recorded:          {}\n", sweeps);
  std::print(stream, "    round trips:              {}\n", roundTrips);
  if (roundTrips > 0uz)
  {
    std::print(stream, "    mean round-trip time:     {:.1f} sweeps ({:.0f} cycles)\n", meanRoundTripSweeps(),
               meanRoundTripSweeps() * static_cast<double>(swapEvery));
    std::print(stream, "    round trips per replica:  {:.3f} (per {} cycles)\n",
               static_cast<double>(roundTrips) / static_cast<double>(std::max(1uz, numberOfReplicas)),
               sweeps * swapEvery);
  }
  else
  {
    std::print(stream, "    mean round-trip time:     n/a (no round trip completed; the run is too short for\n");
    std::print(stream, "                              configurations to traverse the ladder, or a pair blocks it)\n");
  }
  std::print(stream, "\n");

  std::print(stream, "    replica    temperature [K]    f(T) up-fraction    labeled sweeps\n");
  std::print(stream, "    -----------------------------------------------------------------\n");
  for (std::size_t replica = 0; replica < numberOfReplicas; ++replica)
  {
    const std::optional<double> f = upFraction(replica);
    if (f.has_value())
    {
      std::print(stream, "    {:4d}       {:10.4f}         {:8.4f}            {:9d}\n", replica, temperatures[replica],
                 f.value(), upVisitsPerReplica[replica] + downVisitsPerReplica[replica]);
    }
    else
    {
      std::print(stream, "    {:4d}       {:10.4f}              n/a            {:9d}\n", replica,
                 temperatures[replica], 0uz);
    }
  }
  std::print(stream, "\n");

  std::print(stream, "    walker (initial replica)    round trips    current replica\n");
  std::print(stream, "    ---------------------------------------------------------\n");
  std::vector<std::size_t> replicaOfWalker(numberOfReplicas);
  for (std::size_t replica = 0; replica < numberOfReplicas; ++replica)
  {
    replicaOfWalker[walkerAtReplica[replica]] = replica;
  }
  for (std::size_t walker = 0; walker < numberOfReplicas; ++walker)
  {
    std::print(stream, "    {:6d}                      {:9d}        {:6d}\n", walker, roundTripsPerWalker[walker],
               replicaOfWalker[walker]);
  }
  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json ReplicaRoundTrips::jsonStatistics() const
{
  nlohmann::json j;
  j["sweeps"] = sweeps;
  j["roundTrips"] = roundTrips;
  j["meanRoundTripSweeps"] = meanRoundTripSweeps();
  j["roundTripsPerWalker"] = roundTripsPerWalker;
  j["walkerAtReplica"] = walkerAtReplica;

  std::vector<double> upFractions;
  upFractions.reserve(numberOfReplicas);
  for (std::size_t replica = 0; replica < numberOfReplicas; ++replica)
  {
    upFractions.push_back(upFraction(replica).value_or(std::numeric_limits<double>::quiet_NaN()));
  }
  j["upFraction"] = upFractions;
  j["upVisitsPerReplica"] = upVisitsPerReplica;
  j["downVisitsPerReplica"] = downVisitsPerReplica;
  return j;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ReplicaRoundTrips& r)
{
  archive << r.versionNumber;

  archive << r.numberOfReplicas;
  archive << r.sweeps;
  archive << r.walkerAtReplica;
  archive << r.direction;
  archive << r.leftBottomSweep;
  archive << r.roundTrips;
  archive << r.roundTripSweepsTotal;
  archive << r.roundTripsPerWalker;
  archive << r.upVisitsPerReplica;
  archive << r.downVisitsPerReplica;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ReplicaRoundTrips& r)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > r.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ReplicaRoundTrips' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> r.numberOfReplicas;
  archive >> r.sweeps;
  archive >> r.walkerAtReplica;
  archive >> r.direction;
  archive >> r.leftBottomSweep;
  archive >> r.roundTrips;
  archive >> r.roundTripSweepsTotal;
  archive >> r.roundTripsPerWalker;
  archive >> r.upVisitsPerReplica;
  archive >> r.downVisitsPerReplica;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error("ReplicaRoundTrips: error in binary restart\n");
  }

  return archive;
}
