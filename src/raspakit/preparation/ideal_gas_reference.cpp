module;

module ideal_gas_reference;

import std;

import threadpool;
import component;
import atom;
import cbmc;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import simulationbox;
import mc_moves_move_types;
import mc_moves_probabilities;
import units;

namespace
{

// Grows of the isolated molecule at one temperature: ln of the weight of each (the Rosenbluth weight, with the
// Ewald compensation of the intramolecular pairs taken back out as in 'System::equilibratedIdealGasConformation')
// and its intramolecular energy, of which generation density times weight is exp(-beta U) / Z_b. A grow that
// fails has the weight zero: it is counted in 'attempts' only.
struct IdealGasGrows
{
  double beta{0.0};
  std::size_t attempts{0};
  std::vector<double> logWeight{};
  std::vector<double> energy{};
  bool converged{false};

  // ln sum_j W_j exp(-betaShift U_j)
  double logSum(double betaShift) const
  {
    double largest = -1e300;
    for (std::size_t j = 0; j < logWeight.size(); ++j)
    {
      largest = std::max(largest, logWeight[j] - betaShift * energy[j]);
    }
    double sum = 0.0;
    for (std::size_t j = 0; j < logWeight.size(); ++j)
    {
      sum += std::exp(logWeight[j] - betaShift * energy[j] - largest);
    }
    return largest + std::log(sum);
  }

  double logMeanWeight() const { return logSum(0.0) - std::log(static_cast<double>(attempts)); }

  // the standard error of the mean weight over the mean weight (the grows are independent of each other)
  double relativeError() const
  {
    if (logWeight.empty() || attempts < 2uz) return std::numeric_limits<double>::infinity();
    const double largest = *std::ranges::max_element(logWeight);
    double sum = 0.0, sumOfSquares = 0.0;
    for (const double value : logWeight)
    {
      const double weight = std::exp(value - largest);
      sum += weight;
      sumOfSquares += weight * weight;
    }
    const double n = static_cast<double>(attempts);
    const double mean = sum / n;
    const double variance = std::max(0.0, sumOfSquares / n - mean * mean);
    return std::sqrt(variance / (n - 1.0)) / mean;
  }

  // the spread of the energy in the Boltzmann distribution (the grows weighted with W)
  double energySpread() const
  {
    const double largest = *std::ranges::max_element(logWeight);
    double norm = 0.0, mean = 0.0, meanSquare = 0.0;
    for (std::size_t j = 0; j < logWeight.size(); ++j)
    {
      const double weight = std::exp(logWeight[j] - largest);
      norm += weight;
      mean += weight * energy[j];
      meanSquare += weight * energy[j] * energy[j];
    }
    mean /= norm;
    return std::sqrt(std::max(0.0, meanSquare / norm - mean * mean));
  }
};

// the widest step of a free-energy perturbation, in standard deviations of beta U: wider, a temperature is put
// in between
constexpr double idealGasMaximumPerturbation = 1.0;
constexpr std::size_t idealGasMaximumBisections = 8;

// 'body(i)' for i in [0, count), on the threads of the pool and the calling thread. The first exception is thrown
// after every thread is done.
void forEachChunk(std::size_t count, const std::function<void(std::size_t)>& body)
{
  if (count == 0uz) return;
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  const std::size_t helpers = pool.getThreadingType() == ThreadPool::ThreadingType::ThreadPool
                                  ? std::min(pool.getThreadCount(), count - 1uz)
                                  : 0uz;
  if (helpers == 0uz)
  {
    for (std::size_t i = 0; i < count; ++i) body(i);
    return;
  }

  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  auto drain = [&]()
  {
    for (;;)
    {
      const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
      if (i >= count || failed.load(std::memory_order_relaxed)) return;
      try
      {
        body(i);
      }
      catch (...)
      {
        failed.store(true, std::memory_order_relaxed);
        throw;
      }
    }
  };

  std::vector<std::future<void>> pending;
  pending.reserve(helpers);
  for (std::size_t i = 0; i != helpers; ++i) pending.push_back(pool.enqueue(drain));

  std::exception_ptr error{};
  try
  {
    drain();
  }
  catch (...)
  {
    error = std::current_exception();
  }
  for (std::future<void>& task : pending)
  {
    try
    {
      task.get();
    }
    catch (...)
    {
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
}

// Grows of the isolated molecule of component 'componentId' at the temperature of 'system' (whose growth plans and
// conformation reservoirs are those of that temperature).
//
// Every chunk has a generator of its own, seeded with 'seed' and the number of the chunk, and a copy of the
// component of its own (the growth fills caches of the component as it goes): what a chunk gives does not depend on
// the thread that makes it nor on the chunks made before it.
IdealGasGrows growIdealGas(const System& system, std::size_t componentId, std::size_t seed,
                           const Preparation::IdealGasSettings& settings)
{
  const CBMC::GrowContext context = system.makeIdealGasGrowContext();

  // the growth plans and the tables they are grown with are made once, in the component the chunks copy
  Component prototype = system.components[componentId];
  {
    RandomNumber random(seed);
    [[maybe_unused]] const std::optional<CBMC::GrowResult> growData =
        CBMC::growNewMolecule(random, context, prototype, {.componentId = componentId, .moleculeId = 0uz});
  }

  struct Chunk
  {
    std::vector<double> logWeight{};
    std::vector<double> energy{};
  };

  // the tolerance that was asked for the component ('IdealGasRosenbluthWeightTolerance') comes before the default
  const double relativeTolerance = prototype.idealGasRosenbluthWeightTolerance.value_or(settings.relativeTolerance);
  const std::size_t growsPerChunk = std::max(1uz, settings.growsPerChunk);
  const std::size_t maximum = std::max(settings.minimumNumberOfGrows, settings.maximumNumberOfGrows);

  IdealGasGrows grows;
  grows.beta = system.beta;

  std::size_t target = std::max(growsPerChunk, settings.minimumNumberOfGrows);
  std::size_t numberOfChunks = 0;
  for (;;)
  {
    const std::size_t wanted = (target + growsPerChunk - 1uz) / growsPerChunk;
    std::vector<Chunk> chunks(wanted - numberOfChunks);
    forEachChunk(chunks.size(), [&](std::size_t i)
    {
      Component component = prototype;
      RandomNumber random(seed * 1000003uz + numberOfChunks + i + 1uz);
      Chunk& chunk = chunks[i];
      chunk.logWeight.reserve(growsPerChunk);
      chunk.energy.reserve(growsPerChunk);
      for (std::size_t attempt = 0; attempt < growsPerChunk; ++attempt)
      {
        const std::optional<CBMC::GrowResult> growData =
            CBMC::growNewMolecule(random, context, component, {.componentId = componentId, .moleculeId = 0uz});
        if (!growData.has_value()) continue;
        const double compensation = system.ewaldIntraCoulombCompensationEnergy(componentId, growData->atoms);
        chunk.logWeight.push_back(growData->logRosenbluthWeight + system.beta * compensation);
        chunk.energy.push_back(component.intraMolecularPotentials
                                   .computeInternalEnergies(system.forceField, system.simulationBox, growData->atoms)
                                   .potentialEnergy() -
                               compensation);
      }
    });
    for (const Chunk& chunk : chunks)
    {
      grows.logWeight.insert(grows.logWeight.end(), chunk.logWeight.begin(), chunk.logWeight.end());
      grows.energy.insert(grows.energy.end(), chunk.energy.begin(), chunk.energy.end());
    }
    numberOfChunks = wanted;
    grows.attempts = numberOfChunks * growsPerChunk;

    if (grows.logWeight.empty())
    {
      throw std::runtime_error(std::format(
          "[Ideal-gas reference]: no isolated molecule of component '{}' could be grown at {} K in {} attempts\n",
          prototype.name, system.temperature, grows.attempts));
    }
    if (grows.relativeError() <= relativeTolerance)
    {
      grows.converged = true;
      break;
    }
    if (target >= maximum) break;
    target = std::min(2uz * target, maximum);
  }
  return grows;
}

// brings 'scratch' to 'temperature': the growth plans and the conformation reservoirs (the seeds of the ring
// closures) are those of the temperature
void setTemperature(System& scratch, double temperature)
{
  if (scratch.temperature == temperature) return;
  scratch.temperature = temperature;
  scratch.beta = 1.0 / (Units::KB * temperature);
  scratch.buildConformationReservoirs();
}

// ln q(beta_b) / q(beta_a) and the disagreement of its two estimates
std::pair<double, double> logPartitionFunctionRatio(System& scratch, std::size_t componentId, const IdealGasGrows& a,
                                                    const IdealGasGrows& b, std::size_t& seed, std::size_t depth,
                                                    const Preparation::IdealGasSettings& settings)
{
  const double betaShift = b.beta - a.beta;
  if (depth < idealGasMaximumBisections &&
      std::abs(betaShift) * std::max(a.energySpread(), b.energySpread()) > idealGasMaximumPerturbation)
  {
    setTemperature(scratch, 1.0 / (Units::KB * 0.5 * (a.beta + b.beta)));
    const IdealGasGrows middle = growIdealGas(scratch, componentId, seed++, settings);
    const std::pair<double, double> lower =
        logPartitionFunctionRatio(scratch, componentId, a, middle, seed, depth + 1, settings);
    const std::pair<double, double> upper =
        logPartitionFunctionRatio(scratch, componentId, middle, b, seed, depth + 1, settings);
    return {lower.first + upper.first, lower.second + upper.second};
  }
  const double forward = a.logSum(betaShift) - a.logSum(0.0);
  const double backward = b.logSum(0.0) - b.logSum(-betaShift);
  return {0.5 * (forward + backward), std::abs(forward - backward)};
}

Preparation::IdealGasReference makeReference(double temperature, const IdealGasGrows& grows)
{
  return {.temperature = temperature,
          .logRosenbluthWeight = grows.logMeanWeight(),
          .relativeError = grows.relativeError(),
          .numberOfGrows = grows.attempts,
          .converged = grows.converged};
}

void setComputedWeight(Component& component, const Preparation::IdealGasReference& reference)
{
  component.idealGasRosenbluthWeight = std::exp(reference.logRosenbluthWeight);
  component.idealGasRosenbluthWeightEstimate =
      Component::IdealGasRosenbluthWeightEstimate{.temperature = reference.temperature,
                                                  .relativeError = reference.relativeError,
                                                  .numberOfGrows = reference.numberOfGrows,
                                                  .converged = reference.converged};
}

// a weight that was computed at another temperature than that of the system is none
void dropStaleWeights(std::span<System> systems)
{
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      if (component.idealGasRosenbluthWeightEstimate.has_value() &&
          component.idealGasRosenbluthWeightEstimate->temperature != system.temperature)
      {
        component.idealGasRosenbluthWeight.reset();
        component.idealGasRosenbluthWeightEstimate.reset();
      }
    }
  }
}

bool hasMove(const MCMoveProbabilities& probabilities, std::initializer_list<Move::Types> moves)
{
  return std::ranges::any_of(moves, [&](Move::Types move) { return probabilities.getProbability(move) > 0.0; });
}

// The moves that take other components along, whichever component (or the system) they were asked for. The pair
// and the group swap grow their molecules without the name saying so, and divide by the weight as well.
bool hasSharedGrowingMove(const System& system)
{
  const std::initializer_list<Move::Types> shared = {
      Move::Types::IdentityChangeCBMC, Move::Types::PairSwap,         Move::Types::PairSwapCBMC,
      Move::Types::PairSwapCBCFCMC,    Move::Types::GroupSwap,        Move::Types::GroupSwapCBMC,
      Move::Types::GroupSwapCBCFCMC,   Move::Types::ReactionCBMC,     Move::Types::ReactionConventionalCBCFCMC,
      Move::Types::ReactionCBCFCMC};
  return hasMove(system.mc_moves_probabilities, shared) ||
         std::ranges::any_of(system.components,
                             [&](const Component& other) { return hasMove(other.mc_moves_probabilities, shared); });
}

}  // namespace

bool Preparation::isFlexibleOpenComponent(const Component& component)
{
  return component.swappable && !component.rigid && component.atoms.size() > 1uz;
}

bool Preparation::insertionsDivideByIdealGasRosenbluthWeight(const System& system, std::size_t componentId)
{
  // the moves of the component itself
  if (hasMove(system.components[componentId].mc_moves_probabilities,
              {Move::Types::SwapCBMC, Move::Types::SwapMultiTrialCBMC, Move::Types::SwapCBCFCMC}))
  {
    return true;
  }

  return hasSharedGrowingMove(system);
}

bool Preparation::needsIdealGasRosenbluthWeight(const System& system, std::size_t componentId)
{
  const Component& component = system.components[componentId];
  if (component.rigid || component.atoms.size() <= 1uz) return false;

  // the moves of the component itself
  if (hasMove(component.mc_moves_probabilities,
              {Move::Types::SwapCBMC, Move::Types::SwapMultiTrialCBMC, Move::Types::SwapCBCFCMC, Move::Types::Widom}))
  {
    return true;
  }

  return hasSharedGrowingMove(system);
}

Preparation::IdealGasReference Preparation::idealGasRosenbluthWeight(const System& system, std::size_t componentId,
                                                                     const IdealGasSettings& settings)
{
  return makeReference(system.temperature, growIdealGas(system, componentId, 1500uz, settings));
}

std::vector<Preparation::IdealGasReference> Preparation::idealGasReferences(const System& system,
                                                                            std::size_t componentId,
                                                                            std::span<const double> temperatures,
                                                                            const IdealGasSettings& settings)
{
  std::vector<double> sorted(temperatures.begin(), temperatures.end());
  std::ranges::sort(sorted);
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

  std::vector<IdealGasReference> references(sorted.size());
  if (!sorted.empty())
  {
    // a copy: the temperature, the growth plans and the reservoirs of 'system' are left as they are
    System scratch = system;
    std::size_t seed = 1500uz;
    setTemperature(scratch, sorted.front());
    IdealGasGrows previous = growIdealGas(scratch, componentId, seed++, settings);
    references.front() = makeReference(sorted.front(), previous);
    for (std::size_t i = 1; i < sorted.size(); ++i)
    {
      setTemperature(scratch, sorted[i]);
      IdealGasGrows current = growIdealGas(scratch, componentId, seed++, settings);
      const auto [logRatio, disagreement] =
          logPartitionFunctionRatio(scratch, componentId, previous, current, seed, 0uz, settings);
      references[i] = makeReference(sorted[i], current);
      references[i].logPartitionFunction = references[i - 1].logPartitionFunction + logRatio;
      references[i].logPartitionFunctionError = references[i - 1].logPartitionFunctionError + disagreement;
      previous = std::move(current);
    }
  }

  std::vector<IdealGasReference> result;
  result.reserve(temperatures.size());
  for (const double temperature : temperatures)
  {
    result.push_back(references[static_cast<std::size_t>(
        std::distance(sorted.begin(), std::ranges::lower_bound(sorted, temperature)))]);
  }
  return result;
}

void Preparation::prepareIdealGasRosenbluthWeights(std::span<System> systems, const IdealGasSettings& settings)
{
  dropStaleWeights(systems);
  for (System& system : systems)
  {
    for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
    {
      Component& component = system.components[componentId];
      if (component.idealGasRosenbluthWeight.has_value()) continue;
      if (!needsIdealGasRosenbluthWeight(system, componentId)) continue;
      setComputedWeight(component, idealGasRosenbluthWeight(system, componentId, settings));
    }
  }
}

std::vector<std::vector<double>> Preparation::prepareIdealGasReferences(std::span<System> systems,
                                                                       std::span<const double> ladderTemperatures,
                                                                       std::span<const double> otherTemperatures,
                                                                       const IdealGasSettings& settings)
{
  const std::size_t numberOfComponents = systems.empty() ? 0uz : systems.front().components.size();
  std::vector<std::vector<double>> others(numberOfComponents, std::vector<double>(otherTemperatures.size(), 0.0));
  if (systems.empty()) return others;

  dropStaleWeights(systems);

  // the weight of each replica, from the list that was given for the ladder
  for (System& system : systems)
  {
    for (Component& component : system.components)
    {
      if (component.idealGasRosenbluthWeightLadder.empty()) continue;
      if (component.idealGasRosenbluthWeightLadder.size() != ladderTemperatures.size())
      {
        throw std::runtime_error(std::format(
            "[Ideal-gas reference]: component '{}' has {} values of 'IdealGasRosenbluthWeight' for a ladder of {} "
            "temperatures; give one per temperature (or a single number for a rigid molecule)\n",
            component.name, component.idealGasRosenbluthWeightLadder.size(), ladderTemperatures.size()));
      }
      const auto it = std::ranges::find(ladderTemperatures, system.temperature);
      if (it == ladderTemperatures.end()) continue;
      const double weight = component.idealGasRosenbluthWeightLadder[static_cast<std::size_t>(
          std::distance(ladderTemperatures.begin(), it))];
      if (!(weight > 0.0))
      {
        throw std::runtime_error(std::format(
            "[Ideal-gas reference]: the 'IdealGasRosenbluthWeight' of component '{}' at {} K must be positive\n",
            component.name, system.temperature));
      }
      component.idealGasRosenbluthWeight = weight;
      component.idealGasRosenbluthWeightEstimate.reset();
    }
  }

  std::vector<double> temperatures;
  for (const System& system : systems) temperatures.push_back(system.temperature);
  const std::size_t numberOfSystems = temperatures.size();
  temperatures.insert(temperatures.end(), otherTemperatures.begin(), otherTemperatures.end());
  const bool oneTemperature =
      std::ranges::all_of(temperatures, [&](double T) { return T == temperatures.front(); });

  for (std::size_t componentId = 0; !oneTemperature && componentId < numberOfComponents; ++componentId)
  {
    const Component& reference = systems.front().components[componentId];
    if (!isFlexibleOpenComponent(reference)) continue;
    // one number, copied to every replica of the ladder, is the weight of one temperature
    if (!ladderTemperatures.empty() && reference.idealGasRosenbluthWeightLadder.empty() &&
        reference.idealGasRosenbluthWeight.has_value() && !reference.idealGasRosenbluthWeightEstimate.has_value() &&
        std::ranges::any_of(systems, [&](const System& system) { return system.temperature != temperatures.front(); }))
    {
      throw std::runtime_error(std::format(
          "[Ideal-gas reference]: the 'IdealGasRosenbluthWeight' of the flexible component '{}' is one number for "
          "a ladder of temperatures; it depends on the temperature: give a list with one value per temperature, or "
          "leave it out to have it computed\n",
          reference.name));
    }

    const std::vector<IdealGasReference> references =
        idealGasReferences(systems.front(), componentId, temperatures, settings);
    for (std::size_t i = 0; i < numberOfSystems; ++i)
    {
      Component& component = systems[i].components[componentId];
      if (!component.idealGasRosenbluthWeight.has_value())
      {
        setComputedWeight(component, references[i]);
      }
      // ln(W_used Z_b) with Z_b = q / <W>: what the insertions of this replica divide beta f by. The insertions
      // that draw the conformation instead of growing it do not divide by the weight and sample beta f / q,
      // whatever weight was given.
      component.logInsertionNormalization = references[i].logPartitionFunction;
      if (insertionsDivideByIdealGasRosenbluthWeight(systems[i], componentId))
      {
        component.logInsertionNormalization.value() +=
            std::log(component.idealGasRosenbluthWeight.value()) - references[i].logRosenbluthWeight;
      }
    }
    for (std::size_t i = 0; i < otherTemperatures.size(); ++i)
    {
      others[componentId][i] = references[numberOfSystems + i].logPartitionFunction;
    }
  }

  // What is still without a weight: every component on a ladder of one temperature, and the components that are
  // not exchanged with a reservoir but probed or converted. The systems are replicas: those of one temperature
  // share the weight.
  for (std::size_t componentId = 0; componentId < numberOfComponents; ++componentId)
  {
    std::map<double, IdealGasReference> computed;
    for (System& system : systems)
    {
      Component& component = system.components[componentId];
      if (component.idealGasRosenbluthWeight.has_value()) continue;
      if (!needsIdealGasRosenbluthWeight(system, componentId)) continue;
      auto it = computed.find(system.temperature);
      if (it == computed.end())
      {
        it = computed.emplace(system.temperature, idealGasRosenbluthWeight(system, componentId, settings)).first;
      }
      setComputedWeight(component, it->second);
    }
  }
  return others;
}
