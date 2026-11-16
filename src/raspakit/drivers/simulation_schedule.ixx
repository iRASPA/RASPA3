module;

export module simulation_schedule;

import std;

/**
 * \brief Run-control settings shared by the Monte Carlo and Molecular Dynamics drivers.
 *
 * SimulationSchedule groups the cycle counts of the individual simulation stages together with the
 * intervals at which periodic bookkeeping actions are performed. It is passed to the simulation
 * driver constructors as a single parameter object.
 */
export struct SimulationSchedule
{
  std::size_t numberOfProductionCycles{0};                   ///< Number of production cycles.
  std::size_t numberOfPreInitializationCycles{0};  ///< Number of pre-initialization cycles.
  std::size_t numberOfInitializationCycles{0};     ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles{0};      ///< Number of equilibration cycles.

  std::size_t printEvery{5000};               ///< Frequency of printing status reports.
  std::size_t writeBinaryRestartEvery{5000};  ///< Frequency of writing binary restart files.
  std::size_t rescaleWangLandauEvery{5000};   ///< Frequency of rescaling Wang-Landau factors.
  std::size_t optimizeMCMovesEvery{5000};     ///< Frequency of optimizing MC moves.
};

/**
 * \brief The number of Monte Carlo steps of a production cycle, fixed before production starts.
 *
 * A cycle of max(20, N) steps with the live molecule count N makes the length of a cycle depend on the
 * state it starts from. Properties are sampled once per cycle, so in an open ensemble the states of
 * high N are then sampled less often per step than the states of low N: the cycle-boundary chain
 * P^K(x) with a state-dependent K(x) does not leave the Boltzmann distribution invariant, and a
 * replica exchange between such chains does not either. With a state-independent K it does.
 *
 * The stages before production keep the live count (nothing is sampled there, and a box that starts
 * empty must be allowed to fill). The second half of the equilibration records the molecule count;
 * production then runs max(20, <N>) steps per component per cycle, <N> being that average (or the
 * count at the start of production when there was no equilibration).
 */
export struct ProductionCycleLength
{
  std::size_t moleculeCountSum{0};  ///< Sum of the molecule counts sampled during the equilibration.
  std::size_t numberOfSamples{0};   ///< Number of sampled molecule counts.
  std::size_t stepsPerCycle{0};     ///< The fixed number of steps of a production cycle.
  bool fixed{false};                ///< Whether 'stepsPerCycle' has been set.

  /// The number of steps of a cycle outside production: follows the live molecule count.
  static std::size_t adaptive(std::size_t numberOfMolecules, std::size_t numberOfComponents)
  {
    return std::max(numberOfMolecules, 20uz) * numberOfComponents;
  }

  void clear() { *this = ProductionCycleLength{}; }

  void sample(std::size_t numberOfMolecules)
  {
    moleculeCountSum += numberOfMolecules;
    ++numberOfSamples;
  }

  /// Fixes the cycle length from the recorded average, or from 'currentNumberOfMolecules' without samples.
  void fix(std::size_t currentNumberOfMolecules, std::size_t numberOfComponents)
  {
    const std::size_t mean = (numberOfSamples != 0uz)
                                 ? (moleculeCountSum + numberOfSamples / 2uz) / numberOfSamples
                                 : currentNumberOfMolecules;
    stepsPerCycle = adaptive(mean, numberOfComponents);
    fixed = true;
  }
};
