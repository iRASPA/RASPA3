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
