module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_probabilities_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;

/**
 * \brief Holds and manages the probabilities of various Monte Carlo moves in particle simulations.
 *
 * The MCMoveProbabilitiesParticles struct encapsulates the probabilities for different types of
 * Monte Carlo moves used in particle simulations, such as translations, rotations, volume changes,
 * swap moves, and more. It provides methods to normalize these probabilities so that they sum to 1,
 * and maintains accumulated probabilities for efficient selection of moves during simulations.
 */
export struct MCMoveProbabilitiesParticles
{
  uint64_t versionNumber{1};  ///< Version number for serialization purposes.

  bool operator==(MCMoveProbabilitiesParticles const &) const = default;

  /**
   * \brief Constructs an MCMoveProbabilitiesParticles object with specified probabilities.
   *
   * Initializes the move probabilities for various Monte Carlo moves. The probabilities
   * are normalized so that their sum equals 1.
   *
   * \param translationProbability Probability of performing a translation.
   * \param randomTranslationProbability Probability of performing a random translation.
   * \param rotationProbability Probability of performing a rotation.
   * \param randomRotationProbability Probability of performing a random rotation.
   * \param volumeChangeProbability Probability of performing a volume change.
   * \param reinsertionCBMCProbability Probability of performing a reinsertion CBMC.
   * \param identityChangeCBMCProbability Probability of performing an identity change CBMC.
   * \param swapProbability Probability of performing a swap.
   * \param swapCBMCProbability Probability of performing a swap CBMC.
   * \param swapCFCMCProbability Probability of performing a swap CFCMC.
   * \param swapCBCFCMCProbability Probability of performing a swap CBCFCMC.
   * \param gibbsVolumeChangeProbability Probability of performing a Gibbs volume change.
   * \param gibbsSwapCBMCProbability Probability of performing a Gibbs swap CBMC.
   * \param gibbsSwapCFCMCProbability Probability of performing a Gibbs swap CFCMC.
   * \param gibbsSwapCBCFCMCProbability Probability of performing a Gibbs swap CBCFCMC.
   * \param widomProbability Probability of performing a Widom insertion.
   * \param widomCFCMCProbability Probability of performing a Widom CFCMC.
   * \param widomCBCFCMCProbability Probability of performing a Widom CBCFCMC.
   * \param parallelTemperingProbability Probability of performing a parallel tempering.
   * \param hybridMCProbability Probability of performing a hybrid MC move.
   */
  MCMoveProbabilitiesParticles(double translationProbability = 0.0, double randomTranslationProbability = 0.0,
                               double rotationProbability = 0.0, double randomRotationProbability = 0.0,
                               double volumeChangeProbability = 0.0, double reinsertionCBMCProbability = 0.0,
                               double identityChangeCBMCProbability = 0.0, double swapProbability = 0.0,
                               double swapCBMCProbability = 0.0, double swapCFCMCProbability = 0.0,
                               double swapCBCFCMCProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                               double gibbsSwapCBMCProbability = 0.0, double gibbsSwapCFCMCProbability = 0.0,
                               double gibbsSwapCBCFCMCProbability = 0.0, double widomProbability = 0.0,
                               double widomCFCMCProbability = 0.0, double widomCBCFCMCProbability = 0.0,
                               double parallelTemperingProbability = 0.0, double hybridMCProbability = 0.0);

  double translationProbability;         ///< Probability of performing a translation.
  double randomTranslationProbability;   ///< Probability of performing a random translation.
  double rotationProbability;            ///< Probability of performing a rotation.
  double randomRotationProbability;      ///< Probability of performing a random rotation.
  double volumeChangeProbability;        ///< Probability of performing a volume change.
  double reinsertionCBMCProbability;     ///< Probability of performing a reinsertion CBMC.
  double identityChangeCBMCProbability;  ///< Probability of performing an identity change CBMC.
  double swapProbability;                ///< Probability of performing a swap.
  double swapCBMCProbability;            ///< Probability of performing a swap CBMC.
  double swapCFCMCProbability;           ///< Probability of performing a swap CFCMC.
  double swapCBCFCMCProbability;         ///< Probability of performing a swap CBCFCMC.
  double gibbsVolumeChangeProbability;   ///< Probability of performing a Gibbs volume change.
  double gibbsSwapCBMCProbability;       ///< Probability of performing a Gibbs swap CBMC.
  double gibbsSwapCFCMCProbability;      ///< Probability of performing a Gibbs swap CFCMC.
  double gibbsSwapCBCFCMCProbability;    ///< Probability of performing a Gibbs swap CBCFCMC.
  double widomProbability;               ///< Probability of performing a Widom insertion.
  double widomCFCMCProbability;          ///< Probability of performing a Widom CFCMC.
  double widomCBCFCMCProbability;        ///< Probability of performing a Widom CBCFCMC.
  double parallelTemperingProbability;   ///< Probability of performing a parallel tempering.
  double hybridMCProbability;            ///< Probability of performing a Hybrid MC.

  double accumulatedTranslationProbability{0.0};         ///< Normalized probability to perform translation.
  double accumulatedRandomTranslationProbability{0.0};   ///< Normalized probability to perform random translation.
  double accumulatedRotationProbability{0.0};            ///< Normalized probability to perform rotation.
  double accumulatedRandomRotationProbability{0.0};      ///< Normalized probability to perform random rotation.
  double accumulatedVolumeChangeProbability{0.0};        ///< Normalized probability to perform volume change.
  double accumulatedReinsertionCBMCProbability{0.0};     ///< Normalized probability to perform reinsertion CBMC.
  double accumulatedIdentityChangeCBMCProbability{0.0};  ///< Normalized probability to perform identity change CBMC.
  double accumulatedSwapProbability{0.0};                ///< Normalized probability to perform swap.
  double accumulatedSwapCBMCProbability{0.0};            ///< Normalized probability to perform swap CBMC.
  double accumulatedSwapCFCMCProbability{0.0};           ///< Normalized probability to perform swap CFCMC.
  double accumulatedSwapCBCFCMCProbability{0.0};         ///< Normalized probability to perform swap CBCFCMC.
  double accumulatedGibbsVolumeChangeProbability{0.0};   ///< Normalized probability to perform Gibbs volume change.
  double accumulatedGibbsSwapCBMCProbability{0.0};       ///< Normalized probability to perform Gibbs swap CBMC.
  double accumulatedGibbsSwapCFCMCProbability{0.0};      ///< Normalized probability to perform Gibbs swap CFCMC.
  double accumulatedGibbsSwapCBCFCMCProbability{0.0};    ///< Normalized probability to perform Gibbs swap CBCFCMC.
  double accumulatedWidomProbability{0.0};               ///< Normalized probability to perform Widom insertion.
  double accumulatedWidomCFCMCProbability{0.0};          ///< Normalized probability to perform Widom CFCMC.
  double accumulatedWidomCBCFCMCProbability{0.0};        ///< Normalized probability to perform Widom CBCFCMC.
  double accumulatedParallelTemperingProbability{0.0};   ///< Normalized probability to perform parallel tempering.
  double accumulatedHybridMCProbability{0.0};            ///< Normalized probability to perform hybrid MC.

  /**
   * \brief Normalizes the move probabilities so that they sum to 1.
   *
   * Adjusts the move probabilities so that their total sum equals 1. Also computes
   * the accumulated probabilities for efficient move selection during simulations.
   */
  void normalizeMoveProbabilities();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesParticles &p);
};
