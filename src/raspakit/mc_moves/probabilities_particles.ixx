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
   * \param translationProbability Probability of performing a translation move.
   * \param randomTranslationProbability Probability of performing a random translation move.
   * \param rotationProbability Probability of performing a rotation move.
   * \param randomRotationProbability Probability of performing a random rotation move.
   * \param volumeChangeProbability Probability of performing a volume change move.
   * \param reinsertionCBMCProbability Probability of performing a reinsertion CBMC move.
   * \param identityChangeCBMCProbability Probability of performing an identity change CBMC move.
   * \param swapProbability Probability of performing a swap move.
   * \param swapCBMCProbability Probability of performing a swap CBMC move.
   * \param swapCFCMCProbability Probability of performing a swap CFCMC move.
   * \param swapCBCFCMCProbability Probability of performing a swap CBCFCMC move.
   * \param gibbsVolumeChangeProbability Probability of performing a Gibbs volume change move.
   * \param gibbsSwapCBMCProbability Probability of performing a Gibbs swap CBMC move.
   * \param gibbsSwapCFCMCProbability Probability of performing a Gibbs swap CFCMC move.
   * \param gibbsSwapCBCFCMCProbability Probability of performing a Gibbs swap CBCFCMC move.
   * \param widomProbability Probability of performing a Widom insertion move.
   * \param widomCFCMCProbability Probability of performing a Widom CFCMC move.
   * \param widomCBCFCMCProbability Probability of performing a Widom CBCFCMC move.
   * \param parallelTemperingProbability Probability of performing a parallel tempering move.
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
                               double parallelTemperingProbability = 0.0);

  double translationProbability;         ///< Probability of performing a translation move. ///<<
  double randomTranslationProbability;   ///< Probability of performing a random translation move. ///<<
  double rotationProbability;            ///< Probability of performing a rotation move. ///<<
  double randomRotationProbability;      ///< Probability of performing a random rotation move. ///<<
  double volumeChangeProbability;        ///< Probability of performing a volume change move. ///<<
  double reinsertionCBMCProbability;     ///< Probability of performing a reinsertion CBMC move. ///<<
  double identityChangeCBMCProbability;  ///< Probability of performing an identity change CBMC move. ///<<
  double swapProbability;                ///< Probability of performing a swap move. ///<<
  double swapCBMCProbability;            ///< Probability of performing a swap CBMC move. ///<<
  double swapCFCMCProbability;           ///< Probability of performing a swap CFCMC move. ///<<
  double swapCBCFCMCProbability;         ///< Probability of performing a swap CBCFCMC move. ///<<
  double gibbsVolumeChangeProbability;   ///< Probability of performing a Gibbs volume change move. ///<<
  double gibbsSwapCBMCProbability;       ///< Probability of performing a Gibbs swap CBMC move. ///<<
  double gibbsSwapCFCMCProbability;      ///< Probability of performing a Gibbs swap CFCMC move. ///<<
  double gibbsSwapCBCFCMCProbability;    ///< Probability of performing a Gibbs swap CBCFCMC move. ///<<
  double widomProbability;               ///< Probability of performing a Widom insertion move. ///<<
  double widomCFCMCProbability;          ///< Probability of performing a Widom CFCMC move. ///<<
  double widomCBCFCMCProbability;        ///< Probability of performing a Widom CBCFCMC move. ///<<
  double parallelTemperingProbability;   ///< Probability of performing a parallel tempering move. ///<<

  double accumulatedTranslationProbability{0.0};  ///< Accumulated probability up to translation move. ///<<
  double accumulatedRandomTranslationProbability{
      0.0};                                           ///< Accumulated probability up to random translation move. ///<<
  double accumulatedRotationProbability{0.0};         ///< Accumulated probability up to rotation move. ///<<
  double accumulatedRandomRotationProbability{0.0};   ///< Accumulated probability up to random rotation move. ///<<
  double accumulatedVolumeChangeProbability{0.0};     ///< Accumulated probability up to volume change move. ///<<
  double accumulatedReinsertionCBMCProbability{0.0};  ///< Accumulated probability up to reinsertion CBMC move. ///<<
  double accumulatedIdentityChangeCBMCProbability{
      0.0};                                       ///< Accumulated probability up to identity change CBMC move. ///<<
  double accumulatedSwapProbability{0.0};         ///< Accumulated probability up to swap move. ///<<
  double accumulatedSwapCBMCProbability{0.0};     ///< Accumulated probability up to swap CBMC move. ///<<
  double accumulatedSwapCFCMCProbability{0.0};    ///< Accumulated probability up to swap CFCMC move. ///<<
  double accumulatedSwapCBCFCMCProbability{0.0};  ///< Accumulated probability up to swap CBCFCMC move. ///<<
  double accumulatedGibbsVolumeChangeProbability{
      0.0};                                          ///< Accumulated probability up to Gibbs volume change move. ///<<
  double accumulatedGibbsSwapCBMCProbability{0.0};   ///< Accumulated probability up to Gibbs swap CBMC move. ///<<
  double accumulatedGibbsSwapCFCMCProbability{0.0};  ///< Accumulated probability up to Gibbs swap CFCMC move. ///<<
  double accumulatedGibbsSwapCBCFCMCProbability{0.0};  ///< Accumulated probability up to Gibbs swap CBCFCMC move. ///<<
  double accumulatedWidomProbability{0.0};             ///< Accumulated probability up to Widom insertion move. ///<<
  double accumulatedWidomCFCMCProbability{0.0};        ///< Accumulated probability up to Widom CFCMC move. ///<<
  double accumulatedWidomCBCFCMCProbability{0.0};      ///< Accumulated probability up to Widom CBCFCMC move. ///<<
  double accumulatedParallelTemperingProbability{
      0.0};  ///< Accumulated probability up to parallel tempering move. ///<<

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
