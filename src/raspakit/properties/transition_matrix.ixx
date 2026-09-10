module;

export module transition_matrix;

import std;

import archive;
import double3;

/**
 * \brief Represents the transition matrix for TMMC simulations.
 *
 * The TransitionMatrix struct encapsulates the data and methods required for performing
 * Transition Matrix Monte Carlo (TMMC) simulations. It maintains the collection matrix,
 * bias factors, natural logarithm of probability distributions, and histograms to
 * compute macrostates and adjust biases during the simulation. It provides methods
 * to initialize data structures, update matrices and histograms, calculate bias factors,
 * adjust biases, clear statistics, and write statistics to files.
 */
export struct TransitionMatrix
{
  std::uint64_t versionNumber{3};  ///< Version number for serialization compatibility.

  bool operator==(TransitionMatrix const&) const = default;

  std::vector<double3> cmatrix;        ///< Collection matrix: x=deletion, y=no change, z=insertion.
  std::vector<double> bias;            ///< Bias factors for each macrostate.
  std::vector<double> lnpi;            ///< Natural logarithm of the probability distribution.
  std::vector<double> forward_lnpi;    ///< Forward ln(pi) for debugging purposes.
  std::vector<double> reverse_lnpi;    ///< Reverse ln(pi) for debugging purposes.
  std::vector<std::size_t> histogram;  ///< Histogram of macrostate visits.

  std::size_t numberOfSteps = {0};        ///< Number of steps performed in the TMMC simulation.
  std::size_t minMacrostate = {0};        ///< Minimum molecule count (window bound).
  std::size_t maxMacrostate = {100};      ///< Maximum molecule count (window bound).
  /// CFCMC lambda bins along the flattened (N, λ) chain. 1 keeps the historical 1D N-only matrix.
  std::size_t numberOfLambdaBins = {1};
  /// Lambda bin of the selected component *before* the current move (for updateMatrix).
  std::size_t currentLambdaBin = {0};

  [[nodiscard]] bool lambdaChain() const { return numberOfLambdaBins > 1uz; }
  [[nodiscard]] std::size_t lambdaBinCount() const { return std::max(1uz, numberOfLambdaBins); }
  [[nodiscard]] std::size_t lastLambdaBin() const { return lambdaBinCount() - 1uz; }
  [[nodiscard]] std::size_t numberOfMoleculeStates() const { return maxMacrostate - minMacrostate + 1uz; }
  [[nodiscard]] std::size_t numberOfChainStates() const { return numberOfMoleculeStates() * lambdaBinCount(); }

  /// Flatten (N, k) onto the 1D chain: index = (N - minN) * nλ + k. Out of range returns numberOfChainStates().
  [[nodiscard]] std::size_t chainIndex(std::size_t N, std::size_t lambdaBin) const
  {
    if (N < minMacrostate || N > maxMacrostate) return numberOfChainStates();
    const std::size_t nLambda = lambdaBinCount();
    const std::size_t k = lambdaChain() ? lambdaBin : 0uz;
    if (k >= nLambda) return numberOfChainStates();
    return (N - minMacrostate) * nLambda + k;
  }
  std::size_t updateTMEvery = {1000000};  ///< Number of steps between bias updates.
  std::size_t numberOfUpdates = {0};      ///< Number of times the bias has been updated.

  bool doTMMC = {false};   ///< Flag indicating whether to perform TMMC simulation.
  bool useBias = {false};  ///< Flag indicating whether to use bias for changing macrostates.

  /// File name for the statistics output; parallel drivers set a unique name per walker.
  std::string statisticsFileName{"tmmc/tmmc_statistics.txt"};
  bool useTMBias = {true};                   ///< Flag indicating whether to use Transition Matrix bias.
  bool rejectOutOfBound = {true};            ///< Flag indicating whether to reject moves outside macrostate bounds.
  bool rezeroAfterInitialization = {false};  ///< Flag indicating whether to reset statistics after initialization.

  /// Wang-Landau flattening of the macrostate bias.
  ///
  /// The transition-matrix bias can only be sharpened where the collection matrix already holds counts, and
  /// `adjustBias` flattens everything outside the visited range to the value at its edge. A walker that
  /// reaches equilibrium loading therefore sees no gradient back down and stops exploring: at 77 K in a
  /// zeolite ln Pi spans several hundred over the window, which is not a barrier a walk crosses unaided.
  /// Wang-Landau penalises whichever state is occupied, so the bias grows in the states already seen and
  /// pushes the walker into the ones it has not, which is what makes the walk over N ergodic. The estimate
  /// itself is unaffected: the collection matrix records the unbiased acceptance probabilities, so ln Pi
  /// stays an unbiased estimator no matter what bias drove the sampling.
  std::vector<double> wangLandauHistogram;  ///< Visits per macrostate since the last flatness test.
  double wangLandauFactor = {1.0};          ///< The modification factor ln f, halved on each flat histogram.
  double wangLandauFactorFloor = {1.0e-6};  ///< Below this ln f the bias is effectively frozen.
  double wangLandauFlatness = {0.8};        ///< A histogram is flat when its least entry reaches this * mean.
  std::size_t wangLandauCheckEvery = {10000};  ///< Visits between flatness tests.
  std::size_t wangLandauVisits = {0};          ///< Visits since the last flatness test.
  bool useWangLandau = {false};                ///< Flag indicating whether Wang-Landau owns the bias.

  /**
   * \brief Initializes the transition matrix and related data structures.
   *
   * Sets up the collection matrix, bias vectors, and histograms based on the specified
   * macrostate range. Initializes bias factors and probability distributions to default values.
   * Should be called before starting the TMMC simulation.
   */
  void initialize();

  /**
   * \brief Updates the collection matrix with acceptance probabilities.
   *
   * Records the probabilities of transitions between macrostates in the collection matrix.
   * The acceptance probabilities vector Pacc contains probabilities for deletion (x),
   * no change (y), and insertion (z) moves. A zero vector is a no-op (used by 1D-N TMMC
   * to skip CFCMC lambda hops). With numberOfLambdaBins > 1 the matrix is the flattened
   * (N, λ) chain and lambda hops are real ±1 transitions; the 2-argument overload records
   * at currentLambdaBin.
   *
   * \param Pacc A double3 vector containing acceptance probabilities for deletion (x),
   *             no change (y), and insertion (z) moves.
   * \param oldN The molecule count before the move.
   */
  void updateMatrix(double3 Pacc, std::size_t oldN);
  void updateMatrix(double3 Pacc, std::size_t oldN, std::size_t oldLambdaBin);

  /**
   * \brief Updates the histogram of macrostate visits.
   *
   * Increments the histogram count for the given macrostate N, tracking how often each
   * macrostate is visited during the simulation.
   *
   * \param N The current molecule count.
   * \param lambdaBin The current lambda bin (ignored when numberOfLambdaBins is 1).
   */
  void updateHistogram(std::size_t N, std::size_t lambdaBin = 0);

  /**
   * \brief Penalises the occupied macrostate and runs the Wang-Landau schedule.
   *
   * Subtracts the current modification factor from the bias of macrostate N and counts the visit. Once
   * `wangLandauCheckEvery` visits have accumulated the histogram is tested for flatness, and a flat one
   * halves the modification factor and starts a new stage. Does nothing unless `useWangLandau` is set.
   *
   * \param N The current molecule count.
   * \param lambdaBin The current lambda bin (ignored when numberOfLambdaBins is 1).
   */
  void visitWangLandau(std::size_t N, std::size_t lambdaBin = 0);

  /**
   * \brief Calculates the bias factor between two macrostates.
   *
   * Computes the bias factor used to adjust the acceptance probability of moves
   * between macrostates. The bias factor is calculated as the exponential of the
   * difference in bias between the new and old macrostate.
   *
   * \param newN The molecule count after the move.
   * \param oldN The molecule count before the move.
   * \return The bias factor for the transition from oldN to newN.
   *
   * The two-argument form infers lambda bins: same N uses currentLambdaBin; insertion
   * is (N, last) → (N+1, 0); deletion is (N, 0) → (N−1, last). Lambda hops must call
   * the four-argument form.
   */
  double biasFactor(std::size_t newN, std::size_t oldN);
  double biasFactor(std::size_t newN, std::size_t oldN, std::size_t newLambdaBin, std::size_t oldLambdaBin);

  /**
   * \brief Adjusts the bias factors based on collected statistics.
   *
   * Updates the bias vector and the natural logarithm of the probability distribution
   * (lnpi) using the current state of the collection matrix and histogram. This method
   * recalculates biases to improve sampling efficiency.
   */
  void adjustBias();

  /**
   * \brief Hands the bias from Wang-Landau to the transition matrix.
   *
   * Clears `useWangLandau` and immediately rebuilds the bias as -ln Pi from the collection
   * matrix (the exact flattening bias), without waiting for the next `updateTMEvery`
   * boundary. Call once exploration is done, e.g. at the start of production.
   */
  void switchToTMBias();

  /**
   * \brief Recomputes ln Pi from the collection matrix and re-derives the bias.
   *
   * The shared core of adjustBias() and switchToTMBias(): runs the detailed-balance
   * recursion over the visited range, flattens ln Pi outside it, normalizes, and writes
   * bias = -ln Pi unless Wang-Landau owns the bias. Does nothing before the first visit.
   */
  void recomputeLnPiAndBias();

  /**
   * \brief Clears the collection matrix and resets counters.
   *
   * Resets the collection matrix, histogram, bias factors, and related counters to
   * their initial values. Typically used after initialization cycles to reset statistics.
   * Only has an effect when rezeroAfterInitialization is set.
   */
  void clearCMatrix();

  /**
   * \brief Resets the collection matrix and visit histogram but keeps the bias.
   *
   * Unlike clearCMatrix() the bias (including the Wang-Landau state) survives, so an
   * ongoing biased walk continues undisturbed while the statistics start from zero.
   * Used at production start to drop the pre-production samples, which were taken on
   * configurations that had not yet relaxed into the deep adsorption sites.
   */
  void clearStatisticsKeepBias();

  /**
   * \brief Writes the transition matrix statistics to a file.
   *
   * Outputs the current state of the collection matrix, bias factors, probability
   * distributions, and histogram data to a text file for analysis and debugging.
   */
  void writeStatistics();

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const TransitionMatrix& m);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, TransitionMatrix& m);
};
