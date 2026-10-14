module;

export module dctst_free_energy_settings;

import std;

import double3;

/**
 * \brief Settings of the dcTST free-energy driver ('DynamicallyCorrectedTSTFreeEnergy').
 *
 * Read from the 'dcTST' object of the system block. The driver computes the free-energy profile F(q) of
 * one tagged component along the barrier normal, with the precision at the barrier top that the
 * transition-state rate k^TST = <|qdot|>/2 P(q*) needs (an error delta in beta F(q*) is a relative error
 * delta in the diffusion coefficient).
 *
 * Stages:
 *  1. Initialization: unconstrained MC.
 *  2. Pilot: unconstrained MC; the histogram of the reaction coordinate of all molecules of the tagged
 *     component (converged in the cage interior, the unbiased reference state of the reweighting) and
 *     the joint cage-occupancy distribution P(n_A, n_B).
 *  3. Windows: one molecule of the tagged component is sampled under a bias w_k(q) = 1/2 kappa (q - q_k)^2
 *     + w_0(q) per window k (w_0 an optional flattening bias from the pilot), all other molecules free.
 *     The windows are distributed over the replicas (threads).
 *  4. Analysis: the pilot histogram and the window histograms are combined by binned multistate
 *     reweighting (WHAM, equivalent to MBAR for a bias that is a function of q only) on a fine grid; the
 *     barrier top is fitted by a parabola so that P(q*) is the density at the dividing plane and not a
 *     bin average; errors from blocks.
 *
 * The geometry keys are identical to those of the transmission driver so that the 'dcTST' block can be
 * copied between the two inputs.
 */
export struct DcTSTFreeEnergySettings
{
  // --- geometry (shared with the transmission driver) ----------------------------------------------
  std::string taggedComponentName{};
  std::size_t taggedComponent{0};
  double3 barrierPosition{0.5, 0.5, 0.5};
  double3 barrierNormal{1.0, 0.0, 0.0};
  double latticeDistance{0.0};
  std::optional<double3> cageCenter{};
  /// Occupancy-resolved mode ('CageOccupancies': [n_A, n_B]): the profile over cage A conditioned on n_A
  /// molecules in A (the tagged one included) and n_B in B.
  std::optional<std::pair<std::size_t, std::size_t>> cageOccupancies{};

  // --- stages ---------------------------------------------------------------------------------------
  /// Unconstrained MC cycles of the pilot stage ('NumberOfPilotCycles').
  std::size_t numberOfPilotCycles{0};
  /// Biased MC cycles per window, production ('NumberOfWindowCycles').
  std::size_t numberOfWindowCycles{0};
  /// Biased MC cycles per window before the production ('NumberOfWindowEquilibrationCycles').
  std::size_t numberOfWindowEquilibrationCycles{0};
  /// Sample the histograms every this many cycles ('SampleEvery').
  std::size_t sampleEvery{1};

  // --- windows --------------------------------------------------------------------------------------
  /// Number of windows ('NumberOfWindows'); 0 selects one window per 'WindowSpacing' Angstrom.
  std::size_t numberOfWindows{0};
  /// Spacing of the automatic window layout [A] ('WindowSpacing').
  double windowSpacing{0.5};
  /// Range [q_begin, q_end] of the window centers [A] ('WindowRange'); default the full period
  /// [-lambda/2, lambda/2) (q = 0 is the barrier, q = +-lambda/2 the cage center).
  std::optional<std::pair<double, double>> windowRange{};
  /// Explicit window centers [A] ('WindowCenters'); overrides the automatic layout.
  std::vector<double> windowCenters{};
  /// Stiffness beta kappa of the harmonic bias [1/A^2] ('WindowStiffness'); 0 selects 4 / spacing^2
  /// (standard deviation of the biased distribution equal to half the spacing).
  double windowStiffness{0.0};
  /// Add the flattening bias w_0(q) = -beta F_pilot(q) (smoothed, extrapolated into unsampled bins) to every
  /// window so that the sampling inside a window is nearly uniform ('FlatteningBias').
  bool flatteningBias{true};
  /// Fold the profile with q -> -q (a cage symmetric about the window plane): the windows on the negative
  /// side are dropped and both signs of q are histogrammed ('MirrorSymmetry').
  bool mirrorSymmetry{false};

  // --- grid and analysis ----------------------------------------------------------------------------
  /// Width of the fine bins of the reweighting [A] ('FineBinWidth').
  double fineBinWidth{0.01};
  /// Bins of the reported profile over one period ('NumberOfFreeEnergyBins').
  std::size_t numberOfFreeEnergyBins{100};
  /// Half-width of the parabolic fit of the barrier top [A] ('BarrierFitRange').
  double barrierFitRange{0.3};
  /// Convergence of the WHAM iteration ('WhamTolerance', 'WhamMaximumIterations').
  double whamTolerance{1.0e-9};
  std::size_t whamMaximumIterations{200000};
  /// Write the per-window fine histograms (for an external MBAR cross-check) ('WriteWindowHistograms').
  bool writeWindowHistograms{true};

  // --- tagged-molecule moves ------------------------------------------------------------------------
  /// Initial maximum displacement of the biased translation [A], adapted during the run
  /// ('MaximumTranslationTagged').
  double maximumTranslationTagged{0.5};
  /// Probability that a move of the tagged molecule is a CBMC reinsertion within the window support
  /// ('ReinsertionTaggedProbability').
  double reinsertionTaggedProbability{0.3};
  /// Fraction of the Monte Carlo steps of the window stage that move the tagged molecule (the other steps
  /// select a random molecule as usual, which may again be the tagged one). With N molecules a random
  /// selection moves the tagged molecule only once per N steps; a fixed, configuration-independent selection
  /// probability keeps detailed balance and gives the window histogram many more decorrelated samples
  /// ('TaggedMoveFraction').
  double taggedMoveFraction{0.5};

  /// Independent replicas (threads) ('NumberOfReplicas', defaults to the number of threads). The pilot
  /// runs on every replica; the windows are distributed over the replicas.
  std::size_t numberOfReplicas{1};
};
