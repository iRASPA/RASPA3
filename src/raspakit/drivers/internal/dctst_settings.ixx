module;

export module dctst_settings;

import std;

import double3;

/**
 * \brief Settings of the dynamically corrected transition state theory (dcTST) driver.
 *
 * Read from the 'dcTST' object of the system block of the input file. The driver computes, for one
 * tagged component, the self and collective transmission functions kappa_S(t) and kappa_C(t) along the
 * reaction coordinate q (the projection of the component's reaction bead onto the barrier normal).
 * Combined with the probability density P(q*) at the dividing plane from the free-energy driver
 * ('DynamicallyCorrectedTSTFreeEnergy') it reports the transition-state-theory hopping rate k^TST, the
 * self-diffusion coefficient D_S = kappa_S k^TST lambda^2 and the collective (transport) diffusion
 * coefficient D_C = kappa_C k^TST lambda^2, where lambda is the cage-to-cage distance.
 *
 * The dividing surface is the plane through 'barrierPosition' (fractional coordinates) with unit
 * normal 'barrierNormal' (Cartesian). The periodic images of the plane at multiples of
 * 'latticeDistance' along the normal bound the cages; cage index c = floor(q / lambda).
 */
export struct DcTSTSettings
{
  /**
   * \brief Stochastic bath acting on the molecules during the trajectories.
   *
   * In a rigid framework at low loading a molecule that crosses the window keeps its excess energy:
   * nothing thermalizes it in the next cage, so it runs on to a second window or bounces back through
   * the first one, and kappa(t) has no plateau (Kramers' energy-diffusion limit). The bath stands in for
   * the lattice vibrations of a real framework. Its rate should lie between the time scale of the
   * recrossing dynamics in the window (~0.1-0.3 ps) and the cage residence time (1/k), where kappa is
   * insensitive to it.
   */
  enum class TrajectoryThermostat : std::size_t
  {
    None = 0,      ///< NVE trajectories (diagnostic mode at low loading).
    Andersen = 1,  ///< Per molecule: with probability rate*dt per step the full molecular velocity is redrawn.
    Langevin = 2   ///< Ornstein-Uhlenbeck friction + noise on the molecular velocities, friction = rate.
  };

  /// Name of the component whose hopping is studied ('TaggedComponent').
  std::string taggedComponentName{};
  /// Index of the tagged component, resolved by the input reader.
  std::size_t taggedComponent{0};

  /// Position of the barrier plane in fractional coordinates of the simulation box ('BarrierPosition').
  double3 barrierPosition{0.5, 0.5, 0.5};
  /// Normal of the barrier plane in Cartesian coordinates, normalized by the driver ('BarrierNormal').
  double3 barrierNormal{1.0, 0.0, 0.0};
  /// Distance between equivalent barrier planes along the normal (the cage-to-cage distance, in
  /// Angstrom) ('LatticeDistance').
  double latticeDistance{0.0};
  /// Centre of a cage in fractional coordinates of the framework *unit cell* ('CageCenter'); it fixes the
  /// cage lattice of the occupancy-resolved and event-terminated modes (one cage per unit cell). Must lie
  /// half a lattice distance from a barrier plane along the normal. When absent the driver uses the barrier
  /// point shifted by lambda/2 along the normal, which puts the in-plane centre at the in-plane position of
  /// 'BarrierPosition' -- correct only when that supercell fraction happens to fall on a cage centre (e.g. an
  /// odd number of unit cells in the plane for LTA).
  std::optional<double3> cageCenter{};

  /// Monte Carlo cycles of the transmission stage (tagged molecule constrained to the barrier plane,
  /// trajectories launched periodically) ('NumberOfTransmissionCycles').
  std::size_t numberOfTransmissionCycles{0};
  /// Monte Carlo cycles of constrained equilibration before the first trajectory is launched
  /// ('NumberOfConstrainedEquilibrationCycles').
  std::size_t numberOfConstrainedEquilibrationCycles{0};

  /// The free-energy profile is computed by the separate driver 'DynamicallyCorrectedTSTFreeEnergy'. Its
  /// JSON output ('FreeEnergyInput': path of the output_*.dctst_free_energy.json file) supplies the
  /// probability density P(q*) of the reaction bead at the dividing plane, from which the rates and the
  /// diffusion coefficients follow; alternatively P(q*) and its error are given directly
  /// ('BarrierProbability', 'BarrierProbabilityError' [1/Angstrom]). Without either, only the transmission
  /// coefficients are reported.
  std::optional<std::string> freeEnergyInput{};
  std::optional<double> barrierProbability{};
  double barrierProbabilityError{0.0};

  /// Length of every trajectory in picoseconds ('TrajectoryTime').
  double trajectoryTime{10.0};
  /// Launch a trajectory pair (forward and time-reversed) every this many constrained MC cycles
  /// ('LaunchTrajectoryEvery').
  std::size_t launchTrajectoryEvery{10};
  /// Number of independent Maxwell-Boltzmann velocity draws (trajectory pairs) per launched
  /// configuration ('TrajectoriesPerConfiguration').
  std::size_t trajectoriesPerConfiguration{1};
  /// Store the transmission functions every this many integration steps ('SampleTrajectoryEvery').
  std::size_t sampleTrajectoryEvery{10};

  /// Largest cage shell (in cage-index distance from the pair of cages A and B) for which the
  /// collective transmission function is accumulated separately; the sum over all molecules is
  /// always reported too ('MaximumCollectiveShell').
  std::size_t maximumCollectiveShell{2};

  /// Time window [begin, end] (ps) over which the plateau values of the transmission functions are
  /// averaged; defaults to the second half of the trajectory ('PlateauWindow').
  std::optional<std::pair<double, double>> plateauWindow{};

  /// Number of independent replicas (threads); every replica runs the full protocol with its own
  /// random-number stream and the results are averaged ('NumberOfReplicas', defaults to the number
  /// of threads).
  std::size_t numberOfReplicas{1};

  /// Maximum in-plane displacement of the tagged molecule (Angstrom), adapted during the run
  /// ('MaximumTranslationInPlane').
  double maximumTranslationInPlane{1.0};
  /// Probability that a Monte Carlo move of the tagged molecule is an in-plane reinsertion (the
  /// remainder are in-plane translations) ('ReinsertionInPlaneProbability').
  double reinsertionInPlaneProbability{0.5};

  /// Bath acting on the molecules during the trajectories ('TrajectoryThermostat': "None", "Andersen",
  /// "Langevin"). Default Andersen.
  TrajectoryThermostat trajectoryThermostat{TrajectoryThermostat::Andersen};
  /// Rate of the bath in 1/ps: the collision frequency per molecule (Andersen) or the friction
  /// coefficient (Langevin) ('BathRate').
  double bathRate{1.0};

  /// Distance |q| from the barrier (Angstrom) at which the tagged molecule is considered committed to a
  /// cage for the 'committed' estimator of kappa_S (its contribution is frozen at +-1/2 from then on,
  /// the RASPA2 convention of terminating the trajectory); 0 selects lambda/2, the cage center
  /// ('CommitDistance').
  double commitDistance{0.0};

  /// Write the trajectory diagnostics (energy drift per trajectory) ('WriteTrajectoryDiagnostics').
  bool writeTrajectoryDiagnostics{false};

  /// Write one line per completed trajectory pair to 'dctst/pairs_{T}_{P}.r{k}.txt' (per replica, flushed after
  /// every pair), from which every estimator can be rebuilt as a ratio of sums and runs can be pooled exactly
  /// ('WritePairRecords', default true).
  bool writePairRecords{true};
  /// Interval in seconds at which every replica writes a JSON snapshot of its accumulators to
  /// 'dctst/accumulator_{T}_{P}.r{k}.json' (atomically replaced; the transmission functions of a run that was
  /// killed can be recovered from it) ('SnapshotInterval', default 600; 0 disables).
  double snapshotInterval{600.0};

  /**
   * \brief How a trajectory ends ('TrajectoryTermination').
   *
   * With fixed-time trajectories of a few picoseconds the collective transmission function equals the
   * committed self one by construction: the correlated response to a hop (refilling of the hole in A,
   * departure of the excess from B, back-hop of the tagged molecule) takes place on the time scale of
   * the hops themselves, and the far-field hops of the other molecules are pure noise. The event-terminated
   * estimator follows only the local collective coordinate (the axial window crossings in the column of
   * cages A', A, B, B') and ends every trajectory at the 'terminationEvents'-th local hop event after the
   * tagged molecule has committed; 'trajectoryTime' is then the cap.
   */
  enum class TrajectoryTermination : std::size_t
  {
    FixedTime = 0,  ///< Every trajectory runs for 'trajectoryTime'.
    Event = 1       ///< Terminate at the n-th local axial hop event after commitment ('trajectoryTime' is the cap).
  };
  TrajectoryTermination trajectoryTermination{TrajectoryTermination::FixedTime};
  /// Number of local axial window crossings, counted after the commitment of the tagged molecule, at
  /// which an event-terminated trajectory ends ('TerminationEvents').
  std::size_t terminationEvents{1};

  /**
   * \brief Occupancy-resolved mode ('CageOccupancies': [n_A, n_B]).
   *
   * Fixes the occupancies of the two cages that share the window of the tagged molecule: n_A molecules
   * in cage A (the tagged molecule is one of them) and n_B molecules in cage B, both counted before the
   * hop A -> B and over all adsorbate components. In the transmission stage moves of other molecules
   * that would change either occupancy are rejected. The matching free-energy run (the same
   * 'CageOccupancies' in 'DynamicallyCorrectedTSTFreeEnergy') supplies the conditioned P(q*|n_A, n_B)
   * and the joint distribution P(n_A, n_B), so that the pair of runs yields the occupancy-resolved
   * hopping rate k(n_A -> n_B) = kappa(n_A, n_B) k^TST(n_A, n_B). Cages are the unit cells of the
   * framework lattice centered on the cage centers (one cage per unit cell).
   */
  std::optional<std::pair<std::size_t, std::size_t>> cageOccupancies{};
};
