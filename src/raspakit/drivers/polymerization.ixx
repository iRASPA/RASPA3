module;

export module polymerization;

import std;

import archive;
import input_reader;
import molecular_dynamics_spatial_decomposition;
import polymerization_settings;

/**
 * \brief Stochastic free-radical chain-growth polymerization with the spatial-decomposition MD engine.
 *
 * Selected with 'SimulationType' : 'Polymerization'. The driver is the MolecularDynamicsSpatialDecomposition
 * driver with a reaction step: the Monte Carlo (pre-)initialization and the MD equilibration run on the
 * unreacted monomers; the production stage is the curing stage, in which every 'ScanEvery' MD steps the
 * reactive sites ('PolymerizationSites' of the molecule files) are scanned and reactions are accepted with
 * their probabilities (see PolymerizationSettings):
 *
 *   initiation   M -> R*                  (monomer site becomes a radical)
 *   propagation  R* + M -> R-R*           (bond radical tail - monomer head; the monomer is the new radical)
 *   termination  R* + R* -> R-R           (bond tail - tail, combination)
 *
 * A reaction retypes the head and the tail of the reacted site to the product pseudo-atoms and adds the
 * bond with its bends, torsions and non-bonded exclusions as an overlay on the component topology (see
 * PolymerLinkTable); the new terms are ramped in over 'RelaxationSteps' steps. The molecules of the system
 * never change: the polymer is the graph of bonds between them, written at the end to
 * 'output/polymer_topology.s<i>.json' (and every reaction to 'output/polymerization_bonds.s<i>.txt').
 *
 * Restrictions: the resident (fully on-device) integrator and the device bonded kernels are switched off (the
 * overlay is evaluated on the host; the device pair and mesh kernels stay in use); 'ReactiveSites' /
 * 'CrossLinkBonds' (the equilibrium cross-link moves) are a separate mechanism and are not combined with it.
 */
export struct Polymerization
{
  Polymerization(InputReader& reader);

  std::uint64_t versionNumber{1};

  MolecularDynamicsSpatialDecomposition md;
  PolymerizationSettings settings{};

  /// Per system: whether the reactions have stopped (target conversion reached, or nothing left to react).
  std::vector<std::uint8_t> reactionsFinished{};
  /// Per system: the number of reaction scans performed.
  std::vector<std::size_t> numberOfScans{};
  /// Per system: the production cycle at which the reactions were activated.
  std::vector<std::size_t> activationCycle{};

  std::vector<std::ofstream> bondLogs{};  ///< Not part of the restart (reopened in append mode).

  void run();
  /// Reopens the output files of a restarted run (the MD driver's streams and the bond logs, in append mode).
  void createOutputFiles();
  void openBondLogs();

  /// The hook of the MD driver: called after every MD step of a system.
  void afterStep(std::size_t systemId);
  /// The status block of the MD driver's reports.
  std::string status(std::size_t systemId);

  /// Activates the table of a system (site layout, initial radicals) at the start of the reactive stage.
  void activate(std::size_t systemId);
  /// One reaction scan of a system; returns the number of accepted reactions.
  std::size_t scan(std::size_t systemId);
  /// After accepted reactions: tells the engine, refreshes the forces and the energies of the system.
  void topologyChanged(std::size_t systemId);
  void writeSettings();
  void writeTopology();

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Polymerization& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Polymerization& p);
};
