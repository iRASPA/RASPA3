module;

export module interactions_neighbourhood;

import std;

import double3;
import atom;
import running_energy;
import simulationbox;
import forcefield;
import mc_cell_list;

export namespace Interactions
{
/**
 * \brief The background atoms a batch of Monte Carlo trial positions interacts with, gathered once into a
 * structure-of-arrays copy with the periodic images resolved.
 *
 * The trial positions of one Monte Carlo decision share their surroundings: the k trial directions of a
 * CBMC growth step lie within a bond length of the same anchor, the trial orientations of a rigid
 * reinsertion within the molecule radius of its first bead, the new and the old positions of a translation
 * within the maximum displacement of each other. Instead of visiting the cells around every trial set
 * separately and fetching the 64-byte Atom records pair by pair, the atoms within 'cutOff' of every point of
 * the batch's bounding sphere (centre, batchRadius) are gathered once: their positions as periodic images
 * placed relative to the centre, and their type, charge, scaling factors, group id and molecule id in parallel
 * arrays. The trial energies are then dense loops over these arrays ('computeNeighbourhoodEnergy'), without
 * a minimum-image operation per pair and without touching the Atom array -- and, for a plain Lennard-Jones /
 * Ewald force field between fully coupled atoms, in SIMD lanes (the same portable 'SIMD::Lanes' as the cluster
 * kernel of the spatial-decomposition MD engine).
 *
 * Images: when the gather radius (batchRadius + cutOff) is below half the smallest perpendicular width of
 * the box, the image nearest to the centre is the minimum image of every pair within the cut-off of any trial
 * position in the batch ('nearestImageUnique'). In a smaller box every periodic image within the gather
 * radius is gathered (several records per atom), which is exact as long as the cut-off itself stays within
 * half the box (at most one image of an atom within the cut-off of a trial position) and the gather radius
 * stays below the box width (the images enumerated are the nearest one and its neighbours across the box);
 * either way the kernels need no minimum image ('imagesResolved'). A cut-off beyond half the box, or a batch
 * so wide that the gather radius reaches the box width, leaves the minimum image to the kernel, pair by pair
 * (one record per atom, its nearest image to the centre).
 *
 * The gather goes through the cell list of the system when it is enabled and valid (the cells within the
 * gather radius), and over all atoms otherwise. For trial positions spread over the whole box (the first
 * beads of a CBMC growth in a box too small for the cell list) 'gatherAll' copies every atom unshifted
 * instead; such a neighbourhood covers any trial position and the kernels apply the minimum image per pair.
 * A neighbourhood is a reusable scratch buffer: the vectors keep their capacity between gathers. The arrays
 * are padded to a multiple of 'padding' entries with records far from the centre carrying the reserved
 * 'paddingMoleculeId', which the kernels mask out, so that a lane kernel can run over whole blocks.
 *
 * Background: the same gathers and kernels serve the framework atoms (a rigid framework is a static
 * background; 'System::frameworkCellList' is a cell list over it built once). A framework neighbourhood
 * ('Background::Framework') records every atom under the reserved 'frameworkMoleculeId' (the framework atoms
 * carry molecule id 0, which a molecule may carry too, and framework-molecule pairs are never excluded), and
 * the kernels book its energies as framework-molecule terms.
 */
struct Neighbourhood
{
  static constexpr std::size_t padding = 8;
  static constexpr std::uint32_t paddingMoleculeId = std::numeric_limits<std::uint32_t>::max();
  static constexpr std::uint32_t frameworkMoleculeId = std::numeric_limits<std::uint32_t>::max() - 1;

  enum class Background : std::uint8_t
  {
    Molecules = 0,
    Framework = 1,
  };
  Background background{Background::Molecules};

  double3 centre{};
  double batchRadius{0.0};  ///< Every trial position of the batch lies within this distance of the centre.
  double cutOff{0.0};       ///< The pair cut-off the gather covers (the largest cut-off a query may use).
  bool nearestImageUnique{true};
  bool imagesResolved{true};  ///< The records are the images the pairs use: no minimum image in the kernel.
  bool fullyCoupled{true};    ///< Every gathered atom has scaling 1 and no dU/dlambda group (no fractional atoms).
  bool allAtoms{false};       ///< Every atom is present ('gatherAll'): any trial position is covered.
  std::size_t numberOfCandidates{0};  ///< Records visited by the last gather (statistics).

  std::vector<double> x{}, y{}, z{};  ///< Image positions (see the class comment), padded.
  std::vector<double> charge{}, scalingVDW{}, scalingCoulomb{};
  std::vector<std::uint32_t> type{};
  std::vector<std::uint8_t> groupId{};
  std::vector<std::uint32_t> moleculeId{};
  std::vector<double> moleculeIdValue{};  ///< 'moleculeId' as doubles (exact), for the lane compare.

  /// The Lennard-Jones parameters of one trial-atom type against every record (4 epsilon, sigma^6 and the
  /// shift; zeros and sigma^6 = 1 for a pair without van der Waals interaction), laid out along the records so
  /// that the lane kernel streams through them. Derived from the live force field on the first query of a
  /// type after a gather and kept for the following queries of the same gather (the trial sets of a batch).
  struct TypeRow
  {
    std::uint32_t type{};
    std::vector<double> epsilon4{}, sigma6{}, shift{};
  };
  [[nodiscard]] const TypeRow& lennardJonesRow(const ForceField& forceField, std::uint32_t type) const;

  /// Number of gathered records (without the padding).
  [[nodiscard]] std::size_t size() const { return count; }
  /// Length of the arrays (a multiple of 'padding').
  [[nodiscard]] std::size_t paddedSize() const { return x.size(); }

  /**
   * \brief Gathers the atoms of 'moleculeAtoms' within batchRadius + cutOff of 'centre' (all periodic images
   * when the radius exceeds half the box, see the class comment), through 'cellList' when it is enabled and
   * valid (it must then be built over 'moleculeAtoms'), over all atoms otherwise. The atoms of 'skipMolecule'
   * are left out.
   */
  void gather(const SimulationBox& simulationBox, const MCCellList* cellList, std::span<const Atom> moleculeAtoms,
              const double3& centre, double batchRadius, double cutOff,
              std::optional<std::uint32_t> skipMolecule = std::nullopt);

  /**
   * \brief Copies every atom of 'moleculeAtoms' (but those of 'skipMolecule') unshifted: a neighbourhood that
   * covers any trial position ('allAtoms'), with the minimum image left to the kernels ('imagesResolved'
   * false). For queries spread over the whole box when no cell list is usable.
   */
  void gatherAll(std::span<const Atom> moleculeAtoms, double cutOff,
                 std::optional<std::uint32_t> skipMolecule = std::nullopt);

  /// 'gather' over the framework atoms (see the class comment on 'Background'); 'cellList', when usable, must
  /// have been built over 'frameworkAtoms'.
  void gatherFramework(const SimulationBox& simulationBox, const MCCellList* cellList,
                       std::span<const Atom> frameworkAtoms, const double3& centre, double batchRadius, double cutOff);

  /// 'gatherAll' over the framework atoms.
  void gatherAllFramework(std::span<const Atom> frameworkAtoms, double cutOff);

  /// Whether every atom of 'atoms' lies within the batch radius of the centre and the cut-offs do not
  /// exceed the gathered one: the neighbourhood then holds every background atom the query can interact with.
  [[nodiscard]] bool covers(std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb, bool useCharge) const;

  /// The bounding sphere (centroid and the largest distance to it) of the positions of several atom sets.
  [[nodiscard]] static std::pair<double3, double> boundingSphere(std::span<const std::span<const Atom>> sets);
  [[nodiscard]] static std::pair<double3, double> boundingSphere(std::span<const Atom> atoms)
  {
    return boundingSphere(std::span<const std::span<const Atom>>(&atoms, 1));
  }

 private:
  std::size_t count{0};
  mutable std::vector<TypeRow> rows{};  ///< 'lennardJonesRow' cache, cleared by the gathers.
  void append(const Atom& atom, const double3& image, std::uint32_t recordedMoleculeId);
  void clear();
  void pad();
  void gatherImpl(const SimulationBox& simulationBox, const MCCellList* cellList, std::span<const Atom> atoms,
                  const double3& centre, double batchRadius, double cutOff, std::optional<std::uint32_t> skipMolecule,
                  Background kind);
  void gatherAllImpl(std::span<const Atom> atoms, double cutOff, std::optional<std::uint32_t> skipMolecule,
                     Background kind);
  /// Marks the cached Lennard-Jones rows stale (their buffers are kept).
  void invalidateRows();
  /// The block gather out of a usable cell list (see 'gatherImpl'); 'skipped' is the molecule id left out.
  void gatherFromCells(const SimulationBox& simulationBox, const MCCellList& cellList, double radius, bool framework,
                       std::uint32_t skipped);
};

/**
 * \brief Whether 'computeNeighbourhoodEnergy' evaluates 'atoms' against 'neighbourhood' with the lane kernel:
 * fully coupled atoms on both sides, plain (possibly shifted) 12-6 Lennard-Jones or no van der
 * Waals interaction between the trial atom types and every type, and Ewald electrostatics read from the force
 * field's real-space table (or no charges). Exposed for the tests and the status report.
 */
[[nodiscard]] bool neighbourhoodFastKernelApplies(const ForceField& forceField, const Neighbourhood& neighbourhood,
                                                  std::span<const Atom> atoms) noexcept;

/**
 * \brief Van der Waals and real-space Coulomb energy of 'atoms' against the atoms of the neighbourhood at
 * the given cut-offs, scaled by 'sign' (+1 for a trial configuration, -1 for the configuration it replaces).
 *
 * Pairs of the same molecule id do not interact (a framework neighbourhood records a reserved id, so every
 * framework-molecule pair counts). The sums are booked as molecule-molecule or framework-molecule terms after
 * the neighbourhood's 'background'. With 'checkOverlap' a VDW pair energy above the force
 * field's overlap criterion abandons the evaluation (std::nullopt). The pair arithmetic is that of
 * 'Interactions::evaluatePair' with the background atom as atom A (dr = posA - posB), so the sums equal those
 * of the brute-force routines to rounding (the lane kernel reorders the sums and reads the Lennard-Jones
 * parameters as 4 epsilon and sigma^6, the generic kernel matches the summation order of the brute force).
 */
[[nodiscard]] std::optional<RunningEnergy> computeNeighbourhoodEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const Neighbourhood& neighbourhood,
    std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap) noexcept;

/**
 * \brief Energy of a single atom against the records of a usable cell list, evaluated straight out of the cells.
 *
 * The same result as a gather of batch radius zero around the atom (at the larger of the two cut-offs, the
 * atoms of 'skipMolecule' left out) followed by 'computeNeighbourhoodEnergy' -- without the gather: for a lone
 * trial position (a first bead of a CBMC growth, which shares no neighbourhood with the other first-bead
 * trials) the filtered copy of the stencil's records costs more than the pair loop over them. The cells within
 * the cut-off are visited and every bucket is evaluated in place, in SIMD lanes when the lane kernel applies
 * (see 'neighbourhoodFastKernelApplies'; a list with partially coupled atoms goes through the generic pair
 * arithmetic). The caller guarantees 'cellList' is enabled, valid and built over the atoms the energy is
 * taken against; 'background' selects how the sums are booked (a framework list records every pair).
 */
[[nodiscard]] std::optional<RunningEnergy> computeCellListAtomEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const MCCellList& cellList, const Atom& atom,
    double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap, std::optional<std::uint32_t> skipMolecule,
    Neighbourhood::Background background) noexcept;

/**
 * \brief The total intermolecular energy of 'moleculeAtoms' (the atoms the cell list is built over), every pair
 * once, evaluated straight out of the cells.
 *
 * The full-energy counterpart of 'computeCellListAtomEnergy': each atom is taken against the records of its
 * stencil with a higher atom index (the pairs within a molecule excluded), in SIMD lanes when the lane kernel
 * applies and through the generic pair arithmetic (soft-core scaling, dU/dlambda) otherwise. The same pairs,
 * cut-offs and pair terms as the pair enumeration of 'MCCellList::forEachPairOnce' with the scalar kernel. The
 * caller guarantees 'cellList' is enabled, valid and built over exactly 'moleculeAtoms'.
 */
[[nodiscard]] RunningEnergy computeCellListInterMolecularEnergy(const ForceField& forceField,
                                                                const SimulationBox& simulationBox,
                                                                const MCCellList& cellList,
                                                                std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief The total intermolecular energy of 'moleculeAtoms' without a cell list (a box too small for one, or the
 * scaled positions of a volume move), every pair once.
 *
 * The atoms are copied once into a structure-of-arrays neighbourhood ('Neighbourhood::gatherAll', in atom order)
 * and every atom is taken against the records behind its molecule: the minimum-image squared distances streamed
 * over them, the records within the cut-off compacted, then the same second pass as 'computeCellListAtomEnergy'
 * (SIMD lanes when the lane kernel applies, the generic pair arithmetic otherwise). The same pairs, cut-offs and
 * pair terms as the scalar triangular pair loop.
 */
[[nodiscard]] RunningEnergy computeAllPairsInterMolecularEnergy(const ForceField& forceField,
                                                                const SimulationBox& simulationBox,
                                                                std::span<const Atom> moleculeAtoms) noexcept;
}  // namespace Interactions
