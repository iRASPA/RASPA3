module;

export module amber_prmtop_reader;

import std;

import json;

/**
 * \brief AMBER -> RASPA: reads a 'prmtop' topology (plus, optionally, the matching 'inpcrd' / 'rst7'
 * coordinates) and produces RASPA's force_field.json, one component JSON per distinct molecule, a
 * simulation.json skeleton, and a restart JSON that seeds the coordinates.
 *
 * What is converted:
 *   - ATOM_TYPE_INDEX / AMBER_ATOM_TYPE / MASS / ATOMIC_NUMBER -> pseudo-atoms, one per distinct
 *     (Lennard-Jones type, type name, mass) combination; LENNARD_JONES_ACOEF/BCOEF -> Lennard-Jones self
 *     interactions (epsilon = B^2 / 4A, sigma = (A/B)^(1/6)) and, where the AMBER pair table deviates from
 *     Lorentz-Berthelot mixing, binary interactions. CHARGE (q * 18.2223; CHAMBER: q * sqrt(332.0716)) ->
 *     charges, per atom.
 *   - BOND_* -> HARMONIC bonds (p_0 = 2 k), ANGLE_* -> HARMONIC bends (p_0 = 2 k), DIHEDRAL_* -> one torsion
 *     per bonded quadruple: the sum of all terms with a phase of 0 or 180 degrees and periodicity <= 5 becomes a
 *     POLYNOMIAL in cos(phi) (exact, including the constant), terms with other phases or higher periodicities
 *     are appended as CVFF entries in 'ImproperTorsions'; AMBER impropers (negative fourth index) -> CVFF
 *     'ImproperTorsions'. SCEE/SCNB -> Intra14ChargeChargeScalingValue / Intra14VanDerWaalsScalingValue.
 *   - Molecules (connected components of the bond graph) are compared by types, charges and parametrised
 *     topology; each distinct molecule becomes a component. Residues named in 'rigidResidues' (TIP3P water by
 *     default) and single atoms become rigid components: water takes the ideal geometry from its equilibrium
 *     bond lengths and angle, so it is integrated as a rigid body (centre of mass + quaternion).
 *   - CMAP_COUNT / CMAP_RESOLUTION / CMAP_PARAMETER_nn / CMAP_INDEX (ff19SB; also the CHARMM_CMAP_* names of
 *     CHAMBER topologies) -> 'CMAPs' maps in force_field.json (kcal/mol -> K) and 'CMAPTorsions' terms in the
 *     component definitions.
 *   - BOX_DIMENSIONS / the inpcrd box line -> the simulation box; without a box a vacuum box is built around the
 *     molecules with direct Coulomb and a cut-off spanning every pair.
 *
 * The coordinates (inpcrd) are written to 'restart.json' as per-component position lists, referenced from
 * simulation.json via 'RestartFileName' so RASPA starts from the AMBER configuration (rigid molecules are fitted
 * to their reference geometry).
 *
 * Not converted: 10-12 hydrogen-bond terms, polarizabilities, per-atom radii/screening.
 */
export namespace AMBER
{
/// The raw sections of a prmtop file (AMBER7 format, '%FLAG' / '%FORMAT'). Indices are as in the file (1-based
/// where AMBER uses 1-based, coordinate-array offsets 3*(i-1) for the bond/angle/dihedral lists).
struct Prmtop
{
  std::string title{};
  std::vector<std::int64_t> pointers{};
  std::vector<std::string> atomNames{};
  std::vector<double> charges{};  ///< in units of e (file value / 18.2223; CHAMBER: / sqrt(332.0716))
  std::vector<std::int64_t> atomicNumbers{};
  std::vector<double> masses{};
  std::vector<std::int64_t> atomTypeIndex{};
  std::vector<std::int64_t> numberExcludedAtoms{};
  std::vector<std::int64_t> nonbondedParmIndex{};
  std::vector<std::string> residueLabels{};
  std::vector<std::int64_t> residuePointer{};
  std::vector<double> bondForceConstant{};
  std::vector<double> bondEquilValue{};
  std::vector<double> angleForceConstant{};
  std::vector<double> angleEquilValue{};
  std::vector<double> dihedralForceConstant{};
  std::vector<double> dihedralPeriodicity{};
  std::vector<double> dihedralPhase{};
  std::vector<double> sceeScaleFactor{};
  std::vector<double> scnbScaleFactor{};
  std::vector<double> lennardJonesACoef{};
  std::vector<double> lennardJonesBCoef{};
  std::vector<double> lennardJones14ACoef{};  ///< CHAMBER: the Lennard-Jones coefficients of the 1-4 pairs
  std::vector<double> lennardJones14BCoef{};
  std::vector<std::int64_t> bondsIncHydrogen{};
  std::vector<std::int64_t> bondsWithoutHydrogen{};
  std::vector<std::int64_t> anglesIncHydrogen{};
  std::vector<std::int64_t> anglesWithoutHydrogen{};
  std::vector<std::int64_t> dihedralsIncHydrogen{};
  std::vector<std::int64_t> dihedralsWithoutHydrogen{};
  std::vector<std::int64_t> excludedAtomsList{};
  std::vector<std::string> amberAtomTypes{};
  std::vector<std::int64_t> solventPointers{};
  std::vector<std::int64_t> atomsPerMolecule{};
  std::vector<double> boxDimensions{};  ///< beta, a, b, c (when IFBOX > 0)
  std::vector<std::int64_t> cmapResolution{};          ///< per CMAP map: grid points per angle
  std::vector<std::vector<double>> cmapParameters{};  ///< per CMAP map: resolution^2 energies [kcal/mol]
  std::vector<std::int64_t> cmapIndex{};               ///< per CMAP term: 5 atoms (1-based) and the map (1-based)
  std::set<std::string> flags{};        ///< every %FLAG seen

  std::size_t numberOfAtoms() const { return pointers.empty() ? 0 : static_cast<std::size_t>(pointers[0]); }
  std::size_t numberOfTypes() const { return pointers.size() > 1 ? static_cast<std::size_t>(pointers[1]) : 0; }
  bool hasBox() const { return pointers.size() > 27 && pointers[27] > 0; }
};

/// Coordinates from an ASCII inpcrd / restrt (rst7) file.
struct Coordinates
{
  std::string title{};
  std::vector<std::array<double, 3>> positions{};
  std::vector<std::array<double, 3>> velocities{};
  std::optional<std::array<double, 6>> box{};  ///< a, b, c, alpha, beta, gamma
};

struct ReadOptions
{
  /// Residue labels whose molecules are modelled as rigid bodies (centre of mass + quaternion).
  std::set<std::string> rigidResidues{"WAT", "HOH", "TIP3", "TP3", "SOL"};
  /// Van der Waals and real-space Coulomb cut-off in Angstrom for periodic systems.
  double cutOff{10.0};
  /// The van der Waals truncation written to force_field.json: "truncated" (AMBER / OpenMM default), "shifted",
  /// "switched" (OpenMM switching function), or "force-switched" (CHARMM vfswitch).
  std::string truncationMethod{"truncated"};
  /// Where the switching starts [Angstrom] for the switched truncation methods; 0: 2 Angstrom below the cut-off.
  double switchingDistance{0.0};
  /// Temperature written to simulation.json.
  double temperature{300.0};
};

struct ReadComponent
{
  std::string name{};
  nlohmann::json definition{};
  std::size_t count{};
  std::size_t atomsPerMolecule{};
  bool rigid{};
  std::vector<std::size_t> firstAtoms{};  ///< the first AMBER atom index (0-based) of every molecule
};

struct ReadResult
{
  nlohmann::json forceField{};
  std::vector<ReadComponent> components{};
  nlohmann::json simulation{};
  nlohmann::json restart{};  ///< 'SimulationBox' + per-component positions (empty without coordinates)
  std::array<double, 3> boxLengths{};
  std::array<double, 3> boxAngles{90.0, 90.0, 90.0};
  std::vector<std::array<double, 3>> positions{};  ///< all atoms in file order, Angstrom
  std::vector<std::string> warnings{};
};

Prmtop parsePrmtop(const std::filesystem::path& file);
Coordinates readCoordinates(const std::filesystem::path& file, std::size_t numberOfAtoms);

ReadResult readPrmtop(const std::filesystem::path& prmtopFile,
                      const std::optional<std::filesystem::path>& coordinateFile, const ReadOptions& options = {});

/// Writes force_field.json, <component>.json files, simulation.json and restart.json into 'directory'.
void writeRaspaInput(const ReadResult& result, const std::filesystem::path& directory);
}  // namespace AMBER
