module;

export module lammps_topology;

import std;

import double3;
import int3;
import atom;
import atom_dynamics;
import molecule;
import component;
import simulationbox;
import forcefield;
import framework;
import polymer_links;
import lammps_styles;

/**
 * \brief Intermediate representation of a RASPA system in LAMMPS terms.
 *
 * The builder walks the components, framework and force field once and produces everything the file
 * writers need: deduplicated bonded types (one LAMMPS type per distinct (style, coefficients)), the
 * topology lists in LAMMPS atom ids, wrapped coordinates with image flags, per-atom rigid-fragment ids,
 * pair coefficients per type pair, table bodies for forms without an analytic LAMMPS style, and the
 * settings that only the input script can carry (styles, special_bonds, cut-offs, kspace, constraints).
 *
 * A polymer overlay (the bonds formed by a polymerization, see PolymerLinkTable) is merged into the
 * topology when given: the polymer bonds with their bends and torsions are added to the lists, and the
 * bonded terms of the monomers whose parameters changed with the retyping take the product parameters.
 * LAMMPS derives the 1-2 / 1-3 exclusions and the 1-4 scaling across the new bonds from the bond list
 * itself (special_bonds), so the pair corrections of the overlay need no entries. The force ramps of
 * fresh bonds and the stereocentre restraints have no LAMMPS equivalent: the terms are written at full
 * strength and the restraints are dropped (both reported as warnings).
 */
export namespace LAMMPS
{
struct ExportOptions
{
  TableGrid grid{};
  std::string dataFile{"raspa.data"};
  std::string tableFile{"raspa.table"};
  std::string pairListFile{"raspa.pairs"};
  std::string comment{};  ///< Appended to the first (comment) line of the data file.
};

/// Distinct LAMMPS types of one bonded class, in order of first appearance.
struct TypeTable
{
  std::vector<Term> types{};
  std::map<std::string, std::size_t> index{};

  /// Adds (or finds) the term and returns its 1-based LAMMPS type id.
  std::size_t add(const Term &term);
  std::vector<std::string> styles() const;
  bool hybrid() const { return styles().size() > 1; }
  bool empty() const { return types.empty(); }
  bool hasStyle(std::string_view style) const;
};

struct BondedEntry
{
  std::size_t type{};                     ///< 1-based LAMMPS type id
  std::array<std::size_t, 4> atoms{};     ///< 1-based LAMMPS atom ids (trailing entries unused)
};

struct AtomEntry
{
  std::size_t id{};        ///< 1-based LAMMPS atom id
  std::size_t molecule{};  ///< 1-based LAMMPS molecule id
  std::size_t type{};      ///< 1-based LAMMPS atom type
  double charge{};
  double3 position{};      ///< wrapped, Angstrom
  int3 image{};
  double3 velocity{};      ///< Angstrom / fs
  std::size_t fragment{};  ///< rigid-body id (0: not part of a rigid body)
  bool fixed{false};       ///< frozen in the laboratory frame
};

struct PairEntry
{
  std::size_t i{}, j{};  ///< 1-based atom types, i <= j
  PairTerm term{};
  std::string tableKeyword{};
};

struct MoleculeGroup
{
  std::string name{};
  std::size_t firstMolecule{}, lastMolecule{};  ///< inclusive 1-based LAMMPS molecule ids
  bool rigid{false};
  bool hasBonds{false};
};

struct SpecialBonds
{
  bool present{false};   ///< at least one flexible molecule with bonds
  bool uniform{true};    ///< every component follows the 1-2/1-3 excluded, 1-4 scaled, rest full pattern
  double vdw14{0.0};
  double coul14{0.0};
};

struct Topology
{
  ExportOptions options{};
  std::vector<std::string> warnings{};

  // box (Angstrom)
  double3 lengths{};
  double xy{}, xz{}, yz{};
  bool triclinic{false};

  std::vector<std::pair<std::string, double>> masses{};  ///< per atom type: name, mass
  std::vector<AtomEntry> atoms{};

  TypeTable bonds{}, angles{}, dihedrals{}, impropers{};
  std::vector<BondedEntry> bondList{}, angleList{}, dihedralList{}, improperList{};

  std::vector<PairEntry> pairs{};
  std::vector<std::pair<std::string, std::string>> tables{};  ///< keyword, body
  std::vector<std::string> pairList{};                        ///< 'pair_style list' lines (1-4 pairs with non-uniform scaling)

  // force-field settings
  bool useCharge{true};
  ForceField::ChargeMethod chargeMethod{ForceField::ChargeMethod::Ewald};
  double cutOffVDW{12.0};
  double cutOffFrameworkVDW{12.0};
  double cutOffCoulomb{12.0};
  double ewaldPrecision{1e-6};
  double alpha{0.0};
  bool anyTail{false}, anyShift{false}, allShift{false};
  ForceField::MixingRule mixingRule{ForceField::MixingRule::Lorentz_Berthelot};
  SpecialBonds special{};

  // constraints and rigid bodies
  std::vector<std::size_t> shakeBondTypes{}, shakeAngleTypes{};
  std::size_t numberOfFragments{0};
  bool hasFixedAtoms{false};
  bool hasFramework{false};
  std::vector<MoleculeGroup> groups{};

  // polymer overlay (zero without one)
  std::size_t polymerBonds{0};      ///< Polymer bonds merged into the bond list.
  std::size_t polymerOverrides{0};  ///< Monomer terms that took the product parameters.
  std::size_t polymerBrokenTerms{0};  ///< Monomer bonds, angles and dihedrals removed by broken bonds.

  std::vector<std::string> pairStyles() const;
  bool pairCoefficientsInDataFile() const;  ///< single analytic VDW style: PairIJ Coeffs can live in the data file
};

/// 'polymerLinks': the polymer overlay to merge (its derived term lists must be current, see
/// PolymerLinkTable::rebuild); none to export the component topology only.
Topology buildTopology(std::span<const Component> components, std::span<const Atom> atomData,
                       std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                       const SimulationBox &simulationBox, const ForceField &forceField,
                       std::span<const std::size_t> numberOfIntegerMoleculesPerComponent,
                       const std::optional<Framework> &framework, ExportOptions options = {},
                       const PolymerLinkTable *polymerLinks = nullptr);
}  // namespace LAMMPS
