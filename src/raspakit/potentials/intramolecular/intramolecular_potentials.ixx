module;

export module intra_molecular_potentials;

import std;

import archive;
import double3x3;
import atom;
import atom_dynamics;
import chiral_center;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import cmap_potential;
import van_der_waals_potential;
import coulomb_potential;
import running_energy;
import forcefield;
import vdwparameters;
import simulationbox;
import intra_molecular_exclusions;

export namespace Potentials
{
/**
 * \brief The data to materialise an implicit non-bonded pair of a molecule as an explicit pair term.
 *
 * The implicit pairs (see IntraMolecularPotentials) are evaluated with the force field; the CBMC growth plans,
 * which are built without the force field, need explicit VanDerWaalsPotential / CoulombPotential terms per growth
 * step (their pair potentials and charges feed the lookahead guide). This table captures what those terms need:
 * the force field's pair potential (VDWParameters, any functional form) per pair of atom types of the molecule and
 * the charge per atom. It is linear in the molecule size.
 */
struct ImplicitPairParameters
{
  std::uint64_t versionNumber{3};  ///< Version number for serialization.

  std::vector<std::uint32_t> typeIndex{};           ///< Per atom: the index of its type in 'parameters'.
  std::vector<double> charge{};                     ///< Per atom: the charge.
  std::size_t numberOfTypes{0};                     ///< The number of distinct atom types of the molecule.
  std::vector<VDWParameters> parameters{};          ///< [a * numberOfTypes + b]: the force-field pair potential.
  std::vector<VDWParameters> parameters14{};        ///< [a * numberOfTypes + b]: the pair potential of the 1-4 pairs.
  bool useCharge{false};                            ///< Whether the force field uses charges (Coulomb terms exist).

  bool operator==(const ImplicitPairParameters &) const = default;

  [[nodiscard]] bool empty() const { return typeIndex.empty(); }

  /// The explicit van der Waals term of the pair (A, B) with the given scaling (a 1-4-parameter pair gets the
  /// 1-4 pair potential).
  [[nodiscard]] VanDerWaalsPotential vanDerWaalsTerm(std::size_t A, std::size_t B, double scaling, bool pair14) const;
  /// The explicit Coulomb term of the pair (A, B) with the given scaling.
  [[nodiscard]] CoulombPotential coulombTerm(std::size_t A, std::size_t B, double scaling) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ImplicitPairParameters &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ImplicitPairParameters &p);
};

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ImplicitPairParameters &p);
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ImplicitPairParameters &p);

/**
 * \brief The intramolecular potentials of a component (or of one growth step, or of a flexible framework).
 *
 * The bonded terms carry their own parameters. The non-bonded pairs come in two forms:
 *
 *  - the implicit pairs of a whole molecule: every pair of atoms that is not excluded ('exclusions', see
 *    IntraMolecularExclusions) interacts, the scaled (1-4) pairs with their factors, every other pair at full
 *    strength. They are enumerated, never stored, so a component costs memory linear in its size;
 *  - the explicit pair lists 'vanDerWaals' and 'coulombs': the pair terms of a CBMC growth step (materialised from
 *    the implicit pairs of the component by 'filteredInteractions') and of a flexible framework (which lists its
 *    pairs with periodic image shifts).
 *
 * A component has only implicit pairs, a growth step or framework only explicit ones; the evaluation routines
 * visit both ('forEachVanDerWaalsPair', 'forEachCoulombPair'). Either way the energies are evaluated with the
 * regular force-field pair potentials (the same truncated/shifted van der Waals potential and the same Ewald /
 * shifted Coulomb potential as between molecules, under the minimum-image convention) times the pair scaling, which
 * is why the evaluation routines take the force field and the box. The (epsilon, sigma) and charges stored in the
 * explicit pair entries (and in 'implicitParameters' for materialising them) only serve the CBMC lookahead guide (a
 * positive bias that needs no exactness) and the printing.
 */
struct IntraMolecularPotentials
{
  std::uint64_t versionNumber{3};  ///< Version number for serialization.

  IntraMolecularExclusions exclusions{};        ///< The exclusion topology defining the implicit non-bonded pairs.
  ImplicitPairParameters implicitParameters{};  ///< The parameters to materialise the implicit pairs.

  std::vector<ChiralCenter> chiralCenters{};               ///< List of chiral centers in the component.
  std::vector<BondPotential> bonds{};                      ///< List of bond potentials.
  std::vector<UreyBradleyPotential> ureyBradleys{};        ///< List of Urey-Bradley potentials.
  std::vector<BendPotential> bends{};                      ///< List of bend potentials.
  std::vector<InversionBendPotential> inversionBends{};    ///< List of inversion-bend potentials.
  std::vector<OutOfPlaneBendPotential> outOfPlaneBends{};  ///< List of out-of-plane-bend potentials.
  std::vector<TorsionPotential> torsions{};                ///< List of torsion potentials.
  std::vector<TorsionPotential> improperTorsions{};        ///< List of improper torsion potentials.
  std::vector<BondBondPotential> bondBonds{};              ///< List of bond-bond potentials.
  std::vector<BondBendPotential> bondBends{};              ///< List of bond-bend potentials.
  std::vector<BondTorsionPotential> bondTorsions{};        ///< List of bond-torsion potentials.
  std::vector<BendBendPotential> bendBends{};              ///< List of bend-bend potentials.
  std::vector<BendTorsionPotential> bendTorsions{};        ///< List of bend-torsion potentials.
  std::vector<CMAPMap> cmapMaps{};                         ///< The CMAP correction maps the 'cmaps' terms refer to.
  std::vector<CMAPPotential> cmaps{};                      ///< List of CMAP (phi, psi) correction terms.
  std::vector<VanDerWaalsPotential> vanDerWaals{};         ///< The explicit van der Waals pair terms.
  std::vector<CoulombPotential> coulombs{};                ///< The explicit Coulomb pair terms.

  /// Visits every non-bonded van der Waals pair as 'function(A, B, scaling, pair14)': the explicit terms, then
  /// the implicit pairs ('pair14': the pair uses the 1-4 pair parameters of the force field).
  template <typename Function>
  void forEachVanDerWaalsPair(Function &&function) const
  {
    for (const VanDerWaalsPotential &pair : vanDerWaals)
    {
      function(pair.identifiers[0], pair.identifiers[1], pair.scaling, pair.pair14);
    }
    exclusions.forEachNonExcludedPair([&](std::size_t A, std::size_t B, double scalingVDW, double, bool pair14)
                                      { function(A, B, scalingVDW, pair14); });
  }

  /// Visits every non-bonded Coulomb pair as 'function(A, B, scaling)': the explicit terms, then the implicit
  /// pairs (when the force field the component was built with uses charges).
  template <typename Function>
  void forEachCoulombPair(Function &&function) const
  {
    for (const CoulombPotential &pair : coulombs)
    {
      function(pair.identifiers[0], pair.identifiers[1], pair.scaling);
    }
    if (!implicitParameters.useCharge) return;
    exclusions.forEachNonExcludedPair([&](std::size_t A, std::size_t B, double, double scalingCoulomb, bool)
                                      { function(A, B, scalingCoulomb); });
  }

  /// Visits every non-bonded van der Waals pair as an explicit term, 'function(const VanDerWaalsPotential &)';
  /// the implicit pairs are materialised from 'implicitParameters'.
  template <typename Function>
  void forEachVanDerWaalsTerm(Function &&function) const
  {
    for (const VanDerWaalsPotential &pair : vanDerWaals) function(pair);
    exclusions.forEachNonExcludedPair([&](std::size_t A, std::size_t B, double scalingVDW, double, bool pair14)
                                      { function(implicitParameters.vanDerWaalsTerm(A, B, scalingVDW, pair14)); });
  }

  /// Visits every non-bonded Coulomb pair as an explicit term, 'function(const CoulombPotential &)'.
  template <typename Function>
  void forEachCoulombTerm(Function &&function) const
  {
    for (const CoulombPotential &pair : coulombs) function(pair);
    if (!implicitParameters.useCharge) return;
    exclusions.forEachNonExcludedPair([&](std::size_t A, std::size_t B, double, double scalingCoulomb, bool)
                                      { function(implicitParameters.coulombTerm(A, B, scalingCoulomb)); });
  }

  [[nodiscard]] std::size_t numberOfVanDerWaalsPairs() const
  {
    return vanDerWaals.size() + exclusions.numberOfNonExcludedPairs();
  }
  [[nodiscard]] std::size_t numberOfCoulombPairs() const
  {
    return coulombs.size() + (implicitParameters.useCharge ? exclusions.numberOfNonExcludedPairs() : 0);
  }
  [[nodiscard]] bool hasNonBondedPairs() const { return numberOfVanDerWaalsPairs() + numberOfCoulombPairs() > 0; }

  std::optional<BondPotential> findBondPotential(std::size_t A, std::size_t B) const;

  double calculateBondSmallMCEnergies(const std::span<Atom> atoms) const;
  double calculateBendSmallMCEnergies(const std::span<Atom> atoms) const;

  double calculateTorsionEnergies(const std::span<Atom> atoms) const;

  RunningEnergy computeInternalEnergies(const ForceField &forceField, const SimulationBox &simulationBox,
                                        const std::span<const Atom> atoms) const;

  /**
   * \brief Computes the internal interactions that are not sampled during CBMC / recoil-growth.
   *
   * During growing/retracing only the bond, bend, and torsion potentials are sampled, and the
   * intramolecular van der Waals and Coulomb interactions enter through the selection of the
   * beads. This routine computes all remaining internal interactions: Urey-Bradley,
   * inversion-bend, out-of-plane-bend, improper torsion, and the cross-terms (bond-bond,
   * bond-bend, bond-torsion, bend-bend, bend-torsion). The Boltzmann factor of this energy is
   * used to correct the Rosenbluth weight of the grown/retraced chain.
   */
  RunningEnergy computeInternalEnergiesNotSampledDuringGrowth(const std::span<const Atom> atoms) const;

  /**
   * \brief Scales the bonded intramolecular Hamiltonian by 'lambda' (solute tempering).
   *
   * The bonded terms (bond, bend, torsion, ... and the cross terms) have their energy parameters
   * multiplied by lambda. The intramolecular van der Waals and Coulomb pairs are evaluated with the
   * force field and the atomic charges, which the solute tempering scales itself (the force field's
   * solute-solute epsilons by lambda, the solute charges by sqrt(lambda)); only the (epsilon, sigma)
   * and charges stored for the CBMC lookahead guide are rescaled here to keep its bias in step.
   */
  void scaleEnergy(double lambda);

  RunningEnergy computeInternalBondEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalUreyBradleyEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalInversionBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalOutOfPlaneBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalImproperTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondBondEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalCMAPEnergies(const std::span<const Atom> atoms) const;

  /// Appends a CMAP term over the atoms (A, B, C, D, E) with the map named 'mapName' from 'cmapMaps'; throws when
  /// the map is unknown.
  void addCMAP(const std::array<std::size_t, 5> &identifiers, const std::string &mapName);

  /**
   * \brief Intramolecular van der Waals energy of the listed pairs: the regular force-field pair potential (with its
   * cutoff and shift) under the minimum-image convention, times the pair scaling.
   */
  RunningEnergy computeInternalIntraVanDerWaalsEnergies(const ForceField &forceField,
                                                        const SimulationBox &simulationBox,
                                                        const std::span<const Atom> atoms) const;

  /**
   * \brief Intramolecular Coulomb energy of the listed pairs (see Potentials::intraMolecularCoulomb): the bare
   * Coulomb interaction scaled by the pair scaling plus the exclusion kernel of the charge method, inside the
   * Coulomb cutoff, under the minimum-image convention. The lambda derivative of the kernel term goes to
   * 'dudlambdaEwald'.
   */
  RunningEnergy computeInternalIntraCoulombEnergies(const ForceField &forceField, const SimulationBox &simulationBox,
                                                    const std::span<const Atom> atoms) const;

  /**
   * \brief Intramolecular van der Waals and Coulomb energy of the beads.
   *
   * Convenience sum of computeInternalIntraVanDerWaalsEnergies and computeInternalIntraCoulombEnergies.
   * These are the intramolecular non-bonded terms that enter the CBMC / recoil-growth bead selection.
   */
  RunningEnergy computeInternalIntraVanDerWaalsAndCoulombEnergies(const ForceField &forceField,
                                                                  const SimulationBox &simulationBox,
                                                                  const std::span<const Atom> atoms) const;

  RunningEnergy computeInternalGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                        std::span<const Atom> atoms, std::span<AtomDynamics> dynamics) const;

  std::pair<RunningEnergy, double3x3> computeInternalStrainDerivative(const ForceField &forceField,
                                                                      const SimulationBox &simulationBox,
                                                                      std::span<const Atom> atoms,
                                                                      std::span<AtomDynamics> dynamics) const;

  /// The bonded terms only (bonds, bends, torsions, cross terms; no van der Waals / Coulomb pairs): for engines
  /// that evaluate the non-excluded pairs of a molecule in their pair loops.
  RunningEnergy computeInternalBondedGradient(const SimulationBox &simulationBox, std::span<const Atom> atoms,
                                              std::span<AtomDynamics> dynamics) const;
  std::pair<RunningEnergy, double3x3> computeInternalBondedStrainDerivative(const SimulationBox &simulationBox,
                                                                            std::span<const Atom> atoms,
                                                                            std::span<AtomDynamics> dynamics) const;

  /// The van der Waals / Coulomb pairs only ('vanDerWaals', 'coulombs').
  std::pair<RunningEnergy, double3x3> computeInternalNonBondedStrainDerivative(const ForceField &forceField,
                                                                               const SimulationBox &simulationBox,
                                                                               std::span<const Atom> atoms,
                                                                               std::span<AtomDynamics> dynamics) const;

  /**
   * \brief The terms of one growth step: every term whose atoms are all placed or to be placed, with at least one
   * to be placed. The non-bonded pairs of the step are explicit terms (the implicit pairs of the molecule are
   * materialised through 'implicitParameters'); the result has no implicit pairs of its own.
   */
  Potentials::IntraMolecularPotentials filteredInteractions(std::size_t numberOfBeads,
                                                            const std::span<std::size_t> beadsAlreadyPlaced,
                                                            const std::span<std::size_t> beadsToBePlaced) const;

  std::string printStatus() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive,
                                            const Potentials::IntraMolecularPotentials& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Potentials::IntraMolecularPotentials& p);
};

  Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive,
                                          const Potentials::IntraMolecularPotentials& p);
  Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Potentials::IntraMolecularPotentials& p);
}  // namespace Potentials

/*
export template <>
struct std::formatter<Potentials::IntraMolecularPotentials>: std::formatter<string_view>
{
  auto format(const Potentials::IntraMolecularPotentials& p, std::format_context& ctx) const
  {
    std::string temp{};
    std::format_to(std::back_inserter(temp), "(bonds: {}-{})", p.bonds);
    return std::formatter<string_view>::format(temp, ctx);
  }
};
*/
