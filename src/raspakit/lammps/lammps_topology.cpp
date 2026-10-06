module;

module lammps_topology;

import std;

import double3;
import double3x3;
import int3;
import units;
import atom;
import atom_dynamics;
import molecule;
import component;
import simulationbox;
import forcefield;
import vdwparameters;
import pseudo_atom;
import framework;
import fragment;
import fragment_graph;
import connectivity_table;
import intra_molecular_potentials;
import intra_molecular_exclusions;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import torsion_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import lammps_styles;

namespace
{
using namespace LAMMPS;

template <typename Range>
std::string parameterKey(std::string_view prefix, std::size_t type, const Range &parameters)
{
  std::string key = std::format("{}|{}", prefix, type);
  for (double p : parameters) key += std::format("|{:.12g}", p);
  return key;
}

/// Tables are shared between identical potentials; the keyword is the class prefix plus a running index.
struct TableRegistry
{
  std::map<std::string, std::string> keywords{};
  std::size_t counter{0};

  std::string keywordFor(Topology &topology, const std::string &key, std::string_view prefix,
                         const std::function<std::string(std::string_view)> &body)
  {
    auto it = keywords.find(key);
    if (it != keywords.end()) return it->second;
    std::string keyword = std::format("{}{}", prefix, ++counter);
    keywords.emplace(key, keyword);
    topology.tables.emplace_back(keyword, body(keyword));
    return keyword;
  }
};

std::array<std::size_t, 4> sortedKey(std::span<const std::size_t> ids)
{
  std::array<std::size_t, 4> key{0, 0, 0, 0};
  for (std::size_t i = 0; i < ids.size() && i < 4; ++i) key[i] = ids[i];
  return key;
}

/// Undirected key for a triple (A,B,C) ~ (C,B,A) and a quadruple (A,B,C,D) ~ (D,C,B,A).
std::array<std::size_t, 4> canonicalKey(std::span<const std::size_t> ids)
{
  std::array<std::size_t, 4> key = sortedKey(ids);
  std::array<std::size_t, 4> reversed = key;
  std::reverse(reversed.begin(), reversed.begin() + static_cast<std::ptrdiff_t>(ids.size()));
  return std::min(key, reversed);
}

/// Adds the bonded terms of one molecule instance (or of the framework) to the topology.
///
/// 'globalId' maps a template atom index onto the 1-based LAMMPS atom id of this instance. Types are
/// deduplicated through the topology's type tables, so calling this once per instance is cheap.
void addBondedTerms(Topology &topology, TableRegistry &tables, const Potentials::IntraMolecularPotentials &intra,
                    const std::function<std::size_t(std::size_t)> &globalId, std::string_view label,
                    bool warn)
{
  auto warning = [&](std::string message)
  {
    if (warn) topology.warnings.push_back(std::format("{}: {}", label, message));
  };

  auto ids2 = [&](std::span<const std::size_t> ids, const std::array<std::size_t, 4> &order = {0, 1, 2, 3})
  {
    std::array<std::size_t, 4> out{0, 0, 0, 0};
    for (std::size_t i = 0; i < ids.size(); ++i) out[i] = globalId(ids[order[i]]);
    return out;
  };

  // ---- bonds --------------------------------------------------------------------------------------
  for (const BondPotential &bond : intra.bonds)
  {
    Term term = bondTerm(bond);
    if (term.style == "table")
    {
      term.tableKeyword = tables.keywordFor(
          topology, parameterKey("B", static_cast<std::size_t>(bond.type), bond.parameters), "B",
          [&](std::string_view keyword) { return bondTableBody(bond, keyword, topology.options.grid); });
    }
    const std::size_t type = topology.bonds.add(term);
    if (bond.type == BondType::Fixed &&
        std::find(topology.shakeBondTypes.begin(), topology.shakeBondTypes.end(), type) ==
            topology.shakeBondTypes.end())
    {
      topology.shakeBondTypes.push_back(type);
    }
    topology.bondList.push_back({type, ids2(bond.identifiers)});
  }

  // ---- cross terms indexed by their central triple / quadruple -----------------------------------------
  std::map<std::array<std::size_t, 4>, const BondBondPotential *> bondBondByTriple{};
  std::map<std::array<std::size_t, 4>, const BondBendPotential *> bondBendByTriple{};
  std::map<std::array<std::size_t, 4>, const BendTorsionPotential *> bendTorsionByQuad{};
  for (const BondBondPotential &term : intra.bondBonds)
  {
    bondBondByTriple[canonicalKey(std::span<const std::size_t>(term.identifiers))] = &term;
  }
  for (const BondBendPotential &term : intra.bondBends)
  {
    bondBendByTriple[canonicalKey(std::span<const std::size_t>(term.identifiers).first(3))] = &term;
  }
  for (const BendTorsionPotential &term : intra.bendTorsions)
  {
    bendTorsionByQuad[canonicalKey(std::span<const std::size_t>(term.identifiers))] = &term;
  }
  std::set<const void *> consumedCross{};

  // ---- bends (with class2 promotion when a CFF bond-bond / bond-bend cross term sits on the triple) ----
  for (const BendPotential &bend : intra.bends)
  {
    std::span<const std::size_t> ids(bend.identifiers);
    const std::array<std::size_t, 4> key = canonicalKey(ids);
    const BondBondPotential *bondBond = bondBondByTriple.contains(key) ? bondBondByTriple.at(key) : nullptr;
    const BondBendPotential *bondBend = bondBendByTriple.contains(key) ? bondBendByTriple.at(key) : nullptr;

    Term term = bendTerm(bend);
    std::array<std::size_t, 3> order{0, 1, 2};
    if (bondBond || bondBend)
    {
      std::optional<std::string> class2 = bendClass2(bend);
      std::optional<std::string> bb = bondBond ? bondBondClass2(*bondBond) : std::optional<std::string>{};
      std::optional<std::string> ba =
          bondBend ? bondBendClass2(*bondBend, bend.parameters[1]) : std::optional<std::string>{};
      if (class2.has_value() && (!bondBond || bb.has_value()) && (!bondBend || ba.has_value()))
      {
        // the cross term fixes the A-B-C orientation (r1 belongs to the first bond)
        const std::array<std::size_t, 3> crossIds = bondBond
                                                        ? bondBond->identifiers
                                                        : std::array<std::size_t, 3>{bondBend->identifiers[0],
                                                                                     bondBend->identifiers[1],
                                                                                     bondBend->identifiers[2]};
        if (crossIds[0] != bend.identifiers[0]) order = {2, 1, 0};
        term = Term{};
        term.style = "class2";
        term.coefficients = *class2;
        term.extra["bb"] = bb.value_or("0.0 0.0 0.0");
        term.extra["ba"] = ba.value_or("0.0 0.0 0.0 0.0");
        if (bondBend && std::abs(bondBend->parameters[0] - bend.parameters[1]) > 1e-8)
        {
          term.exact = false;
          term.note = "class2 uses the angle's theta0 for the bond-angle cross term (RASPA value differs)";
        }
        if (bondBond) consumedCross.insert(bondBond);
        if (bondBend) consumedCross.insert(bondBend);
      }
      else
      {
        warning("bond-bond / bond-bend cross term on a bend that is not of class2 form; cross term dropped");
      }
    }
    if (term.style == "table")
    {
      term.tableKeyword = tables.keywordFor(
          topology, parameterKey("A", static_cast<std::size_t>(bend.type), bend.parameters), "A",
          [&](std::string_view keyword) { return angleTableBody(bend, keyword, topology.options.grid); });
    }
    const std::size_t type = topology.angles.add(term);
    if (bend.type == BendType::Fixed &&
        std::find(topology.shakeAngleTypes.begin(), topology.shakeAngleTypes.end(), type) ==
            topology.shakeAngleTypes.end())
    {
      topology.shakeAngleTypes.push_back(type);
    }
    topology.angleList.push_back({type, ids2(ids, {order[0], order[1], order[2], 3})});
  }

  // ---- Urey-Bradley: an angle of style charmm with K = 0 on the triple A-B-C ---------------------------
  for (const UreyBradleyPotential &ureyBradley : intra.ureyBradleys)
  {
    const std::size_t A = ureyBradley.identifiers[0], C = ureyBradley.identifiers[1];
    std::optional<std::size_t> B{};
    for (const BendPotential &bend : intra.bends)
    {
      if ((bend.identifiers[0] == A && bend.identifiers[2] == C) ||
          (bend.identifiers[0] == C && bend.identifiers[2] == A))
      {
        B = bend.identifiers[1];
        break;
      }
    }
    if (!B.has_value())
    {
      warning(std::format("Urey-Bradley {}-{} has no matching bend; dropped", A, C));
      continue;
    }
    Term term = ureyBradleyTerm(ureyBradley);
    const std::size_t type = topology.angles.add(term);
    topology.angleList.push_back({type, {globalId(A), globalId(*B), globalId(C), 0}});
  }

  // ---- torsions (class2 promotion when a CFF bend-torsion cross term sits on the quadruple) ------------
  for (const TorsionPotential &torsion : intra.torsions)
  {
    std::span<const std::size_t> ids(torsion.identifiers);
    const std::array<std::size_t, 4> key = canonicalKey(ids);
    const BendTorsionPotential *bendTorsion = bendTorsionByQuad.contains(key) ? bendTorsionByQuad.at(key) : nullptr;

    Term term = torsionTerm(torsion);
    std::array<std::size_t, 4> order{0, 1, 2, 3};
    if (bendTorsion)
    {
      std::optional<std::string> class2 = torsionClass2(torsion);
      std::optional<std::string> aat = bendTorsionClass2(*bendTorsion);
      if (class2.has_value() && aat.has_value())
      {
        if (bendTorsion->identifiers[0] != torsion.identifiers[0]) order = {3, 2, 1, 0};
        term = Term{};
        term.style = "class2";
        term.coefficients = *class2;
        term.extra["mbt"] = "0.0 0.0 0.0 0.0";
        term.extra["ebt"] = "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0";
        term.extra["at"] = "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0";
        term.extra["aat"] = *aat;
        term.extra["bb13"] = "0.0 0.0 0.0";
        consumedCross.insert(bendTorsion);
      }
      else
      {
        warning("bend-torsion cross term on a torsion that is not of class2 (CFF) form; cross term dropped");
      }
    }
    if (term.style == "table")
    {
      term.tableKeyword = tables.keywordFor(
          topology, parameterKey("D", static_cast<std::size_t>(torsion.type), torsion.parameters), "D",
          [&](std::string_view keyword) { return dihedralTableBody(torsion, keyword, topology.options.grid); });
    }
    const std::size_t type = topology.dihedrals.add(term);
    topology.dihedralList.push_back({type, ids2(ids, order)});
  }

  // ---- impropers: improper torsions, inversion bends, out-of-plane bends, CFF bend-bend ---------------
  for (const TorsionPotential &improper : intra.improperTorsions)
  {
    for (const Term &term : improperTorsionTerms(improper))
    {
      const std::size_t type = topology.impropers.add(term);
      topology.improperList.push_back({type, ids2(improper.identifiers)});
    }
  }
  for (const InversionBendPotential &inversion : intra.inversionBends)
  {
    Term term = inversionBendTerm(inversion);
    const std::size_t type = topology.impropers.add(term);
    topology.improperList.push_back({type, ids2(inversion.identifiers, term.atomOrder)});
  }
  for (const OutOfPlaneBendPotential &outOfPlane : intra.outOfPlaneBends)
  {
    Term term = outOfPlaneBendTerm(outOfPlane);
    const std::size_t type = topology.impropers.add(term);
    topology.improperList.push_back({type, ids2(outOfPlane.identifiers)});
  }
  for (const BendBendPotential &bendBend : intra.bendBends)
  {
    std::optional<std::string> aa = bendBendClass2(bendBend);
    Term term{};
    if (aa.has_value())
    {
      term.style = "class2";
      term.coefficients = "0.0 0.0";
      term.extra["aa"] = *aa;
    }
    else
    {
      term.style = "zero";
      term.exact = false;
      term.note = "bend-bend cross term without class2 form dropped";
    }
    const std::size_t type = topology.impropers.add(term);
    topology.improperList.push_back({type, ids2(bendBend.identifiers)});
  }

  // ---- what could not be mapped ---------------------------------------------------------------------
  for (const BondBondPotential &term : intra.bondBonds)
  {
    if (!consumedCross.contains(&term)) warning("bond-bond cross term without matching bend dropped");
  }
  for (const BondBendPotential &term : intra.bondBends)
  {
    if (!consumedCross.contains(&term)) warning("bond-bend cross term without matching bend dropped");
  }
  for (const BendTorsionPotential &term : intra.bendTorsions)
  {
    if (!consumedCross.contains(&term)) warning("bend-torsion cross term without matching torsion dropped");
  }
  if (!intra.bondTorsions.empty()) warning("bond-torsion (MM3) cross terms have no LAMMPS equivalent; dropped");
  if (!intra.cmaps.empty())
    warning(std::format("{} CMAP terms are not exported (LAMMPS needs 'fix cmap' with a separate map file); dropped",
                        intra.cmaps.size()));
}

/// Bond-graph separations up to 'maximum' from every atom (BFS).
std::vector<std::vector<std::size_t>> bondSeparations(const ConnectivityTable &connectivity, std::size_t n,
                                                      std::size_t maximum)
{
  std::vector<std::vector<std::size_t>> result(n, std::vector<std::size_t>(n, maximum + 1));
  for (std::size_t start = 0; start < n; ++start)
  {
    std::vector<std::size_t> &row = result[start];
    row[start] = 0;
    std::deque<std::size_t> queue{start};
    while (!queue.empty())
    {
      const std::size_t current = queue.front();
      queue.pop_front();
      if (row[current] >= maximum) continue;
      for (std::size_t neighbor : connectivity.findAllNeighbors(current))
      {
        if (row[neighbor] > row[current] + 1)
        {
          row[neighbor] = row[current] + 1;
          queue.push_back(neighbor);
        }
      }
    }
  }
  return result;
}

struct ComponentSpecial
{
  bool present{false};
  bool standard{true};
  double vdw14{0.0}, coul14{0.0};
  std::vector<IntraMolecularExclusions::ScaledPair> nonStandard14{};  ///< 1-4 pairs with a non-uniform VDW scaling
};

/// Compares a component's intramolecular non-bonded model with LAMMPS's special_bonds pattern.
///
/// RASPA excludes 1-2, 1-3 and same-rigid-fragment pairs, scales the 1-4 pairs with the component's 1-4 scaling
/// (or a per-pair override) and computes every other pair at full strength. LAMMPS excludes 1-2 and 1-3 via
/// special_bonds, scales 1-4 with a single global factor and computes everything else at full strength. The
/// deviations are all in the scaled pairs of the exclusion topology.
ComponentSpecial analyseSpecial(const Component &component, std::vector<std::string> &warnings)
{
  ComponentSpecial result{};
  const Potentials::IntraMolecularPotentials &intra = component.intraMolecularPotentials;
  if (component.rigid || intra.bonds.empty()) return result;
  result.present = true;

  const std::size_t n = component.atoms.size();
  const std::vector<std::vector<std::size_t>> separation = bondSeparations(component.connectivityTable, n, 4);
  const bool charged =
      std::any_of(component.atoms.begin(), component.atoms.end(), [](const Atom &a) { return std::abs(a.charge) > 0.0; });

  // every 1-4 pair carries the component's 1-4 scaling unless overridden (an override equal to (1, 1) is not a
  // scaled pair: then the 1-4 scaling is non-uniform as well)
  const std::vector<std::array<std::size_t, 2>> pairs14 = IntraMolecularExclusions::pairs14(component.connectivityTable);
  const double vdw14 = component.intra14VanDerWaalsScaling;
  const double coul14 = component.intra14ChargeChargeScaling;
  bool inconsistent = false;

  for (const IntraMolecularExclusions::ScaledPair &pair : intra.exclusions.scaledPairs)
  {
    const std::size_t A = pair.atomA;
    const std::size_t B = pair.atomB;
    const std::size_t s = separation[A][B];
    // a pair with its own 1-4 pair parameters is not a special_bonds scaling of the regular pair
    if (pair.pair14) result.standard = false;
    if (s == 3)
    {
      if (std::abs(vdw14 - pair.scalingVDW) > 1e-12) result.standard = false;
      if (charged && std::abs(coul14 - pair.scalingCoulomb) > 1e-12)
      {
        inconsistent = true;
        warnings.push_back(std::format("{}: non-uniform 1-4 Coulomb scaling; special_bonds is global", component.name));
      }
    }
    else
    {
      if (std::abs(pair.scalingVDW - 1.0) > 1e-12)
      {
        inconsistent = true;
        warnings.push_back(std::format("{}: VDW pair {}-{} beyond 1-4 has scaling {}; LAMMPS applies the full pair",
                                       component.name, A, B, pair.scalingVDW));
      }
      if (charged && std::abs(pair.scalingCoulomb - 1.0) > 1e-12) inconsistent = true;
    }
  }
  // 1-4 pairs that are not scaled pairs interact at full strength: non-uniform unless the 1-4 scaling is one
  std::size_t unscaled14 = 0;
  for (const std::array<std::size_t, 2> &pair : pairs14)
  {
    if (intra.exclusions.isExcluded(pair[0], pair[1])) continue;
    const IntraMolecularExclusions::ScaledPair scaling = intra.exclusions.scalingOf(pair[0], pair[1]);
    if (scaling.scalingVDW == 1.0 && scaling.scalingCoulomb == 1.0 && !scaling.pair14) ++unscaled14;
  }
  if (unscaled14 > 0 && (std::abs(vdw14 - 1.0) > 1e-12 || (charged && std::abs(coul14 - 1.0) > 1e-12)))
  {
    result.standard = false;
  }

  // pairs RASPA excludes because they sit inside one rigid fragment but that LAMMPS will compute
  std::size_t rigidExcluded = 0;
  for (const std::array<std::uint32_t, 2> &pair : component.intraMolecularPotentials.exclusions.pairs)
  {
    if (separation[pair[0]][pair[1]] >= 3) ++rigidExcluded;
  }
  if (rigidExcluded > 0)
  {
    warnings.push_back(std::format("{}: {} intramolecular pairs beyond 1-3 are excluded in RASPA (same rigid "
                                   "fragment) but will be computed by LAMMPS{}",
                                   component.name, rigidExcluded,
                                   charged ? " (VDW and Coulomb)" : ""));
  }
  if (inconsistent) result.standard = false;

  result.vdw14 = vdw14;
  result.coul14 = coul14;
  if (!result.standard)
  {
    for (const IntraMolecularExclusions::ScaledPair &pair : intra.exclusions.scaledPairs)
    {
      if (separation[pair.atomA][pair.atomB] == 3) result.nonStandard14.push_back(pair);
    }
  }
  return result;
}

std::string sanitize(std::string_view name)
{
  std::string out{};
  for (char c : name) out += (std::isalnum(static_cast<unsigned char>(c)) ? c : '_');
  if (out.empty() || std::isdigit(static_cast<unsigned char>(out.front()))) out = "c" + out;
  return out;
}
}  // namespace

namespace LAMMPS
{
std::size_t TypeTable::add(const Term &term)
{
  const std::string key = term.key();
  auto it = index.find(key);
  if (it != index.end()) return it->second;
  types.push_back(term);
  index.emplace(key, types.size());
  return types.size();
}

std::vector<std::string> TypeTable::styles() const
{
  std::vector<std::string> result{};
  for (const Term &term : types)
  {
    if (std::find(result.begin(), result.end(), term.style) == result.end()) result.push_back(term.style);
  }
  return result;
}

bool TypeTable::hasStyle(std::string_view style) const
{
  return std::any_of(types.begin(), types.end(), [&](const Term &t) { return t.style == style; });
}

std::vector<std::string> Topology::pairStyles() const
{
  std::vector<std::string> result{};
  for (const PairEntry &pair : pairs)
  {
    if (pair.term.style == "zero") continue;
    if (std::find(result.begin(), result.end(), pair.term.style) == result.end()) result.push_back(pair.term.style);
  }
  return result;
}

bool Topology::pairCoefficientsInDataFile() const
{
  const std::vector<std::string> styles = pairStyles();
  const bool anyZero = std::any_of(pairs.begin(), pairs.end(), [](const PairEntry &p) { return p.term.style == "zero"; });
  return styles.size() == 1 && styles.front() != "table" && !anyZero;
}

Topology buildTopology(std::span<const Component> components, std::span<const Atom> atomData,
                       std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                       const SimulationBox &simulationBox, const ForceField &forceField,
                       std::span<const std::size_t> numberOfIntegerMoleculesPerComponent,
                       const std::optional<Framework> &framework, ExportOptions options)
{
  Topology topology{};
  topology.options = std::move(options);
  TableRegistry tables{};

  const double toAngstrom = Units::LengthConversionFactor * 1e10;
  const double toAngstromPerFemtosecond = Units::VelocityConversionFactor * 1e-5;

  // ---- box --------------------------------------------------------------------------------------------
  {
    const double3 lengths = simulationBox.lengths() * toAngstrom;
    const double3 angles = simulationBox.angles();
    const double xy = lengths.y * std::cos(angles.z);
    const double xz = lengths.z * std::cos(angles.y);
    const double yz = (lengths.y * lengths.z * std::cos(angles.x) - xy * xz) / (lengths.y * std::sin(angles.z));
    topology.triclinic = (std::abs(xy) > 1e-10) || (std::abs(xz) > 1e-10) || (std::abs(yz) > 1e-10);
    topology.lengths = lengths;
    if (topology.triclinic)
    {
      topology.xy = xy;
      topology.xz = xz;
      topology.yz = yz;
      topology.lengths.y = std::sqrt(lengths.y * lengths.y - xy * xy);
      topology.lengths.z = std::sqrt(lengths.z * lengths.z - xz * xz - yz * yz);
    }
  }

  // ---- masses -----------------------------------------------------------------------------------------
  for (const PseudoAtom &pseudoAtom : forceField.pseudoAtoms)
  {
    topology.masses.emplace_back(pseudoAtom.name, pseudoAtom.mass);
  }

  // ---- force-field settings ---------------------------------------------------------------------------
  topology.useCharge = forceField.useCharge;
  topology.mixingRule = forceField.mixingRule;
  topology.chargeMethod = forceField.chargeMethod;
  topology.cutOffVDW = forceField.cutOffMoleculeVDW;
  topology.cutOffFrameworkVDW = forceField.cutOffFrameworkVDW;
  topology.cutOffCoulomb = forceField.cutOffCoulomb;
  topology.ewaldPrecision = forceField.EwaldPrecision;
  topology.alpha = forceField.EwaldAlpha;
  topology.anyTail = std::any_of(forceField.tailCorrections.begin(), forceField.tailCorrections.end(),
                                 [](bool b) { return b; });
  topology.anyShift = std::any_of(forceField.shiftPotentials.begin(), forceField.shiftPotentials.end(),
                                  [](bool b) { return b; });
  topology.allShift = !forceField.shiftPotentials.empty() &&
                      std::all_of(forceField.shiftPotentials.begin(), forceField.shiftPotentials.end(),
                                  [](bool b) { return b; });
  if (topology.anyShift && !topology.allShift)
  {
    topology.warnings.push_back("only some pair potentials are shifted in RASPA; LAMMPS pair_modify shift is global");
  }
  {
    const double3 widths = simulationBox.perpendicularWidths();
    const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
    const double largestCutOff = std::max({topology.cutOffVDW, topology.cutOffFrameworkVDW, topology.cutOffCoulomb});
    if (largestCutOff > halfWidth)
    {
      topology.warnings.push_back(std::format(
          "cut-off {:.3g} A exceeds half the smallest perpendicular box width ({:.3g} A): LAMMPS sums over all "
          "periodic images (including a molecule's own images) whereas RASPA uses the minimum image and unwrapped "
          "intramolecular distances; energies will differ unless the box is enlarged",
          largestCutOff, halfWidth));
    }
  }
  topology.hasFramework = framework.has_value();

  // ---- pair coefficients ------------------------------------------------------------------------------
  for (std::size_t i = 0; i < forceField.numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = i; j < forceField.numberOfPseudoAtoms; ++j)
    {
      PairEntry entry{i + 1, j + 1, pairTerm(forceField, i, j), ""};
      if (entry.term.style == "table")
      {
        entry.tableKeyword = std::format("P{}_{}", i + 1, j + 1);
        topology.tables.emplace_back(
            entry.tableKeyword,
            pairTableBody(forceField, i, j, forceField.cutOffMoleculeVDW, entry.tableKeyword, topology.options.grid));
      }
      topology.pairs.push_back(std::move(entry));
    }
  }

  // ---- atoms: wrapping helper -------------------------------------------------------------------------
  auto makeAtom = [&](const Atom &atom, std::size_t id, std::size_t molecule, double3 velocity)
  {
    AtomEntry entry{};
    entry.id = id;
    entry.molecule = molecule;
    entry.type = static_cast<std::size_t>(atom.type) + 1;
    entry.charge = atom.charge;
    const double3 s = simulationBox.inverseCell * atom.position;
    const int3 image{static_cast<std::int32_t>(std::floor(s.x)), static_cast<std::int32_t>(std::floor(s.y)),
                     static_cast<std::int32_t>(std::floor(s.z))};
    const double3 shift = simulationBox.cell * double3{static_cast<double>(image.x), static_cast<double>(image.y),
                                                       static_cast<double>(image.z)};
    entry.position = (atom.position - shift) * toAngstrom;
    entry.image = image;
    entry.velocity = velocity * toAngstromPerFemtosecond;
    return entry;
  };

  std::size_t atomId = 0;
  std::size_t moleculeId = 0;

  // ---- framework --------------------------------------------------------------------------------------
  const std::size_t numberOfFrameworkAtoms = framework.has_value() ? framework->atoms.size() : 0uz;
  if (framework.has_value() && numberOfFrameworkAtoms > 0)
  {
    ++moleculeId;
    const std::size_t firstId = atomId + 1;
    std::vector<std::size_t> fragmentOfAtom(numberOfFrameworkAtoms, 0);
    std::vector<bool> fixedAtom(numberOfFrameworkAtoms, framework->rigid);
    if (!framework->groups.empty())
    {
      std::fill(fixedAtom.begin(), fixedAtom.end(), false);
      for (const FrameworkGroup &group : framework->groups)
      {
        if (group.type == FrameworkGroupType::Fixed)
        {
          for (std::size_t a : group.atoms) fixedAtom[a] = true;
        }
        else if (group.type == FrameworkGroupType::Rigid && group.atoms.size() >= 2)
        {
          const std::size_t fragment = ++topology.numberOfFragments;
          for (std::size_t a : group.atoms) fragmentOfAtom[a] = fragment;
        }
      }
    }
    for (std::size_t a = 0; a < numberOfFrameworkAtoms && a < atomData.size(); ++a)
    {
      const double3 velocity = (a < atomDynamics.size()) ? atomDynamics[a].velocity : double3{};
      AtomEntry entry = makeAtom(atomData[a], ++atomId, moleculeId, velocity);
      entry.fragment = fragmentOfAtom[a];
      entry.fixed = fixedAtom[a];
      if (entry.fixed) topology.hasFixedAtoms = true;
      topology.atoms.push_back(entry);
    }
    if (!framework->rigid || !framework->intraMolecularPotentials.bonds.empty())
    {
      addBondedTerms(
          topology, tables, framework->intraMolecularPotentials, [&](std::size_t local) { return firstId + local; },
          "framework", true);
    }
    topology.groups.push_back({"framework", moleculeId, moleculeId, false,
                               !framework->intraMolecularPotentials.bonds.empty()});
  }

  // ---- molecules --------------------------------------------------------------------------------------
  // atomData holds the framework first, then the components in order, each molecule contiguous.
  std::size_t cursor = numberOfFrameworkAtoms;
  std::size_t moleculeDataOffset = 0;
  std::vector<ComponentSpecial> specials{};
  std::vector<std::pair<std::size_t, std::size_t>> exportedMolecules{};  // component index, first LAMMPS atom id
  for (std::size_t c = 0; c < components.size(); ++c)
  {
    const Component &component = components[c];
    const std::size_t atomsPerMolecule = component.atoms.size();
    const std::size_t integerMolecules = (c < numberOfIntegerMoleculesPerComponent.size())
                                             ? numberOfIntegerMoleculesPerComponent[c]
                                             : 0uz;

    // rigid-body fragments of the template
    std::vector<std::size_t> fragmentIndexOfAtom(atomsPerMolecule, std::numeric_limits<std::size_t>::max());
    std::size_t fragmentsPerMolecule = 0;
    if (component.rigid && atomsPerMolecule >= 2)
    {
      std::fill(fragmentIndexOfAtom.begin(), fragmentIndexOfAtom.end(), 0uz);
      fragmentsPerMolecule = 1;
    }
    else
    {
      for (const Fragment &fragment : component.fragmentGraph.fragments)
      {
        if (fragment.atoms.size() < 2) continue;
        for (std::size_t a : fragment.atoms)
        {
          if (a < atomsPerMolecule) fragmentIndexOfAtom[a] = fragmentsPerMolecule;
        }
        ++fragmentsPerMolecule;
      }
    }

    const std::size_t firstMolecule = moleculeId + 1;
    std::size_t moleculesSeen = 0;
    while (cursor < atomData.size() && atomData[cursor].componentId == c)
    {
      // one molecule: consecutive atoms with the same moleculeId
      const std::uint32_t localMolecule = atomData[cursor].moleculeId;
      std::size_t end = cursor;
      while (end < atomData.size() && atomData[end].componentId == c && atomData[end].moleculeId == localMolecule)
      {
        ++end;
      }
      ++moleculesSeen;
      if (localMolecule < integerMolecules && (end - cursor) == atomsPerMolecule)
      {
        ++moleculeId;
        const std::size_t firstId = atomId + 1;
        exportedMolecules.emplace_back(c, firstId);
        const std::size_t fragmentBase = topology.numberOfFragments;
        topology.numberOfFragments += fragmentsPerMolecule;
        const std::size_t moleculeIndex = moleculeDataOffset + localMolecule;
        for (std::size_t a = cursor; a < end; ++a)
        {
          const std::size_t local = a - cursor;
          double3 velocity{};
          if (component.rigid && moleculeIndex < moleculeData.size())
          {
            velocity = moleculeData[moleculeIndex].velocity;
          }
          else if (a < atomDynamics.size())
          {
            velocity = atomDynamics[a].velocity;
          }
          AtomEntry entry = makeAtom(atomData[a], ++atomId, moleculeId, velocity);
          if (fragmentIndexOfAtom[local] != std::numeric_limits<std::size_t>::max())
          {
            entry.fragment = fragmentBase + fragmentIndexOfAtom[local] + 1;
          }
          topology.atoms.push_back(entry);
        }
        if (!component.rigid)
        {
          addBondedTerms(
              topology, tables, component.intraMolecularPotentials,
              [&](std::size_t local) { return firstId + local; }, component.name, localMolecule == 0);
        }
      }
      else if (localMolecule >= integerMolecules)
      {
        // fractional (CFCMC) molecule: not part of the exported system
      }
      else
      {
        topology.warnings.push_back(std::format("{}: molecule {} has {} atoms instead of {}; skipped",
                                                component.name, localMolecule, end - cursor, atomsPerMolecule));
      }
      cursor = end;
    }
    moleculeDataOffset += moleculesSeen;
    if (moleculeId >= firstMolecule)
    {
      topology.groups.push_back({sanitize(component.name), firstMolecule, moleculeId, component.rigid,
                                 !component.intraMolecularPotentials.bonds.empty()});
    }
    specials.push_back(analyseSpecial(component, topology.warnings));
  }
  if (cursor < atomData.size())
  {
    topology.warnings.push_back(std::format("{} atoms after the last component were not exported",
                                            atomData.size() - cursor));
  }

  // ---- special_bonds ----------------------------------------------------------------------------------
  {
    std::optional<double> vdw14{}, coul14{};
    bool uniform = true;
    for (const ComponentSpecial &s : specials)
    {
      if (!s.present) continue;
      topology.special.present = true;
      if (!s.standard) uniform = false;
      if (!vdw14.has_value())
      {
        vdw14 = s.vdw14;
        coul14 = s.coul14;
      }
      else if (std::abs(*vdw14 - s.vdw14) > 1e-12 || std::abs(*coul14 - s.coul14) > 1e-12)
      {
        uniform = false;
        topology.warnings.push_back("components use different 1-4 scalings; special_bonds is global");
      }
    }
    topology.special.uniform = uniform;
    topology.special.vdw14 = vdw14.value_or(0.0);
    topology.special.coul14 = coul14.value_or(0.0);

    if (!uniform)
    {
      // 1-4 VDW pairs go into a 'pair_style list' file with their own scaling; special_bonds lj becomes 0
      topology.special.vdw14 = 0.0;
      std::set<std::string> reported{};
      for (const auto &[c, firstId] : exportedMolecules)
      {
        const Component &component = components[c];
        for (const IntraMolecularExclusions::ScaledPair &pair : specials[c].nonStandard14)
        {
          const std::size_t typeA = component.atoms[pair.atomA].type;
          const std::size_t typeB = component.atoms[pair.atomB].type;
          const VDWParameters &vdw = forceField.pair(typeA, typeB, pair.pair14);
          // (the potential-switched form is plain Lennard-Jones at the 1-4 distances)
          if (vdw.type != VDWParameters::Type::LennardJones && vdw.type != VDWParameters::Type::LennardJonesSwitched)
          {
            if (reported.insert(component.name).second)
            {
              topology.warnings.push_back(std::format(
                  "{}: scaled 1-4 pair of non-LJ type cannot go into pair_style list; dropped", component.name));
            }
            continue;
          }
          topology.pairList.push_back(std::format(
              "{} {} lj126 {} {} {}", firstId + pair.atomA, firstId + pair.atomB,
              formatValue(pair.scalingVDW * vdw.parameters.x * Units::EnergyToKCalPerMol),
              formatValue(vdw.parameters.y), formatValue(forceField.cutOffMoleculeVDW)));
        }
      }
    }
  }

  return topology;
}
}  // namespace LAMMPS
