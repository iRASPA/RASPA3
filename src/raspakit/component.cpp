module;

#if !defined(_WIN32)
#include <assert.h>
#endif

module component;

import std;

import int3;
import simd_quatd;
import double3;
import double3x3;
import randomnumbers;
import archive;
import json;
import units;
import skposcarparser;
import characterset;
import stringutils;
import skparser;
import skstructure;
import skatom;
import skcell;
import skspacegroup;
import forcefield;
import atom;
import property_lambda_probability_histogram;
import property_widom;
import simulationbox;
import cif_reader;
import move_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import mc_moves_cputime;
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
import van_der_waals_potential;
import coulomb_potential;
import intra_molecular_potentials;
import chiral_center;
import cbmc_growth_plan;
import cbmc_flexible_base;
import cbmc_bridge_closure;
import vdwparameters;
import blocking_pockets;

namespace
{
std::vector<std::array<std::size_t, 2>> readExplicitIntraPairs(
    const nlohmann::basic_json<nlohmann::raspa_map>& parsedData, std::string_view key,
    std::size_t numberOfAtoms)
{
  const auto& value = parsedData.at(std::string{key});
  if (!value.is_array())
  {
    throw std::runtime_error(std::format("[Component reader]: {} must be an array of atom-index pairs\n", key));
  }

  std::vector<std::array<std::size_t, 2>> pairs{};
  pairs.reserve(value.size());
  for (const auto& item : value)
  {
    if (!item.is_array() || item.size() != 2 || !item[0].is_number_unsigned() ||
        !item[1].is_number_unsigned())
    {
      throw std::runtime_error(
          std::format("[Component reader]: {} item {} must contain two unsigned atom indices\n", key, item.dump()));
    }

    std::array<std::size_t, 2> pair{item[0].get<std::size_t>(), item[1].get<std::size_t>()};
    if (pair[0] >= numberOfAtoms || pair[1] >= numberOfAtoms || pair[0] == pair[1])
    {
      throw std::runtime_error(
          std::format("[Component reader]: {} item {} contains an invalid atom index\n", key, item.dump()));
    }
    if (pair[1] < pair[0]) std::swap(pair[0], pair[1]);
    if (std::ranges::find(pairs, pair) != pairs.end())
    {
      throw std::runtime_error(
          std::format("[Component reader]: {} contains duplicate pair [{}, {}]\n", key, pair[0], pair[1]));
    }
    pairs.push_back(pair);
  }
  return pairs;
}
}  // namespace

// default constructor, needed for binary restart-file
Component::Component() {}

// create Component in 'inputreader.cpp'
Component::Component(Component::Type type, std::size_t componentId, const ForceField &forceField,
                     const std::string &componentName, std::optional<const std::string> fileName,
                     std::size_t numberOfBlocks, std::size_t numberOfLambdaBins,
                     const MCMoveProbabilities &particleProbabilities, std::optional<double> fugacityCoefficient,
                     bool thermodynamicIntegration) noexcept(false)
    : type(type),
      name(componentName),
      filenameData(fileName),
      fugacityCoefficient(fugacityCoefficient),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      lambdaPairSwap(numberOfBlocks, numberOfLambdaBins),
      lambdaPairSwapCB(numberOfBlocks, numberOfLambdaBins),
      lambdaGroupSwap(numberOfBlocks, numberOfLambdaBins),
      lambdaGroupSwapCB(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbabilities),
      averageRosenbluthWeights(numberOfBlocks),
      averageGibbsRosenbluthWeights(numberOfBlocks)
{
  if (filenameData.has_value())
  {
    readComponent(componentId, forceField, filenameData.value());
  }
  lambdaGC.computeDUdlambda = thermodynamicIntegration;

  cbmcMoveStatistics = std::vector<CBMC::InternalMoveStatistics>(definedAtoms.size());
}

// create programmatically an 'adsorbate' component
Component::Component(const ForceField &forceField, std::string componentName, double T_c,
                     double P_c, double w, std::vector<Atom> atomList, const ConnectivityTable &connectivityTable,
                     const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::size_t numberOfBlocks,
                     std::size_t numberOfLambdaBins, const MCMoveProbabilities &particleProbabilities,
                     std::optional<double> fugacityCoefficient, bool thermodynamicIntegration,
                     std::vector<double4> blockingPockets) noexcept(false)
    : type(Type::Adsorbate),
      name(componentName),
      criticalTemperature(T_c),
      criticalPressure(P_c),
      acentricFactor(w),
      fugacityCoefficient(fugacityCoefficient),
      connectivityTable(connectivityTable),
      intraMolecularPotentials(intraMolecularPotentials),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      lambdaPairSwap(numberOfBlocks, numberOfLambdaBins),
      lambdaPairSwapCB(numberOfBlocks, numberOfLambdaBins),
      lambdaGroupSwap(numberOfBlocks, numberOfLambdaBins),
      lambdaGroupSwapCB(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbabilities),
      cbmcMoveStatistics(atomList.size()),
      averageRosenbluthWeights(numberOfBlocks),
      averageGibbsRosenbluthWeights(numberOfBlocks),
      blockingPockets(blockingPockets)
{
  totalMass = 0.0;
  netCharge = 0.0;
  for (const Atom &atom : atomList)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    double mass = forceField.pseudoAtoms[atomType].mass;
    totalMass += mass;
    netCharge += atom.charge;
    definedAtoms.push_back({atom, mass});
  }

  // A default-constructed connectivity table (no bonds) must still cover all atoms so the fragment
  // graph can index them.
  if (this->connectivityTable.numberOfBeads != definedAtoms.size())
  {
    this->connectivityTable = ConnectivityTable(definedAtoms.size());
  }

  // A programmatic component is flexible when it carries bonded interactions; otherwise it is a
  // single rigid body covering all atoms. No rigid sub-bodies are defined through this path.
  rigid = intraMolecularPotentials.bonds.empty();

  computeRigidProperties();
  lambdaGC.computeDUdlambda = thermodynamicIntegration;

  std::vector<std::vector<std::size_t>> rigidBodies{};
  if (rigid)
  {
    std::vector<std::size_t> wholeMolecule(definedAtoms.size());
    std::iota(wholeMolecule.begin(), wholeMolecule.end(), std::size_t{0});
    rigidBodies.push_back(std::move(wholeMolecule));
  }
  buildFragmentGraph(rigidBodies);
}

// read the component from the molecule-file
void Component::readComponent(std::size_t componentId, const ForceField &forceField, const std::string &fileName)
{
  const std::string defaultMoleculeFileName = addExtension(fileName, ".json");

  std::string moleculeFileName = defaultMoleculeFileName;
  if (!std::filesystem::exists(std::filesystem::path{moleculeFileName}))
  {
    const char *env_p = std::getenv("RASPA_DIR");
    if (env_p)
    {
      moleculeFileName = env_p + std::string("/") + defaultMoleculeFileName;
    }
  }

  if (!std::filesystem::exists(moleculeFileName))
  {
    throw std::runtime_error(std::format("[Component reader]: File '{}' not found\n", moleculeFileName));
  }

  std::filesystem::path moleculePathfile = std::filesystem::path(moleculeFileName);
  std::ifstream moleculeStream{moleculePathfile};
  if (!moleculeStream)
  {
    throw std::runtime_error(
        std::format("[Component reader] File '{}' exists, but error opening file\n", moleculeFileName));
  }

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  try
  {
    parsed_data = nlohmann::json::parse(moleculeStream);
  }
  catch (nlohmann::json::parse_error &ex)
  {
    throw std::runtime_error(std::format("[Component reader]: Parse error of file {} at byte {}\n{}\n",
                                         std::format("{}.json", fileName), ex.byte, ex.what()));
  }

  try
  {
    criticalTemperature = parsed_data.value("CriticalTemperature", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'CriticalTemperature' listed as {} must be floating point number\n{}\n",
                    parsed_data["CriticalTemperature"].dump(), ex.what()));
  }

  try
  {
    criticalPressure = parsed_data.value("CriticalPressure", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'CriticalPressure' listed as {} must be floating point number\n{}\n",
                    parsed_data["CriticalPressure"].dump(), ex.what()));
  }

  try
  {
    acentricFactor = parsed_data.value("AcentricFactor", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'AcentricFactor' listed as {} must be floating point number\n{}\n",
                    parsed_data["AcentricFactor"].dump(), ex.what()));
  }

  if (parsed_data.contains("BlockingPockets"))
  {
    BlockingPockets::Specification specification = BlockingPockets::parse(parsed_data["BlockingPockets"]);
    blockingPockets.insert(blockingPockets.end(), specification.pockets.begin(), specification.pockets.end());
    automaticBlockingPockets = automaticBlockingPockets || specification.automatic;
  }

  std::size_t jsonNumberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();

  definedAtoms.clear();
  definedAtoms.reserve(jsonNumberOfPseudoAtoms);

  if (jsonNumberOfPseudoAtoms == 0)
  {
    throw std::runtime_error(std::format("[Component reader]: key 'PseudoAtoms' empty]\n"));
  }

  if (!parsed_data.contains("PseudoAtoms"))
  {
    throw std::runtime_error(
        std::format("[Component reader]: No pseudo-atoms found [keyword 'PseudoAtoms' missing]\n"));
  }

  for (auto &[_, item] : parsed_data["PseudoAtoms"].items())
  {
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    if (item.size() != 2)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with two elements, "
                      "the pseudo-atom-name and an array with the x,y,z positions\n",
                      item.dump()));
    }

    if (!item[0].is_string())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an string (the name of the pseudo-atom)\n", item[0].dump()));
    }
    std::string pseudoAtomName = item[0].get<std::string>();

    // find atom-type based on read 'atomTypeString'
    std::optional<std::size_t> index = forceField.findPseudoAtom(pseudoAtomName);
    if (!index.has_value())
    {
      throw std::runtime_error(
          std::format("[Component reader]: unknown pseudo-atom '{}', please lookup type in in 'force_field.json'\n",
                      pseudoAtomName));
    }
    std::size_t pseudoAtomType = index.value();

    if (!item[1].is_array())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array (with the positions)\n", item[1].dump()));
    }

    if (item[1].size() != 3)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with three elements, "
                      "the x,y,z positions\n",
                      item[1].dump()));
    }

    std::vector<double> position{};
    try
    {
      position = item[1].is_array() ? item[1].get<std::vector<double>>() : std::vector<double>{};
    }
    catch (nlohmann::json::exception &ex)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be array of three floating point numbers \n{}\n",
                      item[1].dump(), ex.what()));
    }

    double mass = forceField.pseudoAtoms[pseudoAtomType].mass;
    double charge = forceField.pseudoAtoms[pseudoAtomType].charge;
    double scaling = 1.0;

    definedAtoms.push_back(
        {Atom(double3(position[0], position[1], position[2]), charge, scaling, 0,
              static_cast<std::uint16_t>(pseudoAtomType), static_cast<std::uint8_t>(componentId), 0, 0),
         mass});
  }

  if (parsed_data.contains("StartingBead"))
  {
    if (!parsed_data["StartingBead"].is_number_integer())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an integer\n", parsed_data["StartingBead"].dump()));
    }

    std::int64_t starting_bead = parsed_data["StartingBead"].get<std::int64_t>();
    if (starting_bead >= 0)
    {
      startingBead = static_cast<std::size_t>(starting_bead);
    }
  }

  totalMass = 0.0;
  netCharge = 0.0;
  for (auto [atom, mass] : definedAtoms)
  {
    totalMass += mass;
    netCharge += atom.charge;
  }

  connectivityTable = readConnectivityTable(definedAtoms.size(), parsed_data);

  // New molecule schema: a molecule is described by its 'Connectivity' and an optional list of
  // 'RigidBodies' (each a list of atom indices that move as one rigid unit). Every atom not listed in
  // a rigid body is a single-atom (flexible) fragment. A molecule with no bonds and no explicitly
  // listed rigid body is a single rigid body covering all atoms (the common rigid-adsorbate case).
  std::vector<std::vector<std::size_t>> rigidBodies = readRigidBodies(parsed_data);
  bool hasBonds = !connectivityTable.findAllBonds().empty();
  if (rigidBodies.empty() && !hasBonds)
  {
    std::vector<std::size_t> wholeMolecule(definedAtoms.size());
    std::iota(wholeMolecule.begin(), wholeMolecule.end(), std::size_t{0});
    rigidBodies.push_back(std::move(wholeMolecule));
  }

  std::vector<std::vector<std::size_t>> partition =
      FragmentGraph::partitionAtoms(definedAtoms.size(), rigidBodies);

  // Fully rigid molecule when the whole molecule collapses into a single fragment; otherwise the
  // molecule is flexible or semi-flexible and its internal potentials are sampled during growth.
  rigid = partition.size() == 1;

  computeRigidProperties();

  if (!rigid)
  {
    // Read the intramolecular potentials before building the fragment graph filter is applied: the
    // bond/bend/torsion readers use 'fragmentGraph' to skip interactions that lie entirely inside a
    // single rigid fragment (their geometry is fixed), so build the graph first.
    buildFragmentGraph(rigidBodies);

    intraMolecularPotentials.bonds = readBondPotentials(forceField, parsed_data);
    intraMolecularPotentials.bends = readBendPotentials(forceField, parsed_data);
    intraMolecularPotentials.torsions = readTorsionPotentials(forceField, parsed_data);
    intraMolecularPotentials.ureyBradleys = readUreyBradleyPotentials(parsed_data);
    intraMolecularPotentials.inversionBends = readInversionBendPotentials(parsed_data);
    intraMolecularPotentials.outOfPlaneBends = readOutOfPlaneBendPotentials(parsed_data);
    intraMolecularPotentials.improperTorsions = readImproperTorsionPotentials(parsed_data);
    intraMolecularPotentials.bondBonds = readBondBondPotentials(parsed_data);
    intraMolecularPotentials.bondBends = readBondBendPotentials(parsed_data);
    intraMolecularPotentials.bondTorsions = readBondTorsionPotentials(parsed_data);
    intraMolecularPotentials.bendBends = readBendBendPotentials(parsed_data);
    intraMolecularPotentials.bendTorsions = readBendTorsionPotentials(parsed_data);
    intraMolecularPotentials.vanDerWaals = readVanDerWaalsPotentials(forceField, parsed_data);
    intraMolecularPotentials.coulombs = readCoulombPotentials(forceField, parsed_data);
    intraMolecularPotentials.chiralCenters = readChiralCenters(parsed_data);

    partialReinsertionFixedAtoms = readPartialReinsertionFixedAtoms(parsed_data);
    repeatUnits = readRepeatUnits(parsed_data);
    endToEndAtoms = determineEndToEndAtoms(parsed_data);

    // Whether an atom lies in a cyclic cluster (a flexible ring): its fragment maps to a cluster id.
    auto cyclicClusterOf = [&](std::size_t atom) -> std::optional<std::size_t>
    { return fragmentGraph.fragmentCyclicClusterIds[fragmentGraph.atomFragmentIds[atom]]; };

    // Parse-time validation: inside a ring, the reference coordinates of a Fixed bond are
    // load-bearing. Acyclic growth regenerates a fixed bond at its exact length every placement, but
    // ring-closure growth takes the ring geometry from the reference conformation and afterwards only
    // applies constraint-preserving rotations -- which preserve whatever length the reference has.
    // Since a Fixed bond carries no energy, a violated reference length is never corrected and every
    // sampled configuration silently inherits it, so reject it here instead.
    for (const BondPotential &bond : intraMolecularPotentials.bonds)
    {
      if (bond.type != BondType::Fixed) continue;
      std::size_t atomA = bond.identifiers[0];
      std::size_t atomB = bond.identifiers[1];
      // Only bonds inside one cyclic cluster matter; two nullopts compare equal, so exclude those first.
      if (!cyclicClusterOf(atomA).has_value() || cyclicClusterOf(atomA) != cyclicClusterOf(atomB)) continue;
      double referenceLength = (atoms[atomA].position - atoms[atomB].position).length();
      if (std::abs(referenceLength - bond.parameters[0]) > 1e-3)
      {
        throw std::runtime_error(std::format(
            "[Component reader]: component '{}': FIXED ring bond {}-{} has length {:.6f} Å in the reference "
            "coordinates but declares a fixed length of {:g} Å. Ring geometry is taken from the reference "
            "coordinates and a FIXED bond is never re-imposed, so the reference must satisfy it; adjust the "
            "'pseudoAtoms' coordinates.\n",
            name, atomA, atomB, referenceLength, bond.parameters[0]));
      }
    }

    // Parse-time diagnostic: a flexible ring atom pinned by three or more Fixed bonds. Fixed bonds
    // are holonomic constraints carrying no energy, so the ring-closure Monte Carlo may only move
    // such an atom by rotations about its fixed neighbours; with one or two fixed bonds a sphere or
    // circle of freedom remains, but with three or more no single-atom rotation exists. The
    // constraints stay satisfied, but the atom is frozen at its as-grown position relative to the
    // ring and its conformational degree of freedom is never sampled (this typically occurs at
    // fused-ring junctions where every ring bond is declared FIXED).
    std::vector<std::size_t> fixedBondCounts(definedAtoms.size(), 0uz);
    for (const BondPotential &bond : intraMolecularPotentials.bonds)
    {
      if (bond.type != BondType::Fixed) continue;
      ++fixedBondCounts[bond.identifiers[0]];
      ++fixedBondCounts[bond.identifiers[1]];
    }
    for (const std::vector<std::size_t> &cluster : fragmentGraph.cyclicClusters)
    {
      for (std::size_t atom : cluster)
      {
        if (fragmentGraph.fragments[fragmentGraph.atomFragmentIds[atom]].isRigidBody()) continue;
        if (fixedBondCounts[atom] < 3) continue;
        std::print(std::cerr,
                   "[Component reader]: warning: component '{}': flexible ring atom {} has {} FIXED bonds; its "
                   "position is fully constrained, so ring Monte-Carlo cannot move it and its conformation is not "
                   "sampled. Consider stiff HARMONIC bonds or a rigid body instead.\n",
                   name, atom, fixedBondCounts[atom]);
      }
    }
  }
  else
  {
    buildFragmentGraph(rigidBodies);
  }

  // Warm the growth-plan cache for the placed sets used during the simulation: full growth from the
  // starting bead and every partial-reinsertion fixed set.
  if (definedAtoms.size() > 1)
  {
    growthPlan({startingBead});
    for (const std::vector<std::size_t> &fixedAtoms : partialReinsertionFixedAtoms)
    {
      growthPlan(fixedAtoms);
    }

    // Parse-time diagnostic: growth steps that can only be grown with the engine's default geometry
    // (a bonded pair without a bond potential, a rigid junction without its bend). Reported once per
    // distinct message, since the partial plans repeat the steps of the full plan.
    std::set<std::string> reported{};
    for (const auto &[placed, plan] : growthPlanCache)
    {
      for (const std::string &warning : CBMC::growthPlanDefaultGeometryWarnings(plan))
      {
        if (!reported.insert(warning).second) continue;
        std::print(std::cerr, "[Component reader]: warning: component '{}': {}\n", name, warning);
      }
    }
  }
}

std::vector<std::vector<std::size_t>> Component::readRigidBodies(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<std::vector<std::size_t>> rigidBodies{};

  if (!parsed_data.contains("RigidBodies"))
  {
    return rigidBodies;
  }

  if (!parsed_data["RigidBodies"].is_array())
  {
    throw std::runtime_error(std::format("[Component reader]: 'RigidBodies' must be an array of atom-index lists\n"));
  }

  std::size_t numberOfBeads = definedAtoms.size();
  for (auto &[_, item] : parsed_data["RigidBodies"].items())
  {
    if (!item.is_array())
    {
      throw std::runtime_error(
          std::format("[Component reader]: each entry of 'RigidBodies' must be an array of atom indices\n"));
    }

    std::vector<std::size_t> bodyAtoms = item.get<std::vector<std::size_t>>();
    if (bodyAtoms.size() < 2)
    {
      throw std::runtime_error(
          std::format("[Component reader]: a rigid body must contain at least two atoms\n"));
    }

    for (std::size_t atom : bodyAtoms)
    {
      if (atom >= numberOfBeads)
      {
        throw std::runtime_error(std::format(
            "[Component reader]: rigid-body atom index {} out of range (molecule has {} atoms)\n", atom,
            numberOfBeads));
      }
    }

    // A rigid body must be a connected subgraph of the connectivity table (its atoms move together
    // as one body, so they must be bonded into one piece).
    if (!connectivityTable.checkIsConnectedSubgraph(bodyAtoms))
    {
      std::stringstream result{};
      std::copy(bodyAtoms.begin(), bodyAtoms.end(), std::ostream_iterator<std::size_t>(result, " "));
      throw std::runtime_error(
          std::format("[Component reader]: rigid body ({}) is not a connected subgraph\n", result.str()));
    }

    rigidBodies.push_back(std::move(bodyAtoms));
  }

  return rigidBodies;
}

void Component::buildFragmentGraph(const std::vector<std::vector<std::size_t>> &rigidBodies)
{
  // The topology changes, so any cached growth plans and the data derived from them are stale.
  growthPlanCache.clear();
  baseCouplingConstantsMemo.clear();
  closureGuideMemo.clear();
  recoilReferenceStepEnergiesCache.clear();

  std::size_t numberOfBeads = definedAtoms.size();

  std::vector<std::vector<std::size_t>> partition = FragmentGraph::partitionAtoms(numberOfBeads, rigidBodies);

  // Reference geometry and per-atom masses used to compute the rigid-body body frames. 'atoms' is the
  // reference geometry (already in the whole-molecule principal frame) the CBMC growth places
  // fragments from, so the body-fixed positions match the offsets it produces.
  std::vector<double3> referencePositions(numberOfBeads);
  std::vector<double> masses(numberOfBeads);
  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    referencePositions[i] = atoms[i].position;
    masses[i] = definedAtoms[i].second;
  }

  fragmentGraph.build(connectivityTable, partition, startingBead, referencePositions, masses);

  // A semi-flexible molecule is integrated as a set of rigid-body fragments (each contributing three
  // center-of-mass translations plus its rotational degrees of freedom) together with per-atom
  // Cartesian degrees of freedom for every single-atom (flexible) fragment.
  if (fragmentGraph.isSemiFlexible())
  {
    std::size_t rigidFragmentAtoms{};
    std::size_t rotationalDof{};
    for (const Fragment &fragment : fragmentGraph.fragments)
    {
      if (!fragment.isRigidBody()) continue;
      rigidFragmentAtoms += fragment.atoms.size();
      rotationalDof += fragment.rotationalDegreesOfFreedom;
    }
    std::size_t flexibleAtoms = numberOfBeads - rigidFragmentAtoms;
    translationalDegreesOfFreedom = 3 * flexibleAtoms + 3 * fragmentGraph.numberOfRigidFragments();
    rotationalDegreesOfFreedom = rotationalDof;
  }
}

const std::vector<CBMC::GrowStep> &Component::growthPlan(const std::vector<std::size_t> &beadsAlreadyPlaced) const
{
  auto it = growthPlanCache.find(beadsAlreadyPlaced);
  if (it == growthPlanCache.end())
  {
    it = growthPlanCache
             .emplace(beadsAlreadyPlaced, CBMC::buildGrowthPlan(connectivityTable, fragmentGraph,
                                                                intraMolecularPotentials, beadsAlreadyPlaced))
             .first;
    if (growthPlanBeta.has_value())
    {
      CBMC::prepareBaseCouplingConstants(growthPlanBeta.value(), atoms.size(), it->second, baseCouplingConstantsMemo);
      CBMC::prepareClosureGuides(growthPlanBeta.value(), it->second, closureGuideMemo);
    }
  }
  return it->second;
}

void Component::prepareGrowthPlans(double beta) const
{
  if (growthPlanBeta.has_value() && growthPlanBeta.value() == beta) return;

  growthPlanBeta = beta;
  baseCouplingConstantsMemo.clear();
  closureGuideMemo.clear();
  for (auto &[beadsAlreadyPlaced, plan] : growthPlanCache)
  {
    CBMC::prepareBaseCouplingConstants(beta, atoms.size(), plan, baseCouplingConstantsMemo);
    CBMC::prepareClosureGuides(beta, plan, closureGuideMemo);
  }
}

void Component::setRecoilReferenceConformations(std::vector<std::vector<Atom>> conformations)
{
  recoilReferenceConformations = std::move(conformations);
  recoilReferenceStepEnergiesCache.clear();
}

const std::vector<double> &Component::recoilReferenceStepEnergies(
    const std::vector<std::size_t> &beadsAlreadyPlaced) const
{
  auto it = recoilReferenceStepEnergiesCache.find(beadsAlreadyPlaced);
  if (it != recoilReferenceStepEnergiesCache.end()) return it->second;

  const std::vector<CBMC::GrowStep> &plan = growthPlan(beadsAlreadyPlaced);
  std::vector<double> reference(plan.size());
  for (std::size_t seg = 0; seg != plan.size(); ++seg)
  {
    const Potentials::IntraMolecularPotentials &intra = plan[seg].intra;
    double referenceEnergy = 0.0;
    if (recoilReferenceConformations.empty())
    {
      // No reference conformations built: fall back to the component's declared geometry.
      referenceEnergy = intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(atoms).potentialEnergy();
    }
    else
    {
      // The MAXIMUM over the equilibrated conformations: any strain a valid chain exhibits at this
      // step is acceptable there (see the recoil-growth openness test for the rationale).
      for (const std::vector<Atom> &conformation : recoilReferenceConformations)
      {
        referenceEnergy = std::max(
            referenceEnergy, intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(conformation).potentialEnergy());
      }
    }
    reference[seg] = std::max(0.0, referenceEnergy);
  }
  return recoilReferenceStepEnergiesCache.emplace(beadsAlreadyPlaced, std::move(reference)).first->second;
}

const std::vector<Component::PivotBond> &Component::pivotBonds() const
{
  if (pivotBondsCache.has_value()) return pivotBondsCache.value();

  std::vector<PivotBond> result;
  const std::size_t numberOfBeads = connectivityTable.numberOfBeads;

  // Collects the atoms reachable from 'root' without crossing the bond (root, excluded). The two
  // axis atoms are not part of the result (they lie on the rotation axis and do not move). Returns
  // nullopt when 'excluded' is reached through another path, i.e. the bond is part of a ring and
  // there is no well-defined part to rotate.
  auto partition = [&](std::size_t root, std::size_t excluded) -> std::optional<std::vector<std::size_t>>
  {
    std::vector<bool> visited(numberOfBeads, false);
    visited[root] = true;
    std::vector<std::size_t> stack{root};
    std::vector<std::size_t> collected;
    while (!stack.empty())
    {
      std::size_t current = stack.back();
      stack.pop_back();
      for (std::size_t neighbor : connectivityTable.findAllNeighbors(current))
      {
        if (neighbor == excluded)
        {
          if (current == root) continue;  // the pivot bond itself
          return std::nullopt;            // reached the other axis atom through a ring
        }
        if (!visited[neighbor])
        {
          visited[neighbor] = true;
          stack.push_back(neighbor);
          collected.push_back(neighbor);
        }
      }
    }
    return collected;
  };

  for (const std::array<std::size_t, 2> &bond : connectivityTable.findAllBonds())
  {
    // A bond between two atoms of the same rigid fragment cannot be a pivot axis: rotating about it
    // would change the internal geometry of the rigid body. Rigid fragments that lie entirely on
    // one side of a pivot bond are fine (they rotate rigidly, which preserves their geometry).
    const std::vector<std::size_t> &atomFragmentIds = fragmentGraph.atomFragmentIds;
    if (!atomFragmentIds.empty() && atomFragmentIds[bond[0]] == atomFragmentIds[bond[1]] &&
        fragmentGraph.fragments[atomFragmentIds[bond[0]]].isRigidBody())
    {
      continue;
    }

    std::optional<std::vector<std::size_t>> sideOfSecond = partition(bond[1], bond[0]);
    if (!sideOfSecond.has_value()) continue;  // ring bond

    // Rotate the smaller of the two parts; the choice is deterministic given the bond, so the
    // pivot-move proposal stays symmetric.
    std::size_t sizeOfFirstSide = numberOfBeads - 2 - sideOfSecond->size();
    std::vector<std::size_t> rotatedAtoms;
    if (sideOfSecond->size() <= sizeOfFirstSide)
    {
      rotatedAtoms = std::move(sideOfSecond.value());
    }
    else
    {
      rotatedAtoms = std::move(partition(bond[0], bond[1]).value());
    }
    if (rotatedAtoms.empty()) continue;  // terminal bond, nothing to rotate

    result.push_back(PivotBond{bond, std::move(rotatedAtoms)});
  }

  pivotBondsCache = std::move(result);
  return pivotBondsCache.value();
}

const std::vector<Component::CrankshaftUnit> &Component::crankshaftUnits() const
{
  if (crankshaftUnitsCache.has_value()) return crankshaftUnitsCache.value();

  std::vector<CrankshaftUnit> result;
  const std::size_t numberOfBeads = connectivityTable.numberOfBeads;
  const std::vector<std::size_t> &atomFragmentIds = fragmentGraph.atomFragmentIds;

  // Enumerate all anchor pairs (a, b). Removing both anchors from the molecular graph splits it
  // into connected components; every component that is bonded to both anchors lies 'between' them
  // and is a candidate segment for rotation about the a-b axis. Since the anchors are fixed and
  // every path from the segment to the rest of the molecule passes through an anchor, a rigid
  // rotation of the segment preserves all bond lengths (distances to on-axis points are invariant).
  // This construction handles branched molecules (side chains attached to the segment are part of
  // its component and rotate along) and flexible rings (both arcs between two ring atoms are valid
  // segments), which the pivot move cannot sample.
  for (std::size_t anchorA = 0; anchorA + 1 < numberOfBeads; ++anchorA)
  {
    for (std::size_t anchorB = anchorA + 1; anchorB < numberOfBeads; ++anchorB)
    {
      std::vector<bool> visited(numberOfBeads, false);
      visited[anchorA] = true;
      visited[anchorB] = true;

      for (std::size_t seed = 0; seed < numberOfBeads; ++seed)
      {
        if (visited[seed]) continue;

        // Collect the connected component containing 'seed' in the graph without the anchors, and
        // record whether it is bonded to each anchor.
        std::vector<std::size_t> segment{seed};
        std::vector<std::size_t> stack{seed};
        visited[seed] = true;
        bool bondedToA = false;
        bool bondedToB = false;
        while (!stack.empty())
        {
          std::size_t current = stack.back();
          stack.pop_back();
          for (std::size_t neighbor : connectivityTable.findAllNeighbors(current))
          {
            if (neighbor == anchorA)
            {
              bondedToA = true;
            }
            else if (neighbor == anchorB)
            {
              bondedToB = true;
            }
            else if (!visited[neighbor])
            {
              visited[neighbor] = true;
              stack.push_back(neighbor);
              segment.push_back(neighbor);
            }
          }
        }

        if (!bondedToA || !bondedToB) continue;
        if (segment.size() > crankshaftMaxSegmentSize) continue;

        // Every rigid fragment must lie entirely inside or entirely outside the rotated segment;
        // a fragment straddling the boundary (including one containing an anchor) would deform.
        std::vector<bool> inSegment(numberOfBeads, false);
        for (std::size_t atom : segment) inSegment[atom] = true;
        bool deformsRigidFragment = false;
        for (std::size_t atom : segment)
        {
          if (atomFragmentIds.empty()) break;
          const auto &fragment = fragmentGraph.fragments[atomFragmentIds[atom]];
          if (!fragment.isRigidBody()) continue;
          for (std::size_t fragmentAtom : fragment.atoms)
          {
            if (!inSegment[fragmentAtom])
            {
              deformsRigidFragment = true;
              break;
            }
          }
          if (deformsRigidFragment) break;
        }
        if (deformsRigidFragment) continue;

        std::sort(segment.begin(), segment.end());
        result.push_back(CrankshaftUnit{{anchorA, anchorB}, std::move(segment)});
      }
    }
  }

  crankshaftUnitsCache = std::move(result);
  return crankshaftUnitsCache.value();
}

bool Component::isInsideRigidFragment(std::span<const std::size_t> ids) const
{
  return fragmentGraph.isInsideRigidFragment(ids);
}

std::optional<std::size_t> Component::rigidFragmentContaining(std::size_t bead) const
{
  if (fragmentGraph.atomFragmentIds.empty() || bead >= fragmentGraph.atomFragmentIds.size())
  {
    return std::nullopt;
  }
  std::size_t f = fragmentGraph.atomFragmentIds[bead];
  if (f < fragmentGraph.fragments.size() && fragmentGraph.fragments[f].isRigidBody())
  {
    return f;
  }
  return std::nullopt;
}

bool Component::isSemiFlexible() const { return fragmentGraph.isSemiFlexible(); }

std::size_t Component::numberOfRigidFragments() const { return fragmentGraph.numberOfRigidFragments(); }

void Component::regenerateFragmentAtoms(const GroupState &state, std::size_t fragmentIndex,
                                        std::span<Atom> moleculeAtoms) const
{
  fragmentGraph.fragments[fragmentIndex].regenerateAtoms(state, moleculeAtoms);
}

GroupState Component::deriveFragmentState(std::size_t fragmentIndex, std::span<const Atom> moleculeAtoms) const
{
  return fragmentGraph.fragments[fragmentIndex].deriveState(moleculeAtoms);
}

void Component::computeRigidProperties()
{
  double3 com{};
  double total_mass{};
  for (auto &[atom, mass] : definedAtoms)
  {
    com += mass * atom.position;
    total_mass += mass;
  }
  com = com / total_mass;

  double3x3 inertiaTensor{};
  for (auto &[atom, mass] : definedAtoms)
  {
    double3 dr = atom.position - com;
    inertiaTensor.ax += mass * dr.x * dr.x;
    inertiaTensor.bx += mass * dr.y * dr.x;
    inertiaTensor.cx += mass * dr.z * dr.x;
    inertiaTensor.ay += mass * dr.x * dr.y;
    inertiaTensor.by += mass * dr.y * dr.y;
    inertiaTensor.cy += mass * dr.z * dr.y;
    inertiaTensor.az += mass * dr.x * dr.z;
    inertiaTensor.bz += mass * dr.y * dr.z;
    inertiaTensor.cz += mass * dr.z * dr.z;
  }

  // the local body frame is taken to be that in which the rotational inertia tensor is diagonal
  double3 eigenvalues{};
  double3x3 eigenvectors{};
  inertiaTensor.EigenSystemSymmetric(eigenvalues, eigenvectors);

  inertiaVector = double3{0.0, 0.0, 0.0};
  atoms.clear();
  for (auto [atom, mass] : definedAtoms)
  {
    double3 dr = atom.position - com;
    double3 pos = eigenvectors.transpose() * dr;
    if (std::abs(pos.x) < 1e-8) pos.x = 0.0;
    if (std::abs(pos.y) < 1e-8) pos.y = 0.0;
    if (std::abs(pos.z) < 1e-8) pos.z = 0.0;

    inertiaVector.x += mass * (pos.y * pos.y + pos.z * pos.z);
    inertiaVector.y += mass * (pos.x * pos.x + pos.z * pos.z);
    inertiaVector.z += mass * (pos.x * pos.x + pos.y * pos.y);

    // correct the position
    atom.position = pos;
    atoms.push_back(atom);
  }

  // set axis system: Ixx >= Iyy >= Izz
  double rot_xyz = std::max({inertiaVector.x, inertiaVector.y, inertiaVector.z});
  if (rot_xyz >= inertiaVector.x)
  {
    if (inertiaVector.y >= rot_xyz)
    {
      for (auto &atom : atoms)
      {
        double temp = atom.position.x;
        atom.position.x = atom.position.y;
        atom.position.y = -temp;
      }
      inertiaVector.y = inertiaVector.x;
      inertiaVector.x = rot_xyz;
    }
    else if (inertiaVector.z >= rot_xyz)
    {
      for (auto &atom : atoms)
      {
        double temp = atom.position.x;
        atom.position.x = atom.position.z;
        atom.position.z = -temp;
      }
      inertiaVector.z = inertiaVector.x;
      inertiaVector.x = rot_xyz;
    }
  }
  if (inertiaVector.z > inertiaVector.y)
  {
    for (auto &atom : atoms)
    {
      double temp = atom.position.y;
      atom.position.y = atom.position.z;
      atom.position.z = -temp;
    }
    double temp = inertiaVector.z;
    inertiaVector.z = inertiaVector.y;
    inertiaVector.y = temp;
  }

  double rotlim = std::max(1.0e-2, inertiaVector.x + inertiaVector.y + inertiaVector.z) * 1.0e-5;

  inverseInertiaVector.x = (inertiaVector.x < rotlim) ? 0.0 : 1.0 / inertiaVector.x;
  inverseInertiaVector.y = (inertiaVector.y < rotlim) ? 0.0 : 1.0 / inertiaVector.y;
  inverseInertiaVector.z = (inertiaVector.z < rotlim) ? 0.0 : 1.0 / inertiaVector.z;

  double rotall = inertiaVector.x + inertiaVector.y + inertiaVector.z > 1.0e-5
                      ? inertiaVector.x + inertiaVector.y + inertiaVector.z
                      : 1.0;

  std::size_t index = 0;
  if (inertiaVector.x / rotall < 1.0e-5) ++index;
  if (inertiaVector.y / rotall < 1.0e-5) ++index;
  if (inertiaVector.z / rotall < 1.0e-5) ++index;

  if(rigid)
  {
    translationalDegreesOfFreedom = 3;
    rotationalDegreesOfFreedom = 3;
    if (inertiaVector.x / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
    if (inertiaVector.y / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
    if (inertiaVector.z / rotall < 1.0e-5) --rotationalDegreesOfFreedom;

    shapeType = Component::Shape{index};
  }
  else
  {
    translationalDegreesOfFreedom = 3 * definedAtoms.size();
    rotationalDegreesOfFreedom = 0;
    shapeType = Component::Shape::NonLinear;
  }
}

std::vector<Atom> Component::rotatePositions(const simd_quatd &q) const
{
  double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);
  std::vector<Atom> rotatedAtoms{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    Atom a = atoms[i];
    a.position = atoms[startingBead].position + rotationMatrix * (atoms[i].position - atoms[startingBead].position);
    rotatedAtoms.push_back(a);
  }
  return rotatedAtoms;
}

Molecule Component::createMoleculeRecord(std::span<const Atom> moleculeAtoms) const
{
  double3 com{};
  for (std::size_t i = 0; i != moleculeAtoms.size(); ++i)
  {
    com += definedAtoms[i].second * moleculeAtoms[i].position;
  }
  com = com / totalMass;

  simd_quatd orientation(0.0, 0.0, 0.0, 1.0);
  if (rigid && moleculeAtoms.size() > 1)
  {
    // 'atoms' is the COM-centered principal-frame reference geometry, so the Procrustes fit returns
    // R with p_i - com = R * atoms[i].position. The stored quaternion follows the rigid-molecule
    // convention p = com + buildRotationMatrixInverse(q) * atoms[i].position, so
    // buildRotationMatrixInverse(q) = R and q is the quaternion of R^T.
    std::vector<double3> referencePositions(moleculeAtoms.size());
    std::vector<double3> labPositions(moleculeAtoms.size());
    for (std::size_t i = 0; i != moleculeAtoms.size(); ++i)
    {
      referencePositions[i] = atoms[i].position;
      labPositions[i] = moleculeAtoms[i].position;
    }
    double3x3 rotation =
        double3x3::computeRotationMatrix(double3(0.0, 0.0, 0.0), referencePositions, com, labPositions);
    double3x3 rotationTranspose = rotation.transpose();
    orientation = rotationTranspose.quaternion().normalized();
  }

  return Molecule(com, orientation, totalMass, static_cast<std::size_t>(moleculeAtoms.front().componentId),
                  moleculeAtoms.size());
}

double3 Component::computeCenterOfMass(std::span<Atom> atom_list) const
{
  std::vector<std::pair<Atom, double>> a{definedAtoms};
  for (std::size_t i = 0; i != a.size(); ++i)
  {
    a[i].first.position = atom_list[i].position;
  }
  double3 com{};
  double total_mass{};
  for (const auto &[atom, mass] : a)
  {
    com += mass * atom.position;
    total_mass += mass;
  }
  return com / total_mass;
}

std::string Component::printStatus(std::size_t componentId, const ForceField &forceField, double inputPressure) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n\n", componentId, name);

  if (type == Component::Type::Adsorbate)
  {
    std::print(stream, "    Type: Adsorbate\n\n");
  }
  else if (type == Component::Type::Cation)
  {
    std::print(stream, "    Type: Cation\n\n");
  }

  std::print(stream, "    Critical temperature:  {} [{}]\n", criticalTemperature, Units::unitOfTemperatureString);
  std::print(stream, "    Critical pressure:     {} [{}]\n", criticalPressure, Units::unitOfPressureString);
  std::print(stream, "    Acentric factor:       {} [-]\n\n", acentricFactor);

  if(swappable)
  {
    std::print(stream, "    Mol-fraction:                 {} [-]\n", molFraction);
    std::print(stream << std::boolalpha, "    Swappable:                   {}\n\n", swappable);
    std::print(stream, "    Mass:                         {} [-]\n", totalMass);
    std::print(stream, "    Fugacity:                     {} [Pa]\n",
               molFraction * fugacityCoefficient.value_or(1.0) * inputPressure);
    if (fugacityCoefficient.has_value())
    {
      std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient.value());
    }
    std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
    std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
    std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);
  }
  else
  {
    std::print(stream << std::boolalpha, "    Swappable:                   {}\n\n", swappable);
  }

  std::print(stream, "    Number Of Atoms:              {}\n", atoms.size());
  std::print(stream, "    CBMC starting bead:           {}\n", startingBead);
  if (forceField.useRecoilGrowth && atoms.size() >= 2 && !rigid)
  {
    // Which openness reference the recoil-growth open/closed test measures against (see
    // 'recoilReferenceStepEnergies'); the fallback is valid but loses the strain compensation.
    if (recoilReferenceConformations.empty())
    {
      std::print(stream, "    Recoil-growth reference:      declared geometry (no ideal-gas reference conformations)\n");
    }
    else
    {
      std::print(stream, "    Recoil-growth reference:      {} ideal-gas conformations\n",
                 recoilReferenceConformations.size());
    }
  }
  if (idealGasRosenbluthWeight.has_value())
  {
    std::print(stream, "    Ideal gas Rosenbluth weight:  {:10.8f}\n", idealGasRosenbluthWeight.value());
  }
  if (crossSection.has_value())
  {
    std::print(stream, "    Cross-section (BET):          {:10.4f} [Å²]\n", crossSection.value());
  }
  if (liquidVolume.has_value())
  {
    std::print(stream, "    Liquid volume (Gurvich):      {:10.4f} [Å³/molecule]\n", liquidVolume.value());
  }
  if (saturationPressure.has_value())
  {
    std::print(stream, "    Saturation pressure P0:       {:10.4f} [Pa]\n", saturationPressure.value());
  }
  if (lnPartitionFunction != 0.0)
  {
    std::print(stream, "    Ln partition function:        {:10.4f} [ln(q/V) with q/V in A^-3]\n", lnPartitionFunction);
  }
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    std::size_t atomType = static_cast<std::size_t>(atoms[i].type);
    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomTypeString,
               atoms[i].position.x, atoms[i].position.y, atoms[i].position.z, atoms[i].charge);
  }
  std::print(stream, "    Diagonalized inertia-vector:      {:10.8f} {:10.8f} {:10.8f}\n", inertiaVector.x,
             inertiaVector.y, inertiaVector.z);
  std::print(stream, "    Translational degrees of freedom: {}\n", translationalDegreesOfFreedom);
  std::print(stream, "    Rotational degrees of freedom:    {}\n", rotationalDegreesOfFreedom);
  switch (shapeType)
  {
    case Shape::NonLinear:
      std::print(stream, "    Shape of the molecule is: 'Non-linear'\n");
      break;
    case Shape::Linear:
      std::print(stream, "    Shape of the molecule is: 'Linear'\n");
      break;
    case Shape::Point:
      std::print(stream, "    Shape of the molecule is: 'Point-particle'\n");
      break;
  }
  std::print(stream, "    Net-charge:      {:12.8f} [e]\n", netCharge);
  std::print(stream, "\n");

  const std::vector<double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
  std::print(stream, "    Move probabilities:\n");
  for (std::size_t i = 0; i < normalizedProbabilities.size(); ++i)
  {
    std::print(stream, "    {:<30} {:8.6f} [-]\n", Move::moveNames[i] + ":", normalizedProbabilities[i]);
  }
  std::print(stream, "\n");

  std::print(stream, "    number of blocking-pockets: {}{}\n", blockingPockets.size(),
             automaticBlockingPockets ? " (still to be computed from the framework)" : "");
  for (std::size_t i = 0; i < blockingPockets.size(); ++i)
  {
    std::print(stream, "        fractional s_x,s_y,s_z: {},{},{} radius: {}\n", blockingPockets[i].x,
               blockingPockets[i].y, blockingPockets[i].z, blockingPockets[i].w);
  }
  std::print(stream, "\n");

  std::print(stream, "    molecule is modelled as '{}'\n\n", rigid ? "rigid" : "flexible");

  if (!rigid)
  {
    {
      std::print(stream, "    Fragments:\n");
      std::print(stream, "        number of fragments:         {}\n", fragmentGraph.fragments.size());
      std::print(stream, "        rigid-body fragments:        {}\n", fragmentGraph.numberOfRigidFragments());
      std::print(stream, "        closure bonds (rings):       {}\n", fragmentGraph.closureBonds.size());
      std::print(stream << std::boolalpha, "        semi-flexible molecule:      {}\n", isSemiFlexible());
      for (std::size_t i = 0; i < fragmentGraph.fragments.size(); ++i)
      {
        const Fragment &fragment = fragmentGraph.fragments[i];
        const std::string fragmentType = fragment.isRigidBody() ? "Rigid" : "Flexible";

        std::ostringstream atomList;
        std::ostringstream atomTypeList;
        for (std::size_t j = 0; j < fragment.atoms.size(); ++j)
        {
          if (j != 0)
          {
            std::print(atomList, " ");
            std::print(atomTypeList, " ");
          }
          const std::size_t atomIndex = fragment.atoms[j];
          std::print(atomList, "{}", atomIndex);
          const std::size_t atomType = static_cast<std::size_t>(atoms[atomIndex].type);
          std::print(atomTypeList, "{}", forceField.pseudoAtoms[atomType].name);
        }

        std::print(stream, "        fragment {:3d}: {:8}  atoms: [{}]  types: [{}]\n", i, fragmentType,
                   atomList.str(), atomTypeList.str());
        if (fragment.isRigidBody())
        {
          std::print(stream, "                   mass: {:10.5f}  rotational DOF: {}\n", fragment.mass,
                     fragment.rotationalDegreesOfFreedom);
        }
      }
      std::print(stream, "\n");
    }

    std::print(stream, "    connectivity:\n");
    std::print(stream, "{}\n", connectivityTable.print("        "));

    std::size_t number_of_total_bonds = connectivityTable.findAllBonds().size();
    std::print(stream, "    number of bond potentials: {} (out of {})\n", intraMolecularPotentials.bonds.size(),
               number_of_total_bonds);
    for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
    {
      std::print(stream, "        {}", intraMolecularPotentials.bonds[i].print());
    }
    std::print(stream, "\n");

    if (!intraMolecularPotentials.ureyBradleys.empty())
    {
      std::print(stream, "    number of Urey-Bradley potentials: {}\n", intraMolecularPotentials.ureyBradleys.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.ureyBradleys.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.ureyBradleys[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bends.empty())
    {
      std::size_t number_of_total_bends = connectivityTable.findAllBends().size();
      std::print(stream, "    number of bend potentials: {} (out of {})\n", intraMolecularPotentials.bends.size(),
                 number_of_total_bends);
      for (std::size_t i = 0; i < intraMolecularPotentials.bends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.inversionBends.empty())
    {
      std::print(stream, "    number of inversion-bend potentials: {}\n",
                 intraMolecularPotentials.inversionBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.inversionBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.inversionBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.outOfPlaneBends.empty())
    {
      std::print(stream, "    number of out-of-plane bend potentials: {}\n",
                 intraMolecularPotentials.outOfPlaneBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.outOfPlaneBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.outOfPlaneBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.torsions.empty())
    {
      std::size_t number_of_total_torsions = connectivityTable.findAllTorsions().size();
      std::print(stream, "    number of torsion potentials: {} (out of {})\n", intraMolecularPotentials.torsions.size(),
                 number_of_total_torsions);
      for (std::size_t i = 0; i < intraMolecularPotentials.torsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.torsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.improperTorsions.empty())
    {
      std::print(stream, "    number of improper-torsion potentials: {}\n",
                 intraMolecularPotentials.improperTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.improperTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.improperTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondBonds.empty())
    {
      std::print(stream, "    number of bond-bond potentials: {}\n", intraMolecularPotentials.bondBonds.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondBonds.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondBonds[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondBends.empty())
    {
      std::print(stream, "    number of bond-bend potentials: {}\n", intraMolecularPotentials.bondBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondTorsions.empty())
    {
      std::print(stream, "    number of bond-torsion potentials: {}\n", intraMolecularPotentials.bondTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bendBends.empty())
    {
      std::print(stream, "    number of bend-bend potentials: {}\n", intraMolecularPotentials.bendBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bendBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bendBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bendTorsions.empty())
    {
      std::print(stream, "    number of bend-torsion potentials: {}\n", intraMolecularPotentials.bendTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bendTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bendTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.vanDerWaals.empty())
    {
      std::print(stream, "    number of Van der Waals potentials: {}\n", intraMolecularPotentials.vanDerWaals.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.vanDerWaals.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.vanDerWaals[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.coulombs.empty())
    {
      std::print(stream, "    number of coulomb potentials: {}\n", intraMolecularPotentials.coulombs.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.coulombs.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.coulombs[i].print());
      }
      std::print(stream, "\n");
    }

    if (!partialReinsertionFixedAtoms.empty())
    {
      std::print(stream, "    number of partial-reinsertion moves: {}\n", partialReinsertionFixedAtoms.size());
      for (std::size_t i = 0; i < partialReinsertionFixedAtoms.size(); ++i)
      {
        std::stringstream result{};
        std::copy(partialReinsertionFixedAtoms[i].begin(), partialReinsertionFixedAtoms[i].end(), 
                  std::ostream_iterator<std::size_t>(result, " "));
        std::print(stream, "        fixed atoms: {}\n", result.str());
      }
      std::print(stream, "\n");
    }
    if (!identityChanges.empty())
    {
      std::stringstream result{};
      std::copy(identityChanges.begin(), identityChanges.end(), std::ostream_iterator<std::size_t>(result, " "));
      std::print(stream, "    identity changes: {}\n", result.str());
    }
    if (!gibbsIdentityChanges.empty())
    {
      std::stringstream result{};
      std::copy(gibbsIdentityChanges.begin(), gibbsIdentityChanges.end(),
                std::ostream_iterator<std::size_t>(result, " "));
      std::print(stream, "    gibbs identity changes: {}\n", result.str());
    }
  }

  return stream.str();
}

nlohmann::json Component::jsonStatus() const
{
  nlohmann::json status;

  status["name"] = name;
  status["criticalTemperature"] = criticalTemperature;
  status["criticalPressure"] = criticalPressure;
  status["acentricFactor"] = acentricFactor;
  status["molFraction"] = molFraction;
  status["swappable"] = swappable;
  status["mass"] = totalMass;
  status["n_atoms"] = atoms.size();
  status["cbmcStartingBead"] = startingBead;
  status["diagonalizedInertiaVector"] = inertiaVector;
  status["translationalDOF"] = translationalDegreesOfFreedom;
  status["rotationalDOF"] = rotationalDegreesOfFreedom;
  status["netCharge"] = netCharge;

  if (fugacityCoefficient.has_value())
  {
    status["fugacityCoefficient"] = fugacityCoefficient.value();
  }

  if (lnPartitionFunction != 0.0)
  {
    status["lnPartitionFunction"] = lnPartitionFunction;
  }

  switch (shapeType)
  {
    case Shape::NonLinear:
      status["shape"] = "nonlinear";
      break;
    case Shape::Linear:
      status["shape"] = "linear";
      break;
    case Shape::Point:
      status["shape"] = "point-particle";
      break;
  }

  nlohmann::json moves;
  std::vector<double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
  status["moveProbabilities"] = normalizedProbabilities;

  status["n_bonds"] = intraMolecularPotentials.bonds.size();
  std::vector<std::string> bondTypes(intraMolecularPotentials.bonds.size());
  for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
  {
    bondTypes[i] = intraMolecularPotentials.bonds[i].print();
  }
  status["bondTypes"] = bondTypes;
  return status;
}

std::vector<Atom> Component::copiedAtoms(std::span<Atom> molecule) const
{
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].position = molecule[i].position - molecule[startingBead].position;
  }
  return copied_atoms;
}

std::pair<Molecule, std::vector<Atom>> Component::translate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                            double3 displacement) const
{
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  Molecule trialMolecule = molecule;

  if (rigid)
  {
    simd_quatd q = trialMolecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    trialMolecule.centerOfMassPosition += displacement;
    double3 com = trialMolecule.centerOfMassPosition;
    for (std::size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    trialMolecule.centerOfMassPosition += displacement;
    for (Atom &atom : trialAtoms)
    {
      atom.position += displacement;
    }
  }

  return {trialMolecule, trialAtoms};
}

std::pair<Molecule, std::vector<Atom>> Component::rotate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                         simd_quatd rotation) const
{
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  Molecule trialMolecule = molecule;

  if (rigid)
  {
    simd_quatd q = rotation * trialMolecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    trialMolecule.orientation = q;
    double3 com = trialMolecule.centerOfMassPosition;
    for (std::size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    double3x3 rotationMatrix = double3x3(rotation);
    for (std::size_t i = 0; i < trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = rotationMatrix * (molecule_atoms[i].position - molecule_atoms[startingBead].position) +
                               molecule_atoms[startingBead].position;
    }

    // A flexible molecule is rotated rigidly about its starting bead, not its center of mass, so the
    // center-of-mass moves. Recompute it from the rotated atoms to keep the stored center-of-mass consistent
    // with the atom positions (it is otherwise stale after an accepted rotation).
    trialMolecule.centerOfMassPosition = computeCenterOfMass(trialAtoms);
  }

  return {trialMolecule, trialAtoms};
}


ConnectivityTable Component::readConnectivityTable(std::size_t size,
                                                   const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  ConnectivityTable table(size);

  if (parsed_data.contains("Connectivity"))
  {
    for (auto &[_, item] : parsed_data["Connectivity"].items())
    {
      if (!item.is_array())
      {
        throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
      }

      if (item.size() != 2)
      {
        throw std::runtime_error(
            std::format("[Component reader]: item {} must be an array with two unsigned integers.\n", item.dump()));
      }

      try
      {
        std::vector<std::size_t> identifiers = item.get<std::vector<std::size_t>>();

        table[identifiers[0], identifiers[1]] = true;
        table[identifiers[1], identifiers[0]] = true;
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in connectivities ({}): {}\n", item.dump(), e.what()));
      }
    }
  }
  return table;
}

std::vector<BondPotential> Component::readBondPotentials(const ForceField &forceField,
                                                         const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<BondPotential> bond_potentials{};

  std::vector<std::array<std::size_t, 2>> found_bonds = connectivityTable.findAllBonds();

  // Connectivity edges entirely inside a rigid group carry no bond potential: the geometry of a
  // rigid group is fixed by its reference atoms and is never sampled.
  std::erase_if(found_bonds, [&](const std::array<std::size_t, 2> &ids) { return isInsideRigidFragment(ids); });

  if (parsed_data.contains("Bonds"))
  {
    intraMolecularPotentials.bonds.reserve(found_bonds.size());

    // fill in each of the bonds
    for (const std::array<std::size_t, 2> found_bond : found_bonds)
    {
      // check which bond-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool bond_identified{ false };

      for (auto &[_, item] : parsed_data["Bonds"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the two identifiers of the bond, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);

            if (!(type_A.has_value() && type_B.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {}\n", identifiers[0], identifiers[1]));
            }

            if ((atoms[found_bond[0]].type == type_A.value() && atoms[found_bond[1]].type == type_B.value()) ||
                (atoms[found_bond[0]].type == type_B.value() && atoms[found_bond[1]].type == type_A.value()))
            {
              if (bond_identified)
              {
                throw std::runtime_error(std::format("Error in Bond-potential: ambiguity in define bond-potential (bond {}-{} matches on multiple definitions)\n",
                      found_bond[0], found_bond[1]));
              }

              bond_identified = true;

              BondPotential bond_potential =
                  BondPotential({found_bond[0], found_bond[1]}, BondPotential::definitionForString.at(potential_name),
                                potential_parameters);
              bond_potentials.push_back(bond_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            if ((found_bond[0] == A && found_bond[1] == B) ||
                (found_bond[0] == B && found_bond[1] == A))
            {
              if (bond_identified)
              {
                throw std::runtime_error(std::format("Error in Bond-potential: ambiguity in define bond-potential (bond {}-{} matches on multiple definitions)\n",
                      found_bond[0], found_bond[1]));
              }

              bond_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> bond =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              BondPotential bond_potential = BondPotential(
                  {bond[0], bond[1]}, BondPotential::definitionForString.at(potential_name), potential_parameters);
              bond_potentials.push_back(bond_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Bond-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_bonds.size() != bond_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Bond-potential: not all bond-potentials defined ({} missing)\n",
                                         found_bonds.size() - bond_potentials.size()));
  }

  return bond_potentials;
}

std::vector<BendPotential> Component::readBendPotentials(const ForceField &forceField,
                                                         const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<BendPotential> bend_potentials{};

  std::vector<std::array<std::size_t, 3>> found_bends = connectivityTable.findAllBends();

  // Bends entirely inside a rigid group are never sampled; their angle is fixed by the rigid geometry.
  std::erase_if(found_bends, [&](const std::array<std::size_t, 3> &ids) { return isInsideRigidFragment(ids); });

  if (parsed_data.contains("Bends"))
  {
    intraMolecularPotentials.bends.reserve(found_bends.size());

    // fill in each of the bends
    for (const std::array<std::size_t, 3> found_bend : found_bends)
    {
      // check which bend-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool bend_identified{ false };

      for (auto &[_, item] : parsed_data["Bends"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the three identifiers of the bend, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string() && item[0][2].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);
            std::optional<std::size_t> type_C = forceField.findPseudoAtom(identifiers[2]);

            if (!(type_A.has_value() && type_B.has_value() && type_C.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {} {}\n", identifiers[0], identifiers[1], identifiers[2]));
            }

            if ((atoms[found_bend[0]].type == type_A.value() && atoms[found_bend[1]].type == type_B.value() && atoms[found_bend[2]].type == type_C.value()) ||
                (atoms[found_bend[0]].type == type_C.value() && atoms[found_bend[1]].type == type_B.value() && atoms[found_bend[2]].type == type_A.value()))
            {
              if (bend_identified)
              {
                throw std::runtime_error(std::format("Error in Bend-potential: ambiguity in define bend-potential (bend {}-{}-{} matches on multiple definitions)\n",
                      found_bend[0], found_bend[1], found_bend[2]));
              }

              bend_identified = true;

              BendPotential bend_potential =
                  BendPotential({found_bend[0], found_bend[1], found_bend[2]}, BendPotential::definitionForString.at(potential_name),
                                potential_parameters);
              bend_potentials.push_back(bend_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned() && item[0][2].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            std::size_t C = item[0][2].get<std::size_t>();
            if ((found_bend[0] == A && found_bend[1] == B && found_bend[2] == C) ||
                (found_bend[0] == C && found_bend[1] == B && found_bend[2] == A))
            {
              if (bend_identified)
              {
                throw std::runtime_error(std::format("Error in Bend-potential: ambiguity in define bend-potential (bend {}-{}-{} matches on multiple definitions)\n",
                      found_bend[0], found_bend[1], found_bend[2]));
              }

              bend_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> bend =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              BendPotential bend_potential = BendPotential(
                  {bend[0], bend[1], bend[2]}, BendPotential::definitionForString.at(potential_name), potential_parameters);
              bend_potentials.push_back(bend_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Bend-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_bends.size() != bend_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Bend-potential: not all bend-potentials defined ({} missing)\n",
                                         found_bends.size() - bend_potentials.size()));
  }

  return bend_potentials;
}

std::vector<TorsionPotential> Component::readTorsionPotentials(
    const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<TorsionPotential> torsion_potentials{};

  std::vector<std::array<std::size_t, 4>> found_torsions = connectivityTable.findAllTorsions();

  // Torsions entirely inside a rigid group are never sampled; their angle is fixed by the rigid geometry.
  std::erase_if(found_torsions, [&](const std::array<std::size_t, 4> &ids) { return isInsideRigidFragment(ids); });

  if (parsed_data.contains("Torsions"))
  {
    intraMolecularPotentials.torsions.reserve(found_torsions.size());

    // fill in each of the torsions
    for (const std::array<std::size_t, 4> found_torsion : found_torsions)
    {
      // check which torsion-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool torsion_identified{ false };

      for (auto &[_, item] : parsed_data["Torsions"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the four identifiers of the torsion, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string() && item[0][2].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);
            std::optional<std::size_t> type_C = forceField.findPseudoAtom(identifiers[2]);
            std::optional<std::size_t> type_D = forceField.findPseudoAtom(identifiers[3]);

            if (!(type_A.has_value() && type_B.has_value() && type_C.has_value() && type_D.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {} {} {}\n", identifiers[0], identifiers[1], identifiers[2], identifiers[3]));
            }

            if ((atoms[found_torsion[0]].type == type_A.value() && atoms[found_torsion[1]].type == type_B.value() && 
                 atoms[found_torsion[2]].type == type_C.value() && atoms[found_torsion[3]].type == type_D.value()) ||
                (atoms[found_torsion[0]].type == type_D.value() && atoms[found_torsion[1]].type == type_C.value() && 
                 atoms[found_torsion[2]].type == type_B.value() && atoms[found_torsion[3]].type == type_A.value()))
            {
              if (torsion_identified)
              {
                throw std::runtime_error(std::format("Error in Torsion-potential: ambiguity in define torsion-potential (torsion {}-{}-{}-{} matches on multiple definitions)\n",
                      found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]));
              }

              torsion_identified = true;

              TorsionPotential torsion_potential =
                  TorsionPotential({found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]}, TorsionPotential::definitionForString.at(potential_name),
                                potential_parameters);
              torsion_potentials.push_back(torsion_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned() && item[0][2].is_number_unsigned() && item[0][3].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            std::size_t C = item[0][2].get<std::size_t>();
            std::size_t D = item[0][3].get<std::size_t>();
            if ((found_torsion[0] == A && found_torsion[1] == B && found_torsion[2] == C && found_torsion[3] == D) ||
                (found_torsion[0] == D && found_torsion[1] == C && found_torsion[2] == B && found_torsion[3] == A))
            {
              if (torsion_identified)
              {
                throw std::runtime_error(std::format("Error in Torsion-potential: ambiguity in define torsion-potential (torsion {}-{}-{}-{} matches on multiple definitions)\n",
                      found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]));
              }

              torsion_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> torsion =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              TorsionPotential torsion_potential = TorsionPotential(
                  {torsion[0], torsion[1], torsion[2], torsion[3]}, TorsionPotential::definitionForString.at(potential_name), potential_parameters);
              torsion_potentials.push_back(torsion_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Torsion-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_torsions.size() != torsion_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Torsion-potential: not all torsion-potentials defined ({} missing)\n",
                                         found_torsions.size() - torsion_potentials.size()));
  }

  return torsion_potentials;
}

// Generic reader for the internal potentials that are explicitly listed in the molecule-file with
// atom indices (Urey-Bradley, inversion bends, improper torsions, and the cross terms). Each JSON
// entry has the form [[id_0, ..., id_{N-1}], "TYPE", [p_0, p_1, ...]].
template <typename PotentialType, std::size_t N, typename TypeEnum>
static std::vector<PotentialType> readExplicitInternalPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data, const std::string &key,
    const std::map<std::string, TypeEnum, caseInsensitiveComparator> &definitionForString)
{
  std::vector<PotentialType> potentials{};

  if (!parsed_data.contains(key)) return potentials;

  for (auto &[_, item] : parsed_data[key].items())
  {
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    if (item.size() != 3)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with three elements, "
                      "(1) an array with the {} identifiers, (2) the potential type, "
                      "and (3) an array containing the potential parameters\n",
                      item.dump(), N));
    }

    try
    {
      std::vector<std::size_t> identifiers = item[0].get<std::vector<std::size_t>>();
      if (identifiers.size() != N)
      {
        throw std::runtime_error(std::format("expected {} identifiers, got {}", N, identifiers.size()));
      }

      std::string potential_name = item[1].get<std::string>();
      if (!definitionForString.contains(potential_name))
      {
        throw std::runtime_error(std::format("unknown potential type '{}'", potential_name));
      }

      std::vector<double> potential_parameters =
          item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

      std::array<std::size_t, N> ids{};
      std::copy(identifiers.begin(), identifiers.end(), ids.begin());

      potentials.push_back(PotentialType(ids, definitionForString.at(potential_name), potential_parameters));
    }
    catch (std::exception const &e)
    {
      throw std::runtime_error(std::format("Error in {}-potential ({}): {}\n", key, item.dump(), e.what()));
    }
  }

  return potentials;
}

std::vector<UreyBradleyPotential> Component::readUreyBradleyPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<UreyBradleyPotential, 2>(parsed_data, "UreyBradleys",
                                                                 UreyBradleyPotential::definitionForString);
}

std::vector<InversionBendPotential> Component::readInversionBendPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<InversionBendPotential, 4>(parsed_data, "InversionBends",
                                                                   InversionBendPotential::definitionForString);
}

std::vector<OutOfPlaneBendPotential> Component::readOutOfPlaneBendPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<OutOfPlaneBendPotential, 4>(parsed_data, "OutOfPlaneBends",
                                                                    OutOfPlaneBendPotential::definitionForString);
}

std::vector<TorsionPotential> Component::readImproperTorsionPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<TorsionPotential, 4>(parsed_data, "ImproperTorsions",
                                                             TorsionPotential::definitionForString);
}

std::vector<BondBondPotential> Component::readBondBondPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<BondBondPotential, 3>(parsed_data, "BondBonds",
                                                              BondBondPotential::definitionForString);
}

std::vector<BondBendPotential> Component::readBondBendPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<BondBendPotential, 4>(parsed_data, "BondBends",
                                                              BondBendPotential::definitionForString);
}

std::vector<BondTorsionPotential> Component::readBondTorsionPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<BondTorsionPotential, 4>(parsed_data, "BondTorsions",
                                                                 BondTorsionPotential::definitionForString);
}

std::vector<BendBendPotential> Component::readBendBendPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<BendBendPotential, 4>(parsed_data, "BendBends",
                                                              BendBendPotential::definitionForString);
}

std::vector<BendTorsionPotential> Component::readBendTorsionPotentials(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  return readExplicitInternalPotentials<BendTorsionPotential, 4>(parsed_data, "BendTorsions",
                                                                 BendTorsionPotential::definitionForString);
}

std::vector<VanDerWaalsPotential> Component::readVanDerWaalsPotentials(
    const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<VanDerWaalsPotential> van_der_waals_potentials{};

  if (parsed_data.contains("IntraVanDerWaalsPairs") &&
      parsed_data.contains("Intra14VanDerWaalsScalingValue"))
  {
    throw std::runtime_error(
        "[Component reader]: IntraVanDerWaalsPairs and Intra14VanDerWaalsScalingValue are mutually exclusive\n");
  }

  if (parsed_data.contains("Intra14VanDerWaalsScalingValue"))
  {
    if(parsed_data["Intra14VanDerWaalsScalingValue"].is_number_float())
    {
      double scaling = parsed_data["Intra14VanDerWaalsScalingValue"].get<double>();

      if(scaling > 0.0)
      {
        std::vector<std::array<std::size_t, 4>> found_14_van_der_waals = connectivityTable.findAllTorsions();
        for (std::array<std::size_t, 4> &found_14_van_der_waal : found_14_van_der_waals)
        {
          std::size_t A = found_14_van_der_waal[0];
          std::size_t B = found_14_van_der_waal[3];

          // pairs inside the same rigid group have a fixed distance; their energy is a constant
          if (isInsideRigidFragment(std::array<std::size_t, 2>{A, B})) continue;

          std::size_t typeA = static_cast<std::size_t>(atoms[A].type);
          std::size_t typeB = static_cast<std::size_t>(atoms[B].type);
          
          [[maybe_unused]] VDWParameters::Type potentialType = forceField(typeA, typeB).type;
          double4 parameters = forceField(typeA, typeB).parameters;
          //double shift = forceField(typeA, typeB).shift;

          // FIX: unit conversion
          VanDerWaalsPotential potential = VanDerWaalsPotential(
              {A, B}, VanDerWaalsType::LennardJones,
              {parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, scaling);

          van_der_waals_potentials.push_back(potential);
        }
      }
    }
  }

  std::vector<std::array<std::size_t, 2>> found_van_der_waals =
      parsed_data.contains("IntraVanDerWaalsPairs")
          ? readExplicitIntraPairs(parsed_data, "IntraVanDerWaalsPairs", atoms.size())
          : connectivityTable.findAllVanDerWaals();

  for (std::array<std::size_t, 2> &found_van_der_waal : found_van_der_waals)
  {
    std::size_t A = found_van_der_waal[0];
    std::size_t B = found_van_der_waal[1];

    // pairs inside the same rigid group have a fixed distance; their energy is a constant
    if (isInsideRigidFragment(found_van_der_waal)) continue;

    std::size_t typeA = static_cast<std::size_t>(atoms[A].type);
    std::size_t typeB = static_cast<std::size_t>(atoms[B].type);

    [[maybe_unused]] VDWParameters::Type potentialType = forceField(typeA, typeB).type;
    double4 parameters = forceField(typeA, typeB).parameters;

    // FIX: unit conversion
    VanDerWaalsPotential potential = VanDerWaalsPotential(
        {A, B}, VanDerWaalsType::LennardJones,
        {parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, 1.0);

    van_der_waals_potentials.push_back(potential);
  }

  return van_der_waals_potentials;
}

std::vector<CoulombPotential> Component::readCoulombPotentials(
    const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<CoulombPotential> coulomb_potentials{};

  if(!forceField.useCharge) return coulomb_potentials;

  if (parsed_data.contains("IntraCoulombPairs") &&
      parsed_data.contains("Intra14ChargeChargeScalingValue"))
  {
    throw std::runtime_error(
        "[Component reader]: IntraCoulombPairs and Intra14ChargeChargeScalingValue are mutually exclusive\n");
  }

  if (parsed_data.contains("Intra14ChargeChargeScalingValue"))
  {
    if(parsed_data["Intra14ChargeChargeScalingValue"].is_number_float())
    {
      double scaling = parsed_data["Intra14ChargeChargeScalingValue"].get<double>();

      if(scaling > 0.0)
      {
        std::vector<std::array<std::size_t, 4>> found_14_coulombs = connectivityTable.findAllTorsions();
        for (std::array<std::size_t, 4> &found_14_coulomb : found_14_coulombs)
        {
          std::size_t A = found_14_coulomb[0];
          std::size_t B = found_14_coulomb[3];

          if (isInsideRigidFragment(std::array<std::size_t, 2>{A, B})) continue;

          double chargeA = atoms[A].charge;
          double chargeB = atoms[B].charge;
          
          CoulombPotential potential = CoulombPotential({A, B}, CoulombType::Coulomb, chargeA, chargeB, scaling);

          coulomb_potentials.push_back(potential);
        }
      }
    }
  }

  std::vector<std::array<std::size_t, 2>> found_coulombs =
      parsed_data.contains("IntraCoulombPairs")
          ? readExplicitIntraPairs(parsed_data, "IntraCoulombPairs", atoms.size())
          : connectivityTable.findAllVanDerWaals();

  for (std::array<std::size_t, 2> &found_coulomb : found_coulombs)
  {
    std::size_t A = found_coulomb[0];
    std::size_t B = found_coulomb[1];

    if (isInsideRigidFragment(found_coulomb)) continue;

    double chargeA = atoms[A].charge;
    double chargeB = atoms[B].charge;

    CoulombPotential potential = CoulombPotential({A, B}, CoulombType::Coulomb, chargeA, chargeB, 1.0);

    coulomb_potentials.push_back(potential);
  }

  return coulomb_potentials;
}

std::vector<std::vector<std::size_t>> Component::readRepeatUnits(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  if (!parsed_data.contains("RepeatUnits")) return {};

  std::vector<std::vector<std::size_t>> units{};
  try
  {
    units = parsed_data["RepeatUnits"].get<std::vector<std::vector<std::size_t>>>();
  }
  catch (std::exception const &e)
  {
    throw std::runtime_error(
        std::format("[Component reader]: 'RepeatUnits' must be an array of atom-index lists: {}\n", e.what()));
  }

  const std::size_t numberOfBeads = definedAtoms.size();
  const std::size_t numberOfUnits = units.size();
  if (numberOfUnits < 2)
  {
    throw std::runtime_error(std::format("[Component reader]: 'RepeatUnits' needs at least two repeat units\n"));
  }
  const std::size_t unitSize = units.front().size();
  if (unitSize == 0)
  {
    throw std::runtime_error(std::format("[Component reader]: 'RepeatUnits' entries must not be empty\n"));
  }

  // The units must form an ordered partition of the atoms into equally sized blocks.
  std::vector<std::size_t> unitOf(numberOfBeads, 0);
  std::vector<std::size_t> offsetOf(numberOfBeads, 0);
  std::vector<bool> seen(numberOfBeads, false);
  for (std::size_t k = 0; k != numberOfUnits; ++k)
  {
    if (units[k].size() != unitSize)
    {
      throw std::runtime_error(std::format("[Component reader]: all 'RepeatUnits' must have the same size\n"));
    }
    for (std::size_t j = 0; j != unitSize; ++j)
    {
      std::size_t atom = units[k][j];
      if (atom >= numberOfBeads || seen[atom])
      {
        throw std::runtime_error(
            std::format("[Component reader]: 'RepeatUnits' must partition the atom indices (atom {} invalid or "
                        "listed twice)\n",
                        atom));
      }
      seen[atom] = true;
      unitOf[atom] = k;
      offsetOf[atom] = j;
    }
  }
  if (numberOfUnits * unitSize != numberOfBeads)
  {
    throw std::runtime_error(std::format("[Component reader]: 'RepeatUnits' must cover every atom of the molecule\n"));
  }

  // Shift periodicity of the chemistry: atom j of every unit must have the same pseudo-atom type
  // and charge. Reptation relabels physical unit k+1 into the slots of unit k, so a chain with
  // distinct end groups (caps, initiator fragments) is chemically not shift-periodic and cannot
  // reptate by relabeling.
  for (std::size_t k = 1; k != numberOfUnits; ++k)
  {
    for (std::size_t j = 0; j != unitSize; ++j)
    {
      const Atom &current = definedAtoms[units[k][j]].first;
      const Atom &previous = definedAtoms[units[k - 1][j]].first;
      if (current.type != previous.type || current.charge != previous.charge)
      {
        throw std::runtime_error(
            std::format("[Component reader]: 'RepeatUnits' are not shift-periodic: atom {} of unit {} differs in "
                        "type or charge from atom {} of unit {}\n",
                        units[k][j], k, units[k - 1][j], k - 1));
      }
    }
  }

  auto shiftDown = [&](std::size_t atom) -> std::size_t { return units[unitOf[atom] - 1][offsetOf[atom]]; };

  // Canonical form of an identifier tuple: a term and its reverse describe the same interaction.
  auto normalized = [](std::vector<std::size_t> ids) -> std::vector<std::size_t>
  {
    std::vector<std::size_t> reversed(ids.rbegin(), ids.rend());
    return std::min(ids, reversed);
  };

  // Shift periodicity of the topology: the set of terms living entirely in units 0..M-2 must equal
  // the set of terms in units 1..M-1 shifted down by one unit. This single set equality covers the
  // interior periodicity and the equality of the head and tail connections at once. (Potential
  // parameters are generated from pseudo-atom type patterns, so identical types plus identical
  // index topology imply identical parameters.)
  auto checkPeriodicTerms = [&](const auto &identifierLists, std::string_view name)
  {
    std::set<std::vector<std::size_t>> lowTerms;
    std::set<std::vector<std::size_t>> shiftedTerms;
    for (const auto &identifiers : identifierLists)
    {
      std::vector<std::size_t> ids(identifiers.begin(), identifiers.end());
      if (std::ranges::all_of(ids, [&](std::size_t atom) { return unitOf[atom] <= numberOfUnits - 2; }))
      {
        lowTerms.insert(normalized(ids));
      }
      if (std::ranges::all_of(ids, [&](std::size_t atom) { return unitOf[atom] >= 1; }))
      {
        std::vector<std::size_t> shifted(ids.size());
        std::ranges::transform(ids, shifted.begin(), shiftDown);
        shiftedTerms.insert(normalized(shifted));
      }
    }
    if (lowTerms != shiftedTerms)
    {
      throw std::runtime_error(std::format(
          "[Component reader]: 'RepeatUnits' are not shift-periodic: the {} terms do not map onto themselves "
          "under the one-unit shift\n",
          name));
    }
  };

  checkPeriodicTerms(connectivityTable.findAllBonds(), "connectivity");
  checkPeriodicTerms(intraMolecularPotentials.bonds |
                         std::views::transform([](const BondPotential &p) { return p.identifiers; }),
                     "bond");
  checkPeriodicTerms(intraMolecularPotentials.bends |
                         std::views::transform([](const BendPotential &p) { return p.identifiers; }),
                     "bend");
  checkPeriodicTerms(intraMolecularPotentials.torsions |
                         std::views::transform([](const TorsionPotential &p) { return p.identifiers; }),
                     "torsion");
  checkPeriodicTerms(intraMolecularPotentials.vanDerWaals |
                         std::views::transform([](const VanDerWaalsPotential &p) { return p.identifiers; }),
                     "intramolecular van-der-Waals");
  checkPeriodicTerms(intraMolecularPotentials.coulombs |
                         std::views::transform([](const CoulombPotential &p) { return p.identifiers; }),
                     "intramolecular Coulomb");

  // Rigid fragments must lie inside a single repeat unit, map onto rigid fragments under the shift,
  // and be congruent to their shifted counterparts (the relabeled body must fit the slot's rigid
  // reference geometry exactly).
  for (const auto &fragment : fragmentGraph.fragments)
  {
    if (!fragment.isRigidBody()) continue;
    std::size_t unit = unitOf[fragment.atoms.front()];
    for (std::size_t atom : fragment.atoms)
    {
      if (unitOf[atom] != unit)
      {
        throw std::runtime_error(std::format(
            "[Component reader]: 'RepeatUnits' invalid: rigid fragment straddles two repeat units (atom {})\n",
            atom));
      }
    }
    if (unit == 0) continue;

    const auto &shiftedFragment = fragmentGraph.fragments[fragmentGraph.atomFragmentIds[shiftDown(fragment.atoms.front())]];
    std::set<std::size_t> shiftedAtoms{};
    for (std::size_t atom : fragment.atoms) shiftedAtoms.insert(shiftDown(atom));
    if (std::set<std::size_t>(shiftedFragment.atoms.begin(), shiftedFragment.atoms.end()) != shiftedAtoms)
    {
      throw std::runtime_error(std::format(
          "[Component reader]: 'RepeatUnits' are not shift-periodic: rigid fragments do not map onto each other\n"));
    }
    for (std::size_t a = 0; a != fragment.atoms.size(); ++a)
    {
      for (std::size_t b = a + 1; b != fragment.atoms.size(); ++b)
      {
        double distance =
            (definedAtoms[fragment.atoms[a]].first.position - definedAtoms[fragment.atoms[b]].first.position).length();
        double shiftedDistance = (definedAtoms[shiftDown(fragment.atoms[a])].first.position -
                                  definedAtoms[shiftDown(fragment.atoms[b])].first.position)
                                     .length();
        if (std::abs(distance - shiftedDistance) > 1e-6)
        {
          throw std::runtime_error(
              std::format("[Component reader]: 'RepeatUnits' are not shift-periodic: rigid fragment reference "
                          "geometries differ between units\n"));
        }
      }
    }
  }

  // The chain that remains after removing an end unit must be connected, so that the CBMC growth of
  // the arriving unit has a well-defined anchored plan.
  std::array<std::vector<std::size_t>, 2> endPlacedSets{};
  for (std::size_t end = 0; end != 2; ++end)
  {
    const std::vector<std::size_t> &endUnit = (end == 0) ? units.front() : units.back();
    std::vector<bool> inEndUnit(numberOfBeads, false);
    for (std::size_t atom : endUnit) inEndUnit[atom] = true;
    for (std::size_t atom = 0; atom != numberOfBeads; ++atom)
    {
      if (!inEndUnit[atom]) endPlacedSets[end].push_back(atom);
    }
    if (!connectivityTable.checkIsConnectedSubgraph(endPlacedSets[end]))
    {
      throw std::runtime_error(std::format(
          "[Component reader]: 'RepeatUnits' invalid: removing an end unit disconnects the molecule\n"));
    }
  }

  // Reptation pairs the Rosenbluth weight of a tail-plan growth against a head-plan retrace (and
  // vice versa). Flexible acyclic steps draw their trial conformations from the exact bonded-
  // Boltzmann base sampler, so that pairing satisfies detailed balance for ANY pair of plans: no
  // congruence between the two end plans is required. Rigid-body and ring-closure steps still use
  // approximate internal-MC samplers whose deviations depend on the step structure; those deviations
  // cancel only when grow and retrace share the plan, so repeat units whose end plans contain such
  // steps must have congruent head and tail plans (same step sequence: step kind, rigid-body flag,
  // number and types of grown beads, and per-step counts of the filtered potential terms).
  {
    const std::vector<CBMC::GrowStep> &headPlan = growthPlan(endPlacedSets[0]);
    const std::vector<CBMC::GrowStep> &tailPlan = growthPlan(endPlacedSets[1]);
    auto usesApproximateSampler = [](const CBMC::GrowStep &step)
    { return step.rigidBody || step.kind == CBMC::GrowStep::Kind::CloseRing; };
    bool needsCongruence = std::any_of(headPlan.begin(), headPlan.end(), usesApproximateSampler) ||
                           std::any_of(tailPlan.begin(), tailPlan.end(), usesApproximateSampler);
    if (needsCongruence)
    {
      auto stepSignature = [&](const CBMC::GrowStep &step)
      {
        std::vector<std::size_t> grownTypes{};
        for (std::size_t bead : step.nextBeads) grownTypes.push_back(definedAtoms[bead].first.type);
        std::sort(grownTypes.begin(), grownTypes.end());
        return std::make_tuple(step.kind, step.rigidBody, grownTypes, step.intra.bonds.size(), step.intra.bends.size(),
                               step.intra.torsions.size(), step.intra.vanDerWaals.size(), step.intra.coulombs.size());
      };
      bool congruent = headPlan.size() == tailPlan.size();
      for (std::size_t i = 0; congruent && i != headPlan.size(); ++i)
      {
        congruent = stepSignature(headPlan[i]) == stepSignature(tailPlan[i]);
      }
      if (!congruent)
      {
        throw std::runtime_error(std::format(
            "[Component reader]: 'RepeatUnits' invalid for reptation: the repeat unit contains rigid-body or "
            "ring-closure growth steps and its head and tail growth plans are not congruent. Those steps use "
            "approximate internal-MC trial samplers whose deviations only cancel between grow and retrace on the "
            "same plan; pairing them across non-congruent plans (as reptation does) would bias the sampled "
            "intramolecular distributions\n"));
      }
    }
  }

  return units;
}

std::optional<std::array<std::size_t, 2>> Component::determineEndToEndAtoms(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data) const
{
  const std::size_t numberOfBeads = definedAtoms.size();

  // Explicit specification always wins.
  if (parsed_data.contains("EndToEndAtoms"))
  {
    std::array<std::size_t, 2> pair{};
    try
    {
      pair = parsed_data["EndToEndAtoms"].get<std::array<std::size_t, 2>>();
    }
    catch (std::exception const &e)
    {
      throw std::runtime_error(
          std::format("[Component reader]: 'EndToEndAtoms' must be an array of two atom indices: {}\n", e.what()));
    }
    if (pair[0] >= numberOfBeads || pair[1] >= numberOfBeads || pair[0] == pair[1])
    {
      throw std::runtime_error(std::format(
          "[Component reader]: 'EndToEndAtoms' [{}, {}] must be two distinct atom indices below {}\n", pair[0],
          pair[1], numberOfBeads));
    }
    return pair;
  }

  if (numberOfBeads < 2) return std::nullopt;

  // A chain declared with 'RepeatUnits': the ends are the backbone entry bead of the first unit and
  // the backbone exit bead of the last unit. The units are shift-periodic, so the (unique up to
  // slots) bond crossing the first junction identifies both backbone slots: the slot bonding
  // FORWARD into the next unit (exit) and the slot receiving that bond (entry). This deliberately
  // skips side-chain atoms, which can form the longest graph path but are not the polymer ends.
  if (repeatUnits.size() >= 2)
  {
    const std::vector<std::size_t> &firstUnit = repeatUnits.front();
    const std::vector<std::size_t> &secondUnit = repeatUnits[1];
    for (std::size_t slotA = 0; slotA != firstUnit.size(); ++slotA)
    {
      for (std::size_t slotB = 0; slotB != secondUnit.size(); ++slotB)
      {
        if (connectivityTable[firstUnit[slotA], secondUnit[slotB]])
        {
          // slotA bonds forward (exit), slotB receives (entry).
          return std::array<std::size_t, 2>{firstUnit[slotB], repeatUnits.back()[slotA]};
        }
      }
    }
    // Validated repeat units are always connected; not reached.
  }

  // General molecule: the endpoints of the topological diameter of the bond graph (the two termini
  // of a linear chain; for branched molecules an explicit 'EndToEndAtoms' may be preferable).
  return connectivityTable.graphDiameterEndpoints();
}

std::vector<std::vector<std::size_t>> Component::readPartialReinsertionFixedAtoms(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<std::vector<std::size_t>> config_moves{};

  if (parsed_data.contains("Partial-reinsertion"))
  {
    auto item = parsed_data["Partial-reinsertion"];
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    try
    {
      config_moves = parsed_data["Partial-reinsertion"].get<std::vector<std::vector<std::size_t>>>();
    }
    catch (std::exception const &e)
    {
      throw std::runtime_error(std::format("Error in defined partial reinsertion ({}): {}\n", item.dump(), e.what()));
    }
  }

  // A fixed set need not be connected: fixed atoms on both sides of a regrown part make it a
  // fixed-endpoint regrowth (the interior segment is grown from one side and closed onto the other,
  // see cbmc_bridge_closure). It must be a proper, non-empty set of valid atoms that leaves something
  // to regrow; whether the regrown part can be grown (no rigid body or ring closed onto a second
  // fixed atom, no bead bonded to three fixed atoms) is decided when its growth plan is built below.
  for (const std::vector<std::size_t> &config_move : config_moves)
  {
    std::stringstream result{};
    std::copy(config_move.begin(), config_move.end(), std::ostream_iterator<std::size_t>(result, " "));

    if (config_move.empty())
    {
      throw std::runtime_error("Error in defined partial reinsertion (an empty fixed set)\n");
    }
    for (std::size_t atom : config_move)
    {
      if (atom >= connectivityTable.numberOfBeads)
      {
        throw std::runtime_error(std::format(
            "Error in defined partial reinsertion ({}): atom index {} out of range (molecule has {} atoms)\n",
            result.str(), atom, connectivityTable.numberOfBeads));
      }
    }
    if (std::set<std::size_t>(config_move.begin(), config_move.end()).size() != config_move.size())
    {
      throw std::runtime_error(
          std::format("Error in defined partial reinsertion ({}): an atom is listed twice\n", result.str()));
    }
    if (config_move.size() >= connectivityTable.numberOfBeads)
    {
      throw std::runtime_error(
          std::format("Error in defined partial reinsertion ({}): every atom is fixed, nothing to regrow\n",
                      result.str()));
    }
  }

  for (const std::vector<std::size_t> &config_move : config_moves)
  {
    std::vector<std::size_t> beads_already_placed = config_move;

    do
    {
      // The growth is fragment-aware: each rigid body is grown as a single unit and the flexible
      // beads follow the fragment-graph spanning tree, so a partial-reinsertion fixed set must allow
      // the growth to complete one fragment at a time.
      std::optional<std::vector<std::size_t>> nextBeads =
          fragmentGraph.nextGrowBeads(beads_already_placed, connectivityTable);
      if (nextBeads.has_value())
      {
        beads_already_placed.insert(beads_already_placed.end(), nextBeads->begin(), nextBeads->end());
      }
      else
      {
        std::stringstream result{};
        std::copy(config_move.begin(), config_move.end(),
                  std::ostream_iterator<std::size_t>(result, " "));

        throw std::runtime_error(
            std::format("Error in defined partial reinsertion\n"
                        "{} does not satisfy requirement that all branches need to be grown at the same time\n",
                        result.str()));
      }
    } while (beads_already_placed.size() < connectivityTable.numberOfBeads);
  }

  return config_moves;
}

std::vector<ChiralCenter> Component::readChiralCenters(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<ChiralCenter> chiral_centers{};

  if (!parsed_data.contains("ChiralCenters"))
  {
    return chiral_centers;
  }

  if (!parsed_data["ChiralCenters"].is_array())
  {
    throw std::runtime_error(
        std::format("[Component reader]: 'ChiralCenters' must be an array of atom-index quadruples\n"));
  }

  std::size_t numberOfBeads = definedAtoms.size();
  for (auto &[_, item] : parsed_data["ChiralCenters"].items())
  {
    std::array<std::size_t, 4> ids{};
    try
    {
      ids = item.get<std::array<std::size_t, 4>>();
    }
    catch (std::exception const &e)
    {
      throw std::runtime_error(std::format(
          "[Component reader]: each entry of 'ChiralCenters' must be an array of four atom indices "
          "[center, neighbor1, neighbor2, neighbor3] ({}): {}\n",
          item.dump(), e.what()));
    }

    for (std::size_t id : ids)
    {
      if (id >= numberOfBeads)
      {
        throw std::runtime_error(
            std::format("[Component reader]: chiral-center atom index {} out of range (molecule has {} atoms)\n", id,
                        numberOfBeads));
      }
    }
    std::array<std::size_t, 4> sorted = ids;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
    {
      throw std::runtime_error(std::format(
          "[Component reader]: the four atom indices of a chiral center must be distinct ({})\n", item.dump()));
    }

    // The handedness follows from the reference geometry: the signed volume of the tetrahedron
    // spanned by the three neighbors relative to the center. A (near-)planar center carries no
    // handedness and cannot be preserved during growth.
    double3 p0 = atoms[ids[0]].position;
    double3 d1 = atoms[ids[1]].position - p0;
    double3 d2 = atoms[ids[2]].position - p0;
    double3 d3 = atoms[ids[3]].position - p0;
    double signedVolume = double3::dot(d1, double3::cross(d2, d3));
    if (std::abs(signedVolume) < 1e-6)
    {
      throw std::runtime_error(std::format(
          "[Component reader]: chiral center ({}) is planar in the reference geometry, its parity is undefined\n",
          item.dump()));
    }

    ChiralCenter center{};
    center.type = signedVolume > 0.0 ? ChiralCenter::Chirality::R_Chiral : ChiralCenter::Chirality::S_Chiral;
    center.ids = ids;
    chiral_centers.push_back(center);
  }

  return chiral_centers;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Component &c)
{
  archive << c.versionNumber;

  archive << c.type;

  archive << c.name;
  archive << c.filenameData;
  archive << c.filename;

  archive << c.blockingPockets;
  archive << c.automaticBlockingPockets;

  archive << c.rigid;
  archive << c.translationalDegreesOfFreedom;
  archive << c.rotationalDegreesOfFreedom;

  archive << c.criticalTemperature;
  archive << c.criticalPressure;
  archive << c.acentricFactor;
  archive << c.molFraction;
  archive << c.swappable;
  archive << c.partialPressure;

  archive << c.totalMass;
  archive << c.fugacityCoefficient;
  archive << c.amountOfExcessMolecules;
  archive << c.bulkFluidDensity;
  archive << c.compressibility;

  archive << c.idealGasRosenbluthWeight;
  archive << c.idealGasEnergy;
  archive << c.crossSection;
  archive << c.liquidVolume;
  archive << c.saturationPressure;

  archive << c.netCharge;
  archive << c.startingBead;
  archive << c.pivotRandomizationFraction;
  archive << c.crankshaftRandomizationFraction;
  archive << c.crankshaftMaxSegmentSize;
  archive << c.repeatUnits;
  archive << c.endToEndAtoms;
  archive << c.definedAtoms;

  archive << c.inertiaVector;
  archive << c.inverseInertiaVector;
  archive << c.shapeType;
  archive << c.atoms;

  archive << c.connectivityTable;
  archive << c.intraMolecularPotentials;
  archive << c.fragmentGraph;
  archive << c.partialReinsertionFixedAtoms;
  archive << c.identityChanges;
  archive << c.gibbsIdentityChanges;
  archive << c.identitySwitches;

  archive << c.initialNumberOfMolecules;

  archive << c.lambdaGC;
  archive << c.lambdaGibbs;
  archive << c.hasFractionalMolecule;
  archive << c.fixedLambdaBin;
  archive << c.fixedLambdaCoordinate;

  archive << c.mc_moves_probabilities;
  archive << c.mc_moves_statistics;
  archive << c.mc_moves_cputime;

  archive << c.cbmcMoveStatistics;

  archive << c.averageRosenbluthWeights;
  archive << c.averageGibbsRosenbluthWeights;

  archive << c.lnPartitionFunction;

  archive << c.pairComponentId;
  archive << c.maximumPairDistance;

  archive << c.groupComponentIds;
  archive << c.maximumGroupDistance;

  archive << c.tetheredProtonHopSiteGroups;

  archive << c.lambdaPairSwap;
  archive << c.lambdaPairSwapCB;

  archive << c.lambdaGroupSwap;
  archive << c.lambdaGroupSwapCB;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Component &c)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  // The version-5 layout is incompatible with older archives (scratch/duplicated fields were
  // dropped), so a mismatched version cannot be read field-by-field and must be rejected outright.
  if (versionNumber != c.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format(
        "Invalid version {} reading 'Component' (expected {}) at line {} in file {}\n", versionNumber,
        c.versionNumber, location.line(), location.file_name()));
  }

  archive >> c.type;

  archive >> c.name;
  archive >> c.filenameData;
  archive >> c.filename;

  archive >> c.blockingPockets;
  archive >> c.automaticBlockingPockets;

  archive >> c.rigid;
  archive >> c.translationalDegreesOfFreedom;
  archive >> c.rotationalDegreesOfFreedom;

  archive >> c.criticalTemperature;
  archive >> c.criticalPressure;
  archive >> c.acentricFactor;
  archive >> c.molFraction;
  archive >> c.swappable;
  archive >> c.partialPressure;

  archive >> c.totalMass;
  archive >> c.fugacityCoefficient;
  archive >> c.amountOfExcessMolecules;
  archive >> c.bulkFluidDensity;
  archive >> c.compressibility;

  archive >> c.idealGasRosenbluthWeight;
  archive >> c.idealGasEnergy;
  archive >> c.crossSection;
  archive >> c.liquidVolume;
  archive >> c.saturationPressure;

  archive >> c.netCharge;
  archive >> c.startingBead;
  archive >> c.pivotRandomizationFraction;
  archive >> c.crankshaftRandomizationFraction;
  archive >> c.crankshaftMaxSegmentSize;
  archive >> c.repeatUnits;
  archive >> c.endToEndAtoms;
  archive >> c.definedAtoms;

  archive >> c.inertiaVector;
  archive >> c.inverseInertiaVector;
  archive >> c.shapeType;
  archive >> c.atoms;

  archive >> c.connectivityTable;
  archive >> c.intraMolecularPotentials;
  archive >> c.fragmentGraph;
  archive >> c.partialReinsertionFixedAtoms;
  archive >> c.identityChanges;
  archive >> c.gibbsIdentityChanges;
  // 'identitySwitches' was added in version 2; older restart files do not contain it.
  if (versionNumber >= 2)
  {
    archive >> c.identitySwitches;
  }

  archive >> c.initialNumberOfMolecules;

  archive >> c.lambdaGC;
  archive >> c.lambdaGibbs;
  archive >> c.hasFractionalMolecule;
  archive >> c.fixedLambdaBin;
  archive >> c.fixedLambdaCoordinate;

  archive >> c.mc_moves_probabilities;
  archive >> c.mc_moves_statistics;
  archive >> c.mc_moves_cputime;

  archive >> c.cbmcMoveStatistics;

  archive >> c.averageRosenbluthWeights;
  archive >> c.averageGibbsRosenbluthWeights;

  archive >> c.lnPartitionFunction;

  archive >> c.pairComponentId;
  archive >> c.maximumPairDistance;

  archive >> c.groupComponentIds;
  archive >> c.maximumGroupDistance;

  archive >> c.tetheredProtonHopSiteGroups;

  archive >> c.lambdaPairSwap;
  archive >> c.lambdaPairSwapCB;

  archive >> c.lambdaGroupSwap;
  archive >> c.lambdaGroupSwapCB;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Component: Error in binary restart\n"));
  }
#endif

  return archive;
}

std::string Component::repr() const { return std::string("Component test"); }

Component Component::makeMethane(const ForceField &forceField, std::size_t id)
{
  std::optional<std::size_t> type_ch4 = forceField.findPseudoAtom("CH4");
  if (!type_ch4.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "CH4"));
  }

  return Component(forceField, "methane", 190.564, 45599200, 0.01142, 
                   {Atom({0, 0, 0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(type_ch4.value()), static_cast<std::uint8_t>(id), false, false)},
                   {}, {}, 5, 21);
}

Component Component::makeCO2(const ForceField &forceField, std::size_t id, bool useCharges)
{
  const double qC = useCharges ? 0.6512 : 0.0;
  const double qO = useCharges ? -0.3256 : 0.0;

  std::optional<std::size_t> type_c_co2 = forceField.findPseudoAtom("C_co2");
  if (!type_c_co2.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "C_co2"));
  }

  std::optional<std::size_t> type_o_co2 = forceField.findPseudoAtom("O_co2");
  if (!type_o_co2.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "O_co2"));
  }

  return Component(
      forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {Atom({0, 0,  1.149},  qO, 1.0, 0, static_cast<std::uint16_t>(type_o_co2.value()), static_cast<std::uint8_t>(id), false, false), 
       Atom({0, 0,  0.000},  qC, 1.0, 0, static_cast<std::uint16_t>(type_c_co2.value()), static_cast<std::uint8_t>(id), false, false),
       Atom({0, 0, -1.149},  qO, 1.0, 0, static_cast<std::uint16_t>(type_o_co2.value()), static_cast<std::uint8_t>(id), false, false)},
      {}, {}, 5, 21);
}

Component Component::makeN2(const ForceField &forceField, std::size_t id, bool useCharges)
{
  std::optional<std::size_t> type_n_n2 = forceField.findPseudoAtom("N_n2");
  if (!type_n_n2.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "N_n2"));
  }

  std::optional<std::size_t> type_n_com = forceField.findPseudoAtom("N_com");
  if (!type_n_com.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "N_com"));
  }

  const double qN = useCharges ? forceField.pseudoAtoms[type_n_n2.value()].charge : 0.0;
  const double qCom = useCharges ? forceField.pseudoAtoms[type_n_com.value()].charge : 0.0;

  Component component(
      forceField, "N2", 126.192, 3395800.0, 0.0372,
      {Atom({0, 0,  0.55}, qN, 1.0, 0, static_cast<std::uint16_t>(type_n_n2.value()), static_cast<std::uint8_t>(id), false, false),
       Atom({0, 0,  0.00}, qCom, 1.0, 0, static_cast<std::uint16_t>(type_n_com.value()), static_cast<std::uint8_t>(id), false, false),
       Atom({0, 0, -0.55}, qN, 1.0, 0, static_cast<std::uint16_t>(type_n_n2.value()), static_cast<std::uint8_t>(id), false, false)},
      {}, {}, 5, 21);
  component.crossSection = 16.2;
  component.liquidVolume = 57.7;
  component.saturationPressure = 101325.0;
  return component;
}

Component Component::makeSphericalProbe(const ForceField &forceField, std::size_t id, std::string_view typeName)
{
  const std::string name(typeName);
  std::optional<std::size_t> type = forceField.findPseudoAtom(name);
  if (!type.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", name));
  }

  Component component(
      forceField, name, 126.192, 3395800.0, 0.0372,
      {Atom({0, 0, 0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(type.value()), static_cast<std::uint8_t>(id), false,
            false)},
      {}, {}, 5, 21);
  if (!(component.totalMass > 0.0))
  {
    component.totalMass = 28.0134;
    if (!component.definedAtoms.empty()) component.definedAtoms.front().second = 28.0134;
  }
  return component;
}

Component Component::makeWater(const ForceField &forceField, std::size_t id, bool useCharges)
{
  const double qh = useCharges ? 0.241 : 0.0;
  const double ql = useCharges ? -0.241 : 0.0;

  std::optional<std::size_t> type_ow = forceField.findPseudoAtom("Ow");
  if (!type_ow.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Ow"));
  }

  std::optional<std::size_t> type_hw = forceField.findPseudoAtom("Hw");
  if (!type_hw.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Hw"));
  }

  std::optional<std::size_t> type_lw = forceField.findPseudoAtom("Lw");
  if (!type_lw.has_value())
  {
    throw std::runtime_error(
        std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", "Lw"));
  }

  return Component(
      forceField, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, static_cast<std::uint16_t>(type_ow.value()), static_cast<std::uint8_t>(id), false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), qh, 1.0, 0, static_cast<std::uint16_t>(type_hw.value()), static_cast<std::uint8_t>(id), false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), qh, 1.0, 0, static_cast<std::uint16_t>(type_hw.value()), static_cast<std::uint8_t>(id), false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), ql, 1.0, 0, static_cast<std::uint16_t>(type_lw.value()), static_cast<std::uint8_t>(id), false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), ql, 1.0, 0, static_cast<std::uint16_t>(type_lw.value()), static_cast<std::uint8_t>(id), false, false)},
      {}, {}, 5, 21);
}

Component Component::makeIon(const ForceField &forceField, std::size_t id, std::string_view name, std::size_t type, double q)
{
  return Component(forceField, std::string{name}, 0.0, 0.0, 0.0,
                   {Atom({0, 0, 0}, q, 1.0, 0, static_cast<std::uint16_t>(type), static_cast<std::uint8_t>(id), false, false)}, {},
                   {}, 5, 21);
}
