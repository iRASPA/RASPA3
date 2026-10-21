module;

module dctst_free_energy;

import std;

import stringutils;
import hardware_info;
import system;
import framework;
import randomnumbers;
import input_reader;
import component;
import averages;
import units;
import simulationbox;
import forcefield;
import running_energy;
import atom;
import molecule;
import int3;
import double3;
import double3x3;
import simd_quatd;
import cbmc;
import mc_moves;
import mc_moves_move_types;
import mc_moves_translation;
import mc_moves_random_translation;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_partial_reinsertion;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import json;
import dctst_free_energy_settings;
import dctst_barrier;
import dctst_wham;

// =================================================================================================
// Accumulators
// =================================================================================================

DcTSTWindowAccumulator::DcTSTWindowAccumulator(std::size_t numberOfBlocks, std::size_t numberOfBins,
                                               double maximumTranslation)
    : numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      histogram(numberOfBlocks, std::vector<double>(numberOfBins, 0.0)),
      count(numberOfBlocks, 0.0),
      sumX(numberOfBlocks, 0.0),
      sumXSquared(numberOfBlocks, 0.0),
      sumBias(numberOfBlocks, 0.0),
      maximumTranslation(maximumTranslation),
      barrierOccupancyPairs((DcTSTPilotAccumulator::maximumOccupancy + 1uz) * (DcTSTPilotAccumulator::maximumOccupancy + 1uz),
                            0.0)
{
}

void DcTSTWindowAccumulator::optimize(double targetAcceptance)
{
  if (translationTrials > 0)
  {
    const double acceptance = static_cast<double>(translationAccepted) / static_cast<double>(translationTrials);
    const double scaling = std::clamp(acceptance / targetAcceptance, 0.5, 1.5);
    maximumTranslation = std::clamp(maximumTranslation * scaling, 0.01, 5.0);
  }
  translationTotalTrials += translationTrials;
  translationTotalAccepted += translationAccepted;
  translationTrials = 0;
  translationAccepted = 0;
  reinsertionTotalTrials += reinsertionTrials;
  reinsertionTotalAccepted += reinsertionAccepted;
  reinsertionTrials = 0;
  reinsertionAccepted = 0;
}

DcTSTPilotAccumulator::DcTSTPilotAccumulator(std::size_t numberOfBlocks, std::size_t numberOfBins,
                                             std::size_t numberOfConditionedBins)
    : numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      histogram(numberOfBlocks, std::vector<double>(numberOfBins, 0.0)),
      occupancyPairs(numberOfBlocks, std::vector<double>((maximumOccupancy + 1uz) * (maximumOccupancy + 1uz), 0.0)),
      conditionedHistogram(numberOfBlocks, std::vector<double>(numberOfConditionedBins, 0.0))
{
}

// =================================================================================================
// Helpers
// =================================================================================================

namespace
{
/// Mean and 2 sigma error of the mean of block values.
std::pair<double, double> meanAndError(std::span<const double> values)
{
  if (values.empty()) return {0.0, 0.0};
  const double n = static_cast<double>(values.size());
  const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
  if (values.size() < 2) return {mean, 0.0};
  double sumSquares{0.0};
  for (double v : values) sumSquares += (v - mean) * (v - mean);
  return {mean, 2.0 * std::sqrt(sumSquares / (n - 1.0) / n)};
}

/// A bias large enough to forbid a bin without producing infinities under -ffast-math.
constexpr double forbiddenBias = 1e30;

/// The inverse effective mass of the reaction bead along the unit normal n for the current orientation of a
/// molecule: the velocity of the bead along n is Gaussian with variance kT times this value. For a rigid
/// molecule v_bead.n = v_com.n + omega_body.(r_body x n_body) with independent Gaussian v_com (variance kT/M)
/// and omega_body (variances kT/I_i), so 1/m_eff = 1/M + sum_i (r_body x n_body)_i^2 / I_i; for a flexible
/// molecule the bead is an atom of mass m_bead.
double beadInverseEffectiveMass(const Component& component, const Molecule& molecule, std::size_t bead,
                                const double3& normal)
{
  if (!component.rigid) return 1.0 / component.definedAtoms[bead].second;
  double inverseMass = 1.0 / molecule.mass;
  if (molecule.numberOfAtoms > 1)
  {
    simd_quatd orientation = molecule.orientation;
    const double3 normalBody = orientation.inverse() * normal;
    const double3 lever = double3::cross(component.atoms[bead].position, normalBody);
    const double3 invI = component.inverseInertiaVector;
    inverseMass += lever.x * lever.x * invI.x + lever.y * lever.y * invI.y + lever.z * lever.z * invI.z;
  }
  return inverseMass;
}

DcTSTGrid makeGrid(double begin, double end, double requestedWidth, double barrier, double normalizationBegin,
                   double normalizationEnd)
{
  DcTSTGrid grid;
  grid.begin = begin;
  grid.end = end;
  grid.numberOfBins = std::max(2uz, static_cast<std::size_t>(std::llround((end - begin) / requestedWidth)));
  grid.width = (end - begin) / static_cast<double>(grid.numberOfBins);
  grid.barrier = barrier;
  grid.normalizationBegin = normalizationBegin;
  grid.normalizationEnd = normalizationEnd;
  return grid;
}

/// Weighted least-squares fit y = a + b t + c t^2; returns {a, b, c} and whether the system was regular.
struct QuadraticFit
{
  double a{0.0};
  double b{0.0};
  double c{0.0};
  bool valid{false};
  std::size_t points{0};
};

QuadraticFit fitQuadratic(std::span<const double> t, std::span<const double> y, std::span<const double> weight)
{
  QuadraticFit fit;
  fit.points = t.size();
  if (t.size() < 3) return fit;
  // normal equations: sum w t^(i+j) p_j = sum w t^i y
  std::array<double, 5> s{};
  std::array<double, 3> r{};
  for (std::size_t n = 0; n < t.size(); ++n)
  {
    double tp = 1.0;
    for (std::size_t i = 0; i < 5; ++i)
    {
      s[i] += weight[n] * tp;
      if (i < 3) r[i] += weight[n] * tp * y[n];
      tp *= t[n];
    }
  }
  const double m[3][3] = {{s[0], s[1], s[2]}, {s[1], s[2], s[3]}, {s[2], s[3], s[4]}};
  const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                     m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  if (std::abs(det) < 1e-300) return fit;
  auto solveColumn = [&](std::size_t column)
  {
    double mc[3][3];
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j) mc[i][j] = (j == column) ? r[i] : m[i][j];
    return (mc[0][0] * (mc[1][1] * mc[2][2] - mc[1][2] * mc[2][1]) - mc[0][1] * (mc[1][0] * mc[2][2] - mc[1][2] * mc[2][0]) +
            mc[0][2] * (mc[1][0] * mc[2][1] - mc[1][1] * mc[2][0])) /
           det;
  };
  fit.a = solveColumn(0);
  fit.b = solveColumn(1);
  fit.c = solveColumn(2);
  fit.valid = true;
  return fit;
}
}  // namespace

// =================================================================================================
// Construction and setup
// =================================================================================================

DynamicallyCorrectedTSTFreeEnergy::DynamicallyCorrectedTSTFreeEnergy(InputReader& reader)
    : random(reader.randomSeed),
      settings(reader.dctstFreeEnergySettings),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      numberOfBlocks(reader.numberOfBlocks),
      numberOfReplicas(std::max(1uz, reader.dctstFreeEnergySettings.numberOfReplicas))
{
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  occupancyMode = settings.cageOccupancies.has_value();

  barrier = DcTSTBarrier(templateSystem.simulationBox, settings.barrierPosition, settings.barrierNormal,
                         settings.latticeDistance);

  double3 lattice_origin = barrier.point + 0.5 * barrier.lambda * barrier.normal;
  if (settings.cageCenter.has_value())
  {
    if (!templateSystem.framework.has_value())
    {
      throw std::runtime_error("dcTST: 'CageCenter' requires a framework (it is a unit-cell fractional position)\n");
    }
    lattice_origin = templateSystem.framework->simulationBox.cell * settings.cageCenter.value();
    const double offset =
        std::abs(barrier.reducedCoordinate(barrier.projection(templateSystem.simulationBox, lattice_origin)));
    if (std::abs(offset - 0.5 * barrier.lambda) > 1e-3 * barrier.lambda)
    {
      throw std::runtime_error(
          std::format("dcTST: 'CageCenter' must lie half a 'LatticeDistance' from the barrier plane along the normal "
                      "(it lies {:.4f} [A] from the nearest plane, expected {:.4f})\n",
                      offset, 0.5 * barrier.lambda));
    }
  }
  if (templateSystem.framework.has_value())
  {
    cageLattice = DcTSTCageLattice(templateSystem.framework->simulationBox,
                                   templateSystem.framework->numberOfUnitCells, lattice_origin);
  }
  else
  {
    cageLattice = DcTSTCageLattice(templateSystem.simulationBox, int3(1, 1, 1), lattice_origin);
  }

  // grids: the period for the pilot; the period or the frame of cage A (extended a quarter cage into B so
  // that the barrier top can be fitted from both sides) for the reweighting
  const double lambda = barrier.lambda;
  pilotGrid = makeGrid(-0.5 * lambda, 0.5 * lambda, settings.fineBinWidth, 0.0, -0.5 * lambda, 0.5 * lambda);
  if (occupancyMode)
  {
    grid = makeGrid(-0.5 * lambda, 0.75 * lambda, settings.fineBinWidth, 0.5 * lambda, -0.5 * lambda, 0.5 * lambda);
  }
  else
  {
    grid = pilotGrid;
  }

  templateSystem.thermostat = std::nullopt;
  templateSystem.thermobarostat = std::nullopt;

  systems.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  randoms.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    randoms.emplace_back(random.seed + replicaId + 1);
  }

  pilot.assign(numberOfReplicas, DcTSTPilotAccumulator(numberOfBlocks, pilotGrid.numberOfBins, grid.numberOfBins));
  windowStates.assign(numberOfReplicas, DcTSTWindowState{});
  flatteningBias.assign(grid.numberOfBins, 0.0);
}

void DynamicallyCorrectedTSTFreeEnergy::run()
{
  setup();
  runInitialization();
  runPilot();
  buildWindows();
  runWindows();
  output();
}

void DynamicallyCorrectedTSTFreeEnergy::setup()
{
  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  std::filesystem::create_directories("dctst");

  const System& front = systems.front();
  stream.open(std::format("output/output_{}_{}.dctst_free_energy.txt", front.temperature, front.input_pressure),
              std::ios::out);
  outputJsonFileName =
      std::format("output/output_{}_{}.dctst_free_energy.json", front.temperature, front.input_pressure);

  std::print(stream, "{}", front.writeOutputHeader());
  std::print(stream, "Random seed: {}\n\n", random.seed);
  std::print(stream, "{}\n", HardwareInfo::writeInfo());
  std::print(stream, "{}", Units::printStatus());

  const Component& tagged = front.components[settings.taggedComponent];
  std::print(stream, "Dynamically corrected transition state theory (dcTST): free-energy profile\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Tagged component:                        {} ({})\n", tagged.name, settings.taggedComponent);
  std::print(stream, "Reaction bead (starting bead):           {}\n", tagged.startingBead);
  std::print(stream, "Barrier position (fractional):           {} {} {}\n", settings.barrierPosition.x,
             settings.barrierPosition.y, settings.barrierPosition.z);
  std::print(stream, "Barrier point (Cartesian):               {:.5f} {:.5f} {:.5f} [A]\n", barrier.point.x,
             barrier.point.y, barrier.point.z);
  std::print(stream, "Barrier normal (Cartesian, unit):        {:.5f} {:.5f} {:.5f}\n", barrier.normal.x,
             barrier.normal.y, barrier.normal.z);
  std::print(stream, "Lattice (cage-to-cage) distance lambda:  {} [A]\n", barrier.lambda);
  std::print(stream, "Cages along the normal in the box:       {:.3f}\n",
             double3::dot(double3::abs(barrier.normal), front.simulationBox.perpendicularWidths()) / barrier.lambda);
  {
    const double3 f = cageLattice.inverseCell * cageLattice.origin;
    const double3 reduced(f.x - std::floor(f.x), f.y - std::floor(f.y), f.z - std::floor(f.z));
    std::print(stream, "Cage center (unit-cell fractional):      {:.4f} {:.4f} {:.4f}{}\n", reduced.x, reduced.y,
               reduced.z,
               settings.cageCenter.has_value() ? " ('CageCenter')"
                                               : " (derived from 'BarrierPosition'; set 'CageCenter' to be explicit)");
  }
  std::print(stream, "Number of replicas (threads):            {}\n", numberOfReplicas);
  std::print(stream, "Number of blocks:                        {}\n", numberOfBlocks);
  std::print(stream, "Initialization cycles (unconstrained):   {}\n", numberOfInitializationCycles);
  std::print(stream, "Pilot cycles (unconstrained):            {} per replica\n", settings.numberOfPilotCycles);
  std::print(stream, "Window cycles:                           {} equilibration + {} production per window\n",
             settings.numberOfWindowEquilibrationCycles, settings.numberOfWindowCycles);
  std::print(stream, "Histogram sampling:                      every {} cycle(s)\n", settings.sampleEvery);
  std::print(stream, "Fine grid:                               {} bins of {:.5f} [A] over [{:.4f}, {:.4f})\n",
             grid.numberOfBins, grid.width, grid.begin, grid.end);
  std::print(stream, "Reported profile:                        {} bins\n", settings.numberOfFreeEnergyBins);
  std::print(stream, "Barrier fit half-range:                  {} [A]\n", settings.barrierFitRange);
  std::print(stream, "Flattening bias from the pilot:          {}\n", settings.flatteningBias ? "yes" : "no");
  std::print(stream, "Mirror symmetry q -> -q:                 {}\n", settings.mirrorSymmetry ? "yes" : "no");
  std::print(stream, "Tagged translation (initial maximum):    {} [A]\n", settings.maximumTranslationTagged);
  std::print(stream, "Tagged reinsertion probability:          {}\n", settings.reinsertionTaggedProbability);
  std::print(stream, "Tagged move fraction (window stage):     {}\n", settings.taggedMoveFraction);
  std::print(stream, "WHAM tolerance / maximum iterations:     {:.1e} / {}\n", settings.whamTolerance,
             settings.whamMaximumIterations);
  if (occupancyMode)
  {
    std::print(stream,
               "Occupancy-resolved mode:                 n_A = {} (tagged molecule included), n_B = {}; the profile is "
               "over cage A (x = 0 its center) extended a quarter cage into B (window to B at x = {:.4f})\n",
               settings.cageOccupancies->first, settings.cageOccupancies->second, 0.5 * barrier.lambda);
  }
  std::print(stream, "\n");
  std::flush(stream);

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["initialConditions"] = front.jsonSystemStatus();
  outputJson["initialization"]["components"] = front.jsonComponentStatus();
  outputJson["dcTST"]["taggedComponent"] = tagged.name;
  outputJson["dcTST"]["reactionBead"] = tagged.startingBead;
  outputJson["dcTST"]["barrierPosition"] = {settings.barrierPosition.x, settings.barrierPosition.y,
                                            settings.barrierPosition.z};
  outputJson["dcTST"]["barrierNormal"] = {barrier.normal.x, barrier.normal.y, barrier.normal.z};
  outputJson["dcTST"]["latticeDistance"] = barrier.lambda;
  outputJson["dcTST"]["numberOfReplicas"] = numberOfReplicas;
  outputJson["dcTST"]["temperature"] = front.temperature;
  outputJson["dcTST"]["numberOfMolecules"] = front.numberOfMolecules();
  if (occupancyMode)
  {
    outputJson["dcTST"]["cageOccupancies"] = {settings.cageOccupancies->first, settings.cageOccupancies->second};
  }

  replicaStreams.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(std::format("output/output_{}_{}.dctst_free_energy.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                std::ios::out);
    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
    std::print(replicaStream, "{}", system.writeOutputHeader());
    std::print(replicaStream, "dcTST free energy: replica {} of {}\n", replicaId, numberOfReplicas);
    std::print(replicaStream, "Random seed of this replica: {}\n\n", randoms[replicaId].seed);
    std::print(replicaStream, "{}", system.writeSystemStatus());
    std::print(replicaStream, "{}", system.forceField.printPseudoAtomStatus());
    std::print(replicaStream, "{}", system.forceField.printForceFieldStatus());
    std::print(replicaStream, "{}", system.writeComponentStatus());
    std::print(replicaStream, "{}", system.writeNumberOfPseudoAtoms());
    std::flush(replicaStream);
  }

  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t replicaId = 1; replicaId < systems.size(); ++replicaId)
  {
    systems[replicaId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[replicaId].interpolationGrids = systems.front().interpolationGrids;
  }
}

std::size_t DynamicallyCorrectedTSTFreeEnergy::taggedComponentOffset(const System& system) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < settings.taggedComponent; ++i)
  {
    index += system.numberOfMoleculesPerComponent[i];
  }
  return index;
}

// =================================================================================================
// Coordinates and biases
// =================================================================================================

double DynamicallyCorrectedTSTFreeEnergy::coordinate(const System& system, const DcTSTWindowState& state,
                                                     double3 position) const
{
  if (occupancyMode)
  {
    return double3::dot(system.simulationBox.applyPeriodicBoundaryConditions(position - cageLattice.center(state.cageA)),
                        barrier.normal);
  }
  return barrier.reducedCoordinate(barrier.projection(system.simulationBox, position));
}

double DynamicallyCorrectedTSTFreeEnergy::foldedCoordinate(double x) const
{
  return (settings.mirrorSymmetry && !occupancyMode) ? std::abs(x) : x;
}

double DynamicallyCorrectedTSTFreeEnergy::windowBias(std::size_t windowIndex, double x) const
{
  if (!grid.contains(x)) return forbiddenBias;
  return windows[windowIndex].bias[grid.bin(x)];
}

std::vector<std::size_t> DynamicallyCorrectedTSTFreeEnergy::countOccupancies(const System& system) const
{
  std::vector<std::size_t> occupancy(cageLattice.numberOfCages(), 0uz);
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  for (const Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    ++occupancy[cageLattice.cageOf(atoms[molecule.atomIndex + component.startingBead].position)];
  }
  return occupancy;
}

std::optional<double3> DynamicallyCorrectedTSTFreeEnergy::pointInCageFrame(RandomNumber& rng, const System& system,
                                                                           const DcTSTWindowState& state,
                                                                           double x) const
{
  // a random point of the box shifted along the normal to coordinate x in the frame of cage A; it is kept
  // when it lies laterally within cage A (the cage of its projection onto the plane through the center of A)
  const double3 centerA = cageLattice.center(state.cageA);
  const double3 r = system.simulationBox.randomPosition(rng);
  const double q = double3::dot(system.simulationBox.applyPeriodicBoundaryConditions(r - centerA), barrier.normal);
  const double3 point = r + (x - q) * barrier.normal;
  if (cageLattice.cageOf(point - x * barrier.normal) != state.cageA) return std::nullopt;
  return point;
}

// =================================================================================================
// Stages
// =================================================================================================

void DynamicallyCorrectedTSTFreeEnergy::runInitialization()
{
  const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  std::print(stream, "\nInitialization stage: {} MC cycles on {} replicas/threads\n", numberOfInitializationCycles,
             numberOfReplicas);
  std::flush(stream);
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId]()
          {
            System& system = systems[replicaId];
            std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
            system.precomputeTotalRigidEnergy();
            system.runningEnergies = system.computeTotalEnergies();
            std::print(replicaStream, "\nInitialization stage\n{}\n\n",
                       system.runningEnergies.printMC("Recomputed from scratch"));
            for (std::size_t cycle = 0; cycle != numberOfInitializationCycles; ++cycle)
            {
              performUnconstrainedCycle(replicaId, SimulationStage::Initialization);
              if (cycle % optimizeMCMovesEvery == 0uz) system.optimizeMCMoves();
              if (cycle % printEvery == 0uz)
              {
                std::print(replicaStream, "Cycle {} of {} (initialization)\n{}\n", cycle, numberOfInitializationCycles,
                           system.runningEnergies.printMC(""));
                std::flush(replicaStream);
              }
            }
          });
    }
  }
  const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalInitializationTime += (t2 - t1);
  totalSimulationTime += (t2 - t1);
}

void DynamicallyCorrectedTSTFreeEnergy::runPilot()
{
  const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  std::print(stream, "\nPilot stage: {} MC cycles on {} replicas/threads\n", settings.numberOfPilotCycles,
             numberOfReplicas);
  std::flush(stream);
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId]()
          {
            System& system = systems[replicaId];
            std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
            system.mc_moves_statistics.clearMoveStatistics();
            system.mc_moves_cputime.clearTimingStatistics();
            for (Component& component : system.components)
            {
              component.mc_moves_statistics.clearMoveStatistics();
              component.mc_moves_cputime.clearTimingStatistics();
            }
            system.precomputeTotalRigidEnergy();
            system.runningEnergies = system.computeTotalEnergies();
            std::print(replicaStream, "\nPilot stage\n{}\n\n", system.runningEnergies.printMC("Recomputed from scratch"));

            const std::size_t numberOfCycles = settings.numberOfPilotCycles;
            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfCycles));
            for (std::size_t cycle = 0; cycle != numberOfCycles; ++cycle)
            {
              estimation.setCurrentSample(cycle);
              performUnconstrainedCycle(replicaId, SimulationStage::Pilot);
              if (cycle % settings.sampleEvery == 0uz) samplePilot(replicaId, estimation.currentBin);
              if (cycle % optimizeMCMovesEvery == 0uz) system.optimizeMCMoves();
              if (cycle % printEvery == 0uz)
              {
                std::print(replicaStream, "Cycle {} of {} (pilot)\n{}\n", cycle, numberOfCycles,
                           system.runningEnergies.printMC(""));
                std::flush(replicaStream);
              }
            }
            RunningEnergy recomputed = system.computeTotalEnergies();
            std::print(replicaStream, "\n{}\n{}\n", recomputed.printMCDiff(system.runningEnergies),
                       system.writeMCMoveStatistics());
            std::flush(replicaStream);
          });
    }
  }
  const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalPilotTime += (t2 - t1);
  totalSimulationTime += (t2 - t1);
}

void DynamicallyCorrectedTSTFreeEnergy::buildWindows()
{
  const double lambda = barrier.lambda;
  const bool mirror = settings.mirrorSymmetry && !occupancyMode;
  const bool periodic = !occupancyMode && !mirror && !settings.windowRange.has_value();

  // --- flattening bias from the pilot histogram ---------------------------------------------------
  if (settings.flatteningBias && settings.numberOfPilotCycles > 0uz)
  {
    const std::vector<double> histogram = pooledPilotHistogram();
    const double total = std::accumulate(histogram.begin(), histogram.end(), 0.0);
    if (total > 0.0)
    {
      std::vector<double> f(grid.numberOfBins, 0.0);
      std::vector<bool> defined(grid.numberOfBins, false);
      for (std::size_t b = 0; b < grid.numberOfBins; ++b)
      {
        if (histogram[b] > 0.0)
        {
          f[b] = -std::log(histogram[b] / total);
          defined[b] = true;
        }
      }
      // running average over +-0.1 A (periodic in the unconditioned mode)
      const std::size_t half = std::max(1uz, static_cast<std::size_t>(std::llround(0.1 / grid.width)));
      std::vector<double> smooth(grid.numberOfBins, 0.0);
      std::vector<bool> smoothDefined(grid.numberOfBins, false);
      double maximum = -std::numeric_limits<double>::max();
      for (std::size_t b = 0; b < grid.numberOfBins; ++b)
      {
        double sum{0.0};
        std::size_t n{0};
        for (std::int64_t d = -static_cast<std::int64_t>(half); d <= static_cast<std::int64_t>(half); ++d)
        {
          std::int64_t i = static_cast<std::int64_t>(b) + d;
          if (!occupancyMode)
          {
            const std::int64_t nb = static_cast<std::int64_t>(grid.numberOfBins);
            i = ((i % nb) + nb) % nb;
          }
          else if (i < 0 || i >= static_cast<std::int64_t>(grid.numberOfBins))
          {
            continue;
          }
          if (defined[static_cast<std::size_t>(i)])
          {
            sum += f[static_cast<std::size_t>(i)];
            ++n;
          }
        }
        if (n > 0)
        {
          smooth[b] = sum / static_cast<double>(n);
          smoothDefined[b] = true;
          maximum = std::max(maximum, smooth[b]);
        }
      }
      double minimum = std::numeric_limits<double>::max();
      for (std::size_t b = 0; b < grid.numberOfBins; ++b)
      {
        flatteningBias[b] = -(smoothDefined[b] ? smooth[b] : maximum);
        minimum = std::min(minimum, flatteningBias[b]);
      }
      for (double& w : flatteningBias) w -= minimum;
    }
  }

  // --- window centers -----------------------------------------------------------------------------
  std::vector<double> centers = settings.windowCenters;
  double spacing = settings.windowSpacing;
  if (centers.empty())
  {
    double begin = -0.5 * lambda;
    double end = 0.5 * lambda;
    if (occupancyMode)
    {
      begin = -0.5 * lambda;
      end = 0.5 * lambda + 0.125 * lambda;
    }
    if (mirror)
    {
      begin = 0.0;
      end = 0.5 * lambda;
    }
    if (settings.windowRange.has_value())
    {
      begin = settings.windowRange->first;
      end = settings.windowRange->second;
    }
    const double length = end - begin;
    std::size_t K = settings.numberOfWindows;
    if (K == 0uz) K = std::max(2uz, static_cast<std::size_t>(std::ceil(length / settings.windowSpacing)));
    if (periodic)
    {
      // K windows over the period, one of them on the barrier
      spacing = lambda / static_cast<double>(K);
      for (std::size_t k = 0; k < K; ++k)
      {
        centers.push_back(barrier.reducedCoordinate(static_cast<double>(k) * spacing));
      }
    }
    else
    {
      K = std::max(2uz, K);
      spacing = length / static_cast<double>(K - 1uz);
      for (std::size_t k = 0; k < K; ++k)
      {
        centers.push_back(begin + static_cast<double>(k) * spacing);
      }
    }
  }
  else if (centers.size() > 1)
  {
    std::vector<double> sorted = centers;
    std::sort(sorted.begin(), sorted.end());
    spacing = sorted[1] - sorted[0];
    for (std::size_t k = 1; k + 1 < sorted.size(); ++k) spacing = std::min(spacing, sorted[k + 1] - sorted[k]);
  }
  const double stiffness = (settings.windowStiffness > 0.0) ? settings.windowStiffness : 4.0 / (spacing * spacing);

  windows.clear();
  for (std::size_t k = 0; k < centers.size(); ++k)
  {
    DcTSTWindow window;
    window.center = centers[k];
    window.stiffness = stiffness;
    window.supportHalfWidth = 4.0 / std::sqrt(stiffness);
    window.replica = k % numberOfReplicas;
    window.source = k;
    window.bias.resize(grid.numberOfBins);
    for (std::size_t b = 0; b < grid.numberOfBins; ++b)
    {
      double delta = foldedCoordinate(grid.center(b)) - window.center;
      if (!occupancyMode && !mirror) delta = barrier.reducedCoordinate(delta);
      window.bias[b] = 0.5 * stiffness * delta * delta + flatteningBias[b];
    }
    windows.push_back(std::move(window));
  }
  // more replicas than centers: repeat the centers round-robin so that no replica idles; a repeated window is
  // an independent sampling of the same state (its own accumulator, a separate state of the reweighting with
  // the identical bias, which the reweighting handles exactly)
  const std::size_t numberOfCenters = windows.size();
  for (std::size_t r = numberOfCenters; r < numberOfReplicas && numberOfCenters > 0uz; ++r)
  {
    DcTSTWindow window = windows[r % numberOfCenters];
    window.replica = r;
    window.source = r % numberOfCenters;
    windows.push_back(std::move(window));
  }
  windowData.assign(windows.size(),
                    DcTSTWindowAccumulator(numberOfBlocks, grid.numberOfBins, settings.maximumTranslationTagged));

  std::print(stream, "\nWindows: {} windows, stiffness beta kappa = {:.4f} [1/A^2] (sigma = {:.4f} [A]), spacing {:.4f} [A]\n",
             windows.size(), stiffness, 1.0 / std::sqrt(stiffness), spacing);
  for (std::size_t k = 0; k < windows.size(); ++k)
  {
    std::print(stream, "    window {:3}: center {:9.4f} [A], replica {}{}\n", k, windows[k].center, windows[k].replica,
               windows[k].source == k ? "" : std::format(" (repeat of window {})", windows[k].source));
  }
  std::flush(stream);
}

void DynamicallyCorrectedTSTFreeEnergy::runWindows()
{
  const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  std::print(stream, "\nWindow stage: {} windows on {} replicas/threads, {} + {} cycles per window\n", windows.size(),
             numberOfReplicas, settings.numberOfWindowEquilibrationCycles, settings.numberOfWindowCycles);
  std::flush(stream);
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId]()
          {
            System& system = systems[replicaId];
            std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
            const Component& component = system.components[settings.taggedComponent];

            for (std::size_t windowIndex = 0; windowIndex < windows.size(); ++windowIndex)
            {
              if (windows[windowIndex].replica != replicaId) continue;
              DcTSTWindowAccumulator& accumulator = windowData[windowIndex];
              const DcTSTWindowState& state = windowStates[replicaId];

              startWindow(replicaId, windowIndex);
              system.mc_moves_statistics.clearMoveStatistics();
              for (Component& c : system.components) c.mc_moves_statistics.clearMoveStatistics();
              system.precomputeTotalRigidEnergy();
              system.runningEnergies = system.computeTotalEnergies();
              std::print(replicaStream, "\nWindow {} (center {:.4f} [A])\n{}\n\n", windowIndex,
                         windows[windowIndex].center, system.runningEnergies.printMC("Recomputed from scratch"));

              const std::size_t equilibration = settings.numberOfWindowEquilibrationCycles;
              const std::size_t production = settings.numberOfWindowCycles;
              const std::size_t total = equilibration + production;
              BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, production));
              for (std::size_t cycle = 0; cycle != total; ++cycle)
              {
                const bool inProduction = cycle >= equilibration;
                if (inProduction) estimation.setCurrentSample(cycle - equilibration);
                performWindowCycle(replicaId, windowIndex);
                if (inProduction && (cycle - equilibration) % settings.sampleEvery == 0uz)
                {
                  sampleWindow(replicaId, windowIndex, estimation.currentBin);
                }
                if (cycle % optimizeMCMovesEvery == 0uz)
                {
                  system.optimizeMCMoves();
                  accumulator.optimize();
                }
                if (cycle % printEvery == 0uz)
                {
                  const std::span<const Atom> atoms =
                      system.spanOfMolecule(settings.taggedComponent, state.taggedLocalIndex);
                  const double x = coordinate(system, state, atoms[component.startingBead].position);
                  const std::size_t tt = accumulator.translationTotalTrials + accumulator.translationTrials;
                  const std::size_t ta = accumulator.translationTotalAccepted + accumulator.translationAccepted;
                  const std::size_t rt = accumulator.reinsertionTotalTrials + accumulator.reinsertionTrials;
                  const std::size_t ra = accumulator.reinsertionTotalAccepted + accumulator.reinsertionAccepted;
                  std::print(replicaStream, "Window {} cycle {} of {} ({}): x = {:.4f} [A], bias {:.3f} [kT]\n",
                             windowIndex, cycle, total, inProduction ? "production" : "equilibration", x,
                             windowBias(windowIndex, x));
                  std::print(replicaStream, "{}", system.runningEnergies.printMC(""));
                  std::print(replicaStream,
                             "    tagged translation: {}/{} accepted ({:.2f}%), maximum displacement {:.4f} [A]; "
                             "reinsertion: {}/{} accepted ({:.2f}%)\n",
                             ta, tt, 100.0 * static_cast<double>(ta) / static_cast<double>(std::max(1uz, tt)),
                             accumulator.maximumTranslation, ra, rt,
                             100.0 * static_cast<double>(ra) / static_cast<double>(std::max(1uz, rt)));
                  if (occupancyMode)
                  {
                    std::print(replicaStream,
                               "    occupancy constraint: cage A {} (n = {}), cage B {} (n = {}); {} cage changes, {} "
                               "reverted, {} tagged moves out of range\n",
                               state.cageA, state.occupancy[state.cageA] + 1uz, state.cageB,
                               state.occupancy[state.cageB], accumulator.cageChanges, accumulator.rejectedByOccupancy,
                               accumulator.rejectedByRange);
                  }
                  std::print(replicaStream, "\n");
                  std::flush(replicaStream);
                }
              }
              {
                const std::span<const Atom> atoms =
                    system.spanOfMolecule(settings.taggedComponent, state.taggedLocalIndex);
                accumulator.endCoordinate = coordinate(system, state, atoms[component.startingBead].position);
              }
              RunningEnergy recomputed = system.computeTotalEnergies();
              std::print(replicaStream, "\n{}\n", recomputed.printMCDiff(system.runningEnergies));
              std::flush(replicaStream);
            }
          });
    }
  }
  const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalWindowTime += (t2 - t1);
  totalSimulationTime += (t2 - t1);
}

// =================================================================================================
// Monte Carlo
// =================================================================================================

void DynamicallyCorrectedTSTFreeEnergy::performUnconstrainedCycle(std::size_t replicaId, SimulationStage stage)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  std::size_t fractionalMoleculeSystem = 0uz;
  const std::size_t numberOfSteps = std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();
  for (std::size_t step = 0; step != numberOfSteps; ++step)
  {
    const std::size_t selectedComponent = system.randomComponent(rng);
    if (stage == SimulationStage::Initialization)
    {
      MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
    }
    else
    {
      MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem, 0uz);
    }
  }
}

void DynamicallyCorrectedTSTFreeEnergy::startWindow(std::size_t replicaId, std::size_t windowIndex)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTWindowState& state = windowStates[replicaId];
  DcTSTWindowAccumulator& accumulator = windowData[windowIndex];
  std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
  const DcTSTWindow& window = windows[windowIndex];

  const std::size_t selectedComponent = settings.taggedComponent;
  const Component& component = system.components[selectedComponent];
  const std::size_t numberOfTagged = system.numberOfMoleculesPerComponent[selectedComponent];
  if (numberOfTagged == 0uz)
  {
    throw std::runtime_error(
        std::format("dcTST: the tagged component '{}' has no molecules (replica {})\n", component.name, replicaId));
  }

  // the molecule of the tagged component nearest to the window center. In the occupancy-resolved mode the
  // coordinate is taken relative to the molecule's own cage (in [-lambda/2, lambda/2)); a window beyond the
  // barrier (center >= lambda/2, in cage B) wants a molecule at center - lambda relative to its own cage, which
  // then becomes cage B and the cage at -normal of it cage A.
  const bool windowInCageB = occupancyMode && window.center >= 0.5 * barrier.lambda;
  std::size_t best{0};
  double bestDistance = std::numeric_limits<double>::max();
  for (std::size_t m = 0; m < numberOfTagged; ++m)
  {
    const std::span<const Atom> atoms = system.spanOfMolecule(selectedComponent, m);
    const double3 position = atoms[component.startingBead].position;
    double x, target;
    if (occupancyMode)
    {
      const std::size_t ownCage = cageLattice.cageOf(position);
      x = double3::dot(system.simulationBox.applyPeriodicBoundaryConditions(position - cageLattice.center(ownCage)),
                       barrier.normal);
      target = windowInCageB ? window.center - barrier.lambda : window.center;
    }
    else
    {
      x = foldedCoordinate(barrier.reducedCoordinate(barrier.projection(system.simulationBox, position)));
      target = window.center;
    }
    double delta = x - target;
    if (!occupancyMode && !(settings.mirrorSymmetry)) delta = barrier.reducedCoordinate(delta);
    if (std::abs(delta) < bestDistance)
    {
      bestDistance = std::abs(delta);
      best = m;
    }
  }
  state.taggedLocalIndex = best;

  if (!occupancyMode)
  {
    const std::span<const Atom> atoms = system.spanOfMolecule(selectedComponent, best);
    accumulator.startCoordinate = coordinate(system, state, atoms[component.startingBead].position);
    std::print(replicaStream, "Window {}: tagged molecule {} of component '{}' at x = {:.4f} [A] ({:.4f} from the center)\n",
               windowIndex, best, component.name, accumulator.startCoordinate, bestDistance);
    std::flush(replicaStream);
    return;
  }

  // occupancy-resolved mode: cage A is the cage of the tagged molecule, B its neighbor along +normal; move
  // other molecules in or out until the requested occupancies hold
  const auto [nA, nB] = settings.cageOccupancies.value();
  const std::size_t taggedIndex = taggedComponentOffset(system) + best;
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t numberOfCages = cageLattice.numberOfCages();
  if (nA - 1uz + nB > numberOfMolecules - 1uz)
  {
    throw std::runtime_error(
        std::format("dcTST: 'CageOccupancies' [{}, {}] needs {} molecules besides the tagged one, the system has {}\n",
                    nA, nB, nA - 1uz + nB, numberOfMolecules - 1uz));
  }
  {
    const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    state.cageOfMolecule.assign(numberOfMolecules, 0uz);
    state.occupancy.assign(numberOfCages, 0uz);
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      const Molecule& molecule = system.moleculeData[j];
      const Component& c = system.components[molecule.componentId];
      state.cageOfMolecule[j] = cageLattice.cageOf(atoms[molecule.atomIndex + c.startingBead].position);
      if (j != taggedIndex) ++state.occupancy[state.cageOfMolecule[j]];
    }
    if (windowInCageB)
    {
      state.cageB = state.cageOfMolecule[taggedIndex];
      state.cageA = cageLattice.neighbor(state.cageB, -barrier.lambda * barrier.normal);
    }
    else
    {
      state.cageA = state.cageOfMolecule[taggedIndex];
      state.cageB = cageLattice.neighbor(state.cageA, barrier.lambda * barrier.normal);
    }
  }
  if (state.cageA == state.cageB)
  {
    throw std::runtime_error(
        "dcTST: cages A and B coincide; the occupancy-resolved mode needs at least two unit cells along the barrier "
        "normal\n");
  }
  const bool hasOtherCages = numberOfCages > 2uz;

  const CBMC::GrowContext context = system.makeGrowContext(CBMC::CutOffMode::Full);
  auto relocate = [&](std::size_t j, std::size_t targetCage)
  {
    Molecule& molecule = system.moleculeData[j];
    const Component& c = system.components[molecule.componentId];
    std::size_t offset{0};
    for (std::size_t k = 0; k < molecule.componentId; ++k) offset += system.numberOfMoleculesPerComponent[k];
    std::span<Atom> moleculeAtoms = system.spanOfMolecule(molecule.componentId, j - offset);
    for (std::size_t attempt = 0; attempt < 100000uz; ++attempt)
    {
      const double3 point = cageLattice.randomPointInCage(rng, targetCage);
      std::optional<CBMC::GrowResult> growData = CBMC::growNewMolecule(
          rng, context, c,
          CBMC::NewMoleculeIdentity{.componentId = molecule.componentId,
                                    .moleculeId = static_cast<std::size_t>(moleculeAtoms.front().moleculeId)},
          CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = point});
      if (!growData) continue;
      if (system.insideBlockedPockets(c, growData->atoms)) continue;
      if (cageLattice.cageOf(growData->atoms[c.startingBead].position) != targetCage) continue;
      for (std::size_t i = 0; i < moleculeAtoms.size(); ++i) moleculeAtoms[i].position = growData->atoms[i].position;
      Molecule record = growData->molecule;
      record.atomIndex = molecule.atomIndex;
      record.numberOfAtoms = molecule.numberOfAtoms;
      molecule = record;
      system.invalidateCellList();
      --state.occupancy[state.cageOfMolecule[j]];
      ++state.occupancy[targetCage];
      state.cageOfMolecule[j] = targetCage;
      return;
    }
    throw std::runtime_error(
        std::format("dcTST: could not relocate a molecule into cage {} (replica {})\n", targetCage, replicaId));
  };
  auto randomMoleculeWhere = [&](auto&& predicate) -> std::optional<std::size_t>
  {
    std::vector<std::size_t> candidates;
    for (std::size_t j = 0; j < numberOfMolecules; ++j)
    {
      if (j != taggedIndex && predicate(state.cageOfMolecule[j])) candidates.push_back(j);
    }
    if (candidates.empty()) return std::nullopt;
    return candidates[std::min(candidates.size() - 1uz,
                               static_cast<std::size_t>(rng.uniform() * static_cast<double>(candidates.size())))];
  };
  auto randomOtherCage = [&]()
  {
    for (;;)
    {
      const std::size_t cage =
          std::min(numberOfCages - 1uz, static_cast<std::size_t>(rng.uniform() * static_cast<double>(numberOfCages)));
      if (cage != state.cageA && cage != state.cageB) return cage;
    }
  };
  const std::size_t initialA = state.occupancy[state.cageA];
  const std::size_t initialB = state.occupancy[state.cageB];
  std::size_t relocations{0};
  for (const auto [cage, target] : {std::pair{state.cageA, nA - 1uz}, std::pair{state.cageB, nB}})
  {
    while (state.occupancy[cage] > target)
    {
      if (!hasOtherCages) throw std::runtime_error("dcTST: no cage besides A and B to move the excess molecules to\n");
      relocate(randomMoleculeWhere([&](std::size_t c) { return c == cage; }).value(), randomOtherCage());
      ++relocations;
    }
  }
  for (const auto [cage, target] : {std::pair{state.cageA, nA - 1uz}, std::pair{state.cageB, nB}})
  {
    while (state.occupancy[cage] < target)
    {
      const std::optional<std::size_t> j =
          randomMoleculeWhere([&](std::size_t c) { return c != state.cageA && c != state.cageB; });
      if (!j.has_value()) throw std::runtime_error("dcTST: not enough molecules outside cages A and B to fill them\n");
      relocate(j.value(), cage);
      ++relocations;
    }
  }
  {
    const std::span<const Atom> atoms = system.spanOfMolecule(selectedComponent, best);
    accumulator.startCoordinate = coordinate(system, state, atoms[component.startingBead].position);
  }
  std::print(replicaStream,
             "Window {}: tagged molecule {} at x = {:.4f} [A]; occupancy constraint: cage A = {} (had {} other "
             "molecules, now {}), cage B = {} (had {}, now {}); {} molecules relocated\n",
             windowIndex, best, accumulator.startCoordinate, state.cageA, initialA, state.occupancy[state.cageA],
             state.cageB, initialB, state.occupancy[state.cageB], relocations);
  std::flush(replicaStream);
}

void DynamicallyCorrectedTSTFreeEnergy::performWindowCycle(std::size_t replicaId, std::size_t windowIndex)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTWindowState& state = windowStates[replicaId];
  DcTSTWindowAccumulator& accumulator = windowData[windowIndex];

  const std::size_t numberOfSteps = std::max(system.numberOfMolecules(), 20uz) * system.numerOfAdsorbateComponents();
  for (std::size_t step = 0; step != numberOfSteps; ++step)
  {
    // a fixed fraction of the steps moves the tagged molecule, the rest a random molecule (which may again be
    // the tagged one); the selection probabilities are configuration independent so detailed balance holds
    std::size_t selectedComponent, selectedMolecule;
    if (rng.uniform() < settings.taggedMoveFraction)
    {
      selectedComponent = settings.taggedComponent;
      selectedMolecule = state.taggedLocalIndex;
    }
    else
    {
      selectedComponent = system.randomComponent(rng);
      const std::size_t numberOfMolecules = system.numberOfMoleculesPerComponent[selectedComponent];
      if (numberOfMolecules == 0uz) continue;
      selectedMolecule = std::min(numberOfMolecules - 1uz,
                                  static_cast<std::size_t>(rng.uniform() * static_cast<double>(numberOfMolecules)));
    }

    std::optional<RunningEnergy> energyDifference;
    if (selectedComponent == settings.taggedComponent && selectedMolecule == state.taggedLocalIndex)
    {
      energyDifference = (rng.uniform() < settings.reinsertionTaggedProbability)
                             ? biasedReinsertionMove(replicaId, windowIndex)
                             : biasedTranslationMove(replicaId, windowIndex);
    }
    else
    {
      Component& component = system.components[selectedComponent];
      const Move::Types moveType = component.mc_moves_probabilities.sample(rng);
      component.mc_moves_statistics.addAllCounts(moveType);

      const std::size_t moleculeIndex =
          occupancyMode ? system.moleculeIndexOfComponent(selectedComponent, selectedMolecule) : 0uz;
      std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
      std::vector<Atom> oldAtoms;
      Molecule oldMolecule;
      if (occupancyMode)
      {
        oldAtoms.assign(moleculeAtoms.begin(), moleculeAtoms.end());
        oldMolecule = system.moleculeData[moleculeIndex];
      }

      switch (moveType)
      {
        case Move::Types::RandomTranslation:
          energyDifference = MC_Moves::randomTranslationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::Rotation:
          energyDifference = MC_Moves::rotationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::RandomRotation:
          energyDifference = MC_Moves::randomRotationMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::ReinsertionCBMC:
          energyDifference = MC_Moves::reinsertionMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::PartialReinsertionCBMC:
          energyDifference = MC_Moves::partialReinsertionMove(rng, system, selectedComponent, selectedMolecule);
          break;
        case Move::Types::Translation:
        default:
          energyDifference = MC_Moves::translationMove(rng, system, selectedComponent, selectedMolecule);
          break;
      }

      if (occupancyMode && energyDifference.has_value())
      {
        const std::size_t oldCage = state.cageOfMolecule[moleculeIndex];
        const std::size_t newCage = cageLattice.cageOf(moleculeAtoms[component.startingBead].position);
        if (newCage != oldCage)
        {
          ++accumulator.cageChanges;
          if (oldCage == state.cageA || oldCage == state.cageB || newCage == state.cageA || newCage == state.cageB)
          {
            if (system.forceField.useCharge && system.forceField.usesEwaldFourier())
            {
              Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, system.components, oldAtoms, moleculeAtoms);
              Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
            }
            std::copy(oldAtoms.cbegin(), oldAtoms.cend(), moleculeAtoms.begin());
            system.cellListAtomsMoved(moleculeAtoms);
            system.moleculeData[moleculeIndex] = oldMolecule;
            energyDifference = std::nullopt;
            ++accumulator.rejectedByOccupancy;
          }
          else
          {
            --state.occupancy[oldCage];
            ++state.occupancy[newCage];
            state.cageOfMolecule[moleculeIndex] = newCage;
          }
        }
      }
    }
    if (energyDifference)
    {
      system.runningEnergies += energyDifference.value();
    }
  }
}

std::optional<RunningEnergy> DynamicallyCorrectedTSTFreeEnergy::biasedTranslationMove(std::size_t replicaId,
                                                                                      std::size_t windowIndex)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTWindowState& state = windowStates[replicaId];
  DcTSTWindowAccumulator& statistics = windowData[windowIndex];

  const std::size_t selectedComponent = settings.taggedComponent;
  Component& component = system.components[selectedComponent];
  std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, state.taggedLocalIndex);
  Molecule& molecule = system.moleculeData[taggedComponentOffset(system) + state.taggedLocalIndex];

  ++statistics.translationTrials;

  const double3 displacement(statistics.maximumTranslation * 2.0 * (rng.uniform() - 0.5),
                             statistics.maximumTranslation * 2.0 * (rng.uniform() - 0.5),
                             statistics.maximumTranslation * 2.0 * (rng.uniform() - 0.5));
  std::pair<Molecule, std::vector<Atom>> trialMolecule = component.translate(molecule, moleculeAtoms, displacement);

  if (system.insideBlockedPockets(component, trialMolecule.second)) return std::nullopt;

  const double xOld = coordinate(system, state, moleculeAtoms[component.startingBead].position);
  const double xNew = coordinate(system, state, trialMolecule.second[component.startingBead].position);
  if (!grid.contains(xNew))
  {
    ++statistics.rejectedByRange;
    return std::nullopt;
  }
  if (occupancyMode)
  {
    // laterally the bead must stay in the column of cages A and B
    const std::size_t cage = cageLattice.cageOf(trialMolecule.second[component.startingBead].position);
    if (cage != ((xNew < 0.5 * barrier.lambda) ? state.cageA : state.cageB))
    {
      ++statistics.rejectedByRange;
      return std::nullopt;
    }
  }
  const double biasDifference = windowBias(windowIndex, xNew) - windowBias(windowIndex, xOld);

  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialMolecule.second, moleculeAtoms);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, system.frameworkCellList(),
      system.spanOfFrameworkAtoms(), trialMolecule.second, moleculeAtoms);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(), trialMolecule.second,
      moleculeAtoms);
  if (!interMolecule.has_value()) return std::nullopt;

  const RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, trialMolecule.second, moleculeAtoms);

  const RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

  if (rng.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() - biasDifference))
  {
    ++statistics.translationAccepted;
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), moleculeAtoms.begin());
    system.cellListAtomsMoved(moleculeAtoms);
    molecule = trialMolecule.first;
    return energyDifference;
  }
  return std::nullopt;
}

std::optional<RunningEnergy> DynamicallyCorrectedTSTFreeEnergy::biasedReinsertionMove(std::size_t replicaId,
                                                                                      std::size_t windowIndex)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];
  DcTSTWindowState& state = windowStates[replicaId];
  DcTSTWindowAccumulator& statistics = windowData[windowIndex];
  const DcTSTWindow& window = windows[windowIndex];

  const std::size_t selectedComponent = settings.taggedComponent;
  const Component& component = system.components[selectedComponent];
  std::span<Atom> moleculeAtoms = system.spanOfMolecule(selectedComponent, state.taggedLocalIndex);
  Molecule& molecule = system.moleculeData[taggedComponentOffset(system) + state.taggedLocalIndex];

  ++statistics.reinsertionTrials;

  // the reaction bead is placed at a uniformly random coordinate of the window support (clipped to the grid)
  // and a uniformly random in-plane position; the old bead must lie in the same support for the move to be
  // reversible
  const bool mirror = settings.mirrorSymmetry && !occupancyMode;
  const double low = std::max(grid.begin, window.center - window.supportHalfWidth);
  const double high = std::min(grid.end, window.center + window.supportHalfWidth);
  if (high <= low) return std::nullopt;
  const double xOld = coordinate(system, state, moleculeAtoms[component.startingBead].position);
  const double xOldFolded = foldedCoordinate(xOld);
  const bool oldInside = (!occupancyMode && !mirror)
                             ? std::abs(barrier.reducedCoordinate(xOld - window.center)) <= window.supportHalfWidth
                             : (xOldFolded >= low && xOldFolded < high);
  if (!oldInside) return std::nullopt;

  double xNew = low + (high - low) * rng.uniform();
  if (mirror && rng.uniform() < 0.5) xNew = -xNew;
  if (!occupancyMode) xNew = barrier.reducedCoordinate(xNew);
  if (!grid.contains(xNew)) return std::nullopt;

  double3 point;
  if (occupancyMode)
  {
    const std::optional<double3> p = pointInCageFrame(rng, system, state, xNew);
    if (!p.has_value()) return std::nullopt;
    point = p.value();
  }
  else
  {
    point = barrier.randomPointOnPlane(rng, system.simulationBox, xNew);
  }
  const double biasDifference = windowBias(windowIndex, xNew) - windowBias(windowIndex, xOld);

  const CBMC::GrowContext context = system.makeGrowContext();
  std::optional<CBMC::GrowResult> growData = CBMC::growNewMolecule(
      rng, context, component,
      CBMC::NewMoleculeIdentity{.componentId = selectedComponent,
                                .moleculeId = static_cast<std::size_t>(moleculeAtoms.front().moleculeId)},
      CBMC::GrowRequest{.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = point});
  if (!growData) return std::nullopt;
  if (system.insideBlockedPockets(component, growData->atoms)) return std::nullopt;
  if (occupancyMode)
  {
    const std::size_t cage = cageLattice.cageOf(growData->atoms[component.startingBead].position);
    if (cage != ((xNew < 0.5 * barrier.lambda) ? state.cageA : state.cageB))
    {
      ++statistics.rejectedByRange;
      return std::nullopt;
    }
  }

  const CBMC::RetraceResult retraceData =
      CBMC::retraceMolecule(rng, context, component, moleculeAtoms, {.firstBead = CBMC::FirstBeadScheme::Pinned});

  const RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, growData->atoms, moleculeAtoms);

  const double logAcceptance = -system.beta * energyFourierDifference.potentialEnergy() +
                               growData->logRosenbluthWeight - retraceData.logRosenbluthWeight - biasDifference;

  if (rng.uniform() < std::exp(logAcceptance))
  {
    ++statistics.reinsertionAccepted;
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    for (std::size_t i = 0; i < moleculeAtoms.size(); ++i) moleculeAtoms[i].position = growData->atoms[i].position;
    system.cellListAtomsMoved(moleculeAtoms);
    Molecule record = growData->molecule;
    record.atomIndex = molecule.atomIndex;
    record.numberOfAtoms = molecule.numberOfAtoms;
    molecule = record;
    return (growData->energies - retraceData.energies) + energyFourierDifference;
  }
  return std::nullopt;
}

// =================================================================================================
// Sampling
// =================================================================================================

void DynamicallyCorrectedTSTFreeEnergy::samplePilot(std::size_t replicaId, std::size_t block)
{
  System& system = systems[replicaId];
  DcTSTPilotAccumulator& accumulator = pilot[replicaId];
  const std::size_t selectedComponent = settings.taggedComponent;
  const Component& component = system.components[selectedComponent];
  const bool mirror = settings.mirrorSymmetry && !occupancyMode;

  const std::vector<std::size_t> occupancy = countOccupancies(system);
  const double3 windowShift = 0.5 * barrier.lambda * barrier.normal;
  const double kT = Units::KB * system.temperature;

  for (std::size_t m = 0; m < system.numberOfMoleculesPerComponent[selectedComponent]; ++m)
  {
    const std::span<const Atom> atoms = system.spanOfMolecule(selectedComponent, m);
    const double3 position = atoms[component.startingBead].position;

    // thermal speed of the bead along the normal for this orientation
    {
      const Molecule& molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, m)];
      const double inverseMass = beadInverseEffectiveMass(component, molecule, component.startingBead, barrier.normal);
      accumulator.beadSpeedSum += std::sqrt(2.0 * kT * inverseMass / std::numbers::pi);
      accumulator.beadSpeedCount += 1.0;
    }
    const double x = barrier.reducedCoordinate(barrier.projection(system.simulationBox, position));
    if (mirror)
    {
      accumulator.histogram[block][pilotGrid.bin(x)] += 0.5;
      accumulator.histogram[block][pilotGrid.bin(-x)] += 0.5;
    }
    else
    {
      accumulator.histogram[block][pilotGrid.bin(x)] += 1.0;
    }

    const std::size_t cage = cageLattice.cageOf(position);
    const std::size_t neighbor = cageLattice.neighbor(cage, 2.0 * windowShift);
    const std::size_t nA = occupancy[cage];
    const std::size_t nB = occupancy[neighbor];
    accumulator.occupancyPairs[block][DcTSTPilotAccumulator::occupancyIndex(nA, nB)] += 1.0;
    if (occupancyMode && nA == settings.cageOccupancies->first && nB == settings.cageOccupancies->second)
    {
      const double xc = double3::dot(
          system.simulationBox.applyPeriodicBoundaryConditions(position - cageLattice.center(cage)), barrier.normal);
      if (grid.contains(xc)) accumulator.conditionedHistogram[block][grid.bin(xc)] += 1.0;
    }
  }
}

void DynamicallyCorrectedTSTFreeEnergy::sampleWindow(std::size_t replicaId, std::size_t windowIndex, std::size_t block)
{
  System& system = systems[replicaId];
  const DcTSTWindowState& state = windowStates[replicaId];
  DcTSTWindowAccumulator& accumulator = windowData[windowIndex];
  const Component& component = system.components[settings.taggedComponent];
  const bool mirror = settings.mirrorSymmetry && !occupancyMode;

  const std::span<const Atom> atoms = system.spanOfMolecule(settings.taggedComponent, state.taggedLocalIndex);
  const double x = coordinate(system, state, atoms[component.startingBead].position);
  if (!grid.contains(x)) return;
  if (mirror)
  {
    accumulator.histogram[block][grid.bin(x)] += 0.5;
    accumulator.histogram[block][grid.bin(-x)] += 0.5;
  }
  else
  {
    accumulator.histogram[block][grid.bin(x)] += 1.0;
  }
  const double xf = foldedCoordinate(x);
  accumulator.count[block] += 1.0;
  accumulator.sumX[block] += xf;
  accumulator.sumXSquared[block] += xf * xf;
  accumulator.sumBias[block] += windowBias(windowIndex, x);

  // occupancies of the two cages sharing the window while the tagged molecule is on the dividing plane
  if (!occupancyMode && std::abs(x) < DcTSTWindowAccumulator::barrierOccupancyHalfWidth)
  {
    const double3 position = atoms[component.startingBead].position;
    const std::vector<std::size_t> occupancy = countOccupancies(system);
    const double3 windowShift = 0.5 * barrier.lambda * barrier.normal;
    const std::size_t behind = cageLattice.cageOf(position - windowShift);
    const std::size_t ahead = cageLattice.cageOf(position + windowShift);
    const std::size_t own = cageLattice.cageOf(position);
    // the other molecules only; the tagged molecule is then counted in whichever cage plays A
    std::size_t mBehind = occupancy[behind] - ((own == behind) ? 1uz : 0uz);
    std::size_t mAhead = occupancy[ahead] - ((own == ahead) ? 1uz : 0uz);
    accumulator.barrierOccupancyPairs[DcTSTPilotAccumulator::occupancyIndex(mBehind + 1uz, mAhead)] += 0.5;
    accumulator.barrierOccupancyPairs[DcTSTPilotAccumulator::occupancyIndex(mAhead + 1uz, mBehind)] += 0.5;
    accumulator.barrierOccupancySamples += 1.0;
  }
}

double DynamicallyCorrectedTSTFreeEnergy::beadThermalSpeed() const
{
  double sum{0.0};
  double count{0.0};
  for (const DcTSTPilotAccumulator& accumulator : pilot)
  {
    sum += accumulator.beadSpeedSum;
    count += accumulator.beadSpeedCount;
  }
  if (count > 0.0) return sum / count;
  const System& system = systems.front();
  const Component& tagged = system.components[settings.taggedComponent];
  return std::sqrt(2.0 * Units::KB * system.temperature / tagged.totalMass / std::numbers::pi);
}

// =================================================================================================
// Analysis
// =================================================================================================

std::vector<double> DynamicallyCorrectedTSTFreeEnergy::pooledPilotHistogram(std::optional<std::size_t> block) const
{
  std::vector<double> histogram(grid.numberOfBins, 0.0);
  for (const DcTSTPilotAccumulator& accumulator : pilot)
  {
    for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
    {
      if (block.has_value() && b != block.value()) continue;
      const std::vector<double>& source = occupancyMode ? accumulator.conditionedHistogram[b] : accumulator.histogram[b];
      for (std::size_t i = 0; i < grid.numberOfBins; ++i) histogram[i] += source[i];
    }
  }
  return histogram;
}

DynamicallyCorrectedTSTFreeEnergy::ProfileResult DynamicallyCorrectedTSTFreeEnergy::analyzeProfile() const
{
  ProfileResult result;
  const std::size_t numberOfBins = grid.numberOfBins;
  const std::size_t numberOfStates = 1uz + windows.size();

  // reduced potentials: the pilot is unbiased but, in the occupancy-resolved mode, confined to cage A
  std::vector<std::vector<double>> reducedPotential(numberOfStates, std::vector<double>(numberOfBins, 0.0));
  if (occupancyMode)
  {
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (grid.center(b) >= 0.5 * barrier.lambda) reducedPotential[0][b] = forbiddenBias;
    }
  }
  for (std::size_t k = 0; k < windows.size(); ++k) reducedPotential[k + 1uz] = windows[k].bias;

  // the counts of a set of blocks (all blocks for an empty set)
  auto countsOfBlocks = [&](const std::vector<std::size_t>& blocks)
  {
    std::vector<std::vector<double>> counts(numberOfStates, std::vector<double>(numberOfBins, 0.0));
    if (blocks.empty())
    {
      counts[0] = pooledPilotHistogram(std::nullopt);
    }
    else
    {
      for (std::size_t j : blocks)
      {
        const std::vector<double> more = pooledPilotHistogram(j);
        for (std::size_t b = 0; b < numberOfBins; ++b) counts[0][b] += more[b];
      }
    }
    for (std::size_t k = 0; k < windows.size(); ++k)
    {
      const DcTSTWindowAccumulator& accumulator = windowData[k];
      for (std::size_t j = 0; j < accumulator.numberOfBlocks; ++j)
      {
        if (!blocks.empty() && std::find(blocks.begin(), blocks.end(), j) == blocks.end()) continue;
        for (std::size_t b = 0; b < numberOfBins; ++b) counts[k + 1uz][b] += accumulator.histogram[j][b];
      }
    }
    return counts;
  };
  auto countsOf = [&](std::optional<std::size_t> block)
  {
    return block.has_value() ? countsOfBlocks({block.value()}) : countsOfBlocks({});
  };

  // --- a solution: normalized density, free energy, barrier estimates ------------------------------
  struct Solution
  {
    bool available{false};
    std::vector<double> density;  ///< P per fine bin (normalized over the cage), 0 where undefined.
    std::vector<bool> defined;
    double barrierProbabilityBin{0.0};
    QuadraticFit fit;
    double barrierProbability{0.0};
    double barrierOffset{0.0};
    std::vector<double> freeEnergies;
    std::size_t iterations{0};
    double residual{0.0};
    bool converged{false};
  };
  auto solveFor = [&](const std::vector<std::vector<double>>& counts, const std::vector<double>& initial,
                      const std::vector<double>& fitWeights) -> Solution
  {
    Solution solution;
    double total{0.0};
    for (const std::vector<double>& c : counts) total += std::accumulate(c.begin(), c.end(), 0.0);
    if (total <= 0.0) return solution;
    const DcTSTWham::Result wham =
        DcTSTWham::solve(counts, reducedPotential, settings.whamTolerance, settings.whamMaximumIterations, initial);
    solution.freeEnergies = wham.freeEnergies;
    solution.iterations = wham.iterations;
    solution.residual = wham.residual;
    solution.converged = wham.converged;
    solution.density.assign(numberOfBins, 0.0);
    solution.defined.assign(numberOfBins, false);
    double largest = DcTSTWham::logZero;
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (wham.logDensity[b] > DcTSTWham::logZeroThreshold) largest = std::max(largest, wham.logDensity[b]);
    }
    if (largest <= DcTSTWham::logZeroThreshold) return solution;
    double norm{0.0};
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (wham.logDensity[b] <= DcTSTWham::logZeroThreshold) continue;
      solution.defined[b] = true;
      solution.density[b] = std::exp(wham.logDensity[b] - largest);
      const double x = grid.center(b);
      if (x >= grid.normalizationBegin && x < grid.normalizationEnd) norm += solution.density[b] * grid.width;
    }
    if (norm <= 0.0) return solution;
    for (double& p : solution.density) p /= norm;
    solution.available = true;

    // the bin estimate: the mean density of the bins whose centers are within 1.5 bin widths of the plane (the
    // bin on the plane and its two neighbors; the plane is on a bin center for an odd number of bins and on a
    // bin edge otherwise, 1.5 avoids the equality case)
    {
      double sum{0.0};
      std::size_t n{0};
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        if (solution.defined[b] && std::abs(grid.center(b) - grid.barrier) < 1.5 * grid.width)
        {
          sum += solution.density[b];
          ++n;
        }
      }
      solution.barrierProbabilityBin = (n > 0) ? sum / static_cast<double>(n) : 0.0;
    }
    // the parabolic fit of ln P over |x - x*| <= fit range
    {
      std::vector<double> t, y, w;
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        const double dx = grid.center(b) - grid.barrier;
        if (!solution.defined[b] || std::abs(dx) > settings.barrierFitRange) continue;
        t.push_back(dx);
        y.push_back(std::log(solution.density[b]));
        w.push_back(fitWeights.empty() ? 1.0 : fitWeights[b]);
      }
      solution.fit = fitQuadratic(t, y, w);
      if (solution.fit.valid)
      {
        solution.barrierProbability = std::exp(solution.fit.a);
        // ln P = a + b t + c t^2 has its minimum (the maximum of beta F) at t = -b / (2c) when c > 0
        solution.barrierOffset = (solution.fit.c > 0.0) ? -solution.fit.b / (2.0 * solution.fit.c) : 0.0;
      }
      else
      {
        solution.barrierProbability = solution.barrierProbabilityBin;
      }
    }
    return solution;
  };

  // --- full solution, first unweighted to get the block errors per bin, then weighted --------------
  const std::vector<std::vector<double>> allCounts = countsOf(std::nullopt);
  Solution full = solveFor(allCounts, {}, {});
  if (!full.available) return result;
  // the convergence of the cold start (the later solves are warm-started from it and converge immediately)
  const std::size_t coldIterations = full.iterations;
  const double coldResidual = full.residual;
  const bool coldConverged = full.converged;

  std::vector<Solution> blockSolutions;
  for (std::size_t j = 0; j < numberOfBlocks; ++j)
  {
    Solution s = solveFor(countsOf(j), full.freeEnergies, {});
    if (s.available) blockSolutions.push_back(std::move(s));
  }
  // per-bin errors of beta F from the blocks (relative to each block's own minimum region: use ln P directly,
  // the normalization is common)
  std::vector<double> fineError(numberOfBins, 0.0);
  std::vector<double> fitWeights(numberOfBins, 1.0);
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    std::vector<double> values;
    for (const Solution& s : blockSolutions)
    {
      if (s.defined[b]) values.push_back(-std::log(s.density[b]));
    }
    fineError[b] = meanAndError(values).second;
    if (values.size() >= 2 && fineError[b] > 0.0) fitWeights[b] = 1.0 / (fineError[b] * fineError[b]);
  }
  full = solveFor(allCounts, full.freeEnergies, fitWeights);
  // re-solve the blocks with the same weights for the fit
  blockSolutions.clear();
  for (std::size_t j = 0; j < numberOfBlocks; ++j)
  {
    Solution s = solveFor(countsOf(j), full.freeEnergies, fitWeights);
    if (s.available) blockSolutions.push_back(std::move(s));
  }

  result.available = true;
  result.converged = coldConverged && full.converged;
  result.iterations = coldIterations;
  result.residual = coldResidual;
  result.stateFreeEnergies = full.freeEnergies;

  // fine profile
  result.fineFreeEnergy.assign(numberOfBins, 0.0);
  result.fineDefined = full.defined;
  result.fineFreeEnergyError = fineError;
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    if (full.defined[b])
    {
      result.fineFreeEnergy[b] = -std::log(full.density[b]);
      minimum = std::min(minimum, result.fineFreeEnergy[b]);
    }
  }
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    if (full.defined[b]) result.fineFreeEnergy[b] -= minimum;
  }

  // barrier estimates and their block errors
  result.barrierProbabilityBin = full.barrierProbabilityBin;
  result.barrierProbability = full.barrierProbability;
  result.fitValid = full.fit.valid;
  result.fitPoints = full.fit.points;
  result.barrierOffset = full.barrierOffset;
  result.curvature = -2.0 * full.fit.c;  // beta F'' = -(ln P)''
  result.barrierFreeEnergy = -std::log(std::max(full.barrierProbability, 1e-300)) - minimum;
  result.barrierFreeEnergyBin = -std::log(std::max(full.barrierProbabilityBin, 1e-300)) - minimum;
  {
    std::vector<double> pBin, pFit, fBin, fFit, offsets;
    for (const Solution& s : blockSolutions)
    {
      pBin.push_back(s.barrierProbabilityBin);
      pFit.push_back(s.barrierProbability);
      fFit.push_back(-std::log(std::max(s.barrierProbability, 1e-300)));
      offsets.push_back(s.barrierOffset);
    }
    result.barrierProbabilityBinError = meanAndError(pBin).second;
    result.barrierProbabilityError = meanAndError(pFit).second;
    result.barrierFreeEnergyError = meanAndError(fFit).second;
    result.barrierOffsetError = meanAndError(offsets).second;
  }
  // half the number of blocks: merge pairs of consecutive blocks (an odd last block joins the last pair)
  {
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t j = 0; j < numberOfBlocks; ++j)
    {
      if (j % 2 == 0 && !(j + 1 == numberOfBlocks && !groups.empty())) groups.emplace_back();
      groups.back().push_back(j);
    }
    std::vector<double> pFit;
    for (const std::vector<std::size_t>& group : groups)
    {
      Solution s = solveFor(countsOfBlocks(group), full.freeEnergies, fitWeights);
      if (s.available) pFit.push_back(s.barrierProbability);
    }
    result.barrierProbabilityErrorHalfBlocks = meanAndError(pFit).second;
    result.numberOfHalfBlocks = pFit.size();
  }

  // reporting grid
  const std::size_t numberOfReportBins = settings.numberOfFreeEnergyBins;
  result.q.assign(numberOfReportBins, 0.0);
  result.freeEnergy.assign(numberOfReportBins, 0.0);
  result.freeEnergyError.assign(numberOfReportBins, 0.0);
  result.defined.assign(numberOfReportBins, false);
  result.pilotFreeEnergy.assign(numberOfReportBins, 0.0);
  result.pilotDefined.assign(numberOfReportBins, false);
  const double reportWidth = (grid.end - grid.begin) / static_cast<double>(numberOfReportBins);
  auto reportBin = [&](double x) -> std::size_t
  {
    if (!occupancyMode) return barrier.bin(x, numberOfReportBins);
    return std::min(numberOfReportBins - 1uz, static_cast<std::size_t>(std::floor((x - grid.begin) / reportWidth)));
  };
  for (std::size_t i = 0; i < numberOfReportBins; ++i)
  {
    result.q[i] = occupancyMode ? grid.begin + (static_cast<double>(i) + 0.5) * reportWidth
                                : barrier.binCenter(i, numberOfReportBins);
  }
  {
    std::vector<double> density(numberOfReportBins, 0.0);
    std::vector<double> errorSum(numberOfReportBins, 0.0);
    std::vector<std::size_t> errorCount(numberOfReportBins, 0uz);
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (!full.defined[b]) continue;
      const std::size_t i = reportBin(grid.center(b));
      density[i] += full.density[b] * grid.width;
      errorSum[i] += fineError[b];
      ++errorCount[i];
    }
    double reportMinimum = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < numberOfReportBins; ++i)
    {
      if (density[i] > 0.0)
      {
        result.freeEnergy[i] = -std::log(density[i] / reportWidth);
        result.defined[i] = true;
        result.freeEnergyError[i] = errorSum[i] / static_cast<double>(std::max(1uz, errorCount[i]));
        reportMinimum = std::min(reportMinimum, result.freeEnergy[i]);
      }
    }
    for (std::size_t i = 0; i < numberOfReportBins; ++i)
    {
      if (result.defined[i]) result.freeEnergy[i] -= reportMinimum;
    }
  }
  // the pilot alone
  {
    const std::vector<double> histogram = allCounts[0];
    double norm{0.0};
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      const double x = grid.center(b);
      if (x >= grid.normalizationBegin && x < grid.normalizationEnd) norm += histogram[b] * grid.width;
    }
    if (norm > 0.0)
    {
      std::vector<double> density(numberOfReportBins, 0.0);
      double sum{0.0};
      std::size_t n{0};
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        density[reportBin(grid.center(b))] += histogram[b] / norm * grid.width;
        if (std::abs(grid.center(b) - grid.barrier) < 1.5 * grid.width)
        {
          sum += histogram[b] / norm;
          ++n;
        }
      }
      double pilotMinimum = std::numeric_limits<double>::max();
      for (std::size_t i = 0; i < numberOfReportBins; ++i)
      {
        if (density[i] > 0.0)
        {
          result.pilotFreeEnergy[i] = -std::log(density[i] / reportWidth);
          result.pilotDefined[i] = true;
          pilotMinimum = std::min(pilotMinimum, result.pilotFreeEnergy[i]);
        }
      }
      for (std::size_t i = 0; i < numberOfReportBins; ++i)
      {
        if (result.pilotDefined[i]) result.pilotFreeEnergy[i] -= pilotMinimum;
      }
      result.pilotBarrierProbability = (n > 0) ? sum / static_cast<double>(n) : 0.0;
      result.pilotBarrierDefined = result.pilotBarrierProbability > 0.0;
      // block error of the pilot barrier probability
      std::vector<double> values;
      for (std::size_t j = 0; j < numberOfBlocks; ++j)
      {
        const std::vector<double> h = pooledPilotHistogram(j);
        double bnorm{0.0}, bsum{0.0};
        std::size_t bn{0};
        for (std::size_t b = 0; b < numberOfBins; ++b)
        {
          const double x = grid.center(b);
          if (x >= grid.normalizationBegin && x < grid.normalizationEnd) bnorm += h[b] * grid.width;
          if (std::abs(x - grid.barrier) < 1.5 * grid.width)
          {
            bsum += h[b];
            ++bn;
          }
        }
        if (bnorm > 0.0 && bn > 0) values.push_back(bsum / static_cast<double>(bn) / bnorm);
      }
      result.pilotBarrierProbabilityError = meanAndError(values).second;
    }
  }

  // window diagnostics: mean, sigma, counts, mean bias, the overlap of consecutive windows and the full
  // overlap matrix of all states (pilot included)
  {
    std::vector<std::size_t> order(windows.size());
    std::iota(order.begin(), order.end(), 0uz);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return windows[a].center < windows[b].center; });
    result.windowMean.assign(windows.size(), 0.0);
    result.windowSigma.assign(windows.size(), 0.0);
    result.windowCounts.assign(windows.size(), 0.0);
    result.windowMeanBias.assign(windows.size(), 0.0);
    std::vector<std::vector<double>> normalized(numberOfStates);
    for (std::size_t s = 0; s < numberOfStates; ++s)
    {
      normalized[s] = allCounts[s];
      const double total = std::accumulate(normalized[s].begin(), normalized[s].end(), 0.0);
      if (total > 0.0)
        for (double& v : normalized[s]) v /= total;
    }
    for (std::size_t k = 0; k < windows.size(); ++k)
    {
      const DcTSTWindowAccumulator& accumulator = windowData[k];
      double n{0.0}, sx{0.0}, sxx{0.0}, sw{0.0};
      for (std::size_t j = 0; j < accumulator.numberOfBlocks; ++j)
      {
        n += accumulator.count[j];
        sx += accumulator.sumX[j];
        sxx += accumulator.sumXSquared[j];
        sw += accumulator.sumBias[j];
      }
      result.windowCounts[k] = n;
      if (n > 0.0)
      {
        result.windowMean[k] = sx / n;
        result.windowSigma[k] = std::sqrt(std::max(0.0, sxx / n - (sx / n) * (sx / n)));
        result.windowMeanBias[k] = sw / n;
      }
    }
    result.overlapMatrix.assign(numberOfStates * numberOfStates, 0.0);
    for (std::size_t i = 0; i < numberOfStates; ++i)
    {
      for (std::size_t j = i; j < numberOfStates; ++j)
      {
        double overlap{0.0};
        for (std::size_t b = 0; b < numberOfBins; ++b) overlap += std::min(normalized[i][b], normalized[j][b]);
        result.overlapMatrix[i * numberOfStates + j] = overlap;
        result.overlapMatrix[j * numberOfStates + i] = overlap;
      }
    }
    // overlap with the next distinct center (repeated windows share a center and are skipped)
    result.overlap.assign(windows.size(), 0.0);
    for (std::size_t i = 0; i < order.size(); ++i)
    {
      std::size_t j = i + 1;
      while (j < order.size() && windows[order[j]].center <= windows[order[i]].center) ++j;
      if (j == order.size()) continue;
      result.overlap[order[i]] = result.overlapMatrix[(order[i] + 1uz) * numberOfStates + (order[j] + 1uz)];
    }
  }

  // consistency of the pilot with the reweighted profile where the pilot is well sampled (the reweighting is
  // dominated by the pilot there, so a mismatch beyond the statistics signals a biased window or a wrong bias
  // table rather than noise)
  {
    const std::vector<double>& pilotHistogram = allCounts[0];
    double pilotNorm{0.0};
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      const double x = grid.center(b);
      if (x >= grid.normalizationBegin && x < grid.normalizationEnd) pilotNorm += pilotHistogram[b] * grid.width;
    }
    std::vector<double> differences;
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (pilotNorm <= 0.0 || !full.defined[b]) continue;
      if (pilotHistogram[b] < ProfileResult::pilotConsistencyMinimumCounts) continue;
      differences.push_back(-std::log(pilotHistogram[b] / pilotNorm) + std::log(full.density[b]));
    }
    if (!differences.empty())
    {
      const double mean = std::accumulate(differences.begin(), differences.end(), 0.0) /
                          static_cast<double>(differences.size());
      double sumSquares{0.0};
      for (double d : differences) sumSquares += (d - mean) * (d - mean);
      result.pilotConsistencyRms = std::sqrt(sumSquares / static_cast<double>(differences.size()));
      result.pilotConsistencyBins = differences.size();
    }
  }

  return result;
}

DynamicallyCorrectedTSTFreeEnergy::OccupancyDistribution DynamicallyCorrectedTSTFreeEnergy::analyzeOccupancies() const
{
  OccupancyDistribution result;
  result.size = DcTSTPilotAccumulator::maximumOccupancy + 1uz;
  result.joint.assign(result.size * result.size, 0.0);
  result.marginalA.assign(result.size, 0.0);
  result.marginalB.assign(result.size, 0.0);

  std::vector<double> blockRequested;
  const std::optional<std::size_t> requestedIndex =
      occupancyMode ? std::optional<std::size_t>(DcTSTPilotAccumulator::occupancyIndex(
                          settings.cageOccupancies->first, settings.cageOccupancies->second))
                    : std::nullopt;
  for (const DcTSTPilotAccumulator& accumulator : pilot)
  {
    for (std::size_t b = 0; b < accumulator.numberOfBlocks; ++b)
    {
      double blockTotal{0.0};
      for (std::size_t i = 0; i < result.joint.size(); ++i)
      {
        result.joint[i] += accumulator.occupancyPairs[b][i];
        blockTotal += accumulator.occupancyPairs[b][i];
      }
      if (blockTotal > 0.0 && requestedIndex.has_value())
      {
        blockRequested.push_back(accumulator.occupancyPairs[b][requestedIndex.value()] / blockTotal);
      }
    }
  }
  result.total = std::accumulate(result.joint.begin(), result.joint.end(), 0.0);
  if (result.total <= 0.0) return result;
  result.available = true;
  for (std::size_t nA = 0; nA < result.size; ++nA)
  {
    for (std::size_t nB = 0; nB < result.size; ++nB)
    {
      double& p = result.joint[nA * result.size + nB];
      p /= result.total;
      result.marginalA[nA] += p;
      result.marginalB[nB] += p;
      result.meanA += static_cast<double>(nA) * p;
      result.meanB += static_cast<double>(nB) * p;
    }
  }
  if (requestedIndex.has_value())
  {
    result.requestedProbability = result.joint[requestedIndex.value()];
    result.requestedProbabilityError = meanAndError(blockRequested).second;
  }

  // P(n_A, n_B | q*) from the window stage
  result.jointAtBarrier.assign(result.size * result.size, 0.0);
  for (const DcTSTWindowAccumulator& accumulator : windowData)
  {
    for (std::size_t i = 0; i < result.jointAtBarrier.size(); ++i)
    {
      result.jointAtBarrier[i] += accumulator.barrierOccupancyPairs[i];
    }
    result.totalAtBarrier += accumulator.barrierOccupancySamples;
  }
  if (result.totalAtBarrier > 0.0)
  {
    result.availableAtBarrier = true;
    for (double& p : result.jointAtBarrier) p /= result.totalAtBarrier;
  }
  return result;
}

// =================================================================================================
// Output
// =================================================================================================

void DynamicallyCorrectedTSTFreeEnergy::output()
{
  const System& system = systems.front();
  const Component& tagged = system.components[settings.taggedComponent];
  const ProfileResult profile = analyzeProfile();
  const OccupancyDistribution occupancies = analyzeOccupancies();

  const std::string suffix = std::format("{}_{}", system.temperature, system.input_pressure);
  const std::string label = occupancyMode ? "conditioned_" : "";

  // --- occupancy distribution file ---------------------------------------------------------------
  if (occupancies.available)
  {
    std::ofstream file(std::format("dctst/cage_occupancy_distribution_{}.txt", suffix));
    std::print(file,
               "# dcTST joint distribution P(n_A, n_B) of the occupancy n_A of the cage of a molecule of component "
               "'{}' (itself included)\n# and the occupancy n_B of the neighboring cage across its window towards "
               "+normal (all components counted); {} samples\n",
               tagged.name, occupancies.total);
    std::print(file, "# <n_A> = {:.4f} (loading + 1 in the mean-field limit), <n_B> = {:.4f}\n", occupancies.meanA,
               occupancies.meanB);
    std::print(file, "# column 1: n_A\n# column 2: n_B\n# column 3: P(n_A, n_B)\n# column 4: P(n_A)\n# column 5: P(n_B)\n");
    if (occupancies.availableAtBarrier)
    {
      std::print(file,
                 "# column 6: P(n_A, n_B | q*), the occupancies while the tagged molecule is within {} A of the dividing "
                 "plane ({} window-stage samples)\n# column 7: k^TST(n_A -> n_B) / k^TST = P(n_A, n_B | q*) / P(n_A, n_B)\n",
                 DcTSTWindowAccumulator::barrierOccupancyHalfWidth, occupancies.totalAtBarrier);
    }
    for (std::size_t nA = 0; nA < occupancies.size; ++nA)
    {
      for (std::size_t nB = 0; nB < occupancies.size; ++nB)
      {
        const std::size_t index = nA * occupancies.size + nB;
        const double p = occupancies.joint[index];
        const double pBarrier = occupancies.availableAtBarrier ? occupancies.jointAtBarrier[index] : 0.0;
        if (p <= 0.0 && pBarrier <= 0.0) continue;
        std::print(file, "{:4} {:4} {:14.8f} {:14.8f} {:14.8f}", nA, nB, p, occupancies.marginalA[nA],
                   occupancies.marginalB[nB]);
        if (occupancies.availableAtBarrier)
        {
          std::print(file, " {:14.8f} {:14.6f}", pBarrier, (p > 0.0) ? pBarrier / p : 0.0);
        }
        std::print(file, "\n");
      }
    }
  }

  // --- profile files ------------------------------------------------------------------------------
  if (profile.available)
  {
    {
      std::ofstream file(std::format("dctst/{}free_energy_profile_{}.txt", label, suffix));
      std::print(file, "# dcTST free-energy profile of component '{}' along the barrier normal (multistate reweighting "
                       "of the pilot histogram and {} umbrella windows)\n",
                 tagged.name, windows.size());
      if (occupancyMode)
        std::print(file, "# x is measured from the center of cage A, the window to cage B is at x = {:.4f} [A]; "
                         "conditioned on n_A = {}, n_B = {}\n",
                   0.5 * barrier.lambda, settings.cageOccupancies->first, settings.cageOccupancies->second);
      else
        std::print(file, "# bin 0 is centered on the barrier (q = 0); one period [-lambda/2, lambda/2), lambda = {} [A]\n",
                   barrier.lambda);
      std::print(file, "# beta F is shifted so that its minimum is zero; errors are 2 sigma over the blocks; unsampled bins 'inf'\n");
      std::print(file, "# column 1: q [A]\n# column 2: beta F(q) reweighted\n# column 3: error\n# column 4: beta F(q) pilot "
                       "histogram only\n");
      for (std::size_t i = 0; i < profile.q.size(); ++i)
      {
        std::print(file, "{:12.6f} {} {:12.6f} {}\n", profile.q[i],
                   profile.defined[i] ? std::format("{:14.6f}", profile.freeEnergy[i]) : std::format("{:>14}", "inf"),
                   profile.freeEnergyError[i],
                   profile.pilotDefined[i] ? std::format("{:14.6f}", profile.pilotFreeEnergy[i])
                                           : std::format("{:>14}", "inf"));
      }
    }
    {
      std::ofstream file(std::format("dctst/{}free_energy_profile_fine_{}.txt", label, suffix));
      std::print(file, "# dcTST free-energy profile of component '{}' on the fine grid of the reweighting ({} bins of {:.5f} A)\n",
                 tagged.name, grid.numberOfBins, grid.width);
      std::print(file, "# column 1: x [A]\n# column 2: beta F(x)\n# column 3: error (2 sigma over the blocks)\n# column 4: "
                       "flattening bias w_0(x) [kT]\n");
      for (std::size_t b = 0; b < grid.numberOfBins; ++b)
      {
        std::print(file, "{:12.6f} {} {:12.6f} {:12.6f}\n", grid.center(b),
                   profile.fineDefined[b] ? std::format("{:14.6f}", profile.fineFreeEnergy[b])
                                          : std::format("{:>14}", "inf"),
                   profile.fineFreeEnergyError[b], flatteningBias[b]);
      }
    }
    if (settings.writeWindowHistograms)
    {
      std::ofstream file(std::format("dctst/{}window_histograms_{}.txt", label, suffix));
      std::print(file, "# dcTST umbrella windows: fine histograms of the tagged molecule (counts, all blocks) and the "
                       "reduced bias per window [kT]\n");
      std::print(file, "# window centers [A]:");
      for (const DcTSTWindow& w : windows) std::print(file, " {:.4f}", w.center);
      std::print(file, "\n# stiffness beta kappa [1/A^2]: {:.6f}\n", windows.empty() ? 0.0 : windows.front().stiffness);
      std::print(file, "# barrier x* [A]: {:.6f}\n# normalization range [A]: {:.6f} {:.6f}\n# fine bin width [A]: {:.6f}\n",
                 grid.barrier, grid.normalizationBegin, grid.normalizationEnd, grid.width);
      std::print(file, "# fit range |x - x*| <= [A]: {:.6f}\n", settings.barrierFitRange);
      std::print(file, "# pilot forbidden above [A]: {}\n",
                 occupancyMode ? std::format("{:.6f}", grid.barrier) : std::string("none"));
      std::print(file, "# column 1: x [A]\n# column 2: pilot histogram (all molecules)\n");
      std::print(file, "# columns 3..{}: counts of window 0..{}\n# columns {}..{}: bias of window 0..{}\n",
                 2 + windows.size(), windows.size() - 1, 3 + windows.size(), 2 + 2 * windows.size(), windows.size() - 1);
      const std::vector<double> pilotHistogram = pooledPilotHistogram();
      for (std::size_t b = 0; b < grid.numberOfBins; ++b)
      {
        std::print(file, "{:12.6f} {:14.1f}", grid.center(b), pilotHistogram[b]);
        for (std::size_t k = 0; k < windows.size(); ++k)
        {
          double counts{0.0};
          for (std::size_t j = 0; j < windowData[k].numberOfBlocks; ++j) counts += windowData[k].histogram[j][b];
          std::print(file, " {:12.1f}", counts);
        }
        for (std::size_t k = 0; k < windows.size(); ++k) std::print(file, " {:12.5f}", windows[k].bias[b]);
        std::print(file, "\n");
      }
    }
  }

  // --- summary ----------------------------------------------------------------------------------
  std::print(stream, "\n\nResults\n");
  std::print(stream, "===============================================================================\n\n");

  const double kBT_over_m = Units::KB * system.temperature / tagged.totalMass;
  const double thermalSpeed = std::sqrt(2.0 * kBT_over_m / std::numbers::pi);
  const double meanSpeed = beadThermalSpeed();
  const double lambdaSquared = barrier.lambda * barrier.lambda;
  const double conversion = Units::DiffusionConversionFactor;

  if (occupancies.available)
  {
    std::print(stream, "Cage occupancies (pilot stage, {} samples):\n", occupancies.total);
    std::print(stream, "    <n_A> (cage of a molecule, itself included): {:.4f}\n", occupancies.meanA);
    std::print(stream, "    <n_B> (neighbor across its window):           {:.4f}\n", occupancies.meanB);
    std::print(stream, "    P(n_A, n_B) written to dctst/cage_occupancy_distribution_{}.txt\n", suffix);
    if (occupancyMode)
    {
      std::print(stream, "    P(n_A = {}, n_B = {}):                      {:.6f} +/- {:.6f}\n",
                 settings.cageOccupancies->first, settings.cageOccupancies->second, occupancies.requestedProbability,
                 occupancies.requestedProbabilityError);
    }
    if (occupancies.availableAtBarrier)
    {
      std::print(stream,
                 "    occupancies with the tagged molecule within {} A of the plane ({} window-stage samples), the\n"
                 "    ratio P(n_A, n_B | q*) / P(n_A, n_B) = k^TST(n_A -> n_B) / k^TST (pairs with P(n_A, n_B) >= 0.005):\n",
                 DcTSTWindowAccumulator::barrierOccupancyHalfWidth, occupancies.totalAtBarrier);
      std::print(stream, "    {:>4} {:>4} {:>12} {:>14} {:>10} {:>8}\n", "n_A", "n_B", "P(n_A,n_B)", "P(n_A,n_B|q*)",
                 "ratio", "samples");
      double consistency{0.0};
      for (std::size_t nA = 0; nA < occupancies.size; ++nA)
      {
        for (std::size_t nB = 0; nB < occupancies.size; ++nB)
        {
          const std::size_t index = nA * occupancies.size + nB;
          const double p = occupancies.joint[index];
          const double pBarrier = occupancies.jointAtBarrier[index];
          consistency += pBarrier;
          if (p < 0.005) continue;
          std::print(stream, "    {:>4} {:>4} {:>12.6f} {:>14.6f} {:>10.4f} {:>8.0f}\n", nA, nB, p, pBarrier, pBarrier / p,
                     pBarrier * occupancies.totalAtBarrier);
        }
      }
      std::print(stream, "    (sum_(n_A,n_B) P(n_A,n_B) k^TST(n_A -> n_B) / k^TST = {:.4f}, exactly 1 by construction)\n",
                 consistency);
    }
    std::print(stream, "\n");
  }

  if (!profile.available)
  {
    std::print(stream, "Free-energy profile: no samples\n\n");
  }
  else
  {
    std::print(stream,
               "Windows ({}; mean and sigma of the sampled coordinate, counts, mean bias <w_k>, overlap with the next "
               "distinct center):\n",
               windows.size());
    std::print(stream, "    {:>4} {:>10} {:>10} {:>10} {:>12} {:>10} {:>10} {:>10} {:>10}\n", "k", "center", "<x>",
               "sigma", "counts", "<w_k>", "overlap", "g_k", "acc.tr/re");
    for (std::size_t k = 0; k < windows.size(); ++k)
    {
      const DcTSTWindowAccumulator& a = windowData[k];
      const double accTr = static_cast<double>(a.translationTotalAccepted + a.translationAccepted) /
                           static_cast<double>(std::max(1uz, a.translationTotalTrials + a.translationTrials));
      const double accRe = static_cast<double>(a.reinsertionTotalAccepted + a.reinsertionAccepted) /
                           static_cast<double>(std::max(1uz, a.reinsertionTotalTrials + a.reinsertionTrials));
      std::print(stream, "    {:>4} {:>10.4f} {:>10.4f} {:>10.4f} {:>12.0f} {:>10.4f} {:>10.4f} {:>10.4f} {:>4.2f}/{:<4.2f}{}\n",
                 k, windows[k].center, profile.windowMean[k], profile.windowSigma[k], profile.windowCounts[k],
                 profile.windowMeanBias[k], profile.overlap[k], profile.stateFreeEnergies[k + 1uz], accTr, accRe,
                 windows[k].source == k ? "" : std::format("  (repeat of {})", windows[k].source));
    }
    double minimumOverlap = 1.0;
    for (std::size_t k = 0; k < windows.size(); ++k)
    {
      bool last = true;
      for (std::size_t j = 0; j < windows.size(); ++j)
        if (windows[j].center > windows[k].center) last = false;
      if (!last) minimumOverlap = std::min(minimumOverlap, profile.overlap[k]);
    }
    std::print(stream, "    minimum overlap of consecutive windows: {:.4f}{}\n\n", minimumOverlap,
               minimumOverlap < 0.05 ? "  WARNING: poor overlap, add windows or lower the stiffness" : "");

    const std::size_t numberOfStates = windows.size() + 1uz;
    std::print(stream, "Overlap matrix of the histograms (state 0 = pilot, state k = window k-1):\n    {:>5}", "");
    for (std::size_t j = 0; j < numberOfStates; ++j) std::print(stream, " {:>5}", j);
    std::print(stream, "\n");
    for (std::size_t i = 0; i < numberOfStates; ++i)
    {
      std::print(stream, "    {:>5}", i);
      for (std::size_t j = 0; j < numberOfStates; ++j)
      {
        std::print(stream, " {:>5.2f}", profile.overlapMatrix[i * numberOfStates + j]);
      }
      std::print(stream, "\n");
    }
    std::print(stream, "\n");

    std::print(stream, "Multistate reweighting (WHAM/binned MBAR, {} states x {} fine bins):\n", numberOfStates,
               grid.numberOfBins);
    std::print(stream, "    converged:                           {} ({} iterations, residual {:.2e})\n",
               profile.converged ? "yes" : "NO", profile.iterations, profile.residual);
    std::print(stream, "    g_0 (pilot) = 0 by gauge; window free energies above\n");
    if (profile.pilotConsistencyBins > 0)
    {
      std::print(stream,
                 "    pilot vs reweighted profile:         rms {:.4f} kT over {} fine bins with >= {:.0f} pilot counts\n",
                 profile.pilotConsistencyRms, profile.pilotConsistencyBins, ProfileResult::pilotConsistencyMinimumCounts);
    }
    std::print(stream, "\n");

    std::print(stream, "Free-energy profile (reweighted):\n");
    std::print(stream, "    beta F(q*) - beta F(min), top fit:    {:.5f} +/- {:.5f}\n", profile.barrierFreeEnergy,
               profile.barrierFreeEnergyError);
    std::print(stream, "    beta F(q*) - beta F(min), bin:        {:.5f}\n", profile.barrierFreeEnergyBin);
    std::print(stream, "    P(q*) at the plane, top fit:          {:.6e} +/- {:.3e} [1/A] ({} points, |x - x*| <= {} A)\n",
               profile.barrierProbability, profile.barrierProbabilityError, profile.fitPoints, settings.barrierFitRange);
    std::print(stream, "    P(q*) at the plane, bin:              {:.6e} +/- {:.3e} [1/A]\n",
               profile.barrierProbabilityBin, profile.barrierProbabilityBinError);
    std::print(stream, "    error with half the blocks:           {:.3e} ({} merged blocks; should agree with the full-block "
                       "error{})\n",
               profile.barrierProbabilityErrorHalfBlocks, profile.numberOfHalfBlocks,
               profile.numberOfHalfBlocks < 3 ? ", but is itself uncertain with fewer than 3 merged blocks" : "");
    if (profile.fitValid)
    {
      std::print(stream, "    fitted maximum offset from the plane: {:.4f} +/- {:.4f} [A]{}\n", profile.barrierOffset,
                 profile.barrierOffsetError,
                 (profile.curvature > 0.0)
                     ? "  WARNING: the fit is a minimum, not a maximum (window site or wrong plane)"
                     : "");
      std::print(stream, "    curvature beta F'' at the maximum:    {:.4f} [1/A^2] (negative at a barrier)\n",
                 profile.curvature);
    }
    if (profile.pilotBarrierDefined)
    {
      std::print(stream, "    P(q*) from the pilot histogram alone: {:.6e} +/- {:.3e} [1/A]\n",
                 profile.pilotBarrierProbability, profile.pilotBarrierProbabilityError);
    }
    else
    {
      std::print(stream, "    P(q*) from the pilot histogram alone: not sampled\n");
    }
    std::print(stream, "\n");

    const double kTST = 0.5 * meanSpeed * profile.barrierProbability;
    const double kTSTError = 0.5 * meanSpeed * profile.barrierProbabilityError;
    std::print(stream, "Transition-state rate:\n");
    std::print(stream, "    <|qdot|> of the reaction bead:       {:.6f} [A/ps] (pilot orientations; sqrt(2 kT / pi M) of the "
                       "molecular mass: {:.6f})\n",
               meanSpeed, thermalSpeed);
    std::print(stream, "    k^TST = <|qdot|>/2 P(q*):            {:.6e} +/- {:.2e} [1/ps]\n", kTST, kTSTError);
    std::print(stream, "    k^TST lambda^2 (D for kappa = 1):    {:.6e} +/- {:.2e} [m^2/s]\n",
               kTST * lambdaSquared * conversion, kTSTError * lambdaSquared * conversion);
    std::print(stream,
               "    (the transmission driver combines P(q*) with its own measured <|qdot|> of the reaction bead)\n\n");

    nlohmann::json& json = outputJson["results"];
    json["barrierProbability"] = profile.barrierProbability;
    json["barrierProbabilityError"] = profile.barrierProbabilityError;
    json["barrierProbabilityErrorHalfBlocks"] = profile.barrierProbabilityErrorHalfBlocks;
    json["barrierProbabilityBin"] = profile.barrierProbabilityBin;
    json["barrierProbabilityBinError"] = profile.barrierProbabilityBinError;
    json["barrierFreeEnergy"] = profile.barrierFreeEnergy;
    json["barrierFreeEnergyError"] = profile.barrierFreeEnergyError;
    json["barrierFreeEnergyBin"] = profile.barrierFreeEnergyBin;
    json["barrierOffset"] = profile.barrierOffset;
    json["barrierOffsetError"] = profile.barrierOffsetError;
    json["curvature"] = profile.curvature;
    json["fitValid"] = profile.fitValid;
    json["fitPoints"] = profile.fitPoints;
    json["pilotBarrierProbability"] = profile.pilotBarrierProbability;
    json["pilotBarrierProbabilityError"] = profile.pilotBarrierProbabilityError;
    json["thermalSpeed"] = thermalSpeed;
    json["meanSpeed"] = meanSpeed;
    json["kTST"] = kTST;
    json["kTSTError"] = kTSTError;
    json["whamConverged"] = profile.converged;
    json["whamIterations"] = profile.iterations;
    json["minimumWindowOverlap"] = minimumOverlap;
    json["pilotConsistencyRms"] = profile.pilotConsistencyRms;
    json["pilotConsistencyBins"] = profile.pilotConsistencyBins;
    for (std::size_t i = 0; i < numberOfStates; ++i)
    {
      json["overlapMatrix"].push_back(std::vector<double>(profile.overlapMatrix.begin() + static_cast<std::ptrdiff_t>(i * numberOfStates),
                                                          profile.overlapMatrix.begin() + static_cast<std::ptrdiff_t>((i + 1uz) * numberOfStates)));
    }
    for (std::size_t k = 0; k < windows.size(); ++k)
    {
      nlohmann::json w;
      w["center"] = windows[k].center;
      w["stiffness"] = windows[k].stiffness;
      w["replica"] = windows[k].replica;
      w["source"] = windows[k].source;
      w["mean"] = profile.windowMean[k];
      w["sigma"] = profile.windowSigma[k];
      w["counts"] = profile.windowCounts[k];
      w["meanBias"] = profile.windowMeanBias[k];
      w["overlapWithNext"] = profile.overlap[k];
      w["freeEnergy"] = profile.stateFreeEnergies[k + 1uz];
      json["windows"].push_back(w);
    }
    for (std::size_t i = 0; i < profile.q.size(); ++i)
    {
      json["profile"]["q"].push_back(profile.q[i]);
      json["profile"]["freeEnergy"].push_back(profile.defined[i] ? profile.freeEnergy[i] : -1.0);
      json["profile"]["freeEnergyError"].push_back(profile.freeEnergyError[i]);
    }
    if (occupancyMode)
    {
      json["conditioned"] = {{"cageOccupancies", {settings.cageOccupancies->first, settings.cageOccupancies->second}},
                             {"barrierProbability", profile.barrierProbability},
                             {"barrierProbabilityError", profile.barrierProbabilityError}};
    }
  }

  if (occupancies.available)
  {
    nlohmann::json& json = outputJson["results"]["cageOccupancies"];
    json["meanA"] = occupancies.meanA;
    json["meanB"] = occupancies.meanB;
    json["samples"] = occupancies.total;
    for (std::size_t nA = 0; nA < occupancies.size; ++nA)
    {
      for (std::size_t nB = 0; nB < occupancies.size; ++nB)
      {
        const double p = occupancies.joint[nA * occupancies.size + nB];
        if (p > 0.0) json["joint"][std::format("{},{}", nA, nB)] = p;
      }
    }
    json["samplesAtBarrier"] = occupancies.totalAtBarrier;
    json["barrierHalfWidth"] = DcTSTWindowAccumulator::barrierOccupancyHalfWidth;
    if (occupancies.availableAtBarrier)
    {
      json["jointAtBarrier"] = nlohmann::json::object();
      for (std::size_t nA = 0; nA < occupancies.size; ++nA)
      {
        for (std::size_t nB = 0; nB < occupancies.size; ++nB)
        {
          const double p = occupancies.jointAtBarrier[nA * occupancies.size + nB];
          if (p > 0.0) json["jointAtBarrier"][std::format("{},{}", nA, nB)] = p;
        }
      }
    }
    if (occupancyMode)
    {
      json["requested"] = {settings.cageOccupancies->first, settings.cageOccupancies->second};
      json["requestedProbability"] = occupancies.requestedProbability;
      json["requestedProbabilityError"] = occupancies.requestedProbabilityError;
    }
  }

  std::print(stream, "Timings\n");
  std::print(stream, "    initialization:                      {:.3f} [s]\n", totalInitializationTime.count());
  std::print(stream, "    pilot:                               {:.3f} [s]\n", totalPilotTime.count());
  std::print(stream, "    windows:                             {:.3f} [s]\n", totalWindowTime.count());
  std::print(stream, "    total:                               {:.3f} [s]\n", totalSimulationTime.count());
  std::flush(stream);

  outputJson["timings"]["total"] = totalSimulationTime.count();
  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}
