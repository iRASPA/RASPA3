#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import intra_molecular_exclusions;
import molecule;
import simd_quatd;
import system;
import simulationbox;
import running_energy;
import energy_status;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import bond_bond_potential;
import integrators_update;
import integrators_compute;
import thermostat;
import randomnumbers;
import spatial_decomposition_settings;
import spatial_decomposition_cell_list;
import spatial_decomposition_pppm;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_worker_team;
import spatial_decomposition_force_engine;
import force_engine;
import spatial_decomposition_device_step;

namespace
{
ForceField makeWaterForceField(ForceField::ChargeMethod chargeMethod = ForceField::ChargeMethod::Ewald)
{
  ForceField forceField = ForceField(
      {PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true), PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true)},
      {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot, 9.0,
      9.0, 9.0, false, false, true);
  forceField.chargeMethod = chargeMethod;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(14, 14, 14);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 196;
  return forceField;
}

Component makeWater(const ForceField& forceField)
{
  return Component(forceField, "H2O", 304.1282, 7377300.0, 0.22394,
                   {Atom(double3(0.00000, -0.06461, 0.00000), -0.84760, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false),
                    Atom(double3(-0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false)},
                   {}, {}, 5, 21);
}

// Places the molecules of the system on a jittered lattice with random orientations (no overlaps) and gives the
// atoms of a flexible component small random displacements around their component geometry.
void randomizeConfiguration(System& system, RandomNumber& random)
{
  const std::size_t numberOfMolecules = system.moleculeData.size();
  const std::size_t perAxis = static_cast<std::size_t>(std::ceil(std::cbrt(static_cast<double>(numberOfMolecules))));
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::size_t index = 0;
  for (std::size_t m = 0; m < numberOfMolecules; ++m)
  {
    const std::size_t ix = index % perAxis;
    const std::size_t iy = (index / perAxis) % perAxis;
    const std::size_t iz = index / (perAxis * perAxis);
    ++index;
    const double3 fractional(
        (static_cast<double>(ix) + 0.5 + 0.3 * (random.uniform() - 0.5)) / static_cast<double>(perAxis),
        (static_cast<double>(iy) + 0.5 + 0.3 * (random.uniform() - 0.5)) / static_cast<double>(perAxis),
        (static_cast<double>(iz) + 0.5 + 0.3 * (random.uniform() - 0.5)) / static_cast<double>(perAxis));
    // an offset so that some molecules straddle the box boundary
    const double3 center = system.simulationBox.cell * (fractional - double3(0.3, 0.2, 0.1));

    Molecule& molecule = system.moleculeData[m];
    const Component& component = system.components[molecule.componentId];
    molecule.centerOfMassPosition = center;
    molecule.orientation = random.randomSimdQuatd();
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
    for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
    {
      double3 local = component.atoms[k].position;
      if (!component.rigid)
      {
        local += 0.1 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
      }
      atoms[molecule.atomIndex + k].position = center + rotation * local;
    }
  }
}

RunningEnergy referenceGradients(System& system)
{
  return Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
}

std::vector<double3> gradientsOf(const System& system)
{
  std::vector<double3> result;
  for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) result.push_back(dynamics.gradient);
  return result;
}

double rmsDifference(const std::vector<double3>& a, const std::vector<double3>& b)
{
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += double3::dot(a[i] - b[i], a[i] - b[i]);
  return std::sqrt(sum / static_cast<double>(a.size()));
}

double rmsNorm(const std::vector<double3>& a)
{
  double sum = 0.0;
  for (const double3& g : a) sum += double3::dot(g, g);
  return std::sqrt(sum / static_cast<double>(a.size()));
}

double maxAbsDifference(const double3x3& a, const double3x3& b)
{
  return std::max({std::abs(a.ax - b.ax), std::abs(a.ay - b.ay), std::abs(a.az - b.az), std::abs(a.bx - b.bx),
                   std::abs(a.by - b.by), std::abs(a.bz - b.bz), std::abs(a.cx - b.cx), std::abs(a.cy - b.cy),
                   std::abs(a.cz - b.cz)});
}

double maxAbs(const double3x3& a)
{
  return std::max({std::abs(a.ax), std::abs(a.ay), std::abs(a.az), std::abs(a.bx), std::abs(a.by), std::abs(a.bz),
                   std::abs(a.cx), std::abs(a.cy), std::abs(a.cz)});
}

SpatialDecompositionSettings settingsFor(std::size_t threads, double skin = 1.0, double meshSpacing = 0.5)
{
  SpatialDecompositionSettings settings;
  settings.numberOfThreads = threads;
  settings.verletSkin = skin;
  settings.meshSpacing = meshSpacing;
  settings.interpolationOrder = 5;
  return settings;
}
}  // namespace

TEST(spatial_decomposition, bspline_weights_partition_of_unity)
{
  for (std::size_t order = 3; order <= 7; ++order)
  {
    for (double w : {0.0, 0.1, 0.37, 0.5, 0.83, 0.999})
    {
      std::array<double, 8> weights{}, derivatives{};
      PPPM::bsplineWeights(order, w, std::span<double>(weights.data(), order),
                           std::span<double>(derivatives.data(), order));
      double sumW = 0.0, sumD = 0.0;
      for (std::size_t j = 0; j < order; ++j)
      {
        EXPECT_GE(weights[j], -1e-15);
        sumW += weights[j];
        sumD += derivatives[j];
      }
      EXPECT_NEAR(sumW, 1.0, 1e-12);
      EXPECT_NEAR(sumD, 0.0, 1e-12);

      // derivative by finite differences
      const double h = 1e-6;
      std::array<double, 8> plus{}, minus{}, dummy{};
      if (w > h && w < 1.0 - h)
      {
        PPPM::bsplineWeights(order, w + h, std::span<double>(plus.data(), order),
                             std::span<double>(dummy.data(), order));
        PPPM::bsplineWeights(order, w - h, std::span<double>(minus.data(), order),
                             std::span<double>(dummy.data(), order));
        for (std::size_t j = 0; j < order; ++j)
        {
          EXPECT_NEAR(derivatives[j], (plus[j] - minus[j]) / (2.0 * h), 1e-6);
        }
      }
    }
  }
  EXPECT_EQ(PPPM::nextFFTFriendly(31), 32u);
  EXPECT_EQ(PPPM::nextFFTFriendly(33), 36u);
  EXPECT_EQ(PPPM::nextFFTFriendly(97), 100u);
}

TEST(spatial_decomposition, ewald_real_space_table_matches_erfc)
{
  const double alpha = 0.265;
  const double cutoff = 12.0;
  EwaldRealSpaceTable table;
  table.build(alpha, cutoff);
  double maximumValueError = 0.0;
  double maximumDerivativeError = 0.0;
  for (std::size_t k = 0; k <= 20000; ++k)
  {
    const double r = 0.3 + (cutoff - 0.3) * static_cast<double>(k) / 20000.0;
    const double rr = r * r;
    double value, derivative, exactValue, exactDerivative;
    table.evaluate(rr, value, derivative);
    table.exact(rr, exactValue, exactDerivative);
    maximumValueError = std::max(maximumValueError, std::abs(value - exactValue) / std::abs(exactValue));
    maximumDerivativeError =
        std::max(maximumDerivativeError, std::abs(derivative - exactDerivative) / std::abs(exactDerivative));
  }
  EXPECT_LT(maximumValueError, 1e-9);
  EXPECT_LT(maximumDerivativeError, 1e-7);
}

TEST(spatial_decomposition, worker_team_runs_phases_and_propagates_exceptions)
{
  WorkerTeam team(4);
  std::vector<std::size_t> counter(4, 0);
  std::atomic<std::size_t> total{0};
  team.run(
      [&](std::size_t thread)
      {
        team.phase([&] { counter[thread] += 1; });
        team.phase([&] { total.fetch_add(counter[thread]); });
      });
  EXPECT_EQ(total.load(), 4u);

  EXPECT_THROW(team.run(
                   [&](std::size_t thread)
                   {
                     team.phase(
                         [&]
                         {
                           if (thread == 2) throw std::runtime_error("boom");
                         });
                     team.phase([&] { counter[thread] += 1; });
                   }),
               std::runtime_error);

  // the team is usable again afterwards
  total = 0;
  team.run([&](std::size_t) { total.fetch_add(1); });
  EXPECT_EQ(total.load(), 4u);

  WorkerTeam single(1);
  std::size_t calls = 0;
  single.run(
      [&](std::size_t thread)
      {
        EXPECT_EQ(thread, 0u);
        single.phase([&] { ++calls; });
        single.phase([&] { ++calls; });
      });
  EXPECT_EQ(calls, 2u);
}

namespace
{
// Checks the per-domain neighbour lists of a cell list against the brute-force pair set of a box: every pair once,
// consistent ghost images and shifts, balanced ownership. Atoms are placed outside the box on purpose (unwrapped
// positions, several box translations in the shifts).
// `ambiguous` receives the number of ambiguous stencil entries (neighbour cells reached through two offsets).
void checkCellListAgainstBruteForce(const SimulationBox& box, std::size_t numberOfAtoms, double cutoff, double skin,
                                    std::size_t& ambiguous)
{
  ambiguous = 0;
  RandomNumber random(12345);
  std::vector<Atom> atoms;
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    const double3 fractional(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
    Atom atom(box.cell * fractional, 0.0, 1.0, static_cast<std::uint32_t>(i / 3), 0, 0, false, false);
    atoms.push_back(atom);
  }

  const double listCutoffSquared = (cutoff + skin) * (cutoff + skin);

  std::set<std::pair<std::uint32_t, std::uint32_t>> bruteForce;
  for (std::uint32_t i = 0; i < atoms.size(); ++i)
  {
    for (std::uint32_t j = i + 1; j < atoms.size(); ++j)
    {
      if (atoms[i].moleculeId == atoms[j].moleculeId) continue;
      double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
      if (double3::dot(dr, dr) < listCutoffSquared) bruteForce.insert({i, j});
    }
  }

  for (std::size_t threads : {1uz, 2uz, 4uz, 6uz, 8uz, 12uz})
  {
    CellList cellList;
    cellList.setup(box, cutoff, skin, threads, std::nullopt);
    cellList.bin(box, atoms);
    for (std::size_t d = 0; d < threads; ++d) cellList.buildLists(d, box);
    ambiguous = 0;
    for (const CellList::StencilEntry& entry : cellList.stencilList) ambiguous += entry.ambiguous;

    // every pair of the system is stored exactly once, as an owned neighbour or as a ghost image whose periodic
    // shift places it within the list cutoff of the owned atom (no minimum-image operation in the kernel)
    std::set<std::pair<std::uint32_t, std::uint32_t>> all;
    std::size_t imagePairs = 0;
    std::size_t smallest = std::numeric_limits<std::size_t>::max();
    std::size_t largest = 0;
    for (std::size_t d = 0; d < threads; ++d)
    {
      const CellList::DomainLists& domain = cellList.domains[d];
      const std::size_t owned = domain.ownedAtoms.size();
      smallest = std::min(smallest, owned);
      largest = std::max(largest, owned);
      ASSERT_EQ(domain.neighbourStart.size(), owned + 1);
      ASSERT_EQ(domain.imageAtom.size(), domain.imageOwner.size());
      ASSERT_EQ(domain.imageAtom.size(), domain.imageShift.size());
      ASSERT_EQ(domain.positions.size(), domain.numberOfLocalAtoms());
      for (std::size_t k = 0; k < owned; ++k)
      {
        const std::uint32_t sortedI = domain.ownedAtoms[k];
        const std::uint32_t i = cellList.sortedToOriginal[sortedI];
        EXPECT_EQ(cellList.ownerOfAtom[sortedI], d);
        const double3 ri = atoms[i].position;
        EXPECT_EQ(domain.positions[k].x, ri.x);
        for (std::uint32_t n = domain.neighbourStart[k]; n < domain.neighbourStart[k + 1]; ++n)
        {
          const std::uint32_t local = domain.neighbourList[n];
          ASSERT_LT(local, domain.numberOfLocalAtoms());
          std::uint32_t sortedJ;
          double3 rj;
          if (local < owned)
          {
            sortedJ = domain.ownedAtoms[local];
            rj = atoms[cellList.sortedToOriginal[sortedJ]].position;
          }
          else
          {
            const std::uint32_t slot = local - static_cast<std::uint32_t>(owned);
            sortedJ = domain.imageAtom[slot];
            EXPECT_EQ(domain.imageOwner[slot], cellList.ownerOfAtom[sortedJ]);
            rj = atoms[cellList.sortedToOriginal[sortedJ]].position +
                 CellList::shiftVector(box, domain.imageShift[slot]);
            ++imagePairs;
          }
          // the gathered image position is the shifted position, and it is the interacting image
          const double3 gathered(domain.positions[local].x, domain.positions[local].y, domain.positions[local].z);
          EXPECT_LT((gathered - rj).length(), 1e-9);
          const double3 dr = ri - rj;
          EXPECT_LT(double3::dot(dr, dr), listCutoffSquared);
          const std::uint32_t j = cellList.sortedToOriginal[sortedJ];
          EXPECT_NE(i, j);
          const auto pair = std::make_pair(std::min(i, j), std::max(i, j));
          EXPECT_TRUE(all.insert(pair).second) << "pair stored twice";
        }
      }
      // the reduction table lists every image slot exactly once, under its owner
      std::size_t listedSlots = 0;
      ASSERT_EQ(domain.imagesByOwner.size(), threads);
      for (std::size_t owner = 0; owner < threads; ++owner)
      {
        for (const std::uint32_t slot : domain.imagesByOwner[owner])
        {
          EXPECT_EQ(domain.imageOwner[slot], owner);
          ++listedSlots;
        }
      }
      EXPECT_EQ(listedSlots, domain.imageAtom.size());
    }
    if (threads > 1) EXPECT_GT(imagePairs, 0u);
    EXPECT_EQ(cellList.totalPairs(), bruteForce.size());
    EXPECT_EQ(all, bruteForce) << "threads " << threads;
    // the sub-domain cuts balance the atom count
    EXPECT_LE(largest - smallest, 4u) << "threads " << threads;
  }
}
}  // namespace

TEST(spatial_decomposition, cell_list_matches_brute_force_pairs_triclinic)
{
  const SimulationBox box(28.0, 31.0, 26.0, 100.0 * (std::numbers::pi / 180.0), 95.0 * (std::numbers::pi / 180.0),
                          75.0 * (std::numbers::pi / 180.0));
  std::size_t ambiguous = 0;
  checkCellListAgainstBruteForce(box, 1500, 8.0, 1.0, ambiguous);
  EXPECT_EQ(ambiguous, 0u);
}

TEST(spatial_decomposition, cell_list_matches_brute_force_pairs_small_box)
{
  // cutoff + skin exactly half the box: 6 cells per axis with a 7-wide stencil, so neighbour cells are reached
  // through two offsets (the minimum-image fallback of the cell-pair build)
  const SimulationBox box(18.0, 18.0, 18.0);
  std::size_t ambiguous = 0;
  checkCellListAgainstBruteForce(box, 600, 8.0, 1.0, ambiguous);
  EXPECT_GT(ambiguous, 0u);
}

TEST(spatial_decomposition, engine_matches_exact_ewald_rigid_water)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(7);
  randomizeConfiguration(system, random);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;

  SpatialDecompositionForceEngine engine(settingsFor(1));
  engine.initialize(system);
  EXPECT_TRUE(engine.usesFastKernel());  // plain Lennard-Jones + Ewald, fully coupled atoms
  const RunningEnergy energy = engine.computeGradients(system, true);
  const std::vector<double3> engineGradient = gradientsOf(system);

  // real-space terms: the specialised kernel with the tabulated erfc agrees with the exact kernels to rounding
  EXPECT_NEAR(energy.moleculeMoleculeVDW, reference.moleculeMoleculeVDW,
              1e-8 * std::abs(reference.moleculeMoleculeVDW));
  EXPECT_NEAR(energy.moleculeMoleculeCharge, reference.moleculeMoleculeCharge,
              1e-8 * std::abs(reference.moleculeMoleculeCharge));
  EXPECT_NEAR(energy.ewald_self, reference.ewald_self, 1e-10 * std::abs(reference.ewald_self));
  EXPECT_NEAR(energy.ewald_exclusion, reference.ewald_exclusion, 1e-10 * std::abs(reference.ewald_exclusion));

  // the mesh approximates the (converged) direct k-space sum
  EXPECT_NEAR(energy.ewald_fourier, reference.ewald_fourier, 1e-5 * std::abs(reference.ewald_fourier));
  EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()));

  const double rms = rmsNorm(referenceGradient);
  EXPECT_LT(rmsDifference(engineGradient, referenceGradient), 1e-4 * rms);

  const double3x3 enginePressure = engine.molecularPressureTensor();
  EXPECT_LT(maxAbsDifference(enginePressure, referencePressure), 1e-4 * maxAbs(referencePressure));
}

TEST(spatial_decomposition, force_engine_is_a_movable_value)
{
  static_assert(std::is_nothrow_move_constructible_v<ForceEngine> && std::is_nothrow_move_assignable_v<ForceEngine>);
  static_assert(!std::is_copy_constructible_v<ForceEngine>);

  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(7);
  randomizeConfiguration(system, random);

  // an initialized engine (worker threads, FFTW plans and meshes, cell list) is moved into a vector that
  // reallocates, as the driver's std::vector<ForceEngine> does; the moved engine must give the same results
  ForceEngine original(settingsFor(4));
  original.initialize(system);
  const RunningEnergy before = original.computeGradients(system, true);
  const std::vector<double3> gradientBefore = gradientsOf(system);
  const double3x3 pressureBefore = original.molecularPressureTensor();

  std::vector<ForceEngine> engines;
  engines.push_back(std::move(original));
  engines.emplace_back(settingsFor(1));  // reallocation moves the first engine again
  ForceEngine& moved = engines.front();
  EXPECT_TRUE(moved.initialized());
  EXPECT_EQ(moved.numberOfThreads(), 4uz);

  const RunningEnergy after = moved.computeGradients(system, true);
  EXPECT_EQ(after.potentialEnergy(), before.potentialEnergy());
  EXPECT_EQ(after.ewald_fourier, before.ewald_fourier);
  EXPECT_LT(rmsDifference(gradientsOf(system), gradientBefore), 1e-12 * rmsNorm(gradientBefore));
  EXPECT_LT(maxAbsDifference(moved.molecularPressureTensor(), pressureBefore), 1e-12 * maxAbs(pressureBefore));
  EXPECT_EQ(moved.timings().steps, 2uz);

  ForceEngine assigned(settingsFor(1));
  assigned = std::move(engines.front());
  EXPECT_EQ(assigned.computeGradients(system, true).potentialEnergy(), before.potentialEnergy());
}

TEST(spatial_decomposition, engine_threads_agree_with_serial_rigid_water)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(11);
  randomizeConfiguration(system, random);

  SpatialDecompositionForceEngine serial(settingsFor(1));
  serial.initialize(system);
  const RunningEnergy serialEnergy = serial.computeGradients(system, true);
  const std::vector<double3> serialGradient = gradientsOf(system);
  const double3x3 serialPressure = serial.molecularPressureTensor();
  const double rms = rmsNorm(serialGradient);

  std::vector<double3> originalPositions;
  for (const Atom& atom : system.spanOfMoleculeAtoms()) originalPositions.push_back(atom.position);

  for (std::size_t threads : {2uz, 3uz, 4uz, 8uz})
  {
    {
      std::span<Atom> atoms = system.spanOfMoleculeAtoms();
      for (std::size_t i = 0; i < atoms.size(); ++i) atoms[i].position = originalPositions[i];
    }
    SpatialDecompositionForceEngine parallel(settingsFor(threads));
    parallel.initialize(system);
    const RunningEnergy energy = parallel.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    EXPECT_NEAR(energy.moleculeMoleculeVDW, serialEnergy.moleculeMoleculeVDW,
                1e-9 * std::abs(serialEnergy.moleculeMoleculeVDW));
    EXPECT_NEAR(energy.moleculeMoleculeCharge, serialEnergy.moleculeMoleculeCharge,
                1e-9 * std::abs(serialEnergy.moleculeMoleculeCharge));
    EXPECT_NEAR(energy.ewald_fourier, serialEnergy.ewald_fourier, 1e-9 * std::abs(serialEnergy.ewald_fourier));
    EXPECT_NEAR(energy.ewald_exclusion, serialEnergy.ewald_exclusion, 1e-10 * std::abs(serialEnergy.ewald_exclusion));
    EXPECT_NEAR(energy.potentialEnergy(), serialEnergy.potentialEnergy(),
                1e-9 * std::abs(serialEnergy.potentialEnergy()));
    EXPECT_LT(rmsDifference(gradient, serialGradient), 1e-9 * rms) << "threads " << threads;
    EXPECT_LT(maxAbsDifference(parallel.molecularPressureTensor(), serialPressure), 1e-9 * maxAbs(serialPressure));

    // a second evaluation after small displacements reuses the lists (no rebuild) and must still agree
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    for (Atom& atom : atoms)
    {
      atom.position += 0.05 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    }
    const RunningEnergy movedParallel = parallel.computeGradients(system, false);
    const std::vector<double3> movedGradient = gradientsOf(system);
    const RunningEnergy movedSerial = serial.computeGradients(system, false);
    const std::vector<double3> movedSerialGradient = gradientsOf(system);
    EXPECT_NEAR(movedParallel.potentialEnergy(), movedSerial.potentialEnergy(),
                1e-9 * std::abs(movedSerial.potentialEnergy()));
    EXPECT_LT(rmsDifference(movedGradient, movedSerialGradient), 1e-9 * rmsNorm(movedSerialGradient));
    EXPECT_GT(parallel.timings().steps, parallel.timings().rebuilds);
  }
}

namespace
{
struct KernelTolerances
{
  double vdw;       ///< relative, molecule-molecule Lennard-Jones energy
  double charge;    ///< relative, real-space Coulomb energy
  double total;     ///< relative, potential energy
  double gradient;  ///< relative rms
  double pressure;  ///< relative max-abs
  /// relative, Ewald Fourier energy and self + exclusion correction (1e-10 unless the mesh / the molecular terms
  /// are evaluated on the device)
  double fourier{1e-10};
  double correction{1e-10};
};

/// Compares the engine with \p candidate settings against the engine with the default (scalar double) kernel on
/// a water box, at 1 and 4 threads, including a second evaluation after small displacements (lists reused).
void expectKernelAgreesWithScalar(const SpatialDecompositionSettings& candidate, const KernelTolerances& tolerance)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(13);
  randomizeConfiguration(system, random);

  std::vector<double3> originalPositions;
  for (const Atom& atom : system.spanOfMoleculeAtoms()) originalPositions.push_back(atom.position);

  for (std::size_t threads : {1uz, 4uz})
  {
    {
      std::span<Atom> atoms = system.spanOfMoleculeAtoms();
      for (std::size_t i = 0; i < atoms.size(); ++i) atoms[i].position = originalPositions[i];
    }
    SpatialDecompositionForceEngine reference(settingsFor(threads));
    reference.initialize(system);
    EXPECT_TRUE(reference.usesFastKernel());
    const RunningEnergy referenceEnergy = reference.computeGradients(system, true);
    const std::vector<double3> referenceGradient = gradientsOf(system);
    const double3x3 referencePressure = reference.molecularPressureTensor();
    const double rms = rmsNorm(referenceGradient);

    SpatialDecompositionSettings settings = candidate;
    settings.numberOfThreads = threads;
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    EXPECT_TRUE(engine.usesFastKernel());
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);

    EXPECT_NEAR(energy.moleculeMoleculeVDW, referenceEnergy.moleculeMoleculeVDW,
                tolerance.vdw * std::abs(referenceEnergy.moleculeMoleculeVDW))
        << "threads " << threads;
    EXPECT_NEAR(energy.moleculeMoleculeCharge, referenceEnergy.moleculeMoleculeCharge,
                tolerance.charge * std::abs(referenceEnergy.moleculeMoleculeCharge))
        << "threads " << threads;
    // the mesh, self and exclusion terms are evaluated in double by the same code in both engines, unless on the
    // device (the self and exclusion energies cancel strongly: their sum is compared)
    EXPECT_NEAR(energy.ewald_fourier, referenceEnergy.ewald_fourier,
                tolerance.fourier * std::abs(referenceEnergy.ewald_fourier))
        << "threads " << threads;
    EXPECT_NEAR(energy.ewald_self + energy.ewald_exclusion,
                referenceEnergy.ewald_self + referenceEnergy.ewald_exclusion,
                tolerance.correction * std::abs(referenceEnergy.ewald_self + referenceEnergy.ewald_exclusion))
        << "threads " << threads;
    EXPECT_NEAR(energy.potentialEnergy(), referenceEnergy.potentialEnergy(),
                tolerance.total * std::abs(referenceEnergy.potentialEnergy()))
        << "threads " << threads;
    EXPECT_LT(rmsDifference(gradient, referenceGradient), tolerance.gradient * rms) << "threads " << threads;
    EXPECT_LT(maxAbsDifference(engine.molecularPressureTensor(), referencePressure),
              tolerance.pressure * maxAbs(referencePressure))
        << "threads " << threads;

    // after small displacements (lists reused, positions refreshed relative to the sub-domain origin)
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    for (Atom& atom : atoms)
    {
      atom.position += 0.05 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    }
    const RunningEnergy movedEngine = engine.computeGradients(system, false);
    const std::vector<double3> movedGradient = gradientsOf(system);
    const RunningEnergy movedReference = reference.computeGradients(system, false);
    const std::vector<double3> movedReferenceGradient = gradientsOf(system);
    EXPECT_NEAR(movedEngine.potentialEnergy(), movedReference.potentialEnergy(),
                tolerance.total * std::abs(movedReference.potentialEnergy()))
        << "threads " << threads;
    EXPECT_LT(rmsDifference(movedGradient, movedReferenceGradient),
              tolerance.gradient * rmsNorm(movedReferenceGradient))
        << "threads " << threads;
    EXPECT_GT(engine.timings().steps, engine.timings().rebuilds);

    // larger displacements: beyond half the prune skin (the pruned list is rebuilt from the Verlet list) but
    // within half the Verlet skin (no rebuild)
    for (Atom& atom : atoms)
    {
      atom.position += 0.4 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    }
    const RunningEnergy prunedEngine = engine.computeGradients(system, false);
    const std::vector<double3> prunedGradient = gradientsOf(system);
    const RunningEnergy prunedReference = reference.computeGradients(system, false);
    const std::vector<double3> prunedReferenceGradient = gradientsOf(system);
    EXPECT_NEAR(prunedEngine.potentialEnergy(), prunedReference.potentialEnergy(),
                tolerance.total * std::abs(prunedReference.potentialEnergy()))
        << "threads " << threads;
    EXPECT_LT(rmsDifference(prunedGradient, prunedReferenceGradient),
              tolerance.gradient * rmsNorm(prunedReferenceGradient))
        << "threads " << threads;
    EXPECT_EQ(engine.timings().rebuilds, 1uz);
  }
}
}  // namespace

TEST(spatial_decomposition, engine_cluster_kernel_double_matches_scalar_rigid_water)
{
  // the same arithmetic in a different order: agreement to rounding
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.clusterKernelForDouble = true;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-10, .charge = 1e-9, .total = 1e-9, .gradient = 1e-9, .pressure = 1e-9});
}

TEST(spatial_decomposition, engine_mixed_precision_agrees_with_double_rigid_water)
{
  // Single-precision pair geometry (positions within ~25 Angstrom of the sub-domain origin, so about 2e-6
  // Angstrom) and closed-form erfc (absolute error 1.5e-7), double accumulation. Measured: Lennard-Jones energy
  // 2e-8, real-space Coulomb energy 3e-5 (the erfc error is systematic and the water charges cancel strongly),
  // total energy 3e-6, gradient rms 8e-7, pressure 1.4e-6 relative.
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairPrecision = PairPrecision::Mixed;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-6, .charge = 1e-4, .total = 1e-5, .gradient = 1e-5, .pressure = 1e-5});
}

/// The device backends: every device test runs on each available one (skipped where the device is missing).
class SpatialDecompositionDevice : public testing::TestWithParam<PairDevice>
{
};
INSTANTIATE_TEST_SUITE_P(spatial_decomposition, SpatialDecompositionDevice,
                         testing::Values(PairDevice::OpenCL, PairDevice::Metal, PairDevice::CUDA),
                         [](const testing::TestParamInfo<PairDevice>& info) { return pairDeviceName(info.param); });

TEST_P(SpatialDecompositionDevice, engine_device_pair_kernel_agrees_with_double_rigid_water)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // the device kernel: single-precision pair geometry with the minimum image of wrapped positions, closed-form
  // erfc, single-precision accumulation of the per-atom forces and per-cluster partial sums
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = pairDevice;
  settings.deviceMesh = false;
  settings.deviceBonded = false;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-6, .charge = 1e-4, .total = 1e-5, .gradient = 1e-5, .pressure = 1e-5});
}

TEST_P(SpatialDecompositionDevice, engine_device_pair_kernel_without_pruning_rigid_water)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // no pruning: the lane lists hold the whole outer list (cutoff + Verlet skin), compacted once per build
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = pairDevice;
  settings.pruneSkin = 0.0;
  settings.deviceMesh = false;
  settings.deviceBonded = false;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-6, .charge = 1e-4, .total = 1e-5, .gradient = 1e-5, .pressure = 1e-5});
}

TEST_P(SpatialDecompositionDevice, engine_device_mesh_and_molecular_terms_agree_with_double_rigid_water)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // the complete device step: pairs, PPPM (fixed-point spreading, single-precision FFT and influence function,
  // gather interpolation) and the molecular terms (self + exclusion without their cancellation, virial correction).
  // Measured: Fourier energy 1e-6, self + exclusion sum 2e-7, total energy 3.4e-6, gradient rms 2.5e-6, pressure
  // 1e-6 relative.
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = pairDevice;
  expectKernelAgreesWithScalar(settings, {.vdw = 1e-6,
                                          .charge = 1e-4,
                                          .total = 1e-5,
                                          .gradient = 1e-5,
                                          .pressure = 1e-5,
                                          .fourier = 1e-5,
                                          .correction = 1e-5});
}

TEST_P(SpatialDecompositionDevice, engine_device_pair_kernel_small_grid_rigid_water)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // cutoff + skin = half the box: the device grid has 2 cells per axis, every neighbour cell is reached through
  // several stencil offsets and the list build takes the minimum image (the outer list then holds every pair)
  SpatialDecompositionSettings settings = settingsFor(1, 6.0);
  settings.pairDevice = pairDevice;
  expectKernelAgreesWithScalar(settings, {.vdw = 1e-6,
                                          .charge = 1e-4,
                                          .total = 1e-5,
                                          .gradient = 1e-5,
                                          .pressure = 1e-5,
                                          .fourier = 1e-5,
                                          .correction = 1e-5});
}

TEST(spatial_decomposition, engine_matches_exact_ewald_flexible_chains_triclinic)
{
  ForceField forceField = ForceField(
      {{"CH2", false, 14.02658, 0.0, 0.0, 8, false}, {"O", false, 15.9994, 0.0, 0.0, 8, false}},
      {{56.0, 3.96}, {80.0, 3.1}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(14, 14, 14);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 196;

  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                                    BondPotential({1, 2}, BondType::Harmonic, {96500.0, 1.43})};
  intraMolecularPotentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0})};

  Component chain = Component(
      forceField, "chain", 500.0, 4871800.0, 0.0993,
      {Atom({0.0, 0.0, 0.0}, 0.3, 1.0, 0, 0, 0, false, false), Atom({1.54, 0.0, 0.0}, 0.2, 1.0, 0, 0, 0, false, false),
       Atom({2.1, 1.3, 0.0}, -0.5, 1.0, 0, 1, 0, false, false)},
      connectivityTable, intraMolecularPotentials, 5, 21);

  const SimulationBox box(30.0, 31.0, 29.0, 100.0 * (std::numbers::pi / 180.0), 95.0 * (std::numbers::pi / 180.0),
                          75.0 * (std::numbers::pi / 180.0));
  System system = System(forceField, box, false, 300.0, 1e5, 1.0, {}, {chain}, {}, {300}, 5);
  RandomNumber random(3);
  randomizeConfiguration(system, random);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;
  EXPECT_NE(reference.bond, 0.0);
  EXPECT_NE(reference.bend, 0.0);

  SpatialDecompositionForceEngine serial(settingsFor(1, 1.5, 0.5));
  serial.initialize(system);
  const RunningEnergy energy = serial.computeGradients(system, true);
  const std::vector<double3> serialGradient = gradientsOf(system);

  EXPECT_NEAR(energy.bond, reference.bond, 1e-10 * std::abs(reference.bond));
  EXPECT_NEAR(energy.bend, reference.bend, 1e-10 * std::abs(reference.bend));
  EXPECT_NEAR(energy.moleculeMoleculeVDW, reference.moleculeMoleculeVDW,
              1e-8 * std::abs(reference.moleculeMoleculeVDW));
  EXPECT_NEAR(energy.moleculeMoleculeCharge, reference.moleculeMoleculeCharge,
              1e-8 * std::abs(reference.moleculeMoleculeCharge));
  EXPECT_NEAR(energy.ewald_fourier, reference.ewald_fourier, 1e-5 * std::abs(reference.ewald_fourier));
  EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()));
  EXPECT_LT(rmsDifference(serialGradient, referenceGradient), 1e-4 * rmsNorm(referenceGradient));
  EXPECT_LT(maxAbsDifference(serial.molecularPressureTensor(), referencePressure), 1e-4 * maxAbs(referencePressure));

  for (std::size_t threads : {2uz, 4uz, 8uz})
  {
    SpatialDecompositionForceEngine parallel(settingsFor(threads, 1.5, 0.5));
    parallel.initialize(system);
    const RunningEnergy parallelEnergy = parallel.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    EXPECT_NEAR(parallelEnergy.potentialEnergy(), energy.potentialEnergy(), 1e-9 * std::abs(energy.potentialEnergy()));
    EXPECT_NEAR(parallelEnergy.bond, energy.bond, 1e-10 * std::abs(energy.bond));
    EXPECT_LT(rmsDifference(gradient, serialGradient), 1e-9 * rmsNorm(serialGradient)) << "threads " << threads;
    EXPECT_LT(maxAbsDifference(parallel.molecularPressureTensor(), serial.molecularPressureTensor()),
              1e-9 * maxAbs(serial.molecularPressureTensor()));
  }

  // the device pairs in the triclinic cell (general minimum image in the list build, pruning and pair kernel)
  for (const PairDevice pairDevice : {PairDevice::OpenCL, PairDevice::Metal, PairDevice::CUDA})
  {
    if (!DeviceStep::available(pairDevice)) continue;
    SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
    settings.pairDevice = pairDevice;
    SpatialDecompositionForceEngine device(settings);
    device.initialize(system);
    const RunningEnergy deviceEnergy = device.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    // the random configuration overlaps: the VDW energy is ~1e9 and the per-lane single-precision accumulators
    // of the pair kernel round at ~100 per addition, so the backends differ from the host (and from each other,
    // by their contraction patterns) at the 1e-6 level
    EXPECT_NEAR(deviceEnergy.moleculeMoleculeVDW, energy.moleculeMoleculeVDW,
                5e-6 * std::abs(energy.moleculeMoleculeVDW));
    EXPECT_NEAR(deviceEnergy.potentialEnergy(), energy.potentialEnergy(), 1e-5 * std::abs(energy.potentialEnergy()));
    // the stiff bonds on the device: single precision of the positions relative to the molecule
    EXPECT_NEAR(deviceEnergy.bond, energy.bond, 1e-5 * std::abs(energy.bond));
    EXPECT_NEAR(deviceEnergy.bend, energy.bend, 1e-5 * std::abs(energy.bend));
    EXPECT_LT(rmsDifference(gradient, serialGradient), 1e-5 * rmsNorm(serialGradient));
    EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), serial.molecularPressureTensor()),
              1e-5 * maxAbs(serial.molecularPressureTensor()));
  }
}

namespace
{
/// A charged chain of 'numberOfBeads' beads with harmonic bonds and bends, TraPPE torsions, intramolecular
/// Lennard-Jones and Coulomb pairs (the first 1-4 pair scaled by 0.5, the other 1-4 pairs excluded, everything
/// beyond 1-4 at full strength) and, optionally, a bond-bond cross term (which the device kernels do not cover).
System makeChainSystem(bool withBondBond, RandomNumber& random, std::size_t numberOfBeads = 4)
{
  ForceField forceField = ForceField(
      {{"CH3", false, 15.03452, 0.0, 0.0, 6, false}, {"CH2", false, 14.02658, 0.0, 0.0, 6, false}},
      {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  // the 40 Angstrom box needs more wave vectors than the 30 Angstrom boxes above for a converged reference
  // (exp(-k^2 / 4 alpha^2) ~ 1e-5 at k = 2 pi 14 / 40 with alpha = 0.32, ~1e-10 at 20 wave vectors)
  forceField.numberOfWaveVectors = int3(20, 20, 20);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 400;

  ConnectivityTable connectivityTable(numberOfBeads);
  for (std::size_t i = 0; i + 1 < numberOfBeads; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }
  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i + 1 < numberOfBeads; ++i)
  {
    potentials.bonds.push_back(BondPotential({i, i + 1}, BondType::Harmonic, {96500.0, 1.54}));
  }
  for (std::size_t i = 0; i + 2 < numberOfBeads; ++i)
  {
    potentials.bends.push_back(BendPotential({i, i + 1, i + 2}, BendType::Harmonic, {62500.0, 114.0}));
  }
  for (std::size_t i = 0; i + 3 < numberOfBeads; ++i)
  {
    potentials.torsions.push_back(
        TorsionPotential({i, i + 1, i + 2, i + 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32}));
  }
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VDWParameters::Type::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.25, 0.25, 0.5)};
  if (withBondBond)
  {
    potentials.bondBonds = {BondBondPotential({0, 1, 2}, BondBondType::CFF, {5000.0, 1.54, 1.54})};
  }

  // an all-trans zig-zag in the xy-plane; the end beads are CH3 (+0.25), the inner beads CH2 (-0.25)
  std::vector<Atom> beads{};
  for (std::size_t i = 0; i < numberOfBeads; ++i)
  {
    const bool end = (i == 0 || i + 1 == numberOfBeads);
    const double x = 1.29 * static_cast<double>(i) - 0.645 * static_cast<double>(numberOfBeads - 1);
    const double y = (i % 2 == 0) ? -0.42 : 0.42;
    beads.push_back(Atom({x, y, 0.0}, end ? 0.25 : -0.25, 1.0, 0, end ? 0 : 1, 0, false, false));
  }
  Component chain =
      Component(forceField, "chain", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  // 4 x 4 x 4 molecules on a lattice of about 10 Angstrom: no close contacts (the fallback test compares the
  // device pairs with the host pairs, so the pair energies must not be dominated by overlaps)
  System system = System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain}, {}, {64}, 5);
  randomizeConfiguration(system, random);
  return system;
}

// A compact chain of n^3 beads on a boustrophedon path through an n x n x n lattice (3.8 Angstrom spacing, like a
// coarse-grained protein backbone), with harmonic bonds and bends, TraPPE torsions, alternating charges and the
// 1-4 pairs scaled by one half: a single molecule with far more atoms than the 256 of the old device slot packing.
System makeLatticeChainSystem(RandomNumber& random, std::size_t n, std::size_t numberOfMolecules)
{
  ForceField forceField = ForceField(
      {{"CH3", false, 15.03452, 0.0, 0.0, 6, false}, {"CH2", false, 14.02658, 0.0, 0.0, 6, false}},
      {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(20, 20, 20);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 400;

  const std::size_t numberOfBeads = n * n * n;
  ConnectivityTable connectivityTable(numberOfBeads);
  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i + 1 < numberOfBeads; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
    potentials.bonds.push_back(BondPotential({i, i + 1}, BondType::Harmonic, {96500.0, 3.8}));
  }
  for (std::size_t i = 0; i + 2 < numberOfBeads; ++i)
  {
    potentials.bends.push_back(BendPotential({i, i + 1, i + 2}, BendType::Harmonic, {62500.0, 100.0}));
  }
  for (std::size_t i = 0; i + 3 < numberOfBeads; ++i)
  {
    potentials.torsions.push_back(
        TorsionPotential({i, i + 1, i + 2, i + 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32}));
    // the 1-4 pairs at half strength (the component's default 1-4 scaling is zero)
    potentials.vanDerWaals.push_back(
        VanDerWaalsPotential({i, i + 3}, VDWParameters::Type::LennardJones, {46.0, 3.95}, 0.5));
    potentials.coulombs.push_back(CoulombPotential({i, i + 3}, CoulombType::Coulomb, 0.25, 0.25, 0.5));
  }

  std::vector<Atom> beads{};
  const double spacing = 3.8;
  const double offset = 0.5 * spacing * static_cast<double>(n - 1);
  for (std::size_t z = 0; z < n; ++z)
  {
    for (std::size_t row = 0; row < n; ++row)
    {
      const std::size_t y = (z % 2 == 0) ? row : n - 1 - row;
      for (std::size_t column = 0; column < n; ++column)
      {
        const std::size_t x = ((z * n + row) % 2 == 0) ? column : n - 1 - column;
        const std::size_t i = beads.size();
        const bool end = (i == 0 || i + 1 == numberOfBeads);
        beads.push_back(Atom({spacing * static_cast<double>(x) - offset, spacing * static_cast<double>(y) - offset,
                              spacing * static_cast<double>(z) - offset},
                             (i % 2 == 0) ? 0.25 : -0.25, 1.0, 0, end ? 0 : 1, 0, false, false));
      }
    }
  }
  Component chain =
      Component(forceField, "lattice-chain", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  // the molecules are inserted at the component geometry (a CBMC growth of a 343-bead chain would take minutes)
  // and placed by randomizeConfiguration
  std::vector<double3> positions{};
  for (std::size_t m = 0; m < numberOfMolecules; ++m)
  {
    for (const Atom& bead : beads) positions.push_back(bead.position);
  }
  System system = System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain},
                         {positions}, {0}, 5);
  randomizeConfiguration(system, random);
  return system;
}
}  // namespace

// The same-molecule pairs of a flexible chain in the cell lists: a six-bead chain has excluded 1-2 / 1-3 pairs, a
// scaled 1-4 pair, excluded (scaling 0) 1-4 pairs and full-strength 1-5 / 1-6 pairs. The lists must hold exactly
// the full-strength same-molecule pairs within the list cutoff, and the engine (all host kernels) must reproduce
// the exact (all-pairs Ewald) energies, gradients and molecular pressure.
TEST(spatial_decomposition, cell_list_lists_the_non_excluded_same_molecule_pairs)
{
  RandomNumber random(17);
  System system = makeChainSystem(false, random, 6);
  const Component& chain = system.components[0];
  ASSERT_EQ(chain.intraMolecularPotentials.exclusions.pairs.size(), 9uz);       // 5 bonds + 4 bends
  ASSERT_EQ(chain.intraMolecularPotentials.exclusions.scaledPairs.size(), 3uz);  // the 1-4 pairs
  EXPECT_TRUE(chain.intraMolecularPotentials.exclusions.isExcludedFromPairList(0, 3));
  EXPECT_TRUE(chain.intraMolecularPotentials.exclusions.isExcludedFromPairList(1, 4));
  EXPECT_FALSE(chain.intraMolecularPotentials.exclusions.isExcludedFromPairList(0, 4));
  EXPECT_FALSE(chain.intraMolecularPotentials.exclusions.isExcludedFromPairList(0, 5));
  EXPECT_FALSE(chain.intraMolecularPotentials.exclusions.isExcluded(0, 3));

  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  const double cutoff = 9.0, skin = 1.5;
  const double listCutoffSquared = (cutoff + skin) * (cutoff + skin);
  for (std::size_t threads : {1uz, 4uz})
  {
    CellList cellList;
    cellList.setup(box, cutoff, skin, threads, std::nullopt);
    cellList.bin(box, atoms, system.components);
    for (std::size_t d = 0; d < threads; ++d) cellList.buildLists(d, box);

    std::set<std::pair<std::uint32_t, std::uint32_t>> sameMolecule;
    for (std::size_t d = 0; d < threads; ++d)
    {
      const CellList::DomainLists& domain = cellList.domains[d];
      const std::size_t owned = domain.ownedAtoms.size();
      for (std::size_t k = 0; k < owned; ++k)
      {
        const std::uint32_t i = cellList.sortedToOriginal[domain.ownedAtoms[k]];
        for (std::uint32_t n = domain.neighbourStart[k]; n < domain.neighbourStart[k + 1]; ++n)
        {
          const std::uint32_t local = domain.neighbourList[n];
          const std::uint32_t sortedJ = local < owned ? domain.ownedAtoms[local] : domain.imageAtom[local - owned];
          const std::uint32_t j = cellList.sortedToOriginal[sortedJ];
          if (atoms[i].moleculeId != atoms[j].moleculeId) continue;
          EXPECT_TRUE(sameMolecule.insert({std::min(i, j), std::max(i, j)}).second);
        }
      }
    }

    std::set<std::pair<std::uint32_t, std::uint32_t>> expected;
    for (const Molecule& molecule : system.moleculeData)
    {
      for (std::uint32_t a = 0; a < molecule.numberOfAtoms; ++a)
      {
        for (std::uint32_t b = a + 1; b < molecule.numberOfAtoms; ++b)
        {
          if (chain.intraMolecularPotentials.exclusions.isExcludedFromPairList(a, b)) continue;
          const std::uint32_t i = static_cast<std::uint32_t>(molecule.atomIndex) + a;
          const std::uint32_t j = static_cast<std::uint32_t>(molecule.atomIndex) + b;
          const double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
          if (double3::dot(dr, dr) < listCutoffSquared) expected.insert({i, j});
        }
      }
    }
    EXPECT_EQ(expected.size(), 3uz * system.moleculeData.size());  // (0,4), (1,5), (0,5) per chain
    EXPECT_EQ(sameMolecule, expected) << "threads " << threads;

    // without components: no same-molecule pair at all
    CellList plain;
    plain.setup(box, cutoff, skin, threads, std::nullopt);
    plain.bin(box, atoms);
    for (std::size_t d = 0; d < threads; ++d) plain.buildLists(d, box);
    EXPECT_EQ(plain.totalPairs() + expected.size(), cellList.totalPairs());
  }
}

TEST(spatial_decomposition, engine_matches_exact_ewald_long_flexible_chains)
{
  RandomNumber random(19);
  System system = makeChainSystem(false, random, 6);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;
  EXPECT_NE(reference.intraVDW, 0.0);
  EXPECT_NE(reference.intraCoul, 0.0);
  EXPECT_NE(reference.torsion, 0.0);

  // scalar double, cluster double, cluster mixed
  std::vector<SpatialDecompositionSettings> candidates;
  for (std::size_t threads : {1uz, 4uz})
  {
    candidates.push_back(settingsFor(threads, 1.5, 0.5));
    SpatialDecompositionSettings cluster = settingsFor(threads, 1.5, 0.5);
    cluster.clusterKernelForDouble = true;
    candidates.push_back(cluster);
    SpatialDecompositionSettings mixed = settingsFor(threads, 1.5, 0.5);
    mixed.pairPrecision = PairPrecision::Mixed;
    candidates.push_back(mixed);
  }
  for (const SpatialDecompositionSettings& settings : candidates)
  {
    const bool mixed = settings.pairPrecision == PairPrecision::Mixed;
    const double tolerance = mixed ? 1e-5 : 1e-8;
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    const std::string label = std::format("threads {} cluster {} mixed {}", settings.numberOfThreads,
                                          settings.clusterKernelForDouble, mixed);

    // the engine evaluates the unscaled same-molecule pairs in its pair loops (booked with the molecule-molecule
    // pairs); only the scaled 1-4 pairs stay in the intramolecular slots
    EXPECT_NEAR(energy.intraVDW + energy.moleculeMoleculeVDW, reference.intraVDW + reference.moleculeMoleculeVDW,
                tolerance * std::abs(reference.intraVDW + reference.moleculeMoleculeVDW))
        << label;
    EXPECT_NEAR(energy.intraCoul + energy.moleculeMoleculeCharge,
                reference.intraCoul + reference.moleculeMoleculeCharge,
                tolerance * std::abs(reference.intraCoul + reference.moleculeMoleculeCharge))
        << label;
    EXPECT_NE(energy.intraVDW, 0.0) << label;
    EXPECT_NE(energy.intraVDW, reference.intraVDW) << label;
    EXPECT_NEAR(energy.bond, reference.bond, 1e-10 * std::abs(reference.bond)) << label;
    EXPECT_NEAR(energy.torsion, reference.torsion, 1e-10 * std::abs(reference.torsion)) << label;
    EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()))
        << label;
    EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-4 * rmsNorm(referenceGradient)) << label;
    EXPECT_LT(maxAbsDifference(engine.molecularPressureTensor(), referencePressure), 1e-4 * maxAbs(referencePressure))
        << label;
  }
}

TEST_P(SpatialDecompositionDevice, engine_device_molecular_terms_with_torsions)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  RandomNumber random(11);
  System system = makeChainSystem(false, random);

  // the same device pair kernel and mesh, with the molecular terms (bonds, bends, torsions, exclusions, virial
  // correction) on the host (double) and on the device (single precision, positions relative to the molecule)
  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  EXPECT_TRUE(host.usesDeviceMesh());
  EXPECT_FALSE(host.usesDeviceBonded());
  const RunningEnergy hostEnergy = host.computeGradients(system, true);
  const std::vector<double3> hostGradient = gradientsOf(system);
  EXPECT_NE(hostEnergy.torsion, 0.0);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  EXPECT_TRUE(device.usesDeviceMesh());
  EXPECT_TRUE(device.usesDeviceBonded());
  const RunningEnergy energy = device.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);

  EXPECT_NEAR(energy.bond, hostEnergy.bond, 1e-5 * std::abs(hostEnergy.bond));
  EXPECT_NEAR(energy.bend, hostEnergy.bend, 1e-5 * std::abs(hostEnergy.bend));
  EXPECT_NEAR(energy.torsion, hostEnergy.torsion, 1e-5 * std::abs(hostEnergy.torsion));
  EXPECT_NE(hostEnergy.intraVDW, 0.0);
  EXPECT_NE(hostEnergy.intraCoul, 0.0);
  EXPECT_NEAR(energy.intraVDW, hostEnergy.intraVDW, 1e-5 * std::abs(hostEnergy.intraVDW));
  EXPECT_NEAR(energy.intraCoul, hostEnergy.intraCoul, 1e-5 * std::abs(hostEnergy.intraCoul));
  EXPECT_NEAR(energy.ewald_self + energy.ewald_exclusion, hostEnergy.ewald_self + hostEnergy.ewald_exclusion,
              1e-5 * std::abs(hostEnergy.ewald_self + hostEnergy.ewald_exclusion));
  EXPECT_NEAR(energy.potentialEnergy(), hostEnergy.potentialEnergy(), 1e-5 * std::abs(hostEnergy.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, hostGradient), 1e-5 * rmsNorm(hostGradient));
  EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), host.molecularPressureTensor()),
            1e-5 * maxAbs(host.molecularPressureTensor()));

  // deterministic: the same configuration evaluates to the same bits (fixed-point spreading, fixed reductions)
  const double3x3 pressure = device.molecularPressureTensor();
  const RunningEnergy again = device.computeGradients(system, true);
  EXPECT_EQ(again.potentialEnergy(), energy.potentialEnergy());
  EXPECT_EQ(again.ewald_fourier, energy.ewald_fourier);
  EXPECT_EQ(again.torsion, energy.torsion);
  EXPECT_EQ(rmsDifference(gradientsOf(system), gradient), 0.0);
  EXPECT_EQ(maxAbsDifference(device.molecularPressureTensor(), pressure), 0.0);
}

TEST_P(SpatialDecompositionDevice, engine_device_molecular_terms_long_chain)
{
  // six beads: the device exclusion lists (1-2, 1-3), a scaled 1-4 pair, excluded 1-4 pairs (scaling 0) and
  // full-strength 1-5 / 1-6 Lennard-Jones and Coulomb pairs
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  RandomNumber random(13);
  System system = makeChainSystem(false, random, 6);
  ASSERT_EQ(system.components[0].intraMolecularPotentials.exclusions.pairs.size(), 9uz);
  ASSERT_EQ(system.components[0].intraMolecularPotentials.numberOfVanDerWaalsPairs(), 6uz);

  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  EXPECT_FALSE(host.usesDeviceBonded());
  const RunningEnergy hostEnergy = host.computeGradients(system, true);
  const std::vector<double3> hostGradient = gradientsOf(system);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  EXPECT_TRUE(device.usesDeviceBonded());
  const RunningEnergy energy = device.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);

  // the scaled 1-4 pairs in the intramolecular slots, the full-strength same-molecule pairs with the
  // molecule-molecule pairs (device lists and host lists alike)
  EXPECT_NE(hostEnergy.intraVDW, 0.0);
  EXPECT_NE(hostEnergy.intraCoul, 0.0);
  EXPECT_NEAR(energy.intraVDW, hostEnergy.intraVDW, 1e-5 * std::abs(hostEnergy.intraVDW));
  EXPECT_NEAR(energy.intraCoul, hostEnergy.intraCoul, 1e-5 * std::abs(hostEnergy.intraCoul));
  EXPECT_NEAR(energy.moleculeMoleculeVDW, hostEnergy.moleculeMoleculeVDW,
              1e-5 * std::abs(hostEnergy.moleculeMoleculeVDW));
  EXPECT_NEAR(energy.moleculeMoleculeCharge, hostEnergy.moleculeMoleculeCharge,
              1e-5 * std::abs(hostEnergy.moleculeMoleculeCharge));
  EXPECT_NEAR(energy.ewald_self + energy.ewald_exclusion, hostEnergy.ewald_self + hostEnergy.ewald_exclusion,
              1e-5 * std::abs(hostEnergy.ewald_self + hostEnergy.ewald_exclusion));
  EXPECT_NEAR(energy.potentialEnergy(), hostEnergy.potentialEnergy(), 1e-5 * std::abs(hostEnergy.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, hostGradient), 1e-5 * rmsNorm(hostGradient));
  EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), host.molecularPressureTensor()),
            1e-5 * maxAbs(host.molecularPressureTensor()));

  // and against the exact all-pairs Ewald evaluation
  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-4 * rmsNorm(referenceGradient));
}

TEST_P(SpatialDecompositionDevice, engine_device_molecular_terms_large_molecules)
{
  // two 343-atom chains: the device bonded path has no molecule-size limit (slots find their molecule through
  // the slot's system atom), and the virial correction uses the per-molecule centers of the center kernel
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  RandomNumber random(23);
  System system = makeLatticeChainSystem(random, 7, 2);
  ASSERT_EQ(system.components[0].atoms.size(), 343uz);
  ASSERT_EQ(system.components[0].intraMolecularPotentials.exclusions.scaledPairs.size(), 340uz);

  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  EXPECT_FALSE(host.usesDeviceBonded());
  const RunningEnergy hostEnergy = host.computeGradients(system, true);
  const std::vector<double3> hostGradient = gradientsOf(system);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  EXPECT_TRUE(device.usesDeviceBonded());
  const RunningEnergy energy = device.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);

  EXPECT_NE(hostEnergy.torsion, 0.0);
  EXPECT_NE(hostEnergy.intraVDW, 0.0);
  EXPECT_NE(hostEnergy.intraCoul, 0.0);
  EXPECT_NEAR(energy.bond, hostEnergy.bond, 1e-5 * std::abs(hostEnergy.bond));
  EXPECT_NEAR(energy.bend, hostEnergy.bend, 1e-5 * std::abs(hostEnergy.bend));
  EXPECT_NEAR(energy.torsion, hostEnergy.torsion, 1e-5 * std::abs(hostEnergy.torsion));
  EXPECT_NEAR(energy.intraVDW, hostEnergy.intraVDW, 1e-5 * std::abs(hostEnergy.intraVDW));
  EXPECT_NEAR(energy.intraCoul, hostEnergy.intraCoul, 1e-5 * std::abs(hostEnergy.intraCoul));
  EXPECT_NEAR(energy.ewald_self + energy.ewald_exclusion, hostEnergy.ewald_self + hostEnergy.ewald_exclusion,
              1e-5 * std::abs(hostEnergy.ewald_self + hostEnergy.ewald_exclusion));
  EXPECT_NEAR(energy.potentialEnergy(), hostEnergy.potentialEnergy(), 1e-5 * std::abs(hostEnergy.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, hostGradient), 1e-5 * rmsNorm(hostGradient));
  EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), host.molecularPressureTensor()),
            1e-5 * maxAbs(host.molecularPressureTensor()));

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-4 * rmsNorm(referenceGradient));
}

TEST(spatial_decomposition, engine_splits_large_molecules_over_threads)
{
  // two 343-atom chains (1700 atoms + terms each): the host bonded work slices them into atom ranges and term
  // ranges per kind, so that the energies, gradients and the molecular pressure come out the same for any
  // number of threads, and the same as the whole-molecule reference
  RandomNumber random(29);
  System system = makeLatticeChainSystem(random, 7, 2);
  ASSERT_EQ(system.components[0].atoms.size(), 343uz);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;
  EXPECT_NE(reference.torsion, 0.0);
  EXPECT_NE(reference.intraVDW, 0.0);
  EXPECT_NE(reference.intraCoul, 0.0);
  EXPECT_NE(reference.ewald_exclusion, 0.0);

  std::optional<RunningEnergy> firstEnergy;
  std::vector<double3> firstGradient;
  double3x3 firstPressure{};
  for (std::size_t threads : {1uz, 3uz, 8uz})
  {
    SpatialDecompositionForceEngine engine(settingsFor(threads, 1.5, 0.5));
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    const double3x3 pressure = engine.molecularPressureTensor();
    const std::string label = std::format("threads {}", threads);

    EXPECT_NEAR(energy.bond, reference.bond, 1e-10 * std::abs(reference.bond)) << label;
    EXPECT_NEAR(energy.bend, reference.bend, 1e-10 * std::abs(reference.bend)) << label;
    EXPECT_NEAR(energy.torsion, reference.torsion, 1e-10 * std::abs(reference.torsion)) << label;
    // the unscaled same-molecule pairs are booked with the molecule-molecule pairs by the engine
    EXPECT_NEAR(energy.intraVDW + energy.moleculeMoleculeVDW, reference.intraVDW + reference.moleculeMoleculeVDW,
                1e-8 * std::abs(reference.intraVDW + reference.moleculeMoleculeVDW))
        << label;
    EXPECT_NEAR(energy.intraCoul + energy.moleculeMoleculeCharge,
                reference.intraCoul + reference.moleculeMoleculeCharge,
                1e-8 * std::abs(reference.intraCoul + reference.moleculeMoleculeCharge))
        << label;
    EXPECT_NE(energy.intraVDW, 0.0) << label;
    EXPECT_NEAR(energy.ewald_self + energy.ewald_exclusion, reference.ewald_self + reference.ewald_exclusion,
                1e-10 * std::abs(reference.ewald_self + reference.ewald_exclusion))
        << label;
    EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 2e-5 * std::abs(reference.potentialEnergy()))
        << label;
    EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-4 * rmsNorm(referenceGradient)) << label;
    EXPECT_LT(maxAbsDifference(pressure, referencePressure), 1e-4 * maxAbs(referencePressure)) << label;

    // the same answer for every thread count (up to the summation order of the pair domains)
    if (!firstEnergy)
    {
      firstEnergy = energy;
      firstGradient = gradient;
      firstPressure = pressure;
    }
    else
    {
      EXPECT_NEAR(energy.potentialEnergy(), firstEnergy->potentialEnergy(),
                  1e-9 * std::abs(firstEnergy->potentialEnergy()))
          << label;
      EXPECT_LT(rmsDifference(gradient, firstGradient), 1e-9 * rmsNorm(firstGradient)) << label;
      EXPECT_LT(maxAbsDifference(pressure, firstPressure), 1e-9 * maxAbs(firstPressure)) << label;
    }
  }
}

TEST_P(SpatialDecompositionDevice, engine_device_molecular_terms_fall_back_to_the_host)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  RandomNumber random(11);
  System system = makeChainSystem(true, random);

  SpatialDecompositionForceEngine reference(settingsFor(4, 1.5, 0.5));
  reference.initialize(system);
  const RunningEnergy referenceEnergy = reference.computeGradients(system, true);
  const std::vector<double3> referenceGradient = gradientsOf(system);

  // the bond-bond cross term keeps the molecular terms on the host; the pairs and the mesh stay on the device
  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  EXPECT_TRUE(device.usesDeviceMesh());
  EXPECT_FALSE(device.usesDeviceBonded());
  const RunningEnergy energy = device.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);

  EXPECT_EQ(energy.bond, referenceEnergy.bond);
  EXPECT_EQ(energy.bend, referenceEnergy.bend);
  EXPECT_EQ(energy.torsion, referenceEnergy.torsion);
  EXPECT_EQ(energy.intraVDW, referenceEnergy.intraVDW);
  EXPECT_EQ(energy.bondBond, referenceEnergy.bondBond);
  EXPECT_NE(energy.bondBond, 0.0);
  EXPECT_NEAR(energy.potentialEnergy(), referenceEnergy.potentialEnergy(),
              1e-5 * std::abs(referenceEnergy.potentialEnergy()));
  EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-5 * rmsNorm(referenceGradient));
  EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), reference.molecularPressureTensor()),
            1e-5 * maxAbs(reference.molecularPressureTensor()));
}

TEST(spatial_decomposition, engine_matches_exact_damped_shifted_force)
{
  ForceField forceField = makeWaterForceField(ForceField::ChargeMethod::DampedShiftedForce);
  forceField.EwaldAlpha = 0.2;
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(5);
  randomizeConfiguration(system, random);

  const RunningEnergy reference = referenceGradients(system);
  const std::vector<double3> referenceGradient = gradientsOf(system);
  const double3x3 referencePressure = system.computeMolecularPressure().second;

  for (std::size_t threads : {1uz, 4uz})
  {
    SpatialDecompositionForceEngine engine(settingsFor(threads));
    engine.initialize(system);
    EXPECT_FALSE(engine.usesFastKernel());  // the shifted real-space electrostatics use the generic kernels
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    EXPECT_NEAR(energy.moleculeMoleculeCharge, reference.moleculeMoleculeCharge,
                1e-9 * std::abs(reference.moleculeMoleculeCharge));
    EXPECT_NEAR(energy.ewald_self, reference.ewald_self, 1e-10 * std::abs(reference.ewald_self));
    EXPECT_NEAR(energy.ewald_exclusion, reference.ewald_exclusion, 1e-9 * std::abs(reference.ewald_exclusion));
    EXPECT_NEAR(energy.ewald_fourier, 0.0, 1e-12);
    EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), 1e-9 * std::abs(reference.potentialEnergy()));
    EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-9 * rmsNorm(referenceGradient));
    EXPECT_LT(maxAbsDifference(engine.molecularPressureTensor(), referencePressure), 1e-8 * maxAbs(referencePressure));
  }
}

TEST(spatial_decomposition, engine_rejects_unsupported_systems)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {10}, 5);
  std::string reason;
  EXPECT_TRUE(SpatialDecompositionForceEngine::supports(system, reason)) << reason;

  system.hasExternalField = true;
  EXPECT_FALSE(SpatialDecompositionForceEngine::supports(system, reason));
  EXPECT_EQ(reason, "an external field");
  system.hasExternalField = false;

  system.forceField.omitInterInteractions = true;
  EXPECT_FALSE(SpatialDecompositionForceEngine::supports(system, reason));
  system.forceField.omitInterInteractions = false;

  // too many threads for the box (30 A, list cutoff 10 A: at most 12 x 12 x 12 cells of a quarter cutoff)
  SpatialDecompositionForceEngine engine(settingsFor(2048, 1.0));
  EXPECT_THROW(engine.initialize(system), std::runtime_error);

  // skin too large for the minimum-image lists
  SpatialDecompositionForceEngine wide(settingsFor(1, 7.0));
  EXPECT_THROW(wide.initialize(system), std::runtime_error);
}

namespace
{
/// The host velocity-Verlet step of the MD driver with the engine forces (Nosé–Hoover at both ends when the
/// system has a thermostat).
RunningEnergy hostVelocityVerlet(System& system, SpatialDecompositionForceEngine& engine)
{
  const auto thermostat = [&]
  {
    if (!system.thermostat.has_value()) return;
    const double translational = Integrators::computeTranslationalKineticEnergy(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
        system.spanOfGroupData(), system.spanOfFrameworkGroupData());
    const double rotational = Integrators::computeRotationalKineticEnergy(
        system.moleculeData, system.components, system.spanOfGroupData(), system.framework,
        system.spanOfFrameworkGroupData());
    const std::pair<double, double> scaling = system.thermostat->NoseHooverNVT(translational, rotational);
    Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                 system.components, scaling, system.framework, system.spanOfFrameworkDynamics(),
                                 system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  };
  thermostat();
  Integrators::updateVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                system.components, system.timeStep, system.framework, system.spanOfFrameworkAtoms(),
                                system.spanOfFrameworkDynamics(), &system.forceField, system.spanOfGroupData(),
                                system.spanOfFrameworkGroupData());
  Integrators::updatePositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                               system.components, system.timeStep, system.framework, system.spanOfFrameworkAtoms(),
                               system.spanOfFrameworkDynamics(), system.spanOfGroupData(),
                               system.spanOfFrameworkGroupData());
  Integrators::noSquishFreeRotorOrderTwo(system.moleculeData, system.components, system.timeStep,
                                         system.spanOfGroupData(), system.framework, system.spanOfFrameworkGroupData());
  Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components,
                                        system.spanOfGroupData(), system.framework, system.spanOfFrameworkAtoms(),
                                        system.spanOfFrameworkGroupData());
  RunningEnergy energies = engine.computeGradients(system, true);
  Integrators::updateCenterOfMassAndQuaternionGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.spanOfGroupData(), system.framework, system.spanOfFrameworkDynamics(), system.spanOfFrameworkGroupData());
  Integrators::updateVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                                system.components, system.timeStep, system.framework, system.spanOfFrameworkAtoms(),
                                system.spanOfFrameworkDynamics(), &system.forceField, system.spanOfGroupData(),
                                system.spanOfFrameworkGroupData());
  thermostat();
  energies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  energies.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                  system.framework, system.spanOfFrameworkGroupData());
  if (system.thermostat.has_value()) energies.NoseHooverEnergy = system.thermostat->getEnergy();
  return energies;
}

/// Forces and molecular gradients of the start of an MD stage (the driver's recomputeGradients).
void startState(System& system, SpatialDecompositionForceEngine& engine)
{
  engine.computeGradients(system, true);
  Integrators::updateCenterOfMassAndQuaternionGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.spanOfGroupData(), system.framework, system.spanOfFrameworkDynamics(), system.spanOfFrameworkGroupData());
}

void giveVelocities(System& system, std::size_t seed, bool thermostat)
{
  RandomNumber random(seed);
  Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                    system.spanOfMoleculeDynamics(), system.components, system.temperature);
  Integrators::removeCenterOfMassVelocityDrift(system.moleculeData, system.spanOfMoleculeAtoms(),
                                               system.spanOfMoleculeDynamics(), system.components);
  if (thermostat)
  {
    system.setThermostat(Thermostat(3, 1, 0.15));
    system.thermostat->initialize(random);
  }
}

struct ResidentTolerances
{
  double position;  ///< absolute, Angstrom, max over atoms
  double velocity;  ///< relative to the rms velocity, max over atoms (flexible) / molecules (rigid)
  double energy;    ///< relative, kinetic and potential energies of every step
  double drift;     ///< relative, conserved-energy drift of the resident run (NVE)
};

/// Runs `steps` steps with the host integrator (device forces, double integration) and with the resident
/// integrator (device forces, df64 integration on the device) from the same start and compares the trajectories,
/// the reported energies and, without thermostat, the conserved-energy drift. The resident state is downloaded
/// to the host every `downloadEvery` steps (the download must not perturb the trajectory).
void expectResidentMatchesHost(const std::function<System(std::size_t)>& makeSystem, PairDevice pairDevice,
                               bool thermostat, std::size_t steps, double timeStep, double skin,
                               const ResidentTolerances& tolerance, std::size_t downloadEvery = 10)
{
  System hostSystem = makeSystem(7);
  System deviceSystem = makeSystem(7);
  hostSystem.timeStep = timeStep;
  deviceSystem.timeStep = timeStep;
  giveVelocities(hostSystem, 21, thermostat);
  giveVelocities(deviceSystem, 21, thermostat);

  SpatialDecompositionSettings hostSettings = settingsFor(4, skin, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.resident = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(hostSystem);
  EXPECT_FALSE(host.usesResident());
  startState(hostSystem, host);

  SpatialDecompositionSettings deviceSettings = settingsFor(4, skin, 0.5);
  deviceSettings.pairDevice = pairDevice;
  deviceSettings.resident = true;
  SpatialDecompositionForceEngine device(deviceSettings);
  device.initialize(deviceSystem);
  ASSERT_TRUE(device.usesResident()) << device.writeStatus();
  startState(deviceSystem, device);

  const bool rigid = hostSystem.components[0].rigid;
  double rmsVelocity = 0.0;
  {
    std::size_t count = 0;
    if (rigid)
    {
      for (const Molecule& molecule : hostSystem.moleculeData)
      {
        rmsVelocity += double3::dot(molecule.velocity, molecule.velocity);
        ++count;
      }
    }
    else
    {
      for (const AtomDynamics& dynamics : hostSystem.spanOfMoleculeDynamics())
      {
        rmsVelocity += double3::dot(dynamics.velocity, dynamics.velocity);
        ++count;
      }
    }
    rmsVelocity = std::sqrt(rmsVelocity / static_cast<double>(count));
  }

  double hostReference = 0.0;
  double deviceReference = 0.0;
  double hostDrift = 0.0;
  double deviceDrift = 0.0;
  double maxPositionError = 0.0;
  double maxVelocityError = 0.0;
  for (std::size_t step = 1; step <= steps; ++step)
  {
    const RunningEnergy hostEnergy = hostVelocityVerlet(hostSystem, host);
    const RunningEnergy deviceEnergy = device.residentVelocityVerlet(deviceSystem);

    EXPECT_NEAR(deviceEnergy.potentialEnergy(), hostEnergy.potentialEnergy(),
                tolerance.energy * std::abs(hostEnergy.potentialEnergy()))
        << "step " << step;
    EXPECT_NEAR(deviceEnergy.translationalKineticEnergy, hostEnergy.translationalKineticEnergy,
                tolerance.energy * std::abs(hostEnergy.translationalKineticEnergy))
        << "step " << step;
    if (rigid)
    {
      EXPECT_NEAR(deviceEnergy.rotationalKineticEnergy, hostEnergy.rotationalKineticEnergy,
                  tolerance.energy * std::abs(hostEnergy.rotationalKineticEnergy))
          << "step " << step;
    }
    else
    {
      EXPECT_EQ(deviceEnergy.rotationalKineticEnergy, 0.0);
    }
    if (thermostat)
    {
      EXPECT_NEAR(deviceEnergy.NoseHooverEnergy, hostEnergy.NoseHooverEnergy,
                  tolerance.energy * std::max(1.0, std::abs(hostEnergy.NoseHooverEnergy)))
          << "step " << step;
    }

    if (step == 1)
    {
      hostReference = hostEnergy.conservedEnergy();
      deviceReference = deviceEnergy.conservedEnergy();
    }
    hostDrift = std::max(hostDrift, std::abs((hostEnergy.conservedEnergy() - hostReference) / hostReference));
    deviceDrift = std::max(deviceDrift, std::abs((deviceEnergy.conservedEnergy() - deviceReference) / deviceReference));

    if (step % downloadEvery == 0 || step == steps)
    {
      device.downloadResidentState(deviceSystem);
      std::span<const Atom> hostAtoms = hostSystem.spanOfMoleculeAtoms();
      std::span<const Atom> deviceAtoms = deviceSystem.spanOfMoleculeAtoms();
      for (std::size_t i = 0; i < hostAtoms.size(); ++i)
      {
        const double3 difference = deviceAtoms[i].position - hostAtoms[i].position;
        maxPositionError = std::max(maxPositionError, std::sqrt(double3::dot(difference, difference)));
      }
      if (rigid)
      {
        for (std::size_t m = 0; m < hostSystem.moleculeData.size(); ++m)
        {
          const double3 difference =
              deviceSystem.moleculeData[m].velocity - hostSystem.moleculeData[m].velocity;
          maxVelocityError = std::max(maxVelocityError, std::sqrt(double3::dot(difference, difference)));
        }
      }
      else
      {
        std::span<const AtomDynamics> hostDynamics = hostSystem.spanOfMoleculeDynamics();
        std::span<const AtomDynamics> deviceDynamics = deviceSystem.spanOfMoleculeDynamics();
        for (std::size_t i = 0; i < hostDynamics.size(); ++i)
        {
          const double3 difference = deviceDynamics[i].velocity - hostDynamics[i].velocity;
          maxVelocityError = std::max(maxVelocityError, std::sqrt(double3::dot(difference, difference)));
        }
      }
    }
  }
  EXPECT_LT(maxPositionError, tolerance.position);
  EXPECT_LT(maxVelocityError, tolerance.velocity * rmsVelocity);
  if (!thermostat)
  {
    // the resident integrator conserves the energy as well as the host integrator does (the drift of the host
    // run is that of the single-precision forces and the time step, not of the integration arithmetic)
    EXPECT_LE(deviceDrift, 1.1 * hostDrift + tolerance.drift) << "host drift " << hostDrift;
  }
  // the resident path integrates on the device: a single host force evaluation at the start, and the rebuild
  // decisions from the device displacement check agree with the host's
  EXPECT_EQ(device.timings().steps, steps + 1);
  EXPECT_GE(device.timings().rebuilds, 2uz);
  EXPECT_NEAR(static_cast<double>(device.timings().rebuilds), static_cast<double>(host.timings().rebuilds), 1.0);
}

System makeWaterSystem(std::size_t seed)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {343}, 5);
  RandomNumber random(seed);
  randomizeConfiguration(system, random);
  return system;
}
}  // namespace

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_rigid_water_nve)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // 100 steps of 1 fs with a 0.5 A skin: several neighbour-list rebuilds from the device displacement check (the
  // trajectories of the two integrators separate exponentially from the rounding of the single-precision
  // positions; the per-step energies agree to the force precision)
  expectResidentMatchesHost(makeWaterSystem, pairDevice, false, 100, 0.001, 0.5,
                            {.position = 1e-3, .velocity = 1e-3, .energy = 1e-4, .drift = 1e-5});
}

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_rigid_water_nvt)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  expectResidentMatchesHost(makeWaterSystem, pairDevice, true, 60, 0.002, 0.3,
                            {.position = 1e-3, .velocity = 1e-3, .energy = 1e-4, .drift = 1e-5});
}

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_flexible_chains_nve)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  const auto makeSystem = [](std::size_t seed)
  {
    RandomNumber random(seed);
    return makeChainSystem(false, random);
  };
  expectResidentMatchesHost(makeSystem, pairDevice, false, 100, 0.001, 0.3,
                            {.position = 1e-4, .velocity = 1e-4, .energy = 1e-4, .drift = 1e-5});
}

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_flexible_chains_nvt)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  const auto makeSystem = [](std::size_t seed)
  {
    RandomNumber random(seed);
    return makeChainSystem(false, random);
  };
  expectResidentMatchesHost(makeSystem, pairDevice, true, 60, 0.001, 0.3,
                            {.position = 1e-4, .velocity = 1e-4, .energy = 1e-4, .drift = 1e-5});
}

TEST_P(SpatialDecompositionDevice, resident_integrator_resumes_after_host_evaluation)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // a host force evaluation in between (status report, restart) invalidates the device copy: the next resident
  // step re-uploads the host state and continues the same trajectory
  System uninterrupted = makeWaterSystem(3);
  System interrupted = makeWaterSystem(3);
  uninterrupted.timeStep = 0.002;
  interrupted.timeStep = 0.002;
  giveVelocities(uninterrupted, 5, true);
  giveVelocities(interrupted, 5, true);

  SpatialDecompositionSettings settings = settingsFor(4, 1.0, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine first(settings);
  first.initialize(uninterrupted);
  startState(uninterrupted, first);
  SpatialDecompositionForceEngine second(settings);
  second.initialize(interrupted);
  startState(interrupted, second);
  ASSERT_TRUE(first.usesResident());
  ASSERT_TRUE(second.usesResident());

  for (std::size_t step = 0; step < 10; ++step)
  {
    first.residentVelocityVerlet(uninterrupted);
    second.residentVelocityVerlet(interrupted);
  }
  second.downloadResidentState(interrupted);
  const RunningEnergy hostEvaluation = second.computeGradients(interrupted, true);
  Integrators::updateCenterOfMassAndQuaternionGradients(
      interrupted.moleculeData, interrupted.spanOfMoleculeAtoms(), interrupted.spanOfMoleculeDynamics(),
      interrupted.components, interrupted.spanOfGroupData(), interrupted.framework,
      interrupted.spanOfFrameworkDynamics(), interrupted.spanOfFrameworkGroupData());
  for (std::size_t step = 0; step < 10; ++step)
  {
    const RunningEnergy a = first.residentVelocityVerlet(uninterrupted);
    const RunningEnergy b = second.residentVelocityVerlet(interrupted);
    EXPECT_NEAR(a.potentialEnergy(), b.potentialEnergy(), 1e-5 * std::abs(a.potentialEnergy())) << "step " << step;
    EXPECT_NEAR(a.translationalKineticEnergy, b.translationalKineticEnergy, 1e-5 * a.translationalKineticEnergy)
        << "step " << step;
    EXPECT_NEAR(a.rotationalKineticEnergy, b.rotationalKineticEnergy, 1e-5 * a.rotationalKineticEnergy)
        << "step " << step;
  }
  // the host evaluation saw the state of step 10
  EXPECT_NE(hostEvaluation.potentialEnergy(), 0.0);
  first.downloadResidentState(uninterrupted);
  second.downloadResidentState(interrupted);
  std::span<const Atom> a = uninterrupted.spanOfMoleculeAtoms();
  std::span<const Atom> b = interrupted.spanOfMoleculeAtoms();
  double maxDifference = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    const double3 difference = a[i].position - b[i].position;
    maxDifference = std::max(maxDifference, std::sqrt(double3::dot(difference, difference)));
  }
  EXPECT_LT(maxDifference, 1e-5);
}
