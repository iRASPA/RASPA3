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
import randomnumbers;
import spatial_decomposition_settings;
import spatial_decomposition_cell_list;
import spatial_decomposition_pppm;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_worker_team;
import spatial_decomposition_force_engine;
import force_engine;
import spatial_decomposition_opencl_pair_kernel;

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

TEST(spatial_decomposition, engine_opencl_pair_kernel_agrees_with_double_rigid_water)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  // the device kernel: single-precision pair geometry with the minimum image of wrapped positions, closed-form
  // erfc, single-precision accumulation of the per-atom forces and per-cluster partial sums
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = PairDevice::OpenCL;
  settings.deviceMesh = false;
  settings.deviceBonded = false;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-6, .charge = 1e-4, .total = 1e-5, .gradient = 1e-5, .pressure = 1e-5});
}

TEST(spatial_decomposition, engine_opencl_pair_kernel_without_pruning_rigid_water)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  // no pruning: the lane lists hold the whole outer list (cutoff + Verlet skin), compacted once per build
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = PairDevice::OpenCL;
  settings.pruneSkin = 0.0;
  settings.deviceMesh = false;
  settings.deviceBonded = false;
  expectKernelAgreesWithScalar(settings,
                               {.vdw = 1e-6, .charge = 1e-4, .total = 1e-5, .gradient = 1e-5, .pressure = 1e-5});
}

TEST(spatial_decomposition, engine_opencl_mesh_and_molecular_terms_agree_with_double_rigid_water)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  // the complete device step: pairs, PPPM (fixed-point spreading, single-precision FFT and influence function,
  // gather interpolation) and the molecular terms (self + exclusion without their cancellation, virial correction).
  // Measured: Fourier energy 1e-6, self + exclusion sum 2e-7, total energy 3.4e-6, gradient rms 2.5e-6, pressure
  // 1e-6 relative.
  SpatialDecompositionSettings settings = settingsFor(1);
  settings.pairDevice = PairDevice::OpenCL;
  expectKernelAgreesWithScalar(settings, {.vdw = 1e-6,
                                          .charge = 1e-4,
                                          .total = 1e-5,
                                          .gradient = 1e-5,
                                          .pressure = 1e-5,
                                          .fourier = 1e-5,
                                          .correction = 1e-5});
}

TEST(spatial_decomposition, engine_opencl_pair_kernel_small_grid_rigid_water)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  // cutoff + skin = half the box: the device grid has 2 cells per axis, every neighbour cell is reached through
  // several stencil offsets and the list build takes the minimum image (the outer list then holds every pair)
  SpatialDecompositionSettings settings = settingsFor(1, 6.0);
  settings.pairDevice = PairDevice::OpenCL;
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
  if (OpenCLPairKernel::available())
  {
    SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
    settings.pairDevice = PairDevice::OpenCL;
    SpatialDecompositionForceEngine device(settings);
    device.initialize(system);
    const RunningEnergy deviceEnergy = device.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    EXPECT_NEAR(deviceEnergy.moleculeMoleculeVDW, energy.moleculeMoleculeVDW,
                1e-6 * std::abs(energy.moleculeMoleculeVDW));
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
/// A charged four-bead chain with harmonic bonds and bends, a TraPPE torsion, intramolecular 1-4 Lennard-Jones and
/// Coulomb pairs and, optionally, a bond-bond cross term (which the device kernels do not cover).
System makeChainSystem(bool withBondBond, RandomNumber& random)
{
  ForceField forceField = ForceField(
      {{"CH3", false, 15.03452, 0.0, 0.0, 6, false}, {"CH2", false, 14.02658, 0.0, 0.0, 6, false}},
      {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(14, 14, 14);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 196;

  ConnectivityTable connectivityTable(4);
  for (std::size_t i = 0; i + 1 < 4; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({1, 2}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({2, 3}, BondType::Harmonic, {96500.0, 1.54})};
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0}),
                      BendPotential({1, 2, 3}, BendType::Harmonic, {62500.0, 114.0})};
  potentials.torsions = {TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32})};
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VanDerWaalsType::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.25, 0.25, 0.5)};
  if (withBondBond)
  {
    potentials.bondBonds = {BondBondPotential({0, 1, 2}, BondBondType::CFF, {5000.0, 1.54, 1.54})};
  }

  Component chain = Component(forceField, "butane", 425.0, 3796000.0, 0.199,
                              {Atom({-1.85, -0.7, -0.15}, 0.25, 1.0, 0, 0, 0, false, false),
                               Atom({-0.31, -0.7, -0.15}, -0.25, 1.0, 0, 1, 0, false, false),
                               Atom({0.32, 0.71, -0.15}, -0.25, 1.0, 0, 1, 0, false, false),
                               Atom({1.86, 0.71, 0.15}, 0.25, 1.0, 0, 0, 0, false, false)},
                              connectivityTable, potentials, 5, 21);
  // 4 x 4 x 4 molecules on a lattice of about 10 Angstrom: no close contacts (the fallback test compares the
  // device pairs with the host pairs, so the pair energies must not be dominated by overlaps)
  System system = System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain}, {}, {64}, 5);
  randomizeConfiguration(system, random);
  return system;
}
}  // namespace

TEST(spatial_decomposition, engine_opencl_molecular_terms_with_torsions)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  RandomNumber random(11);
  System system = makeChainSystem(false, random);

  // the same device pair kernel and mesh, with the molecular terms (bonds, bends, torsions, exclusions, virial
  // correction) on the host (double) and on the device (single precision, positions relative to the molecule)
  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = PairDevice::OpenCL;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  EXPECT_TRUE(host.usesDeviceMesh());
  EXPECT_FALSE(host.usesDeviceBonded());
  const RunningEnergy hostEnergy = host.computeGradients(system, true);
  const std::vector<double3> hostGradient = gradientsOf(system);
  EXPECT_NE(hostEnergy.torsion, 0.0);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = PairDevice::OpenCL;
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

TEST(spatial_decomposition, engine_opencl_molecular_terms_fall_back_to_the_host)
{
  if (!OpenCLPairKernel::available()) GTEST_SKIP() << "no OpenCL device";
  RandomNumber random(11);
  System system = makeChainSystem(true, random);

  SpatialDecompositionForceEngine reference(settingsFor(4, 1.5, 0.5));
  reference.initialize(system);
  const RunningEnergy referenceEnergy = reference.computeGradients(system, true);
  const std::vector<double3> referenceGradient = gradientsOf(system);

  // the bond-bond cross term keeps the molecular terms on the host; the pairs and the mesh stay on the device
  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = PairDevice::OpenCL;
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
