#include <gtest/gtest.h>

import std;

import ewald_structure_factor;

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
import mc_moves_hybridmc;
import mc_moves_move_types;
import polymer_links;
import reactive_group;
import polymer_reactions;

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

void DOUBLE3_EXPECT_NEAR(double3 a, double3 b, double tol);

namespace
{
std::vector<double3> positionsOf(const System& system)
{
  std::vector<double3> result;
  for (const Atom& atom : system.spanOfMoleculeAtoms()) result.push_back(atom.position);
  return result;
}

// A dilute lattice of waters (no overlaps, so that the short trajectories stay well-behaved) with the running
// energies of the system set as the Monte Carlo drivers do before the moves.
System makeHybridMCWaterSystem(ForceField& forceField, Component& water, RandomNumber& random)
{
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {216}, 5);
  randomizeConfiguration(system, random);
  system.runningEnergies = system.computeTotalEnergies();
  system.numberOfHybridMCSteps = 5;
  system.forceEngineSettings = settingsFor(2);
  return system;
}
}  // namespace

TEST(spatial_decomposition, hybrid_mc_engine_path_follows_the_exact_trajectory_rigid_water)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  RandomNumber configuration(11);
  System exact = makeHybridMCWaterSystem(forceField, water, configuration);
  System engine = exact;
  exact.hybridMCUseForceEngine = false;
  engine.hybridMCUseForceEngine = true;
  exact.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5e-4);
  engine.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 5e-4);

  // the same random stream: the same Boltzmann velocities and the same acceptance draw
  RandomNumber randomExact(5);
  RandomNumber randomEngine(5);
  const std::optional<RunningEnergy> acceptedExact = MC_Moves::hybridMCMove(randomExact, exact);
  const std::optional<RunningEnergy> acceptedEngine = MC_Moves::hybridMCMove(randomEngine, engine);
  ASSERT_TRUE(acceptedExact.has_value());
  ASSERT_TRUE(acceptedEngine.has_value());
  ASSERT_TRUE(engine.hybridMCForceEngine != nullptr);
  EXPECT_TRUE(exact.hybridMCForceEngine == nullptr);

  // the engine forces (mesh Ewald, tabulated erfc) steer the same trajectory to within their accuracy
  const std::vector<double3> exactPositions = positionsOf(exact);
  const std::vector<double3> enginePositions = positionsOf(engine);
  EXPECT_LT(rmsDifference(enginePositions, exactPositions), 1e-5);

  // the returned energies are the exact energies of the new configuration (what the caller stores as the running
  // energies), tail correction included, kinetic energies of the trajectory attached
  const RunningEnergy recomputed = engine.computeTotalEnergies();
  EXPECT_NEAR(acceptedEngine->potentialEnergy(), recomputed.potentialEnergy(),
              1e-10 * std::abs(recomputed.potentialEnergy()));
  EXPECT_NEAR(acceptedEngine->moleculeMoleculeVDW, recomputed.moleculeMoleculeVDW,
              1e-10 * std::abs(recomputed.moleculeMoleculeVDW));
  EXPECT_NEAR(acceptedEngine->ewald_fourier, recomputed.ewald_fourier, 1e-10 * std::abs(recomputed.ewald_fourier));
  EXPECT_NEAR(acceptedEngine->tail, recomputed.tail, 1e-10 * std::abs(recomputed.tail));
  EXPECT_GT(acceptedEngine->translationalKineticEnergy, 0.0);
  EXPECT_GT(acceptedEngine->rotationalKineticEnergy, 0.0);
  EXPECT_NEAR(acceptedEngine->potentialEnergy(), acceptedExact->potentialEnergy(),
              1e-6 * std::abs(acceptedExact->potentialEnergy()));
  EXPECT_DOUBLE_EQ(engine.timeStep, 5e-4);

  // the stored Ewald structure factors belong to the new configuration: a subsequent move must see no difference
  const EwaldStructureFactor stored = engine.storedEik;
  engine.computeTotalEnergies();
  ASSERT_EQ(stored.size(), engine.storedEik.size());
  for (std::size_t k = 0; k < stored.size(); ++k)
  {
    EXPECT_NEAR(std::abs(stored.sum[k] - engine.storedEik.sum[k]), 0.0, 1e-9);
  }
}

TEST(spatial_decomposition, hybrid_mc_engine_path_restores_the_state_on_rejection)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  RandomNumber configuration(13);
  System system = makeHybridMCWaterSystem(forceField, water, configuration);
  // a time step far too large: the energy drift makes the rejection (almost) certain
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 0.05);

  const std::vector<double3> initialPositions = positionsOf(system);
  const std::vector<Molecule> initialMolecules = system.moleculeData;
  const EwaldStructureFactor initialStored(system.storedEik);
  const RunningEnergy initialEnergies = system.runningEnergies;

  RandomNumber random(17);
  std::optional<RunningEnergy> accepted;
  std::size_t rejected = 0;
  for (std::size_t attempt = 0; attempt < 5 && !accepted.has_value(); ++attempt)
  {
    accepted = MC_Moves::hybridMCMove(random, system);
    if (accepted.has_value()) break;
    ++rejected;
    const std::vector<double3> positions = positionsOf(system);
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
      DOUBLE3_EXPECT_NEAR(positions[i], initialPositions[i], 0.0);
    }
    for (std::size_t m = 0; m < initialMolecules.size(); ++m)
    {
      DOUBLE3_EXPECT_NEAR(system.moleculeData[m].centerOfMassPosition, initialMolecules[m].centerOfMassPosition, 0.0);
      EXPECT_EQ(system.moleculeData[m].orientation.r, initialMolecules[m].orientation.r);
    }
    ASSERT_EQ(system.storedEik.size(), initialStored.size());
    for (std::size_t k = 0; k < initialStored.size(); ++k)
    {
      EXPECT_EQ(system.storedEik.sum[k], initialStored.sum[k]);
    }
    EXPECT_FALSE(system.thermostat.has_value());
    EXPECT_DOUBLE_EQ(system.timeStep, 0.0005);
    // the restored configuration has the running energies it had
    const RunningEnergy recomputed = system.computeTotalEnergies();
    EXPECT_NEAR(recomputed.potentialEnergy(), initialEnergies.potentialEnergy(),
                1e-10 * std::abs(initialEnergies.potentialEnergy()));
  }
  EXPECT_GT(rejected, 0uz);
}

TEST(spatial_decomposition, hybrid_mc_engine_path_flexible_chains_matches_the_exact_energy)
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

  System system =
      System(forceField, SimulationBox(32.0, 32.0, 32.0), false, 300.0, 1e5, 1.0, {}, {chain}, {}, {125}, 5);
  RandomNumber configuration(3);
  randomizeConfiguration(system, configuration);
  system.runningEnergies = system.computeTotalEnergies();
  system.numberOfHybridMCSteps = 4;
  system.forceEngineSettings = settingsFor(4, 1.5, 0.5);
  system.mc_moves_statistics.setMaxChange(Move::Types::HybridMC, 2e-4);

  RandomNumber random(29);
  std::optional<RunningEnergy> accepted;
  for (std::size_t attempt = 0; attempt < 10 && !accepted.has_value(); ++attempt)
  {
    accepted = MC_Moves::hybridMCMove(random, system);
  }
  ASSERT_TRUE(accepted.has_value());
  ASSERT_TRUE(system.hybridMCForceEngine != nullptr);
  const RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(accepted->potentialEnergy(), recomputed.potentialEnergy(),
              1e-10 * std::abs(recomputed.potentialEnergy()));
  EXPECT_NEAR(accepted->bond, recomputed.bond, 1e-10 * std::abs(recomputed.bond));
  EXPECT_NEAR(accepted->bend, recomputed.bend, 1e-10 * std::abs(recomputed.bend));
  EXPECT_NE(accepted->bond, 0.0);
}

namespace
{
/// Moves molecule 'b' (of a single-component system of five-bead molecules) rigidly so that its atom 'atomB' sits
/// at the bond distance from atom 'atomA' of molecule 'a' along the open direction of 'atomA' (away from its
/// bonded neighbours), with 'b' rotated so that the open direction of 'atomB' points back at 'atomA': a chain-like
/// junction. The spin of 'b' about the new bond is the one that keeps its other beads farthest from every other
/// atom of the system (the rest of 'a' and the other molecules), so that the configuration has no overlaps and can
/// be integrated.
void approachMolecules(System& system, std::size_t a, std::size_t atomA, std::size_t b, std::size_t atomB,
                       double distance = 1.6, std::optional<double3> approachDirection = std::nullopt)
{
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const ConnectivityTable& connectivity = system.components[0].connectivityTable;
  const std::size_t n = system.components[0].atoms.size();
  const auto openDirection = [&](std::size_t m, std::size_t atom)
  {
    double3 sum{};
    for (std::size_t k = 0; k < n; ++k)
    {
      if (!connectivity[atom, k]) continue;
      sum += double3::normalize(atoms[n * m + k].position - atoms[n * m + atom].position);
    }
    return double3::normalize(-sum);
  };
  const double3 direction = approachDirection.value_or(openDirection(a, atomA));
  // Rodrigues: rotate 'b' about its atom 'atomB' so that openDirection(b, atomB) becomes -direction; an atom with
  // a single neighbour (a head) is tilted so that the junction bend is 110 degrees instead of collinear (the
  // dihedrals through a collinear bend are undefined)
  const double3 u = openDirection(b, atomB);
  double3 v = -direction;
  std::size_t neighbours = 0;
  for (std::size_t k = 0; k < n; ++k) neighbours += connectivity[atomB, k] ? 1 : 0;
  if (neighbours == 1)
  {
    const double3 trial = std::abs(direction.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
    const double3 perpendicular = double3::normalize(double3::cross(direction, trial));
    const double angle = 110.0 * std::numbers::pi / 180.0;
    // the bond of atomB points at cos(angle) along -direction, sin(angle) along the perpendicular
    v = -(std::cos(angle) * (-direction) + std::sin(angle) * perpendicular);
  }
  const double3 axis = double3::cross(u, v);
  const double sine = std::sqrt(double3::dot(axis, axis));
  const double cosine = double3::dot(u, v);
  const double3 pivot = atoms[n * b + atomB].position;
  for (std::size_t k = 0; k < n; ++k)
  {
    const double3 r = atoms[n * b + k].position - pivot;
    double3 rotated = r;
    if (sine > 1e-12)
    {
      const double3 unit = axis * (1.0 / sine);
      rotated = r * cosine + double3::cross(unit, r) * sine + unit * (double3::dot(unit, r) * (1.0 - cosine));
    }
    else if (cosine < 0.0)
    {
      rotated = -r;
    }
    atoms[n * b + k].position = pivot + rotated;
  }
  const double3 shift = atoms[n * a + atomA].position + distance * direction - pivot;
  for (std::size_t k = 0; k < n; ++k) atoms[n * b + k].position += shift;
  system.moleculeData[b].centerOfMassPosition += shift;

  // the spin about the bond axis: the angle with the largest minimum distance of b's beads to the other atoms
  const double3 axisPoint = atoms[n * b + atomB].position;
  const auto spun = [&](std::size_t k, double angle)
  {
    const double3 r = atoms[n * b + k].position - axisPoint;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return axisPoint + r * c + double3::cross(direction, r) * s + direction * (double3::dot(direction, r) * (1.0 - c));
  };
  double bestAngle = 0.0;
  double bestDistance = -1.0;
  for (std::size_t step = 0; step < 72; ++step)
  {
    const double angle = static_cast<double>(step) * (2.0 * std::numbers::pi / 72.0);
    double minimum = std::numeric_limits<double>::max();
    for (std::size_t k = 0; k < n; ++k)
    {
      if (k == atomB) continue;
      const double3 position = spun(k, angle);
      for (std::size_t j = 0; j < atoms.size(); ++j)
      {
        // not b itself, nor the bonded partner (the 1-3 distances across the new bond are set by the geometry)
        if (j / n == b || j == n * a + atomA) continue;
        const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(position - atoms[j].position);
        minimum = std::min(minimum, double3::dot(dr, dr));
      }
    }
    if (minimum > bestDistance)
    {
      bestDistance = minimum;
      bestAngle = angle;
    }
  }
  for (std::size_t k = 0; k < n; ++k) atoms[n * b + k].position = spun(k, bestAngle);
}

/// The free-radical chain growth of a 'vinyl' group (roles head, tail, side): initiation (monomer -> radical, the
/// head and tail retyped to the backbone types), propagation (radical tail - monomer head, the monomer becomes
/// the radical, its head and tail retyped; a stereocentre on the radical's tail whose next centre is the monomer)
/// and termination by combination (tail - tail, a stereocentre on both tails).
std::vector<ReactionRule> vinylChainGrowth()
{
  ReactionRule initiation{};
  initiation.name = "initiation";
  initiation.reactants = {{.label = "A", .type = "vinyl", .state = "monomer"}};
  initiation.products = {{.state = "radical", .retype = {{"head", "CH2_bb"}, {"tail", "CH_bb"}}}};

  ReactionRule propagation{};
  propagation.name = "propagation";
  propagation.reactants = {{.label = "R", .type = "vinyl", .state = "radical"},
                           {.label = "M", .type = "vinyl", .state = "monomer"}};
  propagation.products = {{.state = "dead"}, {.state = "radical", .retype = {{"head", "CH2_bb"}, {"tail", "CH_bb"}}}};
  propagation.capture = std::array<RoleRef, 2>{RoleRef{"R", "tail"}, RoleRef{"M", "head"}};
  propagation.formBonds = {{RoleRef{"R", "tail"}, RoleRef{"M", "head"}}};
  propagation.stereo = {{.centre = {"R", "tail"},
                         .substituents = {RoleRef{"R", "head"}, RoleRef{"M", "head"}, RoleRef{"R", "side"}},
                         .next = "M"}};
  propagation.probability = 0.5;

  ReactionRule termination{};
  termination.name = "termination";
  termination.reactants = {{.label = "A", .type = "vinyl", .state = "radical"},
                           {.label = "B", .type = "vinyl", .state = "radical"}};
  termination.products = {{.state = "dead"}, {.state = "dead"}};
  termination.capture = std::array<RoleRef, 2>{RoleRef{"A", "tail"}, RoleRef{"B", "tail"}};
  termination.formBonds = {{RoleRef{"A", "tail"}, RoleRef{"B", "tail"}}};
  termination.stereo = {
      {.centre = {"A", "tail"}, .substituents = {RoleRef{"A", "head"}, RoleRef{"B", "tail"}, RoleRef{"A", "side"}}},
      {.centre = {"B", "tail"}, .substituents = {RoleRef{"B", "head"}, RoleRef{"A", "tail"}, RoleRef{"B", "side"}}}};
  return {initiation, propagation, termination};
}

/// A branched five-bead monomer (a united-atom vinyl: head 0 = CH3, tail 1 = CH2 with the side group 4, and the
/// tail's chain 2 - 3) whose head and tail take the backbone types CH2_bb / CH_bb when the group reacts. Every
/// polymer-bond template exists (all type combinations, with parameters that depend on the number of backbone
/// types in the term, so that the retyped-term overrides are non-zero), and the component's 1-4 pairs are scaled
/// by one half with the first pair containing the head (a retyped atom of a scaled pair).
System makePolymerSystem(std::size_t seed, std::size_t dimers, std::size_t terminations, std::uint64_t currentStep,
                         std::size_t relaxationSteps, double chiralForceConstant, double bondDistance = 1.6)
{
  ForceField forceField = ForceField(
      {{"CH3", false, 15.03452, 0.25, 0.0, 6, false},
       {"CH2", false, 14.02658, -0.25, 0.0, 6, false},
       {"CH2_bb", false, 14.02658, 0.35, 0.0, 6, false},
       {"CH_bb", false, 13.01864, -0.35, 0.0, 6, false}},
      {{98.0, 3.75}, {46.0, 3.95}, {60.0, 3.85}, {20.0, 4.30}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0,
      9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(20, 20, 20);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 400;

  ConnectivityTable connectivityTable(5);
  for (const auto [a, b] : {std::pair{0uz, 1uz}, std::pair{1uz, 2uz}, std::pair{2uz, 3uz}, std::pair{1uz, 4uz}})
  {
    connectivityTable[a, b] = true;
    connectivityTable[b, a] = true;
  }
  Potentials::IntraMolecularPotentials potentials{};
  for (const auto [a, b] : {std::pair{0uz, 1uz}, std::pair{1uz, 2uz}, std::pair{2uz, 3uz}, std::pair{1uz, 4uz}})
  {
    potentials.bonds.push_back(BondPotential({a, b}, BondType::Harmonic, {96500.0, 1.54}));
  }
  for (const auto [a, b, c] : {std::tuple{0uz, 1uz, 2uz}, std::tuple{0uz, 1uz, 4uz}, std::tuple{2uz, 1uz, 4uz},
                               std::tuple{1uz, 2uz, 3uz}})
  {
    potentials.bends.push_back(BendPotential({a, b, c}, BendType::Harmonic, {62500.0, 112.0}));
  }
  potentials.torsions = {
      TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32}),
      TorsionPotential({4, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32})};
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VDWParameters::Type::LennardJones, {98.0, 3.75}, 0.5),
                            VanDerWaalsPotential({4, 3}, VDWParameters::Type::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.25, 0.25, 0.5),
                         CoulombPotential({4, 3}, CoulombType::Coulomb, 0.0, 0.25, 0.5)};

  const std::vector<Atom> beads{Atom({-1.26, -0.89, 0.0}, 0.25, 1.0, 0, 0, 0, false, false),
                                Atom({0.0, 0.0, 0.0}, -0.25, 1.0, 0, 1, 0, false, false),
                                Atom({1.26, -0.89, 0.0}, -0.25, 1.0, 0, 1, 0, false, false),
                                Atom({2.52, 0.0, 0.0}, 0.25, 1.0, 0, 0, 0, false, false),
                                Atom({0.0, 0.77, 1.33}, 0.0, 1.0, 0, 0, 0, false, false)};
  Component monomer =
      Component(forceField, "vinyl", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  monomer.intra14VanDerWaalsScaling = 0.5;
  monomer.intra14ChargeChargeScaling = 0.5;
  monomer.reactiveGroups.push_back(
      ReactiveGroup{.type = "vinyl", .atoms = {{"head", 0}, {"tail", 1}, {"side", 4}}, .initialState = "monomer"});

  System system =
      System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {monomer}, {}, {64}, 5);
  RandomNumber random(seed);
  randomizeConfiguration(system, random);

  // the templates: every combination of the four types (canonical orientation only), backbone-dependent parameters
  const std::array<std::string, 4> names{"CH3", "CH2", "CH2_bb", "CH_bb"};
  const auto backbone = [](const std::string& name) { return name.ends_with("_bb") ? 1.0 : 0.0; };
  PolymerLinkTable& table = system.polymerLinks;
  for (std::size_t a = 0; a < 4; ++a)
  {
    for (std::size_t b = a; b < 4; ++b)
    {
      const double n = backbone(names[a]) + backbone(names[b]);
      table.templates.bonds.push_back(
          {.types = {names[a], names[b]}, .type = BondType::Harmonic, .parameters = {96500.0 + 5000.0 * n, 1.54 - 0.01 * n}});
    }
  }
  for (std::size_t a = 0; a < 4; ++a)
  {
    for (std::size_t b = 0; b < 4; ++b)
    {
      for (std::size_t c = a; c < 4; ++c)
      {
        const double n = backbone(names[a]) + backbone(names[b]) + backbone(names[c]);
        table.templates.bends.push_back({.types = {names[a], names[b], names[c]},
                                         .type = BendType::Harmonic,
                                         .parameters = {62500.0 + 2000.0 * n, 112.0 - 1.5 * n}});
      }
    }
  }
  for (std::size_t a = 0; a < 4; ++a)
  {
    for (std::size_t b = 0; b < 4; ++b)
    {
      for (std::size_t c = 0; c < 4; ++c)
      {
        for (std::size_t d = 0; d < 4; ++d)
        {
          const std::array<std::size_t, 4> key{a, b, c, d};
          const std::array<std::size_t, 4> reversedKey{d, c, b, a};
          if (reversedKey < key) continue;
          const double n = backbone(names[a]) + backbone(names[b]) + backbone(names[c]) + backbone(names[d]);
          table.templates.torsions.push_back({.types = {names[a], names[b], names[c], names[d]},
                                              .type = TorsionType::TraPPE,
                                              .parameters = {0.0, 355.03 + 10.0 * n, -68.19, 791.32 + 20.0 * n}});
        }
      }
    }
  }
  table.templates.resolve(forceField);
  table.scheme.rules = vinylChainGrowth();
  table.scheme.resolve(system.components, forceField);
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.relaxationSteps = relaxationSteps;
  table.currentStep = currentStep;
  table.chiralForceConstant = chiralForceConstant * Units::KelvinToEnergy;
  table.chiralReferenceVolume = 0.6;

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  // one group per molecule: the group index is the molecule index
  const auto apply = [&](std::size_t rule, std::initializer_list<std::size_t> groups,
                         std::initializer_list<std::int8_t> chiralities, std::uint64_t step)
  {
    const std::vector<std::int8_t> next(chiralities.size(), 0);
    table.apply(rule, std::span<const std::size_t>(groups.begin(), groups.size()),
                std::span<const std::int8_t>(chiralities.begin(), chiralities.size()), next, forceField,
                system.components, atoms, std::span<Molecule>(system.moleculeData), system.simulationBox, step);
  };

  // the first 2 * dimers molecules pair up by propagation (tail of the radical to the head of the monomer), the
  // next 2 * terminations by two initiations and a termination (tail to tail); the first half of the bonds formed
  // at step 0, the second half at step 'currentStep' / 2 (different ramp factors when a ramp is active)
  std::size_t molecule = 0;
  std::size_t formed = 0;
  const std::size_t total = dimers + terminations;
  const auto stepOf = [&](std::size_t index) { return index < total / 2 ? 0 : currentStep / 2; };
  for (std::size_t n = 0; n < dimers; ++n, molecule += 2)
  {
    const std::size_t a = molecule;
    const std::size_t b = molecule + 1;
    approachMolecules(system, a, 1, b, 0, bondDistance);
    apply(0, {a}, {}, stepOf(formed));
    const std::int8_t chirality = (n % 3 == 0) ? 0 : ((n % 3 == 1) ? 1 : -1);
    apply(1, {a, b}, {chirality}, stepOf(formed));
    ++formed;
  }
  for (std::size_t n = 0; n < terminations; ++n, molecule += 2)
  {
    const std::size_t a = molecule;
    const std::size_t b = molecule + 1;
    approachMolecules(system, a, 1, b, 1, bondDistance);
    apply(0, {a}, {}, stepOf(formed));
    apply(0, {b}, {}, stepOf(formed));
    apply(2, {a, b}, {static_cast<std::int8_t>((n % 2 == 0) ? 1 : -1), static_cast<std::int8_t>((n % 2 == 0) ? -1 : 1)},
          stepOf(formed));
    ++formed;
  }
  table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, atoms);
  return system;
}

// ---- ring opening (bond breaking) --------------------------------------------------------------------------------

/// The pseudo-atoms of the epoxy-amine model: a monomer carrying a three-membered ring (C_ep 0, CH2_ep 1, O_ep 2)
/// with a linker (CH2 3) to a primary amine (N 4); the ring opening retypes the attacked carbon and the oxygen
/// (C_op, O_op) and the amine (NH). The product types carry the index >= 5.
const std::array<std::string, 8> ringTypeNames{"C_ep", "CH2_ep", "O_ep", "CH2", "N", "C_op", "O_op", "NH"};

ForceField ringForceField()
{
  ForceField forceField = ForceField(
      {{"C_ep", false, 13.01864, 0.20, 0.0, 6, false},
       {"CH2_ep", false, 14.02658, 0.20, 0.0, 6, false},
       {"O_ep", false, 15.9994, -0.40, 0.0, 8, false},
       {"CH2", false, 14.02658, 0.30, 0.0, 6, false},
       {"N", false, 16.0226, -0.30, 0.0, 7, false},
       {"C_op", false, 13.01864, 0.25, 0.0, 6, false},
       {"O_op", false, 17.0073, -0.45, 0.0, 8, false},
       {"NH", false, 15.0146, -0.30, 0.0, 7, false}},
      {{50.0, 3.90}, {46.0, 3.95}, {60.0, 3.00}, {46.0, 3.95}, {80.0, 3.30}, {55.0, 3.85}, {90.0, 3.05}, {85.0, 3.35}},
      ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(20, 20, 20);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 400;
  return forceField;
}

/// The bonded parameters of the product templates: a formula in the number of product types in the term (so that
/// every retyped term changes its parameters).
double ringProductCount(std::span<const std::size_t> types)
{
  double n = 0.0;
  for (std::size_t t : types) n += t >= 5 ? 1.0 : 0.0;
  return n;
}
std::vector<double> ringTemplateBond(std::span<const std::size_t> types)
{
  const double n = ringProductCount(types);
  return {96500.0 + 5000.0 * n, 1.47 - 0.01 * n};
}
std::vector<double> ringTemplateBend(std::span<const std::size_t> types)
{
  const double n = ringProductCount(types);
  return {62500.0 + 2000.0 * n, 112.0 - 1.5 * n};
}
std::vector<double> ringTemplateTorsion(std::span<const std::size_t> types)
{
  const double n = ringProductCount(types);
  return {0.0, 355.03 + 10.0 * n, -68.19, 791.32 + 20.0 * n};
}

/// The monomer: ring bonds 0-1, 1-2, 2-0, the linker 0-3 and the amine 3-4; harmonic bonds and bends, TraPPE
/// torsions, the 1-4 pairs (1,4) and (2,4) scaled by one half.
const std::vector<std::pair<std::size_t, std::size_t>> ringMonomerBonds{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {3, 4}};
const std::vector<std::array<std::size_t, 3>> ringMonomerBends{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {3, 0, 1}, {3, 0, 2}, {0, 3, 4}};
const std::vector<std::array<std::size_t, 4>> ringMonomerTorsions{{4, 3, 0, 1}, {4, 3, 0, 2}, {3, 0, 1, 2}, {3, 0, 2, 1}};
std::vector<double> ringMonomerBondParameters(std::size_t a, std::size_t b)
{
  const std::array<std::size_t, 2> key{std::min(a, b), std::max(a, b)};
  if (key == std::array<std::size_t, 2>{0, 3}) return {96500.0, 1.52};
  if (key == std::array<std::size_t, 2>{3, 4}) return {96500.0, 1.47};
  return {96500.0, 1.47};
}
std::vector<double> ringMonomerBendParameters(const std::array<std::size_t, 3>& bend)
{
  const bool ring = bend[0] <= 2 && bend[1] <= 2 && bend[2] <= 2;
  return {62500.0, ring ? 60.0 : (bend[1] == 3 ? 112.0 : 118.0)};
}
std::vector<double> ringMonomerTorsionParameters() { return {0.0, 355.03, -68.19, 791.32}; }

Component ringMonomer(const ForceField& forceField)
{
  ConnectivityTable connectivityTable(5);
  Potentials::IntraMolecularPotentials potentials{};
  for (const auto [a, b] : ringMonomerBonds)
  {
    connectivityTable[a, b] = true;
    connectivityTable[b, a] = true;
    potentials.bonds.push_back(BondPotential({a, b}, BondType::Harmonic, ringMonomerBondParameters(a, b)));
  }
  for (const std::array<std::size_t, 3>& bend : ringMonomerBends)
  {
    potentials.bends.push_back(BendPotential(bend, BendType::Harmonic, ringMonomerBendParameters(bend)));
  }
  for (const std::array<std::size_t, 4>& torsion : ringMonomerTorsions)
  {
    potentials.torsions.push_back(TorsionPotential(torsion, TorsionType::TraPPE, ringMonomerTorsionParameters()));
  }
  potentials.vanDerWaals = {VanDerWaalsPotential({1, 4}, VDWParameters::Type::LennardJones, {46.0, 3.95}, 0.5),
                            VanDerWaalsPotential({2, 4}, VDWParameters::Type::LennardJones, {60.0, 3.00}, 0.5)};
  potentials.coulombs = {CoulombPotential({1, 4}, CoulombType::Coulomb, 0.20, -0.30, 0.5),
                         CoulombPotential({2, 4}, CoulombType::Coulomb, -0.40, -0.30, 0.5)};

  const double3 linker = 1.52 * double3::normalize(double3(-0.45, -0.35, 0.82));
  const double3 amine = linker + 1.47 * double3::normalize(double3(-0.9, 0.3, 0.3));
  const std::vector<Atom> beads{Atom({0.0, 0.0, 0.0}, 0.20, 1.0, 0, 0, 0, false, false),
                                Atom({1.47, 0.0, 0.0}, 0.20, 1.0, 0, 1, 0, false, false),
                                Atom({0.735, 1.273, 0.0}, -0.40, 1.0, 0, 2, 0, false, false),
                                Atom(linker, 0.30, 1.0, 0, 3, 0, false, false),
                                Atom(amine, -0.30, 1.0, 0, 4, 0, false, false)};
  Component monomer =
      Component(forceField, "epoxyamine", 400.0, 3000000.0, 0.2, beads, connectivityTable, potentials, 5, 21);
  monomer.intra14VanDerWaalsScaling = 0.5;
  monomer.intra14ChargeChargeScaling = 0.5;
  monomer.reactiveGroups.push_back(
      ReactiveGroup{.type = "epoxide", .atoms = {{"C1", 0}, {"C2", 1}, {"O", 2}}, .initialState = "closed"});
  monomer.reactiveGroups.push_back(ReactiveGroup{.type = "amine", .atoms = {{"N", 4}}, .initialState = "primary"});
  return monomer;
}

/// The ring opening: the amine attacks the ring carbon; the C1-N bond forms, the C1-O bond breaks, C1, O and N
/// are retyped.
std::vector<ReactionRule> ringOpening()
{
  ReactionRule opening{};
  opening.name = "opening";
  opening.reactants = {{.label = "E", .type = "epoxide", .state = "closed"},
                       {.label = "A", .type = "amine", .state = "primary"}};
  opening.products = {{.state = "opened", .retype = {{"C1", "C_op"}, {"O", "O_op"}}},
                      {.state = "secondary", .retype = {{"N", "NH"}}}};
  opening.capture = std::array<RoleRef, 2>{RoleRef{"E", "C1"}, RoleRef{"A", "N"}};
  opening.formBonds = {{RoleRef{"E", "C1"}, RoleRef{"A", "N"}}};
  opening.breakBonds = {{RoleRef{"E", "C1"}, RoleRef{"E", "O"}}};
  // the second addition to the same amine (no retyping: the test of the ramp continuity compares the energy with
  // and without this reaction at its step)
  ReactionRule second{};
  second.name = "second_opening";
  second.reactants = {{.label = "E", .type = "epoxide", .state = "closed"},
                      {.label = "A", .type = "amine", .state = "secondary"}};
  second.products = {{.state = "opened"}, {.state = "tertiary"}};
  second.capture = std::array<RoleRef, 2>{RoleRef{"E", "C1"}, RoleRef{"A", "N"}};
  second.formBonds = {{RoleRef{"E", "C1"}, RoleRef{"A", "N"}}};
  second.breakBonds = {{RoleRef{"E", "C1"}, RoleRef{"E", "O"}}};
  return {opening, second};
}

/// 'molecules' epoxy-amine monomers; the first 2 * 'reactions' pair up (the amine of the second attacks the ring
/// of the first), at step 0 for the first half and at 'currentStep' / 2 for the second. With 'secondAdditions' > 0
/// the amines of the first 'secondAdditions' pairs (already secondary) also attack the rings of the molecules
/// 2 * 'reactions', 2 * 'reactions' + 1, ... at 'secondStep' ('currentStep' / 2 by default; the second rule, no
/// retyping).
System makeRingOpeningSystem(std::size_t seed, std::size_t molecules, std::size_t reactions, std::uint64_t currentStep,
                             std::size_t relaxationSteps, double bondDistance = 1.5, double boxLength = 40.0,
                             std::size_t secondAdditions = 0, std::optional<std::uint64_t> secondStep = std::nullopt)
{
  const ForceField forceField = ringForceField();
  const Component monomer = ringMonomer(forceField);
  System system = System(forceField, SimulationBox(boxLength, boxLength - 1.0, boxLength + 1.0), false, 300.0, 1e5,
                         1.0, {}, {monomer}, {}, {molecules}, 5);
  RandomNumber random(seed);
  randomizeConfiguration(system, random);

  PolymerLinkTable& table = system.polymerLinks;
  const std::size_t count = ringTypeNames.size();
  for (std::size_t a = 0; a < count; ++a)
  {
    for (std::size_t b = a; b < count; ++b)
    {
      const std::array<std::size_t, 2> key{a, b};
      table.templates.bonds.push_back(
          {.types = {ringTypeNames[a], ringTypeNames[b]}, .type = BondType::Harmonic, .parameters = ringTemplateBond(key)});
    }
  }
  for (std::size_t a = 0; a < count; ++a)
  {
    for (std::size_t b = 0; b < count; ++b)
    {
      for (std::size_t c = a; c < count; ++c)
      {
        const std::array<std::size_t, 3> key{a, b, c};
        table.templates.bends.push_back({.types = {ringTypeNames[a], ringTypeNames[b], ringTypeNames[c]},
                                         .type = BendType::Harmonic,
                                         .parameters = ringTemplateBend(key)});
      }
    }
  }
  for (std::size_t a = 0; a < count; ++a)
  {
    for (std::size_t b = 0; b < count; ++b)
    {
      for (std::size_t c = 0; c < count; ++c)
      {
        for (std::size_t d = 0; d < count; ++d)
        {
          const std::array<std::size_t, 4> key{a, b, c, d};
          const std::array<std::size_t, 4> reversedKey{d, c, b, a};
          if (reversedKey < key) continue;
          table.templates.torsions.push_back(
              {.types = {ringTypeNames[a], ringTypeNames[b], ringTypeNames[c], ringTypeNames[d]},
               .type = TorsionType::TraPPE,
               .parameters = ringTemplateTorsion(key)});
        }
      }
    }
  }
  table.templates.resolve(forceField);
  table.scheme.rules = ringOpening();
  table.scheme.resolve(system.components, forceField);
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.relaxationSteps = relaxationSteps;
  table.currentStep = currentStep;

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t n = 0; n < reactions; ++n)
  {
    const std::size_t a = 2 * n;
    const std::size_t b = 2 * n + 1;
    approachMolecules(system, a, 0, b, 4, bondDistance);
    const std::array<std::size_t, 2> groups{table.groupIndex(0, a, 0), table.groupIndex(0, b, 1)};
    table.apply(0, groups, {}, {}, forceField, system.components, atoms, std::span<Molecule>(system.moleculeData),
                system.simulationBox, n < reactions / 2 ? 0 : currentStep / 2);
  }
  const std::size_t n = monomer.atoms.size();
  for (std::size_t s = 0; s < secondAdditions; ++s)
  {
    const std::size_t b = 2 * s + 1;        // the amine, bonded to C1 of molecule 2 s
    const std::size_t c = 2 * reactions + s;  // the ring it attacks
    // away from both the CH2 of the amine and its first partner
    const double3 nitrogen = atoms[n * b + 4].position;
    const double3 direction = double3::normalize(-(double3::normalize(atoms[n * b + 3].position - nitrogen) +
                                                   double3::normalize(atoms[n * (2 * s) + 0].position - nitrogen)));
    approachMolecules(system, b, 4, c, 0, bondDistance, direction);
    const std::array<std::size_t, 2> groups{table.groupIndex(0, c, 0), table.groupIndex(0, b, 1)};
    table.apply(1, groups, {}, {}, forceField, system.components, atoms, std::span<Molecule>(system.moleculeData),
                system.simulationBox, secondStep.value_or(currentStep / 2));
  }
  table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, atoms);
  return system;
}

/// The reacted pair of 'makeRingOpeningSystem' as one molecule of a product component with the opened topology:
/// the atoms 0-4 of the attacked monomer (C1 and O retyped, the C1-O bond gone), 5-9 of the attacking one (N
/// retyped), the bond 0-9. The terms that contain a retyped atom or the new bond take the template parameters,
/// the others keep the monomer's; the 1-4 pairs of the new graph are scaled by one half.
System makeRingOpeningProduct(const System& reacted)
{
  const ForceField forceField = ringForceField();
  const std::size_t n = 10;
  const std::vector<std::size_t> types{5, 1, 6, 3, 4, 0, 1, 2, 3, 7};
  const std::set<std::size_t> retyped{0, 2, 9};
  std::vector<std::pair<std::size_t, std::size_t>> bonds{{0, 1}, {1, 2}, {0, 3}, {3, 4},
                                                         {5, 6}, {6, 7}, {7, 5}, {5, 8}, {8, 9}, {0, 9}};
  ConnectivityTable connectivity(n);
  for (const auto [a, b] : bonds)
  {
    connectivity[a, b] = true;
    connectivity[b, a] = true;
  }
  const auto neighbours = [&](std::size_t atom)
  {
    std::vector<std::size_t> result{};
    for (std::size_t k = 0; k < n; ++k)
    {
      if (connectivity[atom, k]) result.push_back(k);
    }
    return result;
  };
  const auto local = [](std::size_t atom) { return atom % 5; };
  // a term takes the template parameters when it contains a retyped atom or runs through the new bond (atoms of
  // both monomers)
  const auto fromTemplate = [&](std::span<const std::size_t> term)
  {
    return std::ranges::any_of(term, [&](std::size_t a) { return retyped.contains(a) || a / 5 != term.front() / 5; });
  };
  const auto typesOf = [&](std::span<const std::size_t> term)
  {
    std::vector<std::size_t> result{};
    for (std::size_t atom : term) result.push_back(types[atom]);
    return result;
  };

  Potentials::IntraMolecularPotentials potentials{};
  for (const auto [a, b] : bonds)
  {
    const std::array<std::size_t, 2> term{a, b};
    potentials.bonds.push_back(BondPotential(term, BondType::Harmonic,
                                             fromTemplate(term) ? ringTemplateBond(typesOf(term))
                                                                : ringMonomerBondParameters(local(a), local(b))));
  }
  for (std::size_t b = 0; b < n; ++b)
  {
    const std::vector<std::size_t> nb = neighbours(b);
    for (std::size_t i = 0; i < nb.size(); ++i)
    {
      for (std::size_t j = i + 1; j < nb.size(); ++j)
      {
        const std::array<std::size_t, 3> term{nb[i], b, nb[j]};
        std::vector<double> parameters{};
        if (fromTemplate(term))
        {
          parameters = ringTemplateBend(typesOf(term));
        }
        else
        {
          // the monomer lists the bend in one orientation; the parameters depend on the local atoms only
          const std::array<std::size_t, 3> localTerm{local(term[0]), local(term[1]), local(term[2])};
          parameters = ringMonomerBendParameters(localTerm);
        }
        potentials.bends.push_back(BendPotential(term, BendType::Harmonic, parameters));
      }
    }
  }
  for (const auto [b, c] : bonds)
  {
    for (std::size_t a : neighbours(b))
    {
      if (a == c) continue;
      for (std::size_t d : neighbours(c))
      {
        if (d == b || d == a) continue;
        const std::array<std::size_t, 4> term{a, b, c, d};
        potentials.torsions.push_back(TorsionPotential(
            term, TorsionType::TraPPE,
            fromTemplate(term) ? ringTemplateTorsion(typesOf(term)) : ringMonomerTorsionParameters()));
      }
    }
  }
  // (the explicit pair terms of a programmatic component fix the scaling of their pairs; the parameters come from
  // the force field)
  for (const std::array<std::size_t, 2>& pair : IntraMolecularExclusions::pairs14(connectivity))
  {
    potentials.vanDerWaals.push_back(VanDerWaalsPotential(pair, VDWParameters::Type::LennardJones, {50.0, 3.9}, 0.5));
    potentials.coulombs.push_back(CoulombPotential(pair, CoulombType::Coulomb, forceField.pseudoAtoms[types[pair[0]]].charge,
                                                   forceField.pseudoAtoms[types[pair[1]]].charge, 0.5));
  }

  std::span<const Atom> reactedAtoms = reacted.spanOfMoleculeAtoms();
  std::vector<Atom> beads{};
  for (std::size_t k = 0; k < n; ++k)
  {
    Atom atom = reactedAtoms[k];
    atom.type = static_cast<std::uint16_t>(types[k]);
    atom.charge = forceField.pseudoAtoms[types[k]].charge;
    atom.moleculeId = 0;
    beads.push_back(atom);
  }
  Component product = Component(forceField, "product", 400.0, 3000000.0, 0.2, beads, connectivity, potentials, 5, 21);
  product.intra14VanDerWaalsScaling = 0.5;
  product.intra14ChargeChargeScaling = 0.5;

  System system = System(forceField, reacted.simulationBox, false, 300.0, 1e5, 1.0, {}, {product}, {}, {1}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  double3 com{};
  double mass = 0.0;
  for (std::size_t k = 0; k < n; ++k)
  {
    atoms[k].position = reactedAtoms[k].position;
    const double m = forceField.pseudoAtoms[types[k]].mass;
    com += m * atoms[k].position;
    mass += m;
  }
  system.moleculeData[0].centerOfMassPosition = com / mass;
  return system;
}

void expectPolymerOverlayMatchesHost(PairDevice pairDevice, System& system, const char* label)
{
  // the same device pair kernel and mesh; the molecular terms and the polymer overlay on the host (double) and
  // on the device (single precision; the overlay from the wrapped slot positions)
  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  EXPECT_FALSE(host.usesDeviceBonded());
  EXPECT_FALSE(host.usesDevicePolymerOverlay());
  const RunningEnergy hostEnergy = host.computeGradients(system, true);
  const std::vector<double3> hostGradient = gradientsOf(system);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  EXPECT_TRUE(device.usesDeviceBonded());
  ASSERT_TRUE(device.usesDevicePolymerOverlay()) << device.writeStatus();
  const RunningEnergy energy = device.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);

  const auto near = [&](double a, double b, double scale, const char* what)
  {
    EXPECT_FALSE(std::isnan(a)) << label << ": " << what;
    EXPECT_NEAR(a, b, 1e-5 * std::max(std::abs(scale), 1e-3)) << label << ": " << what;
  };
  EXPECT_NE(hostEnergy.crossLink, 0.0) << label;
  near(energy.crossLink, hostEnergy.crossLink, hostEnergy.crossLink, "cross-link");
  near(energy.bond, hostEnergy.bond, hostEnergy.bond, "bond");
  near(energy.bend, hostEnergy.bend, hostEnergy.bend, "bend");
  near(energy.torsion, hostEnergy.torsion, hostEnergy.torsion, "torsion");
  near(energy.intraVDW, hostEnergy.intraVDW, hostEnergy.intraVDW, "intra VDW");
  near(energy.intraCoul, hostEnergy.intraCoul, hostEnergy.intraCoul, "intra Coulomb");
  near(energy.moleculeMoleculeVDW, hostEnergy.moleculeMoleculeVDW, hostEnergy.moleculeMoleculeVDW, "VDW");
  near(energy.moleculeMoleculeCharge, hostEnergy.moleculeMoleculeCharge, hostEnergy.moleculeMoleculeCharge,
       "Coulomb");
  near(energy.ewald_self + energy.ewald_exclusion, hostEnergy.ewald_self + hostEnergy.ewald_exclusion,
       hostEnergy.ewald_self + hostEnergy.ewald_exclusion, "Ewald exclusion");
  near(energy.potentialEnergy(), hostEnergy.potentialEnergy(), hostEnergy.potentialEnergy(), "potential");
  EXPECT_LT(rmsDifference(gradient, hostGradient), 1e-5 * rmsNorm(hostGradient)) << label;
  EXPECT_LT(maxAbsDifference(device.molecularPressureTensor(), host.molecularPressureTensor()),
            1e-5 * maxAbs(host.molecularPressureTensor()))
      << label;

  // deterministic
  const RunningEnergy again = device.computeGradients(system, true);
  EXPECT_EQ(again.potentialEnergy(), energy.potentialEnergy()) << label;
  EXPECT_EQ(again.crossLink, energy.crossLink) << label;
  EXPECT_EQ(rmsDifference(gradientsOf(system), gradient), 0.0) << label;
}
}  // namespace

TEST_P(SpatialDecompositionDevice, engine_device_polymer_overlay_matches_host)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";

  // all ramps complete: the bonds, bends and torsions across the polymer bonds, the overrides of the retyped
  // monomer terms, the pair corrections (excluded 1-2 / 1-3, scaled 1-4) and the stereo restraints at full strength
  {
    System system = makePolymerSystem(7, 20, 6, 1000, 100, 2000.0);
    ASSERT_EQ(system.polymerLinks.numberOfBonds(), 26uz);
    EXPECT_FALSE(system.polymerLinks.rampsActive());
    EXPECT_FALSE(system.polymerLinks.stereoCentres.empty());
    expectPolymerOverlayMatchesHost(pairDevice, system, "complete ramps");
  }
  // mid-ramp: half of the bonds at lambda = 0.6, the other half at lambda = 0.1 (the retype corrections stay at 1)
  {
    System system = makePolymerSystem(11, 20, 6, 60, 100, 2000.0);
    EXPECT_TRUE(system.polymerLinks.rampsActive());
    expectPolymerOverlayMatchesHost(pairDevice, system, "active ramps");
  }
  // no ramp, no stereo restraint (the centres are listed but the restraint is off)
  {
    System system = makePolymerSystem(13, 12, 0, 0, 0, 0.0);
    EXPECT_EQ(system.polymerLinks.chiralForceConstant, 0.0);
    expectPolymerOverlayMatchesHost(pairDevice, system, "no ramp");
  }
}

TEST_P(SpatialDecompositionDevice, engine_device_polymer_overlay_follows_the_ramp_and_new_bonds)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  System system = makePolymerSystem(17, 16, 4, 20, 100, 2000.0);

  SpatialDecompositionSettings hostSettings = settingsFor(4, 1.5, 0.5);
  hostSettings.pairDevice = pairDevice;
  hostSettings.deviceBonded = false;
  SpatialDecompositionForceEngine host(hostSettings);
  host.initialize(system);
  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = pairDevice;
  SpatialDecompositionForceEngine device(settings);
  device.initialize(system);
  ASSERT_TRUE(device.usesDevicePolymerOverlay()) << device.writeStatus();

  const auto compare = [&](const char* label)
  {
    const RunningEnergy hostEnergy = host.computeGradients(system, true);
    const std::vector<double3> hostGradient = gradientsOf(system);
    const RunningEnergy energy = device.computeGradients(system, true);
    EXPECT_NEAR(energy.crossLink, hostEnergy.crossLink, 1e-5 * std::abs(hostEnergy.crossLink)) << label;
    EXPECT_NEAR(energy.potentialEnergy(), hostEnergy.potentialEnergy(), 1e-5 * std::abs(hostEnergy.potentialEnergy()))
        << label;
    EXPECT_LT(rmsDifference(gradientsOf(system), hostGradient), 1e-5 * rmsNorm(hostGradient)) << label;
    return hostEnergy.crossLink;
  };

  // the ramp advances with the step the driver sets (the lambda buffer is rewritten, the topology stays)
  const double start = compare("step 20");
  system.polymerLinks.currentStep = 70;
  const double later = compare("step 70");
  EXPECT_GT(std::abs(later), std::abs(start));
  system.polymerLinks.currentStep = 500;
  compare("step 500");

  // a reaction: two more molecules react (retyping, a new bond), the driver's topologyChanged() on both engines
  PolymerLinkTable& table = system.polymerLinks;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::size_t a = 40;
  const std::size_t b = 41;
  approachMolecules(system, a, 1, b, 0);
  const std::array<std::size_t, 1> initiated{a};
  table.apply(0, initiated, {}, {}, system.forceField, system.components, atoms,
              std::span<Molecule>(system.moleculeData), system.simulationBox, 500);
  const std::array<std::size_t, 2> pair{a, b};
  const std::array<std::int8_t, 1> chirality{1};
  const std::array<std::int8_t, 1> next{0};
  table.apply(1, pair, chirality, next, system.forceField, system.components, atoms,
              std::span<Molecule>(system.moleculeData), system.simulationBox, 500);
  table.rebuild(system.forceField, system.components, system.numberOfMoleculesPerComponent, atoms);
  host.topologyChanged(system);
  device.topologyChanged(system);
  system.polymerLinks.currentStep = 560;
  compare("after the reaction");
  system.polymerLinks.currentStep = 2000;
  compare("after the reaction, ramp complete");
}

// The pairs a polymer bond excludes (1-2 and 1-3 through the bond) are left out of the pair lists and the overlay
// supplies their whole ramped value, lambda value(wanted) + (1 - lambda) value(counted); with the pairs in the
// lists and the overlay adding lambda (value(wanted) - value(counted)) the total is the same in exact arithmetic
// (the difference is the single-precision cancellation of the 1e8 K/A repulsion of the bonded pairs, which gave
// the device pair kernel a net force per bond). Host engine, double precision, with and without ramps.
TEST(spatial_decomposition, engine_polymer_unlisted_pairs_match_the_listed_formulation)
{
  for (const std::uint64_t currentStep : {1000uz, 60uz})
  {
    System system = makePolymerSystem(7, 20, 6, currentStep, 100, 2000.0);
    PolymerLinkTable& table = system.polymerLinks;
    ASSERT_FALSE(table.unlisted.empty());
    std::size_t unlistedCount = 0;
    for (const PolymerLinkTable::PairCorrection& pair : table.pairCorrections)
    {
      if (pair.listed) continue;
      ++unlistedCount;
      // the first entry of the chain of a pair the lists held at full strength that the topology excludes
      EXPECT_EQ(pair.counted, PolymerLinkTable::PairState::Full);
      EXPECT_TRUE(table.unlisted.contains(pair.atoms[0], pair.atoms[1]));
      EXPECT_TRUE(table.unlisted.contains(pair.atoms[1], pair.atoms[0]));
    }
    EXPECT_EQ(unlistedCount, table.unlisted.partner.size() / 2);
    // the bond itself is among them
    for (const PolymerLinkTable::LinkBondTerm& bond : table.linkBonds)
    {
      EXPECT_TRUE(table.unlisted.contains(bond.atoms[0], bond.atoms[1]));
    }

    SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
    settings.pairDevice = PairDevice::CPU;
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);
    const std::vector<double3> gradient = gradientsOf(system);
    const double3x3 pressure = engine.molecularPressureTensor();

    // the listed formulation: every pair in the lists, the overlay the difference
    for (PolymerLinkTable::PairCorrection& pair : table.pairCorrections) pair.listed = true;
    table.unlisted = PairListExclusions{};
    SpatialDecompositionForceEngine reference(settings);
    reference.initialize(system);
    const RunningEnergy referenceEnergy = reference.computeGradients(system, true);
    const std::vector<double3> referenceGradient = gradientsOf(system);

    const auto near = [&](double a, double b, const char* what)
    { EXPECT_NEAR(a, b, 1e-9 * std::max(std::abs(b), 1.0)) << what << " at step " << currentStep; };
    near(energy.potentialEnergy(), referenceEnergy.potentialEnergy(), "potential");
    near(energy.moleculeMoleculeVDW, referenceEnergy.moleculeMoleculeVDW, "VDW");
    near(energy.moleculeMoleculeCharge, referenceEnergy.moleculeMoleculeCharge, "Coulomb");
    near(energy.ewald_exclusion, referenceEnergy.ewald_exclusion, "Ewald exclusion");
    near(energy.crossLink, referenceEnergy.crossLink, "cross-link");
    EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-9 * rmsNorm(referenceGradient)) << currentStep;
    EXPECT_LT(maxAbsDifference(pressure, reference.molecularPressureTensor()),
              1e-9 * maxAbs(reference.molecularPressureTensor()))
        << currentStep;
  }
}

// A broken bond ('BreakBonds' of a reaction rule): the bond and the bends and torsions of the molecule through it
// are overridden to zero, the connectivity the overlay follows loses the bond, and the pairs the molecule's
// exclusions withhold that the new topology no longer excludes receive their interaction. With the ramp complete
// the reacted pair has to have the energy and the gradients of one molecule with the product topology.
TEST(spatial_decomposition, engine_polymer_broken_bond_matches_the_product_topology)
{
  System system = makeRingOpeningSystem(5, 2, 1, 1000, 100);
  PolymerLinkTable& table = system.polymerLinks;
  ASSERT_EQ(table.numberOfBonds(), 1uz);
  ASSERT_EQ(table.numberOfBrokenBonds(), 1uz);
  EXPECT_FALSE(table.rampsActive());
  EXPECT_EQ(table.brokenBonds[0].a.atomIndex, 0u);
  EXPECT_EQ(table.brokenBonds[0].b.atomIndex, 2u);

  // the connectivity: C1 (atom 0) lost the oxygen (2) and gained the amine nitrogen of the other molecule (9)
  const std::vector<std::uint32_t> neighbours = table.neighboursOf(0, system.components);
  EXPECT_FALSE(std::ranges::contains(neighbours, 2u));
  EXPECT_TRUE(std::ranges::contains(neighbours, 1u));
  EXPECT_TRUE(std::ranges::contains(neighbours, 3u));
  EXPECT_TRUE(std::ranges::contains(neighbours, 9u));
  EXPECT_FALSE(std::ranges::contains(table.neighboursOf(2, system.components), 0u));

  // the terms switched off: the bond, the bends 1-2-0, 2-0-1, 3-0-2 and the torsions 3-0-2-1, 4-3-0-2
  std::size_t removedBonds = 0, removedBends = 0, removedTorsions = 0;
  for (const PolymerLinkTable::OverrideBondTerm& term : table.overrideBonds)
  {
    if (term.replacement.type != BondType::None) continue;
    ++removedBonds;
    EXPECT_EQ(std::min(term.atoms[0], term.atoms[1]), 0u);
    EXPECT_EQ(std::max(term.atoms[0], term.atoms[1]), 2u);
  }
  for (const PolymerLinkTable::OverrideBendTerm& term : table.overrideBends)
  {
    if (term.replacement.type == BendType::Fixed) ++removedBends;
  }
  for (const PolymerLinkTable::OverrideTorsionTerm& term : table.overrideTorsions)
  {
    if (term.replacement.type == TorsionType::Fixed) ++removedTorsions;
  }
  EXPECT_EQ(removedBonds, 1uz);
  EXPECT_EQ(removedBends, 3uz);
  EXPECT_EQ(removedTorsions, 2uz);
  // no bend or torsion of the new bond runs through the broken one (the oxygen is reached through C2 only:
  // the torsion O-C2-C1-N exists, O-C1-N does not)
  const auto throughBroken = [](std::span<const std::uint32_t> atoms)
  {
    for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
    {
      if ((atoms[i] == 0u && atoms[i + 1] == 2u) || (atoms[i] == 2u && atoms[i + 1] == 0u)) return true;
    }
    return false;
  };
  bool oxygenTorsion = false;
  for (const PolymerLinkTable::LinkBendTerm& term : table.linkBends)
  {
    EXPECT_FALSE(throughBroken(term.atoms)) << "bend through the broken bond";
    EXPECT_FALSE(std::ranges::contains(term.atoms, 2u)) << "bend O-C1-N";
  }
  for (const PolymerLinkTable::LinkTorsionTerm& term : table.linkTorsions)
  {
    EXPECT_FALSE(throughBroken(term.atoms)) << "torsion through the broken bond";
    if (std::ranges::contains(term.atoms, 2u)) oxygenTorsion = true;
  }
  EXPECT_TRUE(oxygenTorsion);

  // the released pairs: the linker - oxygen (3, 2) was 1-3 (excluded) and is 1-4 now; the oxygen - nitrogen (2, 4)
  // was 1-4 (scaled) and is 1-5 now; C1 - O (0, 2) is 1-3 through C2 and stays excluded
  bool linkerOxygen = false, oxygenNitrogen = false;
  for (const PolymerLinkTable::PairCorrection& pair : table.pairCorrections)
  {
    EXPECT_NE(pair.atoms, (std::array<std::uint32_t, 2>{0u, 2u}));
    if (pair.atoms == std::array<std::uint32_t, 2>{2u, 3u})
    {
      linkerOxygen = true;
      EXPECT_EQ(pair.counted, PolymerLinkTable::PairState::Excluded);
      EXPECT_EQ(pair.wanted, PolymerLinkTable::PairState::Scaled);
      EXPECT_TRUE(pair.listed);
    }
    if (pair.atoms == std::array<std::uint32_t, 2>{2u, 4u})
    {
      oxygenNitrogen = true;
      EXPECT_EQ(pair.counted, PolymerLinkTable::PairState::Scaled);
      EXPECT_EQ(pair.wanted, PolymerLinkTable::PairState::Full);
      EXPECT_TRUE(pair.listed);
    }
  }
  EXPECT_TRUE(linkerOxygen);
  EXPECT_TRUE(oxygenNitrogen);

  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = PairDevice::CPU;
  SpatialDecompositionForceEngine engine(settings);
  engine.initialize(system);
  const RunningEnergy energy = engine.computeGradients(system, true);
  const std::vector<double3> gradient = gradientsOf(system);
  EXPECT_NE(energy.crossLink, 0.0);

  System product = makeRingOpeningProduct(system);
  SpatialDecompositionForceEngine reference(settings);
  reference.initialize(product);
  const RunningEnergy referenceEnergy = reference.computeGradients(product, true);
  const std::vector<double3> referenceGradient = gradientsOf(product);
  EXPECT_EQ(referenceEnergy.crossLink, 0.0);
  EXPECT_NEAR(energy.potentialEnergy(), referenceEnergy.potentialEnergy(),
              1e-8 * std::max(std::abs(referenceEnergy.potentialEnergy()), 1.0));
  EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-8 * rmsNorm(referenceGradient));

  // the overlay is not a no-op: without the released pairs and the removed terms the energy differs
  System untouched = makeRingOpeningSystem(5, 2, 1, 1000, 100);
  PolymerLinkTable& partial = untouched.polymerLinks;
  std::erase_if(partial.overrideBonds, [](const auto& t) { return t.replacement.type == BondType::None; });
  std::erase_if(partial.overrideBends, [](const auto& t) { return t.replacement.type == BendType::Fixed; });
  std::erase_if(partial.overrideTorsions, [](const auto& t) { return t.replacement.type == TorsionType::Fixed; });
  std::erase_if(partial.pairCorrections, [](const auto& p) { return p.counted == PolymerLinkTable::PairState::Excluded; });
  SpatialDecompositionForceEngine partialEngine(settings);
  partialEngine.initialize(untouched);
  const RunningEnergy partialEnergy = partialEngine.computeGradients(untouched, true);
  EXPECT_GT(std::abs(partialEnergy.potentialEnergy() - referenceEnergy.potentialEnergy()),
            1e-3 * std::abs(referenceEnergy.potentialEnergy()));
}

// A second bond at an atom of an earlier bond: the pairs of the earlier bond are within three bonds of the new one
// as well, but their treatment does not change, so they keep the ramp of the earlier bond (a restarted ramp would
// switch the full repulsion of two atoms at the bond length back on). The pair corrections form chains in the
// order of the events, and at the step of the second reaction (its ramp at zero) the energy and the forces are
// those of the system with the first reaction only.
TEST(spatial_decomposition, engine_polymer_second_addition_keeps_the_ramps_of_the_first_bond)
{
  // molecules 0-1 and 2-3 pair up (at step 0 and 30); molecule 1's amine (atom 9), bonded to C1 of molecule 0
  // (atom 0) at step 0, also bonds to C1 of molecule 4 (atom 20) at step 30; the ramps take 100 steps
  System system = makeRingOpeningSystem(5, 5, 2, 60, 100, 1.5, 40.0, 1, 30);
  PolymerLinkTable& table = system.polymerLinks;
  ASSERT_EQ(table.numberOfBonds(), 3uz);
  ASSERT_EQ(table.numberOfBrokenBonds(), 3uz);
  EXPECT_EQ(table.bonds[0].formedAtStep, 0u);
  EXPECT_EQ(table.bonds[2].formedAtStep, 30u);
  EXPECT_NEAR(table.rampFactor(0), 0.6, 1e-12);
  EXPECT_NEAR(table.rampFactor(30), 0.3, 1e-12);

  std::map<std::array<std::uint32_t, 2>, std::vector<PolymerLinkTable::PairCorrection>> chains{};
  for (const PolymerLinkTable::PairCorrection& pair : table.pairCorrections) chains[pair.atoms].push_back(pair);
  for (const auto& [key, chain] : chains)
  {
    for (std::size_t i = 0; i < chain.size(); ++i)
    {
      // an entry changes the treatment, continues the previous one, and only the first may be the unlisted form
      EXPECT_FALSE(chain[i].wanted == chain[i].counted &&
                   chain[i].wantedScaling.scalingVDW == chain[i].countedScaling.scalingVDW &&
                   chain[i].wantedScaling.scalingCoulomb == chain[i].countedScaling.scalingCoulomb)
          << "pair " << key[0] << "-" << key[1];
      if (i > 0)
      {
        EXPECT_GT(chain[i].step, chain[i - 1].step) << "pair " << key[0] << "-" << key[1];
        EXPECT_EQ(chain[i].counted, chain[i - 1].wanted) << "pair " << key[0] << "-" << key[1];
        EXPECT_TRUE(chain[i].listed) << "pair " << key[0] << "-" << key[1];
      }
    }
  }
  const auto chainOf = [&](std::uint32_t a, std::uint32_t b) -> const std::vector<PolymerLinkTable::PairCorrection>&
  { return chains[std::array<std::uint32_t, 2>{a, b}]; };
  // the first bond's own pair and its 1-3 pairs keep the step of the first bond
  ASSERT_EQ(chainOf(0u, 9u).size(), 1uz);
  EXPECT_EQ(chainOf(0u, 9u)[0].step, 0u);
  EXPECT_EQ(chainOf(0u, 9u)[0].wanted, PolymerLinkTable::PairState::Excluded);
  EXPECT_FALSE(chainOf(0u, 9u)[0].listed);
  ASSERT_EQ(chainOf(1u, 9u).size(), 1uz);  // C2 of molecule 0 - N
  EXPECT_EQ(chainOf(1u, 9u)[0].step, 0u);
  ASSERT_EQ(chainOf(0u, 8u).size(), 1uz);  // C1 of molecule 0 - CH2 of the amine
  EXPECT_EQ(chainOf(0u, 8u)[0].step, 0u);
  // the second bond's pairs have its step, among them C1 - N - C1 (1-3 through the nitrogen)
  ASSERT_EQ(chainOf(9u, 20u).size(), 1uz);
  EXPECT_EQ(chainOf(9u, 20u)[0].step, 30u);
  EXPECT_FALSE(chainOf(9u, 20u)[0].listed);
  ASSERT_EQ(chainOf(0u, 20u).size(), 1uz);
  EXPECT_EQ(chainOf(0u, 20u)[0].step, 30u);
  EXPECT_EQ(chainOf(0u, 20u)[0].wanted, PolymerLinkTable::PairState::Excluded);
  EXPECT_EQ(chainOf(0u, 20u)[0].counted, PolymerLinkTable::PairState::Full);

  // the overrides of the terms of molecule 1 that contain the nitrogen (retyped at step 0; the second addition
  // changes the state of its group again but not its type) keep the step of the first reaction
  ASSERT_EQ(table.retypings.size(), 6uz);  // C1 and O of molecules 0 and 2, N of molecules 1 and 3
  std::size_t nitrogenOverrides = 0;
  for (const PolymerLinkTable::OverrideBendTerm& term : table.overrideBends)
  {
    if (!std::ranges::contains(term.atoms, 9u)) continue;
    ++nitrogenOverrides;
    EXPECT_EQ(term.step, 0u);
  }
  EXPECT_GT(nitrogenOverrides, 0uz);

  // at the step of the second addition its ramp is zero: the same energy and forces as without it (the same
  // positions: the fifth molecule is placed the same way, and the second rule does not retype)
  System both = makeRingOpeningSystem(5, 5, 2, 30, 100, 1.5, 40.0, 1, 30);
  System without = makeRingOpeningSystem(5, 5, 2, 30, 100);
  {
    std::span<const Atom> placed = both.spanOfMoleculeAtoms();
    std::span<Atom> atoms = without.spanOfMoleculeAtoms();
    ASSERT_EQ(placed.size(), atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      atoms[i].position = placed[i].position;
      EXPECT_EQ(atoms[i].type, placed[i].type);
    }
    without.moleculeData[4].centerOfMassPosition = both.moleculeData[4].centerOfMassPosition;
  }
  EXPECT_EQ(both.polymerLinks.rampFactor(30), 0.0);
  SpatialDecompositionSettings settings = settingsFor(4, 1.5, 0.5);
  settings.pairDevice = PairDevice::CPU;
  SpatialDecompositionForceEngine engine(settings);
  engine.initialize(both);
  const RunningEnergy energy = engine.computeGradients(both, true);
  const std::vector<double3> gradient = gradientsOf(both);
  SpatialDecompositionForceEngine reference(settings);
  reference.initialize(without);
  const RunningEnergy referenceEnergy = reference.computeGradients(without, true);
  const std::vector<double3> referenceGradient = gradientsOf(without);
  EXPECT_NEAR(energy.potentialEnergy(), referenceEnergy.potentialEnergy(),
              1e-8 * std::max(std::abs(referenceEnergy.potentialEnergy()), 1.0));
  EXPECT_LT(rmsDifference(gradient, referenceGradient), 1e-8 * rmsNorm(referenceGradient));
}

TEST(spatial_decomposition, polymer_ramp_work_is_kept_out_of_the_conserved_energy_drift)
{
  // one bond formed at step 50 with a 100-step ramp, the table at step 100; the table's step is moved by
  // hand to play the bookkeeping of the MD driver through the ramp, its end, and a later reaction
  System system = makeRingOpeningSystem(5, 2, 1, 100, 100);
  PolymerLinkTable& table = system.polymerLinks;
  ASSERT_TRUE(table.activated);
  ASSERT_EQ(table.bonds[0].formedAtStep, 50u);
  ASSERT_TRUE(table.rampsActive());

  system.runningEnergies = RunningEnergy{};
  system.runningEnergies.moleculeMoleculeVDW = -5.0e4;
  system.runningEnergies.translationalKineticEnergy = 2.0e4;
  system.resetConservedEnergyBookkeeping();
  const double start = system.conservedEnergy;
  EXPECT_DOUBLE_EQ(system.referenceEnergy, start);
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 0.0);
  EXPECT_TRUE(system.polymerRampWasActive);

  // a step inside the ramp: the whole change is ramp work, the drift stays zero
  system.runningEnergies.bond += 300.0;
  system.accountConservedEnergy();
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 300.0);
  EXPECT_DOUBLE_EQ(system.conservedEnergy - system.referenceEnergy, 0.0);
  EXPECT_DOUBLE_EQ(system.accumulatedDrift, 0.0);

  // the step that ends the ramp (lambda reaches one in its force evaluation) still counts as ramp work
  table.currentStep = 150;
  ASSERT_FALSE(table.rampsActive());
  system.runningEnergies.bond += 50.0;
  system.accountConservedEnergy();
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 350.0);
  EXPECT_DOUBLE_EQ(system.conservedEnergy - system.referenceEnergy, 0.0);
  EXPECT_FALSE(system.polymerRampWasActive);

  // steps without a ramp: the change is integration drift
  table.currentStep = 151;
  system.runningEnergies.bond += 2.0;
  system.accountConservedEnergy();
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 350.0);
  EXPECT_DOUBLE_EQ(system.conservedEnergy - system.referenceEnergy, 2.0);
  EXPECT_DOUBLE_EQ(system.accumulatedDrift, std::abs(2.0 / system.referenceEnergy));
  const double driftSoFar = system.accumulatedDrift;

  // a reaction: the jump of the potential of the new topology is ramp work; the drift carries over
  system.runningEnergies.bond += 1.0e4;
  system.accountTopologyChange();
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 350.0 + 1.0e4);
  EXPECT_DOUBLE_EQ(system.polymerReactionJumps, 1.0e4);
  EXPECT_DOUBLE_EQ(system.conservedEnergy - system.referenceEnergy, 2.0);
  EXPECT_DOUBLE_EQ(system.accumulatedDrift, driftSoFar);
  EXPECT_TRUE(system.polymerRampWasActive);
  EXPECT_DOUBLE_EQ(system.conservedEnergy, start + 350.0 + 2.0 + 1.0e4);

  // the first step after the reaction (the table has no new ramp here, the flag covers it)
  system.runningEnergies.bond -= 100.0;
  system.accountConservedEnergy();
  EXPECT_DOUBLE_EQ(system.polymerRampWork, 250.0 + 1.0e4);
  EXPECT_DOUBLE_EQ(system.conservedEnergy - system.referenceEnergy, 2.0);
  EXPECT_FALSE(system.polymerRampWasActive);

  // a system without polymer links: every change is drift
  System plain = makeRingOpeningSystem(5, 2, 0, 100, 100);
  plain.runningEnergies = RunningEnergy{};
  plain.runningEnergies.moleculeMoleculeVDW = -5.0e4;
  plain.resetConservedEnergyBookkeeping();
  plain.runningEnergies.bond += 300.0;
  plain.accountConservedEnergy();
  EXPECT_DOUBLE_EQ(plain.polymerRampWork, 0.0);
  EXPECT_DOUBLE_EQ(plain.conservedEnergy - plain.referenceEnergy, 300.0);
}

TEST_P(SpatialDecompositionDevice, engine_device_polymer_overlay_matches_host_with_broken_bonds)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  {
    System system = makeRingOpeningSystem(7, 24, 8, 1000, 100);
    ASSERT_EQ(system.polymerLinks.numberOfBrokenBonds(), 8uz);
    EXPECT_FALSE(system.polymerLinks.rampsActive());
    expectPolymerOverlayMatchesHost(pairDevice, system, "broken bonds, complete ramps");
  }
  {
    System system = makeRingOpeningSystem(11, 24, 8, 60, 100);
    EXPECT_TRUE(system.polymerLinks.rampsActive());
    expectPolymerOverlayMatchesHost(pairDevice, system, "broken bonds, active ramps");
  }
  {
    // second additions: two bonds at one nitrogen, chains of pair corrections
    System system = makeRingOpeningSystem(13, 24, 8, 60, 100, 1.5, 40.0, 4);
    ASSERT_EQ(system.polymerLinks.numberOfBonds(), 12uz);
    EXPECT_TRUE(system.polymerLinks.rampsActive());
    expectPolymerOverlayMatchesHost(pairDevice, system, "second additions, active ramps");
  }
}

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_polymer_overlay_nve)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  // the overlay on the device in both runs (host integration with device forces against the resident integrator);
  // a mid-ramp state that stays constant (the driver advances the step). Mid-ramp the pairs across the new bonds
  // still carry (1 - lambda) of their full-strength repulsion, so the bonds start near the capture distance of the
  // driver rather than at the bond length (where that repulsion is ~1e7 K per pair, not integrable)
  const auto makeSystem = [](std::size_t seed) { return makePolymerSystem(seed, 20, 6, 60, 100, 2000.0, 3.2); };
  expectResidentMatchesHost(makeSystem, pairDevice, false, 60, 0.0005, 0.3,
                            {.position = 1e-4, .velocity = 1e-4, .energy = 1e-4, .drift = 1e-5});
}

TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_polymer_overlay_nvt)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  const auto makeSystem = [](std::size_t seed) { return makePolymerSystem(seed, 20, 6, 1000, 100, 2000.0); };
  expectResidentMatchesHost(makeSystem, pairDevice, true, 40, 0.0005, 0.3,
                            {.position = 1e-4, .velocity = 1e-4, .energy = 1e-4, .drift = 1e-5});
}

// The resident integrator with broken bonds and second additions (chains of pair corrections and overrides), with
// the ramps active: the removed ring terms, the released pairs and the retyped terms on the device.
TEST_P(SpatialDecompositionDevice, resident_integrator_matches_host_polymer_overlay_with_broken_bonds)
{
  const PairDevice pairDevice = GetParam();
  if (!DeviceStep::available(pairDevice)) GTEST_SKIP() << "no " << pairDeviceName(pairDevice) << " device";
  const auto makeSystem = [](std::size_t seed) { return makeRingOpeningSystem(seed, 24, 8, 60, 100, 3.2, 40.0, 4); };
  {
    System system = makeSystem(7);
    ASSERT_EQ(system.polymerLinks.numberOfBrokenBonds(), 12uz);
    EXPECT_TRUE(system.polymerLinks.rampsActive());
  }
  expectResidentMatchesHost(makeSystem, pairDevice, true, 40, 0.0005, 0.3,
                            {.position = 1e-4, .velocity = 1e-4, .energy = 1e-4, .drift = 1e-5});
}
