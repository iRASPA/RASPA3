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
import framework;
import component;
import molecule;
import system;
import simulationbox;
import running_energy;
import energy_status;
import energy_status_inter;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import randomnumbers;
import mc_cell_list;
import mc_moves;
import mc_moves_move_types;
import mc_moves_probabilities;
import cbmc_external_energy;
import cbmc_grow_context;
import cbmc_results;
import interactions_neighbourhood;
import interactions_pair_kernel;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

// The Monte Carlo cell list replaces the O(N) loops of the single-molecule energy differences and of the CBMC
// trial energies by a sphere gather over the fine cells (cut-off / 3 wide) around a trial set, into a
// structure-of-arrays neighbourhood shared by all trial sets of a batch. These tests check that
//   - the grid rule (fine cells, collapsed directions, disabled list) is as documented,
//   - the buckets and back-references are consistent after builds and incremental updates,
//   - the sphere gather is complete and visits every atom once, with the right image, also in short, elongated
//     and triclinic boxes, and the pair enumeration covers every pair once,
//   - the cell-list and neighbourhood energy routines reproduce the brute-force routines (also with a skipped
//     molecule, and in a box too small for the list),
//   - the list stays consistent with the atom positions over many accepted MC moves of all kinds and the running
//     energies stay in step with the recomputed total energies,
//   - the Ewald difference restricted to the moved atoms equals the full-molecule difference.

namespace
{
ForceField makeWaterForceField(ForceField::ChargeMethod chargeMethod = ForceField::ChargeMethod::Ewald,
                               double cutOff = 9.0)
{
  ForceField forceField = ForceField(
      {PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true), PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true)},
      {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot,
      cutOff, cutOff, cutOff, false, false, true);
  forceField.chargeMethod = chargeMethod;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(10, 10, 10);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 100;
  return forceField;
}

Component makeWater(const ForceField& forceField, const MCMoveProbabilities& probabilities = MCMoveProbabilities())
{
  return Component(forceField, "H2O", 304.1282, 7377300.0, 0.22394,
                   {Atom(double3(0.00000, -0.06461, 0.00000), -0.84760, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false),
                    Atom(double3(-0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false)},
                   {}, {}, 5, 21, probabilities);
}

// Places the molecules of the system on a jittered lattice with random orientations (no overlaps).
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

ForceField makeChainForceField()
{
  ForceField forceField = ForceField(
      {{"CH3", false, 15.03452, 0.0, 0.0, 6, false}, {"CH2", false, 14.02658, 0.0, 0.0, 6, false}},
      {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(10, 10, 10);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 100;
  return forceField;
}

// A flexible, partially charged butane with the intramolecular terms the polymer moves need (bonds, bends, a
// torsion and the 1-4 pair exclusions), optionally with Monte Carlo move probabilities.
Component makeChain(const ForceField& forceField, const MCMoveProbabilities& probabilities = MCMoveProbabilities())
{
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
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VDWParameters::Type::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.25, 0.25, 0.5)};

  return Component(forceField, "butane", 425.0, 3796000.0, 0.199,
                   {Atom({-1.85, -0.7, -0.15}, 0.25, 1.0, 0, 0, 0, false, false),
                    Atom({-0.31, -0.7, -0.15}, -0.25, 1.0, 0, 1, 0, false, false),
                    Atom({0.32, 0.71, -0.15}, -0.25, 1.0, 0, 1, 0, false, false),
                    Atom({1.86, 0.71, 0.15}, 0.25, 1.0, 0, 0, 0, false, false)},
                   connectivityTable, potentials, 5, 21, probabilities);
}

// The Ewald structure factor and the running energies of the current configuration, as the Monte Carlo driver
// prepares them before the first move.
void prepareForMonteCarlo(System& system)
{
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
}

void expectSameEnergies(const std::optional<RunningEnergy>& a, const std::optional<RunningEnergy>& b)
{
  ASSERT_EQ(a.has_value(), b.has_value());
  if (!a.has_value()) return;
  const double scale = std::max(1.0, std::abs(a->moleculeMoleculeVDW) + std::abs(a->moleculeMoleculeCharge));
  EXPECT_NEAR(a->moleculeMoleculeVDW, b->moleculeMoleculeVDW, 1e-10 * scale);
  EXPECT_NEAR(a->moleculeMoleculeCharge, b->moleculeMoleculeCharge, 1e-10 * scale);
  EXPECT_NEAR(a->potentialEnergy(), b->potentialEnergy(), 1e-10 * scale);
}

// Trial positions of the molecule 'moleculeIndex' translated by 'shift' (for a rigid water that is a legal move).
std::vector<Atom> translatedMolecule(const System& system, std::size_t moleculeIndex, const double3& shift)
{
  const Molecule& molecule = system.moleculeData[moleculeIndex];
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> trial(atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex),
                          atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex + molecule.numberOfAtoms));
  for (Atom& atom : trial) atom.position += shift;
  return trial;
}

// The framework-molecule terms of two evaluations agree (the charge term to the precision of the Ewald real-space
// table against the closed form).
void expectSameFrameworkEnergies(const std::optional<RunningEnergy>& a, const std::optional<RunningEnergy>& b)
{
  ASSERT_EQ(a.has_value(), b.has_value());
  if (!a.has_value()) return;
  const double scaleVDW = std::max(1.0, std::abs(a->frameworkMoleculeVDW));
  const double scaleCharge = std::max(1.0, std::abs(a->frameworkMoleculeCharge));
  EXPECT_NEAR(a->frameworkMoleculeVDW, b->frameworkMoleculeVDW, 1e-10 * scaleVDW);
  EXPECT_NEAR(a->frameworkMoleculeCharge, b->frameworkMoleculeCharge, 1e-8 * scaleCharge);
  EXPECT_DOUBLE_EQ(a->moleculeMoleculeVDW, b->moleculeMoleculeVDW);
  EXPECT_DOUBLE_EQ(a->moleculeMoleculeCharge, b->moleculeMoleculeCharge);
}

// A rigidly displaced and rotated copy of molecule 'moleculeIndex' (about its first atom).
std::vector<Atom> displacedMolecule(const System& system, RandomNumber& random, std::size_t moleculeIndex,
                                    double maximumDisplacement)
{
  const Molecule& molecule = system.moleculeData[moleculeIndex];
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> trial(atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex),
                          atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex + molecule.numberOfAtoms));
  const double3 shift = maximumDisplacement * double3(2.0 * random.uniform() - 1.0, 2.0 * random.uniform() - 1.0,
                                                      2.0 * random.uniform() - 1.0);
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
  const double3 pivot = trial.front().position;
  for (Atom& atom : trial) atom.position = pivot + shift + rotation * (atom.position - pivot);
  return trial;
}
}  // namespace

namespace
{
// The grid rule of the cell list: floor(3 w / cutOff) cells along a direction of perpendicular width w, collapsed
// to one when a query of the cut-off radius would visit every cell of the direction anyway.
int expectedCells(double width, double cutOff)
{
  const int n = static_cast<int>(std::floor(MCCellList::maximumSubdivision * width / cutOff));
  if (n <= 1) return 1;
  const int range = static_cast<int>(std::ceil(cutOff * n / width));
  return 2 * range + 1 < n ? n : 1;
}
}  // namespace

TEST(mc_cell_list, grid_selection_and_disabled_for_small_boxes)
{
  // 40 A box, 9 A cut-off: floor(3 * 40 / 9) = 13 cells per direction, a query's stencil is 7 wide
  EXPECT_EQ(MCCellList::gridForSubdivision(SimulationBox(40.0, 40.0, 40.0), 9.0, 3), int3(13, 13, 13));
  EXPECT_TRUE(MCCellList::wouldBeEnabled(SimulationBox(40.0, 40.0, 40.0), 9.0));

  // 75.75 A box, 14 A cut-off (the pHDDA-10 case): 16x16x16
  EXPECT_EQ(MCCellList::gridForSubdivision(SimulationBox(75.75, 75.75, 75.75), 14.0, 3), int3(16, 16, 16));

  // 36 A, 12 A cut-off: 9 cells, a stencil of 7: still a reduction
  EXPECT_EQ(MCCellList::gridForSubdivision(SimulationBox(36.0, 36.0, 36.0), 12.0, 3), int3(9, 9, 9));

  // 30 A, 12 A cut-off: 7 cells, the stencil of 7 visits all of them, every direction collapses, list disabled
  EXPECT_EQ(MCCellList::gridForSubdivision(SimulationBox(30.0, 30.0, 30.0), 12.0, 3), int3(1, 1, 1));
  EXPECT_FALSE(MCCellList::wouldBeEnabled(SimulationBox(30.0, 30.0, 30.0), 12.0));

  // the water case of the examples: a 24.83 A box with a 12.41 A cut-off is disabled
  EXPECT_FALSE(MCCellList::wouldBeEnabled(SimulationBox(24.83, 24.83, 24.83), 12.41));

  // an elongated box: only the long direction is divided
  EXPECT_EQ(MCCellList::gridForSubdivision(SimulationBox(22.0, 22.0, 100.0), 9.0, 3), int3(1, 1, 33));
  EXPECT_TRUE(MCCellList::wouldBeEnabled(SimulationBox(22.0, 22.0, 100.0), 9.0));

  // a huge box is capped (the cells grow beyond cutOff / 3, which stays correct)
  const int3 huge = MCCellList::gridForSubdivision(SimulationBox(600.0, 600.0, 600.0), 9.0, 3);
  EXPECT_LE(static_cast<std::size_t>(huge.x) * static_cast<std::size_t>(huge.y) * static_cast<std::size_t>(huge.z),
            32768uz);
  EXPECT_GT(huge.x, 8);

  // a triclinic box uses the perpendicular widths, not the vector lengths
  const SimulationBox triclinic(40.0, 40.0, 40.0, 60.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0,
                                90.0 * std::numbers::pi / 180.0);
  const double3 widths = triclinic.perpendicularWidths();
  const int3 grid = MCCellList::gridForSubdivision(triclinic, 9.0, 3);
  for (std::size_t d = 0; d < 3; ++d)
  {
    EXPECT_EQ(grid[static_cast<int>(d)], expectedCells(widths[d], 9.0));
  }
  EXPECT_LT(widths.y, 40.0);
}

TEST(mc_cell_list, subdivision_follows_the_atom_density)
{
  // 48 A box, 12 A cut-off: 12^3 cells at 3 per cut-off, 8^3 at 2, 4^3 at 1
  const SimulationBox box(48.0, 48.0, 48.0);
  EXPECT_EQ(MCCellList::gridForSubdivision(box, 12.0, 3), int3(12, 12, 12));
  EXPECT_EQ(MCCellList::gridForSubdivision(box, 12.0, 2), int3(8, 8, 8));
  EXPECT_EQ(MCCellList::gridForSubdivision(box, 12.0, 1), int3(4, 4, 4));

  // dense (at least 4 atoms per cell of the finest grid): the finest grid
  EXPECT_EQ(MCCellList::gridFor(box, 12.0, 8000), (std::pair<int3, int>{int3(12, 12, 12), 3}));
  // 2800 atoms: 1.6 per fine cell, 5.5 per cell at 2 per cut-off
  EXPECT_EQ(MCCellList::gridFor(box, 12.0, 2800), (std::pair<int3, int>{int3(8, 8, 8), 2}));
  // 300 atoms: 4.7 per cell at 1 per cut-off
  EXPECT_EQ(MCCellList::gridFor(box, 12.0, 300), (std::pair<int3, int>{int3(4, 4, 4), 1}));
  // too few atoms for any grid: the coarsest one that divides a direction
  EXPECT_EQ(MCCellList::gridFor(box, 12.0, 10), (std::pair<int3, int>{int3(4, 4, 4), 1}));

  // 36 A box, 12 A cut-off: the grid at 1 per cut-off (3 cells) gives no reduction, so the coarsest enabled
  // grid is the one at 2 per cut-off
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(36.0, 36.0, 36.0), 12.0, 10), (std::pair<int3, int>{int3(6, 6, 6), 2}));
  // no grid divides a 30 A box
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(30.0, 30.0, 30.0), 12.0, 10000).first, int3(1, 1, 1));

  // the list built over a dilute configuration reports the coarser subdivision
  std::vector<Atom> atoms(2800);
  RandomNumber random(3);
  for (Atom& atom : atoms)
  {
    atom.position = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
  }
  MCCellList list;
  list.build(atoms, box, 12.0);
  EXPECT_TRUE(list.enabled);
  EXPECT_EQ(list.subdivision, 2);
  EXPECT_EQ(list.numberOfCells, int3(8, 8, 8));
}

namespace
{
// Checks 'forEachRecordWithin' against the definition: every atom whose nearest image to 'centre' lies within
// 'radius' is visited exactly once, with that image; no other atom is visited.
void expectCompleteGather(const MCCellList& list, const SimulationBox& box, std::span<const Atom> atoms,
                          const double3& centre, double radius)
{
  std::map<std::uint32_t, double3> visited;
  list.forEachRecordWithin(box, centre, radius,
                           [&](const MCCellList::Record& record, const double3& image)
                           {
                             EXPECT_FALSE(visited.contains(record.atomIndex))
                                 << "atom " << record.atomIndex << " visited twice";
                             visited[record.atomIndex] = image;
                           });
  for (std::uint32_t i = 0; i < atoms.size(); ++i)
  {
    const double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - centre);
    const bool inside = double3::dot(dr, dr) <= radius * radius;
    ASSERT_EQ(visited.contains(i), inside) << "atom " << i << (inside ? " missed" : " visited although outside");
    if (inside)
    {
      // the image handed out is a periodic image of the atom, and it is the one nearest to the centre
      const double3 image = visited[i];
      const double3 offset = image - atoms[i].position;
      const double3 lattice = box.inverseCell * offset;
      EXPECT_NEAR(lattice.x, std::round(lattice.x), 1e-9);
      EXPECT_NEAR(lattice.y, std::round(lattice.y), 1e-9);
      EXPECT_NEAR(lattice.z, std::round(lattice.z), 1e-9);
      EXPECT_NEAR((image - centre).length_squared(), double3::dot(dr, dr), 1e-9);
    }
  }
}
}  // namespace

TEST(mc_cell_list, build_and_incremental_updates_are_consistent)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {512}, 5);
  RandomNumber random(11);
  randomizeConfiguration(system, random);

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  MCCellList list;
  list.build(atoms, system.simulationBox, 9.0);
  EXPECT_TRUE(list.valid);
  EXPECT_TRUE(list.enabled);
  // 1536 atoms are too few for 13^3 or 8^3 cells: the grid at 1 cell per cut-off is used
  EXPECT_EQ(list.subdivision, 1);
  EXPECT_EQ(list.numberOfCells, int3(4, 4, 4));
  EXPECT_EQ(list.numberOfAtoms, atoms.size());
  EXPECT_TRUE(list.orthorhombic);
  EXPECT_TRUE(list.verify(atoms, system.simulationBox));
  EXPECT_TRUE(list.isCurrent(atoms, system.simulationBox, 9.0));
  EXPECT_FALSE(list.isCurrent(atoms, system.simulationBox, 10.0));
  EXPECT_FALSE(list.isCurrent(atoms, SimulationBox(41.0, 40.0, 40.0), 9.0));

  // the stencil of a cut-off query is 3 cells wide: 27 of the 64 cells
  EXPECT_EQ(list.stencilRange(0, 9.0), 1);
  EXPECT_NEAR(list.averageNeighbourhoodSize(), 27.0 / 64.0 * static_cast<double>(atoms.size()),
              0.05 * static_cast<double>(atoms.size()));

  // the sphere gather around a position (also one far outside the box) is complete and visits every atom once,
  // for radii below and above the list's cut-off
  for (std::size_t trial = 0; trial < 50; ++trial)
  {
    const double3 position =
        system.simulationBox.cell *
        double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
    expectCompleteGather(list, system.simulationBox, atoms, position, trial % 2 == 0 ? 9.0 : 11.5);
  }

  // move molecules around (also across cells and across the periodic boundary) and apply the updates
  for (std::size_t step = 0; step < 2000; ++step)
  {
    const std::size_t m = static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
    const Molecule& molecule = system.moleculeData[m];
    const double3 shift = 15.0 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k) atoms[molecule.atomIndex + k].position += shift;
    list.updateAtoms(atoms, molecule.atomIndex, molecule.numberOfAtoms);
    if (!list.valid)
    {
      // a bucket overflowed: the owner rebuilds
      list.build(atoms, system.simulationBox, 9.0);
    }
    ASSERT_TRUE(list.verify(atoms, system.simulationBox)) << "after update " << step;
  }
  EXPECT_GT(list.numberOfAtomUpdates, 0uz);

  // a position far outside the box is wrapped into the grid
  atoms[0].position += system.simulationBox.cell * double3(7.0, -3.0, 2.0);
  list.updateAtoms(atoms, 0, 1);
  EXPECT_TRUE(list.verify(atoms, system.simulationBox));
}

// Regression: a CIF cell with 90-degree angles, built through the angle constructor and the supercell scaling,
// used to carry off-diagonal elements of about 1e-15 (std::cos(pi/2) is not 0). A wrap through that matrix
// moved an in-box coordinate by one ulp per application, so a slot written at one time no longer matched a
// fresh wrap of the same atom and 'verify' aborted the run ("Monte Carlo cell list inconsistent").
TEST(mc_cell_list, right_angle_boxes_are_exactly_diagonal_and_in_box_positions_are_stored_unchanged)
{
  const double right = 90.0 * Units::DegreesToRadians;
  const SimulationBox unitCell(22.3869, 22.3869, 22.3869, right, right, right);
  const SimulationBox box = unitCell.scaled(int3(2, 2, 2));
  EXPECT_EQ(box.type, SimulationBox::Type::Rectangular);
  EXPECT_EQ(box.cell.ay, 0.0);
  EXPECT_EQ(box.cell.az, 0.0);
  EXPECT_EQ(box.cell.bx, 0.0);
  EXPECT_EQ(box.cell.bz, 0.0);
  EXPECT_EQ(box.cell.cx, 0.0);
  EXPECT_EQ(box.cell.cy, 0.0);
  EXPECT_EQ(box.inverseCell.bx, 0.0);
  EXPECT_EQ(box.inverseCell.cy, 0.0);
  EXPECT_DOUBLE_EQ(box.cell.ax, 2.0 * 22.3869);

  // two right angles of a monoclinic cell are exact as well
  const SimulationBox monoclinic(10.0, 12.0, 14.0, right, 100.0 * Units::DegreesToRadians, right);
  EXPECT_EQ(monoclinic.type, SimulationBox::Type::Triclinic);
  EXPECT_EQ(monoclinic.cell.bx, 0.0);
  EXPECT_EQ(monoclinic.cell.cy, 0.0);
  EXPECT_NE(monoclinic.cell.cx, 0.0);

  RandomNumber random(7);
  std::vector<Atom> atoms;
  for (std::size_t i = 0; i < 600; ++i)
  {
    Atom atom(box.cell * double3(random.uniform(), random.uniform(), random.uniform()), 0.0, 1.0, 0, 0, 0, false,
              false);
    atom.moleculeId = static_cast<std::uint32_t>(i);
    atoms.push_back(atom);
  }
  // the atom of the reported failure
  atoms[2].position = double3(39.721948621946979, 43.451152171399592, 43.892688375709476);

  MCCellList list;
  list.build(atoms, box, 12.0);
  EXPECT_TRUE(list.orthorhombic);
  EXPECT_EQ(list.numberOfCells, int3(7, 7, 7));
  EXPECT_TRUE(list.verify(atoms, box));

  // an in-box position is stored bitwise, and the wrap is idempotent also for positions outside the box
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    const std::size_t slot = static_cast<std::size_t>(list.cellOfAtom[i]) * list.capacity + list.slotOfAtom[i];
    EXPECT_EQ(list.slotX[slot], atoms[i].position.x);
    EXPECT_EQ(list.slotY[slot], atoms[i].position.y);
    EXPECT_EQ(list.slotZ[slot], atoms[i].position.z);
    const double3 outside = atoms[i].position + box.cell * double3(3.0, -2.0, 1.0);
    const double3 once = list.wrappedPosition(outside);
    EXPECT_EQ(list.wrappedPosition(once), once);
    EXPECT_EQ(list.cellIndexOf(outside), list.cellIndexOf(atoms[i].position));
  }

  // a stale slot that is off by a rounding error is tolerated, a stale slot that is really elsewhere is not
  atoms[2].position.x = std::nextafter(atoms[2].position.x, 0.0);
  atoms[2].position.z = std::nextafter(atoms[2].position.z, 100.0);
  EXPECT_TRUE(list.verify(atoms, box));
  atoms[2].position.x += 1e-6;
  EXPECT_FALSE(list.verify(atoms, box));
}

namespace
{
std::vector<Atom> randomAtoms(const SimulationBox& box, std::size_t count, RandomNumber& random)
{
  std::vector<Atom> atoms;
  atoms.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    // positions also outside the primary box, as unwrapped molecule atoms are
    Atom atom(
        box.cell * double3(random.uniform() * 2.0 - 0.5, random.uniform() * 2.0 - 0.5, random.uniform() * 2.0 - 0.5),
        0.0, 1.0, 0, 0, 0, false, false);
    atom.moleculeId = static_cast<std::uint32_t>(i);
    atoms.push_back(atom);
  }
  return atoms;
}
}  // namespace

TEST(mc_cell_list, gather_is_complete_in_short_elongated_and_triclinic_boxes)
{
  RandomNumber random(41);

  // 36 x 36 x 100 A, 12 A cut-off, 3000 atoms: 6 x 6 x 16 cells (2 per cut-off; too dilute for 9 x 9 x 25).
  // A radius of 12 fits the stencil (5 of 6 cells); 14 makes the stencil wrap around the short axes
  // (minimum-image placement); 20 exceeds half the short widths.
  {
    const SimulationBox box(36.0, 36.0, 100.0);
    const std::vector<Atom> atoms = randomAtoms(box, 3000, random);
    MCCellList list;
    list.build(atoms, box, 12.0);
    ASSERT_TRUE(list.enabled);
    EXPECT_EQ(list.subdivision, 2);
    EXPECT_EQ(list.numberOfCells, int3(6, 6, 16));
    for (std::size_t trial = 0; trial < 30; ++trial)
    {
      const double3 centre =
          box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      expectCompleteGather(list, box, atoms, centre, 12.0);
      expectCompleteGather(list, box, atoms, centre, 14.0);
      expectCompleteGather(list, box, atoms, centre, 20.0);
    }
  }

  // 22 x 22 x 100 A, 9 A cut-off: only the long direction is divided, the short ones are visited whole
  {
    const SimulationBox box(22.0, 22.0, 100.0);
    const std::vector<Atom> atoms = randomAtoms(box, 2000, random);
    MCCellList list;
    list.build(atoms, box, 9.0);
    ASSERT_TRUE(list.enabled);
    EXPECT_EQ(list.numberOfCells, int3(1, 1, 33));
    for (std::size_t trial = 0; trial < 30; ++trial)
    {
      const double3 centre =
          box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      expectCompleteGather(list, box, atoms, centre, 9.0);
      expectCompleteGather(list, box, atoms, centre, 10.5);
    }
  }

  // a triclinic box: the per-cell shifts are lattice vectors, the cell distance bound is the face-normal gap
  {
    const SimulationBox box(40.0, 40.0, 40.0, 60.0 * std::numbers::pi / 180.0, 80.0 * std::numbers::pi / 180.0,
                            70.0 * std::numbers::pi / 180.0);
    const std::vector<Atom> atoms = randomAtoms(box, 3000, random);
    MCCellList list;
    list.build(atoms, box, 9.0);
    ASSERT_TRUE(list.enabled);
    EXPECT_FALSE(list.orthorhombic);
    for (std::size_t trial = 0; trial < 30; ++trial)
    {
      const double3 centre =
          box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      expectCompleteGather(list, box, atoms, centre, 9.0);
      expectCompleteGather(list, box, atoms, centre, 13.0);
    }
  }
}

namespace
{
// Checks the all-images gather against the definition: every periodic image (lattice shifts up to two cells) of
// every atom within 'radius' of 'centre' is visited exactly once.
void expectCompleteImageGather(const MCCellList& list, const SimulationBox& box, std::span<const Atom> atoms,
                               const double3& centre, double radius)
{
  std::multiset<std::uint32_t> visited;
  std::vector<double3> images;
  list.forEachRecordWithin(
      box, centre, radius,
      [&](const MCCellList::Record& record, const double3& image)
      {
        visited.insert(record.atomIndex);
        images.push_back(image);
        EXPECT_LE((image - centre).length_squared(), radius * radius * (1.0 + 1e-12));
      },
      /*allImages=*/true);
  std::size_t expected = 0;
  for (std::uint32_t i = 0; i < atoms.size(); ++i)
  {
    std::size_t imagesOfAtom = 0;
    // the atoms (randomAtoms) and the centres may lie a box length outside the box: enumerate generously
    for (int kx = -4; kx <= 4; ++kx)
    {
      for (int ky = -4; ky <= 4; ++ky)
      {
        for (int kz = -4; kz <= 4; ++kz)
        {
          const double3 image = atoms[i].position + box.cell * double3(kx, ky, kz);
          if ((image - centre).length_squared() <= radius * radius) ++imagesOfAtom;
        }
      }
    }
    EXPECT_EQ(visited.count(i), imagesOfAtom) << "atom " << i;
    expected += imagesOfAtom;
  }
  EXPECT_EQ(images.size(), expected);
  // no two visits hand out the same image
  std::sort(images.begin(), images.end(),
            [](const double3& a, const double3& b) { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); });
  for (std::size_t k = 1; k < images.size(); ++k)
  {
    EXPECT_GT((images[k] - images[k - 1]).length_squared(), 1e-12);
  }
}
}  // namespace

TEST(mc_cell_list, all_images_gather_is_complete_beyond_half_the_box)
{
  RandomNumber random(53);

  // a box where one direction is divided and the others collapse: 36 x 36 x 100 A at a 12 A cut-off with a
  // gather radius beyond 18 A reaches the images across the short axes
  {
    const SimulationBox box(36.0, 36.0, 100.0);
    const std::vector<Atom> atoms = randomAtoms(box, 3000, random);
    MCCellList list;
    list.build(atoms, box, 12.0);
    ASSERT_TRUE(list.enabled);
    for (std::size_t trial = 0; trial < 20; ++trial)
    {
      const double3 centre =
          box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      expectCompleteImageGather(list, box, atoms, centre, 19.5);
      expectCompleteImageGather(list, box, atoms, centre, 24.0);
      // a radius at and beyond the box length (a batch spanning the box, as a random translation makes): the
      // images on both sides and further are within reach
      expectCompleteImageGather(list, box, atoms, centre, 36.5);
      expectCompleteImageGather(list, box, atoms, centre, 50.0);
    }
  }

  // a triclinic box
  {
    const SimulationBox box(30.0, 30.0, 30.0, 60.0 * std::numbers::pi / 180.0, 80.0 * std::numbers::pi / 180.0,
                            70.0 * std::numbers::pi / 180.0);
    const std::vector<Atom> atoms = randomAtoms(box, 2000, random);
    MCCellList list;
    list.build(atoms, box, 9.0);
    ASSERT_TRUE(list.enabled);
    for (std::size_t trial = 0; trial < 20; ++trial)
    {
      const double3 centre =
          box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      expectCompleteImageGather(list, box, atoms, centre, 14.0);
      expectCompleteImageGather(list, box, atoms, centre, 40.0);
    }
  }

  // the all-atom fallback (no cell list) through the neighbourhood gather in the small water box: with a batch
  // within the box every image within the gather radius is a record; with a batch spanning the box (a random
  // translation: new and old positions a box apart) the gather radius reaches the box width, the images
  // across the box are no longer the only further ones, and the gather falls back to one record per atom (its
  // nearest image) with the minimum image left to the kernels
  {
    const SimulationBox box(24.83, 24.83, 24.83);
    const std::vector<Atom> atoms = randomAtoms(box, 1500, random);
    for (std::size_t trial = 0; trial < 20; ++trial)
    {
      const double3 centre = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
      const bool spanning = trial % 2 == 1;
      const double batchRadius = spanning ? 13.0 : 1.1;
      const double radius = batchRadius + 12.41;
      Interactions::Neighbourhood neighbourhood;
      neighbourhood.gather(box, nullptr, atoms, centre, batchRadius, 12.41);
      EXPECT_FALSE(neighbourhood.nearestImageUnique);
      EXPECT_EQ(neighbourhood.imagesResolved, !spanning);
      std::size_t expected = 0;
      for (const Atom& atom : atoms)
      {
        if (spanning)
        {
          const double3 d = box.applyPeriodicBoundaryConditions(atom.position - centre);
          if (d.length_squared() <= radius * radius) ++expected;
          continue;
        }
        for (int kx = -4; kx <= 4; ++kx)
        {
          for (int ky = -4; ky <= 4; ++ky)
          {
            for (int kz = -4; kz <= 4; ++kz)
            {
              const double3 image = atom.position + box.cell * double3(kx, ky, kz);
              if ((image - centre).length_squared() <= radius * radius) ++expected;
            }
          }
        }
      }
      EXPECT_EQ(neighbourhood.size(), expected);
      for (std::size_t k = 0; k < neighbourhood.size(); ++k)
      {
        const double3 image(neighbourhood.x[k], neighbourhood.y[k], neighbourhood.z[k]);
        EXPECT_LE((image - centre).length_squared(), radius * radius * (1.0 + 1e-12));
      }
      // the padding lies far away
      for (std::size_t k = neighbourhood.size(); k < neighbourhood.paddedSize(); ++k)
      {
        EXPECT_GT(std::abs(neighbourhood.x[k] - centre.x), 1e5);
        EXPECT_EQ(neighbourhood.moleculeId[k], std::numeric_limits<std::uint32_t>::max());
      }
    }
  }
}

TEST(mc_cell_list, pair_enumeration_covers_every_pair_within_the_cut_off_once)
{
  RandomNumber random(43);
  for (const SimulationBox& box :
       {SimulationBox(40.0, 40.0, 40.0), SimulationBox(36.0, 36.0, 100.0), SimulationBox(22.0, 22.0, 100.0)})
  {
    const std::vector<Atom> atoms = randomAtoms(box, 1500, random);
    MCCellList list;
    list.build(atoms, box, 9.0);
    ASSERT_TRUE(list.enabled);

    std::set<std::pair<std::uint32_t, std::uint32_t>> pairs;
    list.forEachPairOnce(
        [&](const MCCellList::Record& a, const MCCellList::Record& b)
        {
          const auto key = std::minmax(a.atomIndex, b.atomIndex);
          EXPECT_TRUE(pairs.insert({key.first, key.second}).second) << "pair visited twice";
        });
    for (std::uint32_t i = 0; i < atoms.size(); ++i)
    {
      for (std::uint32_t j = i + 1; j < atoms.size(); ++j)
      {
        const double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
        if (double3::dot(dr, dr) < 81.0) EXPECT_TRUE(pairs.contains({i, j})) << "pair " << i << " " << j << " missed";
      }
    }
  }
}

TEST(mc_cell_list, energy_difference_matches_brute_force_rigid_water)
{
  for (ForceField::ChargeMethod chargeMethod :
       {ForceField::ChargeMethod::Ewald, ForceField::ChargeMethod::DampedShiftedForce})
  {
    ForceField forceField = makeWaterForceField(chargeMethod);
    Component water = makeWater(forceField);
    System system =
        System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {700}, 5);
    RandomNumber random(3);
    randomizeConfiguration(system, random);

    const MCCellList& list = system.cellList();
    EXPECT_TRUE(list.enabled);
    EXPECT_TRUE(system.verifyCellList());

    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    std::size_t overlaps = 0;
    for (std::size_t trial = 0; trial < 200; ++trial)
    {
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      const Molecule& molecule = system.moleculeData[m];
      // small, medium and large displacements: within the cell, into a neighbouring cell, anywhere
      const double scale = trial % 3 == 0 ? 0.5 : (trial % 3 == 1 ? 6.0 : 40.0);
      const double3 shift = scale * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
      const std::vector<Atom> trialAtoms = translatedMolecule(system, m, shift);
      std::span<const Atom> oldAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);

      const std::optional<RunningEnergy> bruteForce = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, atoms, trialAtoms, oldAtoms);
      const std::optional<RunningEnergy> withList = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, list, atoms, trialAtoms, oldAtoms);
      expectSameEnergies(bruteForce, withList);
      if (!bruteForce.has_value()) ++overlaps;
    }
    // the large displacements produce overlaps now and then, both routines must report them identically
    EXPECT_GT(overlaps, 0uz);
  }
}

TEST(mc_cell_list, cbmc_trial_energy_matches_brute_force)
{
  ForceField forceField = makeChainForceField();
  Component chain = makeChain(forceField);
  System system =
      System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain}, {}, {300}, 5);
  RandomNumber random(5);
  randomizeConfiguration(system, random);

  const MCCellList& list = system.cellList();
  EXPECT_TRUE(list.enabled);
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();

  for (std::size_t trial = 0; trial < 200; ++trial)
  {
    // trial beads of a (re)growing molecule at random places, several per call as CBMC does
    const std::size_t grown =
        static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
    std::vector<Atom> trialAtoms;
    const std::size_t count = 1 + static_cast<std::size_t>(random.uniform() * 4.0);
    for (std::size_t k = 0; k < count; ++k)
    {
      Atom atom = atoms[system.moleculeData[grown].atomIndex + (k % 4)];
      atom.position = system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      trialAtoms.push_back(atom);
    }

    // the full cut-offs and the reduced CBMC cut-offs (which must not exceed those of the list)
    for (double cutOff : {9.0, 6.5})
    {
      const std::optional<RunningEnergy> bruteForce = CBMC::computeInterMolecularEnergy(
          system.forceField, system.simulationBox, atoms, cutOff, cutOff, trialAtoms, grown);
      const std::optional<RunningEnergy> withList = CBMC::computeInterMolecularEnergy(
          system.forceField, system.simulationBox, list, atoms, cutOff, cutOff, trialAtoms, grown);
      expectSameEnergies(bruteForce, withList);

      // without a skipped molecule (a growing new molecule)
      const std::optional<RunningEnergy> bruteForceAll =
          CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atoms, cutOff, cutOff, trialAtoms);
      const std::optional<RunningEnergy> withListAll = CBMC::computeInterMolecularEnergy(
          system.forceField, system.simulationBox, list, atoms, cutOff, cutOff, trialAtoms);
      expectSameEnergies(bruteForceAll, withListAll);
    }
  }
}

TEST(mc_cell_list, batched_neighbourhood_matches_brute_force_for_every_trial_set)
{
  // a large box (cell list enabled) and the small water box of the examples (cell list disabled, the gather
  // runs over all atoms and the kernel applies the minimum image per pair)
  for (const SimulationBox& box : {SimulationBox(40.0, 40.0, 40.0), SimulationBox(24.83, 24.83, 24.83)})
  {
    ForceField forceField = makeWaterForceField(ForceField::ChargeMethod::Ewald, 12.0);
    Component water = makeWater(forceField);
    const std::size_t count = box.cell.ax > 30.0 ? 700 : 170;
    System system = System(forceField, box, false, 300.0, 1e5, 1.0, {}, {water}, {}, {count}, 5);
    RandomNumber random(47);
    randomizeConfiguration(system, random);

    const MCCellList& list = system.cellList();
    EXPECT_EQ(list.enabled, box.cell.ax > 30.0);
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();

    for (std::size_t round = 0; round < 40; ++round)
    {
      // a batch of trial orientations of one molecule around a common first bead, as a rigid regrow produces
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      const Molecule& molecule = system.moleculeData[m];
      const double3 anchor = system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      std::vector<std::vector<Atom>> trialSets;
      for (std::size_t k = 0; k < 8; ++k)
      {
        std::vector<Atom> set(atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex),
                              atoms.begin() + static_cast<std::ptrdiff_t>(molecule.atomIndex + molecule.numberOfAtoms));
        const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
        for (std::size_t a = 0; a < set.size(); ++a)
        {
          set[a].position = anchor + rotation * water.atoms[a].position;
        }
        trialSets.push_back(std::move(set));
      }
      std::vector<std::span<const Atom>> views(trialSets.begin(), trialSets.end());

      Interactions::Neighbourhood neighbourhood;
      const CBMC::GrowContext context = system.makeGrowContext();
      CBMC::gatherNeighbourhood(context, neighbourhood, views);
      EXPECT_EQ(neighbourhood.nearestImageUnique, box.cell.ax > 30.0);
      // in the small box the gather radius exceeds half the box and all images within it are gathered (see
      // all_images_gather_is_complete_beyond_half_the_box); in the large box every atom appears once. Either way
      // the kernels need no minimum image.
      EXPECT_TRUE(neighbourhood.imagesResolved);
      if (box.cell.ax > 30.0)
      {
        std::map<std::uint32_t, std::size_t> entriesPerMolecule;
        for (std::size_t k = 0; k < neighbourhood.size(); ++k) ++entriesPerMolecule[neighbourhood.moleculeId[k]];
        EXPECT_TRUE(std::ranges::all_of(entriesPerMolecule, [](const auto& e) { return e.second <= 3; }));
      }
      EXPECT_EQ(neighbourhood.paddedSize() % Interactions::Neighbourhood::padding, 0uz);
      EXPECT_TRUE(neighbourhood.fullyCoupled);
      EXPECT_LT(neighbourhood.batchRadius, 1.5);
      // water (plain Lennard-Jones, Ewald through the table) qualifies for the lane kernel
      EXPECT_TRUE(Interactions::neighbourhoodFastKernelApplies(system.forceField, neighbourhood, views[0]));

      for (std::span<const Atom> set : views)
      {
        ASSERT_TRUE(neighbourhood.covers(set, 12.0, 12.0, true));
        EXPECT_FALSE(neighbourhood.covers(set, 12.5, 12.0, true));
        for (double cutOff : {12.0, 7.0})
        {
          const std::optional<RunningEnergy> bruteForce =
              CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atoms, cutOff, cutOff, set);
          const std::optional<RunningEnergy> batched = CBMC::computeInterMolecularEnergy(
              system.forceField, system.simulationBox, neighbourhood, cutOff, cutOff, set);
          expectSameEnergies(bruteForce, batched);
        }
      }

      // the skipped molecule (identity change) is left out of the gather
      const std::size_t other =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      Interactions::Neighbourhood skipping;
      CBMC::gatherNeighbourhood(context.withSkippedMolecule(other), skipping, views);
      for (std::uint32_t id : skipping.moleculeId) EXPECT_NE(id, static_cast<std::uint32_t>(other));
      const std::optional<RunningEnergy> bruteForceSkipped = CBMC::computeInterMolecularEnergy(
          system.forceField, system.simulationBox, atoms, 12.0, 12.0, views[0], other);
      const std::optional<RunningEnergy> batchedSkipped =
          CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, skipping, 12.0, 12.0, views[0]);
      expectSameEnergies(bruteForceSkipped, batchedSkipped);
    }
  }
}

TEST(mc_cell_list, lane_kernel_matches_generic_kernel_and_brute_force)
{
  // shifted Lennard-Jones (a non-zero shift per pair) with Ewald in the small water box (all images gathered)
  // and in a large box (cell list)
  for (const SimulationBox& box : {SimulationBox(24.83, 24.83, 24.83), SimulationBox(40.0, 40.0, 40.0)})
  {
    ForceField forceField =
        ForceField({PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true),
                    PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true)},
                   {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)},
                   ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, /*shifted=*/true, false, true);
    forceField.automaticEwald = false;
    forceField.EwaldAlpha = 0.3;
    forceField.numberOfWaveVectors = int3(8, 8, 8);
    forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
    forceField.reciprocalIntegerCutOffSquared = 64;
    forceField.updateEwaldRealSpaceTable();
    EXPECT_NE(forceField(0, 0).shift, 0.0);

    Component water = makeWater(forceField);
    const std::size_t count = box.cell.ax > 30.0 ? 700 : 170;
    System system = System(forceField, box, false, 300.0, 1e5, 1.0, {}, {water}, {}, {count}, 5);
    RandomNumber random(59);
    randomizeConfiguration(system, random);
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();

    // the generic kernel through the exact erfc: the same force field with the table switched off
    ForceField exactForceField = system.forceField;
    exactForceField.useEwaldRealSpaceTable = false;
    exactForceField.updateEwaldRealSpaceTable();

    std::size_t overlaps = 0;
    for (std::size_t round = 0; round < 60; ++round)
    {
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      // near an existing molecule now and then, so that overlaps occur
      const double3 anchor =
          round % 4 == 0 ? atoms[system.moleculeData[m].atomIndex].position +
                               0.8 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5)
                         : system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      std::vector<Atom> set = translatedMolecule(system, m, double3(0.0, 0.0, 0.0));
      const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
      for (std::size_t a = 0; a < set.size(); ++a) set[a].position = anchor + rotation * water.atoms[a].position;
      // a new id: the trial molecule interacts with every background molecule, the anchor molecule included
      for (Atom& atom : set) atom.moleculeId = static_cast<std::uint32_t>(system.moleculeData.size() + 7);

      const auto [centre, radius] = Interactions::Neighbourhood::boundingSphere(set);
      Interactions::Neighbourhood neighbourhood;
      neighbourhood.gather(system.simulationBox, &system.cellList(), atoms, centre, radius, 12.0);
      EXPECT_TRUE(neighbourhood.imagesResolved);
      ASSERT_TRUE(Interactions::neighbourhoodFastKernelApplies(system.forceField, neighbourhood, set));
      EXPECT_FALSE(Interactions::neighbourhoodFastKernelApplies(exactForceField, neighbourhood, set));

      for (double cutOff : {12.0, 8.0})
      {
        const std::optional<RunningEnergy> lane = Interactions::computeNeighbourhoodEnergy(
            system.forceField, system.simulationBox, neighbourhood, set, cutOff, cutOff, 1.0, true);
        const std::optional<RunningEnergy> generic = Interactions::computeNeighbourhoodEnergy(
            exactForceField, system.simulationBox, neighbourhood, set, cutOff, cutOff, 1.0, true);
        const std::optional<RunningEnergy> bruteForce =
            CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atoms, cutOff, cutOff, set);
        expectSameEnergies(bruteForce, lane);
        ASSERT_EQ(lane.has_value(), generic.has_value());
        if (lane.has_value())
        {
          // the table reproduces erfc(alpha r) / r to 1e-9 relative
          const double scale =
              std::max(1.0, std::abs(lane->moleculeMoleculeVDW) + std::abs(lane->moleculeMoleculeCharge));
          EXPECT_NEAR(lane->moleculeMoleculeVDW, generic->moleculeMoleculeVDW, 1e-10 * scale);
          EXPECT_NEAR(lane->moleculeMoleculeCharge, generic->moleculeMoleculeCharge, 1e-8 * scale);
        }
        else
        {
          ++overlaps;
        }
        // the sign and the old-configuration mode (no overlap check) of the energy differences
        const std::optional<RunningEnergy> negated = Interactions::computeNeighbourhoodEnergy(
            system.forceField, system.simulationBox, neighbourhood, set, cutOff, cutOff, -1.0, false);
        ASSERT_TRUE(negated.has_value());
        if (lane.has_value())
        {
          EXPECT_DOUBLE_EQ(negated->moleculeMoleculeVDW, -lane->moleculeMoleculeVDW);
          EXPECT_DOUBLE_EQ(negated->moleculeMoleculeCharge, -lane->moleculeMoleculeCharge);
        }
      }
    }
    EXPECT_GT(overlaps, 0uz);

    // a fractional background molecule: the lane kernel declines, the generic one carries the soft-core terms
    {
      const Molecule& fractional = system.moleculeData[3];
      for (std::size_t k = 0; k < fractional.numberOfAtoms; ++k)
      {
        atoms[fractional.atomIndex + k].scalingVDW = 0.5;
        atoms[fractional.atomIndex + k].scalingCoulomb = 0.25;
        atoms[fractional.atomIndex + k].groupId = 1;
      }
      system.invalidateCellList();
      std::vector<Atom> set = translatedMolecule(system, 3, double3(0.0, 0.0, 0.0));
      for (Atom& atom : set) atom.position += double3(2.5, -1.0, 0.7);
      for (Atom& atom : set) atom.moleculeId = static_cast<std::uint32_t>(system.moleculeData.size() + 7);
      const auto [centre, radius] = Interactions::Neighbourhood::boundingSphere(set);
      Interactions::Neighbourhood neighbourhood;
      neighbourhood.gather(system.simulationBox, &system.cellList(), atoms, centre, radius, 12.0);
      EXPECT_FALSE(neighbourhood.fullyCoupled);
      EXPECT_FALSE(Interactions::neighbourhoodFastKernelApplies(system.forceField, neighbourhood, set));
      const std::optional<RunningEnergy> generic = Interactions::computeNeighbourhoodEnergy(
          system.forceField, system.simulationBox, neighbourhood, set, 12.0, 12.0, 1.0, true);
      const std::optional<RunningEnergy> bruteForce =
          CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atoms, 12.0, 12.0, set);
      expectSameEnergies(bruteForce, generic);
      ASSERT_TRUE(generic.has_value());
      EXPECT_NE(generic->dudlambdaVDW[0], 0.0);
      EXPECT_NEAR(generic->dudlambdaVDW[0], bruteForce->dudlambdaVDW[0], 1e-10 * std::abs(bruteForce->dudlambdaVDW[0]));
      EXPECT_NEAR(generic->dudlambdaCharge[0], bruteForce->dudlambdaCharge[0],
                  1e-10 * std::abs(bruteForce->dudlambdaCharge[0]));
    }
  }
}

TEST(mc_cell_list, all_atom_neighbourhood_with_minimum_image_matches_brute_force)
{
  // first beads at random places in boxes too small for the cell list: one unshifted copy of all atoms, the
  // minimum image per lane (rectangular and triclinic), the lane and the generic kernel against the brute force
  for (const SimulationBox& box : {SimulationBox(24.83, 24.83, 24.83),
                                   SimulationBox(27.0, 26.0, 28.0, 75.0 * std::numbers::pi / 180.0,
                                                 85.0 * std::numbers::pi / 180.0, 80.0 * std::numbers::pi / 180.0)})
  {
    ForceField forceField = makeWaterForceField(ForceField::ChargeMethod::Ewald, 12.0);
    Component water = makeWater(forceField);
    System system = System(forceField, box, false, 300.0, 1e5, 1.0, {}, {water}, {}, {170}, 5);
    RandomNumber random(61);
    randomizeConfiguration(system, random);
    EXPECT_FALSE(system.cellList().enabled);
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();

    ForceField exactForceField = system.forceField;
    exactForceField.useEwaldRealSpaceTable = false;
    exactForceField.updateEwaldRealSpaceTable();

    Interactions::Neighbourhood all;
    all.gatherAll(atoms, 12.0);
    EXPECT_FALSE(all.imagesResolved);
    EXPECT_EQ(all.size(), atoms.size());
    EXPECT_EQ(all.paddedSize() % Interactions::Neighbourhood::padding, 0uz);
    for (std::size_t k = all.size(); k < all.paddedSize(); ++k)
    {
      EXPECT_EQ(all.moleculeId[k], Interactions::Neighbourhood::paddingMoleculeId);
    }

    std::size_t overlaps = 0;
    for (std::size_t round = 0; round < 150; ++round)
    {
      // a water oxygen (or a hydrogen) anywhere in the box, also outside the box and close to other molecules
      Atom bead = atoms[random.uniform() < 0.5 ? 0 : 1];
      bead.moleculeId = static_cast<std::uint32_t>(system.moleculeData.size() + 3);
      if (round % 3 == 0)
      {
        const std::size_t m =
            static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
        bead.position = atoms[system.moleculeData[m].atomIndex].position +
                        1.2 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
      }
      else
      {
        bead.position = box.cell * double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0,
                                           random.uniform() * 3.0 - 1.0);
      }
      const std::span<const Atom> set(&bead, 1);
      ASSERT_TRUE(all.covers(set, 12.0, 12.0, true));
      ASSERT_TRUE(Interactions::neighbourhoodFastKernelApplies(system.forceField, all, set));

      const std::optional<RunningEnergy> lane =
          Interactions::computeNeighbourhoodEnergy(system.forceField, box, all, set, 12.0, 12.0, 1.0, true);
      const std::optional<RunningEnergy> generic =
          Interactions::computeNeighbourhoodEnergy(exactForceField, box, all, set, 12.0, 12.0, 1.0, true);
      const std::optional<RunningEnergy> bruteForce =
          CBMC::computeInterMolecularEnergy(system.forceField, box, atoms, 12.0, 12.0, set);
      expectSameEnergies(bruteForce, lane);
      // the generic kernel through the exact erfc: equal up to the interpolation error of the table
      ASSERT_EQ(generic.has_value(), lane.has_value());
      if (lane.has_value())
      {
        const double scale =
            std::max(1.0, std::abs(lane->moleculeMoleculeVDW) + std::abs(lane->moleculeMoleculeCharge));
        EXPECT_NEAR(lane->moleculeMoleculeVDW, generic->moleculeMoleculeVDW, 1e-10 * scale);
        EXPECT_NEAR(lane->moleculeMoleculeCharge, generic->moleculeMoleculeCharge, 1e-8 * scale);
      }
      else
      {
        ++overlaps;
      }
    }
    EXPECT_GT(overlaps, 0uz);

    // the first-bead batch of the CBMC growth takes this route with the system's (disabled) cell list and must
    // agree with the brute-force loop of a context without a list
    const CBMC::GrowContext context = system.makeGrowContext();
    ASSERT_NE(context.cellList, nullptr);
    const CBMC::GrowContext bruteForceContext = context.withCellList(nullptr);
    for (std::size_t round = 0; round < 10; ++round)
    {
      std::vector<Atom> trialPositions;
      for (std::size_t k = 0; k < 10; ++k)
      {
        Atom bead = atoms[0];
        bead.moleculeId = static_cast<std::uint32_t>(system.moleculeData.size() + 3);
        bead.position = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
        trialPositions.push_back(bead);
      }
      const std::vector<CBMC::FirstBeadTrial> shared =
          CBMC::computeExternalNonOverlappingEnergies(context, water, trialPositions);
      const std::vector<CBMC::FirstBeadTrial> reference =
          CBMC::computeExternalNonOverlappingEnergies(bruteForceContext, water, trialPositions);
      ASSERT_EQ(shared.size(), reference.size());
      for (std::size_t k = 0; k < shared.size(); ++k)
      {
        EXPECT_EQ(shared[k].position.position.x, reference[k].position.position.x);
        expectSameEnergies(std::optional{reference[k].energy}, std::optional{shared[k].energy});
      }
    }
  }
}

TEST(mc_cell_list, ewald_difference_of_moved_atoms_matches_full_molecule)
{
  for (ForceField::ChargeMethod chargeMethod :
       {ForceField::ChargeMethod::Ewald, ForceField::ChargeMethod::DampedShiftedForce})
  {
    ForceField forceField = makeChainForceField();
    forceField.chargeMethod = chargeMethod;
    Component chain = makeChain(forceField);
    System system =
        System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain}, {}, {64}, 5);
    RandomNumber random(9);
    randomizeConfiguration(system, random);
    prepareForMonteCarlo(system);

    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t trial = 0; trial < 40; ++trial)
    {
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      const Molecule& molecule = system.moleculeData[m];
      std::span<const Atom> oldAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      std::vector<Atom> trialAtoms(oldAtoms.begin(), oldAtoms.end());

      // one bead (bead displacement / flip) or two adjacent beads (crankshaft)
      std::vector<std::size_t> moved;
      if (trial % 2 == 0)
        moved = {static_cast<std::size_t>(random.uniform() * 4.0)};
      else
        moved = {1, 2};
      for (std::size_t index : moved)
      {
        trialAtoms[index].position +=
            0.6 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
      }
      EXPECT_EQ(Interactions::movedAtomIndices(trialAtoms, oldAtoms), moved);

      const RunningEnergy full = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
          system.simulationBox, system.components, trialAtoms, oldAtoms);
      const RunningEnergy partial = Interactions::energyDifferenceEwaldFourierMovedAtoms(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
          system.simulationBox, system.components, trialAtoms, oldAtoms, moved);

      const double scale = std::max(1.0, std::abs(full.potentialEnergy()));
      EXPECT_NEAR(partial.ewald_fourier, full.ewald_fourier, 1e-9 * scale);
      EXPECT_NEAR(partial.ewald_exclusion, full.ewald_exclusion, 1e-9 * scale);
      EXPECT_NEAR(partial.ewald_self, full.ewald_self, 1e-9 * scale);
      EXPECT_NEAR(partial.potentialEnergy(), full.potentialEnergy(), 1e-9 * scale);
      if (chargeMethod == ForceField::ChargeMethod::Ewald)
      {
        // the moved beads change the structure factor, a non-trivial difference is being compared
        EXPECT_NE(full.ewald_fourier, 0.0);
      }
    }
  }
}

namespace
{
void expectSameStatus(EnergyStatus a, EnergyStatus b, std::size_t numberOfComponents, double tolerance)
{
  for (std::size_t compA = 0; compA < numberOfComponents; ++compA)
  {
    for (std::size_t compB = 0; compB < numberOfComponents; ++compB)
    {
      const EnergyInter& x = a.componentEnergy(compA, compB);
      const EnergyInter& y = b.componentEnergy(compA, compB);
      EXPECT_NEAR(x.VanDerWaals.energy, y.VanDerWaals.energy, tolerance) << compA << " " << compB;
      EXPECT_NEAR(x.CoulombicReal.energy, y.CoulombicReal.energy, tolerance) << compA << " " << compB;
      EXPECT_NEAR(x.VanDerWaalsTailCorrection.energy, y.VanDerWaalsTailCorrection.energy, tolerance)
          << compA << " " << compB;
    }
  }
}

double maxAbsDifference(const double3x3& a, const double3x3& b)
{
  return std::max({std::abs(a.ax - b.ax), std::abs(a.ay - b.ay), std::abs(a.az - b.az), std::abs(a.bx - b.bx),
                   std::abs(a.by - b.by), std::abs(a.bz - b.bz), std::abs(a.cx - b.cx), std::abs(a.cy - b.cy),
                   std::abs(a.cz - b.cz)});
}
}  // namespace

TEST(mc_cell_list, full_energy_and_strain_derivative_match_brute_force)
{
  // two components (rigid water and a flexible, charged chain), so that the per-component booking of the
  // tail correction and of the pair energies is checked as well
  ForceField forceField =
      ForceField({PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true),
                  PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true),
                  {"CH3", false, 15.03452, 0.0, 0.0, 6, false},
                  {"CH2", false, 14.02658, 0.0, 0.0, 6, false}},
                 {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0), {98.0, 3.75}, {46.0, 3.95}},
                 ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, /*tailCorrections=*/true, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.32;
  forceField.numberOfWaveVectors = int3(10, 10, 10);
  forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
  forceField.reciprocalIntegerCutOffSquared = 100;

  Component water = makeWater(forceField);
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
  Component chain = Component(forceField, "butane", 425.0, 3796000.0, 0.199,
                              {Atom({-1.85, -0.7, -0.15}, 0.25, 1.0, 0, 2, 1, false, false),
                               Atom({-0.31, -0.7, -0.15}, -0.25, 1.0, 0, 3, 1, false, false),
                               Atom({0.32, 0.71, -0.15}, -0.25, 1.0, 0, 3, 1, false, false),
                               Atom({1.86, 0.71, 0.15}, 0.25, 1.0, 0, 2, 1, false, false)},
                              connectivityTable, potentials, 5, 21);

  System system = System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e5, 1.0, {}, {water, chain}, {},
                         {300, 150}, 5);
  RandomNumber random(17);
  randomizeConfiguration(system, random);

  const MCCellList& list = system.cellList();
  ASSERT_TRUE(list.enabled);
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();

  // energy
  const RunningEnergy bruteForce = Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, atoms);
  const RunningEnergy withList =
      Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, list, atoms);
  const double scale = std::abs(bruteForce.moleculeMoleculeVDW) + std::abs(bruteForce.moleculeMoleculeCharge);
  EXPECT_NE(bruteForce.moleculeMoleculeVDW, 0.0);
  EXPECT_NE(bruteForce.moleculeMoleculeCharge, 0.0);
  EXPECT_NEAR(withList.moleculeMoleculeVDW, bruteForce.moleculeMoleculeVDW, 1e-10 * scale);
  EXPECT_NEAR(withList.moleculeMoleculeCharge, bruteForce.moleculeMoleculeCharge, 1e-10 * scale);

  // strain derivative, gradients and per-component energy status
  std::vector<AtomDynamics> dynamicsBrute(atoms.size()), dynamicsList(atoms.size());
  auto [statusBrute, tensorBrute] = Interactions::computeInterMolecularEnergyStrainDerivative(
      forceField, system.components, system.simulationBox, atoms, dynamicsBrute);
  auto [statusList, tensorList] = Interactions::computeInterMolecularEnergyStrainDerivative(
      forceField, system.components, system.simulationBox, list, atoms, dynamicsList);
  expectSameStatus(statusBrute, statusList, 2, 1e-10 * scale);
  EXPECT_NE(statusBrute.componentEnergy(0, 1).VanDerWaalsTailCorrection.energy, 0.0);
  EXPECT_LT(maxAbsDifference(tensorBrute, tensorList), 1e-9 * scale);
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    EXPECT_LT((dynamicsBrute[i].gradient - dynamicsList[i].gradient).length(),
              1e-9 * std::max(1.0, dynamicsBrute[i].gradient.length()));
  }

  // the fused polarization gather (field and its strain response per atom)
  std::vector<double3> fieldBrute(atoms.size()), fieldList(atoms.size());
  std::vector<std::array<double3, 9>> fieldStrainBrute(atoms.size()), fieldStrainList(atoms.size());
  std::vector<double3> comOffset(atoms.size(), double3(0.1, -0.2, 0.3));
  std::vector<double> polarizability(atoms.size(), 1.0);
  Interactions::PolarizationFieldStrain gatherBrute{fieldBrute, fieldStrainBrute, comOffset, polarizability};
  Interactions::PolarizationFieldStrain gatherList{fieldList, fieldStrainList, comOffset, polarizability};
  std::fill(dynamicsBrute.begin(), dynamicsBrute.end(), AtomDynamics{});
  std::fill(dynamicsList.begin(), dynamicsList.end(), AtomDynamics{});
  auto [statusBruteP, tensorBruteP] = Interactions::computeInterMolecularEnergyStrainDerivative(
      forceField, system.components, system.simulationBox, atoms, dynamicsBrute, &gatherBrute);
  auto [statusListP, tensorListP] = Interactions::computeInterMolecularEnergyStrainDerivative(
      forceField, system.components, system.simulationBox, list, atoms, dynamicsList, &gatherList);
  expectSameStatus(statusBruteP, statusListP, 2, 1e-10 * scale);
  EXPECT_LT(maxAbsDifference(tensorBruteP, tensorListP), 1e-9 * scale);
  double fieldScale = 0.0;
  for (const double3& f : fieldBrute) fieldScale = std::max(fieldScale, f.length());
  EXPECT_GT(fieldScale, 0.0);
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    EXPECT_LT((fieldBrute[i] - fieldList[i]).length(), 1e-9 * fieldScale);
    for (std::size_t k = 0; k < 9; ++k)
    {
      EXPECT_LT((fieldStrainBrute[i][k] - fieldStrainList[i][k]).length(), 1e-9 * fieldScale);
    }
  }
}

namespace
{
void expectSameTotal(const RunningEnergy& a, const RunningEnergy& b)
{
  const double scale = std::abs(a.moleculeMoleculeVDW) + std::abs(a.moleculeMoleculeCharge);
  EXPECT_NE(a.moleculeMoleculeVDW, 0.0);
  EXPECT_NE(a.moleculeMoleculeCharge, 0.0);
  EXPECT_NEAR(b.moleculeMoleculeVDW, a.moleculeMoleculeVDW, 1e-10 * scale);
  EXPECT_NEAR(b.moleculeMoleculeCharge, a.moleculeMoleculeCharge, 1e-10 * scale);
  for (std::size_t g = 0; g < maximumNumberOfDUDlambdaGroups; ++g)
  {
    EXPECT_NEAR(b.dudlambdaVDW[g], a.dudlambdaVDW[g], 1e-10 * std::max(1.0, std::abs(a.dudlambdaVDW[g])));
    EXPECT_NEAR(b.dudlambdaCharge[g], a.dudlambdaCharge[g], 1e-10 * std::max(1.0, std::abs(a.dudlambdaCharge[g])));
  }
}
}  // namespace

TEST(mc_cell_list, total_energy_soa_paths_match_the_scalar_pair_loop)
{
  ForceField forceField = makeWaterForceField();
  // the real-space table at the alpha of the test: the lane kernel reads it
  forceField.updateEwaldRealSpaceTable();
  Component water = makeWater(forceField);
  RandomNumber random(23);

  // a box with a cell list (the sum out of the cells) and the same configuration without one (the all-pairs sum
  // over a structure-of-arrays copy, minimum image per pair), both against the scalar triangular loop
  {
    System system =
        System(forceField, SimulationBox(36.0, 36.0, 36.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {500}, 5);
    randomizeConfiguration(system, random);
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    const MCCellList& list = system.cellList();
    ASSERT_TRUE(list.enabled);

    const RunningEnergy reference =
        Interactions::computeInterMolecularEnergyScalar(forceField, system.simulationBox, atoms);
    expectSameTotal(reference,
                    Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, list, atoms));
    expectSameTotal(reference,
                    Interactions::computeAllPairsInterMolecularEnergy(forceField, system.simulationBox, atoms));
    expectSameTotal(reference, Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, atoms));

    // a fractional molecule: the generic pair arithmetic (soft-core scaling, dU/dlambda) on both paths
    std::span<Atom> mutableAtoms = system.spanOfMoleculeAtoms();
    for (std::size_t k = 0; k < 3; ++k)
    {
      mutableAtoms[3 * 7 + k].scalingVDW = 0.6;
      mutableAtoms[3 * 7 + k].scalingCoulomb = 0.2;
      mutableAtoms[3 * 7 + k].groupId = 1;
    }
    const MCCellList& fractionalList = system.rebuiltCellList();
    ASSERT_TRUE(fractionalList.enabled);
    EXPECT_GT(fractionalList.numberOfPartiallyCoupledAtoms, 0uz);
    const RunningEnergy fractionalReference =
        Interactions::computeInterMolecularEnergyScalar(forceField, system.simulationBox, atoms);
    EXPECT_NE(fractionalReference.dudlambdaVDW[0], 0.0);
    EXPECT_NE(fractionalReference.dudlambdaCharge[0], 0.0);
    expectSameTotal(fractionalReference,
                    Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, fractionalList, atoms));
    expectSameTotal(fractionalReference,
                    Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, atoms));
  }

  // a box too small for the cell list at this cut-off: the all-pairs sum behind the cell-list overload
  {
    System system =
        System(forceField, SimulationBox(22.0, 22.0, 22.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {300}, 5);
    randomizeConfiguration(system, random);
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    const MCCellList& list = system.cellList();
    ASSERT_FALSE(list.enabled);
    const RunningEnergy reference =
        Interactions::computeInterMolecularEnergyScalar(forceField, system.simulationBox, atoms);
    expectSameTotal(reference,
                    Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, list, atoms));
    expectSameTotal(reference, Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, atoms));
  }

  // a triclinic box: the general minimum image
  {
    System system = System(forceField,
                           SimulationBox(30.0, 32.0, 34.0, 80.0 * std::numbers::pi / 180.0,
                                         95.0 * std::numbers::pi / 180.0, 100.0 * std::numbers::pi / 180.0),
                           false, 300.0, 1e5, 1.0, {}, {water}, {}, {350}, 5);
    randomizeConfiguration(system, random);
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    const RunningEnergy reference =
        Interactions::computeInterMolecularEnergyScalar(forceField, system.simulationBox, atoms);
    const MCCellList& list = system.cellList();
    expectSameTotal(reference,
                    Interactions::computeInterMolecularEnergy(forceField, system.simulationBox, list, atoms));
    expectSameTotal(reference,
                    Interactions::computeAllPairsInterMolecularEnergy(forceField, system.simulationBox, atoms));
  }
}

TEST(mc_cell_list, total_energies_and_pressure_consistent_after_unmaintained_moves)
{
  ForceField forceField = makeWaterForceField();
  Component water = makeWater(forceField);
  System system =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {500}, 5);
  RandomNumber random(29);
  randomizeConfiguration(system, random);
  prepareForMonteCarlo(system);
  const RunningEnergy reference = system.runningEnergies;
  const double3x3 referencePressure = system.computeMolecularPressure().second;

  // move every molecule behind the back of the cell list (as an integrator would) and move them back: the
  // full recomputations rebuild the list and must not pick up the stale buckets
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<Atom> saved(atoms.begin(), atoms.end());
  for (Atom& atom : atoms) atom.position += double3(13.0, -7.0, 21.0);
  const RunningEnergy shifted = system.computeTotalEnergies();
  // a rigid translation of everything leaves the energy unchanged
  EXPECT_NEAR(shifted.potentialEnergy(), reference.potentialEnergy(), 1e-8 * std::abs(reference.potentialEnergy()));
  for (std::size_t i = 0; i < atoms.size(); ++i) atoms[i].position = saved[i].position;
  const RunningEnergy restored = system.computeTotalEnergies();
  EXPECT_NEAR(restored.moleculeMoleculeVDW, reference.moleculeMoleculeVDW,
              1e-10 * std::abs(reference.moleculeMoleculeVDW));
  EXPECT_NEAR(restored.moleculeMoleculeCharge, reference.moleculeMoleculeCharge,
              1e-10 * std::abs(reference.moleculeMoleculeCharge));
  EXPECT_NEAR(restored.tail, reference.tail, 1e-10 * std::abs(reference.tail));
  EXPECT_LT(maxAbsDifference(system.computeMolecularPressure().second, referencePressure),
            1e-9 * std::abs(referencePressure.trace()));
  EXPECT_TRUE(system.verifyCellList());
}

namespace
{
// Runs 'steps' production moves on the system and checks, at the end, that the incrementally maintained cell list
// still describes the atom positions and that the running energies match a full recomputation.
void runMovesAndCheck(System& system, RandomNumber& random, std::size_t steps)
{
  prepareForMonteCarlo(system);
  const double scale = std::max(1.0, std::abs(system.runningEnergies.potentialEnergy()));

  std::size_t fractionalMoleculeSystem = 0;
  std::size_t accepted = 0;
  const std::size_t buildsBefore = system.cellList().numberOfBuilds;
  for (std::size_t step = 0; step < steps; ++step)
  {
    const std::size_t selectedComponent = system.randomComponent(random);
    const RunningEnergy before = system.runningEnergies;
    MC_Moves::performRandomMoveProduction(random, system, system, selectedComponent, fractionalMoleculeSystem, 0);
    if (before.potentialEnergy() != system.runningEnergies.potentialEnergy()) ++accepted;
    if (step % 97 == 0) ASSERT_TRUE(system.verifyCellList()) << "cell list inconsistent after step " << step;
  }
  EXPECT_GT(accepted, 0uz);
  EXPECT_TRUE(system.verifyCellList());
  // the list was maintained incrementally, not rebuilt after every move
  EXPECT_LT(system.cellList().numberOfBuilds - buildsBefore, steps / 10 + 1);

  const RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1e-6 * scale);
  EXPECT_NEAR(system.runningEnergies.moleculeMoleculeVDW, recomputed.moleculeMoleculeVDW, 1e-6 * scale);
  EXPECT_NEAR(system.runningEnergies.moleculeMoleculeCharge, recomputed.moleculeMoleculeCharge, 1e-6 * scale);
  EXPECT_NEAR(system.runningEnergies.frameworkMoleculeVDW, recomputed.frameworkMoleculeVDW, 1e-6 * scale);
  EXPECT_NEAR(system.runningEnergies.frameworkMoleculeCharge, recomputed.frameworkMoleculeCharge, 1e-6 * scale);
  EXPECT_NEAR(system.runningEnergies.ewald_fourier, recomputed.ewald_fourier, 1e-6 * scale);
}
}  // namespace

TEST(mc_cell_list, running_energies_stay_consistent_rigid_water_moves)
{
  ForceField forceField = makeWaterForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::RandomTranslation, 0.5);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::RandomRotation, 0.5);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::Swap, 0.5);
  probabilities.setProbability(Move::Types::SwapCBMC, 0.5);
  Component water = makeWater(forceField, probabilities);
  System system =
      System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 400.0, 1e5, 1.0, {}, {water}, {}, {400}, 5);
  RandomNumber random(21);
  randomizeConfiguration(system, random);
  EXPECT_TRUE(system.cellList().enabled);

  runMovesAndCheck(system, random, 3000);
}

TEST(mc_cell_list, running_energies_stay_consistent_polymer_moves)
{
  ForceField forceField = makeChainForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::BeadDisplacement, 2.0);
  probabilities.setProbability(Move::Types::Crankshaft, 1.0);
  probabilities.setProbability(Move::Types::Pivot, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  probabilities.setProbability(Move::Types::PartialReinsertionCBMC, 1.0);
  Component chain = makeChain(forceField, probabilities);
  System system =
      System(forceField, SimulationBox(40.0, 39.0, 41.0), false, 400.0, 1e5, 1.0, {}, {chain}, {}, {120}, 5);
  RandomNumber random(23);
  randomizeConfiguration(system, random);
  EXPECT_TRUE(system.cellList().enabled);

  runMovesAndCheck(system, random, 3000);
}

TEST(mc_cell_list, partially_divided_in_a_framework_box_and_moves_still_consistent)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities;

  // the constructor grows the initial molecules into the pores with CBMC (a lattice would overlap the framework)
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {co2}, {}, {40}, 5);
  RandomNumber random(31);

  // 40.2 x 39.8 x 26.8 A with a 12 A cut-off and 120 atoms: 2 cells per cut-off (6 x 6 along a and b; the
  // finer 10 x 9 grid is too dilute), the short c direction (visited whole by the stencil) collapses
  EXPECT_TRUE(MCCellList::wouldBeEnabled(system.simulationBox, system.cellListCutOff()));
  EXPECT_EQ(system.cellList().subdivision, 2);
  EXPECT_EQ(system.cellList().numberOfCells, int3(6, 6, 1));
  EXPECT_TRUE(system.verifyCellList());

  // the gathers place the records of the collapsed direction with the minimum image; the moves exercise them
  EXPECT_TRUE(system.frameworkCellList().enabled);
  runMovesAndCheck(system, random, 600);
}

TEST(mc_cell_list, small_framework_box_moves_still_consistent)
{
  // ITQ-29 1x1x1 (11.9 A): neither list is enabled, the framework and the molecules are gathered from all atoms
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities = probabilities;
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {co2}, {}, {4}, 5);
  RandomNumber random(41);
  EXPECT_FALSE(system.cellList().enabled);
  EXPECT_FALSE(system.frameworkCellList().enabled);

  runMovesAndCheck(system, random, 600);
}

namespace
{
// CO2 (Ewald) in MFI 2x2x2 (the framework cell list is 10 x 9 x 1, the collapsed direction gathers with the minimum
// image) or in ITQ-29 1x1x1 (11.9 A: no cell list at all, the framework gathers run over all framework atoms and the
// kernels apply the minimum image per pair).
System makeFrameworkSystem(bool large, std::size_t numberOfMolecules)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = large ? Framework::makeMFI(forceField, int3(2, 2, 2)) : Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, true);
  return System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {co2}, {}, {numberOfMolecules}, 5);
}
}  // namespace

TEST(mc_cell_list, framework_energy_difference_matches_brute_force)
{
  for (bool large : {true, false})
  {
    System system = makeFrameworkSystem(large, large ? 40 : 4);
    RandomNumber random(53);
    const MCCellList& list = system.frameworkCellList();
    EXPECT_EQ(list.enabled, large);
    if (large)
    {
      EXPECT_EQ(list.numberOfCells, int3(10, 9, 1));
      EXPECT_EQ(list.numberOfAtoms, system.spanOfFrameworkAtoms().size());
      EXPECT_TRUE(list.verify(system.spanOfFrameworkAtoms(), system.simulationBox));
    }
    const std::size_t buildsBefore = list.numberOfBuilds;

    std::size_t overlaps = 0;
    for (std::size_t round = 0; round < 120; ++round)
    {
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      const Molecule& molecule = system.moleculeData[m];
      std::span<const Atom> oldAtoms = system.spanOfMoleculeAtoms().subspan(molecule.atomIndex, molecule.numberOfAtoms);
      // small displacements (accepted-move sized) and large ones that run into the framework now and then
      const std::vector<Atom> newAtoms = displacedMolecule(system, random, m, round % 3 == 0 ? 3.0 : 0.5);

      const std::optional<RunningEnergy> bruteForce = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
      const std::optional<RunningEnergy> withList = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.frameworkCellList(), system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
      expectSameFrameworkEnergies(bruteForce, withList);
      if (!bruteForce.has_value()) ++overlaps;

      // insertion and deletion shapes (an empty side)
      expectSameFrameworkEnergies(
          Interactions::computeFrameworkMoleculeEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.interpolationGrids, system.framework,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, {}),
          Interactions::computeFrameworkMoleculeEnergyDifference(
              system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
              system.frameworkCellList(), system.spanOfFrameworkAtoms(), newAtoms, {}));
      expectSameFrameworkEnergies(
          Interactions::computeFrameworkMoleculeEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.interpolationGrids, system.framework,
                                                                 system.spanOfFrameworkAtoms(), {}, oldAtoms),
          Interactions::computeFrameworkMoleculeEnergyDifference(
              system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
              system.frameworkCellList(), system.spanOfFrameworkAtoms(), {}, oldAtoms));
    }
    EXPECT_GT(overlaps, 0uz);
    // the framework never moves: the list is built once
    EXPECT_EQ(list.numberOfBuilds, buildsBefore);
  }
}

TEST(mc_cell_list, single_atom_cell_kernel_matches_gather_and_brute_force)
{
  using Background = Interactions::Neighbourhood::Background;

  // molecules: shifted Lennard-Jones with Ewald in a large water box (the lane path), then with a fractional
  // background molecule (the generic path)
  {
    ForceField forceField =
        ForceField({PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true),
                    PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true)},
                   {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)},
                   ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, /*shifted=*/true, false, true);
    forceField.automaticEwald = false;
    forceField.EwaldAlpha = 0.3;
    forceField.numberOfWaveVectors = int3(8, 8, 8);
    forceField.reciprocalCutOffSquared = std::numeric_limits<double>::max();
    forceField.reciprocalIntegerCutOffSquared = 64;
    forceField.updateEwaldRealSpaceTable();
    Component water = makeWater(forceField);
    System system =
        System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 300.0, 1e5, 1.0, {}, {water}, {}, {700}, 5);
    RandomNumber random(67);
    randomizeConfiguration(system, random);
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    ASSERT_TRUE(system.cellList().enabled);

    for (bool fractional : {false, true})
    {
      if (fractional)
      {
        const Molecule& molecule = system.moleculeData[5];
        for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
        {
          atoms[molecule.atomIndex + k].scalingVDW = 0.5;
          atoms[molecule.atomIndex + k].scalingCoulomb = 0.25;
          atoms[molecule.atomIndex + k].groupId = 1;
        }
        system.invalidateCellList();
      }
      const MCCellList& list = system.cellList();
      EXPECT_EQ(list.numberOfPartiallyCoupledAtoms, fractional ? 3uz : 0uz);

      std::size_t overlaps = 0;
      for (std::size_t round = 0; round < 150; ++round)
      {
        const std::size_t m =
            static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
        // an oxygen or a hydrogen as a first bead, now and then close to an existing molecule (overlaps), and
        // also outside the primary box (the kernel wraps it)
        Atom trial = atoms[system.moleculeData[m].atomIndex + (round % 3 == 0 ? 1 : 0)];
        trial.position =
            round % 4 == 0
                ? trial.position + 0.9 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5)
                : system.simulationBox.cell *
                      double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
        trial.moleculeId = static_cast<std::uint32_t>(system.moleculeData.size() + 7);
        const std::span<const Atom> set(&trial, 1);
        // every other round the molecule 'm' is the one being regrown: left out of the background
        const std::optional<std::uint32_t> skipped =
            round % 2 == 0 ? std::optional{static_cast<std::uint32_t>(m)} : std::nullopt;
        const std::optional<std::size_t> skippedSize =
            skipped.has_value() ? std::optional{static_cast<std::size_t>(skipped.value())} : std::nullopt;

        for (double cutOff : {12.0, 7.5})
        {
          const std::optional<RunningEnergy> direct =
              Interactions::computeCellListAtomEnergy(system.forceField, system.simulationBox, list, trial, cutOff,
                                                      cutOff, 1.0, true, skipped, Background::Molecules);
          Interactions::Neighbourhood neighbourhood;
          neighbourhood.gather(system.simulationBox, &list, atoms, trial.position, 0.0, cutOff, skipped);
          const std::optional<RunningEnergy> gathered = Interactions::computeNeighbourhoodEnergy(
              system.forceField, system.simulationBox, neighbourhood, set, cutOff, cutOff, 1.0, true);
          const std::optional<RunningEnergy> bruteForce = CBMC::computeInterMolecularEnergy(
              system.forceField, system.simulationBox, atoms, cutOff, cutOff, set, skippedSize);
          expectSameEnergies(bruteForce, direct);
          expectSameEnergies(gathered, direct);
          if (!direct.has_value())
          {
            ++overlaps;
            continue;
          }
          if (fractional)
          {
            ASSERT_TRUE(bruteForce.has_value());
            EXPECT_NEAR(direct->dudlambdaVDW[0], bruteForce->dudlambdaVDW[0],
                        1e-10 * std::max(1.0, std::abs(bruteForce->dudlambdaVDW[0])));
            EXPECT_NEAR(direct->dudlambdaCharge[0], bruteForce->dudlambdaCharge[0],
                        1e-10 * std::max(1.0, std::abs(bruteForce->dudlambdaCharge[0])));
          }
          // the old-configuration mode: negated, no overlap check
          const std::optional<RunningEnergy> negated =
              Interactions::computeCellListAtomEnergy(system.forceField, system.simulationBox, list, trial, cutOff,
                                                      cutOff, -1.0, false, skipped, Background::Molecules);
          ASSERT_TRUE(negated.has_value());
          EXPECT_DOUBLE_EQ(negated->moleculeMoleculeVDW, -direct->moleculeMoleculeVDW);
          EXPECT_DOUBLE_EQ(negated->moleculeMoleculeCharge, -direct->moleculeMoleculeCharge);
        }
      }
      EXPECT_GT(overlaps, 0uz);
    }
  }

  // the framework background (a list with a collapsed direction): a CO2 oxygen or carbon as a first bead
  {
    System system = makeFrameworkSystem(true, 20);
    RandomNumber random(71);
    const MCCellList& list = system.frameworkCellList();
    ASSERT_TRUE(list.enabled);
    std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    std::size_t overlaps = 0;
    for (std::size_t round = 0; round < 150; ++round)
    {
      const std::size_t m =
          static_cast<std::size_t>(random.uniform() * static_cast<double>(system.moleculeData.size()));
      Atom trial = atoms[system.moleculeData[m].atomIndex + (round % 2)];
      trial.position =
          round % 3 == 0
              ? trial.position + 2.0 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5)
              : system.simulationBox.cell *
                    double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
      const std::span<const Atom> set(&trial, 1);
      for (double cutOff : {12.0, 9.0})
      {
        const std::optional<RunningEnergy> direct =
            Interactions::computeCellListAtomEnergy(system.forceField, system.simulationBox, list, trial, cutOff,
                                                    cutOff, 1.0, true, std::nullopt, Background::Framework);
        Interactions::Neighbourhood neighbourhood;
        neighbourhood.gatherFramework(system.simulationBox, &list, frameworkAtoms, trial.position, 0.0, cutOff);
        const std::optional<RunningEnergy> gathered = Interactions::computeNeighbourhoodEnergy(
            system.forceField, system.simulationBox, neighbourhood, set, cutOff, cutOff, 1.0, true);
        expectSameFrameworkEnergies(gathered, direct);
        // against the brute force over the framework atoms (the pair sums; this system has no grids)
        const std::optional<RunningEnergy> bruteForce =
            CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                 system.framework, frameworkAtoms, cutOff, cutOff, set);
        expectSameFrameworkEnergies(bruteForce, direct);
        if (!direct.has_value()) ++overlaps;
      }
    }
    EXPECT_GT(overlaps, 0uz);
  }
}

TEST(mc_cell_list, cbmc_framework_energy_matches_brute_force_and_batches)
{
  for (bool large : {true, false})
  {
    System system = makeFrameworkSystem(large, large ? 40 : 4);
    RandomNumber random(59);
    const Component& co2 = system.components[0];
    const CBMC::GrowContext context = system.makeGrowContext();
    ASSERT_NE(context.frameworkCellList, nullptr);
    EXPECT_EQ(context.frameworkCellList->enabled, large);
    const CBMC::GrowContext reference = context.withFrameworkCellList(nullptr);

    std::size_t overlaps = 0;
    for (std::size_t round = 0; round < 30; ++round)
    {
      // a batch of trial orientations around a common first bead, as a rigid regrow produces
      const double3 anchor = system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      std::vector<std::vector<Atom>> trialSets;
      for (std::size_t k = 0; k < 8; ++k)
      {
        std::vector<Atom> set(co2.atoms);
        const double3x3 rotation = double3x3::buildRotationMatrixInverse(random.randomSimdQuatd());
        for (std::size_t a = 0; a < set.size(); ++a)
        {
          set[a].position = anchor + rotation * co2.atoms[a].position;
          set[a].moleculeId = 7;
        }
        trialSets.push_back(std::move(set));
      }
      std::vector<std::span<const Atom>> views(trialSets.begin(), trialSets.end());

      Interactions::Neighbourhood frameworkNeighbourhood;
      ASSERT_TRUE(CBMC::gatherFrameworkNeighbourhood(context, frameworkNeighbourhood, views));
      EXPECT_EQ(frameworkNeighbourhood.background, Interactions::Neighbourhood::Background::Framework);
      EXPECT_EQ(frameworkNeighbourhood.paddedSize() % Interactions::Neighbourhood::padding, 0uz);
      EXPECT_TRUE(frameworkNeighbourhood.fullyCoupled);
      for (std::uint32_t id : frameworkNeighbourhood.moleculeId)
      {
        EXPECT_TRUE(id == Interactions::Neighbourhood::frameworkMoleculeId ||
                    id == Interactions::Neighbourhood::paddingMoleculeId);
      }
      EXPECT_TRUE(Interactions::neighbourhoodFastKernelApplies(system.forceField, frameworkNeighbourhood, views[0]));
      const CBMC::GrowContext batch = context.withFrameworkNeighbourhood(&frameworkNeighbourhood);

      for (std::span<const Atom> set : views)
      {
        ASSERT_TRUE(frameworkNeighbourhood.covers(set, context.cutOffFrameworkVDW, context.cutOffCoulomb, true));
        const std::optional<RunningEnergy> bruteForce = CBMC::computeFrameworkMoleculeEnergy(reference, set);
        const std::optional<RunningEnergy> gathered = CBMC::computeFrameworkMoleculeEnergy(context, set);
        const std::optional<RunningEnergy> batched = CBMC::computeFrameworkMoleculeEnergy(batch, set);
        expectSameFrameworkEnergies(bruteForce, gathered);
        expectSameFrameworkEnergies(bruteForce, batched);
        if (!bruteForce.has_value()) ++overlaps;
        // a single atom of the set (a first bead, or the brute-force path for one atom without a usable list)
        expectSameFrameworkEnergies(CBMC::computeFrameworkMoleculeEnergy(reference, set.subspan(0, 1)),
                                    CBMC::computeFrameworkMoleculeEnergy(context, set.subspan(0, 1)));
      }

      // the first-bead batch (spread over the whole box): the complete external energies agree with the
      // brute-force framework loop, through the per-trial gathers (large box) or one copy of all framework
      // atoms shared by the batch (small box)
      std::vector<Atom> firstBeads;
      for (std::size_t k = 0; k < 10; ++k)
      {
        Atom bead = co2.atoms[co2.startingBead];
        bead.position = system.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
        bead.moleculeId = 7;
        firstBeads.push_back(bead);
      }
      const std::vector<CBMC::FirstBeadTrial> expected =
          CBMC::computeExternalNonOverlappingEnergies(reference, co2, firstBeads);
      const std::vector<CBMC::FirstBeadTrial> actual =
          CBMC::computeExternalNonOverlappingEnergies(context, co2, firstBeads);
      ASSERT_EQ(expected.size(), actual.size());
      for (std::size_t k = 0; k < expected.size(); ++k)
      {
        EXPECT_EQ(expected[k].position.position.x, actual[k].position.position.x);
        expectSameFrameworkEnergies(expected[k].energy, actual[k].energy);
      }
    }
    EXPECT_GT(overlaps, 0uz);
  }
}

TEST(mc_cell_list, disabled_in_a_small_box_and_moves_still_consistent)
{
  // the water box of the examples: 24.83 A with a 12 A cut-off, every direction collapses
  ForceField forceField = makeWaterForceField(ForceField::ChargeMethod::Ewald, 12.0);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  // the random translation displaces over up to a whole box length: the batch of its new and old positions
  // spans the box and the gather falls back to the nearest image per atom with the minimum image per pair
  probabilities.setProbability(Move::Types::RandomTranslation, 1.0);
  probabilities.setProbability(Move::Types::RandomRotation, 0.5);
  Component water = makeWater(forceField, probabilities);
  System system =
      System(forceField, SimulationBox(24.83, 24.83, 24.83), false, 400.0, 1e5, 1.0, {}, {water}, {}, {150}, 5);
  RandomNumber random(37);
  randomizeConfiguration(system, random);
  EXPECT_FALSE(system.cellList().enabled);
  EXPECT_TRUE(system.verifyCellList());

  // the all-atom gather with the per-pair minimum image is exercised by the moves
  runMovesAndCheck(system, random, 2000);
}
