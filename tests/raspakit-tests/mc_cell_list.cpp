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
import interactions_pair_kernel;
import interactions_intermolecular;
import interactions_ewald;

// The Monte Carlo cell list replaces the O(N) loops of the single-molecule energy differences and of the CBMC
// trial energies by a visit of the 27 cells around every trial position. These tests check that
//   - the buckets and back-references are consistent after builds and incremental updates,
//   - the cell-list energy routines reproduce the brute-force routines,
//   - the list stays consistent with the atom positions over many accepted MC moves of all kinds and the running
//     energies stay in step with the recomputed total energies,
//   - the Ewald difference restricted to the moved atoms equals the full-molecule difference,
//   - the list is disabled for a box with fewer than four cells in every direction.

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
}  // namespace

TEST(mc_cell_list, grid_selection_and_disabled_for_small_boxes)
{
  // 40 A box, 9 A cut-off: floor(40/9) = 4 cells per direction
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(40.0, 40.0, 40.0), 9.0), int3(4, 4, 4));
  EXPECT_TRUE(MCCellList::wouldBeEnabled(SimulationBox(40.0, 40.0, 40.0), 9.0));

  // 75.75 A box, 14 A cut-off (the pHDDA-10 case): 5x5x5
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(75.75, 75.75, 75.75), 14.0), int3(5, 5, 5));

  // 3 cells of 12 A fit in 36 A: no reduction, every direction collapses, list disabled
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(36.0, 36.0, 36.0), 12.0), int3(1, 1, 1));
  EXPECT_FALSE(MCCellList::wouldBeEnabled(SimulationBox(36.0, 36.0, 36.0), 12.0));

  // an elongated box: only the long direction is divided
  EXPECT_EQ(MCCellList::gridFor(SimulationBox(30.0, 30.0, 100.0), 9.0), int3(1, 1, 11));
  EXPECT_TRUE(MCCellList::wouldBeEnabled(SimulationBox(30.0, 30.0, 100.0), 9.0));

  // a triclinic box uses the perpendicular widths, not the vector lengths
  const SimulationBox triclinic(40.0, 40.0, 40.0, 60.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0,
                                90.0 * std::numbers::pi / 180.0);
  const double3 widths = triclinic.perpendicularWidths();
  const int3 grid = MCCellList::gridFor(triclinic, 9.0);
  for (std::size_t d = 0; d < 3; ++d)
  {
    const int expected = static_cast<int>(std::floor(widths[d] / 9.0));
    EXPECT_EQ(grid[d], expected >= 4 ? expected : 1);
  }
}

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
  EXPECT_EQ(list.numberOfCells, int3(4, 4, 4));
  EXPECT_EQ(list.numberOfAtoms, atoms.size());
  EXPECT_TRUE(list.verify(atoms, system.simulationBox));
  EXPECT_TRUE(list.isCurrent(atoms, system.simulationBox, 9.0));
  EXPECT_FALSE(list.isCurrent(atoms, system.simulationBox, 10.0));
  EXPECT_FALSE(list.isCurrent(atoms, SimulationBox(41.0, 40.0, 40.0), 9.0));

  // a query visits 27 of the 64 cells
  EXPECT_NEAR(list.averageNeighbourhoodSize(), 27.0 / 64.0 * static_cast<double>(atoms.size()),
              0.05 * static_cast<double>(atoms.size()));

  // the 27-cell neighbourhood of a position contains every atom within the cut-off of that position
  for (std::size_t trial = 0; trial < 50; ++trial)
  {
    const double3 position =
        system.simulationBox.cell *
        double3(random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0, random.uniform() * 3.0 - 1.0);
    std::set<std::uint32_t> visited;
    list.forEachNeighbourRecord(position, [&](const MCCellList::Record& record) { visited.insert(record.atomIndex); });
    for (std::uint32_t i = 0; i < atoms.size(); ++i)
    {
      const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(atoms[i].position - position);
      if (double3::dot(dr, dr) < 81.0) EXPECT_TRUE(visited.contains(i)) << "atom " << i << " missed";
    }
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

TEST(mc_cell_list, disabled_in_a_small_framework_box_and_moves_still_consistent)
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

  // 40 x 40 x 27 A with a 12 A cut-off: 3 x 3 x 2 cells, no direction reaches four cells
  EXPECT_FALSE(MCCellList::wouldBeEnabled(system.simulationBox, system.cellListCutOff()));
  EXPECT_FALSE(system.cellList().enabled);
  EXPECT_TRUE(system.verifyCellList());

  // the brute-force fallback of the cell-list routines is exercised by the moves
  runMovesAndCheck(system, random, 600);
}
