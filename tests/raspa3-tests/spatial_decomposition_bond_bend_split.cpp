#include <gtest/gtest.h>

import std;

import double3;
import atom;
import atom_dynamics;
import forcefield;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bond_bend_potential;
import simulationbox;
import system;
import running_energy;
import integrators_update;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;

TEST(SPATIAL_DECOMPOSITION_BOND_BEND_SPLIT, split_molecule_bond_bend_gradients_match_the_whole_molecule)
{
  // A 343-bead chain (atoms plus terms well above the split threshold) is sliced into term ranges per kind. A
  // bond-bend term has three atoms and three gradients although its identifier array has four entries: the slot
  // layout must use three per term, or the gradients are read misaligned and a fourth slot lands on atom 0.
  ForceField forceField({{"CH2", false, 14.02658, 0.0, 0.0, 6, false}}, {{46.0, 3.95}},
                        ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, false);

  const std::size_t n = 7;
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
    potentials.bondBends.push_back(
        BondBendPotential({i, i + 1, i + 2}, BondBendType::CVFF, {100.0, 2000.0, 3.8, 1500.0, 3.7}));
  }

  // a snake through a 7x7x7 lattice, perturbed so that the bond lengths and angles are off their reference values
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
        const double i = static_cast<double>(beads.size());
        beads.push_back(Atom({spacing * static_cast<double>(x) - offset + 0.2 * std::sin(1.3 * i),
                              spacing * static_cast<double>(y) - offset + 0.2 * std::cos(0.7 * i),
                              spacing * static_cast<double>(z) - offset + 0.2 * std::sin(2.1 * i)},
                             0.0, 1.0, 0, 0, 0, false, false));
      }
    }
  }
  Component chain(forceField, "chain", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  std::vector<double3> positions{};
  for (const Atom& bead : beads) positions.push_back(bead.position + double3(20.0, 19.5, 20.5));
  System system(forceField, SimulationBox(40.0, 39.0, 41.0), false, 300.0, 1e5, 1.0, {}, {chain}, {positions}, {0}, 5);

  const RunningEnergy reference = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
  std::vector<double3> referenceGradient{};
  for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) referenceGradient.push_back(dynamics.gradient);
  ASSERT_NE(reference.bondBend, 0.0);

  for (std::size_t threads : {1uz, 4uz})
  {
    SpatialDecompositionSettings settings;
    settings.numberOfThreads = threads;
    settings.verletSkin = 1.5;
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);

    EXPECT_NEAR(energy.bondBend, reference.bondBend, 1e-10 * std::abs(reference.bondBend)) << "threads " << threads;
    std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
    for (std::size_t i = 0; i < dynamics.size(); ++i)
    {
      const double3 difference = dynamics[i].gradient - referenceGradient[i];
      EXPECT_LT(std::sqrt(double3::dot(difference, difference)), 1e-6 * (1.0 + referenceGradient[i].length()))
          << "threads " << threads << ", atom " << i;
    }
  }
}
