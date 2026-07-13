#include <gtest/gtest.h>

import std;

import double3;
import atom;
import atom_dynamics;
import forcefield;
import component;
import system;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_intramolecular;

namespace
{
double bondEnergy(const Potentials::IntraMolecularPotentials &potentials, const double3 &posA, const double3 &posB)
{
  return std::get<0>(potentials.bonds.front().potentialEnergyGradientStrain(posA, posB));
}

double bondEnergyAtAtoms(std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &potentials)
{
  return bondEnergy(potentials, atoms[0].position, atoms[1].position);
}

void perturbAtom(std::span<Atom> atoms, std::size_t atom, std::size_t axis, double delta)
{
  (&atoms[atom].position.x)[axis] += delta;
}
}  // namespace

TEST(minimization_hessian_bond, harmonic_ethane_flexible_analytic_matches_finite_difference)
{
  const double delta = 1e-5;
  const double tolerance = 1e-2;

  ForceField forceField = ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 8, false},
                                      {"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                                      {"CH2", false, 14.02658, 0.0, 0.0, 8, false},
                                      {"CH", false, 13.01864, 0.0, 0.0, 8, false},
                                      {"C", false, 12.0, 0.0, 0.0, 8, false}},
                                     {{158.5, 3.72}, {108.0, 3.76}, {56.0, 3.96}, {17.0, 4.67}, {0.8, 6.38}},
                                     ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);

  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component component = Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                                  {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false),
                                   Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1},
                         5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularBondHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                 dynamics);

  struct DofLabel
  {
    std::size_t atom;
    std::size_t axis;
  };

  std::array<DofLabel, 6> labels{};
  for (std::size_t dof = 0; dof < 6; ++dof)
  {
    bool found = false;
    for (std::size_t atom = 0; atom < 2 && !found; ++atom)
    {
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        if (layout.flexibleAtomDof(0, atom, static_cast<MinimizationDofAxis>(axis)) == dof)
        {
          labels[dof] = {atom, axis};
          found = true;
          break;
        }
      }
    }
    ASSERT_TRUE(found);
  }

  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t col = 0; col < 6; ++col)
    {
      const DofLabel rowLabel = labels[row];
      const DofLabel colLabel = labels[col];
      double numerical{};

      if (row == col)
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        const double ePlus = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMinus = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        numerical = (ePlus - 2.0 * bondEnergyAtAtoms(atoms, intraMolecularPotentials) + eMinus) / (delta * delta);
      }
      else
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, delta);
        const double ePP = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -2.0 * delta);
        const double ePM = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMM = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, 2.0 * delta);
        const double eMP = bondEnergyAtAtoms(atoms, intraMolecularPotentials);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -delta);
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }

      EXPECT_NEAR(numerical, hessian(row, col), tolerance) << "row=" << row << " col=" << col;
    }
  }
}
