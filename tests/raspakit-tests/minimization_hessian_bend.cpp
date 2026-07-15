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
import bend_potential;
import bond_potential;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_intramolecular;

namespace
{
double bendEnergy(const Potentials::IntraMolecularPotentials &potentials, std::span<const Atom> atoms)
{
  double energy = 0.0;
  for (const BendPotential &bend : potentials.bends)
  {
    const std::size_t A = bend.identifiers[0];
    const std::size_t B = bend.identifiers[1];
    const std::size_t C = bend.identifiers[2];
    energy += bend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, std::nullopt);
  }
  return energy;
}

void perturbAtom(std::span<Atom> atoms, std::size_t atom, std::size_t axis, double delta)
{
  (&atoms[atom].position.x)[axis] += delta;
}
}  // namespace

TEST(minimization_hessian_bend, harmonic_propane_flexible_analytic_matches_finite_difference)
{
  // delta ~ eps^(1/4) balances truncation against roundoff for a second-derivative central
  // difference; 1e-5 puts the roundoff noise (~ eps * k / delta^2) above the tolerance.
  const double delta = 1e-4;
  // The finite-difference reference carries O(delta^2) truncation and roundoff noise;
  // the analytic entries are exact, so the tolerance reflects the FD error, not the Hessian.
  const double relativeTolerance = 5e-3;

  ForceField forceField = ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 8, false},
                                      {"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                                      {"CH2", false, 14.02658, 0.0, 0.0, 8, false},
                                      {"CH", false, 13.01864, 0.0, 0.0, 8, false},
                                      {"C", false, 12.0, 0.0, 0.0, 8, false}},
                                     {{158.5, 3.72}, {108.0, 3.76}, {56.0, 3.96}, {17.0, 4.67}, {0.8, 6.38}},
                                     ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);

  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                                    BondPotential({1, 2}, BondType::Harmonic, {96500.0, 1.54})};
  intraMolecularPotentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0})};

  Component component = Component(forceField, "propane", 369.825, 4247660.0, 0.1524,
                                  {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false),
                                   Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({3.08, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  // Deterministic geometry a few degrees away from the equilibrium bend angle. Growing the
  // molecule instead (initialNumberOfMolecules = 1) uses a time-seeded RNG, which makes the
  // stored energy -- and hence the finite-difference noise floor -- vary from run to run.
  const std::vector<std::vector<double3>> initialPositions = {
      {double3(0.0, 0.0, 0.0), double3(1.54, 0.0, 0.0), double3(2.10, 1.47, 0.10)}};

  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component},
                         initialPositions, {0}, 5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 9u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularBendHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                 dynamics);

  struct DofLabel
  {
    std::size_t atom;
    std::size_t axis;
  };

  std::array<DofLabel, 9> labels{};
  for (std::size_t dof = 0; dof < 9; ++dof)
  {
    bool found = false;
    for (std::size_t atom = 0; atom < 3 && !found; ++atom)
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

  for (std::size_t row = 0; row < 9; ++row)
  {
    for (std::size_t col = 0; col < 9; ++col)
    {
      const DofLabel rowLabel = labels[row];
      const DofLabel colLabel = labels[col];
      double numerical{};

      if (row == col)
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        const double ePlus = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMinus = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        numerical = (ePlus - 2.0 * bendEnergy(intraMolecularPotentials, atoms) + eMinus) / (delta * delta);
      }
      else
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, delta);
        const double ePP = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -2.0 * delta);
        const double ePM = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMM = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, 2.0 * delta);
        const double eMP = bendEnergy(intraMolecularPotentials, atoms);
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -delta);
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }

      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }
}
