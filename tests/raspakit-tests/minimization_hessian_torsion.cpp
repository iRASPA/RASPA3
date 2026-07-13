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
import torsion_potential;
import urey_bradley_potential;
import van_der_waals_potential;
import coulomb_potential;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_intramolecular;

namespace
{
ForceField makeTestForceField()
{
  return ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 8, false},
                     {"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                     {"CH2", false, 14.02658, 0.0, 0.0, 8, false},
                     {"CH", false, 13.01864, 0.0, 0.0, 8, false},
                     {"C", false, 12.0, 0.0, 0.0, 8, false}},
                    {{158.5, 3.72}, {108.0, 3.76}, {56.0, 3.96}, {17.0, 4.67}, {0.8, 6.38}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
}

void perturbAtom(std::span<Atom> atoms, std::size_t atom, std::size_t axis, double delta)
{
  (&atoms[atom].position.x)[axis] += delta;
}

struct DofLabel
{
  std::size_t atom;
  std::size_t axis;
};

std::vector<DofLabel> buildDofLabels(const MinimizationDofLayout &layout, std::size_t numAtoms)
{
  std::vector<DofLabel> labels(layout.numDofs());
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    bool found = false;
    for (std::size_t atom = 0; atom < numAtoms && !found; ++atom)
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
    EXPECT_TRUE(found);
  }
  return labels;
}

void compareHessianToFiniteDifference(std::span<Atom> atoms, const std::vector<DofLabel> &labels,
                                      const GeneralizedHessian &hessian,
                                      const std::function<double()> &energyFunction, double delta,
                                      double relativeTolerance)
{
  const std::size_t numDofs = labels.size();
  for (std::size_t row = 0; row < numDofs; ++row)
  {
    for (std::size_t col = 0; col < numDofs; ++col)
    {
      const DofLabel rowLabel = labels[row];
      const DofLabel colLabel = labels[col];
      double numerical{};

      if (row == col)
      {
        const double eZero = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        const double ePlus = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMinus = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        numerical = (ePlus - 2.0 * eZero + eMinus) / (delta * delta);
      }
      else
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, delta);
        const double ePP = energyFunction();
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -2.0 * delta);
        const double ePM = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMM = energyFunction();
        perturbAtom(atoms, colLabel.atom, colLabel.axis, 2.0 * delta);
        const double eMP = energyFunction();
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

System makeButaneSystem(const ForceField &forceField, TorsionType torsionType, std::vector<double> parameters)
{
  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.torsions = {TorsionPotential({0, 1, 2, 3}, torsionType, std::move(parameters))};

  Component component = Component(forceField, "butane", 425.125, 3796000.0, 0.201,
                                  {Atom({-0.5, 1.4, 0.3}, 0.0, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({2.1, 1.2, 1.0}, 0.0, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  return System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1}, 5);
}

void runTorsionFiniteDifferenceTest(TorsionType torsionType, std::vector<double> parameters)
{
  const double delta = 1e-5;
  const double relativeTolerance = 5e-3;

  ForceField forceField = makeTestForceField();
  System system = makeButaneSystem(forceField, torsionType, std::move(parameters));

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularTorsionHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                    dynamics);

  const Potentials::IntraMolecularPotentials &potentials =
      system.components[0].intraMolecularPotentials;
  auto energyFunction = [&]() -> double
  {
    double energy = 0.0;
    for (const TorsionPotential &torsion : potentials.torsions)
    {
      energy += torsion.calculateEnergy(atoms[torsion.identifiers[0]].position, atoms[torsion.identifiers[1]].position,
                                        atoms[torsion.identifiers[2]].position, atoms[torsion.identifiers[3]].position);
    }
    return energy;
  };

  const std::vector<DofLabel> labels = buildDofLabels(layout, 4);
  compareHessianToFiniteDifference(atoms, labels, hessian, energyFunction, delta, relativeTolerance);
}
}  // namespace

TEST(minimization_hessian_torsion, trappe_butane_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32});
}

TEST(minimization_hessian_torsion, harmonic_dihedral_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::Harmonic, {5000.0, 65.0});
}

TEST(minimization_hessian_torsion, harmonic_cosine_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::HarmonicCosine, {800.0, 60.0});
}

TEST(minimization_hessian_intramolecular_pairs, urey_bradley_vdw_coulomb_match_finite_difference)
{
  const double delta = 1e-5;
  const double relativeTolerance = 5e-3;

  ForceField forceField = makeTestForceField();

  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.ureyBradleys = {UreyBradleyPotential({0, 2}, UreyBradleyType::Harmonic, {40000.0, 2.6})};
  intraMolecularPotentials.vanDerWaals = {
      VanDerWaalsPotential({0, 2}, VanDerWaalsType::LennardJones, {120.0, 2.5}, 1.0)};
  intraMolecularPotentials.coulombs = {CoulombPotential({0, 2}, CoulombType::Coulomb, 0.4, -0.35, 1.0)};

  Component component = Component(forceField, "triatomic", 369.825, 4247660.0, 0.1524,
                                  {Atom({-0.3, 1.45, 0.2}, 0.4, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({1.5, -0.2, 0.4}, -0.35, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1},
                         5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 9u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularUreyBradleyHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                        dynamics);
  Interactions::computeIntraMolecularVanDerWaalsHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                        dynamics);
  Interactions::computeIntraMolecularCoulombHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                    dynamics);

  const Potentials::IntraMolecularPotentials &potentials = system.components[0].intraMolecularPotentials;
  auto energyFunction = [&]() -> double
  {
    double energy = 0.0;
    for (const UreyBradleyPotential &ureyBradley : potentials.ureyBradleys)
    {
      energy += ureyBradley.calculateEnergy(atoms[ureyBradley.identifiers[0]].position,
                                            atoms[ureyBradley.identifiers[1]].position);
    }
    for (const VanDerWaalsPotential &vanDerWaals : potentials.vanDerWaals)
    {
      energy += vanDerWaals.calculateEnergy(atoms[vanDerWaals.identifiers[0]].position,
                                            atoms[vanDerWaals.identifiers[1]].position);
    }
    for (const CoulombPotential &coulomb : potentials.coulombs)
    {
      energy += coulomb.calculateEnergy(atoms[coulomb.identifiers[0]].position, atoms[coulomb.identifiers[1]].position);
    }
    return energy;
  };

  const std::vector<DofLabel> labels = buildDofLabels(layout, 3);
  compareHessianToFiniteDifference(atoms, labels, hessian, energyFunction, delta, relativeTolerance);
}
