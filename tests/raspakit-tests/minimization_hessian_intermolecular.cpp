#include <gtest/gtest.h>

import std;

import double3;
import atom;
import atom_dynamics;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_intermolecular;
import interactions_intermolecular;

namespace
{
void useSecondOrderTaylorShiftedLennardJones(ForceField &forceField)
{
  for (VDWParameters &parameters : forceField.data)
  {
    if (parameters.type == VDWParameters::Type::LennardJones)
    {
      parameters.type = VDWParameters::Type::LennardJonesSecondOrderTaylorShifted;
    }
  }
  forceField.preComputeDerivedParameters();
  forceField.preComputePotentialShift();
}
}  // namespace

TEST(minimization_hessian_intermolecular, flexible_methane_pair_vdw_matches_finite_difference)
{
  const double delta = 1e-5;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;

  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(8.0, 10.0, 10.0);
  atoms[1].position = double3(9.4, 10.0, 10.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  Interactions::computeInterMolecularHessian(system.forceField, system.simulationBox, system.moleculeData,
                                             system.components, system.spanOfMoleculeAtoms(), layout, hessian,
                                             dynamics);

  auto interEnergy = [&]()
  {
    return Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atoms).moleculeMoleculeVDW;
  };

  for (std::size_t atom = 0; atom < atoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double &coordinate = (&atoms[atom].position.x)[axis];
      coordinate += delta;
      const double plus = interEnergy();
      coordinate -= 2.0 * delta;
      const double minus = interEnergy();
      coordinate += delta;
      const double numericalGradient = (plus - minus) / (2.0 * delta);
      EXPECT_NEAR((&dynamics[atom].gradient.x)[axis], numericalGradient,
                  relativeTolerance * std::max(1.0, std::abs(numericalGradient)));
    }
  }

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    for (std::size_t col = 0; col < layout.numDofs(); ++col)
    {
      const std::size_t moleculeA = row / 3;
      const std::size_t moleculeB = col / 3;
      if (moleculeA == moleculeB)
      {
        continue;
      }

      const std::size_t axisA = row % 3;
      const std::size_t axisB = col % 3;

      double &coordA = (&atoms[moleculeA].position.x)[axisA];
      double &coordB = (&atoms[moleculeB].position.x)[axisB];

      coordA += delta;
      coordB += delta;
      const double ePP = interEnergy();
      coordB -= 2.0 * delta;
      const double ePM = interEnergy();
      coordA -= 2.0 * delta;
      const double eMM = interEnergy();
      coordB += 2.0 * delta;
      const double eMP = interEnergy();
      coordA += delta;
      coordB -= delta;

      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }
}
