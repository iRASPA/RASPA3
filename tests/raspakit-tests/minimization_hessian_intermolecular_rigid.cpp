#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import generalized_hessian;
import running_energy;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import interactions_hessian_intermolecular;
import interactions_intermolecular;

namespace
{
void setRigidMoleculePositions(System &system, std::size_t moleculeIndex)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  const Component &component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  const double3 com = molecule.centerOfMassPosition;
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        com + rotation * component.atoms[localAtom].position;
  }
}

void perturbRigidCom(System &system, std::size_t moleculeIndex, std::size_t axis, double delta)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  (&molecule.centerOfMassPosition.x)[axis] += delta;
  setRigidMoleculePositions(system, moleculeIndex);
}

simd_quatd quatFromRotationVector(const double3 &omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

struct RigidState
{
  double3 com;
  simd_quatd orientation;
};

// Chart matching the analytic rigid kinematics: orientation displacements premultiply
// the base orientation with the exponential map of the rotation vector.
void applyDisplacedState(System &system, std::size_t moleculeIndex, const RigidState &base,
                         std::span<const double> displacement)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  molecule.centerOfMassPosition = base.com + double3(displacement[0], displacement[1], displacement[2]);
  const double3 omega(displacement[3], displacement[4], displacement[5]);
  molecule.orientation = (quatFromRotationVector(omega) * base.orientation).normalized();
  setRigidMoleculePositions(system, moleculeIndex);
}

// Same chart with the isotropic exp(epsilon) strain applied: the center of mass scales with the
// cell while the internal geometry (offsets and orientation) is unaffected.
void applyDisplacedStrainedState(System &system, std::size_t moleculeIndex, const RigidState &base,
                                 std::span<const double> displacement, double strainFactor)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  molecule.centerOfMassPosition =
      strainFactor * (base.com + double3(displacement[0], displacement[1], displacement[2]));
  const double3 omega(displacement[3], displacement[4], displacement[5]);
  molecule.orientation = (quatFromRotationVector(omega) * base.orientation).normalized();
  setRigidMoleculePositions(system, moleculeIndex);
}

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

TEST(minimization_hessian_intermolecular, rigid_methane_pair_com_vdw_matches_finite_difference)
{
  const double delta = 1e-5;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component methane = Component::makeMethane(forceField, 0);

  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);

  system.moleculeData[0].centerOfMassPosition = double3(8.0, 10.0, 10.0);
  system.moleculeData[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  system.moleculeData[1].centerOfMassPosition = double3(9.4, 10.0, 10.0);
  system.moleculeData[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  setRigidMoleculePositions(system, 0);
  setRigidMoleculePositions(system, 1);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeInterMolecularHessian(system, layout, hessian);

  auto interEnergy = [&]()
  {
    return Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox,
                                                     system.spanOfMoleculeAtoms())
        .moleculeMoleculeVDW;
  };

  for (std::size_t moleculeA = 0; moleculeA < 2; ++moleculeA)
  {
    for (std::size_t moleculeB = 0; moleculeB < 2; ++moleculeB)
    {
      if (moleculeA == moleculeB)
      {
        continue;
      }

      for (std::size_t axisA = 0; axisA < 3; ++axisA)
      {
        for (std::size_t axisB = 0; axisB < 3; ++axisB)
        {
          const auto dofA = layout.rigidMoleculeDof(moleculeA, static_cast<RigidDof>(axisA));
          const auto dofB = layout.rigidMoleculeDof(moleculeB, static_cast<RigidDof>(axisB));
          ASSERT_TRUE(dofA);
          ASSERT_TRUE(dofB);

          perturbRigidCom(system, moleculeA, axisA, delta);
          perturbRigidCom(system, moleculeB, axisB, delta);
          const double ePP = interEnergy();
          perturbRigidCom(system, moleculeB, axisB, -2.0 * delta);
          const double ePM = interEnergy();
          perturbRigidCom(system, moleculeA, axisA, -2.0 * delta);
          const double eMM = interEnergy();
          perturbRigidCom(system, moleculeB, axisB, 2.0 * delta);
          const double eMP = interEnergy();
          perturbRigidCom(system, moleculeA, axisA, delta);
          perturbRigidCom(system, moleculeB, axisB, -delta);

          const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
          const double analytic = hessian(*dofA, *dofB);
          const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
          EXPECT_NEAR(numerical, analytic, relativeTolerance * scale)
              << "molA=" << moleculeA << " molB=" << moleculeB << " axisA=" << axisA << " axisB=" << axisB;
        }
      }
    }
  }
}

TEST(minimization_hessian_intermolecular, rigid_co2_pair_com_orientation_vdw_matches_finite_difference)
{
  // delta balances roundoff (energy cancellation over 4*delta^2) against O(delta^2) truncation.
  const double delta = 1e-4;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component co2 = Component::makeCO2(forceField, 0, false);

  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);

  const std::array<RigidState, 2> baseStates = {
      RigidState{double3(8.0, 10.0, 10.0),
                 simd_quatd::fromAxisAngle(0.3, double3(1.0, 2.0, 3.0) / std::sqrt(14.0)).normalized()},
      RigidState{double3(11.5, 10.7, 10.4),
                 simd_quatd::fromAxisAngle(-0.7, double3(0.5, -1.0, 2.0) / std::sqrt(5.25)).normalized()}};

  for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
  {
    system.moleculeData[moleculeIndex].centerOfMassPosition = baseStates[moleculeIndex].com;
    system.moleculeData[moleculeIndex].orientation = baseStates[moleculeIndex].orientation;
    setRigidMoleculePositions(system, moleculeIndex);
  }

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeInterMolecularHessian(system, layout, hessian);

  // Global DOF index -> (molecule, local index within its 6-vector displacement).
  std::array<std::pair<std::size_t, std::size_t>, 12> dofMap{};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
  {
    for (std::size_t local = 0; local < 6; ++local)
    {
      const auto dof = layout.rigidMoleculeDof(moleculeIndex, static_cast<RigidDof>(local));
      ASSERT_TRUE(dof);
      dofMap[*dof] = {moleculeIndex, local};
    }
  }

  auto energyAtDisplacement = [&](std::span<const double> displacement)
  {
    for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
    {
      applyDisplacedState(system, moleculeIndex, baseStates[moleculeIndex],
                          displacement.subspan(moleculeIndex * 6, 6));
    }
    const RunningEnergy energy = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox,
                                                                           system.spanOfMoleculeAtoms());
    return energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge;
  };

  std::array<double, 12> displacement{};
  for (std::size_t row = 0; row < 12; ++row)
  {
    for (std::size_t col = 0; col < 12; ++col)
    {
      const std::size_t rowSlot = dofMap[row].first * 6 + dofMap[row].second;
      const std::size_t colSlot = dofMap[col].first * 6 + dofMap[col].second;

      displacement.fill(0.0);
      displacement[rowSlot] += delta;
      displacement[colSlot] += delta;
      const double ePP = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[rowSlot] += delta;
      displacement[colSlot] -= delta;
      const double ePM = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[rowSlot] -= delta;
      displacement[colSlot] += delta;
      const double eMP = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[rowSlot] -= delta;
      displacement[colSlot] -= delta;
      const double eMM = energyAtDisplacement(displacement);

      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }

  // Restore base state.
  displacement.fill(0.0);
  energyAtDisplacement(displacement);
}

TEST(minimization_hessian_intermolecular, rigid_co2_pair_strain_matches_finite_difference)
{
  const double delta = 1e-4;
  const double strainStep = 2e-4;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component co2 = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);

  const std::array<RigidState, 2> baseStates = {
      RigidState{double3(12.0, 15.0, 15.0),
                 simd_quatd::fromAxisAngle(0.3, double3(1.0, 2.0, 3.0) / std::sqrt(14.0)).normalized()},
      RigidState{double3(15.5, 15.7, 15.4),
                 simd_quatd::fromAxisAngle(-0.7, double3(0.5, -1.0, 2.0) / std::sqrt(5.25)).normalized()}};

  std::array<double, 12> displacement{};
  auto applyState = [&](double strainExponent)
  {
    const double factor = std::exp(strainExponent);
    for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
    {
      applyDisplacedStrainedState(system, moleculeIndex, baseStates[moleculeIndex],
                                  std::span<const double>(displacement).subspan(moleculeIndex * 6, 6), factor);
    }
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  Interactions::computeInterMolecularHessian(system, layout, hessian);

  std::array<std::size_t, 12> dofToSlot{};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
  {
    for (std::size_t local = 0; local < 6; ++local)
    {
      const auto dof = layout.rigidMoleculeDof(moleculeIndex, static_cast<RigidDof>(local));
      ASSERT_TRUE(dof);
      dofToSlot[*dof] = moleculeIndex * 6 + local;
    }
  }

  const double3 boxLengths(30.0, 30.0, 30.0);
  auto energyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * boxLengths.x, factor * boxLengths.y, factor * boxLengths.z);
    const RunningEnergy energy =
        Interactions::computeInterMolecularEnergy(system.forceField, strainedBox, system.spanOfMoleculeAtoms());
    return energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge;
  };

  // strain-strain block.
  {
    displacement.fill(0.0);
    auto finiteDifference = [&](double step)
    {
      const double eZero = energyAt(0.0);
      const double ePlus = energyAt(step);
      const double eMinus = energyAt(-step);
      return (ePlus - 2.0 * eZero + eMinus) / (step * step);
    };
    const double numerical = (4.0 * finiteDifference(strainStep) - finiteDifference(2.0 * strainStep)) / 3.0;
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks over all rigid DOFs.
  for (std::size_t dof = 0; dof < 12; ++dof)
  {
    const std::size_t slot = dofToSlot[dof];

    auto finiteDifference = [&](double positionStep, double strainDelta)
    {
      displacement.fill(0.0);
      displacement[slot] += positionStep;
      const double ePP = energyAt(strainDelta);
      const double ePM = energyAt(-strainDelta);
      displacement.fill(0.0);
      displacement[slot] -= positionStep;
      const double eMP = energyAt(strainDelta);
      const double eMM = energyAt(-strainDelta);
      return (ePP - ePM - eMP + eMM) / (4.0 * positionStep * strainDelta);
    };
    const double numerical =
        (4.0 * finiteDifference(delta, strainStep) - finiteDifference(2.0 * delta, 2.0 * strainStep)) / 3.0;
    const double analytic = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "position-strain dof=" << dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
}

TEST(minimization_hessian_intermolecular, rigid_co2_generalized_gradient_matches_finite_difference)
{
  constexpr double delta = 2e-6;
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  forceField.useCharge = false;
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system =
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);

  const std::array<RigidState, 2> baseStates = {
      RigidState{double3(12.0, 15.0, 15.0),
                 simd_quatd::fromAxisAngle(0.3, double3(1.0, 2.0, 3.0) / std::sqrt(14.0)).normalized()},
      RigidState{double3(15.5, 15.7, 15.4),
                 simd_quatd::fromAxisAngle(-0.7, double3(0.5, -1.0, 2.0) / std::sqrt(5.25)).normalized()}};
  for (std::size_t molecule = 0; molecule < 2; ++molecule)
  {
    applyDisplacedState(system, molecule, baseStates[molecule], std::array<double, 6>{});
  }

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs());
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, derivatives);

  for (std::size_t molecule = 0; molecule < 2; ++molecule)
  {
    for (std::size_t slot = 0; slot < 6; ++slot)
    {
      std::array<double, 6> displacement{};
      displacement[slot] = delta;
      applyDisplacedState(system, molecule, baseStates[molecule], displacement);
      const double plus =
          Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox,
                                                    system.spanOfMoleculeAtoms())
              .moleculeMoleculeVDW;
      displacement[slot] = -delta;
      applyDisplacedState(system, molecule, baseStates[molecule], displacement);
      const double minus =
          Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox,
                                                    system.spanOfMoleculeAtoms())
              .moleculeMoleculeVDW;
      applyDisplacedState(system, molecule, baseStates[molecule], std::array<double, 6>{});

      const auto dof = layout.rigidMoleculeDof(molecule, static_cast<RigidDof>(slot));
      ASSERT_TRUE(dof);
      const double numerical = (plus - minus) / (2.0 * delta);
      EXPECT_NEAR(gradient[*dof], numerical, 1e-5 * std::max(1.0, std::abs(numerical)))
          << "molecule=" << molecule << " slot=" << slot;
    }
  }
}

TEST(minimization_hessian_intermolecular, mixed_flexible_rigid_co2_pair_strain_matches_finite_difference)
{
  const double delta = 2e-4;
  const double strainStep = 2e-4;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component co2Rigid = Component::makeCO2(forceField, 0, true);
  Component co2Flexible = Component::makeCO2(forceField, 1, true);
  co2Flexible.rigid = false;

  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
                         {co2Rigid, co2Flexible}, {}, {1, 1}, 5);

  // Molecule 0 is rigid, molecule 1 is flexible (slightly bent to avoid degenerate geometry).
  const RigidState rigidBase{double3(12.0, 15.0, 15.0),
                             simd_quatd::fromAxisAngle(0.3, double3(1.0, 2.0, 3.0) / std::sqrt(14.0)).normalized()};
  const std::array<double3, 3> flexibleBase = {double3(15.5, 15.9, 15.2), double3(16.4, 15.6, 15.0),
                                               double3(17.3, 16.0, 14.7)};

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::size_t flexibleAtomBase = system.moleculeData[1].atomIndex;

  // Slots 0..5: rigid COM and orientation; slots 6..14: flexible atom displacements.
  std::array<double, 15> displacement{};
  auto applyState = [&](double strainExponent)
  {
    const double factor = std::exp(strainExponent);
    applyDisplacedStrainedState(system, 0, rigidBase, std::span<const double>(displacement).subspan(0, 6), factor);
    for (std::size_t localAtom = 0; localAtom < 3; ++localAtom)
    {
      const double3 offset(displacement[6 + localAtom * 3 + 0], displacement[6 + localAtom * 3 + 1],
                           displacement[6 + localAtom * 3 + 2]);
      atoms[flexibleAtomBase + localAtom].position = factor * (flexibleBase[localAtom] + offset);
    }
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 15u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  Interactions::computeInterMolecularHessian(system, layout, hessian);

  std::array<std::size_t, 15> dofToSlot{};
  for (std::size_t local = 0; local < 6; ++local)
  {
    const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(local));
    ASSERT_TRUE(dof);
    dofToSlot[*dof] = local;
  }
  for (std::size_t localAtom = 0; localAtom < 3; ++localAtom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const auto dof = layout.flexibleAtomDof(1, localAtom, static_cast<MinimizationDofAxis>(axis));
      ASSERT_TRUE(dof);
      dofToSlot[*dof] = 6 + localAtom * 3 + axis;
    }
  }

  const double3 boxLengths(30.0, 30.0, 30.0);
  auto energyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * boxLengths.x, factor * boxLengths.y, factor * boxLengths.z);
    const RunningEnergy energy =
        Interactions::computeInterMolecularEnergy(system.forceField, strainedBox, system.spanOfMoleculeAtoms());
    return energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge;
  };

  // Position-position block (covers the mixed flexible-rigid pair path).
  for (std::size_t row = 0; row < 15; ++row)
  {
    for (std::size_t col = 0; col < 15; ++col)
    {
      const std::size_t rowSlot = dofToSlot[row];
      const std::size_t colSlot = dofToSlot[col];

      auto finiteDifference = [&](double step)
      {
        displacement.fill(0.0);
        displacement[rowSlot] += step;
        displacement[colSlot] += step;
        const double ePP = energyAt(0.0);

        displacement.fill(0.0);
        displacement[rowSlot] += step;
        displacement[colSlot] -= step;
        const double ePM = energyAt(0.0);

        displacement.fill(0.0);
        displacement[rowSlot] -= step;
        displacement[colSlot] += step;
        const double eMP = energyAt(0.0);

        displacement.fill(0.0);
        displacement[rowSlot] -= step;
        displacement[colSlot] -= step;
        const double eMM = energyAt(0.0);
        return (ePP - ePM - eMP + eMM) / (4.0 * step * step);
      };
      const double numerical = (4.0 * finiteDifference(delta) - finiteDifference(2.0 * delta)) / 3.0;
      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }

  // strain-strain block.
  {
    displacement.fill(0.0);
    auto finiteDifference = [&](double step)
    {
      const double eZero = energyAt(0.0);
      const double ePlus = energyAt(step);
      const double eMinus = energyAt(-step);
      return (ePlus - 2.0 * eZero + eMinus) / (step * step);
    };
    const double numerical = (4.0 * finiteDifference(strainStep) - finiteDifference(2.0 * strainStep)) / 3.0;
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks over all DOFs.
  for (std::size_t dof = 0; dof < 15; ++dof)
  {
    const std::size_t slot = dofToSlot[dof];

    auto finiteDifference = [&](double positionStep, double strainDelta)
    {
      displacement.fill(0.0);
      displacement[slot] += positionStep;
      const double ePP = energyAt(strainDelta);
      const double ePM = energyAt(-strainDelta);
      displacement.fill(0.0);
      displacement[slot] -= positionStep;
      const double eMP = energyAt(strainDelta);
      const double eMM = energyAt(-strainDelta);
      return (ePP - ePM - eMP + eMM) / (4.0 * positionStep * strainDelta);
    };
    const double numerical =
        (4.0 * finiteDifference(delta, strainStep) - finiteDifference(2.0 * delta, 2.0 * strainStep)) / 3.0;
    const double analytic = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "position-strain dof=" << dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
}
