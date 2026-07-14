#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import system;
import simulationbox;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import interactions_hessian_framework_molecule;
import interactions_framework_molecule;

namespace
{
void setRigidMoleculePositions(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  const double3 com = molecule.centerOfMassPosition;
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        com + rotation * component.atoms[localAtom].position;
  }
}

simd_quatd quatFromRotationVector(const double3& omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

void useSecondOrderTaylorShiftedLennardJones(ForceField& forceField)
{
  for (VDWParameters& parameters : forceField.data)
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

TEST(minimization_hessian_framework, flexible_framework_molecule_cross_block_matches_finite_difference)
{
  ForceField forceField({{"A", false, 12.0, 0.0, 0.0, 6, false}}, {{120.0, 3.4}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, true, false, false);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const SimulationBox box(20.0, 20.0, 20.0);
  const Atom frameworkAtom({0.25, 0.25, 0.25}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "flexible-framework", box, 1, {frameworkAtom}, {frameworkAtom}, {1, 1, 1});
  framework.rigid = false;

  ConnectivityTable connectivity(1);
  Component component(forceField, "single-site", 100.0, 1.0e6, 0.1,
                      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, connectivity,
                      Potentials::IntraMolecularPotentials{}, 0, 0);
  component.rigid = false;
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {component}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfMoleculeAtoms()[0].position = {8.8, 5.2, 5.1};

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size());
  ASSERT_EQ(layout.numDofs(), 6u);
  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeFrameworkMoleculeHessian(system, layout, hessian, system.spanOfMoleculeDynamics(),
                                                system.spanOfFrameworkDynamics());

  const auto energy = [&]()
  {
    const RunningEnergy value = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
    return value.frameworkMoleculeVDW + value.frameworkMoleculeCharge;
  };
  constexpr double delta = 2.0e-4;
  constexpr double tolerance = 5.0e-2;
  std::array<Atom*, 2> atoms{&system.spanOfFrameworkAtoms()[0], &system.spanOfMoleculeAtoms()[0]};
  const double referenceEnergy = energy();
  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t column = 0; column < 6; ++column)
    {
      Atom& atomA = *atoms[row / 3];
      Atom& atomB = *atoms[column / 3];
      double numerical{};
      if (row == column)
      {
        (&atomA.position.x)[row % 3] += delta;
        const double plus = energy();
        (&atomA.position.x)[row % 3] -= 2.0 * delta;
        const double minus = energy();
        (&atomA.position.x)[row % 3] += delta;
        numerical = (plus - 2.0 * referenceEnergy + minus) / (delta * delta);
      }
      else
      {
        (&atomA.position.x)[row % 3] += delta;
        (&atomB.position.x)[column % 3] += delta;
        const double ePP = energy();
        (&atomB.position.x)[column % 3] -= 2.0 * delta;
        const double ePM = energy();
        (&atomA.position.x)[row % 3] -= 2.0 * delta;
        const double eMM = energy();
        (&atomB.position.x)[column % 3] += 2.0 * delta;
        const double eMP = energy();
        (&atomA.position.x)[row % 3] += delta;
        (&atomB.position.x)[column % 3] -= delta;
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }
      EXPECT_NEAR(hessian(row, column), numerical, tolerance) << "row=" << row << " column=" << column;
    }
  }
}

TEST(minimization_hessian_framework, rigid_co2_in_itq29_vdw_matches_finite_difference)
{
  const double delta = 4e-4;
  const double relativeTolerance = 1e-5;
  // Roundoff floor of the FD reference: the total framework energy is large, so the
  // four-point cancellation leaves noise of order eps * E / (4 delta^2).
  const double absoluteTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  // 2x2x2 unit cells so the box exceeds twice the cutoff (minimum image stays consistent).
  Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component co2 = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 baseCom(5.93355, 7.93355, 5.93355);
  const simd_quatd baseOrientation =
      simd_quatd::fromAxisAngle(0.4, double3(1.0, -2.0, 1.5) / std::sqrt(7.25)).normalized();

  system.moleculeData[0].centerOfMassPosition = baseCom;
  system.moleculeData[0].orientation = baseOrientation;
  setRigidMoleculePositions(system, 0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeFrameworkMoleculeHessian(system, layout, hessian);

  auto energyAtDisplacement = [&](std::span<const double> displacement)
  {
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = baseCom + double3(displacement[0], displacement[1], displacement[2]);
    const double3 omega(displacement[3], displacement[4], displacement[5]);
    molecule.orientation = (quatFromRotationVector(omega) * baseOrientation).normalized();
    setRigidMoleculePositions(system, 0);

    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  std::array<double, 6> displacement{};
  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t col = 0; col < 6; ++col)
    {
      const auto dofRow = layout.rigidMoleculeDof(0, static_cast<RigidDof>(row));
      const auto dofCol = layout.rigidMoleculeDof(0, static_cast<RigidDof>(col));
      ASSERT_TRUE(dofRow);
      ASSERT_TRUE(dofCol);

      auto finiteDifference = [&](double step)
      {
        displacement.fill(0.0);
        displacement[row] += step;
        displacement[col] += step;
        const double ePP = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] += step;
        displacement[col] -= step;
        const double ePM = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] -= step;
        displacement[col] += step;
        const double eMP = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] -= step;
        displacement[col] -= step;
        const double eMM = energyAtDisplacement(displacement);
        return (ePP - ePM - eMP + eMM) / (4.0 * step * step);
      };
      const double numerical = (4.0 * finiteDifference(delta) - finiteDifference(2.0 * delta)) / 3.0;
      const double analytic = hessian(*dofRow, *dofCol);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale)
          << "row=" << row << " col=" << col;
    }
  }
}

TEST(minimization_hessian_framework, rigid_co2_in_itq29_strain_matches_finite_difference)
{
  const double strainStep = 1.2e-3;
  const double relativeTolerance = 1e-5;
  // Roundoff floor: the total framework-molecule energy is large, so the finite-difference
  // cancellation leaves absolute noise.
  const double absoluteTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component co2 = Component::makeCO2(forceField, 0, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 baseCom(5.93355, 7.93355, 5.93355);
  const simd_quatd baseOrientation =
      simd_quatd::fromAxisAngle(0.4, double3(1.0, -2.0, 1.5) / std::sqrt(7.25)).normalized();

  // Base framework positions and cell; under strain both scale with exp(epsilon), matching the
  // affine convention of the analytic strain blocks.
  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<double3> baseFrameworkPositions(frameworkAtoms.size());
  for (std::size_t i = 0; i < frameworkAtoms.size(); ++i)
  {
    baseFrameworkPositions[i] = frameworkAtoms[i].position;
  }
  const SimulationBox baseBox = system.simulationBox;

  std::array<double, 6> displacement{};
  auto applyState = [&](double strainExponent)
  {
    const double factor = std::exp(strainExponent);
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = factor * (baseCom + double3(displacement[0], displacement[1], displacement[2]));
    const double3 omega(displacement[3], displacement[4], displacement[5]);
    molecule.orientation = (quatFromRotationVector(omega) * baseOrientation).normalized();
    setRigidMoleculePositions(system, 0);
    for (std::size_t i = 0; i < frameworkAtoms.size(); ++i)
    {
      frameworkAtoms[i].position = factor * baseFrameworkPositions[i];
    }
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  const RunningEnergy analyticEnergy = Interactions::computeFrameworkMoleculeHessian(system, layout, hessian);

  auto energyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * baseBox.cell.ax, factor * baseBox.cell.by, factor * baseBox.cell.cz);
    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, strainedBox, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  // Consistency: the Hessian routine and the FD evaluator must see the same base state.
  EXPECT_NEAR(analyticEnergy.frameworkMoleculeVDW + analyticEnergy.frameworkMoleculeCharge, energyAt(0.0), 1e-8)
      << "base-state energy mismatch between analytic Hessian and FD evaluator";

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
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks over the rigid COM and orientation DOFs.
  auto gradientAt = [&](double strainExponent)
  {
    displacement.fill(0.0);
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    system.simulationBox = SimulationBox(factor * baseBox.cell.ax, factor * baseBox.cell.by, factor * baseBox.cell.cz);

    std::vector<double> gradient(layout.numDofs());
    GeneralizedHessian scratchHessian(layout.numDofs(), 0);
    DerivativeCapabilities capabilities{.energy = false, .gradient = true, .hessianPositionPosition = true};
    DerivativeResults derivatives{.gradient = gradient, .hessian = scratchHessian};
    evaluateDerivatives(system, layout, capabilities, derivatives);
    for (std::size_t slot = 0; slot < 3; ++slot)
    {
      gradient[slot] *= factor;
    }
    return gradient;
  };
  constexpr double gradientStep = 5e-5;
  const std::vector<double> gradientPlus = gradientAt(gradientStep);
  const std::vector<double> gradientMinus = gradientAt(-gradientStep);
  for (std::size_t slot = 0; slot < 6; ++slot)
  {
    const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(slot));
    ASSERT_TRUE(dof);

    const double numerical = (gradientPlus[*dof] - gradientMinus[*dof]) / (2.0 * gradientStep);
    const double analytic = hessian.positionStrain()[*dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale) << "position-strain dof=" << *dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
  system.simulationBox = baseBox;
}

TEST(minimization_hessian_framework, single_framework_atom_rigid_co2_strain_step_convergence)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const Atom frameworkSite(double3(0.25, 0.25, 0.25), 0.0, 1.0, 0, 3, 0, false, false);
  Framework framework(forceField, "single-site", SimulationBox(40.0, 40.0, 40.0), 1, {frameworkSite}, {frameworkSite},
                      int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 frameworkBase(10.0, 10.0, 10.0);
  const double3 comBase(14.0, 12.0, 11.0);
  const simd_quatd orientationBase =
      simd_quatd::fromAxisAngle(0.6, double3(1.0, -2.0, 0.7) / std::sqrt(5.49)).normalized();
  system.spanOfFrameworkAtoms()[0].position = frameworkBase;

  std::array<double, 6> displacement{};
  auto applyState = [&](double epsilon)
  {
    const double scale = std::exp(epsilon);
    system.spanOfFrameworkAtoms()[0].position = scale * frameworkBase;
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = scale * (comBase + double3(displacement[0], displacement[1], displacement[2]));
    molecule.orientation =
        (quatFromRotationVector(double3(displacement[3], displacement[4], displacement[5])) * orientationBase)
            .normalized();
    setRigidMoleculePositions(system, 0);
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);
  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  Interactions::computeFrameworkMoleculeHessian(system, layout, hessian, dynamics);

  constexpr double gradientStep = 1e-5;
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t atom = 0; atom < moleculeAtoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double& coordinate = (&moleculeAtoms[atom].position.x)[axis];
      coordinate += gradientStep;
      const RunningEnergy plus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate -= 2.0 * gradientStep;
      const RunningEnergy minus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate += gradientStep;
      const double numerical = ((plus.frameworkMoleculeVDW + plus.frameworkMoleculeCharge) -
                                (minus.frameworkMoleculeVDW + minus.frameworkMoleculeCharge)) /
                               (2.0 * gradientStep);
      EXPECT_NEAR((&dynamics[atom].gradient.x)[axis], numerical, 1e-5 * std::max(1.0, std::abs(numerical)));
    }
  }

  auto energyAt = [&](double epsilon)
  {
    applyState(epsilon);
    const double scale = std::exp(epsilon);
    const SimulationBox box(40.0 * scale, 40.0 * scale, 40.0 * scale);
    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, box, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  // A step sweep distinguishes a bad analytic expression from FD truncation/roundoff.
  // Each estimate must converge to the same analytic value over the stable central-difference range.
  constexpr std::array<double, 4> strainSteps = {5e-4, 4e-4, 3e-4, 2e-4};
  constexpr double positionStep = 2e-5;
  for (double strainStep : strainSteps)
  {
    displacement.fill(0.0);
    const double numericalStrain =
        (energyAt(strainStep) - 2.0 * energyAt(0.0) + energyAt(-strainStep)) / (strainStep * strainStep);
    EXPECT_NEAR(numericalStrain, hessian.strainStrain()[0],
                1e-5 * std::max({1.0, std::abs(numericalStrain), std::abs(hessian.strainStrain()[0])}))
        << "strainStep=" << strainStep;

    for (std::size_t slot = 0; slot < 6; ++slot)
    {
      const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(slot));
      ASSERT_TRUE(dof);
      displacement.fill(0.0);
      displacement[slot] = positionStep;
      const double ePP = energyAt(strainStep);
      const double ePM = energyAt(-strainStep);
      displacement[slot] = -positionStep;
      const double eMP = energyAt(strainStep);
      const double eMM = energyAt(-strainStep);
      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * positionStep * strainStep);
      const double analytic = hessian.positionStrain()[*dof];
      EXPECT_NEAR(numerical, analytic, 1e-5 * std::max({1.0, std::abs(numerical), std::abs(analytic)}))
          << "slot=" << slot << " strainStep=" << strainStep;
    }
  }
}
